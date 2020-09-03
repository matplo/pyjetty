#!/usr/bin/env python3

from pyjetty.mputils import MPBase, pwarning, pinfo, perror, treewriter, jet_analysis
import random
import uproot
import pandas as pd
import numpy as np
import fastjet as fj
import fjext
import os
import tqdm
import argparse


def unique_fname(fn):
	counter = 0
	outfn = fn.replace('.root', '_{}.root'.format(counter))
	while os.path.exists(outfn):
		counter += 1
		outfn = fn.replace('.root', '_{}.root'.format(counter))
	return outfn


class HFAIO(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(d0_tree_name='PWGHF_TreeCreator/tree_D0', 
								 track_tree_name='PWGHF_TreeCreator/tree_Particle',
								 event_tree_name='PWGHF_TreeCreator/tree_event_char',
								 enable_jet=True,
								 enable_d0=True,
								 offset_parts=0,
								 output_prefix='./HFAIO',
								 input_file = None)
		super(HFAIO, self).__init__(**kwargs)
		self.analyses = []
		self.df_grouped = None
		self.df_events = None

		# temp output
		out_file_1 = unique_fname(self.output_prefix + '_djet_tout.root')
		self.tw = treewriter.RTreeWriter(name = 'd0j', file_name = out_file_1)
		out_file_2 = unique_fname(self.output_prefix + '_djet_correl_tout.root')
		self.twjc = treewriter.RTreeWriter(name = 'd0jc', file_name = out_file_2)

		if self.input_file:
			self.process_file(self.input_file)

	def __del__(self):
		self.tw.write_and_close()
		self.twjc.write_and_close()

	def pd_tree(self, path, tname, squery=None):
		try:
			tree = uproot.open(path)[tname]
		except:
			pwarning('error getting', tname, 'from file:', path)
			return None
		if not tree:
			perror('Tree {} not found in file {}'.format(tname, path))
			return None
		df = tree.pandas.df()
		if squery:
			#df.query(squery, inplace=True)
			df = df.query(squery)
			df.reset_index(drop=True)
		return df

	def process_file(self, fname):
		_ev_cuts = "is_ev_rej == 0 & abs(z_vtx_reco) < 10."
		self.event_df = self.pd_tree(path=fname, tname=self.event_tree_name, squery=_ev_cuts)
		if self.event_df is None:
			return False
		pinfo('events from', fname, len(self.event_df.index))

		_d0cuts_base = "(pt_cand > 4.0 & pt_prong0 > 0.5 & pt_prong1 > 0.5 & abs(eta_cand) < 0.8) & "
		_d0cuts_kpi = _d0cuts_base
		_d0cuts_kpi += "((abs(nsigTPC_Pi_0) < 3. & (abs(nsigTOF_Pi_0) < 3. | nsigTOF_Pi_0 < -900) & abs(nsigTPC_K_1) < 3. & (abs(nsigTOF_K_1) < 3. | nsigTOF_K_0 < -900)) | "
		_d0cuts_kpi += "(abs(nsigTPC_Pi_1) < 3. & (abs(nsigTOF_Pi_1) < 3. | nsigTOF_Pi_1 < -900) & abs(nsigTPC_K_0) < 3. & (abs(nsigTOF_K_0) < 3. | nsigTOF_K_0 < -900)))"

		self.d0_df = self.pd_tree(path=fname, tname=self.d0_tree_name, squery=_d0cuts_kpi)
		if self.d0_df is None:
			return False
		pinfo('d0s from', fname, len(self.d0_df.index))
		# pinfo(list(self.event_df))
		if 'ev_id_ext' in list(self.event_df):
			self.d0ev_df = pd.merge(self.d0_df, self.event_df, on=['run_number', 'ev_id', 'ev_id_ext'])
		else:
			self.d0ev_df = pd.merge(self.d0_df, self.event_df, on=['run_number', 'ev_id'])
		self.d0ev_df.query(_ev_cuts, inplace=True)
		self.d0ev_df_grouped = self.d0ev_df.groupby(['run_number','ev_id'])
		pinfo('d0s after event cuts from ', fname, len(self.d0ev_df.index))
		pinfo('N d0 groups after event cuts from ', fname, len(self.d0ev_df_grouped))

		self.track_df = self.pd_tree(path=fname, tname=self.track_tree_name)
		# self.track_df = _track_df.groupby(['run_number','ev_id'])
		if self.track_df is None:
			return False
		pinfo('tracks from', fname, len(self.track_df.index))

		# event based processing - not efficient for D0 analysis
		# self.pbar.close()
		# 	self.event_df.apply(self.process_event, axis=1)
		# with tqdm.tqdm(total=len(self.event_df.index)) as self.pbar:

		# d0 based processing
		# with tqdm.tqdm(total=len(self.d0ev_df.index)) as self.pbar:
		#	self.d0ev_df.apply(self.process_d0s, axis=1)
		with tqdm.tqdm(total=len(self.d0ev_df_grouped)) as self.pbar:
			_tmp = self.d0ev_df_grouped.apply(self.process_d0s)
		self.pbar.close()

		self.event_df = None
		self.d0_df = None
		self.d0ev_df = None
		self.d0ev_df_grouped = None
		self.track_df = None

	def d0_jet_correl(self, jets, _d0s, _d0_imass_list):
		_filled = False
		for j in jets:
			jc = j.constituents()
			for ic in range(len(jc)):
				c = jc[ic]
				if c.user_index() >= self._user_index_offset:
					_d0_index = c.user_index() - self._user_index_offset 
					_d0 = _d0s[_d0_index]
					self.twjc.fill_branches(jet = j, d0 = _d0, dR = j.delta_R(_d0), minv = _d0_imass_list[_d0_index], m = _d0.m(), z = _d0.perp() / j.perp())
					_filled = True
		if _filled:
			self.twjc.fill_tree()

	def process_d0s(self, df):
		self.pbar.update(1)
		_n_d0s = len(df)
		if _n_d0s < 1:
			return
		# pinfo(df)
		if 'ev_id_ext' in list(self.event_df):
			# _ev_query = "run_number == {} & ev_id == {} & ev_id_ext == {}".format(df['run_number'], df['ev_id'], df['ev_id_ext'])
			_ev_query = "run_number == {} & ev_id == {} & ev_id_ext == {}".format(df['run_number'].values[0], df['ev_id'].values[0], df['ev_id_ext'].values[0])
		else:
			_ev_query = "run_number == {} & ev_id == {}".format(df['run_number'].values[0], df['ev_id'].values[0])
		_df_tracks = self.track_df.query(_ev_query)
		_df_tracks.reset_index(drop=True)
		_parts = fjext.vectorize_pt_eta_phi(_df_tracks['ParticlePt'].values, _df_tracks['ParticleEta'].values, _df_tracks['ParticlePhi'].values)
		self._user_index_offset = 10000
		# _d0s = fjext.vectorize_pt_eta_phi([df['pt_cand']], [df['eta_cand']], [df['phi_cand']], self._user_index_offset)
		_d0s = fjext.vectorize_pt_eta_phi(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values, self._user_index_offset)
		_d0s_gh = [p * 1.e-6 for p in _d0s]

		_parts_and_ds = _parts
		_tmp = [_parts_and_ds.push_back(p) for p in _d0s_gh]
		# pinfo('n parts = ', len(_parts_and_ds))

		ja = jet_analysis.JetAnalysis(jet_R = 0.4, particle_eta_max=0.9, jet_pt_min=10.0)
		ja.analyze_event(_parts_and_ds)

		_d0_imass_list = df['inv_mass'].values.tolist()
		# _d0_imass_list = [df['inv_mass']]
		self.tw.fill_branches(dpsj = _d0s)
		self.tw.fill_branches(dpsjgh = _d0s_gh)
		self.tw.fill_branches(minv = _d0_imass_list)
		self.tw.fill_branches(jets = ja.jets_as_psj_vector())
		self.tw.fill_tree()

		self.d0_jet_correl(ja.jets, _d0s, _d0_imass_list)

		return True

	def process_event(self, df):
		self.pbar.update(1)
		if 'ev_id_ext' in list(self.event_df):
			_ev_query = "run_number == {} & ev_id == {} & ev_id_ext == {}".format(df['run_number'], df['ev_id'], df['ev_id_ext'])
		else:
			_ev_query = "run_number == {} & ev_id == {}".format(df['run_number'], df['ev_id'])
		_df_d0 = self.d0_df.query(_ev_query)
		_df_d0.reset_index(drop=True)
		_n_d0s = len(_df_d0.index)
		if _n_d0s < 1:
			return
		_ev_query = ""
		_df_tracks = self.track_df.query(_ev_query)
		_df_tracks.reset_index(drop=True)
		_parts = fjext.vectorize_pt_eta_phi(_df_tracks['ParticlePt'].values, _df_tracks['ParticleEta'].values, _df_tracks['ParticlePhi'].values)
		self._user_index_offset = 10000
		_d0s = fjext.vectorize_pt_eta_phi(_df_d0['pt_cand'].values, _df_d0['eta_cand'].values, _df_d0['phi_cand'].values, self._user_index_offset)
		_d0s_gh = [p * 1.e-6 for p in _d0s]

		_parts_and_ds = _parts
		_tmp = [_parts_and_ds.push_back(p) for p in _d0s_gh]
		# pinfo('n parts = ', len(_parts_and_ds))

		ja = jet_analysis.JetAnalysis(jet_R = 0.6, particle_eta_max=0.9, jet_pt_min=5.0)
		ja.analyze_event(_parts_and_ds)

		_d0_imass_list = _df_d0['inv_mass'].values.tolist()
		self.tw.fill_branches(dpsj = _d0s)
		self.tw.fill_branches(dpsjgh = _d0s_gh)
		self.tw.fill_branches(minv = _d0_imass_list)
		self.tw.fill_branches(jets = ja.jets_as_psj_vector())
		self.tw.fill_tree()

		self.d0_jet_correl(ja.jets, _d0s, _d0_imass_list)

def process_files(fname):
	pinfo('reading file list from', fname)
	with open(fname) as f:
		flist = f.readlines()
	pinfo('number of files', len(flist))
	for ifn, fn in enumerate(flist):
		pinfo('file', ifn, 'of', len(flist))
		HFAIO(output_file='./hfaio_rfile_{}'.format(ifn), input_file=fn.strip('\n'))


def main():
	parser = argparse.ArgumentParser(description='D0 analysis on alice data', prog=os.path.basename(__file__))
	parser.add_argument('-f', '--flist', help='single root file or a file with a list of files to process', type=str, default=None, required=True)
	parser.add_argument('-n', '--nfiles', help='max n files to process', type=int, default=0, required=False)
	parser.add_argument('-o', '--output', help="prefix output file names", type=str, default='./hfaio_rfile')
	args = parser.parse_args()

	if '.root' in args.flist:
		HFAIO(output_prefix=args.output, input_file=args.flist)
	else:
		process_files(args.flist)		

if __name__ == '__main__':
	main()