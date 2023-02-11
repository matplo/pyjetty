#!/usr/bin/env python

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import fastjet as fj
import fjext
import fjcontrib
import fjtools

import pythia8
import pythiafjext
import pythiaext
from heppy.pythiautils import configuration as pyconf

# from tqdm.notebook import tqdm
from tqdm import tqdm
import argparse
import os
import sys
from pyjetty.mputils.generic_object import GenericObject
from pyjetty.mputils.csubtractor import CEventSubtractor
from pyjetty.mputils.jpickle import JetPickleIO

# from pyjetty.mputils.data_io import DataBackgroundIO
try:
	import uproot3 as uproot
	uproot_version = 3
except:
	try:
		import uproot
		uproot_version = 4 # assumed
	except:
		pass
import pandas as pd


def get_args_from_settings(ssettings):
    sys.argv=[' '] + ssettings.split()
    parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly')
    pyconf.add_standard_pythia_args(parser)
    parser.add_argument('--output', default="test_ang_ue.root", type=str)
    parser.add_argument('--user-seed', help='pythia seed', default=1111, type=int)
    args = parser.parse_args()
    return args


def analyze_event(bg_parts, jet, args):
	mpt = fjtools.matched_pt(j1, j2)
	jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
	jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
	if self.debug_level > 2:
		print('WTA jet definition is:', jet_def_wta)
	reclusterer_wta =  fjcontrib.Recluster(jet_def_wta)

	sd = fjcontrib.SoftDrop(0, args.sd_zcut, args.jet_R0)

	jet_wta = reclusterer_wta.result(jet)
	jet_groomed = 0#jet_groomed_lund.pair()


class DataEvent(GenericObject):
	def __init__(self, **kwargs):
		super(DataEvent, self).__init__(**kwargs)


class DataReader(GenericObject):
	_default_track_tree_name = 'PWGHF_TreeCreator/tree_Particle'
	_default_event_tree_name = 'PWGHF_TreeCreator/tree_event_char'
	_default_event_cols_select = ['run_number', 'ev_id', 'z_vtx_reco', 'is_ev_rej', 'perc_v0m', 'n_tracks']
	_default_user_index_offset = 10000
	_defaults = {
		'track_tree_name' : _default_track_tree_name,
		'event_tree_name' : _default_event_tree_name,
		'default_event_cols_select' : _default_event_cols_select,
		'user_index_offset' : _default_user_index_offset
	}

	def set_defaults(self):
		for d in DataReader._defaults:
			if self.__getattr__(d) is None:
				self.__setattr__(d, DataReader._defaults[d])
    
	def __init__(self, **kwargs):
		self.set_select_columns_pbpb()
		super(DataReader, self).__init__(**kwargs)
		self.set_defaults()
		if self.args:
			if type(self.args) == dict:
				self.configure_from_dict(self.args, ignore_none=True)
			else:
				self.configure_from_dict(self.args.__dict__)
		self.read_file_list()
		self.current_file_number = 0
		self.current_event_in_file = 0
		self.current_event = 0
  
	def read_file_list(self):
		self.list_of_files = []
		if self.file_list:
			with open(self.file_list) as f:
				self.list_of_files = [fn.strip() for fn in f.readlines()]

	def set_select_columns_pbpb(self, selection=None):
		if selection is not None:
			self.selected_event_columns = selection
		else:
			self.selected_event_columns=['run_number', 'ev_id', 'z_vtx_reco', 'is_ev_rej', 'perc_v0m', 'n_tracks']

	def get_pandas(self, fname, treename, selected_columns = None):
		if uproot_version < 4:
			_df = uproot.open(fname)[treename].pandas.df(selected_columns)
		else:
			_df = uproot.open(fname)[treename].arrays(selected_columns, library='pd')
		_df.reset_index(drop=True)
		if 'is_ev_rej' in _df.columns:
			_df_ok = _df.query('is_ev_rej == 0')
			_df_ok.reset_index(drop=True)
			return _df_ok
		return _df

	def build_events_obj(self, df):
		event = DataEvent()
		event.particles = fjext.vectorize_pt_eta_phi(df['ParticlePt'].values, df['ParticleEta'].values, df['ParticlePhi'].values,
                                               		self.user_index_offset)
		for c in df.columns:
			if 'Particle' not in c:
				event.__setattr__(c, df[c].values[0])
		return True

	def build_events(self, df):	
		event = {}
		event['particles'] = fjext.vectorize_pt_eta_phi(df['ParticlePt'].values, df['ParticleEta'].values, df['ParticlePhi'].values, 
                                                  		self.user_index_offset)
		for c in df.columns:
			if 'Particle' not in c:
				event[c] = df[c].values[0]
		self.events.append(event)
		return True
		
	def next_file(self):
		if self.current_file_number >= len(self.list_of_files):
			return False
		_fname = self.list_of_files[self.current_file_number]
		print('[i] loading file', _fname, file=sys.stderr)

		self.event_df = self.get_pandas(_fname, self.event_tree_name, self.default_event_cols_select)
		# print(f' - event_df \n{self.event_df}')
		self.track_df = self.get_pandas(_fname, self.track_tree_name)
		# print(f' - track_df \n{self.track_df}')
		self.merged_df = pd.merge(self.track_df, self.event_df, on=['run_number', 'ev_id'])
		# print(f' - merged_df \n{self.merged_df}')
		self.grouped_df = self.track_df.groupby(['run_number','ev_id'])
		# print(f' - grouped_df \n{self.grouped_df}')

		self.events = []
		self.events_df = self.grouped_df.apply(self.build_events)
		# print(f' - events_df', type(self.events_df), f'\n{self.events_df}')
		# print(f' - events_df \n{self.events_df.describe()}')

		if len(self.events_df) < 1:
			return self.next_file()
		if len(self.grouped_df) < 1:
			return self.next_file()

		# self.events = []
		# for e in self.events_df:
		# 	self.events.append(e)
		self.n_event_in_file = len(self.events)
		self.current_event_in_file = 0
		self.current_file_number += 1
    
		return True

	def next_event(self):
		if self.events_df is None:
			if self.next_file() is False:
				return False
		if self.current_event_in_file >= self.n_event_in_file:
			if self.next_file():
				return self.next_event()
			return False
		# print (f'event in file {self.current_event_in_file}', file=sys.stderr)
		# _rv = self.events_df[self.current_event_in_file]
		self.current_event_in_file += 1
		self.current_event += 1	
		return self.events[self.current_event_in_file-1]


def print_jet(j):
	print(f'[i] jet {j.pt()} {j.phi()} {j.eta()}')
	for c in fj.sorted_by_pt(j.constituents()):
		print(f'    {c.pt()} {c.phi()} {c.eta()} {c.user_index()}')
 

class Analysis(GenericObject):
	_default_file_output = __file__.replace('.py', '_output.root')
	_defaults = {
		'output' 			: _default_file_output,
		'jet_R0' 			: 0.2
	}

	def set_defaults(self):
		for d in DataReader._defaults:
			if self.__getattr__(d) is None:
				self.__setattr__(d, DataReader._defaults[d])

	def __init__(self, **kwargs):
		super(Analysis, self).__init__(**kwargs)
		self.set_defaults()
		if self.args:
			if type(self.args) == dict:
				self.configure_from_dict(self.args, ignore_none=True)
			else:
				self.configure_from_dict(self.args.__dict__)

		self.max_eta_hadron=0.9

		self.embedding = self.datalist
 
		if self.embedding:
			self.data 	= DataReader(file_list=self.datalist)

		self.cs 	= CEventSubtractor(alpha=self.alpha, 
                              			max_distance=self.dRmax, 
                                 		max_eta=self.max_eta_hadron, 
										bge_rho_grid_size=0.25, max_pt_correct=100)

		self.parts_selector_h = fj.SelectorAbsEtaMax(self.max_eta_hadron)
		self.jet_selector = fj.SelectorPtMin(self.jet_pt_min) & fj.SelectorPtMax(self.jet_pt_max) & fj.SelectorAbsEtaMax(self.max_eta_hadron - 1.05 * self.jet_R0)

		# set up our jet definition and a jet selector
		self.jet_def 		= fj.JetDefinition(fj.antikt_algorithm, self.jet_R0)
		self.jet_def_CAR0 	= fj.JetDefinition(fj.cambridge_aachen_algorithm, self.jet_R0)
		self.jet_def_CARinf = fj.JetDefinition(fj.cambridge_aachen_algorithm, self.jet_R0 * 10.)

		self.jet_def_CAR0_WTA 	= fj.JetDefinition(fj.cambridge_aachen_algorithm, self.jet_R0)
		self.jet_def_CAR0_WTA.set_recombination_scheme(fj.WTA_pt_scheme)
		self.jet_def_CARinf_WTA = fj.JetDefinition(fj.cambridge_aachen_algorithm, self.jet_R0 * 10.)
		self.jet_def_CARinf_WTA.set_recombination_scheme(fj.WTA_pt_scheme)

		self.sd0 = fjcontrib.SoftDrop(0, 0.1, self.jet_R0)

		self.fout = ROOT.TFile(self.output, 'recreate')
		self.fout.cd()
		self.tn = ROOT.TNtuple('jadiffs', 'jadiffs', 'pt:dR_CAR0:dR_CARinf:dR_CA_R0_inf:dR_CA_R0_inf_WTA:dR_CA_R0_inf_SD:dR_ues_CA_R0_inf:dR_ues_CA_R0_inf_WTA:dR_ues_CA_R0_inf_SD')
		self.tn_check = ROOT.TNtuple('tn_check', 'tn_check', 'pt:delta_pt:dR_CAR0_WTA:dR_CAR0_SD:dR_CAR0_WTA_SD:dR_CARinf_WTA:dR_CARinf_SD:dR_CARinf_WTA_SD')

		if self.dump:
			self.pickle_std = JetPickleIO(fname='{}_std.pkl'.format(self.output.replace('.root', '')))
			self.pickle_caR0SD = JetPickleIO(fname='{}_caR0SD.pkl'.format(self.output.replace('.root', '')))
			self.pickle_caRinfSD = JetPickleIO(fname='{}_caRinfSD.pkl'.format(self.output.replace('.root', '')))


	def run(self, pythia):
		self.parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, False)
		self.parts_pythia_h_selected = self.parts_selector_h(self.parts_pythia_h)
		self.jets_h = fj.sorted_by_pt(self.jet_selector(self.jet_def(self.parts_pythia_h_selected)))

		if len(self.jets_h) < 1:
			return False

		if self.embedding:
			e = self.data.next_event()
		# print('number of particles in the UE', len(e['particles']))
		# _py_ue = [p for p in e['particles']]
		# _ = [_py_ue.append(p) for p in self.parts_pythia_h_selected]
		_py_ue = fj.vectorPJ()
		if self.embedding:
			_ = [ _py_ue.push_back(p) for p in e['particles']]
		_ = [ _py_ue.push_back(p) for p in self.parts_pythia_h_selected]
		self.jets_h_w_ue = fj.sorted_by_pt(self.jet_def(_py_ue))

		_py_ue_subtr = self.cs.process_event(_py_ue)
		self.jets_h_w_ue_subtr = fj.sorted_by_pt(self.jet_def(_py_ue_subtr))

		# print_jet(self.jets_h[0])
		# print_jet(self.jets_h_w_ue[0])
		# print_jet(self.jets_h_w_ue_subtr[0])

		for j in self.jets_h:
			j_ue 		= None
			j_ue_subtr 	= None
			for _j_ue in self.jets_h_w_ue:
				mpt_ue = fjtools.matched_pt(_j_ue, j)
				if mpt_ue > 0.5:
					j_ue = _j_ue
			for _j_ue_subtr in self.jets_h_w_ue_subtr:
				mpt_ue_subtr = fjtools.matched_pt(_j_ue_subtr, j)
				if mpt_ue_subtr > 0.5:
					j_ue_subtr = _j_ue_subtr
			if j_ue and j_ue_subtr:
				j_CAR0 				= fj.sorted_by_pt(self.jet_def_CAR0		(j.constituents()))[0]
				j_CARinf 			= fj.sorted_by_pt(self.jet_def_CARinf	(j.constituents()))[0]
				j_ue_CAR0 			= fj.sorted_by_pt(self.jet_def_CAR0		(j_ue.constituents()))[0]
				j_ue_CARinf 		= fj.sorted_by_pt(self.jet_def_CARinf	(j_ue.constituents()))[0]
				j_ue_subtr_CAR0 	= fj.sorted_by_pt(self.jet_def_CAR0		(j_ue_subtr.constituents()))[0]
				j_ue_subtr_CARinf 	= fj.sorted_by_pt(self.jet_def_CARinf	(j_ue_subtr.constituents()))[0]

				j_CAR0_WTA 				= fj.sorted_by_pt(self.jet_def_CAR0_WTA		(j.constituents()))[0]
				j_CARinf_WTA 			= fj.sorted_by_pt(self.jet_def_CARinf_WTA	(j.constituents()))[0]
				j_ue_CAR0_WTA 			= fj.sorted_by_pt(self.jet_def_CAR0_WTA		(j_ue.constituents()))[0]
				j_ue_CARinf_WTA 		= fj.sorted_by_pt(self.jet_def_CARinf_WTA	(j_ue.constituents()))[0]
				j_ue_subtr_CAR0_WTA 	= fj.sorted_by_pt(self.jet_def_CAR0_WTA		(j_ue_subtr.constituents()))[0]
				j_ue_subtr_CARinf_WTA 	= fj.sorted_by_pt(self.jet_def_CARinf_WTA	(j_ue_subtr.constituents()))[0]

				j_CAR0_SD = self.sd0.result(j_CAR0)
				j_CARinf_SD = self.sd0.result(j_CARinf)
				j_ue_CAR0_SD = self.sd0.result(j_ue_CAR0)
				j_ue_CARinf_SD = self.sd0.result(j_ue_CARinf)
				j_ue_subtr_CAR0_SD = self.sd0.result(j_ue_subtr_CAR0)
				j_ue_subtr_CARinf_SD = self.sd0.result(j_ue_subtr_CARinf)

				# self.tn.Fill(j.pt(), 
                #  			 j.delta_R(j_CAR0), 
                #      		 j.delta_R(j_CARinf),
                #  			 j_CAR0.delta_R(j_CARinf),
                #      		 j_ue_CAR0.delta_R(j_ue_CARinf),
                #         	 j_ue_subtr_CAR0.delta_R(j_ue_subtr_CARinf))

				self.tn.Fill(j.pt(), 
							 j.delta_R(j_CAR0),
							 j.delta_R(j_CARinf),

                			 j_CAR0.delta_R(j_CARinf),
                			 j_CAR0_WTA.delta_R(j_CARinf_WTA),
                			 j_CAR0_SD.delta_R(j_CARinf_SD),

                			 j_ue_subtr_CAR0.delta_R(j_ue_subtr_CARinf),
                			 j_ue_subtr_CAR0_WTA.delta_R(j_ue_subtr_CARinf_WTA),
                			 j_ue_subtr_CAR0_SD.delta_R(j_ue_subtr_CARinf_SD)
				)
    
				self.tn_check.Fill(j.pt(), j_ue_subtr.pt() - j.pt(),
									j.delta_R(j_CAR0_WTA), j.delta_R(j_CAR0_SD), j_CAR0_WTA.delta_R(j_CAR0_SD),
									j.delta_R(j_CARinf_WTA), j.delta_R(j_CARinf_SD), j_CARinf_WTA.delta_R(j_CARinf_SD)
                )

				if self.dump:
					if j.delta_R(j_CAR0_SD) > 0.05:
						print('adding a jet to pickle')
						self.pickle_std.add_jet(j)
						self.pickle_caR0SD.add_jet(j_CAR0_SD)
						self.pickle_caRinfSD.add_jet(j_CARinf_SD)

				# print_jet(j)
				# print_jet(j_ue)
				# print_jet(j_ue_subtr)

				# run the C/A with R>>R0 and R==R0
				# now run the jet axis algo (SD, WTA)
				# plot the residual dR's
				# dR vs Delta(dR)/dR - to get the percentile diffs

		return True

	def finalize(self):
		self.fout.Write()
		self.fout.Close()
		print ('[i] written ', self.fout.GetName(), file=sys.stderr)
		self.pickle_std.write_to_file()
		self.pickle_caR0SD.write_to_file()
		self.pickle_caRinfSD.write_to_file()


 
def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	parser.add_argument('--output', default=Analysis._default_file_output, type=str)
	parser.add_argument('--datalist', help='run through a file list', default='', type=str)
	parser.add_argument('--alpha', default=0, type=float)
	parser.add_argument('--dRmax', default=0.25, type=float)
	parser.add_argument('--sd-zcut', default=0.1, type=float)
	parser.add_argument('--jet-R0', default=0.2, type=float)
	parser.add_argument('--overwrite', help="overwrite output", default=False, action='store_true')
	parser.add_argument('--dump', help="dump some jets to pickle", default=False, action='store_true')
	parser.add_argument('--jet-pt-min', help='remove jets below the cut', default=20., type=float)
	parser.add_argument('--jet-pt-max', help='remove jets above the cut', default=40., type=float)
	parser.add_argument('--nev', help='number of events to run', default=100, type=int)
	args = parser.parse_args()
	print(args)
	print(args.datalist)
	
	an = Analysis(args=args)
	print(an)

	fj.ClusterSequence.print_banner()

	mycfg = []
	ssettings = f"--py-ecm 5000 --user-seed=100000 --nev {args.nev} --py-hardQCD --py-pthatmin {args.jet_pt_min}"
	args = get_args_from_settings(ssettings)
	pythia_hard = pyconf.create_and_init_pythia_from_args(args, mycfg)
  
	for n in tqdm(range(args.nev)):
		while 1:	
			if not pythia_hard.next():
				continue
			if an.run(pythia_hard):
				break
			
  
	pythia_hard.stat()

	an.finalize()


if __name__ == '__main__':
	main()
