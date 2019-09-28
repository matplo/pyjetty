#!/usr/bin/env python

import mputils
import csubtractor
import treewriter

import uproot
import pandas as pd
import numpy as np

import fastjet as fj
import fjcontrib
import fjext

import ROOT
ROOT.gROOT.SetBatch(True)

import joblib
import argparse
import os

class MPJetAnalysisFileIO(mputils.MPBase):
	def __init__(self, **kwargs):
		self.configure_constants(file_input = None, tree_name='tree_Particle')
		self.event_number = 0
		self.event_tree_name = 'PWGHF_TreeCreator/tree_event_char'
		super(MPJetAnalysisFileIO, self).__init__(**kwargs)
		print(self)
		self.reset_dfs()

	def reset_dfs(self):
		self.track_df_orig = None
		self.track_df = None
		self.track_df_grouped = None
		self.df_fjparticles = None
		self.event_tree = None
		self.event_df_orig = None
		self.event_df = None
		self.track_tree = None
		self.track_tree_name = None

	def get_fjparticles(self, df_tracks):
		# Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of pt,eta,phi
		fj_particles = fjext.vectorize_pt_eta_phi(df_tracks['ParticlePt'].values, df_tracks['ParticleEta'].values, df_tracks['ParticlePhi'].values)  
		return fj_particles

	def load_file(self, file_input, tree_name='tree_Particle'):
		self.file_input = file_input
		self.tree_name = tree_name
		print('[i]', self.__class__.__name__, self.file_input, self.tree_name)
		self.reset_dfs()
		# Load event tree into dataframe, and apply event selection
		self.event_tree = uproot.open(file_input)[self.event_tree_name]
		if not self.event_tree:
			print('[e] Tree {} not found in file {}'.format(self.tree_event_name, file_input))
		self.event_df_orig = self.event_tree.pandas.df(['run_number', 'ev_id', 'z_vtx_reco','is_ev_rej'])
		self.event_df_orig.reset_index(drop=True)
		self.event_df = self.event_df_orig.query('is_ev_rej == 0')
		self.event_df.reset_index(drop=True)
		# Load track tree into dataframe
		self.track_tree_name = 'PWGHF_TreeCreator/{}'.format(tree_name)
		self.track_tree = uproot.open(file_input)[self.track_tree_name]
		if not self.track_tree:
			print('[e] Tree {} not found in file {}'.format(tree_name, file_input))
		self.track_df_orig = self.track_tree.pandas.df()
		# Merge event info into track tree
		self.track_df = pd.merge(self.track_df_orig, self.event_df, on=['run_number', 'ev_id'])
		# (i) Group the track dataframe by event
		#     track_df_grouped is a DataFrameGroupBy object with one track dataframe per event
		self.track_df_grouped = self.track_df.groupby(['run_number','ev_id'])
		# (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet particles
		self.df_fjparticles = self.track_df_grouped.apply(self.get_fjparticles)
		self.event_number = self.event_number + len(self.df_fjparticles)


class MPJetAnalysis(mputils.MPBase):
	def __init__(self, **kwargs):
		self.configure_constants(	particle_eta_max=0.9,
									particle_pt_max=100,  # reject jets with constituent with pT > max_particle_pt
									jet_R=0.4, 
									jet_algorithm=fj.antikt_algorithm, 
									jet_pt_min=-200,
									bge_rho_grid_size=-1, # no background subtraction for -1
									sd_z_cuts=[0.1],
									sd_betas=[0, 2],
									show_progress=False)
		self.jet_eta_max = self.particle_eta_max - self.jet_R * 1.05
		self.jet_def = fj.JetDefinition(self.jet_algorithm, self.jet_R)
		self.jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(self.jet_eta_max)
		self.particle_selector = fj.SelectorPtMin(0.15) & fj.SelectorAbsRapMax(self.particle_eta_max)
		self.jet_area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(self.particle_eta_max))
		super(MPJetAnalysis, self).__init__(**kwargs)
		fj.ClusterSequence.print_banner()
		print()
		print(self)
		print(' . particle selector', self.particle_selector)
		print(' . jet selector', self.jet_selector)
		print(' . jet definition:', self.jet_def)
		# define softdrops
		self.sds = []
		for beta in self.sd_betas:
			for zcut in self.sd_z_cuts:
				_sd = fjcontrib.SoftDrop(beta, zcut, self.jet_R)
				self.sds.append(_sd)
		for _sd in self.sds:
			print(' . sd:', _sd.description())

	def analyze_file_data(self, df_fjparticles, bg_estimator=None, csub=None, embed_parts=None,
							jet_root_output=None):
		# Use list comprehension to do jet-finding and fill histograms
		result = [self.analyze_jets_from_particles(fj_particles, bg_estimator, csub, embed_parts, jet_root_output) for fj_particles in df_fjparticles]

	def analyze_jets_from_particles(self, fj_particles, bg_estimator=None, csub=None, embed_parts=None,
									jet_root_output=None):
		if self.show_progress: mputils.CursorSpin()
		# assume this is done already (when creating the ntuples):
		# self.fj_particles_selected = self.particle_selector(fj_particles)
		if csub:
			self.fj_particles_selected = csub.process_event(fj_particles)
		else:
			self.fj_particles_selected = fj_particles
		self.cs = fj.ClusterSequenceArea(self.fj_particles_selected, self.jet_def, self.jet_area_def)
		# self.cs = fj.ClusterSequence(self.fj_particles_selected, self.jet_def)
		self.jets = fj.sorted_by_pt(self.cs.inclusive_jets())
		self.rho = 0
		if bg_estimator:
			bg_estimator.set_particles(fj_particles);
			self.rho = bg_estimator.rho()
		if self.jets[0].pt() - self.rho * self.jets[0].area() < self.jet_pt_min:
			return False
		self.jets_selected = self.jet_selector(self.jets)
		# remove jets containing bad (high pT) particles
		self.jets_detok = list(filter(lambda j: max([t.pt() for t in j.constituents()]) < self.particle_pt_max, self.jets_selected))
		# print(len(self.jets_selected), len(self.jets_detok))
		self.sd_jets = []
		self.sd_jets_info = []
		for isd, sd in enumerate(self.sds):
			self.sd_jets.append([])
			self.sd_jets_info.append([])
		for j in self.jets_detok:
			for isd, sd in enumerate(self.sds):
				jet_sd = sd.result(j)
				self.sd_jets[isd].append(jet_sd)
				self.sd_jets_info[isd].append(fjcontrib.get_SD_jet_info(jet_sd))

		if jet_root_output:
			jet_root_output.fill_branch('jet', [j for j in self.jets_detok])
			for isd, sd in enumerate(self.sds):
				bname = 'jet_sd{}pt'.format(self.sd_betas[isd])
				jet_root_output.fill_branch(bname, [j.pt() for j in self.sd_jets[isd]])
				bname = 'jet_sd{}zg'.format(self.sd_betas[isd])
				jet_root_output.fill_branch(bname, [j.z for j in self.sd_jets_info[isd]])
				bname = 'jet_sd{}Rg'.format(self.sd_betas[isd])
				jet_root_output.fill_branch(bname, [j.dR for j in self.sd_jets_info[isd]])
				bname = 'jet_sd{}thetag'.format(self.sd_betas[isd])
				jet_root_output.fill_branch(bname, [j.dR/self.jet_R for j in self.sd_jets_info[isd]])
			jet_root_output.fill_tree()
		return True

	def df(self):
		dfout = pd.DataFrame({	'pt' : [j.pt() for j in self.jets_detok],
								'a' : [j.area() for j in self.jets_detok],
								'eta' : [j.eta() for j in self.jets_detok],
								'phi' : [j.phi() for j in self.jets_detok]})
		for isd, sd in enumerate(self.sds):
			dfout['sd{}pt'.format(self.sd_betas[isd])] = [j.pt() for j in self.sd_jets[isd]]
			dfout['sd{}zg'.format(self.sd_betas[isd])] = [j.z for j in self.sd_jets_info[isd]]
			dfout['sd{}Rg'.format(self.sd_betas[isd])] = [j.dR for j in self.sd_jets_info[isd]]
			dfout['sd{}thetag'.format(self.sd_betas[isd])] = [j.dR/self.jet_R for j in self.sd_jets_info[isd]]

		dfout['ptsub'] = [j.pt() - self.rho * j.area() for j in self.jets_detok]
		return dfout


def analyze_file_list(file_inputs=[], output_prefix='jets', tree_name='tree_Particle'):
	jet_tw = treewriter.RTreeWriter(tree_name='t', file_name='{}_std.root'.format(output_prefix))

	eventIO = MPJetAnalysisFileIO(file_input=None, tree_name=tree_name)
	# cs008 = csubtractor.CEventSubtractor(alpha=0, max_distance=0.8, max_eta=0.9, bge_rho_grid_size=0.25, max_pt_correct=100)
	cs004 = csubtractor.CEventSubtractor(alpha=0, max_distance=0.4, max_eta=0.9, bge_rho_grid_size=0.25, max_pt_correct=100)
	cs404 = csubtractor.CEventSubtractor(alpha=4, max_distance=0.4, max_eta=0.9, bge_rho_grid_size=0.25, max_pt_correct=100)
	# bg_estimator = fj.GridMedianBackgroundEstimator(0.9, 0.25)
	# print('[i]', bg_estimator.description())
	an = MPJetAnalysis()

	output_std = None
	output_cs004 = None
	output_cs404 = None
	for file_input in file_inputs:
		eventIO.load_file(file_input)
		# rv = an.analyze_file_data(df_fjparticles=eventIO.df_fjparticles, bg_estimator=None, csub=cs008, embed_parts=None)
		_rv = an.analyze_file_data(df_fjparticles=eventIO.df_fjparticles, bg_estimator=None, csub=cs004, embed_parts=None)
		if output_cs004 is None:
			output_cs004 = an.df()
		else:
			output_cs004.append(an.df(), ignore_index = True)

		_rv = an.analyze_file_data(df_fjparticles=eventIO.df_fjparticles, bg_estimator=None, csub=cs404, embed_parts=None)
		if output_cs404 is None:
			output_cs404 = an.df()
		else:
			output_cs404.append(an.df(), ignore_index = True)

		_rv = an.analyze_file_data(df_fjparticles=eventIO.df_fjparticles, 
									bg_estimator=cs404.bge_rho, 
									csub=None, 
									embed_parts=None,
									jet_root_output=jet_tw)
		if output_std is None:
			output_std = an.df()
		else:
			output_std.append(an.df(), ignore_index = True)

	print('[i] analyzed events:', eventIO.event_number)

	joblib.dump(output_std, '{}_std.pd'.format(output_prefix))
	joblib.dump(output_cs004, '{}_cs004.pd'.format(output_prefix))
	joblib.dump(output_cs404, '{}_cs404.pd'.format(output_prefix))

	jet_tw.write_and_close()

def main():
	parser = argparse.ArgumentParser(description='analyze PbPb data', prog=os.path.basename(__file__))
	parser.add_argument('--input', required=True, default="", type=str)
	parser.add_argument('--outprefix', default="matPbPb_out_jets", type=str)

	args = parser.parse_args()
	
	analyze_file_list([args.input], args.outprefix)

if __name__ == '__main__':
	main()

