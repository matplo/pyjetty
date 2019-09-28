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

class DataEvent(object):
	def __init__(self, particles, run_number, ev_id):
		self.particles = particles
		self.run_number = run_number
		self.ev_id = ev_id

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
		self.df_events = None
		self.event_tree = None
		self.event_df_orig = None
		self.event_df = None
		self.track_tree = None
		self.track_tree_name = None

	def get_event(self, df_tracks):
		# Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of pt,eta,phi
		event = DataEvent([], -1, -1)
		event.particles = fjext.vectorize_pt_eta_phi(df_tracks['ParticlePt'].values, df_tracks['ParticleEta'].values, df_tracks['ParticlePhi'].values)  
		if len(event.particles) > 0:
			event.run_number = float(df_tracks['run_number'].values[0])
			event.ev_id = float(df_tracks['ev_id'].values[0])
		else:
			event.run_number = -1
			event.ev_id = -1
		return event

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
		self.df_events = self.track_df_grouped.apply(self.get_event)
		self.event_number = self.event_number + len(self.df_events)


class MPAnalysisDriver(mputils.MPBase):
	def __init__(self, **kwargs):
		self.configure_constants(	analyses_list=[], 
									file_list=[],
									tree_name='tree_Particle')
		super(MPAnalysisDriver, self).__init__(**kwargs)
		self.eventIO = MPJetAnalysisFileIO(file_input=None, tree_name=self.tree_name)
		print(self)

	def run(self):
		for file_input in self.file_list:
			if os.path.isfile(file_input):
				print('[i] processing', file_input)
				self.eventIO.load_file(file_input)
				r = [self.process_event(event) for event in self.eventIO.df_events]
			else:
				print('[w] skip - not a file? ', file_input)

	def process_event(self, event):
		for an in self.analyses_list:
			an.process_event(event)


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
									show_progress=False,
									name="MPJetAnalysis",
									output_prefix='./')
		self.jet_eta_max = self.particle_eta_max - self.jet_R * 1.05
		self.jet_def = fj.JetDefinition(self.jet_algorithm, self.jet_R)
		self.jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(self.jet_eta_max)
		self.particle_selector = fj.SelectorPtMin(0.15) & fj.SelectorAbsRapMax(self.particle_eta_max)
		self.jet_area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(self.particle_eta_max))
		self.bg_estimator = None
		self.constit_subtractor = None
		self.jet_output = None
		self.event_output = None
		self.fout = None
		super(MPJetAnalysis, self).__init__(**kwargs)
		self.filename_output ='{}{}.root'.format(self.output_prefix, self.name)
		if self.fout is None:
			self.fout = ROOT.TFile(self.filename_output, 'recreate')
		if self.jet_output is None:
			self.jet_output = treewriter.RTreeWriter(tree_name='tjet', fout=self.fout)
		if self.event_output is None:
			self.event_output = treewriter.RTreeWriter(tree_name='tev', fout=self.fout)
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

	def __del__(self):
		print('[i] writing {}'.format(self.fout.GetName()))
		self.fout.Write()
		self.fout.Purge()
		self.fout.Close()

	def process_event(self, event):
		# assume this is done already (when creating the ntuples):
		# self.fj_particles_selected = self.particle_selector(fj_particles)
		if self.constit_subtractor:
			self.fj_particles_selected = self.constit_subtractor.process_event(event.particles)
		else:
			self.fj_particles_selected = event.particles
		self.cs = fj.ClusterSequenceArea(self.fj_particles_selected, self.jet_def, self.jet_area_def)
		# self.cs = fj.ClusterSequence(self.fj_particles_selected, self.jet_def)
		self.jets = fj.sorted_by_pt(self.cs.inclusive_jets())
		self.rho = 0
		if self.bg_estimator:
			self.bg_estimator.set_particles(self.fj_particles_selected);
			self.rho = self.bg_estimator.rho()
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

		if self.event_output:
			self.event_output.fill_branch('ev_id', event.ev_id)
			self.event_output.fill_branch('run_number', event.run_number)
			self.event_output.fill_branch('rho', self.rho)
			self.event_output.fill_tree()

		if self.jet_output:
			self.jet_output.fill_branch('ev_id', event.ev_id)
			self.jet_output.fill_branch('run_number', event.run_number)
			self.jet_output.fill_branch('rho', self.rho)
			self.jet_output.fill_branch('jet', [j for j in self.jets_detok])
			self.jet_output.fill_branch('jet_ptsub', [j.pt() - (self.rho * j.area()) for j in self.jets_detok])
			for isd, sd in enumerate(self.sds):
				bname = 'jet_sd{}pt'.format(self.sd_betas[isd])
				self.jet_output.fill_branch(bname, [j.pt() for j in self.sd_jets[isd]])
				bname = 'jet_sd{}zg'.format(self.sd_betas[isd])
				self.jet_output.fill_branch(bname, [j.z for j in self.sd_jets_info[isd]])
				bname = 'jet_sd{}Rg'.format(self.sd_betas[isd])
				self.jet_output.fill_branch(bname, [j.dR for j in self.sd_jets_info[isd]])
				bname = 'jet_sd{}thetag'.format(self.sd_betas[isd])
				self.jet_output.fill_branch(bname, [j.dR/self.jet_R for j in self.sd_jets_info[isd]])
			self.jet_output.fill_tree()
		return True


def analyze_file_list(file_inputs=[], output_prefix='rg', tree_name='tree_Particle'):

	anl = []

	bg_rho_range = fj.SelectorAbsRapMax(0.9)
	bg_jet_def = fj.JetDefinition(fj.kt_algorithm, 0.4)
	bg_area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(0.9))
	bg_estimator = fj.JetMedianBackgroundEstimator(bg_rho_range, bg_jet_def, bg_area_def)
	print('[i]', bg_estimator.description())
	an_std = MPJetAnalysis(output_prefix=output_prefix, name='std', bg_estimator=bg_estimator)
	anl.append(an_std)

	g_bg_estimator = fj.GridMedianBackgroundEstimator(0.9, 0.4)
	print('[i]', g_bg_estimator.description())
	an_gbg = MPJetAnalysis(output_prefix=output_prefix, name='bgb', bg_estimator=g_bg_estimator)
	anl.append(an_gbg)

	cs004 = csubtractor.CEventSubtractor(alpha=0, max_distance=0.4, max_eta=0.9, bge_rho_grid_size=0.25, max_pt_correct=100)
	an_cs_004 = MPJetAnalysis(output_prefix=output_prefix, name='cs004', bg_estimator=cs004.bge_rho, constit_subtractor=cs004)
	anl.append(an_cs_004)

	cs404 = csubtractor.CEventSubtractor(alpha=4, max_distance=0.4, max_eta=0.9, bge_rho_grid_size=0.25, max_pt_correct=100)
	an_cs_404 = MPJetAnalysis(output_prefix=output_prefix, name='cs404', bg_estimator=cs404.bge_rho, constit_subtractor=cs404)
	anl.append(an_cs_404)

	ad = MPAnalysisDriver(file_list=file_inputs, analyses_list = anl, tree_name=tree_name)
	ad.run()

	print()
	print('    ---')
	print()
	print(ad)	
	print('[i] done.')

def main():
	parser = argparse.ArgumentParser(description='analyze PbPb data', prog=os.path.basename(__file__))
	parser.add_argument('--input', default="", type=str)
	parser.add_argument('--outprefix', default="matPbPb_out_jets", type=str)

	args = parser.parse_args()

	filename, file_extension = os.path.splitext(args.input)

	if file_extension == '.root':
		analyze_file_list([args.input], args.outprefix)
	else:
		if args.input_file:
			clines = []
			with open(args.input_file) as f:
				clines = [l.strip('\n') for l in f.readlines()]
			analyze_file_list(clines, args.outprefix)

if __name__ == '__main__':
	main()

