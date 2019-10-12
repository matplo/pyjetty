#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os
import numpy as np

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

import numpy as np
import array
import copy
import random
import uproot
import pandas as pd

from pyjetty.mputils import logbins
from pyjetty.mputils import MPBase
from pyjetty.mputils import BoltzmannEvent
from pyjetty.mputils import CEventSubtractor
from pyjetty.mputils import RTreeWriter
from alice_efficiency import AliceChargedParticleEfficiency

import ROOT
ROOT.gROOT.SetBatch(True)

class DataEvent(object):
	def __init__(self, particles, run_number, ev_id):
		self.particles = particles
		self.run_number = run_number
		self.ev_id = ev_id

class MPJetAnalysisFileIO(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(file_input = None, tree_name='tree_Particle')
		self.event_number = 0
		self.event_tree_name = 'PWGHF_TreeCreator/tree_event_char'
		super(MPJetAnalysisFileIO, self).__init__(**kwargs)
		self.reset_dfs()
		if self.file_input:
			self.load_file(self.file_input, self.tree_name)

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
		self.reset_dfs()
		# Load event tree into dataframe, and apply event selection
		self.event_tree = uproot.open(file_input)[self.event_tree_name]
		if not self.event_tree:
			print('[e] Tree {} not found in file {}'.format(self.event_tree_name, file_input))
			return False
		self.event_df_orig = self.event_tree.pandas.df(['run_number', 'ev_id', 'z_vtx_reco','is_ev_rej'])
		self.event_df_orig.reset_index(drop=True)
		self.event_df = self.event_df_orig.query('is_ev_rej == 0')
		self.event_df.reset_index(drop=True)
		# Load track tree into dataframe
		self.track_tree_name = 'PWGHF_TreeCreator/{}'.format(tree_name)
		self.track_tree = uproot.open(file_input)[self.track_tree_name]
		if not self.track_tree:
			print('[e] Tree {} not found in file {}'.format(tree_name, file_input))
			return False
		self.track_df_orig = self.track_tree.pandas.df()
		# Merge event info into track tree
		self.track_df = pd.merge(self.track_df_orig, self.event_df, on=['run_number', 'ev_id'])
		# (i) Group the track dataframe by event
		#     track_df_grouped is a DataFrameGroupBy object with one track dataframe per event
		self.track_df_grouped = self.track_df.groupby(['run_number','ev_id'])
		# (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet particles
		self.df_events = self.track_df_grouped.apply(self.get_event)
		self.event_number = self.event_number + len(self.df_events)
		return True


class DataBackground(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(file_list='PbPb_file_list.txt')	
		super(DataBackground, self).__init__(**kwargs)
		self.current_event_in_file = 0
		self.file_io = None
		self.list_of_files = []
		self.read_file_list()

	def set_file_list(self, sfile):
		self.file_list = sfile
		self.read_file_list()

	def read_file_list(self):
		self.list_of_files = []
		with open(self.file_list) as f:
			self.list_of_files = [fn.strip() for fn in f.readlines()]

	def open_file(self):
		self.file_io = None
		self.current_event_in_file = 0
		afile = random.choice(self.list_of_files)
		print('[i] opening data file', afile)
		self.file_io = MPJetAnalysisFileIO(file_input=afile)
		print('    number of events', len(self.file_io.df_events))

	def load_event(self):
		if self.file_io is None:
			self.open_file()
		if self.file_io is None:
			print('[e] unable to load the background file')
		if self.current_event_in_file >= len(self.file_io.df_events):
			self.current_event_in_file = 0
			self.file_io = None
			return self.load_event()
		event = self.file_io.df_events[self.current_event_in_file]
		_tmp = [p.set_user_index(10000+ip) for ip,p in enumerate(event.particles)]
		# print('loaded event:', self.current_event_in_file)
		self.current_event_in_file = self.current_event_in_file + 1
		self.particles = event.particles
		return self.particles

class JetAnalysis(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(jet_R=0.4, jet_algorithm=fj.antikt_algorithm, particle_eta_max=0.9, particles=None)
		super(JetAnalysis, self).__init__(**kwargs)

		self.bg_rho_range = fj.SelectorAbsEtaMax(self.particle_eta_max)
		self.bg_jet_def = fj.JetDefinition(fj.kt_algorithm, self.jet_R)
		self.bg_area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(self.particle_eta_max))
		self.bg_estimator = fj.JetMedianBackgroundEstimator(self.bg_rho_range, self.bg_jet_def, self.bg_area_def)
		self.rho = 0

		self.jet_eta_max = self.particle_eta_max - self.jet_R * 1.05
		self.jet_def = fj.JetDefinition(self.jet_algorithm, self.jet_R)
		self.jet_selector = fj.SelectorPtMin(0.0) & fj.SelectorAbsEtaMax(self.jet_eta_max)
		self.jet_area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(self.particle_eta_max))

		if self.particles:
			self.analyze_event(self.particles)

	def analyze_event(self, parts):
		if len(parts) < 1:
			self.rho = 0.0
			self.cs = None
			self.jets = []
			self.corr_jet_pt = []
		else:
			self.bg_estimator.set_particles(parts)
			self.rho = self.bg_estimator.rho()
			self.cs = fj.ClusterSequenceArea(parts, self.jet_def, self.jet_area_def)
			self.jets = fj.sorted_by_pt(self.jet_selector(self.cs.inclusive_jets()))
			self.corr_jet_pt = [j.pt() - j.area() * self.rho for j in self.jets]


def matched_pt_constituent(c0, j1):
	_pt = [c.pt() for i, c in enumerate(j1.constituents()) if c.user_index() == c0.user_index()]
	if len(_pt) == 1:
		return _pt[0]
	if len(_pt) > 1:
		print('[e] more than one match of constituents?')
		return _pt[0]
	return 0

def g_matched_pt_constituent(c0, j1):
	g = (c.pt() for i, c in enumerate(j1.constituents()) if c.user_index() == c0.user_index())
	return next(g)

def matched_pt(j0, j1):
	# pts = [g_matched_pt_constituent(c0, j1) for c0 in j0.constituents()]
	pts = [c.pt() 
			for c0 in j0.constituents()
			for i, c in enumerate(j1.constituents()) if c.user_index() == c0.user_index()]
	idxs = [c.user_index()
			for c0 in j0.constituents()
			for i, c in enumerate(j1.constituents()) if c.user_index() == c0.user_index()]
	# if len(idxs) > 0:
	# 	idxs_v = [c.user_index() for c in j1.constituents()]
	# 	print('idxs_vacuum:', sorted(idxs_v))
	# 	print('idxs_match :', sorted(idxs))
	pt_sum = sum(pts)
	return pt_sum / j1.pt()

# 			tn.Fill(i, 
# 					sjet.pt(), sjet.phi(), sjet.eta(), 
# 					j.pt(), j.phi(), j.eta(),
# 					j.delta_R(sjet), j.pt() - sjet.pt(),
# 					sd_info_signal_jet.dR, sd_info_j.dR, sd_info_j.dR - sd_info_signal_jet.dR,
# 					sd_info_signal_jet.z, sd_info_j.z, sd_info_j.z - sd_info_signal_jet.z)

def fill_tree(signal_jet, emb_jet, tw, sd, rho, iev=None, weight=None):
	if matched_pt(emb_jet, signal_jet) <= 0.5:
		return None

	if iev:
		tw.fill_branch('ev_id', iev)
	if weight:
		tw.fill_branch('weight', weight)

	sd_signal_jet = sd.result(signal_jet)
	sd_info_signal_jet = fjcontrib.get_SD_jet_info(sd_signal_jet)
	sd_emb_jet = sd.result(emb_jet)
	sd_info_emb_jet = fjcontrib.get_SD_jet_info(sd_emb_jet)

	tw.fill_branch('rho', rho)

	tw.fill_branch('j', signal_jet)
	tw.fill_branch('sd_j', sd_signal_jet)
	tw.fill_branch('sd_j_z', sd_info_signal_jet.z)
	tw.fill_branch('sd_j_dR', sd_info_signal_jet.dR)

	tw.fill_branch('ej', emb_jet)
	tw.fill_branch('ej_ptc', emb_jet.pt() - emb_jet.area() * rho)
	tw.fill_branch('sd_ej', sd_emb_jet)
	tw.fill_branch('sd_ej_cpt', sd_emb_jet.pt() - sd_emb_jet.area() * rho)
	tw.fill_branch('sd_ej_z', sd_info_emb_jet.z)
	tw.fill_branch('sd_ej_dR', sd_info_emb_jet.dR)

	p1 = fj.PseudoJet()
	p2 = fj.PseudoJet()
	has_parents_signal = sd_signal_jet.has_parents(p1, p2)
	# print('signal_jet:', has_parents, len(p1.constituents()), len(p2.constituents()))
	tw.fill_branch('j_p1', p1)
	tw.fill_branch('j_p2', p2)


	pe1 = fj.PseudoJet()
	pe2 = fj.PseudoJet()
	has_parents_emb = sd_emb_jet.has_parents(pe1, pe2)
	tw.fill_branch('ej_p1', pe1)
	tw.fill_branch('ej_p2', pe2)
	if has_parents_emb:
		tw.fill_branch('ej_p1_ptc', pe1.pt() - pe1.area() * rho)
		tw.fill_branch('ej_p2_ptc', pe2.pt() - pe2.area() * rho)
	else:
		tw.fill_branch('ej_p1_ptc', -1000)
		tw.fill_branch('ej_p2_ptc', -1000)

	mpt1 = -1.0 # not passed SD
	mpt2 = -1.0 # not passed SD

	if has_parents_signal and has_parents_emb:
		mpt1 = matched_pt(pe1, p1)
		mpt2 = matched_pt(pe2, p2)
	tw.fill_branch('mpt1', mpt1)
	tw.fill_branch('mpt2', mpt2)

		# print('signal_jet:', has_parents, len(pe1.constituents()), len(pe2.constituents()))
		# print('emb_jets', has_parents, len(pe1.constituents()), len(pe2.constituents()))

	# for c in pe2.constituents():
	# 	cp1 = fj.PseudoJet()
	# 	cp2 = fj.PseudoJet()
	# 	print(' - ', c.has_parents(cp1, cp2))

	#tw.fill_branch('jsd', sd_j)
	#tw.fill_branch('jm', ej)
	tw.fill_tree()
	return emb_jet

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--ignore-mycfg', help="ignore some settings hardcoded here", default=False, action='store_true')
	parser.add_argument('--output', default="output.root", type=str)
	parser.add_argument('--alpha', default=0, type=float)
	parser.add_argument('--dRmax', default=0.0, type=float)
	parser.add_argument('--zcut', default=0.1, type=float)
	parser.add_argument('--overwrite', help="overwrite output", default=False, action='store_true')
	parser.add_argument('--embed', help='run embedding from a file list', default='', type=str)
	parser.add_argument('--efficiency', help='apply charged particle efficiency', default=False, action='store_true')
	parser.add_argument('--benchmark', help='benchmark pthat setting - 80 GeV', default=False, action='store_true')
	args = parser.parse_args()

	if args.output == 'output.root':
		args.output = 'output_alpha_{}_dRmax_{}_SDzcut_{}.root'.format(args.alpha, args.dRmax, args.zcut)
		if args.py_seed >= 0:
			args.output = 'output_alpha_{}_dRmax_{}_SDzcut_{}_seed_{}.root'.format(args.alpha, args.dRmax, args.zcut, args.py_seed)
		if args.embed:
			args.output = args.output.replace('.root', '_emb.root')
		if args.efficiency:
			args.output = args.output.replace('.root', '_effi.root')

	if os.path.isfile(args.output):
		if not args.overwrite:
			print('[i] output', args.output, 'exists - use --overwrite to do just that...')
			return

	print(args)

	# alice specific
	max_eta = 0.9

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	print(jet_def)

	mycfg = []
	if args.benchmark:
		mycfg = ['PhaseSpace:pThatMin = 80', 'PhaseSpace:pThatMax = -1']
		jet_selector = fj.SelectorPtMin(80.0) & fj.SelectorPtMax(100.0) & fj.SelectorAbsEtaMax(max_eta - 1.05 * jet_R0)
		# jet_selector_cs = fj.SelectorPtMin(50.0) & fj.SelectorAbsEtaMax(max_eta - 1.05 * jet_R0)
	else:
		args.py_biaspow = 4
		args.py_biasref = 10
		jet_selector = fj.SelectorPtMin(5) & fj.SelectorAbsEtaMax(max_eta - 1.05 * jet_R0)
		# jet_selector_cs = fj.SelectorPtMin(50.0) & fj.SelectorAbsEtaMax(max_eta - 1.05 * jet_R0)

	if args.ignore_mycfg:
		mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return

	sd_zcut = args.zcut
	sd = fjcontrib.SoftDrop(0, sd_zcut, jet_R0)

	ja = JetAnalysis(jet_R=jet_R0, jet_algorithm=fj.antikt_algorithm, particle_eta_max=max_eta)

	be = None
	embd = None
	if len(args.embed) > 0:
		embd = DataBackground(file_list=args.embed)
		print(embd)
	else:
		be = BoltzmannEvent(mean_pt=0.6, multiplicity=2000 * max_eta * 2, max_eta=max_eta, max_pt=100)
		print(be)

	cs = None
	if args.dRmax > 0:
		cs = CEventSubtractor(alpha=args.alpha, max_distance=args.dRmax, max_eta=max_eta, bge_rho_grid_size=0.25, max_pt_correct=100)
		print(cs)

	parts_selector = fj.SelectorAbsEtaMax(max_eta)

	if args.nev < 1:
		args.nev = 1

	outf = ROOT.TFile(args.output, 'recreate')
	outf.cd()
	t = ROOT.TTree('t', 't')
	tw = RTreeWriter(tree=t)

	# effi_pp = AliceChargedParticleEfficiency(csystem='pp')
	effi_PbPb = None
	if args.efficiency:
		effi_PbPb = AliceChargedParticleEfficiency(csystem='PbPb')
		print(effi_PbPb)

	### EVENT LOOP STARTS HERE
	for iev in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue

		parts_pythia = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged])
		parts_gen = parts_selector(parts_pythia)
		if effi_PbPb:
			parts = effi_PbPb.apply_efficiency(parts_gen)
		else:
			parts = parts_gen

		signal_jets = fj.sorted_by_pt(jet_selector(jet_def(parts)))
		if len(signal_jets) < 1:
			continue

		for sjet in signal_jets:
			if embd:
				bg_parts = embd.load_event()
			else:
				bg_parts = be.generate(offset=10000)
			full_event = bg_parts
			tmp = [full_event.push_back(psj) for psj in sjet.constituents()]
			if cs:
				cs_parts = cs.process_event(full_event)
				rho = cs.bge_rho.rho()
				ja.analyze_event(cs_parts)
			else:
				ja.analyze_event(full_event)
				rho = ja.rho
			tmp = [fill_tree(sjet, ej, tw, sd, rho, iev, pythia.info.sigmaGen()) for ej in ja.jets]

	pythia.stat()
	outf.Write()
	outf.Close()
	print('[i] written', outf.GetName())


if __name__ == '__main__':
	main()
