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

from pyjetty.mputils import MPBase
from pyjetty.mputils import BoltzmannEvent
from pyjetty.mputils import CEventSubtractor
from pyjetty.mputils import RTreeWriter

import ROOT
ROOT.gROOT.SetBatch(True)


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

def fill_tree(signal_jet, emb_jet, tw, sd, rho, iev=None):
	if matched_pt(emb_jet, signal_jet) <= 0.5:
		return None

	if iev:
		tw.fill_branch('ev_id', iev)

	sd_signal_jet = sd.result(signal_jet)
	sd_info_signal_jet = fjcontrib.get_SD_jet_info(sd_signal_jet)
	sd_emb_jet = sd.result(emb_jet)
	sd_info_emb_jet = fjcontrib.get_SD_jet_info(sd_emb_jet)

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
	parser.add_argument('--py-seed', help='pythia seed', default=-1, type=int)

	args = parser.parse_args()

	if args.output == 'output.root':
		args.output = 'output_alpha_{}_dRmax_{}_SDzcut_{}.root'.format(args.alpha, args.dRmax, args.zcut)
		if args.py_seed >= 0:
			args.output = 'output_alpha_{}_dRmax_{}_SDzcut_{}_seed_{}.root'.format(args.alpha, args.dRmax, args.zcut, args.py_seed)


	if os.path.isfile(args.output):
		if not args.overwrite:
			print('[i] output', args.output, 'exists - use --overwrite to do just that...')
			return

	print(args)

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(80.0) & fj.SelectorPtMax(100.0) & fj.SelectorAbsEtaMax(1 - 1.05 * jet_R0)
	jet_selector_cs = fj.SelectorPtMin(50.0) & fj.SelectorAbsEtaMax(1 - 1.05 * jet_R0)
	# jet_selector = fj.SelectorPtMin(10.0) & fj.SelectorPtMax(20.0) & fj.SelectorAbsEtaMax(1 - 1.05 * jet_R0)
	# jet_selector_cs = fj.SelectorPtMin(0.0) & fj.SelectorAbsEtaMax(1 - 1.05 * jet_R0)
	print(jet_def)

	mycfg = ['PhaseSpace:pThatMin = 80', 'PhaseSpace:pThatMax = -1']
	if args.py_seed >= 0:
		mycfg.append('Random:setSeed=on')
		mycfg.append('Random:seed={}'.format(args.py_seed))

	# mycfg = []
	if args.ignore_mycfg:
		mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return

	sd_zcut = args.zcut
	sd = fjcontrib.SoftDrop(0, sd_zcut, jet_R0)

	max_eta = 1
	be = BoltzmannEvent(mean_pt=0.6, multiplicity=2000 * max_eta * 2, max_eta=max_eta, max_pt=100)
	ja = JetAnalysis(jet_R=jet_R0, jet_algorithm=fj.antikt_algorithm, particle_eta_max=max_eta)
	parts_selector = fj.SelectorAbsEtaMax(max_eta)

	if args.nev < 1:
		args.nev = 1

	outf = ROOT.TFile(args.output, 'recreate')
	outf.cd()
	t = ROOT.TTree('t', 't')
	tw = RTreeWriter(tree=t)

	for iev in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue

		parts_pythia = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged])
		parts = parts_selector(parts_pythia)

		signal_jets = fj.sorted_by_pt(jet_selector(jet_def(parts)))
		if len(signal_jets) < 1:
			continue
		sjet = signal_jets[0]

		bg_parts = be.generate(offset=10000)
		full_event = bg_parts

		lc = [full_event.push_back(psj) for psj in sjet.constituents()]
		# idxs = [psj.user_index() for psj in sjet.constituents()]
		# print('pythia jet:', idxs)

		# cs_signal_jets = fj.sorted_by_pt(jet_selector_cs(jet_def(cs_parts)))

		if args.dRmax > 0:
			cs = CEventSubtractor(alpha=args.alpha, max_distance=args.dRmax, max_eta=max_eta, bge_rho_grid_size=0.25, max_pt_correct=100)
			cs_parts = cs.process_event(full_event)
			rho = cs.bge_rho.rho()
			ja.analyze_event(cs_parts)
		else:
			ja.analyze_event(full_event)
			rho = ja.rho

		r = [fill_tree(sjet, ej, tw, sd, rho, iev) for ej in ja.jets]

		# emb_jets = fj.sorted_by_pt(jet_selector_cs(jet_def(full_event)))
		# r = [fill_tree(sjet, ej, tw, sd, rho, iev) for ej in emb_jets]

		# max_eta_part = max([j.eta() for j in full_event])
		# print ('number of particles', len(full_event), max_eta_part)
		# mean_pt = sum([j.pt() for j in bg_parts]) / len(bg_parts)
		# print ('mean pt in bg', mean_pt)
		# print ('number of CS particles', len(cs_parts))


	pythia.stat()
	outf.Write()
	outf.Close()


if __name__ == '__main__':
	main()
