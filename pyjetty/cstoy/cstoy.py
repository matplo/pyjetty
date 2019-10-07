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

from pyjetty.mputils import BoltzmannEvent
from pyjetty.mputils import CEventSubtractor
from pyjetty.mputils import RTreeWriter

import ROOT
ROOT.gROOT.SetBatch(True)

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
	# print(idxs)
	pt_sum = sum(pts)
	return pt_sum / j1.pt()

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--ignore-mycfg', help="ignore some settings hardcoded here", default=False, action='store_true')
	parser.add_argument('--output', default="output.root", type=str)
	parser.add_argument('--alpha', default=0, type=float)
	parser.add_argument('--dRmax', default=0.25, type=float)
	parts.add_argument('--zcut', default=0.1, type=float)
	parser.add_argument('--overwrite', help="overwrite output", default=False, action='store_true')

	args = parser.parse_args()

	if args.output == 'output.root':
		args.output = 'output_alpha_{}_dRmax_{}_SDzcut_{}.root'.format(args.alpha, args.dRmax, args.zcut)

	if os.path.isfile(args.output):
		if not args.overwrite:
			print('[i] output', args.output, 'exists - use --overwrite to do just that...')
			return

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
	cs = CEventSubtractor(alpha=args.alpha, max_distance=args.dRmax, max_eta=max_eta, bge_rho_grid_size=0.25, max_pt_correct=100)
	be = BoltzmannEvent(mean_pt=0.6, multiplicity=1000 * max_eta * 2, max_eta=max_eta, max_pt=100)
	parts_selector = fj.SelectorAbsEtaMax(max_eta)

	if args.nev < 100:
		args.nev = 100

	outf = ROOT.TFile(args.output, 'recreate')
	outf.cd()
	tn = ROOT.TNtuple('tn', 'tn', 'n:pt:phi:eta:cspt:csphi:cseta:dR:dpt:rg:csrg:z:csz:dRg:dzg')
	hpt = ROOT.TH1F('hpt', 'hpt', 40, 00, 160)
	hptcs = ROOT.TH1F('hptcs', 'hptcs', 40, 00, 160)
	hdpt = ROOT.TH1F('hdpt', 'hdpt', 40, -50, 50)
	hrg = ROOT.TH1F('hrg', 'hrg', 40, -1.1, 0.9)
	hrgcs = ROOT.TH1F('hrgcs', 'hrgcs', 40, -1.1, 0.9)
	hdrg = ROOT.TH1F('hdrg', 'hdrg', 40, -1.1, 0.9)
	hdz = ROOT.TH1F('hdz', 'hdz', 40, -1.1, 0.9)
	hdzz = ROOT.TH2F('hdzz', 'hdzz', 40, -1.1, 0.9, 40, -1.1, 0.9)
	hdphi = ROOT.TH1F('hdphi', 'hdphi', 90, -ROOT.TMath.Pi(), ROOT.TMath.Pi())

	t = ROOT.TTree('t', 't')
	tw = RTreeWriter(tree=t)

	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue

		parts_pythia = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal])
		parts = parts_selector(parts_pythia)
		signal_jets = fj.sorted_by_pt(jet_selector(jet_def(parts)))

		if len(signal_jets) < 1:
			continue

		bg_parts = be.generate(offset=10000)

		full_event = bg_parts
		sjet = signal_jets[0]
		lc = [full_event.push_back(psj) for psj in sjet.constituents()]
		# idxs = [psj.user_index() for psj in sjet.constituents()]
		# print('pythia jet:', idxs)
		cs_parts = cs.process_event(full_event)
		cs_signal_jets = fj.sorted_by_pt(jet_selector_cs(jet_def(cs_parts)))

		emb_jets = fj.sorted_by_pt(jet_selector_cs(jet_def(full_event)))

		# max_eta_part = max([j.eta() for j in full_event])
		# print ('number of particles', len(full_event), max_eta_part)
		# mean_pt = sum([j.pt() for j in bg_parts]) / len(bg_parts)
		# print ('mean pt in bg', mean_pt)
		# print ('number of CS particles', len(cs_parts))

		sd_signal_jet = sd.result(sjet)
		sd_info_signal_jet = fjcontrib.get_SD_jet_info(sd_signal_jet)

		# for j in cs_signal_jets:
		for j in emb_jets:
			if matched_pt(sjet, j) <= 0.5:
				continue
			sd_j = sd.result(j)
			sd_info_j = fjcontrib.get_SD_jet_info(sd_j)
			rho = cs.bge_rho.rho()
			tn.Fill(i, 
					sjet.pt(), sjet.phi(), sjet.eta(), 
					j.pt(), j.phi(), j.eta(),
					j.delta_R(sjet), j.pt() - sjet.pt(),
					sd_info_signal_jet.dR, sd_info_j.dR, sd_info_j.dR - sd_info_signal_jet.dR,
					sd_info_signal_jet.z, sd_info_j.z, sd_info_j.z - sd_info_signal_jet.z)
			hpt.Fill(sjet.pt())
			hptcs.Fill(j.pt())
			hdpt.Fill(j.pt() - sjet.pt())
			hrg.Fill(sd_info_signal_jet.dR)
			hrgcs.Fill(sd_info_j.dR)
			hdphi.Fill(sjet.delta_phi_to(j))

			if sd_info_j.dR > 0 and sd_info_signal_jet.dR > 0:
				hdrg.Fill(sd_info_j.dR - sd_info_signal_jet.dR)
			if sd_info_j.z > 0 and sd_info_signal_jet.z > 0:
				hdz.Fill(sd_info_j.z - sd_info_signal_jet.z)
				hdzz.Fill(sd_info_j.z, sd_info_j.z - sd_info_signal_jet.z)

		# for i,j in enumerate(signal_jets):
		# 	j_sd = sd.result(j)
		# 	sd_info = fjcontrib.get_SD_jet_info(j_sd)
		# 	# print("  |-> SD jet params z={0:10.3f} dR={1:10.3f} mu={2:10.3f}".format(sd_info.z, sd_info.dR, sd_info.mu))

	pythia.stat()
	outf.Write()
	outf.Close()


if __name__ == '__main__':
	main()
