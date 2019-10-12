#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os
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
from data_io import DataBackgroundIO
from jet_analysis import fill_tree_matched, JetAnalysis

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

import ROOT
ROOT.gROOT.SetBatch(True)


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
		embd = DataBackgroundIO(file_list=args.embed)
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
			tmp = [fill_tree_matched(sjet, ej, tw, sd, rho, iev, pythia.info.sigmaGen()) for ej in ja.jets]

	pythia.stat()
	outf.Write()
	outf.Close()
	print('[i] written', outf.GetName())


if __name__ == '__main__':
	main()
