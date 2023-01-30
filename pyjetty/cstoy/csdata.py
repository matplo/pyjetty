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
import time

from pyjetty.mputils.mputils import logbins
from pyjetty.mputils.mputils import MPBase
from pyjetty.mputils.boltzmann import BoltzmannEvent
from pyjetty.mputils.csubtractor import CEventSubtractor
from pyjetty.mputils.treewriter import RTreeWriter

from alice_efficiency import AliceChargedParticleEfficiency
from pyjetty.mputils.data_io import DataIO
from pyjetty.mputils.jet_analysis import fill_tree_data, JetAnalysisWithRho

import ROOT
ROOT.gROOT.SetBatch(True)


def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	parser.add_argument('--output', default="output.root", type=str)
	parser.add_argument('datalist', help='run through a file list', default='', type=str)
	parser.add_argument('--alpha', default=0, type=float)
	parser.add_argument('--dRmax', default=0.0, type=float)
	parser.add_argument('--zcut', default=0.1, type=float)
	parser.add_argument('--overwrite', help="overwrite output", default=False, action='store_true')
	parser.add_argument('--benchmark', help='benchmark pthat setting - 80 GeV', default=False, action='store_true')
	parser.add_argument('--jetptcut', help='remove jets below the cut', default=-100, type=float)
	parser.add_argument('--nev', help='number of events to run', default=0, type=int)
	args = parser.parse_args()

	if args.output == 'output.root':
		args.output = 'output_data_alpha_{}_dRmax_{}_SDzcut_{}.root'.format(args.alpha, args.dRmax, args.zcut)
		if args.jetptcut > -100:
			args.output = 'output_data_alpha_{}_dRmax_{}_SDzcut_{}_jpt_{}.root'.format(args.alpha, args.dRmax, args.zcut, args.jetptcut)

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

	if args.benchmark:
		jet_selector = fj.SelectorPtMin(80.0) & fj.SelectorPtMax(100.0) & fj.SelectorAbsEtaMax(max_eta - 1.05 * jet_R0)
		# jet_selector_cs = fj.SelectorPtMin(50.0) & fj.SelectorAbsEtaMax(max_eta - 1.05 * jet_R0)
	else:
		jet_selector = fj.SelectorAbsEtaMax(max_eta - 1.05 * jet_R0)

	sd_zcut = args.zcut
	sd = fjcontrib.SoftDrop(0, sd_zcut, jet_R0)
	ja = JetAnalysisWithRho(jet_R=jet_R0, jet_algorithm=fj.antikt_algorithm, particle_eta_max=max_eta)
	data = DataIO(file_list=args.datalist)
	print(data)
	cs = None
	if args.dRmax > 0:
		cs = CEventSubtractor(alpha=args.alpha, max_distance=args.dRmax, max_eta=max_eta, bge_rho_grid_size=0.25, max_pt_correct=100)
		print(cs)
	parts_selector = fj.SelectorAbsEtaMax(max_eta)

	outf = ROOT.TFile(args.output, 'recreate')
	outf.cd()
	t = ROOT.TTree('t', 't')
	tw = RTreeWriter(tree=t)

	# need to change this for data to drive...
	delta_t = 1e-6
	start_t = time.time()
	iev = 0
	while data.load_event():
		iev = iev + 1
		_data_parts = data.particles
		if cs:
			cs_parts = cs.process_event(_data_parts)
			rho = cs.bge_rho.rho()
			ja.analyze_event(cs_parts)
		else:
			ja.analyze_event(_data_parts)
			rho = ja.rho
		tmp = [fill_tree_data(j, tw, sd, rho, iev, 1.) for j in ja.jets if j.pt() > args.jetptcut]
		if iev % 1000 == 0:
			delta_t = time.time()-start_t
			print('[i] processing event', iev, ' - ev/sec = ', iev/delta_t, 'elapsed = ', delta_t)
		if args.nev > 0:
			if iev > args.nev:
				break

	print('[i] processed events', iev, ' - ev/sec = ', iev/delta_t, 'elapsed = ', delta_t)
	outf.Write()
	outf.Close()
	print('[i] written', outf.GetName())


if __name__ == '__main__':
	main()
