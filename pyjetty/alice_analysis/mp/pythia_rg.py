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

import ROOT
from pyjetty.mputils import RTreeWriter

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--output', default=None, type=str)
	parser.add_argument('--jetR', default=0.4, type=float)
	parser.add_argument('--dry', default=False, action='store_true')
	parser.add_argument('--overwrite', default=False, action='store_true')
	parser.add_argument('--charged', default=False, action='store_true')
	parser.add_argument('--runid', default=0, type=int)
	parser.add_argument('--RreclusterR0', default=True, help="use the same R for soft drop as R for jet finding - default is TRUE", action='store_true')
	args = parser.parse_args()

	if args.output is None:
		hlevel = 'on_'
		if 'HadronLevel:all' in args.pythiaopts:
			hlevel = args.pythiaopts.split('HadronLevel:all=')[1].split(' ')[0].lower()
			if hlevel == 'on': 
				hlevel = 'on_'
		isrlevel = '__ISR'
		if args.py_noISR:
			isrlevel = 'noISR'
		mpilevel = '__MPI'
		if args.py_noMPI:
			mpilevel = 'noMPI'
		args.output = "pythia_rg_jetR{}_pthatmin{}_runid{}_hadr{}_{}_{}.root".format(args.jetR, args.py_pthatmin, args.runid, hlevel.upper(), isrlevel, mpilevel)
		if args.charged:
			args.output = args.output.replace('.root', '_charged.root')
	print('[i] output goes to:', args.output)
	if args.dry:
		return
	if os.path.exists(args.output):
		print ('[w] output exists - skip.')
		if args.overwrite:
			print ('[w] output exists - overwrite ON')
		else:
			return

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	jet_R0 = args.jetR
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(20) & fj.SelectorAbsEtaMax(0.9 - jet_R0)
	print(jet_def)
	sd_jet_def_recluster = None
	sd_reclusterer = None
	if args.RreclusterR0:
		sd_jet_def_recluster = fj.JetDefinition(fj.cambridge_algorithm, jet_R0)
		sd_reclusterer = fjcontrib.Recluster(sd_jet_def_recluster)
		print('[i] SD reclusterer', sd_reclusterer.description())
	sd0 = fjcontrib.SoftDrop(0, 0.1, jet_R0)
	if sd_reclusterer:
		sd0.set_reclustering(True, sd_reclusterer)	
	sd1 = fjcontrib.SoftDrop(1, 0.1, jet_R0)
	if sd_reclusterer:
		sd1.set_reclustering(True, sd_reclusterer)	
	sd2 = fjcontrib.SoftDrop(2, 0.1, jet_R0)
	if sd_reclusterer:
		sd2.set_reclustering(True, sd_reclusterer)	
	print(sd0.description())
	print(sd1.description())
	print(sd2.description())
	sds = []
	sds.append(sd0)
	sds.append(sd1)
	sds.append(sd2)

	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 1:
		args.nev = 1

	fout = ROOT.TFile(args.output, 'recreate')
	jet_output = RTreeWriter(tree_name='tjet', fout=fout)
	event_output = RTreeWriter(tree_name='tev', fout=fout)

	for iev in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		parts = []
		if args.charged:
			parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged]);
		else:
			parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal]);

		jets = jet_selector(jet_def(parts))

		ev_s = pythia.info.sigmaGen()
		ev_s_err = pythia.info.sigmaErr()
		ev_w = pythia.info.weight()
		ev_code = pythia.info.code()
		pthard = pythia.info.pTHat()

		event_output.fill_branch('ev_id', iev)
		event_output.fill_branch('sigma', ev_s)
		event_output.fill_branch('sigmaErr', ev_s_err)
		event_output.fill_branch('weight', ev_w)
		event_output.fill_branch('code', ev_code)
		event_output.fill_branch('pthard', pthard)
		event_output.fill_tree()

		for i,j in enumerate(jets):
			jet_output.fill_branch('ev_id', iev)
			jet_output.fill_branch('sigma', ev_s)
			jet_output.fill_branch('sigmaErr', ev_s_err)
			jet_output.fill_branch('weight', ev_w)
			jet_output.fill_branch('code', ev_code)
			jet_output.fill_branch('pthard', pthard)
			jet_output.fill_branch('jet', j)
			for isd, sd in enumerate(sds):
				j_sd = sd.result(j)
				sd_info = fjcontrib.get_SD_jet_info(j_sd)
				bname = 'jet_sd{}pt'.format(sd.beta()).replace('.0', '_')
				jet_output.fill_branch(bname, j_sd.pt())
				bname = 'jet_sd{}zg'.format(sd.beta()).replace('.0', '_')
				jet_output.fill_branch(bname, sd_info.z)
				bname = 'jet_sd{}Rg'.format(sd.beta()).replace('.0', '_')
				jet_output.fill_branch(bname, sd_info.dR)
				bname = 'jet_sd{}thetag'.format(sd.beta()).replace('.0', '_')
				jet_output.fill_branch(bname, sd_info.dR/jet_R0)
		jet_output.fill_tree()

	pythia.stat()

	print('[i] writing {}'.format(fout.GetName()))
	fout.Write()
	fout.Purge()
	fout.Close()

if __name__ == '__main__':
	main()