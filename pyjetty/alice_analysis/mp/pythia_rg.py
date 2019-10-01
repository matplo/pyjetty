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
import treewriter

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--output', default=None, type=str)
	parser.add_argument('--jetR', default=0.4, type=float)
	parser.add_argument('--dry', default=False, action='store_true')
	parser.add_argument('--overwrite', default=False, action='store_true')
	args = parser.parse_args()	

	if args.output is None:
		hlevel = 'on_'
		if 'HadronLevel:all' in args.pythiaopts:
			hlevel = args.pythiaopts.split('HadronLevel:all=')[1].split(' ')[0].lower()
			if hlevel == 'on': 
				hlevel = 'on_'
		isrlevel = '__ISR'
		if args.noISR:
			isrlevel = 'noISR'
		mpilevel = '__MPI'
		if args.noMPI:
			mpilevel = 'noMPI'
		args.output = "pythia_rg_jetR{}_pthatmin{}_hadr{}_{}_{}.root".format(args.jetR, args.pthatmin, hlevel.upper(), isrlevel, mpilevel)
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
	sd0 = fjcontrib.SoftDrop(0, 0.1, jet_R0)
	sd1 = fjcontrib.SoftDrop(1, 0.1, jet_R0)
	sd2 = fjcontrib.SoftDrop(2, 0.1, jet_R0)
	print(sd0)
	print(sd2)
	sds = []
	sds.append(sd0)
	sds.append(sd1)
	sds.append(sd2)

	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 100:
		args.nev = 100

	fout = ROOT.TFile(args.output, 'recreate')
	jet_output = treewriter.RTreeWriter(tree_name='tjet', fout=fout)
	event_output = treewriter.RTreeWriter(tree_name='tev', fout=fout)

	for iev in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		# parts_final = pythiafjext.vectorize(pythia, True, -1, 1, False)
		# parts_all = pythiafjext.vectorize_select(pythia, [], False);
		parts_final = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], False);
		# parts_neutral = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kNeutral], False);
		# parts_charged = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], False);
		# parts_charged_visible = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged, pythiafjext.kVisible], False);
		# parts_charged_hadrons = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged, pythiafjext.kHadron], False);
		# print (len(parts_all), len(parts_final), len(parts_charged) + len(parts_neutral), len(parts_charged), len(parts_neutral))

		# parts_hadrons = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kHadron], False);
		# parts_leptons = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kLepton], False);
		# parts_photons = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kPhoton], False);
		# print ('final', len(parts_final), 'h', len(parts_hadrons), 'l', len(parts_leptons), 'gamma', len(parts_photons), 'sum', len(parts_hadrons) + len(parts_leptons) + len(parts_photons))

		parts = parts_final
		jets = jet_selector(jet_def(parts))

		event_output.fill_branch('ev_id', iev)
		event_output.fill_branch('sigma', pythia.info.sigmaGen())
		event_output.fill_branch('sigma_err', pythia.info.sigmaErr())
		event_output.fill_tree()

		for i,j in enumerate(jets):
			jet_output.fill_branch('ev_id', iev)
			jet_output.fill_branch('sigma', pythia.info.sigmaGen())
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