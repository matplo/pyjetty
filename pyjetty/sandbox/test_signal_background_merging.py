#!/usr/bin/env python

import fastjet as fj
import fjcontrib
import fjext
import fjtools
import ecorrel

import tqdm
import argparse
import os
import numpy as np

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from pyjetty.mputils import treewriter

from pyjetty.mputils import BoltzmannEvent

import pythia8
import pythiafjext
import pythiaext
from heppy.pythiautils import configuration as pyconf

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--embed', help='run embedding from a file list', default='', type=str)
	args = parser.parse_args()

	if args.nev < 1:
		args.nev = 1

	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	part_selection = [pythiafjext.kFinal, pythiafjext.kCharged]

	max_eta = 1.
	be = None
	embd = None
	if len(args.embed) > 0:
		embd = DataBackgroundIO(file_list=args.embed)
		print(embd)
	else:
		be = BoltzmannEvent(mean_pt=0.7, multiplicity=200 * max_eta * 2, max_eta=max_eta, max_pt=100)
		print(be)

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.6
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_def_emb = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(10.0) & fj.SelectorAbsEtaMax(1)

	for iev in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		parts = []
		parts = pythiafjext.vectorize_select(pythia, part_selection, 0, False)
		jets = fj.sorted_by_pt(jet_selector(jet_def(parts)))

		if embd:
			bg_parts = embd.load_event(offset=10000)
		else:
			bg_parts = be.generate(offset=10000)

		for j in jets:
			embedded_jet = ecorrel.merge_signal_background_pjvectors(ecorrel.constituents_as_vector(j), bg_parts, 1., 10000)
			print(len(j.constituents()), bg_parts.size(), embedded_jet.size())
			_uidxs_pythia = [_p.user_index() for _p in ecorrel.constituents_as_vector(j)]
			_uidxs = [_p.user_index() for _p in embedded_jet]
			print(_uidxs_pythia, _uidxs)


	pythia.stat()

	# fout.Write()
	# fout.Close()
	# print('[i] written ', fout.GetName())


if __name__ == '__main__':
	main()

