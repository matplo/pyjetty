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

from pyjetty.mputils import RTreeWriter
from pyjetty.mputils import pwarning, pinfo, perror


import ROOT
ROOT.gROOT.SetBatch(True)


def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	_default_output_filename = os.path.basename(__file__).replace(".py", "") + "_output.root"
	parser.add_argument('--output', default=_default_output_filename, type=str)
	parser.add_argument('--debug', default=0, type=int)
	args = parser.parse_args()

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.6
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(2.0) & fj.SelectorAbsEtaMax(2)
	print(jet_def)

	jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, jet_R0)
	lund_gen = fjcontrib.LundGenerator(jet_def_lund)
	print(jet_def_lund)
	print(lund_gen)

	outf = ROOT.TFile(args.output, 'recreate')
	outf.cd()
	t = ROOT.TTree('t', 't')
	tw = RTreeWriter(tree=t)

	# mycfg = ['PhaseSpace:pThatMin = 100']
	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 100:
		args.nev = 100
	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		if args.debug:
			pwarning('-- event', i)
		# parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
		parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)
		if args.debug > 5:
			parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kHadron], 0, True)
		if args.debug > 10:
			parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kAny], 0, True)
		if args.debug > 0:
			for p in parts:
				pypart = pythiafjext.getPythia8Particle(p)
				if pypart.name()[:2] == 'D0':
					pinfo(pypart.name(), pypart.id(), pypart.status(), 'final =?', pypart.isFinal())
		jets = jet_selector(jet_def(parts))

		for j in jets:
			isD0_lead = False
			lead_part = fj.sorted_by_E(j.constituents())[0]
			pypart = pythiafjext.getPythia8Particle(lead_part)
			if args.debug:
				pinfo('leading id is', pypart.id(), pypart.name(), 'jet', j)
			if abs(pypart.id()) == 421:
				# pinfo('leading D0')
				isD0_lead = True
			l = lund_gen.result(j)
			if len(l) > 0:
				tw.fill_branch('Epair', [s.pair().e() for s in l])
				tw.fill_branch('z', [s.z() for s in l])
				tw.fill_branch('kt', [s.kt() for s in l])
				tw.fill_branch('delta', [s.Delta() for s in l])
				tw.fill_branch('D0lead', isD0_lead)
				tw.fill_branch('lead_id', pypart.id())
				tw.fill_tree()
			else:
				if args.debug:
					pwarning("len of a lund is less than 1?", len(l), l)

	pythia.stat()
	outf.Write()
	outf.Close()
	print('[i] written', outf.GetName())


if __name__ == '__main__':
	main()
