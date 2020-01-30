#!/usr/bin/env python

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os
import numpy as np

import ROOT
from pyjetty.mputils import treewriter

import pythia8
import pythiafjext
import pythiaext
from heppy.pythiautils import configuration as pyconf


def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()

	if args.nev < 1:
		args.nev = 1

	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	part_selection = [pythiafjext.kFinal]

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.6
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(10.0) & fj.SelectorAbsEtaMax(1)

	fout = ROOT.TFile('pythia_8jet.root', 'recreate')
	fout.cd()
	hfraction = ROOT.TProfile('hfraction', 'hfraction', 10, 0, 100)
	hdpt = ROOT.TProfile('hdpt', 'hdpt', 10, 0, 100)

	hfraction2D = ROOT.TH2F('hfraction2D', 'hfraction', 10, 0, 100, 20, 0, 1)
	hdpt2D = ROOT.TH2F('hdpt2D', 'hdpt', 10, 0, 100, 20, -1, 0)

	for iev in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		parts = []
		parts = pythiafjext.vectorize_select(pythia, part_selection, 0, True)
		jets = fj.sorted_by_pt(jet_selector(jet_def(parts)))
		for j in jets:
			cs = fj.sorted_by_pt(j.constituents())
			_sum_top_n = 0
			_sum_all = 0
			for i in range(len(cs)):
				_sum_all += cs[i].pt()
				if i < 8:
					_sum_top_n += cs[i].pt()
			_fraction_pt = _sum_top_n / _sum_all
			hfraction.Fill(j.pt(), _fraction_pt)
			hdpt.Fill(j.pt(), (_sum_top_n - _sum_all) / _sum_all)
			hfraction2D.Fill(j.pt(), _fraction_pt)
			hdpt2D.Fill(j.pt(), (_sum_top_n - _sum_all) / _sum_all)

	fg = ROOT.TF1('fg', 'gaus', 0, 1)
	fg.SetParameter(0, 1)
	fg.SetParameter(1, 0.8)
	fg.SetParameter(2, 0.1)
	hfraction2D.FitSlicesY(fg)

	fgdpt = ROOT.TF1('fgdpt', 'gaus', -1, 0)
	fgdpt.SetParameter(0, 1)
	fgdpt.SetParameter(1, -0.2)
	fgdpt.SetParameter(2, 0.1)
	hdpt2D.FitSlicesY(fgdpt)

	fout.Write()
	fout.Close()


if __name__ == '__main__':
	main()