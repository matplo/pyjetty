#!/usr/bin/env python

import fastjet as fj
import fjcontrib
import fjext
import fjtools

import tqdm
import argparse
import os
import numpy as np

import ROOT
from pyjetty.mputils import treewriter

from pyjetty.mputils import BoltzmannEvent

import pythia8
import pythiafjext
import pythiaext
from heppy.pythiautils import configuration as pyconf


def calc_n_lead(j, n=8):
	cs = fj.sorted_by_pt(j.constituents())
	_sum_top_n = 0
	_sum_all = 0
	for i in range(len(cs)):
		_sum_all += cs[i].pt()
		if i < n:
			_sum_top_n += cs[i].pt()
	_fraction_pt = _sum_top_n / _sum_all
	return _sum_all, _sum_top_n, _fraction_pt

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
		be = BoltzmannEvent(mean_pt=0.6, multiplicity=2000 * max_eta * 2, max_eta=max_eta, max_pt=100)
		print(be)

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.6
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_def_emb = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(10.0) & fj.SelectorAbsEtaMax(1)

	fout = ROOT.TFile('pythia_8jet.root', 'recreate')
	fout.cd()
	hfraction = ROOT.TProfile('hfraction', 'hfraction', 10, 0, 100)
	hdpt = ROOT.TProfile('hdpt', 'hdpt', 10, 0, 100)
	hfraction2D = ROOT.TH2F('hfraction2D', 'hfraction', 10, 0, 100, 20, 0, 1)
	hdpt2D = ROOT.TH2F('hdpt2D', 'hdpt', 10, 0, 100, 20, -1, 0)

	hfraction_emb = ROOT.TProfile('hfraction_emb', 'hfraction', 10, 0, 100)
	hdpt_emb = ROOT.TProfile('hdpt_emb', 'hdpt', 10, 0, 100)
	hfraction2D_emb = ROOT.TH2F('hfraction2D_emb', 'hfraction', 10, 0, 100, 20, 0, 1)
	hdpt2D_emb = ROOT.TH2F('hdpt2D_emb', 'hdpt', 10, 0, 100, 20, -1, 0)


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
			_sum_all, _sum_top_n, _fraction_pt	= calc_n_lead(j, 8)
			hfraction.Fill(j.pt(), _fraction_pt)
			hdpt.Fill(j.pt(), (_sum_top_n - _sum_all) / _sum_all)
			hfraction2D.Fill(j.pt(), _fraction_pt)
			hdpt2D.Fill(j.pt(), (_sum_top_n - _sum_all) / _sum_all)

			full_event = fj.vectorPJ()
			tmp = [full_event.push_back(psj) for psj in bg_parts]
			tmp = [full_event.push_back(psj) for psj in j.constituents()]
			embd_jets = fj.sorted_by_pt(jet_selector(jet_def_emb(full_event)))
			for jemb in embd_jets:
				mpt = fjtools.matched_pt(jemb, j)
				if mpt < 0.5:
					continue
				_sum_all_emb, _sum_top_n_emb, _fraction_pt_emb	= calc_n_lead(jemb, 8)
				_fraction_pt_emb = _sum_top_n_emb / _sum_all;
				hfraction_emb.Fill(jemb.pt(), _fraction_pt)
				hdpt_emb.Fill(jemb.pt(), (_sum_top_n_emb - _sum_all) / _sum_all)
				hfraction2D_emb.Fill(jemb.pt(), _fraction_pt)
				hdpt2D_emb.Fill(jemb.pt(), (_sum_top_n_emb - _sum_all) / _sum_all)

	fg = ROOT.TF1('fg', 'gaus', 0, 1)
	fg.SetParameter(0, 1)
	fg.SetParameter(1, 0.8)
	fg.SetParameter(2, 0.1)
	hfraction2D.FitSlicesY(fg)
	hfraction2D_emb.FitSlicesY(fg)

	fgdpt = ROOT.TF1('fgdpt', 'gaus', -1, 0)
	fgdpt.SetParameter(0, 1)
	fgdpt.SetParameter(1, -0.2)
	fgdpt.SetParameter(2, 0.1)
	hdpt2D.FitSlicesY(fgdpt)
	hdpt2D_emb.FitSlicesY(fgdpt)

	fout.Write()
	fout.Close()


if __name__ == '__main__':
	main()