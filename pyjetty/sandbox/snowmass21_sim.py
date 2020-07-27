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

from pyjetty.mputils import logbins, linbins

import ROOT
ROOT.gROOT.SetBatch(True)


def part_int_h2(h2):
	h2i = h2.Clone(h2.GetName() + '_pint')
	h2i.Reset()
	for ipt in range(1, h2.GetNbinsX()+1):
		for ir in range(1, h2.GetNbinsY()+1):
			hint = h2.Integral(ipt, ipt, 1, ir, 'width')
			# hint = h2.Integral(ipt, ipt, 1, ir)
			h2i.SetBinContent(ipt, ir, hint)
			h2i.SetBinError(ipt, ir, ROOT.TMath.Sqrt(hint))
		hptint = h2.Integral(ipt, ipt, 1, h2.GetNbinsY()+1, 'width')
		for ir in range(1, h2.GetNbinsY()+1):
			nval = h2i.GetBinContent(ipt, ir) / hptint
			h2i.SetBinContent(ipt, ir, nval)
			h2i.SetBinError(ipt, ir, 0)
	return h2i

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--ignore-mycfg', help="ignore some settings hardcoded here", default=True, action='store_true')
	parser.add_argument('--part-info', help="attach particle info to psj", default=False, action='store_true')
	parser.add_argument('--output', help="output file name", default='snowmass21_sim.root', type=str)

	args = parser.parse_args()

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R04 = 0.4
	jet_def_R04 = fj.JetDefinition(fj.antikt_algorithm, jet_R04)
	jet_R10 = 1.0
	jet_def_R10 = fj.JetDefinition(fj.antikt_algorithm, jet_R10)
	# select jets
	jet_selector = fj.SelectorPtMin(10.0) & fj.SelectorAbsEtaMax(2.)

	# pythia init
	mycfg = ['PhaseSpace:pThatMin = 100']
	if args.ignore_mycfg:
		mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return

	fout = ROOT.TFile(args.output, 'recreate')
	fout.cd()
	ptbins = logbins(10, 1000, 25)
	irbins = linbins(0, 1, 25)
	nsbins = linbins(0, 25, 25)
	print(ptbins)
	print(irbins)
	hERptR04 = ROOT.TH2D('hERptR04', 'hERptR04;p_{T}^{jet} (GeV/c); r', 25, ptbins, 25, irbins)
	hERptR10 = ROOT.TH2D('hERptR10', 'hERptR10;p_{T}^{jet} (GeV/c); r', 25, ptbins, 25, irbins)
	hR04nsd01pt = ROOT.TH2D('hR04nsd01pt', 'hR04nsd01pt', 25, ptbins, 25, nsbins)
	hR04nsd02pt = ROOT.TH2D('hR04nsd02pt', 'hR04nsd02pt', 25, ptbins, 25, nsbins)
	hR04nlundpt = ROOT.TH2D('hR04nlundpt', 'hR04nlundpt', 25, ptbins, 25, nsbins)
	hR04sd02Rg = ROOT.TH2D('hR04sd02Rg', 'hR04sd02Rg', 25, ptbins, 25, irbins)
	hR04sd02Rg_n = ROOT.TH2D('hR04sd02Rg_n', 'hR04sd02Rg_n', 25, ptbins, 25, irbins)

	hR10nsd01pt = ROOT.TH2D('hR10nsd01pt', 'hR10nsd01pt', 25, ptbins, 25, nsbins)
	hR10nsd02pt = ROOT.TH2D('hR10nsd02pt', 'hR10nsd02pt', 25, ptbins, 25, nsbins)
	hR10nlundpt = ROOT.TH2D('hR10nlundpt', 'hR10nlundpt', 25, ptbins, 25, nsbins)
	hR10sd02Rg = ROOT.TH2D('hR10sd02Rg', 'hR10sd02Rg', 25, ptbins, 25, irbins)
	hR10sd02Rg_n = ROOT.TH2D('hR10sd02Rg_n', 'hR10sd02Rg_n', 25, ptbins, 25, irbins)

	# event loop
	if args.nev < 10:
		args.nev = 100
	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		attach_pythia_particle_info = args.part_info
		parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], attach_pythia_particle_info)
		jets_R04 = jet_selector(jet_def_R04(parts))
		jets_R10 = jet_selector(jet_def_R10(parts))

		gshops = [fjcontrib.GroomerShop(j, 0.4, fj.cambridge_algorithm) for j in jets_R04]
		for ij, jj in enumerate(jets_R04):
			for c in jj.constituents():
				hERptR04.Fill(jj.perp(), c.delta_R(jj), c.perp())
			lund_splits = gshops[ij].lund_splits()
			n_SD_01 = len([s.z() for s in lund_splits if s.z() > 0.1])
			hR04nsd01pt.Fill(jj.perp(), n_SD_01)
			n_SD_02 = len([s.z() for s in lund_splits if s.z() > 0.2])
			hR04nsd02pt.Fill(jj.perp(), n_SD_02)
			n_splits = len(lund_splits)
			hR04nlundpt.Fill(jj.perp(), n_splits)
			sd02 = gshops[ij].soft_drop(0.0, 0.2, 1.)
			hR04sd02Rg.Fill(jj.perp(), sd02.Delta())
			[hR04sd02Rg_n.Fill(jj.perp(), s.Delta()) for s in lund_splits if s.z() > 0.2]

		gshops = [fjcontrib.GroomerShop(j, 1.0, fj.cambridge_algorithm) for j in jets_R10]
		for ij, jj in enumerate(jets_R10):
			for c in jj.constituents():
				hERptR10.Fill(jj.perp(), c.delta_R(jj), c.perp())
			lund_splits = gshops[ij].lund_splits()
			n_SD_01 = len([s.z() for s in lund_splits if s.z() > 0.1])
			hR10nsd01pt.Fill(jj.perp(), n_SD_01)
			n_SD_02 = len([s.z() for s in lund_splits if s.z() > 0.2])
			hR10nsd02pt.Fill(jj.perp(), n_SD_02)
			n_splits = len(lund_splits)
			hR10nlundpt.Fill(jj.perp(), n_splits)
			sd02 = gshops[ij].soft_drop(0.0, 0.2, 1.)
			hR10sd02Rg.Fill(jj.perp(), sd02.Delta())
			[hR10sd02Rg_n.Fill(jj.perp(), s.Delta()) for s in lund_splits if s.z() > 0.2]

	pythia.stat()

	fout.cd()
	for h in [hERptR04, hERptR10]:
		hi = part_int_h2(h)
		hi.Write()

	for h in [hERptR04, hR04nsd01pt, hR04nsd02pt, hR04nlundpt, hR04sd02Rg, hR04sd02Rg_n]:
		tp = h.ProfileX(h.GetName() + '_pfx', 1, -1, 's')
		tp.Write()

	for h in [hERptR10, hR10nsd01pt, hR10nsd02pt, hR10nlundpt, hR10sd02Rg, hR10sd02Rg_n]:
		tp = h.ProfileX(h.GetName() + '_pfx', 1, -1, 's')
		tp.Write()

	fout.Write()

if __name__ == '__main__':
	main()
