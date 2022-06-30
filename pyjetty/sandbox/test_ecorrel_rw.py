#!/usr/bin/env python3

from pyjetty.mputils.mputils import logbins
import operator as op
import itertools as it
import sys
import os
import argparse
from tqdm import tqdm
from heppy.pythiautils import configuration as pyconf
import pythiaext
import pythiafjext
import pythia8
import fjtools
import ecorrel
import fjcontrib
import fjext
import fastjet as fj
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from tqdm import tqdm
import argparse
import os


def get_args_from_settings(ssettings):
	sys.argv = sys.argv + ssettings.split()
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly')
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--output', default="test_ecorrel_rw.root", type=str)
	parser.add_argument('--user-seed', help='pythia seed',
											default=1111, type=int)
	args = parser.parse_args()
	return args


def main():
	mycfg = []
	ssettings = "--py-ecm 7000 --py-pthatmin 500"
	args = get_args_from_settings(ssettings)
	pythia_hard = pyconf.create_and_init_pythia_from_args(args, mycfg)
	jet_R0 = 0.5
	max_eta_jet = 1.9 # open CMS data
	max_eta_hadron = max_eta_jet + jet_R0 * 2.
	parts_selector_h = fj.SelectorAbsEtaMax(max_eta_hadron)
	# jet_selector = fj.SelectorPtMin(500.0) & fj.SelectorPtMax(550.0) & fj.SelectorAbsEtaMax(max_eta_hadron - 1.05 * jet_R0)
	jet_selector = fj.SelectorPtMin(500.0) & fj.SelectorPtMax(550.0) & fj.SelectorAbsEtaMax(max_eta_jet) 
	pfc_selector0 = fj.SelectorPtMin(0.)
	pfc_selector1 = fj.SelectorPtMin(1.)

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	print(jet_def)

	nbins = int(9.)
	lbins = logbins(1.e-3, 1., nbins)
	fout = ROOT.TFile(args.output, 'recreate')
	fout.cd()

	heec_0 = ROOT.TH1F('heec_0', 'heec_0', nbins, lbins)
	he3c_0 = ROOT.TH1F('he3c_0', 'he3c_0', nbins, lbins)
	he4c_0 = ROOT.TH1F('he4c_0', 'he4c_0', nbins, lbins)

	heec_1 = ROOT.TH1F('heec_1', 'heec_1', nbins, lbins)
	he3c_1 = ROOT.TH1F('he3c_1', 'he3c_1', nbins, lbins)
	he4c_1 = ROOT.TH1F('he4c_1', 'he4c_1', nbins, lbins)

	# heec_s0 = ROOT.TH1F('heec_s0', 'heec_s0', nbins, lbins)
	# he3c_s0 = ROOT.TH1F('he3c_s0', 'he3c_s0', nbins, lbins)
	# he4c_s0 = ROOT.TH1F('he4c_s0', 'he4c_s0', nbins, lbins)
	# heec_s1 = ROOT.TH1F('heec_s1', 'heec_s1', nbins, lbins)
	# he3c_s1 = ROOT.TH1F('he3c_s1', 'he3c_s1', nbins, lbins)
	# he4c_s1 = ROOT.TH1F('he4c_s1', 'he4c_s1', nbins, lbins)

	hjpt = ROOT.TH1F('hjpt', 'hjpt', 10, 500, 550)
 
	for n in tqdm(range(args.nev)):
		if not pythia_hard.next():
				continue
		# parts_pythia_h = pythiafjext.vectorize_select(pythia_hard, [pythiafjext.kFinal, pythiafjext.kCharged], 0, True)
		parts_pythia_h = pythiafjext.vectorize_select(pythia_hard, [pythiafjext.kFinal], 0, True)
		parts_pythia_h_selected = parts_selector_h(parts_pythia_h)

		jets_h = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h_selected))) 
		if len(jets_h) < 1:
				continue

		# print('njets:', hjpt.Integral())
		
		for j in jets_h:
			hjpt.Fill(j.perp())
			# note: the EEC,... takes vector<PseudoJet> while PseudoJet::constituents() returns a tuple in python
			# so we use a helper function (in SWIG only basic types typles handled easily...)
			# vconstits = ecorrel.constituents_as_vector(j)
			# eecs_alt = ecorrel.EEC(vconstits, scale=j.perp())   
			# e3cs_alt = ecorrel.E3C(vconstits, scale=j.perp())
			# e4cs_alt = ecorrel.E4C(vconstits, scale=j.perp())

   			# alternative: push constutents to a vector in python
			_v = fj.vectorPJ()
			_ = [_v.push_back(c) for c in j.constituents()]
			# select only charged constituents
			_vc = fj.vectorPJ()
			_ = [_vc.push_back(c) for c in j.constituents()
                            if pythiafjext.getPythia8Particle(c).isCharged()]
			# select only charged constituents with 1 GeV cut
			_vc1 = fj.vectorPJ()
			_ = [_vc1.push_back(c) for c in pfc_selector1(j.constituents())
                            if pythiafjext.getPythia8Particle(c).isCharged()]

			eecs_cont2 = ecorrel.EEC(_v, j.perp())
			eecs_cont3 = ecorrel.E3C(_v, j.perp())
			eecs_cont4 = ecorrel.E4C(_v, j.perp())
   
			_ = [heec_0.Fill(eecs_cont2.rs()[i], eecs_cont2.weights()[i]) for i in range(0, eecs_cont2.rs().size())]
			_ = [he3c_0.Fill(eecs_cont3.rs()[i], eecs_cont3.weights()[i]) for i in range(0, eecs_cont3.rs().size())]
			_ = [he4c_0.Fill(eecs_cont4.rs()[i], eecs_cont4.weights()[i]) for i in range(0, eecs_cont4.rs().size())]

			_v = fj.vectorPJ()
			_ = [_v.push_back(c) for c in pfc_selector1(j.constituents()) if pythiafjext.getPythia8Particle(c).isCharged()]

			eecs_cont2 = ecorrel.EEC(_v, j.perp())
			eecs_cont3 = ecorrel.E3C(_v, j.perp())
			eecs_cont4 = ecorrel.E4C(_v, j.perp())

			_ = [heec_1.Fill(eecs_cont2.rs()[i], eecs_cont2.weights()[i]) for i in range(0, eecs_cont2.rs().size())]
			_ = [he3c_1.Fill(eecs_cont3.rs()[i], eecs_cont3.weights()[i]) for i in range(0, eecs_cont3.rs().size())]
			_ = [he4c_1.Fill(eecs_cont4.rs()[i], eecs_cont4.weights()[i]) for i in range(0, eecs_cont4.rs().size())]

	njets = hjpt.Integral()
	if njets == 0:
		njets = 1.

	fout.cd()
	for h in [heec_0, he3c_0, he4c_0, heec_1, he3c_1, he4c_1]:
		h.Sumw2()
		h.Scale(1./h.Integral())

	for h in [he3c_0, he4c_0]:
		fout.cd()
		hc = h.Clone(h.GetName() + '_ratio_EEC')
		hc.Sumw2()
		hc.Divide(heec_0)  
		hc.Write()

	for h in [he3c_1, he4c_1]:
		fout.cd()
		hc = h.Clone(h.GetName() + '_ratio_EEC')
		hc.Sumw2()
		hc.Divide(heec_1)
		hc.Write()

	pythia_hard.stat()

	fout.Write()
	fout.Close()
	print('[i] written ', fout.GetName())


if __name__ == '__main__':
	main()
