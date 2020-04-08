#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os
import numpy as np

from pyjetty.mputils import *

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

import ROOT

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	# could use --py-seed
	parser.add_argument('--fj-R', help='jet finder R', default=0.8, type=float)
	parser.add_argument('--user-seed', help='pythia seed', default=1111, type=int)
	parser.add_argument('--output', default="pythia_dyngroom_test_output.root", type=str)
	parser.add_argument('--min-jet-pt', help='jet pt selection', default=450., type=float)
	parser.add_argument('--max-jet-pt', help='jet pt selection', default=1000., type=float)
	args = parser.parse_args()

	if args.user_seed < 0:
		args.user_seed = 1111
	pinfo('user seed for pythia', args.user_seed)
	# mycfg = ['PhaseSpace:pThatMin = 100']
	mycfg = ['Random:setSeed=on', 'Random:seed={}'.format(args.user_seed)]
	mycfg.append('HadronLevel:all=off')
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 100:
		args.nev = 100

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = args.fj_R
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	print(jet_def)

	# hadron level - ALICE
	max_eta_hadron = 3.
	pwarning('max eta for particles after hadronization set to', max_eta_hadron)
	parts_selector_h = fj.SelectorAbsEtaMax(max_eta_hadron)
	jet_selector = fj.SelectorPtMin(args.min_jet_pt) & fj.SelectorPtMax(args.max_jet_pt) & fj.SelectorAbsEtaMax(max_eta_hadron - 1.05 * jet_R0)

	max_eta_parton = max_eta_hadron + 2. * jet_R0
	pwarning('max eta for partons set to', max_eta_parton)
	parts_selector_p = fj.SelectorAbsEtaMax(max_eta_parton)

	outf = ROOT.TFile(args.output, 'recreate')
	outf.cd()
	t = ROOT.TTree('t', 't')
	tw = RTreeWriter(tree=t)

	jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
	dy_groomer = fjcontrib.DynamicalGroomer(jet_def_lund)
	print (dy_groomer.description())

	# event loop
	for iev in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue

		# parts_pythia_p = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], add_particle_info = True)
		parts_pythia_p = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)
		parts_pythia_p_selected = parts_selector_p(parts_pythia_p)

		hstatus = pythia.forceHadronLevel()
		if not hstatus:
			pwarning('forceHadronLevel false event', iev)
			continue
		# parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kHadron, pythiafjext.kCharged])
		# parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], add_particle_info = True)
		parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)
		parts_pythia_h_selected = parts_selector_h(parts_pythia_h)

		jets_p = jet_selector(fj.sorted_by_pt(jet_def(parts_pythia_p)))
		jets_h = jet_selector(fj.sorted_by_pt(jet_def(parts_pythia_h)))

		for j in jets_p:
			tw.fill_branch('iev', iev)
			tw.fill_branch('parton_j', j)
			for a in [0.1, 1.0, 2.0]:
				dy_groomed = dy_groomer.result(j, a)
				if dy_groomed.pair().pt() > 0:
					tw.fill_branch('parton_j_logkt_dg{:.1f}'.format(a), ROOT.TMath.Log(dy_groomed.kt()))
					tw.fill_branch('parton_j_log1odR_dg{:.1f}'.format(a), ROOT.TMath.Log(1/dy_groomed.Delta()))
			max_pt_groomed = dy_groomer.max_pt_softer(j)
			if max_pt_groomed.pair().pt() > 0:
				tw.fill_branch('parton_j_logkt_max_pt_softer', ROOT.TMath.Log(max_pt_groomed.kt()))
				tw.fill_branch('parton_j_log1odR_max_pt_softer', ROOT.TMath.Log(1/max_pt_groomed.Delta()))
			max_z_groomed = dy_groomer.max_z(j)
			if max_z_groomed.pair().pt() > 0:
				tw.fill_branch('parton_j_logkt_max_z', ROOT.TMath.Log(max_z_groomed.kt()))
				tw.fill_branch('parton_j_log1odR_max_z', ROOT.TMath.Log(1/max_z_groomed.Delta()))
		for j in jets_h:
			tw.fill_branch('hadron_j', j)
			for a in [0.1, 1.0, 2.0]:
				dy_groomed = dy_groomer.result(j, a)
				if dy_groomed.pair().pt() > 0:
					tw.fill_branch('hadron_j_logkt_dg{:.1f}'.format(a), ROOT.TMath.Log(dy_groomed.kt()))
					tw.fill_branch('hadron_j_log1odR_dg{:.1f}'.format(a), ROOT.TMath.Log(1/dy_groomed.Delta()))
			max_pt_groomed = dy_groomer.max_pt_softer(j)
			if max_pt_groomed.pair().pt() > 0:
				tw.fill_branch('hadron_j_logkt_max_pt_softer', ROOT.TMath.Log(max_pt_groomed.kt()))
				tw.fill_branch('hadron_j_log1odR_max_pt_softer', ROOT.TMath.Log(1/max_pt_groomed.Delta()))
			max_z_groomed = dy_groomer.max_z(j)
			if max_z_groomed.pair().pt() > 0:
				tw.fill_branch('hadron_j_logkt_max_z', ROOT.TMath.Log(max_z_groomed.kt()))
				tw.fill_branch('hadron_j_log1odR_max_z', ROOT.TMath.Log(1/max_z_groomed.Delta()))
		if len(jets_p) > 0 or len(jets_h) > 0:
			tw.fill_tree()

	pythia.stat()
	outf.Write()
	outf.Close()

if __name__ == '__main__':
	main()
