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

def fill_branches(tw, j, dy_groomer, alphas=[], sds=[]):
	tw.fill_branch('j', j)
	for a in alphas:
		dy_groomed = dy_groomer.result(j, a)
		# if dy_groomed.pair().pt() > 0:
		# 	tw.fill_branch('dg_{:.1f}'.format(a), dy_groomed.harder())
		# 	tw.fill_branch('dg_{:.1f}'.format(a), dy_groomed.softer())
		tw.fill_branch('dg_{:.1f}'.format(a), dy_groomed)
	max_pt_groomed = dy_groomer.max_pt_softer(j)
	tw.fill_branch('max_ptsoft', max_pt_groomed)
	max_z_groomed = dy_groomer.max_z(j)
	tw.fill_branch('max_z', max_z_groomed)
	max_kt_groomed = dy_groomer.max_kt(j)
	tw.fill_branch('max_kt', max_kt_groomed)
	max_kappa_groomed = dy_groomer.max_kappa(j)
	tw.fill_branch('max_kappa', max_kappa_groomed)
	max_tf_groomed = dy_groomer.max_tf(j)
	tw.fill_branch('max_tf', max_tf_groomed)
	min_tf_groomed = dy_groomer.min_tf(j)
	tw.fill_branch('min_tf', min_tf_groomed)

	for i,sd in enumerate(sds):
		j_sd = sd.result(j)
		tw.fill_branch('sd{}'.format(i), j_sd)
		sd_info = fjcontrib.get_SD_jet_info(j_sd)
		tw.fill_branch('sd{}_z'.format(i), sd_info.z)
		tw.fill_branch('sd{}_Delta'.format(i), sd_info.dR)
		tw.fill_branch('sd{}_mu'.format(i), sd_info.mu)
		tw.fill_branch('sd{}_kt'.format(i), sd_info.z * j_sd.pt() * sd_info.dR)


def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	# could use --py-seed
	parser.add_argument('--fj-R', help='jet finder R', default=0.8, type=float)
	parser.add_argument('--user-seed', help='pythia seed', default=-1, type=int)
	parser.add_argument('--output', default='{}.root'.format(os.path.basename(__file__)), type=str)
	parser.add_argument('--min-jet-pt', help='jet pt selection', default=50., type=float)
	parser.add_argument('--max-jet-pt', help='jet pt selection', default=1000., type=float)
	args = parser.parse_args()

	if args.user_seed < 0:
		args.user_seed = -1
		mycfg = []
	else:
		pinfo('user seed for pythia', args.user_seed)
		# mycfg = ['PhaseSpace:pThatMin = 100']
		mycfg = ['Random:setSeed=on', 'Random:seed={}'.format(args.user_seed)]

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

	hepmc2output = '{}.hepmc2.dat'.format(args.output.replace('.root', ''))
	pyhepmc2writer = pythiaext.Pythia8HepMC2Wrapper(hepmc2output)

	outf = ROOT.TFile(args.output, 'recreate')
	outf.cd()
	t = ROOT.TTree('t', 't')
	tw = RTreeWriter(tree=t)
	tch = ROOT.TTree('tch', 'tch')
	twch = RTreeWriter(tree=tch)

	jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
	dy_groomer = fjcontrib.DynamicalGroomer(jet_def_lund)
	print (dy_groomer.description())

	sds = []
	sd01 = fjcontrib.SoftDrop(0, 0.1, 1.0)
	sd02 = fjcontrib.SoftDrop(0, 0.2, 1.0)
	sds.append(sd01)
	sds.append(sd02)

	# event loop
	for iev in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue

		parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, False)
		parts_pythia_h_selected = parts_selector_h(parts_pythia_h)
		jets_h = jet_selector(fj.sorted_by_pt(jet_def(parts_pythia_h)))

		[fill_branches(tw, j, dy_groomer, alphas=[0.1, 1.0, 2.0], sds=sds) for j in jets_h]
		if len(jets_h) > 0:
			tw.fill_tree()
			pyhepmc2writer.fillEvent(pythia)

		parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, False)
		parts_pythia_h_selected = parts_selector_h(parts_pythia_h)
		jets_h_ch = jet_selector(fj.sorted_by_pt(jet_def(parts_pythia_h)))

		[fill_branches(twch, j, dy_groomer, alphas=[0.1, 1.0, 2.0], sds=sds) for j in jets_h_ch]
		if len(jets_h_ch) > 0:
			twch.fill_tree()

	pythia.stat()
	outf.Write()
	outf.Close()

if __name__ == '__main__':
	main()
