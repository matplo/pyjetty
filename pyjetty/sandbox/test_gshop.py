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

def groom(jet, jetR, zcut, beta):
	gshop = fjcontrib.GroomerShop(jet)
	return gshop.soft_drop(beta, zcut, jetR)

def groom_copy(jet, jetR, zcut, beta):
	gshop = fjcontrib.GroomerShop(jet)
	return gshop.soft_drop(beta, zcut, jetR).copy()

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--ignore-mycfg', help="ignore some settings hardcoded here", default=False, action='store_true')

	args = parser.parse_args()

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(100.0) & fj.SelectorAbsEtaMax(1)
	print(jet_def)

	all_jets = []

	mycfg = ['PhaseSpace:pThatMin = 100']
	if args.ignore_mycfg:
		mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return
	if args.nev < 10:
		args.nev = 10
	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		attach_pythia_particle_info = True
		parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], attach_pythia_particle_info)
		jets = jet_selector(jet_def(parts))
		if len(jets) < 2:
			continue
		gsetting = {'sd': [0.1, 0]}
		g1 = groom(jets[0], 0.4, 0.1, 0)
		print('  |-> GroomerShop::soft_drop 1 no copy', g1.as_string())
		g2 = groom(jets[1], 0.4, 0.1, 0)
		print('  |-> GroomerShop::soft_drop 2 no copy', g2.as_string())

		g1_copy = groom_copy(jets[0], 0.4, 0.1, 0)
		g2_copy = groom_copy(jets[1], 0.4, 0.1, 0)
		print('  |-> GroomerShop::soft_drop 1 w/ copy', g1_copy.as_string())
		print('  |-> GroomerShop::soft_drop 1 no copy', g1.as_string(), ' * this can be a hard to debug runtime error *')
		print('  |-> GroomerShop::soft_drop 2 w/ copy', g2_copy.as_string())
		print('  |-> GroomerShop::soft_drop 2 no copy', g2.as_string())

		gshops = [fjcontrib.GroomerShop(j) for j in jets]
		for ij, jj in enumerate(jets):
			print('orig jet:', jj, 'from gshop:', gshops[ij].jet(), 'SD:', gshops[ij].soft_drop(0.0, 0.1, 0.4).as_string())
		# print (' xxx ')
		# fjcontrib.setGroomer(jets[0], 0.4)
		# fjcontrib.setGroomer(jets[0], 0.4)
		# fjcontrib.setGroomer(jets[0], 0.4)
		# print (' xxx ', fj.cambridge_algorithm)

		# algo = fj.cambridge_algorithm
		# this is removed from the library - GroomerShopUI 
		# [print(' |-> Groomed SD', fjcontrib.groom(j).soft_drop(0.0, 0.1, 0.4).as_string()) for j in jets]
		# [print(' |-> Groomed max-tf', fjcontrib.groom(j).max_tf().as_string()) for j in jets]

		print('---')

	pythia.stat()

if __name__ == '__main__':
	main()
