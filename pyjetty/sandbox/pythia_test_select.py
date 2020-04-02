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

def print_pythia_particle(p):
	print(' ', p.name(), p.id())

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--charged', default=False, action='store_true')
	parser.add_argument('--neutral', default=False, action='store_true')
	parser.add_argument('--partons', default=False, action='store_true')
	parser.add_argument('--quarks', default=False, action='store_true')
	parser.add_argument('--diquarks', default=False, action='store_true')
	parser.add_argument('--gluons', default=False, action='store_true')
	parser.add_argument('--leptons', default=False, action='store_true')
	parser.add_argument('--Ncharged', default=False, action='store_true')
	parser.add_argument('--Nneutral', default=False, action='store_true')
	parser.add_argument('--Npartons', default=False, action='store_true')
	parser.add_argument('--Nquarks', default=False, action='store_true')
	parser.add_argument('--Ndiquarks', default=False, action='store_true')
	parser.add_argument('--Ngluons', default=False, action='store_true')
	parser.add_argument('--Nleptons', default=False, action='store_true')
	args = parser.parse_args()

	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 1:
		args.nev = 1

	part_selection = [pythiafjext.kFinal]
	if args.charged:
		part_selection.append(pythiafjext.kCharged)
	if args.neutral:
		part_selection.append(pythiafjext.kNeutral)
	if args.partons:
		part_selection.append(pythiafjext.kParton)
	if args.quarks:
		part_selection.append(pythiafjext.kQuark)
	if args.diquarks:
		part_selection.append(pythiafjext.kDiquark)
	if args.gluons:
		part_selection.append(pythiafjext.kGluon)
	if args.leptons:
		part_selection.append(pythiafjext.kLepton)

	if args.Ncharged:
		part_selection.append(-pythiafjext.kCharged)
	if args.Nneutral:
		part_selection.append(-pythiafjext.kNeutral)
	if args.Npartons:
		part_selection.append(-pythiafjext.kParton)
	if args.Nquarks:
		part_selection.append(-pythiafjext.kQuark)
	if args.Ndiquarks:
		part_selection.append(-pythiafjext.kDiquark)
	if args.Ngluons:
		part_selection.append(-pythiafjext.kGluon)
	if args.Nleptons:
		part_selection.append(-pythiafjext.kLepton)

	for iev in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		parts = []
		parts = pythiafjext.vectorize_select(pythia, part_selection, 0, True)
		print('[i] number of particles selected:', len(parts))
		for p in parts:
			pypart = pythiafjext.getPythia8Particle(p)
			print_pythia_particle(pypart)

if __name__ == '__main__':
	main()