#!/usr/bin/env python

from __future__ import print_function

# oder matters - for fastjet to play nice with pythia load fj first
import fastjet
import fjcontrib
import fjext

import tqdm
import argparse
import os

import pythia8
import pythiaext

from heppy.pythiautils import configuration as pyconf


def main():
	parser = argparse.ArgumentParser(description='pythia8 in python', prog=os.path.basename(__file__))
	parser.add_argument('-o', '--output', help='output file name', default='pythia_hepmc.dat', type=str)
	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()	

	hepmc2output = args.output
	pyhepmc2writer = pythiaext.Pythia8HepMC2Wrapper(hepmc2output)

	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 10:
		args.nev = 10
	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		pyhepmc2writer.fillEvent(pythia)
	pythia.stat()

	print("[i] file written: {}".format(hepmc2output))

if __name__ == '__main__':
	main()
