#!/usr/bin/env python

from __future__ import print_function

import fastjet
import fjcontrib
import fjext

import tqdm
import argparse
import os

import pythia8
import pythiaext
# pythia8+hepmc3 writing 
# import pythiahepmc3

from heppy.pythiautils import configuration as pyconf


def main():
	parser = argparse.ArgumentParser(description='pythia8 in python', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()	

	hepmc2output = "pythia_gen_test_hepmc2.dat"
	pyhepmc2writer = pythiaext.Pythia8HepMC2Wrapper(hepmc2output)

	# experimental hepmc3
	# hepmc3output = "pyhia_gen_test_hepmc3.dat"
	# pyhepmc3writer = pythiahepmc3.Pythia8HepMC3Wrapper(hepmc3output)

	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 10:
		args.nev = 10
	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		pyhepmc2writer.fillEvent(pythia)
		# pyhepmc3writer.fillEvent(pythia)
	pythia.stat()
	pythia.settings.writeFile(args.py_cmnd_out)
	print("[i] file written: {}".format(hepmc2output))

if __name__ == '__main__':
	main()
