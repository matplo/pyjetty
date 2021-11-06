#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
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

	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)

	max_eta_hadron=3
	jet_R0 = 0.4
	jet_selector = fj.SelectorPtMin(100.0) & fj.SelectorPtMax(125.0) & fj.SelectorAbsEtaMax(max_eta_hadron - 1.05 * jet_R0)
	parts_selector_h = fj.SelectorAbsEtaMax(max_eta_hadron)

	fj.ClusterSequence.print_banner()
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)

	if args.nev < 10:
		args.nev = 10
	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue


	pythia.stat()
	pythia.settings.writeFile(args.py_cmnd_out)


if __name__ == '__main__':
	main()
