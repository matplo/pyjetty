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

import csevent

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
	jet_selector = fj.SelectorPtMin(100.0) & fj.SelectorAbsEtaMax(1 - 1.05 * jet_R0)
	print(jet_def)

	cs = csevent.CSEventSubtractor(alpha=0, max_distance=0.25, max_eta=1, bge_rho_grid_size=0.1)

	mycfg = ['PhaseSpace:pThatMin = 100', 'PhaseSpace:pThatMax = 110']
	if args.ignore_mycfg:
		mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return
	if args.nev < 100:
		args.nev = 100
	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
		signal_jets = jet_selector(jet_def(parts))
		cs_parts = cs.process_event(parts)
		cs_signal_jets = jet_selector(jet_def(cs_parts))

		sd = fjcontrib.SoftDrop(0, 0.1, 1.0)
		for i,j in enumerate(signal_jets):
			j_sd = sd.result(j)
			sd_info = fjcontrib.get_SD_jet_info(j_sd)
			# print("  |-> SD jet params z={0:10.3f} dR={1:10.3f} mu={2:10.3f}".format(sd_info.z, sd_info.dR, sd_info.mu))

	pythia.stat()



if __name__ == '__main__':
	main()
