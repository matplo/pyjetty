#!/usr/local/bin/python3

import argparse
import os
import tqdm

import fastjet as fj
import pyjetty
import pythia8
from recursivetools import pyrecursivetools as rt
from lundplane import pylundplane as lund

from pythiautils import configuration as pyconf

def main(args):
	if args.nev < 10:
		args.nev = 100
	print(args)
	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(100.0) & fj.SelectorPtMax(200.0) & fj.SelectorAbsEtaMax(1)
	sd = rt.SoftDrop(0, 0.1, 1.0)

	all_jets = []

	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		parts = pyjetty.vectorize(pythia, True, -1, 1, False)
		jets = jet_selector(jet_def(parts))
		all_jets.extend(jets)

	pythia.stat()

	jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
	lund_gen = lund.LundGenerator(jet_def_lund)

	print('making lund diagram for all jets...')
	lunds = [lund_gen.result(j) for j in all_jets]

	print('listing lund plane points... Delta, kt - for {} selected jets'.format(len(all_jets)))
	for l in lunds:
	    print ('- jet pT={0:5.2f} eta={1:5.2f}'.format(l[0].pair().perp(), l[0].pair().eta()))
	    print ('  Deltas={}'.format([s.Delta() for s in l]))
	    print ('  kts={}'.format([s.Delta() for s in l]))
	    print ( )

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='jet reco on alice data', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()	
	main(args)
