#!/usr/bin/env python3

from __future__ import print_function
import tqdm
import argparse
import os
import numpy as np
import sys

import pyhepmc
import particle
# import ROOT if needed
# import ROOT
import math
import array

import hepmc_count_events

import fastjet as fj

def logbins(xmin, xmax, nbins):
        lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
        arr = array.array('f', lspace)
        return arr

def main():
	parser = argparse.ArgumentParser(description='read hepmc and analyze eecs', prog=os.path.basename(__file__))
	parser.add_argument('-i', '--input', help='input file', default='low', type=str, required=True)
	parser.add_argument('--hepmc', help='what format 2 or 3', default=2, type=int)
	parser.add_argument('--nev', help='number of events', default=-1, type=int)
	parser.add_argument('-o','--output', help='root output filename', default='test_output.root', type=str)
	args = parser.parse_args()	


	if args.nev < 0:
		args.nev = hepmc_count_events.get_n_per_file(args.input, args.hepmc)[0]

	###
	# now lets read the HEPMC file and do some jet finding
	if args.hepmc == 3:
		input_hepmc = pyhepmc.io.ReaderAscii(args.input)
	if args.hepmc == 2:
		input_hepmc = pyhepmc.io.ReaderAsciiHepMC2(args.input)

	if input_hepmc.failed():
		print ("[error] unable to read from {}".format(args.input))
		sys.exit(1)

	# jet finder
	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	particle_eta_max = 0.9
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
 
	# change the selectors as needed
	jet_selector = fj.SelectorPtMin(20.0)
	jet_selector = fj.SelectorPtMin(20.0) * fj.SelectorPtMax(500.0) * fj.SelectorAbsEtaMax(particle_eta_max - jet_R0 * 1.05)

	event_hepmc = pyhepmc.GenEvent()
	pbar = tqdm.tqdm(range(args.nev))
	njets = 0
	while not input_hepmc.failed():
		ev = input_hepmc.read_event(event_hepmc)
		if input_hepmc.failed():
			break
		fjparts = fj.vectorPJ()
		for i,p in enumerate(event_hepmc.particles):
			if p.status == 1 and not p.end_vertex:
				if particle.Particle.from_pdgid(p.pid).charge == 0:
					continue
				psj = fj.PseudoJet(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e)
				psj.set_user_index(i)
				fjparts.push_back(psj)

		jets = fj.sorted_by_pt(jet_selector(jet_def(fjparts)))

		njets += len(jets)
		for j in jets:
			for c in j.constituents():
				# your code here
				pass

		pbar.update()
		if pbar.n >= args.nev:
			break

	print('[i] number of jets:', njets)

if __name__ == '__main__':
	main()
