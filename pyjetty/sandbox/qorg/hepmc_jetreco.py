#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os
import numpy as np

import pyhepmc_ng


def find_jets_hepmc(jet_def, jet_selector, hepmc_event):
	fjparts = []
	for i,p in enumerate(hepmc_event.particles):
		if p.status == 1 and not p.end_vertex:
			psj = fj.PseudoJet(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e)
			# psj.set_user_index(i)
			fjparts.append(psj)
	jets = jet_selector(jet_def(fjparts))
	return jets


def main():
	parser = argparse.ArgumentParser(description='pythia8 in python', prog=os.path.basename(__file__))
	parser.add_argument('-i', '--input', help='input file', default='low', type=str, required=True)
	parser.add_argument('--hepmc', help='what format 2 or 3', default=2, type=int)
	parser.add_argument('--nev', help='number of events', default=10, type=int)
	args = parser.parse_args()	

	###
	# now lets read the HEPMC file and do some jet finding
	if args.hepmc == 3:
		input_hepmc = pyhepmc_ng.ReaderAscii(args.input)
	if args.hepmc == 2:
		input_hepmc = pyhepmc_ng.ReaderAsciiHepMC2(args.input)

	if input_hepmc.failed():
		print ("[error] unable to read from {}".format(args.input))
		sys.exit(1)

	# jet finder
	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(100.0) # & fj.SelectorPtMax(200.0) & fj.SelectorAbsEtaMax(3)

	all_jets = []
	event_hepmc = pyhepmc_ng.GenEvent()
	pbar = tqdm.tqdm(range(args.nev))
	while not input_hepmc.failed():
		ev = input_hepmc.read_event(event_hepmc)
		if input_hepmc.failed():
			break
		jets_hepmc = find_jets_hepmc(jet_def, jet_selector, event_hepmc)
		all_jets.extend(jets_hepmc)
		pbar.update()
		if pbar.n >= args.nev:
			break

	jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
	lund_gen = fjcontrib.LundGenerator(jet_def_lund)

	print('[i] making lund diagram for all jets...')
	lunds = [lund_gen.result(j) for j in all_jets]

	print('[i] listing lund plane points... Delta, kt - for {} selected jets'.format(len(all_jets)))
	for l in lunds:
		print ('- jet pT={0:5.2f} eta={1:5.2f}'.format(l[0].pair().perp(), l[0].pair().eta()))
		print ('  Deltas={}'.format([s.Delta() for s in l]))
		print ('  kts={}'.format([s.Delta() for s in l]))
		print ( )

	print('[i] reclustering and using soft drop...')
	jet_def_rc = fj.JetDefinition(fj.cambridge_algorithm, 0.1)
	print('[i] Reclustering:', jet_def_rc)

	rc = fjcontrib.Recluster(jet_def_rc, True)
	sd = fjcontrib.SoftDrop(0, 0.1, 1.0)
	for i,j in enumerate(all_jets):
		j_rc = rc.result(j)
		print()
		print('- [{0:3d}] orig pT={1:10.3f} reclustered pT={2:10.3f}'.format(i, j.perp(), j_rc.perp()))
		j_sd = sd.result(j)
		print('  |-> after soft drop pT={0:10.3f} delta={1:10.3f}'.format(j_sd.perp(), j_sd.perp() - j.perp()))
		sd_info = fjcontrib.get_SD_jet_info(j_sd)
		print("  |-> SD jet params z={0:10.3f} dR={1:10.3f} mu={2:10.3f}".format(sd_info.z, sd_info.dR, sd_info.mu))

if __name__ == '__main__':
	main()
