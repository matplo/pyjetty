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

from pyjetty.mputils import MPBase, pwarning, pinfo, perror, treewriter, jet_analysis

def match_dR(j, partons, drmatch = 0.1):
	mps = [p for p in partons if j.delta_R(p) < drmatch]
	# for p in fj.sorted_by_pt(mps)[0]:
	if len(mps) < 1:
		return 0, -1, -1
	p = fj.sorted_by_pt(mps)[0]
	pyp = pythiafjext.getPythia8Particle(p)
	# print(p, pyp.id(), pyp.isQuark(), pyp.isGluon())
	return pyp.id(), pyp.isQuark(), pyp.isGluon()

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--ignore-mycfg', help="ignore some settings hardcoded here", default=False, action='store_true')
	parser.add_argument('--output', help="output file name", default="test_hjet_parton.root", type=str)
	parser.add_argument('--no-tt', help="do not require TT to accept the event", default=False, action='store_true')
	args = parser.parse_args()

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	hadron_eta_max = 2.0
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(10.0) & fj.SelectorPtMax(100.0) & fj.SelectorAbsEtaMax(hadron_eta_max - jet_R0)
	# jet_selector = fj.SelectorPtMin(40.0) & fj.SelectorPtMax(200.0) &fj.SelectorAbsEtaMax(hadron_eta_max - jet_R0)
	hTT6_selector 	= fj.SelectorPtMin(6) & fj.SelectorPtMax(7) & fj.SelectorAbsEtaMax(hadron_eta_max)
	hTT12_selector 	= fj.SelectorPtMin(12) & fj.SelectorPtMax(50) & fj.SelectorAbsEtaMax(hadron_eta_max)
	hTT20_selector 	= fj.SelectorPtMin(20) & fj.SelectorPtMax(50) & fj.SelectorAbsEtaMax(hadron_eta_max)

	print(jet_def)

	all_jets = []

	# mycfg = ['PhaseSpace:pThatMin = 80']
	# mycfg = ['PhaseSpace:pThatMin = 6']
	# mycfg = ['PhaseSpace:pThatMin = 12']
	# mycfg = ['PhaseSpace:pThatMin = 40']	
	mycfg = []
	if args.ignore_mycfg:
		mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		perror("pythia initialization failed.")
		return

	# tw = treewriter.RTreeWriter(name = 'lsjvsx', file_name = 'leadsj_vs_x.root')
	tw = treewriter.RTreeWriter(name = 'lsjvsx', file_name = args.output)

	zero_psj = fj.PseudoJet(0,0,10,10)

	if args.nev < 100:
		args.nev = 100
	t = tqdm.tqdm(total = args.nev)
	while t.n < args.nev:
		if not pythia.next():
			continue
		# parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
		partons = pythiafjext.vectorize_select(pythia, [pythiafjext.kParton], 0, True)
		# parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, False)
		parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, False)

		hTT6 = zero_psj
		hTT6s = fj.sorted_by_pt(hTT6_selector(parts))
		if len(hTT6s) > 0:
			hTT6 = hTT6s[0]

		hTT12 = zero_psj
		hTT12s = fj.sorted_by_pt(hTT12_selector(parts))
		if len(hTT12s) > 0:
			hTT12 = hTT12s[0]

		hTT20 = zero_psj
		hTT20s = fj.sorted_by_pt(hTT20_selector(parts))
		if len(hTT20s) > 0:
			hTT20 = hTT20s[0]

		if args.no_tt is False:
			if hTT12.perp() < 1 and hTT6.perp() < 1:
				continue

		jets = jet_selector(jet_def(parts))

		# for j in tqdm.tqdm(jets):
		for j in jets:
			t.update(1)
			j_type = match_dR(j, partons, jet_R0 / 2.)
			if j_type[0] is None:
				continue
			tw.fill_branches(			j=j,
                            hTT6=hTT6,
                            hTT12=hTT12,
                            hTT20=hTT20,
                            dphi6=j.delta_phi_to(hTT6),
                            dphi12=j.delta_phi_to(hTT12),
                            dphi20=j.delta_phi_to(hTT20),
                            ppid=j_type[0],
                            pquark=j_type[1],
                            pglue=j_type[2]  # this is redundancy
                    )
			tw.fill_tree()

	t.close()
	pythia.stat()

	tw.write_and_close()


if __name__ == '__main__':
	main()
