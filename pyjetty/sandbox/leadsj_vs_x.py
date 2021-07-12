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
		return None, False, False
	p = fj.sorted_by_pt(mps)[0]
	pyp = pythiafjext.getPythia8Particle(p)
	# print(p, pyp.id(), pyp.isQuark(), pyp.isGluon())
	return pyp.id(), pyp.isQuark(), pyp.isGluon()

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
	jet_selector = fj.SelectorPtMin(80.0) & fj.SelectorPtMax(100.0) &fj.SelectorAbsEtaMax(1)
	# jet_selector = fj.SelectorPtMin(40.0) & fj.SelectorPtMax(200.0) &fj.SelectorAbsEtaMax(1)
	print(jet_def)

	all_jets = []

	mycfg = ['PhaseSpace:pThatMin = 80']
	# mycfg = ['PhaseSpace:pThatMin = 40']	
	if args.ignore_mycfg:
		mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		perror("pythia initialization failed.")
		return

	jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
	lund_gen = fjcontrib.LundGenerator(jet_def_lund)
	print (lund_gen.description())
	dy_groomer = fjcontrib.DynamicalGroomer(jet_def_lund)
	print (dy_groomer.description())
	# sd = fjcontrib.SoftDrop(0, 0.1, 1.0)
	sd = fjcontrib.SoftDrop(0, 0.2, 1.0)
	print (sd)

	# jet_def_rc01 = fj.JetDefinition(fj.cambridge_algorithm, 0.1)
	# jet_def_rc02 = fj.JetDefinition(fj.cambridge_algorithm, 0.2)
	# print (jet_def_rc01)
	# print (jet_def_rc02)
	# rc = fjcontrib.Recluster(jet_def_rc, True)

	jet_def_rc01 = fj.JetDefinition(fj.antikt_algorithm, 0.1)
	jet_def_rc02 = fj.JetDefinition(fj.antikt_algorithm, 0.2)
	print (jet_def_rc01)
	print (jet_def_rc02)
	#rc = fjcontrib.Recluster(jet_def_rc, True)

	# tw = treewriter.RTreeWriter(name = 'lsjvsx', file_name = 'leadsj_vs_x.root')
	tw = treewriter.RTreeWriter(name = 'lsjvsx', file_name = 'leadsj_vs_x_bias80.root')

	if args.nev < 100:
		args.nev = 100
	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		# parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
		partons = pythiafjext.vectorize_select(pythia, [pythiafjext.kParton], 0, True)
		parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, False)
		# parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, False)
		jets = jet_selector(jet_def(parts))

		# for j in tqdm.tqdm(jets):
		for j in jets:
			j_type = match_dR(j, partons, jet_R0 / 2.)
			if j_type[0] is None:
				pwarning('Jet with no parton label')
				continue
			j_sd = sd.result(j)
			sd_info = fjcontrib.get_SD_jet_info(j_sd)
			rc_sjets01 = fj.sorted_by_pt(jet_def_rc01(j.constituents()))
			rc_sjets02 = fj.sorted_by_pt(jet_def_rc02(j.constituents()))
			tw.fill_branches(j 			= j, 
							 lund 		= [ls for ls in lund_gen.result(j)], 
							 dyg1 		= dy_groomer.result(j, 1), 
							 sd 		= j_sd, 
							 sd_z 		= sd_info.z, 
							 sd_mu 		= sd_info.mu, 
							 sd_Delta 	= sd_info.dR, 
							 lsjet01 	= rc_sjets01[0],
							 nsjet01    = len(rc_sjets01),
							 sjet01     = rc_sjets01,
							 lsjet02 	= rc_sjets02[0],
							 nsjet02    = len(rc_sjets02),
							 sjet02     = rc_sjets02,
							 ppid       = j_type[0],
							 pquark     = j_type[1],
							 pglue      = j_type[2] # this is redundancy
							 )
			tw.fill_tree()

	pythia.stat()

	tw.write_and_close()


if __name__ == '__main__':
	main()
