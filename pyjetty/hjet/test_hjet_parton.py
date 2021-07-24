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
	parser.add_argument('--charged', help="analyze only the charged particles of the FS", default=False, action='store_true')
	parser.add_argument('--max-jet-pt', help="maximum jet pT to consider", type=float, default=100.)
	args = parser.parse_args()

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	hadron_eta_max = 2.0
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(10.0) & fj.SelectorPtMax(args.max_jet_pt) & fj.SelectorAbsEtaMax(hadron_eta_max - jet_R0)
	# jet_selector = fj.SelectorPtMin(40.0) & fj.SelectorPtMax(200.0) &fj.SelectorAbsEtaMax(hadron_eta_max - jet_R0)
	hTT6_selector 	= fj.SelectorPtMin(6) & fj.SelectorPtMax(7) & fj.SelectorAbsEtaMax(hadron_eta_max)
	hTT12_selector 	= fj.SelectorPtMin(12) & fj.SelectorPtMax(50) & fj.SelectorAbsEtaMax(hadron_eta_max)
	hTT20_selector 	= fj.SelectorPtMin(20) & fj.SelectorPtMax(50) & fj.SelectorAbsEtaMax(hadron_eta_max)

	pythia_fs_part_selection = [pythiafjext.kFinal]
	if args.charged is True:
		pwarning('running with charged particles in the final state')
		pythia_fs_part_selection.append(pythiafjext.kCharged)
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
	tw = treewriter.RTreeWriter(name = 'lsjvsx', file_name = args.output)

	zero_psj = fj.PseudoJet(0,0,10,10)

	if args.nev < 100:
		args.nev = 100
	t = tqdm.tqdm(total = args.nev)
	while t.n < args.nev:
		if not pythia.next():
			continue

		# information about the leading process
		# print(pythia.info.code(), pythia.info.nameProc(pythia.info.code()))
		# continue
		# parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
		partons = pythiafjext.vectorize_select(pythia, [pythiafjext.kParton], 0, True)
		# parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, False)
		# parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, False)
		parts = pythiafjext.vectorize_select(pythia, pythia_fs_part_selection, 0, False)

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
			if hTT12.perp() < 1 and hTT6.perp() < 1 and hTT20.perp() < 1:
				continue

		jets = jet_selector(jet_def(parts))

		# for j in tqdm.tqdm(jets):
		for j in jets:
			t.update(1)
			j_type = match_dR(j, partons, jet_R0 / 2.)
			if j_type[0] is None:
				continue
			j_sd = sd.result(j)
			sd_info = fjcontrib.get_SD_jet_info(j_sd)
			rc_sjets01 = fj.sorted_by_pt(jet_def_rc01(j.constituents()))
			rc_sjets02 = fj.sorted_by_pt(jet_def_rc02(j.constituents()))
			tw.fill_branches(	j 			= j,
								mult 		= len(parts),

								lund 		= [ls for ls in lund_gen.result(j)], 
								dyg1 		= dy_groomer.result(j, 1), 
								sd 			= j_sd, 
								sd_z 		= sd_info.z, 
								sd_mu 		= sd_info.mu, 
								sd_Delta 	= sd_info.dR, 
								lsjet01 	= rc_sjets01[0],
								nsjet01    	= len(rc_sjets01),
								sjet01     	= rc_sjets01,
								lsjet02 	= rc_sjets02[0],
								nsjet02    	= len(rc_sjets02),
								sjet02     	= rc_sjets02,

	                            hTT6 		= hTT6,
	                            hTT12 		= hTT12,
	                            hTT20 		= hTT20,
	                            dphi6 		= j.delta_phi_to(hTT6),
	                            dphi12 		= j.delta_phi_to(hTT12),
	                            dphi20 		= j.delta_phi_to(hTT20),

	                            ppid 		= j_type[0],
	                            pquark 		= j_type[1],
	                            pglue 		= j_type[2],  # this is redundancy

	                            pycode 		= pythia.info.code(),
	                            pysigmagen  = pythia.info.sigmaGen(),
	                            pysigmaerr  = pythia.info.sigmaErr(),
	                            pyid1       = pythia.info.id1pdf(),
	                            pyid2       = pythia.info.id1pdf(),
	                            pyx1 	    = pythia.info.x1pdf(),
	                            pyx2       	= pythia.info.x2pdf(),
	                            pypdf1      = pythia.info.pdf1(),
	                            pyQfac 		= pythia.info.QFac(),
	                            pyalphaS 	= pythia.info.alphaS(),

	                            pypthat 	= pythia.info.pTHat(),
	                            pymhat 		= pythia.info.mHat()
                    )
			tw.fill_tree()

	t.close()
	pythia.stat()

	tw.write_and_close()


if __name__ == '__main__':
	main()
