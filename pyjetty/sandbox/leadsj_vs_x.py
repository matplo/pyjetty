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
	parser.add_argument('--nw', help="no warn", default=False, action='store_true')
	parser.add_argument('--ignore-mycfg', help="ignore some settings hardcoded here", default=False, action='store_true')
	parser.add_argument('--enable-background', help="enable background calc", default=False, action='store_true')
	parser.add_argument('--output', help="output file name", default='leadsj_vs_x_output.root', type=str)

	# for background
	parser.add_argument('--cent-bin', help="centraility bin 0 is the  0-5 percent most central bin", type=int, default=0)
	parser.add_argument('--seed', help="pr gen seed", type=int, default=1111)
	parser.add_argument('--harmonics', help="set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: uniform dN/dphi no harmonics) , (6 : v1 - v2 , v4 - v5) , (7 : v1 - v3 , v5) , (8 : v1 , v3 - v5) , (9 : v1 only) , (10 : v2 only) , (11 : v3 only) , (12 : v4 only) , (13 : v5 only)", 
						type=int, default=5)
	parser.add_argument('--eta', help="set eta range must be uniform (e.g. abs(eta) < 0.9, which is ALICE TPC fiducial acceptance)",
						type=float, default=0.9)
	parser.add_argument('--qa', help="PrintOutQAHistos", default=False, action='store_true')

	parser.add_argument('--dRmax', default=0.25, type=float)
	parser.add_argument('--alpha', default=0, type=float)


	args = parser.parse_args()

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(args.py_pthatmin) & fj.SelectorPtMax(1000.0) & fj.SelectorAbsEtaMax(args.eta - jet_R0)
	# jet_selector = fj.SelectorPtMin(40.0) & fj.SelectorPtMax(200.0) &fj.SelectorAbsEtaMax(1)
	print(jet_def)

	all_jets = []

	# mycfg = ['PhaseSpace:pThatMin = 80']
	# mycfg = ['PhaseSpace:pThatMin = 40']	
	mycfg = ['']	
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
	sd01 = fjcontrib.SoftDrop(0, 0.1, 1.0)
	print (sd01)
	sd02 = fjcontrib.SoftDrop(0, 0.2, 1.0)
	print (sd02)

	# jet_def_rc01 = fj.JetDefinition(fj.cambridge_algorithm, 0.1)
	# jet_def_rc02 = fj.JetDefinition(fj.cambridge_algorithm, 0.2)
	# print (jet_def_rc01)
	# print (jet_def_rc02)
	# rc = fjcontrib.Recluster(jet_def_rc, True)

	jet_def_rc005 = fj.JetDefinition(fj.antikt_algorithm, 0.05)
	jet_def_rc01 = fj.JetDefinition(fj.antikt_algorithm, 0.1)
	jet_def_rc02 = fj.JetDefinition(fj.antikt_algorithm, 0.2)
	jet_def_rc03 = fj.JetDefinition(fj.antikt_algorithm, 0.3)
	print(jet_def_rc005)
	print(jet_def_rc01)
	print (jet_def_rc02)
	print (jet_def_rc03)
	#rc = fjcontrib.Recluster(jet_def_rc, True)

	# tw = treewriter.RTreeWriter(name = 'lsjvsx', file_name = 'leadsj_vs_x.root')
	tw = treewriter.RTreeWriter(name = 'lsjvsx', file_name = args.output)
	tgbkg = None
	be = None
	if args.enable_background:
		# ROOT.gSystem.Load("libpyjetty_TennGen.dylib")
		# tgbkg = ROOT.TennGen() # //constructor
		# tgbkg.SetCentralityBin(args.cent_bin) # //centraility bin 0 is the  0-5 % most central bin
		# tgbkg.SetRandomSeed(args.seed) # //setting the seed
		# tgbkg.SetHarmonics(args.harmonics) # // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: uniform dN/dphi no harmonics) , (6 : v1 - v2 , v4 - v5) , (7 : v1 - v3 , v5) , (8 : v1 , v3 - v5) , (9 : v1 only) , (10 : v2 only) , (11 : v3 only) , (12 : v4 only) , (13 : v5 only)
		# tgbkg.SetEtaRange(args.eta) # //set eta range must be uniform (e.g. |eta| < 0.9, which is ALICE TPC fiducial acceptance)
		# tgbkg.PrintOutQAHistos(args.qa) #
		# tgbkg.InitializeBackground() #

		from pyjetty.mputils import BoltzmannEvent
		be = BoltzmannEvent(mean_pt=0.7, multiplicity=2000 * args.eta * 2, max_eta=max_eta, max_pt=100)
		print(be)

		from pyjetty.mputils import CEventSubtractor, CSubtractorJetByJet
		cs = CEventSubtractor(alpha=args.alpha, max_distance=args.dRmax, max_eta=args.eta, bge_rho_grid_size=0.25, max_pt_correct=100)
		print(cs)


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
				if args.nw:
					continue
				pwarning('Jet with no parton label')
				continue

			j_sd02 = sd02.result(j)
			sd02_info = fjcontrib.get_SD_jet_info(j_sd02)
			j_sd01 = sd01.result(j)
			sd01_info = fjcontrib.get_SD_jet_info(j_sd01)

			rc_sjets005 = fj.sorted_by_pt(jet_def_rc005(j.constituents()))
			rc_sjets01 = fj.sorted_by_pt(jet_def_rc01(j.constituents()))
			rc_sjets02 = fj.sorted_by_pt(jet_def_rc02(j.constituents()))
			rc_sjets03 = fj.sorted_by_pt(jet_def_rc03(j.constituents()))
			tw.fill_branches(	j 			= j, 
								lund 		= [ls for ls in lund_gen.result(j)], 
								dyg1 		= dy_groomer.result(j, 1), 

								sd01 		= j_sd01, 
								sd01_z 		= sd01_info.z, 
								sd01_mu 	= sd01_info.mu, 
								sd01_Delta 	= sd01_info.dR, 

								sd02 		= j_sd02, 
								sd02_z 		= sd02_info.z, 
								sd02_mu 	= sd02_info.mu, 
								sd02_Delta 	= sd02_info.dR, 

								# breaking compatibility
								# sd 			= j_sd, 
								# sd_z 		= sd_info.z, 
								# sd_mu 		= sd_info.mu, 
								# sd_Delta 	= sd_info.dR, 

								lsjet005 	= rc_sjets005[0],
								nsjet005    	= len(rc_sjets005),
								sjet005     	= rc_sjets005,

								lsjet01 	= rc_sjets01[0],
								nsjet01    	= len(rc_sjets01),
								sjet01     	= rc_sjets01,

								lsjet02 	= rc_sjets02[0],
								nsjet02    	= len(rc_sjets02),
								sjet02     	= rc_sjets02,

								lsjet03 	= rc_sjets03[0],
								nsjet03    	= len(rc_sjets03),
								sjet03     	= rc_sjets03,

								ppid       	= j_type[0],
								pquark     	= j_type[1],
								pglue      	= j_type[2], # this is redundancy

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
			if be:
				bg_parts = be.generate(offset=10000)
				full_event = bg_parts
				tmp = [full_event.push_back(psj) for psj in j.constituents()]
				if cs:
					cs_parts = cs.process_event(full_event)
					rho = cs.bge_rho.rho()
					bg_jets = fj.sorted_by_pt(jet_def(cs_parts))
					for bj in bg_jets:
						if fjtools.matched_pt(bj, j) > 0.5:
							pass

			tw.fill_tree()

	pythia.stat()

	tw.write_and_close()


if __name__ == '__main__':
	main()
