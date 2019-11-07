#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os
import numpy as np
import array
import copy
import random
import uproot
import pandas as pd

from pyjetty.mputils import logbins
from pyjetty.mputils import MPBase
from pyjetty.mputils import BoltzmannEvent
from pyjetty.mputils import CEventSubtractor, CSubtractorJetByJet
from pyjetty.mputils import RTreeWriter
from pyjetty.mputils import fill_tree_matched, fill_tree_data, JetAnalysis, JetAnalysisWithRho, JetAnalysisPerJet
from pyjetty.mputils import DataBackgroundIO

from alice_efficiency import AliceChargedParticleEfficiency

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

import ROOT
ROOT.gROOT.SetBatch(True)


def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--ignore-mycfg', help="ignore some settings hardcoded here", default=False, action='store_true')
	parser.add_argument('--output', default="output.root", type=str)
	parser.add_argument('--alpha', default=0, type=float)
	parser.add_argument('--dRmax', default=0.0, type=float)
	parser.add_argument('--zcut', default=0.1, type=float)
	parser.add_argument('--overwrite', help="overwrite output", default=False, action='store_true')
	parser.add_argument('--embed', help='run embedding from a file list', default='', type=str)
	parser.add_argument('--SDsignal', help='embed only SD signal prongs', default=False, action='store_true')
	parser.add_argument('--SDsignal-single', help='embed only SD signal - only leading prong!', default=False, action='store_true')
	parser.add_argument('--efficiency', help='apply charged particle efficiency', default=False, action='store_true')
	parser.add_argument('--benchmark', help='benchmark pthat setting - 80 GeV', default=False, action='store_true')
	parser.add_argument('--csjet', help='constituent subtration jet-by-jet', default=False, action='store_true')
	args = parser.parse_args()

	if args.output == 'output.root':
		args.output = 'output_alpha_{}_dRmax_{}_SDzcut_{}.root'.format(args.alpha, args.dRmax, args.zcut)
		if args.py_seed >= 0:
			args.output = 'output_alpha_{}_dRmax_{}_SDzcut_{}_seed_{}.root'.format(args.alpha, args.dRmax, args.zcut, args.py_seed)
		if args.embed:
			args.output = args.output.replace('.root', '_emb.root')
		if args.efficiency:
			args.output = args.output.replace('.root', '_effi.root')
		if args.SDsignal:
			args.output = args.output.replace('.root', '_SDsignal.root')
		if args.SDsignal_single:
			args.output = args.output.replace('.root', '_SDsignal_single.root')
		if args.csjet:
			args.output = args.output.replace('.root', '_csjet.root')

	if os.path.isfile(args.output):
		if not args.overwrite:
			print('[i] output', args.output, 'exists - use --overwrite to do just that...')
			return

	print(args)

	# alice specific
	max_eta = 0.9

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	print(jet_def)

	mycfg = []
	if args.benchmark:
		mycfg = ['PhaseSpace:pThatMin = 80', 'PhaseSpace:pThatMax = -1']
		jet_selector = fj.SelectorPtMin(80.0) & fj.SelectorPtMax(100.0) & fj.SelectorAbsEtaMax(max_eta - 1.05 * jet_R0)
		# jet_selector_cs = fj.SelectorPtMin(50.0) & fj.SelectorAbsEtaMax(max_eta - 1.05 * jet_R0)
	else:
		args.py_biaspow = 4
		args.py_biasref = 10
		jet_selector = fj.SelectorPtMin(20) & fj.SelectorAbsEtaMax(max_eta - 1.05 * jet_R0)
		# jet_selector_cs = fj.SelectorPtMin(50.0) & fj.SelectorAbsEtaMax(max_eta - 1.05 * jet_R0)

	if args.ignore_mycfg:
		mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return

	sd_zcut = args.zcut
	sd = fjcontrib.SoftDrop(0, sd_zcut, jet_R0)

	jarho = JetAnalysisWithRho(jet_R=jet_R0, jet_algorithm=fj.antikt_algorithm, particle_eta_max=max_eta)
	ja = JetAnalysis(jet_R=jet_R0, jet_algorithm=fj.antikt_algorithm, particle_eta_max=max_eta)

	be = None
	embd = None
	if len(args.embed) > 0:
		embd = DataBackgroundIO(file_list=args.embed)
		print(embd)
	else:
		be = BoltzmannEvent(mean_pt=0.6, multiplicity=2000 * max_eta * 2, max_eta=max_eta, max_pt=100)
		print(be)

	cs = None
	if args.dRmax > 0:
		cs = CEventSubtractor(alpha=args.alpha, max_distance=args.dRmax, max_eta=max_eta, bge_rho_grid_size=0.25, max_pt_correct=100)
		print(cs)
	csjet = None
	if args.csjet:
		csjet = CSubtractorJetByJet(max_eta=max_eta, bge_rho_grid_size=0.25)
		print(csjet)

	parts_selector = fj.SelectorAbsEtaMax(max_eta)

	if args.nev < 1:
		args.nev = 1

	outf = ROOT.TFile(args.output, 'recreate')
	outf.cd()
	t = ROOT.TTree('t', 't')
	tw = RTreeWriter(tree=t)
	te = ROOT.TTree('te', 'te')
	twe = RTreeWriter(tree=te)

	# effi_pp = AliceChargedParticleEfficiency(csystem='pp')
	effi_PbPb = None
	if args.efficiency:
		effi_PbPb = AliceChargedParticleEfficiency(csystem='PbPb')
		print(effi_PbPb)

	### EVENT LOOP STARTS HERE
	for iev in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue

		parts_pythia = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged])
		parts_gen = parts_selector(parts_pythia)
		if effi_PbPb:
			parts = effi_PbPb.apply_efficiency(parts_gen)
		else:
			parts = parts_gen

		signal_jets = fj.sorted_by_pt(jet_selector(jet_def(parts)))
		if len(signal_jets) < 1:
			continue


		for sjet in signal_jets:
			if args.SDsignal or args.SDsignal_single:
				sd_sjet = sd.result(sjet)
				pe1 = fj.PseudoJet()
				pe2 = fj.PseudoJet()
				has_parents = sd_sjet.has_parents(pe1, pe2)
				if has_parents:
					jparts = fj.vectorPJ()
					pe1.set_user_index(0)
					pe2.set_user_index(1)
					if args.SDsignal_single:
						if pe1.pt() > pe2.pt():
							jparts.push_back(pe1)
						else:
							jparts.push_back(pe2)
					else:
						jparts.push_back(pe1)
						jparts.push_back(pe2)						
					sjets = fj.sorted_by_pt(jet_selector(jet_def(jparts)))
					if len(sjets) == 1:
						sjet = sjets[0]
					else:
						continue
				else:
					continue
			if embd:
				bg_parts = embd.load_event(offset=10000)
				# for p in bg_parts:
				# 	print(p.user_index())
			else:
				bg_parts = be.generate(offset=10000)
				# for p in bg_parts:
				# 	print(p.user_index())
			full_event = bg_parts
			tmp = [full_event.push_back(psj) for psj in sjet.constituents()]
			if cs:
				cs_parts = cs.process_event(full_event)
				rho = cs.bge_rho.rho()
				jarho.analyze_event(cs_parts)							
				tmp = [fill_tree_data(ej, twe, sd, rho, iev, pythia.info.weight(), pythia.info.sigmaGen()) for ej in jarho.jets]
				tmp = [fill_tree_matched(sjet, ej, tw, sd, rho, iev, pythia.info.weight(), pythia.info.sigmaGen()) for ej in jarho.jets]
			else:
				jarho.analyze_event(full_event)
				rho = jarho.rho
				if csjet:
					#_csjet = fjcontrib.ConstituentSubtractor(jarho.bg_estimator)
					# subtr_jets = [_csjet.result(ej) for ej in jarho.jets]
					csjet.set_event_particles(full_event)
					#subtr_jets = [csjet.process_jet(ej) for ej in jarho.jets]
					#print ('jbyj cs', len(subtr_jets), 'from', len(jarho.jets))
					#subtr_jets_wconstits = [_j for _j in subtr_jets if _j.has_constituents()]
					#for _j in subtr_jets_wconstits:
					#	print(len(_j.constituents()))
					subtr_jets_wconstits = csjet.process_jets(jarho.jets)
					japerjet = JetAnalysisPerJet(jet_R=jet_R0, jet_algorithm=fj.antikt_algorithm, particle_eta_max=max_eta, input_jets=subtr_jets_wconstits)
					# for _j in japerjet.jets:
					# 	for _c in _j.constituents():
					# 		if _c.user_index() >= 0:
					# 			print('user index kept?', _c.user_index())
					# 		# else:
					# 		# 	print('user index kept?', _c.user_index(), _c.pt())
					# 	_sd_j = sd.result(_j)
					# https://phab.hepforge.org/source/fastjetsvn/browse/contrib/contribs/RecursiveTools/trunk/Recluster.cc L 270
					# tmp = [fill_tree_matched(sjet, ej, tw, sd, rho, iev, pythia.info.sigmaGen()) for ej in subtr_jets_wcs]
					tmp = [fill_tree_data(ej, twe, sd, rho, iev, pythia.info.weight(), pythia.info.sigmaGen()) for ej in japerjet.jets]
					tmp = [fill_tree_matched(sjet, ej, tw, sd, rho, iev, pythia.info.weight(), pythia.info.sigmaGen()) for ej in japerjet.jets]
				else:
					tmp = [fill_tree_data(ej, twe, sd, rho, iev, pythia.info.weight(), pythia.info.sigmaGen()) for ej in jarho.jets]
					tmp = [fill_tree_matched(sjet, ej, tw, sd, rho, iev, pythia.info.weight(), pythia.info.sigmaGen()) for ej in jarho.jets]
	pythia.stat()
	outf.Write()
	outf.Close()
	print('[i] written', outf.GetName())


if __name__ == '__main__':
	main()
