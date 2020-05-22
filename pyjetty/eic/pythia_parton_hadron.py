#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os
import numpy as np

from pyjetty.mputils import *

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	# could use --py-seed
	parser.add_argument('--user-seed', help='pythia seed', default=1111, type=int)
	parser.add_argument('--output', default="output.root", type=str)
	parser.add_argument('--beta', help='sd beta', default=0, type=float)
	parser.add_argument('--jetR', help='jet radius', default=0.4, type=float)
	args = parser.parse_args()

	if args.user_seed < 0:
		args.user_seed = 1111
	pinfo('user seed for pythia', args.user_seed)
	# mycfg = ['PhaseSpace:pThatMin = 100']
	mycfg = ['Random:setSeed=on', 'Random:seed={}'.format(args.user_seed)]
	mycfg.append('HadronLevel:all=off')
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 100:
		args.nev = 100

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = args.jetR
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	print(jet_def)

	# hadron level - ALICE
	max_eta_hadron = 0.9
	pwarning('max eta for particles after hadronization set to', max_eta_hadron)
	parts_selector_h = fj.SelectorAbsEtaMax(max_eta_hadron)
	jet_selector = fj.SelectorPtMin(20.0) & fj.SelectorPtMax(100.0) & fj.SelectorAbsEtaMax(max_eta_hadron - 1.05 * jet_R0)

	max_eta_parton = max_eta_hadron + 2. * jet_R0
	pwarning('max eta for partons set to', max_eta_parton)
	parts_selector_p = fj.SelectorAbsEtaMax(max_eta_parton)

	outf = ROOT.TFile(args.output.replace('.root', '_beta{}.root'.format(args.beta)), 'recreate')
	outf.cd()
	t = ROOT.TTree('t', 't')
	tw = RTreeWriter(tree=t)

	# event loop
	for iev in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue

		parts_pythia_p = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)
		parts_pythia_p_selected = parts_selector_p(parts_pythia_p)

		hstatus = pythia.forceHadronLevel()
		if not hstatus:
			pwarning('forceHadronLevel false event', iev)
			continue
		# parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kHadron, pythiafjext.kCharged])
		parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)
		parts_pythia_h_selected = parts_selector_h(parts_pythia_h)

		parts_pythia_hch = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, True)
		parts_pythia_hch_selected = parts_selector_h(parts_pythia_hch)

		# pinfo('debug partons...')
		# for p in parts_pythia_p_selected:
		# 	pyp = pythiafjext.getPythia8Particle(p)
		# 	print(pyp.name())
		# pinfo('debug hadrons...')
		# for p in parts_pythia_h_selected:
		# 	pyp = pythiafjext.getPythia8Particle(p)
		# 	print(pyp.name())
		# pinfo('debug ch. hadrons...')
		# for p in parts_pythia_hch_selected:
		# 	pyp = pythiafjext.getPythia8Particle(p)
		# 	print(pyp.name())

		# parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
		jets_p = fj.sorted_by_pt(jet_def(parts_pythia_p))
		jets_h = fj.sorted_by_pt(jet_def(parts_pythia_h))
		jets_ch_h = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_hch)))

		sd = fjcontrib.SoftDrop(args.beta, 0.1, jet_R0)

		for i,jchh in enumerate(jets_ch_h):
			# match hadron (full) jet
			for j,jh in enumerate(jets_h):
				drhh = jchh.delta_R(jh)
				if drhh < jet_R0 / 2.:
					# match parton level jet
					for k,jp in enumerate(jets_p):
						dr = jh.delta_R(jp)
						if dr < jet_R0 / 2.:
							jchh_sd = sd.result(jchh)
							jchh_sd_info = fjcontrib.get_SD_jet_info(jchh_sd)
							jh_sd = sd.result(jh)
							jh_sd_info = fjcontrib.get_SD_jet_info(jh_sd)
							jp_sd = sd.result(jp)
							jp_sd_info = fjcontrib.get_SD_jet_info(jp_sd)

							# pwarning('event', iev)
							# pinfo('matched jets: ch.h:', jchh.pt(), 'h:', jh.pt(), 'p:', jp.pt(), 'dr:', dr)

							tw.fill_branch('iev', iev)
							tw.fill_branch('ch', jchh)
							tw.fill_branch('h', jh)
							tw.fill_branch('p', jp)

							tw.fill_branch('p_zg', jp_sd_info.z)
							tw.fill_branch('p_Rg', jp_sd_info.dR)
							tw.fill_branch('p_thg', jp_sd_info.dR/jet_R0)
							tw.fill_branch('p_mug', jp_sd_info.mu)

							tw.fill_branch('h_zg', jh_sd_info.z)
							tw.fill_branch('h_Rg', jh_sd_info.dR)
							tw.fill_branch('h_thg', jh_sd_info.dR/jet_R0)
							tw.fill_branch('h_mug', jh_sd_info.mu)

							tw.fill_branch('ch_zg', jchh_sd_info.z)
							tw.fill_branch('ch_Rg', jchh_sd_info.dR)
							tw.fill_branch('ch_thg', jchh_sd_info.dR/jet_R0)
							tw.fill_branch('ch_mug', jchh_sd_info.mu)

							tw.fill_tree()
			#print("  |-> SD jet params z={0:10.3f} dR={1:10.3f} mu={2:10.3f}".format(sd_info.z, sd_info.dR, sd_info.mu))

	pythia.stat()
	outf.Write()
	outf.Close()

if __name__ == '__main__':
	main()
