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

def tuple_to_string(tup):
	rets = [str(t) for t in tup]
	return  ' '.join(rets) 

class SplitInfo(object):
	def __init__(self, pythia):
		self.n = 0
		self.pythia = pythia
		self.event = self.pythia.event
		self.sinfo = []

	def find_split(self, i, n=0):
		if n == 0:
			self.sinfo = []
			self.sinfo.append(tuple_to_string(['---', i, self.event[i].name()]))
		ns = ' '.join([' ' for _in in range(n)])
		pyp = self.event[i]
		m1 = self.event[i].mother1()
		m2 = self.event[i].mother2()
		if m1 == m2:
			m = self.event[m1]
			d1 = m.daughter1()
			d2 = m.daughter2()
			if d1 == d2:
				self.sinfo.append(tuple_to_string([ns+'-', m1, '->', d1, '(same daughter)']))
				self.find_split(m1, n+1)
			else:
				tuple_to_string([ns, m1, '->', d1, d2, '(different daughters)'])
		else:
			# self.sinfo.append(tuple_to_string([ns, i, 'different mothers?', m1, m2]))
			self.sinfo.append(tuple_to_string([ns, ' * m1', m1, self.event[m1].name(), '->', i, pyp.name()]))
			self.sinfo.append(tuple_to_string([ns, ' * m2', m2, self.event[m2].name(), '->', i, pyp.name()]))
		if n == 0:
			self.sinfo.append(tuple_to_string(['---', i, self.event[i].name()]))

	def __str__(self):
		return '\n'.join(reversed(self.sinfo)) + '\n'


def parton_splittings_parton_only(pythia, selected_jet, maxR=0.4):
	spinfo = SplitInfo(pythia)
	for p in fj.sorted_by_pt(selected_jet.constituents()):
		pyp = pythiafjext.getPythia8Particle(p)
		if pyp.isParton():
			if p.delta_R(selected_jet) < maxR:
				spinfo.find_split(p.user_index())
				print (spinfo)
		else:
			continue

def parton_splittings(pythia, selected_jet, maxR=0.4):
	spinfo = SplitInfo(pythia)
	for p in fj.sorted_by_pt(selected_jet.constituents()):
		pyp = pythiafjext.getPythia8Particle(p)
		spinfo.find_split(p.user_index())
		print (spinfo)


def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	# could use --py-seed
	parser.add_argument('--user-seed', help='pythia seed', default=1111, type=int)
	parser.add_argument('--output', default="output_ph_residuals.root", type=str)
	parser.add_argument('--beta', help='sd beta', default=0, type=float)
	args = parser.parse_args()

	if args.user_seed < 0:
		args.user_seed = 1111
	pinfo('user seed for pythia', args.user_seed)
	# mycfg = ['PhaseSpace:pThatMin = 100']
	mycfg = ['Random:setSeed=on', 'Random:seed={}'.format(args.user_seed)]
	mycfg.append('HadronLevel:all=off')
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 1:
		args.nev = 1

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.4
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

		parts_pythia_p = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], add_particle_info = True)
		parts_pythia_p_selected = parts_selector_p(parts_pythia_p)
		jets_p = fj.sorted_by_pt(jet_def(parts_pythia_p))

		# kt split on parton level
		if jets_p[0].pt() < 50:
			continue
		parton_splittings(pythia, jets_p[0], jet_R0)
		break

		hstatus = pythia.forceHadronLevel()
		if not hstatus:
			pwarning('forceHadronLevel false event', iev)
			continue
		# parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kHadron, pythiafjext.kCharged])
		parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], add_particle_info = True)
		parts_pythia_h_selected = parts_selector_h(parts_pythia_h)

		parts_pythia_hch = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], add_particle_info = True)
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
		jets_h = fj.sorted_by_pt(jet_def(parts_pythia_h))
		jets_ch_h = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_hch)))

		sd = fjcontrib.SoftDrop(args.beta, 0.1, jet_R0)

		for j,jh in enumerate(jets_h):
			# match parton level jet
			for k,jp in enumerate(jets_p):
				dr = jh.delta_R(jp)
				if dr < jet_R0 / 2.:
					jh_sd = sd.result(jh)
					jh_sd_info = fjcontrib.get_SD_jet_info(jh_sd)
					jp_sd = sd.result(jp)
					jp_sd_info = fjcontrib.get_SD_jet_info(jp_sd)

					tw.fill_branch('iev', iev)
					tw.fill_branch('jh', jh)
					tw.fill_branch('jp', jp)

					tw.fill_branch('jp_zg', jp_sd_info.z)
					tw.fill_branch('jp_Rg', jp_sd_info.dR)
					tw.fill_branch('jp_thg', jp_sd_info.dR/jet_R0)
					tw.fill_branch('jp_mug', jp_sd_info.mu)

					tw.fill_branch('jh_zg', jh_sd_info.z)
					tw.fill_branch('jh_Rg', jh_sd_info.dR)
					tw.fill_branch('jh_thg', jh_sd_info.dR/jet_R0)
					tw.fill_branch('jh_mug', jh_sd_info.mu)

					tw.fill_tree()
			#print("  |-> SD jet params z={0:10.3f} dR={1:10.3f} mu={2:10.3f}".format(sd_info.z, sd_info.dR, sd_info.mu))

	pythia.stat()
	outf.Write()
	outf.Close()

if __name__ == '__main__':
	main()
