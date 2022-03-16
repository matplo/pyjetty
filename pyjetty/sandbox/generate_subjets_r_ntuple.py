#!/usr/bin/env python3

import operator as op
import itertools as it
import sys
import os
import argparse
from tqdm import tqdm
from heppy.pythiautils import configuration as pyconf
import pythiaext
import pythiafjext
import pythia8
import fjtools
import fjcontrib
import fjext
import fastjet as fj
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from tqdm import tqdm
import argparse
import os


def get_args_from_settings(ssettings):
	sys.argv = [' '] + ssettings.split()
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly')
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--output', default="subjets_r_ntuple.root", type=str)
	parser.add_argument('--user-seed', help='pythia seed',
											default=1111, type=int)
	args = parser.parse_args()
	return args


def get_leading_charge_jet(j):
	if j.has_constituents():
		p = fj.sorted_by_pt(j.constituents())[0]
		pyp = pythiafjext.getPythia8Particle(p)
		return pyp.charge()
	return 0


class SubjetCorrelator(object):
	def __init__(self, sjr, algorithm = fj.cambridge_algorithm):
		self.r = sjr
		self.algorithm = algorithm
		self.jet_def = fj.JetDefinition(self.algorithm, self.r)
		self.subjets = []
		self.pairs = []
		self.mixed = []
		_tname = 's_{}_{}'.format(self.r, str(self.algorithm))
		self.tnS = ROOT.TNtuple(_tname, _tname, 'ptj:dtheta:dz:kt:z1:z2:cdzdtheta:qq')
		_tname = 'b_{}_{}'.format(self.r, str(self.algorithm))
		self.tnB = ROOT.TNtuple(_tname, _tname, 'ptj:dtheta:dz:kt:z1:z2:cdzdtheta:qq')
		_tname = 'j_{}_{}'.format(self.r, str(self.algorithm))
		self.tnJ = ROOT.TNtuple(_tname, _tname, 'ptj:eta:phi:nc')

	def process_jets(self, jets):
		for j in tqdm(jets, desc='jets'):
			self.tnJ.Fill(j.perp(), j.eta(), j.phi(), len(j.constituents()))
			for cp in it.combinations(fj.sorted_by_pt(j.constituents()), 2):
				_dtheta = abs(j.delta_R(cp[0]) - j.delta_R(cp[1]))
				_z1 = cp[0].perp()/j.perp()
				_z2 = cp[1].perp()/j.perp()
				_dz = (_z1 - _z2)
				_kt = (cp[0].perp() - cp[1].perp()) / 2.
				_pyp1 = pythiafjext.getPythia8Particle(cp[0])
				_pyp2 = pythiafjext.getPythia8Particle(cp[1])
				_qq = _pyp1.charge() * _pyp2.charge()
				self.tnS.Fill(j.perp(), _dtheta, _dz, _kt, _z1, _z2, _dz*_dtheta, _qq)
		_mixed_pairs_jets = [pair for pair in it.combinations(jets, 2)]
		for m in _mixed_pairs_jets:
			print(m[0].perp(), m[1].perp())
		# sjs = self.jet_def(j.constituents())
		# _ = [self.pairs.append(pair) for pair in it.combinations(sjs, 2)]


def main():
	mycfg = []
	ssettings = "--py-ecm 5000 --user-seed=100000 --nev 1000 --py-pthatmin 100"
	args = get_args_from_settings(ssettings)
	pythia_hard = pyconf.create_and_init_pythia_from_args(args, mycfg)

	max_eta_hadron = 2
	parts_selector_h = fj.SelectorAbsEtaMax(max_eta_hadron)
	jet_R0 = 0.4
	jet_selector = fj.SelectorPtMin(100.0) & fj.SelectorPtMax(120.0) & fj.SelectorAbsEtaMax(max_eta_hadron - 1.05 * jet_R0)
	n_pileup = 1  # 5

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	print(jet_def)
 
	lead_jets = []

	for n in tqdm(range(args.nev)):
		if not pythia_hard.next():
				continue
		parts_pythia_h = pythiafjext.vectorize_select(
				pythia_hard, [pythiafjext.kFinal, pythiafjext.kCharged], 0, True)
		parts_pythia_h_selected = parts_selector_h(parts_pythia_h)

		jets_h = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h_selected))) 
		if len(jets_h) < 1:
				continue

		lead_jets.append(jets_h[0])
  
	pythia_hard.stat()

	#subjet reconstruction definition
	fout = ROOT.TFile(args.output, 'recreate')
	fout.cd()
	correlator = SubjetCorrelator(jet_R0/10.)
	correlator.process_jets(lead_jets)
	fout.Write()
	fout.Close()
	print('[i] written ', fout.GetName())


if __name__ == '__main__':
	main()
