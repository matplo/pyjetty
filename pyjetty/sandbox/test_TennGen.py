#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import os
import tqdm
import argparse

from pyjetty.mputils import pwarning, pinfo, perror, treewriter

import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libpyjetty_TennGen.dylib")

def tlv_from_tmcparticle(p):
	_tlv = ROOT.TLorentzVector()
	_tlv.SetPxPyPzE(p.GetPx(), p.GetPy(), p.GetPz(), p.GetEnergy())
	return _tlv


def main():
	parser = argparse.ArgumentParser(description='test the TennGen', prog=os.path.basename(__file__))
	# use pythia?
	# pyconf.add_standard_pythia_args(parser)
	# parser.add_argument('--ignore-mycfg', help="ignore some settings hardcoded here", default=False, action='store_true')
	parser.add_argument('--cent-bin', help="centraility bin 0 is the  0-5 % most central bin", type=int, default=0)
	parser.add_argument('--seed', help="pr gen seed", type=int, default=1111)
	parser.add_argument('--harmonics', help="set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: uniform dN/dphi no harmonics) , (6 : v1 - v2 , v4 - v5) , (7 : v1 - v3 , v5) , (8 : v1 , v3 - v5) , (9 : v1 only) , (10 : v2 only) , (11 : v3 only) , (12 : v4 only) , (13 : v5 only)", 
						type=int, default=5)
	parser.add_argument('--eta', help="set eta range must be uniform (e.g. |eta| < 0.9, which is ALICE TPC fiducial acceptance)",
						type=float, default=0.9)
	parser.add_argument('--qa', help="PrintOutQAHistos", default=False, action='store_true')

	parser.add_argument('--nev', help="number of events", default=100, type=int)

	args = parser.parse_args()

	tgbkg = ROOT.TennGen() # //constructor
	tgbkg.SetCentralityBin(args.cent_bin) # //centraility bin 0 is the  0-5 % most central bin
	tgbkg.SetRandomSeed(args.seed) # //setting the seed
	tgbkg.SetHarmonics(args.harmonics) # // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: uniform dN/dphi no harmonics) , (6 : v1 - v2 , v4 - v5) , (7 : v1 - v3 , v5) , (8 : v1 , v3 - v5) , (9 : v1 only) , (10 : v2 only) , (11 : v3 only) , (12 : v4 only) , (13 : v5 only)
	tgbkg.SetEtaRange(args.eta) # //set eta range must be uniform (e.g. |eta| < 0.9, which is ALICE TPC fiducial acceptance)
	tgbkg.PrintOutQAHistos(args.qa) #
	tgbkg.InitializeBackground() #

	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(10.0) & fj.SelectorPtMax(1000.0) &fj.SelectorAbsEtaMax(args.eta - jet_R0)
	# jet_selector = fj.SelectorPtMin(40.0) & fj.SelectorPtMax(200.0) &fj.SelectorAbsEtaMax(1)
	print(jet_def)

	tw = treewriter.RTreeWriter(name = 'tparts', file_name = 'test_TennGen.root')

	for iev in range(args.nev):
		bg_tclones = tgbkg.GetBackground()
		tgbkg.GetRandomSeed()
		nParticles = bg_tclones.GetEntries();
		pinfo('event', iev, 'number of parts', nParticles)
		# _parts = { 'pt' : [], 'eta' : [], 'phi' : [], 'kf' : []}
		_parts = [[], [], [], []]
		_ = [[_parts[0].append(p[0].Pt()), _parts[1].append(p[0].Eta()), _parts[2].append(p[0].Phi()), _parts[3].append(p[1])] for p in [[tlv_from_tmcparticle(_p), _p.GetKF()] for _p in bg_tclones if _p.GetEnergy()>0]]
		_fj_parts = fjext.vectorize_pt_eta_phi(_parts[0], _parts[1], _parts[2])
		# [tw.fill_branches(pt = _parts[0], eta = _parts[1], _phi = _parts[2])
		_ = [tw.fill_branches(pt = _pfj.perp(), eta = _pfj.eta(), phi = _pfj.phi()) for _pfj in _fj_parts]
		jets = jet_selector(jet_def(_fj_parts))
		pinfo('number of jets 10 < pT < 1000 GeV/c', len(jets))
		tw.fill_branches(j = jets)
		tw.fill_tree()
		bg_tclones.Clear()

	tgbkg.CloserFunction()

	tw.write_and_close()

if __name__ == '__main__':
	main()
