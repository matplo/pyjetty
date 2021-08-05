#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext
import fjtools

import os
import tqdm
import argparse

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

from pyjetty.mputils import pwarning, pinfo, perror, treewriter

import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libpyjetty_TennGen.dylib")


def tlv_from_tmcparticle(p):
	_tlv = ROOT.TLorentzVector()
	_tlv.SetPxPyPzE(p.GetPx(), p.GetPy(), p.GetPz(), p.GetEnergy())
	return _tlv


def matched_jet(j, refjets):
	for rj in refjets:
		mpt = fjtools.matched_pt(j, rj)
		if mpt > 0.5:
			return rj
	return None


def delta_pt_matched(j, refjets, rho):
	for refj in refjets:
		mpt = fjtools.matched_pt(j, refj)
		if mpt > 0.5:
			return j.perp() - refj.perp() - rho * j.area()
	return -1000.


def main():
	parser = argparse.ArgumentParser(description='test the TennGen', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--ignore-mycfg', help="ignore some settings hardcoded here", default=False, action='store_true')
	parser.add_argument('--cent-bin', help="centraility bin 0 is the  0-5 % most central bin", type=int, default=0)
	parser.add_argument('--seed', help="pr gen seed", type=int, default=1111)
	parser.add_argument('--harmonics', help="set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: uniform dN/dphi no harmonics) , (6 : v1 - v2 , v4 - v5) , (7 : v1 - v3 , v5) , (8 : v1 , v3 - v5) , (9 : v1 only) , (10 : v2 only) , (11 : v3 only) , (12 : v4 only) , (13 : v5 only)", 
						type=int, default=5)
	parser.add_argument('--eta', help="set eta range must be uniform (e.g. |eta| < 0.9, which is ALICE TPC fiducial acceptance)",
						type=float, default=0.9)
	parser.add_argument('--qa', help="PrintOutQAHistos", default=False, action='store_true')

	args = parser.parse_args()

	args.py_pthatmin = 100
	mycfg = ['PhaseSpace:pThatMin = {}'.format(args.py_pthatmin)]
	if args.ignore_mycfg:
		mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		perror("pythia initialization failed.")
		return

	tgbkg = ROOT.TennGen() # //constructor
	tgbkg.SetCentralityBin(args.cent_bin) # //centraility bin 0 is the  0-5 % most central bin
	tgbkg.SetRandomSeed(args.seed) # //setting the seed
	tgbkg.SetHarmonics(args.harmonics) # // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: uniform dN/dphi no harmonics) , (6 : v1 - v2 , v4 - v5) , (7 : v1 - v3 , v5) , (8 : v1 , v3 - v5) , (9 : v1 only) , (10 : v2 only) , (11 : v3 only) , (12 : v4 only) , (13 : v5 only)
	tgbkg.SetEtaRange(args.eta) # //set eta range must be uniform (e.g. |eta| < 0.9, which is ALICE TPC fiducial acceptance)
	tgbkg.PrintOutQAHistos(args.qa) #
	tgbkg.InitializeBackground() #

	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector_pythia = fj.SelectorPtMin(args.py_pthatmin) & fj.SelectorPtMax(1000.0) &fj.SelectorAbsEtaMax(args.eta - jet_R0)
	jet_selector_hybrid = fj.SelectorPtMin(10) & fj.SelectorPtMax(1000.0) &fj.SelectorAbsEtaMax(args.eta - jet_R0)
	# jet_selector = fj.SelectorPtMin(40.0) & fj.SelectorPtMax(200.0) &fj.SelectorAbsEtaMax(1)
	parts_selector = fj.SelectorAbsEtaMax(args.eta)
	print(jet_def)

	tw = treewriter.RTreeWriter(name = 'tparts', file_name = 'test_TennGen.root')

	if args.nev < 100:
		args.nev = 100
	pbar = tqdm.tqdm(total = args.nev)
	while pbar.n < args.nev:
		if not pythia.next():
			continue

		# get pythia particles
		# parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, False)
		_py_fj_parts = parts_selector(pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, False))
		# get jets w/o area determination
		# pythia_jets = jet_selector_pythia(jet_def(_py_fj_parts))

		# with area determination
		jet_area_def = fj.AreaDefinition(fj.active_area, fj.GhostedAreaSpec(args.eta))
		cs = fj.ClusterSequenceArea(_py_fj_parts, jet_def, jet_area_def)
		pythia_jets = jet_selector_pythia(cs.inclusive_jets())

		if len(pythia_jets) < 1:
			continue
		pbar.update(1)
		tw.fill_branches(pyj = pythia_jets)

		# now generate bg
		bg_tclones = tgbkg.GetBackground()
		# tgbkg.GetRandomSeed()
		nParticles = bg_tclones.GetEntries();
		# pinfo('event', pbar.n, 'number of parts', nParticles)
		# _parts = { 'pt' : [], 'eta' : [], 'phi' : [], 'kf' : []}
		_parts = [[], [], [], []]
		_ = [[_parts[0].append(p[0].Pt()), _parts[1].append(p[0].Eta()), _parts[2].append(p[0].Phi()), _parts[3].append(p[1])] for p in [[tlv_from_tmcparticle(_p), _p.GetKF()] for _p in bg_tclones if _p.GetEnergy()>0]]
		_bg_fj_parts = fjext.vectorize_pt_eta_phi(_parts[0], _parts[1], _parts[2], 1000) #bg particles with index > 1000

		# add background and pythia 
		_fj_parts = []
		_ = [_fj_parts.append(_p) for _p in _py_fj_parts]
		_ = [_fj_parts.append(_p) for _p in _bg_fj_parts]

		# stream all particles
		_ = [tw.fill_branches(part_pt = _pfj.perp(), part_eta = _pfj.eta(), part_phi = _pfj.phi(), part_idx=_pfj.user_index()) for _pfj in _fj_parts]

		# find jets in the hybrid event
		# w/o area
		# jets = jet_selector_hybrid(jet_def(_fj_parts))
		# w/area
		cs_hybrid = fj.ClusterSequenceArea(_fj_parts, jet_def, jet_area_def)
		jets = jet_selector_hybrid(cs_hybrid.inclusive_jets())
		# stream jets from the hybrid event
		tw.fill_branches(j = jets)

		# estimate the background
		bg_rho_range = fj.SelectorAbsEtaMax(args.eta * 1.1)
		bg_jet_def = fj.JetDefinition(fj.kt_algorithm, jet_R0)
		bg_area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(args.eta))
		# bg_area_def = fj.AreaDefinition(fj.active_area, fj.GhostedAreaSpec(args.eta)) #active area defunct for bg estim
		bg_estimator = fj.JetMedianBackgroundEstimator(bg_rho_range, bg_jet_def, bg_area_def)
		bg_estimator.set_particles(_fj_parts)
		if len(_fj_parts) < 0:
			perror('no particles in the hybrid event?')
			continue
		rho = bg_estimator.rho()
		sigma = bg_estimator.sigma()
		corr_jet_pt = [j.pt() - j.area() * rho for j in jets]
		# matches = [j.perp(), matched_jet(j, pythia_jets) for j in jets]
		delta_pt = [delta_pt_matched(j, pythia_jets, rho) for j in jets]
		tw.fill_branches(j_corr_pt = corr_jet_pt, dpt = delta_pt)
		tw.fill_branches(rho = rho, rho_sigma = sigma)

		tw.fill_tree()
		bg_tclones.Clear()

	pbar.close()

	tgbkg.CloserFunction()
	tw.write_and_close()

if __name__ == '__main__':
	main()
