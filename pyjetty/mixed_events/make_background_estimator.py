#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import ROOT
ROOT.gROOT.SetBatch(True)

import tqdm
import argparse
import os
import numpy as np
import array
import copy
import random
import uproot
import pandas as pd

from pyjetty.mputils import logbins, MPBase, CEventSubtractor, CSubtractorJetByJet, RTreeWriter, DataIO, UniqueString
from pyjetty.mputils import fill_tree_matched, fill_tree_data, JetAnalysis, JetAnalysisWithRho, JetAnalysisPerJet
from pyjetty.mputils import matched_pt
from pyjetty.mputils import BoltzmannEvent, BoltzmannSubtractor

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

class BackgroundEstimator(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(grid_size_phi=0.1, grid_size_eta=0.1, particle_eta_max=1., input_file=None)
		super(BackgroundEstimator, self).__init__(**kwargs)
		self._sname = self.name
		if self.input_file:
			self.fin = ROOT.TFile(self.input_file)
			self.grid_eta_phi = self.fin.Get(self.name)	
			self.tlv = ROOT.TLorentzVector()
			self.tlvS = ROOT.TLorentzVector()
		else:
			self.nevents = 0
			self._nbinsX = int(2 * self.particle_eta_max/self.grid_size_eta)
			self._nbinsY = int(2 * (ROOT.TMath.Pi())/self.grid_size_phi)
			self.grid_eta_phi = ROOT.TProfile2D(self._sname, self._sname, 
												self._nbinsX, -self.particle_eta_max, self.particle_eta_max, 
												# int(2 * (ROOT.TMath.Pi())/self.grid_size), -ROOT.TMath.Pi(), ROOT.TMath.Pi())
												self._nbinsY, 0., 2. * ROOT.TMath.Pi())
			self.grid_eta_phi.SetDirectory(0)
			self._sname = UniqueString.str('h_grid_bg_estim_tmp')
			self.grid_eta_phi_tmp = ROOT.TProfile2D(self._sname, self._sname, 
													self._nbinsX, -self.particle_eta_max, self.particle_eta_max, 
													# int(2 * (ROOT.TMath.Pi())/self.grid_size), -ROOT.TMath.Pi(), ROOT.TMath.Pi())
													self._nbinsY, 0., 2. * ROOT.TMath.Pi())
			self.grid_eta_phi_tmp.SetDirectory(0)

	def write(self):
		self.grid_eta_phi.Write()

	def fill_grid(self, parts):
		self.grid_eta_phi_tmp.Reset()
		_tmp = [self.grid_eta_phi_tmp.Fill(p.eta(), p.phi(), p.pt()) for p in parts]
		self.grid_eta_phi.Add(self.grid_eta_phi_tmp, 1.)
		self.nevents = self.nevents + 1

	def get_subtracted_vector(self, p):
		self.tlv.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), 0)
		ieta = self.grid_eta_phi.GetXaxis().FindBin(p.eta())
		iphi = self.grid_eta_phi.GetYaxis().FindBin(p.phi())
		spt = self.grid_eta_phi.GetBinContent(ieta, iphi)
		self.tlvS.SetPtEtaPhiM(spt, p.eta(), p.phi(), 0)
		self.tlv = self.tlv - self.tlvS
		return self.tlv

	def subtract_particle(self, p):
		sp = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
		_sv = self.get_subtracted_vector(p)
		sp.reset_PtYPhiM(_sv.Pt(), _sv.Rapidity(), _sv.Phi(), _sv.M())
		return sp

	def subtracted_particles(self, parts):
		sparts = []
		for p in parts:
			self.tlv.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.m())
			ieta = self.grid_eta_phi.GetXaxis().FindBin(p.eta())
			iphi = self.grid_eta_phi.GetYaxis().FindBin(p.phi())
			spt = self.grid_eta_phi.GetBinContent(ieta, iphi)
			if spt <= self.tlv.Pt():
				self.tlvS.SetPtEtaPhiM(spt, p.eta(), p.phi(), p.m())
				self.tlv = self.tlv - self.tlvS
				sp = fj.PseudoJet(self.tlv.Px(), self.tlv.Py(), self.tlv.Pz(), self.tlv.E())
				sparts.append(sp)
		return sparts


def main_make(args):
	if os.path.isfile(args.output):
		if not args.overwrite:
			print('[i] output', args.output, 'exists - use --overwrite to do just that...')
			return

	print(args)

	# alice specific
	max_eta = 0.9

	outf = ROOT.TFile(args.output, 'recreate')
	outf.cd()

	bgestim = BackgroundEstimator(name='bg_estim', grid_size=0.1, particle_eta_max=max_eta)
	print(bgestim)

	be = None
	data = None
	if len(args.data) > 0:
		data = DataIO(file_list=args.data)
		print(data)
	else:
		be = BoltzmannEvent(mean_pt=0.6, multiplicity=2000 * max_eta * 2, max_eta=max_eta, max_pt=100)
		print(be)

	parts_selector = fj.SelectorAbsEtaMax(max_eta)

	if args.nev < 1:
		args.nev = 1

	### EVENT LOOP STARTS HERE
	for iev in tqdm.tqdm(range(args.nev)):
		if data:
			bg_parts = data.load_event(offset=10000)
		else:
			if be:
				bg_parts = be.generate(offset=10000)
		bgestim.fill_grid(bg_parts)

	outf.cd()
	if be:
		be.write()
	bgestim.write()
	outf.Write()
	outf.Close()
	print('[i] written', outf.GetName())


def main_test(args):
	if os.path.isfile(args.output):
		if not args.overwrite:
			print('[i] output', args.output, 'exists - use --overwrite to do just that...')
			return

	print(args)

	# alice specific
	max_eta = 0.9
	jet_R0 = 0.4

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


	parts_selector = fj.SelectorAbsEtaMax(max_eta)

	fout = ROOT.TFile(args.output, 'recreate')
	fout.cd()
	tw = RTreeWriter(tree_name='tw', fout=fout)

	bgestim = BackgroundEstimator(input_file=args.input, name='bg_estim')
	print(bgestim)

	boltzmann_subtractor = BoltzmannSubtractor(max_pt_subtract=10.)
	print(boltzmann_subtractor)

	be = None
	data = None
	if len(args.data) > 0:
		data = DataIO(file_list=args.data)
		print(data)
	else:
		be = BoltzmannEvent(mean_pt=0.6, multiplicity=2000 * max_eta * 2, max_eta=max_eta, max_pt=100)
		print(be)

	parts_selector = fj.SelectorAbsEtaMax(max_eta)

	if args.nev < 1:
		args.nev = 1

	ja = JetAnalysisWithRho(jet_R=jet_R0, jet_algorithm=fj.antikt_algorithm, particle_eta_max=max_eta)
	jas = JetAnalysisWithRho(jet_R=jet_R0, jet_algorithm=fj.antikt_algorithm, particle_eta_max=max_eta)

	### EVENT LOOP STARTS HERE
	for iev in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue

		parts_pythia = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged])
		parts_gen = parts_selector(parts_pythia)
		jas.analyze_event(parts_gen)
		signal_jets = fj.sorted_by_pt(jet_selector(jas.jets))
		if len(signal_jets) < 1:
			continue
		if data:
			bg_parts = data.load_event(offset=10000)
		else:
			if be:
				bg_parts = be.generate(offset=10000)

		for sjet in signal_jets:
			tmp = [bg_parts.push_back(psj) for psj in sjet.constituents()]

		for sjet in signal_jets:
			tw.fill_branch('sjet', sjet)

			# subtr_parts = bgestim.subtracted_particles(bg_parts)
			# ja.analyze_event(subtr_parts)
			# for j in ja.jets:
			# 	if matched_pt(j, sjet) > 0.5:
			# 		tw.fill_branch('sjet_bg', j)
			# 		tw.fill_branch('delta_pt', sjet.pt() - j.pt())

			ja.analyze_event(bg_parts)
			bg_signal_jets_b = fj.sorted_by_pt(ja.jets)
			for j in bg_signal_jets_b:
				if matched_pt(j, sjet) > 0.5:
					tw.fill_branch('sjet_hybrid', j)
					tw.fill_branch('delta_pt_hybrid', sjet.pt() - j.pt())

			subtr_parts_b = boltzmann_subtractor.subtracted_particles(bg_parts)
			ja.analyze_event(subtr_parts_b)
			bg_signal_jets_b = fj.sorted_by_pt(ja.jets)
			for j in bg_signal_jets_b:
				if matched_pt(j, sjet) > 0.5:
					tw.fill_branch('sjet_bg_b', j)
					tw.fill_branch('delta_pt_b', sjet.pt() - j.pt())

			tw.fill_branch('p', bg_parts)
			# tw.fill_branch('sp', subtr_parts)
			tw.fill_branch('sp_b', subtr_parts_b)
			tw.fill_tree()

	fout.cd()
	fout.Write()
	fout.Close()
	print('[i] written', fout.GetName())

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--ignore-mycfg', help="ignore some settings hardcoded here", default=False, action='store_true')
	parser.add_argument('--output', default="background_estimator.root", type=str)
	parser.add_argument('--input', default="background_estimator.root", type=str)
	parser.add_argument('--data', help='data from a file list', default='', type=str)
	parser.add_argument('--overwrite', help="overwrite output", default=False, action='store_true')
	parser.add_argument('--benchmark', help='benchmark pthat setting - 80 GeV', default=False, action='store_true')

	group_ex = parser.add_mutually_exclusive_group()
	group_ex.add_argument('--make', help='make the estimator', default=False, action='store_true')
	group_ex.add_argument('--test', help='test the estimator', default=False, action='store_true')
	args = parser.parse_args()
	if args.make:
		main_make(args)
	if args.test:
		args.output = args.output.replace('.root', '_test.root')
		main_test(args)

