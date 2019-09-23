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

import csevent

import ROOT
ROOT.gROOT.SetBatch(True)

import numpy as np
import array
import copy

def logbins(xmin, xmax, nbins):
	if xmin <= 0:
		xmin = 1e-2
	lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
	arr = array.array('f', lspace)
	return arr


class BoltzmannEvent(object):
	def __init__(self, **kwargs):
		self.configure_constants(mean_pt=0.7, multiplicity=1, max_eta=1, max_pt=100)
		for key, value in kwargs.items():
			self.__setattr__(key, value)
		self.particles = fj.vectorPJ()
		self.funbg = ROOT.TF1("funbg", "2. / [0] * x * TMath::Exp(-(2. / [0]) * x)", 0, self.max_pt, 1);
		self.funbg.SetParameter(0, self.mean_pt)
		self.funbg.SetNpx(1000)
		self.ROOT_random = ROOT.TRandom()
		self.histogram_pt = ROOT.TH1F("BoltzmannEvent_pt", "BoltzmannEvent_pt;p_{T} (GeV/c)", 100, logbins(1e-1, self.max_pt, 100))
		self.histogram_pt.SetDirectory(0)
		self.histogram_eta = ROOT.TH1F("BoltzmannEvent_eta", "BoltzmannEvent_eta;#eta", 100, -self.max_eta, self.max_eta)
		self.histogram_eta.SetDirectory(0)
		self.histogram_phi = ROOT.TH1F("BoltzmannEvent_phi", "BoltzmannEvent_phi;#varphi (rad)", 100, -ROOT.TMath.Pi(), ROOT.TMath.Pi())
		self.histogram_phi.SetDirectory(0)
		self.nEvent = 0
		print (self)

	def configure_constants(self, **kwargs):
		for key, value in kwargs.items():
			self.__setattr__(key, value)

	def __str__(self):
		return "BoltzmannEvent with mean_pt={} max_pt={} multiplicity={} max_eta={}".format(self.mean_pt, self.max_pt, self.multiplicity, self.max_eta)

	def _boltzmann(self, pt):
		return 2. / self.mean_pt * pt * math.exp(-(2. / self.mean_pt) * pt);

	def generate(self, multiplicity=None, offset=0):
		if multiplicity:
			self.multiplicity = multiplicity
		self.particles.clear()
		for n in range(0, self.multiplicity):
			_pt  = self.funbg.GetRandom(0, self.max_pt)
			_eta = self.ROOT_random.Rndm() * self.max_eta * 2. - self.max_eta;
			_phi = self.ROOT_random.Rndm() * ROOT.TMath.Pi() * 2. - ROOT.TMath.Pi();
			_p = fj.PseudoJet()
			_p.reset_PtYPhiM (_pt, _eta, _phi, 0.0)
			_p.set_user_index(n + offset)
			self.particles.push_back(_p)
			self.histogram_pt.Fill(_pt)
			self.histogram_eta.Fill(_eta)
			self.histogram_phi.Fill(_phi)
		self.nEvent = self.nEvent + 1
		return self.particles

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--ignore-mycfg', help="ignore some settings hardcoded here", default=False, action='store_true')
	parser.add_argument('--output', default="output.root", type=str)
	parser.add_argument('--alpha', default=0, type=float)
	parser.add_argument('--dRmax', default=0.25, type=float)
	parser.add_argument('--overwrite', help="overwrite output", default=False, action='store_true')

	args = parser.parse_args()

	if args.output == 'output.root':
		args.output = 'output_alpha_{}_dRmax_{}.root'.format(args.alpha, args.dRmax)

	if os.path.isfile(args.output):
		if not args.overwrite:
			print('[i] output', args.output, 'exists - use --overwrite to do just that...')
			return

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(100.0) & fj.SelectorPtMax(110.0) & fj.SelectorAbsEtaMax(1 - 1.05 * jet_R0)
	jet_selector_cs = fj.SelectorPtMin(50.0) & fj.SelectorAbsEtaMax(1 - 1.05 * jet_R0)
	print(jet_def)

	mycfg = ['PhaseSpace:pThatMin = 100', 'PhaseSpace:pThatMax = 110']
	if args.ignore_mycfg:
		mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		print("[e] pythia initialization failed.")
		return

	sd = fjcontrib.SoftDrop(0, 0.1, 1.0)
	cs = csevent.CSEventSubtractor(alpha=args.alpha, max_distance=args.dRmax, max_eta=1, bge_rho_grid_size=0.1, max_pt_correct=100)
	be = BoltzmannEvent(mean_pt=0.6, multiplicity=2000, max_eta=1, max_pt=100)

	if args.nev < 100:
		args.nev = 100

	outf = ROOT.TFile(args.output, 'recreate')
	outf.cd()
	tn = ROOT.TNtuple('tn', 'tn', 'n:pt:phi:eta:cspt:csphi:cseta:dR:dpt:rg:csrg:z:csz:dRg:dzg')
	hpt = ROOT.TH1F('hpt', 'hpt', 100, 50, 150)
	hptcs = ROOT.TH1F('hptcs', 'hptcs', 100, 50, 150)
	hdpt = ROOT.TH1F('hdpt', 'hdpt', 100, -50, 50)
	hrg = ROOT.TH1F('hrg', 'hrg', 100, -1.1, 0.9)
	hrgcs = ROOT.TH1F('hrgcs', 'hrgcs', 100, -1.1, 0.9)
	hdrg = ROOT.TH1F('hdrg', 'hdrg', 100, -1.1, 0.9)

	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue

		parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
		signal_jets = fj.sorted_by_pt(jet_selector(jet_def(parts)))

		if len(signal_jets) < 1:
			continue

		bg_parts = be.generate()

		full_event = bg_parts
		sjet = signal_jets[0]
		for psj in sjet.constituents(): full_event.push_back(psj)

		cs_parts = cs.process_event(full_event)
		cs_signal_jets = fj.sorted_by_pt(jet_selector_cs(jet_def(cs_parts)))

		sd_signal_jet = sd.result(sjet)
		sd_info_signal_jet = fjcontrib.get_SD_jet_info(sd_signal_jet)

		for j in cs_signal_jets:
			sd_j = sd.result(j)
			sd_info_j = fjcontrib.get_SD_jet_info(sd_j)
			rho = cs.bge_rho.rho()
			tn.Fill(i, 
					sjet.pt(), sjet.phi(), sjet.eta(), 
					j.pt(), j.phi(), j.eta(),
					j.delta_R(sjet), j.pt() - sjet.pt(),
					sd_info_signal_jet.dR, sd_info_j.dR, sd_info_j.dR - sd_info_signal_jet.dR,
					sd_info_signal_jet.z, sd_info_j.z, sd_info_j.z - sd_info_signal_jet.z)
			hpt.Fill(sjet.pt())
			hptcs.Fill(j.pt())
			hdpt.Fill(j.pt() - sjet.pt())
			hrg.Fill(sd_info_signal_jet.dR)
			hrgcs.Fill(sd_info_j.dR)
			if sd_info_j.dR > 0 and sd_info_signal_jet.dR > 0:
				hdrg.Fill(sd_info_j.dR - sd_info_signal_jet.dR)

		# for i,j in enumerate(signal_jets):
		# 	j_sd = sd.result(j)
		# 	sd_info = fjcontrib.get_SD_jet_info(j_sd)
		# 	# print("  |-> SD jet params z={0:10.3f} dR={1:10.3f} mu={2:10.3f}".format(sd_info.z, sd_info.dR, sd_info.mu))

	pythia.stat()
	outf.Write()
	outf.Close()


if __name__ == '__main__':
	main()
