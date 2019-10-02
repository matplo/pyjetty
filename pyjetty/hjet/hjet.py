#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

import ROOT as r
import array

from pyjetty.mputils import MPBase
from pyjetty.mputils import RTreeWriter

def logbins(xmin, xmax, nbins):
		lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
		arr = array.array('f', lspace)
		return arr

class V0Detector(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(particles = None)
		super(V0Detector, self).__init__(**kwargs)
		self.V0A_selector = fj.SelectorEtaMin(2.8) & fj.SelectorEtaMax(5.1)
		self.V0C_selector = fj.SelectorEtaMin(-3.7) & fj.SelectorEtaMax(-1.7)
		self.V0_selector = self.V0A_selector | self.V0C_selector
		self.V0A_mult = 0
		self.V0C_mult = 0
		self.V0_mult = 0
		if self.particles:
			self.analyze_event(self.particles)

	def analyze_event(self, parts):
		self.V0A_mult = len(self.V0A_selector(parts))
		self.V0C_mult = len(self.V0C_selector(parts))
		self.V0_mult = self.V0A_mult + self.V0C_mult
		self.V0_and = ((self.V0A_mult > 0) & (self.V0C_mult > 0))
		self.V0_or = ((self.V0A_mult > 0) | (self.V0C_mult > 0))

class JetAnalysis(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(jet_R=0.4, jet_algorithm=fj.antikt_algorithm, particle_eta_max=0.9, particles=None)
		super(JetAnalysis, self).__init__(**kwargs)

		self.bg_rho_range = fj.SelectorAbsRapMax(self.particle_eta_max)
		self.bg_jet_def = fj.JetDefinition(fj.kt_algorithm, self.jet_R)
		self.bg_area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(self.particle_eta_max))
		self.bg_estimator = fj.JetMedianBackgroundEstimator(self.bg_rho_range, self.bg_jet_def, self.bg_area_def)
		self.rho = 0

		self.jet_eta_max = self.particle_eta_max - self.jet_R * 1.05
		self.jet_def = fj.JetDefinition(self.jet_algorithm, self.jet_R)
		self.jet_selector = fj.SelectorPtMin(0.0) & fj.SelectorAbsRapMax(self.jet_eta_max)
		self.jet_area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(self.particle_eta_max))

		if self.particles:
			self.analyze_event(self.particles)

	def analyze_event(self, parts):
		if len(parts) < 1:
			self.rho = 0.0
			self.cs = None
			self.jets = []
			self.corr_jet_pt = []
		else:
			self.bg_estimator.set_particles(parts)
			self.rho = self.bg_estimator.rho()
			self.cs = fj.ClusterSequenceArea(parts, self.jet_def, self.jet_area_def)
			self.jets = fj.sorted_by_pt(self.cs.inclusive_jets())
			self.corr_jet_pt = [j.pt() - j.area() * self.rho for j in self.jets]

def HJetAnalysis(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(output='hjet.root', trigger_ranges= [[6, 7] , [12, 22], [20, 30] ])
		super(HJetAnalysis, self).__init__(**kwargs)
		fout = r.TFile(self.output, 'RECREATE')
		fout.cd()

def generate():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	parser.add_argument('-o', '--output', help='output file name', default='hjet.root', type=str)
	parser.add_argument('-t', '--trigger', help='trigger pt', default=6., type=float)
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--charged', default=False, action='store_true')
	args = parser.parse_args()	

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()

	# mycfg = ['PhaseSpace:pThatMin = 6']
	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 100:
		args.nev = 100

	fout = r.TFile(args.output, 'RECREATE')
	fout.cd()	
	tntriggers = r.TNtuple('triggers', 'triggers', 'eventid:mV0:hpt')
	tnhjet = r.TNtuple('hjet', 'hjet', 'eventid:mV0:hpt:jetpt:dphi:deta')
	jet_output = RTreeWriter(tree_name='tjet', fout=fout)
	event_output = RTreeWriter(tree_name='tev', fout=fout)

	hmV0M = r.TH1F('hmV0M', 'hmV0M', 1000, 0, 1000)
	hmV0A = r.TH1F('hmV0A', 'hmV0A', 1000, 0, 1000)
	hmV0C = r.TH1F('hmV0C', 'hmV0C', 1000, 0, 1000)

	v0det = V0Detector()
	hjet_a = JetAnalysis()
	jet_particle_eta_max = 0.9
	jet_particle_selector = fj.SelectorAbsRapMax(jet_particle_eta_max)

	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue

		parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged])
		jet_parts = jet_particle_selector(parts)

		jet_lead_particle = None
		jet_lead_particle_pt = 0.0
		if len(parts) < 1:
			continue

		jet_lead_particle = fj.sorted_by_pt(jet_parts)[0]
		jet_lead_particle_pt = jet_lead_particle.pt()

		v0det.analyze_event(parts)
		hjet_a.analyze_event(jet_parts)

		hmV0M.Fill(v0det.V0_mult)
		hmV0A.Fill(v0det.V0A_mult)
		hmV0C.Fill(v0det.V0C_mult)
		tntriggers.Fill(i, v0det.V0_mult, jet_lead_particle_pt)
		event_output.fill_branch('ev_id', i)
		event_output.fill_branch('mV0A', v0det.V0A_mult)
		event_output.fill_branch('mV0C', v0det.V0C_mult)
		event_output.fill_branch('mV0', v0det.V0_mult)
		event_output.fill_branch('lead_pt', jet_lead_particle_pt)
		event_output.fill_branch('sigmaGen', pythia.info.sigmaGen())
		event_output.fill_tree()

		jets = hjet_a.jets
		for i, j in enumerate(jets):
			if jet_lead_particle:
				dphi = jet_lead_particle.delta_phi_to(j)
				deta = jet_lead_particle.eta() - j.eta()
				tnhjet.Fill(i, v0det.V0_mult, jet_lead_particle_pt, j.pt(), dphi, deta)

		dphis = []
		detas = []
		if jet_lead_particle:
			dphis = [jet_lead_particle.delta_phi_to(j) for j in jets]
			detas = [jet_lead_particle.eta() - j.eta() for j in jets]

		jet_output.fill_branch('ev_id', i)
		jet_output.fill_branch('mV0A', v0det.V0A_mult)
		jet_output.fill_branch('mV0C', v0det.V0C_mult)
		jet_output.fill_branch('mV0', v0det.V0_mult)
		jet_output.fill_branch('lead_pt', jet_lead_particle_pt)
		jet_output.fill_branch('sigmaGen', pythia.info.sigmaGen())
		jet_output.fill_branch('jet', jets)
		jet_output.fill_branch('jet_ptcorr', hjet_a.corr_jet_pt)
		jet_output.fill_branch('dphi', dphis)
		jet_output.fill_branch('deta', detas)
		jet_output.fill_tree()

	pythia.stat()

	fout.Write()
	fout.Close()

if __name__ == '__main__':
	generate()
