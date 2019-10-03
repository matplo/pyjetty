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
import random

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

	def fill_branches(self, t, ev_id, weight=None):
		t.fill_branch('V0w', weight)
		t.fill_branch('mV0A', v0det.V0A_mult)
		t.fill_branch('mV0C', v0det.V0C_mult)
		t.fill_branch('mV0', v0det.V0_mult)


class JetAnalysis(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(jet_R=0.4, jet_algorithm=fj.antikt_algorithm, particle_eta_max=0.9, particles=None)
		super(JetAnalysis, self).__init__(**kwargs)

		self.bg_rho_range = fj.SelectorAbsEtaMax(self.particle_eta_max)
		self.bg_jet_def = fj.JetDefinition(fj.kt_algorithm, self.jet_R)
		self.bg_area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(self.particle_eta_max))
		self.bg_estimator = fj.JetMedianBackgroundEstimator(self.bg_rho_range, self.bg_jet_def, self.bg_area_def)
		self.rho = 0

		self.jet_eta_max = self.particle_eta_max - self.jet_R * 1.05
		self.jet_def = fj.JetDefinition(self.jet_algorithm, self.jet_R)
		self.jet_selector = fj.SelectorPtMin(0.0) & fj.SelectorAbsEtaMax(self.jet_eta_max)
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
			self.jets = fj.sorted_by_pt(self.jet_selector(self.cs.inclusive_jets()))
			self.corr_jet_pt = [j.pt() - j.area() * self.rho for j in self.jets]

	def fill_branches(self, t, ev_id, weight):
			if ev_id:
				t.fill_branch('ev_id', ev_id)
			if weight:
				t.fill_branch('jetw', weight)
			t.fill_branch('rho', self.rho)
			t.fill_branch('jet', self.jets)
			t.fill_branch('jet_ptcorr', self.corr_jet_pt)

class HJetTree(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(	fout=None,
									output='hjet_tree.root',
									trigger_range = [6, 7],
									jet_ana = None)
		super(HJetTree, self).__init__(**kwargs)
		if self.fout is None:
			self.fout = r.TFile(self.output, 'RECREATE')
		self.fout.cd()
		self.tree_name='tjet_0{}_{}_{}'.format(int(10*self.jet_ana.jet_R), int(self.trigger_range[0]), int(self.trigger_range[1]))
		self.tree_writer = RTreeWriter(tree_name=self.tree_name, fout=self.fout)

	def analyze_event(self, jet_parts):
		self.trigger_particle = None
		t_selector = fj.SelectorPtRange(self.trigger_range[0], self.trigger_range[1])
		t_candidates = t_selector(jet_parts)
		if len(t_candidates) < 1:
			return False
		self.trigger_particle = random.choice(t_candidates)
		return True

	def fill_branches(self, ev_id, weight):
		if self.jet_ana is None:
			print('[e] no jet ana in hjettree...')
			return False
		if self.trigger_particle:
			if self.trigger_particle:
				self.dphis = [self.trigger_particle.delta_phi_to(j) for j in self.jet_ana.jets]
				self.detas = [self.trigger_particle.eta() - j.eta() for j in self.jet_ana.jets]
			else:
				self.dphis = []
				self.detas = []
			self.tree_writer.fill_branch('ev_id', ev_id)
			self.tree_writer.fill_branch('t_w', weight)			
			self.tree_writer.fill_branch('t_pt',  self.trigger_particle.pt())
			self.tree_writer.fill_branch('t_phi', self.trigger_particle.phi())
			self.tree_writer.fill_branch('t_eta', self.trigger_particle.eta())
			self.tree_writer.fill_branch('dphi', self.dphis)
			self.tree_writer.fill_branch('deta', self.detas)
			self.jet_ana.fill_branches(self.tree_writer, ev_id = None, weight = None)

	def fill_tree(self):
		self.tree_writer.fill_tree()

	def write_and_close_file(self):
		self.fout.Write()
		self.fout.Close()


class HJetAnalysis(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(jet_particle_eta_max = 0.9,
								 output='hjet.root', 
								 trigger_ranges = [ [0, 1e3], [6, 7] , [12, 22], [20, 30] ],
								 jet_Rs = [0.2, 0.4, 0.6])
		super(HJetAnalysis, self).__init__(**kwargs)
		self.fout = r.TFile(self.output, 'RECREATE')
		self.fout.cd()
		self.event_output = RTreeWriter(tree_name='tev', fout=self.fout)

		self.hmV0M = r.TH1F('hmV0M', 'hmV0M', 1000, 0, 1000)
		self.hmV0A = r.TH1F('hmV0A', 'hmV0A', 1000, 0, 1000)
		self.hmV0C = r.TH1F('hmV0C', 'hmV0C', 1000, 0, 1000)

		self.jet_particle_selector = fj.SelectorAbsEtaMax(self.jet_particle_eta_max)

		self.jet_ans = []
		self.hjet_ts = []
		for jR in self.jet_Rs:
			j_ana = JetAnalysis(jet_R=jR, particle_eta_max=self.jet_particle_eta_max)
			self.jet_ans.append(j_ana)
			for tr in self.trigger_ranges:
				hjet_tree = HJetTree(fout=self.fout, trigger_range=tr, jet_ana=j_ana)
				self.hjet_ts.append(hjet_tree)

	def analyze_event(self, ev_id, parts, pythia=None):
		# V0 multiplicity and event props
		v0det = V0Detector(particles=parts)
		self.hmV0M.Fill(v0det.V0_mult)
		self.hmV0A.Fill(v0det.V0A_mult)
		self.hmV0C.Fill(v0det.V0C_mult)
		self.event_output.fill_branch('ev_id', ev_id)
		self.event_output.fill_branch('mV0A', v0det.V0A_mult)
		self.event_output.fill_branch('mV0C', v0det.V0C_mult)
		self.event_output.fill_branch('mV0', v0det.V0_mult)
		if pythia:
			ev_w = pythia.info.sigmaGen()
			pthard = pythia.info.pTHat()
		else:
			ev_w = -1
			pthard = -1
		self.event_output.fill_branch('weight', ev_w)
		self.event_output.fill_branch('pthard', pthard)
		jet_parts = self.jet_particle_selector(parts)
		mTot = len(parts)
		mCB  = len(jet_parts)
		self.event_output.fill_branch('mTot', mTot)
		self.event_output.fill_branch('mCB', mCB)
		self.event_output.fill_tree()

		trigger_present = False
		for hj in self.hjet_ts:
			if hj.analyze_event(jet_parts):
				trigger_present = True

		if trigger_present:
			for ja in self.jet_ans:
				ja.analyze_event(jet_parts)
			for hj in self.hjet_ts:
				if hj.trigger_particle:
					hj.fill_branches(ev_id, ev_w)
					hj.tree_writer.fill_branch('ev_id', ev_id)
					hj.tree_writer.fill_branch('mV0A', v0det.V0A_mult)
					hj.tree_writer.fill_branch('mV0C', v0det.V0C_mult)
					hj.tree_writer.fill_branch('mV0', v0det.V0_mult)
					hj.tree_writer.fill_branch('mTot', mTot)
					hj.tree_writer.fill_branch('mCB', mCB)
					hj.tree_writer.fill_branch('weight', ev_w)
					hj.tree_writer.fill_branch('pthard', pthard)
					hj.fill_tree()

	def finalize(self):
		self.fout.Write()
		self.fout.Close()
		print('[i] written', self.fout.GetName())

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
	if args.nev < 1:
		args.nev = 1

	# jet_particle_eta_max = 0.9
	hjet = HJetAnalysis(jet_particle_eta_max=0.9, output=args.output)

	for iev in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged])
		if len(parts) < 1:
			continue
		hjet.analyze_event(iev, parts, pythia)

	hjet.finalize()

	pythia.stat()

if __name__ == '__main__':
	generate()
