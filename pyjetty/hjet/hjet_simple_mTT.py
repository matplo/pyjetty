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

	def fill_branches(self, t, ev_id = None, weight=None):
		if ev_id:
			t.fill_branch('V0_ev_id', ev_id)
		if weight:
			t.fill_branch('V0_w', weight)
		t.fill_branch('mV0A', self.V0A_mult)
		t.fill_branch('mV0C', self.V0C_mult)
		t.fill_branch('mV0', self.V0_mult)


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

	def fill_branches(self, t, ev_id = None, weight = None):
			if ev_id:
				t.fill_branch('jet_ev_id', ev_id)
			if weight:
				t.fill_branch('jet_w', weight)
			t.fill_branch('rho', self.rho)
			t.fill_branch('jet', self.jets)
			t.fill_branch('jet_ptcorr', self.corr_jet_pt)

class HJet(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(trigger_range = [6, 7], particles=None)
		super(HJet, self).__init__(**kwargs)
		if self.particles:
			self.analyze_event(self.particles)

	def analyze_event(self, jet_parts):
		self.trigger_particle = None
		t_selector = fj.SelectorPtRange(self.trigger_range[0], self.trigger_range[1])
		t_candidates = t_selector(jet_parts)
		if len(t_candidates) < 1:
			return False
		self.trigger_particle = random.choice(t_candidates)
		return True

	def fill_branches(self, t, jets=None, ev_id = None, weight = None):
		if self.trigger_particle:
			if jets:
				self.dphis = [self.trigger_particle.delta_phi_to(j) for j in jets]
				self.detas = [self.trigger_particle.eta() - j.eta() for j in jets]
			else:
				self.dphis = []
				self.detas = []
			if ev_id:
				t.fill_branch('t_ev_id', ev_id)
			if weight:
				t.fill_branch('t_w', weight)			
			t.fill_branch('t',  self.trigger_particle)
			t.fill_branch('t_dphi', self.dphis)
			t.fill_branch('t_deta', self.detas)

###########
# trigger_ranges = [ [0, 1e3], [6, 7] , [12, 22], [20, 30] ],
# jet_Rs = [0.2, 0.4, 0.6]

class HJetSet(object):
	def __init__(self, hjet, hjet_output):
		self.hjet = hjet
		self.hjet_output = hjet_output

	def analyze_event(self, jet_parts):
		self.hjet.analyze_event(jet_parts)
		if self.hjet.trigger_particle:
			return True
		return False

def generate():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	parser.add_argument('-d', '--output-dir', help='output directory name - use with default file name', default='.', type=str)
	parser.add_argument('-o', '--output', help='output file name', default='', type=str)	
	parser.add_argument('--overwrite', default=False, action='store_true')
	parser.add_argument('--jet-R', help='jet finder R', default=0.4, type=float)
	parser.add_argument('--charged', default=False, action='store_true')
	parser.add_argument('--runid', default=0, type=int)
	parser.add_argument('--tranges', default='6-7', help='hadron trigger ranges min-max,min1-max1,...', type=str)
	group = parser.add_mutually_exclusive_group()
	group.add_argument('--inel', default=False, action='store_true')
	group.add_argument('--hard', default=False, action='store_true')
	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()	

	# mycfg = ['PhaseSpace:pThatMin = 6']
	mycfg = []

	if '.root' not in args.output:
		rs = '0{}'.format(int(args.jet_R*10))
		if args.jet_R >= 1:
			rs = '{}'.format(int(args.jet_R))
		args.output = '{}/h_jet_R{}'.format(args.output_dir, rs)
		if args.charged:
			args.output = '{}/h_jet_ch_R{}'.format(args.output_dir, rs)
		if args.tranges:
			args.output += '_tranges_{}'.format(args.tranges.replace(',', "_"))
		if args.runid >= 0:
			args.output += '_runid_{}'.format(args.runid)
			args.py_seed = 1000 + args.runid
		if args.py_noMPI:
			args.output += '_noMPI'
		if args.py_noISR:
			args.output += '_noISR'
		if args.py_noHadron:
			args.output += '_noHadr'
		if args.py_minbias > 0:
			args.output += '_minbias'
		if args.py_pthatmin > 0:
			args.output += '_pthatmin_{}'.format(args.py_pthatmin)
		if args.inel:
			args.output += '_inel'
			mycfg.append("SoftQCD:inelastic = on") # Andreas' recommendation
			mycfg.append("HardQCD:all = off") # Andreas' recommendation
		print (args)
		if args.hard or ((args.py_minbias==False) and (args.inel==False)):
			args.output += '_biasref_{}'.format(args.py_biasref)
			args.output += '_hard'
			mycfg.append("HardQCD:all = on") # Andreas' recommendation
		args.output += '.root'

	if os.path.exists(args.output):
		if not args.overwrite:
			print('[w] output file', args.output, 'exists - skipping.')
			return

	print('[i] output file', args.output)

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()

	#and in the code i do not allow decay of the following particles:
	mycfg.append("310:mayDecay  = off") # //K0s
	mycfg.append("3122:mayDecay = off") # //labda0
	mycfg.append("3112:mayDecay = off") # //sigma-
	mycfg.append("3212:mayDecay = off") # //sigma0
	mycfg.append("3222:mayDecay = off") # //sigma+
	mycfg.append("3312:mayDecay = off") # //xi-
	mycfg.append("3322:mayDecay = off") # //xi+
	mycfg.append("3334:mayDecay = off") # //omega-
	print('[i] additional settings', mycfg)	
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		return
	if args.nev < 1:
		args.nev = 1

	tpairs = []
	st = args.tranges.split(',')
	try:
		for _s in st:
			t0 = float(_s.split('-')[0])
			t1 = float(_s.split('-')[1])
			tpairs.append([t0,t1])
	except:
		print('[e] something is not quite right with trigger ranges... stop here.')

	jet_particle_eta_max = 0.9
	jet_particle_selector = fj.SelectorAbsEtaMax(jet_particle_eta_max)

	v0det = V0Detector()
	j_ana = JetAnalysis(jet_R=args.jet_R, particle_eta_max=jet_particle_eta_max)

	fout = r.TFile(args.output, 'RECREATE')
	fout.cd()
	hmV0M = r.TH1F('hmV0M', 'hmV0M', 1000, 0, 1000)
	hmV0A = r.TH1F('hmV0A', 'hmV0A', 1000, 0, 1000)
	hmV0C = r.TH1F('hmV0C', 'hmV0C', 1000, 0, 1000)
	event_output = RTreeWriter(tree_name='evT', fout=fout)
	jet_output = RTreeWriter(tree_name='jetT', fout=fout)

	hjet_sets = []
	for t in tpairs:
		hjet_output = RTreeWriter(tree_name='hjetT_{}_{}'.format(int(t[0]), int(t[1])), fout=fout)
		hjet = HJet(trigger_range=[t[0], t[1]])
		print(hjet)
		hjset = HJetSet(hjet, hjet_output)
		hjet_sets.append(hjset)

	for ev_id in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue

		if args.charged:
			parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged])
		else:
			parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal])

		event_output.clear()
		jet_output.clear()
		_tmp = [s.hjet_output.clear() for s in hjet_sets]

		jet_parts = jet_particle_selector(parts)

		ev_s = pythia.info.sigmaGen()
		ev_w = pythia.info.weight()
		ev_code = pythia.info.code()
		pthard = pythia.info.pTHat()
		mTot = len(parts)
		mCB  = len(jet_parts)

		v0det.analyze_event(parts)
		event_output.fill_branch('ev_id', ev_id)	
		event_output.fill_branch('weight', ev_w)
		event_output.fill_branch('sigma', ev_s)
		event_output.fill_branch('code', ev_code)
		event_output.fill_branch('pthard', pthard)
		event_output.fill_branch('mTot', mTot)
		event_output.fill_branch('mCB', mCB)
		v0det.fill_branches(event_output)
		event_output.fill_tree()

		hmV0M.Fill(v0det.V0_mult)
		hmV0A.Fill(v0det.V0A_mult)
		hmV0C.Fill(v0det.V0C_mult)

		if len(parts) < 1:
			continue

		if True in [s.analyze_event(jet_parts) for s in hjet_sets]:
			j_ana.analyze_event(jet_parts)
			jet_output.fill_branch('ev_id', ev_id)
			jet_output.fill_branch('weight', ev_w)
			jet_output.fill_branch('sigma', ev_s)
			jet_output.fill_branch('code', ev_code)
			jet_output.fill_branch('pthard', pthard)
			jet_output.fill_branch('mTot', mTot)
			jet_output.fill_branch('mCB', mCB)
			j_ana.fill_branches(jet_output)
			jet_output.fill_tree()

			for s in hjet_sets:
				if s.hjet.trigger_particle:
					s.hjet_output.fill_branch('ev_id', ev_id)
					s.hjet_output.fill_branch('weight', ev_w)
					s.hjet_output.fill_branch('sigma', ev_s)
					s.hjet_output.fill_branch('code', ev_code)
					s.hjet_output.fill_branch('pthard', pthard)
					s.hjet_output.fill_branch('mTot', mTot)
					s.hjet_output.fill_branch('mCB', mCB)
					s.hjet.fill_branches(s.hjet_output, j_ana.jets)
					v0det.fill_branches(s.hjet_output)
					j_ana.fill_branches(s.hjet_output)
					s.hjet_output.fill_tree()

	pythia.stat()

	fout.Write()
	fout.Close()
	print('[i] written', fout.GetName())

if __name__ == '__main__':
	generate()
