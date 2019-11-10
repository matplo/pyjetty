#!/usr/bin/env python

# this is modified csdata.py

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext
import fjtools

import tqdm
import argparse
import os
import numpy as np
import array
import copy
import random
import uproot
import pandas as pd
import time

from pyjetty.mputils import logbins
from pyjetty.mputils import MPBase
from pyjetty.mputils import BoltzmannEvent
from pyjetty.mputils import CEventSubtractor
from pyjetty.mputils import RTreeWriter
from pyjetty.mputils import DataIO, DataBackgroundIO
from pyjetty.mputils import fill_tree_data, JetAnalysis, JetAnalysisWithRho
from pyjetty.mputils import ColorS, pwarning, perror, pinfo

from alice_efficiency import AliceChargedParticleEfficiency

import ROOT
ROOT.gROOT.SetBatch(True)

# make it a class
class Embedding(MPBase):
	def add_arguments_to_parser(parser):
		parser.add_argument('--output', default="output.root", type=str)
		parser.add_argument('datalistAA', help='run through a file list', default='', type=str)
		parser.add_argument('simulationpp', help='run through a file list', default='', type=str)
		parser.add_argument('--jetR', default=0.4, type=float)
		parser.add_argument('--alpha', default=0, type=float)
		parser.add_argument('--dRmax', default=0.25, type=float)
		parser.add_argument('--sd-zcut', default=0.1, type=float)
		parser.add_argument('--overwrite', help="overwrite output", default=False, action='store_true')
		parser.add_argument('--benchmark', help='benchmark pthat setting - 80 GeV', default=False, action='store_true')
		parser.add_argument('--jetptcut', help='remove jets below the cut', default=1.e-3, type=float)
		parser.add_argument('--nev', help='number of events to run', default=0, type=int)
		parser.add_argument('--max-eta', help='max eta for particles', default=0.9)

	def __init__(self, **kwargs):
		self.configure_from_args(tree_name='tree_Particle', tree_name_gen='tree_Particle_gen')
		super(Embedding, self).__init__(**kwargs)
		self.copy_attributes(self.args)
		self.jet_def = fj.JetDefinition(fj.antikt_algorithm, self.jetR)
		if self.benchmark:
			self.jet_selector = fj.SelectorPtMin(80.0) & fj.SelectorPtMax(100.0) & fj.SelectorAbsEtaMax(self.max_eta - 1.05 * self.jetR)
			# jet_selector_cs = fj.SelectorPtMin(50.0) & fj.SelectorAbsEtaMax(max_eta - 1.05 * self.jetR)
		else:
			self.jet_selector = fj.SelectorAbsEtaMax(self.max_eta - 1.05 * self.jetR)

		self.sd = fjcontrib.SoftDrop(0, self.sd_zcut, self.jetR)

		self.ja_part 	= JetAnalysis(	jet_R=self.jetR, jet_algorithm=fj.antikt_algorithm, 
										jet_pt_min=self.jetptcut/2., particle_eta_max=self.max_eta)
		self.ja_det 	= JetAnalysis(	jet_R=self.jetR, jet_algorithm=fj.antikt_algorithm, 
										jet_pt_min=self.jetptcut, particle_eta_max=self.max_eta)
		self.ja_hybrid 	= JetAnalysis(	jet_R=self.jetR, jet_algorithm=fj.antikt_algorithm, 
										jet_pt_min=self.jetptcut, particle_eta_max=self.max_eta)

		self.dataPbPb 	= DataBackgroundIO(	name='Data PbPb', file_list=self.datalistAA)
		self.det_sim 	= DataIO(	name='Sim Pythia Detector level', file_list=self.simulationpp, random_file_order=False)
		self.part_sim 	= DataIO(	name='Sim Pythia Particle level', file_list=self.simulationpp, random_file_order=False, 
									tree_name='tree_Particle_gen')

		self.cs = None
		if self.dRmax > 0:
			self.cs = CEventSubtractor(	alpha=self.alpha, max_distance=self.dRmax, max_eta=self.max_eta, 
										bge_rho_grid_size=0.25, max_pt_correct=100)
		self.parts_selector = fj.SelectorAbsEtaMax(self.max_eta)

		self.outf = ROOT.TFile(self.output, 'recreate')
		self.outf.cd()
		self.t 		= ROOT.TTree('t', 't')
		self.tw 	= RTreeWriter(tree=self.t, name='Output Tree')

		self.tpp 	= ROOT.TTree('tpp', 'tpp')
		self.twpp 	= RTreeWriter(tree=self.tpp, name='Output Tree pp simulation')

	def run(self):
		# need to change this for data to drive...
		delta_t = 0
		start_t = time.time()
		iev = 1
		# while self.det_sim.load_event() and self.part_sim.load_event():
		while self.det_sim.load_event():
			iev = iev + 1
			if self.nev > 0:
				if iev > self.nev:
					iev = iev - 1
					break
				if iev % 1000 == 0:
					delta_t = time.time()-start_t
					pinfo('processing event', iev, ' - ev/sec =', iev/delta_t, 'elapsed =', delta_t)

			self.part_sim.open_afile(afile=self.det_sim.file_io.file_input)
			if not self.part_sim.load_event_with_loc(self.det_sim.event.run_number, self.det_sim.event.ev_id, 0):
				perror('unable to load partL event run#:', self.det_sim.event.run_number, 'ev_id:', self.det_sim.event.ev_id)
				continue

			if self.det_sim.event.run_number != self.part_sim.event.run_number:
				perror('run# missmatch detL:', self.det_sim.event.run_number, 'partL:', self.part_sim.event.run_number)
				if self.det_sim.event.run_number > self.part_sim.event.run_number:
					pwarning('- trying to recover by advancing partL')
					self.part_sim.load_event()
					if self.det_sim.event.run_number != self.part_sim.event.run_number:
						pwarning('- failed detL:', self.det_sim.event.run_number, 'partL:', self.part_sim.event.run_number)
					else:
						pwarning('- looks good (but will skip an event anyways)', self.det_sim.event.run_number)
				else:
					pwarning('- trying to recover by advancing detL')
					self.det_sim.load_event()					
					if self.det_sim.event.run_number != self.part_sim.event.run_number:
						pwarning('- failed detL:', self.det_sim.event.run_number, 'partL:', self.part_sim.event.run_number)
					else:
						pwarning('- looks good (but will skip an event anyways)', self.det_sim.event.run_number)
				iev = iev - 1
				continue

			if self.det_sim.event.ev_id != self.part_sim.event.ev_id:
				perror('ev_id# missmatch detL:', self.det_sim.event.ev_id, 'partL:',self.part_sim.event.ev_id)
				if self.det_sim.event.ev_id > self.part_sim.event.ev_id:
					pwarning('- trying to recover by advancing partL')
					self.part_sim.load_event()
					if self.det_sim.event.ev_id != self.part_sim.event.ev_id:
						pwarning('- failed detL:', self.det_sim.event.ev_id, 'partL:', self.part_sim.event.ev_id)
					else:
						pwarning('- looks good (but will skip an event anyways)', self.det_sim.event.ev_id)
				else:
					pwarning('- trying to recover by advancing detL')
					self.det_sim.load_event()					
					if self.det_sim.event.ev_id != self.part_sim.event.ev_id:
						pwarning('- failed detL:', self.det_sim.event.ev_id, 'partL:', self.part_sim.event.ev_id)
					else:
						pwarning('- looks good (but will skip an event anyways)', self.det_sim.event.ev_id)
				iev = iev - 1
				continue

			if len(self.part_sim.particles) < 1 or len(self.det_sim.particles) < 1:
				pwarning(iev, 'event skipped N particle parts', len(self.part_sim.particles), 'N detector parts', len(self.det_sim.particles))
				continue

			self.ja_part.analyze_event(self.part_sim.particles)
			_jets_part = self.ja_part.jets
			self.ja_det.analyze_event(self.det_sim.particles)
			_jets_det = self.ja_det.jets

			if len(_jets_part) < 1 or len(_jets_det) < 1:
				# pwarning(iev, 'event skipped N particle jets', len(_jets_part), 'N detector jets', len(_jets_det))
				continue

			_det_part_matches = []
			_n_matches = 0
			for j_det in _jets_det:
				_part_psjv = self.ja_part.jets_as_psj_vector()
				_mactches_pp = fjtools.matched_Reta(j_det, _part_psjv, 0.6 * self.jetR)
				_n_matches = _n_matches + len(_mactches_pp)
				if len(_mactches_pp) > 1:
					pwarning('event:', iev, 'jet pt=', j_det.pt(), 'more than one match in pp jets', [i for i in _mactches_pp])
				if len(_mactches_pp) == 1:
					_det_part_matches.append([ j_det, _part_psjv[_mactches_pp[0]] ])

			if _n_matches < 1 or _n_matches > 1:
				if _n_matches < 1:
					pwarning('event:', iev, '- no matched jets in simulation!?', len(_det_part_matches))
				else:
					pass
			else:
				self.twpp.fill_branch('iev', iev)
				self.twpp.fill_branch('pt_det', [m[0].pt() for m in _det_part_matches])
				self.twpp.fill_branch('pt_part', [m[1].pt() for m in _det_part_matches])
				self.twpp.fill_branch('dpt', [m[0].pt() - m[1].pt() for m in _det_part_matches])
				self.twpp.fill_branch('dR', [m[0].delta_R(m[1]) for m in _det_part_matches])
				self.twpp.fill_tree()

			_offset = 10000
			while _offset < len(self.det_sim.particles):
				_offset = _offset + 1000
				pwarning('increasing bg index offset to', _offset)

			_PbPb_loaded = 0
			while _PbPb_loaded == 0:
				if not self.dataPbPb.load_event(offset=_offset):
					perror('unable to load next PbPb event')
					_PbPb_loaded = -1
				else:
					_hybrid_event = self.dataPbPb.particles
					_nparts_hybrid_no_emb = len(_hybrid_event)
					if len(_hybrid_event) < 1:
						pwarning('hybrid event with no particles! trying another one')
						_PbPb_loaded = 0
					else:
						_PbPb_loaded = 1
			if _PbPb_loaded < 0:
				perror('unable to load PbPb event - permanent - bailing out here.')
				break

			continue

			_tmp = [_hybrid_event.push_back(p) for p in self.det_sim.particles]
			# _nparts_hybrid = len(_hybrid_event)
			# print(_nparts_hybrid_no_emb, _nparts_hybrid)

			if self.cs:
				cs_parts = self.cs.process_event(_hybrid_event)
				rho = self.cs.bge_rho.rho()
				self.ja_hybrid.analyze_event(cs_parts)
			else:
				self.ja_hybrid.analyze_event(_hybrid_event)
			# tmp = [fill_tree_data(j, self.tw, self.sd, self.ja.rho, iev, 1.) for j in self.ja.jets if j.pt() > self.jetptcut]


		delta_t = time.time()-start_t
		pinfo('processed events', iev, ' - ev/sec =', iev/delta_t, 'elapsed =', delta_t)
		self.outf.Write()
		self.outf.Close()
		pinfo('written', self.outf.GetName())


def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	Embedding.add_arguments_to_parser(parser)
	args = parser.parse_args()

	if args.output == 'output.root':
		args.output = 'output_data_emb_CS_alpha_{}_dRmax_{}_SDzcut_{}.root'.format(args.alpha, args.dRmax, args.sd_zcut)
		if args.jetptcut > -100:
			args.output = 'output_data_emb_CS_alpha_{}_dRmax_{}_SDzcut_{}_jpt_{}.root'.format(args.alpha, args.dRmax, args.sd_zcut, args.jetptcut)

	if os.path.isfile(args.output):
		if not args.overwrite:
			print('[i] output', args.output, 'exists - use --overwrite to do just that...')
			return

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()

	embd = Embedding(args=args)
	print(embd)

	embd.run()


if __name__ == '__main__':
	main()
