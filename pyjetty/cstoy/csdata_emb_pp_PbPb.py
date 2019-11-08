#!/usr/bin/env python

# this is modified csdata.py

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

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
from pyjetty.mputils import DataIO
from pyjetty.mputils import fill_tree_data, JetAnalysisWithRho
from pyjetty.mputils import ColorS

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
		parser.add_argument('--jetptcut', help='remove jets below the cut', default=-100, type=float)
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
		self.ja = JetAnalysisWithRho(jet_R=self.jetR, jet_algorithm=fj.antikt_algorithm, particle_eta_max=self.max_eta)

		self.data = DataIO(name='Data PbPb', file_list=self.datalistAA)
		self.det_sim = DataIO(name='Sim Pythia Detector level', file_list=self.simulationpp)
		self.part_sim = DataIO(name='Sim Pythia Particle level', file_list=self.simulationpp, tree_name='tree_Particle_gen')

		self.cs = None
		if self.dRmax > 0:
			self.cs = CEventSubtractor(	alpha=self.alpha, max_distance=self.dRmax, max_eta=self.max_eta, 
										bge_rho_grid_size=0.25, max_pt_correct=100)
		self.parts_selector = fj.SelectorAbsEtaMax(self.max_eta)

		self.outf = ROOT.TFile(self.output, 'recreate')
		self.outf.cd()
		self.t = ROOT.TTree('t', 't')
		self.tw = RTreeWriter(tree=self.t, name='Output Tree')

	def run(self):
		# need to change this for data to drive...
		delta_t = 0
		start_t = time.time()
		iev = 0
		while self.data.load_event():
			iev = iev + 1
			_data_parts = self.data.particles
			if self.cs:
				cs_parts = self.cs.process_event(_data_parts)
				rho = self.cs.bge_rho.rho()
				self.ja.analyze_event(cs_parts)
			else:
				self.ja.analyze_event(_data_parts)
			tmp = [fill_tree_data(j, self.tw, self.sd, self.ja.rho, iev, 1.) for j in self.ja.jets if j.pt() > self.jetptcut]
			if iev % 1000 == 0:
				delta_t = time.time()-start_t
				print(ColorS.green('[i] processing event', iev, ' - ev/sec = ', iev/delta_t, 'elapsed = ', delta_t))
			if self.nev > 0:
				if iev > self.nev:
					break
		delta_t = time.time()-start_t
		print(ColorS.yellow('[i] processed events', iev, ' - ev/sec = ', iev/delta_t, 'elapsed = ', delta_t))
		self.outf.Write()
		self.outf.Close()
		print('[i] written', self.outf.GetName())


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
