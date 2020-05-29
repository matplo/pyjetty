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
from pyjetty.mputils import ColorS, pwarning, perror, pinfo, pdebug

from pyjetty.alice_analysis.process.base import thermal_generator as thg

import ROOT
ROOT.gROOT.SetBatch(True)

class GroomerOutput(object):
	def __init__(self, outputname, enable_aa_trees=False):
		self.outf = ROOT.TFile(outputname, 'recreate')
		self.outf.cd()
		self.t = ROOT.TTree('t', 't')
		self.tw = RTreeWriter(tree=self.t)
		self.taaw = None
		self.tembw = None
		self.tembwp = None
		if enable_aa_trees:
			self.taa = ROOT.TTree('taa', 'taa')
			self.taaw = RTreeWriter(tree=self.taa)
			self.temb = ROOT.TTree('temb', 'temb')
			self.tembw = RTreeWriter(tree=self.temb)
			self.tembp = ROOT.TTree('tembp', 'tembp')
			self.tembwp = RTreeWriter(tree=self.tembp)

		self.jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
		self.dy_groomer = fjcontrib.DynamicalGroomer(self.jet_def_lund)

		self.alphas = [0.1, 1.0, 2.0]

		print (self.dy_groomer.description())

		self.sds = []
		sd01 = fjcontrib.SoftDrop(0, 0.1, 1.0)
		sd02 = fjcontrib.SoftDrop(0, 0.2, 1.0)
		self.sds.append(sd01)
		self.sds.append(sd02)

	def write(self):
		if self.outf:
			self.outf.Write()
			self.outf.Close()
			pinfo('written', self.outf.GetName())

	def __del__(self):
		self.outf = None

	def fill_branches(self, j, syst, **kwargs):
		sprefix = ''
		if syst == 0:
			tw = self.tw
		if syst == 1:
			tw = self.taaw
		if syst == 2:
			tw = self.tembw
		if syst == 3:
			tw = self.tembw
			sprefix = 'e'
		tw.fill_branch(sprefix+'j', j)
		for a in self.alphas:
			dy_groomed = self.dy_groomer.result(j, a)
			tw.fill_branch(sprefix+'dg_{:.1f}'.format(a), dy_groomed)
		max_pt_groomed = self.dy_groomer.max_pt_softer(j)
		tw.fill_branch(sprefix+'max_ptsoft', max_pt_groomed)
		max_z_groomed = self.dy_groomer.max_z(j)
		tw.fill_branch(sprefix+'max_z', max_z_groomed)
		max_kt_groomed = self.dy_groomer.max_kt(j)
		tw.fill_branch(sprefix+'max_kt', max_kt_groomed)
		max_kappa_groomed = self.dy_groomer.max_kappa(j)
		tw.fill_branch(sprefix+'max_kappa', max_kappa_groomed)
		max_tf_groomed = self.dy_groomer.max_tf(j)
		tw.fill_branch(sprefix+'max_tf', max_tf_groomed)
		min_tf_groomed = self.dy_groomer.min_tf(j)
		tw.fill_branch(sprefix+'min_tf', min_tf_groomed)

		for i,sd in enumerate(self.sds):
			j_sd = sd.result(j)
			tw.fill_branch(sprefix+'sd{}'.format(i), j_sd)
			sd_info = fjcontrib.get_SD_jet_info(j_sd)
			tw.fill_branch(sprefix+'sd{}_z'.format(i), sd_info.z)
			tw.fill_branch(sprefix+'sd{}_Delta'.format(i), sd_info.dR)
			tw.fill_branch(sprefix+'sd{}_mu'.format(i), sd_info.mu)
			tw.fill_branch(sprefix+'sd{}_kt'.format(i), sd_info.z * j_sd.pt() * sd_info.dR)

		for s in kwargs:
			tw.fill_branch(sprefix+s, kwargs[s])

		if syst == 2:
			pass
		else:
			tw.fill_tree()

	def fill_groomers_matched(self, tw, grd_pp, grd_emb, sdescr):
		tw.fill_branch(sdescr+'_pp', grd_pp)
		tw.fill_branch(sdescr+'_emb', grd_emb)
		tw.fill_branch(sdescr+'_dRppemb', grd_emb.pair().delta_R(grd_pp.pair()))
		harder_pt_match = -1
		softer_pt_match = -1
		if grd_pp.pair().perp() > 0 and grd_emb.pair().perp() > 0:
			harder_pt_match = fjtools.matched_pt(grd_emb.harder(), grd_pp.harder())
			softer_pt_match = fjtools.matched_pt(grd_emb.softer(), grd_pp.softer())
		tw.fill_branch(sdescr+'_mpt1', harder_pt_match)
		tw.fill_branch(sdescr+'_mpt2', softer_pt_match)

	def fill_branches_prong_matching(self, j_pp, j_emb, **kwargs):
		if fjtools.matched_pt(j_emb, j_pp) < 0.5:
			return
		tw = self.tembwp
		tw.fill_branch('j_pp', j_pp)
		for a in self.alphas:
			dy_groomed_pp = self.dy_groomer.result(j_pp, a)
			dy_groomed_emb = self.dy_groomer.result(j_emb, a)
			self.fill_groomers_matched(tw, self.dy_groomer.result(j_pp, a), self.dy_groomer.result(j_emb, a), 'dg_{:.1f}'.format(a))

		self.fill_groomers_matched(tw, self.dy_groomer.max_pt_softer(j_pp), self.dy_groomer.max_pt_softer(j_emb), 'max_ptsoft')
		self.fill_groomers_matched(tw, self.dy_groomer.max_z(j_pp), self.dy_groomer.max_z(j_emb), 'max_z')
		self.fill_groomers_matched(tw, self.dy_groomer.max_kt(j_pp), self.dy_groomer.max_kt(j_emb), 'max_kt')
		self.fill_groomers_matched(tw, self.dy_groomer.max_kappa(j_pp), self.dy_groomer.max_kappa(j_emb), 'max_kappa')
		self.fill_groomers_matched(tw, self.dy_groomer.max_tf(j_pp), self.dy_groomer.max_tf(j_emb), 'max_tf')
		self.fill_groomers_matched(tw, self.dy_groomer.min_tf(j_pp), self.dy_groomer.min_tf(j_emb), 'min_tf')

		for i,sd in enumerate(self.sds):
			j_sd_pp = sd.result(j_pp)
			tw.fill_branch('sd{}_pp'.format(i), j_sd_pp)
			sd_info = fjcontrib.get_SD_jet_info(j_sd_pp)
			tw.fill_branch('sd{}_z_pp'.format(i), sd_info.z)
			tw.fill_branch('sd{}_Delta_pp'.format(i), sd_info.dR)
			tw.fill_branch('sd{}_mu_pp'.format(i), sd_info.mu)
			tw.fill_branch('sd{}_kt_pp'.format(i), sd_info.z * j_sd_pp.pt() * sd_info.dR)

			j_sd_emb = sd.result(j_emb)
			tw.fill_branch('sd{}_emb'.format(i), j_sd_emb)
			sd_info = fjcontrib.get_SD_jet_info(j_sd_emb)
			tw.fill_branch('sd{}_z_emb'.format(i), sd_info.z)
			tw.fill_branch('sd{}_Delta_emb'.format(i), sd_info.dR)
			tw.fill_branch('sd{}_mu_emb'.format(i), sd_info.mu)
			tw.fill_branch('sd{}_kt_emb'.format(i), sd_info.z * j_sd_pp.pt() * sd_info.dR)

			# now fill the prong matching
			p1 = fj.PseudoJet()
			p2 = fj.PseudoJet()
			has_parents_pp = j_sd_pp.has_parents(p1, p2)
			tw.fill_branch('sd{}_p1_pp'.format(i), p1)
			tw.fill_branch('sd{}_p2_pp'.format(i), p2)

			pe1 = fj.PseudoJet()
			pe2 = fj.PseudoJet()
			has_parents_emb = j_sd_emb.has_parents(pe1, pe2)
			tw.fill_branch('sd{}_p1_emb'.format(i), pe1)
			tw.fill_branch('sd{}_p2_emb'.format(i), pe2)

			# if has_parents_emb:
			# 	tw.fill_branch('ej_p1_ptc', pe1.pt() - pe1.area() * rho)
			# 	tw.fill_branch('ej_p2_ptc', pe2.pt() - pe2.area() * rho)
			# else:
			# 	tw.fill_branch('ej_p1_ptc', -1000)
			# 	tw.fill_branch('ej_p2_ptc', -1000)

			mpt1 = -1.0 # not passed SD
			mpt2 = -1.0 # not passed SD

			if has_parents_pp and has_parents_emb:
				mpt1 = fjtools.matched_pt(pe1, p1)
				mpt2 = fjtools.matched_pt(pe2, p2)
			tw.fill_branch('sd{}_mpt1'.format(i), mpt1)
			tw.fill_branch('sd{}_mpt2'.format(i), mpt2)


		for s in kwargs:
			tw.fill_branch(s, kwargs[s])

		tw.fill_tree()


npart_cents='''
90pc@6,
80pc@15,
70pc@31,
60pc@56,
50pc@90,
40pc@133,
30pc@186,
20pc@248,
10pc@325,
'''

def main():
	parser = argparse.ArgumentParser(description='test groomers', prog=os.path.basename(__file__))
	parser.add_argument('-o', '--output-filename', default="output.root", type=str)
	parser.add_argument('datalistpp', help='run through a file list', default='', type=str)
	parser.add_argument('--datalistAA', help='run through a file list - embedding mode', default='', type=str)
	parser.add_argument('--jetR', default=0.4, type=float)
	parser.add_argument('--alpha', default=0, type=float)
	parser.add_argument('--dRmax', default=0.25, type=float)
	parser.add_argument('--overwrite', help="overwrite output", default=False, action='store_true')
	parser.add_argument('--jetptcut', help='remove jets below the cut', default=50., type=float)
	parser.add_argument('--nev', help='number of events to run', default=0, type=int)
	parser.add_argument('--max-eta', help='max eta for particles', default=0.9)
	parser.add_argument('--npart-cut', help='npart cut on centrality low,high hint:' + npart_cents, default='325,450', type=str)

	args = parser.parse_args()

	try:
		npart_min = int(args.npart_cut.split(',')[0])
		npart_max = int(args.npart_cut.split(',')[1])
	except:
		perror('unable to parse npart centrality selection - two integer numbers with a coma in-between needed - specified:', args.npart_cut)
		return 1

	# initialize constituent subtractor
	cs = None
	if args.dRmax > 0:
		cs = CEventSubtractor(	alpha=args.alpha, max_distance=args.dRmax, max_eta=args.max_eta, 
								bge_rho_grid_size=0.25, max_pt_correct=100)

	pp_data 	= DataIO(	name='Sim Pythia Detector level', file_list=args.datalistpp, random_file_order=False, tree_name='tree_Particle_gen')
	ja_pp 	= JetAnalysis(jet_R=args.jetR, jet_algorithm=fj.antikt_algorithm, jet_pt_min=50., particle_eta_max=args.max_eta)

	if args.datalistAA:
		aa_data = DataBackgroundIO(	name='PbPb', file_list=args.datalistAA, tree_name='tree_Particle_gen')
		ja_emb 	= JetAnalysis(	jet_R=args.jetR, jet_algorithm=fj.antikt_algorithm, jet_pt_min=50., particle_eta_max=args.max_eta)
		ja_aa 	= JetAnalysis(	jet_R=args.jetR, jet_algorithm=fj.antikt_algorithm, jet_pt_min=50., particle_eta_max=args.max_eta)


	dndeta_selector = fj.SelectorAbsEtaMax(1.)

	# tg = thg.ThermalGenerator()
	print(cs)

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()

	gout = GroomerOutput(args.output_filename, enable_aa_trees=bool(args.datalistAA))

	delta_t = 0
	start_t = time.time()
	iev = 1
	while pp_data.load_event(offset=0):
		iev = iev + 1
		if args.nev > 0:
			if iev > args.nev:
				iev = iev - 1
				break
		if iev % 1000 == 0:
			delta_t = time.time() - start_t
			pinfo('processing event', iev, ' - ev/sec =', iev/delta_t, 'elapsed =', delta_t)

		# find jets on detector level
		if len(pp_data.particles) < 1:
			pwarning(iev, 'pp event skipped N parts', len(pp_data.particles))
			continue
		ja_pp.analyze_event(pp_data.particles)
		if len(ja_pp.jets) < 1:
			continue

		# pinfo('n particles', len(pp_data.particles))
		dndeta0 = dndeta_selector(pp_data.particles)
		[gout.fill_branches(j, syst=0, dndeta=len(dndeta0)/2.) for j in ja_pp.jets]
		# pinfo('n jets', len(ja_pp.jets))

		if args.datalistAA:
			while True:
				aa_loaded = aa_data.load_event(offset=10000)
				if aa_data.event.npart < npart_min or aa_data.event.npart >= npart_max:
					continue
				else:
					if len(aa_data.particles) < 1:
						pwarning(iev, 'AA event skipped N parts', len(aa_data.particles))
						continue
					else:
						break
			if aa_loaded:
				ja_aa.analyze_event(aa_data.particles)
				dndeta1 = dndeta_selector(aa_data.particles)
				if len(ja_aa.jets) > 0:
					[gout.fill_branches(j, syst=1, dndeta=len(dndeta1)/2.) for j in ja_aa.jets]
				else:
					# pwarning('no jets in AA event?', len(ja_aa.jets), 'while dndeta=', len(dndeta1)/2.)
					pass
				emb_event = fj.vectorPJ()
				[emb_event.push_back(p) for p in pp_data.particles]
				[emb_event.push_back(p) for p in aa_data.particles]
				rho = 0
				if cs:
					cs_parts = cs.process_event(emb_event)
					rho = cs.bge_rho.rho()
					ja_emb.analyze_event(cs_parts)
				else:
					ja_emb.analyze_event(emb_event)
				# matches = [[jpp, jemb] for jpp in ja_pp.jets for jemb in ja_emb.jets if fjtools.matched_pt(jemb, jpp) > 0.5]
				# for mj in matches:
				# 	gout.fill_branches(mj[0], syst=2, dndeta=len(dndeta1)/2., rho=rho)
				# 	gout.fill_branches(mj[1], syst=3)
				[gout.fill_branches_prong_matching(j_pp, j_emb, dndeta=len(dndeta1)/2., rho=rho) for j_pp in ja_pp.jets for j_emb in ja_emb.jets]


	delta_t = time.time()-start_t
	pinfo('processed events', iev, ' - ev/sec =', iev/delta_t, 'elapsed =', delta_t)
	gout.write()

if __name__ == '__main__':
	main()
