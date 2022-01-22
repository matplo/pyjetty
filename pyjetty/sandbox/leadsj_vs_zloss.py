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

from pyjetty.mputils import MPBase, pwarning, pinfo, perror, treewriter, jet_analysis
from pyjetty.mputils import logbins

import ROOT

def match_dR(j, partons, drmatch = 0.1):
	mps = [p for p in partons if j.delta_R(p) < drmatch]
	# for p in fj.sorted_by_pt(mps)[0]:
	if len(mps) < 1:
		return None, False, False
	p = fj.sorted_by_pt(mps)[0]
	pyp = pythiafjext.getPythia8Particle(p)
	# print(p, pyp.id(), pyp.isQuark(), pyp.isGluon())
	return pyp.id(), pyp.isQuark(), pyp.isGluon()

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--nw', help="no warn", default=True, action='store_true')
	parser.add_argument('--ignore-mycfg', help="ignore some settings hardcoded here", default=False, action='store_true')
	parser.add_argument('--enable-background', help="enable background calc", default=False, action='store_true')
	parser.add_argument('--output', help="output file name", default='leadsj_vs_zloss.root', type=str)
	parser.add_argument('--jetptmin', help="jet pt minimum", default=-1, type=float)
	parser.add_argument('--jetptmax', help="jet pt maximum", default=1e6, type=float)
	parser.add_argument('--eta', help="jet eta max", default=2.4, type=float)

	args = parser.parse_args()

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(args.py_pthatmin) & fj.SelectorPtMax(1000.0) & fj.SelectorAbsEtaMax(args.eta - jet_R0)
	mycfg = []
	if args.jetptmin > 0:
		mycfg = ['PhaseSpace:pThatMin = {}'.format(args.jetptmin)]
		jet_selector = fj.SelectorPtMin(args.jetptmin) & fj.SelectorPtMax(args.jetptmax) &fj.SelectorAbsEtaMax(args.eta - jet_R0)
	print(jet_def)

	if args.ignore_mycfg:
		mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		perror("pythia initialization failed.")
		return

	nbins = 20
	# sjrs = [0.001 + x * 0.04 for x in range(0, nbins)]
	sjrs = logbins(0.001, jet_R0, nbins)
	print(sjrs)
	print('log(1/r) :', [ROOT.TMath.Log(1/r) for r in sjrs])
	sjdefs = dict()
	for sjr in sjrs:
		_jet_def = fj.JetDefinition(fj.antikt_algorithm, sjr)
		sjdefs[sjr] = _jet_def

	# tw = treewriter.RTreeWriter(name = 'lsjvsx', file_name = 'leadsj_vs_x.root')
	tw = treewriter.RTreeWriter(name = 'lsjvszloss', file_name = args.output)
	tw.fout.cd()
	h_zloss_r_q = dict()
	h_zloss_r_g = dict()
	for sjr in sjrs:
		sname = 'h_zloss_glue_{}'.format(sjr)
		_h_zloss_r_g = ROOT.TH1F(sname, sname, len(sjrs), 0., 1.)
		h_zloss_r_g[sjr] = _h_zloss_r_g
		sname = 'h_zloss_quark_{}'.format(sjr)
		_h_zloss_r_q = ROOT.TH1F(sname, sname, len(sjrs), 0., 1.)
		h_zloss_r_q[sjr] = _h_zloss_r_q

	sname = 'prof_zloss_vs_r_glue'
	prof_g = ROOT.TProfile(sname, sname, nbins, 0, jet_R0)
	lbins = logbins(ROOT.TMath.Log(1./jet_R0), ROOT.TMath.Log(1./sjrs[0]), nbins)
	print('lbins:', lbins)
	prof_g_log = ROOT.TProfile(sname+'_log', sname+'_log', nbins, lbins)
	sname = 'prof_zloss_vs_r_quark'
	prof_q = ROOT.TProfile(sname, sname, nbins, 0, jet_R0)
	prof_q_log = ROOT.TProfile(sname+'_log', sname+'_log', nbins, lbins)
	# prof_q_log = ROOT.TProfile(sname+'_log', sname+'_log', nbins, ROOT.TMath.Log(1./jet_R0), ROOT.TMath.Log(1./sjrs[0]))

	sname = 'h2_zloss_vs_r_glue'
	h2_zloss_r_g = ROOT.TH2F(sname, sname, nbins, 0., jet_R0, len(sjrs), 0., 1.)
	sname = 'h2_zloss_vs_r_quark'
	h2_zloss_r_q = ROOT.TH2F(sname, sname, nbins, 0., jet_R0, len(sjrs), 0., 1.)

	# loop

	if args.nev < 100:
		args.nev = 100
	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		# parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
		partons = pythiafjext.vectorize_select(pythia, [pythiafjext.kParton], 0, True)
		parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, False)
		# parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, False)
		jets = jet_selector(jet_def(parts))

		# for j in tqdm.tqdm(jets):
		for j in jets:
			j_type = match_dR(j, partons, jet_R0 / 2.)
			if j_type[0] is None:
				if args.nw:
					continue
				pwarning('Jet with no parton label')
				continue

			tw.fill_branch(	"j", j)
			for sjr in sjrs:
				rc_jets = fj.sorted_by_pt(sjdefs[sjr](j.constituents()))
				tw.fill_branch( "sjr{}".format(sjr), rc_jets[0])
				zloss = 1. - rc_jets[0].perp() / j.perp()
				tw.fill_branch( "sjr{}_zloss".format(sjr), zloss)
				tw.fill_branch("ppid", j_type[0])
				tw.fill_branch("pquark", j_type[1])
				tw.fill_branch("pglue", j_type[2])

				if j_type[1]:
					h_zloss_r_q[sjr].Fill(zloss)
					h2_zloss_r_q.Fill(sjr, zloss)
					prof_q.Fill(sjr, zloss)
					prof_q_log.Fill(ROOT.TMath.Log(1./sjr), zloss)
				if j_type[2]:
					h_zloss_r_g[sjr].Fill(zloss)
					h2_zloss_r_g.Fill(sjr, zloss)
					prof_g.Fill(sjr, zloss)
					prof_g_log.Fill(ROOT.TMath.Log(1./sjr), zloss)

			tw.fill_tree()

	pythia.stat()
	tw.write_and_close()


if __name__ == '__main__':
	main()
