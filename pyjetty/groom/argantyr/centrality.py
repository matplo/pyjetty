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

# centrality bins with h_dndeta_lores_antyr_npart
# 90pc@6,
# 80pc@15,
# 70pc@31,
# 60pc@56,
# 50pc@90,
# 40pc@133,
# 30pc@186,
# 20pc@248,
# 10pc@325,
# centrality bins with h_dndeta_lores_antyr_nch
nch_cents='''
90pc@277,
80pc@555,
70pc@1158,
60pc@2272,
50pc@3945,
40pc@6239,
30pc@9293,
20pc@13245,
10pc@18467,
'''

def main():
	parser = argparse.ArgumentParser(description='test groomers', prog=os.path.basename(__file__))
	parser.add_argument('-o', '--output-filename', default="centrality_output.root", type=str)
	parser.add_argument('datalist', help='run through a file list', default='', type=str)
	parser.add_argument('--overwrite', help="overwrite output", default=False, action='store_true')
	parser.add_argument('--nev', help='number of events to run', default=0, type=int)
	parser.add_argument('--max-eta', help='max eta for particles', default=0.9)
	parser.add_argument('--thermal', help='enable thermal generator',action='store_true', default=False)
	parser.add_argument('--thermal-default', help='enable thermal generator',action='store_true', default=False)
	parser.add_argument('--particles', help='stream particles',action='store_true', default=False)
	parser.add_argument('--npart-cut', help='npart cut on centrality low,high hint:' + npart_cents, default='325,450', type=str)
	parser.add_argument('--nch-cut', help='nch cut on centrality low,high hint:' + nch_cents, default='18467,50000', type=str)

	args = parser.parse_args()

	try:
		npart_min = int(args.npart_cut.split(',')[0])
		npart_max = int(args.npart_cut.split(',')[1])
	except:
		perror('unable to parse npart centrality selection - two integer numbers with a coma in-between needed - specified:', args.npart_cut)
		return 1

	try:
		nch_min = int(args.nch_cut.split(',')[0])
		nch_max = int(args.nch_cut.split(',')[1])
	except:
		perror('unable to parse nch centrality selection - two integer numbers with a coma in-between needed - specified:', args.nch_cut)
		return 1

	outf = ROOT.TFile(args.output_filename, 'recreate')
	outf.cd()
	t = ROOT.TTree('t', 't')
	tw = RTreeWriter(tree=t)
	hpt_antyr = ROOT.TH1F('hpt_antyr', 'hpt_antyr', 100, 0, 100)
	hpt_antyr_c = ROOT.TH1F('hpt_antyr_c', 'hpt_antyr_c', 100, 0, 100)
	hpt_therm = ROOT.TH1F('hpt_therm', 'hpt_therm', 100, 0, 100)
	hpt_therm_c = ROOT.TH1F('hpt_therm_c', 'hpt_therm_c', 100, 0, 100)

	data = DataIO(name='Sim Pythia Detector level', file_list=args.datalist, random_file_order=False, tree_name='tree_Particle_gen')
	dndeta_selector = fj.SelectorAbsEtaMax(abs(args.max_eta)) & fj.SelectorPtMin(0.15)

	tg_default = None
	if args.thermal_default:
		tg_default = thg.ThermalGenerator()
		print(tg_default)

	tg_central = None
	if args.thermal:
		tg_central = thg.ThermalGenerator(beta=0.5, N_avg=3000, sigma_N=500)
		print(tg_central)

	delta_t = 0
	start_t = time.time()
	iev = 1
	while data.load_event(offset=0):
		iev = iev + 1
		if args.nev > 0:
			if iev > args.nev:
				iev = iev - 1
				break
		if iev % 1000 == 0:
			delta_t = time.time() - start_t
			pinfo('processing event', iev, ' - ev/sec =', iev/delta_t, 'elapsed =', delta_t)

		# find jets on detector level
		if len(data.particles) < 1:
			pwarning(iev, 'pp event skipped N parts', len(data.particles))
			continue

		# print(data.event)

		dndeta0_parts = dndeta_selector(data.particles)
		dndeta0 = len(dndeta0_parts)/(abs(args.max_eta*2.))
		[hpt_antyr.Fill(p.perp()) for p in dndeta0_parts]
		if args.particles:
			tw.fill_branches(dndeta=dndeta0, p=data.particles)
		else:
			tw.fill_branches(dndeta=dndeta0)
		tw.fill_branches_attribs(data.event, ['sigma', 'npart', 'nch', 'nchfwd', 'nchselect'], prefix='antyr_')

		if data.event.npart < npart_min or data.event.npart >= npart_max:
			tw.fill_branches(cent10npart=0)
		else:
			tw.fill_branches(cent10npart=1)
			[hpt_antyr_c.Fill(p.perp()) for p in dndeta0_parts]

		if data.event.nch < nch_min or data.event.nch >= nch_max:
			tw.fill_branches(cent10nch=0)
		else:
			tw.fill_branches(cent10nch=1)

		if tg_default:
			thg_particles = tg_default.load_event()
			dndetathg_default = dndeta_selector(thg_particles)
			if args.particles:
				tw.fill_branches(dndeta_thg_0=len(dndetathg_default)/(abs(args.max_eta*2.)), p_thg_0=thg_particles)
			else:
				tw.fill_branches(dndeta_thg_0=len(dndetathg_default)/(abs(args.max_eta*2.)))

		if tg_central:
			thg_parts_central = tg_central.load_event()
			dndetathg_central = dndeta_selector(thg_parts_central)
			[hpt_therm_c.Fill(p.perp()) for p in dndetathg_central]
			if args.particles:
				tw.fill_branches(dndeta_thg_c=len(dndetathg_central)/(abs(args.max_eta*2.)), p_thg_c=thg_parts_central)
			else:
				tw.fill_branches(dndeta_thg_c=len(dndetathg_central)/(abs(args.max_eta*2.)))

		tw.fill_tree()

	delta_t = time.time()-start_t
	pinfo('processed events', iev, ' - ev/sec =', iev/delta_t, 'elapsed =', delta_t)

	outf.Write()
	outf.Close()


if __name__ == '__main__':
	main()
