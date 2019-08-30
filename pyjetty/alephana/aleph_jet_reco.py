#!/usr/bin/env python

import aleph
import fastjet as fj
import fjext
import fjcontrib
from tqdm import tqdm
import sys
import numpy as nd
import pandas as pd
import time
import aleph_utils
import joblib
from pyjetty.fjutils.fjevent import *
import os

class RT(object):
	df_columns = ['evid', 'ptleadjet', 'dphi']
	def __init__(self):
		self.df = []
	def process_event(self, parts, leadjet, id=-1):
		_parts = fjext.vectorize_px_py_pz_e(parts['px'].values, parts['py'].values, parts['pz'].values, parts['e'].values)
		# self.df = self.df.append(pd.DataFrame([ [ id, leadjet.pt(), p.delta_phi_to(leadjet) ] for p in _parts ], columns=self.df_columns))
		# dphis = [p.delta_phi_to(leadjet) for p in _parts]
		self.df.append(pd.DataFrame([ [ id, leadjet.pt(), p.delta_phi_to(leadjet) ] for p in _parts ], columns=self.df_columns))
	def get_pandas(self):
		return pd.concat(self.df, ignore_index=True)

def main():
	aleph_file="/Volumes/two/data/aleph/LEP1Data1992_recons_aftercut-001.aleph"
	if len(sys.argv) > 1:
		aleph_file = sys.argv[1]
	if not os.path.exists(aleph_file):
		print('[error] {} file does not exists'.format(aleph_file))
		return 0

	fj.ClusterSequence.print_banner()
	print()

	all_jets = []
	reader = aleph.Reader(aleph_file)
	nev = aleph_utils.get_n_events(aleph_file)
	pbar = tqdm(total=nev)
	partSelector = KineSelectorFactory(absetamax=2)
	jetSelector = KineSelectorFactory(absetamax=2, ptmin=10)

	fjev = FJEvent(	R=0.4, algorithm=fj.antikt_algorithm,
					jet_selector=jetSelector,
					particle_selector=partSelector)
	rtanalysis = RT()
	while reader.read_next_event():
		e = reader.get_event()
		vparts = e.get_particles_vdoubles()
		df = pd.DataFrame(vparts, columns=aleph.Particle.descr())
		parts = fjext.vectorize_px_py_pz_e(df['px'].values, df['py'].values, df['pz'].values, df['e'].values)
		fjev.run_jet_finder_csaa(particles=parts, evid=e.get_header().n())
		njets = len(fjev.inclusive_jets)
		if njets > 0:
			# print("njets = {} df:{}".format(njets, len(_jets_df)))
			_jets_df = fjev.jets_df
			# R_T analysis on charged tracks
			df_charged = df.loc[df['pwflag'] != 0 ]
			rtanalysis.process_event(df_charged, fjev.leading_pt_jet(), e.get_header().n())
		pbar.update()
		if pbar.n > 1000:
			break;
	pbar.close()

	all_jets = fjev.get_pandas()
	joblib.dump(all_jets, 'test_jets.joblib')
	joblib.dump(rtanalysis.get_pandas(), 'test_rt.joblib')

if __name__ == '__main__':
	main()