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


def main():
	aleph_file="/Volumes/two/data/aleph/LEP1Data1992_recons_aftercut-001.aleph"
	if len(sys.argv) > 1:
		aleph_file = sys.argv[1]
	if not os.path.exists(aleph_file):
		print('[error] {} file does not exists'.format(aleph_file))
		return 0

	# aleph.dump(aleph_file, 2, False);
	# aleph.dump(aleph_file, -1, True);

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()

	all_jets = []
	reader = aleph.Reader(aleph_file)
	nev = aleph_utils.get_n_events(aleph_file)
	pbar = tqdm(total=nev)
	partSelector = KineSelectorFactory(absetamax=2)
	jetSelector = KineSelectorFactory(absetamax=2, ptmin=10)

	all_jets = pd.DataFrame(columns=FJEvent.df_columns)
	while reader.read_next_event():
		e = reader.get_event()
		vparts = e.get_particles_vdoubles()
		df = pd.DataFrame(vparts, columns=aleph.Particle.descr())
		parts = fjext.vectorize_px_py_pz_e(df['px'].values, df['py'].values, df['pz'].values, df['e'].values)
		fjev = FJEvent(	particles=parts, id=e.get_header().n(), R=0.4, algorithm=fj.antikt_algorithm,
						jet_selector=jetSelector,
						particle_selector=partSelector)
		fjev.run_jet_finder_csaa()
		inclusive_jets = fjev.inclusive_jets
		njets = len(inclusive_jets)
		if njets > 0:
			# print("njets = {} df:{}".format(njets, len(_jets_df)))
			_jets_df = fjev.jets_df
			all_jets = all_jets.append(_jets_df, ignore_index=True)
		pbar.update()
		if pbar.n > 1000:
			break;
	pbar.close()

	joblib.dump(all_jets, 'test_jets.joblib')

if __name__ == '__main__':
	main()