#!/usr/bin/env python

import os
import argparse
import uproot
import ROOT as r

import fastjet as fj
import pythia8
from recursivetools import pyrecursivetools as rt
from lundplane import pylundplane as lund
from pythiafjtools import pypythiafjtools as pyfj
from mptools import pymptools as mp

import pandas as pd
import joblib

def fj_parts_from_tracks(tracks):
	fjparts = []
	_pts  = tracks['ParticlePt']
	_etas = tracks['ParticleEta']
	_phis = tracks['ParticlePhi']
	for index, row in tracks.iterrows():
		lv = r.Math.PtEtaPhiMVector(row['ParticlePt'], row['ParticleEta'], row['ParticlePhi'], 0.0)
		psj = fj.PseudoJet(lv.Px(), lv.Py(), lv.Pz(), lv.E())
		psj.set_user_index(index)
		fjparts.append(psj)
	return fjparts

# MP note: read more about uproot unpacking at https://github.com/scikit-hep/uproot#filling-pandas-dataframes
# MP note: there may be more efficient way to do this...

def main(args):
	fname = args.fname
	file = uproot.open(fname)
	all_ttrees = dict(file.allitems(filterclass=lambda cls: issubclass(cls, uproot.tree.TTreeMethods)))
	tracks = all_ttrees[b'PWGHF_TreeCreator/tree_Particle;1']
	pds_trks = tracks.pandas.df() # entrystop=10)
	events = all_ttrees[b'PWGHF_TreeCreator/tree_event_char;1']
	pds_evs = events.pandas.df()

	# print the banner first
	fj.ClusterSequence.print_banner()
	jet_R0 = args.jetR
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(0.0) & fj.SelectorPtMax(1000.0) & fj.SelectorAbsEtaMax(1)
	print(jet_def)
	print()

	e_jets = pd.DataFrame(columns=['evid', 'pt', 'eta', 'phi'])

	for i, e in pds_evs.iterrows():
		iev_id = int(e['ev_id'])
		_ts = pds_trks.loc[pds_trks['ev_id'] == iev_id]
		_tpsj = fj_parts_from_tracks(_ts)
		_jets = jet_selector(jet_def(_tpsj))
		# _jets_a = [[iev_id, j.perp(), j.eta(), j.phi()] for j in _jets]
		# _jets_a = pd.DataFrame(np.array([[iev_id, j.perp(), j.eta(), j.phi()] for j in _jets]), columns=['evid', 'pt', 'eta', 'phi'])
		_jets_a = pd.DataFrame([[iev_id, j.perp(), j.eta(), j.phi()] for j in _jets], columns=['evid', 'pt', 'eta', 'phi'])
		# , columns=['evid, pt, eta, phi']
		e_jets = e_jets.append(_jets_a, ignore_index=True)
		print('event', i, 'number of parts', len(_tpsj), 'number of jets', len(_jets))
		print(_jets_a.describe())

	print(e_jets.describe())
	joblib.dump(e_jets, args.output)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='jet reco on alice data', prog=os.path.basename(__file__))
	parser.add_argument('-n', '--nevents', help='number of events', default=1000, type=int)
	parser.add_argument('-f', '--fname', help='input file name', type=str, default=None, required=True)
	parser.add_argument('-R', '--jetR', help='jet radius', default=0.4, type=float)
	parser.add_argument('-o', '--output', help='output file name', default='{}_output.joblib'.format(os.path.basename(__file__)), type=str)
	args = parser.parse_args()	
	main(args)
	#fname = '/Users/ploskon/data/HFtree_trains/13-06-2019/488_20190613-0256/unmerged/child_1/0001/AnalysisResults.root'
