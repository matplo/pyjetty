#!/usr/bin/env python3
import sys
import fastjet as fj
import pythia8
from recursivetools import pyrecursivetools as rt
from lundplane import pylundplane as lund
from pythiafjtools import pypythiafjtools as pyfj
from mptools import pymptools as mp
from tqdm import tqdm
import argparse
import os
import ROOT as r
import numpy as np
import array
import pyhepmc_ng
import joblib

def logbins(xmin, xmax, nbins):
	lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
	arr = array.array('f', lspace)
	return arr
	# return lspace

def main_jets(args):
	outfname = '{}_output.root'.format(args.read).replace('.dat', '')
	print("output file: {}".format(outfname))
	foutput = r.TFile(outfname, "recreate")
	foutput.cd()
	lbins = logbins(1., 100, 10)
	hjetpt = r.TH1F('hjetpt', 'hjetpt', 10, lbins)

	input = pyhepmc_ng.ReaderAsciiHepMC2(args.read)
	if input.failed():
		print ("[error] unable to read from {}".format(args.read))
		return

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	jet_selector = fj.SelectorPtMin(0.0) & fj.SelectorPtMax(200.0) & fj.SelectorAbsEtaMax(3)

	all_jets = []
	event = pyhepmc_ng.GenEvent()
	pbar = tqdm(range(args.nevents))
	while not input.failed():
		e = input.read_event(event)
		if input.failed():
			break
		fjparts = []
		for i,p in enumerate(event.particles):
			if p.status == 1 and not p.end_vertex:
				psj = fj.PseudoJet(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e)
				psj.set_user_index(i)
				fjparts.append(psj)
		jets = jet_selector(jet_def(fjparts))
		all_jets.append([ [j.pt(), j.eta()] for j in jets])
		pbar.update()
		for j in jets:
			hjetpt.Fill(j.perp())
		if pbar.n >= args.nevents:
			break
	foutput.Write()
	foutput.Close()
	joblib.dump(all_jets, outfname.replace(".root", ".joblib"))

def main_parts(args):
	outfname = '{}_output.root'.format(args.read).replace('.dat', '')
	print("output file: {}".format(outfname))
	foutput = r.TFile(outfname, "recreate")
	foutput.cd()
	lbins = logbins(1., 100, 10)
	hpt = r.TH1F('hpt', 'hpt', 10, lbins)

	input = pyhepmc_ng.ReaderAsciiHepMC2(args.read)
	if input.failed():
		print ("[error] unable to read from {}".format(args.read))
		return

	# print the banner first
	all_parts = []
	event = pyhepmc_ng.GenEvent()
	pbar = tqdm(range(args.nevents))
	while not input.failed():
		e = input.read_event(event)
		if input.failed():
			break
		fjparts = []
		for i,p in enumerate(event.particles):
			if p.status == 1 and not p.end_vertex:
				psj = fj.PseudoJet(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e)
				psj.set_user_index(i)
				fjparts.append(psj)
		jets = jet_selector(jet_def(fjparts))
		all_jets.append([ [j.pt(), j.eta()] for j in jets])
		pbar.update()
		for j in jets:
			hjetpt.Fill(j.perp())
		if pbar.n >= args.nevents:
			break
	foutput.Write()
	foutput.Close()
	joblib.dump(all_jets, outfname.replace(".root", ".joblib"))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='jets with different with EPPS16nlo_CT14nlo_Pb208',
	                                 prog=os.path.basename(__file__))
	parser.add_argument('-n', '--nevents', help='number of events', default=1000, type=int)
	parser.add_argument('-g', '--generate', help='generate - no writing', action='store_true')
	parser.add_argument('-w', '--write', help='generate and write hepmc file', type=str, default='')
	parser.add_argument('-r', '--read', help='read hepmc file', type=str, default=None, required=True)
	#	sconfig_pythia = get_pythia_config(args.ecm, args.pthatmin, args.biaspow, args.biasref, args.epps16set)
	parser.add_argument('--ecm', help='sqrt(s) GeV', default=5000., type=float)
	parser.add_argument('--pthatmin', help='minimum hat{pT}', default=-1, type=float)
	parser.add_argument('--biaspow', help='power of the bias (hard)', default=4, type=float)
	parser.add_argument('--biasref', help='reference pT for the bias', default='50', type=float)
	parser.add_argument('--epps16set', help='set number of EPPS16nlo_CT14nlo_Pb208', default=None, type=int)
	parser.add_argument('--allset', help='run for all sets in EPPS16nlo_CT14nlo_Pb208', default=False, action='store_true')
	parser.add_argument('--noue', help="no underlying event - equivalend to no ISR and MPIs set to off", default=False, action='store_true')
	parser.add_argument('--fj', help="run fastjet", default=False, action='store_true')
	args = parser.parse_args()
	if args.allset and args.epps16set is None:
		main(args)
		_base_outputname = args.write
		for args.epps16set in range(0, 41):
			if _base_outputname:
				args.write = _base_outputname.rstrip('.dat') + '.dat'
				args.write = args.write.replace('.dat', '_EPPS16_set{0}.dat'.format(args.epps16set))
			main(args)
	else:
		main(args)
