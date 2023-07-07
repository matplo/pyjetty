#!/usr/bin/env python

from __future__ import print_function
from ast import parse
from flavor_tagger import FlavorTaggerUtil, FT_print_psj, FT_psj_to_str

import fastjet as fj
import fjcontrib
import fjext

import array
import tqdm
import argparse
import os
import sys

import pythia8
import pythiaext
import pythiafjext
# pythia8+hepmc3 writing
# import pythiahepmc3

import ecorrel

from heppy.pythiautils import configuration as pyconf
from pyjetty.mputils.mputils import logbins

import importlib
root_avail = importlib.util.find_spec("ROOT")
if root_avail is not None:
	import ROOT
	print('[i] ROOT at', ROOT.__file__)
else:
    print('[e] no ROOT? this script writes ROOT files...')
    sys.exit(1)
ROOT.PyConfig.IgnoreCommandLineOptions = True

# utility function to check how the script was called - exclusive process?
def parton_type_from_args(args):
	sqorg = 'any'
	if args.py_hardQCDgluons:
		sqorg = 'glue'
	if args.py_hardQCDquarks:
		sqorg = 'quark'
	return sqorg

def format_output_file(fname, nfile, args):
	if nfile < 0:
		foutname = '{}_{}{}'.format(os.path.splitext(fname)[0], parton_type_from_args(args), os.path.splitext(fname)[1])
	else:
		foutname = '{}_{}_{}{}'.format(os.path.splitext(fname)[0], parton_type_from_args(args), nfile, os.path.splitext(fname)[1])
	return foutname

class EECHistogramOutput(object):
	def __init__(self, filename):
		self.fout = ROOT.TFile(filename, 'recreate')
		self.fout.cd()
		self.h_eec = {}
		self.h_pt = {}

	def __del__(self):
		for h in self.h_eec:
			njets = self.h_pt[h].GetEntries()
			if njets < 1.:
				njets = 1.
			self.h_eec[h].Scale(1./njets, "width")
		self.fout.Write()
		self.fout.Close()
		print('[i] written', self.fout.GetName())

	def get_flavor_histograms(self, fid):
		try:
			self.h_eec[fid]
		except:
			self.fout.cd()
			self.h_eec[fid] = ROOT.TH1F(f'h_eec_{fid}', 'h_eec_fidx;R_{L};dN/dR_{L}'.replace('fidx', str(fid)), 20, logbins(0.001, 1., 20))
			self.h_pt[fid] = ROOT.TH1F(f'h_pt_{fid}', 'h_pt_fidx;p_{T}^{jet};dN/dp_{T}'.replace('fidx', str(fid)), 1000, 0, 1000)
		return self.h_eec[fid], self.h_pt[fid]

def main():
	parser = argparse.ArgumentParser(
		description='pythia8 in python', prog=os.path.basename(__file__))
	parser.add_argument('--jet-R', help='specify jet R for flavor tagging', default=0.4, type=float)
	parser.add_argument('--jet-ptmin', help='minimum pT cut on jets', default=20., type=float)
	parser.add_argument('--jet-ptmax', help='maximum pT cut on jets', default=1e4, type=float)
	parser.add_argument('--jet-absymax', help='maximum abs(rapidity) cut on jets', default=2., type=float)
	parser.add_argument('-o', '--output', help='output root file', default='eec_flavor_jets.root', type=str)
	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()

	# select particles for jet finding
	part_selector = fj.SelectorAbsRapMax(args.jet_absymax + args.jet_R + 0.1)
	part_selector_p = fj.SelectorAbsRapMax(args.jet_absymax + args.jet_R + 0.1)
	# selector for particles above some pT threshold
	pfc_selector1 = fj.SelectorPtMin(1.)

	jet_selector = fj.SelectorPtMin(args.jet_ptmin) * fj.SelectorPtMax(args.jet_ptmax) * fj.SelectorAbsRapMax(args.jet_absymax)

	# setup the event generator
	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 10:
		args.nev = 10

	# create the output
	output = EECHistogramOutput(args.output)

	# here we will count the jets not events
	event_number = 0
	pbar = tqdm.tqdm(range(args.nev))
	n_total_jets_accepted = 0
	while n_total_jets_accepted < args.nev:
		if not pythia.next():
			continue
		event_number = event_number + 1

		n_jets_accepted = 0 # in this event
  
		# get the final state particles as a vector of pseudojets
		parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kVisible], 0, True)
		# cut on rapidity
		parts_pythia_h = part_selector(parts_pythia_h)
		if len(parts_pythia_h) < 1 :
				continue

		# example how to use flavor tagger
		# note tagging goes with the highest pT parton within the jet as caught by the reco using partons as ghosts...
		fT = FlavorTaggerUtil(pythia)
		parts_pythia_p = pythiafjext.vectorize_select(pythia, [pythiafjext.kParton], 0, True)
		parts_pythia_p = part_selector_p(parts_pythia_p)
		if len(parts_pythia_p) < 1:
				continue

		# get the hadron level + parton level ghosted
		parts = fT.merged_particle_vector(parts_pythia_h)
  
		# run the jet finder
		jet_def = fj.JetDefinition(fj.antikt_algorithm, args.jet_R)
		jets = fj.sorted_by_pt(jet_def(parts))

		# select jets 
		jets_selected = jet_selector(fj.sorted_by_pt(jets))
		if len(jets_selected) < 1:
				continue
		# do something with the jets passing the kinematic cuts
		for j in jets_selected:
			# make sure the jet is not only made of parton ghosts (very rarely but may happen)
			if fT.is_jet_purely_parton_ghost(j) is True:
				continue
			# debug statement
			# FT_print_psj(j, '-         jet', show_constituents=False)
			# now get some flavor information
			j_flavor_psj = fT.flavor_tag_psj(j)
			# if we got a tag - should be almost always the case (unless out of ordinary jet with no partons within)
			if j_flavor_psj:
				# example usage commented out:
				# j_flavor_pythia_part = fT.flavor_tag_pythia_particle(j)
				j_flavor = fT.flavor_tag_id(j)
				# or j_flavor = j_flavor_pythia_part.id() 
				# FT_print_psj(j_flavor_psj, '  lead parton')
				dR = j.delta_R(j_flavor_psj)
				dphi = j.delta_phi_to(j_flavor_psj)
				# debug statement
				# print('  ', 'jet - leading parton deltaR={:.3f}\tdphi={:.3f}'.format(dR, dphi))
				# tighten the parton-jet correlation
				if dR > 0.2:
					continue
				# check if run in special process selected mode - reject or take only gluons
				if args.py_hardQCDgluons and j_flavor != 21:
					continue
				if args.py_hardQCDquarks and j_flavor == 21:
					continue
				if args.py_hardQCDcharm and j_flavor == 21:
					continue
				if args.py_hardQCDbeauty and j_flavor == 21:
					continue
				if args.py_hardQCDuds and j_flavor == 21:
					continue				
				# take absolute parton id - treat q and qbar the same...
				j_flavor = abs(j_flavor)
	   			# select only charged constituents with 1 GeV cut
				_vc1 = fj.vectorPJ()
				_ = [_vc1.push_back(c) for c in pfc_selector1(j.constituents())
            	                if pythiafjext.getPythia8Particle(c).isCharged()]
				# n-point correlator with charged particles pt > 1
				cb1 = ecorrel.CorrelatorBuilder(_vc1, j.perp(), 2) # up to 2-point correlator
				# fill output...
				n_jets_accepted += 1
				h_eec, h_pt = output.get_flavor_histograms(j_flavor)
				h_pt.Fill(j.perp())
				if cb1.correlator(2).rs().size() > 0:
					h_eec.FillN(cb1.correlator(2).rs().size(),
								array.array('d', cb1.correlator(2).rs()), 
								array.array('d', cb1.correlator(2).weights()))
		if n_jets_accepted < 1:
			continue
		n_total_jets_accepted += n_jets_accepted
		pbar.update(n_jets_accepted)
	pbar.close()

if __name__ == '__main__':
	main()
