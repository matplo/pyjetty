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

def create_and_init_pythia(config_strings=[]):
	pythia = pythia8.Pythia()
	for s in config_strings:
		pythia.readString(s)
	for extra_s in ["Next:numberShowEvent = 0", "Next:numberShowInfo = 0", "Next:numberShowProcess = 0", "Next:numberCount = 0"]:
		pythia.readString(extra_s)
	if pythia.init():
		return pythia
	return None


def get_pythia_config(args, procsel=[]):
	# sconfig_pythia = [	"Beams:eCM = {}".format(args.ecm)]
	E_protons = 4000.
	if args.ecm == 'high':
		E_protons = 6500.
	E_lead = E_protons / (208. / 82.)
	sconfig_pythia = [	"Beams:eA = {}".format(E_protons), "Beams:eB = {}".format(E_lead), "Beams:frameType = 2" ];
	for s in procsel:
		sconfig_pythia.append(s)
	if args.pthatmin < 0:
		_extra = [	"PhaseSpace:bias2Selection=on",
					"PhaseSpace:bias2SelectionPow={}".format(args.biaspow),
					"PhaseSpace:bias2SelectionRef={}".format(args.biasref)]
		for _es in  _extra:
			sconfig_pythia.append(_es)
	else:
		sconfig_pythia.append("PhaseSpace:pTHatMin = {}".format(args.pthatmin))
	if args.noue:
		sconfig_pythia.append("PartonLevel:ISR = off")
		sconfig_pythia.append("PartonLevel:MPI = off")
	if args.epps16set is not None:
		if args.epps16set >= 0:
			if abs(args.epps16set) > 40:
				print("[e] bad epps16set {}".format(args.epps16set))
				return None
			else:
				# sconfig_pythia.append("PDF:pSet = LHAPDF6:EPPS16nlo_CT14nlo_Pb208/{0}".format(args.epps16set))
				sconfig_pythia.append("PDF:pSetB = LHAPDF6:EPPS16nlo_CT14nlo_Pb208/{0}".format(args.epps16set))
	return sconfig_pythia


def main(args):
	nevents = args.nevents
	procsel = []
	if args.photon:
		procsel.append("PromptPhoton:all = on")
	if args.charm:
		procsel.append("HardQCD:hardccbar = on")
	if len(procsel) < 1:
		procsel.append("HardQCD:all = on")
	sconfig_pythia = get_pythia_config(args, procsel)

	if args.write:
		pythia = create_and_init_pythia(sconfig_pythia)
		if not pythia:
			return
		pyhepmcwriter = mp.Pythia8HepMCWrapper(args.write)
		for iEvent in tqdm(range(nevents), 'event'):
			if not pythia.next(): continue
			pyhepmcwriter.fillEvent(pythia)
		print("[i] done writing to {}".format(args.write))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='jets with different with EPPS16nlo_CT14nlo_Pb208',
	                                 prog=os.path.basename(__file__))
	parser.add_argument('-n', '--nevents', help='number of events', default=1000, type=int)
	parser.add_argument('-g', '--generate', help='generate - no writing', action='store_true')
	parser.add_argument('-w', '--write', help='generate and write hepmc file', type=str, default='pPb_output.dat')

	#	sconfig_pythia = get_pythia_config(args.ecm, args.pthatmin, args.biaspow, args.biasref, args.epps16set)
	parser.add_argument('--ecm', help='low or high sqrt(s) GeV', default='low', type=str)
	parser.add_argument('--pthatmin', help='minimum hat{pT}', default=-1, type=float)
	parser.add_argument('--biaspow', help='power of the bias (hard)', default=4, type=float)
	parser.add_argument('--biasref', help='reference pT for the bias', default='50', type=float)
	parser.add_argument('--epps16set', help='set number of EPPS16nlo_CT14nlo_Pb208', default=None, type=int)
	parser.add_argument('--allset', help='run for all sets in EPPS16nlo_CT14nlo_Pb208', default=False, action='store_true')
	parser.add_argument('--noue', help="no underlying event - equivalend to no ISR and MPIs set to off", default=False, action='store_true')
	parser.add_argument('--charm', help="enable hardccbar", default=False, action='store_true')
	parser.add_argument('--photon', help="enable prompt photon production",  default=False, action='store_true')
	args = parser.parse_args()

	args.write = os.path.splitext(args.write)[0] + '.dat'

	if args.allset and args.epps16set is None:
		main(args)
		_base_outputname = args.write
		for args.epps16set in range(0, 41):
			if _base_outputname:
				args.write = os.path.splitext(_base_outputname)[0] + '_EPPS16_set{0}.dat'.format(args.epps16set)
			main(args)
	else:
		if args.epps16set is not None:
			if args.epps16set >= 0:
				args.write = os.path.splitext(args.write)[0] + '_EPPS16_set{0}.dat'.format(args.epps16set)
		main(args)
