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


def logbins(xmin, xmax, nbins):
	lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
	arr = array.array('f', lspace)
	return arr
	# return lspace


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
	E_protons = 3500.
	if args.ecm == 'high':
		E_protons = 6500.
	E_lead = E_protons / (208. / 82.)
	if args.ecm == 'high':
		E_lead = 5020.
	if args.ecm == 'low':
		E_lead = 2760.
	# sconfig_pythia = [	"Beams:eA = {}".format(E_protons), "Beams:eB = {}".format(E_lead), "Beams:frameType = 2" ];
	# sconfig_pythia = [	"Beams:eA = {}".format(E_lead), "Beams:eB = {}".format(E_lead), "Beams:frameType = 2" ];
	sconfig_pythia = [	"Beams:eCM = {}".format(E_lead)]

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
				sconfig_pythia.append("PDF:pSet = LHAPDF6:EPPS16nlo_CT14nlo_Pb208/{0}".format(args.epps16set))
				# sconfig_pythia.append("PDF:pSetB = LHAPDF6:EPPS16nlo_CT14nlo_Pb208/{0}".format(args.epps16set))
	return sconfig_pythia


def get_final_daughters(event, idx):
	retv = []
	daughters = event.daughterList(idx)
	for d in daughters:
		if event[d].isFinal():
			retv.append(d)
		else:
			_ds = get_final_daughters(event, d)
			if len(_ds):
				retv.extend(_ds)
	return retv

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

	pythia = None
	if args.write or args.generate:
		pythia = create_and_init_pythia(sconfig_pythia)
	if not pythia:
		print("[error] bad or no pythia create and init. stop here.")
		return

	if args.write:
		pyhepmcwriter = mp.Pythia8HepMCWrapper(args.write)
		for iEvent in tqdm(range(nevents), 'event'):
			if not pythia.next(): continue
			pyhepmcwriter.fillEvent(pythia)
		print("[i] done writing to {}".format(args.write))

	if args.generate:
		outfname = '{}_output.root'.format(args.generate).replace('.dat', '')
		print("output file: {}".format(outfname))
		foutput = r.TFile(outfname, "recreate")
		foutput.cd()
		lbins = logbins(1., 100, 10)
		hjetpt = r.TH1F('hjetpt', 'hjetpt', 10, lbins)
		hevent = r.TH1F('hevent', 'hevent', 10, 0, 10);
		hxsection = r.TNtuple('xsection', 'xsection', 'xsec:xsecErr')
		tnpart = r.TNtuple('tnpart', 'tnpart', 'pt:eta:phi:id5:id6')
		for iEvent in tqdm(range(nevents), 'event'):
			if not pythia.next(): continue
			hevent.Fill(0, pythia.info.weight())
			hxsection.Fill(pythia.info.sigmaGen(), pythia.info.sigmaErr())
			if not args.charm and not args.photon:
				for ip in range(pythia.event.size()):
					if abs(pythia.event[ip].id()) == 111:
						tnpart.Fill(pythia.event[ip].pT(), pythia.event[ip].eta(), pythia.event[ip].phi(),
						            pythia.event[5].id(), pythia.event[6].id())
						if pythia.event[ip].eta() > 3 and pythia.event[ip].eta() < 6:
							hjetpt.Fill(pythia.event[ip].pT())
			if args.charm:
				for ip in range(pythia.event.size()):
					if abs(pythia.event[ip].id()) == 421:
						tnpart.Fill(pythia.event[ip].pT(), pythia.event[ip].eta(), pythia.event[ip].phi(),
						            pythia.event[5].id(), pythia.event[6].id())
						if pythia.event[ip].eta() > 3 and pythia.event[ip].eta() < 6:
							hjetpt.Fill(pythia.event[ip].pT())
			if args.photon:
				for idx in range(5,7):
					if pythia.event[idx].id() == 22:
						finalds = get_final_daughters(pythia.event, idx)
						for fidx in finalds:
							if pythia.event[fidx].id() == 22:
								tnpart.Fill(pythia.event[fidx].pT(), pythia.event[fidx].eta(), pythia.event[fidx].phi(),
								            pythia.event[5].id(), pythia.event[6].id())
								if pythia.event[fidx].eta() > 3 and pythia.event[fidx].eta() < 6:
									hjetpt.Fill(pythia.event[fidx].pT())
		w = pythia.info.weightSum()
		xs = pythia.info.sigmaGen()
		tnpart.SetWeight(xs/w)
		hjetpt.Scale(xs/w)
		foutput.Write()
		foutput.Close()
	if pythia:
		pythia.stat()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='jets with different with EPPS16nlo_CT14nlo_Pb208',
	                                 prog=os.path.basename(__file__))
	parser.add_argument('-n', '--nevents', help='number of events', default=1000, type=int)
	# parser.add_argument('-g', '--generate', help='generate - no writing', action='store_true')
	parser.add_argument('-g', '--generate', help='generate - on-the-fly - write root file', type=str, default=None)
	parser.add_argument('-w', '--write', help='generate and write hepmc file', type=str, default=None)

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

	if args.write:
		args.write    = os.path.splitext(args.write)[0] + '.dat'
	if args.generate:
		args.generate = os.path.splitext(args.generate)[0] + '.dat'

	_base_outputname = args.write
	_base_outputname_gen = args.generate

	if args.allset and args.epps16set is None:
		main(args)
		for args.epps16set in range(0, 41):
			if _base_outputname:
				args.write = os.path.splitext(_base_outputname)[0] + '_EPPS16_set{0}.dat'.format(args.epps16set)
			if _base_outputname_gen:
				args.generate = os.path.splitext(_base_outputname_gen)[0] + '_EPPS16_set{0}.dat'.format(args.epps16set)
			main(args)
	else:
		if args.epps16set is not None:
			if args.epps16set >= 0:
				if args.write:
					args.write = os.path.splitext(args.write)[0] + '_EPPS16_set{0}.dat'.format(args.epps16set)
				if args.generate:
					args.generate = os.path.splitext(args.generate)[0] + '_EPPS16_set{0}.dat'.format(args.epps16set)
		main(args)
