#!/usr/bin/env python

from __future__ import print_function
from ast import arg, parse
from flavor_tagger import FlavorTaggerUtil, FT_print_psj, FT_psj_to_str

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os

import pythia8
import pythiaext
import pythiafjext
# pythia8+hepmc3 writing
# import pythiahepmc3

from heppy.pythiautils import configuration as pyconf

from read_hepmc2_jets_ml_v2 import JetFileOutput, UniqueRunNumber
import ROOT

def run(args):
        
	if args.py_seed <= 0:
		print('[e] has to speficy --py-seed ... bailing out')
		return False

	nev_nperfile = args.nperfile

	part_selector = fj.SelectorAbsRapMax(4.)
	part_selector_p = fj.SelectorAbsRapMax(4.)
	jet_selector = fj.SelectorPtMin(args.jet_ptmin) * fj.SelectorPtMax(args.jet_ptmax) * fj.SelectorAbsRapMax(2.5)

	z_selector = fj.SelectorPtMin(args.Z_ptmin) * fj.SelectorPtMax(args.Z_ptmax) * fj.SelectorAbsRapMax(2.5)

	jet_def_akt = fj.JetDefinition ( fj.antikt_algorithm, args.jet_R)

	fout = JetFileOutput(args=args)

	run_number = 1 # this number should change
	urn = UniqueRunNumber()
	run_number = urn.unique()
 
	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 10:
		args.nev = 10

	if args.Zjet is False:
		del args.mDT
		del args.mDTzcut
		# del args.ZjetR
	else:
		if args.py_cmnd is None:
			print('[e] you have to specify the command file for Zjet generation with --py-cmnd <file> - sorry, no shortcut ;-(')
			return

	if args.py_bias is False:
		del args.py_biaspow
		del args.py_biasref
   
	run_number = args.py_seed
	event_number = 0
	pbar = tqdm.tqdm(range(args.nev))
	# for i in pbar:
	n_total_jets_accepted = 0
	jets_accepted = []
	while n_total_jets_accepted < args.nev:
		if not pythia.next():
			continue
		event_number = event_number + 1
   
		# parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kVisible], 0, True)
		parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)
		parts_pythia_h = part_selector(parts_pythia_h)
		if len(parts_pythia_h) < 1 :
				continue

		# quark or gluon tagging
		if args.Zjet is False:
			fT = FlavorTaggerUtil(pythia)
			parts_pythia_p = pythiafjext.vectorize_select(pythia, [pythiafjext.kParton], 0, True)
			parts_pythia_p = part_selector_p(parts_pythia_p)
			if len(parts_pythia_p) < 1:
					continue
			parts = fT.merged_particle_vector(parts_pythia_h)
			# jet_def = fj.JetDefinition(fj.antikt_algorithm, args.jet_R)
			jets = fj.sorted_by_pt(jet_def_akt(parts))

			jets_selected = jet_selector(fj.sorted_by_pt(jets))
			if len(jets_selected) < 1:
				continue
	
			n_jets_accepted = 0
			for j in jets_selected:
				if fT.is_jet_purely_parton_ghost(j) is False:
					# FT_print_psj(j, '-         jet', show_constituents=False)
					j_flavor_psj = fT.flavor_tag_psj(j)
					if j_flavor_psj:
						j_flavor = fT.flavor_tag_id(j)
						dR = j.delta_R(j_flavor_psj)
						dphi = j.delta_phi_to(j_flavor_psj)
						if dR > 0.2:
							continue
						if args.py_hardQCDgluons and j_flavor != 21:
							continue
						if args.py_hardQCDquarks and j_flavor == 21:
							continue
						jets_accepted.append([j, j_flavor, j_flavor_psj])
						n_jets_accepted += 1
					else:
						pass

			if n_jets_accepted < 1:
				continue

			if len(jets_accepted) >= args.nperfile:
				fout.write_file_partons(jets_accepted, pythia)
				jets_accepted.clear()

			n_total_jets_accepted += n_jets_accepted
			pbar.update(n_jets_accepted)
			continue

		# Zjet generation
		if args.Zjet is True:
			n_jets_accepted = 0
			parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)
			# take all particles and find the Z
			parts_pythia_all = pythiafjext.vectorize_select(pythia, [pythiafjext.kAny], 0, True)
			zs = [p for p in parts_pythia_all if pythiafjext.getPythia8Particle(p).id() == 23]
			# if len(zs) > 1:
			# 	print('[debug] more than one Z in the stack', len(zs), zs[0].perp(), zs[0].phi(), zs[1].perp(), zs[1].phi())
			if len(zs) < 1:
				print('[w] no Z in the event?')
				continue
			# print('[i] {} Z(s) found:'.format(event_number), len(zs))
			# now we got the Z
			psjZ = zs[0]
			# find anti-kT jets
			akt_jets = fj.sorted_by_pt(jet_def_akt(parts_pythia_h))
			# Z should be close to either of two hardest jets (in R space).
			if (len(akt_jets) < 2):
				print('[i] should have two jets - found just one')
				continue
			Zjet0R = psjZ.delta_R(akt_jets[0])
			Zjet1R = psjZ.delta_R(akt_jets[1])
			ZjetDR = None
			Zjet = None
			if Zjet0R < Zjet1R:
				Zjet = akt_jets[0]
				ZjetDR = Zjet0R
			else:
				Zjet = akt_jets[1]
				ZjetDR = Zjet1R
			if Zjet is None:
				print('[w] unable to match Z and and any jet within R < 1.', Zjet0R, Zjet1R, ZjetDR)
				continue
			# stream the anti-kT jets
			if len(jet_selector([Zjet])) < 1:
				continue
			jets_accepted.append([Zjet, psjZ])
			n_jets_accepted = 1

			if n_jets_accepted < 1:
				continue

			if len(jets_accepted) >= args.nperfile:
				fout.write_file_Zjets(jets_accepted, pythia)
				jets_accepted.clear()

			n_total_jets_accepted += n_jets_accepted
			pbar.update(n_jets_accepted)
			continue

		break

	if len(jets_accepted) > 0:
		if args.Zjet is False:
			fout.write_file_partons(jets_accepted, pythia=pythia)
			jets_accepted.clear()
		else:
			fout.write_file_Zjets(jets_accepted, pythia=pythia)
			jets_accepted.clear()

	pbar.close()
	pythia.stat()
	pythia.settings.writeFile(fout.cmnd_file_name)

def main():
	parser = argparse.ArgumentParser(description='pythia8 in python', prog=os.path.basename(__file__))
	parser.add_argument('--hepmc', help='write the hepmc files', default=False, action='store_true')
	parser.add_argument('--jet-R', help='specify jet R for flavor tagging', default=0.8, type=float)
	parser.add_argument('-c', '--nperfile', help='nevents per output file', default=1000, type=int)
	parser.add_argument('--jet-ptmin', help='minimum pT cut on jets', default=500., type=float)
	parser.add_argument('--jet-ptmax', help='maximum pT cut on jets', default=550., type=float)
	parser.add_argument('--Zjet', help='force Zjet generation - will read pythia_gen_qorg_Zjet_master.cmnd from current directory', default=False, action='store_true')
	parser.add_argument('--mDT', help='mDT clean up with z_cut=0.04 (default) on the Zjets', default=True, action='store_true')
	parser.add_argument('--mDTzcut', help='mDT clean up with z_cut=0.04 (default) on the Zjets', default=0.04, type=float)

	parser.add_argument('--Z-ptmin', help='minimum pT cut on Z', default=50., type=float)
	parser.add_argument('--Z-ptmax', help='maximum pT cut on Z', default=550., type=float)

	parser.add_argument('-i', '--input', help='hepmc file input to analyze', default=None, type=str)
	# parser.add_argument('--nev', help='number of events to process', default=1000, type=int)
	parser.add_argument('-l', '--list', help='treat input file as a file containing list of files to process', action='store_true', default=False)
	parser.add_argument('-g', '--debug', help='print some extras', default=False, action='store_true')
	parser.add_argument('-o', '--output', help='overwrite default naming of the output', default='', type=str)
	parser.add_argument('-t', '--threads', help='mutlithread', default=1, type=int)
	parser.add_argument('--slurm', help='get run number from slurm batch job id', action='store_true', default=False)

	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()

	run(args)

if __name__ == "__main__":
	main()