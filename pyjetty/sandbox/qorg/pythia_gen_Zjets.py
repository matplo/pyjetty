#!/usr/bin/env python

from __future__ import print_function

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

import importlib
root_avail = importlib.util.find_spec("ROOT")
if root_avail is not None:
	import ROOT


def parton_type_from_args(args):
	sqorg = 'any'
	if args.py_hardQCDgluons:
		sqorg = 'glue'
	if args.py_hardQCDquarks:
		sqorg = 'quark'
	if args.Zjet:
		sqorg = 'Zjet'
	return sqorg


def format_output_file(fname, nfile, args):
	if nfile < 0:
		foutname = '{}_{}{}'.format(os.path.splitext(
			fname)[0], parton_type_from_args(args), os.path.splitext(fname)[1])
	else:
		foutname = '{}_{}_{}{}'.format(os.path.splitext(
			fname)[0], parton_type_from_args(args), nfile, os.path.splitext(fname)[1])
	return foutname


def main():
	parser = argparse.ArgumentParser(
		description='pythia8 in python', prog=os.path.basename(__file__))
	parser.add_argument('--hepmc', help='write the hepmc files',
	                    default=False, action='store_true')
	parser.add_argument('--ml', help='write the ml format root files',
	                    default=False, action='store_true')
	parser.add_argument(
		'-c', '--nperfile', help='nevents per output file', default=1000, type=int)
	parser.add_argument('--Zjet', help='force Zjet generation - will read pythia_gen_qorg_Zjet_master.cmnd from current directory', default=False, action='store_true')
	parser.add_argument('--mDT', help='mDT clean up with z_cut=0.04 (default) on the Zjets', default=False, action='store_true')
	parser.add_argument('--mDTzcut', help='mDT clean up with z_cut=0.04 (default) on the Zjets', default=0.04, type=float)
	parser.add_argument('--ZjetR', help='specify the radius of the anti-kT jet for Z-match - default is R=1.0', default=1.0, type=float)
	parser.add_argument('--jet-ptmin', help='minimum pT cut on jets', default=5., type=float)
	parser.add_argument('--jet-ptmax', help='maximum pT cut on jets', default=1e4, type=float)
	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()

	if args.hepmc:
		hepmc2_output_base = "pythia_gen_qorgorZjet.hepmc"
		hepmc2_output_name = None
		hepmc2_writer = None
		hepmc2_fileno = 0

	if args.ml:
		if root_avail is None:
			print('[error] unable to load ROOT so --ml is defunct')
			args.ml = False
		ml_root_output_base = "pythia_gen_qorgorZjet.root"
		ml_root_output_name = None
		ml_root_file = None
		ml_root_ntuple_parts = None
		ml_root_ntuple_ev = None
		ml_root_fileno = 0

	if args.hepmc is False and args.ml is False:
		print('[error] one of the --hepmc or --ml is required')
		print('		   --help for more options')
		return

	if args.Zjet:
		if args.py_cmnd is None:
			print('[e] you have to specify the command file for Zjet generation with --py-cmnd <file> - sorry ;-(')
			return

	jet_def_akt = fj.JetDefinition ( fj.antikt_algorithm, args.ZjetR);
	jet_def_ca = fj.JetDefinition ( fj.cambridge_algorithm, fj.JetDefinition.max_allowable_R)
	jet_selector = fj.SelectorPtMin(args.jet_ptmin) * fj.SelectorPtMax(args.jet_ptmax) * fj.SelectorAbsRapMax(1.7)

	mDT = None
	if args.mDT:
		mDT = fjcontrib.ModifiedMassDropTagger(args.mDTzcut)
		if args.hepmc:
			hepmc2_output_base = hepmc2_output_base.replace('.hepmc', '_mDT{}.hepmc'.format(args.mDTzcut))
		if args.ml:
			ml_root_output_base = ml_root_output_base.replace('.root', '_mDT{}.root'.format(args.mDTzcut))


	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 10:
		args.nev = 10
	run_number = args.py_seed
	event_number = 0
	pbar = tqdm.tqdm(range(args.nev))
	# for i in pbar:
	n_total_jets_accepted = 0
	while n_total_jets_accepted < args.nev:
		if not pythia.next():
			continue
		event_number = event_number + 1
		if n_total_jets_accepted % args.nperfile == 0:
			if args.hepmc:
				if hepmc2_writer:
					del hepmc2_writer
				hepmc2_output_name = format_output_file(
					hepmc2_output_base, hepmc2_fileno, args)
				# print('[i] new file', hepmc2_output_name)
				hepmc2_writer = pythiaext.Pythia8HepMC2Wrapper(hepmc2_output_name)
				hepmc2_fileno = hepmc2_fileno + 1
			if args.ml:
				if ml_root_file:
					if ml_root_ntuple_parts:
						ml_root_ntuple_parts.Write()
						ml_root_ntuple_ev.Write()
					ml_root_file.Write()
					ml_root_file.Purge()
					ml_root_file.Close()
					ml_root_file = None
					ml_root_ntuple_ev = None
					ml_root_ntuple_zjet = None
					ml_root_ntuple_parts = None
					ml_root_ntuple_parts_mDT = None
				ml_root_output_name = format_output_file(ml_root_output_base, ml_root_fileno, args)
				ml_root_file = ROOT.TFile(ml_root_output_name, 'recreate')
				ml_root_ntuple_parts = ROOT.TNtuple('tree_Particle_gen', 'particles from PYTHIA8 - {} jets'.format(parton_type_from_args(args)),
                                        'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID')
				ml_root_ntuple_parts_mDT = ROOT.TNtuple('tree_Particle_gen_mDT{}'.format(args.mDTzcut), 'particles from PYTHIA8 - {} jets'.format(parton_type_from_args(args)),
                                        'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID')
				ml_root_ntuple_ev = ROOT.TNtuple('tree_Event_gen', 'event info from PYTHIA8 - {} jets'.format(parton_type_from_args(args)),
                                     'run_number:ev_id:xsec:code')
				ml_root_ntuple_zjet = ROOT.TNtuple('tree_Zjet_gen', 'Zjet kinematics','run_number:ev_id:ZjetPt:ZjetEta:ZjetPhi:ParticlePt:ParticleEta:ParticlePhi:ParticlePID:dR')
				ml_root_fileno = ml_root_fileno + 1

		if args.hepmc:
			hepmc2_writer.fillEvent(pythia)

		if args.ml:
			parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)
			ml_root_ntuple_ev.Fill(run_number, event_number, pythia.info.sigmaGen(), pythia.info.code())
			if args.Zjet:
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
				ml_root_ntuple_zjet.Fill(run_number, event_number, 
										Zjet.pt(), Zjet.eta(), Zjet.phi(), 
										psjZ.pt(), psjZ.eta(), psjZ.phi(), 
										pythiafjext.getPythia8Particle(psjZ).id(),
										ZjetDR)
				for p in Zjet.constituents():
					ml_root_ntuple_parts.Fill(run_number, event_number, p.pt(), p.eta(), p.phi(), pythiafjext.getPythia8Particle(p).id())
				if mDT:
					# cleanup with mDT and stream the resulting jet
					ca_jets = fj.sorted_by_pt(jet_def_ca(Zjet.constituents()))
					if (len(ca_jets) > 1):
						print('[w] CA reclustering returns more than 1 subjet (?)')
					else:
						# Use modified mass-drop tagger to clean up jet.
						taggedJet = mDT.result(ca_jets[0])
						if taggedJet.has_constituents():
							for p in taggedJet.constituents():
								ml_root_ntuple_parts_mDT.Fill(run_number, event_number, p.pt(), p.eta(), p.phi(), pythiafjext.getPythia8Particle(p).id())
					pass
				n_total_jets_accepted += 1
			else:
				# kept for ref for the q or glue
				for p in parts_pythia_h:
					ml_root_ntuple_parts.Fill(run_number, event_number, p.pt(
					), p.eta(), p.phi(), pythiafjext.getPythia8Particle(p).id())
			pbar.update()

	pythia.stat()
	pythia.settings.writeFile(
		format_output_file('pythia_gen_qorgorZjet.cmnd', -1, args))

	if args.ml:
		if ml_root_file:
			if ml_root_ntuple_parts:
				ml_root_ntuple_parts.Write()
				ml_root_ntuple_ev.Write()
			ml_root_file.Write()
			ml_root_file.Purge()
			ml_root_file.Close()


if __name__ == '__main__':
	main()
