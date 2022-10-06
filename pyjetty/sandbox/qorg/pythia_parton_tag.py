#!/usr/bin/env python

from __future__ import print_function
from ast import parse
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
	return sqorg


def format_output_file(fname, nfile, args):
	if nfile < 0:
		foutname = '{}_{}{}'.format(os.path.splitext(fname)[0], parton_type_from_args(args), os.path.splitext(fname)[1])
	else:
		foutname = '{}_{}_{}{}'.format(os.path.splitext(fname)[0], parton_type_from_args(args), nfile, os.path.splitext(fname)[1])
	return foutname


def main():
	parser = argparse.ArgumentParser(
		description='pythia8 in python', prog=os.path.basename(__file__))
	parser.add_argument('--hepmc', help='write the hepmc files',
	                    default=False, action='store_true')
	parser.add_argument('--ml', help='write the ml format root files',
	                    default=False, action='store_true')
	parser.add_argument('--jet-R', help='specify jet R for flavor tagging', default=0.4, type=float)
	parser.add_argument(
		'-c', '--nperfile', help='nevents per output file', default=1000, type=int)
	parser.add_argument('--jet-ptmin', help='minimum pT cut on jets', default=5., type=float)
	parser.add_argument('--jet-ptmax', help='maximum pT cut on jets', default=1e4, type=float)
	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()

	if args.hepmc:
		hepmc2_output_base = "pythia_parton_tag.hepmc"
		hepmc2_output_name = None
		hepmc2_writer = None
		hepmc2_fileno = 0

	if args.ml:
		if root_avail is None:
			print('[error] unable to load ROOT so --ml is defunct')
			args.ml = False
		ml_root_output_base = "pythia_parton_tag.root"
		ml_root_output_name = None
		ml_root_file = None
		ml_root_ntuple_parts = None
		ml_root_ntuple_ev = None
		ml_root_fileno = 0

	if args.hepmc is False and args.ml is False:
		print('[error] one of the --hepmc or --ml is required')
		print('		   --help for more options')
		return

	part_selector = fj.SelectorAbsRapMax(4.)
	part_selector_p = fj.SelectorAbsRapMax(4.)
	jet_selector = fj.SelectorPtMin(args.jet_ptmin) * fj.SelectorPtMax(args.jet_ptmax) * fj.SelectorAbsRapMax(1.7)
 
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
		# if pbar.n == 0 or pbar.n % args.nperfile == 0:
		if n_total_jets_accepted % args.nperfile == 0:
			if args.hepmc:
				if hepmc2_writer:
					del hepmc2_writer
				hepmc2_output_name = format_output_file(hepmc2_output_base, hepmc2_fileno, args)
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
					ml_root_ntuple_parts = None
				ml_root_output_name = format_output_file(ml_root_output_base, ml_root_fileno, args)
				ml_root_file = ROOT.TFile(ml_root_output_name, 'recreate')
				ml_root_ntuple_parts = ROOT.TNtuple('tree_Particle_gen', 'particles from PYTHIA8 - {} jets'.format(parton_type_from_args(args)),
                                        'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID')
				ml_root_ntuple_ev = ROOT.TNtuple('tree_Event_gen', 'event info from PYTHIA8 - {} jets'.format(parton_type_from_args(args)),
                                     'run_number:ev_id:xsec:code:partonID')
				ml_root_fileno = ml_root_fileno + 1

		if args.hepmc:
			hepmc2_writer.fillEvent(pythia)

		n_jets_accepted = 0
		if args.ml:
			parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kVisible], 0, True)
			parts_pythia_h = part_selector(parts_pythia_h)
			if len(parts_pythia_h) < 1 :
					continue

			# stream all particles in the event:
			# ml_root_ntuple_ev.Fill(run_number, event_number,
			#                        pythia.info.sigmaGen(), pythia.info.code())
			# for p in parts_pythia_h:
			# 	ml_root_ntuple_parts.Fill(run_number, event_number, p.pt(), p.eta(), p.phi(), pythiafjext.getPythia8Particle(p).id())

			# example how to use flavor tagging
			# note tagging goes with the highest pT parton within the jet as caught by the reco using partons as ghosts...
			fT = FlavorTaggerUtil(pythia)
			parts_pythia_p = pythiafjext.vectorize_select(pythia, [pythiafjext.kParton], 0, True)
			parts_pythia_p = part_selector_p(parts_pythia_p)
			if len(parts_pythia_p) < 1:
					continue
			parts = fT.merged_particle_vector(parts_pythia_h)
			jet_def = fj.JetDefinition(fj.antikt_algorithm, args.jet_R)
			jets = fj.sorted_by_pt(jet_def(parts))
			# print_pt_cut = pythia.info.pTHat() * 0.8
			# print('* event', event_number, 'this event pThat=', pythia.info.pTHat(), 'showing jets pT >', print_pt_cut)
			# print('* event', event_number, 'this event pThat=', pythia.info.pTHat(), 'showing jets {} < pT < {}'.format(args.jet_ptmin, args.jet_ptmax))
			jets_selected = jet_selector(fj.sorted_by_pt(jets))
			if len(jets_selected) < 1:
					continue
			for j in jets_selected:
				if fT.is_jet_purely_parton_ghost(j) is False:
					# FT_print_psj(j, '-         jet', show_constituents=False)
					j_flavor_psj = fT.flavor_tag_psj(j)
					if j_flavor_psj:
						# example usage commented out:
						# j_flavor_pythia_part = fT.flavor_tag_pythia_particle(j)
						j_flavor = fT.flavor_tag_id(j)
						# or j_flavor = j_flavor_pythia_part.id() 
						# FT_print_psj(j_flavor_psj, '  lead parton')
						dR = j.delta_R(j_flavor_psj)
						dphi = j.delta_phi_to(j_flavor_psj)
						# print('  ', 'jet - leading parton deltaR={:.3f}\tdphi={:.3f}'.format(dR, dphi))
						# print()
						if dR > 0.2:
							continue
						if args.py_hardQCDgluons and j_flavor != 21:
							continue
						if args.py_hardQCDquarks and j_flavor == 21:
							continue
						ev_number_stream = args.nperfile * n_jets_accepted + event_number
						ml_root_ntuple_ev.Fill(run_number, ev_number_stream, pythia.info.sigmaGen(), pythia.info.code(), j_flavor)
						for p in j.constituents():
							if pythiafjext.getPythia8Particle(p).isParton() is False:
								ml_root_ntuple_parts.Fill(run_number, ev_number_stream, p.pt(), p.eta(), p.phi(), pythiafjext.getPythia8Particle(p).id())
						n_jets_accepted += 1
						# print('n jets accepted:', n_jets_accepted)
					else:
						# print('   no parton match'.format(j.perp()))
						pass
			if n_jets_accepted < 1:
					continue
			n_total_jets_accepted += n_jets_accepted
			pbar.update(n_jets_accepted)

	pbar.close()
	pythia.stat()
	pythia.settings.writeFile(
		format_output_file('pythia_test_tagger.cmnd', -1, args))

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
