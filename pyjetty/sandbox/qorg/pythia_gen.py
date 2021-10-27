#!/usr/bin/env python

from __future__ import print_function

import fastjet
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
	parser = argparse.ArgumentParser(description='pythia8 in python', prog=os.path.basename(__file__))
	parser.add_argument('--hepmc', help='write the hepmc files', default=False, action='store_true')
	parser.add_argument('--ml', help='write the ml format root files', default=False, action='store_true')
	parser.add_argument('-c', '--nperfile', help='nevents per output file', default=1000, type=int)
	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()

	if args.hepmc:
		hepmc2_output_base = "pythia_gen_qor.hepmc"
		hepmc2_output_name = None
		hepmc2_writer = None
		hepmc2_fileno = 0

	if args.ml:
		if root_avail is None:
			print('[error] unable to load ROOT so --ml is defunct')
			args.ml = False
		ml_root_output_base = "pythia_gen_qor.root"
		ml_root_output_name = None
		ml_root_file = None
		ml_root_ntuple_parts = None
		ml_root_ntuple_ev = None
		ml_root_fileno = 0

	if args.hepmc is False and args.ml is False:
		print('[error] one of the --hepmc or --ml is required')
		print('		   --help for more options')
		return

	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 10:
		args.nev = 10
	run_number = args.py_seed
	event_number = 0
	pbar = tqdm.tqdm(range(args.nev))
	for i in pbar:
		if not pythia.next():
			pbar.update(-1)
			continue
		event_number = event_number + 1
		if i == 0 or i % args.nperfile == 0:
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
					ml_root_ntuple = None
				ml_root_output_name = format_output_file(ml_root_output_base, ml_root_fileno, args)
				ml_root_file = ROOT.TFile(ml_root_output_name, 'recreate')
				ml_root_ntuple_parts = ROOT.TNtuple('tree_Particle_gen', 'particles from PYTHIA8 - {} jets'.format(parton_type_from_args(args)), 
								'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID')
				ml_root_ntuple_ev = ROOT.TNtuple('tree_Event_gen', 'event info from PYTHIA8 - {} jets'.format(parton_type_from_args(args)), 
								'run_number:ev_id:xsec:code')
				ml_root_fileno = ml_root_fileno + 1

		if args.hepmc:
			hepmc2_writer.fillEvent(pythia)

		if args.ml:
			parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)
			ml_root_ntuple_ev.Fill(run_number, event_number, pythia.info.sigmaGen(), pythia.info.code())
			for p in parts_pythia_h:
				ml_root_ntuple_parts.Fill(run_number, event_number, p.pt(), p.eta(), p.phi(), pythiafjext.getPythia8Particle(p).id())

	pythia.stat()
	pythia.settings.writeFile(format_output_file('pythia_gen_qor.cmnd', -1, args))

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
