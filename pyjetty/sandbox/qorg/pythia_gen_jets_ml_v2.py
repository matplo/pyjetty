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

import importlib
root_avail = importlib.util.find_spec("ROOT")
if root_avail is not None:
	import ROOT


class JetFileOutput(object):
	def __init__(self, **kwargs):
		self.args = {}
		if 'args' in kwargs:
			self.args = kwargs['args']
		if type(self.args) is argparse.Namespace:	
			self.args = vars(self.args)
		print(self.args)
		self.sname = []
		for opt in self.args:
			if self.args[opt]:
				self.sname.append(opt)
				if type(self.args[opt]) in [float, int]:
					self.sname.append(str(self.args[opt]))
		self.sname.append('.cmnd')
		self.cmnd_file_name = '_'.join(self.sname)
		self.nfile = 0
  
		self.jet_def_ca = fj.JetDefinition(fj.cambridge_algorithm, fj.JetDefinition.max_allowable_R)
		self.mDT = None
		try: 
			if self.args['mDT']:
				self.mDT = fjcontrib.ModifiedMassDropTagger(self.args['mDTzcut'])
		except:
			pass
  
	def write_file_partons(self, jets, pythia):
		root_file_name = '_'.join(self.sname).replace('.cmnd', '_{}.root'.format(self.nfile))
		root_file = ROOT.TFile(root_file_name, 'recreate')
		root_ntuple_parts = ROOT.TNtuple('tree_Particle_gen', 'particles from PYTHIA8 - jets',
								'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID:ParticleMass')
		root_ntuple_ev = ROOT.TNtuple('tree_Event_gen', 'event info from PYTHIA8 - jets',
								'run_number:ev_id:xsec:code:partonID')
		run_number = self.args['py_seed'] + self.nfile
		if self.args['py_seed'] < 0:
			run_number = self.nfile
		for i, _j in enumerate(jets):
			ev_number_stream = i
			j = _j[0]
			j_flavor = _j[1]
			root_ntuple_ev.Fill(run_number, ev_number_stream, pythia.info.sigmaGen(), pythia.info.code(), j_flavor)
			for p in j.constituents():
				if pythiafjext.getPythia8Particle(p).isParton() is False:
					root_ntuple_parts.Fill(run_number, ev_number_stream, p.pt(), p.eta(), p.phi(), pythiafjext.getPythia8Particle(p).id(), p.m())
		self.nfile += 1
		root_file.Write()
		root_file.Close()

	def write_file_Zjets(self, jets, pythia):
		root_file_name = '_'.join(self.sname).replace('.cmnd', '_{}.root'.format(self.nfile))
		root_file = ROOT.TFile(root_file_name, 'recreate')
		root_ntuple_parts = ROOT.TNtuple('tree_Particle_gen', 'particles from PYTHIA8 - jets', 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID:ParticleMass')
		root_ntuple_parts_mDT = ROOT.TNtuple('tree_Particle_gen_mDT{}'.format(self.args['mDTzcut']), 'particles from PYTHIA8 - jets', 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID:ParticleMass')
		root_ntuple_ev = ROOT.TNtuple('tree_Event_gen', 'event info from PYTHIA8 - jets', 'run_number:ev_id:xsec:code')
		root_ntuple_zjet = ROOT.TNtuple('tree_Zjet_gen', 'Zjet kinematics','run_number:ev_id:ZjetPt:ZjetEta:ZjetPhi:ParticlePt:ParticleEta:ParticlePhi:ParticlePID:dR')
		run_number = self.args['py_seed'] + self.nfile
		if self.args['py_seed'] < 0:
			run_number = self.nfile
		for i, _j in enumerate(jets):
			ev_number_stream = i
			Zjet = _j[0]
			psjZ = _j[1]
			ZjetDR = psjZ.delta_R(Zjet)
			root_ntuple_ev.Fill(run_number, ev_number_stream, pythia.info.sigmaGen(), pythia.info.code())
			root_ntuple_zjet.Fill(run_number, ev_number_stream, Zjet.pt(), Zjet.eta(), Zjet.phi(), psjZ.pt(), psjZ.eta(), psjZ.phi(), pythiafjext.getPythia8Particle(psjZ).id(), ZjetDR)
			for p in Zjet.constituents():
				root_ntuple_parts.Fill(run_number, ev_number_stream, p.pt(), p.eta(), p.phi(), pythiafjext.getPythia8Particle(p).id(), p.m())

			# now mDT jets
			ca_jets = fj.sorted_by_pt(self.jet_def_ca(Zjet.constituents()))
			if (len(ca_jets) > 1):
				print('[w] CA reclustering returns more than 1 subjet (?)')
			else:
				# Use modified mass-drop tagger to clean up jet.
				taggedJet = self.mDT.result(ca_jets[0])
				if taggedJet.has_constituents():
					for p in taggedJet.constituents():
						root_ntuple_parts_mDT.Fill(run_number, ev_number_stream, p.pt(), p.eta(), p.phi(), pythiafjext.getPythia8Particle(p).id(), p.m())
		self.nfile += 1
		root_file.Write()
		root_file.Close()


def main():
	parser = argparse.ArgumentParser(
		description='pythia8 in python', prog=os.path.basename(__file__))
	parser.add_argument('--hepmc', help='write the hepmc files',
	                    default=False, action='store_true')
	parser.add_argument('--jet-R', help='specify jet R for flavor tagging', default=0.8, type=float)
	parser.add_argument('-c', '--nperfile', help='nevents per output file', default=1000, type=int)
	parser.add_argument('--jet-ptmin', help='minimum pT cut on jets', default=500., type=float)
	parser.add_argument('--jet-ptmax', help='maximum pT cut on jets', default=550., type=float)
	parser.add_argument('--Zjet', help='force Zjet generation - will read pythia_gen_qorg_Zjet_master.cmnd from current directory', default=False, action='store_true')
	parser.add_argument('--mDT', help='mDT clean up with z_cut=0.04 (default) on the Zjets', default=False, action='store_true')
	parser.add_argument('--mDTzcut', help='mDT clean up with z_cut=0.04 (default) on the Zjets', default=0.04, type=float)
	parser.add_argument('--ZjetR', help='specify the radius of the anti-kT jet for Z-match - default is R=1.0', default=1.0, type=float)

	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()

	part_selector = fj.SelectorAbsRapMax(4.)
	part_selector_p = fj.SelectorAbsRapMax(4.)
	jet_selector = fj.SelectorPtMin(args.jet_ptmin) * fj.SelectorPtMax(args.jet_ptmax) * fj.SelectorAbsRapMax(1.7)
 
	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 10:
		args.nev = 10

	jet_def_akt = fj.JetDefinition ( fj.antikt_algorithm, args.ZjetR)
	if args.Zjet is False:
		del args.mDT
		del args.mDTzcut
		del args.ZjetR
	else:
		if args.py_cmnd is None:
			print('[e] you have to specify the command file for Zjet generation with --py-cmnd <file> - sorry, no shortcut ;-(')
			return

	if args.py_bias is False:
		del args.py_biaspow
		del args.py_biasref
  
	fout = JetFileOutput(args=args)
 
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
   
		parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kVisible], 0, True)
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
			jet_def = fj.JetDefinition(fj.antikt_algorithm, args.jet_R)
			jets = fj.sorted_by_pt(jet_def(parts))

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
						jets_accepted.append([j, j_flavor])
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

	pbar.close()
	pythia.stat()
	pythia.settings.writeFile(fout.cmnd_file_name)

if __name__ == "__main__":
	main()