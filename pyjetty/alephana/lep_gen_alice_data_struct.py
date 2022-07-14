#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
# import fjcontrib
# import fjext

import tqdm
import argparse
import os
import numpy as np

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

from pyjetty.mputils import RTreeWriter
from pyjetty.mputils import pwarning, pinfo, perror


import ROOT
ROOT.gROOT.SetBatch(True)


# helper function
def get_pid(p):
	py8p = pythiafjext.getPythia8Particle(p)
	pid = 0
	if py8p:
		pid = py8p.id()
	return pid

# note: initialization for hadronic events LEP1 - main06.cc

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	_default_output_filename = os.path.basename(__file__).replace(".py", "") + "_output.root"
	parser.add_argument('--output', default=_default_output_filename, type=str)
	parser.add_argument('--user-seed', help='pythia seed', default=1111, type=int)
	parser.add_argument('--debug', default=0, type=int)
	args = parser.parse_args()

	# jets
	# print the banner first
	# fj.ClusterSequence.print_banner()
	# print()
	# # set up our jet definition and a jet selector
	jet_R0 = 1.0
	# jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	# print(jet_def)

	# acceptance
	# hadron level acceptamce
	max_eta_hadron = 10.
	from pyjetty.mputils import pwarning
	pwarning('max eta for particles after hadronization set to', max_eta_hadron)
	parts_selector_h = fj.SelectorAbsEtaMax(max_eta_hadron)
	jet_selector = fj.SelectorPtMin(1.0) & fj.SelectorPtMax(100.0) & fj.SelectorAbsEtaMax(max_eta_hadron - 1.05 * jet_R0)
	# parton level acceptamce
	max_eta_parton = max_eta_hadron + 3. * jet_R0
	pwarning('max eta for partons set to', max_eta_parton)
	parts_selector_p = fj.SelectorAbsEtaMax(max_eta_parton)

	# initialize pythia
	# mZ = pythia8.Pythia().particleData.m0(23)
	mZ = 91.188
	beams_eCM ="Beams:eCM={}".format(mZ)
	mycfg = [	"PDF:lepton = off", # Allow no substructure in e+- beams: normal for corrected LEP data.
				"WeakSingleBoson:ffbar2gmZ = on", # Process selection.
				"23:onMode = off", # Switch off all Z0 decays and then switch back on those to quarks.
				"23:onIfAny = 1 2 3 4 5", 
				"Beams:idA =  11", 
				"Beams:idB = -11", 
				beams_eCM, # LEP1 initialization at Z0 mass.
				"HadronLevel:all=off", # parton level first
				"PhaseSpace:bias2Selection=off"] # this is ON by default in pyconf - not OK for these settings
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if not pythia:
		return

	# initialize ROOT output
	outf = ROOT.TFile(args.output, 'recreate')
	outf.cd()
	tdf = ROOT.TDirectoryFile('PWGHF_TreeCreator', 'PWGHF_TreeCreator')
	tdf.cd()
	t_p = ROOT.TNtuple('tree_Particle_P', 'tree_Particle_P', 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticleID:ParticleIDabs:ParticleCharge:isGluon:isQuark:mult')
	# t_p = ROOT.TNtuple('tree_Particle_gen', 'tree_Particle_gen', 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticleID:ParticleIDabs:ParticleCharge')
	t_h = ROOT.TNtuple('tree_Particle_H', 'tree_Particle_H', 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticleID:ParticleIDabs:ParticleCharge:isHadron:isLepton:isVisible:mult')
	t_e = ROOT.TNtuple('tree_event_char', 'tree_event_char', 'run_number:ev_id:z_vtx_reco:is_ev_rej:multP:multH')

	if args.nev < 100:
		args.nev = 100

	run_number = args.user_seed

	# main loop
	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		if args.debug:
			pwarning('-- event', i)

		#select particles
		parts_pythia_p = fj.sorted_by_pt(pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True))
		parts_pythia_p_selected = parts_selector_p(parts_pythia_p)

		# hadronize
		hstatus = pythia.forceHadronLevel()
		if not hstatus:
			pwarning('forceHadronLevel false event', iev)
			continue

		parts_pythia_h = fj.sorted_by_pt(pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True))
		parts_pythia_h_selected = parts_selector_h(parts_pythia_h)

		# charged hadrons/particles only
		# parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kHadron, pythiafjext.kCharged])
		# parts_pythia_hch = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, True)
		# parts_pythia_hch_selected = parts_selector_h(parts_pythia_hch)


		# stream to trees
		ev_id = i
		multP = len(parts_pythia_p)
		multH = len(parts_pythia_h)
		t_e.Fill(run_number, ev_id, 0, 0, multP, multH)
		_tmp = [
				t_p.Fill(
					run_number, ev_id, p.perp(), p.eta(), p.phi(), 
					pythiafjext.getPythia8Particle(p).id(), 
					pythiafjext.getPythia8Particle(p).idAbs(), 
					pythiafjext.getPythia8Particle(p).charge(),
					pythiafjext.getPythia8Particle(p).isGluon(),
					pythiafjext.getPythia8Particle(p).isQuark(),
					multP
					)
				for p in parts_pythia_p
				]

		_tmp = [
				t_h.Fill(
					run_number, ev_id, p.perp(), p.eta(), p.phi(), 
					pythiafjext.getPythia8Particle(p).id(), 
					pythiafjext.getPythia8Particle(p).idAbs(), 
					pythiafjext.getPythia8Particle(p).charge(),
					pythiafjext.getPythia8Particle(p).isHadron(),
					pythiafjext.getPythia8Particle(p).isLepton(),
					pythiafjext.getPythia8Particle(p).isVisible(),
					multH
					)
				for p in parts_pythia_h
				]

	pythia.stat()

	outf.Write()
	outf.Close()

if __name__ == '__main__':
	main()