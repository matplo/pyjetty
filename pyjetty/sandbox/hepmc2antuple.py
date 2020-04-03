#!/usr/bin/env python

from __future__ import print_function

import os
import argparse
import tqdm

import pyhepmc_ng

import ROOT
ROOT.gROOT.SetBatch(True)

from pyjetty.mputils import *

# jit improves execution time by 18% - tested with jetty pythia8 events
# BUT produces invalid root file
# from numba import jit
# @jit

def fill_event(run_number, ev_id, event_hepmc, tw_e, tw_p, pdg):
	tw_e.fill_branch('run_number', run_number)
	tw_e.fill_branch('ev_id', ev_id)
	tw_e.fill_branch('z_vtx_reco', 0)
	tw_e.fill_branch('is_ev_rej', 0)
	tw_e.fill_tree()

	for part in event_hepmc.particles:
		if part.status == 1 and not part.end_vertex and pdg.GetParticle(part.pid).Charge() != 0:
			# print(pdg.GetParticle(p.pid).GetName())
			# tlv = ROOT.TLorentzVector()
			# tlv.SetPxPyPzE(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e)
			tw_p.fill_branch('run_number', run_number)
			tw_p.fill_branch('ev_id', ev_id)
			tw_p.fill_branch('ParticlePt', part.momentum.pt())
			tw_p.fill_branch('ParticleEta', part.momentum.eta())
			tw_p.fill_branch('ParticlePhi', part.momentum.phi())
			tw_p.fill_branch('ParticlePID', part.pid)
			tw_p.fill_tree()


def main():
	parser = argparse.ArgumentParser(description='hepmc to ALICE Ntuple format', prog=os.path.basename(__file__))
	parser.add_argument('-i', '--input', help='input file', default='', type=str, required=True)
	parser.add_argument('-o', '--output', help='output root file', default='', type=str, required=True)
	parser.add_argument('--as-data', help='write as data - tree naming convention', action='store_true', default=False)
	parser.add_argument('--hepmc', help='what format 2 or 3', default=2, type=int)
	parser.add_argument('--nev', help='number of events', default=-1, type=int)
	args = parser.parse_args()	

	if args.hepmc == 3:
		input_hepmc = pyhepmc_ng.ReaderAscii(args.input)
	if args.hepmc == 2:
		input_hepmc = pyhepmc_ng.ReaderAsciiHepMC2(args.input)

	if input_hepmc.failed():
		print ("[error] unable to read from {}".format(args.input))
		sys.exit(1)

	outf = ROOT.TFile(args.output, 'recreate')
	outf.cd()
	tdf = ROOT.TDirectoryFile('PWGHF_TreeCreator', 'PWGHF_TreeCreator')
	tdf.cd()
	if args.as_data:
		t_p = ROOT.TTree('tree_Particle', 'tree_Particle')
	else:
		t_p = ROOT.TTree('tree_Particle_gen', 'tree_Particle_gen')
	t_e = ROOT.TTree('tree_event_char', 'tree_event_char')
	tw_p = RTreeWriter(tree=t_p)
	tw_e = RTreeWriter(tree=t_e)

	# run number will be a double - file size in MB
	run_number = os.path.getsize(args.input) / 1.e6
	ev_id = 0

	# unfortunately pyhepmc_ng does not provide the table
	# pdt = pyhepmc_ng.ParticleDataTable()
	# use ROOT instead
	pdg = ROOT.TDatabasePDG()

	event_hepmc = pyhepmc_ng.GenEvent()

	if args.nev > 0:
		pbar = tqdm.tqdm(range(args.nev))
	else:
		pbar = tqdm.tqdm()

	while not input_hepmc.failed():
		ev = input_hepmc.read_event(event_hepmc)
		if input_hepmc.failed():
			break

		fill_event(run_number, ev_id, event_hepmc, tw_e, tw_p, pdg)

		ev_id = ev_id + 1
		pbar.update()
		if args.nev > 0 and ev_id > args.nev:
			break

	outf.Write()
	outf.Close()

if __name__ == '__main__':
	main()