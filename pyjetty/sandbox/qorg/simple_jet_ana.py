#!/usr/bin/env python

from symbol import file_input
from yaml import parse
from pyjetty.mputils import DataIO
import argparse
import os
import ROOT
import fastjet as fj
def main():
	parser = argparse.ArgumentParser(
	description='analysis of jets from an ntuple', prog=os.path.basename(__file__))
	parser.add_argument('-i', '--input', help='input file', type=str, required=True)
	args = parser.parse_args()
	print(args)
 
	rfout = ROOT.TFile(args.input.replace('.txt', '_jets.root'), 'recreate')
	rfout.cd()
	tnjets = ROOT.TNtuple('tnjets', 'tnjets', 'pt:eta:y:phi:mult:tag:lead_part_pid')

	fj.ClusterSequence.print_banner()
	jet_def = fj.JetDefinition ( fj.antikt_algorithm, 0.8 )

	dfio = DataIO(file_list=args.input, tree_name='tree_Particle_gen', event_tree_name='tree_Event_gen', is_data=False)
	number_of_events = 0
	for e in dfio.next_event(): 
		number_of_events += 1
		jets = fj.sorted_by_pt(jet_def(e.particles))
		for j in jets:
			mult = len(j.constituents())
			try:
				tag = e.partonID
			except:
				tag = 0
			lead_pid = 0
			tnjets.Fill(j.perp(), j.eta(), j.rap(), j.phi(), mult, tag, lead_pid)
	print('[i] total number of events analyzed', number_of_events)

	rfout.Write()
	rfout.Close()
	print('[i] written: ', rfout.GetName())

if __name__ == "__main__":
	main()