#!/usr/bin/env python

from symbol import file_input
# from yaml import parse
from pyjetty.mputils.data_io import DataIO
import argparse
import os
import ROOT
import fastjet as fj
import fjcontrib

def main():
	parser = argparse.ArgumentParser(
	description='analysis of jets from an ntuple', prog=os.path.basename(__file__))
	parser.add_argument('-i', '--input', help='input file', type=str, required=True)
	parser.add_argument('-o', '--output', help='output filename - optional', type=str, required=False)
	parser.add_argument('--ptree', help='particle tree to analyze', type=str, default='tree_Particle_gen')
	parser.add_argument('--etree', help='event tree to analyze', type=str, default='tree_Event_gen')
	parser.add_argument('--jet-R', help='jet R', type=float, default=1.0)
	args = parser.parse_args()
	print(args)

	outfname =  args.input.replace('.txt', '_jets.root')
	if args.output:
		outfname = args.output
	rfout = ROOT.TFile(outfname, 'recreate')
	rfout.cd()
	tnjets = ROOT.TNtuple('tnjets', 'tnjets', 'pt:eta:y:phi:m:mult:tag:lead_part_pid:tau21:n')
	tnstray = ROOT.TNtuple('tnstray', 'tnstray', 'pt:eta:y:dphi:deta:dy:dR:z')

	fj.ClusterSequence.print_banner()
	jet_def = fj.JetDefinition ( fj.antikt_algorithm, args.jet_R)

	beta = 1.0
	nSub1_beta1 = fjcontrib.Nsubjettiness(1,   fjcontrib.OnePass_WTA_KT_Axes(), fjcontrib.UnnormalizedMeasure(beta))
	nSub2_beta1 = fjcontrib.Nsubjettiness(2,   fjcontrib.OnePass_WTA_KT_Axes(), fjcontrib.UnnormalizedMeasure(beta))

	dfio = DataIO(file_list=args.input, tree_name=args.ptree, event_tree_name=args.etree, is_data=False)
	number_of_events = 0
	for e in dfio.next_event(): 
		number_of_events += 1
		jets = fj.sorted_by_pt(jet_def(e.particles))
		for i, j in enumerate(jets):
			mult = len(j.constituents())
			try:
				tag = e.partonID
			except:
				tag = 0
			lead_pid = 0
			tau1 = nSub1_beta1.result(j)
			tau2 = nSub2_beta1.result(j)
			tau_ratio = -1
			if tau1 != 0:
				tau_ratio = tau2/tau1
			tnjets.Fill(j.perp(), j.eta(), j.rap(), j.phi(), j.m(), mult, tag, lead_pid, tau_ratio, i)
			if i > 0:
				for c in j.constituents():
					dphi = jets[0].delta_phi_to(c)
					deta = jets[0].eta() - c.eta()
					dy = jets[0].rap() - c.rap()
					dR = jets[0].delta_R(c)
					tnstray.Fill(c.pt(), c.eta(), c.rap(), dphi, deta, dy, dR, c.pt() / jets[0].pt())

	print('[i] total number of events analyzed', number_of_events)

	rfout.Write()
	rfout.Close()
	print('[i] written: ', rfout.GetName())

if __name__ == "__main__":
	main()