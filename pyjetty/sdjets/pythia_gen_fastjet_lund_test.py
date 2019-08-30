#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os
import numpy as np

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

import ROOT as r
import array

def logbins(xmin, xmax, nbins):
		lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
		arr = array.array('f', lspace)
		return arr


def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	parser.add_argument('-o', '--output', help='output file name', default='pythia_lund_sd.root', type=str)
	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()	

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_R0 = 0.4
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)

	jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
	lund_gen = fjcontrib.LundGenerator(jet_def_lund)
	sd = fjcontrib.SoftDrop(0, 0.1, 1.0)

	jet_selector = fj.SelectorPtMin(10) & fj.SelectorAbsEtaMax(1)
	print(jet_def)

	fout = r.TFile(args.output, 'RECREATE')
	fout.cd()
	tnjets = r.TNtuple('jets', 'jets', 'iev:nj:pt:eta:phi:f5:f6:f5dR:f6dR:proc:xsec')
	tnjets_sd = r.TNtuple('jets_sd', 'jets_sd', 'iev:nj:pt:f5:f6:f5dR:f6dR:proc:xsec:pt_sd:z:dR:mu')
	tnlund = r.TNtuple('lund', 'lund', "iev:nj:pt:f5:f6:f5dR:f6dR:ns:pts:delta:kt")

	nbins = 20
	lbins = logbins(10, 500, nbins)
	hjpt = r.TH1F('hjpt', 'hjpt', nbins, lbins)
	hjptq = r.TH1F('hjptq', 'hjptq', nbins, lbins)
	hjptg = r.TH1F('hjptg', 'hjptg', nbins, lbins)

	hjptsd = r.TH1F('hjptsd', 'hjptsd', nbins, lbins)
	hjptsdq = r.TH1F('hjptsdq', 'hjptsdq', nbins, lbins)
	hjptsdg = r.TH1F('hjptsdg', 'hjptsdg', nbins, lbins)

	hjptsdkt10 	= r.TH1F('hjptsdkt10', 'hjptsdkt10', nbins, lbins)
	hjptsdqkt10 = r.TH1F('hjptsdqkt10', 'hjptsdqkt10', nbins, lbins)
	hjptsdgkt10 = r.TH1F('hjptsdgkt10', 'hjptsdgkt10', nbins, lbins)

	hjptsdkt30 	= r.TH1F('hjptsdkt30', 'hjptsdkt30', nbins, lbins)
	hjptsdqkt30 = r.TH1F('hjptsdqkt30', 'hjptsdqkt30', nbins, lbins)
	hjptsdgkt30 = r.TH1F('hjptsdgkt30', 'hjptsdgkt30', nbins, lbins)

	hjptktsd = r.TH2F('hjptktsd', 'hjptktsd', 50, 0, 500, 50, 0, 500)

	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
	if args.nev < 100:
		args.nev = 100
	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
		parts = pythiafjext.vectorize(pythia, True, -2, 2, False)
		jets = fj.sorted_by_pt(jet_selector(jet_def(parts)))
		f5 = pythia.event[5]
		f6 = pythia.event[6]
		f5psj = fj.PseudoJet(f5.px(), f5.py(), f5.pz(), f5.e())
		f6psj = fj.PseudoJet(f6.px(), f6.py(), f6.pz(), f6.e())
		iproc = pythia.info.code()
		xsec = pythia.info.sigmaGen()
		for nj,j in enumerate(jets):
			f5dR = f5psj.delta_R(j)
			f6dR = f6psj.delta_R(j)
			tnjets.Fill(i, nj, j.pt(), j.eta(), j.phi(), f5.id(), f6.id(), f5dR, f6dR, iproc, xsec) 
			lund = lund_gen.result(j)
			for ns, s in enumerate(lund):
				tnlund.Fill(i, nj, j.pt(), f5.id(), f6.id(), f5dR, f6dR, ns, s.pair().perp(), s.Delta(), s.kt())
			j_sd = sd.result(j)
			sd_info = fjcontrib.get_SD_jet_info(j_sd)
			tnjets_sd.Fill(	i, nj, j.pt(), f5.id(), f6.id(), f5dR, f6dR, iproc, xsec,
							j_sd.pt(),
							sd_info.z, sd_info.dR, sd_info.mu)
			# fill the histograms
			f = None
			if f5dR < 0.4: f = f5
			if f6dR < 0.4: f = f6
			if f is not None:
				kt = sd_info.z * j_sd.perp()
				hjpt.Fill(j.pt())
				if sd_info.z >= 0: 
					hjptsd.Fill(j.pt())
					if kt > 10: hjptsdkt10.Fill(j.pt())
					if kt > 30: hjptsdkt30.Fill(j.pt())
				if f.id() != 21: 
					hjptq.Fill(j.pt())
					if sd_info.z >= 0: 
						hjptsdq.Fill(j.pt())
						if kt > 10: hjptsdqkt10.Fill(j.pt())
						if kt > 30: hjptsdqkt30.Fill(j.pt())
				else:
					hjptg.Fill(j.pt())
					if sd_info.z >= 0: 
						hjptsdg.Fill(j.pt())
						if kt > 10: hjptsdgkt10.Fill(j.pt())
						if kt > 30: hjptsdgkt30.Fill(j.pt())
	pythia.stat()

	fout.cd()
	hjptq_r = hjptq.Clone('hjptq_r')
	hjptq_r.Sumw2()
	hjptq_r.Divide(hjpt)
	hjptg_r = hjptg.Clone('hjptg_r')
	hjptg_r.Sumw2()
	hjptg_r.Divide(hjpt)

	hjptsdq_r = hjptsdq.Clone('hjptsdq_r')
	hjptsdq_r.Sumw2()
	hjptsdq_r.Divide(hjptsd)
	hjptsdg_r = hjptsdg.Clone('hjptsdg_r')
	hjptsdg_r.Sumw2()
	hjptsdg_r.Divide(hjptsd)

	hjptsdqkt10_r = hjptsdqkt10.Clone('hjptsdqkt10_r')
	hjptsdqkt10_r.Sumw2()
	hjptsdqkt10_r.Divide(hjptsdkt10)
	hjptsdgkt10_r = hjptsdgkt10.Clone('hjptsdgkt10_r')
	hjptsdgkt10_r.Sumw2()
	hjptsdgkt10_r.Divide(hjptsdkt10)

	hjptsdqkt30_r = hjptsdqkt30.Clone('hjptsdqkt30_r')
	hjptsdqkt30_r.Sumw2()
	hjptsdqkt30_r.Divide(hjptsdkt30)
	hjptsdgkt30_r = hjptsdgkt30.Clone('hjptsdgkt30_r')
	hjptsdgkt30_r.Sumw2()
	hjptsdgkt30_r.Divide(hjptsdkt30)

	fout.Write()
	fout.Close()
	# print('making lund diagram for all jets...')
	# print('listing lund plane points... Delta, kt - for {} selected jets'.format(len(all_jets)))
	# for l in lunds:
	# 	print ('- jet pT={0:5.2f} eta={1:5.2f}'.format(l[0].pair().perp(), l[0].pair().eta()))
	# 	print ('  Deltas={}'.format([s.Delta() for s in l]))
	# 	print ('  kts={}'.format([s.Delta() for s in l]))
	# 	print ( )

if __name__ == '__main__':
	main()