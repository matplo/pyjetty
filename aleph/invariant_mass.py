#!/usr/bin/env python3
import pythia8
from mptools import pymptools as mpt
import fastjet as fj
from lundplane import pylundplane as lund
from tqdm import tqdm
import sys
import ROOT as r
import aleph_utils
import argparse
import os


def main(args):
	aleph_file=args.input
	if not os.path.isfile(aleph_file):
		print('[e] input does not exists')
		return
	outfile = 'invariant_mass.root'
	if args.oxo:
		outfile = '{}_mass.root'.format(os.path.basename(aleph_file))
	fout = r.TFile(outfile, 'recreate')
	fout.cd()
	hwflag = r.TH1I('hwflag', 'hwflag', 10, 0, 10)
	hwflagmass 		= r.TH2F('hwflag_mass', 	'hwflag_mass', 		10, 0, 10, 20, 0, 2)
	hmass2pK0 		= r.TH1F('hmass2pK0', 		'hmass2pK0', 		1000, 0, 3.0)
	hmass2p 		= r.TH1F('hmass2p', 		'hmass2p', 			530, 1.2, 6.5)
	hmasspt2p 		= r.TH2F('hmasspt2p', 		'hmasspt2p', 		530, 1.2, 6.5, 100, 0, 10)
	hmass3p 		= r.TH1F('hmass3p', 		'hmass3p', 			530, 1.2, 6.5)
	hmass2pK0Ldphi 	= r.TH2F('hmass2pK0Ldphi', 	'hmass2pK0Ldphi', 	100, 0, 3.0, 36, 0, r.TMath.Pi())
	hntracks 		= r.TH1F('hntracks', 		'hntracks', 		100, 0, 100)
	tracks 	= r.TNtuple('tr', 'tr', 'q:pt:phi:eta:pwflag')
	tK0c 	= r.TNtuple('tK0c', 'tK0c', 'lid:lq:lpt:dphi:k0pt:k0m')
	_p1 = r.TLorentzVector()
	_p2 = r.TLorentzVector()
	_p3 = r.TLorentzVector()
	nev = aleph_utils.get_n_events(aleph_file)
	if int(args.nevents) > -1:
		if args.nevents < nev:
			nev = int(args.nevents)
	pbar = tqdm(total=nev)
	reader = mpt.Reader(aleph_file)
	while reader.read_next_event():
		pbar.update()
		e = reader.get_event()
		parts = []
		inparts = e.get_particles()
		nparts = len(inparts)

		leptons = []
		for i in range(nparts):
			l = inparts[i]
			if l.pwflag() == aleph_utils.AlephWFLAG.ALEPH_CHARGED_LEPTONS1 or l.pwflag() == aleph_utils.AlephWFLAG.ALEPH_CHARGED_LEPTONS2:
				_l = r.TLorentzVector()
				_l.SetPxPyPzE(l.px(), l.py(), l.pz(), l.e())
				leptons.append([_l, l])
		ntracks = 0
		for ip1 in range(nparts):
			p1 = inparts[ip1]
			hwflag.Fill(p1.pwflag())
			hwflagmass.Fill(p1.pwflag(), p1.m())
			if p1.pwflag() != aleph_utils.AlephWFLAG.ALEPH_CHARGED_TRACK:
				continue
			_p1.SetPxPyPzE(p1.px(), p1.py(), p1.pz(), p1.e())
			# tracks 	= t.TNtuple('tr', 'tr', 'q:pt:phi:eta:pwflag')
			tracks.Fill(p1.q(), _p1.Pt(), _p1.Phi(), _p1.Eta(), p1.pwflag())
			ntracks += 1
			# print("compare mass = {} ? {}".format(p1.m(), _p1.M()))
			for ip2 in range(ip1+1, nparts):
				p2 = inparts[ip2]
				if p2.pwflag() != aleph_utils.AlephWFLAG.ALEPH_CHARGED_TRACK:
					continue
				_p2.SetPxPyPzE(p2.px(), p2.py(), p2.pz(), p2.e())
				sump = _p1 + _p2
				hmass2p.Fill(sump.M())
				pK0 = None
				if p1.q() != p2.q():
					hmass2pK0.Fill(sump.M())
					hmasspt2p.Fill(sump.M(), sump.Pt())
					pK0 = sump
					for _l in leptons:
						if _l[0].Pt() > 2:
							hmass2pK0Ldphi.Fill(sump.M(), abs(sump.DeltaPhi(_l[0])))
						# tK0c 	= r.TNtuple('tK0c', 'tK0c', 'lid:lq:lpt:dphi:k0pt:k0m')
						tK0c.Fill(_l[1].pwflag(), _l[1].q(), _l[0].Pt(), abs(sump.DeltaPhi(_l[0])), sump.Pt(), sump.M())
				if pK0:
					for ip3 in range(ip2+1, nparts):
						p3 = inparts[ip3]
						if p3.pwflag() != aleph_utils.AlephWFLAG.ALEPH_CHARGED_TRACK:
							continue
						sump3 = pK0 + _p3
						hmass3p.Fill(sump3.M())
		hntracks.Fill(ntracks)
		if pbar.n >= nev:
			break
	pbar.close()
	fout.Write()
	print('[i] written', fout.GetName())

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='some analysis on aleph data',
	                                 prog=os.path.basename(__file__))
	parser.add_argument('-i', '--input', help='root file with set0', type=str, required=True)
	parser.add_argument('-n', '--nevents', help='user number of events', type=int, required=False, default=-1)
	parser.add_argument('--oxo', help='follow input structure for output', required=False, default=False, action="store_true")
	args = parser.parse_args()
	main(args)
