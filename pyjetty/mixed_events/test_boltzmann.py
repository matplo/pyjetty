#!/usr/bin/env python

import sys
import tqdm
import ROOT
from pyjetty.mputils import BoltzmannEvent, BoltzmannSubtractor

import fastjet as fj
import fjtools

def test1():
	max_eta=0.9
	# be = BoltzmannEvent(mean_pt=0.6, multiplicity=2000 * max_eta * 2, max_eta=max_eta, max_pt=100)
	be = fjtools.BoltzmannBackground(0.73, 0.0, 5.)
	bg_parts = be.generate(20000, max_eta=1., offset=10000)

	# funbg = ROOT.TF1("funbg", "[1] * 2. / [0] * x * TMath::Exp(-(2. / [0]) * x)", 0, 3, 2)
	funbg = ROOT.TF1("funbg", "[1] * x * TMath::Exp(-(x / ([0]/2.)))", 0, 3, 2)
	funbg.SetParameter(0, 0.73)
	funbg.SetParameter(1, 1.)
	#funbg.Draw()

	# funbgX = ROOT.TF1("funbg", "[1] * 2. / [0] * x * TMath::Exp(-(2. / [0]) * x)", 0, 3, 2)
	funbgX = ROOT.TF1("funbgX", "[1] * x * TMath::Exp(-(x / ([0]/2.)))", 0, 3, 2)
	funbgX.SetParameter(0, 0.3)
	funbgX.SetParameter(1, 1.)

	h = ROOT.TH1F("h", "h", 400, -10, 10)
	hns = ROOT.TH1F("hns", "hns", 400, -10, 10)
	#for i in range(1000):
	  # h.Fill(funbgX.GetRandom())
	total_pt = sum([p.pt() for p in bg_parts])
	for p in bg_parts:
	  h.Fill(p.pt())
	  hns.Fill(p.pt())
	print(h.Integral(), total_pt)
	h.Scale(1./h.GetEntries())

	tc = ROOT.TCanvas()
	tc.Divide(1,2)
	tc.cd(1)
	h.Draw()
	h.Fit(funbg, "RMN")
	funbg.Draw("same")
	ROOT.gPad.SetLogy()
	ROOT.gPad.SetGridx()
	ROOT.gPad.SetGridy()
	ROOT.gPad.Update()

	hc = ROOT.TH1F("hc", "hc", 400, -10, 10)
	for p in bg_parts:
	  fraction = funbg.Eval(p.pt()) * total_pt
	  hc.Fill(p.pt() - fraction)
	# hc.Scale(1./h.GetEntries())
	hc.SetLineColor(ROOT.kRed)
	tc.cd(2)
	hns.Draw()
	hc.Draw("same")
	ROOT.gPad.SetLogy()

	print(funbg.Eval(1))
	print(funbg.Eval(4))

	ROOT.gPad.SetGridx()
	ROOT.gPad.SetGridy()
	ROOT.gPad.Update()
	ROOT.gApplication.Run()

def test2():
	max_eta=0.9
	# be = BoltzmannEvent(mean_pt=0.6, multiplicity=2000 * max_eta * 2, max_eta=max_eta, max_pt=100)
	be = fjtools.BoltzmannBackground(0.73, 0.0, 5.)

	bs = BoltzmannSubtractor()

	nevents = 100

	h = ROOT.TH1F("h", "h", 400, -10, 10)
	hs = ROOT.TH1F("hs", "hs", 400, -10, 10)
	hs2 = ROOT.TH1F("hs2", "hs2", 400, -10, 10)

	hmult = ROOT.TH1F('hmult', 'hmult', nevents, 0, nevents)
	hfpar0 = ROOT.TH1F("hfpar0", "hfpar0", nevents, 0, nevents)
	hfpar1 = ROOT.TH1F("hfpar1", "hfpar1", nevents, 0, nevents)

	for i in tqdm.tqdm(range(nevents)):
		#bg_parts = be.generate(offset=10000)
		bg_parts = be.generate(20000, max_eta=1., offset=10000)
		for p in bg_parts:
			h.Fill(p.pt())
		subtr_parts_b = bs.subtracted_particles(bg_parts)
		for p in subtr_parts_b:
			hs.Fill(p.pt())
		subtr_parts_b2 = be.subtract(bg_parts)
		for p in subtr_parts_b2:
			hs2.Fill(p.pt())
		hmult.Fill(i, len(bg_parts))
		hfpar0.Fill(i, bs.funbg.GetParameter(0))
		hfpar1.Fill(i, bs.funbg.GetParameter(1))

	tc = ROOT.TCanvas('tc', 'tc', 800, 800)
	tc.Divide(2,2)
	tc.cd(1)
	nentries = hmult.GetEntries()
	h.Scale(1./nentries)
	h.Draw()
	hs.Scale(1./nentries)
	hs.SetLineColor(ROOT.kRed)
	hs.Draw("same")
	hs2.Scale(1./nentries)
	hs2.SetLineColor(ROOT.kBlue)
	hs2.Draw("same")
	ROOT.gPad.SetLogy()
	ROOT.gPad.SetGridx()
	ROOT.gPad.SetGridy()

	tc.cd(2)
	hmult.Draw('hist')

	tc.cd(3)
	hfpar0.SetMinimum(0)
	hfpar0.Draw('hist')
	tc.cd(4)
	hfpar1.SetMinimum(0)
	hfpar1.Draw('hist')

	tc.Update()
	ROOT.gApplication.Run()

def main():
	if len(sys.argv) < 2:
		return
	if sys.argv[1] == '1':
		test1()
	if sys.argv[1] == '2':
		test2()

if __name__ == '__main__':
	main()

# so what to do:
# - fit boltzmann to every event - mean pT and multiplicity
# - cut offs are:
# -- what to take for mean mult
# -- range of the boltzmann fit
# precision is number of bins in the histogram to fit
# subtract the particle according to the boltmann value - probability...


# better algo:
# - put particles into a grid
# - fit the grid
# - for every grid - know the amount of BG to subtract - random from the boltzmann 
# - for every cell remove energy particle by particle until all subtracted - start with softest
