#!/usr/bin/env python

import ROOT

from pyjetty.mputils import BoltzmannEvent
max_eta=0.9
be = BoltzmannEvent(mean_pt=0.6, multiplicity=2000 * max_eta * 2, max_eta=max_eta, max_pt=100)
bg_parts = be.generate(offset=10000)

funbg = ROOT.TF1("funbg", "[1] * 2. / [0] * x * TMath::Exp(-(2. / [0]) * x)", 0, 3, 2)
funbg.SetParameter(0, 0.73)
funbg.SetParameter(1, 1.)
#funbg.Draw()

funbgX = ROOT.TF1("funbg", "[1] * 2. / [0] * x * TMath::Exp(-(2. / [0]) * x)", 0, 3, 2)
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
