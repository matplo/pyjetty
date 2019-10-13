#!/usr/bin/env python
# ==============================================================================
#  File and Version Information:
#       $Id$
#
#  Description:
#       Simple example usage of the RooUnfold package using toy MC.
#
#  Author: Tim Adye <T.J.Adye@rl.ac.uk>
#
# ==============================================================================

import sys
method = "bayes"
if len(sys.argv) > 1: method = sys.argv[1]

from ROOT import gRandom, TH1, TH1D, TCanvas, cout
import ROOT
ROOT.gSystem.Load('$ROOUNFOLDDIR/libRooUnfold')
# ==============================================================================
#  Gaussian smearing, systematic translation, and variable inefficiency
# ==============================================================================

def smear(xt):
  xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  #  efficiency
  x= gRandom.Rndm();
  if x>xeff: return None;
  xsmear= gRandom.Gaus(-2.5,0.2);     #  bias and smear
  return xt+xsmear;

# ==============================================================================
#  Example Unfolding
# ==============================================================================

ROOT.gROOT.SetBatch(True)

response= ROOT.RooUnfoldResponse (40, -10.0, 10.0);

#  Train with a Breit-Wigner, mean 0.3 and width 2.5.
for i in range(100000):
  xt= gRandom.BreitWigner (0.3, 2.5);
  x= smear (xt);
  if x!=None:
    response.Fill (x, xt);
  else:
    response.Miss (xt);

hTrue= TH1D ("true", "Test Truth",    40, -10.0, 10.0);
hMeas= TH1D ("meas", "Test Measured", 40, -10.0, 10.0);
#  Test with a Gaussian, mean 0 and width 2.
for i in range(10000):
  xt= gRandom.Gaus (0.0, 2.0)
  x= smear (xt);
  hTrue.Fill(xt);
  if x!=None: hMeas.Fill(x);

if   method == "bayes":
  unfold= ROOT.RooUnfoldBayes   (response, hMeas, 4);    #  OR
elif method == "svd":
  unfold= ROOT.RooUnfoldSvd     (response, hMeas, 20);   #  OR
elif method == "tunfold":
  unfold= ROOT.RooUnfoldTUnfold (response, hMeas);       #  OR
elif method == "ids":
  unfold= ROOT.RooUnfoldIds     (response, hMeas, 3);    #  OR
else:
  print("Unknown method:",method)
  sys.exit(1)

hReco= unfold.Hreco();

unfold.PrintTable (cout, hTrue);

canvas = ROOT.TCanvas("RooUnfold",method)

hReco.Draw();
hMeas.Draw("SAME");
hTrue.SetLineColor(8);
hTrue.Draw("SAME");

canvas.SaveAs("RooUnfold.pdf")
