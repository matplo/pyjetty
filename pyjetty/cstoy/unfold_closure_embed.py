#!/usr/bin/env python
# ==============================================================================
#  File and Version Information:
#       $Id$
#
#  Description:
#       Simple example usage of the RooUnfold package using toy MC.
#
#  Author: Tim Adye <T.J.Adye@rl.ac.uk>
#  modified by matplo - example loop on tree w/ RTreeReader and unfold
#  note: two files hardcoded - one for response matrix (Response) and different 
# (largely stat indep - different background/smear) for cimparison (Test)
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

from pyjetty.mputils import RTreeReader

### TRAIN
response= ROOT.RooUnfoldResponse (20, -20, 180)
tr = RTreeReader( tree_name='t', 
          branches = ['j_pt', 'ej_pt'],
          file_name='$HOME/devel/pyjetty/pyjetty/cstoy/output_alpha_0_dRmax_0.4_SDzcut_0.2_emb.root')
hTrueA= TH1D ("trueA", "Actual Truth",    20, -20, 180);
hMeasA= TH1D ("measA", "Actual Measured", 20, -20, 180);
for i in range(tr.tree.GetEntries()):
  tr.tree.GetEntry(i)
  if tr.j_pt.size() > 0 and tr.ej_pt.size() > 0:
    response.Fill(tr.ej_pt[0], tr.j_pt[0])
    hTrueA.Fill(tr.j_pt[0])
    hMeasA.Fill(tr.ej_pt[0])
print('[i] trained on', tr.tree.GetEntries(), 'jets')

### TEST with different PbPb data
hTrue= TH1D ("true", "Test Truth",    20, -20, 180);
hMeas= TH1D ("meas", "Test Measured", 20, -20, 180);
trT = RTreeReader( tree_name='t', 
          branches = ['j_pt', 'ej_pt'],
          file_name='$HOME/devel/pyjetty/pyjetty/cstoy/output_alpha_0_dRmax_0.2_SDzcut_0.2_emb.root')
for i in range(trT.tree.GetEntries()):
  trT.tree.GetEntry(i)
  if trT.j_pt.size() > 0:
    hTrue.Fill(trT.j_pt[0])
  if trT.ej_pt.size() > 0:
    hMeas.Fill(trT.ej_pt[0])
print('[i] test against', trT.tree.GetEntries(), 'jets')

if   method == "bayes":
  unfold= ROOT.RooUnfoldBayes   (response, hMeas, 3);    #  OR
elif method == "svd":
  unfold= ROOT.RooUnfoldSvd     (response, hMeas, 4);   #  OR
elif method == "tunfold":
  unfold= ROOT.RooUnfoldTUnfold (response, hMeas);       #  OR
elif method == "ids":
  unfold= ROOT.RooUnfoldIds     (response, hMeas, 3);    #  OR
else:
  print("Unknown method:",method)
  sys.exit(1)

hReco= unfold.Hreco();

unfold.PrintTable (cout, hTrue);

canvas = ROOT.TCanvas("{}_unfold_clodure_embed".format(method),method)

hReco.Draw();
hMeas.Draw("SAME");
hTrue.SetLineColor(8);
hTrue.Draw("SAME");

canvas.SaveAs("{}_unfold_closure_embed.pdf".format(method))

fout = ROOT.TFile('{}_unfold_closure_embed.root'.format(method), 'recreate')
fout.cd()
hReco.SetStats(0)
hReco.SetTitle('Unfold {}'.format(method))
hReco.Write()
hMeasA.Write()
hTrueA.Write()
hMeas.Write()
hTrue.Write()
fout.Close()

print ('[i] written', fout.GetName())