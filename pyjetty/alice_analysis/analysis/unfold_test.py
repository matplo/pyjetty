#! /usr/bin/env python

import ROOT
import math

# Load the RooUnfold library
#ROOT.gSystem.Load("/home/software/users/ploskon/RooUnfold-1.1.1/libRooUnfold.so")
ROOT.gSystem.Load("/software/users/james/RooUnfold/libRooUnfold.so")
#ROOT.gSystem.Load("libRooUnfold.so")

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

# Suppress a lot of standard output
#ROOT.gErrorIgnoreLevel = ROOT.kWarning

# Set plotting options
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

###########################################################################################
###########################################################################################
def unfold_test():

  nbins_data = 100
  nbins_truth = 10

  print('Create data and response histograms...')
  hData = ROOT.TH1F("hData", "hData", nbins_data, 0, 10)
  for bin in range(1, nbins_data+1):
    hData.SetBinContent(bin, bin)

  hResponse = ROOT.TH2F("hResponse", "hResponse", nbins_data, 0, 10, nbins_truth, 0, 10)
  for binx in range(1, nbins_data+1):
      for biny in range(1, nbins_truth+1):
          if abs( (binx/10) - biny ) < 1:
            hData.SetBinContent(binx, biny, 1)

  hTrue = ROOT.TH1F("hTrue", "hTrue", nbins_truth, 0, 10)
  for bin in range(1, nbins_truth+1):
    hTrue.SetBinContent(bin, bin*1.1)

  response = ROOT.RooUnfoldResponse(0, hTrue, hResponse, "hResponse", "hResponse")

  print('Unfolding...')
  errorType = ROOT.RooUnfold.kCovToy

  #unfold = ROOT.RooUnfoldSvd(response, hData, 1)
  unfold = ROOT.RooUnfoldBayes(response, hData, 2)

  hUnfolded = unfold.Hreco(errorType)
  print('Done.')
  print('Success!')

#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  unfold_test()
