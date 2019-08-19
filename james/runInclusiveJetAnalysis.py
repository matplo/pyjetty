#!/usr/bin/env python3

"""
  Simple analysis script to read a ROOT TTree of track information
  and do jet-finding, and plot some basic histograms.
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse
import math
import time

# Data analysis and plotting
import uproot
import pandas
import numpy as np
import ROOT

# Fastjet via python (from external library fjpydev)
import pyjetty
import fastjet as fj
from recursivetools import pyrecursivetools as rt

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

# Set debug level (0 = no debug info, 1 = some debug info, 2 = all debug info)
debugLevel = 0

#---------------------------------------------------------------
def runInclusiveJetAnalysis(inputFileData, inputFileMC, outputDir, fileFormat):

  start_time = time.time()
  
  # Create output dir
  if not outputDir.endswith("/"):
    outputDir = outputDir + "/"
  if not os.path.exists(outputDir):
    os.makedirs(outputDir)

  # Convert ROOT TTree to pandas dataframe
  # track_df is a dataframe with one row per jet constituent: run_number, ev_id, ParticlePt, ParticleEta, ParticlePhi
  print('--- {} seconds ---'.format(time.time() - start_time))
  print('Convert ROOT trees to pandas dataframes...')
  track_df = load_dataframe(inputFile)

  # Transform the track dataframe into a SeriesGroupBy object of fastjet particles per event
  print('--- {} seconds ---'.format(time.time() - start_time))
  print('Transform the track dataframe into a series object of fastjet particles per event...')

  # (i) Group the track dataframe by event
  #     track_df_grouped is a DataFrameGroupBy object with one track dataframe per event
  track_df_grouped = track_df.groupby(['run_number','ev_id'])

  # (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet particles
  df_fjparticles = track_df_grouped.apply(get_fjparticles)
  
  if debugLevel > 0:
    print(df_fjparticles.dtypes)
    print(df_fjparticles)
  print('--- {} seconds ---'.format(time.time() - start_time))

  # Print number of events
  nEvents = track_df_grouped.size().count()
  print('Number of events: {}'.format(nEvents))
  nTracks = len(track_df.index)
  print('Number of tracks: {}'.format(nTracks))

  # Initialize histogram dictionary
  hDict = initializeHistograms()
  
  # Find jets and fill histograms
  print('Find jets...')
  analyzeEvents(df_fjparticles, hDict, outputDir, fileFormat)

  # Plot histograms
  print('Plot histograms...')
  plotHistograms(hDict, outputDir, fileFormat)

  print('--- {} seconds ---'.format(time.time() - start_time))

#---------------------------------------------------------------
def load_dataframe(inputFile):
  
  # Load event tree into dataframe, and apply event selection
  event_tree_name = 'PWGHF_TreeCreator/tree_event_char'
  event_tree = uproot.open(inputFile)[event_tree_name]
  if not event_tree:
    print('Tree {} not found in file {}'.format('tree_event_char', inputFile))
  event_df_orig = event_tree.pandas.df(['run_number', 'ev_id', 'z_vtx_reco','is_ev_rej'])
  event_df_orig.reset_index(drop=True)
  event_df = event_df_orig.query('is_ev_rej == 0')
  event_df.reset_index(drop=True)

  # Load track tree into dataframe
  track_tree_name = 'PWGHF_TreeCreator/tree_Particle'
  track_tree = uproot.open(inputFile)[track_tree_name]
  if not track_tree:
    print('Tree {} not found in file {}'.format('tree_Particle', inputFile))
  track_df_orig = track_tree.pandas.df()

  # Merge event info into track tree
  track_df = pandas.merge(track_df_orig, event_df, on=['run_number', 'ev_id'])
  return track_df

#---------------------------------------------------------------
def get_fjparticles(df_tracks):
  
  # Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of pt,eta,phi
  fj_particles = pyjetty.vectorize_pt_eta_phi(df_tracks['ParticlePt'].values, df_tracks['ParticleEta'].values, df_tracks['ParticlePhi'].values)
  
  return fj_particles

#---------------------------------------------------------------
def initializeHistograms():
  
  hDict = {}
  
  hJetPt = ROOT.TH1F('hJetPt', 'hJetPt', 200, 0, 200)
  hJetPt.GetXaxis().SetTitle('p_{T,jet}')
  hJetPt.GetYaxis().SetTitle('dN/dp_{T}')
  hDict['hJetPt'] = hJetPt

  hAJ = ROOT.TH2F('hAJ', 'hAJ', 100, 0, 1., 100, 0, 4.)
  hAJ.GetXaxis().SetTitle('A_{J}')
  hAJ.GetYaxis().SetTitle('#Delta #phi')
  hDict['hAJ'] = hAJ

  hZg = ROOT.TH1F('hZg', 'hZg', 100, 0, 1.)
  hZg.GetXaxis().SetTitle('z_{g}')
  hZg.GetYaxis().SetTitle('dN/dz_{g}')
  hDict['hZg'] = hZg

  hRg = ROOT.TH1F('hRg', 'hRg', 100, 0, 1.)
  hRg.GetXaxis().SetTitle('R_{g}')
  hRg.GetYaxis().SetTitle('dN/dR_{g}')
  hDict['hRg'] = hRg

  return hDict

#---------------------------------------------------------------
def analyzeEvents(df_fjparticles, hDict, outputDir, fileFormat):
  
  fj.ClusterSequence.print_banner()
  print()
  
  # Set jet definition and a jet selector
  jetR = 0.4
  jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
  jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
  print('jet definition is:', jet_def)
  print('jet selector is:', jet_selector,'\n')
  
  # Define SoftDrop settings
  beta = 0
  zcut = 0.1
  sd = rt.SoftDrop(beta, zcut, jetR)
  print('SoftDrop groomer is: {}'.format(sd.description()));

  # Use list comprehension to do jet-finding and fill histograms
  result = [analyzeJets(fj_particles, jet_def, jet_selector, sd, hDict) for fj_particles in df_fjparticles]

#---------------------------------------------------------------
def analyzeJets(fj_particles, jet_def, jet_selector, sd, hDict):
  
  # Do jet finding
  cs = fj.ClusterSequence(fj_particles, jet_def)
  jets = fj.sorted_by_pt(cs.inclusive_jets())
  jets_accepted = jet_selector(jets)

  fillJetHistograms(hDict, jets_accepted)

  # Loop through jets and perform SoftDrop grooming
  jets_sd = []
  for jet in jets_accepted:
    jets_sd.append(sd.result(jet))

  fillSoftDropHistograms(hDict, jets_sd)

#---------------------------------------------------------------
def fillJetHistograms(hDict, jets_accepted):

  # Loop through jets, and fill histograms
  for jet in jets_accepted:
    if debugLevel > 1:
      print('jet:')
      print(jet)
      
    hDict['hJetPt'].Fill(jet.pt())

  # Find di-jets and fill histograms
  if len(jets_accepted) > 1:

    pT1 = jets_accepted[0].pt()
    pT2 = jets_accepted[1].pt()
    phi1 = jets_accepted[0].phi()
    phi2 = jets_accepted[1].phi()
  
    AJ = (pT1 - pT2) / (pT1 + pT2)
    deltaPhi = abs(phi1 - phi2)
    
    hDict['hAJ'].Fill(AJ, deltaPhi)

#---------------------------------------------------------------
def fillSoftDropHistograms(hDict, jets_sd):
  
  for jet in jets_sd:
    
    sd_info = rt.get_SD_jet_info(jet)
    zg = sd_info.z
    Rg = sd_info.dR
    
    hDict['hZg'].Fill(zg)
    hDict['hRg'].Fill(Rg)

#---------------------------------------------------------------
def plotHistograms(hDict, outputDir, fileFormat):

  outputFilename = "{}hJetPt{}".format(outputDir, fileFormat)
  plotHist(hDict['hJetPt'], outputFilename, 'width', setLogy=True)

  outputFilename = "{}hAJ{}".format(outputDir, fileFormat)
  plotHist(hDict['hAJ'], outputFilename, 'colz')

  outputFilename = "{}hZg{}".format(outputDir, fileFormat)
  plotHist(hDict['hZg'], outputFilename)

  outputFilename = "{}hRg{}".format(outputDir, fileFormat)
  plotHist(hDict['hRg'], outputFilename)

#---------------------------------------------------------------
def plotHist(h, outputFilename, drawOptions = "", setLogy = False, setLogz = False):
  
  h.SetLineColor(1)
  h.SetLineWidth(1)
  h.SetLineStyle(1)
  
  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd()
  ROOT.gPad.SetLeftMargin(0.15)
  if setLogy:
    c.SetLogy()
  if setLogz:
    c.SetLogz()
  ROOT.gPad.SetLeftMargin(0.15)

  h.Draw(drawOptions)
  c.SaveAs(outputFilename)
  c.Close()

#----------------------------------------------------------------------
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Plot analysis histograms')
  parser.add_argument('-f', '--inputFileData', action='store',
                      type=str, metavar='inputFileData',
                      default='AnalysisResults.root',
                      help='Path of ROOT file containing data TTrees')
  parser.add_argument('-r', '--inputFileMC', action='store',
                      type=str, metavar='inputFileMC',
                      default='AnalysisResults.root',
                      help='Path of ROOT file containing MC TTrees')
  parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./myTestFigures',
                      help='Output directory for QA plots to be written to')
  parser.add_argument('-i', '--imageFormat', action='store',
                      type=str, metavar='imageFormat',
                      default='.pdf',
                      help='Image format to save plots in, e.g. \'.pdf\' or \'.png\'')
  
  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('inputFileData: \'{0}\''.format(args.inputFileData))
  print('inputFileMC: \'{0}\''.format(args.inputFileMC))
  print('ouputDir: \'{0}\"'.format(args.outputDir))
  print('imageFormat: \'{0}\''.format(args.imageFormat))
  print('----------------------------------------------------------------')
  
  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFileData):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFileData))
    sys.exit(0)
  if not os.path.exists(args.inputFileMC):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFileMC))
    sys.exit(0)

runInclusiveJetAnalysis(inputFileData = args.inputFileData, inputFileMC = args.inputFileMC, outputDir = args.outputDir, fileFormat = args.imageFormat)
