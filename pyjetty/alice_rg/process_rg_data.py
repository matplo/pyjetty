#!/usr/bin/env python3

"""
  Analysis script to read a ROOT TTree of track information
  and do jet-finding, and save basic histograms.
  
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

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext

# Analysis utilities
import analysis_utils

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

# Set debug level (0 = no debug info, 1 = some debug info, 2 = all debug info)
debugLevel = 0

#---------------------------------------------------------------
def process_rg_data(inputFile, outputDir):
  
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
  track_df = analysis_utils.load_dataframe(inputFile)

  # Transform the track dataframe into a SeriesGroupBy object of fastjet particles per event
  print('--- {} seconds ---'.format(time.time() - start_time))
  print('Transform the track dataframe into a series object of fastjet particles per event...')

  # (i) Group the track dataframe by event
  #     track_df_grouped is a DataFrameGroupBy object with one track dataframe per event
  track_df_grouped = track_df.groupby(['run_number','ev_id'])

  # (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet particles
  df_fjparticles = track_df_grouped.apply(analysis_utils.get_fjparticles)
  
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
  analyzeEvents(df_fjparticles, hDict, outputDir)

  # Plot histograms
  print('Plot histograms...')
  saveHistograms(hDict, outputDir)

  print('--- {} seconds ---'.format(time.time() - start_time))

#---------------------------------------------------------------
def initializeHistograms():
  
  hDict = {}
  
  hJetPt = ROOT.TH1F('hJetPt', 'hJetPt', 300, 0, 300)
  hJetPt.GetXaxis().SetTitle('p_{T,jet}')
  hJetPt.GetYaxis().SetTitle('dN/dp_{T}')
  hDict['hJetPt'] = hJetPt
  
  hThetaG_JetPt = ROOT.TH2F('hThetaG_JetPt', 'hThetaG_JetPt', 300, 0, 300, 150, 0, 1.5)
  hThetaG_JetPt.GetXaxis().SetTitle('p_{T,jet}')
  hThetaG_JetPt.GetYaxis().SetTitle('#theta_{g}')
  hDict['hThetaG_JetPt'] = hThetaG_JetPt

  hZg_JetPt = ROOT.TH2F('hZg_JetPt', 'hZg_JetPt', 300, 0, 300, 100, 0, 1.)
  hZg_JetPt.GetXaxis().SetTitle('p_{T,jet}')
  hZg_JetPt.GetYaxis().SetTitle('z_{g}')
  hDict['hZg_JetPt'] = hZg_JetPt
  
  hM_JetPt = ROOT.TH2F('hM_JetPt', 'hM_JetPt', 300, 0, 300, 100, 0, 50.)
  hM_JetPt.GetXaxis().SetTitle('p_{T,jet}')
  hM_JetPt.GetYaxis().SetTitle('m_{jet}')
  hDict['hM_JetPt'] = hM_JetPt

  return hDict

#---------------------------------------------------------------
def analyzeEvents(df_fjparticles, hDict, outputDir):
  
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
  sd = fjcontrib.SoftDrop(beta, zcut, jetR)
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

  fillSoftDropHistograms(hDict, jets_sd, jet_def.R())

#---------------------------------------------------------------
def fillJetHistograms(hDict, jets_accepted):

  # Loop through jets, and fill histograms
  for jet in jets_accepted:
    if debugLevel > 1:
      print('jet:')
      print(jet)
      
    hDict['hJetPt'].Fill(jet.pt())
    hDict['hM_JetPt'].Fill(jet.pt(), jet.m())

#---------------------------------------------------------------
def fillSoftDropHistograms(hDict, jets_sd, jetR):
  
  for jet in jets_sd:
    
    sd_info = fjcontrib.get_SD_jet_info(jet)
    zg = sd_info.z
    theta_g = sd_info.dR / jetR
    
    hDict['hZg_JetPt'].Fill(jet.pt(), zg)
    hDict['hThetaG_JetPt'].Fill(jet.pt(), theta_g)

#---------------------------------------------------------------
def saveHistograms(hDict, outputDir):
  
  fout = ROOT.TFile('AnalysisResults.root', 'recreate')
  fout.cd()
  for key, val in hDict.items():
    val.Write()
  fout.Close()

#----------------------------------------------------------------------
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Plot analysis histograms')
  parser.add_argument('-f', '--inputFile', action='store',
                      type=str, metavar='inputFile',
                      default='AnalysisResults.root',
                      help='Path of ROOT file containing TTrees')
  parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./myTestFigures',
                      help='Output directory for QA plots to be written to')
  
  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('inputFile: \'{0}\''.format(args.inputFile))
  print('ouputDir: \'{0}\"'.format(args.outputDir))
  print('----------------------------------------------------------------')
  
  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
    sys.exit(0)

  process_rg_data(inputFile = args.inputFile, outputDir = args.outputDir)
