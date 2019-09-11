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
import yaml

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
def process_rg_data(inputFile, configFile, outputDir):
  
  start_time = time.time()
  
  # Read config file
  with open(configFile, 'r') as stream:
    config = yaml.safe_load(stream)
  
  jetR_list = config['jetR']
  beta_list = config['beta']
  
  # Create output dir
  if not outputDir.endswith("/"):
    outputDir = outputDir + "/"
  if not os.path.exists(outputDir):
    os.makedirs(outputDir)

  # Initialize histogram dictionary
  hDict = initializeHistograms(jetR_list, beta_list)

  # Convert ROOT TTree to pandas dataframe with one row per jet constituent:
  #     run_number, ev_id, ParticlePt, ParticleEta, ParticlePhi
  print('--- {} seconds ---'.format(time.time() - start_time))
  print('Convert ROOT trees to pandas dataframes...')
  track_df = analysis_utils.load_dataframe(inputFile, 'tree_Particle')

  # Transform the track dataframe into a SeriesGroupBy object of fastjet particles per event
  print('--- {} seconds ---'.format(time.time() - start_time))
  df_fjparticles = analysis_utils.group_fjparticles(track_df)
  
  if debugLevel > 0:
    print(df_fjparticles.dtypes)
    print(df_fjparticles)
  print('--- {} seconds ---'.format(time.time() - start_time))

  # Print number of events
  nEvents = len(df_fjparticles.index)
  hDict['hNevents'].Fill(1, nEvents)
  print('Number of events: {}'.format(nEvents))
  nTracks = len(track_df.index)
  print('Number of tracks: {}'.format(nTracks))

  # Find jets and fill histograms
  print('Find jets...')
  analyzeEvents(df_fjparticles, hDict, jetR_list, beta_list, outputDir)

  # Plot histograms
  print('Save histograms...')
  saveHistograms(hDict, outputDir)

  print('--- {} seconds ---'.format(time.time() - start_time))

#---------------------------------------------------------------
def initializeHistograms(jetR_list, beta_list):
  
  hDict = {}
  
  name = 'hNevents'
  hNevents = ROOT.TH1F(name, name, 2, -0.5, 1.5)
  hDict[name] = hNevents

  for jetR in jetR_list:
  
    name = 'hJetPt_R{}'.format(jetR)
    hJetPt = ROOT.TH1F(name, name, 300, 0, 300)
    hJetPt.GetXaxis().SetTitle('p_{T,jet}')
    hJetPt.GetYaxis().SetTitle('dN/dp_{T}')
    hDict[name] = hJetPt
    
    name = 'hM_JetPt_R{}'.format(jetR)
    hM_JetPt = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 50.)
    hM_JetPt.GetXaxis().SetTitle('p_{T,jet}')
    hM_JetPt.GetYaxis().SetTitle('m_{jet}')
    hDict[name] = hM_JetPt
    
    for beta in beta_list:

      name = 'hThetaG_JetPt_R{}_B{}'.format(jetR, beta)
      hThetaG_JetPt = ROOT.TH2F(name, name, 300, 0, 300, 150, 0, 1.5)
      hThetaG_JetPt.GetXaxis().SetTitle('p_{T,jet}')
      hThetaG_JetPt.GetYaxis().SetTitle('#theta_{g}')
      hDict[name] = hThetaG_JetPt

      name = 'hZg_JetPt_R{}_B{}'.format(jetR, beta)
      hZg_JetPt = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 1.)
      hZg_JetPt.GetXaxis().SetTitle('p_{T,jet}')
      hZg_JetPt.GetYaxis().SetTitle('z_{g}')
      hDict[name] = hZg_JetPt

      name = 'hMg_JetPt_R{}_B{}'.format(jetR, beta)
      hMg_JetPt = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 50.)
      hMg_JetPt.GetXaxis().SetTitle('p_{T,jet}')
      hMg_JetPt.GetYaxis().SetTitle('m_{g,jet}')
      hDict[name] = hMg_JetPt

  return hDict

#---------------------------------------------------------------
def analyzeEvents(df_fjparticles, hDict, jetR_list, beta_list, outputDir):
  
  fj.ClusterSequence.print_banner()
  print()
  
  for jetR in jetR_list:
    
    for beta in beta_list:
    
      # Set jet definition and a jet selector
      jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
      jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
      print('jet definition is:', jet_def)
      print('jet selector is:', jet_selector,'\n')
      
      # Define SoftDrop settings
      zcut = 0.1
      sd = fjcontrib.SoftDrop(beta, zcut, jetR)
      print('SoftDrop groomer is: {}'.format(sd.description()));

      # Use list comprehension to do jet-finding and fill histograms
      result = [analyzeJets(fj_particles, jet_def, jet_selector, sd, beta, hDict) for fj_particles in df_fjparticles]

#---------------------------------------------------------------
def analyzeJets(fj_particles, jet_def, jet_selector, sd, beta, hDict):
  
  # Do jet finding
  cs = fj.ClusterSequence(fj_particles, jet_def)
  jets = fj.sorted_by_pt(cs.inclusive_jets())
  jets_accepted = jet_selector(jets)

  # Loop through jets
  jetR = jet_def.R()
  for jet in jets_accepted:
    
    if debugLevel > 1:
      print('jet:')
      print(jet)
    
    # Fill histograms
    fillJetHistograms(hDict, jet, jetR)
    
    # Perform SoftDrop grooming and fill histograms
    jet_sd = sd.result(jet)
    fillSoftDropHistograms(hDict, jet_sd, jet, jetR, beta)

#---------------------------------------------------------------
def fillJetHistograms(hDict, jet, jetR):
  
  jet_pt = jet.pt()

  hDict['hJetPt_R{}'.format(jetR)].Fill(jet_pt)
  hDict['hM_JetPt_R{}'.format(jetR)].Fill(jet_pt, jet.m())

#---------------------------------------------------------------
def fillSoftDropHistograms(hDict, jet_sd, jet, jetR, beta):
  
  jet_pt_ungroomed = jet.pt()
  
  sd_info = fjcontrib.get_SD_jet_info(jet_sd)
  theta_g = sd_info.dR / jetR
  zg = sd_info.z
  mg = jet_sd.m()

  hDict['hThetaG_JetPt_R{}_B{}'.format(jetR, beta)].Fill(jet_pt_ungroomed, theta_g)
  hDict['hZg_JetPt_R{}_B{}'.format(jetR, beta)].Fill(jet_pt_ungroomed, zg)
  hDict['hMg_JetPt_R{}_B{}'.format(jetR, beta)].Fill(jet_pt_ungroomed, mg)

#---------------------------------------------------------------
def saveHistograms(hDict, outputDir):
  
  outputfilename = os.path.join(outputDir, 'AnalysisResults.root')
  fout = ROOT.TFile(outputfilename, 'recreate')
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
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='analysis_config.yaml',
                      help="Path of config file for jetscape analysis")
  parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./TestOutput',
                      help='Output directory for QA plots to be written to')
  
  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('inputFile: \'{0}\''.format(args.inputFile))
  print('configFile: \'{0}\''.format(args.configFile))
  print('ouputDir: \'{0}\"'.format(args.outputDir))
  print('----------------------------------------------------------------')
  
  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
    sys.exit(0)
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  process_rg_data(inputFile = args.inputFile, configFile = args.configFile, outputDir = args.outputDir)
