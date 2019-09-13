#!/usr/bin/env python3

"""
  Analysis script to read a ROOT TTree of MC track information
  and do jet-finding, and save response histograms.
  
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
from array import *
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
  jet_matching_distance = config['jet_matching_distance']
  
  # Create output dir
  if not outputDir.endswith("/"):
    outputDir = outputDir + "/"
  if not os.path.exists(outputDir):
    os.makedirs(outputDir)

  # Initialize histogram dictionary
  hDict = initializeHistograms(jetR_list, beta_list)

  # ------------------------------------------------------------------------
  # Convert ROOT TTree to pandas dataframe with one row per jet constituent:
  #     run_number, ev_id, ParticlePt, ParticleEta, ParticlePhi
  print('--- {} seconds ---'.format(time.time() - start_time))
  print('Convert ROOT trees to pandas dataframes...')
  track_df_det = analysis_utils.load_dataframe(inputFile, 'tree_Particle')

  # Transform the track dataframe into a SeriesGroupBy object of fastjet particles per event
  print('--- {} seconds ---'.format(time.time() - start_time))
  df_fjparticles_det = analysis_utils.group_fjparticles(track_df_det)
  
  if debugLevel > 0:
    print(df_fjparticles_det.dtypes)
    print(df_fjparticles_det)
  print('--- {} seconds ---'.format(time.time() - start_time))

  # ------------------------------------------------------------------------
  # Convert ROOT TTree to pandas dataframe with one row per jet constituent:
  #     run_number, ev_id, ParticlePt, ParticleEta, ParticlePhi
  print('Convert ROOT trees to pandas dataframes...')
  track_df_truth = analysis_utils.load_dataframe(inputFile, 'tree_Particle_gen')

  # Transform the track dataframe into a SeriesGroupBy object of fastjet particles per event
  print('--- {} seconds ---'.format(time.time() - start_time))
  df_fjparticles_truth = analysis_utils.group_fjparticles(track_df_truth)
  
  if debugLevel > 0:
    print(df_fjparticles_truth.dtypes)
    print(df_fjparticles_truth)
  print('--- {} seconds ---'.format(time.time() - start_time))

  # Print number of events
  nEvents = len(df_fjparticles_truth.index)
  hDict['hNevents'].Fill(1, nEvents)
  print('Number of events: {}'.format(nEvents))
  nTracks = len(track_df_det.index)
  print('Number of det tracks: {}'.format(nTracks))
  nTracks = len(track_df_truth.index)
  print('Number of truth tracks: {}'.format(nTracks))
  # ------------------------------------------------------------------------

  # Now merge the two SeriesGroupBy to create a groupby df with [ev_id, run_number, fj_1, fj_2]
  # (Need a structure such that we can iterate event-by-event through both fj_1, fj_2 simultaneously)
  print('Merge det-level and truth-level into a single dataframe grouped by event...')
  df_fjparticles = pandas.concat([df_fjparticles_det, df_fjparticles_truth], axis=1)
  df_fjparticles.columns = ['fj_particles_det', 'fj_particles_truth']
  if debugLevel > 0:
    print(df_fjparticles.dtypes)
    print(df_fjparticles)
  print('--- {} seconds ---'.format(time.time() - start_time))

  # ------------------------------------------------------------------------

  # Find jets and fill histograms
  print('Find jets...')
  analyzeEvents(df_fjparticles, hDict, jetR_list, beta_list, jet_matching_distance, outputDir)

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
  
    name = 'hResponse_JetPt_R{}'.format(jetR)
    hResponse_JetPt = ROOT.TH2F(name, name, 300, 0, 300, 300, 0, 300)
    hResponse_JetPt.GetXaxis().SetTitle('p_{T,det}')
    hResponse_JetPt.GetYaxis().SetTitle('p_{T,truth}')
    hDict[name] = hResponse_JetPt
    
    for beta in beta_list:
      
      name = 'hResponse_JetPt_ThetaG_R{}_B{}'.format(jetR, beta)
      nbins = ([300, 300, 150, 150])
      xmin = ([0, 0, 0, 0])
      xmax = ([300, 300, 1.5, 1.5])
      nbins_array = array('i', nbins)
      xmin_array = array('d', xmin)
      xmax_array = array('d', xmax)
      hResponse_JetPt_ThetaG = ROOT.THnSparseF(name, name, 4, nbins_array, xmin_array, xmax_array)
      hResponse_JetPt_ThetaG.GetAxis(0).SetTitle('p_{T,det}')
      hResponse_JetPt_ThetaG.GetAxis(1).SetTitle('p_{T,truth}')
      hResponse_JetPt_ThetaG.GetAxis(2).SetTitle('#theta_{g,det}')
      hResponse_JetPt_ThetaG.GetAxis(3).SetTitle('#theta_{g,truth}')
      hDict[name] = hResponse_JetPt_ThetaG

  return hDict

#---------------------------------------------------------------
def analyzeEvents(df_fjparticles, hDict, jetR_list, beta_list, jet_matching_distance, outputDir):
  
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
      
      # Then can use list comprehension to iterate over the groupby and do jet-finding
      # simultaneously for fj_1 and fj_2 per event, so that I can match jets -- and fill histograms
      result = [analyzeJets(fj_particles_det, fj_particles_truth, jet_def, jet_selector, sd, beta, jet_matching_distance, hDict) for fj_particles_det, fj_particles_truth in zip(df_fjparticles['fj_particles_det'], df_fjparticles['fj_particles_truth'])]

#---------------------------------------------------------------
def analyzeJets(fj_particles_det, fj_particles_truth, jet_def, jet_selector, sd, beta, jet_matching_distance, hDict):

  # Check that the entries exist appropriately
  # (need to check how this can happen -- but it is only a tiny fraction of events)
  if type(fj_particles_det) != fj.vectorPJ or type(fj_particles_truth) != fj.vectorPJ:
    print('fj_particles type mismatch -- skipping event')
    return

  # Do jet finding
  cs_det = fj.ClusterSequence(fj_particles_det, jet_def)
  jets_det = fj.sorted_by_pt(cs_det.inclusive_jets())
  jets_det_accepted = jet_selector(jets_det)
  
  cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)
  jets_truth = fj.sorted_by_pt(cs_truth.inclusive_jets())
  jets_truth_accepted = jet_selector(jets_truth)

  # Set number of jet matches for each jet in user_index (to ensure unique matches)
  jetR = jet_def.R()
  setNJetMatches(jets_det_accepted, jets_truth_accepted, jet_matching_distance, jetR)
  
  # Loop through jets and fill response if both det and truth jets are unique match
  for jet_det in jets_det_accepted:
    for jet_truth in jets_truth_accepted:
      
      if debugLevel > 0:
        print('deltaR: {}'.format(jet_det.delta_R(jet_truth)))
        print('jet_det matches: {}'.format(jet_det.user_index()))
        print('jet_truth matches: {}'.format(jet_truth.user_index()))
      
      # Check that jets match geometrically
      if jet_det.delta_R(jet_truth) < jet_matching_distance*jetR:

        #Check that match is unique
        if jet_det.user_index() == 1 and jet_truth.user_index() == 1:

          fillResponseHistograms(jet_det, jet_truth, sd, hDict, jetR, beta)

#---------------------------------------------------------------
# Loop through jets and store number of matching candidates in user_index
# (In principle could also store matching candidate in user_info)
#---------------------------------------------------------------
def setNJetMatches(jets_det_accepted, jets_truth_accepted, jet_matching_distance, jetR):

  # Reset user_index to 0
  for jet_det in jets_det_accepted:
    jet_det.set_user_index(0)
  for jet_truth in jets_truth_accepted:
    jet_truth.set_user_index(0)
   
  # Loop through jets and store number of matching candidates in user_index
  for jet_det in jets_det_accepted:
    for jet_truth in jets_truth_accepted:
      
      if jet_det.delta_R(jet_truth) < jet_matching_distance*jetR:
        
        jet_det.set_user_index(jet_det.user_index() + 1)
        jet_truth.set_user_index(jet_truth.user_index() + 1)

#---------------------------------------------------------------
def fillResponseHistograms(jet_det, jet_truth, sd, hDict, jetR, beta):
  
  jet_pt_det_ungroomed = jet_det.pt()
  jet_pt_truth_ungroomed = jet_truth.pt()
  
  theta_g_det = theta_g(jet_det, sd, jetR)
  theta_g_truth = theta_g(jet_truth, sd, jetR)
  
  hDict['hResponse_JetPt_R{}'.format(jetR)].Fill(jet_pt_det_ungroomed, jet_pt_truth_ungroomed)
  
  x = ([jet_pt_det_ungroomed, jet_pt_truth_ungroomed, theta_g_det, theta_g_truth])
  x_array = array('d', x)
  hDict['hResponse_JetPt_ThetaG_R{}_B{}'.format(jetR, beta)].Fill(x_array)

#---------------------------------------------------------------
def theta_g(jet, sd, jetR):
  
  jet_sd = sd.result(jet)
  sd_info = fjcontrib.get_SD_jet_info(jet_sd)
  theta_g = sd_info.dR / jetR
  return theta_g

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
