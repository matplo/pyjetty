#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of MC track information
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
from pyjetty.alice_analysis.process_base import analysis_io
from pyjetty.alice_analysis.process_base import analysis_utils
from pyjetty.alice_analysis.process_base import analysis_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class analysis_rg_mc(analysis_base.analysis_base):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
    super(analysis_rg_mc, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
    self.initialize_config()
  
  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def process_rg_mc(self):
    
    start_time = time.time()
    
    # Initialize configuration
    self.initialize_config()
    
    # ------------------------------------------------------------------------
    
    # Use IO helper class to convert detector-level ROOT TTree into
    # a SeriesGroupBy object of fastjet particles per event
    print('--- {} seconds ---'.format(time.time() - start_time))
    io_det = analysis_io.analysis_io(input_file=self.input_file, track_tree_name='tree_Particle')
    df_fjparticles_det = io_det.load_data(self.reject_tracks_fraction)
    self.nEvents_det = len(df_fjparticles_det.index)
    self.nTracks_det = len(io_det.track_df.index)
    print('--- {} seconds ---'.format(time.time() - start_time))
    
    # ------------------------------------------------------------------------

    # Use IO helper class to convert truth-level ROOT TTree into
    # a SeriesGroupBy object of fastjet particles per event
    io_truth = analysis_io.analysis_io(input_file=self.input_file, track_tree_name='tree_Particle_gen')
    df_fjparticles_truth = io_truth.load_data()
    self.nEvents_truth = len(df_fjparticles_truth.index)
    self.nTracks_truth = len(io_truth.track_df.index)
    print('--- {} seconds ---'.format(time.time() - start_time))
    
    # ------------------------------------------------------------------------

    # Now merge the two SeriesGroupBy to create a groupby df with [ev_id, run_number, fj_1, fj_2]
    # (Need a structure such that we can iterate event-by-event through both fj_1, fj_2 simultaneously)
    print('Merge det-level and truth-level into a single dataframe grouped by event...')
    self.df_fjparticles = pandas.concat([df_fjparticles_det, df_fjparticles_truth], axis=1)
    self.df_fjparticles.columns = ['fj_particles_det', 'fj_particles_truth']
    print('--- {} seconds ---'.format(time.time() - start_time))

    # ------------------------------------------------------------------------

    # Initialize histograms
    self.initializeHistograms()
    print(self)
  
    # Find jets and fill histograms
    print('Find jets...')
    self.analyzeEvents()
    
    # Plot histograms
    print('Save histograms...')
    analysis_base.analysis_base.saveHistograms(self)
    
    print('--- {} seconds ---'.format(time.time() - start_time))
  
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
    
    # Call base class initialization
    analysis_base.analysis_base.initialize_config(self)
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    
    self.jet_matching_distance = config['jet_matching_distance']
    self.reject_tracks_fraction = config['reject_tracks_fraction']
    
    # config['beta'] is a dictionary of dictionaries, where each dict is for a value of beta
    beta_dict = config['beta']
    
    # Retrieve list of beta values
    self.beta_list = list(beta_dict.keys())
    
    # Retrieve histogram binnings for each beta value
    for beta in self.beta_list:
      
      binning_dict = beta_dict[beta]
      pt_bins_det = (binning_dict['pt_bins_det'])
      rg_bins_det = (binning_dict['rg_bins_det'])
      pt_bins_truth = (binning_dict['pt_bins_truth'])
      rg_bins_truth = (binning_dict['rg_bins_truth'])
      
      n_pt_bins_det = len(pt_bins_det) - 1
      setattr(self, 'n_pt_bins_det_B{}'.format(beta), n_pt_bins_det)
      
      n_rg_bins_det = len(rg_bins_det) - 1
      setattr(self, 'n_rg_bins_det_B{}'.format(beta), n_rg_bins_det)

      n_pt_bins_truth = len(pt_bins_truth) - 1
      setattr(self, 'n_pt_bins_truth_B{}'.format(beta), n_pt_bins_truth)

      n_rg_bins_truth = len(rg_bins_truth) - 1
      setattr(self, 'n_rg_bins_truth_B{}'.format(beta), n_rg_bins_truth)
      
      det_pt_bin_array = array('d',pt_bins_det)
      setattr(self, 'det_pt_bin_array_B{}'.format(beta), det_pt_bin_array)
      
      det_rg_bin_array = array('d',rg_bins_det)
      setattr(self, 'det_rg_bin_array_B{}'.format(beta), det_rg_bin_array)
      
      truth_pt_bin_array = array('d',pt_bins_truth)
      setattr(self, 'truth_pt_bin_array_B{}'.format(beta), truth_pt_bin_array)

      truth_rg_bin_array = array('d',rg_bins_truth)
      setattr(self, 'truth_rg_bin_array_B{}'.format(beta), truth_rg_bin_array)
  
  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initializeHistograms(self):
    
    name = 'hNevents'
    h = ROOT.TH1F(name, name, 2, -0.5, 1.5)
    h.Fill(1, self.nEvents_det)
    setattr(self, name, h)

    for jetR in self.jetR_list:
      
      name = 'hJetPt_Truth_R{}'.format(jetR)
      h = ROOT.TH1F(name, name, 300, 0, 300)
      setattr(self, name, h)

      name = 'hJES_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
      setattr(self, name, h)
      
      name = 'hDeltaR_All_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
      setattr(self, name, h)

      name = 'hResponse_JetPt_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 300, 0, 300)
      h.GetXaxis().SetTitle('p_{T,det}')
      h.GetYaxis().SetTitle('p_{T,truth}')
      setattr(self, name, h)

      for beta in self.beta_list:
        
        # Retrieve histogram binnings
        n_pt_bins_det = getattr(self, 'n_pt_bins_det_B{}'.format(beta))
        det_pt_bin_array = getattr(self, 'det_pt_bin_array_B{}'.format(beta))
        n_rg_bins_det = getattr(self, 'n_rg_bins_det_B{}'.format(beta))
        det_rg_bin_array = getattr(self, 'det_rg_bin_array_B{}'.format(beta))
        n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_B{}'.format(beta))
        truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_B{}'.format(beta))
        n_rg_bins_truth = getattr(self, 'n_rg_bins_truth_B{}'.format(beta))
        truth_rg_bin_array = getattr(self, 'truth_rg_bin_array_B{}'.format(beta))
        
        # Create THn of response
        dim = 0;
        title = []
        nbins = []
        min = []
        max = []
        bin_edges = []
        
        title.append('p_{T,det}')
        nbins.append(n_pt_bins_det)
        bin_edges.append(det_pt_bin_array)
        min.append(det_pt_bin_array[0])
        max.append(det_pt_bin_array[-1])
        dim+=1
        
        title.append('p_{T,truth}')
        nbins.append(n_pt_bins_truth)
        bin_edges.append(truth_pt_bin_array)
        min.append(truth_pt_bin_array[0])
        max.append(truth_pt_bin_array[-1])
        dim+=1
        
        title.append('#theta_{g,det}')
        nbins.append(n_rg_bins_det)
        bin_edges.append(det_rg_bin_array)
        min.append(det_rg_bin_array[0])
        max.append(det_rg_bin_array[-1])
        dim+=1
        
        title.append('#theta_{g,truth}')
        nbins.append(n_rg_bins_truth)
        bin_edges.append(truth_rg_bin_array)
        min.append(truth_rg_bin_array[0])
        max.append(truth_rg_bin_array[-1])
        dim+=1
        
        name = 'hResponse_JetPt_ThetaG_R{}_B{}'.format(jetR, beta)
        nbins = (nbins)
        xmin = (min)
        xmax = (max)
        nbins_array = array('i', nbins)
        xmin_array = array('d', xmin)
        xmax_array = array('d', xmax)
        h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
        for i in range(0, dim):
          h.GetAxis(i).SetTitle(title[i])
          h.SetBinEdges(i, bin_edges[i])
        setattr(self, name, h)

  #---------------------------------------------------------------
  # Main function to loop through and analyze events
  #---------------------------------------------------------------
  def analyzeEvents(self):
    
    fj.ClusterSequence.print_banner()
    print()
    
    for jetR in self.jetR_list:
      
      for beta in self.beta_list:
      
        # Set jet definition and a jet selector
        jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
        jet_selector_det = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
        jet_selector_truth_matched = fj.SelectorPtMin(5.0)
        print('jet definition is:', jet_def)
        print('jet selector for det-level is:', jet_selector_det,'\n')
        print('jet selector for truth-level matches is:', jet_selector_truth_matched,'\n')
        
        # Define SoftDrop settings
        zcut = 0.1
        sd = fjcontrib.SoftDrop(beta, zcut, jetR)
        print('SoftDrop groomer is: {}'.format(sd.description()));
        
        # Then can use list comprehension to iterate over the groupby and do jet-finding
        # simultaneously for fj_1 and fj_2 per event, so that I can match jets -- and fill histograms
        result = [self.analyzeJets(fj_particles_det, fj_particles_truth, jet_def, jet_selector_det, jet_selector_truth_matched, sd, beta) for fj_particles_det, fj_particles_truth in zip(self.df_fjparticles['fj_particles_det'], self.df_fjparticles['fj_particles_truth'])]

  #---------------------------------------------------------------
  # Analyze jets of a given event.
  # fj_particles is the list of fastjet pseudojets for a single fixed event.
  #---------------------------------------------------------------
  def analyzeJets(self, fj_particles_det, fj_particles_truth, jet_def, jet_selector_det, jet_selector_truth_matched, sd, beta):

    # Check that the entries exist appropriately
    # (need to check how this can happen -- but it is only a tiny fraction of events)
    if type(fj_particles_det) != fj.vectorPJ or type(fj_particles_truth) != fj.vectorPJ:
      print('fj_particles type mismatch -- skipping event')
      return

    # Do jet finding
    cs_det = fj.ClusterSequence(fj_particles_det, jet_def)
    jets_det = fj.sorted_by_pt(cs_det.inclusive_jets())
    jets_det_selected = jet_selector_det(jets_det)
    
    cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)
    jets_truth = fj.sorted_by_pt(cs_truth.inclusive_jets())
    jets_truth_selected = jet_selector_det(jets_truth)
    jets_truth_selected_matched = jet_selector_truth_matched(jets_truth)

    jetR = jet_def.R()

    # Fill truth-level jet histograms (before matching)
    for jet_truth in jets_truth_selected:
      self.fillJetHistograms(jet_truth, jetR)

    # Set number of jet matches for each jet in user_index (to ensure unique matches)
    self.setNJetMatches(jets_det_selected, jets_truth_selected_matched, jetR)
    
    # Loop through jets and fill response if both det and truth jets are unique match
    for jet_det in jets_det_selected:
      for jet_truth in jets_truth_selected_matched:
        
        # Check additional acceptance criteria
        if not self.utils.is_det_jet_accepted(jet_det):
          continue
        
        if self.debug_level > 0:
          print('deltaR: {}'.format(jet_det.delta_R(jet_truth)))
          print('jet_det matches: {}'.format(jet_det.user_index()))
          print('jet_truth matches: {}'.format(jet_truth.user_index()))
        
        # Check that jets match geometrically
        deltaR = jet_det.delta_R(jet_truth)
        getattr(self, 'hDeltaR_All_R{}'.format(jetR)).Fill(jet_det.pt(), deltaR)
        
        if deltaR < self.jet_matching_distance*jetR:

          # Check that match is unique
          if jet_det.user_index() == 1 and jet_truth.user_index() == 1:

            self.fillResponseHistograms(jet_det, jet_truth, sd, jetR, beta)

  #---------------------------------------------------------------
  # Loop through jets and store number of matching candidates in user_index
  # (In principle could also store matching candidate in user_info)
  #---------------------------------------------------------------
  def setNJetMatches(self, jets_det_selected, jets_truth_selected, jetR):

    # Reset user_index to 0
    for jet_det in jets_det_selected:
      jet_det.set_user_index(0)
    for jet_truth in jets_truth_selected:
      jet_truth.set_user_index(0)
    
    # Loop through jets and store number of matching candidates in user_index
    for jet_det in jets_det_selected:
      for jet_truth in jets_truth_selected:
        
        if jet_det.delta_R(jet_truth) < self.jet_matching_distance*jetR:
          
          jet_det.set_user_index(jet_det.user_index() + 1)
          jet_truth.set_user_index(jet_truth.user_index() + 1)

  #---------------------------------------------------------------
  # Fill jet histograms
  #---------------------------------------------------------------
  def fillJetHistograms(self, jet, jetR):

    getattr(self, 'hJetPt_Truth_R{}'.format(jetR)).Fill(jet.pt())

  #---------------------------------------------------------------
  # Fill response histograms
  #---------------------------------------------------------------
  def fillResponseHistograms(self, jet_det, jet_truth, sd, jetR, beta):
    
    jet_pt_det_ungroomed = jet_det.pt()
    jet_pt_truth_ungroomed = jet_truth.pt()
    
    theta_g_det = self.theta_g(jet_det, sd, jetR)
    theta_g_truth = self.theta_g(jet_truth, sd, jetR)
    
    JES = (jet_pt_det_ungroomed - jet_pt_truth_ungroomed) / jet_pt_truth_ungroomed
    getattr(self, 'hJES_R{}'.format(jetR)).Fill(jet_pt_truth_ungroomed, JES)
    
    getattr(self, 'hResponse_JetPt_R{}'.format(jetR)).Fill(jet_pt_det_ungroomed, jet_pt_truth_ungroomed)
    
    x = ([jet_pt_det_ungroomed, jet_pt_truth_ungroomed, theta_g_det, theta_g_truth])
    x_array = array('d', x)
    getattr(self, 'hResponse_JetPt_ThetaG_R{}_B{}'.format(jetR, beta)).Fill(x_array)

  #---------------------------------------------------------------
  # Compute theta_g
  #---------------------------------------------------------------
  def theta_g(self, jet, sd, jetR):
    
    jet_sd = sd.result(jet)
    sd_info = fjcontrib.get_SD_jet_info(jet_sd)
    theta_g = sd_info.dR / jetR
    return theta_g

##################################################################
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Plot analysis histograms')
  parser.add_argument('-f', '--inputFile', action='store',
                      type=str, metavar='inputFile',
                      default='AnalysisResults.root',
                      help='Path of ROOT file containing TTrees')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='config/analysis_config.yaml',
                      help="Path of config file for analysis")
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

  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
    sys.exit(0)
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = analysis_rg_mc(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_rg_mc()
