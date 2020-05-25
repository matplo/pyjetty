#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of track information
  and do jet-finding, and save basic histograms.
  
  Based on code by James Mulligan (james.mulligan@berkeley.edu)
  Modified by Ezra Lesser (elesser@berkeley.edu)
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
from array import *

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io, process_utils, process_base
from pyjetty.alice_analysis.process.user.ang_pp.helpers import lambda_beta_kappa, pT_bin

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class process_ang_data(process_base.ProcessBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
    super(process_ang_data, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
    self.initialize_config()

  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def process_ang_data(self):
    
    start_time = time.time()

    # Use IO helper class to convert ROOT TTree into a SeriesGroupBy object
    # of fastjet particles per event
    print('--- {} seconds ---'.format(time.time() - start_time))
    io = process_io.ProcessIO(input_file=self.input_file, track_tree_name='tree_Particle')
    self.df_fjparticles = io.load_data()
    self.nEvents = len(self.df_fjparticles.index)
    self.nTracks = len(io.track_df.index)
    print('--- {} seconds ---'.format(time.time() - start_time))

    # Initialize configuration and histograms
    self.initialize_config()
    self.initializeHistograms()
    print(self)

    # Find jets and fill histograms
    print('Find jets...')
    self.analyzeEvents()

    # Plot histograms
    print('Save histograms...')
    process_base.ProcessBase.save_output_objects(self)

    print('--- {} seconds ---'.format(time.time() - start_time))

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
    
    # Call base class initialization
    process_base.ProcessBase.initialize_config(self)
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)

    # Set configuration for analysis
    self.jetR_list = config["jetR"]
    self.beta_list = config["betas"]

    self.n_pt_bins = config["n_pt_bins"]
    self.pt_limits = config["pt_limits"]
    self.pTbins = np.arange(self.pt_limits[0], self.pt_limits[1] + 1, 
                            (self.pt_limits[1] - self.pt_limits[0]) / self.n_pt_bins)
    self.n_lambda_bins = config["n_lambda_bins"]
    self.lambda_limits = config["lambda_limits"]
    self.n_rap_bins = config["n_rap_bins"]
    self.rap_limits = config["rap_limits"]

    # SoftDrop configuration
    self.sd_zcut = config["sd_zcut"]
    self.sd_beta = config["sd_beta"]

    # Configs for each jetR / beta
    self.config_dict = config["ang"]
  
  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initializeHistograms(self):
    
    self.hNevents = ROOT.TH1F('hNevents', 'hNevents', 2, -0.5, 1.5)
    self.hNevents.Fill(1, self.nEvents)
    
    for jetR in self.jetR_list:

      name = 'hJetPt_R%s' % jetR
      h = ROOT.TH1F(name, name, 300, 0, 300)
      h.GetXaxis().SetTitle('p_{T,ch jet}')
      h.GetYaxis().SetTitle('dN/dp_{T}')
      setattr(self, name, h)

      '''
      name = 'hM_JetPt_R%s' % jetR
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 50.)
      h.GetXaxis().SetTitle('p_{T,ch jet}')
      h.GetYaxis().SetTitle('m_{ch jet}')
      setattr(self, name, h)
      '''
      
      for beta in self.beta_list:

        '''
        for i, pTmin in list(enumerate(self.pTbins))[0:-1]:

          # Individual angularity plots, \lambda_{\beta}^{\kappa}
          pTmax = self.pTbins[i+1]
          name = "hLambda_pT%i-%i_%s" % (pTmin, pTmax, config)
          h = ROOT.TH1F(name, name, self.n_lambda_bins, 0, 1.0)
          h.GetXaxis().SetTitle('#lambda_{%s}' % beta)
          h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{%s}}' % beta)
          setattr(self, name, h)

          # Angularities with soft drop
          name = "hLambda_pT%i-%i_%s_SD" % (pTmin, pTmax, config)
          h = ROOT.TH1F(name, name, self.n_lambda_bins, 0, 1.0)
          h.GetXaxis().SetTitle('#lambda_{%s}' % beta)
          h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{%s}}' % beta)
          setattr(self, name, h)
        '''

        # Lambda vs pT plots with fine binning
        name = "h_ang_JetPt_R%s_%s" % (jetR, beta)
        h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1], 
                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
        h.GetXaxis().SetTitle('p_{T,ch jet}')
        h.GetYaxis().SetTitle('#lambda_{%s}' % beta)
        setattr(self, name, h)

        # Lambda vs pT plots with fine binning -- with soft drop
        name = "h_ang_JetPt_R%s_%s_SD" % (jetR, beta)
        h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1], 
                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
        h.GetXaxis().SetTitle('p_{T,ch jet}')
        h.GetYaxis().SetTitle('#lambda_{%s}' % beta)
        setattr(self, name, h)

        # Lambda vs rapidity plots with fine binning
        name = "h_ang_JetRap_R%s_%s" % (jetR, beta)
        h = ROOT.TH2F(name, name, self.n_rap_bins, self.rap_limits[0], self.rap_limits[1], 
                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
        h.GetXaxis().SetTitle('#eta_{ch jet}')
        h.GetYaxis().SetTitle('#lambda_{%s}' % beta)
        setattr(self, name, h)

        # Lambda vs pT plots with fine binning -- with soft drop
        name = "h_ang_JetRap_R%s_%s_SD" % (jetR, beta)
        h = ROOT.TH2F(name, name, self.n_rap_bins, self.rap_limits[0], self.rap_limits[1], 
                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
        h.GetXaxis().SetTitle('#eta_{ch jet}')
        h.GetYaxis().SetTitle('#lambda_{%s}' % beta)
        setattr(self, name, h)

  #---------------------------------------------------------------
  # Main function to loop through and analyze events
  #---------------------------------------------------------------
  def analyzeEvents(self):
    
    fj.ClusterSequence.print_banner()
    print()
    
    for jetR in self.jetR_list:
      
      # Set jet definition and a jet selector
      jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
      jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
      print('jet definition is:', jet_def)
      print('jet selector is:', jet_selector,'\n')
        
      # Define SoftDrop settings
      sd = fjcontrib.SoftDrop(self.sd_beta, self.sd_zcut, jetR)
      print('SoftDrop groomer is: {}'.format(sd.description()));
      
      for beta in self.beta_list:

        # Use list comprehension to do jet-finding and fill histograms
        result = [self.analyzeJets(fj_particles, jet_def, jet_selector, beta, sd)
                  for fj_particles in self.df_fjparticles]

  #---------------------------------------------------------------
  # Analyze jets of a given event.
  # fj_particles is the list of fastjet pseudojets for a single fixed event.
  #---------------------------------------------------------------
  def analyzeJets(self, fj_particles, jet_def, jet_selector, beta, sd):
    
    # Do jet finding
    cs = fj.ClusterSequence(fj_particles, jet_def)
    jets = fj.sorted_by_pt(cs.inclusive_jets())
    jets_selected = jet_selector(jets)

    # Loop through jets
    jetR = jet_def.R()
    for jet in jets_selected:
      
      if self.debug_level > 1:
        print('jet:')
        print(jet)
      
      # Check additional acceptance criteria
      if not self.utils.is_det_jet_accepted(jet):
        continue
      
      # Fill histograms
      self.fillJetHistograms(jet, jetR, beta, sd)


  #---------------------------------------------------------------
  # Fill jet histograms
  #---------------------------------------------------------------
  def fillJetHistograms(self, jet, jetR, beta, sd):

    jet_pt = jet.pt()
    config = ("R%s_B%s" % (jetR, beta)).replace('.', '')

    # jet with softdrop
    jet_sd = sd.result(jet)

    # Calculate observable for this jet
    # just use kappa = 1 for now
    l = lambda_beta_kappa(jet, jetR, beta, 1)
    lsd = lambda_beta_kappa(jet_sd, jetR, beta, 1)

    # Fill plots
    getattr(self, ("h_ang_JetPt_R%s_%s" % (jetR, beta))).Fill(jet_pt, l)
    getattr(self, ("h_ang_JetPt_R%s_%s_SD" % (jetR, beta))).Fill(jet_sd.pt(), lsd)
    getattr(self, ("h_ang_JetRap_R%s_%s" % (jetR, beta))).Fill(jet.rap(), l)
    getattr(self, ("h_ang_JetRap_R%s_%s_SD" % (jetR, beta))).Fill(jet_sd.rap(), lsd)

    '''
    # just use kappa = 1 for now
    l = lambda_beta_kappa(jet, jetR, beta, 1)
    (pTmin, pTmax) = pT_bin(jet_pt, self.pTbins)
    if pTmin > -1e-3:  # pTmin will be -1 if not a valid bin
      getattr(self, ("hLambda_pT%i-%i_%s" % (pTmin, pTmax, config))).Fill(l)
      getattr(self, ("hLambda_JetpT_%s" % config)).Fill(jet_pt, l)

    # jet with softdrop
    jet_sd = sd.result(jet)
    lsd = lambda_beta_kappa(jet_sd, jetR, beta, 1)
    (pTmin, pTmax) = pT_bin(jet_sd.pt(), self.pTbins)
    if pTmin > -1e-3:  # pTmin will be -1 if not a valid bin
      getattr(self, ("hLambda_pT%i-%i_%s_SD" % (pTmin, pTmax, config))).Fill(lsd)
      getattr(self, ("hLambda_JetpT_%s_SD" % config)).Fill(jet_sd.pt(), lsd)
    '''

    getattr(self, 'hJetPt_R%s' % str(jetR)).Fill(jet_pt)
    #getattr(self, 'hM_JetPt_R%s' % str(jetR)).Fill(jet_pt, jet.m())


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
                      default='config/angularity.yaml',
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

  analysis = process_ang_data(input_file=args.inputFile, config_file=args.configFile, 
                              output_dir=args.outputDir)
  analysis.process_ang_data()
