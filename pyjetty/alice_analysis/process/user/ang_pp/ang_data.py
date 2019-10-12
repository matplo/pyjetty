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
from pyjetty.alice_analysis.process.user.ang_pp.helper import deltaR, lambda_beta_kappa

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class process_ang_data(process_base.process_base):

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

    # Use IO helper class to convert ROOT TTree into a SeriesGroupBy object of fastjet particles per event
    print('--- {} seconds ---'.format(time.time() - start_time))
    io = process_io.process_io(input_file=self.input_file, track_tree_name='tree_Particle')
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
    process_base.process_base.saveHistograms(self)

    print('--- {} seconds ---'.format(time.time() - start_time))

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
    
    # Call base class initialization
    process_base.process_base.initialize_config(self)
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)

    # Set configuration for analysis
    self.pTbins = config['pTbins']
    self.jetR_list = config['jetR']
    self.beta_list = config['betas']
    self.n_lambda_bins = config['n_lambda_bins']
  
  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initializeHistograms(self):
    '''
    self.hNevents = ROOT.TH1F('hNevents', 'hNevents', 2, -0.5, 1.5)
    self.hNevents.Fill(1, self.nEvents)
    '''
    for jetR in self.jetR_list:
      '''
      name = 'hJetPt_R{}'.format(jetR)
      h = ROOT.TH1F(name, name, 300, 0, 300)
      h.GetXaxis().SetTitle('p_{T,ch jet}')
      h.GetYaxis().SetTitle('dN/dp_{T}')
      setattr(self, name, h)
      
      name = 'hM_JetPt_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 50.)
      h.GetXaxis().SetTitle('p_{T,ch jet}')
      h.GetYaxis().SetTitle('m_{ch jet}')
      setattr(self, name, h)
      '''
      
      for beta in self.beta_list:

        for i, pTmin in list(enumerate(self.pTbins))[0:-1]:

          # Angularities, \lambda_{\beta}^{\kappa}
          pTmax = self.pTbins[i+1]
          name = "hLambda_pT%i-%i_R%s_#beta%s" % (pTmin, pTmax, jetR, beta)
          h = ROOT.TH1F(name, name, self.n_lambda_bins, 0, 1.0)
          h.GetXaxis().SetTitle('#lambda_{%s}' % beta)
          h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{%s}}' % beta)
          setattr(self, name, h)

          # Angularities with soft drop
          name = "hLambda_pT%i-%i_R%s_#beta%s_SD" % (pTmin, pTmax, jetR, beta)
          h = ROOT.TH1F(name, name, self.n_lambda_bins, 0, 1.0)
          h.GetXaxis().SetTitle('#lambda_{%s}' % beta)
          h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{%s}}' % beta)
          setattr(self, name, h)

          # Lambda vs pT plots to estimate the binning that will be needed
          name = "hLambda_JetpT_R%s_#beta%s" % (jetR, beta)
          h = ROOT.TH2F(name, name, len(self.pTbins) - 1, self.pTbins[0], self.pTbins[-1],
                        self.n_lambda_bins, 0, 1)
          h.GetXaxis().SetTitle('p_{T, jet}')
          h.GetYaxis().SetTitle('#lambda_{%s}' % beta)
          setattr(self, name, h)
          
          # Lambda vs pT plots to estimate the binning that will be needed -- with soft drop
          name = "hLambda_JetpT_R%s_#beta%s_SD" % (jetR, beta)
          h = ROOT.TH2F(name, name, len(self.pTbins) - 1, self.pTbins[0], self.pTbins[-1],
                        self.n_lambda_bins, 0, 1)
          h.GetXaxis().SetTitle('p_{T, jet}')
          h.GetYaxis().SetTitle('#lambda_{%s}' % beta)
          setattr(self, name, h)
          
        '''
        # Retrieve histogram binnings
        n_pt_bins_det = getattr(self, 'n_pt_bins_det_B{}'.format(beta))
        det_pt_bin_array = getattr(self, 'det_pt_bin_array_B{}'.format(beta))
        n_rg_bins_det = getattr(self, 'n_rg_bins_det_B{}'.format(beta))
        det_rg_bin_array = getattr(self, 'det_rg_bin_array_B{}'.format(beta))
        
        name = 'hThetaG_JetPt_R{}_B{}_Rebinned'.format(jetR, beta)
        h = ROOT.TH2F(name, name, n_pt_bins_det, det_pt_bin_array, n_rg_bins_det, det_rg_bin_array)
        h.GetXaxis().SetTitle('p_{T,ch jet}')
        h.GetYaxis().SetTitle('#theta_{g,ch}')
        setattr(self, name, h)

        name = 'hThetaG_JetPt_R{}_B{}'.format(jetR, beta)
        h = ROOT.TH2F(name, name, 300, 0, 300, 150, 0, 1.5)
        h.GetXaxis().SetTitle('p_{T,ch jet}')
        h.GetYaxis().SetTitle('#theta_{g,ch}')
        setattr(self, name, h)

        name = 'hZg_JetPt_R{}_B{}'.format(jetR, beta)
        h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 1.)
        h.GetXaxis().SetTitle('p_{T,ch jet}')
        h.GetYaxis().SetTitle('z_{g,ch}')
        setattr(self, name, h)

        name = 'hMg_JetPt_R{}_B{}'.format(jetR, beta)
        h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 50.)
        h.GetXaxis().SetTitle('p_{T,ch jet}')
        h.GetYaxis().SetTitle('m_{g,ch jet}')
        setattr(self, name, h)
        '''

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
        jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
        print('jet definition is:', jet_def)
        print('jet selector is:', jet_selector,'\n')
        
        # Define SoftDrop settings
        zcut = 0.1
        sd = fjcontrib.SoftDrop(beta, zcut, jetR)
        print('SoftDrop groomer is: {}'.format(sd.description()));

        # Use list comprehension to do jet-finding and fill histograms
        result = [self.analyzeJets(fj_particles, jet_def, jet_selector, sd, beta)
                  for fj_particles in self.df_fjparticles]

  #---------------------------------------------------------------
  # Analyze jets of a given event.
  # fj_particles is the list of fastjet pseudojets for a single fixed event.
  #---------------------------------------------------------------
  def analyzeJets(self, fj_particles, jet_def, jet_selector, sd, beta):
    
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
      
      # Perform SoftDrop grooming and fill histograms
      #jet_sd = sd.result(jet)
      #self.fillSoftDropHistograms(jet_sd, jet, jetR, beta)

  #---------------------------------------------------------------
  # Fill jet histograms
  #---------------------------------------------------------------
  def fillJetHistograms(self, jet, jetR, beta, sd):

    #jet_pt = jet.pt()

    # just use kappa = 1 for now
    l = lambda_beta_kappa(jet, jetR, beta, 1)
    (pTmin, pTmax) = self.pT_bin(jet.pt())

    # jet with softdrop
    jet_sd = sd.result(jet)
    lsd = lambda_beta_kappa(jet_sd, jetR, beta, 1)

    if pTmin:  # pTmin will be None if not a valid bin
      getattr(self, "hLambda_pT%i-%i_R%s_#beta%s" % (pTmin, pTmax, jetR, beta)).Fill(l)
      getattr(self, "hLambda_pT%i-%i_R%s_#beta%s_SD" % (pTmin, pTmax, jetR, beta)).Fill(lsd)
      getattr(self, "hLambda_JetpT_R%s_#beta%s" % (jetR, beta)).Fill(jet.pt(), l)
      getattr(self, "hLambda_JetpT_R%s_#beta%s_SD" % (jetR, beta)).Fill(jet.pt(), lsd)

    '''
    getattr(self, 'hJetPt_R{}'.format(jetR)).Fill(jet_pt)
    getattr(self, 'hM_JetPt_R{}'.format(jetR)).Fill(jet_pt, jet.m())
    '''

  '''
  #---------------------------------------------------------------
  # Fill soft drop histograms
  #---------------------------------------------------------------
    jet_pt_ungroomed = jet.pt()
    
    sd_info = fjcontrib.get_SD_jet_info(jet_sd)
    theta_g = sd_info.dR / jetR
    zg = sd_info.z
    mg = jet_sd.m()
    
    
    getattr(self, 'hThetaG_JetPt_R{}_B{}'.format(jetR, beta)).Fill(jet_pt_ungroomed, theta_g)
    getattr(self, 'hThetaG_JetPt_R{}_B{}_Rebinned'.format(jetR, beta)).Fill(jet_pt_ungroomed, theta_g)
    getattr(self, 'hZg_JetPt_R{}_B{}'.format(jetR, beta)).Fill(jet_pt_ungroomed, zg)
    getattr(self, 'hMg_JetPt_R{}_B{}'.format(jetR, beta)).Fill(jet_pt_ungroomed, mg)
  '''

  #---------------------------------------------------------------
  # Helper function for finding the correct jet pT bin
  #---------------------------------------------------------------
  def pT_bin(self, jet_pT):
    for i, pTmin in list(enumerate(self.pTbins))[0:-1]:
      pTmax = self.pTbins[i+1]
      if pTmin <= jet_pT < pTmax:
        return (pTmin, pTmax)
    return (None, None)
      

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

  analysis = process_ang_data(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_ang_data()
