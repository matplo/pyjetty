#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of track information
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
from array import *

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_utils
from pyjetty.alice_analysis.process.base import process_base
from pyjetty.mputils import treewriter
from pyjetty.mputils import CEventSubtractor

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class process_rg_data(process_base.process_base):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
    super(process_rg_data, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
    self.initialize_config()

  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def process_rg_data(self):
    
    self.start_time = time.time()

    # Use IO helper class to convert ROOT TTree into a SeriesGroupBy object of fastjet particles per event
    print('--- {} seconds ---'.format(time.time() - self.start_time))
    io = process_io.process_io(input_file=self.input_file, track_tree_name='tree_Particle')
    self.df_fjparticles = io.load_data()
    self.nEvents = len(self.df_fjparticles.index)
    self.nTracks = len(io.track_df.index)
    print('--- {} seconds ---'.format(time.time() - self.start_time))
    
    # Initialize configuration and histograms
    self.initialize_config()
    self.initialize_output_objects()
    
    # Create constituent subtractor, if configured
    if self.do_constituent_subtraction:
      self.constituent_subtractor = CEventSubtractor(max_distance=self.max_distance, alpha=self.alpha, max_eta=self.max_eta, bge_rho_grid_size=self.bge_rho_grid_size, max_pt_correct=self.max_pt_correct, ghost_area=self.ghost_area, distance_type=fjcontrib.ConstituentSubtractor.deltaR)
    
    print(self)

    # Find jets and fill histograms
    print('Analyze events...')
    self.analyzeEvents()

    # Plot histograms
    print('Save histograms...')
    process_base.process_base.save_output_objects(self)

    print('--- {} seconds ---'.format(time.time() - self.start_time))

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
    
    # Call base class initialization
    process_base.process_base.initialize_config(self)
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    
    # config['beta'] is a dictionary of dictionaries, where each dict is for a value of beta
    beta_dict = config['beta']
    
    # Retrieve list of beta values
    self.beta_list = list(beta_dict.keys())
  
  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_output_objects(self):
    
    self.hNevents = ROOT.TH1F('hNevents', 'hNevents', 2, -0.5, 1.5)
    self.hNevents.Fill(1, self.nEvents)
    
    self.hTrackEtaPhi = ROOT.TH2F('hTrackEtaPhi', 'hTrackEtaPhi', 200, -1., 1., 628, 0., 6.28)
    self.hTrackPt = ROOT.TH1F('hTrackPt', 'hTrackPt', 300, 0., 300.)

    for jetR in self.jetR_list:
      
      for beta in self.beta_list:

        # Initialize tree to write out event variables
        tree_name = 't_R{}_B{}'.format(jetR, beta)
        t = ROOT.TTree(tree_name, tree_name)        
        tree_writer_name = 'tree_writer_R{}_B{}'.format(jetR, beta)
        tree_writer = treewriter.RTreeWriter(tree=t, tree_name=tree_name, name=tree_writer_name, file_name=None)
        setattr(self, tree_writer_name, tree_writer)

  #---------------------------------------------------------------
  # Main function to loop through and analyze events
  #---------------------------------------------------------------
  def analyzeEvents(self):
    
    # Fill track histograms
    print('--- {} seconds ---'.format(time.time() - self.start_time))
    print('Fill track histograms')
    [[self.fillTrackHistograms(track) for track in fj_particles] for fj_particles in self.df_fjparticles]
    print('--- {} seconds ---'.format(time.time() - self.start_time))
    
    print('Find jets...')
    fj.ClusterSequence.print_banner()
    print()
    
    for jetR in self.jetR_list:
      
      for beta in self.beta_list:
        
        print('--- {} seconds ---'.format(time.time() - self.start_time))

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
        result = [self.analyzeJets(fj_particles, jet_def, jet_selector, sd, beta) for fj_particles in self.df_fjparticles]
        
  #---------------------------------------------------------------
  # Analyze jets of a given event.
  # fj_particles is the list of fastjet pseudojets for a single fixed event.
  #---------------------------------------------------------------
  def analyzeJets(self, fj_particles, jet_def, jet_selector, sd, beta):
    
    # Perform constituent subtraction
    if self.do_constituent_subtraction:
      fj_particles = self.constituent_subtractor.process_event(fj_particles)
      rho = self.constituent_subtractor.bge_rho.rho()
    
    # Do jet finding
    cs = fj.ClusterSequence(fj_particles, jet_def)
    jets = fj.sorted_by_pt(cs.inclusive_jets())
    jets_selected = jet_selector(jets)

    # Loop through jets
    jetR = jet_def.R()
    for jet in jets_selected:
      
      if self.debug_level > 1:
        print('jet: {} with pt={}'.format(jet, jet.pt()))
      
      # Check additional acceptance criteria
      if not self.utils.is_det_jet_accepted(jet):
        continue
      
      # Perform SoftDrop grooming and fill tree
      jet_sd = sd.result(jet)
      self.fill_tree(jet, jet_sd, jetR, beta)

  #---------------------------------------------------------------
  # Fill tree
  #---------------------------------------------------------------
  def fill_tree(self, jet, jet_sd, jetR, beta):
  
    name = 'tree_writer_R{}_B{}'.format(jetR, beta)
    tree_writer = getattr(self, name)
  
    # Jet variables
    jet_pt = jet.pt()
    tree_writer.fill_branch('jet_pt', jet_pt)
    tree_writer.fill_branch('jet_m', jet.m())
    [tree_writer.fill_branch('z', constituent.pt() / jet_pt) for constituent in jet.constituents()]
    
    # SD jet variables
    sd_info = fjcontrib.get_SD_jet_info(jet_sd)
    tree_writer.fill_branch('theta_g', sd_info.dR / jetR)
    tree_writer.fill_branch('zg',sd_info.z)
    tree_writer.fill_branch('mg', jet_sd.m())
  
    tree_writer.fill_tree()

  #---------------------------------------------------------------
  # Fill track histograms.
  #---------------------------------------------------------------
  def fillTrackHistograms(self, track):
    
    self.hTrackEtaPhi.Fill(track.eta(), track.phi())
    self.hTrackPt.Fill(track.pt())

##################################################################
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Process data')
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
                      help='Output directory for output to be written to')
  
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

  analysis = process_rg_data(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_rg_data()
