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
class ProcessData(process_base.ProcessBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
    super(ProcessData, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
    self.initialize_config()

  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def process_data(self):
    
    self.start_time = time.time()

    # Use IO helper class to convert ROOT TTree into a SeriesGroupBy object of fastjet particles per event
    print('--- {} seconds ---'.format(time.time() - self.start_time))
    io = process_io.ProcessIO(input_file=self.input_file, track_tree_name='tree_Particle')
    self.df_fjparticles = io.load_data()
    self.nEvents = len(self.df_fjparticles.index)
    self.nTracks = len(io.track_df.index)
    print('--- {} seconds ---'.format(time.time() - self.start_time))
    
    # Initialize configuration and histograms
    self.initialize_config()
    self.initialize_output_objects()
    
    # Create constituent subtractor, if configured
    if not self.is_pp:
      self.constituent_subtractor = CEventSubtractor(max_distance=self.max_distance, alpha=self.alpha, max_eta=self.max_eta, bge_rho_grid_size=self.bge_rho_grid_size, max_pt_correct=self.max_pt_correct, ghost_area=self.ghost_area, distance_type=fjcontrib.ConstituentSubtractor.deltaR)
    
    print(self)

    # Find jets and fill histograms
    print('Analyze events...')
    self.analyzeEvents()

    # Plot histograms
    print('Save histograms...')
    process_base.ProcessBase.save_output_objects(self)

    print('--- {} seconds ---'.format(time.time() - self.start_time))

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
    
    # Call base class initialization
    process_base.ProcessBase.initialize_config(self)
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    
    # Write tree output (default is to write only histograms)
    self.write_tree_output = config['write_tree_output']
    
    if self.do_constituent_subtraction:
        self.is_pp = False
    else:
        self.is_pp = True
       
    # Check if SoftDrop analysis is activated
    self.do_softdrop = False
    if 'SoftDrop' in config:
        print('zg/rg analysis activated')
        self.do_softdrop = True
        sd_config_dict = config['SoftDrop']
        sd_config_list = list(sd_config_dict.keys())
        self.sd_settings = [[sd_config_dict[name]['zcut'], sd_config_dict[name]['beta']] for name in sd_config_list]
           
    # Check if sub-jet analysis is activated
    self.do_subjets = False
    if 'Subjet' in config:
      print('Subjet analysis activated')
      self.do_subjets = True
      obs_config_dict = config['Subjet']
      obs_config_list = list(obs_config_dict.keys())
      self.subjet_R = [obs_config_dict[name]['subjet_R'] for name in obs_config_list]
      self.subjet_def = {}
      for subjetR in self.subjet_R:
        self.subjet_def[subjetR] = fj.JetDefinition(fj.antikt_algorithm, subjetR)
    
    # Check if jet axis analysis is activated
    self.do_jet_axes = False
    if 'JetAxis' in config:
      print('Jet axis analysis activated')
      self.do_jet_axes = True
      obs_config_dict = config['JetAxis']
      obs_config_list = list(obs_config_dict.keys())
      self.axis_list = [obs_config_dict[name]['axis'] for name in obs_config_list]

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_output_objects(self):
    
    self.hNevents = ROOT.TH1F('hNevents', 'hNevents', 2, -0.5, 1.5)
    self.hNevents.Fill(1, self.nEvents)
    
    self.hTrackEtaPhi = ROOT.TH2F('hTrackEtaPhi', 'hTrackEtaPhi', 200, -1., 1., 628, 0., 6.28)
    self.hTrackPt = ROOT.TH1F('hTrackPt', 'hTrackPt', 300, 0., 300.)
    
    if not self.is_pp:
        self.hRho = ROOT.TH1F('hRho', 'hRho', 1000, 0., 1000.)
        
    for jetR in self.jetR_list:
      
      name = 'hZ_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 1.)
      setattr(self, name, h)
      
      if self.do_softdrop:
          for sd_setting in self.sd_settings:
            
            zcut = sd_setting[0]
            beta = sd_setting[1]
            sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)
            
            if self.write_tree_output:

              # Initialize tree to write out event variables
              tree_name = 't_R{}_{}'.format(self.utils.remove_periods(jetR), sd_label)
              t = ROOT.TTree(tree_name, tree_name)
              tree_writer_name = 'tree_writer_R{}_{}'.format(self.utils.remove_periods(jetR), sd_label)
              tree_writer = treewriter.RTreeWriter(tree=t, tree_name=tree_name, name=tree_writer_name, file_name=None)
              setattr(self, tree_writer_name, tree_writer)

            else:
              
              name = 'hThetaG_JetPt_R{}_{}'.format(jetR, sd_label)
              h = ROOT.TH2F(name, name, 300, 0, 300, 150, 0, 1.5)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle('#theta_{g,ch}')
              setattr(self, name, h)
              
              name = 'hZg_JetPt_R{}_{}'.format(jetR, sd_label)
              h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 1.)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle('z_{g,ch}')
              setattr(self, name, h)
              
      if self.do_subjets:
        for subjetR in self.subjet_R:
        
          name = 'hSubjetZ_JetPt_R{}_{}'.format(jetR, subjetR)
          h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 1.)
          h.GetXaxis().SetTitle('p_{T,ch jet}')
          h.GetYaxis().SetTitle('z_{subjet}')
          setattr(self, name, h)
      
      if self.do_jet_axes:
        for axis in self.axis_list:
        
          name = 'hJetAxisDeltaR_JetPt_Standard_{}_R{}'.format(axis, jetR)
          h = ROOT.TH2F(name, name, 300, 0, 300, 200, 0, jetR)
          h.GetXaxis().SetTitle('p_{T,ch jet}')
          h.GetYaxis().SetTitle('#Delta R')
          setattr(self, name, h)

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

      print('--- {} seconds ---'.format(time.time() - self.start_time))

      # Set jet definition and a jet selector
      jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
      jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
      print('jet definition is:', jet_def)
      print('jet selector is:', jet_selector,'\n')
      
      # Use list comprehension to do jet-finding and fill histograms
      result = [self.analyzeJets(fj_particles, jet_def, jet_selector) for fj_particles in self.df_fjparticles]
        
  #---------------------------------------------------------------
  # Analyze jets of a given event.
  # fj_particles is the list of fastjet pseudojets for a single fixed event.
  #---------------------------------------------------------------
  def analyzeJets(self, fj_particles, jet_def, jet_selector):
    
    # Perform constituent subtraction
    if not self.is_pp:
      fj_particles = self.constituent_subtractor.process_event(fj_particles)
      rho = self.constituent_subtractor.bge_rho.rho()
      getattr(self, 'hRho').Fill(rho)
    
    # Do jet finding
    cs = fj.ClusterSequence(fj_particles, jet_def)
    jets = fj.sorted_by_pt(cs.inclusive_jets())
    jets_selected = jet_selector(jets)

    # Loop through jets and fill non-SD histograms
    jetR = jet_def.R()
    result = [self.analyze_accepted_jets(jet, jetR) for jet in jets_selected]
          
    # Loop through SD settings and fill SD histograms
    if self.do_softdrop:
        result = [[self.analyze_softdrop_jet(sd_setting, jet, jetR) for sd_setting in self.sd_settings] for jet in jets_selected]
    
  #---------------------------------------------------------------
  # Fill histograms
  #---------------------------------------------------------------
  def analyze_accepted_jets(self, jet, jetR):
    
    if self.debug_level > 1:
      print('jet: {} with pt={}'.format(jet, jet.pt()))
    
    # Check additional acceptance criteria
    if not self.utils.is_det_jet_accepted(jet):
      return
      
    # Find subjets
    if self.do_subjets:
      result = [self.analyze_subjets(jet, jetR, subjetR) for subjetR in self.subjet_R]
    
    if not self.write_tree_output:
      self.fill_jet_histograms(jet, jetR)
      
  #---------------------------------------------------------------
  # Fill histograms
  #---------------------------------------------------------------
  def analyze_softdrop_jet(self, sd_setting, jet, jetR):
    
    zcut = sd_setting[0]
    beta = sd_setting[1]
    sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)
    
    sd = fjcontrib.SoftDrop(beta, zcut, jetR)
    jet_def_recluster = fj.JetDefinition(fj.cambridge_algorithm, jetR)
    reclusterer = fjcontrib.Recluster(jet_def_recluster)
    sd.set_reclustering(True, reclusterer)
    setattr(self, 'sd_R{}_{}'.format(jetR, sd_label), sd)
    if self.debug_level > 2:
      print('SoftDrop groomer is: {}'.format(sd.description()))

    #sd = getattr(self, 'sd_R{}_{}'.format(jetR, sd_label))
    
    if self.debug_level > 1:
      print('SD -- jet: {} with pt={}'.format(jet, jet.pt()))
  
    # Check additional acceptance criteria
    if not self.utils.is_det_jet_accepted(jet):
      return

    # Perform SoftDrop grooming
    jet_sd = sd.result(jet)
    
    # Fill tree or histograms (as requested)
    if self.write_tree_output:
      self.fill_tree(jet, jet_sd, jetR, sd_label)
    else:
      self.fill_softdrop_histograms(jet_sd, jet, jetR, sd_label)
      
    # Fill jet axis difference
    if self.do_jet_axes:
      self.fill_jet_axis_histograms(jet, jet_sd, jetR)

  #---------------------------------------------------------------
  # Fill histograms
  #---------------------------------------------------------------
  def fill_jet_histograms(self, jet, jetR):
    
    jet_pt = jet.pt()
    hZ = getattr(self, 'hZ_R{}'.format(jetR))
    for constituent in jet.constituents():
      z = constituent.pt() / jet_pt
      hZ.Fill(jet_pt, z)

  #---------------------------------------------------------------
  # For a given jet, find subjets of a given radius, and fill histograms
  #---------------------------------------------------------------
  def analyze_subjets(self, jet, jetR, subjetR):

    cs_subjet = fj.ClusterSequence(jet.constituents(), self.subjet_def[subjetR])
    subjets = fj.sorted_by_pt(cs_subjet.inclusive_jets())
    for subjet in subjets:
      z = subjet.pt() / jet.pt()
      getattr(self, 'hSubjetZ_JetPt_R{}_{}'.format(jetR, subjetR)).Fill(jet.pt(), z)

  #---------------------------------------------------------------
  # Fill soft drop histograms
  #---------------------------------------------------------------
  def fill_softdrop_histograms(self, jet_sd, jet, jetR, sd_label):
    
    if self.debug_level > 1:
      print('Filling SD histograms...')
    
    jet_pt_ungroomed = jet.pt()

    # SD jet variables
    sd_info = fjcontrib.get_SD_jet_info(jet_sd)
    theta_g = sd_info.dR / jetR
    zg = sd_info.z

    getattr(self, 'hThetaG_JetPt_R{}_{}'.format(jetR, sd_label)).Fill(jet_pt_ungroomed, theta_g)
    getattr(self, 'hZg_JetPt_R{}_{}'.format(jetR, sd_label)).Fill(jet_pt_ungroomed, zg)

    if self.debug_level > 1:
      print('Done.')
          
  #---------------------------------------------------------------
  # Fill jet axis histograms
  #---------------------------------------------------------------
  def fill_jet_axis_histograms(self, jet, jet_sd, jetR):

    for axis in self.axis_list:
      
      if axis == 'SD':
        deltaR = jet.delta_R(jet_sd)
      
      elif axis == 'WTA':
        deltaR = jet.delta_R(self.utils.get_leading_constituent(jet))
      
      getattr(self, 'hJetAxisDeltaR_JetPt_Standard_{}_R{}'.format(axis, jetR)).Fill(jet.pt(), deltaR)
          
  #---------------------------------------------------------------
  # Fill tree
  #---------------------------------------------------------------
  def fill_tree(self, jet, jet_sd, jetR, sd_label):
  
    name = 'tree_writer_R{}_{}'.format(self.utils.remove_periods(jetR), sd_label)
    tree_writer = getattr(self, name)
  
    # Jet variables
    jet_pt = jet.pt()
    tree_writer.fill_branch('jet_pt', jet_pt)
    tree_writer.fill_branch('jet_m', jet.m())
    
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

  analysis = ProcessData(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_data()
