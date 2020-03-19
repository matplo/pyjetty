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
    
    if self.do_constituent_subtraction:
        self.is_pp = False
    else:
        self.is_pp = True
       
    self.observable_list = config['process_observables']
    
    # Create dictionaries to store SD settings and observable settings for each observable
    # Each dictionary entry stores a list of subconfiguration parameters
    #   The SD list stores a list of SD settings [zcut, beta]
    #   The observable list stores a the observable setting, e.g. subjetR
    self.obs_sd_settings = {}
    self.obs_settings = {}
    
    for observable in self.observable_list:
        
      # Fill SD settings
      self.obs_sd_settings[observable] = []
      obs_config_dict = config[observable]
      for config_key, subconfig in obs_config_dict.items():
        if 'SoftDrop' in subconfig:
          sd_dict = obs_config_dict[config_key]['SoftDrop']
          self.obs_sd_settings[observable].append([sd_dict['zcut'], sd_dict['beta']])
        else:
          self.obs_sd_settings[observable].append(None)
      
      # Fill observable settings
      self.obs_settings[observable] = []
      obs_config_list = list(obs_config_dict.keys())
      if observable == 'subjet_z':
        self.obs_settings[observable] = [obs_config_dict[name]['subjet_R'] for name in obs_config_list]
        self.subjet_def = {}
        for subjetR in self.obs_settings[observable]:
          self.subjet_def[subjetR] = fj.JetDefinition(fj.antikt_algorithm, subjetR)
      if observable == 'jet_axis':
        self.obs_settings[observable] = [obs_config_dict[name]['axis'] for name in obs_config_list]
        
    # Construct set of unique SD settings
    self.sd_settings = []
    lists = [self.obs_sd_settings[obs] for obs in self.observable_list]
    for observable in lists:
      for setting in observable:
        if setting not in self.sd_settings and setting != None:
          self.sd_settings.append(setting)
    
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
      
      for observable in self.observable_list:

        if observable == 'theta_g':
        
          for sd_setting in self.obs_sd_settings[observable]:
            sd_label = self.utils.sd_label(sd_setting)
            name = 'h_{}_JetPt_R{}_{}'.format(observable, jetR, sd_label)
            h = ROOT.TH2F(name, name, 300, 0, 300, 150, 0, 1.0)
            h.GetXaxis().SetTitle('p_{T,ch jet}')
            h.GetYaxis().SetTitle('#theta_{g,ch}')
            setattr(self, name, h)
            
        if observable == 'zg':
        
          for sd_setting in self.obs_sd_settings[observable]:
            sd_label = self.utils.sd_label(sd_setting)
            name = 'h_{}_JetPt_R{}_{}'.format(observable, jetR, sd_label)
            h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 0.5)
            h.GetXaxis().SetTitle('p_{T,ch jet}')
            h.GetYaxis().SetTitle('z_{g,ch}')
            setattr(self, name, h)
              
        if observable == 'subjet_z':
        
          for subjetR in self.obs_settings[observable]:
            name = 'h_{}_JetPt_R{}_{}'.format(observable, jetR, subjetR)
            h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 1.)
            h.GetXaxis().SetTitle('p_{T,ch jet}')
            h.GetYaxis().SetTitle('z_{subjet}')
            setattr(self, name, h)
            
        if observable == 'jet_axis':
              
          for i, axes in enumerate(self.obs_settings[observable]):
          
            sd_setting = self.obs_sd_settings[observable][i]
            if sd_setting:
              sd_label = self.utils.sd_label(sd_setting)
            else:
              sd_label = ''
            
            name = 'h_{}_JetPt_R{}_{}{}'.format(observable, axes, sd_label, jetR)
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
    if len(self.sd_settings) > 0:
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
    if 'subjet_z' in self.observable_list:
      result = [self.analyze_subjets(jet, jetR, subjetR) for subjetR in self.obs_settings['subjet_z']]
    
    self.fill_jet_histograms(jet, jetR)
      
  #---------------------------------------------------------------
  # Fill histograms
  #---------------------------------------------------------------
  def analyze_softdrop_jet(self, sd_setting, jet, jetR):
    
    zcut = sd_setting[0]
    beta = sd_setting[1]
    sd_label = self.utils.sd_label(sd_setting)
    
    sd = fjcontrib.SoftDrop(beta, zcut, jetR)
    jet_def_recluster = fj.JetDefinition(fj.cambridge_algorithm, jetR)
    reclusterer = fjcontrib.Recluster(jet_def_recluster)
    sd.set_reclustering(True, reclusterer)
    setattr(self, 'sd_R{}_{}'.format(jetR, sd_label), sd)
    if self.debug_level > 2:
      print('SoftDrop groomer is: {}'.format(sd.description()))
    
    if self.debug_level > 1:
      print('SD -- jet: {} with pt={}'.format(jet, jet.pt()))
  
    # Check additional acceptance criteria
    if not self.utils.is_det_jet_accepted(jet):
      return

    # Perform SoftDrop grooming
    jet_sd = sd.result(jet)
    
    # Fill histograms
    self.fill_softdrop_histograms(jet_sd, jet, jetR, sd_setting, sd_label)
      
    # Fill jet axis difference
    if 'jet_axis' in self.observable_list:

        # Recluster with WTA (with larger jet R)
        jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
        jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
        if self.debug_level > 2:
            print('WTA jet definition is:', jet_def)
        reclusterer_wta =  fjcontrib.Recluster(jet_def_wta)
        jet_wta = reclusterer_wta.result(jet)

        self.fill_jet_axis_histograms(jet, jet_sd, jet_wta, jetR, sd_setting, sd_label)

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
      getattr(self, 'h_subjet_z_JetPt_R{}_{}'.format(jetR, subjetR)).Fill(jet.pt(), z)

  #---------------------------------------------------------------
  # Fill soft drop histograms
  #---------------------------------------------------------------
  def fill_softdrop_histograms(self, jet_sd, jet, jetR, sd_setting, sd_label):
    
    if self.debug_level > 1:
      print('Filling SD histograms...')
    
    jet_pt_ungroomed = jet.pt()

    # SD jet variables
    sd_info = fjcontrib.get_SD_jet_info(jet_sd)
    theta_g = sd_info.dR / jetR
    zg = sd_info.z

    if sd_setting in self.obs_sd_settings['theta_g']:
      getattr(self, 'h_theta_g_JetPt_R{}_{}'.format(jetR, sd_label)).Fill(jet_pt_ungroomed, theta_g)
    if sd_setting in self.obs_sd_settings['zg']:
      getattr(self, 'h_zg_JetPt_R{}_{}'.format(jetR, sd_label)).Fill(jet_pt_ungroomed, zg)

    if self.debug_level > 1:
      print('Done.')
          
  #---------------------------------------------------------------
  # Fill jet axis histograms
  #---------------------------------------------------------------
  def fill_jet_axis_histograms(self, jet, jet_sd, jet_wta, jetR, sd_setting, sd_label):

    for axis in self.obs_settings['jet_axis']:

      if axis == 'Standard_SD':
        if sd_setting in self.obs_sd_settings['jet_axis']:
          deltaR = jet.delta_R(jet_sd)
          getattr(self, 'h_jet_axis_JetPt_R{}_Standard_SD{}'.format(sd_label, jetR)).Fill(jet.pt(), deltaR)
          
      if axis == 'Standard_WTA':
        if sd_setting == self.sd_settings[0]:
          deltaR = jet.delta_R(jet_wta)
          getattr(self, 'h_jet_axis_JetPt_R{}_Standard_WTA'.format(jetR)).Fill(jet.pt(), deltaR)
        
      if axis == 'WTA_SD':
        if sd_setting in self.obs_sd_settings['jet_axis']:
          deltaR = jet_sd.delta_R(jet_wta)
          getattr(self, 'h_jet_axis_JetPt_R{}_WTA_SD{}'.format(sd_label, jetR)).Fill(jet.pt(), deltaR)

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
