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
import random

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext
import fjtools

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_io_emb
from pyjetty.alice_analysis.process.base import process_utils
from pyjetty.alice_analysis.process.base import jet_info
from pyjetty.alice_analysis.process.base import process_base
from pyjetty.alice_analysis.process.base import thermal_generator
from pyjetty.mputils import CEventSubtractor

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class ProcessMC(process_base.ProcessBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
    super(ProcessMC, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
    
    # Initialize configuration
    self.initialize_config()
  
  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def process_mc(self):
    
    self.start_time = time.time()
    
    # ------------------------------------------------------------------------
    
    # Use IO helper class to convert detector-level ROOT TTree into
    # a SeriesGroupBy object of fastjet particles per event
    print('--- {} seconds ---'.format(time.time() - self.start_time))
    if self.fast_simulation:
      tree_dir = ''
    else:
      tree_dir = 'PWGHF_TreeCreator'
    io_det = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
                                  track_tree_name='tree_Particle', use_ev_id_ext=self.use_ev_id_ext)
    df_fjparticles_det = io_det.load_data(reject_tracks_fraction=self.reject_tracks_fraction)
    self.nEvents_det = len(df_fjparticles_det.index)
    self.nTracks_det = len(io_det.track_df.index)
    print('--- {} seconds ---'.format(time.time() - self.start_time))
    
    # ------------------------------------------------------------------------

    # Use IO helper class to convert truth-level ROOT TTree into
    # a SeriesGroupBy object of fastjet particles per event
    io_truth = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
                                    track_tree_name='tree_Particle_gen', use_ev_id_ext=self.use_ev_id_ext)
    df_fjparticles_truth = io_truth.load_data()
    self.nEvents_truth = len(df_fjparticles_truth.index)
    self.nTracks_truth = len(io_truth.track_df.index)
    print('--- {} seconds ---'.format(time.time() - self.start_time))
    
    # ------------------------------------------------------------------------

    # Now merge the two SeriesGroupBy to create a groupby df with [ev_id, run_number, fj_1, fj_2]
    # (Need a structure such that we can iterate event-by-event through both fj_1, fj_2 simultaneously)
    print('Merge det-level and truth-level into a single dataframe grouped by event...')
    self.df_fjparticles = pandas.concat([df_fjparticles_det, df_fjparticles_truth], axis=1)
    self.df_fjparticles.columns = ['fj_particles_det', 'fj_particles_truth']
    print('--- {} seconds ---'.format(time.time() - self.start_time))

    # ------------------------------------------------------------------------
    
    # Set up the Pb-Pb embedding object
    if not self.is_pp and not self.thermal_model:
        self.process_io_emb = process_io_emb.ProcessIO_Emb(self.emb_file_list, track_tree_name='tree_Particle')
    
    # ------------------------------------------------------------------------

    # Initialize histograms
    if not self.dry_run:
      self.initialize_output_objects()
    
    # Create constituent subtractor, if configured
    if self.do_constituent_subtraction:
      self.constituent_subtractor = [CEventSubtractor(max_distance=R_max, alpha=self.alpha, max_eta=self.max_eta, bge_rho_grid_size=self.bge_rho_grid_size, max_pt_correct=self.max_pt_correct, ghost_area=self.ghost_area, distance_type=fjcontrib.ConstituentSubtractor.deltaR) for R_max in self.max_distance]
    
    print(self)
    
    # Find jets and fill histograms
    print('Find jets...')
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
      
    self.fast_simulation = config['fast_simulation']
    self.dry_run = config['dry_run']
    self.skip_deltapt_RC_histograms = True
    
    self.jet_matching_distance = config['jet_matching_distance']
    self.reject_tracks_fraction = config['reject_tracks_fraction']
    if 'mc_fraction_threshold' in config:
      self.mc_fraction_threshold = config['mc_fraction_threshold']
    self.use_ev_id_ext = config['use_ev_id_ext']
    
    if self.do_constituent_subtraction:
        self.is_pp = False
        self.emb_file_list = config['emb_file_list']
        self.main_R_max = config['constituent_subtractor']['main_R_max']
    else:
        self.is_pp = True
        
    if 'thermal_model' in config:
      self.thermal_model = True
      beta = config['thermal_model']['beta']
      N_avg = config['thermal_model']['N_avg']
      sigma_N = config['thermal_model']['sigma_N']
      self.thermal_generator = thermal_generator.ThermalGenerator(N_avg, sigma_N, beta)
    else:
      self.thermal_model = False
     
    self.observable_list = config['process_observables']
    
    # Create dictionaries to store grooming settings and observable settings for each observable
    # Each dictionary entry stores a list of subconfiguration parameters
    #   The observable list stores the observable setting, e.g. subjetR
    #   The grooming list stores a list of SD or DG settings {'sd': [zcut, beta]} or {'dg': [a]}
    self.obs_settings = {}
    self.obs_grooming_settings = {}
    
    for observable in self.observable_list:
      
      # Fill observable settings
      self.obs_settings[observable] = []
      obs_config_dict = config[observable]
      obs_config_dict = config[observable]
      obs_config_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
     
      if observable == 'subjet_z':
        self.obs_settings[observable] = [obs_config_dict[name]['subjet_R'] for name in obs_config_list]
        self.subjet_def = {}
        for subjetR in self.obs_settings[observable]:
          self.subjet_def[subjetR] = fj.JetDefinition(fj.antikt_algorithm, subjetR)
      if observable == 'jet_axis':
        self.obs_settings[observable] = [obs_config_dict[name]['axis'] for name in obs_config_list]
         
      # Fill grooming settings
      self.obs_grooming_settings[observable] = self.utils.grooming_settings(obs_config_dict)
      
    # Construct set of unique grooming settings
    self.grooming_settings = []
    lists_grooming = [self.obs_grooming_settings[obs] for obs in self.observable_list]
    for observable in lists_grooming:
      for setting in observable:
        if setting not in self.grooming_settings and setting != None:
          self.grooming_settings.append(setting)

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_output_objects(self):
    
    self.hNevents = ROOT.TH1F('hNevents', 'hNevents', 2, -0.5, 1.5)
    self.hNevents.Fill(1, self.nEvents_det)
    
    self.hTrackEtaPhi = ROOT.TH2F('hTrackEtaPhi', 'hTrackEtaPhi', 200, -1., 1., 628, 0., 6.28)
    self.hTrackPt = ROOT.TH1F('hTrackPt', 'hTrackPt', 300, 0., 300.)
    
    if not self.is_pp:
      self.hRho =  ROOT.TH1F('hRho', 'hRho', 1000, 0., 1000.)
      
    name = 'hN_MeanPt'
    h = ROOT.TH2F(name, name, 200, 0, 5000, 200, 0., 2.)
    setattr(self, name, h)

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_output_objects_R(self, jetR):
        
      name = 'hJES_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
      setattr(self, name, h)
      
      if self.is_pp:
          name = 'hDeltaR_All_R{}'.format(jetR)
          h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
          setattr(self, name, h)
          
      else:
      
          for R_max in self.max_distance:
          
            name = 'hDeltaPt_emb_R{}_Rmax{}'.format(jetR, R_max)
            h = ROOT.TH2F(name, name, 300, 0, 300, 400, -200., 200.)
            setattr(self, name, h)
            
            if not self.skip_deltapt_RC_histograms:
              name = 'hDeltaPt_RC_beforeCS_R{}_Rmax{}'.format(jetR, R_max)
              h = ROOT.TH1F(name, name, 400, -200., 200.)
              setattr(self, name, h)
              
              name = 'hDeltaPt_RC_afterCS_R{}_Rmax{}'.format(jetR, R_max)
              h = ROOT.TH1F(name, name, 400, -200., 200.)
              setattr(self, name, h)
      
            name = 'hDeltaR_ppdet_pptrue_R{}_Rmax{}'.format(jetR, R_max)
            h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
            setattr(self, name, h)
            
            name = 'hDeltaR_combined_ppdet_R{}_Rmax{}'.format(jetR, R_max)
            h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
            setattr(self, name, h)
          
            if 'subjet_z' in self.observable_list:
            
              for subjetR in self.obs_settings['subjet_z']:
            
                name = 'hDeltaR_combined_ppdet_{}_R{}_{}_Rmax{}'.format('subjet_z', jetR, subjetR, R_max)
                h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
                setattr(self, name, h)
                
                name = 'hDeltaR_ppdet_pptrue_{}_R{}_{}_Rmax{}'.format('subjet_z', jetR, subjetR, R_max)
                h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
                setattr(self, name, h)
              
      name = 'hZ_Truth_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 1.)
      setattr(self, name, h)
      
      name = 'hZ_Det_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 1.)
      setattr(self, name, h)
      
      for observable in self.observable_list:

        if observable == 'theta_g':

          for grooming_setting in self.obs_grooming_settings[observable]:
            if grooming_setting:
              grooming_label = self.utils.grooming_label(grooming_setting)
              
              if self.is_pp:
                self.create_theta_g_histograms(observable, jetR, grooming_label)
              else:
                for R_max in self.max_distance:
                  self.create_theta_g_histograms(observable, jetR, grooming_label, R_max)
                  if R_max == self.main_R_max:
                    self.create_theta_g_histograms(observable, jetR, grooming_label, '{}_matched'.format(R_max))
              
              if self.thermal_model:
                for R_max in self.max_distance:
                  name = 'h_{}_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, grooming_label, R_max)
                  h = ROOT.TH2F(name, name, 200, 0, 200, 100, 0, 1.0)
                  h.GetXaxis().SetTitle('p_{T,ch jet}')
                  h.GetYaxis().SetTitle('#theta_{g,ch}')
                  setattr(self, name, h)
                  
              name = 'h_{}_JetPt_Truth_R{}_{}'.format(observable, jetR, grooming_label)
              h = ROOT.TH2F(name, name, 20, 0, 200, 100, 0, 1.0)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle('#theta_{g,ch}')
              setattr(self, name, h)

        if observable == 'zg':
        
          for grooming_setting in self.obs_grooming_settings[observable]:
            if grooming_setting:
              grooming_label = self.utils.grooming_label(grooming_setting)
              
              if self.is_pp:
                self.create_zg_histograms(observable, jetR, grooming_label)
              else:
                for R_max in self.max_distance:
                  self.create_zg_histograms(observable, jetR, grooming_label, R_max)
                  if R_max == self.main_R_max:
                    self.create_zg_histograms(observable, jetR, grooming_label, '{}_matched'.format(R_max))

              if self.thermal_model:
                for R_max in self.max_distance:
                  name = 'h_{}_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, grooming_label, R_max)
                  h = ROOT.TH2F(name, name, 200, 0, 200, 100, 0, 0.5)
                  h.GetXaxis().SetTitle('p_{T,ch jet}')
                  h.GetYaxis().SetTitle('z_{g,ch}')
                  setattr(self, name, h)
                  
              name = 'h_{}_JetPt_Truth_R{}_{}'.format(observable, jetR, grooming_label)
              h = ROOT.TH2F(name, name, 20, 0, 200, 50, 0, 0.5)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle('z_{g,ch}')
              setattr(self, name, h)
              
        if observable == 'subjet_z':

          for subjetR in self.obs_settings[observable]:
            
            name = 'hDeltaR_All_{}_R{}_{}'.format(observable, jetR, subjetR)
            h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
            setattr(self, name, h)
          
            name = 'hResidual_JetPt_{}_R{}_{}'.format(observable, jetR, subjetR)
            h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1.0, 1.0)
            h.GetXaxis().SetTitle('p_{T,truth}')
            h.GetYaxis().SetTitle('#frac{z_{det}-z_{truth}}{z_{truth}}')
            setattr(self, name, h)
            
            # Create THn of response for subjet z
            dim = 4;
            title = ['p_{T,det}', 'p_{T,truth}', 'z_{det}', 'z_{truth}']
            nbins = [30, 30, 100, 50]
            min = [0., 0., 0., 0.]
            max = [150., 300., 1., 1.]
            name = 'hResponse_JetPt_{}_R{}_{}'.format(observable, jetR, subjetR)
            self.create_thn(name, title, dim, nbins, min, max)
            
        if observable == 'jet_axis':
              
          for i, axes in enumerate(self.obs_settings[observable]):
          
            grooming_setting = self.obs_grooming_settings[observable][i]
            if grooming_setting:
              grooming_label = self.utils.grooming_label(grooming_setting)
            else:
              grooming_label = ''
            
            name = 'hResidual_JetPt_{}_R{}_{}{}'.format(observable, jetR, axes, grooming_label)
            h = ROOT.TH2F(name, name, 300, 0, 300, 100, -1*jetR, jetR)
            h.GetXaxis().SetTitle('p_{T,truth}')
            h.GetYaxis().SetTitle('#frac{#DeltaR_{det}-#DeltaR_{truth}}{#DeltaR_{truth}}')
            setattr(self, name, h)

            # Create THn of response for jet axis deltaR
            dim = 4;
            title = ['p_{T,det}', 'p_{T,truth}', '#DeltaR_{det}', '#DeltaR_{truth}']
            nbins = [30, 30, 80, 40]
            min = [0., 0., 0., 0.]
            max = [150., 300., jetR, jetR]
            name = 'hResponse_JetPt_{}_R{}_{}{}'.format(observable, jetR, axes, grooming_label)
            self.create_thn(name, title, dim, nbins, min, max)

      # Plot some Lund planes
      for grooming_setting in self.grooming_settings:
        if grooming_setting:
          grooming_label = self.utils.grooming_label(grooming_setting)
          if self.is_pp:
            self.create_lund_histograms(jetR, grooming_label)
          else:
            for R_max in self.max_distance:
              self.create_lund_histograms(jetR, grooming_label, R_max)

      # Create prong matching histograms
      for grooming_setting in self.obs_grooming_settings['theta_g']:
        if grooming_setting:
          grooming_label = self.utils.grooming_label(grooming_setting)
          if not self.is_pp:
            self.create_prong_matching_histograms(jetR, grooming_label)
            
  #---------------------------------------------------------------
  # Create Lund plane histograms
  #---------------------------------------------------------------
  def create_lund_histograms(self, jetR, grooming_label, R_max = None):
  
    if R_max:
      suffix = '_Rmax{}'.format(R_max)
    else:
      suffix = ''
  
    name = 'hLundPlane_R{}_{}{}'.format(jetR, grooming_label, suffix)
    h = ROOT.TH2F(name, name, 140, 0, 7, 100, -4., 6.)
    h.GetXaxis().SetTitle('log(1/ #DeltaR)')
    h.GetYaxis().SetTitle('log(k_{t})')
    setattr(self, name, h)

  #---------------------------------------------------------------
  # Create theta_g response histograms
  #---------------------------------------------------------------
  def create_prong_matching_histograms(self, jetR, grooming_label):
  
    prong_list = ['leading', 'subleading']
    match_list = ['leading', 'subleading', 'ungroomed', 'outside']

    for R_max in self.max_distance:
      for prong in prong_list:
        for match in match_list:

          name = 'hProngMatching_{}_{}_JetPt_R{}_{}_Rmax{}'.format(prong, match, jetR, grooming_label, R_max)
          h = ROOT.TH3F(name, name, 20, 0, 200, 15, -0.4, 1.1, 20, 0., 2*jetR)
          h.GetXaxis().SetTitle('p_{T,truth}')
          h.GetYaxis().SetTitle('Prong matching fraction')
          h.GetZaxis().SetTitle('#Delta R_{prong}')
          setattr(self, name, h)
          
          name = 'hProngMatching_{}_{}_JetPtDet_R{}_{}_Rmax{}'.format(prong, match, jetR, grooming_label, R_max)
          h = ROOT.TH3F(name, name, 20, 0, 200, 15, -0.4, 1.1, 20, 0., 2*jetR)
          h.GetXaxis().SetTitle('p_{T,pp-det}')
          h.GetYaxis().SetTitle('Prong matching fraction')
          h.GetZaxis().SetTitle('#Delta R_{prong}')
          setattr(self, name, h)
          
          name = 'hProngMatching_{}_{}_JetPtZ_R{}_{}_Rmax{}'.format(prong, match, jetR, grooming_label, R_max)
          h = ROOT.TH3F(name, name, 20, 0, 200, 15, -0.4, 1.1, 50, -0.5, 0.5)
          h.GetXaxis().SetTitle('p_{T,truth}')
          h.GetYaxis().SetTitle('Prong matching fraction')
          h.GetZaxis().SetTitle('#Delta z_{prong}')
          setattr(self, name, h)

      name = 'hProngMatching_subleading-leading_correlation_JetPtDet_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)
      h = ROOT.TH3F(name, name, 20, 0, 200, 15, -0.4, 1.1, 15, -0.4, 1.1)
      h.GetXaxis().SetTitle('p_{T,pp-det}')
      h.GetYaxis().SetTitle('Prong matching fraction, leading_subleading')
      h.GetZaxis().SetTitle('Prong matching fraction, subleading_leading')
      setattr(self, name, h)

  #---------------------------------------------------------------
  # Create theta_g response histograms
  #---------------------------------------------------------------
  def create_theta_g_histograms(self, observable, jetR, grooming_label, R_max = None):

    if R_max:
      suffix = '_Rmax{}'.format(R_max)
    else:
      suffix = ''
            
    # Create THn of response for theta_g
    dim = 4;
    title = ['p_{T,det}', 'p_{T,truth}', '#theta_{g,det}', '#theta_{g,truth}']
    nbins = [30, 20, 100, 100]
    min = [0., 0., 0., 0.]
    max = [150., 200., 1., 1.]
    name = 'hResponse_JetPt_{}_R{}_{}{}'.format(observable, jetR, grooming_label, suffix)
    self.create_thn(name, title, dim, nbins, min, max)
    
    name = 'hResidual_JetPt_{}_R{}_{}{}'.format(observable, jetR, grooming_label, suffix)
    h = ROOT.TH3F(name, name, 20, 0, 200, 100, 0., 1., 200, -2., 2.)
    h.GetXaxis().SetTitle('p_{T,truth}')
    h.GetYaxis().SetTitle('#theta_{g,truth}')
    h.GetZaxis().SetTitle('#frac{#theta_{g,det}-#theta_{g,truth}}{#theta_{g,truth}}')
    setattr(self, name, h)

  #---------------------------------------------------------------
  # Create theta_g response histograms
  #---------------------------------------------------------------
  def create_zg_histograms(self, observable, jetR, grooming_label, R_max = None):

    if R_max:
      suffix = '_Rmax{}'.format(R_max)
    else:
      suffix = ''
    
    # Create THn of response for z_g
    dim = 4;
    title = ['p_{T,det}', 'p_{T,truth}', 'z_{g,det}', 'z_{g,truth}']
    nbins = [30, 20, 50, 50]
    min = [0., 0., 0., 0.]
    max = [150., 200., 0.5, 0.5]
    name = 'hResponse_JetPt_{}_R{}_{}{}'.format(observable, jetR, grooming_label, suffix)
    self.create_thn(name, title, dim, nbins, min, max)

    name = 'hResidual_JetPt_{}_R{}_{}{}'.format(observable, jetR, grooming_label, suffix)
    h = ROOT.TH3F(name, name, 20, 0, 200, 50, 0., 0.5, 200, -2., 2.)
    h.GetXaxis().SetTitle('p_{T,truth}')
    h.GetYaxis().SetTitle('z_{g,truth}')
    h.GetZaxis().SetTitle('#frac{z_{g,det}-z_{g,truth}}{z_{g,truth}}')
    setattr(self, name, h)

  #---------------------------------------------------------------
  # Main function to loop through and analyze events
  #---------------------------------------------------------------
  def analyzeEvents(self):
    
    # Fill track histograms
    if not self.dry_run:
      [self.fillTrackHistograms(fj_particles_det) for fj_particles_det in self.df_fjparticles['fj_particles_det']]
    
    fj.ClusterSequence.print_banner()
    print()
        
    self.event_number = 0
    
    for jetR in self.jetR_list:
      if not self.dry_run:
        self.initialize_output_objects_R(jetR)
    
    # Then can use list comprehension to iterate over the groupby and do jet-finding
    # simultaneously for fj_1 and fj_2 per event, so that I can match jets -- and fill histograms
    result = [self.analyze_event(fj_particles_det, fj_particles_truth) for fj_particles_det, fj_particles_truth in zip(self.df_fjparticles['fj_particles_det'], self.df_fjparticles['fj_particles_truth'])]
    
    if self.debug_level > 0:
      for attr in dir(self):
        obj = getattr(self, attr)
        print('size of {}: {}'.format(attr, sys.getsizeof(obj)))
        
    print('Save thn...')
    process_base.ProcessBase.save_thn_th3_objects(self)
    
  #---------------------------------------------------------------
  # Fill track histograms.
  #---------------------------------------------------------------
  def fillTrackHistograms(self, fj_particles_det):

    # Check that the entries exist appropriately
    # (need to check how this can happen -- but it is only a tiny fraction of events)
    if type(fj_particles_det) != fj.vectorPJ:
      return
    
    for track in fj_particles_det:
      self.hTrackEtaPhi.Fill(track.eta(), track.phi())
      self.hTrackPt.Fill(track.pt())
      
  #---------------------------------------------------------------
  # Analyze jets of a given event.
  # fj_particles is the list of fastjet pseudojets for a single fixed event.
  #---------------------------------------------------------------
  def analyze_event(self, fj_particles_det, fj_particles_truth):
  
    self.event_number += 1
    if self.event_number > self.event_number_max:
      return
    if self.debug_level > 1:
      print('-------------------------------------------------')
      print('event {}'.format(self.event_number))
      
    # Check that the entries exist appropriately
    # (need to check how this can happen -- but it is only a tiny fraction of events)
    if type(fj_particles_det) != fj.vectorPJ or type(fj_particles_truth) != fj.vectorPJ:
      print('fj_particles type mismatch -- skipping event')
      return
        
    # If Pb-Pb, construct embedded event (do this once, for all jetR)
    if not self.is_pp:
        
        # If thermal model, generate a thermal event and add it to the det-level particle list
        if self.thermal_model:
          fj_particles_combined_beforeCS = self.thermal_generator.load_event()
          
          # Form the combined det-level event
          # The pp-det tracks are each stored with a unique user_index >= 0
          #   (same index in fj_particles_combined and fj_particles_det -- which will be used in prong-matching)
          # The thermal tracks are each stored with a unique user_index < 0
          [fj_particles_combined_beforeCS.push_back(p) for p in fj_particles_det]

        # Main case: Get Pb-Pb event and embed it into the det-level particle list
        else:
          fj_particles_combined_beforeCS = self.process_io_emb.load_event()
              
          # Form the combined det-level event
          # The pp-det tracks are each stored with a unique user_index >= 0
          #   (same index in fj_particles_combined and fj_particles_det -- which will be used in prong-matching)
          # The Pb-Pb tracks are each stored with a unique user_index < 0
          [fj_particles_combined_beforeCS.push_back(p) for p in fj_particles_det]
         
        # Perform constituent subtraction for each R_max
        fj_particles_combined = [self.constituent_subtractor[i].process_event(fj_particles_combined_beforeCS) for i, R_max in enumerate(self.max_distance)]
        
        if self.debug_level > 3:
          print([p.user_index() for p in fj_particles_truth])
          print([p.pt() for p in fj_particles_truth])
          print([p.user_index() for p in fj_particles_det])
          print([p.pt() for p in fj_particles_det])
          print([p.user_index() for p in fj_particles_combined_beforeCS])
          print([p.pt() for p in fj_particles_combined_beforeCS])
          
    if self.dry_run:
      return

    # Loop through jetR, and process event for each R
    for jetR in self.jetR_list:
    
      # Keep track of whether to fill R-independent histograms
      self.fill_R_indep_hists = (jetR == self.jetR_list[0])

      # Set jet definition and a jet selector
      jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
      jet_selector_det = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
      jet_selector_truth_matched = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9)
      if self.debug_level > 2:
        print('')
        print('jet definition is:', jet_def)
        print('jet selector for det-level is:', jet_selector_det)
        print('jet selector for truth-level matches is:', jet_selector_truth_matched)
      
      # Analyze
      if self.is_pp:
      
        # Find pp det and truth jets
        cs_det = fj.ClusterSequence(fj_particles_det, jet_def)
        jets_det_pp = fj.sorted_by_pt(cs_det.inclusive_jets())
        jets_det_pp_selected = jet_selector_det(jets_det_pp)
        
        cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)
        jets_truth = fj.sorted_by_pt(cs_truth.inclusive_jets())
        jets_truth_selected = jet_selector_det(jets_truth)
        jets_truth_selected_matched = jet_selector_truth_matched(jets_truth)
      
        self.analyze_jets(jets_det_pp_selected, jets_truth_selected, jets_truth_selected_matched, jetR)
        
      else:
      
        for i, R_max in enumerate(self.max_distance):
            
          if self.debug_level > 1:
            print('')
            print('R_max: {}'.format(R_max))
            print('Total number of combined particles: {}'.format(len([p.pt() for p in fj_particles_combined_beforeCS])))
            print('After constituent subtraction {}: {}'.format(i, len([p.pt() for p in fj_particles_combined[i]])))
            
          # Keep track of whether to fill R_max-independent histograms
          self.fill_Rmax_indep_hists = (i == 0)
          
          # Perform constituent subtraction on det-level, if applicable
          self.fill_background_histograms(fj_particles_combined_beforeCS, fj_particles_combined[i], jetR, i)
      
          # Do jet finding (re-do each time, to make sure matching info gets reset)
          cs_det = fj.ClusterSequence(fj_particles_det, jet_def)
          jets_det_pp = fj.sorted_by_pt(cs_det.inclusive_jets())
          jets_det_pp_selected = jet_selector_det(jets_det_pp)
          
          cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)
          jets_truth = fj.sorted_by_pt(cs_truth.inclusive_jets())
          jets_truth_selected = jet_selector_det(jets_truth)
          jets_truth_selected_matched = jet_selector_truth_matched(jets_truth)
          
          cs_combined = fj.ClusterSequence(fj_particles_combined[i], jet_def)
          jets_combined = fj.sorted_by_pt(cs_combined.inclusive_jets())
          jets_combined_selected = jet_selector_det(jets_combined)
          
          self.analyze_jets(jets_combined_selected, jets_truth_selected, jets_truth_selected_matched, jetR, jets_det_pp_selected = jets_det_pp_selected, R_max = R_max)

  #---------------------------------------------------------------
  # Analyze jets of a given event.
  #---------------------------------------------------------------
  def analyze_jets(self, jets_det_selected, jets_truth_selected, jets_truth_selected_matched, jetR, jets_det_pp_selected = None, R_max = None):
  
    if self.debug_level > 1:
      print('Number of det-level jets: {}'.format(len(jets_det_selected)))
    
    # Fill det-level jet histograms (before matching)
    for jet_det in jets_det_selected:
      
      # Check additional acceptance criteria
      # skip event if not satisfied -- since first jet in event is highest pt
      if not self.utils.is_det_jet_accepted(jet_det):
        if self.fill_R_indep_hists:
          self.hNevents.Fill(0)
        if self.debug_level > 1:
          print('event rejected due to jet acceptance')
        return
      
      self.fill_det_before_matching(jet_det, jetR, R_max)
  
    # Fill truth-level jet histograms (before matching)
    for jet_truth in jets_truth_selected:
    
      if self.is_pp or self.fill_Rmax_indep_hists:
        self.fill_truth_before_matching(jet_truth, jetR)
  
    # Loop through jets and set jet matching candidates for each jet in user_info
    if self.is_pp:
        [[self.set_matching_candidates(jet_det, jet_truth, jetR, 'hDeltaR_All_R{}'.format(jetR)) for jet_truth in jets_truth_selected_matched] for jet_det in jets_det_selected]
    else:
        # First fill the combined-to-pp matches, then the pp-to-pp matches
        [[self.set_matching_candidates(jet_det_combined, jet_det_pp, jetR, 'hDeltaR_combined_ppdet_R{{}}_Rmax{}'.format(R_max), fill_jet1_matches_only=True) for jet_det_pp in jets_det_pp_selected] for jet_det_combined in jets_det_selected]
        [[self.set_matching_candidates(jet_det_pp, jet_truth, jetR, 'hDeltaR_ppdet_pptrue_R{{}}_Rmax{}'.format(R_max)) for jet_truth in jets_truth_selected_matched] for jet_det_pp in jets_det_pp_selected]
        
    # Loop through jets and set accepted matches
    if self.is_pp:
        hname = 'hJetMatchingQA_R{}'.format(jetR)
        [self.set_matches_pp(jet_det, hname) for jet_det in jets_det_selected]
    else:
        hname = 'hJetMatchingQA_R{}_Rmax{}'.format(jetR, R_max)
        [self.set_matches_AA(jet_det_combined, jetR, hname) for jet_det_combined in jets_det_selected]
    
    # Loop through jets and fill matching histograms
    result = [self.fill_matching_histograms(jet_det, jetR, R_max) for jet_det in jets_det_selected]

    # Loop through jets and fill groomed histograms if both det and truth jets are unique match
    for grooming_setting in self.grooming_settings:
      if self.debug_level > 1:
        print('grooming setting: {}'.format(grooming_setting))
      result = [self.fill_groomed_jet_matches(grooming_setting, jet_det, jetR, R_max) for jet_det in jets_det_selected]

  #---------------------------------------------------------------
  # Fill some background histograms
  #---------------------------------------------------------------
  def fill_background_histograms(self, fj_particles_combined_beforeCS, fj_particles_combined, jetR, i):

    # Fill rho
    rho = self.constituent_subtractor[i].bge_rho.rho()
    if self.fill_R_indep_hists and self.fill_Rmax_indep_hists:
      getattr(self, 'hRho').Fill(rho)
    
    # Fill random cone delta-pt before constituent subtraction
    if not self.skip_deltapt_RC_histograms:
      R_max = self.max_distance[i]
      self.fill_deltapt_RC_histogram(fj_particles_combined_beforeCS, rho, jetR, R_max, before_CS=True)
          
      # Fill random cone delta-pt after constituent subtraction
      self.fill_deltapt_RC_histogram(fj_particles_combined, rho, jetR, R_max, before_CS=False)
    
  #---------------------------------------------------------------
  # Fill delta-pt histogram
  #---------------------------------------------------------------
  def fill_deltapt_RC_histogram(self, fj_particles, rho, jetR, R_max, before_CS=False):
  
    # Choose a random eta-phi in the fiducial acceptance
    phi = random.uniform(0., 2*np.pi)
    eta = random.uniform(-0.9+jetR, 0.9-jetR)
    
    # Loop through tracks and sum pt inside the cone
    pt_sum = 0.
    pt_sum_global = 0.
    for track in fj_particles:
        if self.utils.delta_R(track, eta, phi) < jetR:
            pt_sum += track.pt()
        pt_sum_global += track.pt()
            
    if before_CS:
        delta_pt = pt_sum - rho * np.pi * jetR * jetR
        getattr(self, 'hDeltaPt_RC_beforeCS_R{}_Rmax{}'.format(jetR, R_max)).Fill(delta_pt)
    else:
        delta_pt = pt_sum
        getattr(self, 'hDeltaPt_RC_afterCS_R{}_Rmax{}'.format(jetR, R_max)).Fill(delta_pt)
        
    # Fill mean pt
    if before_CS and self.fill_R_indep_hists and self.fill_Rmax_indep_hists:
      N_tracks = len(fj_particles)
      mean_pt = pt_sum_global/N_tracks
      getattr(self, 'hN_MeanPt').Fill(N_tracks, mean_pt)

  #---------------------------------------------------------------
  # Fill truth jet histograms
  #---------------------------------------------------------------
  def fill_truth_before_matching(self, jet, jetR):
    
    jet_pt = jet.pt()
    for constituent in jet.constituents():
      z = constituent.pt() / jet.pt()
      getattr(self, 'hZ_Truth_R{}'.format(jetR)).Fill(jet.pt(), z)
          
    # Fill 2D histogram of truth (pt, obs)
    hname = 'h_{{}}_JetPt_Truth_R{}_{{}}'.format(jetR)
    self.fill_observable_histograms(jet, jetR, hname)

  #---------------------------------------------------------------
  # Fill det jet histograms
  #---------------------------------------------------------------
  def fill_det_before_matching(self, jet, jetR, R_max):
    
    if self.is_pp or self.fill_Rmax_indep_hists:
      jet_pt = jet.pt()
      for constituent in jet.constituents():
        z = constituent.pt() / jet_pt
        getattr(self, 'hZ_Det_R{}'.format(jetR)).Fill(jet_pt, z)
      
    # Fill groomed histograms
    if self.thermal_model:
      hname = 'h_{{}}_JetPt_R{}_{{}}_Rmax{}'.format(jetR, R_max)
      self.fill_observable_histograms(jet, jetR, hname)

  #---------------------------------------------------------------
  # Fill 2D histogram of (pt, obs)
  #---------------------------------------------------------------
  def fill_observable_histograms(self, jet, jetR, hname):
  
    for grooming_setting in self.grooming_settings:
          
      grooming_label = self.utils.grooming_label(grooming_setting)
    
      # Construct SD groomer, and groom jet
      if 'sd' in grooming_setting:

        zcut = grooming_setting['sd'][0]
        beta = grooming_setting['sd'][1]
        
        # Note: Set custom recluster definition, since by default it uses jetR=max_allowable_R
        sd = fjcontrib.SoftDrop(beta, zcut, jetR)
        jet_def_recluster = fj.JetDefinition(fj.cambridge_algorithm, jetR)
        reclusterer = fjcontrib.Recluster(jet_def_recluster)
        sd.set_reclustering(True, reclusterer)
        if self.debug_level > 2:
          print('SoftDrop groomer is: {}'.format(sd.description()))

        jet_det_sd = sd.result(jet)

      # Construct Dynamical groomer, and groom jet
      if 'dg' in grooming_setting:

        a = grooming_setting['dg'][0]
        
        jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
        dy_groomer = fjcontrib.DynamicalGroomer(jet_def_lund)
        if self.debug_level > 2:
          print('Dynamical groomer is: {}'.format(dy_groomer.description()))
        
        jet_det_dg_lund = dy_groomer.result(jet, a)
        jet_det_dg = jet_det_dg_lund.pair()
         
      # Compute groomed observables
      if 'sd' in grooming_setting:

        # If both SD and DG are specified, first apply DG, then SD
        if 'dg' in grooming_setting:
          if jet_det_dg.has_constituents():
            jet_det_groomed = sd.result(jet_det_dg)
          else:
            return
        else:
          jet_det_groomed = jet_det_sd
        
        theta_g_det = self.theta_g(jet_det_groomed, jetR)
        zg_det = self.zg(jet_det_groomed)

      elif 'dg' in grooming_setting:
        jet_det_groomed = jet_det_dg

        theta_g_det = jet_det_dg_lund.Delta()/jetR
        zg_det = jet_det_dg_lund.z()

      # Fill histograms
      jet_pt = jet.pt()
      getattr(self, hname.format('theta_g', grooming_label)).Fill(jet_pt, theta_g_det)
      getattr(self,  hname.format('zg', grooming_label)).Fill(jet_pt, zg_det)
      
  #---------------------------------------------------------------
  # Loop through jets and fill matching histos
  #---------------------------------------------------------------
  def fill_matching_histograms(self, jet_det, jetR, R_max):
      
    if jet_det.has_user_info():
      jet_truth = jet_det.python_info().match
      
      if jet_truth:
        
        jet_pt_det_ungroomed = jet_det.pt()
        jet_pt_truth_ungroomed = jet_truth.pt()
        if self.is_pp or self.fill_Rmax_indep_hists:
          JES = (jet_pt_det_ungroomed - jet_pt_truth_ungroomed) / jet_pt_truth_ungroomed
          getattr(self, 'hJES_R{}'.format(jetR)).Fill(jet_pt_truth_ungroomed, JES)
        
        if not self.is_pp:
        
            # Get pp-det match, and fill delta-pt histogram
            jet_pp_det = jet_truth.python_info().match
            
            if jet_pp_det:
                jet_pp_det_pt = jet_pp_det.pt()
                delta_pt = (jet_pt_det_ungroomed - jet_pp_det_pt)
                getattr(self, 'hDeltaPt_emb_R{}_Rmax{}'.format(jetR, R_max)).Fill(jet_pt_truth_ungroomed, delta_pt)
                
        if 'subjet_z' in self.observable_list:
        
          if self.is_pp:
            result = [self.fill_subjet_matching_histograms(jet_det, jet_truth, None, jetR, subjetR) for subjetR in self.obs_settings['subjet_z']]
          else:
            result = [self.fill_subjet_matching_histograms(jet_det, jet_truth, jet_pp_det, jetR, subjetR) for subjetR in self.obs_settings['subjet_z']]
        
  #---------------------------------------------------------------
  # Loop through jets and fill matching histos
  #---------------------------------------------------------------
  def fill_subjet_matching_histograms(self, jet_det, jet_truth, jet_pp_det, jetR, subjetR):
  
    # Find subjets
    cs_subjet_det = fj.ClusterSequence(jet_det.constituents(), self.subjet_def[subjetR])
    subjets_det = fj.sorted_by_pt(cs_subjet_det.inclusive_jets())

    cs_subjet_truth = fj.ClusterSequence(jet_truth.constituents(), self.subjet_def[subjetR])
    subjets_truth = fj.sorted_by_pt(cs_subjet_truth.inclusive_jets())
    
    if not self.is_pp:
      cs_subjet_det_pp = fj.ClusterSequence(jet_pp_det.constituents(), self.subjet_def[subjetR])
      subjets_det_pp = fj.sorted_by_pt(cs_subjet_truth.inclusive_jets())

    # Loop through subjets and set subjet matching candidates for each subjet in user_info
    if self.is_pp:
        [[self.set_matching_candidates(subjet_det, subjet_truth, subjetR, 'hDeltaR_All_subjet_z_R{}_{}'.format(jetR, subjetR)) for subjet_truth in subjets_truth] for subjet_det in subjets_det]
    else:
        # First fill the combined-to-pp matches, then the pp-to-pp matches
        [[self.set_matching_candidates(subjet_det_combined, subjet_det_pp, subjetR, 'hDeltaR_combined_ppdet_subjet_z_R{}_{}'.format(jetR, subjetR), fill_jet1_matches_only=True) for subjet_det_pp in subjets_det_pp] for subjet_det_combined in subjets_det]
        [[self.set_matching_candidates(subjet_det_pp, subjet_truth, subjetR, 'hDeltaR_ppdet_pptrue_subjet_z_R{}_{}'.format(jetR, subjetR)) for subjet_truth in subjets_truth] for subjet_det_pp in subjets_det_pp]
      
    # Loop through subjets and set accepted matches
    if self.is_pp:
        [self.set_matches_pp(subjet_det, 'hSubjetMatchingQA_R{}_{}'.format(jetR, subjetR)) for subjet_det in subjets_det]
    else:
        [self.set_matches_AA(subjet_det_combined, subjetR, 'hSubjetMatchingQA_R{}_{}'.format(jetR, subjetR)) for subjet_det_combined in subjets_det]

    # Loop through matches and fill histograms
    for subjet_det in subjets_det:

      if subjet_det.has_user_info():
        subjet_truth = subjet_det.python_info().match
      
        if subjet_truth:
          
          z_det = subjet_det.pt() / jet_det.pt()
          z_truth = subjet_truth.pt() / jet_truth.pt()

          x = ([jet_det.pt(), jet_truth.pt(), z_det, z_truth])
          x_array = array('d', x)
          getattr(self, 'hResponse_JetPt_subjet_z_R{}_{}'.format(jetR, subjetR)).Fill(x_array)

          if z_truth > 1e-5:
            z_resolution = (z_det - z_truth) / z_truth
            getattr(self, 'hResidual_JetPt_{}_R{}_{}'.format('subjet_z', jetR, subjetR)).Fill(jet_truth.pt(), z_resolution)

  #---------------------------------------------------------------
  # Loop through jets and fill response if both det and truth jets are unique match
  #---------------------------------------------------------------
  def fill_groomed_jet_matches(self, grooming_setting, jet_det, jetR, R_max):
       
    grooming_label = self.utils.grooming_label(grooming_setting)
    
    if R_max:
      suffix = '_Rmax{}'.format(R_max)
    else:
      suffix = ''
       
    if jet_det.has_user_info():
      jet_truth = jet_det.python_info().match
     
      if jet_truth:
     
        jet_pt_det_ungroomed = jet_det.pt()
        jet_pt_truth_ungroomed = jet_truth.pt()
        
        if self.debug_level > 2:
          print('**** jet_pt_det_ungroomed: {}, jet_pt_truth_ungroomed: {}'.format(jet_pt_det_ungroomed, jet_pt_truth_ungroomed))
       
        # Construct SD groomer, and groom jet
        if 'sd' in grooming_setting:
         
          zcut = grooming_setting['sd'][0]
          beta = grooming_setting['sd'][1]
            
          # Note: Set custom recluster definition, since by default it uses jetR=max_allowable_R
          sd = fjcontrib.SoftDrop(beta, zcut, jetR)
          jet_def_recluster = fj.JetDefinition(fj.cambridge_algorithm, jetR)
          reclusterer = fjcontrib.Recluster(jet_def_recluster)
          sd.set_reclustering(True, reclusterer)
          if self.debug_level > 2:
            print('SoftDrop groomer is: {}'.format(sd.description()))
 
          jet_det_sd = sd.result(jet_det)
          jet_truth_sd = sd.result(jet_truth)

        # Construct Dynamical groomer, and groom jet
        if 'dg' in grooming_setting:

          a = grooming_setting['dg'][0]
          
          jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
          dy_groomer = fjcontrib.DynamicalGroomer(jet_def_lund)
          if self.debug_level > 2:
            print('Dynamical groomer is: {}'.format(dy_groomer.description()))
          
          jet_det_dg_lund = dy_groomer.result(jet_det, a)
          jet_truth_dg_lund = dy_groomer.result(jet_truth, a)
          
          jet_det_dg = jet_det_dg_lund.pair()
          jet_truth_dg = jet_truth_dg_lund.pair()
            
        # Compute groomed observables
        if 'sd' in grooming_setting:

          # If both SD and DG are specified, first apply DG, then SD
          if 'dg' in grooming_setting:
            if self.debug_level > 2:
              print('both SD and DG applied')
            if jet_det_dg.has_constituents() and jet_truth_dg.has_constituents():
              jet_det_groomed = sd.result(jet_det_dg)
              jet_truth_groomed = sd.result(jet_truth_dg)
            else:
              return
          else:
            if self.debug_level > 2:
              print('SD applied')
            jet_det_groomed = jet_det_sd
            jet_truth_groomed = jet_truth_sd
            
          theta_g_det = self.theta_g(jet_det_groomed, jetR)
          theta_g_truth = self.theta_g(jet_truth_groomed, jetR)
          zg_det = self.zg(jet_det_groomed)
          zg_truth = self.zg(jet_truth_groomed)

        elif 'dg' in grooming_setting:
        
          if self.debug_level > 2:
            print('DG applied')
          
          jet_det_groomed = jet_det_dg
          jet_truth_groomed = jet_truth_dg
          
          theta_g_det = jet_det_dg_lund.Delta()/jetR
          theta_g_truth = jet_truth_dg_lund.Delta()/jetR
          zg_det = jet_det_dg_lund.z()
          zg_truth = jet_truth_dg_lund.z()
        
        # Fill Lund planes
        if 'sd' in grooming_setting:

          lund_coords = self.lund_coordinates_SD(jet_truth_groomed)
          name = 'hLundPlane_R{}_{}{}'.format(jetR, grooming_label, suffix)
          if jet_pt_truth_ungroomed > 100.:
            getattr(self, name).Fill(lund_coords[0], lund_coords[1])
          
          if not self.is_pp:
            if 'dg' in grooming_setting:
              groomer_list = [sd, dy_groomer]
              if grooming_setting in self.obs_grooming_settings['theta_g']:
                prong_match = self.fill_prong_matching_histograms(jet_truth, jet_det, jet_det_groomed, groomer_list, jet_pt_truth_ungroomed, jetR, grooming_setting, grooming_label, R_max, type = 'SD+DG')
            else:
              groomer_list = [sd]
              if grooming_setting in self.obs_grooming_settings['theta_g']:
                prong_match = self.fill_prong_matching_histograms(jet_truth, jet_det, jet_det_groomed, groomer_list, jet_pt_truth_ungroomed, jetR, grooming_setting, grooming_label, R_max, type = 'SD')
      
        elif 'dg' in grooming_setting:

          lund_coords = self.lund_coordinates_DG(jet_truth_dg_lund)
          name = 'hLundPlane_R{}_{}{}'.format(jetR, grooming_label, suffix)
          if jet_pt_truth_ungroomed > 100.:
            getattr(self, name).Fill(lund_coords[0], lund_coords[1])
              
          if not self.is_pp and grooming_setting in self.obs_grooming_settings['theta_g']:
            prong_match = self.fill_prong_matching_histograms(jet_truth, jet_det, jet_det_dg_lund, [dy_groomer], jet_pt_truth_ungroomed, jetR, grooming_setting, grooming_label, R_max, type = 'DG')

        # If PbPb, fill extra RM only for successful prong matches
        if self.is_pp:
          prong_match = False
        
        # Fill histograms
        observable = 'theta_g'
        self.fill_groomed_response(observable, jetR, jet_pt_det_ungroomed, jet_pt_truth_ungroomed, theta_g_det, theta_g_truth, grooming_setting, self.obs_grooming_settings[observable], grooming_label, R_max, prong_match = prong_match)
          
        observable = 'zg'
        self.fill_groomed_response(observable, jetR, jet_pt_det_ungroomed, jet_pt_truth_ungroomed, zg_det, zg_truth, grooming_setting, self.obs_grooming_settings[observable], grooming_label, R_max, prong_match = prong_match)

        # Fill jet axis difference
        if 'jet_axis' in self.observable_list:

          # Recluster with WTA (with larger jet R)
          jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
          jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
          if self.debug_level > 3:
              print('WTA jet definition is:', jet_def_wta)
          reclusterer_wta =  fjcontrib.Recluster(jet_def_wta)
          jet_det_wta = reclusterer_wta.result(jet_det)
          jet_truth_wta = reclusterer_wta.result(jet_truth)
      
          self.fill_jet_axis_response(jet_det, jet_truth, jet_det_groomed, jet_truth_groomed, jet_det_wta, jet_truth_wta, jetR, grooming_setting, grooming_label)

  #---------------------------------------------------------------
  # Fill response histograms
  #---------------------------------------------------------------
  def fill_jet_axis_response(self, jet_det, jet_truth, jet_det_sd, jet_truth_sd, jet_det_wta, jet_truth_wta, jetR, grooming_setting, grooming_label):
  
    for axis in self.obs_settings['jet_axis']:

      if axis == 'Standard_SD':
        if grooming_setting in self.obs_grooming_settings['jet_axis']:

          deltaR_det = jet_det.delta_R(jet_det_sd)
          deltaR_truth = jet_truth.delta_R(jet_truth_sd)
          
          x = ([jet_det.pt(), jet_truth.pt(), deltaR_det, deltaR_truth])
          x_array = array('d', x)
          getattr(self, 'hResponse_JetPt_jet_axis_R{}_Standard_SD{}'.format(jetR, grooming_label)).Fill(x_array)
          
          if deltaR_truth > 1e-5:
            axis_resolution = (deltaR_det - deltaR_truth) / deltaR_truth
            getattr(self, 'hResidual_JetPt_{}_R{}_Standard_SD{}'.format('jet_axis', jetR, grooming_label)).Fill(jet_truth.pt(), axis_resolution)
          
      if axis == 'Standard_WTA':
        if grooming_setting == self.grooming_settings[0]:

          deltaR_det = jet_det.delta_R(jet_det_wta)
          deltaR_truth = jet_truth.delta_R(jet_truth_wta)
          
          x = ([jet_det.pt(), jet_truth.pt(), deltaR_det, deltaR_truth])
          x_array = array('d', x)
          getattr(self, 'hResponse_JetPt_jet_axis_R{}_Standard_WTA'.format(jetR)).Fill(x_array)
          
          if deltaR_truth > 1e-5:
            axis_resolution = (deltaR_det - deltaR_truth) / deltaR_truth
            getattr(self, 'hResidual_JetPt_{}_R{}_Standard_WTA'.format('jet_axis', jetR)).Fill(jet_truth.pt(), axis_resolution)
        
      if axis == 'WTA_SD':
        if grooming_setting in self.obs_grooming_settings['jet_axis']:

          deltaR_det = jet_det_sd.delta_R(jet_det_wta)
          deltaR_truth = jet_truth_sd.delta_R(jet_truth_wta)
          
          x = ([jet_det.pt(), jet_truth.pt(), deltaR_det, deltaR_truth])
          x_array = array('d', x)
          getattr(self, 'hResponse_JetPt_jet_axis_R{}_WTA_SD{}'.format(jetR, grooming_label)).Fill(x_array)
          
          if deltaR_truth > 1e-5:
            axis_resolution = (deltaR_det - deltaR_truth) / deltaR_truth
            getattr(self, 'hResidual_JetPt_{}_R{}_WTA_SD{}'.format('jet_axis', jetR, grooming_label)).Fill(jet_truth.pt(), axis_resolution)
    
  #---------------------------------------------------------------
  # Fill response histograms
  #---------------------------------------------------------------
  def fill_groomed_response(self, observable, jetR, jet_pt_det_ungroomed, jet_pt_truth_ungroomed, obs_det, obs_truth, grooming_setting, obs_grooming_settings, grooming_label, R_max, prong_match = False):
  
    if grooming_setting in obs_grooming_settings:

      x = ([jet_pt_det_ungroomed, jet_pt_truth_ungroomed, obs_det, obs_truth])
      x_array = array('d', x)
      name = 'hResponse_JetPt_{}_R{}_{}'.format(observable, jetR, grooming_label)
      if not self.is_pp:
        name += '_Rmax{}'.format(R_max)
      getattr(self, name).Fill(x_array)
        
      if obs_truth > 1e-5:
        obs_resolution = (obs_det - obs_truth) / obs_truth
        name = 'hResidual_JetPt_{}_R{}_{}'.format(observable, jetR, grooming_label)
        if not self.is_pp:
          name += '_Rmax{}'.format(R_max)
        getattr(self, name).Fill(jet_pt_truth_ungroomed, obs_truth, obs_resolution)
      
      # Fill prong-matched response
      if not self.is_pp and R_max == self.main_R_max:
        if prong_match:
        
          name = 'hResponse_JetPt_{}_R{}_{}_Rmax{}_matched'.format(observable, jetR, grooming_label, R_max)
          getattr(self, name).Fill(x_array)
          
          if obs_truth > 1e-5:
            name = 'hResidual_JetPt_{}_R{}_{}_Rmax{}_matched'.format(observable, jetR, grooming_label, R_max)
            getattr(self, name).Fill(jet_pt_truth_ungroomed, obs_truth, obs_resolution)
          
  #---------------------------------------------------------------
  # Do prong-matching
  #---------------------------------------------------------------
  def fill_prong_matching_histograms(self, jet_truth, jet_det, jet_det_groomed, groomer_list, jet_pt_truth_ungroomed, jetR, grooming_setting, grooming_label, R_max, type = 'SD'):
    
    # Do grooming on pp-det jet, and get prongs
    jet_pp_det = jet_truth.python_info().match

    if 'SD' in type:
    
      if 'DG' in type:
      
        # Assumes groomer_list = [sd, dy_groomer]
        a = grooming_setting['dg'][0]
        jet_pp_det_dg_lund = groomer_list[1].result(jet_pp_det, a)
        jet_pp_det_dg = jet_pp_det_dg_lund.pair()
        jet_pp_det_groomed = groomer_list[0].result(jet_pp_det_dg)
      
      else:
      
        jet_pp_det_groomed = groomer_list[0].result(jet_pp_det)
    
      # SD grooming returns the groomed fastjet::PseudoJet
      # Use the fastjet::PseudoJet::has_parents function which returns the last clustering step
      #   If the jet passed SoftDrop, then its parents are the SoftDrop splitting
      #   If the jet didn't pass SoftDrop, then it will have no parents
      jet_pp_det_prong1 = fj.PseudoJet()
      jet_pp_det_prong2 = fj.PseudoJet()
      has_parents_pp_det = jet_pp_det_groomed.has_parents(jet_pp_det_prong1, jet_pp_det_prong2)
      
      # Get prongs of combined jet
      jet_combined_prong1 = fj.PseudoJet()
      jet_combined_prong2 = fj.PseudoJet()
      has_parents_combined = jet_det_groomed.has_parents(jet_combined_prong1, jet_combined_prong2)

    elif type == 'DG':
      
      a = grooming_setting['dg'][0]
      jet_pp_det_groomed_lund = groomer_list[0].result(jet_pp_det, a)
    
      # Dynamical grooming returns a fjcontrib::LundGenerator
      #   The prongs can be retrieved directly from this object.
      #   If the object exists, then it has passed grooming
      jet_pp_det_prong1 = jet_pp_det_groomed_lund.harder()
      jet_pp_det_prong2 = jet_pp_det_groomed_lund.softer()
      has_parents_pp_det = jet_pp_det_groomed_lund

      # Get prongs of combined jet
      jet_combined_prong1 = jet_det_groomed.harder()
      jet_combined_prong2 = jet_det_groomed.softer()
      has_parents_combined = jet_det_groomed
      
      # Get the fastjet::PseudoJets from the fjcontrib::LundGenerators
      jet_pp_det_groomed = jet_pp_det_groomed_lund.pair()
      jet_det_groomed = jet_det_groomed.pair()
          
    if self.debug_level > 1:

        if jet_pt_truth_ungroomed > 80.:
        
            print('=======================================================')
            print(type)
            print('jet_pt_truth_ungroomed: {}'.format(jet_pt_truth_ungroomed))
            print('jet_pt_pp_det_ungroomed: {}'.format(jet_pp_det.pt()))
            print('jet_pt_pp_det_groomed: {}'.format(jet_pp_det_groomed.pt()))
            print('jet_pt_combined_groomed: {}'.format(jet_det_groomed.pt()))
            print('')
            print('jet_pp_det tracks: {}'.format([track.user_index() for track in jet_pp_det.constituents()]))
            print('         track pt: {}'.format([np.around(track.pt(),2) for track in jet_pp_det.constituents()]))
            print('jet_pp_det_groomed tracks: {}'.format([track.user_index() for track in jet_pp_det_groomed.constituents()]))
            print('                 track pt: {}'.format([np.around(track.pt(),2) for track in jet_pp_det_groomed.constituents()]))
            print('jet_combined groomed tracks: {}'.format([track.user_index() for track in jet_det_groomed.constituents()]))
            print('                   track pt: {}'.format([np.around(track.pt(),2) for track in jet_det_groomed.constituents()]))
            print('jet_combined ungroomed tracks: {}'.format([track.user_index() for track in jet_det.constituents()]))
            print('                     track pt: {}'.format([np.around(track.pt(),2) for track in jet_det.constituents()]))

    # Compute fraction of pt of the pp-det prong tracks that is contained in the combined-jet prong,
    # in order to have a measure of whether the combined-jet prong is the "same" prong as the pp-det prong
    deltaR_prong1 = -1.
    deltaR_prong2 = -1.
    deltaZ = -1.
    if has_parents_pp_det and has_parents_combined:
    
        # Subleading jet pt-matching
        # --------------------------
        # (1) Fraction of pt matched: subleading pp-det in subleading combined
        matched_pt_subleading_subleading = fjtools.matched_pt(jet_combined_prong2, jet_pp_det_prong2)
        
        # (2) Fraction of pt matched: subleading pp-det in leading combined
        matched_pt_subleading_leading = fjtools.matched_pt(jet_combined_prong1, jet_pp_det_prong2)
        
        # (3) Fraction of pt matched: subleading pp-det in ungroomed combined jet
        matched_pt_subleading_groomed = fjtools.matched_pt(jet_det_groomed, jet_pp_det_prong2)
        matched_pt_subleading_ungroomed = fjtools.matched_pt(jet_det, jet_pp_det_prong2)
        matched_pt_subleading_ungroomed_notgroomed = matched_pt_subleading_ungroomed - matched_pt_subleading_groomed
        
        # (4) Fraction of pt matched: subleading pp-det not in ungroomed combined jet
        matched_pt_subleading_outside = 1 - matched_pt_subleading_ungroomed

        # Leading jet pt-matching
        # --------------------------
        # (1) Fraction of pt matched: leading pp-det in subleading combined
        matched_pt_leading_subleading = fjtools.matched_pt(jet_combined_prong2, jet_pp_det_prong1)
        
        # (2) Fraction of pt matched: leading pp-det in leading combined
        matched_pt_leading_leading = fjtools.matched_pt(jet_combined_prong1, jet_pp_det_prong1)

        # (3) Fraction of pt matched: leading pp-det in ungroomed combined jet
        matched_pt_leading_groomed = fjtools.matched_pt(jet_det_groomed, jet_pp_det_prong1)
        matched_pt_leading_ungroomed = fjtools.matched_pt(jet_det, jet_pp_det_prong1)
        matched_pt_leading_ungroomed_notgroomed = matched_pt_leading_ungroomed - matched_pt_leading_groomed
        
        # (4) Fraction of pt matched: leading pp-det not in ungroomed combined jet
        matched_pt_leading_outside = 1 - matched_pt_leading_ungroomed

        # Compute delta-R between pp-det prong and combined prong
        # --------------------------
        deltaR_prong1 = jet_combined_prong1.delta_R(jet_pp_det_prong1)
        deltaR_prong2 = jet_combined_prong2.delta_R(jet_pp_det_prong2)
        deltaZ = self.zg(jet_det_groomed) - self.zg(jet_pp_det_groomed)
        
        if self.debug_level > 1:
        
            if jet_pt_truth_ungroomed > 80.:
            
                print('subleading prong tracks -- combined: {}'.format([track.user_index() for track in jet_combined_prong2.constituents()]))
                print('subleading prong tracks -- pp-det: {}'.format([track.user_index() for track in jet_pp_det_prong2.constituents()]))
                print('leading prong tracks -- combined: {}'.format([track.user_index() for track in jet_combined_prong1.constituents()]))
                print('leading prong tracks -- pp-det: {}'.format([track.user_index() for track in jet_pp_det_prong1.constituents()]))
                print('')
                print('leading_prong_pt: {}'.format(jet_combined_prong1.pt()))
                print('matched_pt_leading_subleading fraction: {}'.format(matched_pt_leading_subleading))
                print('matched_pt_leading_leading fraction: {}'.format(matched_pt_leading_leading))
                print('matched_pt_leading_ungroomed_notgroomed fraction: {}'.format(matched_pt_leading_ungroomed_notgroomed))
                print('matched_pt_leading_outside fraction: {}'.format(matched_pt_leading_outside))
                print('')
                print('subleading_prong_pt: {}'.format(jet_combined_prong2.pt()))
                print('matched_pt_subleading_subleading fraction: {}'.format(matched_pt_subleading_subleading))
                print('matched_pt_subleading_leading fraction: {}'.format(matched_pt_subleading_leading))
                print('matched_pt_subleading_ungroomed_notgroomed fraction: {}'.format(matched_pt_subleading_ungroomed_notgroomed))
                print('matched_pt_subleading_outside fraction: {}'.format(matched_pt_subleading_outside))
                print('')
                print('deltaR_prong1: {}'.format(deltaR_prong1))
                print('deltaR_prong2: {}'.format(deltaR_prong2))

    elif has_parents_pp_det: # pp-det passed grooming, but combined jet failed grooming
        matched_pt_leading_leading = matched_pt_leading_subleading = matched_pt_leading_ungroomed_notgroomed = matched_pt_leading_outside = matched_pt_subleading_leading = matched_pt_subleading_subleading = matched_pt_subleading_ungroomed_notgroomed = matched_pt_subleading_outside = -0.1
        
    elif has_parents_combined: # combined jet passed grooming, but pp-det failed grooming
        matched_pt_leading_leading = matched_pt_leading_subleading = matched_pt_leading_ungroomed_notgroomed = matched_pt_leading_outside = matched_pt_subleading_leading = matched_pt_subleading_subleading = matched_pt_subleading_ungroomed_notgroomed = matched_pt_subleading_outside = -0.2
        
    else: # both pp-det and combined jet failed SoftDrop
        matched_pt_leading_leading = matched_pt_leading_subleading = matched_pt_leading_ungroomed_notgroomed = matched_pt_leading_outside = matched_pt_subleading_leading = matched_pt_subleading_subleading = matched_pt_subleading_ungroomed_notgroomed = matched_pt_subleading_outside = -0.3

    # Leading prong
    getattr(self, 'hProngMatching_leading_leading_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_leading, deltaR_prong1)
    getattr(self, 'hProngMatching_leading_subleading_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_subleading, deltaR_prong1)
    getattr(self, 'hProngMatching_leading_ungroomed_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_ungroomed_notgroomed, deltaR_prong1)
    getattr(self, 'hProngMatching_leading_outside_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_outside, deltaR_prong1)
    
    getattr(self, 'hProngMatching_leading_leading_JetPtDet_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pp_det.pt(), matched_pt_leading_leading, deltaR_prong1)
    getattr(self, 'hProngMatching_leading_subleading_JetPtDet_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pp_det.pt(), matched_pt_leading_subleading, deltaR_prong1)
    getattr(self, 'hProngMatching_leading_ungroomed_JetPtDet_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pp_det.pt(), matched_pt_leading_ungroomed_notgroomed, deltaR_prong1)
    getattr(self, 'hProngMatching_leading_outside_JetPtDet_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pp_det.pt(), matched_pt_leading_outside, deltaR_prong1)
    
    getattr(self, 'hProngMatching_leading_leading_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_leading, deltaZ)
    getattr(self, 'hProngMatching_leading_subleading_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_subleading, deltaZ)
    getattr(self, 'hProngMatching_leading_ungroomed_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_ungroomed_notgroomed, deltaZ)
    getattr(self, 'hProngMatching_leading_outside_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_outside, deltaZ)

    # Subleading prong
    getattr(self, 'hProngMatching_subleading_leading_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_leading, deltaR_prong2)
    getattr(self, 'hProngMatching_subleading_subleading_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_subleading, deltaR_prong2)
    getattr(self, 'hProngMatching_subleading_ungroomed_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_ungroomed_notgroomed, deltaR_prong2)
    getattr(self, 'hProngMatching_subleading_outside_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_outside, deltaR_prong2)

    getattr(self, 'hProngMatching_subleading_leading_JetPtDet_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pp_det.pt(), matched_pt_subleading_leading, deltaR_prong2)
    getattr(self, 'hProngMatching_subleading_subleading_JetPtDet_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pp_det.pt(), matched_pt_subleading_subleading, deltaR_prong2)
    getattr(self, 'hProngMatching_subleading_ungroomed_JetPtDet_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pp_det.pt(), matched_pt_subleading_ungroomed_notgroomed, deltaR_prong2)
    getattr(self, 'hProngMatching_subleading_outside_JetPtDet_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pp_det.pt(), matched_pt_subleading_outside, deltaR_prong2)

    getattr(self, 'hProngMatching_subleading_leading_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_leading, deltaZ)
    getattr(self, 'hProngMatching_subleading_subleading_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_subleading, deltaZ)
    getattr(self, 'hProngMatching_subleading_ungroomed_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_ungroomed_notgroomed, deltaZ)
    getattr(self, 'hProngMatching_subleading_outside_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_outside, deltaZ)
    
    # Plot correlation of matched pt fraction for leading-subleading and subleading-leading
    getattr(self, 'hProngMatching_subleading-leading_correlation_JetPtDet_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pp_det.pt(), matched_pt_leading_subleading, matched_pt_subleading_leading)
    
    subleading_match = (matched_pt_subleading_subleading > 0.5)
    leading_match = (matched_pt_leading_leading > 0.5)
    prong_match = subleading_match and leading_match
    return prong_match
        
  #---------------------------------------------------------------
  # Compute theta_g
  #---------------------------------------------------------------
  def theta_g(self, jet_sd, jetR):
    
    sd_info = fjcontrib.get_SD_jet_info(jet_sd)
    theta_g = sd_info.dR / jetR
    return theta_g

  #---------------------------------------------------------------
  # Compute z_g
  #---------------------------------------------------------------
  def zg(self, jet_sd):
    
    sd_info = fjcontrib.get_SD_jet_info(jet_sd)
    zg = sd_info.z
    return zg

##################################################################
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Process MC')
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

  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
    sys.exit(0)
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = ProcessMC(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_mc()
