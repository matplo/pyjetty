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
import time

# Data analysis and plotting
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns
from array import *
import ROOT
import yaml

# Fastjet via python (from heppy)
import fastjet as fj
import fjcontrib
import fjtools

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_io_emb
from pyjetty.alice_analysis.process.base import jet_info
from pyjetty.alice_analysis.process.base import process_base
from pyjetty.alice_analysis.process.base import thermal_generator
from pyjetty.mputils import CEventSubtractor
from pyjetty.mputils import RTreeWriter

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class ProcessGroomers(process_base.ProcessBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
    super(ProcessGroomers, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
    
    # Initialize configuration
    self.initialize_config()
  
  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def process_groomers(self):
    
    self.start_time = time.time()
    
    # ------------------------------------------------------------------------
    # Use IO helper class to convert truth-level ROOT TTree into
    # a SeriesGroupBy object of fastjet particles per event
    tree_dir = 'PWGHF_TreeCreator'
    io_truth = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
                                    track_tree_name='tree_Particle_gen', use_ev_id_ext=self.use_ev_id_ext)
    self.df_fjparticles = io_truth.load_data()
    self.nEvents = len(self.df_fjparticles.index)
    self.nTracks = len(io_truth.track_df.index)
    print('--- {} seconds ---'.format(time.time() - self.start_time))
    
    # ------------------------------------------------------------------------
    
    # Set up the Pb-Pb embedding object
    if not self.thermal_model:
      self.process_io_emb = process_io_emb.ProcessIO_Emb(self.emb_file_list,
                              track_tree_name='tree_Particle_gen', is_pp = True,
                              use_ev_id_ext = self.use_ev_id_ext, remove_used_file=False)
        
    # ------------------------------------------------------------------------

    # Initialize histograms
    self.initialize_output_objects()
    
    # Create constituent subtractor, if configured
    if self.do_constituent_subtraction:
      self.constituent_subtractor = [CEventSubtractor(max_distance=R_max, alpha=self.alpha, max_eta=self.eta_max, bge_rho_grid_size=self.bge_rho_grid_size, max_pt_correct=self.max_pt_correct, ghost_area=self.ghost_area, distance_type=fjcontrib.ConstituentSubtractor.deltaR) for R_max in self.max_distance]
    
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
      
    self.jet_matching_distance = config['jet_matching_distance']
    self.mc_fraction_threshold = config['mc_fraction_threshold']
    self.prong_matching_threshold = config['prong_matching_threshold']
    self.use_ev_id_ext = config['use_ev_id_ext']
    self.main_R_max = config['constituent_subtractor']['main_R_max']
    self.eta_max = config['eta_max']
    self.plot_diagram =  config['plot_diagram']
    
    if 'thermal_model' in config:
      self.min_background_multiplicity = None
      self.thermal_model = True
      beta = config['thermal_model']['beta']
      N_avg = config['thermal_model']['N_avg']
      sigma_N = config['thermal_model']['sigma_N']
      self.thermal_generator = thermal_generator.ThermalGenerator(N_avg = N_avg, sigma_N = sigma_N,
                                                                  beta = beta, eta_max=self.eta_max)
    else:
      self.thermal_model = False
      self.min_background_multiplicity = config['angantyr']['min_background_multiplicity']
      self.emb_file_list = config['angantyr']['emb_file_list']

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
      obs_config_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
      obs_subconfig_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
      
      self.obs_settings[observable] = self.utils.obs_settings(observable, obs_config_dict, obs_subconfig_list)
      if observable == 'subjet_z':
        self.subjet_def = {}
        for subjetR in self.obs_settings[observable]:
          self.subjet_def[subjetR] = fj.JetDefinition(fj.antikt_algorithm, subjetR)
     
      # Fill grooming settings
      self.obs_grooming_settings[observable] = self.utils.grooming_settings(obs_config_dict)

    # Construct set of unique grooming settings
    self.grooming_settings = []
    lists_grooming = [self.obs_grooming_settings[obs] for obs in self.observable_list]
    for observable in lists_grooming:
      for setting in observable:
        if setting not in self.grooming_settings and setting != None:
          self.grooming_settings.append(setting)
          
    # Set reclustering algorithm
    self.recluster_alg = config['reclustering_algorithm']
    if self.recluster_alg == 'CA':
      self.reclustering_algorithm = fj.cambridge_algorithm
    elif self.recluster_alg == 'KT':
      self.reclustering_algorithm = fj.kt_algorithm
    elif self.recluster_alg == 'AKT':
      self.reclustering_algorithm = fj.antikt_algorithm

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_output_objects(self):
    
    self.hNevents = ROOT.TH1F('hNevents', 'hNevents', 2, -0.5, 1.5)
    if self.event_number_max < self.nEvents:
      self.hNevents.Fill(1, self.event_number_max)
    else:
      self.hNevents.Fill(1, self.nEvents)

    self.hRho =  ROOT.TH1F('hRho', 'hRho', 1000, 0., 1000.)
    self.hMult =  ROOT.TH1F('hMult', 'hMult', 100, 0., 20000.)

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_output_objects_R(self, jetR):
            
      for R_max in self.max_distance:
      
        name = 'hDeltaPt_emb_R{}_Rmax{}'.format(jetR, R_max)
        h = ROOT.TH2F(name, name, 300, 0, 300, 400, -200., 200.)
        setattr(self, name, h)

        name = 'hDeltaR_R{}_Rmax{}'.format(jetR, R_max)
        h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
        setattr(self, name, h)
        
        if 'subjet_z' in self.observable_list:
          for subjetR in self.obs_settings['subjet_z']:
        
            name = 'hDeltaR_{}_R{}_{}_Rmax{}'.format('subjet_z', jetR, subjetR, R_max)
            h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
            setattr(self, name, h)
          
      #  Construct THn for each groomer: (pt, zg, theta_g, tag flag)
      for grooming_setting in self.obs_grooming_settings['theta_g']:
        if grooming_setting:
          grooming_label = self.utils.grooming_label(grooming_setting)
          
          for R_max in self.max_distance:
            
            # THn for combined jets
            dim = 4;
            title = ['p_{T,ch jet}', '#it{z}_{g,ch}', '#theta_{g,ch}', 'flag']
            nbins = [30, 50, 100, 9]
            min = [0., 0., 0., 0.5]
            max = [300., 0.5, 1., 9.5]
            name = 'h_theta_g_zg_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)
            self.create_thn(name, title, dim, nbins, min, max)
          
          # TH3 for truth jets
          name = 'h_theta_g_zg_JetPt_Truth_R{}_{}'.format(jetR, grooming_label)
          h = ROOT.TH3F(name, name, 30, 0, 300, 50, 0, 0.5, 100, 0, 1.0)
          h.GetXaxis().SetTitle('p_{T,ch jet}')
          h.GetYaxis().SetTitle('#it{z}_{g,ch}')
          h.GetZaxis().SetTitle('#theta_{g,ch}')
          setattr(self, name, h)
            
      for observable in self.observable_list:

        if observable in ['kappa', 'tf']:
        
          if observable == 'kappa':
            label = '#kappa_{ch}'
          if observable == 'kappa':
            label = '#it{t}_{f}'
        
          for grooming_setting in self.obs_grooming_settings[observable]:
            if grooming_setting:
              grooming_label = self.utils.grooming_label(grooming_setting)

              # TH3 for combined jets
              for R_max in self.max_distance:
                name = 'h_{}_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, grooming_label, R_max)
                h = ROOT.TH3F(name, name, 30, 0, 300, 50, 0, 0.5, 9, 0.5, 9.5)
                h.GetXaxis().SetTitle('p_{T,ch jet}')
                h.GetYaxis().SetTitle(label)
                h.GetZaxis().SetTitle('flag')
                setattr(self, name, h)
                 
              # TH2 for truth jets
              name = 'h_{}_JetPt_Truth_R{}_{}'.format(observable, jetR, grooming_label)
              h = ROOT.TH2F(name, name, 30, 0, 300, 50, 0, 0.5)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle(label)
              setattr(self, name, h)
          
        if observable == 'subjet_z':
          
          for obs_setting in self.obs_settings[observable]:
        
            for R_max in self.max_distance:
            
              label = '#it{z}'
              name = 'h_{}_fraction_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, obs_setting, R_max)
              h = ROOT.TH3F(name, name, 30, 0, 300, 100, 0, 1.0, 15, -0.4, 1.1)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle(label)
              h.GetZaxis().SetTitle('Prong matching fraction')
              setattr(self, name, h)
              
              name = 'h_{}_flag_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, obs_setting, R_max)
              h = ROOT.TH3F(name, name, 30, 0, 300, 100, 0, 1.0, 9, 0.5, 9.5)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle(label)
              h.GetZaxis().SetTitle('flag')
              setattr(self, name, h)
              
              label = '#it{z}_{leading}'
              name = 'h_{}_fraction_leading_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, obs_setting, R_max)
              h = ROOT.TH3F(name, name, 30, 0, 300, 100, 0, 1.0, 15, -0.4, 1.1)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle(label)
              h.GetZaxis().SetTitle('Prong matching fraction')
              setattr(self, name, h)
              
              name = 'h_{}_flag_leading_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, obs_setting, R_max)
              h = ROOT.TH3F(name, name, 30, 0, 300, 100, 0, 1.0, 9, 0.5, 9.5)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle(label)
              h.GetZaxis().SetTitle('flag')
              setattr(self, name, h)
              
      # Create prong matching histograms
      for grooming_setting in self.obs_grooming_settings['theta_g']:
        if grooming_setting:
          grooming_label = self.utils.grooming_label(grooming_setting)
          self.create_prong_matching_histograms(jetR, grooming_label)
            
      # Create tree to store splitting info for all groomers
      self.fill_tree = False
      if self.fill_tree:
        self.t = ROOT.TTree('t', 't')
        self.tw = RTreeWriter(tree=self.t)
                
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
          h = ROOT.TH3F(name, name, 30, 0, 300, 15, -0.4, 1.1, 100, 0., 1)
          h.GetXaxis().SetTitle('p_{T,truth}')
          h.GetYaxis().SetTitle('Prong matching fraction')
          h.GetZaxis().SetTitle('#Delta R_{prong}/R_{g}')
          setattr(self, name, h)
          
          name = 'hProngMatching_{}_{}_JetPtZ_R{}_{}_Rmax{}'.format(prong, match, jetR, grooming_label, R_max)
          h = ROOT.TH3F(name, name, 30, 0, 300, 15, -0.4, 1.1, 200, -0.5, 0.5)
          h.GetXaxis().SetTitle('p_{T,truth}')
          h.GetYaxis().SetTitle('Prong matching fraction')
          h.GetZaxis().SetTitle('#Delta z_{prong}/z')
          setattr(self, name, h)

      name = 'hProngMatching_subleading-leading_correlation_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)
      h = ROOT.TH3F(name, name, 30, 0, 300, 15, -0.4, 1.1, 15, -0.4, 1.1)
      h.GetXaxis().SetTitle('p_{T,truth}')
      h.GetYaxis().SetTitle('Prong matching fraction, leading_subleading')
      h.GetZaxis().SetTitle('Prong matching fraction, subleading_leading')
      setattr(self, name, h)
      
      name = 'hProngMatching_truth_groomed_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)
      h = ROOT.TH2F(name, name, 30, 0, 300, 15, -0.4, 1.1)
      h.GetXaxis().SetTitle('p_{T,truth}')
      h.GetYaxis().SetTitle('Prong matching fraction, truth_groomed')
      setattr(self, name, h)
      
      name = 'hProngMatching_truth2_groomed_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)
      h = ROOT.TH2F(name, name, 30, 0, 300, 15, -0.4, 1.1)
      h.GetXaxis().SetTitle('p_{T,truth}')
      h.GetYaxis().SetTitle('Prong matching fraction, truth2_groomed')
      setattr(self, name, h)
      
  #---------------------------------------------------------------
  # Main function to loop through and analyze events
  #---------------------------------------------------------------
  def analyzeEvents(self):
    
    fj.ClusterSequence.print_banner()
    print()
    self.event_number = 0
    
    for jetR in self.jetR_list:
      self.initialize_output_objects_R(jetR)
    
    # Then can use list comprehension to iterate over the groupby and do jet-finding
    # simultaneously for fj_1 and fj_2 per event, so that I can match jets -- and fill histograms
    result = [self.analyze_event(fj_particles) for fj_particles in self.df_fjparticles]

    print('Save thn...')
    process_base.ProcessBase.save_thn_th3_objects(self)
    
  #---------------------------------------------------------------
  # Analyze jets of a given event.
  # fj_particles is the list of fastjet pseudojets for a single fixed event.
  #---------------------------------------------------------------
  def analyze_event(self, fj_particles_truth):
  
    self.event_number += 1
    if self.event_number > self.event_number_max:
      return
    if self.debug_level > 1:
      print('-------------------------------------------------')
      print('event {}'.format(self.event_number))
      
    # Check that the entries exist appropriately
    # (need to check how this can happen -- but it is only a tiny fraction of events)
    if type(fj_particles_truth) != fj.vectorPJ:
      print('fj_particles type mismatch -- skipping event')
      return
      
    if len(fj_particles_truth) > 1:
      if np.abs(fj_particles_truth[0].pt() - fj_particles_truth[1].pt()) <  1e-5:
        print('WARNING: Duplicate particles may be present')
        print([p.user_index() for p in fj_particles_truth])
        print([p.pt() for p in fj_particles_truth])
                  
    # Clear event in tree
    if self.fill_tree:
      self.tw.clear()
  
    # If Pb-Pb, construct embedded event (do this once, for all jetR)
        
    # If thermal model, generate a thermal event
    if self.thermal_model:
      fj_particles_combined_beforeCS = self.thermal_generator.load_event()

    # Angantyr: Get Pb-Pb event
    else:
      accept_background_event = False
      while not accept_background_event:
        fj_particles_combined_beforeCS = self.process_io_emb.load_event()
        particles = [p.eta() for p in fj_particles_combined_beforeCS]
        multiplicity = sum(np.abs(i) < 0.9 for i in particles)
        if multiplicity > self.min_background_multiplicity:
          accept_background_event = True
        self.hMult.Fill(multiplicity)
        if self.debug_level > 3:
          print('multiplicity: {}'.format(multiplicity))
          print('accepted: {}'.format(accept_background_event))
          
    # Form the combined event
    # The pp-truth tracks are each stored with a unique user_index >= 0
    #   (same index in fj_particles_combined and fj_particles_truth -- which will be used in prong-matching)
    # The Pb-Pb tracks are each stored with a unique user_index < 0
    [fj_particles_combined_beforeCS.push_back(p) for p in fj_particles_truth]
     
    # Perform constituent subtraction for each R_max
    fj_particles_combined = [self.constituent_subtractor[i].process_event(fj_particles_combined_beforeCS) for i, R_max in enumerate(self.max_distance)]
    
    if self.debug_level > 3:
      print([p.user_index() for p in fj_particles_truth])
      print([p.pt() for p in fj_particles_truth])
      print([p.user_index() for p in fj_particles_combined[0]])
      print([p.pt() for p in fj_particles_combined[0]])
      print([p.user_index() for p in fj_particles_combined_beforeCS])
      print([p.pt() for p in fj_particles_combined_beforeCS])
      
    # Loop through jetR, and process event for each R
    for jetR in self.jetR_list:
    
      # Keep track of whether to fill R-independent histograms
      self.fill_R_indep_hists = (jetR == self.jetR_list[0])

      # Set jet definition and a jet selector
      jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
      jet_selector_det = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(self.eta_max - jetR)
      jet_selector_truth_matched = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(self.eta_max)
      if self.debug_level > 2:
        print('')
        print('jet definition is:', jet_def)
        print('jet selector for det-level is:', jet_selector_det)
        print('jet selector for truth-level matches is:', jet_selector_truth_matched)
      
      # Analyze
      for i, R_max in enumerate(self.max_distance):
          
        if self.debug_level > 1:
          print('')
          print('R_max: {}'.format(R_max))
          print('Total number of combined particles: {}'.format(len([p.pt() for p in fj_particles_combined_beforeCS])))
          print('After constituent subtraction {}: {}'.format(i, len([p.pt() for p in fj_particles_combined[i]])))
          
        # Keep track of whether to fill R_max-independent histograms
        self.fill_Rmax_indep_hists = (i == 0)
        
        # Perform constituent subtraction on det-level, if applicable
        rho = self.constituent_subtractor[i].bge_rho.rho()
        if self.fill_R_indep_hists and self.fill_Rmax_indep_hists:
          getattr(self, 'hRho').Fill(rho)
    
        # Do jet finding (re-do each time, to make sure matching info gets reset)
        cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)
        jets_truth = fj.sorted_by_pt(cs_truth.inclusive_jets())
        jets_truth_selected = jet_selector_det(jets_truth)
        jets_truth_selected_matched = jet_selector_truth_matched(jets_truth)
        
        cs_combined = fj.ClusterSequence(fj_particles_combined[i], jet_def)
        jets_combined = fj.sorted_by_pt(cs_combined.inclusive_jets())
        jets_combined_selected = jet_selector_det(jets_combined)
        
        self.analyze_jets(jets_combined_selected, jets_truth_selected, jets_truth_selected_matched,
                          jetR, R_max = R_max)
    
    if self.fill_tree:
      self.tw.fill_tree()

  #---------------------------------------------------------------
  # Analyze jets of a given event.
  #---------------------------------------------------------------
  def analyze_jets(self, jets_combined_selected, jets_truth_selected,
                   jets_truth_selected_matched, jetR, R_max = None):
  
    if self.debug_level > 1:
      print('Number of det-level jets: {}'.format(len(jets_combined_selected)))
      
    # Loop through jets and set jet matching candidates (based on deltaR) for each jet in user_info
    [[self.set_matching_candidates(jet_combined, jet_truth, jetR,
     'hDeltaR_R{{}}_Rmax{}'.format(R_max)) for jet_truth in jets_truth_selected]
                                           for jet_combined in jets_combined_selected]
        
    # Loop through jets and set accepted matches
    hname = 'hJetMatchingQA_R{}_Rmax{}'.format(jetR, R_max)
    [self.set_matches_AA_truth(jet_combined, jetR, hname) for jet_combined in jets_combined_selected]
    
    # Loop through jets and fill delta-pt histograms
    result = [self.fill_matching_histograms(jet_combined, jetR, R_max) for jet_combined in jets_combined_selected]

    # Loop through jets and fill groomed histograms if both det and truth jets are unique match
    for grooming_setting in self.grooming_settings:
      if self.debug_level > 1:
        print('grooming setting: {}'.format(grooming_setting))
      result = [self.fill_groomed_jet_matches(grooming_setting, jet_combined, i_jet, jetR, R_max) for i_jet, jet_combined in enumerate(jets_combined_selected)]
      
  #---------------------------------------------------------------
  # Loop through jets and fill matching histos
  #---------------------------------------------------------------
  def fill_matching_histograms(self, jet_combined, jetR, R_max):
          
    if jet_combined.has_user_info():
      jet_truth = jet_combined.python_info().match
      if jet_truth:
      
        # Fill delta pt
        delta_pt = (jet_combined.pt() - jet_truth.pt())
        getattr(self, 'hDeltaPt_emb_R{}_Rmax{}'.format(jetR, R_max)).Fill(jet_truth.pt(), delta_pt)
        
        # Fill subjet histograms
        if 'subjet_z' in self.observable_list:
          for subjetR in self.obs_settings['subjet_z']:
            self.fill_subjet_histograms(jet_combined, jet_truth, jetR, subjetR, R_max)

  #---------------------------------------------------------------
  # Fill subjet histograms
  #---------------------------------------------------------------
  def fill_subjet_histograms(self, jet_combined, jet_truth, jetR, subjetR, R_max):
  
    # Find subjets
    cs_subjet_combined = fj.ClusterSequence(jet_combined.constituents(), self.subjet_def[subjetR])
    subjets_combined = fj.sorted_by_pt(cs_subjet_combined.inclusive_jets())

    cs_subjet_truth = fj.ClusterSequence(jet_truth.constituents(), self.subjet_def[subjetR])
    subjets_truth = fj.sorted_by_pt(cs_subjet_truth.inclusive_jets())
    
    # Find leading subjets and fill matched pt
    leading_subjet_combined = self.leading_subjet(subjets_combined)
    leading_subjet_truth = self.leading_subjet(subjets_truth)
    
    z_leading_combined = leading_subjet_combined.pt() / jet_combined.pt()
    matched_pt = fjtools.matched_pt(leading_subjet_combined, leading_subjet_truth)
    
    name = 'h_subjet_z_fraction_leading_JetPt_R{}_{}_Rmax{}'.format(jetR, subjetR, R_max)
    getattr(self, name).Fill(jet_truth.pt(), z_leading_combined, matched_pt)
    
    # Set subjet matches
    # Beware that defining geometrical subjet matches highly constrains them to match by pt

    # Loop through subjets and set jet matching candidates (based on deltaR) for each jet in user_info
    [[self.set_matching_candidates(subjet_combined, subjet_truth, subjetR,
     'hDeltaR_subjet_z_R{}_{{}}_Rmax{}'.format(jetR, R_max)) for subjet_truth in subjets_truth]
                                           for subjet_combined in subjets_combined]
        
    # Loop through subjets and set accepted matches
    hname = 'hJetMatchingQA_R{}_Rmax{}'.format(jetR, R_max)
    [self.set_matches_AA_truth(subjet_combined, subjetR, hname) for subjet_combined in subjets_combined]
        
    # Loop through matches and fill histograms
    for subjet_combined in subjets_combined:

      if subjet_combined.has_user_info():
        subjet_truth = subjet_combined.python_info().match
      
        if subjet_truth:
          
          z_truth = subjet_truth.pt() / jet_truth.pt()
          
          # Compute fraction of pt of truth subjet contained in matched combined subjet
          matched_pt = fjtools.matched_pt(subjet_combined, subjet_truth)

          name = 'h_subjet_z_fraction_JetPt_R{}_{}_Rmax{}'.format(jetR, subjetR, R_max)
          getattr(self, name).Fill(jet_truth.pt(), z_truth, matched_pt)

  #---------------------------------------------------------------
  # Loop through jets and fill response if both det and truth jets are unique match
  #---------------------------------------------------------------
  def leading_subjet(self, subjets):
  
    leading_subjet = None
    for subjet in subjets:
    
      if not leading_subjet:
        leading_subjet = subjet
        
      if subjet.pt() > leading_subjet.pt():
        leading_subjet = subjet

    return leading_subjet

  #---------------------------------------------------------------
  # Loop through jets and fill response if both det and truth jets are unique match
  #---------------------------------------------------------------
  def fill_groomed_jet_matches(self, grooming_setting, jet_combined, i_jet, jetR, R_max):

    grooming_label = self.utils.grooming_label(grooming_setting)
       
    if jet_combined.has_user_info():
      jet_truth = jet_combined.python_info().match
     
      if jet_truth:

        jet_pt_combined_ungroomed = jet_combined.pt()
        jet_pt_truth_ungroomed = jet_truth.pt()
        
        if self.debug_level > 2:
          print('**** jet_pt_combined_ungroomed: {}, jet_pt_truth_ungroomed: {}'.format(jet_pt_combined_ungroomed, jet_pt_truth_ungroomed))
          
        # Groom combined jet
        gshop_combined = fjcontrib.GroomerShop(jet_combined, jetR, self.reclustering_algorithm)
        jet_combined_groomed_lund = self.utils.groom(gshop_combined, grooming_setting, jetR)
        if not jet_combined_groomed_lund:
          return
        
        # Groom truth jet
        gshop_truth = fjcontrib.GroomerShop(jet_truth, jetR, self.reclustering_algorithm)
        jet_truth_groomed_lund = self.utils.groom(gshop_truth, grooming_setting, jetR)
        if not jet_truth_groomed_lund:
          return
           
        # Fill some variables
        theta_g_combined = jet_combined_groomed_lund.Delta()/jetR
        theta_g_truth = jet_truth_groomed_lund.Delta()/jetR
        zg_combined = jet_combined_groomed_lund.z()
        zg_truth = jet_truth_groomed_lund.z()

        if 'kappa' in self.observable_list:
          kappa_combined = self.kappa(zg_combined, theta_g_combined)
          kappa_truth = self.kappa(zg_truth, theta_g_truth)
        if 'tf' in self.observable_list:
          tf_combined = self.tf(zg_combined, theta_g_combined)
          tf_truth = self.tf(zg_truth, theta_g_truth)
                  
        # Fill prong matching histograms
        if grooming_setting in self.obs_grooming_settings['theta_g']:
          prong_match = self.fill_prong_matching_histograms(jet_truth, jet_truth_groomed_lund, jet_combined, jet_combined_groomed_lund, jet_pt_truth_ungroomed, jetR, grooming_setting, grooming_label, R_max)
          
        # Plot diagram
        if self.plot_diagram:
          self.diagram(jet_truth, jet_combined, prong_match, i_jet, grooming_setting, jetR)

        # Fill combined histograms
        hname = 'h_theta_g_zg_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)
        x = ([jet_pt_truth_ungroomed, zg_combined, theta_g_combined, prong_match])
        x_array = array('d', x)
        getattr(self, hname).Fill(x_array)
        
        if 'kappa' in self.observable_list:
          hname = 'h_kappa_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)
          getattr(self, hname).Fill(jet_pt_truth_ungroomed, kappa_combined, prong_match)
        
        if 'tf' in self.observable_list:
          hname = 'h_tf_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)
          getattr(self, hname).Fill(jet_pt_truth_ungroomed, tf_combined, prong_match)
        
        # Fill truth histograms
        if self.fill_Rmax_indep_hists:
        
          hname = 'h_theta_g_zg_JetPt_Truth_R{}_{}'.format(jetR, grooming_label)
          getattr(self, hname).Fill(jet_pt_truth_ungroomed, zg_truth, theta_g_truth)
          
          if 'kappa' in self.observable_list:
            hname = 'h_kappa_JetPt_Truth_R{}_{}'.format(jetR, grooming_label)
            getattr(self, hname).Fill(jet_pt_truth_ungroomed, kappa_truth)
          
          if 'tf' in self.observable_list:
            hname = 'h_tf_JetPt_Truth_R{}_{}'.format(jetR, grooming_label)
            getattr(self, hname).Fill(jet_pt_truth_ungroomed, tf_truth)
            
        # Fill tree
        if self.fill_tree:

          if R_max == self.main_R_max:
          
            if grooming_setting == self.grooming_settings[0]:
              self.tw.fill_branch('R{}_jet_pt_truth_ungroomed'.format(jetR), jet_pt_truth_ungroomed)
              self.tw.fill_branch('R{}_jet_pt_combined_ungroomed'.format(jetR), jet_pt_combined_ungroomed)
              
            label = 'R{}_Rmax{}_{}'.format(jetR, R_max, grooming_label)
            jet_combined_groomed = jet_combined_groomed_lund.pair()
            jet_pt_combined_groomed = jet_combined_groomed.pt()
            self.tw.fill_branch('{}_jet_pt_combined_groomed'.format(label), jet_pt_combined_groomed)
            self.tw.fill_branch('{}_zg_combined'.format(label), zg_combined)
            self.tw.fill_branch('{}_theta_g_combined'.format(label), theta_g_combined)
            self.tw.fill_branch('{}_prong_matching_flag'.format(label), prong_match)

  #---------------------------------------------------------------
  # Plot diagram
  #---------------------------------------------------------------
  def diagram(self, jet_truth, jet_combined, prong_match, i_jet, grooming_setting, jetR):
  
    if jet_truth.pt() < 40 or jet_truth.pt() > 100:
      return
  
    # Groom truth jet, and get list of all Lund splits
    gshop_truth = fjcontrib.GroomerShop(jet_truth, jetR, self.reclustering_algorithm)
    jet_truth_groomed_lund = self.utils.groom(gshop_truth, grooming_setting, jetR)
    jet_truth_lunds = gshop_truth.lund_splits()
    if not jet_truth_lunds:
      return

    # Loop through Lund splits, and draw diagram
    self.single_diagram(jet_truth, jet_truth_lunds, jet_pt=jet_truth.pt(),
                        i_jet=i_jet, label='truth')
  
    # Groom combined jet, and get list of all Lund splits
    gshop_combined = fjcontrib.GroomerShop(jet_combined, jetR, self.reclustering_algorithm)
    jet_combined_groomed_lund = self.utils.groom(gshop_combined, grooming_setting, jetR)
    jet_combined_lunds = gshop_combined.lund_splits()
    if not jet_combined_lunds:
      return
    
    # Loop through Lund splits, and draw diagram
    self.single_diagram(jet_combined, jet_combined_lunds, jet_pt=jet_truth.pt(),
                        i_jet=i_jet, prong_match=prong_match, label='combined')
    
    #    1: subleading
    #    2: leading, swap (>10% of leading in subleading)
    #    3: leading, mis-tag (<10% of leading in subleading)
    #    4: ungroomed
    #    5: outside
    #    6: other (i.e. 50% is not in any of the above)
    #    7: pp-truth passed grooming, but combined jet failed grooming
    #    8: combined jet passed grooming, but pp-truth failed grooming
    #    9: both pp-truth and combined jet failed SoftDrop
  
  #---------------------------------------------------------------
  # Plot diagram
  #---------------------------------------------------------------
  def single_diagram(self, jet, jet_lunds, jet_pt=0., i_jet=-1, prong_match='', label=''):
    
    # Plot settings
    linewidth=3.
    color_pythia = sns.xkcd_rgb['denim blue']
    color_thermal = sns.xkcd_rgb['pale red']
    color_primary = sns.xkcd_rgb['grey']
    ymin = 0.98
    ymax = 1.08
    
    # Loop through primary Lund plane
    delta = 0.7/len(jet_lunds)
    found_split = False
    for i, split in enumerate(jet_lunds):
      pt = split.pair().pt()
      z = split.z()
      dr = split.Delta()
            
      # Draw softer splitting
      length = pt*z / jet_pt
      x = [delta*(i+1), delta*(i+1) + length*np.cos(dr)]
      y = [1, 1 + length*np.sin(dr)]
      plt.plot(x, y, color_pythia, linewidth=linewidth, label=('PYTHIA' if (i==0 and label != 'truth') else '_'))
      
      # Draw fraction of splitting from background
      prong = split.softer()
      prong_pt = prong.pt()
      matched_pt = 0.
      for p in prong.constituents():
        if p.user_index() >= 0:
          matched_pt += p.pt()
      
      length *= 1 - matched_pt/prong_pt
      x = [delta*(i+1), delta*(i+1) + length*np.cos(dr)]
      y = [1, 1 + length*np.sin(dr)]
      if length > 1e-3:
        plt.plot(x, y, color_thermal, linewidth=linewidth, label=('Background' if i==0 else '_'))
          
      # Identify first splitting passing SD condition
      if not found_split:
        if z > 0.1:
          found_split = True
          x_split = [x[0], x[0]]
          y_split = [y[0], y[0]]
    
    # Draw leading branch
    x = [0, 1]
    y = [1, 1]
    plt.plot(x, y, color_primary, linewidth=linewidth+1.2)
      
    if found_split:
      if prong_match in [1, 2, '']:
        plt.plot(x_split, y_split, color_pythia, marker='o', markersize=12)
      elif prong_match in [3,4,5,6]:
        plt.plot(x_split, y_split, color_thermal, marker='o', markersize=12)
      
    # Draw main legend
    pt_label = r'$p_{{\mathrm{{T,jet}}}} = {:.0f} \;\mathrm{{GeV}}/c$'.format(jet_pt)
    reclustering_label = '{} reclustering'.format(self.recluster_alg)
    grooming_label = r'Soft Drop $z_{\mathrm{cut}}=0.1$'
    if label == 'truth':
      title = '{} \n{} \n{} \n{}'.format(r'$\bf{{pp}}$', grooming_label,
                                         reclustering_label, pt_label)
    else:
      title = '{} \n{} \n{} \n{}'.format(r'$\bf{{pp + thermal}}$', grooming_label,
                                       reclustering_label, pt_label)
    first_legend = plt.legend(title = title, title_fontsize=15,
                                 loc='upper right', fontsize=12)
    ax = plt.gca().add_artist(first_legend)

    axes = plt.gca()
    axes.set_ylim([ymin, ymax])
    plt.savefig(os.path.join(self.output_dir, 'diagram_ev{}_jet{}_{}{}.pdf'.format(self.event_number, i_jet, label, prong_match)))

    plt.close('all')
    
  #---------------------------------------------------------------
  # Do prong-matching
  #---------------------------------------------------------------
  def fill_prong_matching_histograms(self, jet_truth, jet_truth_groomed_lund, jet_combined, jet_combined_groomed_lund, jet_pt_truth_ungroomed, jetR, grooming_setting, grooming_label, R_max):
    
    # Dynamical grooming returns a fjcontrib::LundGenerator
    #   The prongs can be retrieved directly from this object.
    #   If the object exists, then it has passed grooming
    jet_truth_prong1 = jet_truth_groomed_lund.harder()
    jet_truth_prong2 = jet_truth_groomed_lund.softer()

    # Get prongs of combined jet
    jet_combined_prong1 = jet_combined_groomed_lund.harder()
    jet_combined_prong2 = jet_combined_groomed_lund.softer()
    
    # Get the fastjet::PseudoJets from the fjcontrib::LundGenerators
    jet_truth_groomed = jet_truth_groomed_lund.pair()
    jet_combined_groomed = jet_combined_groomed_lund.pair()
    
    has_parents_truth = jet_truth_groomed.has_constituents()
    has_parents_combined = jet_combined_groomed.has_constituents()
    
    # Check that groomed truth jet doesn't contain any background tracks
    problem = False
    if has_parents_truth:
      for constituent in jet_truth_groomed.constituents():
        if constituent.user_index() < 0:
          problem = True
    if problem:
      print(grooming_setting)
      print(dir(jet_truth_groomed_lund))
      print('pair: {}'.format(jet_truth_groomed_lund.pair()))
      print('kappa: {}'.format(jet_truth_groomed_lund.kappa()))
      print('groomed constituents: {}'.format([track.user_index() for track in jet_truth_groomed.constituents()]))
      print('jet constituents: {}'.format([track.user_index() for track in jet_truth.constituents()]))
      print('prong1 constituents: {}'.format([track.user_index() for track in jet_truth_prong1.constituents()]))
      print('prong2 constituents: {}'.format([track.user_index() for track in jet_truth_prong2.constituents()]))
          
    if self.debug_level > 1:

        if jet_pt_truth_ungroomed > 80.:
        
            print('=======================================================')
            print(type)
            print('jet_pt_truth_ungroomed: {}'.format(jet_pt_truth_ungroomed))
            print('jet_pt_truth_groomed: {}'.format(jet_truth_groomed.pt()))
            print('jet_pt_combined_groomed: {}'.format(jet_combined_groomed.pt()))
            print('')
            print('jet_truth tracks: {}'.format([track.user_index() for track in jet_truth.constituents()]))
            print('         track pt: {}'.format([np.around(track.pt(),2) for track in jet_truth.constituents()]))
            if jet_truth_groomed.has_constituents():
              print('jet_truth_groomed tracks: {}'.format([track.user_index() for track in jet_truth_groomed.constituents()]))
              print('                 track pt: {}'.format([np.around(track.pt(),2) for track in jet_truth_groomed.constituents()]))
            if jet_combined_groomed.has_constituents():
              print('jet_combined groomed tracks: {}'.format([track.user_index() for track in jet_combined_groomed.constituents()]))
              print('                   track pt: {}'.format([np.around(track.pt(),2) for track in jet_combined_groomed.constituents()]))
            if jet_combined.has_constituents():
              print('jet_combined ungroomed tracks: {}'.format([track.user_index() for track in jet_combined.constituents()]))
              print('                     track pt: {}'.format([np.around(track.pt(),2) for track in jet_combined.constituents()]))

    # Compute fraction of pt of the pp-truth prong tracks that is contained in the combined-jet prong,
    # in order to have a measure of whether the combined-jet prong is the "same" prong as the pp-truth prong
    deltaR_prong1 = -1.
    deltaR_prong2 = -1.
    deltaZ = -1.
    rg_truth = -1.
    zg_truth = -1.
    if has_parents_truth and has_parents_combined:
    
        # Subleading jet pt-matching
        # --------------------------
        # (1) Fraction of pt matched: subleading pp-det in subleading combined
        matched_pt_subleading_subleading = fjtools.matched_pt(jet_combined_prong2, jet_truth_prong2)
        
        # (2) Fraction of pt matched: subleading pp-det in leading combined
        matched_pt_subleading_leading = fjtools.matched_pt(jet_combined_prong1, jet_truth_prong2)
        
        # (3) Fraction of pt matched: subleading pp-det in ungroomed combined jet
        matched_pt_subleading_groomed = fjtools.matched_pt(jet_combined_groomed, jet_truth_prong2)
        matched_pt_subleading_ungroomed = fjtools.matched_pt(jet_combined, jet_truth_prong2)
        matched_pt_subleading_ungroomed_notgroomed = matched_pt_subleading_ungroomed - matched_pt_subleading_groomed
        
        # (4) Fraction of pt matched: subleading pp-det not in ungroomed combined jet
        matched_pt_subleading_outside = 1 - matched_pt_subleading_ungroomed

        # Leading jet pt-matching
        # --------------------------
        # (1) Fraction of pt matched: leading pp-det in subleading combined
        matched_pt_leading_subleading = fjtools.matched_pt(jet_combined_prong2, jet_truth_prong1)
        
        # (2) Fraction of pt matched: leading pp-det in leading combined
        matched_pt_leading_leading = fjtools.matched_pt(jet_combined_prong1, jet_truth_prong1)

        # (3) Fraction of pt matched: leading pp-det in ungroomed combined jet
        matched_pt_leading_groomed = fjtools.matched_pt(jet_combined_groomed, jet_truth_prong1)
        matched_pt_leading_ungroomed = fjtools.matched_pt(jet_combined, jet_truth_prong1)
        matched_pt_leading_ungroomed_notgroomed = matched_pt_leading_ungroomed - matched_pt_leading_groomed
        
        # (4) Fraction of pt matched: leading pp-det not in ungroomed combined jet
        matched_pt_leading_outside = 1 - matched_pt_leading_ungroomed

        # Compute delta-R between pp-det prong and combined prong
        # --------------------------
        deltaR_prong1 = jet_combined_prong1.delta_R(jet_truth_prong1)
        deltaR_prong2 = jet_combined_prong2.delta_R(jet_truth_prong2)
        deltaZ = jet_combined_groomed_lund.z() - jet_truth_groomed_lund.z()
        rg_truth = jet_truth_groomed_lund.Delta()
        zg_truth = jet_truth_groomed_lund.z()
        if self.debug_level > 3:
          if rg_truth < 0.:
            print(rg_truth)
        
        if self.debug_level > 1:
        
            if jet_pt_truth_ungroomed > 80.:
            
                print('subleading prong tracks -- combined: {}'.format([track.user_index() for track in jet_combined_prong2.constituents()]))
                print('subleading prong tracks -- pp-truth: {}'.format([track.user_index() for track in jet_truth_prong2.constituents()]))
                print('leading prong tracks -- combined: {}'.format([track.user_index() for track in jet_combined_prong1.constituents()]))
                print('leading prong tracks -- pp-truth: {}'.format([track.user_index() for track in jet_truth_prong1.constituents()]))
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

    elif has_parents_truth: # pp-truth passed grooming, but combined jet failed grooming
        matched_pt_leading_leading = matched_pt_leading_subleading = matched_pt_leading_ungroomed_notgroomed = matched_pt_leading_outside = matched_pt_subleading_leading = matched_pt_subleading_subleading = matched_pt_subleading_ungroomed_notgroomed = matched_pt_subleading_outside = -0.1
        
    elif has_parents_combined: # combined jet passed grooming, but pp-det failed grooming
        matched_pt_leading_leading = matched_pt_leading_subleading = matched_pt_leading_ungroomed_notgroomed = matched_pt_leading_outside = matched_pt_subleading_leading = matched_pt_subleading_subleading = matched_pt_subleading_ungroomed_notgroomed = matched_pt_subleading_outside = -0.2
        
    else: # both pp-det and combined jet failed SoftDrop
        matched_pt_leading_leading = matched_pt_leading_subleading = matched_pt_leading_ungroomed_notgroomed = matched_pt_leading_outside = matched_pt_subleading_leading = matched_pt_subleading_subleading = matched_pt_subleading_ungroomed_notgroomed = matched_pt_subleading_outside = -0.3

    # Leading prong
    getattr(self, 'hProngMatching_leading_leading_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_leading, deltaR_prong1/rg_truth)
    getattr(self, 'hProngMatching_leading_subleading_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_subleading, deltaR_prong1/rg_truth)
    getattr(self, 'hProngMatching_leading_ungroomed_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_ungroomed_notgroomed, deltaR_prong1/rg_truth)
    getattr(self, 'hProngMatching_leading_outside_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_outside, deltaR_prong1/rg_truth)
    
    getattr(self, 'hProngMatching_leading_leading_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_leading, deltaZ/zg_truth)
    getattr(self, 'hProngMatching_leading_subleading_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_subleading, deltaZ/zg_truth)
    getattr(self, 'hProngMatching_leading_ungroomed_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_ungroomed_notgroomed, deltaZ/zg_truth)
    getattr(self, 'hProngMatching_leading_outside_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_outside, deltaZ/zg_truth)

    # Subleading prong
    getattr(self, 'hProngMatching_subleading_leading_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_leading, deltaR_prong2/rg_truth)
    getattr(self, 'hProngMatching_subleading_subleading_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_subleading, deltaR_prong2/rg_truth)
    getattr(self, 'hProngMatching_subleading_ungroomed_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_ungroomed_notgroomed, deltaR_prong2/rg_truth)
    getattr(self, 'hProngMatching_subleading_outside_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_outside, deltaR_prong2/rg_truth)

    getattr(self, 'hProngMatching_subleading_leading_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_leading, deltaZ/zg_truth)
    getattr(self, 'hProngMatching_subleading_subleading_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_subleading, deltaZ/zg_truth)
    getattr(self, 'hProngMatching_subleading_ungroomed_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_ungroomed_notgroomed, deltaZ/zg_truth)
    getattr(self, 'hProngMatching_subleading_outside_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_outside, deltaZ/zg_truth)
    
    # Plot correlation of matched pt fraction for leading-subleading and subleading-leading
    getattr(self, 'hProngMatching_subleading-leading_correlation_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_truth.pt(), matched_pt_leading_subleading, matched_pt_subleading_leading)
    
    # Plot fraction of groomed pp truth pt that is contained in groomed combined jet
    if has_parents_truth and has_parents_combined:
      matched_pt_truth_groomed = fjtools.matched_pt(jet_combined_groomed, jet_truth_groomed)
      getattr(self, 'hProngMatching_truth_groomed_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_truth_groomed)
      
      matched_pt_truth2_groomed = fjtools.matched_pt(jet_combined_groomed, jet_truth_prong2)
      getattr(self, 'hProngMatching_truth2_groomed_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_truth2_groomed)
    
    # Plot initial fraction of subleading truth prong pt in combined ungroomed prong
    # To do...
     
    #  Return flag based on where >50% of subleading matched pt resides:
    #    1: subleading
    #    2: leading, swap (>10% of leading in subleading)
    #    3: leading, mis-tag (<10% of leading in subleading)
    #    4: ungroomed
    #    5: outside
    #    6: other (i.e. 50% is not in any of the above)
    #    7: pp-truth passed grooming, but combined jet failed grooming
    #    8: combined jet passed grooming, but pp-truth failed grooming
    #    9: both pp-truth and combined jet failed SoftDrop
    if matched_pt_subleading_subleading > self.prong_matching_threshold:
      return 1
    elif matched_pt_subleading_leading > self.prong_matching_threshold:
      if matched_pt_leading_subleading > 0.1:
        return 2
      else:
        return 3
    elif matched_pt_subleading_ungroomed_notgroomed > self.prong_matching_threshold:
      return 4
    elif matched_pt_subleading_outside > self.prong_matching_threshold:
      return 5
    elif matched_pt_leading_leading >= 0.:
      return 6
    elif matched_pt_leading_leading == -0.1:
      return 7
    elif matched_pt_leading_leading == -0.2:
      return 8
    elif matched_pt_leading_leading == -0.3:
      return 9
    else:
      print('Warning -- flag not specified!')
      return -1
    
  #---------------------------------------------------------------
  # Compute kappa
  #---------------------------------------------------------------
  def kappa(self, z, theta):
    return z*theta
  
  #---------------------------------------------------------------
  # Compute tf
  #---------------------------------------------------------------
  def tf(self, z, theta):
    return z*theta*theta

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

  analysis = ProcessGroomers(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_groomers()
