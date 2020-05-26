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

    # Initialize histograms
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
      
    self.jet_matching_distance = config['jet_matching_distance']
    self.mc_fraction_threshold = config['mc_fraction_threshold']
    self.prong_matching_threshold = config['prong_matching_threshold']
    self.use_ev_id_ext = config['use_ev_id_ext']
    self.main_R_max = config['constituent_subtractor']['main_R_max']
    
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
    if self.event_number_max < self.nEvents:
      self.hNevents.Fill(1, self.event_number_max)
    else:
      self.hNevents.Fill(1, self.nEvents)

    self.hRho =  ROOT.TH1F('hRho', 'hRho', 1000, 0., 1000.)

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
        
      #  Construct THn for each groomer: (pt, zg, theta_g, tag flag)
      for grooming_setting in self.obs_grooming_settings['theta_g']:
        if grooming_setting:
          grooming_label = self.utils.grooming_label(grooming_setting)
          
          for R_max in self.max_distance:
            
            # THn for combined jets
            dim = 4;
            title = ['p_{T,ch jet}', '#it{z}_{g,ch}', '#theta_{g,ch}', 'flag']
            nbins = [30, 50, 100, 8]
            min = [0., 0., 0., 0.5]
            max = [300., 0.5, 1., 8.5]
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
                h = ROOT.TH3F(name, name, 30, 0, 300, 50, 0, 0.5, 8, 0.5, 8.5)
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
          h = ROOT.TH3F(name, name, 30, 0, 300, 15, -0.4, 1.1, 20, 0., 2*jetR)
          h.GetXaxis().SetTitle('p_{T,truth}')
          h.GetYaxis().SetTitle('Prong matching fraction')
          h.GetZaxis().SetTitle('#Delta R_{prong}')
          setattr(self, name, h)
          
          name = 'hProngMatching_{}_{}_JetPtZ_R{}_{}_Rmax{}'.format(prong, match, jetR, grooming_label, R_max)
          h = ROOT.TH3F(name, name, 30, 0, 300, 15, -0.4, 1.1, 50, -0.5, 0.5)
          h.GetXaxis().SetTitle('p_{T,truth}')
          h.GetYaxis().SetTitle('Prong matching fraction')
          h.GetZaxis().SetTitle('#Delta z_{prong}')
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
      
    # Clear event in tree
    if self.fill_tree:
      self.tw.clear()
  
    # If Pb-Pb, construct embedded event (do this once, for all jetR)
        
    # If thermal model, generate a thermal event and add it to the det-level particle list
    if self.thermal_model:
      fj_particles_combined_beforeCS = self.thermal_generator.load_event()
      
      # Form the combined event
      # The pp-truth tracks are each stored with a unique user_index >= 0
      #   (same index in fj_particles_combined and fj_particles_truth -- which will be used in prong-matching)
      # The thermal tracks are each stored with a unique user_index < 0
      [fj_particles_combined_beforeCS.push_back(p) for p in fj_particles_truth]

    # Main case: Get Pb-Pb event and embed it into the det-level particle list
    else:
      fj_particles_combined_beforeCS = self.process_io_emb.load_event()
          
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
      print([p.user_index() for p in fj_particles_combined])
      print([p.pt() for p in fj_particles_combined])
      print([p.user_index() for p in fj_particles_combined_beforeCS])
      print([p.pt() for p in fj_particles_combined_beforeCS])

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
      result = [self.fill_groomed_jet_matches(grooming_setting, jet_combined, jetR, R_max) for jet_combined in jets_combined_selected]
      
  #---------------------------------------------------------------
  # Loop through jets and fill matching histos
  #---------------------------------------------------------------
  def fill_matching_histograms(self, jet_combined, jetR, R_max):
          
    if jet_combined.has_user_info():
      jet_truth = jet_combined.python_info().match
      if jet_truth:
        delta_pt = (jet_combined.pt() - jet_truth.pt())
        getattr(self, 'hDeltaPt_emb_R{}_Rmax{}'.format(jetR, R_max)).Fill(jet_truth.pt(), delta_pt)

  #---------------------------------------------------------------
  # Loop through jets and fill response if both det and truth jets are unique match
  #---------------------------------------------------------------
  def fill_groomed_jet_matches(self, grooming_setting, jet_combined, jetR, R_max):
       
    grooming_label = self.utils.grooming_label(grooming_setting)
       
    if jet_combined.has_user_info():
      jet_truth = jet_combined.python_info().match
     
      if jet_truth:
     
        jet_pt_combined_ungroomed = jet_combined.pt()
        jet_pt_truth_ungroomed = jet_truth.pt()
        
        if self.debug_level > 2:
          print('**** jet_pt_combined_ungroomed: {}, jet_pt_truth_ungroomed: {}'.format(jet_pt_combined_ungroomed, jet_pt_truth_ungroomed))
       
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
 
          jet_combined_sd = sd.result(jet_combined)
          jet_truth_sd = sd.result(jet_truth)

        # Construct Dynamical groomer, and groom jet
        if 'dg' in grooming_setting:

          a = grooming_setting['dg'][0]
          
          jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
          dy_groomer = fjcontrib.DynamicalGroomer(jet_def_lund)
          if self.debug_level > 2:
            print('Dynamical groomer is: {}'.format(dy_groomer.description()))
          
          jet_combined_dg_lund = self.utils.dy_groom(dy_groomer, jet_combined, a)
          jet_truth_dg_lund = self.utils.dy_groom(dy_groomer, jet_truth, a)
          
          jet_combined_dg = jet_combined_dg_lund.pair()
          jet_truth_dg = jet_truth_dg_lund.pair()
            
        # Compute groomed observables
        if 'sd' in grooming_setting:

          # If both SD and DG are specified, first apply DG, then SD
          if 'dg' in grooming_setting:
            if self.debug_level > 2:
              print('both SD and DG applied')
            if jet_combined_dg.has_constituents() and jet_truth_dg.has_constituents():
              jet_combined_groomed = sd.result(jet_combined_dg)
              jet_truth_groomed = sd.result(jet_truth_dg)
            else:
              return
          else:
            if self.debug_level > 2:
              print('SD applied')
            jet_combined_groomed = jet_combined_sd
            jet_truth_groomed = jet_truth_sd
            
          theta_g_combined = self.theta_g(jet_combined_groomed, jetR)
          theta_g_truth = self.theta_g(jet_truth_groomed, jetR)
          zg_combined = self.zg(jet_combined_groomed)
          zg_truth = self.zg(jet_truth_groomed)

        elif 'dg' in grooming_setting:
        
          if self.debug_level > 2:
            print('DG applied')
          
          jet_combined_groomed = jet_combined_dg
          jet_truth_groomed = jet_truth_dg
          
          theta_g_combined = jet_combined_dg_lund.Delta()/jetR
          theta_g_truth = jet_truth_dg_lund.Delta()/jetR
          zg_combined = jet_combined_dg_lund.z()
          zg_truth = jet_truth_dg_lund.z()
       
        # Fill some variables
        if 'kappa' in self.observable_list:
          kappa_combined = self.kappa(zg_combined, theta_g_combined)
          kappa_truth = self.kappa(zg_truth, theta_g_truth)
        if 'tf' in self.observable_list:
          tf_combined = self.tf(zg_combined, theta_g_combined)
          tf_truth = self.tf(zg_truth, theta_g_truth)
       
        #  Fill prong-matching histograms
        if 'sd' in grooming_setting:

          if 'dg' in grooming_setting:
            groomer_list = [sd, dy_groomer]
            if grooming_setting in self.obs_grooming_settings['theta_g']:
              prong_match = self.fill_prong_matching_histograms(jet_truth, jet_combined, jet_combined_groomed, groomer_list, jet_pt_truth_ungroomed, jetR, grooming_setting, grooming_label, R_max, type = 'SD+DG')
          else:
            groomer_list = [sd]
            if grooming_setting in self.obs_grooming_settings['theta_g']:
              prong_match = self.fill_prong_matching_histograms(jet_truth, jet_combined, jet_combined_groomed, groomer_list, jet_pt_truth_ungroomed, jetR, grooming_setting, grooming_label, R_max, type = 'SD')
    
        elif 'dg' in grooming_setting:

          if grooming_setting in self.obs_grooming_settings['theta_g']:
            prong_match = self.fill_prong_matching_histograms(jet_truth, jet_combined, jet_combined_dg_lund, [dy_groomer], jet_pt_truth_ungroomed, jetR, grooming_setting, grooming_label, R_max, type = 'DG')

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
            jet_pt_combined_groomed = jet_combined_groomed.pt()
            self.tw.fill_branch('{}_jet_pt_combined_groomed'.format(label), jet_pt_combined_groomed)
            self.tw.fill_branch('{}_zg_combined'.format(label), zg_combined)
            self.tw.fill_branch('{}_theta_g_combined'.format(label), theta_g_combined)
            self.tw.fill_branch('{}_prong_matching_flag'.format(label), prong_match)
          
  #---------------------------------------------------------------
  # Do prong-matching
  #---------------------------------------------------------------
  def fill_prong_matching_histograms(self, jet_truth, jet_combined, jet_combined_groomed, groomer_list, jet_pt_truth_ungroomed, jetR, grooming_setting, grooming_label, R_max, type = 'SD'):
    
    # Do grooming on pp-truth jet, and get prongs
    if 'SD' in type:
    
      if 'DG' in type:
      
        # Assumes groomer_list = [sd, dy_groomer]
        a = grooming_setting['dg'][0]
        jet_truth_dg_lund = self.utils.dy_groom(groomer_list[1], jet_truth, a)
        jet_truth_dg = jet_truth_dg_lund.pair()
        jet_truth_groomed = groomer_list[0].result(jet_truth_dg)
      
      else:
      
        jet_truth_groomed = groomer_list[0].result(jet_truth)
    
      # SD grooming returns the groomed fastjet::PseudoJet
      # Use the fastjet::PseudoJet::has_parents function which returns the last clustering step
      #   If the jet passed SoftDrop, then its parents are the SoftDrop splitting
      #   If the jet didn't pass SoftDrop, then it will have no parents
      jet_truth_prong1 = fj.PseudoJet()
      jet_truth_prong2 = fj.PseudoJet()
      has_parents_truth = jet_truth_groomed.has_parents(jet_truth_prong1, jet_truth_prong2)
      
      # Get prongs of combined jet
      jet_combined_prong1 = fj.PseudoJet()
      jet_combined_prong2 = fj.PseudoJet()
      has_parents_combined = jet_combined_groomed.has_parents(jet_combined_prong1, jet_combined_prong2)

    elif type == 'DG':
      
      a = grooming_setting['dg'][0]
      jet_truth_groomed_lund = self.utils.dy_groom(groomer_list[0], jet_truth, a)
    
      # Dynamical grooming returns a fjcontrib::LundGenerator
      #   The prongs can be retrieved directly from this object.
      #   If the object exists, then it has passed grooming
      jet_truth_prong1 = jet_truth_groomed_lund.harder()
      jet_truth_prong2 = jet_truth_groomed_lund.softer()
      has_parents_truth = jet_truth_groomed_lund

      # Get prongs of combined jet
      jet_combined_prong1 = jet_combined_groomed.harder()
      jet_combined_prong2 = jet_combined_groomed.softer()
      has_parents_combined = jet_combined_groomed
      
      # Get the fastjet::PseudoJets from the fjcontrib::LundGenerators
      jet_truth_groomed = jet_truth_groomed_lund.pair()
      jet_combined_groomed = jet_combined_groomed.pair()
          
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
            print('jet_truth_groomed tracks: {}'.format([track.user_index() for track in jet_truth_groomed.constituents()]))
            print('                 track pt: {}'.format([np.around(track.pt(),2) for track in jet_truth_groomed.constituents()]))
            print('jet_combined groomed tracks: {}'.format([track.user_index() for track in jet_combined_groomed.constituents()]))
            print('                   track pt: {}'.format([np.around(track.pt(),2) for track in jet_combined_groomed.constituents()]))
            print('jet_combined ungroomed tracks: {}'.format([track.user_index() for track in jet_combined.constituents()]))
            print('                     track pt: {}'.format([np.around(track.pt(),2) for track in jet_combined.constituents()]))

    # Compute fraction of pt of the pp-truth prong tracks that is contained in the combined-jet prong,
    # in order to have a measure of whether the combined-jet prong is the "same" prong as the pp-truth prong
    deltaR_prong1 = -1.
    deltaR_prong2 = -1.
    deltaZ = -1.
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
        deltaZ = self.zg(jet_combined_groomed) - self.zg(jet_truth_groomed)
        
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
    getattr(self, 'hProngMatching_leading_leading_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_leading, deltaR_prong1)
    getattr(self, 'hProngMatching_leading_subleading_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_subleading, deltaR_prong1)
    getattr(self, 'hProngMatching_leading_ungroomed_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_ungroomed_notgroomed, deltaR_prong1)
    getattr(self, 'hProngMatching_leading_outside_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_outside, deltaR_prong1)
    
    getattr(self, 'hProngMatching_leading_leading_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_leading, deltaZ)
    getattr(self, 'hProngMatching_leading_subleading_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_subleading, deltaZ)
    getattr(self, 'hProngMatching_leading_ungroomed_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_ungroomed_notgroomed, deltaZ)
    getattr(self, 'hProngMatching_leading_outside_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_leading_outside, deltaZ)

    # Subleading prong
    getattr(self, 'hProngMatching_subleading_leading_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_leading, deltaR_prong2)
    getattr(self, 'hProngMatching_subleading_subleading_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_subleading, deltaR_prong2)
    getattr(self, 'hProngMatching_subleading_ungroomed_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_ungroomed_notgroomed, deltaR_prong2)
    getattr(self, 'hProngMatching_subleading_outside_JetPt_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_outside, deltaR_prong2)

    getattr(self, 'hProngMatching_subleading_leading_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_leading, deltaZ)
    getattr(self, 'hProngMatching_subleading_subleading_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_subleading, deltaZ)
    getattr(self, 'hProngMatching_subleading_ungroomed_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_ungroomed_notgroomed, deltaZ)
    getattr(self, 'hProngMatching_subleading_outside_JetPtZ_R{}_{}_Rmax{}'.format(jetR, grooming_label, R_max)).Fill(jet_pt_truth_ungroomed, matched_pt_subleading_outside, deltaZ)
    
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
    #    2: leading
    #    3: ungroomed
    #    4: outside
    #    5: other (i.e. 50% is not in any of the above)
    #    6: pp-truth passed grooming, but combined jet failed grooming
    #    7: combined jet passed grooming, but pp-truth failed grooming
    #    8: both pp-truth and combined jet failed SoftDrop
    if matched_pt_subleading_subleading > self.prong_matching_threshold:
      return 1
    elif matched_pt_subleading_leading > self.prong_matching_threshold:
      return 2
    elif matched_pt_subleading_ungroomed_notgroomed > self.prong_matching_threshold:
      return 3
    elif matched_pt_subleading_outside > self.prong_matching_threshold:
      return 4
    elif matched_pt_leading_leading >= 0.:
      return 5
    elif matched_pt_leading_leading == -0.1:
      return 6
    elif matched_pt_leading_leading == -0.2:
      return 7
    elif matched_pt_leading_leading == -0.3:
      return 8
    else:
      print('Warning -- flag not specified!')
      return -1
        
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
