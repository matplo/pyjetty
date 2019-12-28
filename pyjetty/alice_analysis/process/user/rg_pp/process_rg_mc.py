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
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_io_emb
from pyjetty.alice_analysis.process.base import process_utils
from pyjetty.alice_analysis.process.base import jet_info
from pyjetty.alice_analysis.process.base import process_base
from pyjetty.mputils import treewriter
from pyjetty.mputils import CEventSubtractor

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class process_rg_mc(process_base.process_base):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
    super(process_rg_mc, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
    
    # Initialize configuration
    self.initialize_config()
  
  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def process_rg_mc(self):
    
    start_time = time.time()
    
    # ------------------------------------------------------------------------
    
    # Use IO helper class to convert detector-level ROOT TTree into
    # a SeriesGroupBy object of fastjet particles per event
    print('--- {} seconds ---'.format(time.time() - start_time))
    io_det = process_io.process_io(input_file=self.input_file, track_tree_name='tree_Particle')
    df_fjparticles_det = io_det.load_data(self.reject_tracks_fraction)
    self.nEvents_det = len(df_fjparticles_det.index)
    self.nTracks_det = len(io_det.track_df.index)
    print('--- {} seconds ---'.format(time.time() - start_time))
    
    # ------------------------------------------------------------------------

    # Use IO helper class to convert truth-level ROOT TTree into
    # a SeriesGroupBy object of fastjet particles per event
    io_truth = process_io.process_io(input_file=self.input_file, track_tree_name='tree_Particle_gen')
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
    
    # Set up the Pb-Pb embedding object
    if not self.is_pp:
        self.process_io_emb = process_io_emb.process_io_emb(self.emb_file_list, track_tree_name='tree_Particle')
    
    # ------------------------------------------------------------------------

    # Initialize histograms
    self.initialize_output_objects()
    
    # Create constituent subtractor, if configured
    if self.do_constituent_subtraction:
      self.constituent_subtractor = CEventSubtractor(max_distance=self.max_distance, alpha=self.alpha, max_eta=self.max_eta, bge_rho_grid_size=self.bge_rho_grid_size, max_pt_correct=self.max_pt_correct, ghost_area=self.ghost_area, distance_type=fjcontrib.ConstituentSubtractor.deltaR)
    
    print(self)
  
    # Find jets and fill histograms
    print('Find jets...')
    self.analyzeEvents()
    
    # Plot histograms
    print('Save histograms...')
    process_base.process_base.save_output_objects(self)
    
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
    
    # Write tree output (default is to write only histograms)
    self.write_tree_output = config['write_tree_output']

    self.jet_matching_distance = config['jet_matching_distance']
    self.mc_fraction_threshold = config['mc_fraction_threshold']
    self.reject_tracks_fraction = config['reject_tracks_fraction']
    
    # Retrieve list of SD grooming settings
    sd_config_dict = config['SoftDrop']
    sd_config_list = list(sd_config_dict.keys())
    self.sd_settings = [[sd_config_dict[name]['zcut'], sd_config_dict[name]['beta']] for name in sd_config_list]

    if self.do_constituent_subtraction:
        self.is_pp = False
        self.emb_file_list = config['emb_file_list']
    else:
        self.is_pp = True

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_output_objects(self):
    
    self.hNevents = ROOT.TH1F('hNevents', 'hNevents', 2, -0.5, 1.5)
    self.hNevents.Fill(1, self.nEvents_det)
    
    self.hTrackEtaPhi = ROOT.TH2F('hTrackEtaPhi', 'hTrackEtaPhi', 200, -1., 1., 628, 0., 6.28)
    self.hTrackPt = ROOT.TH1F('hTrackPt', 'hTrackPt', 300, 0., 300.)

    for jetR in self.jetR_list:
      
      name = 'hJES_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
      setattr(self, name, h)
      
      if self.is_pp:
          name = 'hDeltaR_All_R{}'.format(jetR)
          h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
          setattr(self, name, h)
          
          name = 'hJetMatchingQA_R{}'.format(jetR)
          bin_labels = ['all', 'has_matching_candidate', 'unique_match']
          nbins = len(bin_labels)
          h = ROOT.TH2F(name, name, nbins, 0, nbins, 30, 0., 300.)
          for i in range(1, nbins+1):
            h.GetXaxis().SetBinLabel(i,bin_labels[i-1])
          setattr(self, name, h)
      else:
          name = 'hDeltaR_ppdet_pptrue_R{}'.format(jetR)
          h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
          setattr(self, name, h)
          
          name = 'hDeltaR_combined_ppdet_R{}'.format(jetR)
          h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
          setattr(self, name, h)
          
          name = 'hJetMatchingQA_R{}'.format(jetR)
          bin_labels = ['all', 'has_matching_ppdet_candidate', 'has_matching_pptrue_candidate', 'has_matching_pptrue_unique_candidate', 'mc_fraction', 'deltaR_combined-truth', 'passed_all_cuts']
          nbins = len(bin_labels)
          h = ROOT.TH2F(name, name, nbins, 0, nbins,  30, 0., 300.)
          for i in range(1, nbins+1):
            h.GetXaxis().SetBinLabel(i,bin_labels[i-1])
          setattr(self, name, h)
      
      name = 'hZ_Truth_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 1.)
      setattr(self, name, h)
      
      name = 'hZ_Det_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 1.)
      setattr(self, name, h)
      
      if not self.write_tree_output:
      
        name = 'hJetPt_Truth_R{}'.format(jetR)
        h = ROOT.TH1F(name, name, 300, 0, 300)
        setattr(self, name, h)

      for sd_setting in self.sd_settings:
        
        zcut = sd_setting[0]
        beta = sd_setting[1]
        sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)
        
        name = 'hThetaG_RelativeResidual_JetPt_R{}_{}'.format(jetR, sd_label)
        h = ROOT.TH2F(name, name, 300, 0, 300, 200, -2., 2.)
        h.GetXaxis().SetTitle('p_{T,truth}')
        h.GetYaxis().SetTitle('#frac{#theta_{g,det}-#theta_{g,truth}}{#theta_{g,truth}}')
        setattr(self, name, h)
        
        name = 'hZg_RelativeResidual_JetPt_R{}_{}'.format(jetR, sd_label)
        h = ROOT.TH2F(name, name, 300, 0, 300, 200, -2., 2.)
        h.GetXaxis().SetTitle('p_{T,truth}')
        h.GetYaxis().SetTitle('#frac{z_{g,det}-z_{g,truth}}{z_{g,truth}}')
        setattr(self, name, h)
        
        name = 'hThetaGResidual_JetPt_R{}_{}'.format(jetR, sd_label)
        h = ROOT.TH2F(name, name, 300, 0, 300, int(3*jetR*100), -3*jetR, 3*jetR)
        h.GetXaxis().SetTitle('p_{T,truth}')
        h.GetYaxis().SetTitle('#theta_{g,det}-#theta_{g,truth}')
        setattr(self, name, h)
        
        name = 'hZgResidual_JetPt_R{}_{}'.format(jetR, sd_label)
        h = ROOT.TH2F(name, name, 300, 0, 300, 100, -0.5, 0.5)
        h.GetXaxis().SetTitle('p_{T,truth}')
        h.GetYaxis().SetTitle('z_{g,det}-z_{g,truth}')
        setattr(self, name, h)

        if  self.write_tree_output:

          # Initialize tree to write out event variables
          tree_name = 't_R{}_{}'.format(self.utils.remove_periods(jetR), sd_label)
          t = ROOT.TTree(tree_name, tree_name)
          tree_writer_name = 'tree_writer_R{}_{}'.format(self.utils.remove_periods(jetR), sd_label)
          tree_writer = treewriter.RTreeWriter(tree=t, tree_name=tree_name, name=tree_writer_name, file_name=None)
          setattr(self, tree_writer_name, tree_writer)
            
        else:
          
          # Create THn of response for theta_g
          dim = 4;
          title = ['p_{T,det}', 'p_{T,truth}', '#theta_{g,det}', '#theta_{g,truth}']
          nbins = [100, 60, 130, 26]
          min = [0., 0., 0., 0.]
          max = [100., 300., 1.3, 1.3]

          name = 'hResponse_JetPt_ThetaG_R{}_{}'.format(jetR, sd_label)
          nbins = (nbins)
          xmin = (min)
          xmax = (max)
          nbins_array = array('i', nbins)
          xmin_array = array('d', xmin)
          xmax_array = array('d', xmax)
          h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
          for i in range(0, dim):
            h.GetAxis(i).SetTitle(title[i])
            setattr(self, name, h)
              
          # Create THn of response for z_g
          dim = 4;
          title = ['p_{T,det}', 'p_{T,truth}', 'z_{g,det}', 'z_{g,truth}']
          nbins = [100, 60, 50, 10]
          min = [0., 0., 0., 0.]
          max = [100., 300., 0.5, 0.5]

          name = 'hResponse_JetPt_zg_R{}_{}'.format(jetR, sd_label)
          nbins = (nbins)
          xmin = (min)
          xmax = (max)
          nbins_array = array('i', nbins)
          xmin_array = array('d', xmin)
          xmax_array = array('d', xmax)
          h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
          for i in range(0, dim):
            h.GetAxis(i).SetTitle(title[i])
            setattr(self, name, h)

  #---------------------------------------------------------------
  # Main function to loop through and analyze events
  #---------------------------------------------------------------
  def analyzeEvents(self):
    
    # Fill track histograms
    [self.fillTrackHistograms(fj_particles_det) for fj_particles_det in self.df_fjparticles['fj_particles_det']]
    
    fj.ClusterSequence.print_banner()
    print()
    
    for jetR in self.jetR_list:
      
      # Set jet definition and a jet selector
      jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
      jet_selector_det = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
      jet_selector_truth_matched = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9)
      print('jet definition is:', jet_def)
      print('jet selector for det-level is:', jet_selector_det,'\n')
      print('jet selector for truth-level matches is:', jet_selector_truth_matched,'\n')
      
      # Then can use list comprehension to iterate over the groupby and do jet-finding
      # simultaneously for fj_1 and fj_2 per event, so that I can match jets -- and fill histograms
      result = [self.analyzeJets(fj_particles_det, fj_particles_truth, jet_def, jet_selector_det, jet_selector_truth_matched) for fj_particles_det, fj_particles_truth in zip(self.df_fjparticles['fj_particles_det'], self.df_fjparticles['fj_particles_truth'])]

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
  def analyzeJets(self, fj_particles_det, fj_particles_truth, jet_def, jet_selector_det, jet_selector_truth_matched):
      
    # Check that the entries exist appropriately
    # (need to check how this can happen -- but it is only a tiny fraction of events)
    if type(fj_particles_det) != fj.vectorPJ or type(fj_particles_truth) != fj.vectorPJ:
      print('fj_particles type mismatch -- skipping event')
      return
    
    if not self.is_pp:
        
        # Get Pb-Pb event
        fj_particles_combined = self.process_io_emb.load_event()
        
        # Set user_index=0 to keep track that these are Pb-Pb tracks
        # Note:
        #   fj_particles is a std::vector<fastjet::PseudoJet>
        #   Iterating over std::vector returns *copy* of object -- modifying the copy doesn't change the original object
        #   So we re-assign each entry...
        #   This is highly inefficient...should swigify (iterate by reference in c++) ...
        for i, track in enumerate(fj_particles_combined):
            track.set_user_index(0)
            fj_particles_combined[i] = track
        
        # Set user_index >0 on the pp-det tracks to keep track that these are from pp-det
        # Here we set user_index=1, but other values will later be used during prong-matching
        for i, track in enumerate(fj_particles_det):
            track.set_user_index(1)
            fj_particles_det[i] = track
            
        # Form the combined det-level event
        [fj_particles_combined.push_back(p) for p in fj_particles_det]

        # Perform constituent subtraction on det-level, if applicable
        fj_particles_combined = self.constituent_subtractor.process_event(fj_particles_combined)
        rho = self.constituent_subtractor.bge_rho.rho()
    
        # Do jet finding
        cs_combined = fj.ClusterSequence(fj_particles_combined, jet_def)
        jets_combined = fj.sorted_by_pt(cs_combined.inclusive_jets())
        jets_combined_selected = jet_selector_det(jets_combined)

    # Do jet finding
    cs_det = fj.ClusterSequence(fj_particles_det, jet_def)
    jets_det_pp = fj.sorted_by_pt(cs_det.inclusive_jets())
    jets_det_pp_selected = jet_selector_det(jets_det_pp)
    
    cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)
    jets_truth = fj.sorted_by_pt(cs_truth.inclusive_jets())
    jets_truth_selected = jet_selector_det(jets_truth)
    jets_truth_selected_matched = jet_selector_truth_matched(jets_truth)

    jetR = jet_def.R()
    
    # Set det-level jet list as appropriate
    jets_det_selected = None
    if self.is_pp:
        jets_det_selected = jets_det_pp_selected
    else:
        jets_det_selected = jets_combined_selected
    
    # Fill det-level jet histograms (before matching)
    for jet_det in jets_det_selected:
      
      # Check additional acceptance criteria
      # skip event if not satisfied -- since first jet in event is highest pt
      if not self.utils.is_det_jet_accepted(jet_det):
        self.hNevents.Fill(0)
        return
      
      self.fill_det_before_matching(jet_det, jetR)
  
    # Fill truth-level jet histograms (before matching)
    for jet_truth in jets_truth_selected:
      self.fill_truth_before_matching(jet_truth, jetR)
  
    # Loop through jets and set jet matching candidates for each jet in user_info
    if self.is_pp:
        [[self.set_matching_candidates(jet_det, jet_truth, jetR) for jet_truth in jets_truth_selected_matched] for jet_det in jets_det_selected]
    else:
        # First fill the combined-to-pp matches, then the pp-to-pp matches
        [[self.set_matching_candidates(jet_det_combined, jet_det_pp, jetR, 'hDeltaR_combined_ppdet_R{}', fill_jet1_matches_only=True) for jet_det_pp in jets_det_pp_selected] for jet_det_combined in jets_det_selected]
        [[self.set_matching_candidates(jet_det_pp, jet_truth, jetR, 'hDeltaR_ppdet_pptrue_R{}') for jet_truth in jets_truth_selected_matched] for jet_det_pp in jets_det_pp_selected]
        
    # Loop through jets and set accepted matches
    if self.is_pp:
        [self.set_matches_pp(jet_det, jetR) for jet_det in jets_det_selected]
    else:
        [self.set_matches_AA(jet_det_combined, jetR) for jet_det_combined in jets_det_selected]
    
    # Loop through jets and fill matching histograms
    result = [self.fill_matching_histograms(jet_det, jetR) for jet_det in jets_det_selected]

    # Loop through jets and fill response if both det and truth jets are unique match
    for sd_setting in self.sd_settings:

      # Get tree writer
      tree_writer = None
      if self.write_tree_output:
        name = 'tree_writer_R{}_{}'.format(self.utils.remove_periods(jetR), sd_label)
        tree_writer = getattr(self, name)

      result = [self.fill_jet_matches(sd_setting, jet_det, jetR, tree_writer) for jet_det in jets_det_selected]

      # Fill the tree
      if self.write_tree_output:
        tree_writer.fill_tree()

  #---------------------------------------------------------------
  # Loop through jets and store matching candidates in user_info
  #---------------------------------------------------------------
  def set_matching_candidates(self, jet1, jet2, jetR, hname='hDeltaR_All_R{}', fill_jet1_matches_only=False):
  
    # Fill histogram of matching distance of all candidates
    deltaR = jet1.delta_R(jet2)
    getattr(self, hname.format(jetR)).Fill(jet1.pt(), deltaR)
  
    # Add a matching candidate to the list if it is within the geometrical cut
    if deltaR < self.jet_matching_distance*jetR:
        self.set_jet_info(jet1, jet2, deltaR)
        if not fill_jet1_matches_only:
            self.set_jet_info(jet2, jet1, deltaR)
  
  #---------------------------------------------------------------
  # Set 'jet_match' as a matching candidate in user_info of 'jet'
  #---------------------------------------------------------------
  def set_jet_info(self, jet, jet_match, deltaR):
  
    # Get/create object to store list of matching candidates
    jet_user_info = None
    if jet.has_user_info():
        jet_user_info = jet.python_info()
    else:
        jet_user_info = jet_info.jet_info()
      
    jet_user_info.matching_candidates.append(jet_match)
    if deltaR < jet_user_info.closest_jet_deltaR:
        jet_user_info.closest_jet = jet_match
        jet_user_info.closest_jet_deltaR = deltaR
          
    jet.set_python_info(jet_user_info)

  #---------------------------------------------------------------
  # Set accepted jet matches for pp case
  #---------------------------------------------------------------
  def set_matches_pp(self, jet_det, jetR):

    h = getattr(self, 'hJetMatchingQA_R{}'.format(jetR))
    h.Fill('all', jet_det.pt(), 1)
  
    if jet_det.has_user_info():
        
        jet_info_det = jet_det.python_info()
        h.Fill('has_matching_candidate', jet_det.pt(), 1)
        
        if len(jet_info_det.matching_candidates) == 1:
            jet_truth = jet_info_det.closest_jet
            
            # Check that the match is unique
            if jet_truth.has_user_info():
                jet_info_truth = jet_truth.python_info()
                if len(jet_info_truth.matching_candidates) == 1:
                    
                    # Set accepted match
                    jet_info_det.match = jet_truth
                    jet_det.set_python_info(jet_info_det)
                    h.Fill('unique_match', jet_det.pt(), 1)

  #---------------------------------------------------------------
  # Set accepted jet matches for Pb-Pb case
  #---------------------------------------------------------------
  def set_matches_AA(self, jet_det_combined, jetR):
  
    set_match = False
    
    h = getattr(self, 'hJetMatchingQA_R{}'.format(jetR))
    h.Fill('all', jet_det_combined.pt(), 1)

    # Check if combined jet has a pp-det match (uniqueness not currently enforced)
    if jet_det_combined.has_user_info():
    
        h.Fill('has_matching_ppdet_candidate', jet_det_combined.pt(), 1)

        jet_info_combined = jet_det_combined.python_info()
        jet_pp_det = jet_info_combined.closest_jet
          
        # Check if the pp-det jet has a pp-truth match
        if jet_pp_det.has_user_info():
        
          jet_info_pp_det = jet_pp_det.python_info()
          jet_pp_truth = jet_info_pp_det.closest_jet
          h.Fill('has_matching_pptrue_candidate', jet_det_combined.pt(), 1)
          set_match = True
          
          # Check that pp-det to pp-true match is unique
          if self.is_match_unique(jet_pp_det):
              h.Fill('has_matching_pptrue_unique_candidate', jet_det_combined.pt(), 1)
          else:
            set_match = False
              
          # Check matching distance between combined jet and pp-truth
          if jet_det_combined.delta_R(jet_pp_truth) < self.jet_matching_distance*jetR:
            h.Fill('deltaR_combined-truth', jet_det_combined.pt(), 1)
          else:
            set_match = False
    
          # Check if >50% of the pp-det tracks are in the Pb-Pb det jet
          mc_fraction = self.mc_fraction(jet_pp_det, jet_det_combined)
          if mc_fraction > self.mc_fraction_threshold:
            h.Fill('mc_fraction', jet_det_combined.pt(), 1)
          else:
            set_match = False
               
    # Set accepted match
    if set_match:
        jet_info_combined.match = jet_pp_truth
        jet_det_combined.set_python_info(jet_info_combined)
        
        h.Fill('passed_all_cuts', jet_det_combined.pt(), 1)

  #---------------------------------------------------------------
  # Return pt-fraction of tracks in jet_pp_det that are contained in jet_det_combined
  #---------------------------------------------------------------
  def mc_fraction(self, jet_pp_det, jet_det_combined):
  
    pt_total = jet_pp_det.pt()
    
    pt_contained = 0.
    for track in jet_det_combined.constituents():
        if track.user_index() > 0:
            pt_contained += track.pt()
            
    return pt_contained/pt_total

  #---------------------------------------------------------------
  # Return whether a jet has a unique match
  #---------------------------------------------------------------
  def is_match_unique(self, jet_pp_det):
  
    if jet_pp_det.has_user_info():

        jet_info_pp_det = jet_pp_det.python_info()
        jet_pp_truth = jet_info_pp_det.closest_jet

        if len(jet_info_pp_det.matching_candidates) == 1:

            if jet_pp_truth.has_user_info():
            
                jet_info_pp_truth = jet_pp_truth.python_info()
                if len(jet_info_pp_truth.matching_candidates) == 1:
                  return True
                  
    return False

  #---------------------------------------------------------------
  # Fill truth jet histograms
  #---------------------------------------------------------------
  def fill_truth_before_matching(self, jet, jetR):
    
    jet_pt = jet.pt()
    for constituent in jet.constituents():
      z = constituent.pt() / jet.pt()
      getattr(self, 'hZ_Truth_R{}'.format(jetR)).Fill(jet.pt(), z)
      
    getattr(self, 'hJetPt_Truth_R{}'.format(jetR)).Fill(jet.pt())

  #---------------------------------------------------------------
  # Fill det jet histograms
  #---------------------------------------------------------------
  def fill_det_before_matching(self, jet, jetR):
    
    jet_pt = jet.pt()
    for constituent in jet.constituents():
      z = constituent.pt() / jet_pt
      getattr(self, 'hZ_Det_R{}'.format(jetR)).Fill(jet_pt, z)

  #---------------------------------------------------------------
  # Loop through jets and fill matching histos
  #---------------------------------------------------------------
  def fill_matching_histograms(self, jet_det, jetR):
      
    if jet_det.has_user_info():
      jet_truth = jet_det.python_info().match
      
      if jet_truth:
        
        jet_pt_det_ungroomed = jet_det.pt()
        jet_pt_truth_ungroomed = jet_truth.pt()
        JES = (jet_pt_det_ungroomed - jet_pt_truth_ungroomed) / jet_pt_truth_ungroomed
        getattr(self, 'hJES_R{}'.format(jetR)).Fill(jet_pt_truth_ungroomed, JES)

  #---------------------------------------------------------------
  # Loop through jets and fill response if both det and truth jets are unique match
  #---------------------------------------------------------------
  def fill_jet_matches(self, sd_setting, jet_det, jetR, tree_writer):
            
    zcut = sd_setting[0]
    beta = sd_setting[1]
    sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)
    
    # Define SoftDrop settings
    # Note: Set custom recluster definition, since by default it uses jetR=max_allowable_R
    sd = fjcontrib.SoftDrop(beta, zcut, jetR)
    jet_def_recluster = fj.JetDefinition(fj.cambridge_algorithm, jetR)
    reclusterer = fjcontrib.Recluster(jet_def_recluster)
    sd.set_reclustering(True, reclusterer)
    if self.debug_level > 2:
      print('SoftDrop groomer is: {}'.format(sd.description()))

    # Check additional acceptance criteria
    # skip event if not satisfied -- since first jet in event is highest pt
    if not self.utils.is_det_jet_accepted(jet_det):
      return

    if jet_det.has_user_info():
      jet_truth = jet_det.python_info().match
      
      if jet_truth:
        self.fill_response(tree_writer, jet_det, jet_truth, sd, jetR, sd_label)

  #---------------------------------------------------------------
  # Fill response histograms
  #---------------------------------------------------------------
  def fill_response(self, tree_writer, jet_det, jet_truth, sd, jetR, sd_label):
    
    # Fill tree
    jet_pt_det_ungroomed = jet_det.pt()
    jet_pt_truth_ungroomed = jet_truth.pt()
    
    theta_g_det = self.theta_g(jet_det, sd, jetR)
    theta_g_truth = self.theta_g(jet_truth, sd, jetR)
    
    zg_det = self.zg(jet_det, sd)
    zg_truth = self.zg(jet_truth, sd)
    
    if self.write_tree_output:
      tree_writer.fill_branch('jet_pt_det_ungroomed', jet_pt_det_ungroomed)
      tree_writer.fill_branch('jet_pt_truth_ungroomed', jet_pt_truth_ungroomed)
      tree_writer.fill_branch('theta_g_det', theta_g_det)
      tree_writer.fill_branch('theta_g_truth', theta_g_truth)
      tree_writer.fill_branch('zg_det', zg_det)
      tree_writer.fill_branch('zg_truth', zg_truth)
    else:
      x = ([jet_pt_det_ungroomed, jet_pt_truth_ungroomed, theta_g_det, theta_g_truth])
      x_array = array('d', x)
      getattr(self, 'hResponse_JetPt_ThetaG_R{}_{}'.format(jetR, sd_label)).Fill(x_array)

      x = ([jet_pt_det_ungroomed, jet_pt_truth_ungroomed, zg_det, zg_truth])
      x_array = array('d', x)
      getattr(self, 'hResponse_JetPt_zg_R{}_{}'.format(jetR, sd_label)).Fill(x_array)

    # Fill response histograms
    theta_g_resolution = (theta_g_det - theta_g_truth) / theta_g_truth
    getattr(self, 'hThetaG_RelativeResidual_JetPt_R{}_{}'.format(jetR, sd_label)).Fill(jet_pt_truth_ungroomed, theta_g_resolution)
    
    theta_g_resolution = theta_g_det - theta_g_truth
    getattr(self, 'hThetaGResidual_JetPt_R{}_{}'.format(jetR, sd_label)).Fill(jet_pt_truth_ungroomed, theta_g_resolution)
      
    zg_resolution = (zg_det - zg_truth) / zg_truth
    getattr(self, 'hZg_RelativeResidual_JetPt_R{}_{}'.format(jetR, sd_label)).Fill(jet_pt_truth_ungroomed, zg_resolution)

    zg_resolution = zg_det - zg_truth
    getattr(self, 'hZgResidual_JetPt_R{}_{}'.format(jetR, sd_label)).Fill(jet_pt_truth_ungroomed, zg_resolution)

  #---------------------------------------------------------------
  # Compute theta_g
  #---------------------------------------------------------------
  def theta_g(self, jet, sd, jetR):
    
    jet_sd = sd.result(jet)
    sd_info = fjcontrib.get_SD_jet_info(jet_sd)
    theta_g = sd_info.dR / jetR
    return theta_g

  #---------------------------------------------------------------
  # Compute z_g
  #---------------------------------------------------------------
  def zg(self, jet, sd):
    
    jet_sd = sd.result(jet)
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

  analysis = process_rg_mc(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_rg_mc()
