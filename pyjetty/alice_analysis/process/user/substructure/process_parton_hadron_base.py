#!/usr/bin/env python3

"""
Base class to read a ROOT TTree of track information
and do jet-finding, and save basic histograms.
  
To use this class, the following should be done:

  - Implement a user analysis class inheriting from this one, such as in user/james/process_mc_XX.py
    You should implement the following functions:
      - initialize_user_output_objects_R()
      - fill_observable_histograms()
      - fill_matched_jet_histograms()
    
  - You should include the following histograms:
      - Response matrix: hResponse_JetPt_[obs]_R[R]_[subobs]_[grooming setting]
      - Residual distribution: hResidual_JetPt_[obs]_R[R]_[subobs]_[grooming setting]

  - You also should modify observable-specific functions at the top of common_utils.py
  
Author: Ezra Lesser (elesser@berkeley.edu)
        based on substructure framework by James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import time  # for timing code
import gc    # garbage collection for cleaning memory

# Data analysis and plotting
import pandas
import numpy as np
from array import *
import ROOT
import yaml
import random

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io_parton_hadron
from pyjetty.alice_analysis.process.base import process_base
from pyjetty.alice_analysis.process.base import common_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class ProcessPHBase(process_base.ProcessBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
  
    # Initialize base class
    super(ProcessPHBase, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
    
    # Initialize configuration
    self.initialize_config()
    
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
    
    # Call base class initialization
    process_base.ProcessBase.initialize_config(self)
    
    # Read config file
    config = None
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    
    self.jet_matching_distance = config['jet_matching_distance']

    # Load levels desired for the various RMs.
    # Should be the form of a list of tuples: [("level1", "level2", "MPIon/off"), ...]
    self.RM_levels = config['response_levels']

    # Can be used to simplify histogram generation / initialization
    if 'theory_pt_bins' in config:
      self.pt_bins = config['theory_pt_bins']
    if 'theory_obs_bins' in config:
      self.obs_bins = config['theory_obs_bins']

    # Create dictionaries to store grooming settings and observable settings for each observable
    # Each dictionary entry stores a list of subconfiguration parameters
    #   The observable list stores the observable setting, e.g. subjetR
    #   The grooming list stores a list of grooming settings {'sd': [zcut, beta]} or {'dg': [a]}
    self.observable_list = config['process_observables']
    self.obs_settings = {}
    self.obs_grooming_settings = {}
    for observable in self.observable_list:
    
      obs_config_dict = config[observable]
      obs_config_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
      
      obs_subconfig_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
      self.obs_settings[observable] = self.utils.obs_settings(
        observable, obs_config_dict, obs_subconfig_list)
      self.obs_grooming_settings[observable] = self.utils.grooming_settings(obs_config_dict)
      
    # Construct set of unique grooming settings
    self.grooming_settings = []
    lists_grooming = [self.obs_grooming_settings[obs] for obs in self.observable_list]
    for observable in lists_grooming:
      for setting in observable:
        if setting not in self.grooming_settings and setting != None:
          self.grooming_settings.append(setting)
  
  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def process_mc(self):
    
    self.start_time = time.time()
    
    # ------------------------------------------------------------------------
    
    # Use IO helper class to convert parton- & hadron-level ROOT TTree into
    # SeriesGroupBy objects of fastjet particles per event
    print('--- {} seconds ---'.format(time.time() - self.start_time))

    # Load MPI-off trees into Pandas DataFrames
    io_p_MPIoff = process_io_parton_hadron.ProcessIO(input_file=self.input_file, level='p', MPI='off')
    io_h_MPIoff = process_io_parton_hadron.ProcessIO(input_file=self.input_file, level='h', MPI='off')

    df_fjparticles_p_MPIoff = io_p_MPIoff.load_data()
    df_fjparticles_h_MPIoff = io_h_MPIoff.load_data()
    df_fjparticles_ch_MPIoff = io_h_MPIoff.group_fjparticles(ch_cut=True)

    # Save number of events and cross section information for proper scaling
    self.nEvents_MPIoff = io_p_MPIoff.Nev
    self.xsec_MPIoff = io_p_MPIoff.xsec

    # Explicitly clean up memory since files can be pretty large
    del io_p_MPIoff
    del io_h_MPIoff
    gc.collect()

    # Repeat for MPI-on trees
    io_p_MPIon = process_io_parton_hadron.ProcessIO(input_file=self.input_file, level='p', MPI='on')
    io_h_MPIon = process_io_parton_hadron.ProcessIO(input_file=self.input_file, level='h', MPI='on')

    df_fjparticles_p_MPIon = io_p_MPIon.load_data()
    df_fjparticles_h_MPIon = io_h_MPIon.load_data()
    df_fjparticles_ch_MPIon = io_h_MPIon.group_fjparticles(ch_cut=True)

    self.nEvents_MPIon = io_p_MPIon.Nev
    self.xsec_MPIon = io_p_MPIon.xsec

    del io_p_MPIon
    del io_h_MPIon
    gc.collect()

    print('--- {} seconds ---'.format(time.time() - self.start_time))
    
    # ------------------------------------------------------------------------

    # Now merge the SeriesGroupBy to create a groupby df with [ev_id, run_number, fj_p, fj_h, fj_ch]
    # (Need a structure such that we can iterate event-by-event through both fj_p, fj_h, fj_ch simultaneously)
    print('Merge parton-, hadron-, and charged-level into a single dataframe grouped by event...')

    self.df_fjparticles_MPIoff = pandas.concat(
      [df_fjparticles_p_MPIoff, df_fjparticles_h_MPIoff, df_fjparticles_ch_MPIoff], axis=1)
    self.df_fjparticles_MPIon = pandas.concat(
      [df_fjparticles_p_MPIon, df_fjparticles_h_MPIon, df_fjparticles_ch_MPIon], axis=1)

    self.df_fjparticles_MPIoff.columns = ['fj_particles_p', 'fj_particles_h', 'fj_particles_ch']
    self.df_fjparticles_MPIon.columns = ['fj_particles_p', 'fj_particles_h', 'fj_particles_ch']

    print('--- {} seconds ---'.format(time.time() - self.start_time))
    
    # ------------------------------------------------------------------------

    # Initialize histograms
    self.initialize_output_objects()
    
    print(self)
    
    # Find jets and fill histograms
    for MPI in ["off", "on"]:
      print('Find jets with MPI %s...' % MPI)
      self.analyze_events(MPI)

    print('--- {} seconds ---'.format(time.time() - self.start_time))
    
    # ------------------------------------------------------------------------

    print('Save thn...')
    process_base.ProcessBase.save_thn_th3_objects(self)

    # Plot histograms
    print('Save histograms...')
    process_base.ProcessBase.save_output_objects(self)
    
    print('--- {} seconds ---'.format(time.time() - self.start_time))

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_output_objects(self):

    for MPI in ["on", "off"]:

      setattr(self, "hNevents_MPI%s" % MPI, ROOT.TH1F(
        'hNevents_MPI%s' % MPI, 'hNevents_MPI%s' % MPI, 2, -0.5, 1.5))
      getattr(self, "hNevents_MPI%s" % MPI).Fill(1, getattr(self, "nEvents_MPI%s" % MPI))

      for jetR in self.jetR_list:
        self.initialize_output_objects_R(jetR, MPI)

      for level in ["p", "h", "ch"]:

        name_suffix = "%s_MPI%s" % (level, MPI)

        setattr(self, "hTrackEtaPhi_%s" % name_suffix, ROOT.TH2F(
          'hTrackEtaPhi_%s' % name_suffix, 'hTrackEtaPhi_%s' % name_suffix,
          200, -1., 1., 628, 0., 6.28))

        setattr(self, "hTrackPt_%s" % name_suffix, ROOT.TH1F(
          'hTrackPt_%s' % name_suffix, 'hTrackPt_%s' % name_suffix, 300, 0., 300.))

    #self.hN_MeanPt = ROOT.TH2F('hN_MeanPt', 'hN_MeanPt', 200, 0, 5000, 200, 0., 2.)
    for jetR in self.jetR_list:
      # Call user-specific initialization
      self.initialize_user_output_objects_R(jetR)

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_output_objects_R(self, jetR, MPI):

      name_suffix = "MPI%s_R%s" % (MPI, str(jetR))

      # Base histograms
      # Jet energy scale (between the 3 different levels)
      name = 'hJES_p_h_%s' % name_suffix
      h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
      setattr(self, name, h)

      name = 'hJES_h_ch_%s' % name_suffix
      h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
      setattr(self, name, h)

      name = 'hJES_p_ch_%s' % name_suffix
      h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
      setattr(self, name, h)

      # Delta R (between p-h and h-ch)
      name = 'hDeltaR_All_p_h_%s' % name_suffix
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
      setattr(self, name, h)

      name = 'hDeltaR_All_h_ch_%s' % name_suffix
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
      setattr(self, name, h)

      # Create matching QA histograms
      hname_p = "hJetMatchingQA_p_h_MPI%s_R%s" % (MPI, str(jetR))
      hname_h = "hJetMatchingQA_h_ch_MPI%s_R%s" % (MPI, str(jetR))
      for hname in (hname_p, hname_h):
        bin_labels = ['all', 'has_matching_candidate', 'unique_match']
        nbins = len(bin_labels)
        h = ROOT.TH2F(hname, hname, nbins, 0, nbins, 30, 0., 300.)
        for i in range(1, nbins+1):
          h.GetXaxis().SetBinLabel(i, bin_labels[i-1])
        setattr(self, hname, h)

      # z (for each of the 3 different levels)
      for level in ["p", "h", "ch"]:
        name = 'hZ_%s_%s' % (level, name_suffix)
        h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 1.)
        setattr(self, name, h)

  #---------------------------------------------------------------
  # Main function to loop through and analyze events
  #---------------------------------------------------------------
  def analyze_events(self, MPI="on"):

    df = getattr(self, "df_fjparticles_MPI%s" % MPI)

    # Fill track histograms
    for level in ['p', 'h', 'ch']:
      for fj_particles in df['fj_particles_%s' % level]:
        self.fill_track_histograms(fj_particles, level, MPI)
    
    fj.ClusterSequence.print_banner()
    print()
        
    self.event_number = 0

    # match jets -- and fill histograms
    for fj_particles_p, fj_particles_h, fj_particles_ch in zip(
        df['fj_particles_p'], df['fj_particles_h'], df['fj_particles_ch']):

      self.analyze_event(fj_particles_p, fj_particles_h, fj_particles_ch, MPI)

    if self.debug_level > 0:
      for attr in dir(self):
        obj = getattr(self, attr)
        print('size of {}: {}'.format(attr, sys.getsizeof(obj)))
    
  #---------------------------------------------------------------
  # Fill track histograms.
  #---------------------------------------------------------------
  def fill_track_histograms(self, fj_particles, level, MPI):

    # Check that the entries exist appropriately
    # (need to check if indeed this happens for MC, and how...)
    if type(fj_particles) != fj.vectorPJ:
      return

    for track in fj_particles:
      getattr(self, "hTrackEtaPhi_%s_MPI%s" % (level, MPI)).Fill(track.eta(), track.phi())
      getattr(self, "hTrackPt_%s_MPI%s" % (level, MPI)).Fill(track.pt())
      
  #---------------------------------------------------------------
  # Analyze jets of a given event.
  # fj_particles is the list of fastjet pseudojets for a single fixed event.
  #---------------------------------------------------------------
  def analyze_event(self, fj_particles_p, fj_particles_h, fj_particles_ch, MPI):

    self.event_number += 1
    if self.event_number > self.event_number_max:
      return
    if self.debug_level > 1:
      print('-------------------------------------------------')
      print('event %i' % self.event_number)

    # Check that the entries exist appropriately
    # (need to check how this can happen -- but it is only a tiny fraction of events)
    if type(fj_particles_p) != fj.vectorPJ or type(fj_particles_h) != fj.vectorPJ or \
       type(fj_particles_ch) != fj.vectorPJ:
      print('fj_particles type mismatch -- skipping event')
      return

    if self.debug_level > 0 and len(fj_particles_h) > 1:
      if np.abs(fj_particles_h[0].pt() - fj_particles_h[1].pt()) <  1e-10:
        print('WARNING: Duplicate particles may be present')
        print([p.user_index() for p in fj_particles_h])
        print([p.pt() for p in fj_particles_h])

    # Loop through jetR, and process event for each R
    for jetR in self.jetR_list:
    
      # Keep track of whether to fill R-independent histograms
      self.fill_R_indep_hists = (jetR == self.jetR_list[0])

      # Do not consider tracks with pT below 150 MeV/c at ch-level
      track_selector_ch = fj.SelectorPtMin(0.15)
      # Set jet definition and a jet selector
      jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
      jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
      if self.debug_level > 2:
        print('')
        print('track selector for charged level is', track_selector_ch)
        print('jet definition is:', jet_def)
        print('jet selector is:', jet_selector)

      # Find jets at p, h, and ch level
      cs_p = fj.ClusterSequence(fj_particles_p, jet_def)
      jets_p = fj.sorted_by_pt(cs_p.inclusive_jets())
      jets_p_selected = jet_selector(jets_p)

      cs_h = fj.ClusterSequence(fj_particles_h, jet_def)
      jets_h = fj.sorted_by_pt(cs_h.inclusive_jets())
      jets_h_selected = jet_selector(jets_h)

      #fj_part_ch_selected = track_selector_ch(fj_particles_ch)
      cs_ch = fj.ClusterSequence(track_selector_ch(fj_particles_ch), jet_def)
      jets_ch = fj.sorted_by_pt(cs_ch.inclusive_jets())
      jets_ch_selected = jet_selector(jets_ch)

      self.analyze_jets(jets_p_selected, jets_h_selected, jets_ch_selected, jetR, MPI)

  #---------------------------------------------------------------
  # Analyze jets of a given event.
  #---------------------------------------------------------------
  def analyze_jets(self, jets_p_selected, jets_h_selected, jets_ch_selected, jetR, MPI):

    if self.debug_level > 1:
      print('Number of parton-level jets: %i' % len(jets_p_selected))
      print('Number of hadron-level jets: %i' % len(jets_h_selected))
      print('Number of charged-level jets: %i' % len(jets_ch_selected))

    # Fill parton-, hadron- and charged-level jet histograms (before matching)
    for jet_p in jets_p_selected:
      self.fill_level_before_matching(jet_p, jetR, 'p', MPI)
    for jet_h in jets_h_selected:
      self.fill_level_before_matching(jet_h, jetR, 'h', MPI)
    for jet_ch in jets_h_selected:
      self.fill_level_before_matching(jet_ch, jetR, 'ch', MPI)
  
    # Loop through jets and set jet matching candidates for each jet in user_info
    name_suffix = "MPI%s_R%s" % (MPI, str(jetR))
    for jet_p in jets_p_selected:
      for jet_h in jets_h_selected:
        self.set_matching_candidates(jet_p, 'p', jet_h, 'h', jetR, MPI)
        for jet_ch in jets_ch_selected:
          self.set_matching_candidates(jet_h, 'h', jet_ch, 'ch', jetR, MPI)

    # Loop through jets and set accepted matches
    for jet_h in jets_h_selected:
      self.set_matches(jet_h, jetR, MPI)

    # Loop through jets and fill response histograms if both det and truth jets are unique match
    for jet_h in jets_h_selected:
      self.fill_jet_matches(jet_h, jetR, MPI)

  #---------------------------------------------------------------
  # Compare two jets and store matching candidates in user_info
  #---------------------------------------------------------------
  def set_matching_candidates(self, jet1, level1, jet2, level2, jetR, MPI):

    # Fill histogram of matching distance of all candidates
    deltaR = jet1.delta_R(jet2)
    getattr(self, "hDeltaR_All_%s_%s_MPI%s_R%s" % (
      level1, level2, MPI, str(jetR))).Fill(jet1.pt(), deltaR)
  
    # Add a matching candidate to the list if it is within the geometrical cut
    if deltaR < self.jet_matching_distance * jetR:
      self.set_jet_info(jet1, level1, jet2, level2, deltaR)
      self.set_jet_info(jet2, level2, jet1, level1, deltaR)

  #---------------------------------------------------------------
  # Set 'jet_match' as a matching candidate in user_info of 'jet'
  #---------------------------------------------------------------
  def set_jet_info(self, jet1, level1, jet2, level2, deltaR):
  
    # Get/create object to store list of matching candidates
    jet_user_info = None
    if jet1.has_user_info():
      jet_user_info = jet1.python_info()
    else:
      jet_user_info = JetInfo()

    if 'p' in (level1, level2):  # assume matching parton to hadron
      jet_user_info.matching_candidates_p.append(jet2)
      if deltaR < jet_user_info.closest_jet_deltaR_p:
        jet_user_info.closest_jet_p = jet2
        jet_user_info.closest_jet_deltaR_p = deltaR

    else:  # assume matching hadron to charged
      jet_user_info.matching_candidates_ch.append(jet2)
      if deltaR < jet_user_info.closest_jet_deltaR_ch:
        jet_user_info.closest_jet_ch = jet2
        jet_user_info.closest_jet_deltaR_ch = deltaR
          
    jet1.set_python_info(jet_user_info)

  #---------------------------------------------------------------
  # Set accepted jet matches
  #---------------------------------------------------------------
  def set_matches(self, jet_h, jetR, MPI):

    # For filling matching QA histograms
    hname_p = "hJetMatchingQA_p_h_MPI%s_R%s" % (MPI, str(jetR))
    hname_ch = "hJetMatchingQA_h_ch_MPI%s_R%s" % (MPI, str(jetR))

    h_p = getattr(self, hname_p)
    h_p.Fill('all', jet_h.pt(), 1)
    h_ch = getattr(self, hname_ch)
    h_ch.Fill('all', jet_h.pt(), 1)
  
    if jet_h.has_user_info():
        
      jet_info_h = jet_h.python_info()

      if len(jet_info_h.matching_candidates_p) > 0:
        h_p.Fill('has_matching_candidate', jet_h.pt(), 1)
      if len(jet_info_h.matching_candidates_ch) > 0:
        h_ch.Fill('has_matching_candidate', jet_h.pt(), 1)

      # Match with parton-level jet
      if len(jet_info_h.matching_candidates_p) == 1:
        jet_p = jet_info_h.closest_jet_p
          
        # Check that the match is unique
        if jet_p.has_user_info():
          jet_info_p = jet_p.python_info()
          if len(jet_info_p.matching_candidates_p) == 1:
                  
            # Set accepted match
            jet_info_h.match_p = jet_p
            jet_h.set_python_info(jet_info_h)
            h_p.Fill('unique_match', jet_h.pt(), 1)

      # Match with charged-level jet
      if len(jet_info_h.matching_candidates_ch) == 1:
        jet_ch = jet_info_h.closest_jet_ch
          
        # Check that the match is unique
        if jet_ch.has_user_info():
          jet_info_ch = jet_ch.python_info()
          if len(jet_info_ch.matching_candidates_ch) == 1:
                  
            # Set accepted match
            jet_info_h.match_ch = jet_ch
            jet_h.set_python_info(jet_info_h)
            h_ch.Fill('unique_match', jet_h.pt(), 1)

  #---------------------------------------------------------------
  # Fill det jet histograms
  #---------------------------------------------------------------
  def fill_level_before_matching(self, jet, jetR, level, MPI):

    jet_pt = jet.pt()
    label = '%s_MPI%s_R%s' % (level, MPI, jetR)
    for constituent in jet.constituents():
      z = constituent.pt() / jet_pt
      getattr(self, 'hZ_%s' % label).Fill(jet_pt, z)

    self.fill_unmatched_jet_histograms(jet, jetR, level, MPI)
  
  #---------------------------------------------------------------
  # This function is called once for each jet
  #---------------------------------------------------------------
  def fill_unmatched_jet_histograms(self, jet, jetR, level, MPI):

    # Loop through each jet subconfiguration (i.e. subobservable / grooming setting)
    observable = self.observable_list[0]
    for i in range(len(self.obs_settings[observable])):

      obs_setting = self.obs_settings[observable][i]
      grooming_setting = self.obs_grooming_settings[observable][i]
      obs_label = self.utils.obs_label(obs_setting, grooming_setting)

      # Groom jet, if applicable
      if grooming_setting:
        gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
        jet_groomed_lund = self.utils.groom(gshop, grooming_setting, jetR)
        if not jet_groomed_lund:
          continue
      else:
        jet_groomed_lund = None
        
      # Call user function to fill histograms
      self.fill_observable_histograms(jet, jet_groomed_lund, jetR, level, MPI,
                                      obs_setting, grooming_setting, obs_label, jet.pt())
  
  #---------------------------------------------------------------
  # Loop through jets and call user function to fill matched
  # histos if both det and truth jets are unique match.
  #---------------------------------------------------------------
  def fill_jet_matches(self, jet_h, jetR, MPI):

    suffix = "MPI%s_R%s" % (MPI, str(jetR))

    # Get matched jets
    if jet_h.has_user_info():
      jet_p = jet_h.python_info().match_p
      jet_ch = jet_h.python_info().match_ch

      if jet_p and jet_ch:

        jet_pt_p_ungroomed = jet_p.pt()
        jet_pt_h_ungroomed = jet_h.pt()
        jet_pt_ch_ungroomed = jet_ch.pt()

        JES = (jet_pt_h_ungroomed - jet_pt_p_ungroomed) / jet_pt_p_ungroomed
        getattr(self, 'hJES_p_h_%s' % suffix).Fill(jet_pt_p_ungroomed, JES)
        JES = (jet_pt_ch_ungroomed - jet_pt_h_ungroomed) / jet_pt_h_ungroomed
        getattr(self, 'hJES_h_ch_%s' % suffix).Fill(jet_pt_h_ungroomed, JES)
        JES = (jet_pt_ch_ungroomed - jet_pt_p_ungroomed) / jet_pt_p_ungroomed
        getattr(self, 'hJES_p_ch_%s' % suffix).Fill(jet_pt_p_ungroomed, JES)

        # Loop through each jet subconfiguration (i.e. subobservable / grooming setting)
        observable = self.observable_list[0]
        for i in range(len(self.obs_settings[observable])):

          obs_setting = self.obs_settings[observable][i]
          grooming_setting = self.obs_grooming_settings[observable][i]
          obs_label = self.utils.obs_label(obs_setting, grooming_setting)

          if self.debug_level > 3:
            print('obs_label: {}'.format(obs_label))

          jet_p_groomed_lund = None
          jet_h_groomed_lund = None
          jet_ch_groomed_lund = None

          # Groom jets, if applicable
          if grooming_setting:

            # Groom p jet
            gshop_p = fjcontrib.GroomerShop(jet_p, jetR, self.reclustering_algorithm)
            jet_p_groomed_lund = self.utils.groom(gshop_p, grooming_setting, jetR)
            if not jet_p_groomed_lund:
              continue

            # Groom h jet
            gshop_h = fjcontrib.GroomerShop(jet_h, jetR, self.reclustering_algorithm)
            jet_h_groomed_lund = self.utils.groom(gshop_h, grooming_setting, jetR)
            if not jet_h_groomed_lund:
              continue

            # Groom ch jet
            gshop_ch = fjcontrib.GroomerShop(jet_ch, jetR, self.reclustering_algorithm)
            jet_ch_groomed_lund = self.utils.groom(gshop_ch, grooming_setting, jetR)
            if not jet_ch_groomed_lund:
              continue

          # Call user function to fill histos
          self.fill_matched_jet_histograms(
            jetR, obs_setting, grooming_setting, obs_label, jet_p, jet_p_groomed_lund,
            jet_h, jet_h_groomed_lund, jet_ch, jet_ch_groomed_lund,
            jet_pt_p_ungroomed, jet_pt_h_ungroomed, jet_pt_ch_ungroomed, suffix)

  #---------------------------------------------------------------
  # Initialize MPI scaling histograms -- common utility function
  # NOTE: if used in the config, you can pass in self.pt/obs_bins for simplicity
  #---------------------------------------------------------------
  def init_MPI_scaling_hist(self, observable, obs_name, level, jetR, pt_bins, obs_bins, obs_label):

    for MPI in ["on", "off"]:
      name = "h_JetPt_%s_%s_MPI%s_R%s_%s" % (observable, level, MPI, jetR, obs_label)
      h = ROOT.TH2F(name, name, len(pt_bins)-1, array('d', pt_bins),
                    len(obs_bins)-1, array('d', obs_bins))
      h.GetXaxis().SetTitle("p_{T}^{%s jet}" % level)
      h.GetYaxis().SetTitle("#frac{dN}{d%s_{%s}^{%s}}" % (obs_name, obs_label, level))
      h.Sumw2()
      setattr(self, name, h)

  #---------------------------------------------------------------
  # Initialize MPI scaling histograms -- common utility function
  # NOTE: if used in the config, you can pass in self.pt/obs_bins for simplicity
  #---------------------------------------------------------------
  def fill_MPI_scaling_hist(self, observable, level, MPI, jetR, jet_pt, obs, obs_label):

    getattr(self, "h_JetPt_%s_%s_MPI%s_R%s_%s" % \
            (observable, level, MPI, jetR, obs_label)).Fill(jet_pt, obs)

  #---------------------------------------------------------------
  # Initialize response histogram -- common utility function
  # NOTE: if used in the config, you can pass in self.pt/obs_bins for simplicity
  #---------------------------------------------------------------
  def init_response(self, observable, obs_name, level_1, level_2, MPI,
                    jetR, pt_bins, obs_bins, obs_label):

    # Another set of THn for full hadron folding
    title = ["p_{T}^{%s jet}" % level_2, "p_{T}^{%s jet}" % level_1, 
             "%s_{%s}^{%s}" % (obs_name, obs_label, level_2),
             "%s_{%s}^{%s}" % (obs_name, obs_label, level_1)]

    label = "%s_%s_%s_MPI%s_R%s_%s" % (observable, level_1, level_2, MPI, jetR, obs_label)
    name = "hResponse_JetPt_%s" % label
    nbins_array = array('i', [len(pt_bins)-1, len(pt_bins)-1, len(obs_bins)-1, len(obs_bins)-1])
    xmin_array = array('d', [pt_bins[0], pt_bins[0], obs_bins[0], obs_bins[0]])
    xmax_array = array('d', [pt_bins[-1], pt_bins[-1], obs_bins[-1], obs_bins[-1]])

    # assume 4 dimmensions
    h = ROOT.THnF(name, name, 4, nbins_array, xmin_array, xmax_array)
    for i in range(0, 4):
      h.GetAxis(i).SetTitle(title[i])
      if i == 0 or i == 1:
        h.SetBinEdges(i, array('d', pt_bins))
      else:  # i == 2 or i == 3
        h.SetBinEdges(i, array('d', obs_bins))
    h.Sumw2()
    setattr(self, name, h)

    name = "hResidual_JetPt_%s" % label
    h = ROOT.TH3F(name, name, 20, 0, 200, 100, 0., 1., 200, -2., 2.)
    h.GetXaxis().SetTitle("p_{T, truth}^{%s jet}" % level_1)
    h.GetYaxis().SetTitle("%s_{%s}^{%s}" % (obs_name, obs_label, level_1))
    h.GetZaxis().SetTitle("#frac{%s^{%s}-%s^{%s}}{%s^{%s}}" % \
                          (obs_name, level_2, obs_name, level_1, obs_name, level_1))
    setattr(self, name, h)

  #---------------------------------------------------------------
  # Fill response histograms -- common utility function
  #---------------------------------------------------------------
  def fill_response(self, observable, level_1, level_2, MPI, jetR, jet_pt_level_1, jet_pt_level_2,
                    obs_level_1, obs_level_2, obs_label):

    label = "%s_%s_%s_MPI%s_R%s_%s" % (observable, level_1, level_2, MPI, jetR, obs_label)

    x_array = array('d', [jet_pt_level_2, jet_pt_level_1, obs_level_2, obs_level_1])
    getattr(self, "hResponse_JetPt_"+label).Fill(x_array)
      
    if obs_level_1 > 0:
      obs_resolution = (obs_level_2 - obs_level_1) / obs_level_1
      getattr(self, "hResidual_JetPt_"+label).Fill(jet_pt_level_1, obs_level_1, obs_resolution)

  #---------------------------------------------------------------
  # This function is called once for each jetR
  # You must implement this
  #---------------------------------------------------------------
  def initialize_user_output_objects_R(self, jetR, MPI):
      
    raise NotImplementedError('You must implement initialize_user_output_objects_R()!')

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  # You must implement this
  #---------------------------------------------------------------
  def fill_observable_histograms(self, jet, jet_groomed_lund, jetR, level, MPI,
                                 obs_setting, grooming_setting, obs_label, jet_pt_ungroomed):

    raise NotImplementedError('You must implement fill_observable_histograms()!')

  #---------------------------------------------------------------
  # This function is called once for each matched jet subconfiguration
  # You must implement this
  #---------------------------------------------------------------
  def fill_matched_jet_histograms(
      self, jetR, obs_setting, grooming_setting, obs_label,
      jet_p, jet_p_groomed_lund, jet_h, jet_h_groomed_lund, jet_ch, jet_ch_groomed_lund,
      jet_pt_p_ungroomed, jet_pt_h_ungroomed, jet_pt_ch_ungroomed, suffix):

    raise NotImplementedError('You must implement fill_matched_jet_histograms()!')




################################################################
class JetInfo(common_base.CommonBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(JetInfo, self).__init__(**kwargs)

    # Store the parton-to-hadron matching candidates
    self.matching_candidates_p = []
    self.closest_jet_p = None
    self.closest_jet_deltaR_p = 1000.
    self.match_p = None

    # Store the hadron-to-charged matching candidates
    self.matching_candidates_ch = []
    self.closest_jet_ch = None
    self.closest_jet_deltaR_ch = 1000.
    self.match_ch = None
