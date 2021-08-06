#!/usr/bin/env python3

"""
Base class to read a ROOT TTree of track information
and do jet-finding, and save basic histograms.

To use this class, the following should be done:

  - Implement a user analysis class inheriting from this one, such as in user/james/process_data_XX.py
    You should implement the following functions:
      - initialize_user_output_objects()
      - fill_jet_histograms()

  - The histogram of the data should be named h_[obs]_JetPt_R[R]_[subobs]_[grooming setting]
    The grooming part is optional, and should be labeled e.g. zcut01_B0 — from CommonUtils::grooming_label({'sd':[zcut, beta]})
    For example: h_subjet_z_JetPt_R0.4_0.1
    For example: h_subjet_z_JetPt_R0.4_0.1_zcut01_B0

  - You also should modify observable-specific functions at the top of common_utils.py

Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import time

# Data analysis and plotting
import numpy as np
import ROOT
import yaml

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_base
from pyjetty.mputils import CEventSubtractor

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class ProcessDataBase(process_base.ProcessBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):

    # Initialize base class
    super(ProcessDataBase, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)

    # Initialize configuration
    self.initialize_config()

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
  def process_data(self):

    self.start_time = time.time()

    # Use IO helper class to convert ROOT TTree into a SeriesGroupBy object of
    # fastjet particles per event
    print('--- {} seconds ---'.format(time.time() - self.start_time))
    io = process_io.ProcessIO(input_file=self.input_file, track_tree_name='tree_Particle',
                              is_pp=self.is_pp, use_ev_id_ext=True)
    self.df_fjparticles = io.load_data(m=self.m)
    self.nEvents = len(self.df_fjparticles.index)
    self.nTracks = len(io.track_df.index)
    print('--- {} seconds ---'.format(time.time() - self.start_time))

    # Initialize histograms
    self.initialize_output_objects()

    # Create constituent subtractor, if configured
    if not self.is_pp:
      max_dist_li = self.max_distance if isinstance(self.max_distance, list) else \
                    list(np.unique(np.concatenate(list(self.max_distance.values()))))
      self.constituent_subtractor = { R_max : CEventSubtractor(
        max_distance=R_max, alpha=self.alpha, max_eta=self.max_eta,
        bge_rho_grid_size=self.bge_rho_grid_size, max_pt_correct=self.max_pt_correct,
        ghost_area=self.ghost_area, distance_type=fjcontrib.ConstituentSubtractor.deltaR) \
                                      for R_max in max_dist_li}

    print(self)

    # Find jets and fill histograms
    print('Analyze events...')
    self.analyze_events()

    # Plot histograms
    print('Save histograms...')
    process_base.ProcessBase.save_output_objects(self)

    print('--- {} seconds ---'.format(time.time() - self.start_time))

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_output_objects(self):

    # Initialize user-specific histograms
    self.initialize_user_output_objects()

    # Initialize base histograms
    self.hNevents = ROOT.TH1F('hNevents', 'hNevents', 2, -0.5, 1.5)
    if self.event_number_max < self.nEvents:
      self.hNevents.Fill(1, self.event_number_max)
    else:
      self.hNevents.Fill(1, self.nEvents)

    self.hTrackEtaPhi = ROOT.TH2F('hTrackEtaPhi', 'hTrackEtaPhi', 200, -1., 1., 628, 0., 6.28)
    self.hTrackPt = ROOT.TH1F('hTrackPt', 'hTrackPt', 300, 0., 300.)

    if not self.is_pp:
      self.hRho = ROOT.TH1F('hRho', 'hRho', 1000, 0., 1000.)

    for jetR in self.jetR_list:

      name = 'hZ_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 1.)
      setattr(self, name, h)

  #---------------------------------------------------------------
  # Main function to loop through and analyze events
  #---------------------------------------------------------------
  def analyze_events(self):

    # Fill track histograms
    print('--- {} seconds ---'.format(time.time() - self.start_time))
    print('Fill track histograms')
    [[self.fillTrackHistograms(track) for track in fj_particles] for fj_particles in self.df_fjparticles]
    print('--- {} seconds ---'.format(time.time() - self.start_time))

    print('Find jets...')
    fj.ClusterSequence.print_banner()
    print()
    self.event_number = 0

    # Use list comprehension to do jet-finding and fill histograms
    result = [self.analyze_event(fj_particles) for fj_particles in self.df_fjparticles]

    print('--- {} seconds ---'.format(time.time() - self.start_time))
    print('Save thn...')
    process_base.ProcessBase.save_thn_th3_objects(self)

  #---------------------------------------------------------------
  # Fill track histograms.
  #---------------------------------------------------------------
  def fillTrackHistograms(self, track):

    self.hTrackEtaPhi.Fill(track.eta(), track.phi())
    self.hTrackPt.Fill(track.pt())

  #---------------------------------------------------------------
  # Analyze jets of a given event.
  # fj_particles is the list of fastjet pseudojets for a single fixed event.
  #---------------------------------------------------------------
  def analyze_event(self, fj_particles):

    self.event_number += 1
    if self.event_number > self.event_number_max:
      return
    if self.debug_level > 1:
      print('-------------------------------------------------')
      print('event {}'.format(self.event_number))

    if len(fj_particles) > 1:
      if np.abs(fj_particles[0].pt() - fj_particles[1].pt()) <  1e-10:
        print('WARNING: Duplicate particles may be present')
        print([p.user_index() for p in fj_particles])
        print([p.pt() for p in fj_particles])

    # Perform constituent subtraction for each R_max (do this once, for all jetR)
    if not self.is_pp:
      fj_particles_subtracted = { R_max : cs.process_event(fj_particles) for \
                                  R_max, cs in self.constituent_subtractor.items() }

    # Loop through jetR, and process event for each R
    for jetR in self.jetR_list:

      # Keep track of whether to fill R-independent histograms
      self.fill_R_indep_hists = (jetR == self.jetR_list[0])

      # Set jet definition and a jet selector
      jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
      jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
      if self.debug_level > 2:
        print('jet definition is:', jet_def)
        print('jet selector is:', jet_selector,'\n')

      # Analyze
      if self.is_pp:

        # Do jet finding
        cs = fj.ClusterSequence(fj_particles, jet_def)
        jets = fj.sorted_by_pt(cs.inclusive_jets())
        jets_selected = jet_selector(jets)

        self.analyze_jets(jets_selected, jetR)

      else:

        max_distance = self.max_distance if isinstance(self.max_distance, list) else \
                       self.max_distance[jetR]

        for R_max in max_distance:

          if self.debug_level > 1:
            print('R_max: {}'.format(R_max))

          # Keep track of whether to fill R_max-independent histograms
          self.fill_Rmax_indep_hists = (R_max == max_distance[0])

          # Perform constituent subtraction
          rho = self.constituent_subtractor[R_max].bge_rho.rho()
          if self.fill_R_indep_hists and self.fill_Rmax_indep_hists:
            getattr(self, 'hRho').Fill(rho)

          # Do jet finding (re-do each time, to make sure matching info gets reset)
          cs = fj.ClusterSequence(fj_particles_subtracted[R_max], jet_def)
          jets = fj.sorted_by_pt(cs.inclusive_jets())
          jets_selected = jet_selector(jets)

          self.analyze_jets(jets_selected, jetR, R_max = R_max)

  #---------------------------------------------------------------
  # Analyze jets of a given event.
  #---------------------------------------------------------------
  def analyze_jets(self, jets_selected, jetR, R_max = None):

    # Set suffix for filling histograms
    if R_max:
      suffix = '_Rmax{}'.format(R_max)
    else:
      suffix = ''

    # Loop through jets and call user function to fill histos
    result = [self.analyze_accepted_jet(jet, jetR, suffix) for jet in jets_selected]

  #---------------------------------------------------------------
  # Fill histograms
  #---------------------------------------------------------------
  def analyze_accepted_jet(self, jet, jetR, suffix):

    # Check additional acceptance criteria
    if not self.utils.is_det_jet_accepted(jet):
      return

    # Fill base histograms
    jet_pt_ungroomed = jet.pt()
    if self.is_pp or self.fill_Rmax_indep_hists:

      hZ = getattr(self, 'hZ_R{}'.format(jetR))
      for constituent in jet.constituents():
        z = constituent.pt() / jet_pt_ungroomed
        hZ.Fill(jet_pt_ungroomed, z)

    # Loop through each jet subconfiguration (i.e. subobservable / grooming setting)
    # Note that the subconfigurations are defined by the first observable, if multiple are defined
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
      self.fill_jet_histograms(jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
                               obs_label, jet_pt_ungroomed, suffix)

  #---------------------------------------------------------------
  # This function is called once
  # You must implement this
  #---------------------------------------------------------------
  def initialize_user_output_objects(self):

    raise NotImplementedError('You must implement initialize_user_output_objects()!')

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  # You must implement this
  #---------------------------------------------------------------
  def fill_jet_histograms(self, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
                          obs_label, jet_pt_ungroomed, suffix):

    raise NotImplementedError('You must implement fill_jet_histograms()!')
