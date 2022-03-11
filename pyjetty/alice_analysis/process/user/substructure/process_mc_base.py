#!/usr/bin/env python3

"""
Base class to read a ROOT TTree of track information
and do jet-finding, and save basic histograms.

To use this class, the following should be done:

  - Implement a user analysis class inheriting from this one, such as in user/james/process_mc_XX.py
    You should implement the following functions:
      - (optional) calculate_observable()
      - initialize_user_output_objects_R()
      - fill_observable_histograms()
      - fill_matched_jet_histograms()

  - You should include the following histograms:
      - Response matrix: hResponse_JetPt_[obs]_R[R]_[subobs]_[grooming setting]
      - Residual distribution: hResidual_JetPt_[obs]_R[R]_[subobs]_[grooming setting]

  - You also should modify observable-specific functions at the top of common_utils.py

Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import time

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
import fjtools

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_io_emb
from pyjetty.alice_analysis.process.base import process_base
from pyjetty.alice_analysis.process.base import thermal_generator
from pyjetty.mputils import CEventSubtractor

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class ProcessMCBase(process_base.ProcessBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):

    # Initialize base class
    super(ProcessMCBase, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)

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

    self.fast_simulation = config['fast_simulation']
    if 'jetscape' in config:
        self.jetscape = config['jetscape']
    else:
        self.jetscape = False
    if 'event_plane_angle' in config:
      self.event_plane_range = config['event_plane_angle']
    else:
      self.event_plane_range = None
    if 'matching_systematic' in config:
      self.matching_systematic = config['matching_systematic']
    else:
      self.matching_systematic = False
    self.dry_run = config['dry_run']
    self.skip_deltapt_RC_histograms = True
    self.fill_RM_histograms = True

    self.jet_matching_distance = config['jet_matching_distance']
    self.reject_tracks_fraction = config['reject_tracks_fraction']
    if 'mc_fraction_threshold' in config:
      self.mc_fraction_threshold = config['mc_fraction_threshold']

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

    # Whether or not to require jets to contain a track with some leading track pT
    self.min_leading_track_pT = config["min_leading_track_pT"] if \
      "min_leading_track_pT" in config else None

    # Create dictionaries to store grooming settings and observable settings for each observable
    # Each dictionary entry stores a list of subconfiguration parameters
    #   The observable list stores the observable setting, e.g. subjetR
    #   The grooming list stores a list of grooming settings {'sd': [zcut, beta]} or {'dg': [a]}
    self.observable_list = config['process_observables']
    self.obs_settings = {}
    self.obs_grooming_settings = {}
    self.obs_names = {}
    for observable in self.observable_list:

      obs_config_dict = config[observable]
      obs_config_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]

      obs_subconfig_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
      self.obs_settings[observable] = self.utils.obs_settings(observable, obs_config_dict, obs_subconfig_list)
      self.obs_grooming_settings[observable] = self.utils.grooming_settings(obs_config_dict)

      self.obs_names[observable] = obs_config_dict["common_settings"]["xtitle"]

    # Construct set of unique grooming settings
    self.grooming_settings = []
    lists_grooming = [self.obs_grooming_settings[obs] for obs in self.observable_list]
    for observable in lists_grooming:
      for setting in observable:
        if setting not in self.grooming_settings and setting != None:
          self.grooming_settings.append(setting)

    # Flag for creating delta-observable histograms in Pb-Pb case.
    # You can override this by setting the flag to True in your user class.
    # NOTE: requires implementation of self.calculate_observable() and creation
    #     of user histograms with this form:
    #     'hDeltaObs_%s_emb_R%s_%s%s' % (observable, jetR, obs_label, suffix)
    self.fill_delta_obs = False

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
                                  track_tree_name='tree_Particle', use_ev_id_ext=False,
                                  is_jetscape=self.jetscape, event_plane_range=self.event_plane_range)
    df_fjparticles_det = io_det.load_data(m=self.m, reject_tracks_fraction=self.reject_tracks_fraction)
    self.nEvents_det = len(df_fjparticles_det.index)
    self.nTracks_det = len(io_det.track_df.index)
    print('--- {} seconds ---'.format(time.time() - self.start_time))

    # If jetscape, store also the negative status particles (holes)
    if self.jetscape:
        io_det_holes = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
                                           track_tree_name='tree_Particle', use_ev_id_ext=False,
                                            is_jetscape=self.jetscape, holes=True,
                                            event_plane_range=self.event_plane_range)
        df_fjparticles_det_holes = io_det_holes.load_data(m=self.m, reject_tracks_fraction=self.reject_tracks_fraction)
        self.nEvents_det_holes = len(df_fjparticles_det_holes.index)
        self.nTracks_det_holes = len(io_det_holes.track_df.index)
        print('--- {} seconds ---'.format(time.time() - self.start_time))

    # ------------------------------------------------------------------------

    # Use IO helper class to convert truth-level ROOT TTree into
    # a SeriesGroupBy object of fastjet particles per event
    io_truth = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
                                    track_tree_name='tree_Particle_gen', use_ev_id_ext=False,
                                    is_jetscape=self.jetscape, event_plane_range=self.event_plane_range)
    df_fjparticles_truth = io_truth.load_data(m=self.m)
    self.nEvents_truth = len(df_fjparticles_truth.index)
    self.nTracks_truth = len(io_truth.track_df.index)
    print('--- {} seconds ---'.format(time.time() - self.start_time))

    # If jetscape, store also the negative status particles (holes)
    if self.jetscape:
        io_truth_holes = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
                                              track_tree_name='tree_Particle_gen', use_ev_id_ext=False,
                                              is_jetscape=self.jetscape, holes=True,
                                              event_plane_range=self.event_plane_range)
        df_fjparticles_truth_holes = io_truth_holes.load_data(m=self.m, reject_tracks_fraction=self.reject_tracks_fraction)
        self.nEvents_truth_holes = len(df_fjparticles_truth_holes.index)
        self.nTracks_truth_holes = len(io_truth_holes.track_df.index)
        print('--- {} seconds ---'.format(time.time() - self.start_time))

    # ------------------------------------------------------------------------

    # Now merge the two SeriesGroupBy to create a groupby df with [ev_id, run_number, fj_1, fj_2]
    # (Need a structure such that we can iterate event-by-event through both fj_1, fj_2 simultaneously)
    # In the case of jetscape, we merge also the hole collections fj_3, fj_4
    print('Merge det-level and truth-level into a single dataframe grouped by event...')
    if self.jetscape:
        self.df_fjparticles = pandas.concat([df_fjparticles_det, df_fjparticles_truth, df_fjparticles_det_holes, df_fjparticles_truth_holes], axis=1)
        self.df_fjparticles.columns = ['fj_particles_det', 'fj_particles_truth', 'fj_particles_det_holes', 'fj_particles_truth_holes']
    else:
        self.df_fjparticles = pandas.concat([df_fjparticles_det, df_fjparticles_truth], axis=1)
        self.df_fjparticles.columns = ['fj_particles_det', 'fj_particles_truth']
    print('--- {} seconds ---'.format(time.time() - self.start_time))

    # ------------------------------------------------------------------------

    # Set up the Pb-Pb embedding object
    if not self.is_pp and not self.thermal_model:
        self.process_io_emb = process_io_emb.ProcessIO_Emb(self.emb_file_list, track_tree_name='tree_Particle', m=self.m)

    # ------------------------------------------------------------------------

    # Initialize histograms
    if not self.dry_run:
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
    print('Find jets...')
    self.analyze_events()

    # Plot histograms
    print('Save histograms...')
    process_base.ProcessBase.save_output_objects(self)

    print('--- {} seconds ---'.format(time.time() - self.start_time))

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

    if not self.skip_deltapt_RC_histograms:
      name = 'hN_MeanPt'
      h = ROOT.TH2F(name, name, 200, 0, 5000, 200, 0., 2.)
      setattr(self, name, h)

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_output_objects_R(self, jetR):

      # Call user-specific initialization
      self.initialize_user_output_objects_R(jetR)

      # Base histograms
      if self.is_pp:

          name = 'hJES_R{}'.format(jetR)
          h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
          setattr(self, name, h)

          name = 'hDeltaR_All_R{}'.format(jetR)
          h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 2.)
          setattr(self, name, h)

      else:

          max_distance = self.max_distance if isinstance(self.max_distance, list) else \
                       self.max_distance[jetR]

          for R_max in max_distance:

            name = 'hJES_R{}_Rmax{}'.format(jetR, R_max)
            h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
            setattr(self, name, h)

            name = 'hDeltaPt_emb_R{}_Rmax{}'.format(jetR, R_max)
            h = ROOT.TH2F(name, name, 300, 0, 300, 400, -200., 200.)
            setattr(self, name, h)

            # Random cone histograms
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

      name = 'hZ_Truth_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 1.)
      setattr(self, name, h)

      name = 'hZ_Det_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 1.)
      setattr(self, name, h)

  #---------------------------------------------------------------
  # Main function to loop through and analyze events
  #---------------------------------------------------------------
  def analyze_events(self):

    # Fill track histograms
    if not self.dry_run:
      [self.fill_track_histograms(fj_particles_det) for fj_particles_det in self.df_fjparticles['fj_particles_det']]

    fj.ClusterSequence.print_banner()
    print()

    self.event_number = 0

    for jetR in self.jetR_list:
      if not self.dry_run:
        self.initialize_output_objects_R(jetR)

    # Then can use list comprehension to iterate over the groupby and do jet-finding
    # simultaneously for fj_1 and fj_2 per event, so that I can match jets -- and fill histograms
    if self.jetscape:
        result = [self.analyze_event(
          fj_particles_det, fj_particles_truth, fj_particles_det_holes, fj_particles_truth_holes) \
          for fj_particles_det, fj_particles_truth, fj_particles_det_holes, fj_particles_truth_holes \
          in zip(self.df_fjparticles['fj_particles_det'], self.df_fjparticles['fj_particles_truth'],
                 self.df_fjparticles['fj_particles_det_holes'], self.df_fjparticles['fj_particles_truth_holes'])]
    else:
        result = [self.analyze_event(fj_particles_det, fj_particles_truth) \
          for fj_particles_det, fj_particles_truth in zip(
          self.df_fjparticles['fj_particles_det'], self.df_fjparticles['fj_particles_truth'])]

    if self.debug_level > 0:
      for attr in dir(self):
        obj = getattr(self, attr)
        print('size of {}: {}'.format(attr, sys.getsizeof(obj)))

    print('Save thn...')
    process_base.ProcessBase.save_thn_th3_objects(self)

  #---------------------------------------------------------------
  # Fill track histograms.
  #---------------------------------------------------------------
  def fill_track_histograms(self, fj_particles_det):

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
  def analyze_event(self, fj_particles_det, fj_particles_truth, fj_particles_det_holes=None, fj_particles_truth_holes=None):

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
    if self.jetscape:
      if type(fj_particles_det_holes) != fj.vectorPJ or type(fj_particles_truth_holes) != fj.vectorPJ:
        print('fj_particles_holes type mismatch -- skipping event')
        return

    if len(fj_particles_truth) > 1:
      if np.abs(fj_particles_truth[0].pt() - fj_particles_truth[1].pt()) <  1e-10:
        print('WARNING: Duplicate particles may be present')
        print([p.user_index() for p in fj_particles_truth])
        print([p.pt() for p in fj_particles_truth])

    # If Pb-Pb, construct embedded event (do this once, for all jetR)
    if not self.is_pp:

        # If thermal model, generate a thermal event and add it to the det-level particle list
        if self.thermal_model:
          fj_particles_combined_beforeCS = self.thermal_generator.load_event()

          # Form the combined det-level event
          # The pp-det tracks are each stored with a unique user_index >= 0
          #   (same index in fj_particles_combined and fj_particles_det -- which will be used in prong-matching)
          # The thermal tracks are each stored with a unique user_index < 0
          for p in fj_particles_det:
            fj_particles_combined_beforeCS.push_back(p)

        # Main case: Get Pb-Pb event and embed it into the det-level particle list
        else:
          fj_particles_combined_beforeCS = self.process_io_emb.load_event()

          # Form the combined det-level event
          # The pp-det tracks are each stored with a unique user_index >= 0
          #   (same index in fj_particles_combined and fj_particles_det -- which will be used in prong-matching)
          # The Pb-Pb tracks are each stored with a unique user_index < 0
          for p in fj_particles_det:
            fj_particles_combined_beforeCS.push_back(p)

        # Perform constituent subtraction for each R_max
        fj_particles_combined = { R_max : cs.process_event(fj_particles_combined_beforeCS) for \
                                  R_max, cs in self.constituent_subtractor.items() }

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
        if self.min_leading_track_pT:
          print('*** Requiring minimum leading track pT:', self.min_leading_track_pT)
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

        max_distance = self.max_distance if isinstance(self.max_distance, list) else \
                       self.max_distance[jetR]

        # Keep track of whether to fill R_max-independent histograms
        self.fill_Rmax_indep_hists = True

        for R_max in max_distance:

          if self.debug_level > 1:
            print('')
            print('R_max: {}'.format(R_max))
            print('Total number of combined particles:', len([p.pt() for p in fj_particles_combined_beforeCS]))
            print('After constituent subtraction:', len([p.pt() for p in fj_particles_combined[R_max]]))

          # Perform constituent subtraction on det-level, if applicable
          self.fill_background_histograms(fj_particles_combined_beforeCS, fj_particles_combined[R_max], jetR, R_max)

          # Do jet finding (re-do each time, to make sure matching info gets reset)
          cs_det = fj.ClusterSequence(fj_particles_det, jet_def)
          jets_det_pp = fj.sorted_by_pt(cs_det.inclusive_jets())
          jets_det_pp_selected = jet_selector_det(jets_det_pp)

          cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)
          jets_truth = fj.sorted_by_pt(cs_truth.inclusive_jets())
          jets_truth_selected = jet_selector_det(jets_truth)
          jets_truth_selected_matched = jet_selector_truth_matched(jets_truth)

          cs_combined = fj.ClusterSequence(fj_particles_combined[R_max], jet_def)
          jets_combined = fj.sorted_by_pt(cs_combined.inclusive_jets())
          jets_combined_selected = jet_selector_det(jets_combined)

          self.analyze_jets(jets_combined_selected, jets_truth_selected, jets_truth_selected_matched, jetR,
                            jets_det_pp_selected = jets_det_pp_selected, R_max = R_max,
                            fj_particles_det_holes = fj_particles_det_holes,
                            fj_particles_truth_holes = fj_particles_truth_holes)

          # Don't fill R_max-independent histograms in future loops
          self.fill_Rmax_indep_hists = False

  #---------------------------------------------------------------
  # Analyze jets of a given event.
  #---------------------------------------------------------------
  def analyze_jets(self, jets_det_selected, jets_truth_selected, jets_truth_selected_matched, jetR,
                   jets_det_pp_selected = None, R_max = None,
                   fj_particles_det_holes = None, fj_particles_truth_holes = None):

    if self.debug_level > 1:
      print('Number of det-level jets: {}'.format(len(jets_det_selected)))

    # Fill det-level jet histograms (before matching)
    for jet_det in jets_det_selected:

      # Check additional acceptance criteria
      # skip event if not satisfied -- since first jet in event is highest pt
      if not self.utils.is_det_jet_accepted(jet_det, self.min_leading_track_pT):
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
      for jet_det in jets_det_selected:
        for jet_truth in jets_truth_selected_matched:
          self.set_matching_candidates(jet_det, jet_truth, jetR, 'hDeltaR_All_R{}'.format(jetR))
    else:
      for jet_det_pp in jets_det_pp_selected:
        # First fill the combined-to-pp matches, then the pp-to-pp matches
        for jet_det_combined in jets_det_selected:
          self.set_matching_candidates(jet_det_combined, jet_det_pp, jetR,
            'hDeltaR_combined_ppdet_R{{}}_Rmax{}'.format(R_max), fill_jet1_matches_only=True)
        for jet_truth in jets_truth_selected_matched:
          self.set_matching_candidates(
            jet_det_pp, jet_truth, jetR, 'hDeltaR_ppdet_pptrue_R{{}}_Rmax{}'.format(R_max))

    # Loop through jets and set accepted matches
    if self.is_pp:
      hname = 'hJetMatchingQA_R{}'.format(jetR)
      for jet_det in jets_det_selected:
        self.set_matches_pp(jet_det, hname)
    else:
      hname = 'hJetMatchingQA_R{}_Rmax{}'.format(jetR, R_max)
      for jet_det_combined in jets_det_selected:
        self.set_matches_AA(jet_det_combined, jetR, hname)

    # Loop through jets and fill response histograms if both det and truth jets are unique match
    for jet_det in jets_det_selected:
      self.fill_jet_matches(
        jet_det, jetR, R_max, fj_particles_det_holes, fj_particles_truth_holes)

  #---------------------------------------------------------------
  # Fill some background histograms
  #---------------------------------------------------------------
  def fill_background_histograms(self, fj_particles_combined_beforeCS, fj_particles_combined, jetR, R_max):

    # Fill rho
    rho = self.constituent_subtractor[R_max].bge_rho.rho()
    if self.fill_R_indep_hists and self.fill_Rmax_indep_hists:
      getattr(self, 'hRho').Fill(rho)

    # Fill random cone delta-pt before constituent subtraction
    if not self.skip_deltapt_RC_histograms:
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
    self.fill_unmatched_jet_histograms(jet, jetR, hname)

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
      self.fill_unmatched_jet_histograms(jet, jetR, hname)

  #---------------------------------------------------------------
  # This function is called once for each jet
  #---------------------------------------------------------------
  def fill_unmatched_jet_histograms(self, jet, jetR, hname):

    # Loop through each jet subconfiguration (i.e. subobservable / grooming setting)
    for observable in self.observable_list:
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
        self.fill_observable_histograms(observable, hname, jet, jet_groomed_lund,
            jetR, obs_setting, grooming_setting, obs_label, jet.pt())

  #---------------------------------------------------------------
  # Loop through jets and call user function to fill matched
  # histos if both det and truth jets are unique match.
  #---------------------------------------------------------------
  def fill_jet_matches(self, jet_det, jetR, R_max, fj_particles_det_holes, fj_particles_truth_holes):

    # Set suffix for filling histograms
    if R_max:
      suffix = '_Rmax' + str(R_max)
    else:
      suffix = ''

    # Get matched truth jet
    if jet_det.has_user_info():
      jet_truth = jet_det.python_info().match

      if jet_truth:

        jet_pt_det_ungroomed = jet_det.pt()
        jet_pt_truth_ungroomed = jet_truth.pt()
        JES = (jet_pt_det_ungroomed - jet_pt_truth_ungroomed) / jet_pt_truth_ungroomed
        getattr(self, 'hJES_R' + str(jetR) + suffix).Fill(jet_pt_truth_ungroomed, JES)

        # If Pb-Pb case, we need to keep jet_det, jet_truth, jet_pp_det
        jet_pp_det = None
        if not self.is_pp:

          # Get pp-det jet
          jet_pp_det = jet_truth.python_info().match

          # Fill delta-pt histogram
          if jet_pp_det:
            getattr(self, 'hDeltaPt_emb_R' + str(jetR) + suffix).Fill(
                jet_pt_truth_ungroomed, jet_pt_det_ungroomed - jet_pp_det.pt())

        # Loop through each jet subconfiguration (i.e. subobservable / grooming setting)
        for observable in self.observable_list:
          for i in range(len(self.obs_settings[observable])):

            obs_setting = self.obs_settings[observable][i]
            grooming_setting = self.obs_grooming_settings[observable][i]
            obs_label = self.utils.obs_label(obs_setting, grooming_setting)

            if self.debug_level > 3:
              print('obs_label:', obs_label)

            # Groom jets, if applicable
            if grooming_setting:

              # Groom det jet
              gshop_det = fjcontrib.GroomerShop(jet_det, jetR, self.reclustering_algorithm)
              jet_det_groomed_lund = self.utils.groom(gshop_det, grooming_setting, jetR)
              if not jet_det_groomed_lund:
                continue

              # Groom truth jet
              gshop_truth = fjcontrib.GroomerShop(jet_truth, jetR, self.reclustering_algorithm)
              jet_truth_groomed_lund = self.utils.groom(gshop_truth, grooming_setting, jetR)
              if not jet_truth_groomed_lund:
                continue

            else:

              jet_det_groomed_lund = None
              jet_truth_groomed_lund = None

            # If jetscape, pass the list of holes within R of the jet to the user
            holes_in_det_jet = None
            holes_in_truth_jet = None
            if self.jetscape:
              holes_in_det_jet = [hadron for hadron in fj_particles_det_holes if \
                                  jet_det.delta_R(hadron) < jetR]
              holes_in_truth_jet = [hadron for hadron in fj_particles_truth_holes if \
                                    jet_truth.delta_R(hadron) < jetR]

              # Get the corrected jet pt by subtracting the negative recoils within R
              for hadron in holes_in_det_jet:
                  jet_pt_det_ungroomed -= hadron.pt()

              for hadron in holes_in_truth_jet:
                  jet_pt_truth_ungroomed -= hadron.pt()

            # Fill delta-observable histograms
            if self.fill_delta_obs and not self.is_pp:

              # Get pp-det jet
              jet_pp_det = jet_truth.python_info().match
              if jet_pp_det:

                # Groom jet, if applicable
                jet_pp_det_groomed_lund = None
                if grooming_setting:
                  gshop_pp_det = fjcontrib.GroomerShop(
                      jet_pp_det, jetR, self.reclustering_algorithm)
                  jet_pp_det_groomed_lund = self.utils.groom(
                      gshop_pp_det, grooming_setting, jetR)
                  if not jet_pp_det_groomed_lund:
                    continue

                obs_det = self.calculate_observable(observable, jet_det,
                    jet_det_groomed_lund, jetR, obs_setting, grooming_setting,
                    obs_label, jet_pt_det_ungroomed)

                obs_pp_det = self.calculate_observable(observable, jet_pp_det,
                    jet_pp_det_groomed_lund, jetR, obs_setting, grooming_setting,
                    obs_label, jet_pp_det.pt())

                getattr(self, 'hDeltaObs_%s_emb_R%s_%s%s' % (observable, jetR,
                    obs_label, suffix)).Fill(jet_pt_truth_ungroomed, obs_det - obs_pp_det)

            # Call user function to fill histos
            self.fill_matched_jet_histograms(observable, jet_det, jet_det_groomed_lund,
                jet_truth, jet_truth_groomed_lund, jet_pp_det, jetR, obs_setting,
                grooming_setting, obs_label, jet_pt_det_ungroomed, jet_pt_truth_ungroomed,
                R_max, suffix, holes_in_det_jet=holes_in_det_jet,
                holes_in_truth_jet=holes_in_truth_jet)

  #---------------------------------------------------------------
  # Fill response histograms -- common utility function
  #---------------------------------------------------------------
  def fill_response(self, observable, jetR, jet_pt_det_ungroomed, jet_pt_truth_ungroomed,
                    obs_det, obs_truth, obs_label, R_max, prong_match = False):

    if self.fill_RM_histograms:
      x = ([jet_pt_det_ungroomed, jet_pt_truth_ungroomed, obs_det, obs_truth])
      x_array = array('d', x)
      name = 'hResponse_JetPt_{}_R{}_{}'.format(observable, jetR, obs_label) if \
        len(obs_label) else 'hResponse_JetPt_{}_R{}'.format(observable, jetR)
      if not self.is_pp:
        name += '_Rmax{}'.format(R_max)
      getattr(self, name).Fill(x_array)

    if obs_truth > 1e-5:
      obs_resolution = (obs_det - obs_truth) / obs_truth
      name = 'hResidual_JetPt_{}_R{}_{}'.format(observable, jetR, obs_label) if \
        len(obs_label) else 'hResidual_JetPt_{}_R{}'.format(observable, jetR)
      if not self.is_pp:
        name += '_Rmax{}'.format(R_max)
      getattr(self, name).Fill(jet_pt_truth_ungroomed, obs_truth, obs_resolution)

    # Fill prong-matched response
    if not self.is_pp:
      main_R_max = self.main_R_max if isinstance(self.main_R_max, float) else self.main_R_max[jetR]
      if prong_match and R_max == main_R_max:

        name = 'hResponse_JetPt_{}_R{}_{}_Rmax{}_matched'.format(
          observable, jetR, obs_label, R_max) if len(obs_label) else \
          'hResponse_JetPt_{}_R{}_Rmax{}_matched'.format(observable, jetR, R_max)
        getattr(self, name).Fill(x_array)

        if obs_truth > 1e-5:
          name = 'hResidual_JetPt_{}_R{}_{}_Rmax{}_matched'.format(
            observable, jetR, obs_label, R_max) if len(obs_label) else \
            'hResidual_JetPt_{}_R{}_Rmax{}_matched'.format(observable, jetR, R_max)
          getattr(self, name).Fill(jet_pt_truth_ungroomed, obs_truth, obs_resolution)


  #---------------------------------------------------------------
  # Calculate the observable given a jet
  #---------------------------------------------------------------
  def calculate_observable(self, observable, jet, jet_groomed_lund,
      jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed):

    raise NotImplementedError('You must implement caclulate_observable()!')

  #---------------------------------------------------------------
  # This function is called once for each jetR
  # You must implement this
  #---------------------------------------------------------------
  def initialize_user_output_objects_R(self, jetR):

    raise NotImplementedError('You must implement initialize_user_output_objects_R()!')

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  # You must implement this
  #---------------------------------------------------------------
  def fill_observable_histograms(self, observable, hname, jet, jet_groomed_lund,
      jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed):

    raise NotImplementedError('You must implement fill_observable_histograms()!')

  #---------------------------------------------------------------
  # This function is called once for each matched jet subconfiguration
  # You must implement this
  #---------------------------------------------------------------
  def fill_matched_jet_histograms(self, observable, jet_det, jet_det_groomed_lund,
      jet_truth, jet_truth_groomed_lund, jet_pp_det, jetR, obs_setting,
      grooming_setting, obs_label, jet_pt_det_ungroomed, jet_pt_truth_ungroomed,
      R_max, suffix, **kwargs):

    raise NotImplementedError('You must implement fill_matched_jet_histograms()!')

  #---------------------------------------------------------------
  # Optional testing function for embedding
  # Useful for background subtraction studies
  #---------------------------------------------------------------
  def fill_jet_embedding_test_histograms(self):

    # Override default functionality in user code, if desired
    return
