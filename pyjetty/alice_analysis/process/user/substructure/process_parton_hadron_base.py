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
import math  # for ceil
import time  # for timing code
import sys
import gc    # garbage collection for cleaning memory

# Data analysis and plotting
import pandas
import numpy as np
from array import *
import ROOT
import yaml
import random

# Load pyjetty ROOT utils
ROOT.gSystem.Load('libpyjetty_rutil')

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io_parton_hadron
from pyjetty.alice_analysis.process.base import process_base
from pyjetty.alice_analysis.process.base import common_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

# Helper function
# Turns 2D list into 1D list by concat'ing all sublists
def li_concat(li):
  return [item for sublist in li for item in sublist]

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

    # C++ histogram rebinning functions
    # Don't actually need any of these, but have to init to get other
    #     RUtil functions for some reason...
    self.histutils = ROOT.RUtil.HistUtils()

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

    self.make_th3s = False if 'make_th3s' not in config else config['make_th3s']

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

    # These values set the step size and control memory usage
    self.track_df_step_size = 20  #1e5
    self.max_ev_storage = 5e4

    fj.ClusterSequence.print_banner()

  #---------------------------------------------------------------
  # Initialize empty storage dictionaries to store events for use in next iteration
  #---------------------------------------------------------------
  def init_storage(self, level, MPI):
    setattr(self, "df_fjparticles_storage_%s_MPI%s" % (level, MPI),
            { "run_number": [], "ev_id": [], "fj_particle": [] })
    setattr(self, "run_numbers_storage_%s_MPI%s" % (level, MPI), [])
    setattr(self, "unique_ev_ids_per_run_storage_%s_MPI%s" % (level, MPI), [])

  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def process_mc(self):

    self.start_time = time.time()

    # ------------------------------------------------------------------------

    # Initialize histograms
    self.initialize_output_objects()

    if self.debug_level > 0:
      print(self)

    print('--- {} seconds ---'.format(time.time() - self.start_time))

    # ------------------------------------------------------------------------

    for MPI in ["off", "on"]:

      # Initialize necessary objects for looping
      self.init_tree_readers(MPI)
      getattr(self, "hNevents_MPI%s" % MPI).Fill(1, getattr(self, "nEvents_MPI%s" % MPI))
      for level in ["p", "h", "ch"]:
        self.init_storage(level, MPI)

      # Initialize iteration number & max iteration number
      setattr(self, "df_iter_p_MPI"+MPI, 0)
      max_iter_p = math.ceil(getattr(self, "tree_len_p_MPI"+MPI) / self.track_df_step_size)
      setattr(self, "df_iter_h_MPI"+MPI, 0)
      max_iter_h = math.ceil(getattr(self, "tree_len_h_MPI"+MPI) / self.track_df_step_size)

      # Iterate through the track TTree in steps
      # (This strategy saves memory as compared to loading the entire TTree at once)
      while getattr(self, "df_iter_p_MPI"+MPI) < max_iter_p or \
            getattr(self, "df_iter_h_MPI"+MPI) < max_iter_h or \
            len(getattr(self, "df_fjparticles_storage_p_MPI"+MPI)["run_number"]) > 0 or \
            len(getattr(self, "df_fjparticles_storage_h_MPI"+MPI)["run_number"]) > 0:

        if self.debug_level > 1:
          print(getattr(self, "df_iter_p_MPI"+MPI), getattr(self, "df_iter_h_MPI"+MPI),
                len(getattr(self, "df_fjparticles_storage_p_MPI"+MPI)["run_number"]),
                len(getattr(self, "df_fjparticles_storage_h_MPI"+MPI)["run_number"]))

        print("Load TTrees with MPI %s, iteration %i/%i (p) | %i/%i (h)..." % \
              (MPI, getattr(self, "df_iter_p_MPI"+MPI), max_iter_p,
               getattr(self, "df_iter_h_MPI"+MPI), max_iter_h))
        self.load_trees_to_dict(MPI)
        #print('--- {} seconds ---'.format(time.time() - self.start_time))

        # ------------------------------------------------------------------------

        # Find jets and fill histograms
        print('Find jets with MPI %s...' % MPI)
        self.analyze_events(MPI)

        print('\n--- {} seconds ---'.format(time.time() - self.start_time))

    # ------------------------------------------------------------------------

    print("Scale histograms by appropriate weighting...")
    self.scale_objects()

    print('Save thn...')
    process_base.ProcessBase.save_thn_th3_objects(self)

    # Plot histograms
    print('Save histograms...')
    process_base.ProcessBase.save_output_objects(self)

    print('--- {} seconds ---'.format(time.time() - self.start_time))

  #---------------------------------------------------------------
  # Intialize track tree readers and save sizes for proper looping
  #---------------------------------------------------------------
  def init_tree_readers(self, MPI):
    setattr(self, "io_p_MPI"+MPI,
            process_io_parton_hadron.ProcessIO(input_file=self.input_file, level='p', MPI=MPI))
    setattr(self, "tree_len_p_MPI"+MPI, getattr(self, "io_p_MPI"+MPI).tree_length)

    setattr(self, "io_h_MPI"+MPI,
            process_io_parton_hadron.ProcessIO(input_file=self.input_file, level='h', MPI=MPI))
    setattr(self, "tree_len_h_MPI"+MPI, getattr(self, "io_h_MPI"+MPI).tree_length)

    # Save number of events and cross section information for proper scaling
    setattr(self, "nEvents_MPI"+MPI, getattr(self, "io_p_MPI"+MPI).Nev)
    setattr(self, "xsec_MPI"+MPI, getattr(self, "io_p_MPI"+MPI).xsec)

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_output_objects(self):

    for MPI in ["on", "off"]:

      setattr(self, "hNevents_MPI%s" % MPI, ROOT.TH1F(
        'hNevents_MPI%s' % MPI, 'hNevents_MPI%s' % MPI, 2, -0.5, 1.5))

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
  # Open ROOT TTrees, load to dictionaries, and merge
  #---------------------------------------------------------------
  def load_trees_to_dict(self, MPI):

    ###### level == "p" or "h" ######
    copy_h = True  # Save bool for checking ch case separately
    for level in ["p", "h"]:

      # Copy saved events from storage for new processing step
      self.copy_df_from_storage(level, MPI)

      # Check if there are already "enough" ev_id's loaded, to minimize memory usage
      # or if we have already read all the info from the TTree
      it = getattr(self, "df_iter_%s_MPI%s" % (level, MPI))
      tree_len = getattr(self, "tree_len_%s_MPI%s" % (level, MPI))
      if len(getattr(self, "df_fjparticles_%s_MPI%s" % (level, MPI))["run_number"]) > \
         self.max_ev_storage or it*self.track_df_step_size >= tree_len:
        if level == "h":
          copy_h = False
        continue

      # Load ROOT TTrees into Python objects
      df_to_append = getattr(self, "io_%s_MPI%s" % (level, MPI)).load_data(
        start=int(it*self.track_df_step_size), stop=int(min((it+1)*self.track_df_step_size, tree_len)) )
      #print("df len:", len(df_to_append["run_number"]))
      self.append_df_fjparticles(level, MPI, df_to_append)
      setattr(self, "df_iter_%s_MPI%s" % (level, MPI), it+1)

      # Save run_number / ev_id information for pairing & trimming
      self.append_run_evid(getattr(self, "io_%s_MPI%s" % (level, MPI)), level, MPI)

    ###### level == "ch" ######
    # Do the ch case separately since it is slightly different
    self.copy_df_from_storage("ch", MPI)

    if copy_h:
      # Do the charge cut on the hadron-level tree (don't need to reload twice)
      df_to_append = getattr(self, "io_h_MPI"+str(MPI)).group_fjparticles(ch_cut=True)
      self.append_df_fjparticles("ch", MPI, df_to_append)
      self.append_run_evid(getattr(self, "io_h_MPI"+str(MPI)), "ch", MPI)

    #print('--- {} seconds ---'.format(time.time() - self.start_time))

    # ------------------------------------------------------------------------

    # Move any extra ev's at the end to storage
    self.move_extra_ev_to_storage(MPI)

    # Now merge the dfs to create a single df with [ev_id, run_number, fj_p, fj_h, fj_ch]
    # (Need a structure such that we can iterate event-by-event through
    #     fj_p, fj_h, and fj_ch simultaneously)
    print('Merge parton-, hadron-, and charged-level into a single object grouped by event...')

    setattr(self, "df_fjparticles_MPI"+MPI, self.pair_dictionary(MPI))

    # Clean up memory from unpaired trees
    for level in ["p", "h", "ch"]:
      self.del_unpaired_dfs(level, MPI)
    gc.collect()

  #---------------------------------------------------------------
  # Move storage dataframes to "current" dataframes and reset (empty) storage
  #---------------------------------------------------------------
  def copy_df_from_storage(self, level, MPI):

    setattr(self, "df_fjparticles_%s_MPI%s" % (level, MPI), getattr(
      self, "df_fjparticles_storage_%s_MPI%s" % (level, MPI)).copy())
    setattr(self, "run_numbers_%s_MPI%s" % (level, MPI), getattr(
      self, "run_numbers_storage_%s_MPI%s" % (level, MPI)).copy())
    setattr(self, "unique_ev_ids_per_run_%s_MPI%s" % (level, MPI), getattr(
      self, "unique_ev_ids_per_run_storage_%s_MPI%s" % (level, MPI)).copy())

    # Remove any existing events from storage and create new pointers
    self.init_storage(level, MPI)

  #---------------------------------------------------------------
  def append_df_fjparticles(self, level, MPI, df_to_append):
    name = "df_fjparticles_%s_MPI%s" % (level, MPI)
    # Merge events that were not fully loaded in last iteration
    if len(df_to_append["ev_id"]) and len(getattr(self, name)["ev_id"]) and \
       df_to_append["ev_id"][0] == getattr(self, name)["ev_id"][-1]:
      for key, val in df_to_append.items():
        if "fj" in key:
          getattr(self, name)[key][-1] += df_to_append[key][0]
        df_to_append[key].pop(0)

    print("old ev tree")
    self.print_df(getattr(self, name))
    print("appending")
    self.print_df(df_to_append)

    # Merge dfs
    for key, val in df_to_append.items():
      getattr(self, name)[key] += val


  #---------------------------------------------------------------
  # Function to look inside df in convenient way
  #---------------------------------------------------------------
  def print_df(self, df):
    if len(df["run_number"]) < 1:
      print("empty df")
    print("run_number", "ev_id", "pt\t\t", "eta\t\t", "phi", sep='\t')
    for i, run_number in enumerate(df["run_number"]):
      ev_id = df["ev_id"][i]
      for j, fj_particle in enumerate(df["fj_particle"][i]):
        print("%s\t" % run_number, ev_id, fj_particle.pt(), fj_particle.eta(),
              fj_particle.phi(), sep='\t')


  #---------------------------------------------------------------
  # Save the run_numbers and unique ev_ids as attributes for given io, level
  #---------------------------------------------------------------
  def append_run_evid(self, io, level, MPI):

    nr = "run_numbers_%s_MPI%s" % (level, MPI)
    ne = "unique_ev_ids_per_run_%s_MPI%s" % (level, MPI)

    if self.debug_level > 1:
      print("nr before:", getattr(self, nr))

    # Merge ev info with existing rather than creating additional list per run
    if len(getattr(self, nr)) == 0:
      setattr(self, nr, getattr(self, nr) + list(io.run_numbers[0]))

      if len(getattr(self, ne)):
        if getattr(self, ne)[-1][0][-1] == io.unique_ev_ids_per_run[0][0][0]:
          # Remove duplicate event number which was split across iterations
          getattr(self, ne)[-1][0] = getattr(self, ne)[-1][0][:-1]
        setattr(self, ne, list(getattr(self, ne)[:-1]) + [[np.concatenate((
          getattr(self, ne)[-1][0], io.unique_ev_ids_per_run[0][0]))]] + \
                list(io.unique_ev_ids_per_run)[1:])
      else:
        setattr(self, ne, list(io.unique_ev_ids_per_run))

    elif getattr(self, nr)[-1] == io.run_numbers[0][0]:
      # Ignore repeated entry
      setattr(self, nr, getattr(self, nr) + list(io.run_numbers[0][:-1]))

      if getattr(self, ne)[-1][0][-1] == io.unique_ev_ids_per_run[0][0][0]:
        # Remove duplicate event number which was split across iterations
        getattr(self, ne)[-1][0] = getattr(self, ne)[-1][0][1:]
      setattr(self, ne, list(getattr(self, ne)[:-1]) + [[np.concatenate((
        getattr(self, ne)[-1][0], io.unique_ev_ids_per_run[0][0]))]] + \
              list(io.unique_ev_ids_per_run)[1:])

    else:  # different runs
      setattr(self, nr, getattr(self, nr) + list(io.run_numbers[0]))
      setattr(self, ne, getattr(self, ne) + list(io.unique_ev_ids_per_run))

    if self.debug_level > 1:
      print("nr after:", getattr(self, nr))
    #setattr(self, n, getattr(self, n) + list(io.unique_ev_ids_per_run))

  #---------------------------------------------------------------
  # Delete all unpaired track / run / event DFs from memory
  #---------------------------------------------------------------
  def del_unpaired_dfs(self, level, MPI):
    delattr(self, "df_fjparticles_%s_MPI%s" % (level, MPI))
    delattr(self, "run_numbers_%s_MPI%s" % (level, MPI))
    delattr(self, "unique_ev_ids_per_run_%s_MPI%s" % (level, MPI))

  #---------------------------------------------------------------
  # Find the ev's at the end which have not yet been parsed in one or
  # more df's and move them to storage for future parsing
  #---------------------------------------------------------------
  def move_extra_ev_to_storage(self, MPI):

    # Check to make sure that there are events at each level
    # This could happen if there is only 1 ev left but it is missing at some level
    if not len(getattr(self, "df_fjparticles_p_MPI" + MPI)["run_number"]) or \
       not len(getattr(self, "df_fjparticles_h_MPI" + MPI)["run_number"]) or \
       not len(getattr(self, "df_fjparticles_ch_MPI" + MPI)["run_number"]):
      return

    # Check to see if there are any extra run_numbers
    self.check_and_move_extra_runs(MPI)

    # Now check to see if there are extra events
    self.check_and_move_extra_evs(MPI)

  #---------------------------------------------------------------
  # Check track df for extra runs and move to storage
  #---------------------------------------------------------------
  def check_and_move_extra_runs(self, MPI):

    # Load required info for splicing
    last_run_nums = [getattr(self, "df_fjparticles_p_MPI" + MPI)["run_number"][-1],
                     getattr(self, "df_fjparticles_h_MPI" + MPI)["run_number"][-1],
                     getattr(self, "df_fjparticles_ch_MPI" + MPI)["run_number"][-1]]

    if len(np.unique(last_run_nums)) != 1:

      last_run = min(last_run_nums)  # assuming run_number is sorted in TTree
      for i, level in zip([0, 1, 2], ["p", "h", "ch"]):
        if last_run_nums[i] != last_run:
          extra_i, extra_run_i = self.get_extra_index_run(level, MPI, last_run)
          self.move_df_to_storage_index(level, MPI, extra_i)
          self.move_run_to_storage_index(level, MPI, extra_run_i)

  #---------------------------------------------------------------
  # Return first index in track df where extra runs start
  #---------------------------------------------------------------
  def get_extra_index_run(self, level, MPI, last_run):

    run_numbers = getattr(self, "run_numbers_%s_MPI%s" % (level, MPI))
    unique_ev_ids_per_run = getattr(self, "unique_ev_ids_per_run_%s_MPI%s" % (level, MPI))

    extra_runs = []
    for i in range(len(run_numbers)-1, -1, -1):
      if run_numbers[i] != last_run:
        extra_runs += [run_numbers[i]]
      else:
        break

    extra_run_i = len(run_numbers) - len(extra_runs)
    return sum([len(unique_ev_ids_per_run[i][0]) for i in range(extra_run_i)]), extra_run_i

  #---------------------------------------------------------------
  # Check for extra events in track df and move them to storage
  #---------------------------------------------------------------
  def check_and_move_extra_evs(self, MPI):

    # Load required info for splicing
    last_ev_ids = [getattr(self, "df_fjparticles_p_MPI" + MPI)["ev_id"][-1],
                   getattr(self, "df_fjparticles_h_MPI" + MPI)["ev_id"][-1],
                   getattr(self, "df_fjparticles_ch_MPI" + MPI)["ev_id"][-1]]

    last_ev = min(last_ev_ids)  # assuming ev_id is sorted in TTree
    # Unless there is only 1 ev, and only if it is the last event in the file,
    #     move last ev to storage to make sure we have all its tracks
    if getattr(self, "df_fjparticles_p_MPI" + MPI)["ev_id"][0] == max(last_ev_ids):
      # Look for a case where we have not exhausted the datafile
      for level in ["p", "h"]:
        if getattr(self, "tree_len_%s_MPI%s" % (level, MPI)) <= \
           getattr(self, "df_iter_%s_MPI%s" % (level, MPI)) * self.track_df_step_size:
          last_ev -= 1
          break
    else:
      last_ev -= 1
    if self.debug_level > 1:
      print("last_ev:", last_ev)

    if last_ev != min(last_ev_ids) or len(np.unique(last_ev_ids)) != 1:

      for i, level in zip([0, 1, 2], ["p", "h", "ch"]):
        if last_ev_ids[i] != last_ev:
          extra_i, extra_i_ev = self.get_extra_index_ev(level, MPI, last_ev)
          self.move_df_to_storage_index(level, MPI, extra_i)
          self.move_ev_to_storage_index(level, MPI, extra_i_ev)
          # Update run storage to contain run number of these extra events
          self.add_last_run_to_storage(level, MPI)

  #---------------------------------------------------------------
  # Return first index in track df where extra runs start
  #---------------------------------------------------------------
  def get_extra_index_ev(self, level, MPI, last_ev):

    run_numbers = getattr(self, "run_numbers_%s_MPI%s" % (level, MPI))
    unique_ev_ids_per_run = getattr(self, "unique_ev_ids_per_run_%s_MPI%s" % (level, MPI))

    extra_evs = []
    # We have already moved the extra runs, so only have to consider the last one
    run_i = len(run_numbers)-1
    if self.debug_level > 1:
      print("unique:", unique_ev_ids_per_run[run_i][0])
    for i in range(len(unique_ev_ids_per_run[run_i][0])-1, -1, -1):
      if unique_ev_ids_per_run[run_i][0][i] > last_ev:
        extra_evs += [unique_ev_ids_per_run[run_i][0][i]]
      else:
        break

    extra_ev_i = len(unique_ev_ids_per_run[run_i][0]) - len(extra_evs)
    return sum([len(unique_ev_ids_per_run[i][0]) for i in range(run_i)]) + extra_ev_i, extra_ev_i

  #---------------------------------------------------------------
  # Move entries in track df to storage starting at move_index
  #---------------------------------------------------------------
  def move_df_to_storage_index(self, level, MPI, move_index):

    df_name_storage = "df_fjparticles_storage_%s_MPI%s" % (level, MPI)
    df_name = "df_fjparticles_%s_MPI%s" % (level, MPI)
    for key, val in getattr(self, df_name).items():
      # Move to _beginning_ of storage array (we move extra runs first)
      getattr(self, df_name_storage)[key] = val[move_index:] + getattr(self, df_name_storage)[key]
      getattr(self, df_name)[key] = val[:move_index]

  #---------------------------------------------------------------
  # Move entries in runs/ev_ids to storage starting at move_index
  #---------------------------------------------------------------
  def move_run_to_storage_index(self, level, MPI, move_index):

    names = ["run_numbers_%s_MPI%s" % (level, MPI),
             "unique_ev_ids_per_run_%s_MPI%s" % (level, MPI)]
    storage_names = ["run_numbers_storage_%s_MPI%s" % (level, MPI),
                     "unique_ev_ids_per_run_storage_%s_MPI%s" % (level, MPI)]

    for n, sn in zip(names, storage_names):
      setattr(self, sn, getattr(self, sn) + getattr(self, n)[move_index:])
      setattr(self, n, getattr(self, n)[:move_index])

  #---------------------------------------------------------------
  # Move ev_ids in last run to storage starting at move_index
  #---------------------------------------------------------------
  def move_ev_to_storage_index(self, level, MPI, move_index):

    n = "unique_ev_ids_per_run_%s_MPI%s" % (level, MPI)
    sn = "unique_ev_ids_per_run_storage_%s_MPI%s" % (level, MPI)

    setattr(self, sn, getattr(self, sn) + [[getattr(self, n)[-1][0][move_index:]]])
    setattr(self, n, getattr(self, n)[:-1] + [[getattr(self, n)[-1][0][:move_index]]])
    if self.debug_level > 1:
      print("storage ev id %s:" % level, getattr(self, sn)[-1][0])
      print("ana ev id %s:" % level, getattr(self, n)[-1][0])

  #---------------------------------------------------------------
  # Copy last run number in run_numbers to end of storage list
  #---------------------------------------------------------------
  def add_last_run_to_storage(self, level, MPI):

    n = "run_numbers_%s_MPI%s" % (level, MPI)
    sn = "run_numbers_storage_%s_MPI%s" % (level, MPI)

    setattr(self, sn, list(getattr(self, sn)) + [getattr(self, n)[-1]])
    if self.debug_level > 1:
      print("storage run id:", getattr(self, sn), type(getattr(self, sn)))

  #---------------------------------------------------------------
  # Pair NumPy dictionaries at p, h, and ch levels
  #---------------------------------------------------------------
  def pair_dictionary(self, MPI):

    run_numbers_p = getattr(self, "run_numbers_p_MPI" + MPI)
    run_numbers_h = getattr(self, "run_numbers_h_MPI" + MPI)
    run_numbers_ch = getattr(self, "run_numbers_ch_MPI" + MPI)

    unique_ev_ids_per_run_p = getattr(self, "unique_ev_ids_per_run_p_MPI" + MPI)
    unique_ev_ids_per_run_h = getattr(self, "unique_ev_ids_per_run_h_MPI" + MPI)
    unique_ev_ids_per_run_ch = getattr(self, "unique_ev_ids_per_run_ch_MPI" + MPI)

    df_fjparticles_p = getattr(self, "df_fjparticles_p_MPI" + MPI)
    df_fjparticles_h = getattr(self, "df_fjparticles_h_MPI" + MPI)
    df_fjparticles_ch = getattr(self, "df_fjparticles_ch_MPI" + MPI)

    # Need to figure out which ev_id to save.
    # For some reason, it can be the case that the ev_id does not exist
    #    at one level (e.g. parton), while it does at another (e.g. hadron).
    # Deal with this by simply removing these bad run_numbers / ev_ids.
    # First, check if there are any bad run_numbers
    bad_index_p = None
    bad_runs_p = np.setdiff1d(run_numbers_p, run_numbers_ch)
    if len(bad_runs_p):
      # This is what we want, but can get it more trivially using available info
      #bad_index_p = np.isin(df_fjparticles_p_MPIoff["run_number"], bad_runs_p)
      print("bad_runs_p:", bad_runs_p)
      bad_index_p = li_concat([
        [run_numbers_p[i] in bad_runs_p]*len(unique_ev_ids_per_run_p[i][0]) \
        for i in range(len(run_numbers_p))])
    else:
      bad_index_p = li_concat([[False] * len(unique_ev_ids_per_run_p[i][0]) \
                               for i in range(len(run_numbers_p))])

    bad_index_h = None
    bad_runs_h = np.concatenate((np.setdiff1d(run_numbers_h, run_numbers_p),
                                 np.setdiff1d(run_numbers_h, run_numbers_ch)))
    if len(bad_runs_h):
      print("bad_runs_h:", bad_runs_h)
      bad_index_h = li_concat([
        [run_numbers_h[i] in bad_runs_h]*len(unique_ev_ids_per_run_h[i][0]) \
        for i in range(len(run_numbers_h))])
    else:
      bad_index_h = li_concat([[False] * len(unique_ev_ids_per_run_h[i][0]) \
                               for i in range(len(run_numbers_h))])

    bad_index_ch = None
    bad_runs_ch = bad_runs_p  #np.setdiff1d(run_numbers_ch, run_numbers_p)
    if len(bad_runs_ch):
      print("bad_runs_ch:", bad_runs_ch)
      bad_index_ch = li_concat([
        [run_numbers_ch[i] in bad_runs_ch]*len(unique_ev_ids_per_run_ch[i][0]) \
        for i in range(len(run_numbers_ch))])
    else:
      bad_index_ch = li_concat([[False] * len(unique_ev_ids_per_run_ch[i][0]) \
                                for i in range(len(run_numbers_ch))])

    # Now check for bad ev_ids for each run.
    overall_index_p = 0
    run_index_p = 0
    overall_index_h = 0
    run_index_h = 0
    overall_index_ch = 0
    run_index_ch = 0

    while run_index_p < len(run_numbers_p):
      ev_ids_p = unique_ev_ids_per_run_p[run_index_p][0]
      if run_numbers_p[run_index_p] in bad_runs_p:
        overall_index_p += len(ev_ids_p)
        run_index_p += 1
        continue
      ev_ids_h = unique_ev_ids_per_run_h[run_index_h][0]
      if run_numbers_h[run_index_h] in bad_runs_h:
        overall_index_h += len(ev_ids_h)
        run_index_h += 1
        continue
      ev_ids_ch = unique_ev_ids_per_run_ch[run_index_ch][0]
      if run_numbers_ch[run_index_ch] in bad_runs_ch:
        overall_index_ch += len(ev_ids_ch)
        run_index_ch += 1
        continue

      # Use the fact that ev_ids are sorted for efficiency improvement
      good_ev_id_index_1 = ROOT.RUtil.sorted_match(
        ev_ids_p, len(ev_ids_p), ev_ids_ch, len(ev_ids_ch))
      bad_ev_id_index = np.array([not good_ev_id_index_1[i] for i in range(len(ev_ids_p))])
      bad_index_p[overall_index_p:overall_index_p+len(bad_ev_id_index)] = bad_index_p[
        overall_index_p:overall_index_p+len(bad_ev_id_index)] | bad_ev_id_index

      good_ev_id_index_1 = ROOT.RUtil.sorted_match(
        ev_ids_h, len(ev_ids_h), ev_ids_p, len(ev_ids_p))
      good_ev_id_index_2 = ROOT.RUtil.sorted_match(
        ev_ids_h, len(ev_ids_h), ev_ids_ch, len(ev_ids_ch))
      bad_ev_id_index = np.array(
        [not (good_ev_id_index_1[i] and good_ev_id_index_2[i]) for i in range(len(ev_ids_h))])
      bad_index_h[overall_index_h:overall_index_h+len(bad_ev_id_index)] = bad_index_h[
        overall_index_h:overall_index_h+len(bad_ev_id_index)] | bad_ev_id_index

      good_ev_id_index_1 = ROOT.RUtil.sorted_match(
        ev_ids_ch, len(ev_ids_ch), ev_ids_p, len(ev_ids_p))
      bad_ev_id_index = np.array([not good_ev_id_index_1[i] for i in range(len(ev_ids_ch))])
      bad_index_ch[overall_index_ch:overall_index_ch+len(bad_ev_id_index)] = bad_index_ch[
        overall_index_ch:overall_index_ch+len(bad_ev_id_index)] | bad_ev_id_index

      overall_index_p += len(ev_ids_p)
      run_index_p += 1
      overall_index_h += len(ev_ids_h)
      run_index_h += 1
      overall_index_ch += len(ev_ids_ch)
      run_index_ch += 1

    # DEBUGGING CODE ONLY
    # Test random pairing rate by pairing incorrect events
    #bad_index_p[0] = bad_index_p[1] = True
    #bad_index_h[-1] = bad_index_h[0] = True
    #bad_index_ch[-2] = bad_index_ch[-1] = True

    [df_fjparticles_p["run_number"].pop(i) and df_fjparticles_p["ev_id"].pop(i) and \
     df_fjparticles_p["fj_particle"].pop(i) for i in range(len(bad_index_p)-1, -1, -1) \
     if bad_index_p[i]]
    [df_fjparticles_h["fj_particle"].pop(i) for i in range(len(bad_index_h)-1, -1, -1) \
     if bad_index_h[i]]
    [df_fjparticles_ch["fj_particle"].pop(i) for i in range(len(bad_index_ch)-1, -1, -1) \
     if bad_index_ch[i]]

    df_fjparticles_p["fj_particles_p"] = df_fjparticles_p["fj_particle"]
    df_fjparticles_p.pop("fj_particle")
    df_fjparticles_p["fj_particles_h"] = df_fjparticles_h["fj_particle"]
    df_fjparticles_p["fj_particles_ch"] = df_fjparticles_ch["fj_particle"]

    if self.debug_level > 1:
      print("number of p events:", len(df_fjparticles_p["fj_particles_p"]))
      print("number of h events:", len(df_fjparticles_p["fj_particles_h"]))
      print("number of ch events:", len(df_fjparticles_p["fj_particles_ch"]))

    return df_fjparticles_p

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
  def analyze_events(self, MPI):

    df = getattr(self, "df_fjparticles_MPI%s" % MPI)

    # Fill track histograms
    for level in ['p', 'h', 'ch']:
      for fj_particles in df['fj_particles_%s' % level]:
        self.fill_track_histograms(fj_particles, level, MPI)

    self.event_number = 0

    # match jets -- and fill histograms
    for fj_particles_p, fj_particles_h, fj_particles_ch in zip(
        df['fj_particles_p'], df['fj_particles_h'], df['fj_particles_ch']):

      self.analyze_event(fj_particles_p, fj_particles_h, fj_particles_ch, MPI)

    if self.debug_level > 2:
      for attr in dir(self):
        obj = getattr(self, attr)
        print('size of {}: {}'.format(attr, sys.getsizeof(obj)))

  #---------------------------------------------------------------
  # Fill track histograms.
  #---------------------------------------------------------------
  def fill_track_histograms(self, fj_particles, level, MPI):

    # Check that the entries exist appropriately
    # (need to check if indeed this happens for MC, and how...)
    #if type(fj_particles) != fj.vectorPJ:
    #  return

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
    #if type(fj_particles_p) != fj.vectorPJ or type(fj_particles_h) != fj.vectorPJ or \
    #   type(fj_particles_ch) != fj.vectorPJ:
    #  print('fj_particles type mismatch -- skipping event')
    #  return

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
    for jet_ch in jets_ch_selected:
      self.fill_level_before_matching(jet_ch, jetR, 'ch', MPI)

    # Loop through jets and set jet matching candidates for each jet in user_info
    name_suffix = "MPI%s_R%s" % (MPI, str(jetR))
    for jet_h in jets_h_selected:
      for jet_p in jets_p_selected:
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

        # Check for weird events
        if abs(jet_pt_p_ungroomed - jet_pt_h_ungroomed) > 100:
          print("\nweird event!")
          print("p jet momentum: (%s, %s, %s)" % (jet_p.px(), jet_p.py(), jet_p.pz()), "pT:", jet_p.pt())
          print("h jet momentum: (%s, %s, %s)" % (jet_h.px(), jet_h.py(), jet_h.pz()), "pT:", jet_h.pt())
          print("delta R: %s * jetR" % (jet_p.delta_R(jet_h) / jetR) )
          #exit()

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
  # Apply appropriate scaling for pT-hat bin
  #---------------------------------------------------------------
  def scale_objects(self):

    for attr in dir(self):
      obj = getattr(self, attr)
      types = (ROOT.TH1, ROOT.TH2, ROOT.TH3, ROOT.THnBase)
      if isinstance(obj, types):
        name = obj.GetName()
        if "Nevents" not in name:
          if "MPIon" in name:
            obj.Scale(self.xsec_MPIon / self.nEvents_MPIon)
          else:  # "MPIoff" in name
            obj.Scale(self.xsec_MPIoff / self.nEvents_MPIoff)
          obj.SetNameTitle(name+"Scaled", name+"Scaled")

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

    if self.make_th3s:
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

    if self.make_th3s and obs_level_1 > 0:
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
