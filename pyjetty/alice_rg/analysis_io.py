#!/usr/bin/env python3

"""
  Analysis IO class for jet analysis with track dataframe.
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import math

# Data analysis and plotting
import uproot
import pandas
import numpy as np
import ROOT

# Fastjet via python (from external library fjpydev)
import fastjet as fj
import fjext

# Base class
import analysis_base

class analysis_io(analysis_base.analysis_base):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(analysis_io, self).__init__(**kwargs)
    self.event_tree_name = 'PWGHF_TreeCreator/tree_event_char'
    self.event_columns = ['run_number', 'ev_id', 'z_vtx_reco','is_ev_rej']
    self.reset_dataframes()
    print(self)
  
  #---------------------------------------------------------------
  # Clear dataframes
  #---------------------------------------------------------------
  def reset_dataframes(self):
    self.event_tree = None
    self.event_df_orig = None
    self.event_df = None
    self.track_tree = None
    self.track_df_orig = None
    self.track_df = None
    self.track_df_grouped = None
    self.df_fjparticles = None

  #---------------------------------------------------------------
  # Convert ROOT TTree to pandas dataframe
  # Return merged track+event dataframe from a given input file
  # Returned dataframe has one row per jet constituent:
  #     run_number, ev_id, ParticlePt, ParticleEta, ParticlePhi
  #---------------------------------------------------------------
  def load_dataframe(self, input_file, track_tree_name='tree_Particle'):
    
    self.input_file = input_file
    self.track_tree_name = track_tree_name
    self.reset_dataframes()
    
    # Load event tree into dataframe, and apply event selection
    self.event_tree = uproot.open(self.input_file)[self.event_tree_name]
    if not self.event_tree:
      print('Tree {} not found in file {}'.format(self.event_tree_name, self.input_file))
    self.event_df_orig = self.event_tree.pandas.df(self.event_columns)
    self.event_df_orig.reset_index(drop=True)
    self.event_df = self.event_df_orig.query('is_ev_rej == 0')
    self.event_df.reset_index(drop=True)

    # Load track tree into dataframe
    self.track_tree_name = 'PWGHF_TreeCreator/{}'.format(self.track_tree_name)
    self.track_tree = uproot.open(self.input_file)[self.track_tree_name]
    if not self.track_tree:
      print('Tree {} not found in file {}'.format(self.track_tree_name, self.input_file))
    self.track_df_orig = self.track_tree.pandas.df()

    # Merge event info into track tree
    self.track_df = pandas.merge(self.track_df_orig, self.event_df, on=['run_number', 'ev_id'])
    return self.track_df

  #---------------------------------------------------------------
  # Transform the track dataframe into a SeriesGroupBy object
  # of fastjet particles per event.
  #---------------------------------------------------------------
  def group_fjparticles(self):

    print('Transform the track dataframe into a series object of fastjet particles per event...')

    # (i) Group the track dataframe by event
    #     track_df_grouped is a DataFrameGroupBy object with one track dataframe per event
    self.track_df_grouped = self.track_df.groupby(['run_number','ev_id'])
    
    # (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet particles
    self.df_fjparticles = self.track_df_grouped.apply(self.get_fjparticles)
    
    return self.df_fjparticles

  #---------------------------------------------------------------
  # Return fastjet:PseudoJets from a given track dataframe
  #---------------------------------------------------------------
  def get_fjparticles(self, df_tracks):
    
    # Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of pt,eta,phi
    fj_particles = fjext.vectorize_pt_eta_phi(df_tracks['ParticlePt'].values, df_tracks['ParticleEta'].values, df_tracks['ParticlePhi'].values)
    
    return fj_particles
