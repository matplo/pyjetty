#!/usr/bin/env python3

"""
  Analysis IO class for jet analysis with track dataframe.
  Each instance of the class handles the IO of a *single* track tree.
  
  Authors: James Mulligan
           Mateusz Ploskon
           Ezra Lesser
"""

from __future__ import print_function

import os   # for creating file on output

# Data analysis and plotting
import uproot
import pandas
import numpy as np

# Fastjet via python (from external library fjpydev)
import fastjet as fj
import fjext

# Base class
from pyjetty.alice_analysis.process.base import common_base

################################################################
class process_io(common_base.common_base):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', output_dir='', track_tree_name='tree_Particle', **kwargs):
    super(process_io, self).__init__(**kwargs)
    self.input_file = input_file
    self.output_dir = output_dir
    if self.output_dir[-1] != '/':
      self.output_dir += '/'
    self.track_tree_name = track_tree_name
    self.event_tree_name = 'PWGHF_TreeCreator/tree_event_char'
    self.event_columns = ['run_number', 'ev_id', 'z_vtx_reco', 'is_ev_rej']
    self.reset_dataframes()
  
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
  # Convert ROOT TTree to SeriesGroupBy object of fastjet particles per event.
  # Optionally, remove a certain random fraction of tracks
  #---------------------------------------------------------------
  def load_data(self, reject_tracks_fraction=0., offset_indices=False, group_by_evid=True):
    
    self.reject_tracks_fraction = reject_tracks_fraction
    self.reset_dataframes()

    print('Convert ROOT trees to pandas dataframes...')
    print('    track_tree_name = {}'.format(self.track_tree_name))

    
    self.track_df = self.load_dataframe()
    
    if self.reject_tracks_fraction > 1e-3:
      n_remove = int(reject_tracks_fraction * len(self.track_df.index))
      print('    Removing {} of {} tracks from {}'.format(n_remove, len(self.track_df.index), self.track_tree_name))
      np.random.seed()
      indices_remove = np.random.choice(self.track_df.index, n_remove, replace=False)
      self.track_df.drop(indices_remove, inplace=True)

    df_fjparticles = self.group_fjparticles(offset_indices, group_by_evid)

    return df_fjparticles
  
  #---------------------------------------------------------------
  # Convert ROOT TTree to pandas dataframe
  # Return merged track+event dataframe from a given input file
  # Returned dataframe has one row per jet constituent:
  #     run_number, ev_id, ParticlePt, ParticleEta, ParticlePhi
  #---------------------------------------------------------------
  def load_dataframe(self):
    
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
  # Opposite operation as load_dataframe above. Takes a dataframe
  # with the same formatting and saves to class's output_file.
  #---------------------------------------------------------------
  def save_dataframe(self, df, filename, histograms=[]):

    # Create output directory if it does not already exist
    if not os.path.exists(self.output_dir):
      os.makedirs(self.output_dir)

    # Open output directory and (re)create rootfile
    with uproot.recreate(self.output_dir + filename) as f:

      branchdict = {"run_number": int, "ev_id": int, "ParticlePt": float,
                    "ParticleEta": float, "ParticlePhi": float}
      # Subdirectories not yet implemented in uproot
      #title = "PWGHF_TreeCreator/tree_Particle"
      title = "tree_Particle"
      f[title] = uproot.newtree(branchdict, title=title)
      f[title].extend( { "run_number": df["run_number"], "ev_id": df["ev_id"], 
                         "ParticlePt": df["ParticlePt"], "ParticleEta": df["ParticleEta"],
                         "ParticlePhi": df["ParticlePhi"] } )

      # Write histograms to file too, if any are passed
      for h in histograms:
        f[h.name] = h

  #---------------------------------------------------------------
  # Transform the track dataframe into a SeriesGroupBy object
  # of fastjet particles per event.
  #---------------------------------------------------------------
  def group_fjparticles(self, offset_indices=False, group_by_evid=True):

    if group_by_evid:
      print("Transform the track dataframe into a series object of fastjet particles per event...")

      # (i) Group the track dataframe by event
      #     track_df_grouped is a DataFrameGroupBy object with one track dataframe per event
      self.track_df_grouped = self.track_df.groupby(['run_number','ev_id'])
    
      # (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet particles
      self.df_fjparticles = self.track_df_grouped.apply(self.get_fjparticles, offset_indices=offset_indices)
    
    else:
      print("Transform the track dataframe into a dataframe of fastjet particles per track...")

      # Transform into a DataFrame of fastjet particles
      df = self.track_df
      self.df_fjparticles = pandas.DataFrame( {"run_number": df["run_number"], "ev_id": df["ev_id"], 
                                               "fj_particle": self.get_fjparticles(self.track_df)} )

    return self.df_fjparticles

  #---------------------------------------------------------------
  # Return fastjet:PseudoJets from a given track dataframe
  #---------------------------------------------------------------
  def get_fjparticles(self, df_tracks, offset_indices=False):
    
    # If offset_indices is true, then offset the user_index by a large negative value
    user_index_offset = 0
    if offset_indices:
        user_index_offset = int(-1e6)
    
    # Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of pt,eta,phi
    fj_particles = fjext.vectorize_pt_eta_phi(df_tracks['ParticlePt'].values, df_tracks['ParticleEta'].values, df_tracks['ParticlePhi'].values, user_index_offset)
    
    return fj_particles
