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
import sys

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
class ProcessIO(common_base.CommonBase):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', tree_dir='PWGHF_TreeCreator',
               track_tree_name='tree_Particle', event_tree_name='tree_event_char',
               output_dir='', is_pp=True, min_cent=0., max_cent=10.,
               use_ev_id_ext=False, **kwargs):
    super(ProcessIO, self).__init__(**kwargs)
    self.input_file = input_file
    self.output_dir = output_dir
    self.tree_dir = tree_dir
    if len(tree_dir) and tree_dir[-1] != '/':
      self.tree_dir += '/'
    self.track_tree_name = track_tree_name
    self.event_tree_name = event_tree_name
    if len(output_dir) and output_dir[-1] != '/':
      self.output_dir += '/'
    self.reset_dataframes()
    
    self.use_ev_id_ext = use_ev_id_ext
    if self.use_ev_id_ext:
      self.unique_identifier =  ['run_number', 'ev_id', 'ev_id_ext']
    else:
      self.unique_identifier =  ['run_number', 'ev_id']
    
    self.is_pp = is_pp
    if self.is_pp:
      if self.use_ev_id_ext:
        self.event_columns = ['run_number', 'ev_id', 'ev_id_ext', 'z_vtx_reco', 'is_ev_rej']
      else:
        self.event_columns = ['run_number', 'ev_id', 'z_vtx_reco', 'is_ev_rej']
    else:
      if self.use_ev_id_ext:
        self.event_columns = ['run_number', 'ev_id', 'ev_id_ext', 'z_vtx_reco', 'is_ev_rej', 'centrality']
      else:
        self.event_columns = ['run_number', 'ev_id', 'z_vtx_reco', 'is_ev_rej', 'centrality']
      self.min_centrality = min_cent
      self.max_centrality = max_cent
  
  #---------------------------------------------------------------
  # Clear dataframes
  #---------------------------------------------------------------
  def reset_dataframes(self):
    self.event_df_orig = None
    self.track_df = None
  
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
    event_tree = None
    event_df = None
    event_tree_name = self.tree_dir + self.event_tree_name
    event_tree = uproot.open(self.input_file)[event_tree_name]
    if not event_tree:
      sys.exit('Tree {} not found in file {}'.format(event_tree_name, self.input_file))
    self.event_df_orig = event_tree.pandas.df(self.event_columns)
    self.event_df_orig.reset_index(drop=True)
    if self.is_pp:
      event_criteria = 'is_ev_rej == 0'
    else:
      event_criteria = 'is_ev_rej == 0 and centrality > @self.min_centrality and centrality < @self.max_centrality'
    event_df = self.event_df_orig.query(event_criteria)
    event_df.reset_index(drop=True)

    # Load track tree into dataframe
    track_tree = None
    track_df_orig = None
    track_tree_name = self.tree_dir + self.track_tree_name
    track_tree = uproot.open(self.input_file)[track_tree_name]
    if not track_tree:
      sys.exit('Tree {} not found in file {}'.format(track_tree_name, self.input_file))
    track_df_orig = track_tree.pandas.df()

    # Merge event info into track tree
    self.track_df = pandas.merge(track_df_orig, event_df, on=self.unique_identifier)
    return self.track_df

  #---------------------------------------------------------------
  # Opposite operation as load_dataframe above. Takes a dataframe
  # with the same formatting and saves to class's output_file.
  # histograms is list of tuples: [ ("title", np.histogram), ... ]
  #---------------------------------------------------------------
  def save_dataframe(self, filename, df, df_true=False, histograms=[]):

    # Create output directory if it does not already exist
    if not os.path.exists(self.output_dir):
      os.makedirs(self.output_dir)

    # Open output directory and (re)create rootfile
    with uproot.recreate(self.output_dir + filename) as f:

      if df_true:
        # Create tree with truth particle info
        title = 'tree_Particle_gen'
        branchdict = {"run_number": int, "ev_id": int, "ParticlePt": float,
                      "ParticleEta": float, "ParticlePhi": float}
        print("Length of truth track tree: %i" % len(self.track_df))
        f[title] = uproot.newtree(branchdict, title=title)
        f[title].extend( { "run_number": self.track_df["run_number"],
                           "ev_id": self.track_df["ev_id"], 
                           "ParticlePt": self.track_df["ParticlePt"],
                           "ParticleEta": self.track_df["ParticleEta"],
                           "ParticlePhi": self.track_df["ParticlePhi"] } )

      # Create tree with detector-level particle info
      title = 'tree_Particle'
      branchdict = {"run_number": int, "ev_id": int, "ParticlePt": float,
                    "ParticleEta": float, "ParticlePhi": float}
      print("Length of detector-level track tree: %i" % len(df))
      f[title] = uproot.newtree(branchdict, title=title)
      f[title].extend( { "run_number": df["run_number"], "ev_id": df["ev_id"], 
                         "ParticlePt": df["ParticlePt"], "ParticleEta": df["ParticleEta"],
                         "ParticlePhi": df["ParticlePhi"] } )

      # Create tree with event char
      title = self.event_tree_name
      branchdict = {"is_ev_rej": int, "run_number": int, "ev_id": int, "z_vtx_reco": float}
      f[title] = uproot.newtree(branchdict, title=title)
      f[title].extend( {"is_ev_rej": self.event_df_orig["is_ev_rej"], 
                        "run_number": self.event_df_orig["run_number"], 
                        "ev_id": self.event_df_orig["ev_id"],
                        "z_vtx_reco": self.event_df_orig["z_vtx_reco"] } )
        
      # Write hNevents histogram: number of accepted events at detector level
      f["hNevents"] = ( np.array([ 0, df["ev_id"].nunique() ]), np.array([ -0.5, 0.5, 1.5 ]) )

      # Write histograms to file too, if any are passed
      for title, h in histograms:
        f[title] = h

  #---------------------------------------------------------------
  # Transform the track dataframe into a SeriesGroupBy object
  # of fastjet particles per event.
  #---------------------------------------------------------------
  def group_fjparticles(self, offset_indices=False, group_by_evid=True):

    if group_by_evid:
      print("Transform the track dataframe into a series object of fastjet particles per event...")

      # (i) Group the track dataframe by event
      #     track_df_grouped is a DataFrameGroupBy object with one track dataframe per event
      track_df_grouped = None
      track_df_grouped = self.track_df.groupby(self.unique_identifier)
    
      # (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet particles
      df_fjparticles = None
      df_fjparticles = track_df_grouped.apply(self.get_fjparticles, offset_indices=offset_indices)
    
    else:
      print("Transform the track dataframe into a dataframe of fastjet particles per track...")

      # Transform into a DataFrame of fastjet particles
      df = self.track_df
      df_fjparticles = pandas.DataFrame( {"run_number": df["run_number"], "ev_id": df["ev_id"],
                                               "fj_particle": self.get_fjparticles(self.track_df)} )

    return df_fjparticles

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
