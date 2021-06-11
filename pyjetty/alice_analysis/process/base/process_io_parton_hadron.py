#!/usr/bin/env python3

"""
  Analysis IO class for jet analysis with MC track dataframes at parton/hadron level.
  Each instance of the class handles the IO of a *single* track tree.
  
  Authors: Ezra Lesser
           James Mulligan
           Mateusz Ploskon
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
  # - level is either 'p' or 'h' (parton or hadron)
  # - MPI is either 'on' or 'off'
  #---------------------------------------------------------------
  def __init__(self, input_file='', tree_name_base='tree_Particle_gen',
               level='p', MPI='on', **kwargs):

    super(ProcessIO, self).__init__(**kwargs)

    # Input ROOT files containing generated events with jet observables
    self.input_file = input_file

    self.level = level

    if MPI not in ['on', 'off']:
      raise ValueError("MPI must either be 'on' or 'off'")

    # Name of the trees containing the track information
    self.tree_name = '_'.join((tree_name_base, level, "MPI%s" % MPI))

    self.reset_dataframes()

    # Name of the columns to load from each respective tree
    if level == 'p':
      self.columns = ['run_number', 'ev_id', 'ParticleE',
                      'ParticlePx', 'ParticlePy', 'ParticlePz']
    elif level == 'h':
      self.columns = ['run_number', 'ev_id', 'ParticleE', 'ParticlePx',
                        'ParticlePy', 'ParticlePz', 'is_charged']
    else:
      raise ValueError("Particle level %s not recognized / use either 'p' or 'h'")

    # Set the combination of fields that give a unique event id
    self.unique_identifier =  ['run_number', 'ev_id']

    self.load_xsec_Nev(MPI)

  #---------------------------------------------------------------
  # Clear dataframes
  #---------------------------------------------------------------
  def reset_dataframes(self):
    self.track_df = None
    self.run_numbers = None
    self.unique_ev_ids_per_run = None

  #---------------------------------------------------------------
  # Load cross section and number of events
  #---------------------------------------------------------------
  def load_xsec_Nev(self, MPI):

    # Load the histograms containing cross section and N_ev information
    with uproot.open(self.input_file) as in_f:
      h_xsec = in_f["hxsec_MPI%s" % MPI]
      h_Nev = in_f["hNev_MPI%s" % MPI]

      # Save contents
      self.xsec = h_xsec.values(flow=False)[0]
      self.Nev = h_Nev.values(flow=False)[0]

  #---------------------------------------------------------------
  # Convert ROOT TTree to SeriesGroupBy object of fastjet particles per event.
  # Optionally, remove a certain random fraction of tracks
  #---------------------------------------------------------------
  def load_data(self, group_by_evid=True, ch_cut=False):

    self.reset_dataframes()

    print('    tree_name = {}'.format(self.tree_name))
    self.track_df = self.load_dataframe()

    return self.group_fjparticles(group_by_evid, ch_cut)

  #---------------------------------------------------------------
  # Convert ROOT TTree to pandas dataframe
  # Return merged track+event dataframe from a given input file
  # Returned dataframe has one row per jet constituent:
  #     run_number, ev_id, ParticlePt, ParticleEta, ParticlePhi
  #---------------------------------------------------------------
  def load_dataframe(self):

    df = None
    with uproot.open(self.input_file)[self.tree_name] as tree:
      if not tree:
        raise ValueError('Tree {} not found in file {}'.format(self.tree_name, self.input_file))
      # Pandas DataFrame implementation
      #df = uproot.concatenate(tree, self.columns, library="pd")

      # Try saving memory with numpy implementation
      df = uproot.concatenate(tree, self.columns, library="np")
      # Each value is a 2D array for some reason, so fix that
      for key, value in df.items():
        if key == "run_number" or key == "ev_id" or key == "is_charged":
          df[key] = value[0].astype('int32')
        else:
          df[key] = value[0]

    return df

  #---------------------------------------------------------------
  # Transform the track dataframe into a SeriesGroupBy object
  # of fastjet particles per event.
  #---------------------------------------------------------------
  def group_fjparticles(self, group_by_evid=True, ch_cut=False):

    track_df = self.track_df

    if ch_cut:
      if self.level == 'p':
        raise ValueError("ch_cut cannot be set for parton-level tree")
      else:  # self.level == 'h'
        #track_df = track_df.loc[track_df["is_charged"] == True]
        for key, value in track_df.items():
          if key == "is_charged":
            continue
          track_df[key] = value[track_df["is_charged"] == True]
        track_df.pop("is_charged")  # This info is now redundant

    df_fjparticles = None
    fj_particles = self.get_fjparticles(track_df)

    if group_by_evid:
      ''' Pandas implementation
      # Transform the track dataframe into a series object of fastjet particles per event

      # (i) Group the track dataframe by event
      #     track_df_grouped is a DataFrameGroupBy object with one track dataframe per event
      track_df_grouped = track_df.groupby(self.unique_identifier)
    
      # (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet particles
      df_fjparticles = track_df_grouped.apply(self.get_fjparticles)
      '''
      # First split by run_number, then by ev_id
      self.run_numbers = np.unique(track_df["run_number"], return_index=True)
      ev_ids_per_run = np.split(track_df["ev_id"], self.run_numbers[1][1:])
      self.unique_ev_ids_per_run = [np.unique(ev_ids_per_run[i], return_index=True) \
                                    for i in range(len(ev_ids_per_run))]
      split_points = np.concatenate([self.unique_ev_ids_per_run[i][1] + self.run_numbers[1][i] \
                                     for i in range(len(ev_ids_per_run))])[1:]
      df_fjparticles = {
        "run_number": np.concatenate(
          [[self.run_numbers[0][i]] * len(self.unique_ev_ids_per_run[i][0]) \
           for i in range(len(ev_ids_per_run))]),
        "ev_id": np.concatenate(
          [self.unique_ev_ids_per_run[i][0] for i in range(len(ev_ids_per_run))]),
        "fj_particle": np.split(fj_particles, split_points) }

    else:
      # Transform the track dataframe into a dataframe of fastjet particles per track
      ''' Pandas implementation
      df_fjparticles = pandas.DataFrame( 
        {"run_number": track_df["run_number"], "ev_id": track_df["ev_id"], 
         "fj_particle": self.get_fjparticles(track_df)} )
      '''
      df_fjparticles = {
        "run_number": track_df["run_number"], "ev_id": track_df["ev_id"], 
        "fj_particle": fj_particles }

    #print(df_fjparticles)
    #for i in range(1, len(track_df["ev_id"])):
    #  if track_df["ev_id"][i] - track_df["ev_id"][i-1] > 1:
    #    print(i, track_df["ev_id"][i], track_df["ev_id"][i-1])
    #    exit()
    #print(len(df_fjparticles["run_number"]), len(df_fjparticles["ev_id"]),
    #      len(df_fjparticles["fj_particle"]))
    #print(df_fjparticles["run_number"][-1], df_fjparticles["ev_id"][-1],
    #      df_fjparticles["fj_particle"][-1])
    #exit()
    return df_fjparticles

  #---------------------------------------------------------------
  # Return fastjet:PseudoJets from a given track dataframe
  #---------------------------------------------------------------
  def get_fjparticles(self, df_tracks):

    ''' Pandas implementation
    return fjext.vectorize_px_py_pz_e(df_tracks['ParticlePx'].values, df_tracks['ParticlePy'].values,
                                      df_tracks['ParticlePz'].values, df_tracks['ParticleE'].values)
    '''
    return fjext.vectorize_px_py_pz_e(df_tracks['ParticlePx'], df_tracks['ParticlePy'],
                                      df_tracks['ParticlePz'], df_tracks['ParticleE'])
