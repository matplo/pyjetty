#!/usr/bin/env python3

"""
  Analysis IO class for jet analysis with track dataframe.
  The class stores a list of Pb-Pb files, and keeps track of
  a current file and current event -- and returns the current
  event when requested.
  
  Authors: James Mulligan
           Mateusz Ploskon
"""

from __future__ import print_function

# Data analysis and plotting
import uproot
import pandas
import numpy as np
import random

# Base class
from pyjetty.alice_analysis.process.base import common_base

# Main file IO class
from pyjetty.alice_analysis.process.base import process_io

################################################################
class ProcessIO_Emb(common_base.CommonBase):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, file_list='PbPb_file_list.txt', track_tree_name='tree_Particle',
               min_cent=0., max_cent=10., is_pp = False, use_ev_id_ext = True,
               m=0.1396, remove_used_file=True, **kwargs):
    super(ProcessIO_Emb, self).__init__(**kwargs)
    
    self.file_list = file_list
    self.track_tree_name = track_tree_name
    
    self.list_of_files = []
    with open(self.file_list) as f:
      files = [fn.strip() for fn in f.readlines()]
    list_of_files = list(filter(None, files))
    
    # Choose N random files to keep in the list
    n_files = 1000
    random.shuffle(list_of_files)
    self.list_of_files = list_of_files[0:n_files]

    self.current_file_df = None
    self.current_file_nevents = 0
    self.current_event_index = 0
    
    self.min_centrality = min_cent
    self.max_centrality = max_cent
    self.m = m
    
    self.is_pp = is_pp
    self.use_ev_id_ext = use_ev_id_ext
    self.remove_used_file = remove_used_file
          
    # Initialize by loading a file
    self.load_file()
    
  #---------------------------------------------------------------
  # Return current event in current file, and increment current_event.
  # If there are no more events in the current file, open a new file.
  #---------------------------------------------------------------
  def load_event(self):
      
    if self.current_event_index >= self.current_file_nevents:
        self.load_file()
        
    current_event = self.current_file_df.iloc[self.current_event_index]
    self.current_event_index += 1
    #print('Get Pb-Pb event {}/{}'.format(self.current_event_index, self.current_file_nevents))
    return current_event

  #---------------------------------------------------------------
  # Pick a random file from the file list, load it as the
  # current file as a dataframe, and remove it from the file list.
  #---------------------------------------------------------------
  def load_file(self):
      
    input_file = random.choice(self.list_of_files)
    if self.remove_used_file:
      self.list_of_files.remove(input_file)
    print('Opening Pb-Pb file: {}'.format(input_file))

    io = process_io.ProcessIO(input_file=input_file, track_tree_name=self.track_tree_name,
                              is_pp=self.is_pp, min_cent=self.min_centrality,
                              max_cent=self.max_centrality, use_ev_id_ext=self.use_ev_id_ext)
    self.current_file_df = io.load_data(m=self.m, offset_indices=True)
    self.current_file_nevents = len(self.current_file_df.index)
    self.current_event_index = 0
