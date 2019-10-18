#!/usr/bin/env python3

"""
  Analysis task base class.
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import time

# Data analysis and plotting
import ROOT
import yaml

# Analysis utilities
from pyjetty.alice_analysis.process.base import common_base
from pyjetty.alice_analysis.process.base import process_utils
from pyjetty.mputils import treewriter

################################################################
class process_base(common_base.common_base):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
    super(process_base, self).__init__(**kwargs)
    self.input_file = input_file
    self.config_file = config_file
    self.output_dir = output_dir
    self.debug_level = debug_level # (0 = no debug info, 1 = some debug info, 2 = all debug info)
    
    # Create output dir
    if not self.output_dir.endswith("/"):
      self.output_dir = self.output_dir + "/"
    if not os.path.exists(self.output_dir):
      os.makedirs(self.output_dir, 775)

    # Initialize utils class
    self.utils = process_utils.process_utils()
  
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
  
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    
    self.jetR_list = config['jetR']
    self.debug_level = config['debug_level']

    # Check if constituent subtractor is included, and initialize it if so
    self.do_constituent_subtraction = False
    if 'constituent_subtractor' in config:
      print('Constituent subtractor is enabled.')
      self.do_constituent_subtraction = True
      constituent_subtractor = config['constituent_subtractor']
      
      self.max_distance = constituent_subtractor['max_distance']
      self.alpha = constituent_subtractor['alpha']
      self.max_eta = constituent_subtractor['max_eta']
      self.bge_rho_grid_size = constituent_subtractor['bge_rho_grid_size']
      self.max_pt_correct = constituent_subtractor['max_pt_correct']
      self.ghost_area = constituent_subtractor['ghost_area']
    else:
      print('Constituent subtractor is disabled.')

  #---------------------------------------------------------------
  # Save all histograms
  #---------------------------------------------------------------
  def save_output_objects(self):
    
    outputfilename = os.path.join(self.output_dir, 'AnalysisResults.root')
    fout = ROOT.TFile(outputfilename, 'recreate')
    fout.cd()
    for attr in dir(self):
      
      obj = getattr(self, attr)
      
      # If tree writer object, get the tree it contains
      if isinstance(obj, treewriter.RTreeWriter):
        obj = obj.tree
      
      # Write all ROOT histograms and trees to file
      types = (ROOT.TH1, ROOT.THnBase, ROOT.TTree)
      if isinstance(obj, types):
        obj.Write()
  
    fout.Close()
