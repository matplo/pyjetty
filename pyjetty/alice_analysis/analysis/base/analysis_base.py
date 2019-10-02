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
from pyjetty.alice_analysis.analysis.base import base
from pyjetty.alice_analysis.analysis.base import analysis_utils

################################################################
class analysis_base(base.base):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
    super(analysis_base, self).__init__(**kwargs)
    self.input_file = input_file
    self.config_file = config_file
    self.output_dir = output_dir
    self.debug_level = debug_level # (0 = no debug info, 1 = some debug info, 2 = all debug info)
    
    # Create output dir
    if not self.output_dir.endswith("/"):
      self.output_dir = self.output_dir + "/"
    if not os.path.exists(self.output_dir):
      os.makedirs(self.output_dir)

    # Initialize utils class
    self.utils = analysis_utils.analysis_utils()
  
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
  
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    
    self.jetR_list = config['jetR']
    self.debug_level = config['debug_level']
