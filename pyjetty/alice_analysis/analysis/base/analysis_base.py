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
from pyjetty.alice_analysis.analysis.base import analysis_utils
#import base
#import analysis_utils

################################################################
class AnalysisBase(common_base.CommonBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file_data='', input_file_response='', config_file='', output_dir='', file_format='', **kwargs):
    super(AnalysisBase, self).__init__(**kwargs)
    self.input_file_data = input_file_data
    self.input_file_response = input_file_response
    self.config_file = config_file
    self.output_dir = output_dir
    self.file_format = file_format
    
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
