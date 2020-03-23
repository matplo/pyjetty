#!/usr/bin/env python3

"""
  Common utilities for process and analysis for jet analysis with track dataframe.
  
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
from pyjetty.alice_analysis.process.base import common_base

################################################################
class CommonUtils(common_base.CommonBase):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(CommonUtils, self).__init__(**kwargs)
  
  #---------------------------------------------------------------
  # Get SD settings (i.e. list that stores a list of SD settings [zcut, beta])
  # from observable config block
  #---------------------------------------------------------------
  def sd_settings(self, obs_config_dict):
  
    sd_settings = []
    for config_key, subconfig in obs_config_dict.items():
      if config_key == 'common_settings':
        continue
      if 'SoftDrop' in subconfig:
        sd_dict = obs_config_dict[config_key]['SoftDrop']
        sd_settings.append([sd_dict['zcut'], sd_dict['beta']])
      else:
        sd_settings.append(None)
        
    return sd_settings

  # Get formatted Soft Drop label from sd_setting = [zcut, beta]
  #---------------------------------------------------------------
  def sd_label(self, sd_setting):
  
      zcut = sd_setting[0]
      beta = sd_setting[1]
      sd_label = 'zcut{}_B{}'.format(self.remove_periods(zcut), beta)
      return sd_label

  #---------------------------------------------------------------
  # Remove periods from a label
  #---------------------------------------------------------------
  def remove_periods(self, text):
    
    string = str(text)
    return string.replace('.', '')
