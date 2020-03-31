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
  # Get grooming settings (i.e. list that stores a dict of grooming
  # settings, i.e. {'sd': [zcut, beta]} or {'dg': [a]} from observable config block
  #---------------------------------------------------------------
  def grooming_settings(self, obs_config_dict):
  
    grooming_settings = []
    for config_key, subconfig in obs_config_dict.items():
      if config_key == 'common_settings':
        continue
      if 'SoftDrop' in subconfig:
        grooming_dict = obs_config_dict[config_key]['SoftDrop']
        grooming_settings.append({'sd':[grooming_dict['zcut'], grooming_dict['beta']]})
      elif 'DynamicalGrooming' in subconfig:
        grooming_dict = obs_config_dict[config_key]['DynamicalGrooming']
        grooming_settings.append({'dg':[grooming_dict['a']]})
      else:
        grooming_settings.append(None)
        
    return grooming_settings

  # Get formatted grooming label from grooming_setting,
  # i.e. {'sd': [zcut, beta]} or {'dg': [a]}
  #---------------------------------------------------------------
  def grooming_label(self, grooming_setting):
  
    key, value = list(grooming_setting.items())[0]

    if key == 'sd':
      text = 'zcut{}_B{}'.format(self.remove_periods(value[0]), value[1])
    elif key == 'dg':
      text = 'dg{}'.format(value[0])
    else:
      sys.exit('Unknown grooming type!')
    
    return text

  #---------------------------------------------------------------
  # Remove periods from a label
  #---------------------------------------------------------------
  def remove_periods(self, text):
    
    string = str(text)
    return string.replace('.', '')
