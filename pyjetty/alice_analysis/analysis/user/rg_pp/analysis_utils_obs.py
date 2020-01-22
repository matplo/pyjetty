#!/usr/bin/env python3

"""
  Analysis utilities for Soft Drop jet analysis with track dataframe.
  
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
from array import *
import ROOT

# Base class
from pyjetty.alice_analysis.analysis.base import analysis_utils

################################################################
class analysis_utils_obs(analysis_utils.analysis_utils):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(analysis_utils_obs, self).__init__(**kwargs)

  #---------------------------------------------------------------
  # Get observable label
  #---------------------------------------------------------------
  def obs_label(self, observable, obs_setting):
  
    if observable == 'theta_g':
      zcut = obs_setting[0]
      beta = obs_setting[1]
      obs_label = 'zcut{}_B{}'.format(self.remove_periods(zcut), beta)
      
    elif observable == 'zg':
      zcut = obs_setting[0]
      beta = obs_setting[1]
      obs_label = 'zcut{}_B{}'.format(self.remove_periods(zcut), beta)

    elif observable == 'subjet_z':
      obs_label = '{}'.format(obs_setting)
      
    return obs_label

  #---------------------------------------------------------------
  # Get regularization parameter
  #---------------------------------------------------------------
  def get_reg_param(self, obs_settings, obs_config_list, obs_config_dict, obs_label, observable, jetR):
    
    for i, obs_setting in enumerate(obs_settings):
      
      if self.obs_label(observable, obs_setting) == obs_label:
        
        config_name = obs_config_list[i]
        
        if observable == 'theta_g' or observable == 'zg':
          reg_param = obs_config_dict[config_name]['reg_param'][observable][jetR]
        else:
          reg_param = obs_config_dict[config_name]['reg_param'][jetR]
        #print('reg_param for {} {} jetR={}: {}'.format(obs_label, observable, jetR, reg_param))
          
        return reg_param
      
      else:
        continue
