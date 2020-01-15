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
class analysis_utils_sd(analysis_utils.analysis_utils):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(analysis_utils_sd, self).__init__(**kwargs)

  #---------------------------------------------------------------
  # Get regularization parameter
  #---------------------------------------------------------------
  def get_reg_param(self, sd_settings, sd_config_list, sd_config_dict, sd_label, observable, jetR):
    
    for i, sd_setting in enumerate(sd_settings):
      
      zcut = sd_setting[0]
      beta = sd_setting[1]
      label = 'zcut{}_B{}'.format(self.remove_periods(zcut), beta)
      if label == sd_label:
        
        config_name = sd_config_list[i]
        
        reg_param = sd_config_dict[config_name]['reg_param'][observable][jetR]
        #print('reg_param for {} {} jetR={}: {}'.format(sd_label, observable, jetR, reg_param))
          
        return reg_param
      
      else:
        continue
