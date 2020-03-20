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
class AnalysisUtils_Obs(analysis_utils.AnalysisUtils):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(AnalysisUtils_Obs, self).__init__(**kwargs)

  #---------------------------------------------------------------
  # Get observable settings (i.e. list that stores the observable setting, e.g. subjetR)
  # from observable config block
  #---------------------------------------------------------------
  def obs_settings(self, observable, obs_config_dict, obs_subconfig_list):
  
    if observable == 'subjet_z':
      obs_settings = [obs_config_dict[name]['subjet_R'] for name in obs_subconfig_list]
    elif observable == 'jet_axis':
      obs_settings = [obs_config_dict[name]['axis'] for name in obs_subconfig_list]
    else:
      obs_settings = [None for _ in obs_subconfig_list]
      
    return obs_settings

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
  
  #---------------------------------------------------------------
  # Get formatted label from obs_setting and sd_setting
  #---------------------------------------------------------------
  def obs_label(self, obs_setting, sd_setting):

    obs_label = ''
    if obs_setting:
      obs_label += '{}'.format(obs_setting)
    if sd_setting:
      obs_label += '{}'.format(self.sd_label(sd_setting))
    return obs_label

  #---------------------------------------------------------------
  # Get formatted Soft Drop label from sd_setting = [zcut, beta]
  #---------------------------------------------------------------
  def sd_label(self, sd_setting):
  
      zcut = sd_setting[0]
      beta = sd_setting[1]
      sd_label = 'zcut{}_B{}'.format(self.remove_periods(zcut), beta)
      return sd_label
    
  #---------------------------------------------------------------
  # Get name of response THn
  #---------------------------------------------------------------
  def name_thn(self, observable, jetR, obs_label):
  
    if observable == 'jet_axis':
      return 'hResponse_JetPt_{}_{}_R{}Scaled'.format(observable, obs_label, jetR)
    else:
      return 'hResponse_JetPt_{}_R{}_{}Scaled'.format(observable, jetR, obs_label)

  #---------------------------------------------------------------
  # Get name of response THn, rebinned
  #---------------------------------------------------------------
  def name_thn_rebinned(self, observable, jetR, obs_label):
  
    if observable == 'jet_axis':
      return 'hResponse_JetPt_{}_{}_R{}_rebinned'.format(observable, obs_label, jetR)
    else:
      return 'hResponse_JetPt_{}_R{}_{}_rebinned'.format(observable, jetR, obs_label)
  
  #---------------------------------------------------------------
  # Get name of 2D data histogram
  #---------------------------------------------------------------
  def name_data(self, observable, jetR, obs_label):
  
    if observable == 'jet_axis':
      return 'h_{}_JetPt_{}_R{}'.format(observable, obs_label, jetR)
    else:
      return 'h_{}_JetPt_R{}_{}'.format(observable, jetR, obs_label)
  
  #---------------------------------------------------------------
  # Get name of 2D data histogram, rebinned
  #---------------------------------------------------------------
  def name_data_rebinned(self, observable, jetR, obs_label):
  
    if observable == 'jet_axis':
      return 'h_{}_JetPt_{}_R{}_rebinned'.format(observable, obs_label, jetR)
    else:
      return 'h_{}_JetPt_R{}_{}_rebinned'.format(observable, jetR, obs_label)

  #---------------------------------------------------------------
  # Get regularization parameter
  #---------------------------------------------------------------
  def get_reg_param(self, obs_settings, sd_settings, obs_subconfig_list, obs_config_dict, obs_label, observable, jetR):
    
    for i, _ in enumerate(obs_subconfig_list):
    
      obs_setting = obs_settings[i]
      sd_setting = sd_settings[i]
      
      if self.obs_label(obs_setting, sd_setting) == obs_label:
        
        config_name = obs_subconfig_list[i]
        reg_param = obs_config_dict[config_name]['reg_param'][jetR]
        #print('reg_param for {} {} jetR={}: {}'.format(obs_label, observable, jetR, reg_param))
          
        return reg_param
      
      else:
        continue
