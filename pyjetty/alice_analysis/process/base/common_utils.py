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
  # Get observable settings (i.e. list that stores the observable setting, e.g. subjetR)
  # from observable config block
  #---------------------------------------------------------------
  def obs_settings(self, observable, obs_config_dict, obs_subconfig_list):

    if 'subjet_z' in observable:
      return [obs_config_dict[name]['subjet_R'] for name in obs_subconfig_list]
    elif observable == 'jet_axis':
      return [obs_config_dict[name]['axis'] for name in obs_subconfig_list]
    elif observable == 'ang':
      return [obs_config_dict[name]['alpha'] for name in obs_subconfig_list]
    elif observable == 'theta_g':
      return [obs_config_dict[name]['SoftDrop'] for name in obs_subconfig_list if 'SoftDrop' in obs_config_dict[name]]

    # Else observable not implemented
    return [None for _ in obs_subconfig_list]
    
  #---------------------------------------------------------------
  #---------------------------------------------------------------
  #---------------------------------------------------------------
  # All observable-specific edits should be above here!
  #---------------------------------------------------------------
  #---------------------------------------------------------------
  #---------------------------------------------------------------

  #---------------------------------------------------------------
  # Get label from obs_setting and grooming_setting
  #---------------------------------------------------------------
  def obs_label(self, obs_setting, grooming_setting):

    obs_label = ''
    if obs_setting:
      obs_label += '{}'.format(obs_setting)
      if grooming_setting:
        obs_label += '_'
    if grooming_setting:
      obs_label += '{}'.format(self.grooming_label(grooming_setting))
    return obs_label
  
  #---------------------------------------------------------------
  # Get grooming settings (i.e. list that stores a dict of grooming
  # settings, i.e. {'sd': [zcut, beta]} or {'dg': [a]} from observable config block.
  # Note that one can also include multiple grooming settings
  # e.g. {'sd': [zcut, beta], 'dg': [a]}.
  #---------------------------------------------------------------
  def grooming_settings(self, obs_config_dict):
  
    grooming_settings = []
    for config_key, subconfig in obs_config_dict.items():
    
      grooming_setting_dict = {}
    
      if config_key == 'common_settings':
        continue
      if 'SoftDrop' in subconfig:
        grooming_config_dict = obs_config_dict[config_key]['SoftDrop']
        grooming_setting_dict['sd'] = [grooming_config_dict['zcut'], grooming_config_dict['beta']]
      if 'DynamicalGrooming' in subconfig:
        grooming_config_dict = obs_config_dict[config_key]['DynamicalGrooming']
        grooming_setting_dict['dg'] = [grooming_config_dict['a']]
      
      if grooming_setting_dict:
        grooming_settings.append(grooming_setting_dict)
      else:
        grooming_settings.append(None)
        
    return grooming_settings

  # Get formatted grooming label from grooming_setting,
  # i.e. {'sd': [zcut, beta]} or {'dg': [a]} or {'sd': [zcut, beta], 'dg': [a]}
  #---------------------------------------------------------------
  def grooming_label(self, grooming_setting):

    text = ''
    for key, value in grooming_setting.items():
    
      if key == 'sd':
        if text:
          text += '_'
        text += 'SD_zcut{}_B{}'.format(self.remove_periods(value[0]), value[1])
      if key == 'dg':
        if text:
           text += '_'
        if value[0] in  ['max_pt_soft', 'max_z', 'max_kt', 'max_kappa', 'max_tf', 'min_tf']:
          text += value[0]
        else:
          text += 'DG_a{}'.format(self.remove_periods(value[0]))
      
    if not text:
      sys.exit('CommonUtils::grooming_label: Unknown grooming type!')
    
    return text

  #---------------------------------------------------------------
  # Remove periods from a label
  #---------------------------------------------------------------
  def remove_periods(self, text):
    
    string = str(text)
    return string.replace('.', '')
