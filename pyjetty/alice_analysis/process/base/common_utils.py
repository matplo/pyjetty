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
from silx.io.dictdump import dicttoh5, h5todict

# Fastjet via python (from external library fjpydev)
try:
    import fastjet as fj
    import fjext
except ImportError:
    pass

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
      return [obs_config_dict[name]['beta'] for name in obs_subconfig_list]

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

  #---------------------------------------------------------------
  # Write nested dictionary of ndarray to hdf5 file
  # Note: all keys should be strings
  #---------------------------------------------------------------
  def write_data(self, results, output_dir, filename = 'results.h5'):
      print(f'Writing results to {output_dir}/{filename}...')
      dicttoh5(results, os.path.join(output_dir, filename), overwrite_data=True)
      print('done.')
      print()
  
  #---------------------------------------------------------------
  # Read dictionary of ndarrays from hdf5
  # Note: all keys should be strings
  #---------------------------------------------------------------
  def read_data(self, input_file):
      print(f'Loading results from {input_file}...')
      results = h5todict(input_file)
      print('done.')
      return results