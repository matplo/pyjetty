#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of MC track information
  and do jet-finding, and save response histograms.
  
  Author: Ezra Lesser (elesser@berkeley.edu)
          based on substructure framework by James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse

# Data analysis and plotting
import numpy as np
from array import *

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjext

# Analysis utilities
from pyjetty.alice_analysis.process.user.substructure import process_parton_hadron_base

################################################################
class ProcessPH_ang(process_parton_hadron_base.ProcessPHBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
  
    # Initialize base class
    super(ProcessPH_ang, self).__init__(
      input_file, config_file, output_dir, debug_level, **kwargs)

    self.obs_name = "#it{#lambda}"

    # For angularities, easier to define here than type manually in config file
    self.obs_bins = np.concatenate((np.linspace(0, 0.0009, 10), np.linspace(0.001, 0.009, 9),
                                    np.linspace(0.01, 0.1, 19), np.linspace(0.11, 1., 90)))

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects_R(self, jetR):

    for observable in self.observable_list:

      for i in range(len(self.obs_settings[observable])):

        obs_setting = self.obs_settings[observable][i]
        grooming_setting = self.obs_grooming_settings[observable][i]
        obs_label = self.utils.obs_label(obs_setting, grooming_setting)

        # Initialize 2D histograms for doing the MPI scaling at ch level
        self.init_MPI_scaling_hist(
          observable, self.obs_name, "ch", jetR, self.pt_bins, self.obs_bins, obs_label)

        for (level_1, level_2, MPI) in self.RM_levels:
          # Initialize the RMs for doing the folding
          self.init_response(observable, self.obs_name, level_1, level_2, MPI, jetR,
                             self.pt_bins, self.obs_bins, obs_label)

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  # Fill 2D histogram of (pt, obs)
  #---------------------------------------------------------------
  def fill_observable_histograms(self, jet, jet_groomed_lund, jetR, level, MPI,
      obs_setting, grooming_setting, obs_label, jet_pt_ungroomed):

    if level != "ch":
      return

    # Calculate angularity and fill MPI scaling histogram
    ang = fjext.lambda_beta_kappa(jet, jet_groomed_lund.pair(), obs_setting, 1, jetR) \
          if grooming_setting else fjext.lambda_beta_kappa(jet, obs_setting, 1, jetR)
        
    self.fill_MPI_scaling_hist(
      self.observable_list[0], "ch", MPI, jetR, jet_pt_ungroomed, ang, obs_label)

  #---------------------------------------------------------------
  # Fill matched jet histograms
  #---------------------------------------------------------------
  def fill_matched_jet_histograms(self, jetR, obs_setting, grooming_setting, obs_label,
      jet_p, jet_p_groomed_lund, jet_h, jet_h_groomed_lund, jet_ch, jet_ch_groomed_lund,
      jet_pt_p_ungroomed, jet_pt_h_ungroomed, jet_pt_ch_ungroomed, suffix):

    # Compute observables
    obs = None
    if grooming_setting:
      obs = {
        "p" : fjext.lambda_beta_kappa(jet_p, jet_p_groomed_lund.pair(), obs_setting, 1, jetR),
        "h" : fjext.lambda_beta_kappa(jet_h, jet_h_groomed_lund.pair(), obs_setting, 1, jetR),
        "ch": fjext.lambda_beta_kappa(jet_ch, jet_ch_groomed_lund.pair(), obs_setting, 1, jetR) }
    else:
      obs = {
        "p" : fjext.lambda_beta_kappa(jet_p, obs_setting, 1, jetR),
        "h" : fjext.lambda_beta_kappa(jet_h, obs_setting, 1, jetR),
        "ch": fjext.lambda_beta_kappa(jet_ch, obs_setting, 1, jetR) }

    pt = { "p" : jet_pt_p_ungroomed,
           "h" : jet_pt_h_ungroomed,
           "ch": jet_pt_ch_ungroomed }

    # Fill histograms
    observable = self.observable_list[0]
    for (level_1, level_2, MPI) in self.RM_levels:
      self.fill_response(observable, level_1, level_2, MPI, jetR, pt[level_1], pt[level_2],
                         obs[level_1], obs[level_2], obs_label)


##################################################################
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Process MC')
  parser.add_argument('-f', '--inputFile', action='store',
                      type=str, metavar='inputFile',
                      default='AnalysisResults.root',
                      help='Path of ROOT file containing TTrees')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='config/analysis_config.yaml',
                      help="Path of config file for analysis")
  parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./TestOutput',
                      help='Output directory for output to be written to')
  
  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('inputFile: \'{0}\''.format(args.inputFile))
  print('configFile: \'{0}\''.format(args.configFile))
  print('ouputDir: \'{0}\"'.format(args.outputDir))

  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
    sys.exit(0)
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = ProcessPH_ang(
    input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_mc()
