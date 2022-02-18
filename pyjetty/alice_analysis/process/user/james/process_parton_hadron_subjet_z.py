#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of MC track information
  and do jet-finding, and save response histograms.
  
  Author: James Mulligan
          based on classes by Rey Cruz-Torres and Ezra Lesser
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
import fjcontrib

# Analysis utilities
from pyjetty.alice_analysis.process.user.substructure import process_parton_hadron_base

################################################################
class ProcessPH_subjet_z(process_parton_hadron_base.ProcessPHBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):

    # Initialize base class
    super(ProcessPH_subjet_z, self).__init__(
      input_file, config_file, output_dir, debug_level, **kwargs)

    self.obs_name = "#it{z}_{r}"

    # Define subjet finders (from first observable defined)
    self.subjet_def = {}
    for subjetR in self.obs_settings[self.observable_list[0]]:
      self.subjet_def[subjetR] = fj.JetDefinition(fj.antikt_algorithm, subjetR)

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects_R(self, jetR):

    for observable in self.observable_list:

      for i in range(len(self.obs_settings[observable])):

        obs_setting = self.obs_settings[observable][i] # r
        grooming_setting = self.obs_grooming_settings[observable][i]
        obs_label = self.utils.obs_label(obs_setting, grooming_setting)

        obs_bins = np.linspace(0, 1., 101)

        # Initialize 2D histograms for doing the MPI scaling at ch level
        self.init_MPI_scaling_hist(observable, self.obs_name, "ch", jetR, self.pt_bins, obs_bins, obs_label)

        for (level_1, level_2, MPI) in self.RM_levels:
          # Initialize the RMs for doing the folding
          self.init_response(observable, self.obs_name, level_1, level_2, MPI, jetR, self.pt_bins, obs_bins, obs_label)

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  # Fill 2D histogram of (pt, obs)
  #---------------------------------------------------------------
  def fill_observable_histograms(self, jet, jet_groomed_lund, jetR, level, MPI,
      obs_setting, grooming_setting, obs_label, jet_pt_ungroomed):

    # Only doing the MPI scaling at ch level, so ignore everything else
    if level != "ch":
      return

    # For a given jet, find all inclusive subjets of a given subjet radius
    cs_subjet = fj.ClusterSequence(jet.constituents(), self.subjet_def[obs_setting])
    subjets = fj.sorted_by_pt(cs_subjet.inclusive_jets())
    
    for observable in self.observable_list:

      # Fill inclusive subjets
      if 'inclusive' in observable:
        for subjet in subjets:
          z = self.z_r(subjet, jet_pt_ungroomed)
          self.fill_MPI_scaling_hist(observable, "ch", MPI, jetR, jet_pt_ungroomed, z, obs_label)
                    
      # Fill leading subjets
      elif 'leading' in observable:
        leading_subjet = self.leading_jet(subjets)
        z = self.z_r(leading_subjet, jet_pt_ungroomed)
        self.fill_MPI_scaling_hist(observable, "ch", MPI, jetR, jet_pt_ungroomed, z, obs_label)

      else:
        raise ValueError(f'fill_observable_histograms: {observable} is not a recognized observable')

  #---------------------------------------------------------------
  # Fill matched jet histograms
  #---------------------------------------------------------------
  def fill_matched_jet_histograms(self, jetR, obs_setting, grooming_setting, obs_label,
      jet_p, jet_p_groomed_lund, jet_h, jet_h_groomed_lund, jet_ch, jet_ch_groomed_lund,
      jet_pt_p_ungroomed, jet_pt_h_ungroomed, jet_pt_ch_ungroomed, suffix):

    # Only need to fill full<-->charged response
    MPI = 'off'

    # For a given jet, find all inclusive subjets of a given subjet radius
    cs_subjet_ch = fj.ClusterSequence(jet_ch.constituents(), self.subjet_def[obs_setting])
    subjets_ch = fj.sorted_by_pt(cs_subjet_ch.inclusive_jets())

    cs_subjet_h = fj.ClusterSequence(jet_h.constituents(), self.subjet_def[obs_setting])
    subjets_h = fj.sorted_by_pt(cs_subjet_h.inclusive_jets())

    for observable in self.observable_list:

      # Fill inclusive subjets
      if 'inclusive' in observable:
        for subjet_ch in subjets_ch:
            for subjet_h in subjets_h:
                if subjet_ch.delta_R(subjet_h) < 0.6*obs_setting:
                    z_ch = self.z_r(subjet_ch, jet_pt_ch_ungroomed)
                    z_h = self.z_r(subjet_h, jet_pt_h_ungroomed)
                    self.fill_response(observable, 'h', 'ch', MPI, jetR, jet_pt_h_ungroomed, jet_pt_ch_ungroomed,
                            z_h, z_ch, obs_label)

      # Fill leading subjets
      elif 'leading' in observable:
        leading_subjet_ch = self.leading_jet(subjets_ch)
        leading_subjet_h = self.leading_jet(subjets_h)
        z_ch = self.z_r(leading_subjet_ch, jet_pt_ch_ungroomed)
        z_h = self.z_r(leading_subjet_h, jet_pt_h_ungroomed)

        self.fill_response(observable, 'h', 'ch', MPI, jetR, jet_pt_h_ungroomed, jet_pt_ch_ungroomed,
                            z_h, z_ch, obs_label)

      else:
        raise ValueError(f'fill_matched_jet_histograms: {observable} is not a recognized observable')

  #---------------------------------------------------------------
  # Return leading jet (or subjet)
  #---------------------------------------------------------------
  def leading_jet(self, jets):

    leading_jet = None
    leading_jet_pt = 0.
    for jet in jets:
                      
      jet_pt = jet.pt()

      if not leading_jet:
        leading_jet = jet
        leading_jet_pt = jet_pt
      
      if jet_pt > leading_jet_pt:
        leading_jet = jet
        leading_jet_pt = jet_pt

    return leading_jet

  #---------------------------------------------------------------
  # Compute zr of subjet
  #---------------------------------------------------------------
  def z_r(self, subjet, jet_pt):
    z = subjet.pt() / jet_pt          
    if np.isclose(z, 1.):
        z = 0.999   # If z=1, it will be default be placed in overflow bin -- prevent this
    return z

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

  analysis = ProcessPH_subjet_z(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_mc()