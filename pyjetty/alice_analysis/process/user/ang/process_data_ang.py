#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of track information
  and do jet-finding, and save basic histograms.

  Author: Ezra Lesser (elesser@berkeley.edu) with much code borrowed
          from original script by James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import argparse
import numpy as np
from array import array

# Data analysis and plotting
import ROOT

# Fastjet via python (from external library heppy)
import fjext

# Base class
from pyjetty.alice_analysis.process.user.substructure import process_data_base

################################################################
class ProcessData_ang(process_data_base.ProcessDataBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):

    # Initialize base class
    super(ProcessData_ang, self).__init__(
      input_file, config_file, output_dir, debug_level, **kwargs)

    self.pt_bins = array('d', list(range(5, 305, 5)))
    self.obs_bins_ang = np.concatenate((np.linspace(0, 0.009, 10), np.linspace(0.01, 0.1, 19),
                                        np.linspace(0.11, 0.8, 70)))
    self.obs_bins_mass = array('d', list(range(0, 61, 1)))

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects(self):

    for jetR in self.jetR_list:

      for observable in self.observable_list:
        # Should only be two: observable == "ang" or "mass"
        if observable != "ang" and observable != "mass":
          raise ValueError("Observable %s is not implemented in this script" % observable)

        obs_bins = getattr(self, "obs_bins_" + observable)
        for i in range(len(self.obs_settings[observable])):

          obs_setting = self.obs_settings[observable][i]
          grooming_setting = self.obs_grooming_settings[observable][i]
          obs_label = self.utils.obs_label(obs_setting, grooming_setting)

          if self.is_pp or self.include_no_subtraction:
            name = ('h_%s_JetPt_R%s_%s' % (observable, jetR, obs_label)) if \
              len(obs_label) else ('h_%s_JetPt_R%s' % (observable, jetR))
            h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                          len(obs_bins)-1, obs_bins)
            h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
            h.GetYaxis().SetTitle(self.obs_names[observable])
            setattr(self, name, h)

          if not self.is_pp:
            # Pb-Pb: have several R_max for contituent subtraction
            max_distance = self.max_distance if isinstance(self.max_distance, list) \
                           else self.max_distance[jetR]
            for R_max in max_distance:
              name = ('h_%s_JetPt_R%s_%s_Rmax%s' % (
                observable, jetR, obs_label, R_max)) if len(obs_label) else \
                ('h_%s_JetPt_R%s_Rmax%s' % (observable, jetR, R_max))
              h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                            len(obs_bins)-1, obs_bins)
              h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
              h.GetYaxis().SetTitle(self.obs_names[observable])
              setattr(self, name, h)

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  #---------------------------------------------------------------
  def fill_jet_histograms(self, observable, jet, jet_groomed_lund, jetR, obs_setting,
                          grooming_setting, obs_label, jet_pt_ungroomed, suffix):

    if observable == "ang":
      # Calculate angularity
      ang = fjext.lambda_beta_kappa(jet, jet_groomed_lund.pair(), obs_setting, 1, jetR) \
            if grooming_setting else fjext.lambda_beta_kappa(jet, obs_setting, 1, jetR)

      # Fill histograms
      getattr(self, "h_ang_JetPt_R%s_%s%s" % (jetR, obs_label, suffix)).Fill(
        jet_pt_ungroomed, ang)

    # Only do jet mass stuff once per set of angularity configs
    elif observable == "mass":
      name = 'h_mass_JetPt_R%s_%s%s' % (jetR, obs_label, suffix) if \
        grooming_setting else 'h_mass_JetPt_R%s%s' % (jetR, suffix)
      getattr(self, name).Fill(
        jet_pt_ungroomed, jet_groomed_lund.pair().m() if grooming_setting else jet.m())


##################################################################
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Process data')
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
  print('----------------------------------------------------------------')

  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFile):
    raise ValueError("File \"%s\" does not exist" % args.inputFile)

  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    raise ValueError("File \"%s\" does not exist" % args.configFile)

  analysis = ProcessData_ang(
    input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_data()
