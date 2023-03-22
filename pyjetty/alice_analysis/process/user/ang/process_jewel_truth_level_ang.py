#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of track information
  and do jet-finding, and save basic histograms.
  This specific code is to run over jewel generator data to produce histograms (at truth level) that will then be compared to
  the data. That is, this base is to run over MC, but only at truth level, without response matrices.

  Author: Ezra Lesser (elesser@berkeley.edu)
  Based on code by: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse

# Data analysis and plotting
import ROOT
import yaml
import numpy as np
from array import array

# Fastjet via python (from external library heppy)
import fjext

# Load pyjetty ROOT utils
ROOT.gSystem.Load('libpyjetty_rutil')

# Base class
from pyjetty.alice_analysis.process.user.substructure import process_jewel_generated_base

################################################################
class Process_CurvesFromJewelTracks_ang(process_jewel_generated_base.CurvesFromJewelTracks):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', **kwargs):

    # Initialize base class
    super(Process_CurvesFromJewelTracks_ang, self).__init__(input_file, config_file, output_dir, **kwargs)

    self.pt_bins = array('d', list(range(5, 305, 5)))
    self.obs_bins_ang = np.concatenate((np.linspace(0, 0.009, 10), np.linspace(0.01, 0.1, 19),
                                        np.linspace(0.11, 0.8, 70)))
    self.obs_bins_mass = np.concatenate(
      (np.linspace(0, 0.9, 10), np.linspace(1, 9.8, 45), np.linspace(10, 14.5, 10),
       np.linspace(15, 19, 5), np.linspace(20, 60, 9)))

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects(self, label=''):

    for jetR in self.jetR_list:
      for observable in self.observable_list:

        # Should only be two options: observable == "ang" or "mass"
        if observable != "ang" and observable != "mass":
          raise ValueError("Observable %s is not implemented in this script" % observable)

        obs_bins = getattr(self, "obs_bins_" + observable)
        for i, obs_setting in enumerate(self.obs_settings[observable]):

          grooming_setting = self.obs_grooming_settings[observable][i]
          obs_label = self.utils.obs_label(obs_setting, grooming_setting)

          # Histogram name based on recoil subtraction method
          names = []
          if not self.thermal_subtraction_method or \
              'negative_recombiner' in self.thermal_subtraction_method:
            names.append('h_%s_JetPt_R%s_%s%s' % (observable, jetR, obs_label, label) if \
              len(obs_label) else 'h_%s_JetPt_R%s%s' % (observable, jetR, label))
          elif 'gridsub' in self.thermal_subtraction_method:
            for gridsize in self.gridsizes:
              names.append('h_%s_JetPt_R%s_%s_gridsub_%s%s' % (observable, jetR, obs_label, gridsize, label) \
                if len(obs_label) else 'h_%s_JetPt_R%s_gridsub_%s%s' % (observable, jetR, gridsize, label))
          elif '4momsub' in self.thermal_subtraction_method:
            names.append('h_%s_JetPt_R%s_%s_4momsub%s' % (self.observable, jetR, obs_label, label) \
              if len(obs_label) else 'h_%s_JetPt_R%s_4momsub%s' % (self.observable, jetR, label))
          else:
            raise ValueError("Recoil subtraction method not recognized")

          for name in names:
            h = ROOT.TH2F(name, name, len(self.pt_bins) - 1, self.pt_bins, len(obs_bins) - 1, obs_bins)
            h.GetXaxis().SetTitle("#it{p}_{T}^{ch jet}")
            h.GetYaxis().SetTitle(self.obs_names[observable])
            h.Sumw2()
            setattr(self, name, h)


  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  #---------------------------------------------------------------
  def fill_jet_histograms(
      self, observable, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
      obs_label, jet_pt_ungroomed, suffix=None, label=''):

    check_user_index = False

    if not self.thermal_subtraction_method:
      name = 'h_%s_JetPt_R%s_%s%s' % (observable, jetR, obs_label, label) if \
        len(obs_label) else 'h_%s_JetPt_R%s%s' % (observable, jetR, label)
    elif 'negative_recombiner' in self.thermal_subtraction_method:
      name = 'h_%s_JetPt_R%s_%s%s' % (observable, jetR, obs_label, label) if \
        len(obs_label) else 'h_%s_JetPt_R%s%s' % (observable, jetR, label)
      check_user_index = True
    elif 'gridsub' in self.thermal_subtraction_method:
      name = 'h_%s_JetPt_R%s_%s_gridsub_%s%s' % (observable, jetR, obs_label, suffix, label) \
        if len(obs_label) else 'h_%s_JetPt_R%s_gridsub_%s%s' % (observable, jetR, suffix, label)
    elif '4momsub' in self.thermal_subtraction_method:
      name = 'h_%s_JetPt_R%s_%s_4momsub%s' % (self.observable, jetR, obs_label, label) \
        if len(obs_label) else 'h_%s_JetPt_R%s_4momsub%s' % (self.observable, jetR, label)
    else:
      raise ValueError("Recoil subtraction method not recognized")

    obs = None
    groomed_jet = None

    # Check to make sure that the jet is "real" with positive pT
    if grooming_setting:
      groomed_jet = jet_groomed_lund.pair()
      if groomed_jet.user_index() < 0:
        return
    else:  # no grooming
      if jet.user_index() < 0:
        return

    #######################################################################
    if observable == "ang":
      kappa = 1

      if grooming_setting:
        groomed_jet = jet_groomed_lund.pair()
        obs = fjext.lambda_beta_kappa(jet, groomed_jet, obs_setting, kappa, jetR, check_user_index)

      else:
        obs = fjext.lambda_beta_kappa(jet, obs_setting, kappa, jetR, check_user_index)

    #######################################################################
    elif observable == "mass":
      # m^2 = E^2 - p^2

      if grooming_setting:
        j_groomed = jet_groomed_lund.pair()
        if not j_groomed.has_constituents():
          # Untagged jet -- record underflow value
          obs = -1
        else:
          obs = j_groomed.m()

      else:
        obs = jet.m()


    # Fill histograms
    getattr(self, name).Fill(jet_pt_ungroomed, obs)

##################################################################
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Process Generator For Theory Comparison')
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
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
    sys.exit(0)

  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = Process_CurvesFromJewelTracks_ang(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_gen()
