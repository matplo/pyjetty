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

    self.pt_bins = array('d', list(range(5, 301, 1)))

    # groomed jet splitting momentum fraction & angle
    self.obs_bins_zg = np.linspace(0, 0.5, 101)
    self.obs_bins_thetag = np.linspace(0, 1, 101)

    # ln(1/theta)
    self.obs_bins_lundx = np.concatenate((np.array([-1000, -100, -10, -1]), np.linspace(0, 2, 81), np.array([10, 100, 1000])))
    # ln(kT)
    self.obs_bins_lundy = np.concatenate((np.array([-1000, -100, -10]), np.linspace(-2, 2, 81), np.array([10, 100, 1000])))

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects(self):

    for jetR in self.jetR_list:

      for observable in self.observable_list:
        if observable != "lund":
          raise ValueError("Observable %s is not implemented in this script" % observable)

      for i in range(len(self.obs_settings[observable])):

        obs_setting = self.obs_settings[observable][i]
        grooming_setting = self.obs_grooming_settings[observable][i]
        obs_label = self.utils.obs_label(obs_setting, grooming_setting)

        if self.is_pp or self.include_no_subtraction:
          name = 'h_zg_thetag_JetPt_R%s_%s' % (jetR, obs_label)
          h = ROOT.TH3F(name, name, len(self.pt_bins)-1, self.pt_bins,
                        len(self.obs_bins_zg)-1, self.obs_bins_zg,
                        len(self.obs_bins_thetag)-1, self.obs_bins_thetag)
          h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
          h.GetYaxis().SetTitle('#it{z}_{g}')
          h.GetZaxis().SetTitle('#it{#theta}_{g}')
          setattr(self, name, h)

          name = 'h_lund_JetPt_R%s_%s' % (jetR, obs_label)
          h = ROOT.TH3F(name, name, len(self.pt_bins)-1, self.pt_bins,
                        len(self.obs_bins_lundx)-1, self.obs_bins_lundx,
                        len(self.obs_bins_lundy)-1, self.obs_bins_lundy)
          h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
          h.GetYaxis().SetTitle('ln(1/#it{#theta})')
          h.GetZaxis().SetTitle('ln(#it{k}_{T}/GeV)')
          setattr(self, name, h)

        if not self.is_pp:
          # Pb-Pb: have several R_max for contituent subtraction
          max_distance = self.max_distance if isinstance(self.max_distance, list) \
                         else self.max_distance[jetR]
          for R_max in max_distance:
            name = 'h_zg_thetag_JetPt_R%s_%s_Rmax%s' % (jetR, obs_label, R_max)
            h = ROOT.TH3F(name, name, len(self.pt_bins)-1, self.pt_bins,
                          len(self.obs_bins_zg)-1, self.obs_bins_zg,
                          len(self.obs_bins_thetag)-1, self.obs_bins_thetag)
            h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
            h.GetYaxis().SetTitle('#it{z}_{g}')
            h.GetZaxis().SetTitle('#it{#theta}_{g}')
            setattr(self, name, h)

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  #---------------------------------------------------------------
  def fill_jet_histograms(self, observable, jet, jet_groomed_lund, jetR, obs_setting,
                          grooming_setting, obs_label, jet_pt_ungroomed, suffix):

    # Calculate observable
    zg = jet_groomed_lund.z()
    thetag = jet_groomed_lund.Delta() / jetR
    # Values are already returned negative if jet_groomed_lund is empty
    #if not jet_groomed_lund.pair().has_constituents()

    # Fill histograms
    getattr(self, "h_zg_thetag_JetPt_R%s_%s%s" % (jetR, obs_label, suffix)).Fill(jet_pt_ungroomed, zg, thetag)

    # Fill jet lund planes
    name = "h_lund_JetPt_R%s_%s%s" % (jetR, obs_label, suffix)
    split = jet_groomed_lund
    if split.pair().has_constituents():
      lnth = np.log(jetR / split.Delta())
      lnkt = np.log(split.kt())
      getattr(self, name).Fill(jet_pt_ungroomed, lnth, lnkt)
    else:  # Untagged jet -- add underflow value
      getattr(self, name).Fill(jet_pt_ungroomed, -1e9, -1e9)


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
