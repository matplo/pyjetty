#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of track information
  and do jet-finding, and save basic histograms.

  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse

# Data analysis and plotting
import ROOT
import yaml

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib

# Base class
from pyjetty.alice_analysis.process.user.substructure import process_data_base

################################################################
class ProcessData_theta_g(process_data_base.ProcessDataBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):

    # Initialize base class
    super(ProcessData_theta_g, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects(self):

    for jetR in self.jetR_list:
      for observable in self.observable_list:

        # LaTeX-formatted observable name
        obs_name = self.obs_names[observable]

        # Get the maximum observable value for the histograms
        obs_max = None
        if observable == "theta_g":
          obs_max = 1.0
        elif observable == "zg":
          obs_max = 0.5
        else:
          # No other observables are implemented in this script
          raise ValueError("Observable %s not implemented" % observable)

        for grooming_setting in self.obs_grooming_settings[observable]:
          if grooming_setting:
            grooming_label = self.utils.grooming_label(grooming_setting)
            if self.is_pp:
              name = "h_%s_JetPt_R%s_%s" % (observable, jetR, grooming_label)
              h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 1.0)
              h.GetXaxis().SetTitle("#it{p}_{T}^{ch jet}")
              h.GetYaxis().SetTitle(obs_name)
              setattr(self, name, h)
            else:
              for R_max in self.max_distance:
                name = "h_%s_JetPt_R%s_%s_Rmax%s" % (observable, jetR, grooming_label, R_max)
                h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 1.0)
                h.GetXaxis().SetTitle("#it{p}_{T}^{ch jet}")
                h.GetYaxis().SetTitle(obs_name)
                setattr(self, name, h)


  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  #---------------------------------------------------------------
  def calculate_observable(self, observable, jet, jet_groomed_lund,
      jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed):

    if observable == "theta_g":
      return jet_groomed_lund.Delta() / jetR

    elif observable == "zg":
      return zg = jet_groomed_lund.z()

    # No other observables are implemented in this script
    raise ValueError("Observable %s not implemented" % observable)


  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  #---------------------------------------------------------------
  def fill_jet_histograms(self, observable, jet, jet_groomed_lund, jetR, obs_setting,
                          grooming_setting, obs_label, jet_pt_ungroomed, suffix):

    # Fill histograms
    getattr(self, "h_%s_JetPt_R%s_%s%s" % (observable, jetR, obs_label, suffix)).Fill(
      jet_pt_ungroomed, self.calculate_observable(observable, jet, jet_groomed_lund,
          jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed))


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
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
    sys.exit(0)

  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = ProcessData_theta_g(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_data()
