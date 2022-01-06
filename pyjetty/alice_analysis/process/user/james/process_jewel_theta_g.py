#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of track information
  and do jet-finding, and save basic histograms.
  This specific code is to run over jewel generator data to produce histograms (at truth level) that will then be compared to
  the data. That is, this base is to run over MC, but only at truth level, without response matrices.
  Author: James Mulligan (james.mulligan@berkeley.edu)
  Based on code by: Reynier Cruz Torres (reynier@lbl.gov) 
"""

from __future__ import print_function

# General
import os
import sys
import argparse

# Data analysis and plotting
import ROOT
import numpy as np

# Fastjet via python (from external library heppy)
import fastjet as fj

# Load pyjetty ROOT utils
ROOT.gSystem.Load('libpyjetty_rutil')

# Base class
from pyjetty.alice_analysis.process.user.substructure import process_jewel_generated_base

################################################################
class Process_CurvesFromJewelTracks_theta_g(process_jewel_generated_base.CurvesFromJewelTracks):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', **kwargs):
  
    # Initialize base class
    super(Process_CurvesFromJewelTracks_theta_g, self).__init__(input_file, config_file, output_dir, **kwargs)
    
  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects(self,label=''):
    
    for jetR in self.jetR_list:
        for observable in self.observable_list:
            
            if observable == 'zg':

                if not self.thermal_subtraction_method:
                    name = f'h_{observable}_JetPt_R{jetR}{label}'
                    h = ROOT.TH2F(name, name+';p_{T,ch jet};#it{z}_{g}', 200, 0, 200, 100, 0, 0.5)
                    setattr(self, name, h) 
                elif 'gridsub' in self.thermal_subtraction_method.lower():
                    for gridsize in self.gridsizes:
                        name = f'h_{observable}_JetPt_R{jetR}_gridsub_{gridsize}{label}'
                        h = ROOT.TH2F(name, name+';p_{T,ch jet};#it{z}_{g}', 200, 0, 200, 100, 0, 0.5)
                        setattr(self, name, h) 
                elif '4momsub' in self.thermal_subtraction_method.lower():
                    name = f'h_{observable}_JetPt_R{jetR}_4momsub{label}'
                    h = ROOT.TH2F(name, name+';p_{T,ch jet};#it{z}_{g}', 200, 0, 200, 100, 0, 0.5)
                    setattr(self, name, h) 
            
            elif observable == 'theta_g':

                if not self.thermal_subtraction_method:
                    name = f'h_{observable}_JetPt_R{jetR}{label}'
                    h = ROOT.TH2F(name, name+';p_{T,ch jet};#it{#theta}_{g}', 200, 0, 200, 100, 0, 1.)
                    setattr(self, name, h) 
                elif 'gridsub' in self.thermal_subtraction_method.lower():
                    for gridsize in self.gridsizes:
                        name = f'h_{observable}_JetPt_R{jetR}_gridsub_{gridsize}{label}'
                        h = ROOT.TH2F(name, name+';p_{T,ch jet};#it{#theta}_{g}', 200, 0, 200, 100, 0, 1.)
                        setattr(self, name, h) 
                elif '4momsub' in self.thermal_subtraction_method.lower():
                    name = f'h_{observable}_JetPt_R{jetR}_4momsub{label}'
                    h = ROOT.TH2F(name, name+';p_{T,ch jet};#it{#theta}_{g}', 200, 0, 200, 100, 0, 1.)
                    setattr(self, name, h) 

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  # But note that it is only called for the first observable -- in this case 
  #   it happens to work OK since the grooming settings are the same for zg/tg
  #---------------------------------------------------------------
  def fill_jet_histograms(self, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
                          obs_label, jet_pt_ungroomed, suffix=None, label=''):

    # Get groomed observables from Lund object
    theta_g = jet_groomed_lund.Delta() / jetR
    zg = jet_groomed_lund.z()

    for observable in self.observable_list:

        if not self.thermal_subtraction_method:
            name = f'h_{observable}_JetPt_R{jetR}{label}'
        elif 'gridsub' in self.thermal_subtraction_method.lower():
            name = f'h_{observable}_JetPt_R{jetR}_gridsub_{suffix}{label}'
        elif '4momsub' in self.thermal_subtraction_method.lower():
            name = f'h_{observable}_JetPt_R{jetR}_4momsub{label}'

        if observable == 'zg':
            getattr(self, name).Fill(jet_pt_ungroomed, zg)
        elif observable == 'theta_g':
            getattr(self, name).Fill(jet_pt_ungroomed, theta_g)

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

  analysis = Process_CurvesFromJewelTracks_theta_g(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_gen()