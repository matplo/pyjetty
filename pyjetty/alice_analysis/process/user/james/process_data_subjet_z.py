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
class ProcessData_subjet_z(process_data_base.ProcessDataBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
  
    # Initialize base class
    super(ProcessData_subjet_z, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)

    # User-specific initialization
    self.initialize_user_config()

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_user_config(self):
  
    self.observable = self.observable_list[0]
    
    # Define subjet finders
    self.subjet_def = {}
    for subjetR in self.obs_settings[self.observable]:
      self.subjet_def[subjetR] = fj.JetDefinition(fj.antikt_algorithm, subjetR)
    
  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects(self):

    for jetR in self.jetR_list:
      for subjetR in self.obs_settings[self.observable]:
      
        name = 'h_{}_JetPt_R{}_{}'.format(self.observable, jetR, subjetR)
        h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 1.)
        h.GetXaxis().SetTitle('p_{T,ch jet}')
        h.GetYaxis().SetTitle('z_{subjet}')
        setattr(self, name, h)
      
  #---------------------------------------------------------------
  # Fill histograms
  #---------------------------------------------------------------
  def fill_jet_histograms(self, jet, jetR):
    
    # Find subjets
    result = [self.analyze_subjets(jet, jetR, subjetR) for subjetR in self.obs_settings[self.observable]]

  #---------------------------------------------------------------
  # For a given jet, find subjets of a given radius, and fill histograms
  #---------------------------------------------------------------
  def analyze_subjets(self, jet, jetR, subjetR):

    cs_subjet = fj.ClusterSequence(jet.constituents(), self.subjet_def[subjetR])
    subjets = fj.sorted_by_pt(cs_subjet.inclusive_jets())
    for subjet in subjets:
      z = subjet.pt() / jet.pt()
      getattr(self, 'h_{}_JetPt_R{}_{}'.format(self.observable, jetR, subjetR)).Fill(jet.pt(), z)
  
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

  analysis = ProcessData_subjet_z(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_data()
