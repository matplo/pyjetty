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
class ProcessData_jet_axis(process_data_base.ProcessDataBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
  
    # Initialize base class
    super(ProcessData_jet_axis, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
    
    # User-specific initialization
    self.initialize_user_config()

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_user_config(self):
  
    self.observable = self.observable_list[0]
    
  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects(self):
    
    for jetR in self.jetR_list:
    
      for i, axes in enumerate(self.obs_settings[self.observable]):
      
        grooming_setting = self.obs_grooming_settings[self.observable][i]
        if grooming_setting:
          grooming_label = self.utils.grooming_label(grooming_setting)
        else:
          grooming_label = ''
        
        name = 'h_{}_JetPt_R{}_{}{}'.format(self.observable, jetR, axes, grooming_label)
        h = ROOT.TH2F(name, name, 300, 0, 300, 200, 0, jetR)
        h.GetXaxis().SetTitle('p_{T,ch jet}')
        h.GetYaxis().SetTitle('#Delta R')
        setattr(self, name, h)
      
  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  #---------------------------------------------------------------
  def fill_jet_histograms(self, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
                          obs_label, jet_pt_ungroomed, suffix):
    
    # Recluster with WTA (with larger jet R)
    jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
    jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
    if self.debug_level > 2:
        print('WTA jet definition is:', jet_def_wta)
    reclusterer_wta =  fjcontrib.Recluster(jet_def_wta)
    jet_wta = reclusterer_wta.result(jet)

    jet_groomed = jet_groomed_lund.pair()
    self.fill_jet_axis_histogram(jet, jet_groomed, jet_wta, jetR, obs_setting,
                                 grooming_setting, obs_label)

  #---------------------------------------------------------------
  # Fill jet axis histograms
  #---------------------------------------------------------------
  def fill_jet_axis_histograms(self, jet, jet_groomed, jet_wta, jetR, obs_setting,
                               grooming_setting, obs_label):

    if obs_setting == 'Standard_SD':
      if grooming_setting in self.obs_grooming_settings[self.observable]:
        deltaR = jet.delta_R(jet_groomed)
        getattr(self, 'h_{}_JetPt_R{}_{}'.format(self.observable, jetR, obs_label)).Fill(jet.pt(), deltaR)
        
    if obs_setting == 'Standard_WTA':
      if grooming_setting == self.grooming_settings[0]:
        deltaR = jet.delta_R(jet_wta)
        getattr(self, 'h_{}_JetPt_R{}_{}'.format(self.observable, jetR, obs_label)).Fill(jet.pt(), deltaR)
      
    if obs_setting == 'WTA_SD':
      if grooming_setting in self.obs_grooming_settings[self.observable]:
        deltaR = jet_groomed.delta_R(jet_wta)
        getattr(self, 'h_{}_JetPt_R{}_{}'.format(self.observable, jetR, obs_label)).Fill(jet.pt(), deltaR)

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

  analysis = ProcessData_jet_axis(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_data()
