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

        if observable == 'theta_g':
        
          for grooming_setting in self.obs_grooming_settings[observable]:
            if grooming_setting:
              grooming_label = self.utils.grooming_label(grooming_setting)
              if self.is_pp:
                name = 'h_{}_JetPt_R{}_{}'.format(observable, jetR, grooming_label)
                h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 1.0)
                h.GetXaxis().SetTitle('p_{T,ch jet}')
                h.GetYaxis().SetTitle('#theta_{g,ch}')
                setattr(self, name, h)
              else:
                for R_max in self.max_distance:
                  name = 'h_{}_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, grooming_label, R_max)
                  h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 1.0)
                  h.GetXaxis().SetTitle('p_{T,ch jet}')
                  h.GetYaxis().SetTitle('#theta_{g,ch}')
                  setattr(self, name, h)
            
        if observable == 'zg':
        
          for grooming_setting in self.obs_grooming_settings[observable]:
            if grooming_setting:
              grooming_label = self.utils.grooming_label(grooming_setting)
              if self.is_pp:
                name = 'h_{}_JetPt_R{}_{}'.format(observable, jetR, grooming_label)
                h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 0.5)
                h.GetXaxis().SetTitle('p_{T,ch jet}')
                h.GetYaxis().SetTitle('z_{g,ch}')
                setattr(self, name, h)
              else:
                for R_max in self.max_distance:
                  name = 'h_{}_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, grooming_label, R_max)
                  h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0, 0.5)
                  h.GetXaxis().SetTitle('p_{T,ch jet}')
                  h.GetYaxis().SetTitle('z_{g,ch}')
                  setattr(self, name, h)

  #---------------------------------------------------------------
  # This function is called once for each groomed jet
  #---------------------------------------------------------------
  def fill_groomed_jet_histograms(self, grooming_setting, grooming_label,
                                  jet, jet_groomed, jet_pt_ungroomed, jetR, R_max):
                                      
    # Compute groomed observables
    if 'sd' in grooming_setting:

      sd_info = fjcontrib.get_SD_jet_info(jet_groomed)
      theta_g = sd_info.dR / jetR
      zg = sd_info.z

    elif 'dg' in grooming_setting:
        
      # (https://phab.hepforge.org/source/fastjetsvn/browse/contrib/contribs/LundPlane/tags/1.0.3/LundGenerator.hh)
      jet_dg_lund = jet_groomed
      theta_g = jet_dg_lund.Delta() / jetR
      zg = jet_dg_lund.z()
      
    # Fill histograms
    if R_max:
      suffix = '_Rmax{}'.format(R_max)
    else:
      suffix = ''
    if grooming_setting in self.obs_grooming_settings['theta_g']:
      getattr(self, 'h_theta_g_JetPt_R{}_{}{}'.format(jetR, grooming_label, suffix)).Fill(jet_pt_ungroomed, theta_g)
    if grooming_setting in self.obs_grooming_settings['zg']:
      getattr(self, 'h_zg_JetPt_R{}_{}{}'.format(jetR, grooming_label, suffix)).Fill(jet_pt_ungroomed, zg)

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
