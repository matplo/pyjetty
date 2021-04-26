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
import numpy as np

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
        
    # Define subjet finders (from first observable defined)
    self.subjet_def = {}
    for subjetR in self.obs_settings[self.observable_list[0]]:
      self.subjet_def[subjetR] = fj.JetDefinition(fj.antikt_algorithm, subjetR)
    
  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects(self):

    for jetR in self.jetR_list:
    
      for observable in self.observable_list:
      
        for subjetR in self.obs_settings[observable]:
        
          if (jetR - subjetR) < 1e-3:
            continue
        
          if self.is_pp:

              name = 'h_{}_JetPt_R{}_{}'.format(observable, jetR, subjetR)
              h = ROOT.TH2F(name, name, 200, 0, 200, 100, 0, 1.)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle('z_{r}')
              setattr(self, name, h)
              
              if 'leading' in observable:
                  name = 'h_{}_zconst_R{}_{}_z099_1'.format(observable, jetR, subjetR)
                  h = ROOT.TH2F(name, name, 20, 0, 200, 100, 0, 1.)
                  h.GetXaxis().SetTitle('p_{T,ch jet}')
                  h.GetYaxis().SetTitle('z')
                  setattr(self, name, h)
              
          else:
          
              for R_max in self.max_distance:
              
                  name = 'h_{}_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, subjetR, R_max)
                  h = ROOT.TH2F(name, name, 200, 0, 200, 100, 0, 1.)
                  h.GetXaxis().SetTitle('p_{T,ch jet}')
                  h.GetYaxis().SetTitle('z_{r}')
                  setattr(self, name, h)
                  
                  name = 'h_{}_zconst_R{}_{}_z099_1_Rmax{}'.format(observable, jetR, subjetR, R_max)
                  h = ROOT.TH2F(name, name, 20, 0, 200, 100, 0, 1.)
                  h.GetXaxis().SetTitle('p_{T,ch jet}')
                  h.GetYaxis().SetTitle('z')
                  setattr(self, name, h)
      
  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  #---------------------------------------------------------------
  def fill_jet_histograms(self, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
                          obs_label, jet_pt_ungroomed, suffix):
    
    if (jetR - obs_setting) < 1e-3:
      return
    
    # For a given jet, find all inclusive subjets of a given subjet radius
    cs_subjet = fj.ClusterSequence(jet.constituents(), self.subjet_def[obs_setting])
    subjets = fj.sorted_by_pt(cs_subjet.inclusive_jets())
    
    for observable in self.observable_list:

      # Fill inclusive subjets
      if 'inclusive' in observable:
        for subjet in subjets:
          z = subjet.pt() / jet.pt()
          
          # If z=1, it will be default be placed in overflow bin -- prevent this
          if np.isclose(z, 1.):
            z = 0.999
          
          getattr(self, 'h_{}_JetPt_R{}_{}{}'.format(observable, jetR, obs_setting, suffix)).Fill(jet.pt(), z)
          
      # Fill leading subjets
      if 'leading' in observable:
        leading_subjet = self.utils.leading_jet(subjets)
        z_leading = leading_subjet.pt() / jet.pt()
        
        # If z=1, it will be default be placed in overflow bin -- prevent this
        if np.isclose(z_leading, 1.):
            z_leading = 0.999
            
        getattr(self, 'h_{}_JetPt_R{}_{}{}'.format(observable, jetR, obs_setting, suffix)).Fill(jet.pt(), z_leading)
        
        # Fill z of subjet constituents for z~1 subjets
        if z_leading > 0.99:
            name = 'h_{}_zconst_R{}_{}_z099_1{}'.format(observable, jetR, obs_setting, suffix)
            for p in leading_subjet.constituents():
                z = p.pt() / leading_subjet.pt()
                if np.isclose(z, 1.):
                    z = 0.999
                getattr(self, name).Fill(jet.pt(), z)

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
