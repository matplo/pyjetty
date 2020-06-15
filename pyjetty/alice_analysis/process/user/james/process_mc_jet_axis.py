#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of MC track information
  and do jet-finding, and save response histograms.
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse

# Data analysis and plotting
import numpy as np
import ROOT
import yaml
from array import *

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjtools

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_io_emb
from pyjetty.alice_analysis.process.base import jet_info
from pyjetty.alice_analysis.process.user.substructure import process_mc_base
from pyjetty.alice_analysis.process.base import thermal_generator
from pyjetty.mputils import CEventSubtractor

################################################################
class ProcessMC_jet_axis(process_mc_base.ProcessMCBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
  
    # Initialize base class
    super(ProcessMC_jet_axis, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
    
    self.observable = self.observable_list[0]

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects_R(self, jetR):
      
    for i, axes in enumerate(self.obs_settings[self.observable]):
    
      grooming_setting = self.obs_grooming_settings[self.observable][i]
      if grooming_setting:
        grooming_label = self.utils.grooming_label(grooming_setting)
      else:
        grooming_label = ''
      
      name = 'hResidual_JetPt_{}_R{}_{}{}'.format(self.observable, jetR, axes, grooming_label)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, -1*jetR, jetR)
      h.GetXaxis().SetTitle('p_{T,truth}')
      h.GetYaxis().SetTitle('#frac{#DeltaR_{det}-#DeltaR_{truth}}{#DeltaR_{truth}}')
      setattr(self, name, h)

      # Create THn of response for jet axis deltaR
      dim = 4;
      title = ['p_{T,det}', 'p_{T,truth}', '#DeltaR_{det}', '#DeltaR_{truth}']
      nbins = [30, 30, 80, 40]
      min = [0., 0., 0., 0.]
      max = [150., 300., jetR, jetR]
      name = 'hResponse_JetPt_{}_R{}_{}{}'.format(self.observable, jetR, axes, grooming_label)
      self.create_thn(name, title, dim, nbins, min, max)

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  # Fill 2D histogram of (pt, obs)
  #---------------------------------------------------------------
  def fill_observable_histograms(self, hname, jet, jet_groomed_lund, jetR, obs_setting,
                                 grooming_setting, obs_label, jet_pt_ungroomed):
                                 
    # Recluster with WTA (with larger jet R)
    jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
    jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
    if self.debug_level > 3:
        print('WTA jet definition is:', jet_def_wta)
    reclusterer_wta =  fjcontrib.Recluster(jet_def_wta)
    jet_wta = reclusterer_wta.result(jet)

    # Compute jet axis differences
    jet_groomed = jet_groomed_lund.pair()
    if obs_setting == 'Standard_SD':
        deltaR = jet.delta_R(jet_groomed)
    elif obs_setting == 'Standard_WTA':
        deltaR = jet.delta_R(jet_wta)
    elif obs_setting == 'WTA_SD':
        deltaR = jet_groomed.delta_R(jet_wta)

    getattr(self,  hname.format(self.observable, jetR, obs_label)).Fill(jet_pt_ungroomed, deltaR)
      
  #---------------------------------------------------------------
  # Fill matched jet histograms
  #---------------------------------------------------------------
  def fill_matched_jet_histograms(self, jet_det, jet_det_groomed_lund, jet_truth,
                                  jet_truth_groomed_lund, jet_pp_det, jetR,
                                  obs_setting, grooming_setting, obs_label,
                                  jet_pt_det_ungroomed, jet_pt_truth_ungroomed, R_max, suffix):

    # Recluster with WTA (with larger jet R)
    jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
    jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
    if self.debug_level > 3:
        print('WTA jet definition is:', jet_def_wta)
    reclusterer_wta =  fjcontrib.Recluster(jet_def_wta)
    jet_det_wta = reclusterer_wta.result(jet_det)
    jet_truth_wta = reclusterer_wta.result(jet_truth)

    # Compute jet axis differences
    jet_det_groomed = jet_det_groomed_lund.pair()
    jet_truth_groomed = jet_truth_groomed_lund.pair()
    if obs_setting == 'Standard_SD':
        deltaR_det = jet_det.delta_R(jet_det_groomed)
        deltaR_truth = jet_truth.delta_R(jet_truth_groomed)
    elif obs_setting == 'Standard_WTA':
        deltaR_det = jet_det.delta_R(jet_det_wta)
        deltaR_truth = jet_truth.delta_R(jet_truth_wta)
    elif obs_setting == 'WTA_SD':
        deltaR_det = jet_det_groomed.delta_R(jet_det_wta)
        deltaR_truth = jet_truth_groomed.delta_R(jet_truth_wta)
      
    # Fill response
    self.fill_response(self.observable, jetR, jet_pt_det_ungroomed, jet_pt_truth_ungroomed,
                       deltaR_det, deltaR_truth, obs_label, R_max, prong_match = False)

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

  analysis = ProcessMC_jet_axis(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_mc()
