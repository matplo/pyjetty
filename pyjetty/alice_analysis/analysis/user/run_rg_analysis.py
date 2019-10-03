#! /usr/bin/env python

import sys
import os
import argparse
import itertools
from array import *
import numpy
import math
import ROOT
import yaml

import base
import analysis_utils
import roounfold_rg

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class run_rg_analysis(base.base):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(run_rg_analysis, self).__init__(**kwargs)
    self.config_file = config_file
    self.initialize_config()
    print(self)
  
    # Initialize utils class
    self.utils = analysis_utils.analysis_utils()

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    
    self.do_unfolding = config['do_unfolding']
    self.do_systematics = config['do_systematics']
    self.do_plot_final_result = config['do_plot_final_result']
    
    # config['beta'] is a dictionary of dictionaries, where each dict is for a value of beta
    self.jetR_list = config['jetR']
    beta_dict = config['beta']
    
    # Retrieve list of beta values
    self.beta_list = list(beta_dict.keys())
    
    # Retrieve histogram binnings for each beta value
    for beta in self.beta_list:
      
      binning_dict = beta_dict[beta]
      pt_bins_truth = (binning_dict['pt_bins_truth'])

      n_pt_bins_truth = len(pt_bins_truth) - 1
      setattr(self, 'n_pt_bins_truth_B{}'.format(beta), n_pt_bins_truth)

      truth_pt_bin_array = array('d',pt_bins_truth)
      setattr(self, 'truth_pt_bin_array_B{}'.format(beta), truth_pt_bin_array)
      
    self.reg_param_final = config['reg_param']
    self.min_pt_reported = 20
    self.max_pt_reported = 80
    self.regularizationParamName = 'n_iter'

    # List of possible systematic variations
    # Regularization parameter is automatically included
    self.kMain, self.kUnfoldingRange1, self.kUnfoldingRange2, self.kPrior1, self.kPrior2, self.kTrackEff = range(0, 6)

    # Set which systematics should be performed
    self.systematics_list = [self.kMain, self.kTrackEff]
      
    # Load paths to processing output, to be unfolded
    self.main_data = config['main_data']
    self.main_response = config['main_response']
    
    if self.kTrackEff in self.systematics_list:
      self.trkeff_data = config['trkeff_data']
      self.trkeff_response = config['trkeff_response']
          
    # Create output dirs
    self.file_format = config['file_format']
    self.output_dir = config['output_dir']
    if not self.output_dir.endswith("/"):
      self.output_dir = self.output_dir + "/"
      if not os.path.exists(self.output_dir):
        os.makedirs(self.output_dir)

    self.output_dir_main = os.path.join(self.output_dir, 'main')
    if not os.path.isdir(self.output_dir_main):
      os.makedirs(self.output_dir_main)
        
    if self.kTrackEff in self.systematics_list:
      self.output_dir_trkeff = os.path.join(self.output_dir, 'trkeff')
      if not os.path.isdir(self.output_dir_trkeff):
        os.makedirs(self.output_dir_trkeff)

  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def run_rg_analysis(self):

    if self.do_unfolding:
      self.perform_unfolding()
    
    for jetR in self.jetR_list:
      for beta in self.beta_list:

        if self.do_systematics:
          self.compute_systematics(jetR, beta)

        if self.do_plot_final_result:
          self.plot_final_result(jetR, beta)
      
  #----------------------------------------------------------------------
  def perform_unfolding(self):
    print('Perform unfolding for all systematic variations...')

    # Main result
    analysis_main = roounfold_rg.roounfold_rg(self.main_data, self.main_response, self.config_file, self.output_dir_main, self.file_format)
    analysis_main.roounfold_rg()

    # Tracking efficiency variation
    if self.kTrackEff in self.systematics_list:
      analysis_trkeff = roounfold_rg.roounfold_rg(self.trkeff_data, self.trkeff_response, self.config_file, output_dir_trkeff, self.file_format)
      analysis_trkeff.roounfold_rg()

  #----------------------------------------------------------------------
  def compute_systematics(self, jetR, beta):
    print('Compute systematics...')
  
    self.output_dir_systematics = os.path.join(self.output_dir, 'systematics')
    if not os.path.isdir(self.output_dir_systematics):
      os.makedirs(self.output_dir_systematics)

    # Get main result
    path_main = os.path.join(self.output_dir_main, 'fResult.root')
    fMain = ROOT.TFile(path_main, 'READ')
    
    name = 'hUnfolded_R{}_B{}_{}'.format(jetR, beta, self.reg_param_final)
    hMain = fMain.Get(name)
    hMain.SetDirectory(0)
    setattr(self, name, hMain)

    # Regularization parameter +2
    name = 'hUnfolded_R{}_B{}_{}'.format(jetR, beta, self.reg_param_final+2)
    hRegParam1 = fMain.Get(name)
    hRegParam1.SetDirectory(0)
    setattr(self, name, hRegParam1)
    
    # Regularization parameter -2
    name = 'hUnfolded_R{}_B{}_{}'.format(jetR, beta, self.reg_param_final-2)
    hRegParam2 = fMain.Get(name)
    hRegParam2.SetDirectory(0)
    setattr(self, name, hRegParam2)
    
    # Loop through pt slices, and compute systematics for each 1D theta_g distribution
    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_B{}'.format(beta))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_B{}'.format(beta))
    
    for bin in range(1, n_pt_bins_truth-1):
      min_pt_truth = truth_pt_bin_array[bin]
      max_pt_truth = truth_pt_bin_array[bin+1]
      
      self.compute_theta_systematic(jetR, beta, min_pt_truth, max_pt_truth)
    
  #----------------------------------------------------------------------
  def compute_theta_systematic(self, jetR, beta, min_pt_truth, max_pt_truth):
    
    #------------------------------------
    # Get 1D histograms
    # Normalize by integral, i.e. N_jets,inclusive in this pt-bin (cross-check this)

    # Get main histogram
    hMain2D = getattr(self, 'hUnfolded_R{}_B{}_{}'.format(jetR, beta, self.reg_param_final))
    hMain2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    hMain = hMain2D.ProjectionY()
    integral = hMain.Integral()
    hMain.Scale(1./integral, 'width')
    setattr(self, 'hMain_R{}_B{}_{}-{}'.format(jetR, beta, min_pt_truth, max_pt_truth), hMain)

    # Get reg param +2
    hRegParam1_2D = getattr(self, 'hUnfolded_R{}_B{}_{}'.format(jetR, beta, self.reg_param_final+2))
    hRegParam1_2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    hRegParam1 = hRegParam1_2D.ProjectionY()
    integral = hRegParam1.Integral()
    hRegParam1.Scale(1./integral, 'width')
    setattr(self, 'hRegParam1_R{}_B{}_{}-{}'.format(jetR, beta, min_pt_truth, max_pt_truth), hRegParam1)
    
    # Get reg param -2
    hRegParam2_2D = getattr(self, 'hUnfolded_R{}_B{}_{}'.format(jetR, beta, self.reg_param_final-2))
    hRegParam2_2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    hRegParam2 = hRegParam2_2D.ProjectionY()
    integral = hRegParam2.Integral()
    hRegParam2.Scale(1./integral, 'width')
    setattr(self, 'hRegParam2_R{}_B{}_{}-{}'.format(jetR, beta, min_pt_truth, max_pt_truth), hRegParam2)

    #------------------------------------
    # Compute systematics

    # Reg param +2
    name = 'hSystematic_RegParam1_R{}_B{}_{}-{}'.format(jetR, beta, min_pt_truth, max_pt_truth)
    hSystematic_RegParam1 = hMain.Clone()
    hSystematic_RegParam1.SetName(name)
    hSystematic_RegParam1.Divide(hRegParam1)
    self.change_to_per(hSystematic_RegParam1)
    setattr(self, name, hSystematic_RegParam1)
  
    outputFilename = os.path.join(self.output_dir_systematics, 'hSystematic_RegParam1_R{}_B{}_{}-{}{}'.format(jetR, beta, min_pt_truth, max_pt_truth, self.file_format))
    self.utils.plot_hist(hSystematic_RegParam1, outputFilename, 'P E')
  
    # Reg param -2
    name = 'hSystematic_RegParam2_R{}_B{}_{}-{}'.format(jetR, beta, min_pt_truth, max_pt_truth)
    hSystematic_RegParam2 = hMain.Clone()
    hSystematic_RegParam2.SetName(name)
    hSystematic_RegParam2.Divide(hRegParam2)
    self.change_to_per(hSystematic_RegParam2)
    setattr(self, name, hSystematic_RegParam2)
    
    outputFilename = os.path.join(self.output_dir_systematics, 'hSystematic_RegParam2_R{}_B{}_{}-{}{}'.format(jetR, beta, min_pt_truth, max_pt_truth, self.file_format))
    self.utils.plot_hist(hSystematic_RegParam2, outputFilename, 'P E')

  #----------------------------------------------------------------------
  def plot_final_result(self, jetR, beta):
    print('Plot final results...')

    output_dir_final = os.path.join(self.output_dir, 'final_results')
    if not os.path.isdir(output_dir_final):
      os.makedirs(output_dir_final)

  #----------------------------------------------------------------------
  def change_to_per(self, h):
    
    for bin in range(0, h.GetNbinsX()):
      content = h.GetBinContent(bin)
      content_new = math.fabs(1-content)
      h.SetBinContent(bin, content_new*100)

#----------------------------------------------------------------------
if __name__ == '__main__':

  # Define arguments
  parser = argparse.ArgumentParser(description='Unfold theta_g distribution')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='analysis_config.yaml',
                      help='Path of config file for analysis')

  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('configFile: \'{0}\''.format(args.configFile))
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = run_rg_analysis(config_file = args.configFile)
  analysis.run_rg_analysis()
