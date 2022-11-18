#! /usr/bin/env python

import sys
import os
import argparse
from array import *
import numpy as np
import ROOT
import yaml

from pyjetty.alice_analysis.analysis.user.james import run_analysis_james_base
from pyjetty.alice_analysis.analysis.user.james import plotting_utils_subjet_z

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class RunAnalysisSubjetZ(run_analysis_james_base.RunAnalysisJamesBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(RunAnalysisSubjetZ, self).__init__(config_file, **kwargs)

    print(self)

  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #---------------------------------------------------------------
  def plot_single_result(self, jetR, obs_label, obs_setting, grooming_setting):
    print('Plotting each individual result...')
  
    # Plot final result for each 1D substructure distribution (with PYTHIA)
    self.plot_final_result(jetR, obs_label, obs_setting, grooming_setting)
    
  #----------------------------------------------------------------------
  def plot_final_result(self, jetR, obs_label, obs_setting, grooming_setting):
    print('Plot final results for {}: R = {}, {}'.format(self.observable, jetR, obs_label))

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()
    
    # Loop through pt slices, and plot final result for each 1D theta_g distribution
    for bin in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[bin]
      max_pt_truth = self.pt_bins_reported[bin+1]
      
      self.plot_observable(jetR, obs_label, obs_setting, grooming_setting, min_pt_truth, max_pt_truth, plot_pythia=True)
    
  #---------------------------------------------------------------
  # This function is called once after all subconfigurations have been looped over, for each R
  #---------------------------------------------------------------
  def plot_all_results(self, jetR):
    print('Plotting overlay of all results...')
    
    for i_config, overlay_list in enumerate(self.plot_overlay_list):
    
      if len(overlay_list) > 1:
      
        self.plot_final_result_overlay(i_config, jetR, overlay_list)
  
  #----------------------------------------------------------------------
  # This function is called once after all subconfigurations and jetR have been looped over
  #----------------------------------------------------------------------
  def plot_performance(self):
    
    if not self.do_plot_performance:
      return
    print('Plotting performance plots...')
    
    # Initialize performance plotting class, and plot
    if self.is_pp:
    
      self.plotting_utils = plotting_utils_subjet_z.PlottingUtils(self.output_dir_performance, self.config_file)
      self.plot_single_performance(self.output_dir_performance)
      
    else:
      
      # Plot for each R_max
      for R_max in self.max_distance:
      
        output_dir_performance = os.path.join(self.output_dir_performance, 'Rmax{}'.format(R_max))
        self.plotting_utils = plotting_utils_subjet_z.PlottingUtils(output_dir_performance, self.config_file, R_max = R_max)
        self.plot_single_performance(output_dir_performance, R_max)

        # Plot for thermal model
        if self.do_thermal_closure and R_max == self.R_max:
          
          output_dir_performance = os.path.join(self.output_dir_performance, 'thermal')
          self.plotting_utils = plotting_utils_subjet_z.PlottingUtils(output_dir_performance, self.config_file, R_max = R_max, thermal = True)
          self.plot_single_performance(output_dir_performance, R_max)

  #----------------------------------------------------------------------
  # This function is called once after all subconfigurations and jetR have been looped over
  #----------------------------------------------------------------------
  def plot_single_performance(self, output_dir_performance, R_max = None):
  
    if R_max:
      suffix = '_Rmax{}'.format(R_max)
    else:
      suffix = ''
      
    # Create output subdirectories
    self.create_output_subdir(output_dir_performance, 'jet')
    self.create_output_subdir(output_dir_performance, 'resolution')
    self.create_output_subdir(output_dir_performance, 'residual_pt')
    self.create_output_subdir(output_dir_performance, 'residual_obs')
    self.create_output_subdir(output_dir_performance, 'mc_projections_det')
    self.create_output_subdir(output_dir_performance, 'mc_projections_truth')
    self.create_output_subdir(output_dir_performance, 'truth')
    self.create_output_subdir(output_dir_performance, 'data')
    if 'leading' in self.observable:
        self.create_output_subdir(output_dir_performance, 'z1_crosscheck')
    if self.is_pp and 'inclusive' in self.observable:
        self.create_output_subdir(output_dir_performance, 'subjet_matching_pp')
    if not self.is_pp:
      self.create_output_subdir(output_dir_performance, 'delta_pt')
    
    # Generate performance plots
    for jetR in self.jetR_list:
  
      # Plot some subobservable-independent performance plots
      self.plotting_utils.plot_DeltaR(jetR, self.jet_matching_distance)
      self.plotting_utils.plot_JES(jetR)
      self.plotting_utils.plot_JES_proj(jetR, self.pt_bins_reported)
      self.plotting_utils.plotJER(jetR, self.utils.obs_label(self.obs_settings[0], self.grooming_settings[0]))
      self.plotting_utils.plot_jet_reco_efficiency(jetR, self.utils.obs_label(self.obs_settings[0], self.grooming_settings[0]))
      
      if not self.is_pp:
        self.plotting_utils.plot_delta_pt(jetR, self.pt_bins_reported)
      
      # Plot subobservable-dependent performance plots
      for i, _ in enumerate(self.obs_subconfig_list):

        obs_setting = self.obs_settings[i]
        grooming_setting = self.grooming_settings[i]
        obs_label = self.utils.obs_label(obs_setting, grooming_setting)
        
        if (jetR - obs_setting) < 1e-3:
          continue
    
        self.plotting_utils.plot_subjet_DeltaR(jetR, obs_label, self.jet_matching_distance)
        self.plotting_utils.plot_obs_resolution(jetR, obs_label, self.xtitle, self.pt_bins_reported)
        self.plotting_utils.plot_obs_residual_pt(jetR, obs_label, self.xtitle, self.pt_bins_reported)
        self.plotting_utils.plot_obs_residual_obs(jetR, obs_label, self.xtitle)
        self.plotting_utils.plot_obs_projections(jetR, obs_label, obs_setting, grooming_setting, self.xtitle, self.pt_bins_reported)
        self.plotting_utils.plot_obs_truth(jetR, obs_label, obs_setting, grooming_setting, self.xtitle, self.pt_bins_reported)
        if 'leading' in self.observable:
            self.plotting_utils.plot_z1_crosscheck(jetR, obs_label, obs_setting, grooming_setting, self.xtitle, self.pt_bins_reported)
        if self.is_pp and 'inclusive' in self.observable:
            self.plotting_utils.plot_subjet_matching_pp(jetR, obs_label, obs_setting, grooming_setting, self.xtitle, self.pt_bins_reported)
        
      if not self.is_pp:
      
          # Plot subjet matched pt histograms
          self.prong_match_threshold = 0.5
          min_pt = 80.
          max_pt = 100.
          for i, overlay_list in enumerate(self.plot_overlay_list):

            self.create_output_subdir(output_dir_performance, 'matched_pt_fraction_pt')
            hname = 'h_{}_matched_pt_JetPt_R{}'.format(self.observable, jetR)
            self.plotting_utils.plot_subjet_matching(i, jetR, hname, self.obs_subconfig_list, self.obs_settings, self.grooming_settings, overlay_list, self.prong_match_threshold, thermal=True)
            
            self.create_output_subdir(output_dir_performance, 'prong_matching_deltaR')
            self.create_output_subdir(output_dir_performance, 'prong_matching_deltaZ')
            
            name_prefix = f'h_{self.observable}_matched_pt_deltaZ_JetPt_R{jetR}'
            self.plotting_utils.plot_prong_matching_delta(i, jetR, name_prefix, self.obs_subconfig_list, self.obs_settings, self.grooming_settings, overlay_list, self.prong_match_threshold, min_pt, max_pt, plot_deltaz=True, plot_matched=True)
            self.plotting_utils.plot_prong_matching_delta(i, jetR, name_prefix, self.obs_subconfig_list, self.obs_settings, self.grooming_settings, overlay_list, self.prong_match_threshold, min_pt, max_pt, plot_deltaz=True, plot_matched=False)

            name_prefix = f'h_{self.observable}_matched_pt_deltaR_JetPt_R{jetR}'
            self.plotting_utils.plot_prong_matching_delta(i, jetR, name_prefix, self.obs_subconfig_list, self.obs_settings, self.grooming_settings, overlay_list, self.prong_match_threshold, min_pt, max_pt, plot_deltaz=False, plot_matched=True)
            self.plotting_utils.plot_prong_matching_delta(i, jetR, name_prefix, self.obs_subconfig_list, self.obs_settings, self.grooming_settings, overlay_list, self.prong_match_threshold, min_pt, max_pt, plot_deltaz=False, plot_matched=False)
            
          for i, _ in enumerate(self.obs_subconfig_list):

            obs_setting = self.obs_settings[i]
            obs_label = self.utils.obs_label(obs_setting, grooming_setting)
            
            if (jetR - obs_setting) < 1e-3:
              continue
            
            output_dir_money = os.path.join(output_dir_performance, 'matched_pt_money')
            self.create_output_subdir(output_dir_money, os.path.join(str(jetR), str(obs_setting)))
            self.plotting_utils.plot_subjet_money_plot(self.observable, jetR, R_max, self.prong_match_threshold,
                                                       obs_setting, self.pt_bins_reported,
                                                       output_dir_money, self.ytitle, thermal=False)

#----------------------------------------------------------------------
if __name__ == '__main__':

  # Define arguments
  parser = argparse.ArgumentParser(description='Jet substructure analysis')
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

  analysis = RunAnalysisSubjetZ(config_file = args.configFile)
  analysis.run_analysis()
