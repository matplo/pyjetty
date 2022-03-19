#! /usr/bin/env python

'''
This class steers analysis of generic jet substructure analyses which require 2D unfolding.

To use this class, the following should be done:

  - Produce a histogram of the data, with name h_[obs]_JetPt_R[R]_[subobs]_[grooming setting]
    The grooming part is optional, and should be labeled e.g. zcut01_B0 — from CommonUtils::grooming_label({'sd':[zcut, beta]})
    For example: h_subjet_z_JetPt_R0.4_0.1
    For example: h_subjet_z_JetPt_R0.4_0.1_zcut01_B0

  - Produce a histogram of the response, with name hResponse_JetPt_[obs]_R[R]_[subobs]_[grooming setting]

  - Specify a configuration file, see examples in the config/ directory

  - Implement a user analysis class inheriting from this one, such as in james/run_analysis_theta_g.py
    You should implement the following functions:
      - plot_single_result()
      - plot_all_results()
      - plot_performance() (optional)
    See the example for details.

  - You also should modify a few observable-specific functions at the top of substructure/analysis_utils_obs.py
    and common_utils.py

Then, you should run your analysis with:

  python run_analysis_user.py /path/to/my/config.yaml


Author: James Mulligan (james.mulligan@berkeley.edu)
'''

import sys
import os
import argparse
import itertools
from array import *
import numpy as np
import math
import ROOT
import yaml
import shutil
import hepdata_lib
import pickle

from pyjetty.alice_analysis.analysis.base import common_base
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs
from pyjetty.alice_analysis.analysis.user.substructure import roounfold_obs

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class RunAnalysis(common_base.CommonBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(RunAnalysis, self).__init__(**kwargs)
    self.config_file = config_file

    # Initialize utils class
    self.utils = analysis_utils_obs.AnalysisUtils_Obs()

    # Initialize yaml config
    self.initialize_config()
    
    self.ColorArray = [ROOT.kBlue-4, ROOT.kAzure+7, ROOT.kCyan-2, ROOT.kViolet-8,
                       ROOT.kBlue-6, ROOT.kGreen+3, ROOT.kPink-4, ROOT.kRed-4,
                       ROOT.kOrange-3]
    self.MarkerArray = [20, 21, 22, 23, 33, 34, 24, 25, 26, 32]

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):

    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)

    # Set list of observables
    self.observable = config['analysis_observable']
    self.debug_level = config['debug_level']

    # Set which analysis steps to perform
    self.do_unfolding = config['do_unfolding']
    self.do_systematics = config['do_systematics']
    self.do_plot_final_result = config['do_plot_final_result']
    self.force_rebin_response=config['force_rebin']
    self.do_plot_performance = config['do_plot_performance']
    
    # Flag for whether to do kinematic efficiency correction in RM, or by hand
    if 'use_miss_fake' in config:
        self.use_miss_fake = config['use_miss_fake']
    else:
        self.use_miss_fake = True
    
    # Set whether pp or PbPb
    if 'constituent_subtractor' in config:
        self.is_pp = False
        self.R_max = config['constituent_subtractor']['main_R_max']
    else:
        self.is_pp = True
        self.R_max = None

    # Get the sub-configs to unfold
    self.jetR_list = config['jetR']
    self.obs_config_dict = config[self.observable]
    self.obs_subconfig_list = [name for name in list(self.obs_config_dict.keys()) if 'config' in name ]
    self.grooming_settings = self.utils.grooming_settings(self.obs_config_dict)
    self.obs_settings = self.utils.obs_settings(self.observable, self.obs_config_dict, self.obs_subconfig_list)
    self.obs_labels = [self.utils.obs_label(self.obs_settings[i], self.grooming_settings[i])
                       for i in range(len(self.obs_subconfig_list))]
    self.xtitle = self.obs_config_dict['common_settings']['xtitle']
    self.ytitle = self.obs_config_dict['common_settings']['ytitle']
    self.pt_bins_reported = self.obs_config_dict['common_settings']['pt_bins_reported']

    self.use_max_reg_param = False
    if 'max_reg_param' in self.obs_config_dict['common_settings']:
      self.max_reg_param = self.obs_config_dict['common_settings']['max_reg_param']
      if self.max_reg_param < 3:
        print("ERROR: The minimum number of iterations has been set to 3.",
              "Please set max_reg_param to a value >= 3.")
        raise ValueError()
      self.use_max_reg_param = True
    self.reg_param_name = 'n_iter'

    # Retrieve histogram binnings for each observable setting
    for i, _ in enumerate(self.obs_subconfig_list):

      config_name = self.obs_subconfig_list[i]
      obs_label = self.obs_labels[i]

      pt_bins_truth = (self.obs_config_dict[config_name]['pt_bins_truth'])
      n_pt_bins_truth = len(pt_bins_truth) - 1
      setattr(self, 'n_pt_bins_truth_{}'.format(obs_label), n_pt_bins_truth)
      truth_pt_bin_array = array('d',pt_bins_truth)
      setattr(self, 'truth_pt_bin_array_{}'.format(obs_label), truth_pt_bin_array)

      obs_bins_truth = (self.obs_config_dict[config_name]['obs_bins_truth'])
      n_obs_bins_truth = len(obs_bins_truth) - 1
      setattr(self, 'n_obs_bins_truth_{}'.format(obs_label), n_obs_bins_truth)
      truth_obs_bin_array = array('d',obs_bins_truth)
      setattr(self, 'truth_obs_bin_array_{}'.format(obs_label), truth_obs_bin_array)
      
      pt_bins_det = (self.obs_config_dict[config_name]['pt_bins_det'])
      det_pt_bin_array = array('d',pt_bins_det)
      setattr(self, 'det_pt_bin_array_{}'.format(obs_label), det_pt_bin_array)

    # List of systematic variations to perform
    self.systematics_list = config['systematics_list']

    # Load paths to processing output, to be unfolded
    self.main_data = config['main_data']
    self.main_response = config['main_response']

    if 'trkeff' in self.systematics_list:
      self.trkeff_response = config['trkeff_response']
    if 'delta_matching' in self.systematics_list:
      self.delta_matching_response = config['delta_matching_response']
    if 'fastsim_generator0' in self.systematics_list:
      self.fastsim_response_list = config['fastsim_response']
    if 'random_mass' in self.systematics_list:
      self.randmass_data = config['randmass_data']
      self.randmass_response = config['randmass_response']

    self.prior1_variation_parameter = config['prior1_variation_parameter']
    self.prior2_variation_parameter = config['prior2_variation_parameter']
      
    if 'subtraction1' in self.systematics_list:
      self.R_max1 = config['R_max_variation1']
    if 'subtraction2' in self.systematics_list:
      self.R_max2 = config['R_max_variation2']
      
    if 'thermal_closure' in config:
      self.do_thermal_closure = True
      self.fThermal = config['thermal_closure']
    else:
      self.do_thermal_closure = False

    # Create output dirs
    self.file_format = config['file_format']
    self.output_dir = config['output_dir']
    self.create_output_dirs()

  #---------------------------------------------------------------
  # Create a set of output directories for a given observable
  #---------------------------------------------------------------
  def create_output_dirs(self):

    output_dir = os.path.join(self.output_dir, self.observable)
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

    if self.do_plot_performance:
      self.create_output_subdir(output_dir, 'performance')

    for systematic in self.systematics_list:
      self.create_output_subdir(output_dir, systematic)
      
    if self.do_thermal_closure:
      self.create_output_subdir(output_dir, 'thermal_closure')

    if self.do_systematics:
      output_dir_systematics = self.create_output_subdir(output_dir, 'systematics')
      sys_root_filename = os.path.join(output_dir_systematics, 'fSystematics.root')
      fSystematics = ROOT.TFile(sys_root_filename, 'RECREATE')
      
      self.logfile = os.path.join(output_dir_systematics, 'log.txt')
      open(self.logfile, 'w').close()
      
      self.create_output_subdir(output_dir, 'unfolding_tests')

    if self.do_plot_final_result:
      self.output_dir_systematics = os.path.join(output_dir, 'systematics')
      output_dir_final = self.create_output_subdir(output_dir, 'final_results')
      final_result_root_file = os.path.join(output_dir_final, 'fFinalResults.root')
      if self.do_systematics:
        fFinalResults = ROOT.TFile(final_result_root_file, 'RECREATE')

  #---------------------------------------------------------------
  # Create a single output subdirectory
  #---------------------------------------------------------------
  def create_output_subdir(self, output_dir, name):

    output_subdir = os.path.join(output_dir, name)
    setattr(self, 'output_dir_{}'.format(name), output_subdir)
    if not os.path.isdir(output_subdir):
      os.makedirs(output_subdir)

    return output_subdir

  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def run_analysis(self):

    if self.do_unfolding:
      self.perform_unfolding()

    # Loop through jet radii
    for jetR in self.jetR_list:

      # Loop through subconfigurations to unfold
      for i, subconfig in enumerate(self.obs_subconfig_list):

        obs_setting = self.obs_settings[i]
        grooming_setting = self.grooming_settings[i]
        obs_label = self.obs_labels[i]
        self.get_obs_min_bins(self.obs_config_dict[subconfig], obs_label)
        self.get_obs_max_bins(self.obs_config_dict[subconfig], obs_label)

        # Compute systematics and attach to main results
        if self.do_systematics:
            self.compute_systematics(jetR, obs_label, obs_setting, grooming_setting)

        # Plot result for each individial subconfiguration
        if self.do_plot_final_result and self.do_systematics:
          # You must implement this
          self.plot_single_result(jetR, obs_label, obs_setting, grooming_setting)
          
      # Plot final results for all subconfigs (may need to comment this out if only writing hepdata)
      if self.do_plot_final_result and self.do_systematics:
        self.plot_all_results(jetR) # You must implement this
        
    # Write hepdata submission
    if self.do_plot_final_result:
        self.write_hepdata()

    # Plot additional performance plots
    if self.do_plot_performance:
      self.plot_performance() # You must implement this
        
  #----------------------------------------------------------------------
  def perform_unfolding(self):
    print('Perform unfolding for all systematic variations: {} ...'.format(self.observable))

    for systematic in self.systematics_list:

      output_dir = getattr(self, 'output_dir_{}'.format(systematic))
      data = self.main_data
      response = self.main_response

      main_response_location = os.path.join(getattr(self, 'output_dir_main'), 'response.root')
      rebin_response = self.check_rebin_response(output_dir)

      prior_variation_parameter = 0.
      truncation = False
      binning = False
      R_max =  self.R_max
      prong_matching_response = False

      if systematic == 'trkeff':
        response = self.trkeff_response
      elif systematic == 'prior1':
        prior_variation_parameter = self.prior1_variation_parameter
      elif systematic == 'prior2':
        prior_variation_parameter = self.prior2_variation_parameter
      elif systematic == 'truncation':
        truncation = True
      elif systematic == 'binning':
        binning = True
      elif systematic == 'subtraction1':
        R_max = self.R_max1
      elif systematic == 'subtraction2':
        R_max = self.R_max2
      elif systematic == 'prong_matching':
        prong_matching_response = True
      elif systematic == 'fastsim_generator0':
        response = self.fastsim_response_list[0]
      elif systematic == 'fastsim_generator1':
        response = self.fastsim_response_list[1]
      elif systematic == 'random_mass':
        data = self.randmass_data
        response = self.randmass_response
      elif systematic == 'delta_matching':
        response = self.delta_matching_response

      analysis = roounfold_obs.Roounfold_Obs(
        self.observable, data, response, self.config_file, output_dir, self.file_format,
        rebin_response=rebin_response, prior_variation_parameter=prior_variation_parameter,
        truncation=truncation, binning=binning, R_max=R_max,
        prong_matching_response=prong_matching_response, use_miss_fake=self.use_miss_fake)
      analysis.roounfold_obs()
      
    # Unfold thermal closure test
    if self.do_thermal_closure:
    
      output_dir = getattr(self, 'output_dir_thermal_closure')
      rebin_response = self.check_rebin_response(output_dir)
      
      analysis = roounfold_obs.Roounfold_Obs(
        self.observable, self.fThermal, self.fThermal, self.config_file, output_dir,
        self.file_format, rebin_response=rebin_response, R_max=R_max, thermal_model = True,
        use_miss_fake=self.use_miss_fake)
      analysis.roounfold_obs()

  #----------------------------------------------------------------------
  def check_rebin_response(self, output_dir):

    rebin_response = True
    response_path = os.path.join(output_dir, 'response.root')
    if os.path.exists(response_path):
      rebin_response = False
      if self.force_rebin_response:
        print('Response {} exists -- force re-create...'.format(response_path))
        rebin_response = True
      else:
        print('Response {} exists -- don\'t re-create.'.format(response_path))
    else:
      print('Response {} doesn\'t exist -- create it...'.format(response_path))

    return rebin_response

  #----------------------------------------------------------------------
  def get_obs_max_bins(self, subconfig, obs_label):

    if "obs_max_reported" in subconfig:
      obs_bins_truth = self.truth_bin_array(obs_label)
      obs_max_reported = subconfig["obs_max_reported"]
      setattr(self, "obs_max_reported_{}".format(obs_label), obs_max_reported)

      # Idiot checks
      if len(obs_max_reported) != len(self.pt_bins_reported)-1:
        raise ValueError(
          "Length of obs_max_reported {} should match number of pt_bins_reported {}".format(
          obs_max_reported, self.pt_bins_reported))

      for val in obs_max_reported:
        if not val in obs_bins_truth:
          raise ValueError("Final bin cutoff {} is not in the bin edges list {}".format(
            val, obs_bins_truth))

      obs_max_bins = [obs_bins_truth.index(i) for i in obs_max_reported]

    else:
      obs_max_bins = [None for i in self.pt_bins_reported]

    setattr(self, f"obs_max_bins_{obs_label}", obs_max_bins)

  #----------------------------------------------------------------------
  def get_obs_min_bins(self, subconfig, obs_label):

    if "obs_min_reported" in subconfig:
      obs_bins_truth = self.truth_bin_array(obs_label)
      obs_min_reported = subconfig["obs_min_reported"]
      setattr(self, "obs_min_reported_{}".format(obs_label), obs_min_reported)

      # Idiot checks
      if len(obs_min_reported) != len(self.pt_bins_reported)-1:
        raise ValueError(
          "Length of obs_min_reported {} should match number of pt_bins_reported {}".format(
          obs_min_reported, self.pt_bins_reported))

      for val in obs_min_reported:
        if not val in obs_bins_truth:
          raise ValueError("Final bin cutoff {} is not in the bin edges list {}".format(
            val, obs_bins_truth))

      obs_min_bins = [obs_bins_truth.index(i) for i in obs_min_reported]

    else:
      obs_min_bins = [None for i in self.pt_bins_reported]

    setattr(self, f"obs_min_bins_{obs_label}", obs_min_bins)

  #----------------------------------------------------------------------
  def compute_systematics(self, jetR, obs_label, obs_setting, grooming_setting):

    text = 'Compute systematics for {}: R = {}, subobs = {}, grooming = {}'.format(
            self.observable, jetR, obs_label, grooming_setting)
    with open(self.logfile, 'a') as myfile:
      myfile.write(text + '\n')
    print(text)

    # Select final regularization parameter
    try:
      reg_param_final = self.utils.get_reg_param(
        self.obs_settings, self.grooming_settings, self.obs_subconfig_list,
        self.obs_config_dict, obs_label, jetR)
      min_reg_param = reg_param_final
    except:
      if self.use_max_reg_param:
        reg_param_final = self.max_reg_param
        min_reg_param = 3
      else:
        print("ERROR: No max_reg_param or reg_param set in config file for", obs_label)
        exit(1)

    # First, calculate systematics for all possible reg params
    for reg_param in range(min_reg_param, reg_param_final + 1):

      # Load 2D unfolded results from their files into attributes
      self.load_2D_observables(jetR, obs_label, obs_setting, grooming_setting, reg_param)

      # Loop through pt slices, and compute systematics for each 1D observable distribution
      for i in range(0, len(self.pt_bins_reported) - 1):
        min_pt_truth = self.pt_bins_reported[i]
        max_pt_truth = self.pt_bins_reported[i+1]
        minbin = self.obs_min_bins(obs_label)[i]
        maxbin = self.obs_max_bins(obs_label)[i]

        # Load 1D unfolded results for each pt slice into attributes
        self.load_1D_observables(jetR, obs_label, obs_setting, grooming_setting, reg_param,
                                 min_pt_truth, max_pt_truth, minbin, maxbin)

        # Compute systematics of the 1D distributions for each pt slice
        self.compute_obs_systematic(jetR, obs_label, obs_setting, grooming_setting, reg_param,
                                    min_pt_truth, max_pt_truth, minbin, maxbin, final=False)

    # Now determine the best possible reg param by minimizing uncertainty
    for i in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[i]
      max_pt_truth = self.pt_bins_reported[i+1]
      minbin = self.obs_min_bins(obs_label)[i]
      maxbin = self.obs_max_bins(obs_label)[i]

      try:
        reg_param_final = self.utils.get_reg_param(
          self.obs_settings, self.grooming_settings, self.obs_subconfig_list,
          self.obs_config_dict, obs_label, jetR)
      except:
        if self.use_max_reg_param:
          reg_param_final = self.determine_reg_param_final(
            jetR, obs_label, obs_setting, grooming_setting,
            min_pt_truth, max_pt_truth, minbin, maxbin)
          text = 'Optimal regularization parameter for pT={}-{} determined to be {} = {}.'.format(
            min_pt_truth, max_pt_truth, self.reg_param_name, reg_param_final)
          with open(self.logfile, 'a') as myfile:
            myfile.write(text + '\n')
          print(text)
        else:
          print("ERROR: No max_reg_param or reg_param set in config file for", obs_label)
          exit(1)

      self.compute_obs_systematic(jetR, obs_label, obs_setting, grooming_setting, reg_param_final,
                                  min_pt_truth, max_pt_truth, minbin, maxbin, final=True)
                                  
      # Set SD tagging fraction for final reg parameter
      if grooming_setting and 'sd' in grooming_setting:
        f_tagging_name = 'tagging_fraction_R{}_{}_{}-{}'.format(
          jetR, obs_label, min_pt_truth, max_pt_truth)
        f_tagged = getattr(self, '{}_{}'.format(f_tagging_name, reg_param_final))
        setattr(self, f_tagging_name, f_tagged)
        
      # Copy plots of final reg param, for convenience
      self.copy_unfolding_tests(jetR, obs_label, reg_param_final, min_pt_truth, max_pt_truth)
      
  #----------------------------------------------------------------------
  def copy_unfolding_tests(self, jetR, obs_label, reg_param_final, min_pt, max_pt):

      outputdir_main = getattr(self, 'output_dir_main')
      label = 'R{}_{}'.format(self.utils.remove_periods(jetR), obs_label)

      # Make a folder for each obs_label and pt bin
      outputdir = os.path.join(getattr(self, 'output_dir_unfolding_tests'), label)
      outputdir = os.path.join(outputdir, '{}-{}'.format(min_pt, max_pt))
      if not os.path.exists(outputdir):
        os.makedirs(outputdir)

      # Copy refolding test
      outputdir_test = os.path.join(outputdir_main, 'Test_Refolding')
      old_name = 'hFoldedTruth_pt_{}_{}{}'.format(label, reg_param_final, self.file_format)
      new_name = 'hFoldedTruth_pt_{}_final{}'.format(label, self.file_format)
      shutil.copy(os.path.join(outputdir_test, old_name), os.path.join(outputdir, new_name))
      
      pt_det_bins = getattr(self, 'det_pt_bin_array_{}'.format(obs_label))
      if min_pt < pt_det_bins[0]:
        pt_det_bins = pt_det_bins
      else:
        pt_det_bins = pt_det_bins[pt_det_bins.index(min_pt) : pt_det_bins.index(max_pt) + 1]
      for i in range(len(pt_det_bins)-1):
        min_det_pt = int(pt_det_bins[i])
        max_det_pt = int(pt_det_bins[i+1])
        old_name = 'hFoldedTruth_{}_{}-{}_{}{}'.format(
          label, min_det_pt, max_det_pt, reg_param_final, self.file_format)
        new_name = 'hFoldedTruth_{}_{}-{}_final{}'.format(
          label, min_det_pt, max_det_pt, self.file_format)
        shutil.copy(os.path.join(outputdir_test, old_name), os.path.join(outputdir, new_name))
      
      # Copy statistical closure test
      outputdir_test = os.path.join(outputdir_main, 'Test_StatisticalClosure')
      old_name = 'hClosure_pt_{}_{}{}'.format(label, reg_param_final, self.file_format)
      new_name = 'hStatisticalClosure_pt_{}_final{}'.format(label, self.file_format)
      shutil.copy(os.path.join(outputdir_test, old_name), os.path.join(outputdir, new_name))
      old_name = 'hClosure_{}_{}-{}_{}{}'.format(
        label, min_pt, max_pt, reg_param_final, self.file_format)
      new_name = 'hStatisticalClosure_{}_{}-{}_final{}'.format(
        label, min_pt, max_pt, self.file_format)
      shutil.copy(os.path.join(outputdir_test, old_name), os.path.join(outputdir, new_name))
      
      # Copy shape closure test
      parameter = self.utils.remove_periods(self.prior1_variation_parameter)
      outputdir_test = os.path.join(outputdir_main, 'Test_ShapeClosure{}'.format(parameter))
      old_name = 'hClosure_pt_{}_{}{}'.format(label, reg_param_final, self.file_format)
      new_name = 'hShapeClosure{}_pt_{}_final{}'.format(parameter, label, self.file_format)
      shutil.copy(os.path.join(outputdir_test, old_name), os.path.join(outputdir, new_name))
      old_name = 'hClosure_{}_{}-{}_{}{}'.format(
        label, min_pt, max_pt, reg_param_final, self.file_format)
      new_name = 'hShapeClosure{}_{}_{}-{}_final{}'.format(
        parameter, label, min_pt, max_pt, self.file_format)
      shutil.copy(os.path.join(outputdir_test, old_name), os.path.join(outputdir, new_name))
      
      parameter = self.utils.remove_periods(self.prior2_variation_parameter)
      outputdir_test = os.path.join(outputdir_main, 'Test_ShapeClosure{}'.format(parameter))
      old_name = 'hClosure_pt_{}_{}{}'.format(label, reg_param_final, self.file_format)
      new_name = 'hShapeClosure{}_pt_{}_final{}'.format(parameter, label, self.file_format)
      shutil.copy(os.path.join(outputdir_test, old_name), os.path.join(outputdir, new_name))
      old_name = 'hClosure_{}_{}-{}_{}{}'.format(
        label, min_pt, max_pt, reg_param_final, self.file_format)
      new_name = 'hShapeClosure{}_{}_{}-{}_final{}'.format(
        parameter, label, min_pt, max_pt, self.file_format)
      shutil.copy(os.path.join(outputdir_test, old_name), os.path.join(outputdir, new_name))

  #----------------------------------------------------------------------
  def load_2D_observables(self, jetR, obs_label, obs_setting, grooming_setting, reg_param):

    # Get all other systematic variations, and store as attributes
    for systematic in self.systematics_list:

      output_dir = getattr(self, 'output_dir_{}'.format(systematic))
      path = os.path.join(output_dir, 'fResult_R{}_{}.root'.format(jetR, obs_label))
      f = ROOT.TFile(path, 'READ')
      name = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, reg_param)
      self.retrieve_histo_and_set_attribute(name, f, systematic)

      if systematic == 'main':
        # Get regularization parameter variations, and store as attributes
        name = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, reg_param+2)
        self.retrieve_histo_and_set_attribute(name, f)
        name = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, reg_param-2)
        self.retrieve_histo_and_set_attribute(name, f)

  #----------------------------------------------------------------------
  def retrieve_histo_and_set_attribute(self, name, f, suffix = ''):

    h = f.Get(name)
    h.SetDirectory(0)
    setattr(self, '{}{}'.format(name, suffix), h)

  #----------------------------------------------------------------------
  # Get 1D histograms
  # Normalize by integral, i.e. N_jets,inclusive in this pt-bin
  #----------------------------------------------------------------------
  def load_1D_observables(self, jetR, obs_label, obs_setting, grooming_setting,
                          reg_param, min_pt_truth, max_pt_truth, minbin, maxbin):

    # Get all other systematic variations, and store as attributes
    for systematic in self.systematics_list:

      name2D = 'hUnfolded_{}_R{}_{}_{}{}'.format(self.observable, jetR, obs_label,
                                                 reg_param, systematic)
      name1D = 'h{}_{}_R{}_{}_n{}_{}-{}'.format(systematic, self.observable, jetR,
                                            obs_label, reg_param, min_pt_truth, max_pt_truth)
      store_tagging_fraction = (systematic == 'main')
      self.get_obs_distribution(jetR, obs_label, name2D, name1D, reg_param, grooming_setting,
                                min_pt_truth, max_pt_truth, minbin, maxbin,
                                store_tagging_fraction=store_tagging_fraction)

      if systematic == 'main':
        # Get regularization parameter variations, and store as attributes
        name2D = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, reg_param+2)
        name1D = 'hRegParam1_{}_R{}_{}_n{}_{}-{}'.format(self.observable, jetR, obs_label,
                                                         reg_param, min_pt_truth, max_pt_truth)
        hRegParam1 = self.get_obs_distribution(jetR, obs_label, name2D, name1D, reg_param, grooming_setting,
                                               min_pt_truth, max_pt_truth, minbin, maxbin)

        name2D = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, reg_param-2)
        name1D = 'hRegParam2_{}_R{}_{}_n{}_{}-{}'.format(self.observable, jetR, obs_label,
                                                         reg_param, min_pt_truth, max_pt_truth)
        hRegParam2 = self.get_obs_distribution(jetR, obs_label, name2D, name1D, reg_param, grooming_setting,
                                               min_pt_truth, max_pt_truth, minbin, maxbin)

  #----------------------------------------------------------------------
  # Compute systematics
  #----------------------------------------------------------------------
  def compute_obs_systematic(self, jetR, obs_label, obs_setting, grooming_setting,
                             reg_param, min_pt_truth, max_pt_truth, minbin=None, maxbin=None, final=False):

    # Get main result
    name = 'h{}_{}_R{}_{}_n{}_{}-{}'.format('main', self.observable, jetR,
                                            obs_label, reg_param, min_pt_truth, max_pt_truth)
    hMain = getattr(self, name)
    if final:
      # Also save under name without reg param, won't need this info later
      name = 'h{}_{}_R{}_{}_{}-{}'.format('main', self.observable, jetR,
                                              obs_label, min_pt_truth, max_pt_truth)
      setattr(self, name, hMain)

    # Loop through all systematic variations, and take ratio to main result
    h_list = []
    h_unfolding_list = []
    for systematic in self.systematics_list:

      if final:

        # ** For subjet JEWEL systematic, set a larger reg param so that it can converge
        reg_param_original = reg_param
        if not self.is_pp and 'subjet_z' in self.observable and 'generator1' in systematic:
          reg_param = self.max_reg_param

        h_systematic_ratio = self.retrieve_systematic(systematic, jetR, obs_label,
                                                      reg_param, min_pt_truth, max_pt_truth)

        # ** Reset reg param
        reg_param = reg_param_original

      else:
        h_systematic_ratio = self.construct_systematic(systematic, hMain, jetR, obs_label, obs_setting, grooming_setting,
                                                       reg_param, min_pt_truth, max_pt_truth, minbin, maxbin)
      if h_systematic_ratio:
        if systematic in ['main', 'prior1', 'truncation', 'binning']:
          h_unfolding_list.append(h_systematic_ratio)
        else:
          h_list.append(h_systematic_ratio)
            
    # Construct or retrieve unfolding uncertainty
    name = 'hSystematic_Unfolding_R{}_{}_n{}_{}-{}'.format(
               self.utils.remove_periods(jetR), obs_label,
               reg_param, int(min_pt_truth), int(max_pt_truth))
    if final:
      hSystematic_Unfolding = getattr(self, name)
      name = 'hSystematic_Unfolding_R{}_{}_{}-{}'.format(self.utils.remove_periods(jetR), obs_label,
                                                         int(min_pt_truth), int(max_pt_truth))
      setattr(self, name, hSystematic_Unfolding)
      
    else:
      hSystematic_Unfolding = self.construct_unfolding_uncertainty(h_unfolding_list)
      setattr(self, name, hSystematic_Unfolding)
    h_list.append(hSystematic_Unfolding)

    # Construct total systematic uncertainty: Add all systematic uncertainties in quadrature
    name = 'hSystematic_Total_R{}_{}_n{}_{}-{}'.format(
                    self.utils.remove_periods(jetR), obs_label,
                    reg_param, int(min_pt_truth), int(max_pt_truth))
    hSystematic_Total = self.add_in_quadrature(h_list, new_name=name)
    setattr(self, name, hSystematic_Total)

    # Attach total systematic to main result, and save as an attribute
    name = 'hResult_{}_systotal_R{}_{}_n{}_{}-{}'.format(
                    self.observable, jetR, obs_label, reg_param,
                    int(min_pt_truth), int(max_pt_truth))
    if grooming_setting and maxbin:
      hResult_sys = self.truncate_hist(hMain.Clone(), minbin, maxbin+1, name)
    else:
      hResult_sys = self.truncate_hist(hMain.Clone(), minbin, maxbin, name)
    hResult_sys.SetDirectory(0)
    self.AttachErrToHist(hResult_sys, hSystematic_Total)
    setattr(self, name, hResult_sys)

    if final:
      # Save main result under name without reg param
      name = 'hResult_{}_systotal_R{}_{}_n{}_{}-{}'.format(
        self.observable, jetR, obs_label, reg_param, 
        int(min_pt_truth), int(max_pt_truth))
      hResult_sys = getattr(self, name)
      name = 'hResult_{}_systotal_R{}_{}_{}-{}'.format(
        self.observable, jetR, obs_label,
        int(min_pt_truth), int(max_pt_truth))
      setattr(self, name, hResult_sys)
      
      if self.debug_level > 0:
        name = 'hSystematic_Total_R{}_{}_{}-{}{}'.format(
          self.utils.remove_periods(jetR), obs_label,
          int(min_pt_truth), int(max_pt_truth), self.file_format)
        outputFilename = os.path.join(self.output_dir_systematics, name)
        self.utils.plot_hist(hSystematic_Total, outputFilename, 'P E')

      # Plot systematic uncertainties, and write total systematic to a ROOT file
      self.plot_systematic_uncertainties(
        jetR, obs_label, obs_setting, grooming_setting,
        min_pt_truth, max_pt_truth, minbin, maxbin, h_list, hSystematic_Total)

      # Plot unfolding uncertainties
      self.plot_systematic_uncertainties(
        jetR, obs_label, obs_setting, grooming_setting,
        min_pt_truth, max_pt_truth, minbin, maxbin, h_unfolding_list, hSystematic_Unfolding, suffix='Unfolding')

  #----------------------------------------------------------------------
  # Retrieve a given systematic histogram
  def retrieve_systematic(self, systematic, jetR, obs_label,
                           reg_param, min_pt_truth, max_pt_truth):

    if systematic == 'main':
      sys_label = 'RegParam'
    elif systematic in ['prior1', 'prior2']:
      if systematic == 'prior1':
        sys_label = 'prior'
      else:
        return None
    elif systematic in ['subtraction1', 'subtraction2']:
      if systematic == 'subtraction1':
        sys_label = 'subtraction'
      else:
        return None
    elif systematic in ['fastsim_generator0', 'fastsim_generator1']:
      if systematic == 'fastsim_generator1':
        sys_label = 'generator'
      else:
        return None
    else:
      sys_label = systematic
    if reg_param:
        name = 'hSystematic_{}_{}_R{}_{}_n{}_{}-{}'.format(
                self.observable, sys_label, jetR, obs_label,
                reg_param, min_pt_truth, max_pt_truth)
    else:
        name = 'hSystematic_{}_{}_R{}_{}_{}-{}'.format(
                self.observable, sys_label, jetR, obs_label,
                min_pt_truth, max_pt_truth)

    h_systematic_ratio = getattr(self, name)
    
    # Save final sys uncertainty breakdowns
    name = 'hSystematic_{}_{}_R{}_{}_{}-{}'.format(self.observable, sys_label,
                                   jetR, obs_label, min_pt_truth, max_pt_truth)
    setattr(self, name, h_systematic_ratio)

    if self.debug_level > 0:
      output_dir = getattr(self, 'output_dir_systematics')
      name = 'hSystematic_{}_R{}_{}_{}-{}{}'.format(
                 systematic, self.utils.remove_periods(jetR), obs_label,
                 int(min_pt_truth), int(max_pt_truth), self.file_format)
      outputFilename = os.path.join(output_dir, name)
      self.utils.plot_hist(h_systematic_ratio, outputFilename, 'P E')      

    return h_systematic_ratio

  #----------------------------------------------------------------------
  # Get systematic variation and save percentage difference as attribte
  def construct_systematic(self, systematic, hMain, jetR, obs_label, obs_setting, grooming_setting,
                           reg_param, min_pt_truth, max_pt_truth, minbin, maxbin):
                           
    # Set whether to store signed uncertainty value or absolute value
    # For now, only use signed uncertainty in cases where we don't average/combine multiple sources
    signed = systematic in ['trkeff', 'fastsim_generator1']

    if grooming_setting and maxbin:
      maxbin += 1

    # Combine certain systematics as average or max
    if systematic == 'main':
      h_systematic_ratio = self.construct_systematic_average(
            hMain, 'RegParam', jetR, obs_label, reg_param,
            min_pt_truth, max_pt_truth, minbin, maxbin, takeMaxDev=False)
        
    elif systematic in ['prior1', 'prior2']:
      if systematic == 'prior1':
        h_systematic_ratio = self.construct_systematic_average(
              hMain, 'prior', jetR, obs_label, reg_param,
              min_pt_truth, max_pt_truth, minbin, maxbin, takeMaxDev=True)
      else:
        return None
            
    # For model-depedence, take the difference between two fastsim generators
    elif systematic in ['fastsim_generator0', 'fastsim_generator1']:
      if systematic == 'fastsim_generator1':
        # Get reference generator unfolded result (e.g. pythia fastsim)
        name = 'h{}_{}_R{}_{}_n{}_{}-{}'.format(
              'fastsim_generator0', self.observable, jetR, obs_label,
              reg_param, min_pt_truth, max_pt_truth)
        h_reference = getattr(self, name)
            
        # Take the difference of generator to reference generator (e.g. herwig fastsim)
        h_systematic_ratio = self.construct_systematic_percentage(
              h_reference, 'fastsim_generator1', jetR, obs_label,
              reg_param, min_pt_truth, max_pt_truth, minbin, maxbin, signed=signed)
      else:
        return None
        
    elif systematic in ['subtraction1', 'subtraction2']:
      if systematic == 'subtraction1':
        h_systematic_ratio = self.construct_systematic_average(
              hMain, 'subtraction', jetR, obs_label, reg_param,
              min_pt_truth, max_pt_truth, minbin, maxbin, takeMaxDev=True)
      else:
        return None
        
    elif systematic == 'thermal_closure':
      fname = 'nonclosureR{}_{}_{}_{}-{}.root'.format(jetR, obs_setting, grooming_setting, int(min_pt_truth), int(max_pt_truth))
      outf_name = os.path.join(getattr(self, 'output_dir_thermal_closure'), 'Test_ThermalClosure')
      f_nonclosure = ROOT.TFile(os.path.join(outf_name, fname), 'READ')
      h_systematic_ratio_temp = f_nonclosure.Get('hNonclosureRatio')
      h_systematic_ratio_temp.SetDirectory(0)
      f_nonclosure.Close()
          
      name_ratio = 'hSystematic_{}_{}_R{}_{}_n{}_{}-{}'.format(
            self.observable, systematic, jetR, obs_label,
            reg_param, min_pt_truth, max_pt_truth)
      self.change_to_per(h_systematic_ratio_temp)
      h_systematic_ratio_temp.SetMinimum(0.)
      h_systematic_ratio_temp.SetMaximum(1.5*h_systematic_ratio_temp.GetMaximum())
      h_systematic_ratio = self.truncate_hist(h_systematic_ratio_temp, minbin, maxbin, name_ratio)
      setattr(self, name_ratio, h_systematic_ratio)
        
    else:
      h_systematic_ratio = self.construct_systematic_percentage(
            hMain, systematic, jetR, obs_label, reg_param,
            min_pt_truth, max_pt_truth, minbin, maxbin, signed=signed)
              
    name = 'hSystematic_{}_{}_R{}_{}_n{}_{}-{}'.format(self.observable, systematic, jetR, obs_label, reg_param, min_pt_truth, max_pt_truth)
    setattr(self, name, h_systematic_ratio)
              
    return h_systematic_ratio

  #----------------------------------------------------------------------
  # Get systematic variation and save percentage difference as attribte
  def construct_systematic_percentage(self, hMain, systematic, jetR,
                                      obs_label, reg_param, min_pt_truth,
                                      max_pt_truth, minbin, maxbin, signed=False):

    name = 'h{}_{}_R{}_{}_n{}_{}-{}'.format(systematic, self.observable, jetR, obs_label,
                                            reg_param, min_pt_truth, max_pt_truth)
    h_systematic = getattr(self, name)

    if 'fastsim_generator1' in systematic:
        systematic = 'generator'  # For generator systematic, need to set label here
        
    # Normalization
    #integral = hMain.Integral(2, hMain.GetNbinsX(), 'width')
    #hMain.Scale(1./integral)
    #integral = h_systematic.Integral(2, h_systematic.GetNbinsX(), 'width')
    #h_systematic.Scale(1./integral)
    
    name_ratio = 'hSystematic_{}_{}_R{}_{}_n{}_{}-{}'.format(
      self.observable, systematic, jetR, obs_label,
      reg_param, min_pt_truth, max_pt_truth)
    h_systematic_ratio_temp = hMain.Clone()
    h_systematic_ratio_temp.SetName(name_ratio+'_temp')
    h_systematic_ratio_temp.Divide(h_systematic)

    if self.debug_level > 0:
      print("Printing systematic variation ratio for", systematic)
      for i in range(1, h_systematic.GetNbinsX()+1):
        print("main:", hMain.GetBinContent(i), "-- sys:", h_systematic.GetBinContent(i),
              "-- ratio:", h_systematic_ratio_temp.GetBinContent(i))

    self.change_to_per(h_systematic_ratio_temp, signed=signed)
    h_systematic_ratio = self.truncate_hist(h_systematic_ratio_temp, minbin, maxbin, name_ratio)
    del h_systematic_ratio_temp   # No longer need this -- prevents memory leaks
    setattr(self, name_ratio, h_systematic_ratio)
    return h_systematic_ratio

  #----------------------------------------------------------------------
  def construct_systematic_average(self, hMain, sys_label, jetR, obs_label,
                                  reg_param, min_pt_truth, max_pt_truth,
                                  minbin, maxbin, takeMaxDev=False):
  
    h_systematic_ratio1 = self.construct_systematic_percentage(
      hMain, '{}1'.format(sys_label), jetR, obs_label, reg_param,
      min_pt_truth, max_pt_truth, minbin, maxbin)
    h_systematic_ratio2 = self.construct_systematic_percentage(
      hMain, '{}2'.format(sys_label), jetR, obs_label,
      reg_param, min_pt_truth, max_pt_truth, minbin, maxbin)

    name = 'hSystematic_{}_{}_R{}_{}_n{}_{}-{}'.format(self.observable, sys_label, jetR, obs_label,
                                                       reg_param, min_pt_truth, max_pt_truth)
    h_systematic_ratio = self.build_average(h_systematic_ratio1, h_systematic_ratio2,
                                            takeMaxDev=takeMaxDev)
    h_systematic_ratio.SetName(name)
    setattr(self, name, h_systematic_ratio)
    return h_systematic_ratio

  #----------------------------------------------------------------------
  def get_obs_distribution(self, jetR, obs_label, name2D, name1D, reg_param, grooming_setting,
                           min_pt_truth, max_pt_truth, minbin, maxbin, store_tagging_fraction=False):

    h2D = getattr(self, name2D)
    h2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    h = h2D.ProjectionY() # Better to use ProjectionY('{}_py'.format(h2D.GetName()), 1, h2D.GetNbinsX()) ?
    h.SetName(name1D)
    h.SetDirectory(0)

    if not minbin:
      minbin = 1
    if not maxbin:
      maxbin = h.GetNbinsX()
        
    # Normalize by integral, i.e. N_jets,inclusive in this pt-bin
    
    # First scale by bin width -- then normalize by integral
    # (where integral weights by bin width)
    h.Scale(1., 'width')

    if grooming_setting and 'sd' in grooming_setting:
    
      # If SD, the untagged jets are in the first bin
      n_jets_inclusive = h.Integral(minbin, maxbin, 'width')
      n_jets_tagged = h.Integral(minbin+1, maxbin, 'width')
      
      if store_tagging_fraction:
        f_tagging = n_jets_tagged/n_jets_inclusive
        f_tagging_name = 'tagging_fraction_R{}_{}_{}-{}_{}'.format(
          jetR, obs_label, min_pt_truth, max_pt_truth, reg_param)
        setattr(self, f_tagging_name, f_tagging)
      
    else:
      n_jets_inclusive = h.Integral(minbin, maxbin, 'width')
    
    h.Scale(1./n_jets_inclusive)
    
    setattr(self, name1D, h)
    
    return h

  #----------------------------------------------------------------------
  def plot_systematic_uncertainties(self, jetR, obs_label, obs_setting, grooming_setting,
                                    min_pt_truth, max_pt_truth, minbin, maxbin, h_list, h_total, suffix=''):

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    name = 'cSys_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt_truth, max_pt_truth)
    c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()

    c.cd()
    myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    myPad.SetLeftMargin(0.2)
    myPad.SetTopMargin(0.07)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.13)
    myPad.Draw()
    myPad.cd()

    truth_bin_array = self.truth_bin_array(obs_label)
    if maxbin:
      truth_bin_array = truth_bin_array[0:maxbin+1]
    if minbin:
      truth_bin_array = truth_bin_array[minbin:]
    n_bins_truth = len(truth_bin_array) - 1

    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_bins_truth, truth_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle( getattr(self, 'xtitle') )
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle('Systematic uncertainty (%)')
    myBlankHisto.SetMaximum(2.7*h_total.GetMaximum(50))
    if not suffix=="Unfolding":
      myBlankHisto.SetMinimum(-1.1*h_total.GetMaximum(50))
    else:
      # Unfolding uncertainties do not go below 0
      myBlankHisto.SetMinimum(0)
    myBlankHisto.Draw("E")

    leg = ROOT.TLegend(0.67,0.6,0.8,0.92)
    self.utils.setup_legend(leg,0.04)

    for i, h in enumerate(h_list):
      if h:
        h.SetMarkerStyle(self.MarkerArray[i])
        h.SetMarkerSize(1.5)
        h.SetMarkerColor(self.ColorArray[i])
        h.SetLineColor(self.ColorArray[i])
        h.SetLineStyle(1)
        h.SetLineWidth(2)
        if h.GetMaximum() > h_total.GetMaximum():
          myBlankHisto.SetMaximum(1.7*h.GetMaximum())

        h.DrawCopy('P X0 same')

        legend_label = ''
        for systematic in self.systematics_list:
          if 'Unfolding' in h.GetName():
            legend_label = 'unfolding'
          elif systematic in h.GetName():
            legend_label = systematic
          elif 'RegParam' in h.GetName():
            legend_label = 'reg param'
          elif 'prior' in h.GetName():
            legend_label = 'prior'
          elif 'generator' in h.GetName():
            legend_label = 'generator'
          elif 'subtraction' in h.GetName():
            legend_label = 'subtraction'
        leg.AddEntry(h, legend_label, 'P')

    h_total.SetLineStyle(1)
    h_total.SetLineColor(1)
    h_total.SetLineWidth(2)
    h_total.DrawCopy('same hist')
    leg.AddEntry(h_total, 'Total {}'.format(suffix), 'l')

    leg.Draw()

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth)
    text_latex.DrawLatex(0.3, 0.85, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = '#it{R} = ' + str(jetR)
    text_latex.DrawLatex(0.3, 0.78, text)

    subobs_label = self.utils.formatted_subobs_label(self.observable)
    if subobs_label:
      text = '{} = {}'.format(subobs_label, obs_setting)
      text_latex.DrawLatex(0.3, 0.71, text)

    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting)
      if subobs_label and '#it{#beta}' in subobs_label:
        text.replace('#it{#beta}', '#it{#beta}_{SD}')
      text_latex.DrawLatex(0.3, 0.64, text)

    output_dir = getattr(self, 'output_dir_systematics')
    outputFilename = os.path.join(output_dir, 'hSystematics_R{}_{}_{}-{}{}{}'.format(
      self.utils.remove_periods(jetR), obs_label, int(min_pt_truth),
      int(max_pt_truth), suffix, self.file_format))
    c.SaveAs(outputFilename)
    c.Close()

    if not suffix:
        sys_root_filename = os.path.join(output_dir, 'fSystematics.root')
        fSystematics = ROOT.TFile(sys_root_filename, 'UPDATE')
        h_total.Write()
        for h in h_list:
            h.Write()
        fSystematics.Close()

  #----------------------------------------------------------------------
  # Returns truncated 1D histogram from bins [minbin, ..., maxbin] inclusive
  # (Note this uses 1-indexed bin number to comply with ROOT)
  #----------------------------------------------------------------------
  def truncate_hist(self, h, minbin, maxbin, new_name):
    length = h.GetNbinsX()

    # Check if either minbin or maxbin exist, and if so set bin indices 
    if maxbin == None and minbin == None:
      h.SetNameTitle(new_name, new_name)
      return h
    else:
      if maxbin == None:
        bin_range = range(minbin+1, length+2)
        if minbin >= length:
          raise ValueError(f"Min bin number {minbin} larger or equal to histogram size {length}")
        if minbin < 1:
          raise ValueError(f"Min bin number {minbin} cannot be less than 1")

      elif minbin == None:
        bin_range = range(1, maxbin+2)
        if maxbin >= length:
          raise ValueError(f"Max bin number {maxbin} larger or equal to histogram size {length}")
        if maxbin < 1:
          raise ValueError(f"Max bin number {maxbin} cannot be less than 1")

      bin_edges = array('d', [h.GetXaxis().GetBinLowEdge(i) for i in bin_range])
      return h.Rebin(len(bin_edges)-1, new_name, bin_edges)

  #----------------------------------------------------------------------
  # Add a list of (identically-binned) histograms in quadrature, bin-by-bin
  #----------------------------------------------------------------------
  def add_in_quadrature(self, h_list, new_name=None):

    h_new = h_list[0].Clone()
    if not new_name:
        new_name = '{}_new'.format(h_list[0].GetName())
    h_new.SetName(new_name)

    for i in range(1, h_new.GetNbinsX()+1):

      values_i = [h.GetBinContent(i) for h in h_list]

      new_value_squared = 0.
      for value_i in values_i:
        new_value_squared += value_i*value_i
      new_value = math.sqrt(new_value_squared)
      h_new.SetBinContent(i, new_value)

    return h_new

  #----------------------------------------------------------------------
  # Take stdev of a list of (identically-binned) histograms, bin-by-bin
  #----------------------------------------------------------------------
  def construct_unfolding_uncertainty(self, h_list):

    h_new = h_list[0].Clone()
    h_new.SetName('hCombinedUnfoldingSystematic_{}'.format(h_list[0].GetName()))

    for i in range(1, h_new.GetNbinsX()+1):

      values_i = [h.GetBinContent(i) for h in h_list]

      new_value_squared = 0.
      for value_i in values_i:
        new_value_squared += value_i*value_i
      new_value = math.sqrt(new_value_squared / len(h_list))
      h_new.SetBinContent(i, new_value)

    return h_new

  #----------------------------------------------------------------------
  # Finds the regularization parameter than minimizes statistical and
  # all calculated systematic uncertainties when added in quadrature
  # 
  def determine_reg_param_final(self, jetR, obs_label, obs_setting, grooming_setting,
                                min_pt_truth, max_pt_truth, minbin, maxbin):

    assert self.use_max_reg_param == True

    # Obtain the total systematic uncertainty histograms from memory
    min_reg_param = 3
    reg_params = range(min_reg_param, self.max_reg_param + 1)
    names = ["hSystematic_Unfolding_R{}_{}_n{}_{}-{}".format(
      self.utils.remove_periods(jetR), obs_label, reg_param, int(min_pt_truth),
      int(max_pt_truth)) for reg_param in reg_params]
    hSys_list = [getattr(self, name) for name in names]

    # Obtain the statistical uncertainty histograms from disk
    fdir = os.path.join(self.output_dir_main, "Unfolded_stat_uncert",
                        "fResult_name_R%s_%s.root" % (jetR, obs_label))
    f = ROOT.TFile.Open(fdir, "READ")
    names = ["hUnfolded_%s_stat_uncert_R%s_%s_n%i_Pt%i-%i" % \
             (self.observable, self.utils.remove_periods(jetR), obs_label,
              reg_param, min_pt_truth, max_pt_truth) for reg_param in reg_params]
    hStat_list_untruncated = [f.Get(name) for name in names]
    hStat_list = [self.truncate_hist(h, minbin, maxbin, name+'_trunc')
                  for (h, name) in zip(hStat_list_untruncated, names)]

    # Add statistical and systematic uncertainties in quadrature
    hUncert_list = [self.add_in_quadrature([hSys, hStat]) for (hSys, hStat)
                    in zip(hSys_list, hStat_list)]

    # Sum total uncertainties per observable bin and return lowest value
    hUncert_sum_list = [sum([h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)])
                        for h in hUncert_list]
    index_min = min(range(len(hUncert_sum_list)), key=hUncert_sum_list.__getitem__)
    return index_min + min_reg_param

    f.Close()
    
  #----------------------------------------------------------------------
  # Write HEPData submission
  # You can test it here: https://www.hepdata.net/record/sandbox
  #
  # You will want to edit the tables slightly to add observable-specific info:
  #   units, description, etc.
  #----------------------------------------------------------------------
  def write_hepdata(self):
  
    # Create submission
    self.hepdata_dir = os.path.join(getattr(self, 'output_dir_final_results'), 'hepdata')
    if not os.path.exists(self.hepdata_dir):
      os.makedirs(self.hepdata_dir)

    self.hepdata_submission = hepdata_lib.Submission()
    self.hepdata_dict = {}
  
    # Loop through jet radii
    for jetR in self.jetR_list:

      # Loop through subconfigurations
      for i, subconfig in enumerate(self.obs_subconfig_list):

        obs_setting = self.obs_settings[i]
        grooming_setting = self.grooming_settings[i]
        obs_label = self.obs_labels[i]

        # Loop through pt slices, and add table for each result
        for bin in range(0, len(self.pt_bins_reported) - 1):
          min_pt_truth = self.pt_bins_reported[bin]
          max_pt_truth = self.pt_bins_reported[bin+1]
          
          self.add_hepdata_table(i, jetR, obs_label, obs_setting, grooming_setting,
                                  min_pt_truth, max_pt_truth)
                        
    # Write tables to a dict -- this is useful if create_files fails, since can construct submission then with write_hepdata.py script
    with open(os.path.join(self.hepdata_dir, 'tables.pkl'), 'wb') as f:
      pickle.dump(self.hepdata_dict, f)

    # Write submission files
    self.hepdata_submission.create_files(self.hepdata_dir)

  #----------------------------------------------------------------------
  def add_hepdata_table(self, i, jetR, obs_label, obs_setting, grooming_setting, min_pt, max_pt):

    index = self.set_hepdata_table_index(i, jetR, min_pt, grooming_setting, obs_setting)
    table = hepdata_lib.Table(f'Table {index}')
    table.keywords["reactions"] = ['P P --> jet+X']
    table.keywords["cmenergies"] = ['5020']
    
    # Set observable-specific info in table
    x_label, y_label = self.set_hepdata_table_descriptors(table, jetR, obs_label, obs_setting, grooming_setting, min_pt, max_pt)
 
    # Create readers to read histograms
    final_result_root_filename = os.path.join(getattr(self, 'output_dir_final_results'), 'fFinalResults.root')
    hepdata_reader = hepdata_lib.RootFileReader(final_result_root_filename)
 
    systematics_root_filename = os.path.join(self.output_dir_systematics, 'fSystematics.root')
    hepdata_reader_systematics = hepdata_lib.RootFileReader(systematics_root_filename)
 
    # Define variables
    h_name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt, max_pt)
    if self.observable == 'ang':
        h_name += '_trunc'
    h = hepdata_reader.read_hist_1d(h_name)
 
    n_significant_digits = 3 # Note: use 2 significant digits for uncertainties
    
    x = hepdata_lib.Variable(x_label, is_independent=True, is_binned=True, units='')
    x.digits = n_significant_digits
    x.values = h['x_edges']
 
    y = hepdata_lib.Variable(y_label, is_independent=False, is_binned=False, units='')
    y.digits = n_significant_digits
    y.values = h['y']
    if self.is_pp:
        y.add_qualifier('RE', 'P P --> jet+X')
    else:
        y.add_qualifier('RE', 'Pb Pb --> jet+X')
    y.add_qualifier('SQRT(S)', 5.02, 'TeV')
    y.add_qualifier('ETARAP', '|0.9-R|')
    y.add_qualifier('jet radius', jetR)
    y.add_qualifier('jet method', 'Anti-$k_{T}$')
 
    # Define uncertainties
    stat = hepdata_lib.Uncertainty('stat', is_symmetric=True)
    stat.values = [float('{:.2g}'.format(dy)) for dy in h['dy']]
 
    # Add tables to submission
    table.add_variable(x)
    table.add_variable(y)
    y.add_uncertainty(stat)
 
    # Add unfolding systematic
    name = 'hSystematic_Unfolding_R{}_{}_{}-{}'.format(self.utils.remove_periods(jetR), obs_label,
                                                      int(min_pt), int(max_pt))
    h_sys_unfolding = hepdata_reader_systematics.read_hist_1d(getattr(self, name).GetName())
    sys_unfolding = hepdata_lib.Uncertainty('sys,unfolding', is_symmetric=True)
    sys_unfolding.values = ['{:.2g}%'.format(y) for y in h_sys_unfolding['y']]
    y.add_uncertainty(sys_unfolding)

    # Add systematic uncertainty breakdown
    for systematic in self.systematics_list:
 
        if systematic in ['main', 'prior1', 'truncation', 'binning']:
            continue

        h_sys = self.retrieve_systematic(systematic, jetR, obs_label,
                                         None, min_pt, max_pt)
        if not h_sys:
            continue
            
        if 'generator' in systematic:
          sys_label = 'generator'
        else:
          sys_label = systematic
        
        h_sys = hepdata_reader_systematics.read_hist_1d(h_sys.GetName())
        sys = hepdata_lib.Uncertainty('sys,{}'.format(sys_label), is_symmetric=True)
        sys.values = ['{:.2g}%'.format(y) for y in h_sys['y']]
        y.add_uncertainty(sys)
 
    # Add table to the submission
    self.hepdata_dict[index] = table
    self.hepdata_submission.add_table(table)

  #----------------------------------------------------------------------
  def set_hepdata_table_index(self, i, jetR, min_pt, grooming_setting='', obs_setting=None):
  
    index = 1
    index += i
    if self.observable == 'ang':
        if np.isclose(jetR, 0.4):
            if np.isclose(min_pt, 40.):
                index += 8
            elif np.isclose(min_pt, 60.):
                index += 16
            elif np.isclose(min_pt, 80.):
                index += 24
        elif np.isclose(jetR, 0.2):
            if np.isclose(min_pt, 20.):
                index += 32
            if np.isclose(min_pt, 40.):
                index += 40
            if np.isclose(min_pt, 60.):
                index += 48
            if np.isclose(min_pt, 80.):
                index += 56

    elif 'DG_a01' in self.obs_labels:

        if 'sd' in grooming_setting:
          beta = grooming_setting['sd'][1]

          if self.observable == 'zg':
            index = 1 + beta # Fig 2
          elif self.observable == 'theta_g':
            index = 4 + beta # Fig 3

        elif 'dg' in grooming_setting:
          a = grooming_setting['dg'][0]
                
          if self.observable == 'zg':
            index = 7 + int(a) # Fig 5
          elif self.observable == 'theta_g':
            index = 10 + int(a) # Fig 6

    elif 'subjet_z' in self.observable:
            
      if self.observable == 'leading_subjet_z':
        if np.isclose(obs_setting, 0.1):
          index = 3
        elif np.isclose(obs_setting, 0.2):
          index = 4
          
      elif self.observable == 'inclusive_subjet_z':
        if np.isclose(obs_setting, 0.1):
          index = 1
        elif np.isclose(obs_setting, 0.2):
          index = 2


    return index

  #----------------------------------------------------------------------
  def set_hepdata_table_descriptors(self, table, jetR, obs_label, obs_setting, grooming_setting, min_pt, max_pt):
    
    if self.observable == 'ang':
    
        if grooming_setting and 'sd' in grooming_setting.keys():
        
            if np.isclose(jetR, 0.4):
                table.location = 'Figure 3.'
            elif np.isclose(jetR, 0.2):
                table.location = 'Figure 4.'
                
            table.description = r'Groomed jet angularity $\lambda_{{\alpha,g}}$ for $\alpha = {}$.'.format(obs_setting)
            table.description += '\n'
            table.description += r'${}<p_{{\mathrm{{T}}}}^{{\mathrm{{ch jet}}}}<{}$, Soft Drop $z_{{\mathrm{{cut}}}}=0.2, \beta=0$.'.format(min_pt, max_pt)
            table.description += '\n\nNote: The first bin corresponds to the Soft Drop untagged fraction.'
            
            x_label = r'$\lambda_{\alpha,g}$'
            y_label = r'$\frac{1}{\sigma_{inc}} \frac{d\sigma}{d\lambda_{\alpha,g}}$'
            
        else:
        
            if np.isclose(jetR, 0.4):
                table.location = 'Figure 1.'
            elif np.isclose(jetR, 0.2):
                table.location = 'Figure 2.'

            table.description = r'Jet angularity $\lambda_{{\alpha}}$ for $\alpha = {}$.'.format(obs_setting)
            table.description += '\n'
            table.description += r'${}<p_{{\mathrm{{T}}}}^{{\mathrm{{ch\;jet}}}}<{}$ GeV/$c$.'.format(min_pt, max_pt)
            
            x_label = r'$\lambda_{{\alpha}}$'
            y_label = r'$\frac{1}{\sigma} \frac{d\sigma}{d\lambda_{\alpha}}$'
          
        table.description += '\n\n'
        table.description += r'For the "trkeff" and "generator" systematic uncertainty sources, the signed systematic uncertainty breakdowns ($\pm$ vs. $\mp$), denote correlation across bins (both within this table, and across tables). For the remaining sources ("unfolding", "random_mass") no correlation information is specified ($\pm$ is always used).'
            
        if np.isclose(min_pt, 20.):
            table.location += ' upper left'
        elif np.isclose(min_pt, 40.):
            table.location += ' upper right'
        elif np.isclose(min_pt, 60.):
            table.location += ' lower left'
        elif np.isclose(min_pt, 80.):
            table.location += ' lower right'
            
    elif 'DG_a01' in self.obs_labels:
            
        if self.observable == 'zg':
            table.description = r'Groomed jet momentum splitting fraction $z_{{\mathrm{g}}}$'
            x_label = r'$z_{{\mathrm{g}}}$'
            y_label = r'$\frac{1}{\sigma_{\mathrm{jet}}} \frac{d\sigma}{dz_{{\mathrm{g}}}}$'
            if 'sd' in grooming_setting:
                table.location = 'Figure 2'
            elif 'dg' in grooming_setting:
                table.location = 'Figure 5'
            
        elif self.observable == 'theta_g':
            table.description = r'Groomed jet radius (scaled) $\theta_{{\mathrm{g}}}$'
            x_label = r'$\theta_{{\mathrm{g}}}$'
            y_label = r'$\frac{1}{\sigma_{\mathrm{jet}}} \frac{d\sigma}{d\theta_{{\mathrm{g}}}}$'
            if 'sd' in grooming_setting:
                table.location = 'Figure 3'
            elif 'dg' in grooming_setting:
                table.location = 'Figure 6'
        
        table.description += '\n'
        if 'sd' in grooming_setting:
          zcut = grooming_setting['sd'][0]
          beta = grooming_setting['sd'][1]
          table.description += r'${}<p_{{\mathrm{{T}}}}^{{\mathrm{{ch\;jet}}}}<{}$ GeV/$c$, soft drop $z_{{\mathrm{{cut}}}}={}, \beta={}$.'.format(min_pt, max_pt, zcut, beta)
          table.description += '\n\nNote: The first bin corresponds to the Soft Drop untagged fraction.'
        elif 'dg' in grooming_setting:
          a = grooming_setting['dg'][0]
          table.description += r'${}<p_{{\mathrm{{T}}}}^{{\mathrm{{ch\;jet}}}}<{}$ GeV/$c$, dynamical grooming $a={}$.'.format(min_pt, max_pt, a)

        table.description += '\n\n'
        table.description += r'For the "trkeff" and "generator" systematic uncertainty sources, the signed systematic uncertainty breakdowns ($\pm$ vs. $\mp$), denote correlation across bins (both within this table, and across tables). For the remaining sources ("unfolding") no correlation information is specified ($\pm$ is always used).'

    elif self.observable in ['zg', 'theta_g']: # Note: theta_g PbPb paper hepdata is constructed in plot_raa.py
            
        if self.observable == 'zg':
            table.description = r'Groomed jet momentum splitting fraction $z_{{\mathrm{g}}}$'
            x_label = r'$z_{{\mathrm{g}}}$'
            y_label = r'$\frac{1}{\sigma_{inc}} \frac{d\sigma}{dz_{{\mathrm{g}}}}$'
            if np.isclose(jetR, 0.4):
                table.location = 'Figure 2 (right)'
            elif np.isclose(jetR, 0.2):
                table.location = 'Figure 2 (left)'
            
        elif self.observable == 'theta_g':
            table.description = r'Groomed jet radius (scaled) $\theta_{{\mathrm{g}}}$'
            x_label = r'$\theta_{{\mathrm{g}}}$'
            y_label = r'$\frac{1}{\sigma_{inc}} \frac{d\sigma}{d\theta_{{\mathrm{g}}}}$'
            table.location = 'Figure 4'
            
        table.description += '\n'
        table.description += r'${}<p_{{\mathrm{{T}}}}^{{\mathrm{{ch jet}}}}<{}$, Soft Drop $z_{{\mathrm{{cut}}}}=0.2, \beta=0$.'.format(min_pt, max_pt)
        table.description += '\n\nNote: The first bin corresponds to the Soft Drop untagged fraction.'
        
        table.description += '\n\n'
        if self.is_pp:
            table.description += r'For the "trkeff" and "generator" systematic uncertainty sources, the signed systematic uncertainty breakdowns ($\pm$ vs. $\mp$), denote correlation across bins (both within this table, and across tables for a given centrality). For the remaining sources ("unfolding") no correlation information is specified ($\pm$ is always used).'
        else:
            table.description += r'For the "trkeff" systematic uncertainty sources, the signed systematic uncertainty breakdowns ($\pm$ vs. $\mp$), denote correlation across bins (both within this table, and across tables for a given centrality). For the remaining sources ("unfolding", "subtraction", "thermal_closure") no correlation information is specified ($\pm$ is always used).'
    
    elif 'subjet_z' in self.observable:
            
        if self.observable == 'leading_subjet_z':
            table.description = rf'Leading subjet $z_r$ in pp collisions for $r={obs_setting}$'
            x_label = r'$z_r$'
            y_label = r'$\frac{1}{\sigma_{\mathrm{jet}}} \frac{d\sigma}{dz_r}$'
            table.location = 'Figure 3'
            
        elif self.observable == 'inclusive_subjet_z':
            table.description = rf'Inclusive subjet $z_r$ in pp collisions for $r={obs_setting}$'
            x_label = r'$z_r$'
            y_label = r'$\frac{1}{\sigma_{\mathrm{jet}}} \frac{d\sigma}{dz_r}$'
            table.location = 'Figure 2'
        
        table.description += '\n'
        table.description += rf'${min_pt}<p_{{\mathrm{{T}}}}^{{\mathrm{{ch\;jet}}}}<{max_pt}$ GeV/$c$.'

        table.description += '\n\n'
        table.description += r'For the "trkeff" and "generator" systematic uncertainty sources, the signed systematic uncertainty breakdowns ($\pm$ vs. $\mp$), denote correlation across bins (both within this table, and across tables). For the remaining sources ("unfolding") no correlation information is specified ($\pm$ is always used).'


    return x_label, y_label

  #----------------------------------------------------------------------
  def change_to_per(self, h, signed=False):

    for bin in range(0, h.GetNbinsX()+2):
      content = h.GetBinContent(bin)
      
      if signed:
        content_new = 1-content
      else:
        content_new = math.fabs(1-content)

      h.SetBinContent(bin, content_new*100)

  #----------------------------------------------------------------------
  def build_average(self, h1, h2, takeMaxDev=False):

    h_avg = h1.Clone()
    h_avg.SetName('{}_avg'.format(h1.GetName()))

    for i in range(1, h_avg.GetNbinsX()+1):
      value1 = h1.GetBinContent(i)
      value2 = h2.GetBinContent(i)
      avg =  0.5*(value1 + value2)

      if takeMaxDev:
        if value1 > value2:
          avg = value1
        else:  # value2 > value1:
          avg = value2

      h_avg.SetBinContent(i, avg)

    return h_avg

  #----------------------------------------------------------------------
  def AttachErrToHist(self, h, hPercError):

    #Fill array with lower bin edges of data histogram
    for bin in range(1, h.GetNbinsX()+1):
      content = h.GetBinContent(bin)
      perErr = hPercError.GetBinContent(bin)
      h.SetBinError(bin, content*perErr*0.01)

  #----------------------------------------------------------------------
  # This function is called once for each subconfiguration
  # You must implement this
  #----------------------------------------------------------------------
  def plot_single_result(jetR, obs_label, obs_setting, grooming_setting):

    raise NotImplementedError('You must implement plot_single_result()!')

  #----------------------------------------------------------------------
  # This function is called once after all subconfigurations have been looped over, for each R
  # You must implement this
  #----------------------------------------------------------------------
  def plot_all_results(jetR, obs_label, obs_setting, grooming_setting):

    raise NotImplementedError('You must implement plot_all_results()!')

  #----------------------------------------------------------------------
  # This function is called once after all subconfigurations and jetR have been looped over
  # You must implement this
  #----------------------------------------------------------------------
  def plot_performance(self):

    raise NotImplementedError('You must implement plot_performance()!')

  #---------------------------------------------------------------
  # Get n_bins_truth
  #---------------------------------------------------------------
  def n_bins_truth(self, obs_label):

    return getattr(self, 'n_obs_bins_truth_{}'.format(obs_label))

  #---------------------------------------------------------------
  # Get truth_bin_array
  #---------------------------------------------------------------
  def truth_bin_array(self, obs_label):

    return getattr(self, 'truth_obs_bin_array_{}'.format(obs_label))

  #---------------------------------------------------------------
  # Get max obs bin
  #---------------------------------------------------------------
  def obs_max_bins(self, obs_label):

    return getattr(self, 'obs_max_bins_{}'.format(obs_label))

  #---------------------------------------------------------------
  # Get min obs bin
  #---------------------------------------------------------------
  def obs_min_bins(self, obs_label):

    return getattr(self, 'obs_min_bins_{}'.format(obs_label))

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

  analysis = RunAnalysis(config_file = args.configFile)
  analysis.run_analysis()
