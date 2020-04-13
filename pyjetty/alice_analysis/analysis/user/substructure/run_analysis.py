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

  - Implement a user analysis class inheriting from this one, such as in james/run_analysis_james.py
    You should implement the following functions:
      - plot_single_result()
      - plot_all_results()
      - plot_performance() (optional)
    See the example for details.

  - You also should modify a few observable-specific functions at the top of substructure/analysis_utils_obs.py

Then, you should run your analysis with:

  python run_analysis_user.py /path/to/my/config.yaml


Author: James Mulligan (james.mulligan@berkeley.edu)
'''

import sys
import os
import argparse
import itertools
from array import *
import numpy
import math
import ROOT
import yaml
import shutil

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

    # Get the sub-configs to unfold
    self.jetR_list = config['jetR']
    self.obs_config_dict = config[self.observable]
    self.obs_subconfig_list = [name for name in list(self.obs_config_dict.keys()) if 'config' in name ]
    self.grooming_settings = self.utils.grooming_settings(self.obs_config_dict)
    self.obs_settings = self.utils.obs_settings(self.observable, self.obs_config_dict, self.obs_subconfig_list)
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
      obs_label = self.utils.obs_label(self.obs_settings[i], self.grooming_settings[i])

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

    # List of systematic variations to perform
    self.systematics_list = config['systematics_list']

    # Load paths to processing output, to be unfolded
    self.main_data = config['main_data']
    self.main_response = config['main_response']

    if 'trkeff' in self.systematics_list:
      self.trkeff_response = config['trkeff_response']

    if 'prior1' in self.systematics_list:
      self.prior1_variation_parameter = config['prior1_variation_parameter']
    if 'prior2' in self.systematics_list:
      self.prior2_variation_parameter = config['prior2_variation_parameter']

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

    if self.do_systematics:
      output_dir_systematics = self.create_output_subdir(output_dir, 'systematics')
      sys_root_filename = os.path.join(output_dir_systematics, 'fSystematics.root')
      fSystematics = ROOT.TFile(sys_root_filename, 'RECREATE')

    if self.do_plot_final_result:
      output_dir_final = self.create_output_subdir(output_dir, 'final_results')
      final_result_root_file = os.path.join(output_dir_final, 'fFinalResults.root')
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
      for i, _ in enumerate(self.obs_subconfig_list):

        obs_setting = self.obs_settings[i]
        grooming_setting = self.grooming_settings[i]
        obs_label = self.utils.obs_label(obs_setting, grooming_setting)

        # Compute systematics and attach to main results
        if self.do_systematics:
            self.compute_systematics(jetR, obs_label, obs_setting, grooming_setting)

        # Plot result for each individial subconfiguration
        if self.do_plot_final_result:
          self.plot_single_result(jetR, obs_label, obs_setting, grooming_setting) # You must implement this

      # Plot final results for all subconfigs
      if self.do_plot_final_result:
        self.plot_all_results(jetR) # You must implement this

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

      if systematic == 'trkeff':
        response = self.trkeff_response
      elif systematic == 'prior1':
        prior_variation_parameter = self.prior1_variation_parameter
        shutil.copyfile(main_response_location, os.path.join(output_dir, 'response.root'))
        rebin_response = False
      elif systematic == 'prior2':
        prior_variation_parameter = self.prior2_variation_parameter
        shutil.copyfile(main_response_location, os.path.join(output_dir, 'response.root'))
        rebin_response = False
      elif systematic == 'truncation':
        truncation = True
      elif systematic == 'binning':
        binning = True

      analysis = roounfold_obs.Roounfold_Obs(self.observable, data, response, self.config_file,
                                             output_dir, self.file_format,
                                             rebin_response=rebin_response,
                                             prior_variation_parameter=prior_variation_parameter,
                                             truncation=truncation, binning=binning)
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
  def compute_systematics(self, jetR, obs_label, obs_setting, grooming_setting):

    print( "Compute systematics for %s: R = %s, subobs = %s, grooming = %s" % 
           (self.observable, jetR, obs_label, grooming_setting) )

    # Select final regularization parameter
    if self.use_max_reg_param:
      reg_param_final = self.max_reg_param
      min_reg_param = 3
    else:
      reg_param_final = self.utils.get_reg_param(
        self.obs_settings, self.grooming_settings, self.obs_subconfig_list,
        self.obs_config_dict, obs_label, jetR)
      min_reg_param = reg_param_final

    # First, calculate systematics for all possible reg params
    for reg_param in range(min_reg_param, reg_param_final + 1):

      # Load 2D unfolded results from their files into attributes
      self.load_2D_observables(jetR, obs_label, obs_setting, grooming_setting, reg_param)

      # Loop through pt slices, and compute systematics for each 1D observable distribution
      for i in range(0, len(self.pt_bins_reported) - 1):
        min_pt_truth = self.pt_bins_reported[i]
        max_pt_truth = self.pt_bins_reported[i+1]

        # Load 1D unfolded results for each pt slice into attributes
        self.load_1D_observables(jetR, obs_label, obs_setting, grooming_setting, reg_param,
                                 min_pt_truth, max_pt_truth)

        # Compute systematics of the 1D distributions for each pt slice
        self.compute_obs_systematic(jetR, obs_label, obs_setting, grooming_setting, reg_param,
                                    min_pt_truth, max_pt_truth, final=False)

    # Now determine the best possible reg param by minimizing uncertainty
    for i in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[i]
      max_pt_truth = self.pt_bins_reported[i+1]
      if self.use_max_reg_param:
        reg_param_final = self.determine_reg_param_final(jetR, obs_label, obs_setting, grooming_setting,
                                                         min_pt_truth, max_pt_truth)
      print( "Optimal regularization parameter for pT=%i-%i determined to be %s = %i." % 
             (min_pt_truth, max_pt_truth, self.reg_param_name, reg_param_final) )
      self.compute_obs_systematic(jetR, obs_label, obs_setting, grooming_setting, reg_param_final,
                                  min_pt_truth, max_pt_truth, final=True)

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

        # Get tagging rate histogram, and store as an attribute
        if grooming_setting:
          name = 'hTaggingFractions_R{}_{}'.format(jetR, obs_label)
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
                          reg_param, min_pt_truth, max_pt_truth):

    # Get all other systematic variations, and store as attributes
    for systematic in self.systematics_list:

      name2D = 'hUnfolded_{}_R{}_{}_{}{}'.format(self.observable, jetR, obs_label,
                                                 reg_param, systematic)
      name1D = 'h{}_{}_R{}_{}_n{}_{}-{}'.format(systematic, self.observable, jetR,
                                            obs_label, reg_param, min_pt_truth, max_pt_truth)
      self.get_obs_distribution(jetR, obs_label, name2D, name1D, grooming_setting,
                                min_pt_truth, max_pt_truth)

      if systematic == 'main':
        # Get regularization parameter variations, and store as attributes
        name2D = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, reg_param+2)
        name1D = 'hRegParam1_{}_R{}_{}_n{}_{}-{}'.format(self.observable, jetR, obs_label,
                                                         reg_param, min_pt_truth, max_pt_truth)
        hRegParam1 = self.get_obs_distribution(jetR, obs_label, name2D, name1D, grooming_setting,
                                               min_pt_truth, max_pt_truth)

        name2D = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, reg_param-2)
        name1D = 'hRegParam2_{}_R{}_{}_n{}_{}-{}'.format(self.observable, jetR, obs_label,
                                                         reg_param, min_pt_truth, max_pt_truth)
        hRegParam2 = self.get_obs_distribution(jetR, obs_label, name2D, name1D, grooming_setting,
                                               min_pt_truth, max_pt_truth)

  #----------------------------------------------------------------------
  # Compute systematics
  #----------------------------------------------------------------------
  def compute_obs_systematic(self, jetR, obs_label, obs_setting, grooming_setting,
                             reg_param, min_pt_truth, max_pt_truth, final=False):

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
    for systematic in self.systematics_list:

      if final:
        # Save plot of each systematic to directory
        if systematic != 'main':
          name = 'hSystematic_{}_{}_R{}_{}_n{}_{}-{}'.format(self.observable, systematic, jetR, obs_label,
                                                             reg_param, min_pt_truth, max_pt_truth)
        else:  # systematic == 'main':
          name = 'hSystematic_{}_RegParam_R{}_{}_n{}_{}-{}'.format(self.observable, jetR, obs_label,
                                                                   reg_param, min_pt_truth, max_pt_truth)
        h_systematic_ratio = getattr(self, name)

        output_dir = getattr(self, 'output_dir_systematics')
        name = 'hSystematic_{}_R{}_{}_{}-{}{}'.format(systematic, self.utils.remove_periods(jetR), obs_label,
                                                      int(min_pt_truth), int(max_pt_truth), self.file_format)
        outputFilename = os.path.join(output_dir, name)
        self.utils.plot_hist(h_systematic_ratio, outputFilename, 'P E')

      elif systematic != 'main':
        # Get systematic variation and plot ratio 

        name = 'h{}_{}_R{}_{}_n{}_{}-{}'.format(systematic, self.observable, jetR, obs_label,
                                                reg_param, min_pt_truth, max_pt_truth)
        h_systematic = getattr(self, name)

        name_ratio = 'hSystematic_{}_{}_R{}_{}_n{}_{}-{}'.format(self.observable, systematic, jetR, obs_label,
                                                                 reg_param, min_pt_truth, max_pt_truth)
        h_systematic_ratio = hMain.Clone()
        h_systematic_ratio.SetName(name_ratio)
        h_systematic_ratio.Divide(h_systematic)

        print("Printing systematic variation ratio for", systematic) 
        for i in range(1, h_systematic.GetNbinsX()+1):
          print("main:", hMain.GetBinContent(i), "-- sys:", h_systematic.GetBinContent(i),
                "-- ratio:", h_systematic_ratio.GetBinContent(i))

        self.change_to_per(h_systematic_ratio)
        setattr(self, name_ratio, h_systematic_ratio)

      else:  # systematic == 'main':
        # Do the two reg param variations & average them

        name = 'hRegParam1_{}_R{}_{}_n{}_{}-{}'.format(self.observable, jetR, obs_label,
                                                       reg_param, min_pt_truth, max_pt_truth)
        hSystematic_RegParam1 = getattr(self, name)
        name_ratio = 'hSystematic_{}_RegParam1_R{}_{}_n{}_{}-{}'.format(self.observable, jetR, obs_label,
                                                                        reg_param, min_pt_truth, max_pt_truth)
        hSystematic_RegParam1_ratio = hMain.Clone()
        hSystematic_RegParam1_ratio.SetName(name_ratio)
        hSystematic_RegParam1_ratio.Divide(hSystematic_RegParam1)
        self.change_to_per(hSystematic_RegParam1_ratio)
        setattr(self, name_ratio, hSystematic_RegParam1_ratio)

        name = 'hRegParam2_{}_R{}_{}_n{}_{}-{}'.format(self.observable, jetR, obs_label, reg_param,
                                                       min_pt_truth, max_pt_truth)
        hSystematic_RegParam2 = getattr(self, name)
        name_ratio = 'hSystematic_{}_RegParam2_R{}_{}_n{}_{}-{}'.format(self.observable, jetR, obs_label,
                                                                        reg_param, min_pt_truth, max_pt_truth)
        hSystematic_RegParam2_ratio = hMain.Clone()
        hSystematic_RegParam2_ratio.SetName(name_ratio)
        hSystematic_RegParam2_ratio.Divide(hSystematic_RegParam2)
        self.change_to_per(hSystematic_RegParam2_ratio)
        setattr(self, name_ratio, hSystematic_RegParam2_ratio)

        name = 'hSystematic_{}_RegParam_R{}_{}_n{}_{}-{}'.format(self.observable, jetR, obs_label,
                                                                 reg_param, min_pt_truth, max_pt_truth)
        h_systematic_ratio = self.build_average(hSystematic_RegParam1_ratio, hSystematic_RegParam2_ratio)
        setattr(self, name, h_systematic_ratio)

      h_list.append(h_systematic_ratio)

    if not final:
      # Add uncertainties in quadrature
      hSystematic_Total = self.add_in_quadrature(h_list)
      name = 'hSystematic_Total_R{}_{}_n{}_{}-{}'.format(self.utils.remove_periods(jetR), obs_label,
                                                         reg_param, int(min_pt_truth), int(max_pt_truth))
      setattr(self, name, hSystematic_Total)
      # Attach total systematic to main result, and save as an attribute
      name = 'hResult_{}_systotal_R{}_{}_n{}_{}-{}'.format(self.observable, jetR, obs_label, reg_param,
                                                           int(min_pt_truth), int(max_pt_truth))
      hResult_sys = hMain.Clone()
      hResult_sys.SetName(name)
      hResult_sys.SetDirectory(0)
      self.AttachErrToHist(hResult_sys, hSystematic_Total)
      setattr(self, name, hResult_sys)

    else:
      # Save main result under name without reg param
      name = 'hResult_{}_systotal_R{}_{}_n{}_{}-{}'.format(self.observable, jetR, obs_label, reg_param,
                                                           int(min_pt_truth), int(max_pt_truth))
      hResult_sys = getattr(self, name)
      name = 'hResult_{}_systotal_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label,
                                                           int(min_pt_truth), int(max_pt_truth))
      setattr(self, name, hResult_sys)

      # Get total uncertainty plot from memory and save
      name = 'hSystematic_Total_R{}_{}_n{}_{}-{}'.format(self.utils.remove_periods(jetR), obs_label,
                                                         reg_param, int(min_pt_truth), int(max_pt_truth))
      hSystematic_Total = getattr(self, name)
      name = 'hSystematic_Total_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), obs_label,
                                                       int(min_pt_truth), int(max_pt_truth), self.file_format)
      outputFilename = os.path.join(self.output_dir_systematics, name)
      self.utils.plot_hist(hSystematic_Total, outputFilename, 'P E')

      # Plot systematic uncertainties, and write total systematic to a ROOT file
      self.plot_systematic_uncertainties(jetR, obs_label, obs_setting, grooming_setting,
                                         min_pt_truth, max_pt_truth, h_list, hSystematic_Total)

  #----------------------------------------------------------------------
  def get_obs_distribution(self, jetR, obs_label, name2D, name1D, grooming_setting,
                           min_pt_truth, max_pt_truth):

    h2D = getattr(self, name2D)
    h2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    h = h2D.ProjectionY() # Better to use ProjectionY('{}_py'.format(h2D.GetName()), 1, h2D.GetNbinsX()) ?
    h.SetName(name1D)
    h.SetDirectory(0)

    if grooming_setting:
        name = 'hTaggingFractions_R{}_{}'.format(jetR, obs_label)
        hTaggingFrac = getattr(self, name)
        x = (min_pt_truth + max_pt_truth)/2.
        fraction_tagged = hTaggingFrac.GetBinContent(hTaggingFrac.FindBin(x))
        setattr(self, '{}_fraction_tagged'.format(name1D), fraction_tagged)
    else:
        fraction_tagged = 1.

    n_jets_tagged = h.Integral(1, h.GetNbinsX())
    n_jets_inclusive = n_jets_tagged/fraction_tagged
    h.Scale(1./n_jets_inclusive, 'width')

    setattr(self, name1D, h)

    return h

  #----------------------------------------------------------------------
  def plot_systematic_uncertainties(self, jetR, obs_label, obs_setting, grooming_setting,
                                    min_pt_truth, max_pt_truth, h_list, h_total):

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

    n_bins_truth = self.n_bins_truth(obs_label)
    truth_bin_array = self.truth_bin_array(obs_label)

    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_bins_truth, truth_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle( getattr(self, 'xtitle') )
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle('Systematic uncertainty (%)')
    myBlankHisto.SetMaximum(1.5*h_total.GetMaximum())
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")

    leg = ROOT.TLegend(0.67,0.65,0.8,0.92)
    self.utils.setup_legend(leg,0.04)

    for i, h in enumerate(h_list):
      if h:
        if i == 0:
          h.SetMarkerStyle(20)
        if i == 1:
          h.SetMarkerStyle(21)
        if i == 2:
          h.SetMarkerStyle(22)
        if i == 3:
          h.SetMarkerStyle(23)
        if i == 4:
          h.SetMarkerStyle(33)
        if i == 5:
          h.SetMarkerStyle(34)

        h.SetMarkerSize(1.5)
        h.SetMarkerColor(600-5+i)
        h.SetLineStyle(1)
        h.SetLineWidth(2)
        h.SetLineColor(600-5+i)

        h.DrawCopy('P X0 same')

        leg.AddEntry(h, self.systematics_list[i], 'Pe')

    h_total.SetLineStyle(1)
    h_total.SetLineColor(1)
    h_total.SetLineWidth(2)
    h_total.DrawCopy('same hist')
    leg.AddEntry(h_total, 'Total', 'l')

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
      text_latex.DrawLatex(0.3, 0.64, text)

    output_dir = getattr(self, 'output_dir_systematics')
    outputFilename = os.path.join(output_dir,
                                  'hSystematics_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR),
                                                                       obs_label, int(min_pt_truth),
                                                                       int(max_pt_truth), self.file_format))
    c.SaveAs(outputFilename)
    c.Close()

    sys_root_filename = os.path.join(output_dir, 'fSystematics.root')
    fSystematics = ROOT.TFile(sys_root_filename, 'UPDATE')
    h_total.Write()
    fSystematics.Close()

  #----------------------------------------------------------------------
  # Add a list of (identically-binned) histograms in quadrature, bin-by-bin
  #----------------------------------------------------------------------
  def add_in_quadrature(self, h_list):

    h_new = h_list[0].Clone()
    h_new.SetName('{}_new'.format(h_list[0].GetName()))

    for i in range(1, h_new.GetNbinsX()+1):

      values_i = [h.GetBinContent(i) for h in h_list]

      new_value_squared = 0.
      for value_i in values_i:
        new_value_squared += value_i*value_i
      new_value = math.sqrt(new_value_squared)
      h_new.SetBinContent(i, new_value)

    return h_new

  #----------------------------------------------------------------------
  def determine_reg_param_final(self, jetR, obs_label, obs_setting, grooming_setting,
                                min_pt_truth, max_pt_truth):

    # Obtain the total systematic uncertainty histograms from memory
    min_reg_param = 3

    # Select final regularization parameter
    if self.use_max_reg_param:
      reg_param_final = self.max_reg_param
    else:
      reg_param_final = self.utils.get_reg_param(
        self.obs_settings, self.grooming_settings, self.obs_subconfig_list,
        self.obs_config_dict, obs_label, jetR)

    reg_params = range(min_reg_param, reg_param_final + 1)
    names = ["hSystematic_Total_R{}_{}_n{}_{}-{}".format(self.utils.remove_periods(jetR), obs_label,
                                                         reg_param, int(min_pt_truth), int(max_pt_truth))
             for reg_param in reg_params]
    hSys_list = [getattr(self, name) for name in names]

    # Obtain the statistical uncertainty histograms from disk
    fdir = os.path.join(self.output_dir_main, "Unfolded_stat_uncert",
                        "fResult_name_R%s_%s.root" % (jetR, obs_setting))
    f = ROOT.TFile.Open(fdir, "READ")
    names = ["hUnfolded_%s_stat_uncert_R%s_%s_n%i_Pt%i-%i" % (self.observable, self.utils.remove_periods(jetR),
                                                              obs_setting, reg_param, min_pt_truth, max_pt_truth) 
             for reg_param in reg_params]
    hStat_list = [f.Get(name) for name in names]

    # Add statistical and systematic uncertainties in quadrature
    hUncert_list = [self.add_in_quadrature([hSys, hStat]) for (hSys, hStat) in zip(hSys_list, hStat_list)]

    # Sum total uncertainties per observable bin and return lowest value
    hUncert_sum_list = [sum([h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)]) for h in hUncert_list]
    index_min = min(range(len(hUncert_sum_list)), key=hUncert_sum_list.__getitem__)
    return index_min + min_reg_param

    f.Close()

  #----------------------------------------------------------------------
  def change_to_per(self, h):

    for bin in range(0, h.GetNbinsX()):
      content = h.GetBinContent(bin)
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

    n_bins_truth = getattr(self, 'n_obs_bins_truth_{}'.format(obs_label))
    return n_bins_truth

  #---------------------------------------------------------------
  # Get truth_bin_array
  #---------------------------------------------------------------
  def truth_bin_array(self, obs_label):

    truth_bin_array = getattr(self, 'truth_obs_bin_array_{}'.format(obs_label))
    return truth_bin_array

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
