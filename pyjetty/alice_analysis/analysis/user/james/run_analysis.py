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

from pyjetty.alice_analysis.analysis.base import common_base
import analysis_utils_obs
import roounfold_obs

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

    print(self)
  
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    
    # Set list of observables
    self.observable = config['analysis_observable']
    
    # Set which analysis steps to perform
    self.do_unfolding = config['do_unfolding']
    self.do_systematics = config['do_systematics']
    self.do_plot_final_result = config['do_plot_final_result']
    self.force_rebin_response=config['force_rebin']
    
    # Get the sub-configs to unfold
    self.jetR_list = config['jetR']
    self.obs_config_dict = config[self.observable]
    self.obs_subconfig_list = [name for name in list(self.obs_config_dict.keys()) if 'config' in name ]
    self.sd_settings = self.utils.sd_settings(self.obs_config_dict)
    self.obs_settings = self.utils.obs_settings(self.observable, self.obs_config_dict, self.obs_subconfig_list)
    setattr(self, 'xtitle', self.obs_config_dict['common_settings']['xtitle'])
    setattr(self, 'ytitle', self.obs_config_dict['common_settings']['ytitle'])
    setattr(self, 'xmin', self.obs_config_dict['common_settings']['xmin'])
    setattr(self, 'xmax', self.obs_config_dict['common_settings']['xmax'])
    setattr(self, 'ymax', self.obs_config_dict['common_settings']['ymax'])
    setattr(self, 'ymin_ratio', self.obs_config_dict['common_settings']['ymin_ratio'])
    setattr(self, 'ymax_ratio', self.obs_config_dict['common_settings']['ymax_ratio'])

    # Retrieve histogram binnings for each observable setting
    for i, _ in enumerate(self.obs_subconfig_list):

      config_name = self.obs_subconfig_list[i]
      obs_label = self.utils.obs_label(self.obs_settings[i], self.sd_settings[i])
      
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

    # Theory comparisons
    if 'fPythia' in config:
      self.fPythia_name = config['fPythia']
    if 'fNLL' in config:
      self.fNLL = config['fNLL']
    if 'fNPcorrection_numerator' in config and 'fNPcorrection_denominator' in config:
      self.fNPcorrection_numerator = config['fNPcorrection_numerator']
      self.fNPcorrection_denominator = config['fNPcorrection_denominator']

  #---------------------------------------------------------------
  # Create a set of output directories for a given observable
  #---------------------------------------------------------------
  def create_output_dirs(self):
    
    output_dir = os.path.join(self.output_dir, self.observable)
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

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
        sd_setting = self.sd_settings[i]
        obs_label = self.utils.obs_label(obs_setting, sd_setting)
        reg_param_final = self.utils.get_reg_param(self.obs_settings, self.sd_settings, self.obs_subconfig_list, self.obs_config_dict, obs_label, jetR)

        if self.do_systematics:
          self.compute_systematics(jetR, obs_label, obs_setting, sd_setting, reg_param_final)

        if self.do_plot_final_result:
          self.plot_final_result(jetR, obs_label, obs_setting, sd_setting)
          
          if self.observable == 'theta_g' or self.observable == 'zg':
              self.get_nll_tgraph(jetR, obs_label, obs_setting, sd_setting, 20., 40.)
              self.get_nll_tgraph(jetR, obs_label, obs_setting, sd_setting, 40., 60.)
              self.get_nll_tgraph(jetR, obs_label, obs_setting, sd_setting, 60., 80.)

      # Plot final results for all subconfigs
      if self.do_plot_final_result:
        
        if self.observable == 'theta_g' or self.observable == 'zg':
          self.plot_final_result_overlay(jetR)

          #self.plot_NPcorrection(jetR)
    
  #----------------------------------------------------------------------
  def perform_unfolding(self):
    print('Perform unfolding for all systematic variations: {} ...'.format(self.observable))
    
    for systematic in self.systematics_list:
    
      output_dir = getattr(self, 'output_dir_{}'.format(systematic))
      data = self.main_data
      response = self.main_response
      prior_variation_parameter = 0.
      truncation = False
      binning = False
    
      if systematic == 'trkeff':
        response = self.trkeff_response
      if systematic == 'prior1':
        prior_variation_parameter = self.prior1_variation_parameter
      if systematic == 'prior2':
        prior_variation_parameter = self.prior2_variation_parameter
      if systematic == 'truncation':
        truncation = True
      if systematic == 'binning':
        binning = True
      
      analysis = roounfold_obs.Roounfold_Obs(self.observable, data, response, self.config_file, output_dir, self.file_format, rebin_response=self.check_rebin_response(output_dir), prior_variation_parameter=prior_variation_parameter, truncation=truncation, binning=binning)
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
  def compute_systematics(self, jetR, obs_label, obs_setting, sd_setting, reg_param_final):
    print('Compute systematics for {}: R = {}, subobs = {} , SD = {}...'.format(self.observable, jetR, obs_label, sd_setting))

    # Load 2D unfolded results from their files into attributes
    self.load_2D_observables(jetR, obs_label, obs_setting, sd_setting, reg_param_final)
    
    # Loop through pt slices, and compute systematics for each 1D observable distribution
    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(obs_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(obs_label))
    
    for bin in range(1, n_pt_bins_truth-3):
      min_pt_truth = truth_pt_bin_array[bin]
      max_pt_truth = truth_pt_bin_array[bin+1]
      
      # Load 1D unfolded results for each pt slice into attributes
      self.load_1D_observables(jetR, obs_label, obs_setting, sd_setting, reg_param_final, min_pt_truth, max_pt_truth)
      
      # Compute systematics of the 1D distributions for each pt slice
      self.compute_obs_systematic(jetR, obs_label, obs_setting, sd_setting, reg_param_final, min_pt_truth, max_pt_truth)

  #----------------------------------------------------------------------
  def load_2D_observables(self, jetR, obs_label, obs_setting, sd_setting, reg_param_final):
    
    # Get all other systematic variations, and store as attributes
    for systematic in self.systematics_list:
    
      output_dir = getattr(self, 'output_dir_{}'.format(systematic))
      path = os.path.join(output_dir, 'fResult_R{}_{}.root'.format(jetR, obs_label))
      f = ROOT.TFile(path, 'READ')
      name = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, reg_param_final)
      self.retrieve_histo_and_set_attribute(name, f, systematic)
      
      if systematic == 'main':
        # Get regularization parameter variations, and store as attributes
        name = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, reg_param_final+2)
        self.retrieve_histo_and_set_attribute(name, f)
        name = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, reg_param_final-2)
        self.retrieve_histo_and_set_attribute(name, f)

    # Get tagging rate histogram, and store as an attribute
    if sd_setting:
      name = 'hTaggingFractions_R{}_{}'.format(jetR, obs_label)
      self.retrieve_histo_and_set_attribute(name, fMain)
      
  #----------------------------------------------------------------------
  def retrieve_histo_and_set_attribute(self, name, file, suffix = ''):
  
    h = file.Get(name)
    h.SetDirectory(0)
    setattr(self, '{}{}'.format(name, suffix), h)

  #----------------------------------------------------------------------
  # Get 1D histograms
  # Normalize by integral, i.e. N_jets,inclusive in this pt-bin
  #----------------------------------------------------------------------
  def load_1D_observables(self, jetR, obs_label, obs_setting, sd_setting, reg_param_final, min_pt_truth, max_pt_truth):
  
    # Get all other systematic variations, and store as attributes
    for systematic in self.systematics_list:
  
      name2D = 'hUnfolded_{}_R{}_{}_{}{}'.format(self.observable, jetR, obs_label, reg_param_final, systematic)
      name1D = 'h{}_{}_R{}_{}_{}-{}'.format(systematic, self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
      self.get_obs_distribution(jetR, obs_label, name2D, name1D, sd_setting, min_pt_truth, max_pt_truth)

      if systematic == 'main':
        # Get regularization parameter variations, and store as attributes
        name2D = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, reg_param_final+2)
        name1D = 'hRegParam1_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
        hRegParam1 = self.get_obs_distribution(jetR, obs_label, name2D, name1D, sd_setting, min_pt_truth, max_pt_truth)
        
        name2D = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, reg_param_final-2)
        name1D = 'hRegParam2_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
        hRegParam2 = self.get_obs_distribution(jetR, obs_label, name2D, name1D, sd_setting, min_pt_truth, max_pt_truth)

  #----------------------------------------------------------------------
  # Compute systematics
  #----------------------------------------------------------------------
  def compute_obs_systematic(self, jetR, obs_label, obs_setting, sd_setting, reg_param_final, min_pt_truth, max_pt_truth):
    
    # Get main result
    name = 'h{}_{}_R{}_{}_{}-{}'.format('main', self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    hMain = getattr(self, name)
    
    # Loop through all systematic variations, and take ratio to main result
    h_list = []
    for systematic in self.systematics_list:
    
      if systematic != 'main':
        # Do all variations

        name = 'h{}_{}_R{}_{}_{}-{}'.format(systematic, self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
        name_ratio = 'hSystematic_{}_{}_R{}_{}_{}-{}'.format(self.observable, systematic, jetR, obs_label, min_pt_truth, max_pt_truth)
        h_systematic = getattr(self, name)
        h_systematic_ratio = hMain.Clone()
        h_systematic_ratio.SetName(name_ratio)
        h_systematic_ratio.Divide(h_systematic)
        self.change_to_per(h_systematic_ratio)
        setattr(self, name_ratio, h_systematic_ratio)
          
      if systematic == 'main':
        # Do reg param variations
      
        name = 'hRegParam1_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
        name_ratio = 'hSystematic_{}_RegParam1_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
        hSystematic_RegParam1 = getattr(self, name)
        hSystematic_RegParam1_ratio = hMain.Clone()
        hSystematic_RegParam1_ratio.SetName(name)
        hSystematic_RegParam1_ratio.Divide(hSystematic_RegParam1)
        self.change_to_per(hSystematic_RegParam1_ratio)
        setattr(self, name, hSystematic_RegParam1_ratio)

        name = 'hRegParam2_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
        name_ratio = 'hSystematic_{}_RegParam2_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
        hSystematic_RegParam2 = getattr(self, name)
        hSystematic_RegParam2_ratio = hMain.Clone()
        hSystematic_RegParam2_ratio.SetName(name)
        hSystematic_RegParam2_ratio.Divide(hSystematic_RegParam2)
        self.change_to_per(hSystematic_RegParam2_ratio)
        setattr(self, name, hSystematic_RegParam2_ratio)
    
        name = 'hSystematic_{}_RegParam_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
        h_systematic_ratio = self.build_average(hSystematic_RegParam1, hSystematic_RegParam2)
        setattr(self, name, h_systematic_ratio)
        
      if self.debug_level > 0
        output_dir = getattr(self, 'output_dir_systematics')
        outputFilename = os.path.join(output_dir, 'hSystematic_{}_R{}_{}_{}-{}{}'.format(systematic, self.utils.remove_periods(jetR), obs_label, int(min_pt_truth), int(max_pt_truth), self.file_format))
        self.utils.plot_hist(h_systematic_ratio, outputFilename, 'P E')
        
      h_list.append(h_systematic_ratio)
    
    # Add uncertainties in quadrature
    hSystematic_Total = self.add_in_quadrature(h_list)
    output_dir = getattr(self, 'output_dir_systematics')
    outputFilename = os.path.join(output_dir, 'hSystematic_Total_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), obs_label, int(min_pt_truth), int(max_pt_truth), self.file_format))
    #self.utils.plot_hist(hSystematic_Total, outputFilename, 'P E')

    # Attach total systematic to main result, and save as an attribute
    name = 'hResult_{}_systotal_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    hResult_sys = hMain.Clone()
    hResult_sys.SetName(name)
    hResult_sys.SetDirectory(0)
    self.AttachErrToHist(hResult_sys, hSystematic_Total)
    setattr(self, name, hResult_sys)
      
    # Plot systematic uncertainties, and write total systematic to a ROOT file
    self.plot_systematic_uncertainties(jetR, obs_label, obs_setting, sd_setting, min_pt_truth, max_pt_truth, h_list, hSystematic_Total)
  
  #----------------------------------------------------------------------
  def get_obs_distribution(self, jetR, obs_label, name2D, name1D, sd_setting, min_pt_truth, max_pt_truth):
  
    h2D = getattr(self, name2D)
    h2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    h = h2D.ProjectionY() # Better to use ProjectionY('{}_py'.format(h2D.GetName()), 1, h2D.GetNbinsX()) ?
    h.SetName(name1D)
    h.SetDirectory(0)
    
    if sd_setting:
        name = 'hTaggingFractions_R{}_{}'.format(jetR, obs_label)
        hTaggingFrac = getattr(self, name)
        x = (min_pt_truth + max_pt_truth)/2.
        fraction_tagged =  hTaggingFrac.GetBinContent(hTaggingFrac.FindBin(x))
        setattr(self, '{}_fraction_tagged'.format(name1D), fraction_tagged)
    else:
        fraction_tagged = 1.

    n_jets_tagged = h.Integral(1, h.GetNbinsX())
    n_jets_inclusive = n_jets_tagged/fraction_tagged
    h.Scale(1./n_jets_inclusive, 'width')
  
    setattr(self, name1D, h)

    return h
  
  #----------------------------------------------------------------------
  def plot_systematic_uncertainties(self, jetR, obs_label, obs_setting, sd_setting, min_pt_truth, max_pt_truth, h_list, h_total):
  
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
    
    if sd_setting:
      text = self.utils.formatted_sd_label(sd_setting)
      text_latex.DrawLatex(0.3, 0.64, text)
      
    output_dir = getattr(self, 'output_dir_systematics')
    outputFilename = os.path.join(output_dir, 'hSystematics_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), obs_label, int(min_pt_truth), int(max_pt_truth), self.file_format))
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
  def plot_final_result(self, jetR, obs_label, obs_setting, sd_setting):
    print('Plot final results for {}: R = {}, {} ...'.format(self.observable, jetR, obs_label))

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    # Loop through pt slices, and compute systematics for each 1D theta_g distribution
    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(obs_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(obs_label))
    
    for bin in range(1, n_pt_bins_truth-3):
      min_pt_truth = truth_pt_bin_array[bin]
      max_pt_truth = truth_pt_bin_array[bin+1]
      
      #self.get_NPcorrection(self.observable, jetR, obs_label, obs_setting, min_pt_truth, max_pt_truth)
      
      self.plot_observable(self.observable, jetR, obs_label, obs_setting, sd_setting, min_pt_truth, max_pt_truth, plot_pythia=False)
      #self.plot_observable(self.observable, jetR, obs_label, obs_setting, sd_setting, min_pt_truth, max_pt_truth, plot_pythia=True)

  #----------------------------------------------------------------------
  def get_NPcorrection(self, observable, jetR, obs_label, obs_setting, min_pt_truth, max_pt_truth):
  
    fNumerator = ROOT.TFile(self.fNPcorrection_numerator, 'READ')
    fDenominator = ROOT.TFile(self.fNPcorrection_denominator, 'READ')
    
    n_bins_truth = self.n_bins_truth(obs_label)
    truth_bin_array = self.truth_bin_array(obs_label)
    
    hname = 'histogram_h_{}_B{}_{}-{}'.format(observable, obs_setting[1], int(min_pt_truth), int(max_pt_truth))
    
    hNumerator = fNumerator.Get(hname)
    hNumerator.SetDirectory(0)
    n_jets_inclusive = hNumerator.Integral(0, hNumerator.GetNbinsX()+1)
    n_jets_tagged = hNumerator.Integral(hNumerator.FindBin(truth_bin_array[0]), hNumerator.GetNbinsX())
    fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
    hNumerator.Scale(1./n_jets_inclusive, 'width')
      
    hDenominator = fDenominator.Get(hname)
    hDenominator.SetDirectory(0)
    n_jets_inclusive = hDenominator.Integral(0, hDenominator.GetNbinsX()+1)
    n_jets_tagged = hDenominator.Integral(hDenominator.FindBin(truth_bin_array[0]), hDenominator.GetNbinsX())
    fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
    hDenominator.Scale(1./n_jets_inclusive, 'width')
        
    hNumerator.Divide(hDenominator)
    hNumerator.SetName('hNPcorrection_{}_{}_{}-{}'.format(observable, obs_label, min_pt_truth, max_pt_truth))
    setattr(self, 'hNPcorrection_{}_{}_{}-{}'.format(observable, obs_label, min_pt_truth, max_pt_truth), hNumerator)
      
    #output_dir = getattr(self, 'output_dir_final_results')
    #outputFilename = os.path.join(output_dir, 'hNPcorrection_{}_{}_{}-{}.pdf'.format(observable, obs_label, min_pt_truth, max_pt_truth))
    #self.utils.plot_hist(hNumerator, outputFilename)

  #----------------------------------------------------------------------
  def get_nll_tgraph(self, jetR, obs_label, obs_setting, sd_setting, min_pt_truth, max_pt_truth):

    n_bins_truth = self.n_bins_truth(obs_label)
    truth_bin_array = self.truth_bin_array(obs_label)
    
    if self.observable == 'theta_g':
      path_txt = '/Users/jamesmulligan/Analysis_theta_g/NLL/Rg_value/beta{}/{}_{}.dat'.format(sd_setting[1], int(min_pt_truth), int(max_pt_truth))
    if self.observable == 'zg':
      path_txt = '/Users/jamesmulligan/Analysis_theta_g/NLL/zg_value/beta{}/{}_{}.dat'.format(sd_setting[1], int(min_pt_truth), int(max_pt_truth))
    
    filename = open(path_txt, 'r')

    x_list = []
    center_list = []
    low_list = []
    up_list = []
    for line in filename.readlines():
      line = line.rstrip('\n')
      
      row = line.split(' ')
      x_list.append(float(row[0]))
      low_list.append(float(row[1]))
      center_list.append(float(row[2]))
      up_list.append(float(row[3]))

    #x = array('d', x_list)
    x = numpy.array(x_list)
    x_err = numpy.zeros(n_bins_truth)
    center = numpy.array(center_list)
    low = numpy.subtract(center, numpy.array(low_list))
    up = numpy.subtract(numpy.array(up_list), center)
    
    g = ROOT.TGraphAsymmErrors(n_bins_truth, x, center, x_err, x_err, low, up)
    g.SetName('tgraph_{}_{}_{}-{}'.format(self.observable, obs_label, min_pt_truth, max_pt_truth))
    setattr(self, 'tgraph_NLL_{}_{}_{}-{}'.format(self.observable, obs_label, min_pt_truth, max_pt_truth), g)

  #----------------------------------------------------------------------
  def plot_NPcorrection(self, observable, jetR):

    self.plot_NPcorrection_overlay(observable, jetR, 20., 40.)
    self.plot_NPcorrection_overlay(observable, jetR, 40., 60.)
    self.plot_NPcorrection_overlay(observable, jetR, 60., 80.)
    
  #----------------------------------------------------------------------
  def plot_NPcorrection_overlay(self, observable, jetR, min_pt_truth, max_pt_truth):
    
    name = 'cResult_NPcorrection_R{}_allpt_{}-{}'.format(jetR, min_pt_truth, max_pt_truth)
    c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()
    c.cd()
    pad1 = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    pad1.SetLeftMargin(0.2)
    pad1.SetTopMargin(0.07)
    pad1.SetRightMargin(0.04)
    pad1.SetBottomMargin(0.13)
    pad1.Draw()
    pad1.cd()

    xmin = getattr(self, 'xmin')
    xmax = getattr(self, 'xmax')
    ymax = getattr(self, 'ymax')
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, xmin, xmax)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle( getattr(self, 'xtitle') )
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle('Correction')
    myBlankHisto.SetMaximum(ymax)
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")
    
    line = ROOT.TLine(0,1,xmax,1)
    line.SetLineColor(920+2)
    line.SetLineStyle(2)
    line.Draw()
    
    pad1.cd()
    myLegend = ROOT.TLegend(0.76,0.65,0.9,0.85)
    self.utils.setup_legend(myLegend,0.035)
    
    # Retrieve histogram binnings for each observable setting
    for i, obs_setting in enumerate(self.obs_settings):
      
      obs_label = self.utils.obs_label(obs_setting, sd_setting)
      
      n_obs_bins_truth = self.n_bins_truth(obs_label)
      truth_bin_array = self.truth_bin_array(obs_label)

      if i == 0:
        marker = 20
        color = 600-6
      if i == 1:
        marker = 21
        color = 632-4
      if i == 2:
        marker = 33
        color = 416-2

      h = getattr(self, 'hNPcorrection_{}_{}_{}-{}'.format(observable, obs_label, min_pt_truth, max_pt_truth))
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(marker)
      h.SetMarkerColor(color)
      h.SetLineStyle(1)
      h.SetLineWidth(2)
      h.SetLineColor(color)
      h.DrawCopy('PE X0 same')
      
      if observable == 'theta_g' or observable == 'zg':
        myLegend.AddEntry(h, 'pp #beta={}'.format(obs_setting[1]), 'pe')

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE Simulation'
    #text_latex.DrawLatex(0.25, 0.87, text)
  
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.81, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'PYTHIA8 Monash2013'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.75, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    if observable == 'theta_g' or observable == 'zg':
      text = '#it{R} = ' + str(jetR) + '  |#it{#eta}_{jet}| < 0.5' + '  #it{z}_{cut} = ' + str(obs_setting[0])
      text_latex.SetTextSize(0.045)
      text_latex.DrawLatex(0.25, 0.69, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.63, text)
    
    myLegend.Draw()
    
    name = 'hNPcorrection_{}_R{}_{}-{}{}'.format(observable, self.utils.remove_periods(jetR), int(min_pt_truth), int(max_pt_truth), self.file_format)
    output_dir = getattr(self, 'output_dir_final_results')
    outputFilename = os.path.join(output_dir, name)
    c.SaveAs(outputFilename)
    c.Close()
    
  #----------------------------------------------------------------------
  def plot_final_result_overlay(self, jetR):
  
    # Plot overlay of different beta, for fixed pt bin

    # Plot PYTHIA
    self.plot_observable_overlay_beta(jetR, 20., 40., plot_pythia=True, plot_nll = False, plot_ratio = True)
    self.plot_observable_overlay_beta(jetR, 40., 60., plot_pythia=True, plot_nll = False, plot_ratio = True)
    self.plot_observable_overlay_beta(jetR, 60., 80., plot_pythia=True, plot_nll = False, plot_ratio = True)
  
    # Plot NLL
    self.plot_observable_overlay_beta(jetR, 20., 40., plot_pythia=False, plot_nll = True, plot_ratio = False)
    self.plot_observable_overlay_beta(jetR, 40., 60., plot_pythia=False, plot_nll = True, plot_ratio = False)
    self.plot_observable_overlay_beta(jetR, 60., 80., plot_pythia=False, plot_nll = True, plot_ratio = False)

  #----------------------------------------------------------------------
  def plot_observable_overlay_beta(self, jetR, min_pt_truth, max_pt_truth, plot_pythia=False, plot_nll=False, plot_ratio=False):
    
    name = 'cResult_overlay_R{}_allpt_{}-{}'.format(jetR, min_pt_truth, max_pt_truth)
    if plot_ratio:
      c = ROOT.TCanvas(name, name, 600, 650)
    else:
      c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()
    
    c.cd()
    if plot_ratio:
      pad1 = ROOT.TPad('myPad', 'The pad',0,0.3,1,1)
    else:
      pad1 = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    pad1.SetLeftMargin(0.2)
    pad1.SetTopMargin(0.07)
    pad1.SetRightMargin(0.04)
    pad1.SetBottomMargin(0.13)
    if plot_ratio:
      pad1.SetBottomMargin(0.)
    pad1.Draw()
    pad1.cd()
    
    xmin = getattr(self, 'xmin')
    xmax = getattr(self, 'xmax')
    ymax = getattr(self, 'ymax')
    ymin_ratio = getattr(self, 'ymin_ratio')
    ymax_ratio = getattr(self, 'ymax_ratio')
    xtitle = getattr(self, 'xtitle')
    ytitle = getattr(self, 'ytitle')

    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, xmin, xmax)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle(xtitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle(ytitle)
    myBlankHisto.SetMaximum(ymax)
    myBlankHisto.SetMinimum(0.)
    if plot_ratio:
      myBlankHisto.SetMinimum(2e-4) # Don't draw 0 on top panel
      myBlankHisto.GetYaxis().SetTitleSize(0.065)
      myBlankHisto.GetYaxis().SetTitleOffset(1.4)
      myBlankHisto.GetYaxis().SetLabelSize(0.06)
    myBlankHisto.Draw("E")

    # Plot ratio
    if plot_ratio:
      
      c.cd()
      pad2 = ROOT.TPad("pad2", "pad2", 0, 0.02, 1, 0.3)
      pad2.SetTopMargin(0)
      pad2.SetBottomMargin(0.4)
      pad2.SetLeftMargin(0.2)
      pad2.SetRightMargin(0.04)
      pad2.Draw()
      pad2.cd()
      
      myBlankHisto2 = myBlankHisto.Clone("myBlankHisto_C")
      myBlankHisto2.SetYTitle("#frac{Data}{PYTHIA}")
      myBlankHisto2.SetXTitle(xtitle)
      myBlankHisto2.GetXaxis().SetTitleSize(30)
      myBlankHisto2.GetXaxis().SetTitleFont(43)
      myBlankHisto2.GetXaxis().SetTitleOffset(4.)
      myBlankHisto2.GetXaxis().SetLabelFont(43)
      myBlankHisto2.GetXaxis().SetLabelSize(25)
      myBlankHisto2.GetYaxis().SetTitleSize(20)
      myBlankHisto2.GetYaxis().SetTitleFont(43)
      myBlankHisto2.GetYaxis().SetTitleOffset(2.2)
      myBlankHisto2.GetYaxis().SetLabelFont(43)
      myBlankHisto2.GetYaxis().SetLabelSize(25)
      myBlankHisto2.GetYaxis().SetNdivisions(505)
      myBlankHisto2.GetYaxis().SetRangeUser(ymin_ratio, ymax_ratio)
      myBlankHisto2.Draw()
    
      line = ROOT.TLine(0,1,xmax,1)
      line.SetLineColor(920+2)
      line.SetLineStyle(2)
      line.Draw()

    pad1.cd()
    myLegend = ROOT.TLegend(0.66,0.65,0.8,0.85)
    self.utils.setup_legend(myLegend,0.035)
      
      
    for i, _ in enumerate(self.obs_subconfig_list):

      obs_setting = self.obs_settings[i]
      sd_setting = self.sd_settings[i]
      obs_label = self.utils.obs_label(obs_setting, sd_setting)

      n_obs_bins_truth = self.n_bins_truth(obs_label)
      truth_bin_array = self.truth_bin_array(obs_label)

      if i == 0:
        marker = 20
        marker_pythia = marker+4
        color = 600-6
      if i == 1:
        marker = 21
        marker_pythia = marker+4
        color = 632-4
      if i == 2:
        marker = 33
        marker_pythia = 27
        color = 416-2

      pad1.cd()
      if plot_pythia:
        
        fPythia = ROOT.TFile(self.fPythia_name, 'READ')
        hname = 'histogram_h_{}_B{}_{}-{}'.format(self.observable, sd_setting[1], int(min_pt_truth), int(max_pt_truth))
        hPythia = fPythia.Get(hname)
        hPythia.SetDirectory(0)
        
        n_jets_inclusive = hPythia.Integral(0, hPythia.GetNbinsX()+1)
        n_jets_tagged = hPythia.Integral(hPythia.FindBin(truth_bin_array[0]), hPythia.GetNbinsX())
        fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
        hPythia.Scale(1./n_jets_inclusive, 'width')
        
        plot_errors = False
        if plot_errors:
          hPythia.SetMarkerSize(0)
          hPythia.SetMarkerStyle(0)
          hPythia.SetMarkerColor(color)
          hPythia.SetFillColor(color)
          hPythia.DrawCopy('E3 same')
        else:
          hPythia.SetLineColor(color)
          hPythia.SetLineColorAlpha(color, 0.5)
          hPythia.SetLineWidth(4)
          hPythia.DrawCopy('L hist same')

      if plot_nll:
        
        # Get parton-level prediction
        g = getattr(self, 'tgraph_NLL_{}_{}_{}-{}'.format(self.observable, obs_label, min_pt_truth, max_pt_truth))
        
        # Get correction
        apply_nll_correction = False
        if apply_nll_correction:
          h_correction = getattr(self, 'hNPcorrection_{}_{}_{}-{}'.format(self.observable, obs_label, min_pt_truth, max_pt_truth))
        
          # Apply correction
          self.utils.multiply_tgraph(g, h_correction)
        
        g.SetLineColor(color)
        g.SetLineColorAlpha(color, 0.5)
        g.SetLineWidth(4)
        g.SetFillColor(color)
        g.SetFillColorAlpha(color, 0.5)
        g.Draw('L3 same')
    
      h_sys = getattr(self, 'hResult_{}_systotal_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth))
      h_sys.SetLineColor(0)
      h_sys.SetFillColor(color)
      h_sys.SetFillColorAlpha(color, 0.3)
      h_sys.SetFillStyle(1001)
      h_sys.SetLineWidth(0)
      h_sys.DrawCopy('E2 same')
    
      name = 'hMain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
      fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
      h = getattr(self, name)
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(marker)
      h.SetMarkerColor(color)
      h.SetLineStyle(1)
      h.SetLineWidth(2)
      h.SetLineColor(color)
      h.DrawCopy('PE X0 same')
        
      myLegend.AddEntry(h, 'ALICE pp #beta={}'.format(sd_setting[1]), 'pe')

      if plot_ratio:
        pad2.cd()
        
        hRatioSys = h_sys.Clone()
        hRatioSys.SetName('{}_Ratio'.format(h_sys.GetName()))
        if plot_pythia:
          hRatioSys.Divide(hPythia)
          hRatioSys.SetLineColor(0)
          hRatioSys.SetFillColor(color)
          hRatioSys.SetFillColorAlpha(color, 0.3)
          hRatioSys.SetFillStyle(1001)
          hRatioSys.SetLineWidth(0)
          hRatioSys.DrawCopy('E2 same')
        elif plot_nll:
          gRatioSys = g.Clone()
          gRatioSys.SetName('{}_{}_Ratio'.format(obs_label, g.GetName()))
          self.utils.divide_tgraph(hRatioSys, gRatioSys, combine_errors=True)
          gRatioSys.SetLineColor(0)
          gRatioSys.SetFillColor(color)
          gRatioSys.SetFillColorAlpha(color, 0.3)
          gRatioSys.SetFillStyle(1001)
          gRatioSys.SetLineWidth(0)
          gRatioSys.Draw('L3 same')

        hRatioStat = h.Clone()
        hRatioStat.SetName('{}_Ratio'.format(h.GetName()))
        if plot_pythia:
          hRatioStat.Divide(hPythia)
        elif plot_nll:
          self.utils.divide_tgraph(hRatioStat, g, combine_errors=False)
        hRatioStat.SetMarkerSize(1.5)
        hRatioStat.SetMarkerStyle(marker)
        hRatioStat.SetMarkerColor(color)
        hRatioStat.SetLineStyle(1)
        hRatioStat.SetLineWidth(2)
        hRatioStat.SetLineColor(color)
        hRatioStat.DrawCopy('PE X0 same')
        
    pad1.cd()
    myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    if plot_pythia:
      myLegend.AddEntry(hPythia, 'PYTHIA8 Monash2013', 'l')
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE Preliminary'
    text_latex.DrawLatex(0.25, 0.87, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.81, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'Charged jets   anti-#it{k}_{T}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.75, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    if self.observable == 'theta_g' or self.observable == 'zg':
      text = '#it{R} = ' + str(jetR) + '  |#it{#eta}_{jet}| < 0.5' + '  #it{z}_{cut} = ' + str(sd_setting[0])
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.69, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.63, text)
    
    myLegend.Draw()

    name = 'h_{}_R{}_{}-{}{}'.format(self.observable, self.utils.remove_periods(jetR), int(min_pt_truth), int(max_pt_truth), self.file_format)
    if plot_pythia:
      name = 'h_{}_R{}_{}-{}_Pythia{}'.format(self.observable, self.utils.remove_periods(jetR), int(min_pt_truth), int(max_pt_truth), self.file_format)
    if plot_nll:
      name = 'h_{}_R{}_{}-{}_NLL{}'.format(self.observable, self.utils.remove_periods(jetR), int(min_pt_truth), int(max_pt_truth), self.file_format)
    output_dir = getattr(self, 'output_dir_final_results')
    outputFilename = os.path.join(output_dir, name)
    c.SaveAs(outputFilename)
    c.Close()
  
  #----------------------------------------------------------------------
  def plot_observable(self, observable, jetR, obs_label, obs_setting, sd_setting, min_pt_truth, max_pt_truth, plot_pythia=False):
    
    name = 'cResult_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt_truth, max_pt_truth)
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
    
    ymax = getattr(self, 'ymax')
    xtitle = getattr(self, 'xtitle')
    ytitle = getattr(self, 'ytitle')
    
    n_obs_bins_truth = self.n_bins_truth(obs_label)
    truth_bin_array = self.truth_bin_array(obs_label)
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_obs_bins_truth, truth_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle( xtitle )
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle(ytitle)
    myBlankHisto.SetMaximum(ymax)
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")

    if plot_pythia:
    
      plot_pythia_from_response = True
      plot_pythia_from_mateusz = False
      if plot_pythia_from_response:
        hPythia = self.get_pythia_from_response(observable, jetR, obs_label, min_pt_truth, max_pt_truth)
        n_jets_inclusive = hPythia.Integral(0, hPythia.GetNbinsX()+1)
        n_jets_tagged = hPythia.Integral(hPythia.FindBin(truth_bin_array[0]), hPythia.GetNbinsX())
        fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
        hPythia.Scale(1./n_jets_inclusive, 'width')
        hPythia.SetFillStyle(0)
        hPythia.SetMarkerSize(1.5)
        hPythia.SetMarkerStyle(21)
        hPythia.SetMarkerColor(2)
        hPythia.SetLineColor(2)
        hPythia.SetLineWidth(1)
        hPythia.Draw('E2 same')
    
      if plot_pythia_from_mateusz:
        fPythia_name = '/Users/jamesmulligan/Analysis_theta_g/Pythia_new/pythia.root'
        fPythia = ROOT.TFile(fPythia_name, 'READ')
        hname = 'histogram_h_{}_B{}_{}-{}'.format(observable, obs_setting[1], int(min_pt_truth), int(max_pt_truth))
        hPythia2 = fPythia.Get(hname)

        n_jets_inclusive2 = hPythia2.Integral(0, hPythia2.GetNbinsX()+1)
        n_jets_tagged2 = hPythia2.Integral(hPythia2.FindBin(truth_bin_array[0]), hPythia2.GetNbinsX())
        fraction_tagged_pythia =  n_jets_tagged2/n_jets_inclusive2
        hPythia2.Scale(1./n_jets_inclusive2, 'width')
        hPythia2.SetFillStyle(0)
        hPythia2.SetMarkerSize(1.5)
        hPythia2.SetMarkerStyle(21)
        hPythia2.SetMarkerColor(1)
        hPythia2.SetLineColor(1)
        hPythia2.SetLineWidth(1)
        hPythia2.Draw('E2 same')
    
    color = 600-6
    h_sys = getattr(self, 'hResult_{}_systotal_R{}_{}_{}-{}'.format(observable, jetR, obs_label, min_pt_truth, max_pt_truth))
    h_sys.SetLineColor(0)
    h_sys.SetFillColor(color)
    h_sys.SetFillColorAlpha(color, 0.3)
    h_sys.SetFillStyle(1001)
    h_sys.SetLineWidth(0)
    h_sys.DrawCopy('E2 same')
    
    name = 'hMain_{}_R{}_{}_{}-{}'.format(observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    if observable == 'theta_g' or observable == 'zg':
        fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
    h = getattr(self, name)
    h.SetMarkerSize(1.5)
    h.SetMarkerStyle(20)
    h.SetMarkerColor(color)
    h.SetLineStyle(1)
    h.SetLineWidth(2)
    h.SetLineColor(color)
    h.DrawCopy('PE X0 same')
  
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE Preliminary'
    text_latex.DrawLatex(0.57, 0.87, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.8, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'R = ' + str(jetR) + '   | #eta_{jet}| < 0.5'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.73, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.66, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    if observable == 'theta_g' or observable == 'zg':
      text = '#it{z}_{cut} = ' + str(sd_setting[0]) + '  #beta = ' + str(sd_setting[1])
    elif observable == 'subjet_z':
      text = 'R_{subjet} = ' + str(obs_setting)
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.59, text)
    
    if observable == 'theta_g' or observable == 'zg':
        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        text_latex.SetTextSize(0.04)
        text = '#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged
        text_latex.DrawLatex(0.57, 0.52, text)
    
    if plot_pythia:
      text_latex = ROOT.TLatex()
      text_latex.SetNDC()
      text_latex.SetTextSize(0.04)
      text = ('#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged) + (', #it{f}_{tagged}^{pythia} = %3.3f' % fraction_tagged_pythia)
      text_latex.DrawLatex(0.57, 0.52, text)

    myLegend = ROOT.TLegend(0.25,0.7,0.5,0.85)
    self.utils.setup_legend(myLegend,0.035)
    myLegend.AddEntry(h, 'ALICE pp', 'pe')
    myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    if plot_pythia:
      if plot_pythia_from_response:
        myLegend.AddEntry(hPythia, 'PYTHIA Monash2013', 'pe')
      if plot_pythia_from_mateusz:
        myLegend.AddEntry(hPythia2, 'PYTHIA Monash2013', 'pe')
    myLegend.Draw()

    name = 'hUnfolded_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), obs_label, int(min_pt_truth), int(max_pt_truth), self.file_format)
    if plot_pythia:
      name = 'hUnfolded_R{}_{}_{}-{}_Pythia{}'.format(self.utils.remove_periods(jetR), obs_label, int(min_pt_truth), int(max_pt_truth), self.file_format)
    output_dir = getattr(self, 'output_dir_final_results')
    outputFilename = os.path.join(output_dir, name)
    c.SaveAs(outputFilename)
    c.Close()

    # Write result to ROOT file
    if not plot_pythia:
      final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
      fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
      h.Write()
      fFinalResults.Close()

  #----------------------------------------------------------------------
  def get_pythia_from_response(self, observable, jetR, obs_label, min_pt_truth, max_pt_truth):
  
    output_dir = getattr(self, 'output_dir_main')
    file = os.path.join(output_dir, 'response.root')
    f = ROOT.TFile(file, 'READ')

    if observable == 'theta_g':
      thn_name = 'hResponse_JetPt_ThetaG_R{}_{}_rebinned'.format(jetR, obs_label)
    if observable == 'zg':
      thn_name = 'hResponse_JetPt_zg_R{}_{}_rebinned'.format(jetR, obs_label)
    
    thn = f.Get(thn_name)
    thn.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)

    h = thn.Projection(3)
    h.SetName('hPythia_{}_R{}_{}_{}-{}'.format(observable, jetR, obs_label, min_pt_truth, max_pt_truth))
    h.SetDirectory(0)

    return h

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
        if value1>value2:
          avg = value1
        else:
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
