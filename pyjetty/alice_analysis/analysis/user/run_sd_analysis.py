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
from pyjetty.alice_analysis.analysis.base import analysis_utils
import roounfold_sd

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class run_sd_analysis(common_base.common_base):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(run_sd_analysis, self).__init__(**kwargs)
    self.config_file = config_file
    
    # Initialize utils class
    self.utils = analysis_utils.analysis_utils()
    
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
    self.observables = config['observables']
    
    # Set which analysis steps to perform
    self.do_unfolding = config['do_unfolding']
    self.do_systematics = config['do_systematics']
    self.do_plot_final_result = config['do_plot_final_result']
    
    self.force_rebin_response=config['force_rebin']
    
    # Retrieve list of SD grooming settings
    self.jetR_list = config['jetR']
    self.sd_config_dict = config['SoftDrop']
    self.sd_config_list = list(self.sd_config_dict.keys())
    self.sd_settings = [[self.sd_config_dict[name]['zcut'], self.sd_config_dict[name]['beta']] for name in self.sd_config_list]
    
    for observable in self.observables:
    
      if observable == 'theta_g':
        xtitle = '#theta_{g}'
        ytitle = '#frac{1}{#it{N}_{jets, inc}} #frac{d#it{N}}{d#theta_{g}}'
      if observable == 'zg':
        xtitle = '#it{z}_{g}'
        ytitle = '#frac{1}{#it{N}_{jets, inc}} #frac{d#it{N}}{d#it{z}_{g}}'
      setattr(self, 'xtitle_{}'.format(observable), xtitle)
      setattr(self, 'ytitle_{}'.format(observable), ytitle)

    # Retrieve histogram binnings for each SD setting
    for i, sd_setting in enumerate(self.sd_settings):
        
      zcut = sd_setting[0]
      beta = sd_setting[1]
      sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)
      config_name = self.sd_config_list[i]
      
      pt_bins_truth = (self.sd_config_dict[config_name]['pt_bins_truth'])
      
      n_pt_bins_truth = len(pt_bins_truth) - 1
      setattr(self, 'n_pt_bins_truth_{}'.format(sd_label), n_pt_bins_truth)

      truth_pt_bin_array = array('d',pt_bins_truth)
      setattr(self, 'truth_pt_bin_array_{}'.format(sd_label), truth_pt_bin_array)

      if 'theta_g' in self.observables:
      
        rg_bins_truth = (self.sd_config_dict[config_name]['rg_bins_truth'])
      
        n_rg_bins_truth = len(rg_bins_truth) - 1
        setattr(self, 'n_rg_bins_truth_{}'.format(sd_label), n_rg_bins_truth)
        
        truth_rg_bin_array = array('d',rg_bins_truth)
        setattr(self, 'truth_rg_bin_array_{}'.format(sd_label), truth_rg_bin_array)
      
      if 'zg' in self.observables:
        
        zg_bins_truth = (self.sd_config_dict[config_name]['zg_bins_truth'])

        n_zg_bins_truth = len(zg_bins_truth) - 1
        setattr(self, 'n_zg_bins_truth_{}'.format(sd_label), n_zg_bins_truth)
        
        truth_zg_bin_array = array('d',zg_bins_truth)
        setattr(self, 'truth_zg_bin_array_{}'.format(sd_label), truth_zg_bin_array)
    
    self.min_pt_reported = 20
    self.max_pt_reported = 80
    self.regularizationParamName = 'n_iter'

    # List of possible systematic variations
    # Regularization parameter is automatically included
    self.kMain, self.kUnfoldingRange1, self.kUnfoldingRange2, self.kPrior1, self.kPrior2, self.kTrackEff = range(0, 6)

    # Set which systematics should be performed
    self.systematics_list = [self.kMain, self.kPrior1, self.kTrackEff]
    
    # Load paths to processing output, to be unfolded
    self.main_data = config['main_data']
    self.main_response = config['main_response']
    
    if self.kTrackEff in self.systematics_list:
      self.trkeff_data = config['trkeff_data']
      self.trkeff_response = config['trkeff_response']
          
    # Create output dirs
    self.file_format = config['file_format']
    self.output_dir = config['output_dir']
    for observable in self.observables:
      self.create_output_dirs(observable)

    # Path to pythia results
    self.pythia_20 = config['pythia_20']
    self.pythia_40 = config['pythia_40']
    self.pythia_60 = config['pythia_60']

  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def create_output_dirs(self, observable):
    
    output_dir = os.path.join(self.output_dir, observable)
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

    output_dir_main = os.path.join(output_dir, 'main')
    setattr(self, 'output_dir_main_{}'.format(observable), output_dir_main)
    if not os.path.isdir(output_dir_main):
      os.makedirs(output_dir_main)

    if self.kTrackEff in self.systematics_list:
      output_dir_trkeff = os.path.join(output_dir, 'trkeff')
      setattr(self, 'output_dir_trkeff_{}'.format(observable), output_dir_trkeff)
      if not os.path.isdir(output_dir_trkeff):
        os.makedirs(output_dir_trkeff)
      
    if self.kPrior1 in self.systematics_list:
      output_dir_prior = os.path.join(output_dir, 'prior')
      setattr(self, 'output_dir_prior_{}'.format(observable), output_dir_prior)
      if not os.path.isdir(output_dir_prior):
        os.makedirs(output_dir_prior)

    if self.do_systematics:
      output_dir_systematics = os.path.join(output_dir, 'systematics')
      setattr(self, 'output_dir_systematics_{}'.format(observable), output_dir_systematics)
      if not os.path.isdir(output_dir_systematics):
        os.makedirs(output_dir_systematics)
      
    if self.do_plot_final_result:
      output_dir_final = os.path.join(output_dir, 'final_results')
      setattr(self, 'output_dir_final_{}'.format(observable), output_dir_final)
      if not os.path.isdir(output_dir_final):
        os.makedirs(output_dir_final)
      
  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def run_sd_analysis(self):

    for observable in self.observables:
    
      if self.do_unfolding:
        self.perform_unfolding(observable)
      
      for jetR in self.jetR_list:
        for sd_setting in self.sd_settings:
          
          zcut = sd_setting[0]
          beta = sd_setting[1]
          sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)

          if self.do_systematics:
            self.compute_systematics(observable, jetR, sd_label)

          if self.do_plot_final_result:
            self.plot_final_result(observable, jetR, sd_label, zcut, beta)
      
  #----------------------------------------------------------------------
  def perform_unfolding(self, observable):
    print('Perform unfolding for all systematic variations: {} ...'.format(observable))
    
    # Main result
    output_dir = getattr(self, 'output_dir_main_{}'.format(observable))
    analysis_main = roounfold_sd.roounfold_sd(observable, self.main_data, self.main_response, self.config_file, output_dir, self.file_format, rebin_response=self.check_rebin_response(output_dir))
    analysis_main.roounfold_sd()

    # Tracking efficiency variation
    if self.kTrackEff in self.systematics_list:
      output_dir = getattr(self, 'output_dir_trkeff_{}'.format(observable))
      analysis_trkeff = roounfold_sd.roounfold_sd(observable, self.trkeff_data, self.trkeff_response, self.config_file, output_dir, self.file_format, rebin_response=self.check_rebin_response(output_dir))
      analysis_trkeff.roounfold_sd()

    # Prior variation
    if self.kPrior1 in self.systematics_list:
      output_dir = getattr(self, 'output_dir_prior_{}'.format(observable))
      analysis_prior = roounfold_sd.roounfold_sd(observable, self.main_data, self.main_response, self.config_file, output_dir, self.file_format, rebin_response=self.check_rebin_response(output_dir), power_law_offset=0.5)
      analysis_prior.roounfold_sd()

  #----------------------------------------------------------------------
  def check_rebin_response(self, output_dir):
    
    rebin_response = True
    response_path = os.path.join(output_dir, 'response.root')
    if os.path.exists(response_path):
      rebin_response = False
      if self.force_rebin_response:
        print('Response {} exists -- force re-create...'.format(response_path))
      else:
        print('Response {} exists -- don\'t re-create.'.format(response_path))
    else:
      print('Response {} doesn\'t exist -- create it...'.format(response_path))

    return rebin_response

  #----------------------------------------------------------------------
  def compute_systematics(self, observable, jetR, sd_label):
    print('Compute systematics for {}: R = {}, {} ...'.format(observable, jetR, sd_label))

    # Get main result
    output_dir = getattr(self, 'output_dir_main_{}'.format(observable))
    path_main = os.path.join(output_dir, 'fResult_R{}_{}.root'.format(jetR, sd_label))
    fMain = ROOT.TFile(path_main, 'READ')
    
    reg_param_final = self.utils.get_reg_param(self.sd_settings, self.sd_config_list, self.sd_config_dict, sd_label, observable, jetR)

    name = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final)
    hMain = fMain.Get(name)
    hMain.SetDirectory(0)
    setattr(self, name, hMain)
    
    # Tagging rate histogram
    name = 'hTaggingFractions_R{}_{}'.format(jetR, sd_label)
    hTaggingFractions = fMain.Get(name)
    hTaggingFractions.SetDirectory(0)
    setattr(self, name, hTaggingFractions)

    # Regularization parameter +2
    name = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final+2)
    hRegParam1 = fMain.Get(name)
    hRegParam1.SetDirectory(0)
    setattr(self, name, hRegParam1)
    
    # Regularization parameter -2
    name = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final-2)
    hRegParam2 = fMain.Get(name)
    hRegParam2.SetDirectory(0)
    setattr(self, name, hRegParam2)
    
    # Get trkeff result
    if self.kTrackEff in self.systematics_list:
      output_dir = getattr(self, 'output_dir_trkeff_{}'.format(observable))
      path_trkeff = os.path.join(output_dir, 'fResult_R{}_{}.root'.format(jetR, sd_label))
      fTrkEff = ROOT.TFile(path_trkeff, 'READ')
      name = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final)
      hTrkEff = fTrkEff.Get(name)
      hTrkEff.SetDirectory(0)
      setattr(self, '{}_trkeff'.format(name), hTrkEff)
    
    # Get prior result
    if self.kPrior1 in self.systematics_list:
      output_dir = getattr(self, 'output_dir_prior_{}'.format(observable))
      path_prior1 = os.path.join(output_dir, 'fResult_R{}_{}.root'.format(jetR, sd_label))
      fPrior1 = ROOT.TFile(path_prior1, 'READ')
      name = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final)
      hPrior1 = fPrior1.Get(name)
      hPrior1.SetDirectory(0)
      setattr(self, '{}_prior1'.format(name), hPrior1)
    
    # Loop through pt slices, and compute systematics for each 1D theta_g distribution
    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(sd_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(sd_label))
    
    for bin in range(1, n_pt_bins_truth-3):
      min_pt_truth = truth_pt_bin_array[bin]
      max_pt_truth = truth_pt_bin_array[bin+1]
      
      self.compute_sd_observable_systematic(observable, jetR, sd_label, min_pt_truth, max_pt_truth)

  #----------------------------------------------------------------------
  def get_sd_observable_distribution(self, jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth):
  
    h2D = getattr(self, name2D)
    h2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    h = h2D.ProjectionY() # Better to use ProjectionY('{}_py'.format(h2D.GetName()), 1, h2D.GetNbinsX()) ?
    h.SetDirectory(0)
    
    name = 'hTaggingFractions_R{}_{}'.format(jetR, sd_label)
    hTaggingFrac = getattr(self, name)
    x = (min_pt_truth + max_pt_truth)/2.
    fraction_tagged =  hTaggingFrac.GetBinContent(hTaggingFrac.FindBin(x))
    setattr(self, '{}_fraction_tagged'.format(name1D), fraction_tagged)

    n_jets_tagged = h.Integral(1, h.GetNbinsX())
    n_jets_inclusive = n_jets_tagged/fraction_tagged
    h.Scale(1./n_jets_inclusive, 'width')
  
    setattr(self, name1D, h)

    return h

  #----------------------------------------------------------------------
  def compute_sd_observable_systematic(self, observable, jetR, sd_label, min_pt_truth, max_pt_truth):
    
    #------------------------------------
    # Get 1D histograms
    # Normalize by integral, i.e. N_jets,inclusive in this pt-bin (cross-check this)
    
    reg_param_final = self.utils.get_reg_param(self.sd_settings, self.sd_config_list, self.sd_config_dict, sd_label, observable, jetR)
    
    # Get main histogram
    name2D = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final)
    name1D = 'hMain_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    hMain = self.get_sd_observable_distribution(jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth)

    # Get reg param +2
    name2D = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final+2)
    name1D = 'hRegParam1_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    hRegParam1 = self.get_sd_observable_distribution(jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth)
    
    # Get reg param -2
    name2D = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final-2)
    name1D = 'hRegParam2_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    hRegParam2 = self.get_sd_observable_distribution(jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth)

    # Get trk eff
    if self.kTrackEff in self.systematics_list:
      name2D = 'hUnfolded_{}_R{}_{}_{}_trkeff'.format(observable, jetR, sd_label, reg_param_final)
      name1D = 'hTrkEff_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hTrkEff = self.get_sd_observable_distribution(jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth)
    
    # Get prior1
    if self.kPrior1 in self.systematics_list:
      name2D = 'hUnfolded_{}_R{}_{}_{}_prior1'.format(observable, jetR, sd_label, reg_param_final)
      name1D = 'hPrior1_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hPrior1 = self.get_sd_observable_distribution(jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth)
    
    #------------------------------------
    # Compute systematics

    # Reg param +2
    name = 'hSystematic_{}_RegParam1_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    hSystematic_RegParam1 = hMain.Clone()
    hSystematic_RegParam1.SetName(name)
    hSystematic_RegParam1.Divide(hRegParam1)
    self.change_to_per(hSystematic_RegParam1)
    setattr(self, name, hSystematic_RegParam1)

    # Reg param -2
    name = 'hSystematic_{}_RegParam2_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    hSystematic_RegParam2 = hMain.Clone()
    hSystematic_RegParam2.SetName(name)
    hSystematic_RegParam2.Divide(hRegParam2)
    self.change_to_per(hSystematic_RegParam2)
    setattr(self, name, hSystematic_RegParam2)
    
    name = 'hSystematic_{}_RegParam_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    hSystematic_RegParam = self.build_average(hSystematic_RegParam1, hSystematic_RegParam2)
    setattr(self, name, hSystematic_RegParam)
    
    output_dir = getattr(self, 'output_dir_systematics_{}'.format(observable))
    outputFilename = os.path.join(output_dir, 'hSystematic_RegParam_R{}_{}_{}-{}{}'.format(jetR, sd_label, min_pt_truth, max_pt_truth, self.file_format))
    self.utils.plot_hist(hSystematic_RegParam, outputFilename, 'P E')
    
    # Prior
    hSystematic_Prior1 = None
    if self.kPrior1 in self.systematics_list:
      name = 'hSystematic_{}_Prior1_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hSystematic_Prior1 = hMain.Clone()
      hSystematic_Prior1.SetName(name)
      hSystematic_Prior1.Divide(hPrior1)
      self.change_to_per(hSystematic_Prior1)
      setattr(self, name, hSystematic_Prior1)
      
      output_dir = getattr(self, 'output_dir_systematics_{}'.format(observable))
      outputFilename = os.path.join(output_dir, 'hSystematic_Prior1_R{}_{}_{}-{}{}'.format(jetR, sd_label, min_pt_truth, max_pt_truth, self.file_format))
      self.utils.plot_hist(hSystematic_Prior1, outputFilename, 'P E')
    
    # Add shape uncertainties in quadrature
    hSystematic_Shape = self.add_in_quadrature(hSystematic_RegParam, hSystematic_Prior1)
    
    name = 'hResult_{}_shape_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    hResult_shape = hMain.Clone()
    hResult_shape.SetName(name)
    hResult_shape.SetDirectory(0)
    self.AttachErrToHist(hResult_shape, hSystematic_Shape)
    setattr(self, name, hResult_shape)
    
    # Trk eff
    if self.kTrackEff in self.systematics_list:
      name = 'hSystematic_{}_TrkEff_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hSystematic_TrkEff = hMain.Clone()
      hSystematic_TrkEff.SetName(name)
      hSystematic_TrkEff.Divide(hTrkEff)
      self.change_to_per(hSystematic_TrkEff)
      setattr(self, name, hSystematic_TrkEff)
      
      output_dir = getattr(self, 'output_dir_systematics_{}'.format(observable))
      outputFilename = os.path.join(output_dir, 'hSystematic_TrkEff_R{}_{}_{}-{}{}'.format(jetR, sd_label, min_pt_truth, max_pt_truth, self.file_format))
      self.utils.plot_hist(hSystematic_TrkEff, outputFilename, 'P E')
    
      name = 'hResult_{}_trkeff_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hResult_trkeff = hMain.Clone()
      hResult_trkeff.SetName(name)
      hResult_trkeff.SetDirectory(0)
      self.AttachErrToHist(hResult_trkeff, hSystematic_TrkEff)
      setattr(self, name, hResult_trkeff)

  #----------------------------------------------------------------------
  def add_in_quadrature(self, h1, h2):
  
    if not h2:
      return h1
  
    h_new = h1.Clone()
    h_new.SetName('{}_new'.format(h1.GetName()))
    
    for i in range(1, h_new.GetNbinsX()+1):
      value1 = h1.GetBinContent(i)
      value2 = h2.GetBinContent(i)
      new_value = math.sqrt(value1*value1 + value2*value2)
      
    h_new.SetBinContent(i, new_value)
    
    return h_new
  
  #----------------------------------------------------------------------
  def plot_final_result(self, observable, jetR, sd_label, zcut, beta):
    print('Plot final results for {}: R = {}, {} ...'.format(observable, jetR, sd_label))

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    # Loop through pt slices, and compute systematics for each 1D theta_g distribution
    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(sd_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(sd_label))
  
    for bin in range(1, n_pt_bins_truth-3):
      min_pt_truth = truth_pt_bin_array[bin]
      max_pt_truth = truth_pt_bin_array[bin+1]
      
      self.plot_observable(observable, jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth, plot_pythia=False)
      self.plot_observable(observable, jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth, plot_pythia=True)

  #----------------------------------------------------------------------
  def plot_observable(self, observable, jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth, plot_pythia=False):
  
    name = 'cResult_R{}_{}_{}-{}'.format(jetR, sd_label, min_pt_truth, max_pt_truth)
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
    
    if observable == 'theta_g':
      n_rg_bins_truth = getattr(self, 'n_rg_bins_truth_{}'.format(sd_label))
      truth_rg_bin_array = getattr(self, 'truth_rg_bin_array_{}'.format(sd_label))
      ymax = 4.
    if observable == 'zg':
      n_rg_bins_truth = getattr(self, 'n_zg_bins_truth_{}'.format(sd_label))
      truth_rg_bin_array = getattr(self, 'truth_zg_bin_array_{}'.format(sd_label))
      ymax = 10.
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_rg_bins_truth, truth_rg_bin_array)
    myBlankHisto.SetNdivisions(505)
    xtitle = getattr(self, 'xtitle_{}'.format(observable))
    ytitle = getattr(self, 'ytitle_{}'.format(observable))
    myBlankHisto.SetXTitle(xtitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle(ytitle)
    myBlankHisto.SetMaximum(ymax)
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")

    if plot_pythia and observable == 'theta_g':
      if math.fabs(min_pt_truth - 20) < 1e-3:
        fPythia_name = self.pythia_20
        hname = 'theta_g_beta{}_pt20'.format(beta)
      if math.fabs(min_pt_truth - 40) < 1e-3:
        fPythia_name = self.pythia_40
        hname = 'theta_g_beta{}_pt40'.format(beta)
      if math.fabs(min_pt_truth - 60) < 1e-3:
        fPythia_name = self.pythia_60
        hname = 'theta_g_beta{}_pt60'.format(beta)
      
      fPythia = ROOT.TFile(fPythia_name, 'READ')
      hPythia = fPythia.Get(hname)
      n_jets_inclusive = hPythia.Integral(0, hPythia.GetNbinsX()+1)
      n_jets_tagged = hPythia.Integral(hPythia.FindBin(truth_rg_bin_array[0]), hPythia.GetNbinsX())
      fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
      hPythia.Scale(1./n_jets_inclusive, 'width')
      hPythia.SetFillStyle(0)
      hPythia.SetMarkerSize(1.5)
      hPythia.SetMarkerStyle(23)
      hPythia.SetMarkerColor(1)
      hPythia.SetLineColor(1)
      hPythia.SetLineWidth(1)
      hPythia.Draw('E2 same')
    
    name = 'hMain_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    h = getattr(self, name)
    fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
    color = 600-6
    h.SetMarkerSize(1.5)
    h.SetMarkerStyle(20)
    h.SetMarkerColor(color)
    h.SetLineStyle(1)
    h.SetLineWidth(2)
    h.SetLineColor(color)
    h.DrawCopy('PE X0 same')
    
    h_shape = getattr(self, 'hResult_{}_shape_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth))
    h_shape.SetLineColor(0)
    h_shape.SetFillColor(color)
    h_shape.SetFillColorAlpha(color, 0.3)
    h_shape.SetFillStyle(1001)
    h_shape.Draw('E2 same')
    
    h_corr = None
    if self.kTrackEff in self.systematics_list:
      h_corr = getattr(self, 'hResult_{}_trkeff_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth))
      h_corr.SetFillStyle(0)
      h_corr.SetLineColor(color)
      h_corr.SetMarkerColorAlpha(color, 0)
      h_corr.SetLineWidth(1)
      h_corr.Draw('E2 same')
  
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
    text = '#it{z}_{cut} = ' + str(zcut) + '  #beta = ' + str(beta)
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.59, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text_latex.SetTextSize(0.04)
    text = '#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged
    text_latex.DrawLatex(0.57, 0.52, text)
    
    if plot_pythia and observable == 'theta_g':
      text_latex = ROOT.TLatex()
      text_latex.SetNDC()
      text_latex.SetTextSize(0.04)
      text = ('#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged) + (', #it{f}_{tagged}^{pythia} = %3.3f' % fraction_tagged_pythia)
      text_latex.DrawLatex(0.57, 0.52, text)

    myLegend = ROOT.TLegend(0.25,0.7,0.5,0.85)
    self.utils.setup_legend(myLegend,0.035)
    myLegend.AddEntry(h, 'ALICE pp', 'pe')
    myLegend.AddEntry(h_corr, 'Correlated uncertainty', 'f')
    myLegend.AddEntry(h_shape, 'Shape uncertainty', 'f')
    if plot_pythia and observable == 'theta_g':
      myLegend.AddEntry(hPythia, 'PYTHIA', 'pe')
    myLegend.Draw()

    name = 'hUnfolded_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), sd_label, int(min_pt_truth), int(max_pt_truth), self.file_format)
    if plot_pythia and observable == 'theta_g':
      name = 'hUnfolded_R{}_{}_{}-{}_Pythia{}'.format(self.utils.remove_periods(jetR), sd_label, int(min_pt_truth), int(max_pt_truth), self.file_format)
    output_dir = getattr(self, 'output_dir_final_{}'.format(observable))
    outputFilename = os.path.join(output_dir, name)
    c.SaveAs(outputFilename)
    c.Close()

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

  analysis = run_sd_analysis(config_file = args.configFile)
  analysis.run_sd_analysis()
