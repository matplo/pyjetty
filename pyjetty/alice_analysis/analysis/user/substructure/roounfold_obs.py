#! /usr/bin/env python

import sys
import os
import argparse
import itertools
from array import *
import numpy as np
import ROOT
import yaml

# Analysis utilities
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs
from pyjetty.alice_analysis.analysis.base import analysis_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

# Suppress a lot of standard output
ROOT.gErrorIgnoreLevel = ROOT.kWarning

################################################################
class Roounfold_Obs(analysis_base.AnalysisBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, observable='', input_file_data='', input_file_response='', config_file='',
               output_dir='', file_format='', rebin_response=False, truncation=False,
               binning=False, prior_variation_parameter=0., **kwargs):

    super(Roounfold_Obs, self).__init__(input_file_data, input_file_response, config_file,
                                        output_dir, file_format, **kwargs)
    self.utils = analysis_utils_obs.AnalysisUtils_Obs(observable)

    self.fData = ROOT.TFile(self.input_file_data, 'READ')
    self.fResponse = ROOT.TFile(self.input_file_response, 'READ')
    self.observable = observable
    self.truncation = truncation
    self.binning = binning
    self.prior_variation_parameter = prior_variation_parameter

    self.initialize_config()

    self.get_responses(rebin_response)

    # Create output files to store results
    for jetR in self.jetR_list:
      for i, _ in enumerate(self.obs_subconfig_list):

        obs_setting = self.obs_settings[i]
        grooming_setting = self.grooming_settings[i]
        obs_label = self.utils.obs_label(obs_setting, grooming_setting)

        fResult_name = os.path.join(self.output_dir, 'fResult_R{}_{}.root'.format(jetR, obs_label))
        setattr(self, 'fResult_name_R{}_{}'.format(jetR, obs_label), fResult_name)
        fResult = ROOT.TFile(fResult_name, 'RECREATE')
        fResult.Close()

    # Create output directories for unfolding plots
    self.create_output_dirs()

    print(self)

  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def roounfold_obs(self):

    for jetR in self.jetR_list:

      for i, _ in enumerate(self.obs_subconfig_list):

        obs_setting = self.obs_settings[i]
        grooming_setting = self.grooming_settings[i]
        obs_label = self.utils.obs_label(obs_setting, grooming_setting)
        self.unfold_single_setting(jetR, obs_label, obs_setting, grooming_setting)

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):

    # Call base class initialization
    analysis_base.AnalysisBase.initialize_config(self)

    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)

      # Load the RooUnfold library
      ROOT.gSystem.Load(config['roounfold_path'])

      # Get the sub-configs to unfold
      self.jetR_list = config['jetR']
      self.obs_config_dict = config[self.observable]
      self.obs_subconfig_list = [name for name in list(self.obs_config_dict.keys())
                                 if 'config' in name ]
      self.grooming_settings = self.utils.grooming_settings(self.obs_config_dict)
      self.obs_settings = self.utils.obs_settings(self.observable, self.obs_config_dict,
                                                  self.obs_subconfig_list)
      self.xtitle = self.obs_config_dict['common_settings']['xtitle']
      self.ytitle = self.obs_config_dict['common_settings']['ytitle']
      self.pt_bins_reported = self.obs_config_dict['common_settings']['pt_bins_reported']

      self.use_max_reg_param = False
      if 'max_reg_param' in self.obs_config_dict['common_settings']:
        self.max_reg_param = self.obs_config_dict['common_settings']['max_reg_param']
        self.use_max_reg_param = True

      # Retrieve histogram binnings for each observable setting
      for i, _ in enumerate(self.obs_subconfig_list):

        config_name = self.obs_subconfig_list[i]
        obs_label = self.utils.obs_label(self.obs_settings[i], self.grooming_settings[i])

        pt_det_bins_name = 'pt_bins_det'
        if self.truncation:
          pt_det_bins_name += '_sys_truncation'
        pt_bins_det = (self.obs_config_dict[config_name][pt_det_bins_name])
        pt_bins_truth = (self.obs_config_dict[config_name]['pt_bins_truth'])
        n_pt_bins_det = len(pt_bins_det) - 1
        setattr(self, 'n_pt_bins_det_{}'.format(obs_label), n_pt_bins_det)
        n_pt_bins_truth = len(pt_bins_truth) - 1
        setattr(self, 'n_pt_bins_truth_{}'.format(obs_label), n_pt_bins_truth)
        det_pt_bin_array = array('d',pt_bins_det)
        setattr(self, 'det_pt_bin_array_{}'.format(obs_label), det_pt_bin_array)
        truth_pt_bin_array = array('d',pt_bins_truth)
        setattr(self, 'truth_pt_bin_array_{}'.format(obs_label), truth_pt_bin_array)

        det_bins_name = 'obs_bins_det'
        truth_bins_name = 'obs_bins_truth'
        if self.binning:
          det_bins_name += '_sys_binning'
        obs_bins_det = (self.obs_config_dict[config_name][det_bins_name])
        obs_bins_truth = (self.obs_config_dict[config_name][truth_bins_name])
        n_obs_bins_det = len(obs_bins_det) - 1
        setattr(self, 'n_bins_det_{}'.format(obs_label), n_obs_bins_det)
        n_obs_bins_truth = len(obs_bins_truth) - 1
        setattr(self, 'n_bins_truth_{}'.format(obs_label), n_obs_bins_truth)
        det_obs_bin_array = array('d',obs_bins_det)
        setattr(self, 'det_bin_array_{}'.format(obs_label), det_obs_bin_array)
        truth_obs_bin_array = array('d',obs_bins_truth)
        setattr(self, 'truth_bin_array_{}'.format(obs_label), truth_obs_bin_array)

        for jetR in self.jetR_list:

          name_thn = self.utils.name_thn(self.observable, jetR, obs_label)
          name_thn_rebinned = self.utils.name_thn_rebinned(self.observable, jetR, obs_label)
          name_data = self.utils.name_data(self.observable, jetR, obs_label)
          name_data_rebinned = self.utils.name_data_rebinned(self.observable, jetR, obs_label)

          name_roounfold = 'roounfold_response_R{}_{}'.format(jetR, obs_label)
          setattr(self, 'name_thn_R{}_{}'.format(jetR, obs_label), name_thn)
          setattr(self, 'name_thn_rebinned_R{}_{}'.format(jetR, obs_label), name_thn_rebinned)
          setattr(self, 'name_data_R{}_{}'.format(jetR, obs_label), name_data)
          setattr(self, 'name_data_rebinned_R{}_{}'.format(jetR, obs_label), name_data_rebinned)
          setattr(self, 'name_roounfold_R{}_{}'.format(jetR, obs_label), name_roounfold)

      self.reg_param_name = 'n_iter'
      self.errorType = ROOT.RooUnfold.kCovToy

  #---------------------------------------------------------------
  # Get responses, either from file or manually rebin
  #---------------------------------------------------------------
  def get_responses(self, rebin_response=False):

    response_file_name = os.path.join(self.output_dir, 'response.root')
    if rebin_response:
      f = ROOT.TFile(response_file_name, 'RECREATE')
      f.Close()

    # Rebin response matrix, and create RooUnfoldResponse object
    # THn response matrix is: (pt-det, pt-true, obs-det, obs-true)
    for jetR in self.jetR_list:

      for i, _ in enumerate(self.obs_subconfig_list):

        obs_label = self.utils.obs_label(self.obs_settings[i], self.grooming_settings[i])

        name_thn = getattr(self, 'name_thn_R{}_{}'.format(jetR, obs_label))
        name_thn_rebinned = getattr(self, 'name_thn_rebinned_R{}_{}'.format(jetR, obs_label))
        name_data = getattr(self, 'name_data_R{}_{}'.format(jetR, obs_label))
        name_roounfold = getattr(self, 'name_roounfold_R{}_{}'.format(jetR, obs_label))

        # Retrieve desired binnings
        n_pt_bins_det = getattr(self, 'n_pt_bins_det_{}'.format(obs_label))
        det_pt_bin_array = getattr(self, 'det_pt_bin_array_{}'.format(obs_label))
        n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(obs_label))
        truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(obs_label))

        n_bins_det = getattr(self, 'n_bins_det_{}'.format(obs_label))
        det_bin_array = getattr(self, 'det_bin_array_{}'.format(obs_label))
        n_bins_truth = getattr(self, 'n_bins_truth_{}'.format(obs_label))
        truth_bin_array = getattr(self, 'truth_bin_array_{}'.format(obs_label))

        # Rebin if requested, and write to file
        thn = self.fResponse.Get(name_thn)
        thn.SetName(name_thn)
        setattr(self, name_thn, thn)
        if rebin_response:
          # Create rebinned THn and RooUnfoldResponse with these binnings, and write to file
          label = 'R{}_{}'.format(jetR, obs_label)
          self.utils.rebin_response(response_file_name, thn, name_thn_rebinned, name_roounfold,
                                    label, n_pt_bins_det, det_pt_bin_array, n_bins_det,
                                    det_bin_array, n_pt_bins_truth, truth_pt_bin_array,
                                    n_bins_truth, truth_bin_array, self.observable,
                                    self.prior_variation_parameter)

        # Also re-bin the data histogram
        hData = self.fData.Get(name_data)
        h = self.utils.rebin_data(hData, name_data, n_pt_bins_det, det_pt_bin_array,
                                  n_bins_det, det_bin_array)
        h.SetDirectory(0)
        name = '{}_{}'.format(name_data, 'rebinned')
        setattr(self, name, h)

        # Retrieve responses from file
        f = ROOT.TFile(response_file_name, 'READ')
        thn_rebinned = f.Get(name_thn_rebinned)
        #thn_rebinned.SetDirectory(0)
        roounfold_response = f.Get(name_roounfold)
        roounfold_response.UseOverflow(False)

        setattr(self, name_thn_rebinned, thn_rebinned)
        setattr(self, name_roounfold, roounfold_response)
        f.Close()

  #---------------------------------------------------------------
  # Create a set of output directories for a given observable
  #---------------------------------------------------------------
  def create_output_dirs(self):

    dirs = ['RM', 'Data', 'KinematicEfficiency', 'Unfolded_obs', 'Unfolded_pt',
            'Unfolded_ratio', 'Unfolded_stat_uncert', 'Test_StatisticalClosure',
            'Test_Refolding', 'Correlation_Coefficients']
    for i in dirs:
      output_dir = os.path.join(self.output_dir, i)
      setattr(self, 'output_dir_{}'.format(i), output_dir)
      if not os.path.exists(output_dir):
        os.makedirs(output_dir)

  ###################################################################################################
  # Unfold jet substructure observable for a single setting
  ###################################################################################################
  def unfold_single_setting(self, jetR, obs_label, obs_setting, grooming_setting):

    print('jetR = {}, {}'.format(jetR, obs_label))

    # Plot data 2D histogram
    name = getattr(self, 'name_data_rebinned_R{}_{}'.format(jetR, obs_label))
    hData_PerBin = getattr(self, name)
    hData_PerBin.GetXaxis().SetTitle('#it{p}_{T,jet}')
    hData_PerBin.GetYaxis().SetTitle(self.xtitle)
    output_dir = getattr(self, 'output_dir_Data')
    outf_name = 'hData_R{}_{}{}'.format(self.utils.remove_periods(jetR),
                                        obs_label, self.file_format)
    outf_name = os.path.join(output_dir, outf_name)
    self.utils.plot_hist(hData_PerBin, outf_name, 'colz', False, True)

    # Plot various slices of the response matrix (from the THn)
    self.plot_RM_slices(jetR, obs_label)

    # Plot the kinematic efficiency from the response THn, and save it as an attribute
    self.plot_kinematic_efficiency(jetR, obs_label, obs_setting, grooming_setting)

    # Get MC-det and MC-truth 2D projections for (i) grooming tagging rate
    # bin-by-bin correction, and (ii) unfolding closure test
    name = 'hMC_Det_R{}_{}'.format(jetR, obs_label)
    hMC_Det = self.get_MCdet2D(jetR, obs_label)
    setattr(self, name, hMC_Det)

    name = 'hMC_Truth_R{}_{}'.format(jetR, obs_label)
    hMC_Truth = self.get_MCtruth2D(jetR, obs_label)
    setattr(self, name, hMC_Truth)

    # Compute grooming tagging rate
    if grooming_setting:
      self.compute_tagging_rate(jetR, obs_label)

    # Unfold spectrum
    if hData_PerBin and hMC_Det and hMC_Truth:

      self.unfold_observable(jetR, obs_label, obs_setting, grooming_setting)
      self.unfolding_checks(jetR, obs_label, obs_setting, grooming_setting)

  #################################################################################################
  # Compute grooming tagging rate, based on MC correction
  #################################################################################################
  def compute_tagging_rate(self, jetR, obs_label):

    # Get truth binnings, since we need to rebin the det-level histograms
    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(obs_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(obs_label))
    n_bins_truth = getattr(self, 'n_bins_truth_{}'.format(obs_label))
    truth_bin_array = getattr(self, 'truth_bin_array_{}'.format(obs_label))

    # Create a histogram to store the tagging fractions
    name = 'hTaggingFractions_R{}_{}'.format(jetR, obs_label)
    h = ROOT.TH1F(name, name, n_pt_bins_truth, truth_pt_bin_array)

    # Get relevant histograms
    name = getattr(self, 'name_data_rebinned_R{}_{}'.format(jetR, obs_label))
    hData2D = getattr(self, name).Clone()
    hData2D.SetName('{}_tagging'.format(hData2D.GetName()))
    hData2D_rebinned = self.utils.rebin_data(hData2D, name, n_pt_bins_truth, truth_pt_bin_array,
                                             n_bins_truth, truth_bin_array)

    name = 'hMC_Det_R{}_{}'.format(jetR, obs_label)
    hMC_Det = getattr(self, name).Clone()
    hMC_Det.SetName('{}_tagging'.format(hMC_Det.GetName()))
    hMC_Det_rebinned = self.utils.rebin_data(hMC_Det, name, n_pt_bins_truth, truth_pt_bin_array,
                                             n_bins_truth, truth_bin_array)

    name = 'hMC_Truth_R{}_{}'.format(jetR, obs_label)
    hMC_Truth = getattr(self, name).Clone()
    hMC_Truth.SetName('{}_tagging'.format(hMC_Truth.GetName()))

    # Compute the tagging fraction for each pt bin
    for bin in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[bin]
      max_pt_truth = self.pt_bins_reported[bin+1]

      fraction_tagged = self.utils.tagging_rate(jetR, min_pt_truth, max_pt_truth,
                                                hData2D_rebinned, hMC_Det_rebinned, hMC_Truth)

      x = (min_pt_truth + max_pt_truth)/2.
      bin = h.FindBin(x)
      h.SetBinContent(bin, fraction_tagged)

    fResult_name = getattr(self, 'fResult_name_R{}_{}'.format(jetR, obs_label))
    fResult = ROOT.TFile(fResult_name, 'UPDATE')
    h.Write()
    fResult.Close()

  #################################################################################################
  # Unfold jet spectrum
  #################################################################################################
  def unfold_observable(self, jetR, obs_label, obs_setting, grooming_setting):

    # Get response and data to unfold
    response = getattr(self, 'roounfold_response_R{}_{}'.format(jetR, obs_label))

    name_data = getattr(self, 'name_data_rebinned_R{}_{}'.format(jetR, obs_label))
    hData = getattr(self, name_data)

    # Write unfolded results to a ROOT file
    fResult_name = getattr(self, 'fResult_name_R{}_{}'.format(jetR, obs_label))
    fResult = ROOT.TFile(fResult_name, 'UPDATE')

    # Select final regularization parameter
    if self.use_max_reg_param:
      reg_param_final = self.max_reg_param
    else:
      reg_param_final = self.utils.get_reg_param(
        self.obs_settings, self.grooming_settings, self.obs_subconfig_list,
        self.obs_config_dict, obs_label, jetR)

    # Loop over values of regularization parameter
    for i in range(1, reg_param_final + 3):

      # Set up the Bayesian unfolding object
      unfold_bayes = ROOT.RooUnfoldBayes(response, hData, i)
      #unfoldBayes.SetNToys(1000)

      # Perform the unfolding
      print('Unfolding with {} = {}'.format(self.reg_param_name, i))
      hUnfolded = unfold_bayes.Hreco(self.errorType) # Produces the truth distribution, 
      # with errors, PerBin (will scale by bin width below, after refolding checks)

      # Save unfolded solution as class member
      name = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, i)
      hUnfolded.SetName(name)
      setattr(self, name, hUnfolded)
      hUnfolded.SetDirectory(0)

      # Correct by kinematic efficiency
      hKinematicEfficiency = getattr(self, 'hKinematicEfficiency_R{}_{}'.format(jetR, obs_label))
      hUnfolded.Divide(hKinematicEfficiency)

      # Write result to file
      hUnfolded.Write()

      # Plot Pearson correlation coeffs for each iteration, to get a measure of
      # the correlation between the bins
      covariance_matrix = unfold_bayes.Ereco(self.errorType) # Get the covariance matrix
      self.plot_correlation_coefficients(covariance_matrix, jetR, obs_label, i)

    fResult.Close()
    print('Done unfolding')

    # Plot unfolded results
    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    self.plot_unfolded_observable(jetR, obs_label, obs_setting, grooming_setting)
    self.plot_unfolded_pt(jetR, obs_label, obs_setting, grooming_setting)

  #################################################################################################
  # Plot unfolded observable for various pt slices
  #################################################################################################
  def plot_unfolded_observable(self, jetR, obs_label, obs_setting, grooming_setting):

    for bin in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[bin]
      max_pt_truth = self.pt_bins_reported[bin+1]

      self.plot_observable(jetR, obs_label, obs_setting, grooming_setting,
                           min_pt_truth, max_pt_truth)
      self.plot_observable(jetR, obs_label, obs_setting, grooming_setting,
                           min_pt_truth, max_pt_truth, option = 'ratio')
      self.plot_observable(jetR, obs_label, obs_setting, grooming_setting,
                           min_pt_truth, max_pt_truth, option = 'stat_uncert')

  #################################################################################################
  # Plot observable for a single pt slicee
  #################################################################################################
  def plot_observable(self, jetR, obs_label, obs_setting, grooming_setting,
                      min_pt_truth, max_pt_truth, option = ''):

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    name = 'cResult{}_R{}_{}_{}-{}'.format(option, jetR, obs_label, min_pt_truth, max_pt_truth)
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

    n_bins_truth = getattr(self, 'n_bins_truth_{}'.format(obs_label))
    truth_bin_array = getattr(self, 'truth_bin_array_{}'.format(obs_label))

    leg = ROOT.TLegend(0.75,0.65,0.88,0.92)
    self.utils.setup_legend(leg,0.04)

    # Select final regularization parameter
    if self.use_max_reg_param:
      reg_param_final = self.max_reg_param
    else:
      reg_param_final = self.utils.get_reg_param(
        self.obs_settings, self.grooming_settings, self.obs_subconfig_list,
        self.obs_config_dict, obs_label, jetR)

    for i in range(1, reg_param_final + 3):

      h = self.get_unfolded_result(jetR, obs_label, i, min_pt_truth, max_pt_truth, option)

      if option == 'ratio':
        if i > 1:
          h_previous = self.get_unfolded_result(jetR, obs_label, i-1, min_pt_truth,
                                                max_pt_truth, option)
          h.Divide(h_previous)
        else:
          continue

      elif option == 'stat_uncert':
        h = self.get_unfolded_result_uncertainties(jetR, obs_label, i, min_pt_truth,
                                                   max_pt_truth, option)

      # Save the histogram to the result ROOT file for future access
      if len(option):
        outh_name = 'hUnfolded_{}_{}_R{}_{}_n{}_Pt{}-{}'.format(
          self.observable, option, self.utils.remove_periods(jetR),
          obs_label, i, int(min_pt_truth), int(max_pt_truth) )
      else:
        outh_name = 'hUnfolded_{}_R{}_{}_n{}_Pt{}-{}'.format(
          self.observable, self.utils.remove_periods(jetR),
          obs_label, i, int(min_pt_truth), int(max_pt_truth) )
      h.SetNameTitle(outh_name, outh_name)

      if len(option):
        output_dir = getattr(self, 'output_dir_Unfolded_' + option)
      else:
        output_dir = getattr(self, 'output_dir_Unfolded_obs')
      if output_dir[-1] != '/':
        output_dir += '/'
      fResult_name = output_dir + 'fResult_name_R{}_{}.root'.format(jetR, obs_label)
      fResult = ROOT.TFile(fResult_name, 'UPDATE')
      h.Write()
      fResult.Close()

      # Set different markers for each element on superimposed plot
      if i == 1:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(20)
      elif i == 2:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(21)
      elif i == 3:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(22)
      elif i == 4:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(23)
      elif i == 5:
        h.SetMarkerSize(2)
        h.SetMarkerStyle(33)
      elif i == 6:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(34)
      elif i == 7:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(35)
      elif i == 8:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(36)
      else:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(19)
      h.SetMarkerColor(600-6+i)
      h.SetLineStyle(1)
      h.SetLineWidth(2)
      h.SetLineColor(600-6+i)

      # Doesn't work for some reason...
      #shift = 0.5*h.GetBinWidth(1)
      #h.GetXaxis().SetLimits(truth_bin_array[0]+shift,truth_bin_array[-1]+shift)

      # Create empty histogram for plotting superimposed results of each n_iter
      if i == 1 or (i == 2 and option == 'ratio'):
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram',
                                 n_bins_truth, truth_bin_array)
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(self.xtitle)
        myBlankHisto.GetYaxis().SetTitleOffset(1.5)
        myBlankHisto.SetYTitle(self.ytitle)
        myBlankHisto.SetMaximum(3*h.GetMaximum())
        myBlankHisto.SetMinimum(0.)
        if option == 'ratio':
          myBlankHisto.SetMaximum(1.2)
          myBlankHisto.SetMinimum(0.9)
          myBlankHisto.SetYTitle('#frac{n_{iter}}{(n_{iter}-1)}')
        if option == 'stat_uncert':
          myBlankHisto.SetMaximum(40)
          myBlankHisto.SetMinimum(0.)
          myBlankHisto.SetYTitle('statistical uncertainty (%)')
        myBlankHisto.Draw("E")

      # Plot the current calculation on the superimposed plot
      if option == 'ratio':
        h.DrawCopy('P hist X0 same')
      else:
        h.DrawCopy('PE X0 same')

      label = '{} = {}'.format(self.reg_param_name, i)
      leg.AddEntry(h, label, 'Pe')

    leg.Draw()

    # Draw horizontal line at y = 1
    if option == 'ratio':
      line = ROOT.TLine(truth_bin_array[0], 1, truth_bin_array[-1], 1)
      line.SetLineColor(920+2)
      line.SetLineStyle(2)
      line.SetLineWidth(4)
      line.Draw()

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth)
    text_latex.DrawLatex(0.35, 0.85, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = '#it{R} = ' + str(jetR)
    text_latex.DrawLatex(0.35, 0.78, text)

    subobs_label = self.utils.formatted_subobs_label(self.observable)
    if subobs_label:
      text = '{} = {}'.format(subobs_label, obs_setting)
      text_latex.DrawLatex(0.35, 0.71, text)

    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting)
      text_latex.DrawLatex(0.35, 0.64, text)

    if len(option):
      output_dir = getattr(self, 'output_dir_Unfolded_' + option)
    else:
      output_dir = getattr(self, 'output_dir_Unfolded_obs')

    outf_name = 'hUnfolded{}_{}_R{}_{}_{}-{}{}'.format(option, self.observable,
                                                       self.utils.remove_periods(jetR),
                                                       obs_label, int(min_pt_truth),
                                                       int(max_pt_truth), self.file_format)
    outf_name = os.path.join(output_dir, outf_name)
    c.SaveAs(outf_name)
    c.Close()

  #################################################################################################
  # Get unfolded result in 1D, for fixed slice of pt
  #################################################################################################
  def get_unfolded_result(self, jetR, obs_label, i, min_pt_truth, max_pt_truth, option):

    name = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, i)
    h2D = getattr(self, name)
    h2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    h = h2D.ProjectionY()
    h.SetName('{}_{}'.format(h.GetName(), option))

    # Normalize by integral, i.e. N_jets,inclusive in this pt-bin
    n_jets_inclusive = h.Integral(0, h.GetNbinsX()+1)
    h.Scale(1./n_jets_inclusive, 'width')

    return h

  #################################################################################################
  # Get unfolded result uncertainties in 1D, for fixed slice of pt
  #################################################################################################
  def get_unfolded_result_uncertainties(self, jetR, obs_label, i,
                                        min_pt_truth, max_pt_truth, option):

    h = self.get_unfolded_result(jetR, obs_label, i, min_pt_truth, max_pt_truth, option)

    for bin in range(1, h.GetNbinsX()+1):
      content = h.GetBinContent(bin)
      uncertainty = h.GetBinError(bin)
      if content < 1e-10:
        print('content of {} in bin {} is {}'.format(h.GetName(), bin, content))
      h.SetBinContent(bin, uncertainty/content * 100)
      h.SetBinError(bin, 0)

    return h

  #################################################################################################
  # Get unfolded pt in 1D
  #################################################################################################
  def plot_unfolded_pt(self, jetR, obs_label, obs_setting, grooming_setting):

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    name = 'cResultPt_R{}_{}'.format(jetR, obs_label)
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

    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(obs_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(obs_label))
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram',
                             n_pt_bins_truth, truth_pt_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle('#it{p}_{T, ch jet}')
    myBlankHisto.GetYaxis().SetTitleOffset(2.2)
    myBlankHisto.SetYTitle('#frac{dN}{d#it{p}_{T, ch jet}}')
    myBlankHisto.SetMaximum(5000)
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")

    leg = ROOT.TLegend(0.75,0.65,0.88,0.92)
    self.utils.setup_legend(leg,0.04)

    # Select final regularization parameter
    if self.use_max_reg_param:
      reg_param_final = self.max_reg_param
    else:
      reg_param_final = self.utils.get_reg_param(
        self.obs_settings, self.grooming_settings, self.obs_subconfig_list,
        self.obs_config_dict, obs_label, jetR)

    for i in range(1, reg_param_final + 3):

      name = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, i)
      h2D = getattr(self, name)
      h2D.GetXaxis().SetRangeUser(5., 120.)
      h = h2D.ProjectionX()

      h.Scale(1., 'width')

      if i == 1:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(20)
      elif i == 2:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(21)
      elif i == 3:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(22)
      elif i == 4:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(23)
      elif i == 5:
        h.SetMarkerSize(2)
        h.SetMarkerStyle(33)
      elif i == 6:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(34)
      elif i == 7:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(35)
      elif i == 8:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(36)
      else:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(19)
      h.SetMarkerColor(600-6+i)
      h.SetLineStyle(1)
      h.SetLineWidth(2)
      h.SetLineColor(600-6+i)

      h.DrawCopy('PE X0 same')

      label = '{} = {}'.format(self.reg_param_name, i)
      leg.AddEntry(h, label, 'Pe')

    leg.Draw()

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = '#it{R} = ' + str(jetR)
    text_latex.DrawLatex(0.35, 0.85, text)

    subobs_label = self.utils.formatted_subobs_label(self.observable)
    if subobs_label:
      text = '{} = {}'.format(subobs_label, obs_setting)
      text_latex.DrawLatex(0.35, 0.78, text)

    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting)
      text_latex.DrawLatex(0.35, 0.71, text)

    output_dir = getattr(self, 'output_dir_Unfolded_pt')
    outf_name = 'hUnfoldedPt_{}_R{}_{}{}'.format(self.observable,
                                                 self.utils.remove_periods(jetR),
                                                 obs_label, self.file_format)
    outf_name = os.path.join(output_dir, outf_name)
    c.SaveAs(outf_name)
    c.Close()

  ###################################################################################################
  # Plot kinematic efficiency
  # The kinematic efficiency is the ratio:
  #   Numerator: 2D truth-level projection [pt-true, obs-true] using 
  #              [pt-det, obs-det] cut on det-level
  #   Denominator: 2D truth-level projection [pt-true, obs-true] using no cut on det-level
  ###################################################################################################
  def plot_kinematic_efficiency(self, jetR, obs_label, obs_setting, grooming_setting):

    # (pt-det, pt-true, obs-det, obs-true)
    name_response = getattr(self, 'name_thn_rebinned_R{}_{}'.format(jetR, obs_label))
    response = getattr(self,  name_response)
    hResponse = response.Clone()
    hResponse.SetName('hResponse_KinEff_JetPt_Obs_R{}_{}'.format(jetR, obs_label))

    # Denominator -- by default, under/over-flow bins are included in projection
    hDenominator = hResponse.Projection(3, 1)
    hDenominator.SetName('{}_Denominator'.format(hDenominator.GetName()))

    # Numerator -- cut on det-level input binning
    det_pt_bin_array = getattr(self, 'det_pt_bin_array_{}'.format(obs_label))
    min_pt_det = det_pt_bin_array[0]
    max_pt_det = det_pt_bin_array[-1]
    det_bin_array = getattr(self, 'det_bin_array_{}'.format(obs_label))
    min_obs_det = det_bin_array[0]
    max_obs_det = det_bin_array[-1]
    hResponse.GetAxis(0).SetRangeUser(min_pt_det, max_pt_det)
    hResponse.GetAxis(2).SetRangeUser(min_obs_det, max_obs_det)
    hNumerator = hResponse.Projection(3, 1)
    hNumerator.SetName('{}_Numerator'.format(hNumerator.GetName()))

    hKinematicEfficiency = hNumerator.Clone()
    hKinematicEfficiency.GetYaxis().SetTitle(self.xtitle)
    hKinematicEfficiency.SetName('hKinematicEfficiency_R{}_{}'.format(jetR, obs_label))
    hKinematicEfficiency.Divide(hDenominator)
    for bin in range(0, hKinematicEfficiency.GetNcells()+1):
      hKinematicEfficiency.SetBinError(bin, 0)

    output_dir = getattr(self, 'output_dir_KinematicEfficiency')
    outf_name = 'hKinematicEfficiency2D_R{}_{}{}'.format(self.utils.remove_periods(jetR),
                                                         obs_label, self.file_format)
    outf_name = os.path.join(output_dir, outf_name)
    self.utils.plot_hist(hKinematicEfficiency, outf_name, 'colz')

    # Save kinematic efficiency as class member
    setattr(self, 'hKinematicEfficiency_R{}_{}'.format(jetR, obs_label), hKinematicEfficiency)

    # Plot 1D kinematic efficiency
    self.plot_kinematic_efficiency_projections(hKinematicEfficiency, jetR, obs_label,
                                               obs_setting, grooming_setting)

  #################################################################################################
  # Plot kinematic efficiency projections for various pt slices
  #################################################################################################
  def plot_kinematic_efficiency_projections(self, hKinematicEfficiency2D, jetR,
                                            obs_label, obs_setting, grooming_setting):

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    name = 'cKinEff_R{}_{}'.format(jetR, obs_label)
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

    n_bins_truth = getattr(self, 'n_bins_truth_{}'.format(obs_label))
    truth_bin_array = getattr(self, 'truth_bin_array_{}'.format(obs_label))
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_bins_truth, truth_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle(self.xtitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle('#varepsilon_{kin}')
    myBlankHisto.SetMaximum(1.2)
    myBlankHisto.SetMinimum(0.8)
    myBlankHisto.Draw("E")

    leg = ROOT.TLegend(0.6,0.65,0.72,0.92)
    self.utils.setup_legend(leg,0.04)

    for bin in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[bin]
      max_pt_truth = self.pt_bins_reported[bin+1]

      hKinematicEfficiency2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
      h = hKinematicEfficiency2D.ProjectionY()
      name = 'hKinematicEfficiency_R{}_{}_{}-{}'.format(jetR, obs_label,
                                                        min_pt_truth, max_pt_truth)
      h.SetName(name)

      if bin == 1:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(20)
      elif bin == 2:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(21)
      elif bin == 3:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(22)
      elif bin == 4:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(23)
      elif bin == 5:
        h.SetMarkerSize(2)
        h.SetMarkerStyle(33)
      else:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(19)
      h.SetMarkerColor(600-6+bin)
      h.SetLineStyle(1)
      h.SetLineWidth(2)
      h.SetLineColor(600-6+bin)

      h.DrawCopy('P X0 same')

      label = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
      leg.AddEntry(h, label, 'Pe')

    leg.Draw()

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = '#it{R} = ' + str(jetR)
    text_latex.DrawLatex(0.3, 0.85, text)

    subobs_label = self.utils.formatted_subobs_label(self.observable)
    if subobs_label:
      text = '{} = {}'.format(subobs_label, obs_setting)
      text_latex.DrawLatex(0.3, 0.78, text)

    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting)
      text_latex.DrawLatex(0.3, 0.71, text)

    line = ROOT.TLine(truth_bin_array[0], 1, truth_bin_array[-1], 1)
    line.SetLineColor(920+2)
    line.SetLineStyle(2)
    line.SetLineWidth(4)
    line.Draw()

    output_dir = getattr(self, 'output_dir_KinematicEfficiency')
    outf_name = 'hKinematicEfficiency_R{}_{}{}'.format(self.utils.remove_periods(jetR),
                                                       obs_label, self.file_format)
    outf_name = os.path.join(output_dir, outf_name)
    c.SaveAs(outf_name)
    c.Close()

  #################################################################################################
  # Plot various slices of the response matrix (from the THn)
  #################################################################################################
  def plot_RM_slices(self, jetR, obs_label):

    # (pt-det, pt-true, obs-det, obs-true)
    name_response = getattr(self, 'name_thn_R{}_{}'.format(jetR, obs_label))
    hResponse = getattr(self, name_response)

    # Fix pt-true, and plot the 2D observable response
    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(obs_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(obs_label))
    for bin in range(1, n_pt_bins_truth-1):
      min_pt_truth = truth_pt_bin_array[bin]
      max_pt_truth = truth_pt_bin_array[bin+1]

      self.plot_obs_response(jetR, obs_label, min_pt_truth, max_pt_truth, hResponse)

  #################################################################################################
  # Plot 2D observable response for a fixed range of pt-truth
  #################################################################################################
  def plot_obs_response(self, jetR, obs_label, min_pt_truth, max_pt_truth, hResponse):

    hResponse4D = hResponse.Clone()
    hResponse4D.SetName('{}_{}_{}'.format(hResponse4D.GetName(), min_pt_truth, max_pt_truth))

    hResponse4D.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)

    truth_bin_array = getattr(self, 'truth_bin_array_{}'.format(obs_label))
    hResponse4D.GetAxis(2).SetRangeUser(truth_bin_array[0], truth_bin_array[-1])
    hResponse4D.GetAxis(3).SetRangeUser(truth_bin_array[0], truth_bin_array[-1])

    hResponse_Obs = hResponse4D.Projection(3,2)
    hResponse_Obs.SetName('hResponse_Obs_R{}_{}_{}_{}'.format(self.utils.remove_periods(jetR),
                                                              obs_label, int(min_pt_truth),
                                                              int(max_pt_truth)))

    hResponse_Obs_Normalized = self.utils.normalize_response_matrix(hResponse_Obs)

    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet}^{true} < ' + str(max_pt_truth)

    output_dir = getattr(self, 'output_dir_RM')
    outf_name = '{}{}'.format(hResponse_Obs.GetName(), self.file_format)
    outf_name = os.path.join(output_dir, outf_name)
    self.utils.plot_hist(hResponse_Obs_Normalized, outf_name, 'colz', False, True, text)

  #################################################################################################
  # Plot correlation coefficients
  # Note: At the moment this directly plots the global bin index correlation...which
  #       as far as I can tell is what RooUnfold gives us (it maps the 4D RM into 2D)
  #################################################################################################
  def plot_correlation_coefficients(self, covariance_matrix, jetR, obs_label, i):

    nBinsX = covariance_matrix.GetNrows()
    nBinsY = covariance_matrix.GetNcols()
    
    correlation_coefficient_matrix = ROOT.TH2D('correlation_coefficient_matrix', 'correlation_coefficient_matrix', nBinsX, 0, nBinsX, nBinsY, 0, nBinsY)
    correlation_coefficient_matrix.GetXaxis().SetTitle('bin #')
    correlation_coefficient_matrix.GetYaxis().SetTitle('bin #')

    for xbin in range(0, nBinsX):
      varianceX = covariance_matrix(xbin, xbin)
      sigmaX = np.sqrt(varianceX)

      for ybin in range(0, nBinsY):
        varianceY = covariance_matrix(ybin, ybin)
        sigmaY = np.sqrt(varianceY)
        
        covXY = covariance_matrix(xbin, ybin)
        if sigmaX > 0 and sigmaY > 0:
          Cxy = covXY / (sigmaX * sigmaY)
          correlation_coefficient_matrix.SetBinContent(xbin+1, ybin+1, Cxy)

          if self.debug_level > 2:
            print('sigma x: {}, sigmay: {}'.format(sigmaX, sigmaY))
            print('cov (x,y) = {}'.format(covXY))
            print('Cxy = {}'.format(Cxy))

    output_dir = os.path.join(self.output_dir, 'Correlation_Coefficients')
    outputFilename = os.path.join(output_dir, 'hCorrelationCoefficientMatrix{}{}'.format(i, self.file_format))
    self.utils.plot_hist(correlation_coefficient_matrix, outputFilename, 'colz')


  #################################################################################################
  # Perform refolding test and statistical closure test
  #################################################################################################
  def unfolding_checks(self, jetR, obs_label, obs_setting, grooming_setting):

    # Smear MC truth spectrum according to the error bars on the measured spectrum, for closure test
    name_data = getattr(self, 'name_data_rebinned_R{}_{}'.format(jetR, obs_label))
    hData = getattr(self, name_data)
    hMC_Det = getattr(self, 'hMC_Det_R{}_{}'.format(jetR, obs_label))

    measuredErrors = self.getMeasuredErrors(hData)
    self.smearSpectrum(hMC_Det, measuredErrors)

    # Select final regularization parameter
    if self.use_max_reg_param:
      reg_param_final = self.max_reg_param
    else:
      reg_param_final = self.utils.get_reg_param(
        self.obs_settings, self.grooming_settings, self.obs_subconfig_list,
        self.obs_config_dict, obs_label, jetR)

    # Loop over values of regularization parameter to do unfolding checks
    for i in range(1, reg_param_final + 3):

      # Apply RM to unfolded result, and check that I obtain measured spectrum
      # (simple technical check)
      self.plot_refolding_test(i, jetR, obs_label, obs_setting, grooming_setting)

      # Unfold the smeared det-level result with response, and compare to truth-level MC.
      self.plot_closure_test(i, jetR, obs_label, obs_setting, grooming_setting)

  #################################################################################################
  # Apply RM to unfolded result, and check that I obtain measured spectrum (simple technical check)
  #################################################################################################
  def plot_refolding_test(self, i, jetR, obs_label, obs_setting, grooming_setting):

    hUnfolded = getattr(self, 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, obs_label, i))
    response = getattr(self, 'roounfold_response_R{}_{}'.format(jetR, obs_label))

    hFoldedTruth = response.ApplyToTruth(hUnfolded) # Produces folded distribution PerBin
    # (unfolded spectrum is also PerBin at this point)

    n_pt_bins_det = getattr(self, 'n_pt_bins_det_{}'.format(obs_label))
    det_pt_bin_array = getattr(self, 'det_pt_bin_array_{}'.format(obs_label))
    for bin in range(2, n_pt_bins_det):
      min_pt_det = det_pt_bin_array[bin]
      max_pt_det = det_pt_bin_array[bin+1]

      self.plot_obs_refolded_slice(hFoldedTruth, i, jetR, obs_label, obs_setting,
                                   grooming_setting, min_pt_det, max_pt_det)

    self.plot_pt_refolded_slice(hFoldedTruth, i, jetR, obs_label, obs_setting,
                                grooming_setting, det_pt_bin_array[2], det_pt_bin_array[-1])

  #################################################################################################
  # Plot refolding test, for a given pt slice
  #################################################################################################
  def plot_obs_refolded_slice(self, hFoldedTruth, i, jetR, obs_label, obs_setting,
                              grooming_setting, min_pt_det, max_pt_det):

    hFoldedTruth.GetXaxis().SetRangeUser(min_pt_det, max_pt_det)
    hFolded_obs = hFoldedTruth.ProjectionY()
    hFolded_obs.SetName('hFolded_obs_R{}_{}_{}_{}-{}'.format(jetR, obs_label, i,
                                                             min_pt_det, max_pt_det))

    name_data_rebinned = getattr(self, 'name_data_rebinned_R{}_{}'.format(jetR, obs_label))
    hData_PerBin = getattr(self, name_data_rebinned)
    hData_PerBin.GetXaxis().SetRangeUser(min_pt_det, max_pt_det)
    hData_obs = hData_PerBin.ProjectionY()
    hData_obs.SetName('hData_obs_R{}_{}_{}_{}-{}'.format(jetR, obs_label, i,
                                                         min_pt_det, max_pt_det))

    legendTitle = ''
    h1LegendLabel = 'Folded truth, {} = {}'.format(self.reg_param_name,i)
    h2LegendLabel = 'Measured data'
    ratioYAxisTitle = 'Folded truth / Measured'
    output_dir = getattr(self, 'output_dir_Test_Refolding')
    outf_name = 'hFoldedTruth_R{}_{}_{}-{}_{}{}'.format(self.utils.remove_periods(jetR),
                                                        obs_label, int(min_pt_det),
                                                        int(max_pt_det), i, self.file_format)
    outf_name = os.path.join(output_dir, outf_name)
    self.plot_obs_ratio(hFolded_obs, hData_obs, None, self.ytitle, ratioYAxisTitle,
                        int(min_pt_det), int(max_pt_det), jetR, obs_label, obs_setting,
                        grooming_setting, outf_name, 'width', legendTitle,
                        h1LegendLabel, h2LegendLabel)

  #################################################################################################
  # Plot refolding test, for pt dimension
  #################################################################################################
  def plot_pt_refolded_slice(self, hFoldedTruth, i, jetR, obs_label, obs_setting,
                             grooming_setting, min_pt_det, max_pt_det):

    hFoldedTruth.GetXaxis().SetRangeUser(min_pt_det, max_pt_det)
    hFolded_pt = hFoldedTruth.ProjectionX()
    hFolded_pt.SetName('hFolded_pt_R{}_{}_{}_{}-{}'.format(jetR, obs_label, i,
                                                           min_pt_det, max_pt_det))

    name_data_rebinned = getattr(self, 'name_data_rebinned_R{}_{}'.format(jetR, obs_label))
    hData_PerBin = getattr(self, name_data_rebinned)
    hData_PerBin.GetXaxis().SetRangeUser(min_pt_det, max_pt_det)
    hData_pt = hData_PerBin.ProjectionX()
    hData_pt.SetName('hData_pt_R{}_{}_{}_{}-{}'.format(jetR, obs_label, i,
                                                       min_pt_det, max_pt_det))

    legendTitle = ''
    h1LegendLabel = 'Folded truth, {} = {}'.format(self.reg_param_name,i)
    h2LegendLabel = 'Measured data'
    ratioYAxisTitle = 'Folded truth / Measured'
    output_dir = getattr(self, 'output_dir_Test_Refolding')
    outf_name = 'hFoldedTruth_pt_R{}_{}_{}{}'.format(self.utils.remove_periods(jetR),
                                                     obs_label, i, self.file_format)
    outf_name = os.path.join(output_dir, outf_name)
    self.plot_obs_ratio(hFolded_pt, hData_pt, None, self.ytitle, ratioYAxisTitle, 0, 0,
                        jetR, obs_label, obs_setting, grooming_setting, outf_name,
                        'width', legendTitle, h1LegendLabel, h2LegendLabel)

  #################################################################################################
  # Statistical closure test: Smear data, then unfold and compare to original truth
  #################################################################################################
  def plot_closure_test(self, i, jetR, obs_label, obs_setting, grooming_setting):

    # Unfold smeared det-level spectrum with RM
    response = getattr(self, 'roounfold_response_R{}_{}'.format(jetR, obs_label))
    hMC_Det = getattr(self, 'hMC_Det_R{}_{}'.format(jetR, obs_label))
    hMC_Truth = getattr(self, 'hMC_Truth_R{}_{}'.format(jetR, obs_label))

    unfold2 = ROOT.RooUnfoldBayes(response, hMC_Det, i)
    hUnfolded2 = unfold2.Hreco() # Produces the truth distribution, with errors, PerBin

    for bin in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[bin]
      max_pt_truth = self.pt_bins_reported[bin+1]

      self.plot_obs_closure_slice(hUnfolded2, hMC_Truth, i, jetR, obs_label,
                                  obs_setting, grooming_setting, min_pt_truth, max_pt_truth)

    # Closure test for pt dimension
    self.plot_pt_closure_slice(hUnfolded2, hMC_Truth, i, jetR, obs_label,
                               obs_setting, grooming_setting,
                               self.pt_bins_reported[0], self.pt_bins_reported[-1])

  #################################################################################################
  # Plot closure test, for a given pt slice
  #################################################################################################
  def plot_obs_closure_slice(self, hUnfolded, hMC_Truth, i, jetR, obs_label, obs_setting,
                             grooming_setting, min_pt_truth, max_pt_truth):

    hUnfolded.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    hUnfolded_obs = hUnfolded.ProjectionY()
    hUnfolded_obs.SetName('hUnfolded_obs_R{}_{}_{}_{}-{}'.format(jetR, obs_label, i,
                                                                 min_pt_truth, max_pt_truth))

    hMC_Truth.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    hMCTruth_obs = hMC_Truth.ProjectionY()
    hMCTruth_obs.SetName('hMCTruth_obs_R{}_{}_{}_{}-{}'.format(jetR, obs_label, i,
                                                               min_pt_truth, max_pt_truth))

    legendTitle = ''
    h1LegendLabel = 'Unfolded MC-det, {} = {}'.format(self.reg_param_name,i)
    h2LegendLabel = 'MC-truth'
    ratioYAxisTitle = 'Unfolded MC det / Truth'
    output_dir = getattr(self, 'output_dir_Test_StatisticalClosure')
    outf_name = 'hClosure_R{}_{}_{}-{}_{}{}'.format(self.utils.remove_periods(jetR),
                                                    obs_label, int(min_pt_truth),
                                                    int(max_pt_truth), i, self.file_format)
    outf_name = os.path.join(output_dir, outf_name)
    self.plot_obs_ratio(hUnfolded_obs, hMCTruth_obs, None, self.ytitle,
                        ratioYAxisTitle, min_pt_truth, max_pt_truth, jetR,
                        obs_label, obs_setting, grooming_setting, outf_name,
                        'width', legendTitle, h1LegendLabel, h2LegendLabel)

  #################################################################################################
  # Plot closure test, for pt dimension
  #################################################################################################
  def plot_pt_closure_slice(self, hUnfolded, hMC_Truth, i, jetR, obs_label,
                            obs_setting, grooming_setting, min_pt_truth, max_pt_truth):

    hUnfolded.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    hUnfolded_pt = hUnfolded.ProjectionX()
    hUnfolded_pt.SetName('hUnfolded_pt_R{}_{}_{}'.format(jetR, obs_label, i))

    hMC_Truth.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    hMCTruth_pt = hMC_Truth.ProjectionX()
    hMCTruth_pt.SetName('hMCTruth_pt_R{}_{}_{}_'.format(jetR, obs_label, i))

    legendTitle = ''
    h1LegendLabel = 'Unfolded MC-det, {} = {}'.format(self.reg_param_name,i)
    h2LegendLabel = 'MC-truth'
    ratioYAxisTitle = 'Unfolded MC det / Truth'
    output_dir = getattr(self, 'output_dir_Test_StatisticalClosure')
    outf_name = 'hClosure_pt_R{}_{}_{}{}'.format(self.utils.remove_periods(jetR),
                                                 obs_label, i, self.file_format)
    outf_name = os.path.join(output_dir, outf_name)
    self.plot_obs_ratio(hUnfolded_pt, hMCTruth_pt, None, self.ytitle, ratioYAxisTitle,
                        0, 0, jetR, obs_label, obs_setting, grooming_setting, outf_name,
                        'width', legendTitle, h1LegendLabel, h2LegendLabel)

  #################################################################################################
  # Get errors from measured spectrum, stored as dictionary {bin:error}
  #################################################################################################
  def getMeasuredErrors(self, h):

    dictErrors = {}
    for binx in range(1, h.GetNbinsX()+1):
      for biny in range(1, h.GetNbinsY()+1):
        global_bin = h.GetBin(binx, biny)
        content = h.GetBinContent(global_bin)
        error = h.GetBinError(global_bin)
        if content!=0:
          dictErrors[global_bin] = error/content
        #print('Bin ({},{}) with content {} has error {}'.format(binx, biny, content, error/content))
        else:
          print('    Error: check your histogram ranges - something might be wrong')
          print('    Bin {}, content {}, error {}'.format(global_bin, content, error))
          dictErrors[global_bin] = 0
    return dictErrors

  #################################################################################################
  # Smear spectrum according to the error bars on the measured spectrum
  #################################################################################################
  def smearSpectrum(self, h, measuredErrors):

    # Zero errors in all bins
    for binx in range(0, h.GetNbinsX() + 1):
      for biny in range(0, h.GetNbinsY() + 1):
        global_bin = h.GetBin(binx, biny)
        h.SetBinError(global_bin, 0)

    # Loop through relevant bins and smear content according to
    # errors in measured data, and set new errors
    r = ROOT.TRandom3(0)
    for bin, error in measuredErrors.items():
      content = h.GetBinContent(bin)
      errorNew = error * content
      contentNew = content + r.Gaus(0, errorNew)
      h.SetBinContent(bin, contentNew)
      h.SetBinError(bin, errorNew)
      #print('Setting bin {} to have relative error {}'.format(bin, errorNew/contentNew))

  #################################################################################################
  # Get MC-det 2D projection
  #################################################################################################
  def get_MCdet2D(self, jetR, obs_label):

    # (pt-det, pt-true, obs-det, obs-true)
    name_response = getattr(self, 'name_thn_rebinned_R{}_{}'.format(jetR, obs_label))
    hResponse = getattr(self, name_response)

    hResponse4D = hResponse.Clone()
    hResponse4D.SetName('{}_clone'.format(hResponse4D.GetName()))

    hMC_Det = hResponse4D.Projection(2,0)
    hMC_Det.SetName('hMC_Det_R{}_{}'.format(jetR, obs_label))
    return hMC_Det

  #################################################################################################
  # Get MC-det 2D projection
  #################################################################################################
  def get_MCtruth2D(self, jetR, obs_label):

    # (pt-det, pt-true, obs-det, obs-true)
    name_response = getattr(self, 'name_thn_rebinned_R{}_{}'.format(jetR, obs_label))
    hResponse = getattr(self, name_response)

    hResponse4D = hResponse.Clone()
    hResponse4D.SetName('{}_clone'.format(hResponse4D.GetName()))

    hMC_Truth = hResponse4D.Projection(3,1)
    hMC_Truth.SetName('hMC_Truth_R{}_{}'.format(jetR, obs_label))
    return hMC_Truth

  #################################################################################################
  # Plot spectra and ratio of h (and h3, if supplied) to h2
  #################################################################################################
  def plot_obs_ratio(self, h, h2, h3, yAxisTitle, ratioYAxisTitle, min_pt_det,
                     max_pt_det, jetR, obs_label, obs_setting, grooming_setting,
                     outputFilename, scalingOptions = "", legendTitle = "",
                     hLegendLabel = "", h2LegendLabel = "", h3LegendLabel = "",
                     yRatioMax = 2.2):

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    c = ROOT.TCanvas("c","c: pT",800,850)
    c.cd()
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)
    pad1.SetLeftMargin(0.2)
    pad1.SetRightMargin(0.05)
    pad1.SetTopMargin(0.05)
    if '_pt_' in outputFilename:
      pad1.SetLogy()
    pad1.Draw()
    pad1.cd()

    h.SetLineColor(1)
    h.SetLineWidth(2)
    h.SetLineStyle(1)

    h.Scale(1./h.Integral(), scalingOptions)
    if '_pt_' in outputFilename:
      h.GetYaxis().SetTitle('#frac{d#it{N}}{d#it{p}_{T}}')
    else:
      h.GetYaxis().SetTitle(yAxisTitle)

    h.GetYaxis().SetTitleSize(0.06)
    if '_pt_' in outputFilename:
      h.GetYaxis().SetRangeUser(1e-4, 10.)
    else:
      h.GetYaxis().SetRangeUser(0., 3*h.GetMaximum())

    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(20)
    xAxisTitle = self.xtitle
    h.GetXaxis().SetTitle("")

    h2.SetLineColor(4)
    h2.SetLineWidth(2)
    h2.SetLineStyle(1)
    h2.Scale(1./h2.Integral(), scalingOptions)

    h.Draw("hist same E")
    h2.Draw("hist same E")

    if h3:
      h3.SetLineColor(2)
      h3.SetLineWidth(2)
      h3.SetLineStyle(1)
      h3.Scale(1./h3.Integral(), scalingOptions)
      h3.Draw("hist same")

    c.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.2)
    pad2.SetRightMargin(0.05)
    pad2.Draw()
    pad2.cd()

    # plot ratio h/h2
    hRatio = h.Clone()
    hRatio.Divide(h2)
    hRatio.SetMarkerStyle(21)
    hRatio.SetMarkerColor(1)

    hRatio.GetXaxis().SetTitleSize(30)
    hRatio.GetXaxis().SetTitleFont(43)
    hRatio.GetXaxis().SetTitleOffset(4.)
    hRatio.GetXaxis().SetLabelFont(43)
    hRatio.GetXaxis().SetLabelSize(20)
    if '_pt_' in outputFilename:
      hRatio.GetXaxis().SetTitle('#it{p}_{T,jet}')
    else:
      hRatio.GetXaxis().SetTitle(xAxisTitle)

    hRatio.GetYaxis().SetTitle(ratioYAxisTitle)
    hRatio.GetYaxis().SetTitleSize(20)
    hRatio.GetYaxis().SetTitleFont(43)
    hRatio.GetYaxis().SetTitleOffset(2.2)
    hRatio.GetYaxis().SetLabelFont(43)
    hRatio.GetYaxis().SetLabelSize(20)
    hRatio.GetYaxis().SetNdivisions(505)
    min= hRatio.GetBinContent(hRatio.GetMinimumBin())
    max= hRatio.GetBinContent(hRatio.GetMaximumBin())
    #automatic zoom-in for a very small scatter of the points
    if min>0.5 and max<1.5:
      hRatio.GetYaxis().SetRangeUser(0.5,1.5)
    elif yRatioMax>2:
      hRatio.GetYaxis().SetRangeUser(0,yRatioMax)
    else:
      hRatio.GetYaxis().SetRangeUser(2-yRatioMax,yRatioMax)

    hRatio.Draw("P E")

    # plot ratio h3/h2
    if h3:
      hRatio3 = h3.Clone()
      hRatio3.Divide(h2)
      hRatio3.SetMarkerStyle(21)
      hRatio3.SetMarkerColor(2)
      hRatio3.Draw("P E same")

    pad1.cd()

    leg2 = ROOT.TLegend(0.55,0.7,0.8,0.93,legendTitle)
    leg2.SetFillColor(10)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetTextSize(0.04)
    leg2.AddEntry(h, hLegendLabel, "l")
    if h3:
      leg2.AddEntry(h3, h3LegendLabel, "l")
    if h2:
      leg2.AddEntry(h2, h2LegendLabel, "l")
    leg2.Draw("same")

    if not '_pt_' in outputFilename:
      text_latex = ROOT.TLatex()
      text_latex.SetNDC()
      text = str(min_pt_det) + ' < #it{p}_{T,ch jet} < ' + str(max_pt_det)
      text_latex.DrawLatex(0.25, 0.85, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = '#it{R} = ' + str(jetR)
    text_latex.DrawLatex(0.25, 0.78, text)

    subobs_label = self.utils.formatted_subobs_label(self.observable)
    if subobs_label:
      text = '{} = {}'.format(subobs_label, obs_setting)
      text_latex.DrawLatex(0.25, 0.71, text)

    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting)
      text_latex.DrawLatex(0.25, 0.64, text)

    c.SaveAs(outputFilename)
    c.Close()

#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Unfold observable distribution')
  parser.add_argument('-m', '--observable', action='store',
                      type=str, metavar='inputFileData',
                      default='theta_g',
                      help='Observable to be unfolded')
  parser.add_argument('-d', '--inputFileData', action='store',
                      type=str, metavar='inputFileData',
                      default='AnalysisResults.root',
                      help='Path of AnalysisResults.root file containing spectrum to be unfolded')
  parser.add_argument('-r', '--inputFileResponse', action='store',
                      type=str, metavar='inputFileResponse',
                      default='AnalysisResults.root',
                      help='Path of AnalysisResults.root file containing response matrix')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='analysis_config.yaml',
                      help='Path of config file for analysis')
  parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./unfolding_output/',
                      help='Output directory for plots to be written to')
  parser.add_argument('-i', '--imageFormat', action='store',
                      type=str, metavar='imageFormat',
                      default='.pdf',
                      help='Image format to save plots in, e.g. \".pdf\" or \".png\"')

  # Parse the arguments
  args = parser.parse_args()

  print('Configuring...')
  print('inputFileData: \'{0}\''.format(args.inputFileData))
  print('inputFileResponse: \'{0}\''.format(args.inputFileResponse))
  print('configFile: \'{0}\''.format(args.configFile))
  print('ouputDir: \'{0}\''.format(args.outputDir))
  print('imageFormat: \'{0}\''.format(args.imageFormat))

  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFileData):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFileData))
    sys.exit(0)
  if not os.path.exists(args.inputFileResponse):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFileResponse))
    sys.exit(0)
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = Roounfold_Obs(observable=args.observable, input_file_data = args.inputFileData,
                           input_file_response = args.inputFileResponse,
                           config_file = args.configFile, output_dir = args.outputDir,
                           file_format = args.imageFormat, rebin_response=True,
                           truncation=False, binning=False, prior_variation_parameter=0.)
  analysis.roounfold_obs()
