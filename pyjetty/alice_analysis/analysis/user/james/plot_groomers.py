#!/usr/bin/env python3

"""
  Plotting utilities for jet substructure analysis with track dataframe.
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import math
import yaml
import argparse

# Data analysis and plotting
import numpy as np
from array import *
import ROOT

# Base class
from pyjetty.alice_analysis.analysis.base import common_base
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs
from pyjetty.alice_analysis.analysis.user.james import plotting_utils

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

# Suppress a lot of standard output
ROOT.gErrorIgnoreLevel = ROOT.kWarning

################################################################
class PlotGroomers(common_base.CommonBase):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, output_dir = '.', config_file = '', **kwargs):
    super(PlotGroomers, self).__init__(**kwargs)
    
    # Initialize utils class
    self.utils = analysis_utils_obs.AnalysisUtils_Obs()
    
    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()
    
    # Read config file
    self.config_file = config_file
    with open(config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    
    self.observables = config['process_observables']
    thermal_data = config['main_response']
    self.fMC = ROOT.TFile(thermal_data, 'READ')
    
    self.max_distance = config['constituent_subtractor']['max_distance']
    self.R_max = config['constituent_subtractor']['main_R_max']
    
    self.file_format = config['file_format']
    self.output_dir = config['output_dir']

    self.ColorArray = [ROOT.kBlue-4, ROOT.kAzure+7, ROOT.kCyan-2, ROOT.kViolet-8,
                       ROOT.kBlue-6, ROOT.kGreen+3, ROOT.kPink-4, ROOT.kRed-4,
                       ROOT.kOrange-3]
    self.MarkerArray = [20, 21, 22, 23, 33, 34, 24, 25, 26, 32]
    self.OpenMarkerArray = [24, 25, 26, 32, 27, 28]
    
    print(self)

  #---------------------------------------------------------------
  def init_observable(self, observable):
    
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    
    # Get the sub-configs
    self.jetR_list = config['jetR']
    self.obs_config_dict = config[observable]
    self.obs_subconfig_list = [name for name in list(self.obs_config_dict.keys()) if 'config' in name ]
    self.grooming_settings = self.utils.grooming_settings(self.obs_config_dict)
    self.obs_settings = self.utils.obs_settings(observable, self.obs_config_dict, self.obs_subconfig_list)
    self.obs_labels = [self.utils.obs_label(self.obs_settings[i], self.grooming_settings[i])
                       for i in range(len(self.obs_subconfig_list))]
    self.xtitle = self.obs_config_dict['common_settings']['xtitle']
    self.ytitle = self.obs_config_dict['common_settings']['ytitle']
    self.pt_bins_reported = self.obs_config_dict['common_settings']['pt_bins_reported']
    self.plot_overlay_list = self.obs_config_dict['common_settings']['plot_overlay_list']

    # Create output dirs
    output_dir = os.path.join(self.output_dir, observable)
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

  #---------------------------------------------------------------
  def plot_groomers(self):
  
    # Loop through all observables
    for observable in self.observables:
      self.init_observable(observable)
  
      # Plot for each R_max
      for R_max in self.max_distance:
      
        # Plot for each R
        for jetR in self.jetR_list:
        
          # Create output dir
          output_dir = os.path.join(self.output_dir, observable)
          output_dir = os.path.join(output_dir, 'jetR{}'.format(jetR))
          output_dir = os.path.join(output_dir, 'Rmax{}'.format(R_max))
          self.create_output_subdir(output_dir, 'money_plot')
        
          # Plot money plot for all observables
          for i, _ in enumerate(self.obs_subconfig_list):

            obs_setting = self.obs_settings[i]
            grooming_setting = self.grooming_settings[i]
            obs_label = self.utils.obs_label(obs_setting, grooming_setting)

            #self.plot_money_plot(jetR, obs_label, obs_setting, grooming_setting, self.xtitle, self.pt_bins_reported)

          # Plot performance plots only once
          if observable == 'theta_g':

            # Create output subdirectories
            output_dir = os.path.join(self.output_dir, 'performance')
            output_dir = os.path.join(output_dir, 'jetR{}'.format(jetR))
            output_dir = os.path.join(output_dir, 'Rmax{}'.format(R_max))
            self.create_output_subdir(output_dir, 'delta_pt')
            self.create_output_subdir(output_dir, 'prong_matching_fraction_pt')
            self.create_output_subdir(output_dir, 'prong_matching_deltaR')
            self.create_output_subdir(output_dir, 'prong_matching_deltaZ')
            self.create_output_subdir(output_dir, 'prong_matching_correlation')
            
            self.plotting_utils = plotting_utils.PlottingUtils(output_dir, self.config_file, R_max=R_max,
                                                               thermal = False, groomer_studies = True)
        
            # Plot some subobservable-independent performance plots
            self.plotting_utils.plot_delta_pt(jetR, self.pt_bins_reported)

            # Plot prong matching histograms
            self.prong_match_threshold = 0.5
            min_pt = 80.
            max_pt = 100.
            prong_list = ['leading', 'subleading']
            match_list = ['leading', 'subleading', 'ungroomed', 'outside']
            for i, overlay_list in enumerate(self.plot_overlay_list):
              for prong in prong_list:
                for match in match_list:

                  hname = 'hProngMatching_{}_{}_JetPt_R{}'.format(prong, match, jetR)
                  self.plotting_utils.plot_prong_matching(i, jetR, hname, self.obs_subconfig_list, self.obs_settings, self.grooming_settings, overlay_list, self.prong_match_threshold)
                  self.plotting_utils.plot_prong_matching_delta(i, jetR, hname, self.obs_subconfig_list, self.obs_settings, self.grooming_settings, overlay_list, self.prong_match_threshold, min_pt, max_pt, plot_deltaz=False)

                  if 'subleading' in prong or 'leading' in prong:
                    hname = 'hProngMatching_{}_{}_JetPtZ_R{}'.format(prong, match, jetR)
                    self.plotting_utils.plot_prong_matching_delta(i, jetR, hname, self.obs_subconfig_list, self.obs_settings, self.grooming_settings, overlay_list, self.prong_match_threshold, min_pt, max_pt, plot_deltaz=True)

              hname = 'hProngMatching_subleading-leading_correlation_JetPt_R{}'.format(jetR)
              self.plotting_utils.plot_prong_matching_correlation(i, jetR, hname, self.obs_subconfig_list, self.obs_settings, self.grooming_settings, overlay_list, self.prong_match_threshold)

  #---------------------------------------------------------------
  def plot_obs_projections(self, jetR, R_max, obs_label, obs_setting, grooming_setting, xtitle, pt_bins):
    return
    if not self.fData:
      return
      
    # (pt-det, pt-truth, obs-det, obs-truth)
    name = 'hResponse_JetPt_{}_R{}_{}_Rmax{}'.format(self.observable, jetR, obs_label, R_max)
    hRM_obs = self.fMC.Get(name)
    if hRM_obs.GetSumw2() is 0:
      hRM_obs.Sumw2()
    
    name = 'h_{}_JetPt_R{}_{}_Rmax{}'.format(self.observable, jetR, obs_label, R_max)
    hObs_JetPt = self.fData.Get(name)
    if hObs_JetPt.GetSumw2() is 0:
      hObs_JetPt.Sumw2()

    # Plot 2D statistics in data
    self.plot2D_obs_statistics(hObs_JetPt.Clone(), jetR, obs_label)
      
    # Loop through pt slices, and plot:
    #   (a) MC-det and MC-truth 1D distributions, for fixed pt-truth
    #   (b) MC-det and data 1D distributions, for fixed pt-det
    for i in range(0, len(pt_bins) - 1):
      min_pt_truth = pt_bins[i]
      max_pt_truth = pt_bins[i+1]
      
      self.plot_obs_projection(hRM_obs, hObs_JetPt, jetR, obs_label, obs_setting, grooming_setting, xtitle, min_pt_truth, max_pt_truth, option='truth')
      self.plot_obs_projection(hRM_obs, hObs_JetPt, jetR, obs_label, obs_setting, grooming_setting, xtitle, min_pt_truth, max_pt_truth, option='det')

  #---------------------------------------------------------------
  # If option='truth', plot MC-truth and MC-det projections for fixed pt-true
  # If option='det', plot data and MC-det projections for fixed pt-det
  def plot_obs_projection(self, hRM, hObs_JetPt, jetR, obs_label, obs_setting, grooming_setting, xtitle, min_pt, max_pt, option='truth'):

    ytitle = '#frac{{1}}{{N}} #frac{{dN}}{{d{}}}'.format(xtitle)
    
    if self.observable == 'theta_g':
      rebin_val_mcdet = 5
      rebin_val_mctruth = 5
      rebin_val_data = 5
    elif self.observable == 'zg':
      rebin_val_mcdet = 5
      rebin_val_mctruth = 2
      rebin_val_data = 10
    elif self.observable == 'subjet_z':
      rebin_val_mcdet = 2
      rebin_val_mctruth = 1
      rebin_val_data = 2
    elif self.observable == 'jet_axis':
      rebin_val_mcdet = 2
      rebin_val_mctruth = 1
      rebin_val_data = 5
      
    # Get RM, for a given pt cut
    if option == 'det':
      hRM.GetAxis(0).SetRangeUser(min_pt, max_pt)
    if option == 'truth':
      hRM.GetAxis(1).SetRangeUser(min_pt, max_pt)
    
    # Get histogram of observable at MC-det from RM
    hObs_det = hRM.Projection(2)
    hObs_det.SetName('hObs_det_{}'.format(obs_label))
    hObs_det.GetYaxis().SetTitle(ytitle)
    hObs_det.SetLineColor(2)
    hObs_det.SetLineWidth(2)
    self.scale_by_integral(hObs_det)
    hObs_det.Rebin(rebin_val_mcdet)
    hObs_det.Scale(1., 'width')
    if 'sd' in grooming_setting:
      hObs_det.GetXaxis().SetRange(0, hObs_det.GetNbinsX())

    # Get histogram of observable at MC-truth from RM
    if option == 'truth':
      hObs_truth = hRM.Projection(3)
      hObs_truth.SetName('hObs_truth_{}'.format(obs_label))
      hObs_truth.SetLineColor(4)
      hObs_truth.SetLineWidth(2)
      self.scale_by_integral(hObs_truth)
      hObs_truth.Rebin(rebin_val_mctruth)
      hObs_truth.Scale(1., 'width')
      if 'sd' in grooming_setting:
        hObs_truth.GetXaxis().SetRange(0, hObs_truth.GetNbinsX())
      
    # Get histogram of theta_g in data, for given pt-det cut
    if option == 'det':
      hObs_JetPt.GetXaxis().SetRangeUser(min_pt, max_pt)
      hObs_data = hObs_JetPt.ProjectionY()
      hObs_data.SetMarkerStyle(21)
      hObs_data.SetMarkerSize(1)
      self.scale_by_integral(hObs_data)
      hObs_data.Rebin(rebin_val_data)
      hObs_data.Scale(1., 'width')
      if 'sd' in grooming_setting:
        hObs_data.GetXaxis().SetRange(0, hObs_data.GetNbinsX())

    # Draw histogram
    c = ROOT.TCanvas('c','c: hist',600,450)
    c.cd()

    myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    myPad.SetLeftMargin(0.2)
    myPad.SetTopMargin(0.07)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.13)
    myPad.Draw()
    myPad.cd()
    
    leg = ROOT.TLegend(0.7,0.75,0.85,0.85, "")
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1)
    leg.SetTextSize(0.04)
    
    hObs_det.GetYaxis().SetTitleOffset(1.5)
    hObs_det.SetMaximum(2.5*hObs_det.GetMaximum())
    hObs_det.SetMinimum(0.)

    hObs_det.Draw('hist')
    leg.AddEntry(hObs_det, "MC det", "L")
    if option == 'truth':
      hObs_truth.Draw('hist same')
      leg.AddEntry(hObs_truth, "MC truth", "L")
    elif option == 'det':
      hObs_data.Draw('hist E same')
      leg.AddEntry(hObs_data, "data", "PE")

    leg.Draw("same")
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(0.3, 0.85, text)
    
    if self.is_pp:
      text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    else:
      text = 'Pb-Pb #sqrt{#it{s_{NN}}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.3, 0.79, text)

    if option == 'truth':
      text = str(min_pt) + ' < #it{p}_{T, ch jet}^{truth} < ' + str(max_pt) + ' GeV/#it{c}'
    elif option == 'det':
      text = str(min_pt) + ' < #it{p}_{T, ch jet}^{det} < ' + str(max_pt) + ' GeV/#it{c}'
    text_latex.DrawLatex(0.3, 0.73, text)

    text = '#it{R} = ' + str(jetR) + '   | #eta_{jet}| < 0.5'
    text_latex.DrawLatex(0.3, 0.67, text)
    
    subobs_label = self.formatted_subobs_label(self.observable)
    delta = 0.
    if subobs_label:
      text = '{} = {}'.format(subobs_label, obs_setting)
      text_latex.DrawLatex(0.3, 0.61, text)
      delta = 0.07
      
    if grooming_setting:
      text = self.formatted_grooming_label(grooming_setting)
      text_latex.DrawLatex(0.3, 0.61-delta, text)

    output_filename = os.path.join(self.output_dir, 'mc_projections_{}/h_{}_MC_R{}_{}_{}-{}.pdf'.format(option, self.observable, self.remove_periods(jetR), obs_label, min_pt, max_pt))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  # Create a single output subdirectory
  #---------------------------------------------------------------
  def create_output_subdir(self, output_dir, name):

    output_subdir = os.path.join(output_dir, name)
    setattr(self, 'output_dir_{}'.format(name), output_subdir)
    if not os.path.isdir(output_subdir):
      os.makedirs(output_subdir)

    return output_subdir

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

  analysis = PlotGroomers(config_file = args.configFile)
  analysis.plot_groomers()
