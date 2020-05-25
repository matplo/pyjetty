#! /usr/bin/env python

import sys
import os
import argparse
from array import *
import numpy
import ROOT
import yaml

from pyjetty.alice_analysis.analysis.user.substructure import run_analysis
from pyjetty.alice_analysis.analysis.user.james import plotting_utils

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class RunAnalysisAng(run_analysis.RunAnalysis):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(RunAnalysisAng, self).__init__(config_file, **kwargs)
    
    # Initialize yaml config
    self.initialize_user_config()

    print(self)
  
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_user_config(self):
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
      
    self.figure_approval_status = config['figure_approval_status']
    self.plot_overlay_list = self.obs_config_dict['common_settings']['plot_overlay_list']
    
    self.jet_matching_distance = config['jet_matching_distance']
    
    if 'constituent_subtractor' in config:
        self.is_pp = False
    else:
        self.is_pp = True
    print('is_pp: {}'.format(self.is_pp))

    # Theory comparisons
    if 'fPythia' in config:
      self.fPythia_name = config['fPythia']

  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #---------------------------------------------------------------
  def plot_single_result(self, jetR, obs_label, obs_setting, grooming_setting):
    print('Plotting each individual result...')
  
    # Plot final result for each 1D substructure distribution (with PYTHIA)
    self.plot_final_result(jetR, obs_label, obs_setting, grooming_setting)

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
    print('Plotting performance plots...')
    
    # Initialize performance plotting class
    if self.do_plot_performance:
      self.plotting_utils = plotting_utils.PlottingUtils(self.observable, self.is_pp,
                                                         self.main_data, self.main_response,
                                                         self.output_dir_performance,
                                                         self.figure_approval_status)
      
    # Create output subdirectories
    self.create_output_subdir(self.output_dir_performance, 'jet')
    self.create_output_subdir(self.output_dir_performance, 'resolution')
    self.create_output_subdir(self.output_dir_performance, 'residual')
    self.create_output_subdir(self.output_dir_performance, 'residual_relative')
    self.create_output_subdir(self.output_dir_performance, 'mc_projections_det')
    self.create_output_subdir(self.output_dir_performance, 'mc_projections_truth')
    self.create_output_subdir(self.output_dir_performance, 'statistics')
    self.create_output_subdir(self.output_dir_performance, 'lund')
    if not self.is_pp:
      self.create_output_subdir(self.output_dir_performance, 'prong_matching_fraction_pt')
      self.create_output_subdir(self.output_dir_performance, 'prong_matching_fraction_ptdet')
      self.create_output_subdir(self.output_dir_performance, 'prong_matching_deltaR')
      self.create_output_subdir(self.output_dir_performance, 'prong_matching_deltaZ')
      self.create_output_subdir(self.output_dir_performance, 'prong_matching_correlation')
      self.create_output_subdir(self.output_dir_performance, 'prong_matching_N_z')
    
    # Generate performance plots
    for jetR in self.jetR_list:
  
      # Plot some subobservable-independent performance plots
      self.plotting_utils.plot_DeltaR(jetR, self.jet_matching_distance)
      self.plotting_utils.plot_JES(jetR)
      self.plotting_utils.plot_JES_proj(jetR, self.pt_bins_reported)
      self.plotting_utils.plotJER(jetR, self.utils.obs_label(self.obs_settings[0], 
                                                             self.grooming_settings[0]))
      self.plotting_utils.plot_jet_reco_efficiency(jetR, self.utils.obs_label(self.obs_settings[0],
                                                                              self.grooming_settings[0]))
      
      # Plot subobservable-dependent performance plots
      for i, _ in enumerate(self.obs_subconfig_list):

        obs_setting = self.obs_settings[i]
        grooming_setting = self.grooming_settings[i]
        obs_label = self.utils.obs_label(obs_setting, grooming_setting)
    
        self.plotting_utils.plot_obs_resolution(jetR, obs_label, self.xtitle, self.pt_bins_reported)
        self.plotting_utils.plot_obs_residual(jetR, obs_label, self.xtitle, self.pt_bins_reported)
        self.plotting_utils.plot_obs_residual(jetR, obs_label, self.xtitle, self.pt_bins_reported,
                                              relative=True)
        self.plotting_utils.plot_obs_projections(jetR, obs_label, obs_setting, grooming_setting,
                                                 self.xtitle, self.pt_bins_reported)
        
        if grooming_setting and self.observable != 'jet_axis':
          self.plotting_utils.plot_lund_plane(jetR, obs_label, grooming_setting)

      # Plot prong matching histograms
      if not self.is_pp:
        self.prong_match_threshold = 0.5
        min_pt = 80.
        max_pt = 100.
        prong_list = ['leading', 'subleading']
        match_list = ['leading', 'subleading', 'groomed', 'ungroomed', 'outside']
        for i, overlay_list in enumerate(self.plot_overlay_list):
          for prong in prong_list:
            for match in match_list:

              hname = 'hProngMatching_{}_{}_JetPt_R{}'.format(prong, match, jetR)
              self.plotting_utils.plot_prong_matching(i, jetR, hname, self.obs_subconfig_list,
                                                      self.obs_settings, self.grooming_settings,
                                                      overlay_list, self.prong_match_threshold)
              self.plotting_utils.plot_prong_matching_delta(i, jetR, hname, self.obs_subconfig_list,
                                                            self.obs_settings, self.grooming_settings,
                                                            overlay_list, self.prong_match_threshold,
                                                            min_pt, max_pt, plot_deltaz=False)

              hname = 'hProngMatching_{}_{}_JetPtDet_R{}'.format(prong, match, jetR)
              self.plotting_utils.plot_prong_matching(i, jetR, hname, self.obs_subconfig_list,
                                                      self.obs_settings, self.grooming_settings,
                                                      overlay_list, self.prong_match_threshold)

              if 'subleading' in prong:
                hname = 'hProngMatching_{}_{}_JetPtZ_R{}'.format(prong, match, jetR)
                self.plotting_utils.plot_prong_matching_delta(i, jetR, hname, self.obs_subconfig_list,
                                                              self.obs_settings, self.grooming_settings,
                                                              overlay_list, self.prong_match_threshold,
                                                              min_pt, max_pt, plot_deltaz=True)

          hname = 'hProngMatching_subleading-leading_correlation_JetPtDet_R{}'.format(jetR)
          self.plotting_utils.plot_prong_matching_correlation(jetR, hname, self.obs_subconfig_list,
                                                              self.obs_settings, self.grooming_settings,
                                                              overlay_list, self.prong_match_threshold)

        # Plot subobservable-dependent plots
        for i, _ in enumerate(self.obs_subconfig_list):
          obs_setting = self.obs_settings[i]
          grooming_setting = self.grooming_settings[i]
          obs_label = self.utils.obs_label(obs_setting, grooming_setting)
          self.plotting_utils.plot_prong_N_vs_z(jetR, obs_label, 'tagged')
          self.plotting_utils.plot_prong_N_vs_z(jetR, obs_label, 'untagged')

  #----------------------------------------------------------------------
  def plot_final_result(self, jetR, obs_label, obs_setting, grooming_setting):
    print('Plot final results for {}: R = {}, {}'.format(self.observable, jetR, obs_label))

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    # Loop through pt slices, and plot final result for each 1D theta_g distribution
    for i in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[i]
      max_pt_truth = self.pt_bins_reported[i+1]
      maxbin = self.obs_max_bins(obs_label)[i]
      
      self.plot_observable(jetR, obs_label, obs_setting, grooming_setting,
                           min_pt_truth, max_pt_truth, maxbin, plot_pythia=True)

      if min_pt_truth == 40 and (jetR == 0.2 or jetR == 0.4):
        # Only want to compare to girth with \beta=1
        if obs_label == '1':
          self.plot_obs_comp(jetR, obs_label, obs_setting, grooming_setting,
                             min_pt_truth, max_pt_truth, maxbin)

  #----------------------------------------------------------------------
  def plot_observable(self, jetR, obs_label, obs_setting, grooming_setting,
                      min_pt_truth, max_pt_truth, maxbin, plot_pythia=False):
    
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
    
    xtitle = getattr(self, 'xtitle')
    ytitle = getattr(self, 'ytitle')
    color = 600-6
    
    # Get histograms
    name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label,
                                              min_pt_truth, max_pt_truth)
    if grooming_setting:
      fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
    h = self.truncate_hist(getattr(self, name), maxbin, name+'_trunc')
    h.SetMarkerSize(1.5)
    h.SetMarkerStyle(20)
    h.SetMarkerColor(color)
    h.SetLineStyle(1)
    h.SetLineWidth(2)
    h.SetLineColor(color)
    
    h_sys = getattr(self, 'hResult_{}_systotal_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth))
    h_sys.SetLineColor(0)
    h_sys.SetFillColor(color)
    h_sys.SetFillColorAlpha(color, 0.3)
    h_sys.SetFillStyle(1001)
    h_sys.SetLineWidth(0)
    
    n_obs_bins_truth = self.n_bins_truth(obs_label)
    truth_bin_array = self.truth_bin_array(obs_label)
    if maxbin:
      truth_bin_array = truth_bin_array[0:maxbin+1]
      n_obs_bins_truth = len(truth_bin_array)-1
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_obs_bins_truth, truth_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle(xtitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle(ytitle)
    myBlankHisto.SetMaximum(2.5*h.GetMaximum())
    if self.observable == 'subjet_z' or self.observable == 'jet_axis':
      myBlankHisto.SetMaximum(1.7*h.GetMaximum())
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")

    if plot_pythia:
    
      hPythia, fraction_tagged_pythia = self.pythia_prediction(
        jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin)

      if hPythia:
        hPythia.SetFillStyle(0)
        hPythia.SetMarkerSize(1.5)
        hPythia.SetMarkerStyle(21)
        hPythia.SetMarkerColor(1)
        hPythia.SetLineColor(1)
        hPythia.SetLineWidth(1)
        hPythia.Draw('E2 same')
      else:
        print('No PYTHIA prediction for {} {}'.format(self.observable, obs_label))
        plot_pythia = False
    
    h_sys.DrawCopy('E2 same')
    h.DrawCopy('PE X0 same')
  
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(0.57, 0.87, text)
    
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.8, text)

    text = str(min_pt_truth) + ' < #it{p}_{T}^{jet, ch} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.DrawLatex(0.57, 0.73, text)

    text = '#it{R} = ' + str(jetR) + '   | #eta_{jet}| < %s' % str(0.9 - jetR)
    text_latex.DrawLatex(0.57, 0.66, text)
    
    subobs_label = self.utils.formatted_subobs_label(self.observable)
    delta = 0.
    if subobs_label:
      text = '{} = {}'.format(subobs_label, obs_setting)
      text_latex.DrawLatex(0.57, 0.59, text)
      delta = 0.07
    
    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting)
      text_latex.DrawLatex(0.57, 0.59-delta, text)
      
      text_latex.SetTextSize(0.04)
      text = '#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged
      text_latex.DrawLatex(0.57, 0.52-delta, text)
    
      if plot_pythia:
        text_latex.SetTextSize(0.04)
        text = ('#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged) + (
          ', #it{f}_{tagged}^{pythia} = %3.3f' % fraction_tagged_pythia)
        text_latex.DrawLatex(0.57, 0.52-delta, text)

    myLegend = ROOT.TLegend(0.25, 0.7, 0.45, 0.85)
    self.utils.setup_legend(myLegend,0.035)
    myLegend.AddEntry(h, 'ALICE pp', 'pe')
    myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    if plot_pythia:
      myLegend.AddEntry(hPythia, 'PYTHIA8 Monash2013', 'pe')
    myLegend.Draw()

    name = 'hUnfolded_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), obs_label,
                                             int(min_pt_truth), int(max_pt_truth), self.file_format)
    if plot_pythia:
      name = 'hUnfolded_R{}_{}_{}-{}_Pythia{}'.format(self.utils.remove_periods(jetR),
                                                      obs_label, int(min_pt_truth),
                                                      int(max_pt_truth), self.file_format)
    output_dir = getattr(self, 'output_dir_final_results')
    output_dir_single = output_dir + '/single_results'
    if not os.path.exists(output_dir_single):
      os.mkdir(output_dir_single)
    outputFilename = os.path.join(output_dir_single, name)
    c.SaveAs(outputFilename)
    c.Close()

    # Write result to ROOT file
    final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
    fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
    h.Write()
    h_sys.Write()
    hPythia.Write()
    fFinalResults.Close()

  #----------------------------------------------------------------------
  # Get unfolded data from the previous preliminary result (40-60 GeV/c)
  def get_h_prelim(self, jetR):
    
    if jetR == 0.2:
      xedges = [0., 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.12]
      yvals = [2.892837, 11.9108, 17.3579, 17.65965, 15.25709, 13.00818, 9.1359, 2.471203]
      staterror = [0.09723696, 0.2879163, 0.3482209, 0.3487025, 0.3212577,
                   0.2975396, 0.2503627, 0.06595427]
      syserrorl = [0.1654749, 0.4376057, 0.536369, 0.1987916, 0.3712922,
                   0.3223265, 0.3906225, 0.1588837]
      syserrorh = [0.1673992, 0.4359767, 0.5354239, 0.200042, 0.3804927,
                   0.3368305, 0.3948841, 0.1588432]
    else:  #jetR == 0.4:
      xedges = [0, 0.02, 0.03, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.18, 0.22]
      yvals = [0.6078014, 2.815131, 5.960223, 9.770085, 9.899658, 8.603309,
               6.119539, 4.60788, 2.300467, 0.7015587]
      staterror = [0.0430097, 0.1299623, 0.1900576, 0.1710785, 0.1712931, 0.1581808,
                   0.1320032, 0.1172254, 0.059633, 0.03359087]
      syserrorl = [0.1025423, 0.1138034, 0.04146903, 0.096956, 0.06155705, 0.06077894,
                   0.08091901, 0.06775198, 0.03015912, 0.04115888]
      syserrorh = [0.1024212, 0.1204349, 0.07618093, 0.1491075, 0.07706482, 0.03761498,
                   0.1431532, 0.1033103, 0.02661073, 0.0411509]

    # Scale values by R due to observable definition
    xedges_scaled = array('d', [round(val / jetR, 2) for val in xedges])
    setattr(self, "xedges_prev_prelim_%s" % jetR, xedges_scaled)
    scale_factor = [(xedges[i+1]-xedges[i]) / (xedges_scaled[i+1]-xedges_scaled[i])
                    for i in range(0, len(yvals))]
    yvals_scaled = [yvals[i] * scale_factor[i] for i in range(0, len(yvals))]
    staterror_scaled = [staterror[i] * scale_factor[i] for i in range(0, len(staterror))]
    syserrorl_scaled = [syserrorl[i] * scale_factor[i] for i in range(0, len(syserrorl))]
    syserrorh_scaled = [syserrorh[i] * scale_factor[i] for i in range(0, len(syserrorh))]

    # Create histograms
    name = "ALI−PREL−339374"
    hCompStat = ROOT.TH1D(name, name, len(yvals), xedges_scaled)
    hCompSys = ROOT.TH1D(name+'_sys', name+'_sys', len(yvals), xedges_scaled)
    for i in range(1, len(xedges), 1):
      hCompStat.SetBinContent(i, yvals_scaled[i-1])
      hCompSys.SetBinContent(i, yvals_scaled[i-1])
      hCompStat.SetBinError(i, staterror_scaled[i-1])
      hCompSys.SetBinError(i, max(syserrorl_scaled[i-1], syserrorh_scaled[i-1]))

    return hCompStat, hCompSys

  #----------------------------------------------------------------------
  def plot_obs_comp(self, jetR, obs_label, obs_setting, grooming_setting,
                    min_pt_truth, max_pt_truth, maxbin):
    
    # Scale both distributions by integrals
    scale_by_int = False

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
    
    xtitle = getattr(self, 'xtitle')
    ytitle = getattr(self, 'ytitle')
    color = 600-6
    
    # Get histograms
    name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label,
                                              min_pt_truth, max_pt_truth)
    if grooming_setting:
      fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
    h = self.truncate_hist(getattr(self, name), maxbin, name+'_trunc')
    h.SetMarkerSize(1.5)
    h.SetMarkerStyle(20)
    h.SetMarkerColor(color)
    h.SetLineStyle(1)
    h.SetLineWidth(2)
    h.SetLineColor(color)
    if scale_by_int:
      h.Scale(1/h.Integral())
    
    h_sys = getattr(self, 'hResult_{}_systotal_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth))
    h_sys.SetLineColor(0)
    h_sys.SetFillColor(color)
    h_sys.SetFillColorAlpha(color, 0.3)
    h_sys.SetFillStyle(1001)
    h_sys.SetLineWidth(0)
    if scale_by_int:
      h_sys.Scale(1/h_sys.Integral())
    
    n_obs_bins_truth = self.n_bins_truth(obs_label)
    truth_bin_array = self.truth_bin_array(obs_label)
    if maxbin:
      truth_bin_array = truth_bin_array[0:maxbin+1]
      if jetR == 0.2:
        truth_bin_array[-1] = 0.6
      n_obs_bins_truth = len(truth_bin_array)-1
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_obs_bins_truth, truth_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle(xtitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle(ytitle)
    if scale_by_int:
      myBlankHisto.SetYTitle(ytitle + ' / #int(' + ytitle + ')')
    myBlankHisto.SetMaximum(2*h.GetMaximum())
    if self.observable == 'subjet_z' or self.observable == 'jet_axis':
      myBlankHisto.SetMaximum(1.5*h.GetMaximum())
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")

    hCompStat, hCompSys = self.get_h_prelim(jetR)

    #formatting
    hCompStat.SetFillStyle(0)
    hCompStat.SetMarkerSize(1.5)
    hCompStat.SetMarkerStyle(21)
    hCompStat.SetMarkerColor(1)
    hCompStat.SetLineColor(1)
    hCompStat.SetLineWidth(1)
    hCompSys.SetLineColor(0)
    hCompSys.SetFillColor(1)
    hCompSys.SetFillColorAlpha(1, 0.3)
    hCompSys.SetFillStyle(1001)
    hCompSys.SetLineWidth(0)
    if scale_by_int:
      hCompStat.Scale(1/hCompStat.Integral())
      hCompSys.Scale(1/hCompSys.Integral())

    hCompSys.Draw('E2 same')
    hCompStat.Draw('PE X0 same')
    h_sys.DrawCopy('E2 same')
    h.DrawCopy('PE X0 same')
  
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(0.57, 0.87, text)
    
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.8, text)

    text = str(min_pt_truth) + ' < #it{p}_{T}^{jet, ch} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.DrawLatex(0.57, 0.73, text)

    text = '#it{R} = ' + str(jetR) + '   | #eta_{jet}| < %s' % str(0.9 - jetR)
    text_latex.DrawLatex(0.57, 0.66, text)
    
    subobs_label = self.utils.formatted_subobs_label(self.observable)
    delta = 0.
    if subobs_label:
      text = '{} = {}'.format(subobs_label, obs_setting)
      text_latex.DrawLatex(0.57, 0.59, text)
      delta = 0.07
    
    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting)
      text_latex.DrawLatex(0.57, 0.59-delta, text)
      
      text_latex.SetTextSize(0.04)
      text = '#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged
      text_latex.DrawLatex(0.57, 0.52-delta, text)

    myLegend = ROOT.TLegend(0.25, 0.7, 0.45, 0.85)
    self.utils.setup_legend(myLegend,0.035)
    myLegend.AddEntry(h, 'This measurement', 'pe')
    myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    myLegend.AddEntry(hCompStat, 'ALI-PREL-339374', 'pe')
    myLegend.Draw()

    name = 'hUnfoldedComp_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR),
                                                 obs_label, int(min_pt_truth),
                                                 int(max_pt_truth), self.file_format)
    output_dir = getattr(self, 'output_dir_final_results')
    output_dir_single = output_dir + '/single_results'
    if not os.path.exists(output_dir_single):
      os.mkdir(output_dir_single)
    outputFilename = os.path.join(output_dir_single, name)
    c.SaveAs(outputFilename)
    c.Close()

    # Write result to ROOT file
    final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
    fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
    h.Write()
    h_sys.Write()
    hCompStat.Write()
    hCompSys.Write()
    fFinalResults.Close()

  #----------------------------------------------------------------------
  def pythia_prediction(self, jetR, obs_setting, obs_label, min_pt_truth,
                        max_pt_truth, maxbin, overlay=False):
  
    hPythia = self.get_pythia_from_response(jetR, obs_label, min_pt_truth,
                                            max_pt_truth, maxbin, overlay)
    n_jets_inclusive = hPythia.Integral(0, hPythia.GetNbinsX()+1)
    n_jets_tagged = hPythia.Integral(hPythia.FindBin(
      self.truth_bin_array(obs_label)[0]), hPythia.GetNbinsX())

    fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
    hPythia.Scale(1./n_jets_inclusive, 'width')
      
    return [hPythia, fraction_tagged_pythia]

  #----------------------------------------------------------------------
  def get_pythia_from_response(self, jetR, obs_label, min_pt_truth, max_pt_truth,
                               maxbin, overlay=False):

    output_dir = getattr(self, 'output_dir_main')

    prev_prelim = False
    if overlay and (jetR == 0.2 or jetR == 0.4) and min_pt_truth == 40 and obs_label == '1':
      prev_prelim = True
      # Need to rebin response for the binning used by previous preliminary result
      filepath = os.path.join(output_dir, 'response_prev_prelim.root')

      if not os.path.exists(filepath):
        # Create rebinned THn with these binnings, and write to file
        print("Rebinning response matrix for previous preliminary masurement...")
        name_thn = self.utils.name_thn(self.observable, jetR, obs_label)
        name_thn_rebinned = self.utils.name_thn_rebinned(self.observable, jetR, obs_label)
        name_roounfold = 'roounfold_response_R{}_{}'.format(jetR, obs_label)
        thn = ROOT.TFile(self.main_response, 'READ').Get(name_thn)
        thn.SetName(name_thn)
        label = 'R{}_{}'.format(jetR, obs_label)
        pt_bins_truth = array('d', [5, 20, 40, 60, 80, 100, 150, 200])
        pt_bins_det = array('d', [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 120, 150])
        obs_bins = getattr(self, "xedges_prev_prelim_%s" % jetR)
        self.utils.rebin_response(
          filepath, thn, name_thn_rebinned, name_roounfold, label, len(pt_bins_det)-1,
          pt_bins_det, len(obs_bins)-1, obs_bins, len(pt_bins_truth)-1, pt_bins_truth,
          len(obs_bins)-1, obs_bins, self.observable, do_roounfoldresponse=False)
    else:
      filepath = os.path.join(output_dir, 'response.root')
    f = ROOT.TFile(filepath, 'READ')

    thn_name = 'hResponse_JetPt_{}_R{}_{}_rebinned'.format(self.observable, jetR, obs_label)
    thn = f.Get(thn_name)
    thn.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)

    name = 'hPythia_{}_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    if prev_prelim:
      h = thn.Projection(3)
    else:
      h = self.truncate_hist(thn.Projection(3), maxbin, name)
    h.SetDirectory(0)

    return h

  #----------------------------------------------------------------------
  def plot_final_result_overlay(self, i_config, jetR, overlay_list):
    print('Plotting overlay of {}'.format(overlay_list))

    # Plot overlay of different subconfigs, for fixed pt bin
    for i in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[i]
      max_pt_truth = self.pt_bins_reported[i+1]
      maxbins = [self.obs_max_bins(obs_label)[i] for obs_label in self.obs_labels]

      # Plot PYTHIA
      self.plot_observable_overlay_subconfigs(
        i_config, jetR, overlay_list, min_pt_truth,
        max_pt_truth, maxbins, plot_pythia=True, plot_ratio = True)


  #----------------------------------------------------------------------
  def plot_observable_overlay_subconfigs(self, i_config, jetR, overlay_list, min_pt_truth,
                                         max_pt_truth, maxbins, plot_pythia=False,
                                         plot_nll=False, plot_ratio=False):

    # Flag to plot ratio all on the same scale, 0 to 2.2
    plot_ratio_same_scale = True

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

    myLegend = ROOT.TLegend(0.66, 0.65, 0.8, 0.85)
    self.utils.setup_legend(myLegend, 0.035)
    
    name = 'hmain_{}_R{}_{{}}_{}-{}'.format(self.observable, jetR, min_pt_truth, max_pt_truth)
    ymax = self.get_maximum(name, overlay_list)

    h_list = []
    text_list = []

    for i, subconfig_name in enumerate(self.obs_subconfig_list):
    
      if subconfig_name not in overlay_list:
        continue

      obs_setting = self.obs_settings[i]
      grooming_setting = self.grooming_settings[i]
      obs_label = self.obs_labels[i]
      maxbin = maxbins[i]
      
      if subconfig_name == overlay_list[0]:
        marker = 20
        marker_pythia = marker+4
        color = 1
      elif subconfig_name == overlay_list[1]:
        marker = 21
        marker_pythia = marker+4
        color = 600-6
      elif subconfig_name == overlay_list[2]:
        marker = 22
        marker_pythia = marker+4
        color = 632-4
      else:  # subconfig_name == overlay_list[3]:
        marker = 23
        marker_pythia = 32
        color = 416-2

      name = 'hmain_{}_R{}_{}_{}-{}'.format(
        self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
      if grooming_setting:
        fraction_tagged = getattr(self, name+'_fraction_tagged')

      if (jetR == 0.2 or jetR == 0.4) and min_pt_truth == 40 and obs_label == '1':
        # Use previous preliminary result
        h, h_sys = self.get_h_prelim(jetR)
        # Move error bars to different histogram
        h_sys.SetNameTitle("hSysPrelim%s", "hSysPrelim%s")
        for i in range(1, h.GetNbinsX()+1, 1):
          h.SetBinError(i, 0.)
      else:
        h = self.truncate_hist(getattr(self, name), maxbin, name+'_trunc')
        h_sys = getattr(self, 'hResult_{}_systotal_R{}_{}_{}-{}'.format(
          self.observable, jetR, obs_label, min_pt_truth, max_pt_truth))

      h.SetDirectory(0)
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(marker)
      h.SetMarkerColor(color)
      h.SetLineStyle(1)
      h.SetLineWidth(2)
      h.SetLineColor(color)
      
      h_sys.SetLineColor(0)
      h_sys.SetFillColor(color)
      h_sys.SetFillColorAlpha(color, 0.3)
      h_sys.SetFillStyle(1001)
      h_sys.SetLineWidth(0)
      
      if subconfig_name == overlay_list[0]:

        pad1.cd()
        xtitle = getattr(self, 'xtitle')
        ytitle = getattr(self, 'ytitle')
        xmin = self.obs_config_dict[subconfig_name]['obs_bins_truth'][0]
        xmax = self.obs_config_dict[subconfig_name]['obs_bins_truth'][-1]
        if maxbin:
          if (jetR == 0.2 or jetR == 0.4) and min_pt_truth == 40 and obs_label == '1':
            xmax = getattr(self, "xedges_prev_prelim_%s" % jetR)[-1]
          else:
            xmax = self.obs_config_dict[subconfig_name]['obs_bins_truth'][maxbin]
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, xmin, xmax)
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(xtitle)
        myBlankHisto.GetYaxis().SetTitleOffset(1.3)
        myBlankHisto.SetYTitle(ytitle)
        if jetR == 0.2:
          if min_pt_truth == 20:
            myBlankHisto.SetMaximum(1.1*ymax)
          elif min_pt_truth == 40:
            myBlankHisto.SetMaximum(1.2*ymax)
          elif min_pt_truth == 60:
            myBlankHisto.SetMaximum(1.15*ymax)
          else:
            myBlankHisto.SetMaximum(1.3*ymax)
        elif jetR == 0.4:
          if min_pt_truth == 20:
            myBlankHisto.SetMaximum(1.5*ymax)
          elif min_pt_truth == 40:
            myBlankHisto.SetMaximum(1.4*ymax)
          elif min_pt_truth == 60:
            myBlankHisto.SetMaximum(1.15*ymax)
          else:
            myBlankHisto.SetMaximum(1.5*ymax)
        else:
          myBlankHisto.SetMaximum(1.5*ymax)
        myBlankHisto.SetMinimum(0.)
        if plot_ratio:
          myBlankHisto.SetMinimum(2e-4) # Don't draw 0 on top panel
          myBlankHisto.GetYaxis().SetTitleSize(0.065)
          myBlankHisto.GetYaxis().SetTitleOffset(1.1)
          myBlankHisto.GetYaxis().SetLabelSize(0.06)
        myBlankHisto.Draw('E')
        
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
          if plot_ratio_same_scale:
            if jetR == 0.2:
              myBlankHisto2.GetYaxis().SetRangeUser(0.5, 1.75)
            elif jetR == 0.4:
              myBlankHisto2.GetYaxis().SetRangeUser(0, 1.99)
            else: 
              myBlankHisto2.GetYaxis().SetRangeUser(0, 2.2)
          elif jetR == 0.2:
            if min_pt_truth == 20:
              myBlankHisto2.GetYaxis().SetRangeUser(0.6, 1.75)
            elif min_pt_truth == 40:
              myBlankHisto2.GetYaxis().SetRangeUser(0.78, 1.299)
            elif min_pt_truth == 60:
              myBlankHisto2.GetYaxis().SetRangeUser(0.55, 1.499)
            else: 
              myBlankHisto2.GetYaxis().SetRangeUser(0.5, 1.99)
          elif jetR == 0.4:
            if min_pt_truth == 20:
              myBlankHisto2.GetYaxis().SetRangeUser(0.81, 1.72)
            elif min_pt_truth == 40:
              myBlankHisto2.GetYaxis().SetRangeUser(0.7, 2.1)
            elif min_pt_truth == 60:
              myBlankHisto2.GetYaxis().SetRangeUser(0.75, 1.55)
            else: 
              myBlankHisto2.GetYaxis().SetRangeUser(0.5, 1.99)
          else: 
            myBlankHisto2.GetYaxis().SetRangeUser(0.5, 1.99)
          myBlankHisto2.Draw()
        
          line = ROOT.TLine(0,1,xmax,1)
          line.SetLineColor(920+2)
          line.SetLineStyle(2)
          line.Draw()
      
      if plot_pythia:
        hPythia, fraction_tagged_pythia = self.pythia_prediction(
          jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin, overlay=True)

        plot_errors = False
        if plot_errors:
          hPythia.SetMarkerSize(0)
          hPythia.SetMarkerStyle(0)
          hPythia.SetMarkerColor(color)
          hPythia.SetFillColor(color)
        else:
          hPythia.SetLineColor(color)
          hPythia.SetLineColorAlpha(color, 0.5)
          hPythia.SetLineWidth(4)

      if plot_ratio:
        hRatioSys = h_sys.Clone()
        hRatioSys.SetName('{}_Ratio'.format(h_sys.GetName()))
        if plot_pythia:
          hRatioSys.Divide(hPythia)
        hRatioSys.SetLineColor(0)
        hRatioSys.SetFillColor(color)
        hRatioSys.SetFillColorAlpha(color, 0.3)
        hRatioSys.SetFillStyle(1001)
        hRatioSys.SetLineWidth(0)
        hRatioSys.SetMaximum(1.99)
          
        hRatioStat = h.Clone()
        hRatioStat.SetName('{}_Ratio'.format(h.GetName()))
        if plot_pythia:
          hRatioStat.Divide(hPythia)
        hRatioStat.SetMarkerSize(1.5)
        hRatioStat.SetMarkerStyle(marker)
        hRatioStat.SetMarkerColor(color)
        hRatioStat.SetLineStyle(1)
        hRatioStat.SetLineWidth(2)
        hRatioStat.SetLineColor(color)
        hRatioStat.SetMaximum(1.99)

      pad1.cd()
      if plot_pythia:
        plot_errors = False
        if plot_errors:
          hPythia.DrawCopy('E3 same')
        else:
          hPythia.DrawCopy('L hist same')

      h_sys.DrawCopy('E2 same')
      h.DrawCopy('PE X0 same')
      
      if plot_ratio:
        pad2.cd()
        if plot_pythia:
          hRatioSys.DrawCopy('E2 same')
          hRatioStat.DrawCopy('PE X0 same')

      subobs_label = self.utils.formatted_subobs_label(self.observable)
      text = ''
      if subobs_label:
        text += '{} = {}'.format(subobs_label, obs_setting)
      if grooming_setting:
        text += self.utils.formatted_grooming_label(grooming_setting)
      text_list.append(text)
      h_list.append(h)
        
    pad1.cd()
    for h, text in zip(h_list, text_list):
      myLegend.AddEntry(h, text, 'pe')
    myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    if plot_pythia:
      myLegend.AddEntry(hPythia, 'PYTHIA8 Monash2013', 'l')

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(0.25, 0.87, text)
    
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.81, text)

    text = 'Charged jets   anti-#it{k}_{T}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.75, text)
    
    text = '#it{R} = ' + str(jetR) + '   | #eta_{jet}| < %s' % str(0.9 - jetR)
    text_latex.DrawLatex(0.25, 0.69, text)
    
    text = str(min_pt_truth) + ' < #it{p}_{T}^{jet, ch} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.63, text)
    
    myLegend.Draw()

    name = 'h_{}_R{}_{}-{}_{}{}'.format(self.observable, 
                                        self.utils.remove_periods(jetR), int(min_pt_truth), 
                                        int(max_pt_truth), i_config, self.file_format)
    if plot_pythia:
      name = 'h_{}_R{}_{}-{}_Pythia_{}{}'.format(self.observable, self.utils.remove_periods(jetR),
                                                 int(min_pt_truth), int(max_pt_truth),
                                                 i_config, self.file_format)
    output_dir = getattr(self, 'output_dir_final_results') + '/all_results'
    if not os.path.exists(output_dir):
      os.mkdir(output_dir)
    outputFilename = os.path.join(output_dir, name)
    c.SaveAs(outputFilename)
    c.Close()

  #----------------------------------------------------------------------
  # Return maximum y-value of unfolded results in a subconfig list
  def get_maximum(self, name, overlay_list):
  
    max = 0.
    for i, subconfig_name in enumerate(self.obs_subconfig_list):
    
      if subconfig_name not in overlay_list:
        continue

      obs_setting = self.obs_settings[i]
      grooming_setting = self.grooming_settings[i]
      obs_label = self.utils.obs_label(obs_setting, grooming_setting)
      
      h = getattr(self, name.format(obs_label))
      if h.GetMaximum() > max:
        max = h.GetMaximum()
        
    return max

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

  analysis = RunAnalysisAng(config_file = args.configFile)
  analysis.run_analysis()
