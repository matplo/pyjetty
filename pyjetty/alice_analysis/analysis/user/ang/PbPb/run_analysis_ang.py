#! /usr/bin/env python

""" theory_comp.py
Loads theory comparisons, preforms un/folding, makes plots
Ezra Lesser, 2020 (elesser@berkeley.edu)
"""

import sys
import os
import argparse
from array import *
import numpy as np
import math   # for exp()
import ROOT
ROOT.gSystem.Load("$HEPPY_DIR/external/roounfold/roounfold-current/lib/libRooUnfold.so")
import yaml

# For log y plots which ROOT just decides not to work for
#import matplotlib
#matplotlib.rcParams['text.usetex'] = True   # LaTeX labels
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams["yaxis.labellocation"] = 'top'
plt.rcParams["xaxis.labellocation"] = 'right'

from pyjetty.alice_analysis.analysis.user.substructure import run_analysis
from pyjetty.alice_analysis.analysis.user.ang.PbPb import plotting_utils_ang

# Load pyjetty ROOT utils
ROOT.gSystem.Load('libpyjetty_rutil')

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
#######################  RUN ANALYSIS  #########################
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
    self.plot_overlay_list = \
      self.obs_config_dict['common_settings']['plot_overlay_list']

    self.jet_matching_distance = config['jet_matching_distance']

    self.is_pp = True
    self.results_pp = None
    if 'constituent_subtractor' in config:
      self.is_pp = False
      self.results_pp = config["results_pp"] if "results_pp" in config else None
      self.max_distance = config["constituent_subtractor"]["max_distance"]
    print('is_pp: {}'.format(self.is_pp))

    # Whether or not to use the previous measurement in ratio
    self.use_prev_result = config["use_prev_result"]

    self.histutils = ROOT.RUtil.HistUtils()


  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #---------------------------------------------------------------
  def plot_single_result(self, jetR, obs_label, obs_setting, grooming_setting):
    #print('Plotting each individual result...')

    # Plot final result for each 1D substructure distribution (with PYTHIA)
    self.plot_final_result(jetR, obs_label, obs_setting, grooming_setting)


  #---------------------------------------------------------------
  # This function is called once after all subconfigurations have been looped over, for each R
  #---------------------------------------------------------------
  def plot_all_results(self, jetR):

    print('Plotting overlay of all results...')

    for i_config, overlay_list in enumerate(self.plot_overlay_list):

        #if len(overlay_list) > 1:

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

      self.plotting_utils = plotting_utils_theta_g.PlottingUtils(
        self.output_dir_performance, self.config_file)
      self.plot_single_performance(self.output_dir_performance)

    # Pb-Pb case
    else:

      # Plot for each R_max
      for R_max in self.max_distance:

        output_dir_performance = os.path.join(
          self.output_dir_performance, 'Rmax{}'.format(R_max))
        self.plotting_utils = plotting_utils_ang.PlottingUtils(
          output_dir_performance, self.config_file, R_max = R_max)
        self.plot_single_performance(output_dir_performance, R_max)

        # Plot for thermal model
        if self.do_thermal_closure and R_max == self.R_max:

          output_dir_performance = os.path.join(self.output_dir_performance, 'thermal')
          self.plotting_utils = plotting_utils_ang.PlottingUtils(
            output_dir_performance, self.config_file, R_max = R_max, thermal = True)
          self.plot_single_performance(output_dir_performance, R_max)

    return


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
    if not self.is_pp:
      self.create_output_subdir(output_dir_performance, 'delta_pt')

    # Generate performance plots
    for jetR in self.jetR_list:

      # Plot some subobservable-independent performance plots
      self.plotting_utils.plot_DeltaR(jetR, self.jet_matching_distance)
      self.plotting_utils.plot_JES(jetR)
      self.plotting_utils.plot_JES_proj(jetR, self.pt_bins_reported)
      self.plotting_utils.plotJER(
        jetR, self.utils.obs_label(self.obs_settings[0], self.grooming_settings[0]))
      self.plotting_utils.plot_jet_reco_efficiency(
        jetR, self.utils.obs_label(self.obs_settings[0], self.grooming_settings[0]))

      if not self.is_pp:
        self.plotting_utils.plot_delta_pt(jetR, self.pt_bins_reported)

      # Plot subobservable-dependent performance plots
      for i, _ in enumerate(self.obs_subconfig_list):

        obs_setting = self.obs_settings[i]
        grooming_setting = self.grooming_settings[i]
        obs_label = self.utils.obs_label(obs_setting, grooming_setting)

        self.plotting_utils.plot_obs_resolution(
          jetR, obs_label, self.xtitle, self.pt_bins_reported)
        self.plotting_utils.plot_obs_residual_pt(
          jetR, obs_label, self.xtitle, self.pt_bins_reported)
        self.plotting_utils.plot_obs_residual_obs(jetR, obs_label, self.xtitle)
        self.plotting_utils.plot_obs_projections(
          jetR, obs_label, obs_setting, grooming_setting, self.xtitle, self.pt_bins_reported)
        self.plotting_utils.plot_obs_truth(
          jetR, obs_label, obs_setting, grooming_setting, self.xtitle, self.pt_bins_reported)

    return


  #----------------------------------------------------------------------
  def plot_final_result(self, jetR, obs_label, obs_setting, grooming_setting):
    print('Plot final results for {}: R = {}, {}'.format(self.observable, jetR, obs_label))

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    # Loop through pt slices, and plot final result for each 1D distribution
    for i in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[i]
      max_pt_truth = self.pt_bins_reported[i+1]
      maxbin = self.obs_max_bins(obs_label)[i]

      self.plot_observable(jetR, obs_label, obs_setting, grooming_setting,
                           min_pt_truth, max_pt_truth, maxbin, plot_MC=True)

  #----------------------------------------------------------------------
  def plot_observable(self, jetR, obs_label, obs_setting, grooming_setting,
                      min_pt_truth, max_pt_truth, maxbin, plot_MC=False,
                      plot_theory=False, plot_theory_Fnp=False):

    self.set_logy = False
    #if grooming_setting:
    #  self.set_logy = True

    name = 'cResult_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt_truth, max_pt_truth)
    c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()
    
    c.cd()
    myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    myPad.SetLeftMargin(0.15)
    myPad.SetTopMargin(0.03)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.13)
    if self.set_logy:
      myPad.SetLogy()
    myPad.Draw()
    myPad.cd()
    
    xtitle = getattr(self, 'xtitle')
    ytitle = getattr(self, 'ytitle')
    color = 1   # black for data
    
    # Get histograms
    name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label,
                                              min_pt_truth, max_pt_truth)
    if grooming_setting:
      fraction_tagged = getattr(self, 'tagging_fraction_R{}_{}_{}-{}'.format(
        jetR, obs_label, min_pt_truth, max_pt_truth))
      #fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
      # maxbin+1 in grooming case to account for extra tagging bin
    if grooming_setting and maxbin:
      h = self.truncate_hist(getattr(self, name), None, maxbin+1, name+'_trunc')
    else:
      h = self.truncate_hist(getattr(self, name), None, maxbin, name+'_trunc')
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
    myBlankHisto.GetXaxis().SetTitleOffset(1.02)
    myBlankHisto.GetXaxis().SetTitleSize(0.055)
    myBlankHisto.SetYTitle(ytitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.1)
    myBlankHisto.GetYaxis().SetTitleSize(0.055)
    ymin = 1e-4 if self.set_logy else 0
    myBlankHisto.SetMinimum(ymin)
    maxval = max(2.3*h.GetBinContent(int(0.4*h.GetNbinsX())), 1.7*h.GetMaximum())
    myBlankHisto.SetMaximum(maxval)
    myBlankHisto.Draw("E")

    plot_pythia = False; plot_herwig = False;
    plot_jewel_no_recoils = False; plot_jewel_recoils = False;
    if plot_MC:

      hPythia = None; fraction_tagged_pythia = None;
      maxbin_adj = maxbin + 1 if grooming_setting else maxbin

      hPythia, fraction_tagged_pythia = self.MC_prediction(
        jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin_adj, 'Pythia')
      hHerwig, fraction_tagged_herwig = self.MC_prediction(
        jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin_adj, 'Herwig')
      hJewel_no_recoils, fraction_tagged_jewel_no_recoils = self.MC_prediction(
        jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin_adj, 'JEWEL', recoils=False)
      hJewel_recoils, fraction_tagged_jewel_no_recoils = self.MC_prediction(
        jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin_adj, 'JEWEL', recoils=True)

      if hPythia:
        hPythia.SetFillStyle(0)
        hPythia.SetMarkerSize(1.5)
        hPythia.SetMarkerStyle(21)
        hPythia.SetMarkerColor(600-6)
        hPythia.SetLineColor(600-6)
        hPythia.SetLineWidth(1)
        hPythia.Draw('E2 same')
        plot_pythia = True
      else:
        print('No PYTHIA prediction for %s %s' % (self.observable, obs_label))

      if hHerwig:
        hHerwig.SetFillStyle(0)
        hHerwig.SetMarkerSize(1.5)
        hHerwig.SetMarkerStyle(33)
        hHerwig.SetMarkerColor(46)
        hHerwig.SetLineColor(46)
        hHerwig.SetLineWidth(1)
        hHerwig.Draw('E2 same')
        plot_herwig = True
      else:
        print('No Herwig prediction for %s %s' % (self.observable, obs_label))

      if hJewel_no_recoils:
        hJewel_no_recoils.SetFillStyle(0)
        hJewel_no_recoils.SetMarkerSize(1.5)
        hJewel_no_recoils.SetMarkerStyle(34)
        hJewel_no_recoils.SetMarkerColor(42)
        hJewel_no_recoils.SetLineColor(42)
        hJewel_no_recoils.SetLineWidth(1)
        hJewel_no_recoils.Draw('E2 same')
        plot_jewel_no_recoils = True
      else:
        print('No JEWEL (recoils off) prediction for %s %s' % (self.observable, obs_label))

      if hJewel_recoils:
        hJewel_recoils.SetFillStyle(0)
        hJewel_recoils.SetMarkerSize(1.5)
        hJewel_recoils.SetMarkerStyle(47)
        hJewel_recoils.SetMarkerColor(7)
        hJewel_recoils.SetLineColor(7)
        hJewel_recoils.SetLineWidth(1)
        hJewel_recoils.Draw('E2 same')
        plot_jewel_recoils = True
      else:
        print('No JEWEL (recoils on) prediction for %s %s' % (self.observable, obs_label))

    h_sys.Draw("E2 same")
    h.Draw("PE X0 same")

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text_xval = 0.52 if grooming_setting else 0.55
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(text_xval, 0.92, text)

    text = '0-10% Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(text_xval, 0.85, text)

    text = "anti-#it{k}_{T} jets,   #it{R} = %s" % str(jetR)
    text_latex.DrawLatex(text_xval, 0.78, text)

    if plot_theory and show_parton_theory:
      text = str(min_pt_truth) + ' < #it{p}_{T,jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    else:
      text = str(min_pt_truth) + ' < #it{p}_{T,jet}^{ch} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.DrawLatex(text_xval, 0.71, text)

    text = '| #it{#eta}_{jet}| < %s' % str(0.9 - jetR)
    subobs_label = self.utils.formatted_subobs_label(self.observable)
    if subobs_label:
      text += ',   %s = %s' % (subobs_label, obs_setting)
    delta = 0.07
    text_latex.DrawLatex(text_xval, 0.71-delta, text)

    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting) #.replace("#beta}", "#beta}_{SD}")
      text_latex.DrawLatex(text_xval, 0.71-2*delta, text)

      text_latex.SetTextSize(0.04)
      text = '#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged
      if plot_pythia:
        text += (', #it{f}_{tagged}^{pythia} = %3.3f' % fraction_tagged_pythia)
      if plot_herwig:
        text += (', #it{f}_{tagged}^{herwig} = %3.3f' % fraction_tagged_herwig)
      text_latex.DrawLatex(text_xval, 0.71-3*delta, text)

    miny = 0.62
    #if plot_pythia:
    #  if plot_herwig:
    #    miny = 0.72  #TODO
    #  else:
    #    miny = 0.72
    myLegend = ROOT.TLegend(0.16, miny, text_xval-0.02, 0.96)
    self.utils.setup_legend(myLegend, 0.035)
    myLegend.AddEntry(h, 'ALICE Pb-Pb data', 'pe')
    myLegend.AddEntry(h_sys, 'Syst. uncertainty', 'f')
    if plot_pythia:
      myLegend.AddEntry(hPythia, 'PYTHIA8 Monash2013', 'pe')
    if plot_herwig:
      myLegend.AddEntry(hHerwig, 'Herwig7 default tune', 'pe')
    if plot_jewel_no_recoils:
      myLegend.AddEntry(hJewel_no_recoils, 'JEWEL Pb-Pb (recoils off)', 'pe')
    if plot_jewel_recoils:
      myLegend.AddEntry(hJewel_recoils, 'JEWEL Pb-Pb (recoils on)', 'pe')
    myLegend.Draw()

    name = 'hUnfolded_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), obs_label,
                                             int(min_pt_truth), int(max_pt_truth), self.file_format)
    if plot_MC:
      name = 'hUnfolded_R{}_{}_{}-{}_MC{}'.format(
        self.utils.remove_periods(jetR), obs_label, int(min_pt_truth),
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
    if plot_pythia:
      hPythia.Write()
    if plot_herwig:
      hHerwig.Write()
    if plot_jewel_no_recoils:
      hJewel_no_recoils.Write()
    if plot_jewel_recoils:
      hJewel_recoils.Write()
    fFinalResults.Close()


  #----------------------------------------------------------------------
  # Get unfolded data from the previous measurement
  def get_h_prev_result(self, jetR, obs_label, min_pt_truth, max_pt_truth):

    output_dir = getattr(self, 'output_dir_main')

    f = ROOT.TFile(self.results_pp, 'READ')

    # Retrieve previous result and ensure that it has the proper bin range
    h_prev_data_name = 'hmain_%s_R%s_%s_%s-%s_trunc' % \
      (self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    h_prev_data = f.Get(h_prev_data_name)
    if not h_prev_data:
      raise AttributeError("%s not found in file %s" % (h_prev_data_name, self.results_pp))

    # Rename and steal ownership from ROOT
    name = 'h_prev_result_{}_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    h_prev_data.SetNameTitle(name, name)
    h_prev_data.SetDirectory(0)

    # Retrieve previous result systematics and ensure that it has the proper bin range
    h_prev_sys_up_name = 'hResult_%s_sys_plus_R%s_%s_%s-%s' % \
      (self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    h_prev_sys_up = f.Get(h_prev_sys_up_name)
    if not h_prev_sys_up:
      raise AttributeError("%s not found in file %s" % (h_prev_sys_up_name, self.results_pp))
    name = 'h_prev_sys_up_{}_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    h_prev_sys_up.SetNameTitle(name, name)
    h_prev_sys_up.SetDirectory(0)

    h_prev_sys_down_name = 'hResult_%s_sys_minus_R%s_%s_%s-%s' % \
      (self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    h_prev_sys_down = f.Get(h_prev_sys_down_name)
    if not h_prev_sys_down:
      raise AttributeError("%s not found in file %s" % (h_prev_sys_down_name, self.results_pp))
    name = 'h_prev_sys_down_{}_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    h_prev_sys_down.SetNameTitle(name, name)
    h_prev_sys_down.SetDirectory(0)

    # Set normalization to 1 in all histograms
    h_prev_data.Scale(jetR)
    h_prev_sys_up.Scale(jetR)
    h_prev_sys_down.Scale(jetR)

    return h_prev_data, h_prev_sys_up, h_prev_sys_down


  #----------------------------------------------------------------------
  def MC_prediction(self, jetR, obs_setting, obs_label, min_pt_truth,
                    max_pt_truth, maxbin, MC='Pythia', overlay=False, recoils=False):

    if MC.lower() == 'pythia':
      hMC = self.get_pythia_from_response(jetR, obs_label, min_pt_truth,
                                          max_pt_truth, maxbin, overlay)
    elif MC.lower() == 'herwig':
      hMC = self.get_herwig_from_response(jetR, obs_label, min_pt_truth,
                                          max_pt_truth, maxbin, overlay)
    elif MC.lower() == "jewel":
      hMC = self.get_jewel_from_response(
        jetR, obs_label, min_pt_truth, max_pt_truth, maxbin, overlay, recoils)
    else:
      raise NotImplementedError("MC must be PYTHIA, Herwig, or JEWEL.")

    n_jets_inclusive = hMC.Integral(0, hMC.GetNbinsX()+1)
    n_jets_tagged = hMC.Integral(hMC.FindBin(
      self.truth_bin_array(obs_label)[0]), hMC.GetNbinsX())

    fraction_tagged_MC =  n_jets_tagged/n_jets_inclusive
    hMC.Scale(1./n_jets_inclusive, 'width')

    return [hMC, fraction_tagged_MC]

  #----------------------------------------------------------------------
  def get_pythia_from_response(self, jetR, obs_label, min_pt_truth, max_pt_truth,
                               maxbin, overlay=False):

    output_dir = getattr(self, 'output_dir_main')

    filepath = os.path.join(output_dir, 'response.root')
    f = ROOT.TFile(filepath, 'READ')

    thn_name = 'hResponse_JetPt_{}_R{}_{}_rebinned'.format(self.observable, jetR, obs_label)
    thn = f.Get(thn_name)
    thn.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)

    name = 'hPythia_{}_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    h = self.truncate_hist(thn.Projection(3), None, maxbin, name)
    h.SetDirectory(0)

    return h

  #----------------------------------------------------------------------
  def get_herwig_from_response(self, jetR, obs_label, min_pt_truth, max_pt_truth,
                               maxbin, overlay=False):

    filepath = os.path.join(self.output_dir_fastsim_generator1, 'response.root')
    f = ROOT.TFile(filepath, 'READ')

    thn_name = 'hResponse_JetPt_{}_R{}_{}_rebinned'.format(self.observable, jetR, obs_label)
    thn = f.Get(thn_name)
    thn.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)

    name = 'hHerwig_{}_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    h = self.truncate_hist(thn.Projection(3), None, maxbin, name)
    h.SetDirectory(0)

    return h

  #----------------------------------------------------------------------
  def get_jewel_from_response(self, jetR, obs_label, min_pt_truth, max_pt_truth,
                               maxbin, overlay=False, recoils=False):

    gen_to_use = 2 + int(recoils)
    filepath = os.path.join(
      getattr(self, "output_dir_fastsim_generator%i" % gen_to_use), 'response.root')
    f = ROOT.TFile(filepath, 'READ')

    thn_name = 'hResponse_JetPt_{}_R{}_{}_rebinned'.format(self.observable, jetR, obs_label)
    thn = f.Get(thn_name)
    thn.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)

    name = 'hJewel_recoils_{}_{}_R{}_{}_{}-{}'.format(
      "on" if recoils else "off", self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    h = self.truncate_hist(thn.Projection(3), None, maxbin, name)
    h.SetDirectory(0)

    return h

  #----------------------------------------------------------------------
  def get_pp_data(self, jetR, obs_label, min_pt_truth, max_pt_truth,
                  xbins, overlay=False):

    output_dir = getattr(self, 'output_dir_main')

    f = ROOT.TFile(self.results_pp, 'READ')

    # Retrieve pp data and ensure that it has the proper bin range
    h_pp_data_name = 'hmain_%s_R%s_%s_%s-%s_trunc' % \
      (self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    h_pp_data = f.Get(h_pp_data_name)
    if not h_pp_data:
      raise AttributeError("%s not found in file %s" % (h_pp_data_name, self.results_pp))

    name = 'h_pp_data_{}_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    if (xbins == [h_pp_data.GetBinLowEdge(i) for i in range(1, h_pp_data.GetNbinsX()+2)]):
      h_pp_data.SetNameTitle(name, name)
    else:
      h_pp_data = h_pp_data.Rebin(len(xbins)-1, name, array('d', xbins))
    h_pp_data.SetDirectory(0)

    # Retrieve pp systematics and ensure that it has the proper bin range
    h_pp_sys_name = 'hResult_%s_systotal_R%s_%s_n3_%s-%s' % \
      (self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    h_pp_sys = f.Get(h_pp_sys_name)
    if not h_pp_sys:
      raise AttributeError("%s not found in file %s" % (h_pp_sys_name, self.results_pp))

    name = 'h_pp_sys_{}_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    if (xbins == [h_pp_sys.GetBinLowEdge(i) for i in range(1, h_pp_sys.GetNbinsX()+2)]):
      h_pp_sys.SetNameTitle(name, name)
    else:
      h_pp_sys = h_pp_sys.Rebin(len(xbins)-1, name, array('d', xbins))
    h_pp_sys.SetDirectory(0)

    return h_pp_data, h_pp_sys


  #----------------------------------------------------------------------
  def plot_final_result_overlay(self, i_config, jetR, overlay_list):
    print('Plotting overlay of', overlay_list)

    # Plot overlay of different subconfigs, for fixed pt bin
    for i in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[i]
      max_pt_truth = self.pt_bins_reported[i+1]
      maxbins = [self.obs_max_bins(obs_label)[i] for obs_label in self.obs_labels]

      # Plot PYTHIA comparison plots
      self.plot_observable_overlay_subconfigs(
        i_config, jetR, overlay_list, min_pt_truth,
        max_pt_truth, maxbins, plot_MC=True, MC='PYTHIA', plot_ratio=True)

      if not self.is_pp and self.results_pp:
        # Plot Pb-Pb/pp data comparison plots
        self.plot_observable_overlay_subconfigs(
          i_config, jetR, overlay_list, min_pt_truth,
          max_pt_truth, maxbins, plot_pp_data=True, plot_ratio=True)


  #----------------------------------------------------------------------
  def plot_observable_overlay_subconfigs(self, i_config, jetR, overlay_list, min_pt_truth,
                                         max_pt_truth, maxbins, plot_pp_data=False, plot_MC=False,
                                         MC='PYTHIA', plot_nll=False, plot_ratio=False):

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
    pad1.SetTicks(0,1)
    if plot_ratio:
      pad1.SetBottomMargin(0.)
    # set log y axis if all configs are SD
    setlogy = False
    for name in overlay_list:
      if "SD" not in name:
        break
      elif name == overlay_list[-1]:
        setlogy = True
        pad1.SetLogy()
    pad1.Draw()
    pad1.cd()

    myLegend = ROOT.TLegend(0.62, 0.6, 0.96, 0.91)
    self.utils.setup_legend(myLegend, 0.045)
    myLegend2 = ROOT.TLegend(0.81, 0.788, 1, 0.91)
    self.utils.setup_legend(myLegend2, 0.045)

    name = 'hmain_{}_R{}_{{}}_{}-{}'.format(self.observable, jetR, min_pt_truth, max_pt_truth)
    ymax, ymin = self.get_max_min(name, overlay_list, maxbins)

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
        marker_pythia = 24
        color = 1
      elif subconfig_name == overlay_list[1]:
        marker = 21
        marker_pythia = 25
        color = 600-6
      elif subconfig_name == overlay_list[2]:
        marker = 33
        marker_pythia = 27
        color = 632-4
      else:  # subconfig_name == overlay_list[3]:
        marker = 34
        marker_pythia = 28
        color = 416-2

      name = 'hmain_{}_R{}_{}_{}-{}'.format(
        self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
      if grooming_setting:
        fraction_tagged = getattr(self, 'tagging_fraction_R{}_{}_{}-{}'.format(
          jetR, obs_label, min_pt_truth, max_pt_truth))

      if grooming_setting and maxbin:
        h = self.truncate_hist(getattr(self, name), None, maxbin+1, name+'_trunc')
      else:
        h = self.truncate_hist(getattr(self, name), None, maxbin, name+'_trunc')
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
            if self.use_prev_result:
              myBlankHisto.SetMaximum(1.8*ymax)
            else:
              myBlankHisto.SetMaximum(1.2*ymax)
          elif min_pt_truth == 60:
            myBlankHisto.SetMaximum(1.15*ymax)
          else:
            myBlankHisto.SetMaximum(1.3*ymax)
        elif jetR == 0.4:
          if min_pt_truth == 20:
            myBlankHisto.SetMaximum(1.5*ymax)
          elif min_pt_truth == 40:
            myBlankHisto.SetMaximum(1.35*ymax)
          elif min_pt_truth == 60:
            myBlankHisto.SetMaximum(1.15*ymax)
          else:
            myBlankHisto.SetMaximum(1.5*ymax)
        else:
          myBlankHisto.SetMaximum(1.5*ymax)
        if setlogy:
          myBlankHisto.SetMaximum(5*ymax)
        myBlankHisto.SetMinimum(0.)
        if plot_ratio:
          myBlankHisto.SetMinimum(2e-4) # Don't draw 0 on top panel
          if setlogy:  # the minimum value matters
            myBlankHisto.SetMinimum(2*ymin/3)
          myBlankHisto.GetYaxis().SetTitleSize(0.065)
          myBlankHisto.GetYaxis().SetTitleOffset(1.1)
          myBlankHisto.GetYaxis().SetLabelSize(0.06)
        myBlankHisto.Draw('E')

        # Initialize ratio plot
        if plot_ratio:

          c.cd()
          pad2 = ROOT.TPad("pad2", "pad2", 0, 0.02, 1, 0.3)
          pad2.SetTopMargin(0)
          pad2.SetBottomMargin(0.5)
          pad2.SetLeftMargin(0.2)
          pad2.SetRightMargin(0.04)
          pad2.SetTicks(1,1)
          pad2.Draw()
          pad2.cd()

          myBlankHisto2 = myBlankHisto.Clone("myBlankHisto_C")
          if plot_MC:
            myBlankHisto2.SetYTitle("#frac{Data}{%s}" % MC)
          elif plot_pp_data:
            if self.use_prev_result:
              myBlankHisto2.SetYTitle("#frac{5.02 TeV}{2.76 TeV}")
            else:
              myBlankHisto2.SetYTitle("#frac{Pb-Pb}{pp}")
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
              myBlankHisto2.GetYaxis().SetRangeUser(0.4, 2.2)
            elif jetR == 0.4:
              myBlankHisto2.GetYaxis().SetRangeUser(0.5, 1.9)
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

      hMC = None; fraction_tagged_MC = None;
      h_pp_data = None; h_pp_sys = None;
      if plot_MC:
        if MC.lower() == "pythia":
          if grooming_setting and maxbin:
            hMC, fraction_tagged_MC = self.MC_prediction(
              jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin+1,
              MC='Pythia', overlay=True)
          else:
            hMC, fraction_tagged_MC = self.MC_prediction(
              jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin,
              MC='Pythia', overlay=True)

        elif MC.lower() == "herwig":
          if grooming_setting:
            hMC, fraction_tagged_MC = self.MC_prediction(
              jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin+1,
              'Herwig', overlay=True)
          else:
            hMC, fraction_tagged_MC = self.MC_prediction(
              jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin,
              'Herwig', overlay=True)

        else:
          raise NotImplementedError("MC must be either Pythia or Herwig.")

        plot_errors = False
        if plot_errors:
          hMC.SetMarkerSize(0)
          hMC.SetMarkerStyle(0)
          hMC.SetMarkerColor(color)
          hMC.SetFillColor(color)
        else:
          hMC.SetLineColor(color)
          hMC.SetLineColorAlpha(color, 0.5)
          hMC.SetLineWidth(4)

      elif plot_pp_data:

        if self.use_prev_result:
          # Use previous Pb-Pb result (use pp name for convenience)
          h_pp_data, h_prev_sys_up, h_prev_sys_down = self.get_h_prev_result(
            jetR, obs_label, min_pt_truth, max_pt_truth)
          # Construct TGraphAsymmErrors to store asymmetric up and down uncertainties
          x = array('d', [h_prev_sys_up.GetBinCenter(i) for i in \
                          range(1, h_prev_sys_up.GetNbinsX()+1)])
          y = array('d', [h_prev_sys_up.GetBinContent(i) for i in \
                          range(1, h_prev_sys_up.GetNbinsX()+1)])
          exl = array('d', [x[i-1]-h_prev_sys_up.GetBinLowEdge(i) for i in \
                            range(1, h_prev_sys_up.GetNbinsX()+1)])
          eyl = array('d', [h_prev_sys_down.GetBinError(i) for i in \
                            range(1, h_prev_sys_down.GetNbinsX()+1)])
          exh = array('d', [h_prev_sys_up.GetBinLowEdge(i)-x[i-2] for i in \
                            range(2, h_prev_sys_up.GetNbinsX()+2)])
          eyh = array('d', [h_prev_sys_up.GetBinError(i) for i in \
                            range(1, h_prev_sys_up.GetNbinsX()+1)])
          h_pp_sys = ROOT.TGraphAsymmErrors(len(x), x, y, exl, exh, eyl, eyh)
        else:
          h_pp_data, h_pp_sys = self.get_pp_data(
            jetR, obs_label, min_pt_truth, max_pt_truth,
            [h.GetBinLowEdge(i) for i in range(1, h.GetNbinsX()+2)])

        plot_errors = True
        if plot_errors:
          h_pp_data.SetMarkerSize(1.5)
          h_pp_data.SetMarkerStyle(marker_pythia) #27)
          h_pp_data.SetMarkerColor(color)
          h_pp_data.SetFillColor(color)
          h_pp_data.SetLineStyle(9)
          h_pp_data.SetLineWidth(2)
          h_pp_data.SetLineColor(color)
          h_pp_sys.SetLineColor(0)
          h_pp_sys.SetFillColor(color)
          #h_pp_sys.SetFillColorAlpha(color, 0.8)
          h_pp_sys.SetFillStyle(3004)
          h_pp_sys.SetLineWidth(0)
        else:
          h_pp_data.SetLineColor(color)
          h_pp_data.SetLineColorAlpha(color, 0.5)
          h_pp_data.SetLineWidth(4)

      if plot_ratio:
        if self.use_prev_result and plot_pp_data:
          # Take ratios to get the correct systematic uncertainties
          hRatioSysUp = h_sys.Clone()
          hRatioSysUp.SetName('{}_RatioUp'.format(h_sys.GetName()))
          hRatioSysUp.Divide(h_prev_sys_up)
          hRatioSysDown = h_sys.Clone()
          hRatioSysDown.SetName('{}_RatioDown'.format(h_sys.GetName()))
          hRatioSysDown.Divide(h_prev_sys_down)
          # Construct TGraphAsymmErrors
          x = array('d', [hRatioSysUp.GetBinCenter(i) for i in \
                          range(1, hRatioSysUp.GetNbinsX()+1)])
          y = array('d', [hRatioSysUp.GetBinContent(i) for i in \
                          range(1, hRatioSysUp.GetNbinsX()+1)])
          exl = array('d', [x[i-1]-hRatioSysUp.GetBinLowEdge(i) for i in \
                            range(1, hRatioSysUp.GetNbinsX()+1)])
          eyl = array('d', [hRatioSysDown.GetBinError(i) for i in \
                            range(1, hRatioSysDown.GetNbinsX()+1)])
          exh = array('d', [hRatioSysUp.GetBinLowEdge(i)-x[i-2] for i in \
                            range(2, hRatioSysUp.GetNbinsX()+2)])
          eyh = array('d', [hRatioSysUp.GetBinError(i) for i in \
                            range(1, hRatioSysUp.GetNbinsX()+1)])
          hRatioSys = ROOT.TGraphAsymmErrors(len(x), x, y, exl, exh, eyl, eyh)
        else:
          hRatioSys = h_sys.Clone()
          hRatioSys.SetName('{}_Ratio'.format(h_sys.GetName()))
          if plot_MC:
            hRatioSys.Divide(hMC)
          elif plot_pp_data:
            hRatioSys.Divide(h_pp_sys)
        hRatioSys.SetLineColor(0)
        hRatioSys.SetFillColor(color)
        hRatioSys.SetFillColorAlpha(color, 0.3)
        hRatioSys.SetFillStyle(1001)
        hRatioSys.SetLineWidth(0)
        hRatioSys.SetMaximum(1.99)

        hRatioStat = h.Clone()
        hRatioStat.SetName('{}_Ratio'.format(h.GetName()))
        if plot_MC:
          hRatioStat.Divide(hMC)
        elif plot_pp_data:
          hRatioStat.Divide(h_pp_data)
          #for i in range(1, hRatioStat.GetNbinsX()+1):
          #  new_error = math.sqrt(h.GetBinError(i) ** 2 + h_pp_data.GetBinError(i) ** 2)
          #  hRatioStat.SetBinError(i, new_error)
        hRatioStat.SetMarkerSize(1.5)
        hRatioStat.SetMarkerStyle(marker)
        hRatioStat.SetMarkerColor(color)
        hRatioStat.SetLineStyle(1)
        hRatioStat.SetLineWidth(2)
        hRatioStat.SetLineColor(color)
        hRatioStat.SetMaximum(1.99)

      pad1.cd()
      if plot_MC:
        plot_errors = False
        if plot_errors:
          hMC.DrawCopy('E3 same')
        else:
          hMC.DrawCopy('L hist same')

      elif plot_pp_data:
        plot_errors = True
        if plot_errors:
          # TGraphAsymmErrors doesn't have DrawCopy
          if self.use_prev_result:
            h_pp_sys.Draw('E2 same')
            h_pp_data.Draw('PE X0 same')
          else:
            h_pp_sys.DrawCopy('E2 same')
            h_pp_data.DrawCopy('PE X0 same')
        else:
          h_pp_data.DrawCopy('L hist same')

      h_sys.DrawCopy('E2 same')
      h.DrawCopy('PE X0 same')

      if plot_ratio:
        pad2.cd()
        if plot_MC or plot_pp_data:
          # TGraphAsymmErrors doesn't have DrawCopy
          if self.use_prev_result:
            hRatioSys.Draw('E2 same')
            hRatioStat.Draw('PE X0 same')
          else:
            hRatioSys.DrawCopy('E2 same')
            hRatioStat.DrawCopy('PE X0 same')

      subobs_label = self.utils.formatted_subobs_label(self.observable)
      text = ''
      if subobs_label:
        text += '%s = %s' % (subobs_label, obs_setting)
      text_list.append(text)
      h_list.append(h)

    pad1.cd()
    for i, h, text in zip(range(len(h_list)), h_list, text_list):
      if i < 2:
        myLegend.AddEntry(h, text, 'pe')
      else:
        myLegend2.AddEntry(h, text, 'pe')
    myLegend.AddEntry(h_sys, 'Syst. uncertainty', 'f')
    if plot_MC:
      if MC.lower() == "pythia":
        myLegend.AddEntry(hMC, 'PYTHIA8 Monash2013', 'l')
      elif MC.lower() == "herwig":
        myLegend.AddEntry(hMC, 'Herwig7 Default', 'l')
    elif plot_pp_data:
      if self.use_prev_result:
        myLegend.AddEntry(h_pp_data, 'Pb-Pb @ 2.76 TeV', 'pe')
        if plot_errors:
          myLegend.AddEntry(h_pp_sys, '2.76 TeV syst. uncert.', 'f')
      else:
        myLegend.AddEntry(h_pp_data, 'pp', 'pe')
        if plot_errors:
          myLegend.AddEntry(h_pp_sys, 'pp syst. uncert.', 'f')

    text_xval = 0.27
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(text_xval, 0.87, text)

    text = 'Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(text_xval, 0.81, text)

    text = 'Charged jets   anti-#it{k}_{T}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(text_xval, 0.75, text)

    text = '#it{R} = ' + str(jetR) + ',   | #it{#eta}_{jet}| < %s' % str(0.9 - jetR)
    text_latex.DrawLatex(text_xval, 0.69, text)

    text = str(min_pt_truth) + ' < #it{p}_{T,jet}^{ch} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(text_xval, 0.63, text)

    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting)#.replace("#beta}", "#beta}_{SD}")
      text_latex.DrawLatex(text_xval, 0.57, text)

    myLegend.Draw()
    if len(h_list) > 2:
      myLegend2.Draw()

    name = 'h_{}_R{}_{}-{}_{}{}'.format(self.observable,
                                        self.utils.remove_periods(jetR), int(min_pt_truth),
                                        int(max_pt_truth), i_config, self.file_format)
    if plot_MC:
      name = 'h_{}_R{}_{}-{}_{}_{}{}'.format(self.observable, self.utils.remove_periods(jetR),
                                                 int(min_pt_truth), int(max_pt_truth), MC,
                                                 i_config, self.file_format)
    elif plot_pp_data:
      name = 'h_{}_R{}_{}-{}_ppComp_{}{}'.format(self.observable, self.utils.remove_periods(jetR),
                                                 int(min_pt_truth), int(max_pt_truth),
                                                 i_config, self.file_format)

    output_dir = getattr(self, 'output_dir_final_results')
    if not os.path.exists(os.path.join(output_dir, 'all_results')):
      os.mkdir(os.path.join(output_dir, 'all_results'))
    outputFilename = os.path.join(output_dir, 'all_results', name)
    c.SaveAs(outputFilename)

    # Write result to ROOT file
    final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
    fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
    c.Write()

    c.Close()

  #----------------------------------------------------------------------
  # Return maximum & minimum y-values of unfolded results in a subconfig list
  def get_max_min(self, name, overlay_list, maxbins):

    total_min = 1e10
    total_max = -1e10

    for i, subconfig_name in enumerate(self.obs_subconfig_list):

      if subconfig_name not in overlay_list:
        continue

      obs_setting = self.obs_settings[i]
      grooming_setting = self.grooming_settings[i]
      obs_label = self.utils.obs_label(obs_setting, grooming_setting)
      maxbin = maxbins[i]

      h = getattr(self, name.format(obs_label))
      if 'SD' in obs_label:
        content = [ h.GetBinContent(j) for j in range(2, maxbin+2) ]
      else:
        content = [ h.GetBinContent(j) for j in range(1, maxbin+1) ]

      min_val = min(content)
      if min_val < total_min:
        total_min = min_val
      max_val = max(content)
      if max_val > total_max:
        total_max = max_val


    return (total_max, total_min)


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
