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

    self.theory_predictions = config["theory_predictions"] if \
      "theory_predictions" in config else []
    self.theory_predictions_names = config["theory_predictions_names"] if \
      "theory_predictions_names" in config else []

    # Whether or not to use the previous measurement in ratio
    self.use_prev_result = config["use_prev_result"]

    self.histutils = ROOT.RUtil.HistUtils()

    self.colors = [ROOT.kRed+1, ROOT.kGreen+2, ROOT.kBlue, ROOT.kOrange+2,
        ROOT.kViolet+2, ROOT.kCyan+1, ROOT.kPink+10, ROOT.kGray+1,
        ROOT.kYellow+4, ROOT.kAzure+8, ROOT.kRed]
    #self.colors = self.ColorArray

  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #---------------------------------------------------------------
  def plot_single_result(self, jetR, obs_label, obs_setting, grooming_setting):
    #print('Plotting each individual result...')

    # Plot final result for each 1D substructure distribution (with PYTHIA)
    self.plot_final_result(jetR, obs_label, obs_setting, grooming_setting, draw_ratio=True)


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
  def plot_final_result(self, jetR, obs_label, obs_setting, grooming_setting, draw_ratio=False):
    print('Plot final results for {}: R = {}, {}'.format(self.observable, jetR, obs_label))

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    # Loop through pt slices, and plot final result for each 1D distribution
    for i in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[i]
      max_pt_truth = self.pt_bins_reported[i+1]
      maxbin = self.obs_max_bins(obs_label)[i]

      # Do special plotting for comparison to Run 1 case
      if self.use_prev_result:
        self.plot_observable(
          jetR, obs_label, obs_setting, grooming_setting, min_pt_truth,
          max_pt_truth, maxbin, plot_MC=False, draw_ratio=True)
        return

      self.plot_observable(
        jetR, obs_label, obs_setting, grooming_setting, min_pt_truth,
        max_pt_truth, maxbin, plot_MC=True, draw_ratio=True)

      if self.results_pp:
        self.plot_observable(
          jetR, obs_label, obs_setting, grooming_setting, min_pt_truth,
          max_pt_truth, maxbin, plot_pp_data=True, plot_MC=True, draw_ratio=True)

        self.plot_observable(
          jetR, obs_label, obs_setting, grooming_setting, min_pt_truth, max_pt_truth,
          maxbin, plot_pp_data=True, plot_MC=True, plot_PbPb=False, draw_ratio=True)

  #----------------------------------------------------------------------
  def plot_observable(
    self, jetR, obs_label, obs_setting, grooming_setting, min_pt_truth, max_pt_truth,
    maxbin, plot_pp_data=False, plot_MC=False, plot_PbPb=True, draw_ratio=False):

    self.set_logy = True if (grooming_setting and self.observable == "ang") else False
    make_ratio_plot = True if (plot_MC or plot_pp_data and draw_ratio) else False
    match_data_normalization = True
    plot_pythia_and_herwig = (not plot_PbPb)
    ignore_MC_top_panel = (plot_pp_data and draw_ratio and plot_PbPb)

    name = 'cResult_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt_truth, max_pt_truth)
    c = ROOT.TCanvas(name, name, 600, 900 if make_ratio_plot else 450)
    c.Draw()

    c.cd()
    pad_y_split = 0.5
    if self.observable == "mass" and min_pt_truth == 80:
      pad_left_margin = 0.16
    else:
      pad_left_margin = 0.16 if self.set_logy else 0.15
    myPad = ROOT.TPad('myPad', 'The pad', 0, pad_y_split if make_ratio_plot else 0, 1, 1)
    myPad.SetLeftMargin(pad_left_margin)
    myPad.SetTopMargin(0.03)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0 if make_ratio_plot else 0.13)
    if self.set_logy:
      myPad.SetLogy()
    myPad.SetTicks(1, 1)
    myPad.Draw()
    myPad.cd()

    color = 1   # black for data

    # Get histograms
    name = 'hmain_{}_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    if grooming_setting:
      fraction_tagged = getattr(self, 'tagging_fraction_R{}_{}_{}-{}'.format(
        jetR, obs_label, min_pt_truth, max_pt_truth))
      #fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
      # maxbin+1 in grooming case to account for extra tagging bin
    if grooming_setting and maxbin:
      h = self.truncate_hist(getattr(self, name), None, maxbin+1, (name+'_trunc').replace("__","_"))
    else:
      h = self.truncate_hist(getattr(self, name), None, maxbin, (name+'_trunc').replace("__", "_"))
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
    alpha = obs_label.split("_")[0]
    if not make_ratio_plot:
      myBlankHisto.SetXTitle(self.xtitle.replace("alpha}", "alpha}="+alpha))
      myBlankHisto.GetXaxis().SetTitleOffset(1.02)
      myBlankHisto.GetXaxis().SetTitleSize(0.055)
    myBlankHisto.SetYTitle(self.ytitle.replace("alpha}", "alpha}="+alpha))
    if self.observable == "mass" and min_pt_truth == 80:
      myBlankHisto.GetYaxis().SetTitleOffset(1.35)
    else:
      myBlankHisto.GetYaxis().SetTitleOffset(1.3 if self.set_logy else 1.15)
    myBlankHisto.GetYaxis().SetTitleSize(0.055)

    plot_pythia = False; plot_herwig = False;
    plot_jewel_no_recoils = False; plot_jewel_recoils = False; plot_jewel_pp = False
    plot_jetscape = False; plot_jetscape_pp = False;
    plot_hybrid = False; plot_hybrid_pp = False;
    plot_zhang = False; plot_zhang_pp = False;
    n_pp_models = 0; n_AA_models = 0
    if plot_MC:

      hPythia = None; fraction_tagged_pythia = None;
      hHerwig = None; fraction_tagged_herwig = None;
      hJewel_pp = None; fraction_tagged_jewel_pp = None;
      hJetscape_pp = None; fraction_tagged_jetscape_pp = None;
      hJetscape = None; fraction_tagged_jetscape = None;
      hZhang_pp = None; fraction_tagged_zhang_pp = None;
      hZhang = None; fraction_tagged_zhang = None;
      maxbin_adj = maxbin + 1 if (maxbin != None and grooming_setting) else maxbin

      if plot_pp_data:
        if plot_pythia_and_herwig:
          hPythia, fraction_tagged_pythia = self.MC_prediction(
            jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin_adj, 'Pythia')
          hHerwig, fraction_tagged_herwig = self.MC_prediction(
            jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin_adj, 'Herwig')
        hJewel_pp, fraction_tagged_jewel_pp = self.MC_prediction(
          jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin_adj, 'JEWEL pp')
        hJetscape_pp, fraction_tagged_jetscape_pp = self.MC_prediction(
          jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin_adj, 'JETSCAPE pp')
        if not grooming_setting:  # Disable groomed predictions for now
          hZhang_pp, fraction_tagged_zhang_pp = self.MC_prediction(
            jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin_adj, 'Zhang pp')

      hJewel_no_recoils, fraction_tagged_jewel_no_recoils = self.MC_prediction(
        jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin_adj, 'JEWEL', recoils=False)
      hJewel_recoils, fraction_tagged_jewel_recoils = self.MC_prediction(
        jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin_adj, 'JEWEL', recoils=True)
      hJetscape, fraction_tagged_jetscape = self.MC_prediction(
        jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin_adj, 'JETSCAPE')
      if not grooming_setting:  # Disable groomed predictions for now
        hZhang, fraction_tagged_zhang = self.MC_prediction(
          jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin_adj, 'Zhang')

      hHybridNoElastic_pp, hHybridWithElastic_pp, hHybridNoElastic, hHybridWithElastic = self.get_hybrid(
        jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin_adj)

      if hPythia:
        color_pythia = self.colors[0]
        # Create clone with 0 error for solid line
        hPythia_draw = hPythia.Clone(hPythia.GetName()+"_drawclone")
        for i in range(hPythia_draw.GetNbinsX()+1):
          hPythia_draw.SetBinError(i, 0)
        #hPythia_draw.SetMarkerSize(1.5)
        #hPythia_draw.SetMarkerStyle(21)
        #hPythia_draw.SetMarkerColor(600-6)
        #hPythia_draw.SetFillColorAlpha(color_pythia, 0.5)
        #hPythia_draw.SetFillStyle(1001)
        hPythia_draw.SetLineStyle(1)
        hPythia_draw.SetLineColor(color_pythia)
        hPythia_draw.SetLineWidth(4)
        plot_pythia = True
        n_pp_models += 1
      #else:
      #  print('No PYTHIA prediction for %s %s' % (self.observable, obs_label))

      if hHerwig:
        color_herwig = self.colors[1]
        # Create clone with 0 error for solid line
        hHerwig_draw = hHerwig.Clone(hHerwig.GetName()+"_drawclone")
        for i in range(hHerwig_draw.GetNbinsX()+1):
          hHerwig_draw.SetBinError(i, 0)
        #hHerwig_draw.SetMarkerSize(1.5)
        #hHerwig_draw.SetMarkerStyle(33)
        #hHerwig_draw.SetMarkerColor(2)
        #hHerwig_draw.SetFillColorAlpha(color_herwig, 0.5)
        #hHerwig_draw.SetFillStyle(1001)
        hHerwig_draw.SetLineStyle(7)
        hHerwig_draw.SetLineColor(color_herwig)
        hHerwig_draw.SetLineWidth(4)
        plot_herwig = True
        n_pp_models += 1
      #else:
      #  print('No Herwig prediction for %s %s' % (self.observable, obs_label))

      if hJewel_pp:
        # Create clone with 0 error for solid line
        hJewel_pp_draw = hJewel_pp.Clone(hJewel_pp.GetName()+"_drawclone")
        for i in range(hJewel_pp_draw.GetNbinsX()+1):
          hJewel_pp_draw.SetBinError(i, 0)
        color_jewel_pp = self.colors[2]
        #hJewel_pp_draw.SetMarkerSize(1.5)
        #hJewel_pp_draw.SetMarkerStyle(28)
        #hJewel_pp_draw.SetMarkerColor(8)
        #hJewel_pp_draw.SetFillColorAlpha(color_jewel_pp, 0.5)
        #hJewel_pp_draw.SetFillStyle(1001)
        hJewel_pp_draw.SetLineStyle(8)
        hJewel_pp_draw.SetLineColor(color_jewel_pp)
        hJewel_pp_draw.SetLineWidth(4)
        plot_jewel_pp = True
        n_pp_models += 1
      #else:
      #  print('No JEWEL pp prediction for %s %s' % (self.observable, obs_label))

      if hJetscape_pp:
        # Check that binning is the same as data
        b_d = [h.GetBinLowEdge(i) for i in range(1, h.GetNbinsX()+2)]
        b_j = [hJetscape_pp.GetBinLowEdge(i) for i in range(1, hJetscape_pp.GetNbinsX()+2)]
        if b_j != b_d:
          print("JETSCAPE pp BINS ARE DIFFERENT for %s!\n" % obs_label,
                "*** data: ", b_d, "*** jtsp: ", b_j, sep="")
        else:
          color_jetscape = self.colors[3]
          #hJetscape_pp.SetMarkerSize(1.5)
          #hJetscape_pp.SetMarkerStyle(42)
          #hJetscape_pp.SetMarkerColor(44)
          hJetscape_pp.SetFillColorAlpha(color_jetscape, 0.7)
          hJetscape_pp.SetFillStyle(1001)
          hJetscape_pp.SetLineStyle(1)
          hJetscape_pp.SetLineColor(color_jetscape)
          hJetscape_pp.SetLineWidth(0)
          n_pp_models += 1
          plot_jetscape_pp = True
      #else:
      #  print('No JETSCAPE prediction for %s %s' % (self.observable, obs_label))

      if hZhang_pp:
        # Set 0 error for solid line
        for i in range(hZhang_pp.GetNbinsX()+1):
          hZhang_pp.SetBinError(i, 0)
        color_zhang = self.colors[6]
        #hZhang_pp.SetMarkerSize(1.5)
        #hZhang_pp.SetMarkerStyle(27)
        #hZhang_pp.SetMarkerColor(95)
        #hZhang_pp.SetFillColorAlpha(color_zhang, 0.7)
        #hZhang_pp.SetFillStyle(1001)
        hZhang_pp.SetLineStyle(3)
        hZhang_pp.SetLineColor(color_zhang)
        hZhang_pp.SetLineWidth(4)
        n_pp_models += 1
        plot_zhang_pp = True

      if hJewel_no_recoils:
        color_jewel_no_recoils = self.colors[0]
        # Create clone with 0 error for solid line
        hJewel_no_recoils_draw = hJewel_no_recoils.Clone(hJewel_no_recoils.GetName()+"_drawclone")
        for i in range(hJewel_no_recoils_draw.GetNbinsX()+1):
          hJewel_no_recoils_draw.SetBinError(i, 0)        #hJewel_no_recoils.SetMarkerSize(1.5)
        #hJewel_no_recoils.SetMarkerStyle(34)
        #hJewel_no_recoils.SetMarkerColor(42)
        #hJewel_no_recoils.SetFillColorAlpha(color_jewel_no_recoils, 0.5)
        #hJewel_no_recoils.SetFillStyle(1001)
        hJewel_no_recoils_draw.SetLineStyle(1)
        hJewel_no_recoils_draw.SetLineColor(color_jewel_no_recoils)
        hJewel_no_recoils_draw.SetLineWidth(4)
        plot_jewel_no_recoils = True
        n_AA_models += 1
      #else:
      #  print('No JEWEL (recoils off) prediction for %s %s' % (self.observable, obs_label))

      if hJewel_recoils:
        color_jewel_recoils = self.colors[1]
        # Create clone with 0 error for solid line
        hJewel_recoils_draw = hJewel_recoils.Clone(hJewel_recoils.GetName()+"_drawclone")
        for i in range(hJewel_recoils_draw.GetNbinsX()+1):
          hJewel_recoils_draw.SetBinError(i, 0)        #hJewel_no_recoils.SetMarkerSize(1.5)
        #hJewel_recoils.SetMarkerSize(1.5)
        #hJewel_recoils.SetMarkerStyle(47)
        #hJewel_recoils.SetMarkerColor(4)
        #hJewel_recoils.SetFillColorAlpha(color_jewel_recoils, 0.5)
        #hJewel_recoils.SetFillStyle(1001)
        hJewel_recoils_draw.SetLineStyle(7)
        hJewel_recoils_draw.SetLineColor(color_jewel_recoils)
        hJewel_recoils_draw.SetLineWidth(4)
        plot_jewel_recoils = True
        n_AA_models += 1
      #else:
      #  print('No JEWEL (recoils on) prediction for %s %s' % (self.observable, obs_label))

      if hJetscape:
        # Check that binning is the same as data
        b_d = [h.GetBinLowEdge(i) for i in range(1, h.GetNbinsX()+2)]
        b_j = [hJetscape.GetBinLowEdge(i) for i in range(1, hJetscape.GetNbinsX()+2)]
        if b_j != b_d:
          print("JETSCAPE AA BINS ARE DIFFERENT for %s!\n" % obs_label,
                "*** data: ", b_d, "*** jtsp: ", b_j, sep="")
        else:
          color_jetscape = self.colors[3]
          #hJetscape.SetMarkerSize(1.5)
          #hJetscape.SetMarkerStyle(43)
          #hJetscape.SetMarkerColor(6)
          hJetscape.SetFillColorAlpha(color_jetscape, 0.7)
          hJetscape.SetFillStyle(1001)
          hJetscape.SetLineStyle(1)
          hJetscape.SetLineColor(color_jetscape)
          hJetscape.SetLineWidth(0)
          plot_jetscape = True
          n_AA_models += 1
      #else:
      #  print('No JETSCAPE prediction for %s %s' % (self.observable, obs_label))

      if hZhang:
        # Set 0 error for solid line
        for i in range(hZhang.GetNbinsX()+1):
          hZhang.SetBinError(i, 0)
        color_zhang = self.colors[6]
        #hZhang.SetMarkerSize(1.5)
        #hZhang.SetMarkerStyle(33)
        #hZhang.SetMarkerColor(2)
        #hZhang.SetFillColorAlpha(color_zhang, 0.7)
        #hZhang.SetFillStyle(1001)
        hZhang.SetLineStyle(3)
        hZhang.SetLineColor(color_zhang)
        hZhang.SetLineWidth(4)
        n_AA_models += 1
        plot_zhang = True

      if hHybridNoElastic:
        plot_hybrid = True
        color_hybrid_no_elastic = self.colors[4]
        color_hybrid_elastic = self.colors[5]

        #hHybridNoElastic.SetMarkerSize(1.5)
        #hHybridNoElastic.SetMarkerStyle(35)
        #hHybridNoElastic.SetMarkerColor(7)
        hHybridNoElastic.SetFillColorAlpha(color_hybrid_no_elastic, 0.5)
        hHybridNoElastic.SetFillStyle(1001)
        hHybridNoElastic.SetLineStyle(1)
        hHybridNoElastic.SetLineColor(color_hybrid_no_elastic)
        hHybridNoElastic.SetLineWidth(0)

        #hHybridWithElastic.SetMarkerSize(1.5)
        #hHybridWithElastic.SetMarkerStyle(36)
        #hHybridWithElastic.SetMarkerColor(8)
        hHybridWithElastic.SetFillColorAlpha(color_hybrid_elastic, 0.5)
        hHybridWithElastic.SetFillStyle(1001)
        hHybridWithElastic.SetLineStyle(1)
        hHybridWithElastic.SetLineColor(color_hybrid_elastic)
        hHybridWithElastic.SetLineWidth(0)

        n_AA_models += 2

        if hHybridNoElastic_pp and plot_pp_data:
          plot_hybrid_pp = True
          color_hybrid_pp = color_hybrid_no_elastic

          #hHybridNoElastic_pp.SetMarkerSize(1.5)
          #hHybridNoElastic_pp.SetMarkerStyle(35)
          #hHybridNoElastic_pp.SetMarkerColor(7)
          hHybridNoElastic_pp.SetFillColorAlpha(color_hybrid_pp, 0.5)
          hHybridNoElastic_pp.SetFillStyle(1001)
          hHybridNoElastic_pp.SetLineStyle(8)
          hHybridNoElastic_pp.SetLineColor(color_hybrid_pp)
          hHybridNoElastic_pp.SetLineWidth(0)

          #hHybridWithElastic_pp.SetMarkerSize(1.5)
          #hHybridWithElastic_pp.SetMarkerStyle(36)
          #hHybridWithElastic_pp.SetMarkerColor(8)
          hHybridWithElastic_pp.SetFillColorAlpha(color_hybrid_pp, 0.5)
          hHybridWithElastic_pp.SetFillStyle(1001)
          hHybridWithElastic_pp.SetLineStyle(8)
          hHybridWithElastic_pp.SetLineColor(color_hybrid_pp)
          hHybridWithElastic_pp.SetLineWidth(0)

          n_pp_models += 2

    h_pp_data = None;  h_pp_sys = None
    if plot_pp_data:
      h_pp_data, h_pp_sys = self.get_pp_data(
        jetR, obs_label, min_pt_truth, max_pt_truth,
        [h.GetBinLowEdge(i) for i in range(1, h.GetNbinsX()+2)])

      h_pp_data.SetMarkerSize(1.5)
      h_pp_data.SetMarkerStyle(4) #27)
      h_pp_data.SetMarkerColor(1)
      h_pp_data.SetFillColor(1)
      h_pp_data.SetLineStyle(9)
      h_pp_data.SetLineWidth(2)
      h_pp_data.SetLineColor(1)
      h_pp_sys.SetLineColor(0)
      h_pp_sys.SetLineWidth(0)
      #h_pp_sys.SetFillColor(1)
      h_pp_sys.SetFillColorAlpha(color, 0.8)
      h_pp_sys.SetFillStyle(3004)

    if self.observable != "mass" and not grooming_setting:
      maxval = max(2.3*h.GetBinContent(int(0.4*h.GetNbinsX())), 1.7*h.GetMaximum())
      if plot_jetscape:
        maxval = max(maxval, 1.7*hJetscape.GetMaximum())
      if min_pt_truth == 100:
        alpha = obs_label.split("_")[0]
        if alpha == "1.5":
          maxval *= 1.2
    else:
      maxval = 2.3*max(h.GetBinContent(int(0.4*h.GetNbinsX())), h.GetBinContent(2), h.GetBinContent(3))
    ymin = 1e-3  # Prevent ROOT from drawing 0 on plots
    if self.set_logy:
      maxval *= 5e1
      ymin = 5e-1 * h.GetMinimum()
    myBlankHisto.SetMinimum(ymin)
    myBlankHisto.SetMaximum(maxval)
    myBlankHisto.Draw("E")

    if match_data_normalization:
      m = 2 if grooming_setting else 1
      integral = h.Integral(m, h.GetNbinsX(), "width")
      if plot_jewel_no_recoils:
        hJewel_no_recoils_draw.Scale(integral / hJewel_no_recoils.Integral(m, h.GetNbinsX(), "width"))
        hJewel_no_recoils.Scale(integral / hJewel_no_recoils.Integral(m, h.GetNbinsX(), "width"))
      if plot_jewel_recoils:
        hJewel_recoils_draw.Scale(integral / hJewel_recoils.Integral(m, h.GetNbinsX(), "width"))
        hJewel_recoils.Scale(integral / hJewel_recoils.Integral(m, h.GetNbinsX(), "width"))
      if plot_jetscape:
        hJetscape.Scale(integral / hJetscape.Integral(m, h.GetNbinsX(), "width"))
      if plot_zhang:
        hZhang.Scale(integral / hZhang.Integral(m, h.GetNbinsX(), "width"))
      if plot_hybrid:
        hHybridNoElastic.Scale(integral / hHybridNoElastic.Integral(m, h.GetNbinsX(), "width"))
        hHybridWithElastic.Scale(integral / hHybridWithElastic.Integral(m, h.GetNbinsX(), "width"))

      if plot_pp_data:
        integral = h_pp_data.Integral(m, h_pp_data.GetNbinsX(), "width")
        if plot_pythia:
          hPythia_draw.Scale(integral / hPythia.Integral(m, h.GetNbinsX(), "width"))
          hPythia.Scale(integral / hPythia.Integral(m, h.GetNbinsX(), "width"))
        if plot_herwig:
          hHerwig_draw.Scale(integral / hHerwig.Integral(m, h.GetNbinsX(), "width"))
          hHerwig.Scale(integral / hHerwig.Integral(m, h.GetNbinsX(), "width"))
        if plot_jewel_pp:
          hJewel_pp_draw.Scale(integral / hJewel_pp.Integral(m, h.GetNbinsX(), "width"))
          hJewel_pp.Scale(integral / hJewel_pp.Integral(m, h.GetNbinsX(), "width"))
        if plot_jetscape_pp:
          hJetscape_pp.Scale(integral / hJetscape_pp.Integral(m, h.GetNbinsX(), "width"))
        if plot_zhang_pp:
          hZhang_pp.Scale(integral / hZhang_pp.Integral(m, h.GetNbinsX(), "width"))
        if plot_hybrid:
          hHybridNoElastic_pp.Scale(integral / hHybridNoElastic_pp.Integral(m, h.GetNbinsX(), "width"))
          hHybridWithElastic_pp.Scale(integral / hHybridWithElastic_pp.Integral(m, h.GetNbinsX(), "width"))

    if not ignore_MC_top_panel:
      if not plot_PbPb:
        if plot_pythia:
          hPythia_draw.Draw('L 3 same')
        if plot_herwig:
          hHerwig_draw.Draw('L 3 same')
        if plot_jewel_pp:
          hJewel_pp_draw.Draw('L 3 same')
        if plot_zhang_pp:
          hZhang_pp.Draw('L 3 same')
        if plot_hybrid_pp:
          hHybridNoElastic_pp.Draw('E3 same')
          #hHybridWithElastic_pp.Draw('E2 same')
        if plot_jetscape_pp:
          hJetscape_pp.Draw('E3 same')
      else:
        if plot_jewel_no_recoils:
          hJewel_no_recoils_draw.Draw('L 3 same')
        if plot_jewel_recoils:
          hJewel_recoils_draw.Draw('L 3 same')
        if plot_zhang:
          hZhang.Draw('L 3 same')
        if plot_hybrid:
          hHybridNoElastic.Draw('E3 same')
          hHybridWithElastic.Draw('E3 same')
        if plot_jetscape:
          hJetscape.Draw('E3 same')

    if plot_pp_data:
      h_pp_sys.Draw("E2 same")
      h_pp_data.Draw("PE X0 same")

    if plot_PbPb:
      h_sys.Draw("E2 same")
      h.Draw("PE X0 same")

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    #text_xval = 0.53 if plot_PbPb else 0.56
    text_xval = 0.61
    text_yval = 0.9;  delta_y = 0.065
    if self.observable == "ang" and plot_pp_data and not plot_PbPb:
      text = 'ALICE'
    else:
      text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(text_xval, text_yval, text)
    text_yval -= delta_y

    text_latex.SetTextSize(0.045)
    if plot_PbPb:
      if not plot_pp_data:
        text = '0-10% centrality Pb-Pb'
        text_latex.DrawLatex(text_xval, text_yval, text)
        text_yval -= delta_y
      text = '#sqrt{#it{s}_{NN}} = 5.02 TeV'
      text_latex.DrawLatex(text_xval, text_yval, text)
    else:
      text = "pp #sqrt{#it{s}} = 5.02 TeV"
      text_latex.DrawLatex(text_xval, text_yval, text)
    text_yval -= delta_y

    text = "Ch.-particle anti-#it{k}_{T} jets"
    text_latex.DrawLatex(text_xval, text_yval, text)
    text_yval -= delta_y
    #text_xval += 0.1 if plot_PbPb else 0.07

    text = str(min_pt_truth) + ' < #it{p}_{T}^{ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.DrawLatex(text_xval, text_yval, text)
    text_yval -= delta_y

    text = '| #it{#eta}_{jet}| < %s' % str(0.9 - jetR)
    subobs_label = self.utils.formatted_subobs_label(self.observable)
    #if subobs_label:
    #  text += ',   %s = %s' % (subobs_label, obs_setting)
    text += ',   #it{R} = ' + str(jetR)
    text_latex.DrawLatex(text_xval, text_yval, text)
    text_yval -= delta_y

    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting) #.replace("#beta}", "#beta}_{SD}")
      text_latex.DrawLatex(text_xval, text_yval, text)
      text_yval -= delta_y

      if not match_data_normalization:
        text_latex.SetTextSize(0.04)
        text = ['#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged]
        if plot_pythia:
          text.append('#it{f}_{tagged}^{PYTHIA} = %3.3f' % fraction_tagged_pythia)
        if plot_herwig:
          text.append('#it{f}_{tagged}^{Herwig} = %3.3f' % fraction_tagged_herwig)
        if plot_jewel_pp:
          text.append('#it{f}_{tagged}^{JEWEL pp} = %3.3f' % fraction_tagged_jewel_pp)
        if plot_jewel_no_recoils:
          text.append('#it{f}_{tagged}^{JEWEL} = %3.3f' % fraction_tagged_jewel_no_recoils)
        if plot_jewel_recoils:
          text.append('#it{f}_{tagged}^{JEWEL rec.} = %3.3f' % fraction_tagged_jewel_recoils)
        if plot_jetscape:
          text.append('#it{f}_{tagged}^{JETSCAPE} = %3.3f' % fraction_tagged_jetscape)
        if plot_jetscape:
          text.append('#it{f}_{tagged}^{JETSCAPE pp} = %3.3f' % fraction_tagged_jetscape_pp)

        # Print two of each f_tagged on a line
        for list_i, sublist in enumerate([text[i:i+2] for i in range(0, len(text), 2)]):
          line = ", ".join(sublist)
          text_latex.DrawLatex(text_xval, text_yval, line)
          text_yval -= delta_y

    maxy = 0.94
    #miny = maxy - 0.07 * n_AA_models if plot_PbPb else maxy - 0.07 * n_pp_models
    miny = 0.65 if (plot_PbPb and plot_pp_data) else 0.55
    myLegend = ROOT.TLegend(0.18, miny, text_xval-0.02, maxy)
    self.utils.setup_legend(myLegend, 0.035)
    if plot_PbPb:
      if plot_pp_data:
        myLegend.AddEntry(h, 'ALICE 0-10% Pb-Pb data', 'pe')
      else:
        myLegend.AddEntry(h, 'ALICE Pb-Pb data', 'pe')
      myLegend.AddEntry(h_sys, 'Pb-Pb syst. uncert.', 'f')
    if plot_pp_data:
      myLegend.AddEntry(h_pp_data, 'ALICE pp data', 'pe')
      myLegend.AddEntry(h_pp_sys, 'pp syst. uncert.', 'f')
    if not ignore_MC_top_panel:
      if plot_pythia:
        myLegend.AddEntry(hPythia_draw, 'PYTHIA8 Monash2013', 'l')
      if plot_herwig:
        myLegend.AddEntry(hHerwig_draw, 'Herwig7 default tune', 'l')
      if plot_jewel_pp:
        myLegend.AddEntry(hJewel_pp_draw, 'JEWEL pp', 'l')
      if plot_zhang_pp:
        myLegend.AddEntry(hZhang_pp, 'POWHEG+PYTHIA6', 'l')
      if plot_hybrid_pp:
        myLegend.AddEntry(hHybridNoElastic_pp, 'Hybrid model vacuum', 'f')
        #myLegend.AddEntry(hHybridWithElastic_pp, 'Hybrid model (with elastic) baseline')
      if plot_jetscape_pp:
        myLegend.AddEntry(hJetscape_pp, 'JETSCAPE pp', 'f')
      if plot_PbPb:
        if plot_jewel_no_recoils:
          myLegend.AddEntry(hJewel_no_recoils_draw, 'JEWEL (recoils off)', 'l')
        if plot_jewel_recoils:
          myLegend.AddEntry(hJewel_recoils_draw, 'JEWEL (recoils on)', 'l')
        if plot_zhang:
          myLegend.AddEntry(hZhang, 'Higher-Twist parton #it{E}-loss', 'l')
        if plot_jetscape:
          myLegend.AddEntry(hJetscape, 'JETSCAPE (MATTER+LBT)', 'f')
        if plot_hybrid:
          myLegend.AddEntry(hHybridNoElastic, 'Hybrid model (no elastic)', 'f')
          myLegend.AddEntry(hHybridWithElastic, 'Hybrid model (with elastic)', 'f')
    myLegend.Draw()

    ##########################################################################
    # Make ratio plot if desired

    if make_ratio_plot:
      c.cd()
      pad2 = ROOT.TPad('pad2', 'Ratio pad', 0, 0, 1, pad_y_split)
      pad2.SetLeftMargin(pad_left_margin)
      pad2.SetTopMargin(0)
      pad2.SetRightMargin(0.04)
      pad2.SetBottomMargin(0.13)
      #if self.set_logy:
      #  pad2.SetLogy()
      pad2.SetTicks(1, 1)
      pad2.Draw()
      pad2.cd()

      myBlankHisto2 = ROOT.TH1F('myBlankHisto2','Blank Histogram for Ratio', n_obs_bins_truth, truth_bin_array)
      myBlankHisto2.SetNdivisions(505)
      myBlankHisto2.SetXTitle(self.xtitle.replace("alpha}", "alpha}="+alpha))
      myBlankHisto2.GetXaxis().SetTitleOffset(1.02)
      myBlankHisto2.GetXaxis().SetTitleSize(0.055)
      if plot_pp_data and plot_PbPb:
        myBlankHisto2.SetYTitle("#frac{Pb-Pb}{pp}")
      elif plot_MC:
        myBlankHisto2.SetYTitle("#frac{Theory}{Data}")
      myBlankHisto2.GetYaxis().SetTitleOffset(1.3 if self.set_logy else 1.15)
      myBlankHisto2.GetYaxis().SetTitleSize(0.055)
      if plot_pp_data and plot_PbPb:
        myBlankHisto2.SetMinimum(0.3)
        if self.observable == "ang":
          myBlankHisto2.SetMaximum(2.6)
        else:
          myBlankHisto2.SetMaximum(1.99)
      elif plot_MC:
        if min_pt_truth == 100:
          myBlankHisto2.SetMinimum(0.35)
        else:
          myBlankHisto2.SetMinimum(0.5)
        if plot_pp_data:
          myBlankHisto2.SetMaximum(1.85)
        else:
          if min_pt_truth == 100:
            myBlankHisto2.SetMaximum(1.65)
          else:
            myBlankHisto2.SetMaximum(1.55)
      myBlankHisto2.Draw("E")

      # Draw dashed line at ratio = 1
      line = ROOT.TLine(0, 1, h.GetBinLowEdge(h.GetNbinsX()+1), 1)
      line.SetLineColor(920+2)
      line.SetLineStyle(2)
      line.Draw()

      # pp MC / pp data
      if not plot_PbPb:
        h_pp_data_no_error = h_pp_data.Clone(h_pp_data.GetName()+"_no_error")
        for i in range(1, h_pp_data.GetNbinsX()+1):
          h_pp_data_no_error.SetBinError(i, 0)

        if plot_pythia:
          hPythiaRatio = hPythia.Clone(hPythia.GetName()+"_ratio")
          hPythiaRatio.Divide(h_pp_data_no_error)
          hPythiaRatio.SetMarkerSize(0)
          hPythiaRatio.SetFillColorAlpha(color_pythia, 1)
          hPythiaRatio.SetFillStyle(1001)
          hPythiaRatio.SetLineStyle(1)
          hPythiaRatio.SetLineColor(color_pythia)
          hPythiaRatio.SetLineWidth(4)
          hPythiaRatio.Draw('E3 same')

        if plot_herwig:
          hHerwigRatio = hHerwig.Clone(hHerwig.GetName()+"_ratio")
          hHerwigRatio.Divide(h_pp_data_no_error)
          hHerwigRatio.SetMarkerSize(0)
          hHerwigRatio.SetFillColorAlpha(color_herwig, 0.5)
          hHerwigRatio.SetFillStyle(1001)
          hHerwigRatio.SetLineStyle(1)
          hHerwigRatio.SetLineColor(color_herwig)
          hHerwigRatio.SetLineWidth(4)
          hHerwigRatio.Draw('E3 same')

        if plot_jewel_pp:
          hJewelRatio = hJewel_pp.Clone(hJewel_pp.GetName()+"_ratio")
          hJewelRatio.Divide(h_pp_data_no_error)
          hJewelRatio.SetMarkerSize(0)
          hJewelRatio.SetFillColorAlpha(color_jewel_pp, 0.5)
          hJewelRatio.SetFillStyle(1001)
          hJewelRatio.SetLineStyle(1)
          hJewelRatio.SetLineColor(color_jewel_pp)
          hJewelRatio.SetLineWidth(4)
          hJewelRatio.Draw('E3 same')

        if plot_hybrid_pp:
          hHybridNoElasticRatio = hHybridNoElastic_pp.Clone(hHybridNoElastic_pp.GetName()+"_ratio")
          hHybridNoElasticRatio.Divide(h_pp_data_no_error)
          hHybridNoElasticRatio.SetMarkerSize(0)
          hHybridNoElasticRatio.SetFillColorAlpha(color_hybrid_pp, 0.5)
          hHybridNoElasticRatio.SetFillStyle(1001)
          hHybridNoElasticRatio.SetLineStyle(1)
          hHybridNoElasticRatio.SetLineColor(color_hybrid_pp)
          hHybridNoElasticRatio.SetLineWidth(4)
          hHybridNoElasticRatio.Draw('E3 same')

          #hHybridWithElasticRatio = hHybridWithElastic_pp.Clone(hHybridWithElastic_pp.GetName()+"_ratio")
          #hHybridWithElasticRatio.Divide(h_pp_data)
          #hHybridWithElasticRatio.Draw('E2 same')

        if plot_jetscape_pp:
          hJetscapeRatio = hJetscape_pp.Clone(hJetscape_pp.GetName()+"_ratio")
          hJetscapeRatio.Divide(h_pp_data_no_error)
          hJetscapeRatio.SetMarkerSize(0)
          hJetscapeRatio.SetFillColorAlpha(color_jetscape, 0.7)
          hJetscapeRatio.SetFillStyle(1001)
          hJetscapeRatio.SetLineStyle(1)
          hJetscapeRatio.SetLineColor(color_jetscape)
          hJetscapeRatio.SetLineWidth(4)
          hJetscapeRatio.Draw('E3 same')

        if plot_zhang_pp:
          hZhangRatio = hZhang_pp.Clone(hZhang_pp.GetName()+"_ratio")
          hZhangRatio.Divide(h_pp_data_no_error)
          hZhangRatio.SetMarkerSize(0)
          #hZhangRatio.SetFillColorAlpha(color_zhang, 0.7)
          #hZhangRatio.SetFillStyle(1001)
          hZhangRatio.SetLineStyle(3)
          hZhangRatio.SetLineColor(color_zhang)
          hZhangRatio.SetLineWidth(4)
          hZhangRatio.Draw('L 3 same')

        # Draw data stat + systematic uncertainties around ratio = 1
        hSysRatio = h_pp_sys.Clone(h_pp_sys.GetName()+"_ratio")
        for i in range(0, hSysRatio.GetNbinsX()+1):
          old_content = hSysRatio.GetBinContent(i)
          if old_content != 0:
            combined_error = math.sqrt(
              hSysRatio.GetBinError(i) * hSysRatio.GetBinError(i) + \
              h_pp_data.GetBinError(i) * h_pp_data.GetBinError(i) )
            hSysRatio.SetBinError(i, combined_error/old_content)
          hSysRatio.SetBinContent(i, 1)
        hSysRatio.SetLineColor(0)
        hSysRatio.SetFillColor(color)
        hSysRatio.SetFillColorAlpha(color, 0.1)
        hSysRatio.SetFillStyle(1001)
        hSysRatio.SetLineWidth(0)
        hSysRatio.Draw("E2 same")

      # Pb-Pb/pp data ratios
      elif plot_pp_data:
        h_sys_ppratio = h_sys.Clone(h_sys.GetName()+"_ppratio")
        h_sys_ppratio.Divide(h_pp_sys)
        h_sys_ppratio.Draw('E2 same')

        h_data_ppratio = h.Clone(h.GetName()+"_ppratio")
        h_data_ppratio.Divide(h_pp_data)
        h_data_ppratio.Draw('PE X0 same')

        # Theory ratios
        if plot_jewel_no_recoils:
          hJewelNoRecoilsRatio = hJewel_no_recoils.Clone(hJewel_no_recoils.GetName()+"_ppratio")
          hJewelNoRecoilsRatio.Divide(hJewel_pp)
          hJewelNoRecoilsRatio.SetMarkerSize(0)
          hJewelNoRecoilsRatio.SetFillColorAlpha(color_jewel_no_recoils, 1)
          hJewelNoRecoilsRatio.SetFillStyle(1001)
          hJewelNoRecoilsRatio.SetLineStyle(1)
          hJewelNoRecoilsRatio.SetLineColor(color_jewel_no_recoils)
          hJewelNoRecoilsRatio.SetLineWidth(0)
          hJewelNoRecoilsRatio.Draw('E3 same')

        if plot_jewel_recoils:
          hJewelRecoilsRatio = hJewel_recoils.Clone(hJewel_recoils.GetName()+"_ppratio")
          hJewelRecoilsRatio.Divide(hJewel_pp)
          hJewelRecoilsRatio.SetMarkerSize(0)
          hJewelRecoilsRatio.SetFillColorAlpha(color_jewel_recoils, 0.6)
          hJewelRecoilsRatio.SetFillStyle(1001)
          hJewelRecoilsRatio.SetLineStyle(1)
          hJewelRecoilsRatio.SetLineColor(color_jewel_recoils)
          hJewelRecoilsRatio.SetLineWidth(0)
          hJewelRecoilsRatio.Draw('E3 same')

        if plot_jetscape:
          hJetscapeRatio = hJetscape.Clone(hJetscape.GetName()+"_ppratio")
          hJetscapeRatio.Divide(hJetscape_pp)
          hJetscapeRatio.SetMarkerSize(0)
          hJetscapeRatio.SetFillColorAlpha(color_jetscape, 0.5)
          hJetscapeRatio.SetFillStyle(1001)
          hJetscapeRatio.SetLineStyle(1)
          hJetscapeRatio.SetLineColor(color_jetscape)
          hJetscapeRatio.SetLineWidth(0)
          hJetscapeRatio.Draw('E3 same')

        if plot_zhang:
          hZhangRatio = hZhang.Clone(hZhang.GetName()+"_ppratio")
          hZhangRatio.Divide(hZhang_pp)
          hZhangRatio.SetMarkerSize(0)
          #hZhangRatio.SetFillColorAlpha(color_zhang, 0.7)
          #hZhangRatio.SetFillStyle(1001)
          hZhangRatio.SetLineStyle(3)
          hZhangRatio.SetLineColor(color_zhang)
          hZhangRatio.SetLineWidth(4)
          hZhangRatio.Draw('L 3 same')

        if plot_hybrid:
          hHybridNoElasticRatio = hHybridNoElastic.Clone(hHybridNoElastic.GetName()+"_ppratio")
          hHybridNoElasticRatio.Divide(hHybridNoElastic_pp)
          hHybridNoElasticRatio.SetMarkerSize(0)
          hHybridNoElasticRatio.SetFillColorAlpha(color_hybrid_no_elastic, 0.3)
          hHybridNoElasticRatio.SetFillStyle(1001)
          hHybridNoElasticRatio.SetLineStyle(1)
          hHybridNoElasticRatio.SetLineColor(color_hybrid_no_elastic)
          hHybridNoElasticRatio.SetLineWidth(0)
          hHybridNoElasticRatio.Draw('E3 same')

          hHybridWithElasticRatio = hHybridWithElastic.Clone(hHybridWithElastic.GetName()+"_ppratio")
          hHybridWithElasticRatio.Divide(hHybridWithElastic_pp)
          hHybridWithElasticRatio.SetMarkerSize(0)
          hHybridWithElasticRatio.SetFillColorAlpha(color_hybrid_elastic, 0.4)
          hHybridWithElasticRatio.SetFillStyle(1001)
          hHybridWithElasticRatio.SetLineStyle(1)
          hHybridWithElasticRatio.SetLineColor(color_hybrid_elastic)
          hHybridWithElasticRatio.SetLineWidth(0)
          hHybridWithElasticRatio.Draw('E3 same')

      # Calculate MC ratio plots
      elif plot_MC:
        h_AA_no_error = h.Clone(h.GetName() + "_no_error")
        for i in range(1, h.GetNbinsX() + 1):
          h_AA_no_error.SetBinError(i, 0)

        if plot_jewel_no_recoils:
          hJewelNoRecoilsRatio = hJewel_no_recoils.Clone(hJewel_no_recoils.GetName()+"_ratio")
          hJewelNoRecoilsRatio.Divide(h_AA_no_error)
          hJewelNoRecoilsRatio.SetMarkerSize(0)
          hJewelNoRecoilsRatio.SetFillColorAlpha(color_jewel_no_recoils, 1)
          hJewelNoRecoilsRatio.SetFillStyle(1001)
          hJewelNoRecoilsRatio.SetLineStyle(1)
          hJewelNoRecoilsRatio.SetLineColor(color_jewel_no_recoils)
          hJewelNoRecoilsRatio.SetLineWidth(4)
          hJewelNoRecoilsRatio.Draw('E3 same')

        if plot_jewel_recoils:
          hJewelRecoilsRatio = hJewel_recoils.Clone(hJewel_recoils.GetName()+"_ratio")
          hJewelRecoilsRatio.Divide(h_AA_no_error)
          hJewelRecoilsRatio.SetMarkerSize(0)
          hJewelRecoilsRatio.SetFillColorAlpha(color_jewel_recoils, 0.5)
          hJewelRecoilsRatio.SetFillStyle(1001)
          hJewelRecoilsRatio.SetLineStyle(1)
          hJewelRecoilsRatio.SetLineColor(color_jewel_recoils)
          hJewelRecoilsRatio.SetLineWidth(4)
          hJewelRecoilsRatio.Draw('E3 same')

        if plot_jetscape:
          hJetscapeRatio = hJetscape.Clone(hJetscape.GetName()+"_ratio")
          hJetscapeRatio.Divide(h_AA_no_error)
          hJetscapeRatio.SetMarkerSize(0)
          hJetscapeRatio.SetFillColorAlpha(color_jetscape, 0.7)
          hJetscapeRatio.SetFillStyle(1001)
          hJetscapeRatio.SetLineStyle(1)
          hJetscapeRatio.SetLineColor(color_jetscape)
          hJetscapeRatio.SetLineWidth(4)
          hJetscapeRatio.Draw('E3 same')

        if plot_zhang:
          hZhangRatio = hZhang.Clone(hZhang.GetName()+"_ratio")
          hZhangRatio.Divide(h_AA_no_error)
          hZhangRatio.SetMarkerSize(0)
          #hZhangRatio.SetFillColorAlpha(color_zhang, 0.7)
          #hZhangRatio.SetFillStyle(1001)
          hZhangRatio.SetLineStyle(3)
          hZhangRatio.SetLineColor(color_zhang)
          hZhangRatio.SetLineWidth(4)
          hZhangRatio.Draw('L 3 same')

        if plot_hybrid:
          hHybridNoElasticRatio = hHybridNoElastic.Clone(hHybridNoElastic.GetName()+"_ratio")
          hHybridNoElasticRatio.Divide(h_AA_no_error)
          hHybridNoElasticRatio.SetMarkerSize(0)
          hHybridNoElasticRatio.SetFillColorAlpha(color_hybrid_no_elastic, 0.5)
          hHybridNoElasticRatio.SetFillStyle(1001)
          hHybridNoElasticRatio.SetLineStyle(1)
          hHybridNoElasticRatio.SetLineColor(color_hybrid_no_elastic)
          hHybridNoElasticRatio.SetLineWidth(4)
          hHybridNoElasticRatio.Draw('E3 same')

          hHybridWithElasticRatio = hHybridWithElastic.Clone(hHybridWithElastic.GetName()+"_ratio")
          hHybridWithElasticRatio.Divide(h_AA_no_error)
          hHybridWithElasticRatio.SetMarkerSize(0)
          hHybridWithElasticRatio.SetFillColorAlpha(color_hybrid_elastic, 0.5)
          hHybridWithElasticRatio.SetFillStyle(1001)
          hHybridWithElasticRatio.SetLineStyle(1)
          hHybridWithElasticRatio.SetLineColor(color_hybrid_elastic)
          hHybridWithElasticRatio.SetLineWidth(4)
          hHybridWithElasticRatio.Draw('E3 same')

        # Draw data stat + systematic uncertainties around ratio = 1
        hSysRatio = h_sys.Clone(h_sys.GetName()+"_ratio")
        for i in range(0, hSysRatio.GetNbinsX()+1):
          old_content = hSysRatio.GetBinContent(i)
          if old_content != 0:
            combined_error = math.sqrt(
              h_sys.GetBinError(i) * h_sys.GetBinError(i) + \
              h.GetBinError(i) * h.GetBinError(i) )
            hSysRatio.SetBinError(i, combined_error/old_content)
          hSysRatio.SetBinContent(i, 1)
        hSysRatio.SetLineColor(0)
        hSysRatio.SetFillColor(color)
        hSysRatio.SetFillColorAlpha(color, 0.1)
        hSysRatio.SetFillStyle(1001)
        hSysRatio.SetLineWidth(0)
        hSysRatio.Draw("E2 same")

      if ignore_MC_top_panel:
        maxy = 0.97
        miny = maxy - 0.06 * n_AA_models if plot_PbPb else maxy - 0.06 * n_pp_models
        maxx = 0.98
        minx = 0.53
        if plot_pp_data and plot_PbPb and plot_MC and \
          self.observable == "ang" and min_pt_truth == 80 and "SD" not in obs_label:
          maxx -= 0.22; minx -= 0.22
        myLegend2 = ROOT.TLegend(minx, miny, maxx, maxy)
        self.utils.setup_legend(myLegend2, 0.035)
        if plot_jewel_no_recoils:
          myLegend2.AddEntry(hJewelNoRecoilsRatio, 'JEWEL (recoils off)', 'f')
        if plot_jewel_recoils:
          myLegend2.AddEntry(hJewelRecoilsRatio, 'JEWEL (recoils on)', 'f')
        if plot_jetscape:
          myLegend2.AddEntry(hJetscapeRatio, 'JETSCAPE (MATTER+LBT)', 'f')
        if plot_zhang:
          myLegend2.AddEntry(hZhangRatio, 'Higher-Twist parton #it{E}-loss', 'l')
        if plot_hybrid:
          myLegend2.AddEntry(hHybridNoElasticRatio, 'Hybrid model (no elastic)', 'f')
          myLegend2.AddEntry(hHybridWithElasticRatio, 'Hybrid model (with elastic)', 'f')
        myLegend2.Draw()

    ##########################################################################
    # Save plot to output

    name = 'hUnfolded_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), obs_label,
                                             int(min_pt_truth), int(max_pt_truth), self.file_format)

    if not plot_PbPb:
      name = 'hUnfolded_R{}_{}_{}-{}_pp{}'.format(
        self.utils.remove_periods(jetR), obs_label, int(min_pt_truth),
        int(max_pt_truth), self.file_format)

    elif plot_pp_data:
      name = 'hUnfolded_R{}_{}_{}-{}_ppcomp{}'.format(
        self.utils.remove_periods(jetR), obs_label, int(min_pt_truth),
        int(max_pt_truth), self.file_format)

    elif plot_MC:
      name = 'hUnfolded_R{}_{}_{}-{}_MC{}'.format(
        self.utils.remove_periods(jetR), obs_label, int(min_pt_truth),
        int(max_pt_truth), self.file_format)

    name.replace("__", "_")
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
    if plot_jetscape:
      hJetscape.Write()
    if plot_zhang:
      hZhang.Write()
    fFinalResults.Close()


  #----------------------------------------------------------------------
  # Get unfolded data from the previous measurement
  def get_h_prev_result(self, jetR, obs_label, min_pt_truth, max_pt_truth):

    output_dir = getattr(self, 'output_dir_main')

    f = ROOT.TFile(self.results_pp, 'READ')

    # Retrieve previous result and ensure that it has the proper bin range
    h_prev_data_name = ('hmain_%s_R%s_%s_%s-%s_trunc' % \
      (self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)).replace("__", "_")
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

    scale_width = True
    if 'pythia' in MC.lower():
      hMC = self.get_pythia_from_response(
        jetR, obs_label, min_pt_truth, max_pt_truth, maxbin, overlay)
    elif 'herwig' in MC.lower():
      hMC = self.get_herwig_from_response(
        jetR, obs_label, min_pt_truth, max_pt_truth, maxbin, overlay)
    elif "jewel" in MC.lower():
      hMC = self.get_jewel_from_response(
        jetR, obs_label, min_pt_truth, max_pt_truth, maxbin, overlay, MC, recoils)
    elif "jetscape" in MC.lower():
      hMC = self.get_jetscape_from_response(
        jetR, obs_label, min_pt_truth, max_pt_truth, maxbin, overlay, MC)
    elif "zhang" in MC.lower():
      hMC = self.get_zhang_from_file(
        jetR, obs_label, min_pt_truth, max_pt_truth, maxbin, overlay, MC)
    else:
      raise NotImplementedError("Prediction type %s not recognized." % MC)

    if not hMC:
      return [None, None]

    n_jets_inclusive = hMC.Integral(0, hMC.GetNbinsX()+1, "" if scale_width else "width")
    n_jets_tagged = hMC.Integral(hMC.FindBin(
      self.truth_bin_array(obs_label)[0]), hMC.GetNbinsX(), "" if scale_width else "width")

    fraction_tagged_MC =  n_jets_tagged/n_jets_inclusive
    hMC.Scale(1./n_jets_inclusive, "width" if scale_width else "")

    return [hMC, fraction_tagged_MC]

  #----------------------------------------------------------------------
  def get_pythia_from_response(self, jetR, obs_label, min_pt_truth, max_pt_truth,
                               maxbin, overlay=False):

    # Use direct (unmatched) files instead of projecting fastsim RM
    do_direct_files = False #(len(self.theory_predictions) >= (int(recoils) + 1))

    h = None
    if do_direct_files:  # Read from TH2

      f = ROOT.TFile(self.main_response, 'READ')
      name = "h_%s_JetPt_Truth_R%s_%sScaled" % (self.observable, str(jetR), obs_label)
      th2 = f.Get(name)
      if not th2:
        raise AttributeError("%s not found in %s" % (name, self.main_response))
      if not th2.GetSumw2():
        th2.Sumw2()

      # Set range and binning to be the same as data
      name_data = ('hmain_{}_R{}_{}_{}-{}'.format(
        self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)).replace("__","_")
      h_data = getattr(self, name_data)
      pt_bin_array = array('d', [h_data.GetXaxis().GetBinLowEdge(i) for \
                                 i in range(1, h_data.GetNbinsX()+2)])
      obs_bin_array = array('d', [h_data.GetXaxis().GetBinLowEdge(i) for \
                                  i in range(1, h_data.GetNbinsX()+2)])
      move_underflow = (obs_bin_array[0] < 0)

      th2.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
      h = th2.ProjectionY()

      # Finally, rename and truncate the histogram to the correct size
      name = ('hPythia_{}_R{}_{}_{}-{}'.format(
        self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)).replace("__", "_")
      h_rebin = h.Rebin(len(obs_bin_array)-1, name+"_Rebin", obs_bin_array)
      if move_underflow:
        h_rebin.SetBinContent(1, h.GetBinContent(0))
        h_rebin.SetBinError(1, h.GetBinError(0))
      h = self.truncate_hist(h_rebin, None, maxbin, name)
      h.SetDirectory(0)

    else:  # Get projection of RM
      output_dir = getattr(self, 'output_dir_main')

      filepath = os.path.join(output_dir, 'response.root')
      f = ROOT.TFile(filepath, 'READ')

      thn_name = ('hResponse_JetPt_{}_R{}_{}_rebinned'.format(
        self.observable, jetR, obs_label)).replace("__", "_")
      thn = f.Get(thn_name)
      if not thn:
        raise AttributeError("%s not found in %s" % (thn_name, filepath))
      thn.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)

      name = ('hPythia_{}_R{}_{}_{}-{}'.format(
        self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)).replace("__", "_")
      h = self.truncate_hist(thn.Projection(3), None, maxbin, name)
      h.SetDirectory(0)

    return h

  #----------------------------------------------------------------------
  def get_herwig_from_response(self, jetR, obs_label, min_pt_truth, max_pt_truth,
                               maxbin, overlay=False):

    # Use direct (unmatched) files instead of projecting fastsim RM
    do_direct_files = False #(len(self.theory_predictions) >= (int(recoils) + 1))

    h = None

    if do_direct_files:  # Read from TH2

      f = ROOT.TFile(self.fastsim_response_list[1], 'READ')
      name = "h_%s_JetPt_Truth_R%s_%sScaled" % (self.observable, str(jetR), obs_label)
      th2 = f.Get(name)
      if not th2:
        raise AttributeError("%s not found in %s" % (name, self.main_response))
      if not th2.GetSumw2():
        th2.Sumw2()

      # Set range and binning to be the same as data
      name_data = ('hmain_{}_R{}_{}_{}-{}'.format(
        self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)).replace("__", "_")
      h_data = getattr(self, name_data)
      pt_bin_array = array('d', [h_data.GetXaxis().GetBinLowEdge(i) for \
                                 i in range(1, h_data.GetNbinsX()+2)])
      obs_bin_array = array('d', [h_data.GetXaxis().GetBinLowEdge(i) for \
                                  i in range(1, h_data.GetNbinsX()+2)])
      move_underflow = (obs_bin_array[0] < 0)

      th2.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
      h = th2.ProjectionY()

      # Finally, rename and truncate the histogram to the correct size
      name = ('hHerwig_{}_R{}_{}_{}-{}'.format(
        self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)).replace("__", "_")
      h_rebin = h.Rebin(len(obs_bin_array)-1, name+"_Rebin", obs_bin_array)
      if move_underflow:
        h_rebin.SetBinContent(1, h.GetBinContent(0))
        h_rebin.SetBinError(1, h.GetBinError(0))
      h = self.truncate_hist(h_rebin, None, maxbin, name)
      h.SetDirectory(0)

    else:  # Get projection of RM
      try:
        filepath = os.path.join(self.output_dir_fastsim_generator1, 'response.root')
      except AttributeError:  # No fastsim generator
        return None
      f = ROOT.TFile(filepath, 'READ')

      thn_name = 'hResponse_JetPt_{}_R{}_{}_rebinned'.format(
        self.observable, jetR, obs_label).replace("__", "_")
      thn = f.Get(thn_name)
      if not thn:
        raise AttributeError("%s not found in %s" % (thn_name, filepath))
      thn.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)

      name = ('hHerwig_{}_R{}_{}_{}-{}'.format(
        self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)).replace("__", "_")
      h = self.truncate_hist(thn.Projection(3), None, maxbin, name)
      h.SetDirectory(0)

    return h

  #----------------------------------------------------------------------
  def get_jewel_from_response(self, jetR, obs_label, min_pt_truth, max_pt_truth,
                               maxbin, overlay=False, MC="jewel", recoils=False):

    # Use direct (unmatched) files instead of projecting fastsim RM
    do_direct_files = True #(len(self.theory_predictions) >= (int(recoils) + 1))
    # Copy uncertainties them from the fast sim
    copy_errors = False

    # First, try to get from direct calculation files (no fast sim)
    if do_direct_files:
      # Find the theory file corresponding to this prediction
      type = "recoils on" if recoils else "recoils off"
      if "pp" in MC.lower():
        type = "pp"

      filepath = None
      for i, name in enumerate(self.theory_predictions_names):
        if "jewel" in name.lower() and type in name.lower():
          filepath = self.theory_predictions[i]
          break
      if not filepath:
        return None

      f = ROOT.TFile(filepath, 'READ')
      name = "h_%s_JetPt_R%s_%sScaled" % (self.observable, str(jetR), obs_label) if \
        obs_label else "h_%s_JetPt_R%sScaled" % (self.observable, str(jetR))
      th2 = f.Get(name)
      if not th2:
        raise AttributeError("%s not found in %s" % (name, filepath))
      if not th2.GetSumw2():
        th2.Sumw2()

      # Set range and binning to be the same as data
      name_data = 'hmain_{}_R{}_{}_{}-{}'.format(
        self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
      h_data = getattr(self, name_data)
      pt_bin_array = array('d', [h_data.GetXaxis().GetBinLowEdge(i) for \
                                 i in range(1, h_data.GetNbinsX()+2)])
      obs_bin_array = array('d', [h_data.GetXaxis().GetBinLowEdge(i) for \
                                  i in range(1, h_data.GetNbinsX()+2)])
      move_underflow = (obs_bin_array[0] < 0)
      #th2 = self.histutils.rebin_th2(
      #  th2, th2.GetName()+"_rebin", pt_bin_array, len(pt_bin_array)-1,
      #  obs_bin_array, len(obs_bin_array)-1, move_underflow)

      th2.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
      h = th2.ProjectionY()
      #for i in range(1, h.GetNbinsX()+1):
      #  h.SetBinError(i, 1e-10)
      #copy_errors = True

      # Finally, rename and truncate the histogram to the correct size
      name = 'hJewel_{}_{}_R{}_{}_{}-{}'.format(
        type.replace(' ', '_'), self.observable, jetR, obs_label,
        min_pt_truth, max_pt_truth)
      h_rebin = h.Rebin(len(obs_bin_array)-1, name+"_Rebin", obs_bin_array)
      if move_underflow:
        h_rebin.SetBinContent(1, h.GetBinContent(0))
        h_rebin.SetBinError(1, h.GetBinError(0))
      h = self.truncate_hist(h_rebin, None, maxbin, name)
      #print([h.GetBinLowEdge(i) for i in range(1, h.GetNbinsX()+2)])
      h.SetDirectory(0)

    if not do_direct_files or copy_errors:
      # Second, try to get from the fast simulation RM
      gen_to_use = 2 + int(recoils)
      try:
        filepath = os.path.join(
          self.output_dir, self.observable, "fastsim_generator%i" % gen_to_use, 'response.root')
      except AttributeError:
        # Not doing JEWEL for this case
        return None
      f = ROOT.TFile(filepath, 'READ')

      thn_name = 'hResponse_JetPt_{}_R{}_{}_rebinned'.format(
        self.observable, jetR, obs_label).replace("__", "_")
      thn = f.Get(thn_name)
      if not thn:
        raise AttributeError("%s not found in %s" % (thn_name, filepath))
      thn.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)

      name = 'hJewel_recoils_{}_{}_R{}_{}_{}-{}'.format(
        "on" if recoils else "off", self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
      h_proj = self.truncate_hist(thn.Projection(3), None, maxbin, name)
      if copy_errors:
        for i in range(1, h.GetNbinsX()+1):
          h.SetBinError(i, h_proj.GetBinError(i))
      else:
        h = h_proj
        h.SetDirectory(0)

    return h

  #----------------------------------------------------------------------
  def get_jetscape_from_response(self, jetR, obs_label, min_pt_truth, max_pt_truth,
                                 maxbin, overlay=False, MC="jetscape"):

    # Determine type from MC string arg
    type = "AA"  # Default to AA prediction
    if "pp" in MC.lower():
      type = "pp"
    elif "ratio" in MC.lower():
      type = "ratio"

    i_pred = -1
    for i in range(0, len(self.theory_predictions_names)):
      if "jetscape" in self.theory_predictions_names[i].lower():  #and type in self.theory_predictions_names[i].lower():
          i_pred = i
          break

    if i_pred < 0:
      return None
    filepath = self.theory_predictions[i_pred]
    f = ROOT.TFile(filepath, 'READ')

    #if self.observable != "ang":
    #  return None
    gr = "" if "SD" in obs_label else "un"
    if self.observable == "ang":
      alpha = obs_label.split('_')[0]
      name = "h_chjet_angularity_%sgroomed_alice_R%s_pt%i-%i_alpha%s_5020_0-10_%s" % \
        (gr, jetR, min_pt_truth, max_pt_truth, alpha, type)
    elif self.observable == "mass":
      name = "h_chjet_mass_%sgroomed_alice_R%s_pt%i-%i_5020_0-10_%s" % \
        (gr, jetR, min_pt_truth, max_pt_truth, type)

    h = f.Get(name)
    if not h:
      raise AttributeError("%s not found in %s" % (name, filepath))

    name = 'hJetscape_{}_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    width = []
    if "SD" in obs_label:
      # Add negative bin to the groomed plots
      bins = array('d', [-0.1] + [h.GetBinLowEdge(i) for i in range(1, h.GetNbinsX()+2)])
      h_new = ROOT.TH1F(name+"_new", name+"_new", len(bins)-1, bins)
      integral = h.Integral(1, h.GetNbinsX())
      #h_new.SetBinContent(1, 0)
      #h_new.SetBinError(1, 0)
      for i in range(0, h.GetNbinsX()+2):
        width = 1 if i == 0 else abs(h.GetBinLowEdge(i + 1) - h.GetBinLowEdge(i))
        h_new.SetBinContent(i + 1, width * h.GetBinContent(i))
        h_new.SetBinError(i + 1, width * h.GetBinError(i))
      h = self.truncate_hist(h_new, None, maxbin, name)
      h.Scale(1./integral)
    else:
      for i in range(1, h.GetNbinsX()+1):
        width = abs(h.GetBinLowEdge(i + 1) - h.GetBinLowEdge(i))
        h.SetBinContent(i, width * h.GetBinContent(i))
        h.SetBinError(i, width * h.GetBinError(i))
      h = self.truncate_hist(h, None, maxbin, name)
      h.Scale(1./h.Integral())
    h.SetDirectory(0)

    return h

  #----------------------------------------------------------------------
  def get_zhang_from_file(
      self, jetR, obs_label, min_pt_truth, max_pt_truth, maxbin, overlay, MC):

    # At the moment only have R=0.2 predictions
    if str(jetR) != "0.2" or self.observable != "ang":
      return None

    # Determine and open the correct file
    base_dir = "/rstorage/alice/AnalysisResults/ang/PbPb/zhang_predictions"
    gr = "groomed" if "sd" in obs_label.lower() else "ungroomed"
    type = "pp" if "pp" in MC.lower() else "PbPb"
    filename = type + '-' + gr + "-jets.dat"

    # Fixed predictions for groomed jets
    #if gr[0] == 'g':
    #  filename = "new-" + filename

    alpha = obs_label.split('_')[0]
    contents = []
    with open(os.path.join(base_dir, filename), 'r') as f:
      contents = f.readlines()

    # Find the correct set of data -- we have 4 pT bins in each file
    values = self.trim_contents(contents, min_pt_truth, max_pt_truth, alpha)

    # Get the observable binning from the measured data
    name_data = 'hmain_{}_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth).replace("__", "_")
    h_data = self.truncate_hist(getattr(self, name_data), None, maxbin, name_data+"_trunc")
    obs_bin_array = [h_data.GetXaxis().GetBinLowEdge(i) for i in range(1, h_data.GetNbinsX()+2)]
    if "SD" in obs_label:
      values = [0] + values
    obs_bin_array = array('d', obs_bin_array)

    # Create and fill TH1 with values from file
    name = ("hZhang_%s_R%s_%s_%i-%i" % (type, jetR, obs_label, min_pt_truth, max_pt_truth)).replace("__", "_")
    h = ROOT.TH1F(name, name, len(obs_bin_array)-1, obs_bin_array)
    for i in range(1, h.GetNbinsX()+1):
      width = abs(h.GetBinLowEdge(i + 1) - h.GetBinLowEdge(i))
      h.SetBinContent(i, width * values[i-1])
      h.SetBinError(i, 1e-8)
    h.Scale(1./h.Integral())
    del h_data

    return h

  #----------------------------------------------------------------------
  def trim_contents(self, contents, min_pt_truth, max_pt_truth, alpha):

    range = "(%i-%i)GeV" % (min_pt_truth, max_pt_truth)
    column = ["1", "1.5", "2", "3"].index(alpha)

    recording = False
    data = []
    for line in contents:
      # Determine line where the data of interest starts
      if not recording:
        if range in line:
          recording = True
        continue

      # Stop condition
      if line[0] == '-' and len(data):
        break

      ending = column * 10 + 7
      if len(line) >= ending:
        s = line[column*10:ending]
        if s[0].isnumeric():
          data.append(float(s))

    return data

  #----------------------------------------------------------------------
  def get_hybrid(self, jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin):

    # For now only have R=0.2 predictions
    if str(jetR) != "0.2" or self.observable not in ["ang", "mass"]:
      return None, None, None, None

    # Get the observable binning from the measured data
    name_data = 'hmain_{}_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth) #.replace("__","_")
    h_data = self.truncate_hist(getattr(self, name_data), None, maxbin, name_data+"_trunc")
    obs_bin_array = [h_data.GetXaxis().GetBinLowEdge(i) for i in range(1, h_data.GetNbinsX()+2)]
    obs_bin_array = array('d', obs_bin_array)
    del h_data

    filepath = "/rstorage/alice/AnalysisResults/ang/PbPb/hybrid_%s/" % self.observable

    pT_bin = [40, 60, 80, 100].index(min_pt_truth)
    if self.observable == "ang":
      alpha = ["1", "1.5", "2", "3"].index(obs_label.split('_')[0])
    h_pp_NoElastic = None; h_pp_WithElastic = None; h_PbPb_NoElastic = None; h_PbPb_WithElastic = None
    for elastic in ["No", "With"]:  # Elastic Moliere scattering
      filename = ""
      if self.observable == "ang":
        filename = "010_Angus_%sElastic_WantWake_1_JetR_2_Angus_JetBin_%i_Alpha_%i_Groomed_%i.dat" % \
          (elastic, pT_bin, alpha, int("SD" in obs_label))
      elif self.observable == "mass":
        WantWake = True
        el = elastic if elastic == "No" else ""
        filename = "results_chmass/HYBRID_Hadrons_%sElastic_5020_010_Wake_%i_ChJetMass_JetR0p2_Gro_%i_JetBin_%i.dat" % \
          (el, int(WantWake), int("SD" in obs_label), pT_bin+1)

      # Create and fill TH1 with values from file
      name = "hHybrid_%s_%sElastic_pp_R%s_%s_%i-%i" % \
        (self.observable, elastic, jetR, obs_label, min_pt_truth, max_pt_truth)
      h_pp = ROOT.TH1F(name, name, len(obs_bin_array)-1, obs_bin_array)
      name = "hHybrid_%s_%sElastic_PbPb_R%s_%s_%i-%i" % \
        (self.observable, elastic, jetR, obs_label, min_pt_truth, max_pt_truth)
      h_PbPb = ROOT.TH1F(name, name, len(obs_bin_array)-1, obs_bin_array)

      # Read lines from file and fill histograms
      lines = None
      with open(os.path.join(filepath, filename), 'r') as f:
        delin = ' ' if self.observable == "ang" else '\t'
        lines = [line.split(delin) for line in f.readlines()]

      if self.observable == "ang":
        obs_bin_center   = [float(line[0]) for line in lines]
        upper_band_PbPb  = [float(line[1]) for line in lines]
        lower_band_PbPb  = [float(line[2]) for line in lines]
        vacuum           = [float(line[3]) for line in lines]
        error_vacuum     = [float(line[4]) for line in lines]
        upper_band_ratio = [float(line[5]) for line in lines]
        lower_band_ratio = [float(line[6]) for line in lines]

      elif self.observable == "mass":
        obs_bin_center   = [float(line[0]) for line in lines]
        vacuum           = [float(line[1]) for line in lines]
        error_vacuum     = [float(line[2]) for line in lines]
        upper_band_PbPb  = [float(line[3]) for line in lines]
        lower_band_PbPb  = [float(line[4]) for line in lines]
        upper_band_ratio = [float(line[5]) for line in lines]
        lower_band_ratio = [float(line[6]) for line in lines]

      if "SD" in obs_label:
        obs_bin_center   = [-0.05] + obs_bin_center
        upper_band_PbPb  = [0] + upper_band_PbPb
        lower_band_PbPb  = [0] + lower_band_PbPb
        vacuum           = [0] + vacuum
        error_vacuum     = [0] + error_vacuum
        upper_band_ratio = [0] + upper_band_ratio
        lower_band_ratio = [0] + lower_band_ratio

      # Check that obs_bin_center is the same as obs_bin_array
      expected_obs_bin_center = [(obs_bin_array[i+1] - obs_bin_array[i]) / 2 + obs_bin_array[i] for \
                                 i in range(0, len(obs_bin_array)-1)]
      try:
        diff = [abs(expected_obs_bin_center[i] - obs_bin_center[i]) < 1e-8 for i in range(len(expected_obs_bin_center))]
      except:
        raise IndexError("Expected:", expected_obs_bin_center, "\nFound:", obs_bin_center)
      if False in diff:
        raise ValueError("Expected:", expected_obs_bin_center, "\nFound:", obs_bin_center)

      error_PbPb = [(upper_band_PbPb[i] - lower_band_PbPb[i]) / 2 for i in range(0, h_pp.GetNbinsX())]
      center_PbPb = [lower_band_PbPb[i] + error_PbPb[i] for i in range(0, h_pp.GetNbinsX())]

      # Save to TH1 for pp
      for i in range(0, h_pp.GetNbinsX()):
        h_pp.SetBinContent(i+1, vacuum[i])
        h_pp.SetBinError(i+1, error_vacuum[i])
      for i in range(0, h_PbPb.GetNbinsX()):
        h_PbPb.SetBinContent(i+1, center_PbPb[i])
        h_PbPb.SetBinError(i+1, error_PbPb[i])

      if elastic == "No":
        h_pp_NoElastic = h_pp
        h_PbPb_NoElastic = h_PbPb
      else:  # elastic == "With"
        h_pp_WithElastic = h_pp
        h_PbPb_WithElastic = h_PbPb

    return [h_pp_NoElastic, h_pp_WithElastic, h_PbPb_NoElastic, h_PbPb_WithElastic]

  #----------------------------------------------------------------------
  def get_pp_data(self, jetR, obs_label, min_pt_truth, max_pt_truth,
                  xbins, overlay=False):

    output_dir = getattr(self, 'output_dir_main')

    f = ROOT.TFile(self.results_pp, 'READ')

    # Retrieve pp data and ensure that it has the proper bin range
    h_pp_data_name = ('hmain_%s_R%s_%s_%s-%s_trunc' % \
      (self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)).replace('__', '_')
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
    for n_iter in range(2, 20):
      h_pp_sys_name = 'hResult_%s_systotal_R%s_%s_n%i_%s-%s' % \
        (self.observable, jetR, obs_label, n_iter, min_pt_truth, max_pt_truth)
      h_pp_sys = f.Get(h_pp_sys_name)
      if h_pp_sys:
        break
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

      else:
        # Plot data-only plots
        self.plot_observable_overlay_subconfigs(
          i_config, jetR, overlay_list, min_pt_truth,
          max_pt_truth, maxbins)


      # Plot PYTHIA vs pp comparison
      if min_pt_truth == 40 and not self.use_prev_result:
        # Plot PYTHIA comparison plots
        self.plot_observable_overlay_subconfigs(
          i_config, jetR, overlay_list, min_pt_truth, max_pt_truth, maxbins,
          plot_pp_data=True, plot_MC=True, MC='PYTHIA', plot_ratio=True)


  #----------------------------------------------------------------------
  def plot_observable_overlay_subconfigs(self, i_config, jetR, overlay_list, min_pt_truth,
                                         max_pt_truth, maxbins, plot_pp_data=False, plot_MC=False,
                                         MC='PYTHIA', plot_nll=False, plot_ratio=False):

    # Flag to plot ratio all on the same scale, 0 to 2.2
    plot_ratio_same_scale = True

    # In the 2-ratio case, whether to plot just 1 value of alpha
    single_alpha = (self.observable == "mass") or False
    if single_alpha and plot_ratio and plot_pp_data and plot_MC:
      overlay_list = [i for i in overlay_list if i == "config_R0.2_1"]
    if not overlay_list:
      return

    plot_pp_data = plot_pp_data if self.results_pp else False

    name = 'cResult_overlay_R{}_allpt_{}-{}'.format(jetR, min_pt_truth, max_pt_truth)
    if plot_ratio:
      c = ROOT.TCanvas(name, name, 600, 650 if not (plot_pp_data and plot_MC) else 900)
    else:
      c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()

    c.cd()
    minpad1y = 0.5 if (plot_pp_data and plot_MC) else 0.4
    if plot_ratio:
      pad1 = ROOT.TPad('myPad', 'The pad', 0, minpad1y, 1, 1)
    else:
      pad1 = ROOT.TPad('myPad', 'The pad', 0, 0, 1, 1)
    pad1.SetLeftMargin(0.2)
    pad1.SetTopMargin(0.07)
    pad1.SetRightMargin(0.04)
    pad1.SetBottomMargin(0.13)
    pad1.SetTicks(1, 1)
    if plot_ratio:
      pad1.SetBottomMargin(0.)
    # set log y axis if all configs are SD
    setlogy = False
    for name in overlay_list:
      if "SD" not in name:
        break
      if self.observable != "ang":
        break
      elif name == overlay_list[-1]:
        setlogy = True
    pad1.Draw()
    pad1.cd()

    legend_xmin = 0.6
    legend_ymin = 0.52
    legend_ymax = legend_ymin + 0.31
    if plot_ratio:
      if plot_pp_data and plot_MC:
        if single_alpha:
          legend_xmin = 0.52 #0.63
          legend_ymax = 0.91
          legend_ymin = 0.77
        else:
          legend_xmin = 0.5
          legend_ymin = 0.25
          legend_ymax = legend_ymin + 0.31
    else:
      legend_xmin = 0.61
    myLegend = ROOT.TLegend(legend_xmin, legend_ymin, 0.96, legend_ymax)
    self.utils.setup_legend(myLegend, 0.045)

    legend2_ymin = legend_ymin + 0.17
    legend2_xmin = legend_xmin + 0.2
    legend2_ymax = legend_ymin + 0.31
    if plot_ratio:
      if plot_pp_data:
        if plot_MC:
          if single_alpha:
            # Use for MC & pp labels
            legend2_xmin = legend_xmin #- 0.11
            legend2_ymin = legend_ymin - 0.18
            legend2_ymax = legend_ymin
          else:
            legend2_xmin = legend_xmin + 0.25
            legend2_ymin = legend_ymin + 0.21
        else:
          legend2_xmin = legend_xmin + 0.2
          legend2_ymin = legend_ymin + 0.188
    else:
      legend2_ymin = legend_ymin + 0.105

    myLegend2 = ROOT.TLegend(legend2_xmin, legend2_ymin, 0.96, legend2_ymax)
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
        h = self.truncate_hist(getattr(self, name), None, maxbin+1, (name+'_trunc').replace("__", "_"))
      else:
        h = self.truncate_hist(getattr(self, name), None, maxbin, (name+'_trunc').replace("__", "_"))
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
        xmin = self.obs_config_dict[subconfig_name]['obs_bins_truth'][0]
        xmax = self.obs_config_dict[subconfig_name]['obs_bins_truth'][-1]
        if maxbin:
          xmax = self.obs_config_dict[subconfig_name]['obs_bins_truth'][maxbin]
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, xmin, xmax)
        myBlankHisto.SetNdivisions(505)
        xtitle = self.xtitle
        ytitle = self.ytitle
        if single_alpha and plot_ratio and plot_pp_data and plot_MC:
          xtitle = xtitle.replace("#it{#alpha}", "#it{#alpha}=1")
          ytitle = ytitle.replace("#it{#alpha}", "#it{#alpha}=1")
        myBlankHisto.SetXTitle(xtitle)
        myBlankHisto.GetYaxis().SetTitleOffset(1.1 if plot_ratio else 1.4)
        myBlankHisto.SetYTitle(ytitle)
        if jetR == 0.2:
          if min_pt_truth == 20:
            myBlankHisto.SetMaximum(1.1*ymax)
          elif min_pt_truth == 40:
            if self.observable == "mass" and grooming_setting:
              myBlankHisto.SetMaximum(1.9*ymax)
            elif self.use_prev_result:
              myBlankHisto.SetMaximum(1.8*ymax)
            elif plot_pp_data and plot_MC and single_alpha:
              myBlankHisto.SetMaximum(2*ymax)
            else:
              myBlankHisto.SetMaximum(1.7*ymax)
          elif min_pt_truth == 60:
            if self.observable == "mass":
              myBlankHisto.SetMaximum(1.9*ymax)
            else:
              myBlankHisto.SetMaximum(1.15*ymax)
          elif self.observable == "mass":
              myBlankHisto.SetMaximum(1.9*ymax)
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
          if self.observable == "mass":
            myBlankHisto.SetMaximum(1e2*ymax)
          else:
            myBlankHisto.SetMaximum(5*ymax)
          myBlankHisto.SetMinimum(ymin/2)
        elif plot_ratio:
          myBlankHisto.SetMinimum(2e-4) # Don't draw 0 on top panel
        else:
          myBlankHisto.SetMinimum(0.)
        myBlankHisto.GetYaxis().SetTitleSize(0.065)
        myBlankHisto.GetYaxis().SetTitleOffset(1.2)
        myBlankHisto.GetYaxis().SetLabelSize(0.06)
        myBlankHisto.Draw('E')
        if setlogy:
          pad1.SetLogy()

        # Initialize ratio plot
        if plot_ratio:

          c.cd()
          minpad2y = 0.3 if (plot_pp_data and plot_MC) else 0
          pad2 = ROOT.TPad("pad2", "pad2", 0, minpad2y, 1, minpad1y)
          pad2.SetTopMargin(0)
          pad2.SetBottomMargin(0.3 if not (plot_pp_data and plot_MC) else 0)
          pad2.SetLeftMargin(0.2)
          pad2.SetRightMargin(0.04)
          pad2.SetTicks(1, 1)
          pad2.Draw()
          pad2.cd()

          myBlankHisto2 = myBlankHisto.Clone("myBlankHisto_C")
          if plot_MC:
            if self.is_pp:
              myBlankHisto2.SetYTitle("#frac{Data}{%s}" % MC)
            else:
              myBlankHisto2.SetYTitle("#frac{Pb-Pb}{%s}" % MC)
          elif plot_pp_data:
            if self.use_prev_result:
              myBlankHisto2.SetYTitle("#frac{5.02 TeV}{2.76 TeV}")
            else:
              myBlankHisto2.SetYTitle("#frac{Pb-Pb}{pp}")
          myBlankHisto2.GetYaxis().SetTitleSize(20)
          myBlankHisto2.GetYaxis().SetTitleFont(43)
          myBlankHisto2.GetYaxis().SetTitleOffset(3)
          myBlankHisto2.GetYaxis().SetLabelFont(43)
          myBlankHisto2.GetYaxis().SetLabelSize(25)
          myBlankHisto2.GetYaxis().SetNdivisions(505)
          if not plot_pp_data or not plot_MC:
            myBlankHisto2.GetYaxis().SetTitleOffset(2.5)
            myBlankHisto2.SetXTitle(xtitle)
            myBlankHisto2.GetXaxis().SetTitleSize(30)
            myBlankHisto2.GetXaxis().SetTitleFont(43)
            myBlankHisto2.GetXaxis().SetTitleOffset(2.2)
            myBlankHisto2.GetXaxis().SetLabelFont(43)
            myBlankHisto2.GetXaxis().SetLabelSize(25)
          if plot_ratio_same_scale:
            if jetR == 0.2:
              myBlankHisto2.GetYaxis().SetRangeUser(0.2, 2.4)
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

          line = ROOT.TLine(0, 1, xmax, 1)
          line.SetLineColor(920+2)
          line.SetLineStyle(2)
          line.Draw()

          # Initialize third pad if need to also plot pp ratio
          if plot_pp_data and plot_MC:
            c.cd()
            pad3 = ROOT.TPad("pad3", "pad3", 0, 0, 1, minpad2y)
            pad3.SetTopMargin(0)
            pad3.SetBottomMargin(0.4)
            pad3.SetLeftMargin(0.2)
            pad3.SetRightMargin(0.04)
            pad3.SetTicks(1,1)
            pad3.Draw()
            pad3.cd()

            myBlankHisto3 = myBlankHisto2.Clone("myBlankHisto3")
            myBlankHisto3.SetYTitle("#frac{Pb-Pb}{pp}")
            myBlankHisto3.SetXTitle(xtitle)
            myBlankHisto3.GetXaxis().SetTitleSize(30)
            myBlankHisto3.GetXaxis().SetTitleFont(43)
            myBlankHisto3.GetXaxis().SetTitleOffset(3)
            myBlankHisto3.GetXaxis().SetLabelFont(43)
            myBlankHisto3.GetXaxis().SetLabelSize(25)
            myBlankHisto3.GetYaxis().SetTitleSize(20)
            myBlankHisto3.GetYaxis().SetTitleFont(43)
            myBlankHisto3.GetYaxis().SetTitleOffset(3)
            myBlankHisto3.GetYaxis().SetLabelFont(43)
            myBlankHisto3.GetYaxis().SetLabelSize(25)
            myBlankHisto3.GetYaxis().SetNdivisions(505)
            if plot_ratio_same_scale:
              if jetR == 0.2:
                myBlankHisto3.GetYaxis().SetRangeUser(0.2, 2.4)
              elif jetR == 0.4:
                myBlankHisto3.GetYaxis().SetRangeUser(0.5, 1.9)
              else:
                myBlankHisto3.GetYaxis().SetRangeUser(0, 2.2)
            elif jetR == 0.2:
              if min_pt_truth == 20:
                myBlankHisto3.GetYaxis().SetRangeUser(0.6, 1.75)
              elif min_pt_truth == 40:
                myBlankHisto3.GetYaxis().SetRangeUser(0.78, 1.299)
              elif min_pt_truth == 60:
                myBlankHisto3.GetYaxis().SetRangeUser(0.55, 1.499)
              else:
                myBlankHisto3.GetYaxis().SetRangeUser(0.5, 1.99)
            elif jetR == 0.4:
              if min_pt_truth == 20:
                myBlankHisto3.GetYaxis().SetRangeUser(0.81, 1.72)
              elif min_pt_truth == 40:
                myBlankHisto3.GetYaxis().SetRangeUser(0.7, 2.1)
              elif min_pt_truth == 60:
                myBlankHisto3.GetYaxis().SetRangeUser(0.75, 1.55)
              else:
                myBlankHisto3.GetYaxis().SetRangeUser(0.5, 1.99)
            else:
              myBlankHisto3.GetYaxis().SetRangeUser(0.5, 1.99)
            myBlankHisto3.Draw()

            line2 = line.Clone("line2")
            line2.Draw("same")

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

      if plot_pp_data:

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
          hRatioSys = h_sys.Clone('%s_Ratio' % h_sys.GetName())
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

        hRatioStat = h.Clone('%s_Ratio' % h.GetName())
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

        pad2.cd()
        if plot_MC or plot_pp_data:
          # TGraphAsymmErrors doesn't have DrawCopy
          if self.use_prev_result:
            hRatioSys.Draw('E2 same')
            hRatioStat.Draw('PE X0 same')
          else:
            hRatioSys.DrawCopy('E2 same')
            hRatioStat.DrawCopy('PE X0 same')

        # If both plot_MC and plot_pp_data, need to do pad3
        if plot_MC and plot_pp_data:
          pad3.cd()

          hRatioSys2 = h_sys.Clone('%s_Ratio2' % h_sys.GetName())
          hRatioSys2.Divide(h_pp_sys)
          hRatioSys2.SetLineColor(0)
          hRatioSys2.SetFillColor(color)
          hRatioSys2.SetFillColorAlpha(color, 0.3)
          hRatioSys2.SetFillStyle(1001)
          hRatioSys2.SetLineWidth(0)
          hRatioSys2.SetMaximum(1.99)
          hRatioSys2.DrawCopy('E2 same')

          hRatioStat2 = h.Clone('%s_Ratio2' % h.GetName())
          hRatioStat2.Divide(h_pp_data)
          hRatioStat2.SetMarkerSize(1.5)
          hRatioStat2.SetMarkerStyle(marker)
          hRatioStat2.SetMarkerColor(color)
          hRatioStat2.SetLineStyle(1)
          hRatioStat2.SetLineWidth(2)
          hRatioStat2.SetLineColor(color)
          hRatioStat2.SetMaximum(1.99)
          hRatioStat2.DrawCopy('PE X0 same')

      pad1.cd()
      if plot_MC:
        plot_errors = False
        if plot_errors:
          hMC.DrawCopy('E3 same')
        else:
          hMC.DrawCopy('L hist same')

      if plot_pp_data:
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

      subobs_label = self.utils.formatted_subobs_label(self.observable)
      text = ''
      if subobs_label:
        text += subobs_label
        if obs_setting:
          text += ' = ' + str(obs_setting)
      text_list.append(text)
      h_list.append(h)

    pad1.cd()
    for i, h, text in zip(range(len(h_list)), h_list, text_list):
      if i < 2:
        if single_alpha and plot_ratio and plot_pp_data and plot_MC:
          myLegend.AddEntry(h, "0-10% Pb-Pb data") #text + " (girth)", 'pe')
        else:
          myLegend.AddEntry(h, text, 'pe')
      else:
        myLegend2.AddEntry(h, text, 'pe')
    myLegend.AddEntry(h_sys, 'Pb-Pb syst. uncert.', 'f')
    if plot_pp_data:
      if self.use_prev_result:
        myLegend.AddEntry(h_pp_data, 'Pb-Pb @ 2.76 TeV', 'pe')
        if plot_errors:
          myLegend.AddEntry(h_pp_sys, '2.76 TeV syst. uncert.', 'f')
      elif not (plot_MC and single_alpha):
        myLegend.AddEntry(h_pp_data, 'pp data', 'pe')
        if plot_errors:
          myLegend.AddEntry(h_pp_sys, 'pp syst. uncert.', 'f')
      else:
        myLegend2.AddEntry(h_pp_data, 'pp data', 'pe')
        if plot_errors:
          myLegend2.AddEntry(h_pp_sys, 'pp syst. uncert.', 'f')
        if MC.lower() == "pythia":
          myLegend2.AddEntry(hMC, 'PYTHIA8 Monash2013', 'l')
        elif MC.lower() == "herwig":
          myLegend2.AddEntry(hMC, 'Herwig7 Default', 'l')
    if plot_MC and not (plot_pp_data and single_alpha):
      if MC.lower() == "pythia":
        myLegend.AddEntry(hMC, 'PYTHIA8 Monash2013', 'l')
      elif MC.lower() == "herwig":
        myLegend.AddEntry(hMC, 'Herwig7 Default', 'l')

    text_xval = 0.22 if single_alpha and plot_pp_data and plot_MC else 0.27
    if not plot_ratio:
      text_xval = 0.26
    text_yval = 0.85
    delta_y = 0.065 if single_alpha and plot_pp_data and plot_MC else 0.06
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(text_xval, text_yval, text)
    text_yval -= delta_y

    #if single_alpha and plot_ratio and plot_pp_data and plot_MC:
    text = '#sqrt{#it{s}_{NN}} = 5.02 TeV'
    #else:
    #  text = '0-10% Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(text_xval, text_yval, text)
    text_yval -= delta_y

    text = 'Ch.-particle anti-#it{k}_{T} jets'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(text_xval, text_yval, text)
    text_yval -= delta_y

    text = '#it{R} = ' + str(jetR) + ',   | #it{#eta}_{jet}| < %s' % str(0.9 - jetR)
    text_latex.DrawLatex(text_xval, text_yval, text)
    text_yval -= delta_y

    text = str(min_pt_truth) + ' < #it{p}_{T}^{ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(text_xval, text_yval, text)
    text_yval -= delta_y

    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting) #.replace("#beta}", "#beta}_{SD}")
      text_latex.DrawLatex(text_xval, text_yval, text)
      text_yval -= delta_y

    if not (single_alpha and plot_ratio and plot_pp_data and plot_MC):
      text = "0-10% Pb-Pb data"
      xmin = legend_xmin+0.12 if (plot_ratio and plot_pp_data and plot_MC) else legend_xmin+0.09
      if not plot_ratio:
        xmin -= 0.02
      text_latex.DrawLatex(xmin, legend_ymax+0.02, text)

    myLegend.Draw()
    if len(h_list) > 2 or (single_alpha and plot_MC and plot_pp_data):
      myLegend2.Draw()

    name = 'h_{}_R{}_{}-{}_{}{}'.format(self.observable,
                                        self.utils.remove_periods(jetR), int(min_pt_truth),
                                        int(max_pt_truth), i_config, self.file_format)
    if plot_MC:
      if plot_pp_data:
        name = 'h_{}_R{}_{}-{}_ppComp+{}_{}{}'.format(
          self.observable, self.utils.remove_periods(jetR), int(min_pt_truth),
          int(max_pt_truth), MC, i_config, self.file_format)
      else:
        name = 'h_{}_R{}_{}-{}_{}_{}{}'.format(
          self.observable, self.utils.remove_periods(jetR), int(min_pt_truth),
          int(max_pt_truth), MC, i_config, self.file_format)
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
      minbin = 1
      maxbin_adj = maxbin if (maxbin != None) else h.GetNbinsX()
      if 'SD' in obs_label:
        minbin += 1
        if maxbin != None:
          maxbin_adj += 1
      content = [ h.GetBinContent(j) for j in range(minbin, maxbin_adj+1) ]

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
