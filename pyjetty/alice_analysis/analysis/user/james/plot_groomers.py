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

    self.ColorArray = [ROOT.kViolet-8, ROOT.kAzure-4, ROOT.kTeal-8, ROOT.kOrange+6, ROOT.kOrange-3, ROOT.kRed-7, ROOT.kPink+1, ROOT.kCyan-2, ROOT.kGray]
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
        
          # Plot money plot for all observables
          for i, _ in enumerate(self.obs_subconfig_list):

            obs_setting = self.obs_settings[i]
            grooming_setting = self.grooming_settings[i]
            obs_label = self.utils.obs_label(obs_setting, grooming_setting)
            self.create_output_subdir(output_dir, self.utils.grooming_label(grooming_setting))
            


            self.plot_money_plot(observable, jetR, R_max, obs_label, obs_setting, grooming_setting, output_dir)

          # Plot performance plots only once
          if observable == 'theta_g':

            # Create output subdirectories
            output_dir = os.path.join(self.output_dir, 'performance')
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
  def plot_money_plot(self, observable, jetR, R_max, obs_label, obs_setting, grooming_setting, output_dir):
      
    if observable in ['zg', 'theta_g']:
      # (pt, zg, theta_g, flag)
      #  Flag based on where >50% of subleading matched pt resides:
      #    1: subleading
      #    2: leading
      #    3: ungroomed
      #    4: outside
      #    5: other (i.e. 50% is not in any of the above)
      #    6: pp-truth passed grooming, but combined jet failed grooming
      #    7: combined jet passed grooming, but pp-truth failed grooming
      #    8: both pp-truth and combined jet failed SoftDrop
      name = 'h_theta_g_zg_JetPt_R{}_{}_Rmax{}'.format(jetR, obs_label, R_max)
      h4D = self.fMC.Get(name)
      if h4D.GetSumw2() is 0:
        h4D.Sumw2()
        
      name = 'h_theta_g_zg_JetPt_Truth_R{}_{}'.format(jetR, obs_label)
      h3D_truth = self.fMC.Get(name)
      if h3D_truth.GetSumw2() is 0:
        h3D_truth.Sumw2()
        
      # Loop through pt slices, and plot 1D projection onto observable
      for i in range(0, len(self.pt_bins_reported) - 1):
        min_pt_truth = self.pt_bins_reported[i]
        max_pt_truth = self.pt_bins_reported[i+1]
        
        self.plot_obs_projection(observable, h4D.Clone(), h3D_truth.Clone(), jetR, obs_label, obs_setting,
                                 grooming_setting, min_pt_truth, max_pt_truth, output_dir)

  #---------------------------------------------------------------
  def plot_obs_projection(self, observable, h4D, h3D_truth, jetR, obs_label, obs_setting,
                          grooming_setting, min_pt, max_pt, output_dir):
    
    if observable == 'theta_g':
      xmin = 0.
      xmax = 1.
      axis = 2
      rebin_value = 5
    elif observable == 'zg':
      xmin = 0.
      xmax = 0.5
      axis = 1
      rebin_value = 5

    # Optional: cut on dR
    self.min_theta = 0.

    # Set pt range
    h4D.GetAxis(0).SetRangeUser(min_pt, max_pt)
    
    # Normalize by integral
    #integral = h4D.Projection(axis).Integral()
    #if integral > 0:
    #  h4D.Scale(1./integral)
    #else:
    #  print('Integral is 0, check for problem')
    
    # Draw histogram
    c = ROOT.TCanvas('c','c: hist', 600, 650)
    c.cd()
    
    pad1 = ROOT.TPad('myPad', 'The pad',0,0.3,1,1)
    pad1.SetLeftMargin(0.2)
    pad1.SetTopMargin(0.04)
    pad1.SetRightMargin(0.04)
    pad1.SetBottomMargin(0.)
    pad1.Draw()
    pad1.cd()
    
    leg = ROOT.TLegend(0.65,0.5,0.8,0.8)
    self.utils.setup_legend(leg, 0.04)

    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 20, xmin, xmax)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetMinimum(0.01)
    myBlankHisto.GetXaxis().SetTitle(self.xtitle)
    myBlankHisto.GetXaxis().SetTitleSize(30)
    myBlankHisto.GetXaxis().SetTitleFont(43)
    myBlankHisto.GetXaxis().SetTitleOffset(1.95)
    myBlankHisto.GetXaxis().SetLabelFont(43)
    myBlankHisto.GetXaxis().SetLabelSize(25)
    myBlankHisto.GetYaxis().SetTitleSize(25)
    myBlankHisto.GetYaxis().SetTitleFont(43)
    myBlankHisto.GetYaxis().SetTitleOffset(1.7)
    myBlankHisto.GetYaxis().SetLabelFont(43)
    myBlankHisto.GetYaxis().SetLabelSize(25)
    myBlankHisto.GetYaxis().SetNdivisions(505)
    ytitle = '#frac{{d#it{{N}}}}{{d{}}}'.format(self.xtitle)
    myBlankHisto.GetYaxis().SetTitle(ytitle)
    myBlankHisto.Draw('E')
    
    c.cd()
    pad2 = ROOT.TPad('pad2', 'pad2', 0, 0.02, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.2)
    pad2.SetRightMargin(0.04)
    pad2.SetTicks(0,1)
    pad2.Draw()
    pad2.cd()
    
    myBlankHisto2 = myBlankHisto.Clone('myBlankHisto_C')
    ytitle = '#frac{Embedded}{Truth}'
    myBlankHisto2.SetYTitle(ytitle)
    myBlankHisto2.SetXTitle(self.xtitle)
    myBlankHisto2.GetXaxis().SetTitleSize(30)
    myBlankHisto2.GetXaxis().SetTitleFont(43)
    myBlankHisto2.GetXaxis().SetTitleOffset(3.5)
    myBlankHisto2.GetXaxis().SetLabelFont(43)
    myBlankHisto2.GetXaxis().SetLabelSize(25)
    myBlankHisto2.GetYaxis().SetTitleSize(25)
    myBlankHisto2.GetYaxis().SetTitleFont(43)
    myBlankHisto2.GetYaxis().SetTitleOffset(2.)
    myBlankHisto2.GetYaxis().SetLabelFont(43)
    myBlankHisto2.GetYaxis().SetLabelSize(25)
    myBlankHisto2.GetYaxis().SetNdivisions(505)
    myBlankHisto2.SetMinimum(0.61)
    myBlankHisto2.SetMaximum(1.39)
    myBlankHisto2.Draw('')
    
    h_stack = ROOT.THStack('h_stack', 'stacked')
    h_sum = None
    legend_list = ['subleading', 'leading (swap)', 'leading (mis-tag)', 'ungroomed', 'outside', 'other',
                   'combined fail', 'truth fail', 'both fail']
    
    # Loop over each flag
    for i in range(len(legend_list)):
      flag = i+1
      h4D.GetAxis(3).SetRange(flag, flag)
      h4D.GetAxis(2).SetRangeUser(self.min_theta, 1)

      # Project onto 1D
      h1D = h4D.Projection(axis)
      h1D.SetName('h1D_{}'.format(i))
      h1D.Rebin(rebin_value)
        
      h1D.SetLineColor(self.ColorArray[i])
      h1D.SetFillColor(self.ColorArray[i])
      h1D.SetLineWidth(2)

      h_stack.Add(h1D);
      leg.AddEntry(h1D, legend_list[i], 'F')
      
      if i == 0:
        h_sum = h1D.Clone()
        h_sum.Sumw2()
      else:
        h_sum.Add(h1D.Clone())
    
    # Draw truth histogram
    h3D_truth.GetXaxis().SetRangeUser(min_pt, max_pt)
    h3D_truth.GetZaxis().SetRangeUser(self.min_theta, 1)
    if observable == 'theta_g':
      h1D_truth = h3D_truth.Project3D('z')
    elif observable == 'zg':
      h1D_truth = h3D_truth.Project3D('y')
    h1D_truth.Rebin(rebin_value)
    h1D_truth.SetLineColor(1)
    h1D_truth.SetLineWidth(2)
    myBlankHisto.SetMaximum(2.3*h1D_truth.GetMaximum())
    leg.AddEntry(h1D_truth, 'Truth', 'L')
    
    pad1.cd()
    leg.Draw('same')
    h_stack.Draw('same')
    h1D_truth.Draw('same')
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    
    x = 0.23
    y = 0.9
    text_latex.SetTextSize(0.05)
    text = 'PYTHIA8 embedded in thermal background'
    text_latex.DrawLatex(x, y, text)
        
    text = '#sqrt{#it{s_{#it{NN}}}} = 5.02 TeV'
    text_latex.DrawLatex(x, y-0.06, text)

    text = 'Charged jets   anti-#it{k}_{T}'
    text_latex.DrawLatex(x, y-0.12, text)
    
    text = '#it{R} = ' + str(jetR) + '   | #it{#eta}_{jet}| < 0.5'
    text_latex.DrawLatex(x, y-0.18, text)
   
    text = str(int(min_pt)) + ' < #it{p}_{T, ch jet}^{PYTHIA} < ' + str(int(max_pt)) + ' GeV/#it{c}'
    text_latex.DrawLatex(x, y-0.25, text)
   
    text = self.utils.formatted_grooming_label(grooming_setting)
    text_latex.DrawLatex(x, y-0.32, text)
    
    if self.min_theta > 1e-3:
      text = '#Delta#it{{R}} > {}'.format(self.min_theta*jetR)
      text_latex.DrawLatex(x, y-0.38, text)
    
    pad2.cd()
    h_ratio = h_sum.Clone()
    h_ratio.SetMarkerStyle(21)
    h_ratio.SetMarkerColor(1)
    h_ratio.Divide(h1D_truth)
    if h_ratio.GetMaximum() > myBlankHisto2.GetMaximum():
      myBlankHisto2.SetMaximum(1.2*h_ratio.GetMaximum())
    h_ratio.Draw('PE same')
    
    line = ROOT.TLine(xmin,1,xmax,1)
    line.SetLineColor(920+2)
    line.SetLineStyle(2)
    line.Draw('same')

    output_filename = os.path.join(output_dir, '{}/money_plot_{}_{}-{}_dR{}.pdf'.format(self.utils.grooming_label(grooming_setting),
                                          obs_label, min_pt, max_pt, self.utils.remove_periods(self.min_theta*jetR)))
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
