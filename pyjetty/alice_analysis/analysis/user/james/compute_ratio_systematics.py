#! /usr/bin/env python

'''
This class is intended to compute a ratio of two unfolded distributions and the corresponding 
systematic uncertainties on the ratio, taking advantage of partially cancelling systematic uncertainties.

For example, the ratio of subjet z distributions in 0-10% Pb-Pb for two different subjet radii: r=0.1/r=0.2.

We assume that run_analysis.py has been used to unfold the distributions and all systematic variations.
This class is intended to replace only the systematics computation part -- and only for ratios.

Then, you should run your analysis with:

  python compute_ratio_systematics.py /path/to/my/config.yaml

The configuration should contain, as usual, the list of systematics to be computed,
and the locations of the output directory from the run_analysis.py unfolding.

NOTE: This class is still in draft stage, written specific to subjet z analysis. 
      It does not loop through settings (jetR, grooming, subobs, etc), but computes only for a single ratio.


Author: James Mulligan (james.mulligan@berkeley.edu)
'''

import sys
import os
import argparse
from array import *
import numpy as np
import math
import ROOT
import yaml
import hepdata_lib

from pyjetty.alice_analysis.analysis.base import common_base
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class ComputeRatioSystematics(common_base.CommonBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(ComputeRatioSystematics, self).__init__(**kwargs)
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

    self.observable = config['analysis_observable']
    self.jetR = config['jetR'][0]

    # Load needed info from the two subobservable configurations
    obs_config_dict = config[self.observable]
    self.xtitle = obs_config_dict['common_settings']['xtitle']
    self.ytitle = obs_config_dict['common_settings']['ytitle']
    self.pt_bin = obs_config_dict['common_settings']['pt_bins_reported']

    obs_bins_truth = (obs_config_dict['config1']['obs_bins_truth'])
    self.truth_obs_bin_array = np.array(obs_bins_truth)

    # We will assume the numerator is the first setting, and the denominator is the second
    obs_subconfig_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
    self.obs_settings = self.utils.obs_settings(self.observable, obs_config_dict, obs_subconfig_list)
    self.obs_setting_label = self.utils.formatted_subobs_label(self.observable)

    # List of systematic variations to perform
    self.systematics_list = config['systematics_list']

    # Get the two halves of unfolded data -- use one for r=0.1, and one for r=0.2
    input_dir1 = config['output_dir']
    if 'half2' in input_dir1:
      input_dir1 = input_dir1.replace('half2', 'half1')
    input_dir2 = input_dir1.replace('half1', 'half2')
    self.input_dir1 = os.path.join(input_dir1, self.observable)
    self.input_dir2 = os.path.join(input_dir2, self.observable)

    # Create output dir
    self.output_dir = os.path.join(self.input_dir1, 'ratio_systematics')
    if not os.path.exists(self.output_dir):
      os.makedirs(self.output_dir)

    self.file_format = config['file_format']

    print(self)
    print()

  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def compute_ratio_systematics(self):
    print(f'Computing ratio and systematics for ({self.obs_setting_label}={self.obs_settings[0]})/({self.obs_setting_label}={self.obs_settings[1]})...')
    print()

    #-----------------------------------------
    # First, we need to set the regularization parameter to be used
    # (the unfolding machinery saves all n_iter to file, so that we can 
    #  determine it dynamically using the systematic uncertainties)
    #
    # We will just set this manually for each individual result, 
    # to be taken from the run_analysis.py machinery
    #
    # Note: One might think to determine the optimal value using the *ratio* systematics, 
    #       rather than the individual results. However this may cause confusion since then
    #       the ratio will be slightly different than the reported individual results.
    #       So we will just keep the final values from the individual results.
    self.reg_param_numerator = 3
    self.reg_param_denominator = 3
    self.reg_param_generator = 30
    self.numerator_half = 1
    self.denominator_half = 2
    self.min_zr = 0.7
    self.max_zr = 1.
    self.min_pt = 100
    self.max_pt = 150

    #-----------------------------------------
    # Form the ratio of the main result (and save it to file)
    h_numerator_main = self.load_observable('main', self.jetR, self.obs_settings[0], self.reg_param_numerator, half=self.numerator_half)
    h_denominator_main = self.load_observable('main', self.jetR, self.obs_settings[1], self.reg_param_denominator, half=self.denominator_half)

    h_ratio_main = h_numerator_main.Clone()
    h_ratio_main.SetName('h_ratio_main')
    h_ratio_main.Divide(h_denominator_main)

    f = ROOT.TFile(os.path.join(self.input_dir1, 'final_results/fFinalResults_ratio.root'), 'RECREATE')
    h_ratio_main.Write()
    f.Close()

    #-----------------------------------------
    # First pass: Loop through all variations and form ratio, then take ratio to main ratio
    h_ratio_systematic_dict = {}
    for systematic in self.systematics_list:

      if systematic == 'main': # Reg param uncertainty
          h_ratio_systematic_dict['reg_param+'] = self.compute_systematic_variation(h_ratio_main, 'reg_param+')
          h_ratio_systematic_dict['reg_param-'] = self.compute_systematic_variation(h_ratio_main, 'reg_param-')
      elif systematic in ['fastsim_generator0', 'fastsim_generator1']: # We will add these in quadrature
        if systematic == 'fastsim_generator1':
          #h_numerator_generator0 = self.load_observable('fastsim_generator0', self.jetR, self.obs_settings[0], self.reg_param_numerator, half=self.numerator_half)
          #h_denominator_generator0 = self.load_observable('fastsim_generator0', self.jetR, self.obs_settings[1], self.reg_param_denominator, half=self.denominator_half)
          #h_ratio_generator0 = h_numerator_generator0.Clone()
          #h_ratio_generator0.SetName('h_ratio_generator0')
          #h_ratio_generator0.Divide(h_denominator_generator0)
          #h_ratio_systematic_dict['generator'] = self.compute_systematic_variation(h_ratio_generator0, 'generator')
          h_numerator_systematic = self.load_systematic_percentage('generator', self.jetR, self.obs_settings[0], self.reg_param_generator, self.min_pt, self.max_pt, half=self.numerator_half)
          h_denominator_systematic = self.load_systematic_percentage('generator', self.jetR, self.obs_settings[1], self.reg_param_generator, self.min_pt, self.max_pt, half=self.denominator_half)
          h_ratio_systematic_dict['generator'] = self.add_hlist_in_quadrature([h_numerator_systematic, h_denominator_systematic])
      else:
        h_ratio_systematic_dict[systematic] = self.compute_systematic_variation(h_ratio_main, systematic)

    #-----------------------------------------
    # Second pass: Combine related variations into each final "ratio systematic" source

    # First, combine pairs into averaged/max systematic
    # Note: For all of these, we should not keep track of the "signed" uncertainty, since we take abs value before averaging
    h_ratio_systematic_dict_combined = {}
    for systematic in self.systematics_list:

      if systematic == 'main':
        h_ratio_systematic_dict_combined['reg_param'] = self.hist_average(h_ratio_systematic_dict['reg_param+'], h_ratio_systematic_dict['reg_param-'])
      elif systematic in ['prior1', 'prior2']:
        if systematic == 'prior1':
          h_ratio_systematic_dict_combined['prior'] = self.hist_average(h_ratio_systematic_dict['prior1'], h_ratio_systematic_dict['prior1'])
      elif systematic in ['subtraction1', 'subtraction2']:
        if systematic == 'subtraction1':
          h_ratio_systematic_dict_combined['subtraction'] = self.hist_average(h_ratio_systematic_dict['subtraction1'], h_ratio_systematic_dict['subtraction2'])
      elif 'fastsim_generator' in systematic:
        if systematic == 'fastsim_generator0':
          h_ratio_systematic_dict_combined['generator'] = h_ratio_systematic_dict['generator']
      else:
        h_ratio_systematic_dict_combined[systematic] = h_ratio_systematic_dict[systematic]

    # Next, construct the final individual sources
    h_ratio_systematic_dict_grouped = {}
    h_unfolding_list = []
    for sys_label, systematic_histogram in h_ratio_systematic_dict_combined.items():

      if sys_label in ['reg_param', 'prior', 'truncation', 'binning']:
        h_unfolding_list.append(systematic_histogram)
      else:
        h_ratio_systematic_dict_grouped[sys_label] = systematic_histogram

    h_ratio_systematic_dict_grouped['unfolding'] = self.hist_stdev(h_unfolding_list)

    # At this point, the signed uncertainties are: trkeff 
    # And the unsigned uncertainties are: subtraction, unfolding, generator

    #-----------------------------------------
    # Combine systematics in quadrature and write to file
    h_systematic_total = self.add_in_quadrature(h_ratio_systematic_dict_grouped)

    # Attach total systematic to main result and save under name without reg param
    h_ratio_main_sys = h_ratio_main.Clone()
    h_ratio_main_sys.SetName('h_ratio_main_sys')
    self.attach_uncertainty_to_hist(h_ratio_main_sys, h_systematic_total)

    sys_root_filename = os.path.join(self.output_dir, 'fSystematics.root')
    fSystematics = ROOT.TFile(sys_root_filename, 'RECREATE')
    h_ratio_main_sys.Write()
    fSystematics.Close()

    #-----------------------------------------
    # Plot 

    # Plot systematic uncertainties, and write total systematic to a ROOT file
    h_list = list(h_ratio_systematic_dict_grouped.values())
    self.plot_systematic_uncertainties(h_list, h_systematic_total)

    # Plot unfolding uncertainties
    self.plot_systematic_uncertainties(h_unfolding_list, h_ratio_systematic_dict_grouped['unfolding'], suffix='Unfolding')
    print()     

    # Compare "ratio" systematic to quadrature sum systematic

    #-----------------------------------------
    # TODO: Write hepdata submission
    #self.write_hepdata()     

  #----------------------------------------------------------------------
  def load_observable(self, systematic, jetR, obs_setting, reg_param, half=None):

    if half == 1:
      input_dir = self.input_dir1
    elif half == 2:
      input_dir = self.input_dir2

    # Load 2D distribution from file
    input_dir = os.path.join(input_dir, systematic)
    filename = os.path.join(input_dir, f'fResult_R{self.jetR}_{obs_setting}.root')
    f = ROOT.TFile(filename, 'READ')
    hname = f'hUnfolded_{self.observable}_R{jetR}_{obs_setting}_{reg_param}'
    h2D = f.Get(hname)
    h2D.SetDirectory(0)

    # Project to 1D histogram for selected pt range
    h2D.GetXaxis().SetRangeUser(self.pt_bin[0], self.pt_bin[1])
    h2D.GetYaxis().SetRangeUser(self.min_zr, self.max_zr)
    h = h2D.ProjectionY() # Better to use ProjectionY('{}_py'.format(h2D.GetName()), 1, h2D.GetNbinsX()) ?
    h.SetName(f'h{systematic}_{self.observable}_R{jetR}_{obs_setting}_n{reg_param}_{self.pt_bin[0]}-{self.pt_bin[1]}')
    h.SetDirectory(0)
            
    # Scale by bin width -- then normalize by integral
    # (where integral weights by bin width)
    # Note: Sumw2 is already active on these histograms
    h.Scale(1., 'width')
    n_jets_inclusive = h.Integral(1, h.GetNbinsX(), 'width')
    h.Scale(1./n_jets_inclusive)
        
    return h

  #----------------------------------------------------------------------
  def load_systematic_percentage(self, systematic, jetR, obs_setting, reg_param, min_pt, max_pt, half=None):

    if half == 1:
      input_dir = self.input_dir1
    elif half == 2:
      input_dir = self.input_dir2

    # Load 2D distribution from file
    input_dir = os.path.join(input_dir, 'systematics')
    filename = os.path.join(input_dir, 'fSystematics.root')
    f = ROOT.TFile(filename, 'READ')
    hname = f'hSystematic_{self.observable}_{systematic}_R{jetR}_{obs_setting}_n{reg_param}_{min_pt}-{max_pt}'
    h = f.Get(hname)
    h.SetDirectory(0)
        
    return h

  #----------------------------------------------------------------------
  # Get systematic variation and return histogram of percentage difference
  def compute_systematic_variation(self, h_ratio_main, systematic):

    # First, set some systematic-specific info: hname and reg param

    # Set special cases of histogram name
    if 'reg_param' in systematic:
      systematic_name = 'main' 
    elif 'generator' in systematic:
      systematic_name = 'fastsim_generator1'  
    else:
      systematic_name = systematic

    # Set special cases of reg param
    if 'reg_param' in systematic:
      if systematic == 'reg_param+':
        reg_param_numerator = self.reg_param_numerator + 2
        reg_param_denominator = self.reg_param_denominator + 2
      elif systematic == 'reg_param-':
        reg_param_numerator = self.reg_param_numerator - 2
        reg_param_denominator = self.reg_param_denominator - 2   
    elif 'generator' in systematic: # For subjet JEWEL systematic, set a larger reg param so that it can converge
        reg_param_numerator = self.reg_param_generator
        reg_param_denominator = self.reg_param_generator
    else:
      reg_param_numerator = self.reg_param_numerator
      reg_param_denominator = self.reg_param_numerator

    # Load unfolded variations and form ratio
    h_numerator_systematic = self.load_observable(systematic_name, self.jetR, self.obs_settings[0], reg_param_numerator, half=self.numerator_half)
    h_denominator_systematic = self.load_observable(systematic_name, self.jetR, self.obs_settings[1], reg_param_denominator, half=self.denominator_half)

    h_ratio_systematic = h_numerator_systematic.Clone()
    h_ratio_systematic.SetName(f'h_ratio_{systematic}')
    h_ratio_systematic.Divide(h_denominator_systematic)

    # Divide by main result ratio
    h_ratio_systematic.Divide(h_ratio_main)

    # Convert to percentage uncertainty
    self.change_to_per(h_ratio_systematic, signed=True)

    return h_ratio_systematic

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
  def hist_average(self, h1, h2, signed=False, take_max_dev=False):

    h_avg = h1.Clone()
    h_avg.SetName('{}_avg'.format(h1.GetName()))

    for i in range(1, h_avg.GetNbinsX()+1):
      if signed:
        value1 = h1.GetBinContent(i)
        value2 = h2.GetBinContent(i)
      else:
        value1 = np.abs(h1.GetBinContent(i))
        value2 = np.abs(h2.GetBinContent(i))
      avg =  0.5*(value1 + value2)

      if take_max_dev:
        if value1 > value2:
          avg = value1
        else:  # value2 > value1:
          avg = value2

      h_avg.SetBinContent(i, avg)

    return h_avg

  #----------------------------------------------------------------------
  # Take stdev of a list of (identically-binned) histograms, bin-by-bin
  #----------------------------------------------------------------------
  def hist_stdev(self, h_list):

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
  # Add a list of (identically-binned) histograms in quadrature, bin-by-bin
  #----------------------------------------------------------------------
  def add_in_quadrature(self, h_dict):

    h_new = next(iter(h_dict.values())).Clone()
    h_new.SetName('total_systematic')

    for i in range(1, h_new.GetNbinsX()+1):

      values_i = [h.GetBinContent(i) for key,h in h_dict.items()]

      new_value_squared = 0.
      for value_i in values_i:
        new_value_squared += value_i*value_i
      new_value = math.sqrt(new_value_squared)
      h_new.SetBinContent(i, new_value)

    return h_new

  #----------------------------------------------------------------------
  # Add a list of (identically-binned) histograms in quadrature, bin-by-bin
  #----------------------------------------------------------------------
  def add_hlist_in_quadrature(self, h_list, new_name=None):

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
  def attach_uncertainty_to_hist(self, h, hPercError):

    #Fill array with lower bin edges of data histogram
    for bin in range(1, h.GetNbinsX()+1):
      content = h.GetBinContent(bin)
      perErr = hPercError.GetBinContent(bin)
      h.SetBinError(bin, content*perErr*0.01)

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
  def plot_systematic_uncertainties(self, h_list, h_total, suffix=''):

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    name = f'c{suffix}'
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

    truth_bin_array = self.truth_obs_bin_array[(self.truth_obs_bin_array > self.min_zr-1.e-3)]
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
    text = str(self.pt_bin[0]) + ' < #it{p}_{T, ch jet} < ' + str(self.pt_bin[1])
    text_latex.DrawLatex(0.3, 0.85, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = '#it{R} = ' + str(self.jetR)
    text_latex.DrawLatex(0.3, 0.78, text)

    text = f'({self.obs_setting_label}={self.obs_settings[0]})/({self.obs_setting_label}={self.obs_settings[1]})'
    text_latex.DrawLatex(0.3, 0.71, text)

    outputFilename = os.path.join(self.output_dir, f'hSystematics_ratio_R{self.utils.remove_periods(self.jetR)}_{int(self.pt_bin[0])}-{int(self.pt_bin[1])}{suffix}{self.file_format}')
    c.SaveAs(outputFilename)
    c.Close()

    if not suffix:
        sys_root_filename = os.path.join(self.output_dir, 'fSystematics.root')
        fSystematics = ROOT.TFile(sys_root_filename, 'UPDATE')
        h_total.Write()
        for h in h_list:
            h.Write()
        fSystematics.Close()

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
    self.hepdata_submission = hepdata_lib.Submission()
  
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
      
    # Write submission files
    self.hepdata_submission.create_files(self.hepdata_dir)

  #----------------------------------------------------------------------
  def add_hepdata_table(self, i, jetR, obs_label, obs_setting, grooming_setting, min_pt, max_pt):
  
    index = self.set_hepdata_table_index(i, jetR, min_pt)
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
    self.hepdata_submission.add_table(table)

  #----------------------------------------------------------------------
  def set_hepdata_table_index(self, i, jetR, min_pt):
  
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
                
    return index

  #----------------------------------------------------------------------
  def set_hepdata_table_descriptors(self, table, jetR, obs_label, obs_setting, grooming_setting, min_pt, max_pt):
        
    table.description = r'Groomed jet radius (scaled) $\theta_{{\mathrm{g}}}$'
    x_label = r'$\theta_{{\mathrm{g}}}$'
    y_label = r'$\frac{1}{\sigma_{inc}} \frac{d\sigma}{d\theta_{{\mathrm{g}}}}$'
    table.location = 'Figure 4'
        
    table.description += '\n'
    table.description += r'${}<p_{{\mathrm{{T}}}}^{{\mathrm{{ch jet}}}}<{}$, Soft Drop $z_{{\mathrm{{cut}}}}=0.2, \beta=0$.'.format(min_pt, max_pt)
    table.description += '\n\nNote: The first bin corresponds to the Soft Drop untagged fraction.'
    
    table.description += '\n\n'
    table.description += r'For the "trkeff" systematic uncertainty sources, the signed systematic uncertainty breakdowns ($\pm$ vs. $\mp$), denote correlation across bins (both within this table, and across tables for a given centrality). For the remaining sources ("unfolding", "subtraction", "thermal_closure") no correlation information is specified ($\pm$ is always used).'
    
    return x_label, y_label

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

  print()
  print('Configuring...')
  print('configFile: \'{0}\''.format(args.configFile))
  print()

  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = ComputeRatioSystematics(config_file = args.configFile)
  analysis.compute_ratio_systematics()
