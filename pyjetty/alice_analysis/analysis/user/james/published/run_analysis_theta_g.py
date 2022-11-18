#! /usr/bin/env python

import sys
import os
import argparse
from array import *
import numpy as np
import ROOT
import yaml

from pyjetty.alice_analysis.analysis.user.james import run_analysis_james_base
from pyjetty.alice_analysis.analysis.user.james import plotting_utils_theta_g

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class RunAnalysisThetaG(run_analysis_james_base.RunAnalysisJamesBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(RunAnalysisThetaG, self).__init__(config_file, **kwargs)
    
    print(self)

  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #---------------------------------------------------------------
  def plot_single_result(self, jetR, obs_label, obs_setting, grooming_setting):
    print('Plotting each individual result...')
  
    # Plot final result for each 1D substructure distribution (with PYTHIA)
    self.plot_final_result(jetR, obs_label, obs_setting, grooming_setting)
    
    if self.observable == 'theta_g' or self.observable == 'zg':
    
      if 'sd' in grooming_setting:

        # Construct NP correction, and store as an attribute
        if self.NPcorrection:
          self.construct_NPcorrections(jetR, obs_label)
      
        # Construct TGraph of NLL predictions
        if self.fNLL:
          self.construct_nll_tgraphs(jetR, obs_label, obs_setting, grooming_setting)
  
  #---------------------------------------------------------------
  # This function is called once after all subconfigurations have been looped over, for each R
  #---------------------------------------------------------------
  def plot_all_results(self, jetR):
    print('Plotting overlay of all results...')
    
    for i_config, overlay_list in enumerate(self.plot_overlay_list):
    
      if len(overlay_list) > 1:
      
        self.plot_final_result_overlay(i_config, jetR, overlay_list)

        self.plot_NPcorrections(i_config, jetR, overlay_list)
  
  #----------------------------------------------------------------------
  # This function is called once after all subconfigurations and jetR have been looped over
  #----------------------------------------------------------------------
  def plot_performance(self):
    
    if not self.do_plot_performance:
      return
    print('Plotting performance plots...')
    
    # Initialize performance plotting class, and plot
    if self.is_pp:
    
      self.plotting_utils = plotting_utils_theta_g.PlottingUtils(self.output_dir_performance, self.config_file)
      self.plot_single_performance(self.output_dir_performance)
      
    else:
      
      # Plot for each R_max
      for R_max in self.max_distance:
      
        output_dir_performance = os.path.join(self.output_dir_performance, 'Rmax{}'.format(R_max))
        self.plotting_utils = plotting_utils_theta_g.PlottingUtils(output_dir_performance, self.config_file, R_max = R_max)
        self.plot_single_performance(output_dir_performance, R_max)

        # Plot for thermal model
        if self.do_thermal_closure and R_max == self.R_max:
          
          output_dir_performance = os.path.join(self.output_dir_performance, 'thermal')
          self.plotting_utils = plotting_utils_theta_g.PlottingUtils(output_dir_performance, self.config_file, R_max = R_max, thermal = True)
          self.plot_single_performance(output_dir_performance, R_max)

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
    self.create_output_subdir(output_dir_performance, 'lund')
    if not self.is_pp:
      self.create_output_subdir(output_dir_performance, 'delta_pt')
      self.create_output_subdir(output_dir_performance, 'prong_matching_fraction_pt')
      self.create_output_subdir(output_dir_performance, 'prong_matching_fraction_ptdet')
      self.create_output_subdir(output_dir_performance, 'prong_matching_deltaR')
      self.create_output_subdir(output_dir_performance, 'prong_matching_deltaZ')
      self.create_output_subdir(output_dir_performance, 'prong_matching_correlation')
    
    # Generate performance plots
    for jetR in self.jetR_list:
  
      # Plot some subobservable-independent performance plots
      self.plotting_utils.plot_DeltaR(jetR, self.jet_matching_distance)
      self.plotting_utils.plot_JES(jetR)
      self.plotting_utils.plot_JES_proj(jetR, self.pt_bins_reported)
      self.plotting_utils.plotJER(jetR, self.utils.obs_label(self.obs_settings[0], self.grooming_settings[0]))
      self.plotting_utils.plot_jet_reco_efficiency(jetR, self.utils.obs_label(self.obs_settings[0], self.grooming_settings[0]))
      
      if not self.is_pp:
        self.plotting_utils.plot_delta_pt(jetR, self.pt_bins_reported)
      
      # Plot subobservable-dependent performance plots
      for i, _ in enumerate(self.obs_subconfig_list):

        obs_setting = self.obs_settings[i]
        grooming_setting = self.grooming_settings[i]
        obs_label = self.utils.obs_label(obs_setting, grooming_setting)
    
        self.plotting_utils.plot_obs_resolution(jetR, obs_label, self.xtitle, self.pt_bins_reported)
        self.plotting_utils.plot_obs_residual_pt(jetR, obs_label, self.xtitle, self.pt_bins_reported)
        self.plotting_utils.plot_obs_residual_obs(jetR, obs_label, self.xtitle)
        self.plotting_utils.plot_obs_projections(jetR, obs_label, obs_setting, grooming_setting, self.xtitle, self.pt_bins_reported)
        self.plotting_utils.plot_obs_truth(jetR, obs_label, obs_setting, grooming_setting, self.xtitle, self.pt_bins_reported)
        
        if grooming_setting and self.observable != 'jet_axis':
          self.plotting_utils.plot_lund_plane(jetR, obs_label, grooming_setting)

      # Plot prong matching histograms
      if not self.is_pp:
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

              hname = 'hProngMatching_{}_{}_JetPtDet_R{}'.format(prong, match, jetR)
              self.plotting_utils.plot_prong_matching(i, jetR, hname, self.obs_subconfig_list, self.obs_settings, self.grooming_settings, overlay_list, self.prong_match_threshold)

              if 'subleading' in prong or 'leading' in prong:
                hname = 'hProngMatching_{}_{}_JetPtZ_R{}'.format(prong, match, jetR)
                self.plotting_utils.plot_prong_matching_delta(i, jetR, hname, self.obs_subconfig_list, self.obs_settings, self.grooming_settings, overlay_list, self.prong_match_threshold, min_pt, max_pt, plot_deltaz=True)

          hname = 'hProngMatching_subleading-leading_correlation_JetPtDet_R{}'.format(jetR)
          self.plotting_utils.plot_prong_matching_correlation(i, jetR, hname, self.obs_subconfig_list, self.obs_settings, self.grooming_settings, overlay_list, self.prong_match_threshold)

  #----------------------------------------------------------------------
  def plot_final_result(self, jetR, obs_label, obs_setting, grooming_setting):
    print('Plot final results for {}: R = {}, {}'.format(self.observable, jetR, obs_label))

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()
    
    # Construct histogram of tagging fraction, to write to file
    if 'sd' in grooming_setting:
      name = 'h_tagging_fraction_R{}_{}'.format(jetR, obs_label)
      h_tagging = ROOT.TH1D(name, name, len(self.pt_bins_reported) - 1, array('d', self.pt_bins_reported))

    # Loop through pt slices, and plot final result for each 1D theta_g distribution
    for bin in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[bin]
      max_pt_truth = self.pt_bins_reported[bin+1]
      
      self.plot_observable(jetR, obs_label, obs_setting, grooming_setting, min_pt_truth, max_pt_truth, plot_pythia=True)
      
      # Fill tagging fraction
      if 'sd' in grooming_setting:
        fraction_tagged = getattr(self, 'tagging_fraction_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt_truth, max_pt_truth))
        pt = (min_pt_truth + max_pt_truth)/2
        h_tagging.Fill(pt, fraction_tagged)
      
    # Write tagging fraction to ROOT file
    if 'sd' in grooming_setting:
      output_dir = getattr(self, 'output_dir_final_results')
      final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
      fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
      h_tagging.Write()
      fFinalResults.Close()
   
  #----------------------------------------------------------------------
  def construct_NPcorrections(self, jetR, obs_label):
    
    for bin in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[bin]
      max_pt_truth = self.pt_bins_reported[bin+1]
  
      self.get_NPcorrection(jetR, obs_label, min_pt_truth, max_pt_truth)

  #----------------------------------------------------------------------
  def get_NPcorrection(self, jetR, obs_label, min_pt_truth, max_pt_truth):
  
    fNumerator = ROOT.TFile(self.fNPcorrection_numerator, 'READ')
    fDenominator = ROOT.TFile(self.fNPcorrection_denominator, 'READ')
    hname = 'histogram_h_{}_B{}_{}-{}'.format(self.observable, obs_label[-1], int(min_pt_truth), int(max_pt_truth))

    hNumerator = fNumerator.Get(hname)
    if not hNumerator:
      print('NP prediction does not exist for {}'.format(obs_label))
      return
    hNumerator.SetDirectory(0)
    n_jets_inclusive = hNumerator.Integral(0, hNumerator.GetNbinsX())
    n_jets_tagged = hNumerator.Integral(hNumerator.FindBin(self.truth_bin_array(obs_label)[0]), hNumerator.GetNbinsX())
    fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
    hNumerator.Scale(1./n_jets_inclusive, 'width')
      
    hDenominator = fDenominator.Get(hname)
    hDenominator.SetDirectory(0)
    n_jets_inclusive = hDenominator.Integral(0, hDenominator.GetNbinsX())
    n_jets_tagged = hDenominator.Integral(hDenominator.FindBin(self.truth_bin_array(obs_label)[0]), hDenominator.GetNbinsX())
    fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
    hDenominator.Scale(1./n_jets_inclusive, 'width')
        
    hNumerator.Divide(hDenominator)
    hNumerator.SetName('hNPcorrection_{}_{}_{}-{}'.format(self.observable, obs_label, min_pt_truth, max_pt_truth))
    setattr(self, 'hNPcorrection_{}_{}_{}-{}'.format(self.observable, obs_label, min_pt_truth, max_pt_truth), hNumerator)

  #----------------------------------------------------------------------
  def construct_nll_tgraphs(self, jetR, obs_label, obs_setting, grooming_setting):
    
    for bin in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[bin]
      max_pt_truth = self.pt_bins_reported[bin+1]
      
      self.get_nll_tgraph(jetR, obs_label, obs_setting, grooming_setting, min_pt_truth, max_pt_truth)

  #----------------------------------------------------------------------
  def get_nll_tgraph(self, jetR, obs_label, obs_setting, grooming_setting, min_pt_truth, max_pt_truth):

    n_bins_truth = self.n_bins_truth(obs_label)
    truth_bin_array = self.truth_bin_array(obs_label)
    
    key, value = list(grooming_setting.items())[0]
    beta = value[1]
    if self.observable == 'theta_g':
      path_txt = '/Users/jamesmulligan/Analysis_theta_g/NLL/Rg_value/beta{}/{}_{}.dat'.format(beta, int(min_pt_truth), int(max_pt_truth))
    if self.observable == 'zg':
      path_txt = '/Users/jamesmulligan/Analysis_theta_g/NLL/zg_value/beta{}/{}_{}.dat'.format(beta, int(min_pt_truth), int(max_pt_truth))
    
    if not os.path.exists(path_txt):
      print('NLL prediction does not exist for {}'.format(obs_label))
      return
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
    x = np.array(x_list)
    x_err = np.zeros(n_bins_truth)
    center = np.array(center_list)
    low = np.subtract(center, np.array(low_list))
    up = np.subtract(np.array(up_list), center)
    
    g = ROOT.TGraphAsymmErrors(n_bins_truth, x, center, x_err, x_err, low, up)
    g.SetName('tgraph_{}_{}_{}-{}'.format(self.observable, obs_label, min_pt_truth, max_pt_truth))
    setattr(self, 'tgraph_NLL_{}_{}_{}-{}'.format(self.observable, obs_label, min_pt_truth, max_pt_truth), g)

  #----------------------------------------------------------------------
  def plot_NPcorrections(self, i_config, jetR, overlay_list):
  
    # Loop through pt slices, and plot NP correction for each 1D observable
    for bin in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[bin]
      max_pt_truth = self.pt_bins_reported[bin+1]
      
      self.plot_NPcorrection_overlay(i_config, jetR, overlay_list, min_pt_truth, max_pt_truth)

  #----------------------------------------------------------------------
  def plot_NPcorrection_overlay(self, i_config, jetR, overlay_list, min_pt_truth, max_pt_truth):
    
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
    
    pad1.cd()
    myLegend = ROOT.TLegend(0.66,0.65,0.8,0.85)
    self.utils.setup_legend(myLegend,0.035)
    
    # Retrieve histograms for each observable setting
    for i, subconfig_name in enumerate(self.obs_subconfig_list):
    
      if subconfig_name not in overlay_list:
        continue

      obs_setting = self.obs_settings[i]
      grooming_setting = self.grooming_settings[i]
      obs_label = self.utils.obs_label(obs_setting, grooming_setting)
      
      if subconfig_name == overlay_list[0]:
        marker = 20
        color = 600-6
      if subconfig_name == overlay_list[1]:
        marker = 21
        color = 632-4
      if i > 1 and subconfig_name == overlay_list[2]:
        marker = 33
        color = 416-2

      attr_name = 'hNPcorrection_{}_{}_{}-{}'.format(self.observable, obs_label, min_pt_truth, max_pt_truth)
      if hasattr(self, attr_name):
        h = getattr(self, attr_name)
      else:
        return
        
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(marker)
      h.SetMarkerColor(color)
      h.SetLineStyle(1)
      h.SetLineWidth(2)
      h.SetLineColor(color)
      
      if subconfig_name == overlay_list[0]:
      
        xmin = self.obs_config_dict[subconfig_name]['obs_bins_truth'][0]
        xmax = self.obs_config_dict[subconfig_name]['obs_bins_truth'][-1]
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, xmin, xmax)
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle( getattr(self, 'xtitle') )
        myBlankHisto.GetYaxis().SetTitleOffset(1.5)
        myBlankHisto.SetYTitle('Correction')
        myBlankHisto.SetMaximum(2*h.GetMaximum())
        myBlankHisto.SetMinimum(0.)
        myBlankHisto.Draw("E")
        
        line = ROOT.TLine(0,1,xmax,1)
        line.SetLineColor(920+2)
        line.SetLineStyle(2)
        line.Draw()
      
      h.DrawCopy('PE X0 same')
      
      subobs_label = self.utils.formatted_subobs_label(self.observable)
      text = ''
      if subobs_label:
        text += '{} = {}'.format(subobs_label, obs_setting)
      if grooming_setting:
        text += self.utils.formatted_grooming_label(grooming_setting)
      myLegend.AddEntry(h, '{}'.format(text), 'pe')

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(0.25, 0.87, text)
    
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.81, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'PYTHIA8 Monash2013'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.75, text)
    
    text = '#it{R} = ' + str(jetR) + '   | #eta_{jet}| < 0.5'
    text_latex.DrawLatex(0.25, 0.69, text)
    
    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.63, text)
    
    myLegend.Draw()
    
    name = 'hNPcorrection_{}_R{}_{}-{}_{}{}'.format(self.observable, self.utils.remove_periods(jetR), int(min_pt_truth), int(max_pt_truth), i_config, self.file_format)
    output_dir = getattr(self, 'output_dir_final_results') + '/NPcorrection'
    if not os.path.exists(output_dir):
      os.mkdir(output_dir)
    outputFilename = os.path.join(output_dir, name)
    c.SaveAs(outputFilename)
    c.Close()

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

  analysis = RunAnalysisThetaG(config_file = args.configFile)
  analysis.run_analysis()
