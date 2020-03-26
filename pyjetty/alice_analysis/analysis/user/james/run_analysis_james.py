#! /usr/bin/env python

import sys
import os
import argparse
from array import *
import numpy
import ROOT
import yaml

from pyjetty.alice_analysis.analysis.user.substructure import run_analysis

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class RunAnalysisJames(run_analysis.RunAnalysis):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(RunAnalysisJames, self).__init__(config_file, **kwargs)
    
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

    # Theory comparisons
    if 'fPythia' in config:
      self.fPythia_name = config['fPythia']
    if 'fNLL' in config:
      self.fNLL = config['fNLL']
    if 'fNPcorrection_numerator' in config and 'fNPcorrection_denominator' in config:
      self.fNPcorrection_numerator = config['fNPcorrection_numerator']
      self.fNPcorrection_denominator = config['fNPcorrection_denominator']

  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #---------------------------------------------------------------
  def plot_single_result(self, jetR, obs_label, obs_setting, sd_setting):
  
    self.plot_final_result(jetR, obs_label, obs_setting, sd_setting)
    
    if self.observable == 'theta_g' or self.observable == 'zg':
        self.get_nll_tgraph(jetR, obs_label, obs_setting, sd_setting, 20., 40.)
        self.get_nll_tgraph(jetR, obs_label, obs_setting, sd_setting, 40., 60.)
        self.get_nll_tgraph(jetR, obs_label, obs_setting, sd_setting, 60., 80.)

  
  #---------------------------------------------------------------
  # This function is called once after all subconfigurations have been looped over, for each R
  #---------------------------------------------------------------
  def plot_all_results(self, jetR, obs_label, obs_setting, sd_setting):
  
    if self.observable == 'theta_g' or self.observable == 'zg':
      self.plot_final_result_overlay(jetR)

      #self.plot_NPcorrection(jetR)
  
  #----------------------------------------------------------------------
  # This function is called once after all subconfigurations and jetR have been looped over
  #----------------------------------------------------------------------
  def plot_performance(self):
  
    print('Plotting performance plots...')
  
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
      
      self.plot_observable(jetR, obs_label, obs_setting, sd_setting, min_pt_truth, max_pt_truth, plot_pythia=False)
      #self.plot_observable(jetR, obs_label, obs_setting, sd_setting, min_pt_truth, max_pt_truth, plot_pythia=True)

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
  def plot_observable(self, jetR, obs_label, obs_setting, sd_setting, min_pt_truth, max_pt_truth, plot_pythia=False):
    
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
        hPythia = self.get_pythia_from_response(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
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
        hname = 'histogram_h_{}_B{}_{}-{}'.format(self.observable, obs_setting[1], int(min_pt_truth), int(max_pt_truth))
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
    h_sys = getattr(self, 'hResult_{}_systotal_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth))
    h_sys.SetLineColor(0)
    h_sys.SetFillColor(color)
    h_sys.SetFillColorAlpha(color, 0.3)
    h_sys.SetFillStyle(1001)
    h_sys.SetLineWidth(0)
    h_sys.DrawCopy('E2 same')
    
    name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    if sd_setting:
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
    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.73, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = '#it{R} = ' + str(jetR) + '   | #eta_{jet}| < 0.5'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.66, text)
    
    subobs_label = self.utils.formatted_subobs_label(self.observable)
    if subobs_label:
      text = '{} = {}'.format(subobs_label, obs_setting)
      text_latex.DrawLatex(0.57, 0.59, text)
    
    if sd_setting:
      text = self.utils.formatted_sd_label(sd_setting)
      text_latex.DrawLatex(0.57, 0.52, text)
      
      text_latex = ROOT.TLatex()
      text_latex.SetNDC()
      text_latex.SetTextSize(0.04)
      text = '#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged
      text_latex.DrawLatex(0.57, 0.45, text)
    
      if plot_pythia:
        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        text_latex.SetTextSize(0.04)
        text = ('#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged) + (', #it{f}_{tagged}^{pythia} = %3.3f' % fraction_tagged_pythia)
        text_latex.DrawLatex(0.57, 0.38, text)

    myLegend = ROOT.TLegend(0.27,0.7,0.5,0.85)
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

  analysis = RunAnalysisJames(config_file = args.configFile)
  analysis.run_analysis()
