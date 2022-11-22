#! /usr/bin/env python

import sys
import os
import argparse
from array import *
import numpy as np
import ROOT
import yaml

from pyjetty.alice_analysis.analysis.user.substructure import run_analysis

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class RunAnalysisJamesBase(run_analysis.RunAnalysis):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(RunAnalysisJamesBase, self).__init__(config_file, **kwargs)
    
    # Initialize yaml config
    self.initialize_user_config()
  
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
      self.max_distance = config['constituent_subtractor']['max_distance']

    # Theory comparisons
    if 'fPythia' in config:
      self.fPythia_name = config['fPythia']
      
    if 'fNLL' in config:
      self.fNLL = config['fNLL']
    else:
      self.fNLL = False

    if 'fNPcorrection_numerator' in config and 'fNPcorrection_denominator' in config:
      self.NPcorrection = True
      self.fNPcorrection_numerator = config['fNPcorrection_numerator']
      self.fNPcorrection_denominator = config['fNPcorrection_denominator']
    else:
      self.NPcorrection = False
      
  #----------------------------------------------------------------------
  def plot_observable(self, jetR, obs_label, obs_setting, grooming_setting, min_pt_truth, max_pt_truth, plot_pythia=False):
    
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
    name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    h = getattr(self, name)
    h.SetName(name)
    h.SetMarkerSize(1.5)
    h.SetMarkerStyle(20)
    h.SetMarkerColor(color)
    h.SetLineStyle(1)
    h.SetLineWidth(2)
    h.SetLineColor(color)
    
    name = 'hResult_{}_systotal_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    h_sys = getattr(self, name)
    h_sys.SetName(name)
    h_sys.SetLineColor(0)
    h_sys.SetFillColor(color)
    h_sys.SetFillColorAlpha(color, 0.3)
    h_sys.SetFillStyle(1001)
    h_sys.SetLineWidth(0)
    
    n_obs_bins_truth = self.n_bins_truth(obs_label)
    truth_bin_array = self.truth_bin_array(obs_label)
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_obs_bins_truth, truth_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle(xtitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle(ytitle)
    myBlankHisto.SetMaximum(3*h.GetMaximum())
    if 'subjet_z' in self.observable or self.observable == 'jet_axis':
      myBlankHisto.SetMaximum(1.7*h.GetMaximum())
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")

    if plot_pythia:
    
      hPythia, fraction_tagged_pythia = self.pythia_prediction(jetR, obs_setting, grooming_setting, obs_label, min_pt_truth, max_pt_truth)
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

    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.DrawLatex(0.57, 0.73, text)

    text = '#it{R} = ' + str(jetR) + '   | #eta_{jet}| < 0.5'
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
      
      if 'sd' in grooming_setting:
        fraction_tagged = getattr(self, 'tagging_fraction_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt_truth, max_pt_truth))
        text_latex.SetTextSize(0.04)
        text = '#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged
        text_latex.DrawLatex(0.57, 0.52-delta, text)
    
        if plot_pythia:
          text_latex.SetTextSize(0.04)
          text = ('#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged) + (', #it{f}_{tagged}^{pythia} = %3.3f' % fraction_tagged_pythia)
          text_latex.DrawLatex(0.57, 0.52-delta, text)

    myLegend = ROOT.TLegend(0.25,0.7,0.45,0.85)
    self.utils.setup_legend(myLegend,0.035)
    myLegend.AddEntry(h, 'ALICE pp', 'pe')
    myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    if plot_pythia:
      myLegend.AddEntry(hPythia, 'PYTHIA8 Monash2013', 'pe')
    myLegend.Draw()

    name = 'hUnfolded_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), obs_label, int(min_pt_truth), int(max_pt_truth), self.file_format)
    if plot_pythia:
      name = 'hUnfolded_R{}_{}_{}-{}_Pythia{}'.format(self.utils.remove_periods(jetR), obs_label, int(min_pt_truth), int(max_pt_truth), self.file_format)
    output_dir = getattr(self, 'output_dir_final_results')
    output_dir_single = output_dir + '/single_results'
    if not os.path.exists(output_dir_single):
      os.mkdir(output_dir_single)
    outputFilename = os.path.join(output_dir_single, name)
    c.SaveAs(outputFilename)
    c.Close()

    # Write result to ROOT file
    #final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
    #fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
    #h.Write()
    #h_sys.Write()
    #hPythia.Write()
    #fFinalResults.Close()

  #----------------------------------------------------------------------
  def pythia_prediction(self, jetR, obs_setting, grooming_setting, obs_label, min_pt_truth, max_pt_truth):
  
    plot_pythia_from_response = True
    plot_pythia_from_mateusz = False
    
    if plot_pythia_from_response:
    
      hPythia = self.get_pythia_from_response(jetR, obs_label, min_pt_truth, max_pt_truth)
      
      if grooming_setting and 'sd' in grooming_setting:
      
        # If SD, the untagged jets are in the first bin
        n_jets_inclusive = hPythia.Integral(1, hPythia.GetNbinsX())
        n_jets_tagged = hPythia.Integral(hPythia.FindBin(self.truth_bin_array(obs_label)[0]), hPythia.GetNbinsX()+1)
        
      else:
        n_jets_inclusive = hPythia.Integral(1, hPythia.GetNbinsX())
        n_jets_tagged = hPythia.Integral(hPythia.FindBin(self.truth_bin_array(obs_label)[0]), hPythia.GetNbinsX())

    elif plot_pythia_from_mateusz:
    
      fPythia_name = '/Users/jamesmulligan/Analysis_theta_g/Pythia_new/pythia.root'
      fPythia = ROOT.TFile(fPythia_name, 'READ')
      print(fPythia.ls())
      hname = 'histogram_h_{}_B{}_{}-{}'.format(self.observable, obs_label, int(min_pt_truth), int(max_pt_truth))
      hPythia = fPythia.Get(hname)
      n_jets_inclusive = hPythia.Integral(0, hPythia.GetNbinsX())
      n_jets_tagged = hPythia.Integral(hPythia2.FindBin(self.truth_bin_array(obs_label)[0]), hPythia2.GetNbinsX())
      
    fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
    hPythia.Scale(1./n_jets_inclusive, 'width')
      
    return [hPythia, fraction_tagged_pythia]

  #----------------------------------------------------------------------
  def get_pythia_from_response(self, jetR, obs_label, min_pt_truth, max_pt_truth):
  
    output_dir = getattr(self, 'output_dir_main')
    file = os.path.join(output_dir, 'response.root')
    f = ROOT.TFile(file, 'READ')

    thn_name = 'hResponse_JetPt_{}_R{}_{}_rebinned'.format(self.observable, jetR, obs_label)
    thn = f.Get(thn_name)
    thn.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)

    h = thn.Projection(3)
    h.SetName('hPythia_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth))
    h.SetDirectory(0)
    
    for i in range(1, h.GetNbinsX() + 1):
      h.SetBinError(i, 0)

    return h

  #----------------------------------------------------------------------
  def plot_final_result_overlay(self, i_config, jetR, overlay_list):
    print('Plotting overlay of {}'.format(overlay_list))

    # Plot overlay of different subconfigs, for fixed pt bin
    for bin in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[bin]
      max_pt_truth = self.pt_bins_reported[bin+1]

      # Plot PYTHIA
      self.plot_observable_overlay_subconfigs(i_config, jetR, overlay_list, min_pt_truth, max_pt_truth, plot_pythia=True, plot_ratio = True)

      # Plot NLL
      if self.fNLL:
        self.plot_observable_overlay_subconfigs(i_config, jetR, overlay_list, min_pt_truth, max_pt_truth, plot_nll = True, plot_ratio = False)

  #----------------------------------------------------------------------
  def plot_observable_overlay_subconfigs(self, i_config, jetR, overlay_list, min_pt_truth, max_pt_truth, plot_pythia=False, plot_nll=False, plot_ratio=False):
    
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
    pad1.SetTicks(1,1)
    pad1.Draw()
    pad1.cd()

    if self.observable == 'leading_subjet_z':
      myLegend = ROOT.TLegend(0.45,0.3,0.61,0.54)
    elif self.observable == 'inclusive_subjet_z':
      myLegend = ROOT.TLegend(0.38,0.32,0.54,0.56)
    else:
      myLegend = ROOT.TLegend(0.45,0.33,0.61,0.57)
    self.utils.setup_legend(myLegend,0.05, sep=-0.2)

    if self.observable == 'leading_subjet_z':
      myLegend.SetHeader('Leading anti-#it{k}_{T} subjets')
    elif self.observable == 'inclusive_subjet_z':
      myLegend.SetHeader('Inclusive anti-#it{k}_{T} subjets')
    
    name = 'hmain_{}_R{}_{{}}_{}-{}'.format(self.observable, jetR, min_pt_truth, max_pt_truth)
    ymax = 2*self.get_maximum(name, overlay_list)
    ymin = 2e-4
    ymin_ratio = 0.
    ymax_ratio = 1.99
    if self.observable == 'theta_g':
      ymin_ratio = 0.5
      ymax_ratio = 1.69
      ymax*=1.1
    if self.observable == 'leading_subjet_z':
        ymax = 16.99
        ymin = 1e-3
        ymin_ratio = 0.5
        ymax_ratio = 1.79
    if self.observable == 'inclusive_subjet_z':
        ymax = 2e4
        ymin = 2e-1
        pad1.SetLogy()
        ymin_ratio = 0.5
        ymax_ratio = 1.79
    
    # Get xmin and xmax over all hists
    xmin = 1
    xmax = 0
    for i, subconfig_name in enumerate(self.obs_subconfig_list):
      if subconfig_name not in overlay_list:
        continue
      xmin_temp = self.obs_config_dict[subconfig_name]['obs_bins_truth'][0]
      xmax_temp = self.obs_config_dict[subconfig_name]['obs_bins_truth'][-1]
      if xmin_temp < xmin:
        xmin = xmin_temp
      if xmax_temp > xmax:
        xmax = xmax_temp
      
    # Loop through hists
    for i, subconfig_name in enumerate(self.obs_subconfig_list):
    
      if subconfig_name not in overlay_list:
        continue

      obs_setting = self.obs_settings[i]
      grooming_setting = self.grooming_settings[i]
      obs_label = self.utils.obs_label(obs_setting, grooming_setting)
      
      if subconfig_name == overlay_list[0]:
        marker = 20
        marker_pythia = marker+4
        color = 600-6
      if subconfig_name == overlay_list[1]:
        marker = 21
        marker_pythia = marker+4
        color = 632-4
      if i > 1 and subconfig_name == overlay_list[2]:
        marker = 33
        marker_pythia = 27
        color = 416-2
            
      name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
      h = getattr(self, name)
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(marker)
      h.SetMarkerColor(color)
      h.SetLineStyle(1)
      h.SetLineWidth(2)
      h.SetLineColor(color)
      h.GetXaxis().SetRangeUser(self.obs_config_dict[subconfig_name]['obs_bins_truth'][0], self.obs_config_dict[subconfig_name]['obs_bins_truth'][-1])
      
      h_sys = getattr(self, 'hResult_{}_systotal_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth))
      h_sys.SetLineColor(0)
      h_sys.SetFillColor(color)
      h_sys.SetFillColorAlpha(color, 0.3)
      h_sys.SetFillStyle(1001)
      h_sys.SetLineWidth(0)
      h_sys.GetXaxis().SetRangeUser(self.obs_config_dict[subconfig_name]['obs_bins_truth'][0], self.obs_config_dict[subconfig_name]['obs_bins_truth'][-1])
      
      if subconfig_name == overlay_list[0]:

        pad1.cd()
        xtitle = getattr(self, 'xtitle')
        ytitle = getattr(self, 'ytitle')
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, xmin, xmax)
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.GetXaxis().SetTitleSize(0.085)
        myBlankHisto.SetXTitle(xtitle)
        myBlankHisto.GetYaxis().SetTitleOffset(1.5)
        myBlankHisto.SetYTitle(ytitle)
        myBlankHisto.SetMaximum(ymax)
        myBlankHisto.SetMinimum(ymin)
        if plot_ratio:
          myBlankHisto.SetMinimum(ymin) # Don't draw 0 on top panel
          myBlankHisto.GetYaxis().SetTitleSize(0.075)
          myBlankHisto.GetYaxis().SetTitleOffset(1.2)
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
          pad2.SetTicks(0,1)
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
          myBlankHisto2.GetXaxis().SetTickSize(0.07)
          myBlankHisto2.GetYaxis().SetTitleSize(25)
          myBlankHisto2.GetYaxis().SetTitleFont(43)
          myBlankHisto2.GetYaxis().SetTitleOffset(2.2)
          myBlankHisto2.GetYaxis().SetLabelFont(43)
          myBlankHisto2.GetYaxis().SetLabelSize(25)
          myBlankHisto2.GetYaxis().SetNdivisions(505)
          myBlankHisto2.GetYaxis().SetTickSize(0.035)
          
          myBlankHisto2.GetYaxis().SetRangeUser(ymin_ratio, ymax_ratio)
          myBlankHisto2.Draw()
        
          line = ROOT.TLine(xmin,1,xmax,1)
          line.SetLineColor(920+2)
          line.SetLineStyle(2)
          line.Draw()
      
      if plot_pythia:
      
        hPythia, fraction_tagged_pythia = self.pythia_prediction(jetR, obs_setting, grooming_setting, obs_label, min_pt_truth, max_pt_truth)

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
          
        hPythia.GetXaxis().SetRangeUser(self.obs_config_dict[subconfig_name]['obs_bins_truth'][0], self.obs_config_dict[subconfig_name]['obs_bins_truth'][-1])

      if plot_nll:
        
        # Get parton-level prediction
        attr_name = 'tgraph_NLL_{}_{}_{}-{}'.format(self.observable, obs_label, min_pt_truth, max_pt_truth)
        if hasattr(self, attr_name):
          g = getattr(self, attr_name)
        else:
          return
        
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

      # Scale inclusive subjets to N_jets rather than N_subjets -- fill in by hand for now
      z_min = 0.71
      if self.observable == 'leading_subjet_z':
        
        integral_total = h.Integral(1, h.GetNbinsX(), 'width')
        print(f'integral_total, leading_subjet_z, {obs_label}: {integral_total}')
        h.Scale(1./integral_total)
        h_sys.Scale(1./integral_total)
        
        integral_total_pythia = hPythia.Integral(1, hPythia.GetNbinsX(), 'width')
        print(f'integral_total_pythia, leading_subjet_z, {obs_label}: {integral_total_pythia}')
        hPythia.Scale(1./integral_total_pythia)
        
        integral = h.Integral(h.FindBin(z_min), h.GetNbinsX(), 'width')
        print(f'integral, leading_subjet_z, {obs_label}: {integral}')
        
        integral_pythia = hPythia.Integral(hPythia.FindBin(z_min), hPythia.GetNbinsX(), 'width')
        print(f'integral_pythia, leading_subjet_z, {obs_label}: {integral_pythia}')
        
      if self.observable == 'inclusive_subjet_z':
        
        if np.isclose(float(obs_label), 0.1):
            integral_leading_subjet = 0.7865922387180659 # R=0.4, r=0.1
            integral_leading_subjet_pythia = 0.7441921938539977
        elif np.isclose(float(obs_label), 0.2):
            integral_leading_subjet = 0.9365660517984986 # R=0.4, r=0.2
            integral_leading_subjet_pythia = 0.9296991675820692

        integral = h.Integral(h.FindBin(z_min), h.GetNbinsX(), 'width')
        print(f'integral, inclusive_subjet_z, {obs_label}: {integral}')
        normalization = integral_leading_subjet / integral
        print(f'normalization: {normalization}')
        h.Scale(normalization)
        h_sys.Scale(normalization)

        integral_pythia = hPythia.Integral(hPythia.FindBin(z_min), hPythia.GetNbinsX(), 'width')
        print(f'integral_pythia, inclusive_subjet_z, {obs_label}: {integral_pythia}')
        normalization_pythia = integral_leading_subjet_pythia / integral_pythia
        print(f'normalization_pythia: {normalization_pythia}')
        hPythia.Scale(normalization_pythia)
        
      # Compute <N_subjets>
      n_subjets = h.Integral(1, h.GetNbinsX(), 'width')
      print(f'<N_subjets>, {obs_label}: {n_subjets}')
        
      # Compute z_loss for leading subjets
      # Should come up with a better way to decide bin center
      z_moment = 0.
      if self.observable == 'leading_subjet_z':
        for i in range(1, h.GetNbinsX()+1):
            zr = h.GetBinCenter(i)
            content = h.GetBinContent(i)
            width =  h.GetXaxis().GetBinWidth(i)
            z_moment += zr*content*width
            #print(f'bin: {i} (zr = {zr}, width = {width}): content = {content} -- {zr*content*width}')
        z_loss = 1 - z_moment
        #print(z_moment)
        #print(f'z_loss for r={obs_label}: {1-z_moment}')

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
        elif plot_nll:
          gRatioSys = g.Clone()
          gRatioSys.SetName('{}_{}_Ratio'.format(obs_label, g.GetName()))
          self.utils.divide_tgraph(hRatioSys, gRatioSys, combine_errors=True)
          gRatioSys.SetLineColor(0)
          gRatioSys.SetFillColor(color)
          gRatioSys.SetFillColorAlpha(color, 0.3)
          gRatioSys.SetFillStyle(1001)
          gRatioSys.SetLineWidth(0)
          
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

      pad1.cd()

      if plot_pythia:
        plot_errors = False
        if plot_errors:
          hPythia.DrawCopy('E3 same')
        else:
          hPythia.DrawCopy('L hist same')
          
      if plot_nll:
        g.Draw('L3 same')

      h_sys.DrawCopy('E2 same')
      h.DrawCopy('PE X0 same')
      
      if plot_ratio:
        pad2.cd()
        if plot_pythia:
          hRatioSys.DrawCopy('E2 same')
        elif plot_nll:
          gRatioSys.Draw('L3 same')
        hRatioStat.DrawCopy('PE X0 same')
        
      subobs_label = self.utils.formatted_subobs_label(self.observable)
      text = ''
      if subobs_label:
        text += '{} = {}'.format(subobs_label, obs_setting)
      if grooming_setting:
        text += self.utils.formatted_grooming_label(grooming_setting, verbose=True)
      myLegend.AddEntry(h, '{}'.format(text), 'pe')
      
      if self.observable == 'leading_subjet_z':
          pad1.cd()
          text_latex = ROOT.TLatex()
          text_latex.SetNDC()
          text_latex.SetTextSize(0.05)
          x = 0.3
          y = 0.3 - 0.9*float(obs_setting)
          text = f'#LT#it{{z}}^{{loss}}_{{{subobs_label} = {obs_setting} }}#GT = {np.round(z_loss,2):.2f}'
          text_latex.DrawLatex(x, y, text)
      elif self.observable == 'inclusive_subjet_z':
          pad1.cd()
          text_latex = ROOT.TLatex()
          text_latex.SetNDC()
          text_latex.SetTextSize(0.045)
          x = 0.47
          y = 0.33 - 0.8*float(obs_setting)
          text = f'#LT#it{{N}}^{{subjets }}_{{{subobs_label} = {obs_setting}}}#GT = {np.round(n_subjets,2):.1f}'
          text_latex.DrawLatex(x, y, text)

      # Write result to ROOT file
      final_result_root_filename = os.path.join(getattr(self, 'output_dir_final_results'), 'fFinalResults.root')
      fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
      h.Write()
      h_sys.Write()
      hPythia.Write()
      fFinalResults.Close()
        
    pad1.cd()
    myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    if plot_pythia:
      myLegend.AddEntry(hPythia, 'PYTHIA8 Monash 2013', 'l')
    if plot_nll:
      myLegend.AddEntry(g, 'NLL', 'l')
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    
    text_latex.SetTextSize(0.065)
    x = 0.25
    y = 0.855
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(x, y, text)
    
    text_latex.SetTextSize(0.055)
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.DrawLatex(x, y-0.06, text)

    text = 'Charged-particle anti-#it{k}_{T} jets'
    text_latex.DrawLatex(x, y-0.12, text)
    
    text = '#it{R} = ' + str(jetR) + '   | #it{{#eta}}_{{jet}}| < {}'.format(0.9-jetR)
    text_latex.DrawLatex(x, y-0.18, text)
    
    text = str(min_pt_truth) + ' < #it{p}_{T}^{ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.DrawLatex(x, y-0.26, text)
    
    myLegend.Draw()
    
    if self.observable == 'theta_g':
      rg_axis_tf1 = ROOT.TF1('rg_axis_tf1', 'x', 0, jetR-0.01)
      rg_axis = ROOT.TGaxis(xmin, 2*ymax, xmax, 2*ymax, 'rg_axis_tf1', 505, '- S')
      rg_axis.SetTitle('#it{R}_{g}')
      rg_axis.SetTitleSize(25)
      rg_axis.SetTitleFont(43)
      rg_axis.SetTitleOffset(0.6)
      rg_axis.SetLabelFont(43)
      rg_axis.SetLabelSize(25)
      rg_axis.SetTickSize(0.015)
      rg_axis.SetLabelOffset(0.015)
      rg_axis.Draw()

    name = 'h_{}_R{}_{}-{}_{}{}'.format(self.observable, self.utils.remove_periods(jetR), int(min_pt_truth), int(max_pt_truth), i_config, self.file_format)
    if plot_pythia:
      name = 'h_{}_R{}_{}-{}_Pythia_{}{}'.format(self.observable, self.utils.remove_periods(jetR), int(min_pt_truth), int(max_pt_truth), i_config, self.file_format)
    if plot_nll:
      name = 'h_{}_R{}_{}-{}_NLL_{}{}'.format(self.observable, self.utils.remove_periods(jetR), int(min_pt_truth), int(max_pt_truth), i_config, self.file_format)
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
