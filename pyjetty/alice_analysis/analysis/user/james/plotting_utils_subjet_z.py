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

# Data analysis and plotting
import numpy as np
from array import *
import ROOT

# Base class
from pyjetty.alice_analysis.analysis.user.james import plotting_utils_base

################################################################
class PlottingUtils(plotting_utils_base.PlottingUtilsBase):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, output_dir = '.', config_file = '', R_max = None, thermal = False, groomer_studies = False, **kwargs):
    super(PlottingUtils, self).__init__(output_dir, config_file, R_max, thermal, groomer_studies, **kwargs)
    
    self.xmin = 0.
    self.xmax = 1.
    self.rebin_value = 5
    
    print(self)

  #---------------------------------------------------------------
  def plot_subjet_matching(self, i_overlay, jetR, name_prefix, obs_subconfig_list, obs_settings, grooming_settings, overlay_list, prong_match_threshold, thermal=False):

    c = ROOT.TCanvas('c','c: hist',600,450)
    c.cd()
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.15)
    c.GetPad(0).SetTicks(0,1)
    
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 20, 0., 200.)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.GetXaxis().SetTitleOffset(1.4)
    myBlankHisto.GetYaxis().SetTitleOffset(1.4)
    myBlankHisto.GetXaxis().SetTitleSize(0.05)
    myBlankHisto.GetYaxis().SetTitleSize(0.055)
    if self.groomer_studies:
      myBlankHisto.SetMaximum(1.61)
    else:
      myBlankHisto.SetMaximum(2.01)
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.GetYaxis().SetTitle('Subjet purity')
    myBlankHisto.Draw()
    
    if self.groomer_studies:
      myLegend = ROOT.TLegend(0.42,0.64,0.55,0.83)
      self.setup_legend(myLegend,0.035)
      myLegend2 = ROOT.TLegend(0.76,0.64,0.88,0.83)
      self.setup_legend(myLegend2,0.035)
    else:
      myLegend = ROOT.TLegend(0.66,0.6,0.78,0.7)
      self.setup_legend(myLegend,0.043)
    
    h_list = [] # Store hists in a list, since otherwise it seems I lose the marker information
                # (removed from memory?)

    i_reset = 0
    for i, subconfig_name in enumerate(obs_subconfig_list):
    
      if subconfig_name not in overlay_list:
        continue

      obs_setting = obs_settings[i]
      grooming_setting = grooming_settings[i]
      obs_label = self.obs_label(obs_setting, grooming_setting)
      
      if (jetR - obs_setting) < 1e-3:
        continue
      
      name = '{}_{}{}{}'.format(name_prefix, obs_setting, self.suffix, self.scaled_suffix)
      hFraction_vs_pt = self.fMC.Get(name)
      xtitle = hFraction_vs_pt.GetXaxis().GetTitle()
      #myBlankHisto.GetXaxis().SetTitle(xtitle)
      if self.groomer_studies:
        myBlankHisto.GetXaxis().SetTitle('#it{p}_{T,ch jet}^{PYTHIA} (GeV/#it{c})')
      else:
        myBlankHisto.GetXaxis().SetTitle('#it{p}_{T,ch jet}^{pp-det} (GeV/#it{c})')
      
      epsilon = 1e-5
      min_bin = hFraction_vs_pt.GetZaxis().FindBin(0. + epsilon)
      cut_bin = hFraction_vs_pt.GetZaxis().FindBin(prong_match_threshold + epsilon)
      max_bin = hFraction_vs_pt.GetZaxis().FindBin(1. + epsilon)
      
      name_total = 'hTotal_vs_pt{}_{}_{}'.format(i, jetR, obs_label)
      hTotal_vs_pt = hFraction_vs_pt.ProjectionX(name_total, 0, -1, min_bin, max_bin)
      hTotal_vs_pt.SetName(name_total)
      
      name_matched = 'hMatched_vs_pt{}_{}_{}'.format(i, jetR, obs_label)
      hMatched_vs_pt = hFraction_vs_pt.ProjectionX(name_matched, 0, -1, cut_bin, max_bin)
      hMatched_vs_pt.SetName(name_matched)
      
      name_fraction = 'hMatchedFraction_vs_pt{}_{}_{}'.format(i, jetR, obs_label)
      hMatchedFraction_vs_pt = None
      hMatchedFraction_vs_pt = hMatched_vs_pt.Clone()
      hMatchedFraction_vs_pt.SetName(name_fraction)
      hMatchedFraction_vs_pt.Divide(hTotal_vs_pt)
      hMatchedFraction_vs_pt.SetMarkerStyle(self.MarkerArray[i_reset])
      hMatchedFraction_vs_pt.SetMarkerColor(self.ColorArray[i_reset])
      hMatchedFraction_vs_pt.SetLineColor(self.ColorArray[i_reset])
      i_reset += 1
      if self.thermal:
        hMatchedFraction_vs_pt.SetMarkerColor(self.ColorArray[len(self.ColorArray)-i-2])
        hMatchedFraction_vs_pt.SetLineColor(self.ColorArray[len(self.ColorArray)-i-2])
      hMatchedFraction_vs_pt.DrawCopy('P same')

      if self.groomer_studies:
        if i_reset < 7:
          myLegend.AddEntry(hMatchedFraction_vs_pt, '#it{{r}} = {}'.format(obs_setting), 'P')
      else:
        myLegend.AddEntry(hMatchedFraction_vs_pt, '#it{{r}} = {}'.format(obs_setting), 'P')
              
      h_list.append(hMatchedFraction_vs_pt)

    myLegend.Draw('same')
    if self.groomer_studies:
      myLegend2.Draw('same')

    line = ROOT.TLine(0, 1, 200, 1)
    line.SetLineColor(920+2)
    line.SetLineStyle(2)
    line.SetLineWidth(4)
    line.Draw()
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()

    if self.groomer_studies:
      x = 0.17
      y = 0.85
      text_latex.SetTextSize(0.04)
      text = 'PYTHIA8 Monash 2013 embedded in thermal background'
      text_latex.DrawLatex(x, y, text)
    else:
      x = 0.23
      y = 0.8
      text_latex.SetTextSize(0.055)
      text = self.figure_approval_status
      text_latex.DrawLatex(x, y+0.05, text)
      
      text_latex.SetTextSize(0.04)
      if thermal:
        text = 'PYTHIA8 Monash 2013 embedded in thermal background'
      else:
        text = 'PYTHIA8 Monash 2013 embedded in 0#font[122]{-}10% Pb#font[122]{-}Pb'
      text_latex.DrawLatex(x, y, text)
        
    text = '#sqrt{#it{s_{#it{NN}}}} = 5.02 TeV'
    text_latex.DrawLatex(x, y-0.05, text)

    text = 'Charged jets anti-#it{k}_{T}'
    text_latex.DrawLatex(x, y-0.1, text)
    
    text = '#it{R} = ' + str(jetR) + '   | #it{{#eta}}_{{jet}}| < {:.1f}'.format(self.eta_max - jetR)
    text_latex.DrawLatex(x, y-0.15, text)
    
    if 'leading' in name_prefix:
      text = 'Leading subjets'
    else:
      text = 'Inclusive subjets'
    text_latex.DrawLatex(x, y-0.21, text)
    
    outdir = 'matched_pt_fraction_pt'
    output_filename = os.path.join(self.output_dir, '{}/{}_{}.pdf'.format(outdir, self.remove_periods(name_prefix), i_overlay))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  def plot_subjet_money_plot(self, observable, jetR, R_max, prong_match_threshold, obs_setting,
                             pt_bins_reported, output_dir, ytitle, option='', thermal=False):
  
    name = 'h_{}_matched_pt_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, obs_setting, R_max)
    h3D = self.fMC.Get(name)
    #h3D.Sumw2()
  
    # Loop through pt slices, and plot 1D projection onto observable
    for i in range(0, len(pt_bins_reported) - 1):
      min_pt = pt_bins_reported[i]
      max_pt = pt_bins_reported[i+1]
  
      # Set pt range
      h3D.GetXaxis().SetRangeUser(min_pt, max_pt)
  
      self.plot_obs_projection_subjet(observable, h3D, jetR, prong_match_threshold, obs_setting,
                                      min_pt, max_pt, ytitle, output_dir, thermal=thermal)

  #---------------------------------------------------------------
  def plot_obs_projection_subjet(self, observable, h3D, jetR, prong_match_threshold, obs_setting,
                                 min_pt, max_pt, ytitle, output_dir, thermal=False):
  
    # Reset projections for normalization
  
    # z
    h3D.GetYaxis().SetRangeUser(0., self.xmax)
  
    # matching fraction
    h3D.GetZaxis().SetRangeUser(0., 1.)

    # Get normalization for embeddeed case
    h_normalization = h3D.Project3D('y')
    h_normalization.SetName('h_normalization_emb_{}_{}_{}-{}'.format(observable, obs_setting, min_pt, max_pt))
    N_inclusive_emb = h_normalization.Integral()
    if N_inclusive_emb < 1e-3:
      print('Emb Integral is 0, check for problem')
      return
  
    # Draw histogram
    c = ROOT.TCanvas('c','c: hist', 600, 650)
    c.cd()
  
    pad1 = ROOT.TPad('myPad', 'The pad',0.,0.3,1.,1.)
    pad1.SetLeftMargin(0.2)
    pad1.SetTopMargin(0.04)
    pad1.SetRightMargin(0.04)
    pad1.SetBottomMargin(0.)
    pad1.Draw()
    pad1.cd()
  
    leg = ROOT.TLegend(0.7,0.6,0.85,0.8)
    self.setup_legend(leg, 0.04)
    leg.SetHeader('#splitline{PYTHIA subjet}{reconstructed as:}')
    leg.AddEntry(None, '', '')

    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 20, self.xmin, self.xmax)
    myBlankHisto.SetNdivisions(505)
    #myBlankHisto.SetMinimum(0.01)
    xtitle = h3D.GetYaxis().GetTitle()
    myBlankHisto.GetXaxis().SetTitle(xtitle)
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

    leg_ratio = ROOT.TLegend(0.65,0.77,0.8,0.92)
    self.setup_legend(leg_ratio, 0.1)

    myBlankHisto2 = myBlankHisto.Clone('myBlankHisto_C')
    ytitle = 'match / total'
    myBlankHisto2.SetYTitle(ytitle)
    myBlankHisto2.SetXTitle(xtitle)
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
    myBlankHisto2.SetMinimum(0.01)
    myBlankHisto2.SetMaximum(1.99)
    myBlankHisto2.Draw('')
  
    pad1.cd()
    h_stack = ROOT.THStack('h_stack', 'stacked')
    h_sum = None
    h_sum_tagged = None

    # Loop over each flag
    for i in range(2):
      if i == 0:
        h3D.GetZaxis().SetRangeUser(prong_match_threshold, 1.)
        legend_label = 'match'
      elif i == 1:
        h3D.GetZaxis().SetRangeUser(0., prong_match_threshold)
        legend_label = 'mis-tag'

      # Project onto 1D
      h1D = h3D.Project3D('y')
      h1D.SetName('h1D_{}'.format(i))
      h1D.Rebin(self.rebin_value)
      h1D.Scale(1./N_inclusive_emb, 'width')

      h1D.SetLineColor(self.ColorArray[i])
      h1D.SetFillColor(self.ColorArray[i])
      h1D.SetLineWidth(2)

      h_stack.Add(h1D)
      leg.AddEntry(h1D, legend_label, 'F')

      if h_sum:
        h_sum.Add(h1D)
      else:
        h_sum = h1D.Clone()
        h_sum.SetName('h_sum_{}_{}_{}-{}'.format(observable, obs_setting, min_pt, max_pt))

      if i == 0:
        if h_sum_tagged:
          h_sum_tagged.Add(h1D)
        else:
          h_sum_tagged = h1D.Clone()
          h_sum_tagged.SetName('h_sum_tagged_{}_{}_{}-{}'.format(observable, obs_setting, min_pt, max_pt))
  
    if 0.5*h_sum.GetMaximum() > myBlankHisto.GetMaximum():
      myBlankHisto.SetMaximum(2*h_sum.GetMaximum())
  
    leg.Draw('same')
    h_stack.Draw('same hist')
  
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
  
    x = 0.23
    y = 0.9
    text_latex.SetTextSize(0.05)
    if thermal:
        text = 'PYTHIA8 embedded in thermal background'
    else:
        text = 'PYTHIA8 embedded in 0-10% Pb-Pb'
    text_latex.DrawLatex(x, y, text)
  
    text = '#sqrt{#it{s_{#it{NN}}}} = 5.02 TeV'
    text_latex.DrawLatex(x, y-0.06, text)

    text = 'Charged jets   anti-#it{k}_{T}'
    text_latex.DrawLatex(x, y-0.12, text)
  
    text = '#it{R} = ' + str(jetR) + '   | #it{{#eta}}_{{jet}}| < {:.1f}'.format(self.eta_max - jetR)
    text_latex.DrawLatex(x, y-0.18, text)
  
    text = str(int(min_pt)) + ' < #it{p}_{T, ch jet}^{PYTHIA} < ' + str(int(max_pt)) + ' GeV/#it{c}'
    text_latex.DrawLatex(x, y-0.25, text)
  
    text = 'Leading subjets, #it{{r}} = {}'.format(obs_setting)
    text_latex.DrawLatex(x, y-0.32, text)
  
    # Build ratio of embedded tagging purity
    pad2.cd()
    setattr(self, 'h_ratio_tagged_{}'.format(obs_setting), None)
    h_ratio_tagged = h_sum_tagged.Clone()
    h_ratio_tagged.SetMarkerStyle(self.MarkerArray[1])
    h_ratio_tagged.SetMarkerColor(self.ColorArray[0])
    h_ratio_tagged.Divide(h_sum)
    if h_ratio_tagged.GetMinimum() < myBlankHisto2.GetMinimum():
      myBlankHisto2.SetMinimum(0.8*h_ratio_tagged.GetMinimum())
    h_ratio_tagged.Draw('PE same')
    leg_ratio.AddEntry(h_ratio_tagged, 'Tagging purity', 'P')
    setattr(self, 'h_ratio_tagged_{}_{}-{}'.format(obs_setting, min_pt, max_pt), h_ratio_tagged)

    leg_ratio.Draw('same')
  
    line = ROOT.TLine(self.xmin,1,self.xmax,1)
    line.SetLineColor(920+2)
    line.SetLineStyle(2)
    line.Draw('same')

    output_filename = os.path.join(output_dir, f'money_plot_{min_pt}-{max_pt}_R{jetR}_r{obs_setting}.pdf')
    c.SaveAs(output_filename)
    c.Close()
    
  #---------------------------------------------------------------
  def plot_prong_matching_delta(self, i_overlay, jetR, name_prefix, obs_subconfig_list, obs_settings, grooming_settings, overlay_list, prong_match_threshold, min_pt, max_pt, plot_deltaz=False, plot_matched=True):

    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()
    ROOT.gPad.SetLeftMargin(0.2)
    ROOT.gPad.SetBottomMargin(0.15)
    c.GetPad(0).SetTicks(0,1)
    
    xmin = 0.
    if plot_deltaz:
      xmin = -jetR
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 20, xmin, jetR)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.GetXaxis().SetTitleOffset(1.4)
    myBlankHisto.GetYaxis().SetTitleOffset(1.4)
    myBlankHisto.GetXaxis().SetTitleSize(0.05)
    myBlankHisto.GetYaxis().SetTitleSize(0.05)
    max = 0.01
    myBlankHisto.SetMaximum(2*max)
    myBlankHisto.SetMinimum(0.)
    if plot_deltaz:
        myBlankHisto.GetXaxis().SetTitle('#Delta#it{z}_{subjet}')
        myBlankHisto.GetYaxis().SetTitle('#frac{d#it{N}}{d#Delta #it{z}_{subjet}}')
    else:
        xtitle = '#Delta#it{R}_{PYTHIA, PYTHIA#oplusPb#font[122]{-}Pb}^{leading subjet}'
        myBlankHisto.GetXaxis().SetTitle(xtitle)
        myBlankHisto.GetYaxis().SetTitle('#frac{{d#it{{N}}}}{{d{}}}'.format(xtitle))
    myBlankHisto.Draw()
    
    leg = ROOT.TLegend(0.54,0.58,0.67,0.75)
    self.setup_legend(leg,0.035)
    
    h_list = [] # Store hists in a list, sincge otherwise it seems I lose the marker information
                # (removed from memory?)
    
    i_reset = 0
    for i, subconfig_name in enumerate(obs_subconfig_list):
    
      if subconfig_name not in overlay_list:
        continue

      obs_setting = obs_settings[i]
      grooming_setting = grooming_settings[i]
      obs_label = self.obs_label(obs_setting, grooming_setting)
      
      if subconfig_name == overlay_list[0]:
        marker = 20
      elif subconfig_name == overlay_list[1]:
        marker = 21
      elif i > 1 and subconfig_name == overlay_list[2]:
        marker = 33
      else:
        marker = 34
      
      name = '{}_{}{}{}'.format(name_prefix, obs_label, self.suffix, self.scaled_suffix)
      hFraction_vs_pt = self.fMC.Get(name)

      epsilon = 1e-5
      min_bin = hFraction_vs_pt.GetYaxis().FindBin(0. + epsilon)
      cut_bin = hFraction_vs_pt.GetYaxis().FindBin(prong_match_threshold + epsilon)
      max_bin = hFraction_vs_pt.GetYaxis().FindBin(1. + epsilon)
      min_frac_bin = cut_bin
      max_frac_bin = max_bin

      min_pt_bin = hFraction_vs_pt.GetXaxis().FindBin(min_pt + epsilon)
      max_pt_bin = hFraction_vs_pt.GetXaxis().FindBin(max_pt + epsilon)

      # Get projections of deltaR
      hFraction_vs_pt.GetXaxis().SetRange(min_pt_bin, max_pt_bin)
      if plot_matched:
        hFraction_vs_pt.GetYaxis().SetRange(min_frac_bin, max_frac_bin)
      else:
        hFraction_vs_pt.GetYaxis().SetRange(min_bin, min_frac_bin)
      hUnmatched_vs_pt = hFraction_vs_pt.Project3D('z')
      hUnmatched_vs_pt.SetName('hUnmatched_vs_pt{}_{}_{}'.format(i, jetR, obs_label))
      hUnmatched_vs_pt.SetLineColor(self.ColorArray[i_reset])
      i_reset += 1
      hUnmatched_vs_pt.SetLineWidth(4)
      if self.thermal:
        hUnmatched_vs_pt.SetLineColor(self.ColorArray[len(self.ColorArray)-i-2])
      if max < hUnmatched_vs_pt.GetMaximum():
        max = hUnmatched_vs_pt.GetMaximum()
        if jetR == 0.2:
          myBlankHisto.SetMaximum(2.4*max)
        else:
          myBlankHisto.SetMaximum(1.9*max)
      hUnmatched_vs_pt.Draw('L hist same')
      leg.AddEntry(hUnmatched_vs_pt, f'r = {obs_setting}', 'L')
      h_list.append(hUnmatched_vs_pt)
    
    leg.Draw('same')
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()

    if self.groomer_studies:
      x = 0.23
      y = 0.85
      text_latex.SetTextSize(0.04)
      text = 'PYTHIA8 embedded in thermal background'
      text_latex.DrawLatex(x, y, text)
    else:
      x = 0.23
      y = 0.8
      text_latex.SetTextSize(0.055)
      text = self.figure_approval_status
      text_latex.DrawLatex(x, y+0.05, text)
      
      text_latex.SetTextSize(0.04)
      text = 'PYTHIA8 Monash 2013 embedded in 0#font[122]{-}10% Pb#font[122]{-}Pb'
      text_latex.DrawLatex(x, y, text)
        
    text = '#sqrt{#it{s_{#it{NN}}}} = 5.02 TeV'
    text_latex.DrawLatex(x, y-0.05, text)

    text = 'Charged jets   anti-#it{k}_{T}'
    text_latex.DrawLatex(x, y-0.1, text)
    
    text = '#it{R} = ' + str(jetR) + '   | #it{{#eta}}_{{jet}}| < {:.1f}'.format(self.eta_max - jetR)
    text_latex.DrawLatex(x, y-0.15, text)
    
    text = str(int(min_pt)) + ' < #it{p}_{T, ch jet}^{pp-det} < ' + str(int(max_pt)) + ' GeV/#it{c}'
    text_latex.DrawLatex(x, y-0.2, text)
    
    text_latex.SetTextSize(0.04)
    text = 'PYTHIA leading subjet '
    text_latex.DrawLatex(0.23, 0.5, text)
    if plot_matched:
        text = 'tagged in PYTHIA#oplusPb#font[122]{-}Pb'
    else:
        text = 'mis-tagged in PYTHIA#oplusPb#font[122]{-}Pb'
    text_latex.DrawLatex(0.23, 0.46, text)
    if plot_matched:
        text = 'leading subjet'
        text_latex.DrawLatex(0.23, 0.42, text)

    if plot_matched:
        label = 'matched'
    else:
        label = 'unmatched'

    if plot_deltaz:
        output_filename = os.path.join(self.output_dir, f'prong_matching_deltaZ/{self.remove_periods(name_prefix)}_{label}_{i_overlay}.pdf')
    else:
        output_filename = os.path.join(self.output_dir, f'prong_matching_deltaR/{self.remove_periods(name_prefix)}_{label}_{i_overlay}.pdf')
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  def plot_subjet_DeltaR(self, jetR, subjetR, jet_matching_distance):

    if self.is_pp:
      name = f'hDeltaR_ppdet_pptrue_subjet_z_R{jetR}_{subjetR}{self.scaled_suffix}'
    else:
      name = f'hDeltaR_combined_ppdet_subjet_z_R{jetR}_{subjetR}{self.suffix}{self.scaled_suffix}'

    h = self.fMC.Get(name)
    
    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetRightMargin(0.15)
    c.SetLogz()

    h.GetXaxis().SetTitle('#it{p}_{T,det}^{ch jet}')
    h.GetYaxis().SetTitle('#DeltaR_{match}')
    
    x_max = 200.
    h.GetXaxis().SetRangeUser(0., x_max)
    h.GetYaxis().SetRangeUser(0., 1.)
    h.Draw('colz')
    h.Scale(1e6)
    
    deltaR_max = jet_matching_distance * float(subjetR)
    line = ROOT.TLine(0, deltaR_max, x_max, deltaR_max)
    line.SetLineColor(920+2)
    line.SetLineStyle(2)
    line.SetLineWidth(4)
    line.Draw()
    
    text = f'R = {jetR}, r = {subjetR}'
    textFit = ROOT.TLatex()
    textFit.SetTextSize(0.04)
    textFit.SetNDC()
    textFit.DrawLatex(0.62,0.8,text)

    output_filename = os.path.join(self.output_dir, 'jet/{}.pdf'.format(self.remove_periods(name)))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  def plot_z1_crosscheck(self, jetR, obs_label, obs_setting, grooming_setting,
                         xtitle, pt_bins):
            
    name = 'h_{}_zconst_R{}_{}_z099_1{}'.format(self.observable, jetR, obs_label, self.suffix)
    h2D = self.fData.Get(name)
    
    for i in range(0, len(pt_bins) - 1):
      min_pt_truth = pt_bins[i]
      max_pt_truth = pt_bins[i+1]
      
      self.plot_z1_projection(h2D, jetR, obs_label, obs_setting, grooming_setting, xtitle, min_pt_truth, max_pt_truth)
    
  #---------------------------------------------------------------
  def plot_z1_projection(self, h2D, jetR, obs_label, obs_setting,
                            grooming_setting, xtitle, min_pt, max_pt):
    
    ytitle = '#frac{dN}{dz}'

    # Get histogram of theta_g in data, for given pt-det cut
    h2D.GetXaxis().SetRangeUser(min_pt, max_pt)
    hObs_truth = h2D.ProjectionY()
    hObs_truth.SetMarkerStyle(21)
    hObs_truth.SetMarkerSize(1)
    #hObs_truth.Rebin(5)
    #hObs_truth.Scale(1., 'width')
    if grooming_setting and 'sd' in grooming_setting:
      hObs_truth.GetXaxis().SetRange(0, hObs_truth.GetNbinsX())

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
    
    hObs_truth.GetYaxis().SetTitle(ytitle)
    hObs_truth.GetYaxis().SetTitleOffset(1.5)
    hObs_truth.SetMaximum(2.5*hObs_truth.GetMaximum())
    hObs_truth.SetMinimum(0.)

    hObs_truth.Draw('hist')
    leg.AddEntry(hObs_truth, "Data", "L")
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

    text = str(min_pt) + ' < #it{p}_{T, ch jet} < ' + str(max_pt) + ' GeV/#it{c}'
    text_latex.DrawLatex(0.3, 0.73, text)

    text = '#it{R} = ' + str(jetR) + '   | #it{{#eta}}_{{jet}}| < {:.1f}'.format(self.eta_max - jetR)
    text_latex.DrawLatex(0.3, 0.67, text)
    
    subobs_label = self.formatted_subobs_label(self.observable)
    delta = 0.
    if subobs_label:
      text = '{} = {}'.format(subobs_label, obs_setting)
      text_latex.DrawLatex(0.3, 0.61, text)
      delta = 0.07
      
    text = f'0.99 < {xtitle} < 1.0'
    text_latex.DrawLatex(0.5, 0.5, text)
      
    output_filename = os.path.join(self.output_dir, 'z1_crosscheck/h_{}_z1_crosscheck_R{}_{}_{}-{}.pdf'.format(self.observable, self.remove_periods(jetR), obs_label, min_pt, max_pt))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  # Plot fraction of det-level subjets without a unique match
  # (e.g. for pt=80-120), as a function of z
  def plot_subjet_matching_pp(self, jetR, obs_label, obs_setting, grooming_setting,
                              xtitle, pt_bins_reported):
    
    # (jet_truth.pt(), z_det, successful_match)
    name = f'h_match_fraction_{self.observable}_R{jetR}_{obs_label}{self.scaled_suffix}'
    h3D = self.fMC.Get(name)

    # Select pt range
    h3D.GetXaxis().SetRangeUser(80., 120.)
    
    # Make projection for successful_match = True or False
    h3D_denominator = h3D.Clone(f'h_denominator_{h3D.GetName()}')
    h_denominator = h3D_denominator.Project3D('y')
    
    # Make successful_match=True projection
    h3D_numerator = h3D.Clone(f'h_numerator_{h3D.GetName()}')
    h3D_numerator.GetZaxis().SetRange(2, 2)
    h_numerator = h3D_numerator.Project3D('y')

    h_ratio = h_numerator.Clone()
    h_ratio.SetName(f'h_ratio_{h_numerator.GetName()}')
    h_ratio.Divide(h_denominator)
    
    h_ratio.GetYaxis().SetRangeUser(0, 1.5)
    h_ratio.SetMarkerStyle(21)
    h_ratio.SetMarkerSize(1)
    h_ratio.GetYaxis().SetTitle('Matching fraction')
    text=f'Inclusive subjets R={jetR},r={obs_setting}, 80.<p_{{T,jet}}<120.'

    output_filename = os.path.join(self.output_dir, f'subjet_matching_pp/h_match_fraction_{self.observable}_R{jetR}_{obs_label}{self.suffix}.pdf')
    self.plot_hist(h_ratio, output_filename, text=text)
