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
  def plot_subjet_matching(self, i_overlay, jetR, name_prefix, obs_subconfig_list, obs_settings, grooming_settings, overlay_list, prong_match_threshold):

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
                             pt_bins_reported, output_dir, ytitle, option=''):
  
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
                                      min_pt, max_pt, ytitle, output_dir)

  #---------------------------------------------------------------
  def plot_obs_projection_subjet(self, observable, h3D, jetR, prong_match_threshold, obs_setting,
                                 min_pt, max_pt, ytitle, output_dir):
  
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
    xtitle = h3D.GetXaxis().GetTitle()
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
    ytitle = 'Ratio'
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
    text = 'PYTHIA8 embedded in thermal background'
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

    output_filename = os.path.join(output_dir, '{}/{}/money_plot_{}-{}.pdf'.format(jetR, obs_setting,
                                          min_pt, max_pt))
    c.SaveAs(output_filename)
    c.Close()
