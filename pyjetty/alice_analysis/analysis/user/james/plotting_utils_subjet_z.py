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
          myLegend.AddEntry(hMatchedFraction_vs_pt, '#it{{R}}_{{subjet}} = {}'.format(obs_setting), 'P')
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
    
    outdir = 'prong_matching_fraction_pt'
    if 'leading' in name_prefix:
      outdir = 'prong_matching_fraction_pt_leading'
    output_filename = os.path.join(self.output_dir, '{}/{}_{}.pdf'.format(outdir, self.remove_periods(name_prefix), i_overlay))
    c.SaveAs(output_filename)
    c.Close()
