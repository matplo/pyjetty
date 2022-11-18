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
  def plot_lund_plane(self, jetR, obs_label, grooming_setting):

    name = 'hLundPlane_R{}_{}{}{}'.format(jetR, obs_label, self.suffix, self.scaled_suffix)
    hLund = self.fMC.Get(name)
    
    hLund.GetXaxis().SetRangeUser(np.log(1/jetR), 5)
    hLund.GetYaxis().SetRangeUser(-3., 6.)
    
    text = '#it{p}_{T, ch jet}^{truth} > 100 GeV/c'
    output_filename = os.path.join(self.output_dir, 'lund/hLundPlane_R{}_{}.pdf'.format(jetR, obs_label))
    self.plot_hist(hLund, output_filename, drawOptions = 'colz', text = text)

  #---------------------------------------------------------------
  def plot_prong_matching(self, i_overlay, jetR, name_prefix, obs_subconfig_list, obs_settings, grooming_settings, overlay_list, prong_match_threshold):
  
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
    myBlankHisto.SetMaximum(1.61)
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.GetYaxis().SetTitle('Subleading prong purity')
    myBlankHisto.Draw()
    
    if self.groomer_studies:
      myLegend = ROOT.TLegend(0.42,0.64,0.55,0.83)
      self.setup_legend(myLegend,0.035)
      myLegend2 = ROOT.TLegend(0.76,0.64,0.88,0.83)
      self.setup_legend(myLegend2,0.035)
    else:
      myLegend = ROOT.TLegend(0.42,0.2,0.75,0.4)
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
      
      name = '{}_{}{}{}'.format(name_prefix, obs_label, self.suffix, self.scaled_suffix)
      hFraction_vs_pt = self.fMC.Get(name)
      xtitle = hFraction_vs_pt.GetXaxis().GetTitle()
      #myBlankHisto.GetXaxis().SetTitle(xtitle)
      if self.groomer_studies:
        myBlankHisto.GetXaxis().SetTitle('#it{p}_{T,ch jet}^{PYTHIA} (GeV/#it{c})')
      else:
        myBlankHisto.GetXaxis().SetTitle('#it{p}_{T,ch jet}^{pp-det} (GeV/#it{c})')
      
      epsilon = 1e-5
      min_bin = hFraction_vs_pt.GetYaxis().FindBin(0. + epsilon)
      cut_bin = hFraction_vs_pt.GetYaxis().FindBin(prong_match_threshold + epsilon)
      max_bin = hFraction_vs_pt.GetYaxis().FindBin(1. + epsilon)
      
      name_total = 'hTotal_vs_pt{}_{}_{}'.format(i, jetR, obs_label)
      hTotal_vs_pt = hFraction_vs_pt.ProjectionX(name_total, min_bin, max_bin)
      hTotal_vs_pt.SetName(name_total)
      
      name_matched = 'hMatched_vs_pt{}_{}_{}'.format(i, jetR, obs_label)
      hMatched_vs_pt = hFraction_vs_pt.ProjectionX(name_matched, cut_bin, max_bin)
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
          myLegend.AddEntry(hMatchedFraction_vs_pt, self.formatted_grooming_label(grooming_setting, verbose = True), 'P')
        else:
          myLegend2.AddEntry(hMatchedFraction_vs_pt, self.formatted_grooming_label(grooming_setting, verbose = True), 'P')
      else:
        myLegend.AddEntry(hMatchedFraction_vs_pt, self.formatted_grooming_label(grooming_setting, verbose = not self.groomer_studies), 'P')
              
      h_list.append(hMatchedFraction_vs_pt)

    myLegend.Draw('same')
    if self.groomer_studies:
      #myLegend2.AddEntry(None, '', '')
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
    
    text = '{} reclustering'.format(self.reclustering_algorithm)
    text_latex.DrawLatex(x, y-0.21, text)
    
    outdir = 'prong_matching_fraction_pt'
    if 'JetPtDet' in name_prefix:
      outdir = 'prong_matching_fraction_ptdet'
    output_filename = os.path.join(self.output_dir, '{}/{}_{}.pdf'.format(outdir, self.remove_periods(name_prefix), i_overlay))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  def plot_prong_matching_delta(self, i_overlay, jetR, name_prefix, obs_subconfig_list, obs_settings, grooming_settings, overlay_list, prong_match_threshold, min_pt, max_pt, plot_deltaz=False):

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
    if jetR == 0.2:
      myBlankHisto.GetYaxis().SetTitleOffset(1.4)
    else:
      myBlankHisto.GetYaxis().SetTitleOffset(1.2)
    myBlankHisto.GetXaxis().SetTitleSize(0.05)
    myBlankHisto.GetYaxis().SetTitleSize(0.05)
    max = 0.01
    myBlankHisto.SetMaximum(2*max)
    myBlankHisto.SetMinimum(0.)
    if plot_deltaz:
        myBlankHisto.GetXaxis().SetTitle('#Delta#it{z}_{prong}')
        myBlankHisto.GetYaxis().SetTitle('#frac{d#it{N}}{d#Delta #it{z}_{prong}}')
    else:
        xtitle = '#Delta#it{R}_{PYTHIA, PYTHIA#oplusPb#font[122]{-}Pb}^{subleading prong}'
        myBlankHisto.GetXaxis().SetTitle(xtitle)
        myBlankHisto.GetYaxis().SetTitle('#frac{{d#it{{N}}}}{{d{}}}'.format(xtitle))
    myBlankHisto.Draw()
    
    leg = ROOT.TLegend(0.54,0.58,0.67,0.75)
    self.setup_legend(leg,0.035)
    
    h_list = [] # Store hists in a list, since otherwise it seems I lose the marker information
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
      hFraction_vs_pt.GetYaxis().SetRange(min_frac_bin, max_frac_bin)
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
      leg.AddEntry(hUnmatched_vs_pt, self.formatted_grooming_label(grooming_setting, verbose = not self.groomer_studies), 'L')
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
    text = 'PYTHIA subleading prong '
    text_latex.DrawLatex(0.23, 0.5, text)
    text = 'tagged in PYTHIA#oplusPb#font[122]{-}Pb'
    text_latex.DrawLatex(0.23, 0.46, text)
    text = 'leading prong'
    text_latex.DrawLatex(0.23, 0.42, text)

    if plot_deltaz:
        output_filename = os.path.join(self.output_dir, 'prong_matching_deltaZ/{}_{}.pdf'.format(self.remove_periods(name_prefix), i_overlay))
    else:
        output_filename = os.path.join(self.output_dir, 'prong_matching_deltaR/{}_{}.pdf'.format(self.remove_periods(name_prefix), i_overlay))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  def plot_prong_matching_correlation(self, i_overlay, jetR, hname_prefix, obs_subconfig_list, obs_settings, grooming_settings, overlay_list, prong_match_threshold):

    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.15)

    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 20, 0., 200.)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.GetYaxis().SetTitleOffset(1.2)
    myBlankHisto.SetMaximum(1.2)
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.GetYaxis().SetTitle('Fraction tagged but swapped')
    myBlankHisto.Draw()
    
    leg = ROOT.TLegend(0.55,0.3,0.85,0.5)
    self.setup_legend(leg,0.035)

    h_list = [] # Store hists in a list, since otherwise it seems I lose the marker information
                # (removed from memory?)
    
    i_reset = 0
    for i, subconfig_name in enumerate(obs_subconfig_list):
    
      if subconfig_name not in overlay_list:
        continue

      obs_setting = obs_settings[i]
      grooming_setting = grooming_settings[i]
      obs_label = self.obs_label(obs_setting, grooming_setting)

      name = '{}_{}{}{}'.format(hname_prefix, obs_label, self.suffix, self.scaled_suffix)
      hFraction_vs_pt = self.fMC.Get(name)

      epsilon = 1e-5
      min_bin = hFraction_vs_pt.GetYaxis().FindBin(0. + epsilon)
      cut_bin = hFraction_vs_pt.GetYaxis().FindBin(prong_match_threshold + epsilon)
      max_bin = hFraction_vs_pt.GetYaxis().FindBin(1. + epsilon)

      hTotal_vs_pt = hFraction_vs_pt.ProjectionX('hTotal_vs_pt{}_{}_{}'.format(i, jetR, obs_label), 0, -1, cut_bin, max_bin)

      hMatched_vs_pt = hFraction_vs_pt.ProjectionX('hMatched_vs_pt{}_{}_{}'.format(i, jetR, obs_label), cut_bin, max_bin, cut_bin, max_bin)

      hMatchedFraction_vs_pt = hMatched_vs_pt.Clone()
      hMatchedFraction_vs_pt.SetName('hMatchedFraction_vs_pt{}_{}_{}'.format(i, jetR, obs_label))
      hMatchedFraction_vs_pt.Divide(hTotal_vs_pt)
      hMatchedFraction_vs_pt.SetMarkerStyle(self.MarkerArray[i_reset])
      hMatchedFraction_vs_pt.SetMarkerColor(self.ColorArray[i_reset])
      hMatchedFraction_vs_pt.SetLineColor(self.ColorArray[i_reset])
      i_reset += 1
      hMatchedFraction_vs_pt.Draw('P same')
      leg.AddEntry(hMatchedFraction_vs_pt, self.formatted_grooming_label(grooming_setting, not self.groomer_studies), 'P')
      h_list.append(hMatchedFraction_vs_pt)

    leg.Draw('same')

    output_filename = os.path.join(self.output_dir, 'prong_matching_correlation/{}_{}.pdf'.format(self.remove_periods(hname_prefix), i_overlay))
    c.SaveAs(output_filename)
    c.Close()
