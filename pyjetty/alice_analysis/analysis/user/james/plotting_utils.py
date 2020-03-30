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

# Data analysis and plotting
import numpy as np
from array import *
import ROOT

# Base class
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs

################################################################
class PlottingUtils(analysis_utils_obs.AnalysisUtils_Obs):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, observable='', is_pp=False, fnameData='', fnameMC='', output_dir='.', **kwargs):
    super(PlottingUtils, self).__init__(**kwargs)
    
    self.observable = observable
    self.is_pp = is_pp
    self.fData = ROOT.TFile(fnameData, 'READ')
    self.fMC = ROOT.TFile(fnameMC, 'READ')
    self.output_dir = output_dir
    
    print(self)

  #---------------------------------------------------------------
  def plot_DeltaR(self, jetR, jet_matching_distance):

    if self.is_pp:
      name = 'hDeltaR_All_R{}Scaled'.format(jetR)
    else:
      name = 'hDeltaR_combined_ppdet_R{}Scaled'.format(jetR)
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
    
    deltaR_max = jet_matching_distance * jetR
    line = ROOT.TLine(0, deltaR_max, x_max, deltaR_max)
    line.SetLineColor(920+2)
    line.SetLineStyle(2)
    line.SetLineWidth(4)
    line.Draw()
    
    text = 'R = {}'.format(jetR)
    textFit = ROOT.TLatex()
    textFit.SetTextSize(0.04)
    textFit.SetNDC()
    textFit.DrawLatex(0.72,0.8,text)

    output_filename = os.path.join(self.output_dir, '{}.pdf'.format(self.remove_periods(name)))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  def plot_JES(self, jetR):
    
    name = 'hJES_R{}Scaled'.format(jetR)
    histDeltaJES = self.fMC.Get(name)
    histDeltaJES.GetXaxis().SetTitle("#it{p}_{T}^{gen}")
    histDeltaJES.GetYaxis().SetTitle("#frac{#it{p}_{T}^{det} - #it{p}_{T}^{gen}}{#it{p}_{T}^{gen}}")
    histDeltaJES.GetXaxis().SetRangeUser(0., 200.)
    outputFilename = os.path.join(self.output_dir, "histDeltaJES_R{}.pdf".format(self.remove_periods(jetR)))
    #self.utils.plot_hist(histDeltaJES, outputFilename, "colz", False, True)
    
    histDeltaJES.RebinX(4)
    histDeltaJESprof = histDeltaJES.ProfileX()
    histDeltaJESprof.GetXaxis().SetRangeUser(0., 200.)
    histDeltaJESprof.GetYaxis().SetRangeUser(-0.8, 0.2)
    histDeltaJESprof.GetYaxis().SetTitle("#frac{#it{p}_{T}^{det} - #it{p}_{T}^{gen}}{#it{p}_{T}^{gen}}")
    outputFilename = os.path.join(self.output_dir, "histDeltaJESprof_R{}.pdf".format(self.remove_periods(jetR)))
    histDeltaJESprof.SetMarkerStyle(21)
    histDeltaJESprof.SetMarkerColor(4)
    self.plot_hist(histDeltaJESprof, outputFilename, "P")

  #---------------------------------------------------------------
  # Plot JES shift distribution for various fixed pT-gen
  def plot_JES_proj(self, jetR, pt_bins):
    
    name = 'hJES_R{}Scaled'.format(jetR)
    
    cJES = ROOT.TCanvas('cJES','cJES: hist',600,450)
    cJES.cd()
    cJES.SetBottomMargin(0.2)
    
    leg = ROOT.TLegend(0.55,0.55,0.88,0.85, '')
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1)
    leg.SetTextSize(0.04)
    
    kGreen = 416
    kBlue  = 600
    kCyan  = 432
    kAzure   = 860
    kViolet  = 880
    kMagenta = 616
    kPink    = 900
    ColorArray = [kBlue-4, kAzure+7, kCyan-2,7, kViolet-8, kBlue-6, kGreen+3]
    
    # Loop through pt slices, and plot final result for each 1D theta_g distribution
    for i in range(0, len(pt_bins) - 1):
      min_pt_truth = pt_bins[i]
      max_pt_truth = pt_bins[i+1]
      
      hJESProj = self.getJESshiftProj(name, 'hJESproj1', min_pt_truth, max_pt_truth)
      hJESProj.SetMarkerStyle(20)
      hJESProj.SetMarkerColor(ColorArray[i])
      hJESProj.SetLineColor(ColorArray[i])
      
      if i == 0:
      
        hJESProj.GetXaxis().SetTitleOffset(1.6);
        hJESProj.GetYaxis().SetTitle('Probability density')
        hJESProj.GetXaxis().SetTitle('#frac{#it{p}_{T}^{det} - #it{p}_{T}^{gen}}{#it{p}_{T}^{gen}}')
        
        hJESProj.GetYaxis().SetRangeUser(0, 10.4)
        hJESProj.DrawCopy('P E')
        
      else:
      
        hJESProj.DrawCopy('P E same')
    
      leg.AddEntry(hJESProj, '#it{{p}}_{{T}}^{{gen}} = {}-{} GeV'.format(min_pt_truth, max_pt_truth), 'P')
      
    leg.Draw('same')
    
    outputFilename = os.path.join(self.output_dir, 'histDeltaJESproj_R{}.pdf'.format(self.remove_periods(jetR)))
    cJES.SaveAs(outputFilename)
    cJES.Close()

  #---------------------------------------------------------------
  # Get JES shift distribution for a fixed pT-gen
  def getJESshiftProj(self, name, label, minPt, maxPt):
    
    histDeltaJES = self.fMC.Get(name)
    histDeltaJES.SetName('{}_{}'.format(histDeltaJES.GetName(), label))
    
    histDeltaJES.GetXaxis().SetRangeUser(minPt, maxPt)
    histDeltaJES.GetYaxis().SetRangeUser(-1., 1.)
    h = histDeltaJES.ProjectionY()
    
    integral = h.Integral()
    if integral > 0:
      h.Scale(1./integral, 'width')
    
    return h
