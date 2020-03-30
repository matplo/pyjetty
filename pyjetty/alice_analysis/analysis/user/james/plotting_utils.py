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
    histDeltaJESprof.SetMarkerSize(3)
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

  #---------------------------------------------------------------
  def plotJER(self, jetR, obs_label):
    
    # (pt-det, pt-truth, theta_g-det, theta_g-truth)
    name = 'hResponse_JetPt_{}_R{}_{}Scaled'.format(self.observable, jetR, obs_label)
    hRM_4d = self.fMC.Get(name)
    hRM = hRM_4d.Projection(1,0)
    hRM.SetName('hResponse_JetPt_{}_R{}_{}_Proj'.format(self.observable, jetR, obs_label))
    
    # For each pT^gen, compute the standard deviation of the pT^det distribution
    
    # Get the pT^gen profile, with errors as standard deviation of pT^det distribution
    histPtGenProf = hRM.ProfileY('histPtGenProff', 1, -1, "s")
    
    # Create histo to be used to fill JER values
    nBins = 60
    histJER = ROOT.TH1D('histJER_R{}_{}'.format(jetR, obs_label), 'histJER_R{}_{}'.format(jetR, obs_label), nBins, 0., 300.) # same binning for pT^gen as in task
    
    # Loop through the bins, and fill the JER
    for i in range(0,nBins+1):
      sigma = histPtGenProf.GetBinError(i)
      pTgen = histPtGenProf.GetXaxis().GetBinCenter(i)
      JER = sigma/pTgen
      histJER.SetBinContent(i, JER)
    
    histJER.GetYaxis().SetTitle('#frac{#sigma(#it{p}_{T}^{gen})}{#it{p}_{T}^{gen}}')
    histJER.GetXaxis().SetTitle('#it{p}_{T}^{gen}')
    histJER.GetYaxis().SetRangeUser(-0.01, 0.5)
    histJER.GetXaxis().SetRangeUser(5., 100.)
    outputFilename = os.path.join(self.output_dir, 'histJER_R{}.pdf'.format(self.remove_periods(jetR)))
    histJER.SetMarkerStyle(21)
    histJER.SetMarkerColor(2)
    histJER.SetMarkerSize(3)
    self.plot_hist(histJER, outputFilename, 'hist P')
  
  #---------------------------------------------------------------
  def plot_obs_resolution(self, jetR, obs_label, xtitle, pt_bins):
    
    c_resolution = ROOT.TCanvas('cres_{}'.format(obs_label),'cres_{}: hist'.format(obs_label),600,450)
    c_resolution.cd()
    c_resolution.SetBottomMargin(0.17)
    c_resolution.SetLeftMargin(0.13)
    
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 20, 0., 1.)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.GetYaxis().SetTitleOffset(1.2)
    myBlankHisto.SetMaximum(2.)
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw()
    label_gen = '{}^{{gen}}'.format(xtitle)
    myBlankHisto.GetYaxis().SetTitle('#sigma({}) / {}'.format(label_gen, label_gen))
    myBlankHisto.GetXaxis().SetTitle(label_gen)
    
    myLegend = ROOT.TLegend(0.55,0.55,0.88,0.85)
    self.setup_legend(myLegend,0.035)
    
    kGreen = 416
    kBlue  = 600
    kCyan  = 432
    kAzure   = 860
    kViolet  = 880
    kMagenta = 616
    kPink    = 900
    ColorArray = [kBlue-4, kAzure+7, kCyan-2, 7, kViolet-8, kBlue-6, kGreen+3]
    
    # (pt-det, pt-truth, theta_g-det, theta_g-truth)
    name = 'hResponse_JetPt_{}_R{}_{}Scaled'.format(self.observable, jetR, obs_label)
    hRM_4d = self.fMC.Get(name)
    
    for i in range(0, len(pt_bins) - 1):
    
      min_pt_truth = pt_bins[i]
      max_pt_truth = pt_bins[i+1]
      
      hRM_4d_clone = hRM_4d.Clone()
      hRM_4d_clone.SetName('{}_{}'.format(hRM_4d_clone.GetName(), i))
      hResolution = self.get_resolution(hRM_4d_clone, jetR, obs_label, min_pt_truth, max_pt_truth, 'hResolution_{}'.format(i))
 
      hResolution.SetMarkerColor(ColorArray[i])
      hResolution.SetMarkerStyle(21)
      hResolution.SetLineColor(ColorArray[i])
      hResolution.DrawCopy('P same')

      myLegend.AddEntry(hResolution, '#it{{p}}_{{T}}^{{gen}} = {}-{} GeV'.format(min_pt_truth, max_pt_truth), 'P')

      myLegend.Draw('same')

    outputFilename = os.path.join(self.output_dir, 'hResolution_R{}_{}.pdf'.format(self.remove_periods(jetR), obs_label))
    c_resolution.SaveAs(outputFilename)
    c_resolution.Close()

  #---------------------------------------------------------------
  # Get resolution for a fixed pT-gen
  def get_resolution(self, hRM_4d, jetR, obs_label, minPt, maxPt, label):

    hRM = hRM_4d.GetAxis(1).SetRangeUser(minPt, maxPt)
    hRM = hRM_4d.Projection(3,2)
    hrm_name = 'hResponse_JetPt_{}_R{}_{}_{}-{}_{}Proj'.format(self.observable, jetR, obs_label, minPt, maxPt, label)
    hRM.SetName(hrm_name)
    
    # For each pT^gen, compute the standard deviation of the pT^det distribution
    
    # Get the pT^gen profile, with errors as standard deviation of pT^det distribution
    hprof_name = 'histGenProff_{}_{}-{}'.format(obs_label, minPt, maxPt)
    histGenProf = hRM.ProfileY(hprof_name, 1, -1, 's')
    histGenProf.SetName(hprof_name)
    
    # Create histo to be used to fill JER values
    nBins = 20
    hres_name = 'histResolution_R{}_{}_{}-{}_{}'.format(jetR, obs_label, minPt, maxPt, label)
    histResolution = ROOT.TH1D(hres_name, hres_name, nBins, 0., 1.) # same binning for pT^gen as in task
    
    # Loop through the bins, and fill the JER
    for i in range(0,nBins+1):
      sigma = histGenProf.GetBinError(i)
      obs_gen = histGenProf.GetXaxis().GetBinCenter(i)
      resolution = sigma/obs_gen
      histResolution.SetBinContent(i, resolution)
      
    return histResolution

  #---------------------------------------------------------------
  def plot_obs_residual(self, jetR, obs_label, xtitle, pt_bins):

    name = 'hResidual_JetPt_{}_R{}_{}Scaled'.format(self.observable, jetR, obs_label)
    
    c_residual = ROOT.TCanvas('c','c: hist',600,450)
    c_residual.cd()
    c_residual.SetBottomMargin(0.17)
    
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
    ColorArray = [kBlue-4, kAzure+7, kCyan-2, 7, kViolet-8, kBlue-6, kGreen+3]
    
    # Loop through pt slices, and plot final residual for each 1D distribution
    for i in range(0, len(pt_bins) - 1):
      min_pt_truth = pt_bins[i]
      max_pt_truth = pt_bins[i+1]
      
      hResidual = self.get_residual_proj(name, 'hResidual{}'.format(i), min_pt_truth, max_pt_truth)
      hResidual.SetMarkerStyle(20)
      hResidual.SetMarkerColor(ColorArray[i])
      hResidual.SetLineColor(ColorArray[i])

      if i == 0:
      
        hResidual.GetXaxis().SetTitleOffset(1.6);
        hResidual.GetYaxis().SetTitle('Probability density')
        hResidual.GetYaxis().SetRangeUser(0, 1.2*hResidual.GetMaximum())
        hResidual.DrawCopy('P E')
        
      else:

        hResidual.DrawCopy('P E same')

      leg.AddEntry(hResidual, '#it{{p}}_{{T}}^{{gen}} = {}-{} GeV'.format(min_pt_truth, max_pt_truth), 'P')

    leg.Draw('same')
    
    outputFilename = os.path.join(self.output_dir, 'hResidual_R{}_{}.pdf'.format(self.remove_periods(jetR), obs_label))
    c_residual.SaveAs(outputFilename)
    c_residual.Close()

  #---------------------------------------------------------------
  # Get residual for a fixed pT-gen
  def get_residual_proj(self, name, label, minPt, maxPt):
    
    h_residual_pt = self.fMC.Get(name)
    h_residual_pt.SetName('{}_{}'.format(h_residual_pt.GetName(), label))
    
    h_residual_pt.GetXaxis().SetRangeUser(minPt, maxPt)
    h_residual_pt.GetYaxis().SetRangeUser(-0.5, 0.5)
    h = h_residual_pt.ProjectionY()
    
    integral = h.Integral()
    if integral > 0:
      h.Scale(1./integral, 'width')
    
    return h
