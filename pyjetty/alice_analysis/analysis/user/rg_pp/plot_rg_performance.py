#!/usr/bin/env python3

"""
  Plot MC histograms
  
  Some code is based on plotJetPerformance.py (James Mulligan + Eliane Epple)
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse
import math
from array import *

# Data analysis and plotting
import ROOT
import yaml

# Analysis utilities
from pyjetty.alice_analysis.analysis.base import common_base
from pyjetty.alice_analysis.analysis.base import analysis_utils

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

# Set plotting options
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

################################################################
class plot_rg_performance(common_base.common_base):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, mc_file = '', data_file = '', config_file = '', output_dir = '', **kwargs):
    super(plot_rg_performance, self).__init__(**kwargs)
    self.fMC = ROOT.TFile(mc_file, 'READ')
    self.fData = ROOT.TFile(data_file, 'READ')
    self.config_file = config_file
    self.output_dir = output_dir
    
    # Initialize utils class
    self.utils = analysis_utils.analysis_utils()
    
    # Create output dir
    if not self.output_dir.endswith('/'):
      self.output_dir = self.output_dir + '/'
    if not os.path.exists(self.output_dir):
      os.makedirs(self.output_dir)
    
    # Initialize yaml config
    self.initialize_config()
    
    print(self)

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    
    # Set debug level (0 = no debug info, 1 = some debug info, 2 = all debug info)
    self.debug_level = config['debug_level']
    
    # Retrieve list of SD grooming settings
    self.jetR_list = config['jetR']
    self.jet_matching_distance = config['jet_matching_distance']

    sd_config_dict = config['SoftDrop']
    sd_config_list = list(sd_config_dict.keys())
    self.sd_settings = [[sd_config_dict[name]['zcut'], sd_config_dict[name]['beta']] for name in sd_config_list]
    
    # Retrieve histogram binnings for each SD setting
    for i, sd_setting in enumerate(self.sd_settings):
      
      zcut = sd_setting[0]
      beta = sd_setting[1]
      sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)
      config_name = sd_config_list[i]
      
      pt_bins_truth = (sd_config_dict[config_name]['pt_bins_truth'])
      rg_bins_truth = (sd_config_dict[config_name]['rg_bins_truth'])
      
      n_pt_bins_truth = len(pt_bins_truth) - 1
      setattr(self, 'n_pt_bins_truth_{}'.format(sd_label), n_pt_bins_truth)
      
      truth_pt_bin_array = array('d',pt_bins_truth)
      setattr(self, 'truth_pt_bin_array_{}'.format(sd_label), truth_pt_bin_array)
      
      n_rg_bins_truth = len(rg_bins_truth) - 1
      setattr(self, 'n_rg_bins_truth_{}'.format(sd_label), n_rg_bins_truth)
      
      truth_rg_bin_array = array('d',rg_bins_truth)
      setattr(self, 'truth_rg_bin_array_{}'.format(sd_label), truth_rg_bin_array)
  
  #---------------------------------------------------------------
  def plot_rg_performance(self):
    
    for jetR in self.jetR_list:

      self.plotDeltaR(jetR)
      self.plotJES(jetR)
      self.plotJESproj(jetR)

      for sd_setting in self.sd_settings:
        
        zcut = sd_setting[0]
        beta = sd_setting[1]
        sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)
        
        self.plotJER(jetR, sd_label)
        self.plot_theta_resolution(jetR, sd_label)
        self.plot_theta_residual(jetR, sd_label)
        self.plot_zg_residual(jetR, sd_label)
        self.plotJetRecoEfficiency(jetR, sd_label)
        self.plotRg(jetR, sd_label, zcut, beta)

  #---------------------------------------------------------------
  def plotRg(self, jetR, sd_label, zcut, beta):

    # (pt-det, pt-truth, theta_g-det, theta_g-truth)
    name = 'hResponse_JetPt_ThetaG_R{}_{}Scaled'.format(jetR, sd_label)
    hRM_theta = self.fMC.Get(name)
    
    # (pt-det, pt-truth, theta_g-det, theta_g-truth)
    name = 'hResponse_JetPt_zg_R{}_{}Scaled'.format(jetR, sd_label)
    hRM_zg = self.fMC.Get(name)
    
    if self.fData:
      name = 'hThetaG_JetPt_R{}_{}'.format(jetR, sd_label)
      hThetaG_JetPt = self.fData.Get(name)
      hThetaG_JetPt.Sumw2()
    
      name = 'hZg_JetPt_R{}_{}'.format(jetR, sd_label)
      hZg_JetPt = self.fData.Get(name)
      hZg_JetPt.Sumw2()
    else:
      hThetaG_JetPt = None
      hZg_JetPt = None

    self.plot2D_theta_statistics(hThetaG_JetPt.Clone(), jetR, sd_label)
    self.plot2D_zg_statistics(hZg_JetPt.Clone(), jetR, sd_label)

    self.plotRgProjection(hRM_theta, hThetaG_JetPt, jetR, sd_label, zcut, beta, 20, 40)
    self.plotRgProjection(hRM_theta, hThetaG_JetPt, jetR, sd_label, zcut, beta, 40, 60)
    self.plotRgProjection(hRM_theta, hThetaG_JetPt, jetR, sd_label, zcut, beta, 60, 80)

    self.plotZgProjection(hRM_zg, hZg_JetPt, jetR, sd_label, zcut, beta, 20, 40)
    self.plotZgProjection(hRM_zg, hZg_JetPt, jetR, sd_label, zcut, beta, 40, 60)
    self.plotZgProjection(hRM_zg, hZg_JetPt, jetR, sd_label, zcut, beta, 60, 80)

  #---------------------------------------------------------------
  def plot2D_theta_statistics(self, hThetaG_JetPt, jetR, sd_label):

    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()
    ROOT.gPad.SetLeftMargin(0.15)
    
    hThetaG_JetPt.SetMarkerSize(0.5)
    hThetaG_JetPt.GetYaxis().SetRangeUser(0, 1.)
    hThetaG_JetPt.GetXaxis().SetRangeUser(0, 100)
    hThetaG_JetPt.RebinX(5)
    hThetaG_JetPt.RebinY(5)
    hThetaG_JetPt.Draw('text colz')

    output_filename = os.path.join(self.output_dir, 'h2D_theta_statistics_R{}_{}.pdf'.format(self.utils.remove_periods(jetR), sd_label))
    c.SaveAs(output_filename)
    c.Close()
  
  #---------------------------------------------------------------
  def plot2D_zg_statistics(self, hZg_JetPt, jetR, sd_label):
    
    czg = ROOT.TCanvas("czg","czg: hist",600,450)
    czg.cd()
    ROOT.gPad.SetLeftMargin(0.15)
    
    hZg_JetPt.SetMarkerSize(0.5)
    hZg_JetPt.GetYaxis().SetRangeUser(0, 0.5)
    hZg_JetPt.GetXaxis().SetRangeUser(0, 100)
    hZg_JetPt.RebinX(5)
    hZg_JetPt.RebinY(5)
    hZg_JetPt.Draw('text colz')
    
    output_filename = os.path.join(self.output_dir, 'h2D_zg_statistics_R{}_{}.pdf'.format(self.utils.remove_periods(jetR), sd_label))
    czg.SaveAs(output_filename)
    czg.Close()

  #---------------------------------------------------------------
  def plotRgProjection(self, hRM, hThetaG_JetPt, jetR, sd_label, zcut, beta, min_pt_det, max_pt_det):

    rebin_val_det = 10
    rebin_val_truth = 2
    
    # Get histogram of theta_g in data, for given pt-det cut
    if hThetaG_JetPt:
      hThetaG_JetPt.GetXaxis().SetRange(min_pt_det, max_pt_det)
      hThetaG_data = hThetaG_JetPt.ProjectionY()
      hThetaG_data.GetYaxis().SetTitle('#frac{dN}{d#theta_{g}}')
      hThetaG_data.SetMarkerStyle(21)
      hThetaG_data.SetMarkerSize(1)
      self.utils.scale_by_integral(hThetaG_data)
      hThetaG_data.Rebin(rebin_val_det)

    # Get histograms of theta_g in MC, for a given pt-det cut
    hRM.GetAxis(0).SetRange(min_pt_det, max_pt_det)
    
    hThetaG_det = hRM.Projection(2)
    hThetaG_det.SetName('hThetaG_det')
    hThetaG_det.SetLineColor(2)
    hThetaG_det.SetLineWidth(2)
    self.utils.scale_by_integral(hThetaG_det)
    hThetaG_det.Rebin(rebin_val_det)

    hThetaG_truth = hRM.Projection(3)
    hThetaG_truth.SetName('hThetaG_truth')
    hThetaG_truth.GetYaxis().SetTitle('#frac{dN}{d#theta_{g}}')
    hThetaG_truth.SetLineColor(4)
    hThetaG_truth.SetLineWidth(2)
    self.utils.scale_by_integral(hThetaG_truth)
    hThetaG_truth.Rebin(rebin_val_truth)

    # Draw histogram
    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()

    myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    myPad.SetLeftMargin(0.2)
    myPad.SetTopMargin(0.07)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.13)
    myPad.Draw()
    myPad.cd()
    
    hThetaG_truth.GetXaxis().SetRangeUser(0., 1.)
    hThetaG_truth.GetYaxis().SetTitleOffset(1.5)
    hThetaG_truth.SetMaximum(1.6*hThetaG_truth.GetMaximum())
    hThetaG_truth.SetMinimum(0.)

    #if hThetaG_JetPt:
    #hThetaG_data.Draw('E P')
    hThetaG_truth.Draw('hist E same')
    hThetaG_det.Draw('hist E same')
    
    leg = ROOT.TLegend(0.65,0.75,0.85,0.85, "")
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1)
    leg.SetTextSize(0.04)
    #if hThetaG_JetPt:
    #leg.AddEntry(hThetaG_data, "data", "P")
    leg.AddEntry(hThetaG_det, "MC det", "L")
    leg.AddEntry(hThetaG_truth, "MC truth", "L")
    leg.Draw("same")

    text = 'ALICE Simulation'
    textFit = ROOT.TLatex()
    textFit.SetNDC()
    textFit.DrawLatex(0.3,0.85,text)

    text = 'R = {}'.format(jetR)
    textFit = ROOT.TLatex()
    textFit.SetTextSize(0.04)
    textFit.SetNDC()
    textFit.DrawLatex(0.3,0.8,text)
    
    text = '#it{z}_{cut} = ' + str(zcut) + '    #beta = ' + str(beta)
    textFit = ROOT.TLatex()
    textFit.SetTextSize(0.04)
    textFit.SetNDC()
    textFit.DrawLatex(0.3,0.75,text)
    
    text = 'p_{T, ch jet}^{det} = ' + '{}-{} GeV/c'.format(min_pt_det, max_pt_det)
    textFit = ROOT.TLatex()
    textFit.SetTextSize(0.035)
    textFit.SetNDC()
    textFit.DrawLatex(0.3,0.7,text)
    
    output_filename = os.path.join(self.output_dir, 'hTheta_g_MC_R{}_{}_{}-{}.pdf'.format(self.utils.remove_periods(jetR), sd_label, min_pt_det, max_pt_det))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  def plotZgProjection(self, hRM, hZg_JetPt, jetR, sd_label, zcut, beta, min_pt_det, max_pt_det):
  
    rebin_val_det = 5
    rebin_val_truth = 1
    
    # Get histogram of theta_g in data, for given pt-det cut
    if hZg_JetPt:
      hZg_JetPt.GetXaxis().SetRange(min_pt_det, max_pt_det)
      hZg_data = hZg_JetPt.ProjectionY()
      hZg_data.GetYaxis().SetTitle('#frac{dN}{dz_{g}}')
      hZg_data.SetMarkerStyle(21)
      hZg_data.SetMarkerSize(1)
      self.utils.scale_by_integral(hZg_data)
      hZg_data.Rebin(rebin_val_det)
  
    # Get histograms of theta_g in MC, for a given pt-det cut
    hRM.GetAxis(0).SetRange(min_pt_det, max_pt_det)
    
    hZg_det = hRM.Projection(2)
    hZg_det.SetName('hZg_det')
    hZg_det.SetLineColor(2)
    hZg_det.SetLineWidth(2)
    self.utils.scale_by_integral(hZg_det)
    hZg_det.Rebin(rebin_val_det)
    
    hZg_truth = hRM.Projection(3)
    hZg_truth.SetName('hZg_truth')
    hZg_truth.GetYaxis().SetTitle('#frac{dN}{dz_{g}}')
    hZg_truth.SetLineColor(4)
    hZg_truth.SetLineWidth(2)
    self.utils.scale_by_integral(hZg_truth)
    hZg_truth.Rebin(rebin_val_truth)
    
    # Draw histogram
    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()
    
    myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    myPad.SetLeftMargin(0.2)
    myPad.SetTopMargin(0.07)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.13)
    myPad.Draw()
    myPad.cd()
    
    hZg_truth.GetXaxis().SetRangeUser(0., 1.)
    hZg_truth.GetYaxis().SetTitleOffset(1.5)
    hZg_truth.SetMaximum(1.6*hZg_truth.GetMaximum())
    hZg_truth.SetMinimum(0.)
    
    #if hZg_JetPt:
    #hZg_data.Draw('E P')
    hZg_truth.Draw('hist E same')
    hZg_det.Draw('hist E same')
    
    leg = ROOT.TLegend(0.65,0.75,0.85,0.85, "")
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1)
    leg.SetTextSize(0.04)
    #if hZg_JetPt:
    #leg.AddEntry(hZg_data, "data", "P")
    leg.AddEntry(hZg_det, "MC det", "L")
    leg.AddEntry(hZg_truth, "MC truth", "L")
    leg.Draw("same")
    
    text = 'ALICE Simulation'
    textFit = ROOT.TLatex()
    textFit.SetNDC()
    textFit.DrawLatex(0.3,0.85,text)
    
    text = 'R = {}'.format(jetR)
    textFit = ROOT.TLatex()
    textFit.SetTextSize(0.04)
    textFit.SetNDC()
    textFit.DrawLatex(0.3,0.8,text)
    
    text = '#it{z}_{cut} = ' + str(zcut) + '    #beta = ' + str(beta)
    textFit = ROOT.TLatex()
    textFit.SetTextSize(0.04)
    textFit.SetNDC()
    textFit.DrawLatex(0.3,0.75,text)
    
    text = 'p_{T, ch jet}^{det} = ' + '{}-{} GeV/c'.format(min_pt_det, max_pt_det)
    textFit = ROOT.TLatex()
    textFit.SetTextSize(0.035)
    textFit.SetNDC()
    textFit.DrawLatex(0.3,0.7,text)
    
    output_filename = os.path.join(self.output_dir, 'hZg_MC_R{}_{}_{}-{}.pdf'.format(self.utils.remove_periods(jetR), sd_label, min_pt_det, max_pt_det))
    c.SaveAs(output_filename)
    c.Close()
  
  #---------------------------------------------------------------
  def plotDeltaR(self, jetR):

    name = 'hDeltaR_All_R{}Scaled'.format(jetR)
    h = self.fMC.Get(name)
    
    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()
    ROOT.gPad.SetLeftMargin(0.15)
    c.SetLogz()

    h.GetXaxis().SetTitle('#it{p}_{T,det}^{ch jet}')
    h.GetYaxis().SetTitle('#DeltaR_{match}')
    
    x_max = 200.
    h.GetXaxis().SetRangeUser(0., x_max)
    h.GetYaxis().SetRangeUser(0., 1.)
    h.Draw('colz')
    
    deltaR_max = self.jet_matching_distance * jetR
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

    output_filename = os.path.join(self.output_dir, '{}.pdf'.format(self.utils.remove_periods(name)))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  def plotJER(self, jetR, sd_label):
    
    # (pt-det, pt-truth, theta_g-det, theta_g-truth)
    name = 'hResponse_JetPt_ThetaG_R{}_{}Scaled'.format(jetR, sd_label)
    hRM_4d = self.fMC.Get(name)
    hRM = hRM_4d.Projection(1,0)
    hRM.SetName('hResponse_JetPt_ThetaG_R{}_{}_Proj'.format(jetR, sd_label))
    
    # For each pT^gen, compute the standard deviation of the pT^det distribution
    
    # Get the pT^gen profile, with errors as standard deviation of pT^det distribution
    histPtGenProf = hRM.ProfileY("histPtGenProff", 1, -1, "s")
    
    # Create histo to be used to fill JER values
    nBins = 60
    histJER = ROOT.TH1D('histJER_R{}_{}'.format(jetR, sd_label), 'histJER_R{}_{}'.format(jetR, sd_label), nBins, 0., 300.) # same binning for pT^gen as in task
    
    # Loop through the bins, and fill the JER
    for bin in range(0,nBins+1):
      sigma = histPtGenProf.GetBinError(bin)
      pTgen = histPtGenProf.GetXaxis().GetBinCenter(bin)
      JER = sigma/pTgen
      histJER.SetBinContent(bin, JER)
    
    histJER.GetYaxis().SetTitle("#frac{#sigma(#it{p}_{T}^{gen})}{#it{p}_{T}^{gen}}")
    histJER.GetXaxis().SetTitle("#it{p}_{T}^{gen}")
    histJER.GetYaxis().SetRangeUser(-0.01, 0.5)
    histJER.GetXaxis().SetRangeUser(5., 100.)
    outputFilename = os.path.join(self.output_dir, 'histJER_R{}_{}.pdf'.format(self.utils.remove_periods(jetR), sd_label))
    histJER.SetMarkerStyle(21)
    histJER.SetMarkerColor(2)
    self.utils.plot_hist(histJER, outputFilename, "hist P")

  #---------------------------------------------------------------
  def plot_theta_resolution(self, jetR, sd_label):
    
    # (pt-det, pt-truth, theta_g-det, theta_g-truth)
    name = 'hResponse_JetPt_ThetaG_R{}_{}Scaled'.format(jetR, sd_label)
    hRM_4d = self.fMC.Get(name)
    hRM = hRM_4d.GetAxis(1).SetRangeUser(40, 60)
    hRM = hRM_4d.Projection(3,2)
    hRM.SetName('hResponse_JetPt_ThetaG_R{}_{}_Proj_theta'.format(jetR, sd_label))
    
    # For each pT^gen, compute the standard deviation of the pT^det distribution
    
    # Get the pT^gen profile, with errors as standard deviation of pT^det distribution
    histThetaGenProf = hRM.ProfileY("histThetaGenProff", 1, -1, "s")
    
    # Create histo to be used to fill JER values
    nBins = 20
    histThetaResolution = ROOT.TH1D('histThetaResolution_R{}_{}'.format(jetR, sd_label), 'histThetaResolution_R{}_{}'.format(jetR, sd_label), nBins, 0., 1.) # same binning for pT^gen as in task
    
    # Loop through the bins, and fill the JER
    for bin in range(0,nBins+1):
      sigma = histThetaGenProf.GetBinError(bin)
      theta_gen = histThetaGenProf.GetXaxis().GetBinCenter(bin)
      resolution = sigma/theta_gen
      histThetaResolution.SetBinContent(bin, resolution)
    
      histThetaResolution.GetYaxis().SetTitle("#frac{#sigma(#theta_{g}^{gen})}{#theta_{g}^{gen}}")
      histThetaResolution.GetXaxis().SetTitle("#theta_{g}^{gen}")
      histThetaResolution.GetYaxis().SetRangeUser(-0.01, 1.)
      histThetaResolution.GetXaxis().SetRangeUser(5., 100.)
      outputFilename = os.path.join(self.output_dir, 'histThetaResolution_R{}_{}.pdf'.format(self.utils.remove_periods(jetR), sd_label))
      histThetaResolution.SetMarkerStyle(21)
      histThetaResolution.SetMarkerColor(2)
      self.utils.plot_hist(histThetaResolution, outputFilename, "hist P")

  #---------------------------------------------------------------
  def plot_theta_residual(self, jetR, sd_label):

    name = 'hThetaGResidual_JetPt_R{}_{}Scaled'.format(jetR, sd_label)
    
    c_theta_residual = ROOT.TCanvas("cJES","cJES: hist",600,450)
    c_theta_residual.cd()
    c_theta_residual.SetBottomMargin(0.17)
    
    kGreen = 416
    kBlue  = 600
    kCyan  = 432
    kAzure   = 860
    kViolet  = 880
    kMagenta = 616
    kPink    = 900
    ColorArray = [kBlue-4, kAzure+7, kCyan-2, 7, kViolet-8, kBlue-6, kGreen+3]
    
    hThetaResidual1 = self.get_residual_proj(name, 'hThetaResidual1', 20, 40)
    hThetaResidual2 = self.get_residual_proj(name, 'hThetaResidual2', 40, 60)
    hThetaResidual3 = self.get_residual_proj(name, 'hThetaResidual3', 60, 80)
    
    hThetaResidual1.SetMarkerStyle(20)
    hThetaResidual2.SetMarkerStyle(21)
    hThetaResidual3.SetMarkerStyle(22)
    hThetaResidual1.SetMarkerColor(ColorArray[0])
    hThetaResidual2.SetMarkerColor(ColorArray[1])
    hThetaResidual3.SetMarkerColor(ColorArray[2])
    hThetaResidual1.SetLineColor(ColorArray[0])
    hThetaResidual2.SetLineColor(ColorArray[1])
    hThetaResidual3.SetLineColor(ColorArray[2])
    
    hThetaResidual1.GetXaxis().SetTitleOffset(1.6);
    hThetaResidual1.GetYaxis().SetTitle("Probability density")
    
    hThetaResidual1.GetYaxis().SetRangeUser(0, 20.)
    hThetaResidual1.DrawCopy("P E")
    hThetaResidual2.DrawCopy("P E same")
    hThetaResidual3.DrawCopy("P E same")
    
    leg = ROOT.TLegend(0.55,0.55,0.88,0.85, "")
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1)
    leg.SetTextSize(0.04)
    
    leg.AddEntry(hThetaResidual1, "#it{p}_{T}^{gen} = 20-40 GeV", "P")
    leg.AddEntry(hThetaResidual2, "#it{p}_{T}^{gen} = 40-60 GeV", "P")
    leg.AddEntry(hThetaResidual3, "#it{p}_{T}^{gen} = 60-80 GeV", "P")
    leg.Draw("same")
    
    outputFilename = os.path.join(self.output_dir, 'histThetaResidual_R{}_{}.pdf'.format(self.utils.remove_periods(jetR), sd_label))
    c_theta_residual.SaveAs(outputFilename)
    c_theta_residual.Close()

  #---------------------------------------------------------------
  # Get JES shift distribution for a fixed pT-gen
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
  
  #---------------------------------------------------------------
  def plot_zg_residual(self, jetR, sd_label):
    
    name = 'hZgResidual_JetPt_R{}_{}Scaled'.format(jetR, sd_label)
    
    c_zg_residual = ROOT.TCanvas("cJES","cJES: hist",600,450)
    c_zg_residual.cd()
    c_zg_residual.SetBottomMargin(0.17)
    
    kGreen = 416
    kBlue  = 600
    kCyan  = 432
    kAzure   = 860
    kViolet  = 880
    kMagenta = 616
    kPink    = 900
    ColorArray = [kBlue-4, kAzure+7, kCyan-2, 7, kViolet-8, kBlue-6, kGreen+3]
    
    hZgResidual1 = self.get_residual_proj(name, 'hZgResidual1', 20, 40)
    hZgResidual2 = self.get_residual_proj(name, 'hZgResidual2', 40, 60)
    hZgResidual3 = self.get_residual_proj(name, 'hZgResidual3', 60, 80)
    
    hZgResidual1.SetMarkerStyle(20)
    hZgResidual2.SetMarkerStyle(21)
    hZgResidual3.SetMarkerStyle(22)
    hZgResidual1.SetMarkerColor(ColorArray[0])
    hZgResidual2.SetMarkerColor(ColorArray[1])
    hZgResidual3.SetMarkerColor(ColorArray[2])
    hZgResidual1.SetLineColor(ColorArray[0])
    hZgResidual2.SetLineColor(ColorArray[1])
    hZgResidual3.SetLineColor(ColorArray[2])
    
    hZgResidual1.GetXaxis().SetTitleOffset(1.6);
    hZgResidual1.GetYaxis().SetTitle("Probability density")
    
    hZgResidual1.GetYaxis().SetRangeUser(0, 20.)
    hZgResidual1.DrawCopy("P E")
    hZgResidual2.DrawCopy("P E same")
    hZgResidual3.DrawCopy("P E same")
    
    leg = ROOT.TLegend(0.55,0.55,0.88,0.85, "")
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1)
    leg.SetTextSize(0.04)
    
    leg.AddEntry(hZgResidual1, "#it{p}_{T}^{gen} = 20-40 GeV", "P")
    leg.AddEntry(hZgResidual2, "#it{p}_{T}^{gen} = 40-60 GeV", "P")
    leg.AddEntry(hZgResidual3, "#it{p}_{T}^{gen} = 60-80 GeV", "P")
    leg.Draw("same")
    
    outputFilename = os.path.join(self.output_dir, 'histZgResidual_R{}_{}.pdf'.format(self.utils.remove_periods(jetR), sd_label))
    c_zg_residual.SaveAs(outputFilename)
    c_zg_residual.Close()

  #---------------------------------------------------------------
  def plotJES(self, jetR):
    
    name = 'hJES_R{}Scaled'.format(jetR)
    histDeltaJES = self.fMC.Get(name)
    histDeltaJES.GetXaxis().SetTitle("#it{p}_{T}^{gen}")
    histDeltaJES.GetYaxis().SetTitle("#frac{#it{p}_{T}^{det} - #it{p}_{T}^{gen}}{#it{p}_{T}^{gen}}")
    histDeltaJES.GetXaxis().SetRangeUser(0., 200.)
    outputFilename = os.path.join(self.output_dir, "histDeltaJES_R{}.pdf".format(self.utils.remove_periods(jetR)))
    #self.utils.plot_hist(histDeltaJES, outputFilename, "colz", False, True)
    
    histDeltaJES.RebinX(4)
    histDeltaJESprof = histDeltaJES.ProfileX()
    histDeltaJESprof.GetXaxis().SetRangeUser(0., 200.)
    histDeltaJESprof.GetYaxis().SetRangeUser(-0.8, 0.2)
    histDeltaJESprof.GetYaxis().SetTitle("#frac{#it{p}_{T}^{det} - #it{p}_{T}^{gen}}{#it{p}_{T}^{gen}}")
    outputFilename = os.path.join(self.output_dir, "histDeltaJESprof_R{}.pdf".format(self.utils.remove_periods(jetR)))
    histDeltaJESprof.SetMarkerStyle(21)
    histDeltaJESprof.SetMarkerColor(4)
    self.utils.plot_hist(histDeltaJESprof, outputFilename, "P")

  #---------------------------------------------------------------
  # Plot JES shift distribution for various fixed pT-gen
  def plotJESproj(self, jetR):
    
    name = 'hJES_R{}'.format(jetR)
    
    cJES = ROOT.TCanvas("cJES","cJES: hist",600,450)
    cJES.cd()
    cJES.SetBottomMargin(0.17)
    
    kGreen = 416
    kBlue  = 600
    kCyan  = 432
    kAzure   = 860
    kViolet  = 880
    kMagenta = 616
    kPink    = 900
    ColorArray = [kBlue-4, kAzure+7, kCyan-2,7, kViolet-8, kBlue-6, kGreen+3]
    
    hJESProj1 = self.getJESshiftProj(name, "hJESproj1", 20, 30)
    hJESProj2 = self.getJESshiftProj(name, "hJESproj2", 50, 70)
    hJESProj3 = self.getJESshiftProj(name, "hJESproj3", 100, 120)
    hJESProj1b = self.getJESshiftProj(name, "hJESproj1b", 30, 40)
    hJESProj2b = self.getJESshiftProj(name, "hJESproj2b", 40, 50)
    hJESProj3b = self.getJESshiftProj(name, "hJESproj3b", 70, 90)
    
    hJESProj1.SetMarkerStyle(20)
    hJESProj1b.SetMarkerStyle(20)
    hJESProj2.SetMarkerStyle(20)
    hJESProj2b.SetMarkerStyle(20)
    hJESProj3.SetMarkerStyle(20)
    hJESProj3b.SetMarkerStyle(20)
    hJESProj1.SetMarkerColor(ColorArray[0])
    hJESProj2.SetMarkerColor(ColorArray[1])
    hJESProj3.SetMarkerColor(ColorArray[2])
    hJESProj1b.SetMarkerColor(ColorArray[3])
    hJESProj2b.SetMarkerColor(ColorArray[4])
    hJESProj3b.SetMarkerColor(ColorArray[5])
    hJESProj1.SetLineColor(ColorArray[0])
    hJESProj2.SetLineColor(ColorArray[1])
    hJESProj3.SetLineColor(ColorArray[2])
    hJESProj1b.SetLineColor(ColorArray[3])
    hJESProj2b.SetLineColor(ColorArray[4])
    hJESProj3b.SetLineColor(ColorArray[5])
    
    hJESProj1.GetXaxis().SetTitleOffset(1.6);
    hJESProj1.GetYaxis().SetTitle("Probability density")
    
    hJESProj1.GetYaxis().SetRangeUser(0, 10.4)
    hJESProj1.DrawCopy("P E")
    hJESProj2.DrawCopy("P E same")
    hJESProj3.DrawCopy("P E same")
    hJESProj1b.DrawCopy("P E same")
    hJESProj2b.DrawCopy("P E same")
    hJESProj3b.DrawCopy("P E same")
    
    leg = ROOT.TLegend(0.55,0.55,0.88,0.85, "")
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1)
    leg.SetTextSize(0.04)
    
    leg.AddEntry(hJESProj1, "#it{p}_{T}^{gen} = 20-30 GeV", "P")
    leg.AddEntry(hJESProj1b, "#it{p}_{T}^{gen} = 30-40 GeV", "P")
    leg.AddEntry(hJESProj2b, "#it{p}_{T}^{gen} = 40-50 GeV", "P")
    leg.AddEntry(hJESProj2, "#it{p}_{T}^{gen} = 50-70 GeV", "P")
    leg.AddEntry(hJESProj3b, "#it{p}_{T}^{gen} = 70-90 GeV", "P")
    leg.AddEntry(hJESProj3, "#it{p}_{T}^{gen} = 100-120 GeV", "P")
    leg.Draw("same")
    
    outputFilename = os.path.join(self.output_dir, "histDeltaJESproj_R{}.pdf".format(self.utils.remove_periods(jetR)))
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
  def plotJetRecoEfficiency(self, jetR, sd_label):
    
    # For each pT^gen, compute the fraction of matched pT^gen
    
    # First, get the pT^gen spectrum
    name = 'hJetPt_Truth_R{}Scaled'.format(jetR)
    histPtGen = self.fMC.Get(name)
    histPtGen.Rebin(5)
    histPtGen.Scale(1/3.) # TEMP fix since I accidentally filled 3x
    
    # Then, get the pT^gen spectrum for matched jets
    name = 'hResponse_JetPt_ThetaG_R{}_{}Scaled'.format(jetR, sd_label)
    hRM_4d = self.fMC.Get(name)
    hRM = hRM_4d.Projection(1,0)
    hRM.SetName('hResponse_JetPt_ThetaG_R{}_{}_Proj'.format(jetR, sd_label))
    histPtGenMatched = hRM.ProjectionY("_py",1,hRM.GetNbinsX()) #avoid under and overflow bins
    histPtGenMatched.SetName('histPtGenMatched_R{}_{}'.format(jetR, sd_label))
    
    # Compute the ratio
    histEfficiency = histPtGenMatched.Clone()
    histEfficiency.SetName("histEfficiency_{}".format(name))
    histEfficiency.Divide(histPtGenMatched, histPtGen, 1., 1., "B")
    
    histEfficiency.GetXaxis().SetTitle("#it{p}_{T}^{gen}")
    histEfficiency.GetYaxis().SetTitle("Efficiency")
    histEfficiency.GetXaxis().SetRangeUser(0., 100.)
    histEfficiency.GetYaxis().SetRangeUser(0., 1.2)
    histEfficiency.SetMarkerStyle(21)
    histEfficiency.SetMarkerColor(1)
    outputFilename = os.path.join(self.output_dir, '{}_R{}.pdf'.format(self.utils.remove_periods(histEfficiency.GetName()), self.utils.remove_periods(jetR)))
    self.utils.plot_hist(histEfficiency, outputFilename)

#----------------------------------------------------------------------
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Plot analysis histograms')
  parser.add_argument('-r', '--mcFile', action='store',
                      type=str, metavar='mcFile',
                      default='AnalysisResults.root',
                      help='Path of MC ROOT file')
  parser.add_argument('-d', '--dataFile', action='store',
                      type=str, metavar='dataFile',
                      default='AnalysisResults.root',
                      help='Path of data ROOT file')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='analysis_config.yaml',
                      help="Path of config file for jetscape analysis")
  parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./mc_performance_output',
                      help='Output directory for QA plots to be written to')
  
  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('mcFile: \'{0}\''.format(args.mcFile))
  print('dataFile: \'{0}\''.format(args.dataFile))
  print('configFile: \'{0}\''.format(args.configFile))
  print('ouputDir: \'{0}\"'.format(args.outputDir))
  print('----------------------------------------------------------------')
  
  # If invalid mcFile is given, exit
  if not os.path.exists(args.mcFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.mcFile))
    sys.exit(0)

  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = plot_rg_performance(mc_file = args.mcFile, data_file = args.dataFile, config_file = args.configFile, output_dir = args.outputDir)
  analysis.plot_rg_performance()
