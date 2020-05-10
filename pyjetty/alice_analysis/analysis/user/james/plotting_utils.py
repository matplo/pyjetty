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
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs

################################################################
class PlottingUtils(analysis_utils_obs.AnalysisUtils_Obs):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, output_dir = '.', config_file = '', R_max = None, thermal = False, **kwargs):
    super(PlottingUtils, self).__init__(**kwargs)
    
    self.output_dir = output_dir
    self.R_max = R_max
    self.thermal = thermal
    
    # Read config file
    with open(config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    
    self.observable = config['analysis_observable']
    self.is_pp = not 'constituent_subtractor' in config
    self.figure_approval_status = config['figure_approval_status']

    main_data = config['main_data']
    main_response = config['main_response']
    if main_data:
      self.fData = ROOT.TFile(main_data, 'READ')
    else:
      self.fData = None
    self.fMC = ROOT.TFile(main_response, 'READ')
    
    if self.R_max:
      self.suffix = '_Rmax{}'.format(self.R_max)
    else:
      self.suffix = ''
    
    self.ColorArray = [ROOT.kBlue-4, ROOT.kAzure+7, ROOT.kCyan-2, ROOT.kViolet-8,
                       ROOT.kBlue-6, ROOT.kGreen+3, ROOT.kPink-4, ROOT.kRed-4,
                       ROOT.kOrange-3]
    self.MarkerArray = [20, 21, 22, 23, 33, 34, 24, 25, 26, 32]
    
    print(self)

  #---------------------------------------------------------------
  def plot_DeltaR(self, jetR, jet_matching_distance):

    if self.is_pp:
      name = 'hDeltaR_All_R{}{}Scaled'.format(jetR, self.suffix)
    else:
      name = 'hDeltaR_combined_ppdet_R{}{}Scaled'.format(jetR, self.suffix)
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

    output_filename = os.path.join(self.output_dir, 'jet/{}.pdf'.format(self.remove_periods(name)))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  def plot_JES(self, jetR):
    
    name = 'hJES_R{}{}Scaled'.format(jetR, self.suffix)
    histDeltaJES = self.fMC.Get(name)
    histDeltaJES.GetXaxis().SetTitle("#it{p}_{T}^{gen}")
    histDeltaJES.GetYaxis().SetTitle("#frac{#it{p}_{T}^{det} - #it{p}_{T}^{gen}}{#it{p}_{T}^{gen}}")
    histDeltaJES.GetXaxis().SetRangeUser(0., 200.)
    outputFilename = os.path.join(self.output_dir, "histDeltaJES_R{}.pdf".format(self.remove_periods(jetR)))
    #self.utils.plot_hist(histDeltaJES, outputFilename, "colz", False, True)
    
    histDeltaJES.RebinX(4)
    histDeltaJESprof = histDeltaJES.ProfileX()
    histDeltaJESprof.GetXaxis().SetRangeUser(0., 200.)
    histDeltaJESprof.GetYaxis().SetRangeUser(-0.5, 0.5)
    histDeltaJESprof.GetYaxis().SetTitle("#frac{#it{p}_{T}^{det} - #it{p}_{T}^{gen}}{#it{p}_{T}^{gen}}")
    outputFilename = os.path.join(self.output_dir, "jet/histDeltaJESprof_R{}.pdf".format(self.remove_periods(jetR)))
    histDeltaJESprof.SetMarkerStyle(21)
    histDeltaJESprof.SetMarkerColor(4)
    histDeltaJESprof.SetMarkerSize(3)
    self.plot_hist(histDeltaJESprof, outputFilename, "P")

  #---------------------------------------------------------------
  # Plot JES shift distribution for various fixed pT-gen
  def plot_JES_proj(self, jetR, pt_bins):
    
    name = 'hJES_R{}{}Scaled'.format(jetR, self.suffix)
    
    cJES = ROOT.TCanvas('cJES','cJES: hist',600,450)
    cJES.cd()
    cJES.SetBottomMargin(0.2)
    
    leg = ROOT.TLegend(0.55,0.55,0.88,0.85, '')
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1)
    leg.SetTextSize(0.04)

    # Loop through pt slices, and plot final result for each 1D theta_g distribution
    for i in range(0, len(pt_bins) - 1):
      min_pt_truth = pt_bins[i]
      max_pt_truth = pt_bins[i+1]
      
      hJESProj = self.getJESshiftProj(name, 'hJESproj{}'.format(i), min_pt_truth, max_pt_truth)
      hJESProj.SetMarkerStyle(20)
      hJESProj.SetMarkerColor(self.ColorArray[i])
      hJESProj.SetLineColor(self.ColorArray[i])
      
      if i == 0:
      
        hJESProj.GetXaxis().SetTitleOffset(1.6);
        hJESProj.GetYaxis().SetTitle('Probability density')
        hJESProj.GetXaxis().SetTitle('#frac{#it{p}_{T}^{det} - #it{p}_{T}^{gen}}{#it{p}_{T}^{gen}}')
        
        hJESProj.GetYaxis().SetRangeUser(0, 1.3*hJESProj.GetMaximum())
        hJESProj.DrawCopy('P E')
        
      else:
      
        hJESProj.DrawCopy('P E same')
    
      leg.AddEntry(hJESProj, '#it{{p}}_{{T}}^{{gen}} = {}-{} GeV'.format(min_pt_truth, max_pt_truth), 'P')
      
    leg.Draw('same')
    
    outputFilename = os.path.join(self.output_dir, 'jet/histDeltaJESproj_R{}.pdf'.format(self.remove_periods(jetR)))
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
    name = 'hResponse_JetPt_{}_R{}_{}{}Scaled'.format(self.observable, jetR, obs_label, self.suffix)
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
    outputFilename = os.path.join(self.output_dir, 'jet/histJER_R{}.pdf'.format(self.remove_periods(jetR)))
    histJER.SetMarkerStyle(21)
    histJER.SetMarkerColor(2)
    histJER.SetMarkerSize(3)
    self.plot_hist(histJER, outputFilename, 'hist P')
  
  #---------------------------------------------------------------
  def plot_jet_reco_efficiency(self, jetR, obs_label):
    
    # For each pT^gen, compute the fraction of matched pT^gen
    
    # First, get the pT^gen spectrum
    name = 'h_{}_JetPt_Truth_R{}_{}Scaled'.format(self.observable, jetR, obs_label)
    histPtGen = self.fMC.Get(name).ProjectionX()
    
    # Then, get the pT^gen spectrum for matched jets
    name = 'hResponse_JetPt_{}_R{}_{}{}Scaled'.format(self.observable, jetR, obs_label, self.suffix)
    hRM_4d = self.fMC.Get(name)
    hRM = hRM_4d.Projection(1,0)
    hRM.SetName('hResponse_JetPt_{}_R{}_{}_Proj'.format(self.observable, jetR, obs_label))
    histPtGenMatched = hRM.ProjectionY("_py",1,hRM.GetNbinsX()) #avoid under and overflow bins
    histPtGenMatched.SetName('histPtGenMatched_R{}_{}'.format(jetR, obs_label))
    
    # Compute the ratio
    histEfficiency = histPtGenMatched.Clone()
    histEfficiency.SetName('histEfficiency_{}'.format(name))
    histEfficiency.Divide(histPtGenMatched, histPtGen, 1., 1., 'B')
    
    histEfficiency.GetXaxis().SetTitle('#it{p}_{T}^{gen}')
    histEfficiency.GetYaxis().SetTitle('Efficiency')
    histEfficiency.GetXaxis().SetRangeUser(0., 100.)
    histEfficiency.GetYaxis().SetRangeUser(0., 1.2)
    histEfficiency.SetMarkerStyle(21)
    histEfficiency.SetMarkerColor(1)
    histEfficiency.SetMarkerSize(3)
    outputFilename = os.path.join(self.output_dir, 'jet/hJetRecoEfficiency_R{}.pdf'.format(self.remove_periods(jetR)))
    self.plot_hist(histEfficiency, outputFilename)
      
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
    
    # (pt-det, pt-truth, theta_g-det, theta_g-truth)
    name = 'hResponse_JetPt_{}_R{}_{}{}Scaled'.format(self.observable, jetR, obs_label, self.suffix)
    hRM_4d = self.fMC.Get(name)
    
    h_list = [] # Store hists in a list, since otherwise it seems I lose the marker information
                # (removed from memory?)
    
    for i in range(0, len(pt_bins) - 1):
    
      min_pt_truth = pt_bins[i]
      max_pt_truth = pt_bins[i+1]
      
      hRM_4d_clone = hRM_4d.Clone()
      hRM_4d_clone.SetName('{}_{}'.format(hRM_4d_clone.GetName(), i))
      hResolution = self.get_resolution(hRM_4d_clone, jetR, obs_label, min_pt_truth, max_pt_truth, 'hResolution_{}'.format(i))
 
      hResolution.SetMarkerColor(self.ColorArray[i])
      hResolution.SetMarkerStyle(21)
      hResolution.SetLineColor(self.ColorArray[i])
      hResolution.DrawCopy('P same')
      myLegend.AddEntry(hResolution, '#it{{p}}_{{T}}^{{gen}} = {}-{} GeV'.format(min_pt_truth, max_pt_truth), 'P')
      h_list.append(hResolution)

    myLegend.Draw('same')

    outputFilename = os.path.join(self.output_dir, 'resolution/hResolution_R{}_{}.pdf'.format(self.remove_periods(jetR), obs_label))
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
  def plot_obs_residual_pt(self, jetR, obs_label, xtitle, pt_bins):

    if self.observable == 'subjet_z' or self.observable == 'jet_axis':
      return
    else:
      name = 'hResidual_JetPt_{}_R{}_{}{}Scaled'.format(self.observable, jetR, obs_label, self.suffix)
    
    c_residual = ROOT.TCanvas('c','c: hist',600,450)
    c_residual.cd()
    c_residual.SetBottomMargin(0.2)
    
    leg = ROOT.TLegend(0.55,0.55,0.88,0.85, '')
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1)
    leg.SetTextSize(0.04)
    
    # Loop through pt slices, and plot final residual for each 1D distribution
    for i in range(0, len(pt_bins) - 1):
      min_pt_truth = pt_bins[i]
      max_pt_truth = pt_bins[i+1]
      
      hResidual = self.get_residual_proj(name, 'hResidual{}'.format(i), min_pt_truth, max_pt_truth, option='pt')
      hResidual.SetMarkerStyle(self.MarkerArray[i])
      hResidual.SetMarkerColor(self.ColorArray[i])
      hResidual.SetLineColor(self.ColorArray[i])
      
      if i == 0:
      
        hResidual.GetXaxis().SetTitleOffset(1.6);
        hResidual.GetYaxis().SetTitle('Probability density')
        hResidual.GetYaxis().SetRangeUser(0, 2.*hResidual.GetMaximum())
        hResidual.DrawCopy('P E')
        
      else:

        hResidual.DrawCopy('P E same')
    
      leg.AddEntry(hResidual, '#it{{p}}_{{T}}^{{gen}} = {}-{} GeV'.format(min_pt_truth, max_pt_truth), 'P')

    leg.Draw('same')
    
    outputFilename = os.path.join(self.output_dir, 'residual_pt/hResidual_R{}_{}.pdf'.format(self.remove_periods(jetR), obs_label))
    c_residual.SaveAs(outputFilename)
    c_residual.Close()

  #---------------------------------------------------------------
  def plot_obs_residual_obs(self, jetR, obs_label, xtitle):

    name = 'hResidual_JetPt_{}_R{}_{}{}Scaled'.format(self.observable, jetR, obs_label, self.suffix)
    
    c_residual = ROOT.TCanvas('c','c: hist',600,450)
    c_residual.cd()
    c_residual.SetBottomMargin(0.2)
    
    leg = ROOT.TLegend(0.55,0.55,0.88,0.85, '')
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1)
    leg.SetTextSize(0.035)
    
    # Loop through pt slices, and plot final residual for each 1D distribution
    min_pt = 80
    max_pt = 100
    if 'theta' in xtitle:
      obs_true_list = [0., 0.05, 0.1, 0.2, 0.5, 1.]
    else:
      obs_true_list = [0.2, 0.3, 0.4, 0.5]
    for i in range(0, len(obs_true_list) - 1):
      min_obs_truth = obs_true_list[i]
      max_obs_truth = obs_true_list[i+1]
      
      hResidual = self.get_residual_proj(name, 'hResidual{}'.format(i), min_obs_truth, max_obs_truth, option='obs', min_pt=min_pt, max_pt=max_pt)
      hResidual.SetMarkerStyle(self.MarkerArray[i])
      hResidual.SetMarkerColor(self.ColorArray[i])
      hResidual.SetLineColor(self.ColorArray[i])

      if i == 0:
      
        hResidual.GetXaxis().SetTitleOffset(1.6);
        hResidual.GetYaxis().SetTitle('Probability density')
        hResidual.GetYaxis().SetRangeUser(0, 2.*hResidual.GetMaximum())
        hResidual.DrawCopy('P E')
        
      else:

        hResidual.DrawCopy('P E same')

      leg.AddEntry(hResidual, '{} = {}-{}'.format('{}^{{{}}}'.format(xtitle, 'truth'), min_obs_truth, max_obs_truth), 'P')

    leg.Draw('same')
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = '#it{{p}}_{{T}}^{{gen}} = {}-{} GeV'.format(min_pt, max_pt)
    text_latex.DrawLatex(0.2, 0.8, text)
    
    outputFilename = os.path.join(self.output_dir, 'residual_obs/hResidual_R{}_{}.pdf'.format(self.remove_periods(jetR), obs_label))
    c_residual.SaveAs(outputFilename)
    c_residual.Close()

  #---------------------------------------------------------------
  # Get residual for a fixed pT-gen
  def get_residual_proj(self, name, label, min, max, option='pt', min_pt=80., max_pt=100.):
    
    h_residual_pt = self.fMC.Get(name)
    h_residual_pt.SetName('{}_{}'.format(h_residual_pt.GetName(), label))
    
    if option == 'pt':
      h_residual_pt.GetXaxis().SetRangeUser(min, max)
    elif option == 'obs':
      h_residual_pt.GetXaxis().SetRangeUser(min_pt, max_pt)
      h_residual_pt.GetYaxis().SetRangeUser(min, max)
    h_residual_pt.GetZaxis().SetRangeUser(-0.5, 0.5)
    h = h_residual_pt.Project3D('z')
    
    integral = h.Integral()
    if integral > 0:
      h.Scale(1./integral, 'width')
    
    return h

  #---------------------------------------------------------------
  def plot_obs_projections(self, jetR, obs_label, obs_setting, grooming_setting, xtitle, pt_bins):

    if not self.fData:
      return
      
    # (pt-det, pt-truth, obs-det, obs-truth)
    name = 'hResponse_JetPt_{}_R{}_{}{}Scaled'.format(self.observable, jetR, obs_label, self.suffix)
    hRM_obs = self.fMC.Get(name)
    if hRM_obs.GetSumw2() is 0:
      hRM_obs.Sumw2()
    
    if self.thermal:
      name = 'h_{}_JetPt_R{}_{}{}Scaled'.format(self.observable, jetR, obs_label, self.suffix)
    else:
      name = 'h_{}_JetPt_R{}_{}{}'.format(self.observable, jetR, obs_label, self.suffix)
    hObs_JetPt = self.fData.Get(name)
    if hObs_JetPt.GetSumw2() is 0:
      hObs_JetPt.Sumw2()

    # Plot 2D statistics in data
    self.plot2D_obs_statistics(hObs_JetPt.Clone(), jetR, obs_label)
      
    # Loop through pt slices, and plot:
    #   (a) MC-det and MC-truth 1D distributions, for fixed pt-truth
    #   (b) MC-det and data 1D distributions, for fixed pt-det
    for i in range(0, len(pt_bins) - 1):
      min_pt_truth = pt_bins[i]
      max_pt_truth = pt_bins[i+1]
      
      self.plot_obs_projection(hRM_obs, hObs_JetPt, jetR, obs_label, obs_setting, grooming_setting, xtitle, min_pt_truth, max_pt_truth, option='truth')
      self.plot_obs_projection(hRM_obs, hObs_JetPt, jetR, obs_label, obs_setting, grooming_setting, xtitle, min_pt_truth, max_pt_truth, option='det')

  #---------------------------------------------------------------
  def plot_obs_truth(self, jetR, obs_label, obs_setting, grooming_setting,
                     xtitle, pt_bins):
                     
    name = 'h_{}_JetPt_Truth_R{}_{}Scaled'.format(self.observable, jetR, obs_label)
    h2D = self.fMC.Get(name)
    
    for i in range(0, len(pt_bins) - 1):
      min_pt_truth = pt_bins[i]
      max_pt_truth = pt_bins[i+1]
      
      self.plot_truth_projection(h2D, jetR, obs_label, obs_setting, grooming_setting, xtitle, min_pt_truth, max_pt_truth)
    
  #---------------------------------------------------------------
  def plot_truth_projection(self, h2D, jetR, obs_label, obs_setting,
                            grooming_setting, xtitle, min_pt, max_pt):
    
    ytitle = '#frac{{1}}{{N}} #frac{{dN}}{{d{}}}'.format(xtitle)

    # Get histogram of theta_g in data, for given pt-det cut
    h2D.GetXaxis().SetRangeUser(min_pt, max_pt)
    hObs_truth = h2D.ProjectionY()
    hObs_truth.SetMarkerStyle(21)
    hObs_truth.SetMarkerSize(1)
    self.scale_by_integral(hObs_truth)
    hObs_truth.Rebin(5)
    hObs_truth.Scale(1., 'width')
    if 'sd' in grooming_setting:
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
    leg.AddEntry(hObs_truth, "MC truth", "L")
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

    text = str(min_pt) + ' < #it{p}_{T, ch jet}^{truth} < ' + str(max_pt) + ' GeV/#it{c}'
    text_latex.DrawLatex(0.3, 0.73, text)

    text = '#it{R} = ' + str(jetR) + '   | #eta_{jet}| < 0.5'
    text_latex.DrawLatex(0.3, 0.67, text)
    
    subobs_label = self.formatted_subobs_label(self.observable)
    delta = 0.
    if subobs_label:
      text = '{} = {}'.format(subobs_label, obs_setting)
      text_latex.DrawLatex(0.3, 0.61, text)
      delta = 0.07
      
    if grooming_setting:
      text = self.formatted_grooming_label(grooming_setting)
      text_latex.DrawLatex(0.3, 0.61-delta, text)

    output_filename = os.path.join(self.output_dir, 'truth/h_{}_MC_R{}_{}_{}-{}.pdf'.format(self.observable, self.remove_periods(jetR), obs_label, min_pt, max_pt))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  def plot2D_obs_statistics(self, hObs_JetPt, jetR, obs_label):

    c = ROOT.TCanvas('c','c: hist',600,450)
    c.cd()
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.15)
    
    hObs_JetPt.SetMarkerSize(0.5)
    hObs_JetPt.GetYaxis().SetRangeUser(0, 1.)
    hObs_JetPt.GetXaxis().SetRangeUser(0, 100)
    hObs_JetPt.RebinX(5)
    hObs_JetPt.RebinY(5)
    hObs_JetPt.Draw('text colz')

    output_filename = os.path.join(self.output_dir, 'data/h2D_{}_statistics_R{}_{}.pdf'.format(self.observable, self.remove_periods(jetR), obs_label))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  # If option='truth', plot MC-truth and MC-det projections for fixed pt-true
  # If option='det', plot data and MC-det projections for fixed pt-det
  def plot_obs_projection(self, hRM, hObs_JetPt, jetR, obs_label, obs_setting, grooming_setting, xtitle, min_pt, max_pt, option='truth'):

    ytitle = '#frac{{1}}{{N}} #frac{{dN}}{{d{}}}'.format(xtitle)
    
    if self.observable == 'theta_g':
      rebin_val_mcdet = 5
      rebin_val_mctruth = 5
      rebin_val_data = 5
    elif self.observable == 'zg':
      rebin_val_mcdet = 5
      rebin_val_mctruth = 2
      rebin_val_data = 10
    elif self.observable == 'subjet_z':
      rebin_val_mcdet = 2
      rebin_val_mctruth = 1
      rebin_val_data = 2
    elif self.observable == 'jet_axis':
      rebin_val_mcdet = 2
      rebin_val_mctruth = 1
      rebin_val_data = 5
      
    # Get RM, for a given pt cut
    if option == 'det':
      hRM.GetAxis(0).SetRangeUser(min_pt, max_pt)
    if option == 'truth':
      hRM.GetAxis(1).SetRangeUser(min_pt, max_pt)
    
    # Get histogram of observable at MC-det from RM
    hObs_det = hRM.Projection(2)
    hObs_det.SetName('hObs_det_{}'.format(obs_label))
    hObs_det.GetYaxis().SetTitle(ytitle)
    hObs_det.SetLineColor(2)
    hObs_det.SetLineWidth(2)
    self.scale_by_integral(hObs_det)
    hObs_det.Rebin(rebin_val_mcdet)
    hObs_det.Scale(1., 'width')
    if 'sd' in grooming_setting:
      hObs_det.GetXaxis().SetRange(0, hObs_det.GetNbinsX())

    # Get histogram of observable at MC-truth from RM
    if option == 'truth':
      hObs_truth = hRM.Projection(3)
      hObs_truth.SetName('hObs_truth_{}'.format(obs_label))
      hObs_truth.SetLineColor(4)
      hObs_truth.SetLineWidth(2)
      self.scale_by_integral(hObs_truth)
      hObs_truth.Rebin(rebin_val_mctruth)
      hObs_truth.Scale(1., 'width')
      if 'sd' in grooming_setting:
        hObs_truth.GetXaxis().SetRange(0, hObs_truth.GetNbinsX())
      
    # Get histogram of theta_g in data, for given pt-det cut
    if option == 'det':
      hObs_JetPt.GetXaxis().SetRangeUser(min_pt, max_pt)
      hObs_data = hObs_JetPt.ProjectionY()
      hObs_data.SetMarkerStyle(21)
      hObs_data.SetMarkerSize(1)
      self.scale_by_integral(hObs_data)
      hObs_data.Rebin(rebin_val_data)
      hObs_data.Scale(1., 'width')
      if 'sd' in grooming_setting:
        hObs_data.GetXaxis().SetRange(0, hObs_data.GetNbinsX())

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
    
    hObs_det.GetYaxis().SetTitleOffset(1.5)
    hObs_det.SetMaximum(2.5*hObs_det.GetMaximum())
    hObs_det.SetMinimum(0.)

    hObs_det.Draw('hist')
    leg.AddEntry(hObs_det, "MC det", "L")
    if option == 'truth':
      hObs_truth.Draw('hist same')
      leg.AddEntry(hObs_truth, "MC truth", "L")
    elif option == 'det':
      hObs_data.Draw('hist E same')
      leg.AddEntry(hObs_data, "data", "PE")

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

    if option == 'truth':
      text = str(min_pt) + ' < #it{p}_{T, ch jet}^{truth} < ' + str(max_pt) + ' GeV/#it{c}'
    elif option == 'det':
      text = str(min_pt) + ' < #it{p}_{T, ch jet}^{det} < ' + str(max_pt) + ' GeV/#it{c}'
    text_latex.DrawLatex(0.3, 0.73, text)

    text = '#it{R} = ' + str(jetR) + '   | #eta_{jet}| < 0.5'
    text_latex.DrawLatex(0.3, 0.67, text)
    
    subobs_label = self.formatted_subobs_label(self.observable)
    delta = 0.
    if subobs_label:
      text = '{} = {}'.format(subobs_label, obs_setting)
      text_latex.DrawLatex(0.3, 0.61, text)
      delta = 0.07
      
    if grooming_setting:
      text = self.formatted_grooming_label(grooming_setting)
      text_latex.DrawLatex(0.3, 0.61-delta, text)

    output_filename = os.path.join(self.output_dir, 'mc_projections_{}/h_{}_MC_R{}_{}_{}-{}.pdf'.format(option, self.observable, self.remove_periods(jetR), obs_label, min_pt, max_pt))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  def plot_lund_plane(self, jetR, obs_label, grooming_setting):

    name = 'hLundPlane_R{}_{}{}Scaled'.format(jetR, obs_label, self.suffix)
    hLund = self.fMC.Get(name)
    
    hLund.GetXaxis().SetRangeUser(np.log(1/jetR), 5)
    hLund.GetYaxis().SetRangeUser(-3., 6.)
    
    text = '#it{p}_{T, ch jet}^{truth} > 100 GeV/c'
    output_filename = os.path.join(self.output_dir, 'lund/hLundPlane_R{}_{}.pdf'.format(jetR, obs_label))
    self.plot_hist(hLund, output_filename, drawOptions = 'colz', text = text)

  #---------------------------------------------------------------
  def plot_delta_pt(self, jetR, pt_bins):
  
    #self.plot_delta_pt_RC(jetR, 'before')
    #self.plot_delta_pt_RC(jetR, 'after')
    self.plot_delta_pt_emb(jetR, pt_bins)

  #---------------------------------------------------------------
  def plot_delta_pt_RC(self, jetR, CS_label):
  
    name = 'hDeltaPt_RC_{}CS_R{}{}Scaled'.format(CS_label, jetR, self.suffix)
    hDeltaPt = self.fMC.Get(name)
    hDeltaPt.SetMarkerStyle(21)
    hDeltaPt.SetMarkerSize(2)
    hDeltaPt.SetMarkerColor(self.ColorArray[0])
    hDeltaPt.GetYaxis().SetTitle('#frac{dN}{d#delta#it{p}_{T}}')
    if 'after' in CS_label:
      min = 0
      hDeltaPt.GetXaxis().SetTitle('#delta#it{p}_{T} #equiv #it{p}_{T}^{RC}')
    else:
      min = -50
      hDeltaPt.GetXaxis().SetTitle('#delta#it{p}_{T} #equiv #it{p}_{T}^{RC} - #pi#rho#it{R}^{2}')
    hDeltaPt.GetXaxis().SetRangeUser(min, 50)
    hDeltaPt.GetYaxis().SetRangeUser(10, 100*hDeltaPt.GetMaximum())
    
    mean = hDeltaPt.GetMean()
    std_dev = hDeltaPt.GetStdDev()
    text = 'Mean: {:.2f}, #sigma: {:.2f}'.format(mean, std_dev)
    
    output_filename = os.path.join(self.output_dir, 'delta_pt/hDeltaPt_RC_{}CS_R{}.pdf'.format(CS_label, self.remove_periods(jetR)))
    self.plot_hist(hDeltaPt, output_filename, setLogy = True, text = text)
    
  #---------------------------------------------------------------
  def plot_delta_pt_emb(self, jetR, pt_bins):
  
    name = 'hDeltaPt_emb_R{}{}Scaled'.format(jetR, self.suffix)
    
    c = ROOT.TCanvas('c','c: hist',600,450)
    c.cd()
    c.SetBottomMargin(0.17)
    c.SetLeftMargin(0.17)
    
    leg = ROOT.TLegend(0.57,0.53,0.9,0.87, '')
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1)
    leg.SetTextSize(0.04)
    
    # Loop through pt slices, and plot final residual for each 1D distribution
    for i in range(0, len(pt_bins) - 1):
      min_pt_truth = pt_bins[i]
      max_pt_truth = pt_bins[i+1]
      
      hDeltaPt = self.get_delta_pt_proj(name, 'h_delta_pt{}'.format(i), min_pt_truth, max_pt_truth)
      hDeltaPt.SetMarkerStyle(20)
      hDeltaPt.SetMarkerColor(self.ColorArray[i])
      hDeltaPt.SetLineColor(self.ColorArray[i])

      if i == 0:
      
        hDeltaPt.GetXaxis().SetTitleOffset(1.6);
        hDeltaPt.GetYaxis().SetTitleOffset(1.6);
        hDeltaPt.GetXaxis().SetTitle('#delta#it{p}_{T} #equiv #it{p}_{T,jet}^{combined} - #it{p}_{T,jet}^{pp-det}')
        hDeltaPt.GetYaxis().SetTitle('#frac{dN}{d#delta#it{p}_{T}}')
        hDeltaPt.GetXaxis().SetRangeUser(-50, 50)
        hDeltaPt.DrawCopy('P E')
        
      else:

        hDeltaPt.DrawCopy('P E same')

      leg.AddEntry(hDeltaPt, '#it{{p}}_{{T}}^{{gen}} = {}-{} GeV'.format(min_pt_truth, max_pt_truth), 'P')
      
      mean = hDeltaPt.GetMean()
      std_dev = hDeltaPt.GetStdDev()
      text = 'Mean: {:.2f}, #sigma: {:.2f}'.format(mean, std_dev)
      leg.AddEntry(None, text, '')

    leg.Draw('same')
    
    output_filename = os.path.join(self.output_dir, 'delta_pt/hDeltaPt_Emb_R{}.pdf'.format(self.remove_periods(jetR)))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  # Get delta_pt for a fixed pT-gen
  def get_delta_pt_proj(self, name, label, minPt, maxPt):
    
    h_delta_pt = self.fMC.Get(name)
    h_delta_pt.SetName('{}_{}'.format(h_delta_pt.GetName(), label))
    
    h_delta_pt.GetXaxis().SetRangeUser(minPt, maxPt)
    h = h_delta_pt.ProjectionY()
    
    integral = h.Integral()
    if integral > 0:
      h.Scale(1./integral, 'width')
    
    return h

  #---------------------------------------------------------------
  def plot_prong_matching(self, i_overlay, jetR, name_prefix, obs_subconfig_list, obs_settings, grooming_settings, overlay_list, prong_match_threshold):
  
    c = ROOT.TCanvas('c','c: hist',600,450)
    c.cd()
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.15)
    
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 20, 0., 200.)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.GetYaxis().SetTitleOffset(1.2)
    myBlankHisto.SetMaximum(1.2)
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.GetYaxis().SetTitle('Fraction tagged')
    myBlankHisto.Draw()
    
    myLegend = ROOT.TLegend(0.55,0.3,0.85,0.5)
    self.setup_legend(myLegend,0.035)
    
    h_list = [] # Store hists in a list, since otherwise it seems I lose the marker information
                # (removed from memory?)

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
      
      name = '{}_{}{}Scaled'.format(name_prefix, obs_label, self.suffix)
      hFraction_vs_pt = self.fMC.Get(name)
      xtitle = hFraction_vs_pt.GetXaxis().GetTitle()
      myBlankHisto.GetXaxis().SetTitle(xtitle)
      
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
      hMatchedFraction_vs_pt.SetMarkerStyle(marker)
      hMatchedFraction_vs_pt.SetMarkerColor(self.ColorArray[i])
      hMatchedFraction_vs_pt.SetLineColor(self.ColorArray[i])
      hMatchedFraction_vs_pt.DrawCopy('P same')
      myLegend.AddEntry(hMatchedFraction_vs_pt, self.formatted_grooming_label(grooming_setting), 'P')
        
      h_list.append(hMatchedFraction_vs_pt)

    myLegend.Draw('same')

    outdir = 'prong_matching_fraction_pt'
    if 'JetPtDet' in name_prefix:
      outdir = 'prong_matching_fraction_ptdet'
    output_filename = os.path.join(self.output_dir, '{}/{}_{}.pdf'.format(outdir, name_prefix, i_overlay))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  def plot_prong_matching_delta(self, i, jetR, name_prefix, obs_subconfig_list, obs_settings, grooming_settings, overlay_list, prong_match_threshold, min_pt, max_pt, plot_deltaz=False):

    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()
    ROOT.gPad.SetLeftMargin(0.2)
    ROOT.gPad.SetBottomMargin(0.15)
    
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 20, 0., jetR)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.GetYaxis().SetTitleOffset(1.2)
    max = 0.01
    myBlankHisto.SetMaximum(max)
    myBlankHisto.SetMinimum(0.)
    if plot_deltaz:
        myBlankHisto.GetXaxis().SetTitle('#Delta z_{prong}')
        myBlankHisto.GetYaxis().SetTitle('#frac{dN}{d#Delta z_{prong}}')
    else:
        myBlankHisto.GetXaxis().SetTitle('#Delta R_{prong}')
        myBlankHisto.GetYaxis().SetTitle('#frac{dN}{d#Delta R_{prong}}')
    myBlankHisto.Draw()
    
    leg = ROOT.TLegend(0.55,0.3,0.85,0.5)
    self.setup_legend(leg,0.035)
    
    h_list = [] # Store hists in a list, since otherwise it seems I lose the marker information
                # (removed from memory?)
    
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
      
      name = '{}_{}{}Scaled'.format(name_prefix, obs_label, self.suffix)
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
      hUnmatched_vs_pt.SetLineColor(self.ColorArray[i])
      if max < hUnmatched_vs_pt.GetMaximum():
        max = hUnmatched_vs_pt.GetMaximum()
        myBlankHisto.SetMaximum(max)
      hUnmatched_vs_pt.Draw('L hist same')
      leg.AddEntry(hUnmatched_vs_pt, self.formatted_grooming_label(grooming_setting), 'L')
      h_list.append(hUnmatched_vs_pt)

    text = 'p_{{T}} = {} - {} GeV/c'.format(int(min_pt), int(max_pt))
    leg.AddEntry(None, text, '')
    
    leg.Draw('same')

    if plot_deltaz:
        output_filename = os.path.join(self.output_dir, 'prong_matching_deltaZ/{}.pdf'.format(name_prefix))
    else:
        output_filename = os.path.join(self.output_dir, 'prong_matching_deltaR/{}.pdf'.format(name_prefix))
    c.SaveAs(output_filename)
    c.Close()

  #---------------------------------------------------------------
  def plot_prong_matching_correlation(self, jetR, hname_prefix, obs_subconfig_list, obs_settings, grooming_settings, overlay_list, prong_match_threshold):

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

      name = '{}_{}{}Scaled'.format(hname_prefix, obs_label, self.suffix)
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
      hMatchedFraction_vs_pt.SetMarkerStyle(marker)
      hMatchedFraction_vs_pt.SetMarkerColor(self.ColorArray[i])
      hMatchedFraction_vs_pt.SetLineColor(self.ColorArray[i])
      hMatchedFraction_vs_pt.Draw('P same')
      leg.AddEntry(hMatchedFraction_vs_pt, self.formatted_grooming_label(grooming_setting), 'P')
      h_list.append(hMatchedFraction_vs_pt)

    leg.Draw('same')

    output_filename = os.path.join(self.output_dir, 'prong_matching_correlation/{}.pdf'.format(hname_prefix))
    c.SaveAs(output_filename)
    c.Close()
