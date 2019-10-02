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

# Data analysis and plotting
import ROOT
import yaml

# Analysis utilities
from pyjetty.alice_analysis.analysis.base import analysis_utils

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

# Set plotting options
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

# Set debug level (0 = no debug info, 1 = some debug info, 2 = all debug info)
debugLevel = 0

analysis_utils = analysis_utils.analysis_utils()

#---------------------------------------------------------------
def plot_rg_performance(mcFile, dataFile, configFile, outputDir):
  
  fMC = ROOT.TFile(mcFile, 'READ')
  fData = ROOT.TFile(dataFile, 'READ')
  
  # Read config file
  with open(configFile, 'r') as stream:
    config = yaml.safe_load(stream)
  
  jetR_list = config['jetR']
  beta_list = list(config['beta'].keys())
  jet_matching_distance = config['jet_matching_distance']
  
  # Create output dir
  if not outputDir.endswith("/"):
    outputDir = outputDir + "/"
  if not os.path.exists(outputDir):
    os.makedirs(outputDir)

  for jetR in jetR_list:

    plotDeltaR(fMC, jetR, jet_matching_distance, outputDir)
    plotJER(fMC, jetR, outputDir)
    plotJES(fMC, jetR, outputDir)
    plotJESproj(fMC, jetR, outputDir)
    plotJetRecoEfficiency(fMC, jetR, outputDir)

    for beta in beta_list:

      plotRg(fMC, fData, jetR, beta, outputDir)

#---------------------------------------------------------------
def plotRg(fMC, fData, jetR, beta, outputDir):

  # (pt-det, pt-truth, theta_g-det, theta_g-truth)
  name = 'hResponse_JetPt_ThetaG_R{}_B{}'.format(jetR, beta)
  hRM = fMC.Get(name)
  
  if fData:
    name = 'hThetaG_JetPt_R{}_B{}'.format(jetR, beta)
    hThetaG_JetPt = fData.Get(name)
  else:
    hThetaG_JetPt = None

  plot2D_statistics(hThetaG_JetPt.Clone(), jetR, beta, outputDir)

  plotRgProjection(hRM, hThetaG_JetPt, jetR, beta, 20, 40, outputDir)
  plotRgProjection(hRM, hThetaG_JetPt, jetR, beta, 40, 60, outputDir)
  plotRgProjection(hRM, hThetaG_JetPt, jetR, beta, 60, 80, outputDir)
  plotRgProjection(hRM, hThetaG_JetPt, jetR, beta, 80, 100, outputDir)

#---------------------------------------------------------------
def plot2D_statistics(hThetaG_JetPt, jetR, beta, outputDir):

  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd()
  ROOT.gPad.SetLeftMargin(0.15)
  
  hThetaG_JetPt.SetMarkerSize(0.5)
  hThetaG_JetPt.GetXaxis().SetRangeUser(0, 100)
  hThetaG_JetPt.RebinX(5)
  hThetaG_JetPt.RebinY(5)
  hThetaG_JetPt.Draw('text colz')

  output_filename = os.path.join(outputDir, 'h2D_statistics_R{}_B{}.pdf'.format(analysis_utils.remove_periods(jetR), beta))
  c.SaveAs(output_filename)
  c.Close()

#---------------------------------------------------------------
def plotRgProjection(hRM, hThetaG_JetPt, jetR, beta, min_pt_det, max_pt_det, outputDir):

  rebin_val_det = 5
  rebin_val_truth = 1
  
  # Get histogram of theta_g in data, for given pt-det cut
  if hThetaG_JetPt:
    hThetaG_JetPt.GetXaxis().SetRange(min_pt_det, max_pt_det)
    hThetaG_data = hThetaG_JetPt.ProjectionY()
    hThetaG_data.GetYaxis().SetTitle('#frac{dN}{d#theta_{g}}')
    hThetaG_data.SetMarkerStyle(21)
    hThetaG_data.SetMarkerSize(1)
    analysis_utils.scale_by_integral(hThetaG_data)
    hThetaG_data.Rebin(rebin_val_det)

  # Get histograms of theta_g in MC, for a given pt-det cut
  hRM.GetAxis(0).SetRange(min_pt_det, max_pt_det)
  
  hThetaG_det = hRM.Projection(2)
  hThetaG_det.SetName('hThetaG_det')
  hThetaG_det.SetLineColor(2)
  hThetaG_det.SetLineWidth(2)
  analysis_utils.scale_by_integral(hThetaG_det)
  hThetaG_det.Rebin(rebin_val_det)

  hThetaG_truth = hRM.Projection(3)
  hThetaG_truth.SetName('hThetaG_truth')
  hThetaG_truth.GetYaxis().SetTitle('#frac{dN}{d#theta_{g}}')
  hThetaG_truth.SetLineColor(4)
  hThetaG_truth.SetLineWidth(2)
  analysis_utils.scale_by_integral(hThetaG_truth)
  hThetaG_truth.Rebin(rebin_val_truth)

  # Draw histogram
  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd()
  ROOT.gPad.SetLeftMargin(0.15)

  if hThetaG_JetPt:
    hThetaG_data.Draw('E P')
  hThetaG_truth.Draw('hist E same')
  hThetaG_det.Draw('hist E same')
  
  leg = ROOT.TLegend(0.65,0.75,0.85,0.9, "")
  leg.SetFillColor(10)
  leg.SetBorderSize(0)
  leg.SetFillStyle(1)
  leg.SetTextSize(0.04)
  if hThetaG_JetPt:
    leg.AddEntry(hThetaG_data, "data", "P")
  leg.AddEntry(hThetaG_det, "MC det", "L")
  leg.AddEntry(hThetaG_truth, "MC truth", "L")
  leg.Draw("same")
  
  text = 'R = {}'.format(jetR)
  textFit = ROOT.TLatex()
  textFit.SetTextSize(0.04)
  textFit.SetNDC()
  textFit.DrawLatex(0.65,0.7,text)
  
  text = '#beta = {}'.format(beta)
  textFit = ROOT.TLatex()
  textFit.SetTextSize(0.04)
  textFit.SetNDC()
  textFit.DrawLatex(0.65,0.65,text)
  
  text = 'pT,det = {}-{} GeV/c'.format(min_pt_det, max_pt_det)
  textFit = ROOT.TLatex()
  textFit.SetTextSize(0.035)
  textFit.SetNDC()
  textFit.DrawLatex(0.65,0.6,text)
  
  output_filename = os.path.join(outputDir, 'hTheta_g_MC_R{}_B{}_{}-{}.pdf'.format(analysis_utils.remove_periods(jetR), beta, min_pt_det, max_pt_det))
  c.SaveAs(output_filename)
  c.Close()

#---------------------------------------------------------------
def plotDeltaR(f, jetR, jet_matching_distance, outputDir):

  name = 'hDeltaR_All_R{}'.format(jetR)
  h = f.Get(name)
  
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

  output_filename = os.path.join(outputDir, '{}.pdf'.format(analysis_utils.remove_periods(name)))
  c.SaveAs(output_filename)
  c.Close()

#---------------------------------------------------------------
def plotJER(f, jetR, outputDir):
  
  name = 'hResponse_JetPt_R{}'.format(jetR)
  hRM = f.Get(name)
  
  # For each pT^gen, compute the standard deviation of the pT^det distribution
  
  # Get the pT^gen profile, with errors as standard deviation of pT^det distribution
  histPtGenProf = hRM.ProfileY("histPtGenProff", 1, -1, "s")
  
  # Create histo to be used to fill JER values
  nBins = 300
  histJER = ROOT.TH1D("histJER", "histJER", nBins, 0., 300.) # same binning for pT^gen as in task
  
  # Loop through the bins, and fill the JER
  for bin in range(0,nBins+1):
    sigma = histPtGenProf.GetBinError(bin)
    pTgen = histPtGenProf.GetXaxis().GetBinCenter(bin)
    JER = sigma/pTgen
    histJER.SetBinContent(bin, JER)
  
  histJER.GetYaxis().SetTitle("#frac{#sigma(#it{p}_{T}^{gen})}{#it{p}_{T}^{gen}}")
  histJER.GetXaxis().SetTitle("#it{p}_{T}^{gen}")
  histJER.GetYaxis().SetRangeUser(-0.01, 0.5)
  histJER.GetXaxis().SetRangeUser(5., 300.)
  outputFilename = os.path.join(outputDir, "histJER_R{}.pdf".format(analysis_utils.remove_periods(jetR)))
  histJER.SetMarkerStyle(21)
  histJER.SetMarkerColor(2)
  analysis_utils.plotHist(histJER, outputFilename, "hist P")

#---------------------------------------------------------------
def plotJES(f, jetR, outputDir):
  
  name = 'hJES_R{}'.format(jetR)
  histDeltaJES = f.Get(name)
  histDeltaJES.GetXaxis().SetTitle("#it{p}_{T}^{gen}")
  histDeltaJES.GetYaxis().SetTitle("#frac{#it{p}_{T}^{det} - #it{p}_{T}^{gen}}{#it{p}_{T}^{gen}}")
  histDeltaJES.GetXaxis().SetRangeUser(0., 200.)
  outputFilename = os.path.join(outputDir, "histDeltaJES_R{}.pdf".format(analysis_utils.remove_periods(jetR)))
  #analysis_utils.plotHist(histDeltaJES, outputFilename, "colz", False, True)
  
  histDeltaJES.RebinX(4)
  histDeltaJESprof = histDeltaJES.ProfileX()
  histDeltaJESprof.GetXaxis().SetRangeUser(0., 200.)
  histDeltaJESprof.GetYaxis().SetRangeUser(-0.8, 0.2)
  histDeltaJESprof.GetYaxis().SetTitle("#frac{#it{p}_{T}^{det} - #it{p}_{T}^{gen}}{#it{p}_{T}^{gen}}")
  outputFilename = os.path.join(outputDir, "histDeltaJESprof_R{}.pdf".format(analysis_utils.remove_periods(jetR)))
  histDeltaJESprof.SetMarkerStyle(21)
  histDeltaJESprof.SetMarkerColor(4)
  analysis_utils.plotHist(histDeltaJESprof, outputFilename, "P")


#---------------------------------------------------------------
# Plot JES shift distribution for various fixed pT-gen
def plotJESproj(f, jetR, outputDir):
  
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
  
  hJESProj1 = getJESshiftProj(f, name, "hJESproj1", 20, 30)
  hJESProj2 = getJESshiftProj(f, name, "hJESproj2", 50, 70)
  hJESProj3 = getJESshiftProj(f, name, "hJESproj3", 100, 120)
  hJESProj1b = getJESshiftProj(f, name, "hJESproj1b", 30, 40)
  hJESProj2b = getJESshiftProj(f, name, "hJESproj2b", 40, 50)
  hJESProj3b = getJESshiftProj(f, name, "hJESproj3b", 70, 90)
  
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
  
  outputFilename = os.path.join(outputDir, "histDeltaJESproj_R{}.pdf".format(analysis_utils.remove_periods(jetR)))
  cJES.SaveAs(outputFilename)
  cJES.Close()

#---------------------------------------------------------------
# Get JES shift distribution for a fixed pT-gen
def getJESshiftProj(f, name, label, minPt, maxPt):
  
  histDeltaJES = f.Get(name)
  histDeltaJES.SetName('{}_{}'.format(histDeltaJES.GetName(), label))
  
  histDeltaJES.GetXaxis().SetRangeUser(minPt, maxPt)
  histDeltaJES.GetYaxis().SetRangeUser(-1., 1.)
  h = histDeltaJES.ProjectionY()
  
  integral = h.Integral()
  if integral > 0:
    h.Scale(1./integral, 'width')
  
  return h

#---------------------------------------------------------------
def plotJetRecoEfficiency(f, jetR, outputDir):
  
  # For each pT^gen, compute the fraction of matched pT^gen
  
  # First, get the pT^gen spectrum
  name = 'hJetPt_Truth_R{}'.format(jetR)
  histPtGen = f.Get(name)
  
  # Then, get the pT^gen spectrum for matched jets
  name = 'hResponse_JetPt_R{}'.format(jetR)
  hRM = f.Get(name)
  histPtGenMatched = hRM.ProjectionY("_py",1,hRM.GetNbinsX()) #avoid under and overflow bins
  
  # Compute the ratio
  histEfficiency = histPtGenMatched.Clone()
  histEfficiency.SetName("histEfficiency_{}".format(name))
  histEfficiency.Divide(histPtGenMatched, histPtGen, 1., 1., "B")
  
  histEfficiency.GetXaxis().SetTitle("#it{p}_{T}^{gen}")
  histEfficiency.GetYaxis().SetTitle("Efficiency")
  histEfficiency.GetXaxis().SetRangeUser(0., 300.)
  histEfficiency.GetYaxis().SetRangeUser(0.4, 1.2)
  histEfficiency.SetMarkerStyle(21)
  histEfficiency.SetMarkerColor(1)
  outputFilename = os.path.join(outputDir, '{}_R{}.pdf'.format(analysis_utils.remove_periods(histEfficiency.GetName()), analysis_utils.remove_periods(jetR)))
  analysis_utils.plotHist(histEfficiency, outputFilename)

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

  plot_rg_performance(mcFile = args.mcFile, dataFile = args.dataFile, configFile = args.configFile, outputDir = args.outputDir)
