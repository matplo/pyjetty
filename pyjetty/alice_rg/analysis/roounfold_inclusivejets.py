#! /usr/bin/env python

import sys
import os
import argparse
import itertools
from array import *
import numpy
import ROOT

# Load the RooUnfold library
ROOT.gSystem.Load("$ALIBUILD_WORK_DIR/slc7_x86-64/RooUnfold/latest/lib/libRooUnfold.so")

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

# Suppress a lot of standard output
#ROOT.gErrorIgnoreLevel = ROOT.kWarning

# Set plotting options
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

###########################################################################################
###########################################################################################
def roounfold_inclusivejets(input_file_data, input_file_response, output_dir, file_format):

  fData = ROOT.TFile(input_file_data)
  fResponse = ROOT.TFile(input_file_response)
  
  # Create output dir for unfolding histograms and result
  if not os.path.isdir(output_dir):
    os.makedirs(output_dir)
  
  # Read config file
  with open(configFile, 'r') as stream:
    config = yaml.safe_load(stream)
  
  jetR_list = config['jetR']
  beta_list = config['beta']
  
  # Unfolding settings
  reg_param = 4

  #--------------------------------------------------------------

  # Set pT range of input spectrum for unfolding
  min_pt_det = 10
  max_pt_det = 100

  # Set pT range of output spectrum
  min_pt_reported = 20
  max_pt_reported = 100

  # Set pT range of response spectrum
  min_pt_gen = 10
  max_pt_gen = 300

  # Define pT-det and pT-truth binning
  pt_bin_array_truth = ([min_pt_gen, 20, 30, 40, 50, 60, 70, 80, 100, 120, 140, 190, 240, max_pt_gen])
  pt_bin_array_det = ([min_pt_det, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, max_pt_det])

  n_pt_bins_det = len(pt_bin_array_det) - 1
  det_pt_bin_array = array('d',pt_bin_array_det)
  n_pt_bins_truth = len(pt_bin_array_truth) - 1
  pt_truth_bin_array = array('d',pt_bin_array_truth)
  print('n_pt_bins_det: {}'.format(n_pt_bins_det))
  print('n_pt_bins_truth: {}'.format(n_pt_bins_truth))
  
  #--------------------------------------------------------------
  
  # Set pT range of input spectrum for unfolding
  min_rg_det = 0.
  max_rg_det = 1.2
  
  # Set pT range of output spectrum
  min_rg_reported = 20
  max_rg_reported = 100
  
  # Set pT range of response spectrum
  min_rg_gen = 0.
  max_rg_gen = 1.5
  
  # Define pT-det and pT-truth binning
  rg_bin_array_truth = ([min_rg_gen, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, max_rg_gen])
  rg_bin_array_det = ([min_rg_det, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, max_rg_det])
  
  n_rg_bins_det = len(rg_bin_array_det) - 1
  det_rg_bin_array = array('d',rg_bin_array_det)
  n_rg_bins_truth = len(rg_bin_array_truth) - 1
  rg_truth_bin_array = array('d',rg_bin_array_truth)
  print('n_rg_bins_det: {}'.format(n_rg_bins_det))
  print('n_rg_bins_truth: {}'.format(n_rg_bins_truth))
  
  #--------------------------------------------------------------
  
  for jetR in jetR_list:
    for beta in beta_list:
  
  
  
  
  
  
  
  # Define an empty dictionary to store the final unfolded histograms from each label
  hUnfoldedResultDict = {}

  unfoldSingleOutputList(fData, fResponse, hUnfoldedResultDict, reg_param, min_pt_det, max_pt_det, n_bins_det, det_bin_array, min_pt_gen, max_pt_gen, n_bins_truth, truth_bin_array, min_pt_reported, max_pt_reported, output_dir, file_format)

###################################################################################################
# Unfold jet spectrum from a single output list
###################################################################################################
def unfoldSingleOutputList(fData, fResponse, hUnfoldedResultDict, reg_param, min_pt_det, max_pt_det, n_bins_det, det_bin_array, min_pt_gen, max_pt_gen, n_bins_truth, truth_bin_array, min_pt_reported, max_pt_reported, output_dir, file_format):

  # Get N events
  hname_event = 'hNevents'
  hNevent_data = fData.Get(hname_event)
  n_events_data = hNevent_data.GetBinContent(2)
  print('N accepted events in data: {}'.format(n_events_data))
  
  hNevent_response = fResponse.Get(hname_event)
  n_events_response = hNevent_response.GetBinContent(2)/20.
  print('N accepted events in response (avg per bin): {}'.format(n_events_response))

  # Get data jet spectrum
  hname_jetpt_data = 'hJetPt_R0.4'
  hJetSpectrumMeasuredPerBin = getMeasuredSpectrum(fData, hname_jetpt_data, min_pt_det, max_pt_det, n_bins_det, det_bin_array)
  hJetSpectrumMeasuredPerBin.Sumw2()
  outputFilename = os.path.join(output_dir, "hJetSpectrumMeasuredPerBin" + file_format)
  plotHist(hJetSpectrumMeasuredPerBin, outputFilename, "P E", True)

  # Plot the fine-binned response matrix, with normalization of each truth projection bin to 1 (over specified pT ranges)
  hname_response = 'hResponse_JetPt_R0.4Scaled'
  hResponseMatrixFineBinned = getResponseMatrix(fResponse, hname_response, min_pt_gen, max_pt_gen, 0, max_pt_gen, 0, 0, 0, 0, "uncutFineBinned", output_dir)
  hResponseMatrixFineBinned.SetName("hResponseMatrix_FineBinned")
  normalizeResponseMatrix(hResponseMatrixFineBinned, min_pt_det, max_pt_det, min_pt_gen, max_pt_gen, output_dir, file_format)
  
  # Plot the re-binned response matrix, with normalization of each truth projection bin to 1 (over specified pT ranges)
  hResponseMatrixRebinned = getResponseMatrix(fResponse, hname_response, min_pt_det, max_pt_det, min_pt_gen, max_pt_gen, n_bins_det, det_bin_array, n_bins_truth, truth_bin_array, "Rebinned", output_dir)
  hResponseMatrixRebinned.SetName("hResponseMatrix_Rebinned")
  normalizeResponseMatrix(hResponseMatrixRebinned, min_pt_det, max_pt_det, min_pt_gen, max_pt_gen, output_dir, file_format)
  
  # Plot the kinematic efficiency, using the response matrix  maxPtReported
  plotKinematicEfficiency(fResponse, hname_response, min_pt_det, max_pt_det, 0, max_pt_gen, min_pt_reported, max_pt_reported, output_dir, file_format)

  # Get truth-level spectrum (matched) from response matrix projection, before cutting the pT-det
  # range, do not rebin at this point since it will be cut to the range otherwise
  hResponseMatrixUncut = getResponseMatrix(fResponse, hname_response, 0, max_pt_gen, min_pt_gen, max_pt_gen, 0, 0, 0, 0, "uncut", output_dir)
  hJetSpectrumTrueUncutPerBin = hResponseMatrixUncut.ProjectionY()
  # rebin only the projcetion to keep an uncut range (for kinematic efficiency correction)
  hJetSpectrumTrueUncutPerBin = hJetSpectrumTrueUncutPerBin.Rebin(len(truth_bin_array)-1, "{}_NewBinning".format(hJetSpectrumTrueUncutPerBin.GetName()), truth_bin_array)
  hJetSpectrumTrueUncutPerBin.SetName("hJetSpectrumTrueUncutPerBin")
  outputFilename = os.path.join(output_dir, "hJetSpectrumTrueUncutPerBin" + file_format)
  #plotHist(hJetSpectrumTrueUncutPerBin, outputFilename, "hist", True)

  # Get response matrix from response file (Measured, True) to be used for the unfolding,
  # with pT-det range cut to desired range, and re-bin.
  hResponseMatrix = getResponseMatrix(fResponse, hname_response, min_pt_det, max_pt_det, min_pt_gen, max_pt_gen, n_bins_det, det_bin_array, n_bins_truth, truth_bin_array, "", output_dir)

  # Get the truth-level jet spectrum (matched) from response matrix (already re-binned)
  hJetSpectrumTruePerBin = hResponseMatrix.ProjectionY("_py",1,hResponseMatrix.GetNbinsX()) # Do exclude under and overflow bins
  hJetSpectrumTruePerBin.SetName("hJetSpectrumTruePerBin")
  outputFilename = os.path.join(output_dir, "hJetSpectrumTruePerBin" + file_format)
  plotHist(hJetSpectrumTruePerBin, outputFilename, "hist", True)

  # Unfold spectrum
  if hJetSpectrumMeasuredPerBin and hJetSpectrumTruePerBin and hJetSpectrumTrueUncutPerBin and hResponseMatrix:
    
    # Generate the response objects, scale spectra appropriately, and plot some basic checks
    response = prepareToUnfold(hJetSpectrumMeasuredPerBin, hJetSpectrumTruePerBin, hResponseMatrix, hJetSpectrumTrueUncutPerBin, min_pt_det, max_pt_det, n_events_data, output_dir, file_format)
    responseNoKinEff = ROOT.RooUnfoldResponse(0, 0, hResponseMatrix, "hResponseMatrixNoKinEff", "hResponseMatrixNoKinEff")

    unfoldJetSpectrum(hJetSpectrumMeasuredPerBin, hJetSpectrumTruePerBin, response, responseNoKinEff, hUnfoldedResultDict, reg_param, n_bins_truth, truth_bin_array, min_pt_reported, max_pt_reported, output_dir, file_format)

#################################################################################################
# Unfold jet spectrum, SVD
#################################################################################################
def unfoldJetSpectrum(hJetSpectrumMeasuredPerBin, hJetSpectrumTruePerBin, response, responseNoKinEff, hUnfoldedResultDict, regParamFinal, nBinsTruth, truthBinArray, minPtReported, maxPtReported, outputDir, fileFormat):
  
  regularizationParamName = "k"
  
  # Create canvas to superpose iterations on one plot, to examine convergence
  c1 = ROOT.TCanvas("c1","c1: histos",600,450)
  c1.cd()
  c1.SetLogy()
  ROOT.gPad.SetLeftMargin(0.15)
  leg = ROOT.TLegend(0.6,0.5,0.88,0.83,"{} Unfolding".format(type))
  leg.SetFillColor(10)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.04)

  hJetSpectrumUnfoldedPerBin = None
  hJetSpectrumUnfoldedPerGeV = None
  
  # Loop over values of regularization parameter
  for i in range(1,regParamFinal+3):
    
    errorType = ROOT.RooUnfold.kCovToy

    # Set up the SVD unfolding object
    unfoldSVD = ROOT.RooUnfoldSvd(response, hJetSpectrumMeasuredPerBin, i)
    unfoldSVD.SetNToys(1000)
    
    # Perform the unfolding
    hJetSpectrumUnfoldedPerGeV = unfoldSVD.Hreco(errorType) # Produces the truth distribution, with errors, PerBin (will scale by bin width below, after refolding checks)
    hJetSpectrumUnfoldedPerGeV.SetName("hJetSpectrumUnfoldedPerGeV")
    
    # Plot Pearson correlation coeffs for each k, to get a measure of the correlation between the bins
    covarianceMatrix = unfoldSVD.Ereco(errorType) # Get the covariance matrix
    #plotCorrelationCoefficients(covarianceMatrix, i, outputDir, fileFormat)
    
    # Plot d-vector, as a function of k
    if i is 1:
      svdUnfoldObject = unfoldSVD.Impl() # Get TSVDUnfold object
      hDVector = svdUnfoldObject.GetD()
      hDVector.GetXaxis().SetTitle("k")
      hDVector.GetYaxis().SetTitle("d_{k}")
      hDUsed = hDVector.Clone("hDUsed")
      for bin in range(1, hDVector.GetNbinsX() + 1):
        if bin!=regParamFinal:
          hDUsed.SetBinContent(bin,0)
    
      outputFilename = os.path.join(outputDir, "hDVector{}".format(fileFormat))
      plotHist(hDVector, outputFilename, "", True)
    
    if i is regParamFinal:
      hUnfoldedResultDict["Final"] = hJetSpectrumUnfoldedPerGeV

    # Store the spectra near the optimal value of k, for later plotting the ratio
    if i is regParamFinal-1:
      hLowerkResult = hJetSpectrumUnfoldedPerGeV.Clone()
      hLowerkResult.SetName("hLowerkResult")
      hLowerkResult.SetLineColor(1) # Set colors to match the ratio plot
      hLowerkResult.GetXaxis().SetRangeUser(minPtReported, maxPtReported)
    elif i is regParamFinal:
      hMainResult = hJetSpectrumUnfoldedPerGeV.Clone()
      hMainResult.SetName("hMainResult")
      hMainResult.SetLineColor(4)
      hMainResult.GetXaxis().SetRangeUser(minPtReported, maxPtReported)
    elif i is regParamFinal+1:
      hHigherkResult = hJetSpectrumUnfoldedPerGeV.Clone()
      hHigherkResult.SetName("hHigherResult")
      hHigherkResult.SetLineColor(2) # Set colors to match the ratio plot
      hHigherkResult.GetXaxis().SetRangeUser(minPtReported, maxPtReported)
    # This is for the special case of having a very low regularization parameter
    # in this case take the two higher k values as variations
    if regParamFinal==2 and i is regParamFinal+2:
      hLowerkResult = hJetSpectrumUnfoldedPerGeV.Clone()
      hLowerkResult.SetName("hLowerkResult")
      hLowerkResult.SetLineColor(1) # Set colors to match the ratio plot
      hLowerkResult.GetXaxis().SetRangeUser(minPtReported, maxPtReported)
  
  #Plot the spectra comparing only k=+1 and k-1 to the main result
  xRangeMin = minPtReported
  xRangeMax = maxPtReported
  yAxisTitle = "#frac{d#sigma}{dp_{T}} [mb (GeV/c)^{-1}]"
  ratioYAxisTitle = "Ratio to k={}".format(regParamFinal)
  outputFilename = os.path.join(outputDir, "hJetSpectraUnfoldedRatio" + fileFormat)
  legendTitle = "{} Unfolding".format(type)
  hLegendLabel = "k = {}".format(regParamFinal-1)
  h2LegendLabel = "k = {}".format(regParamFinal)
  h3LegendLabel = "k = {}".format(regParamFinal+1)
  if regParamFinal==2:
    hLegendLabel = "k = {}".format(regParamFinal+2)
  if hLowerkResult and hMainResult and hHigherkResult:
    # To get sensible error bars, assume main result has no errors
    for bin in range(1, hMainResult.GetNbinsX() + 1):
      hMainResult.SetBinError(bin, 0)
    plotSpectra(hLowerkResult, hMainResult, hHigherkResult, 1., xRangeMin, xRangeMax, yAxisTitle, ratioYAxisTitle, outputFilename, "", legendTitle, hLegendLabel, h2LegendLabel, h3LegendLabel)

###################################################################################################
# Unfold jet spectrum from a single output list
###################################################################################################
def getMeasuredSpectrum(fData, hname_jetpt_data, min_pt_det, max_pt_det, n_bins_det, det_bin_array):

  jetHistogram = fData.Get(hname_jetpt_data)
  jetHistogram.GetXaxis().SetRangeUser(min_pt_det, max_pt_det)
  
  hJetSpectrumRebinned = jetHistogram.Rebin(n_bins_det, "{}New".format(hname_jetpt_data), det_bin_array)
  return hJetSpectrumRebinned

###################################################################################################
# Get response matrix from response file (Measured, True) and rebin it
###################################################################################################
def getResponseMatrix(fResponse, hname_response, min_pt_det, max_pt_det, min_pt_gen, max_pt_gen, n_bins_det, det_bin_array, n_bins_truth, truth_bin_array, label, output_dir):
  
  hResponseMatrixFineBinned = fResponse.Get(hname_response)
  hResponseMatrixFineBinned.GetYaxis().SetRangeUser(min_pt_gen, max_pt_gen)
  hResponseMatrixFineBinned.GetXaxis().SetRangeUser(min_pt_det, max_pt_det)
  
  # Create a new fine-binned histogram with the appropriate min,max cuts
  # Loop over all bins in fine-binned response matrix, and fill appropriate bin in new response matrix
  # Assume that the bin edges overlap appropriately
  histname = "{}_{}".format(hResponseMatrixFineBinned.GetName(), label)
  title = histname + ";#it{p}_{T,corr}^{det} (GeV/#it{c});#it{p}_{T}^{truth} (GeV/#it{c})"
  hResponseMatrixNew = ROOT.TH2D(histname, title, max_pt_det-min_pt_det, min_pt_det, max_pt_det, max_pt_gen-min_pt_gen, min_pt_gen, max_pt_gen)
  for ibin in range(1, hResponseMatrixFineBinned.GetNbinsX() + 1):
    for jbin in range(1, hResponseMatrixFineBinned.GetNbinsY() + 1):
      
      oldContent = hResponseMatrixFineBinned.GetBinContent(ibin, jbin)
      
      # Find the bin that should be filled in the new histogram, and fill it
      # Need to get (x,y) location from bins (ibin, jbin)
      x = hResponseMatrixNew.GetXaxis().GetBinCenter(ibin)
      y = hResponseMatrixNew.GetYaxis().GetBinCenter(jbin)
      if x > min_pt_det and x < max_pt_det and y > min_pt_gen and y < max_pt_gen:
        hResponseMatrixNew.Fill(x, y, oldContent)

      # Assume 0 errors on response matrix
      #for bin in range(1, hResponseMatrixNew.GetNcells() + 1):
      # hResponseMatrixNew.SetBinError(bin, 0)

  # Re-bin the response matrix, if a binning is provided
  if n_bins_det > 1 and n_bins_truth > 1:
    hResponseMatrixRebinned = rebinResponseMatrix(hResponseMatrixFineBinned, n_bins_det, det_bin_array, n_bins_truth, truth_bin_array)
    return hResponseMatrixRebinned
  else:
    return hResponseMatrixFineBinned

##################################################################################################
# Rebin the response matrix to have variable binning
##################################################################################################
def rebinResponseMatrix(hResponseMatrix, nBinsDet, detBinArray, nBinsTruth, truthBinArray):
  
  histname = "{}NewRebinned".format(hResponseMatrix.GetName())
  title = histname + ";#it{p}_{T,corr}^{det} (GeV/#it{c});#it{p}_{T}^{truth} (GeV/#it{c})"
  hResponseMatrixNew = ROOT.TH2D(histname, title, nBinsDet, detBinArray, nBinsTruth, truthBinArray)
  
  # Loop over all bins in fine-binned response matrix, and fill appropriate bin in new response matrix
  # Assume that the bin edges overlap appropriately
  for ibin in range(1, hResponseMatrix.GetNbinsX() + 1):
    for jbin in range(1, hResponseMatrix.GetNbinsY() + 1):
      
      oldContent = hResponseMatrix.GetBinContent(ibin, jbin)
      
      # Find the bin that should be filled in the new histogram, and fill it
      # Need to get (x,y) location from bins (ibin, jbin)
      x = hResponseMatrix.GetXaxis().GetBinCenter(ibin)
      y = hResponseMatrix.GetYaxis().GetBinCenter(jbin)
      hResponseMatrixNew.Fill(x, y, oldContent)

  # Assume 0 errors on response matrix
  for bin in range(1, hResponseMatrixNew.GetNcells() + 1):
    hResponseMatrixNew.SetBinError(bin, 0)
  
  return hResponseMatrixNew

#################################################################################################
# Prepare to unfold jet spectrum
# Returns RooUnfoldResponse object.
#################################################################################################
def prepareToUnfold(hJetSpectrumMeasuredPerBin, hJetSpectrumTruePerBin, hResponseMatrix, hJetSpectrumTrueUncutPerBin, minPtDet, maxPtDet, nEventsData, outputDir, fileFormat):
  
  visibleMBCrossSection = 50.87 # (mb) V0AND cross section (https://cds.cern.ch/record/2648933)
  vertexEfficiency = 0.95
  hJetSpectrumMeasuredPerBin.Scale(visibleMBCrossSection)
  hJetSpectrumMeasuredPerBin.Scale(vertexEfficiency)
  hJetSpectrumMeasuredPerBin.Scale(1./nEventsData)
  
  # Normalize the response matrix by setting each pT-truth projection to the intended prior distribution
  # Keep response matrix as per-bin probabilities (i.e. don't scale by bin width)
  # Scale also hJetSpectrumTrueUncutPerBin accordingly, soas to preserve the kinematic efficiency (which we will use to create the RooUnfoldResponse object)
  setPrior(hResponseMatrix, hJetSpectrumTrueUncutPerBin, outputDir, fileFormat)
  
  # Set up the RooUnfoldResponse object
  # response = RooUnfoldResponse(0, 0, hResponseMatrix, "hResponseMatrixMain", "hResponseMatrixMain")
  # One can supply the measured and truth distributions is to incorporate fakes and inefficiency
  # For the truth-level, we pass in the projection before the RM was cut to a specific pT-det range, in order to account for the kinematic efficiency
  response = ROOT.RooUnfoldResponse(0, hJetSpectrumTrueUncutPerBin, hResponseMatrix, "hResponseMatrixMain", "hResponseMatrixMain")
  response.UseOverflow(False)
  
  return response

##################################################################################################
# Set prior in repsonse matrix
# RooUnfold takes the prior as the truth-axis projection of the response matrix.
# After normalizing the response matrix to our desired prior, RooUnfold will then
# automatically normalize the response matrix to conserve the number of jets (i.e. each bin of
# truth-axis projection normalized to 1).
##################################################################################################
def setPrior(hResponseMatrix, hJetSpectrumTrueUncutPerBin, outputDir, fileFormat):
  
  priorFolder = os.path.join(outputDir, "PriorProperties/")
  if not os.path.exists(priorFolder):
    os.makedirs(priorFolder)
  
  # Plot response matrix before normalization
  outputFilename = os.path.join(priorFolder, "hResponseMatrixBeforeNormalization" + fileFormat)
  plotHist(hResponseMatrix, outputFilename, "colz")
  
  # Make projection onto pT-true axis (y-axis), and scale appropriately
  hTruthProjectionBefore = hResponseMatrix.ProjectionY("_py",1,hResponseMatrix.GetNbinsX()) # Do exclude under and overflow bins
  hTruthProjectionBefore.SetName("hTruthProjectionBefore")

  #hResponseMatrix has a detector min pT restriction
  outputFilename = os.path.join(priorFolder, "hTruthProjectionBeforeNormalization" + fileFormat)
  plotHist(hTruthProjectionBefore, outputFilename, "hist", True)
  
  #hJetSpectrumTrueUncutPerBin has no detector pT restriction
  outputFilename = os.path.join(priorFolder, "hJetSpectrumTrueUncutPerBin_asPrior{}".format(fileFormat))
  plotHist(hJetSpectrumTrueUncutPerBin, outputFilename,"hist", True)
  
  # Check kinematic efficiency before
  hKinematicEfficiencyBefore = hTruthProjectionBefore.Clone()
  hKinematicEfficiencyBefore.SetName("hKinematicEfficiencyBefore")
  hKinematicEfficiencyBefore.Divide(hTruthProjectionBefore, hJetSpectrumTrueUncutPerBin, 1., 1., "B")
  hKinematicEfficiencyBefore.SetMarkerStyle(21)
  outputFilename = os.path.join(priorFolder, "hKinematicEfficiencyBefore{}".format(fileFormat))
  plotHist(hKinematicEfficiencyBefore, outputFilename, "P E")
  
  # Loop through truth-level bins (of RM and uncut truth spectrum), and normalize each bin (of RM and uncut truth spectrum) to the prior
  nBinsY = hResponseMatrix.GetNbinsY() # pT-gen
  nBinsX = hResponseMatrix.GetNbinsX() # pT-det

  # Plot truth projection after normalization (i.e. prior distribution)
  hTruthProjectionAfter = hResponseMatrix.ProjectionY("_py",1,hResponseMatrix.GetNbinsX()) # Do exclude under and overflow bins
  hTruthProjectionAfter.SetName("hTruthProjectionAfterNormalization")
  outputFilename = os.path.join(priorFolder, "hTruthProjectionAfterNormalization" + fileFormat)
  plotHist(hTruthProjectionAfter, outputFilename, "hist", True)
  
  # Plot response matrix after normalization
  outputFilename = os.path.join(priorFolder, "hResponseMatrixAfterNormalization" + fileFormat)
  plotHist(hResponseMatrix, outputFilename, "colz")
  
  #Prior used for the unfolding
  outputFilename = os.path.join(priorFolder, "hPriorForUnfolding{}".format(fileFormat))
  plotHist(hJetSpectrumTrueUncutPerBin, outputFilename,"hist", True)
  
  # Check kinematic efficiency after normalization
  hKinematicEfficiencyAfter = hTruthProjectionAfter.Clone()
  hKinematicEfficiencyAfter.SetName("hKinematicEfficiencyAfter")
  hKinematicEfficiencyAfter.Divide(hTruthProjectionAfter, hJetSpectrumTrueUncutPerBin, 1., 1., "B")
  hKinematicEfficiencyAfter.SetMarkerStyle(21)
  outputFilename = os.path.join(priorFolder, "hKinematicEfficiencyAfter{}".format(fileFormat))
  plotHist(hKinematicEfficiencyAfter, outputFilename, "P E")

###################################################################################################
# Plot kinematic efficiency (i.e. (pT-truth projection of response matrix with measured pT-det range
# selected) / (pT-truth projection of response matrix with full pT-det range selected)
###################################################################################################
def plotKinematicEfficiency(fResponse, hname_response, min_pt_det, max_pt_det, min_pt_gen, max_pt_gen, min_pt_reported, max_pt_reported, output_dir, file_format):
  
  # Get fine-binned response matrix (Measured, True), with pT-det range cut to desired range, and project truth distribution
  hResponseMatrixCut = getResponseMatrix(fResponse, hname_response, min_pt_det, max_pt_det, 0, max_pt_gen, 0, 0, 0, 0, "Cut", output_dir)
  hJetSpectrumTrueCutPerBin = hResponseMatrixCut.ProjectionY()
  hJetSpectrumTrueCutPerBin.SetName("hJetSpectrumTrueCutPerBin")
  rebinVal = 5
  hJetSpectrumTrueCutPerBin.Rebin(5)
  hJetSpectrumTrueCutPerBin.Scale(1., "width")
  outputFilename = os.path.join(output_dir, "hJetSpectrumTrueCutPerBin.pdf")
  plotHist(hJetSpectrumTrueCutPerBin, outputFilename, "hist E", True)
  
  # Get fine-binned response matrix (Measured, True), with full pT-det range, and project truth distribution
  fMinPt = 0 # Min value of pTdet in response matrix
  hResponseMatrixUncut = getResponseMatrix(fResponse, hname_response, fMinPt, max_pt_gen, 0, max_pt_gen, 0, 0, 0, 0, "Uncut", output_dir)
  hJetSpectrumTrueUncutPerBin = hResponseMatrixUncut.ProjectionY()
  hJetSpectrumTrueUncutPerBin.SetName("hJetSpectrumTrueUncutPerBin_KinEff")
  hJetSpectrumTrueUncutPerBin.Rebin(5)
  hJetSpectrumTrueUncutPerBin.Scale(1., "width")
  outputFilename = os.path.join(output_dir, "hJetSpectrumTrueUncutPerBin.pdf")
  plotHist(hJetSpectrumTrueUncutPerBin, outputFilename, "hist E", True)
  
  # Plot the ratio of the spectra
  hKinematicEfficiency = hJetSpectrumTrueCutPerBin.Clone()
  hKinematicEfficiency.SetName("hKinematicEfficiency")
  hKinematicEfficiency.Divide(hJetSpectrumTrueCutPerBin, hJetSpectrumTrueUncutPerBin, 1., 1., "B")
  
  hKinematicEfficiency.GetXaxis().SetTitle("#it{p}_{T}^{gen}")
  hKinematicEfficiency.GetYaxis().SetTitle("Kinematic Efficiency")
  hKinematicEfficiency.SetMarkerStyle(21)
  hKinematicEfficiency.SetMarkerColor(2)
  
  text = "p_{T}^{det} #in [%d, %d] GeV" % (min_pt_det, max_pt_det)
  outputFilename = os.path.join(output_dir, "hKinematicEfficiency_{}_{}{}".format(min_pt_det, max_pt_det, file_format))
  plotHist(hKinematicEfficiency, outputFilename, "P E", False, False, text)

################################################################################################
# Normalize response matrix
# Normalize the pT-truth projection to 1
################################################################################################
def normalizeResponseMatrix(hResponseMatrix, minPtDet, maxPtDet, minPtGen, maxPtGen, outputDir, fileFormat):
  
  # Plot response matrix before normalization
  outputFilename = os.path.join(outputDir, "hResponseMatrixBeforeNormalization" + fileFormat)
  #plotHist(hResponseMatrix, outputFilename, "colz")
  
  # Make projection onto pT-true axis (y-axis), and scale appropriately
  hTruthProjectionBefore = hResponseMatrix.ProjectionY("_py",1,hResponseMatrix.GetNbinsX()) # Do exclude under and overflow bins
  hTruthProjectionBefore.SetName("hTruthProjectionBefore")
  
  # Loop through truth-level bins, and apply normalization factor to all bins.
  nBinsY = hResponseMatrix.GetNbinsY() # pT-gen
  nBinsX = hResponseMatrix.GetNbinsX() # pT-det
  for truthBin in range(1,nBinsY+1):
    normalizationFactor = hTruthProjectionBefore.GetBinContent(truthBin)
    if normalizationFactor > 0:
      truthBinCenter = hTruthProjectionBefore.GetXaxis().GetBinCenter(truthBin)
      
      for detBin in range(1,nBinsX+1):
        binContent = hResponseMatrix.GetBinContent(detBin, truthBin)
        hResponseMatrix.SetBinContent(detBin, truthBin, binContent/normalizationFactor)

  # Plot response matrix
  outputFilename = os.path.join(outputDir, "{}_{}_{}{}".format(hResponseMatrix.GetName(), minPtDet, maxPtDet, fileFormat))
  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd()
  c.cd().SetLeftMargin(0.15)

  hResponseMatrix.Draw("colz")
  line = ROOT.TLine(minPtDet,0,minPtDet,250)
  line.SetLineColor(0)
  line.SetLineStyle(2)
  line.Draw("same")
  line2 = ROOT.TLine(maxPtDet,0,maxPtDet,250)
  line2.SetLineColor(0)
  line2.SetLineStyle(2)
  line2.Draw("same")
  line3 = ROOT.TLine(0,minPtGen,100,minPtGen)
  line3.SetLineColor(0)
  line3.SetLineStyle(2)
  line3.Draw("same")
  line4 = ROOT.TLine(0,maxPtGen,100,maxPtGen)
  line4.SetLineColor(0)
  line4.SetLineStyle(2)
  line4.Draw("same")
  
  c.SaveAs(outputFilename)
  c.Close()

########################################################################################################
# Plot spectra and ratio of h (and h3, if supplied) to h2    ###########################################
########################################################################################################
def plotSpectra(h, h2, h3, nEvents, xRangeMin, xRangeMax, yAxisTitle, ratioYAxisTitle, outputFilename, scalingOptions = "", legendTitle = "",hLegendLabel = "", h2LegendLabel = "", h3LegendLabel = "", yRatioMax = 2.2):
  
  c = ROOT.TCanvas("c","c: pT",800,850)
  c.cd()
  pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
  pad1.SetBottomMargin(0)
  pad1.SetLeftMargin(0.15)
  pad1.SetRightMargin(0.05)
  pad1.SetTopMargin(0.05)
  pad1.SetLogy()
  pad1.Draw()
  pad1.cd()
  
  h.SetLineColor(1)
  h.SetLineWidth(2)
  h.SetLineStyle(1)
  
  h.Scale(1./nEvents, scalingOptions)
  h.GetYaxis().SetTitle(yAxisTitle)
  h.GetYaxis().SetTitleSize(0.06)
  h.GetXaxis().SetRangeUser(xRangeMin, xRangeMax)
  h.GetYaxis().SetRangeUser(2e-10,2e-3)
  h.GetYaxis().SetLabelFont(43)
  h.GetYaxis().SetLabelSize(20)
  xAxisTitle = h.GetXaxis().GetTitle()
  h.GetXaxis().SetTitle("")
  
  h2.SetLineColor(4)
  h2.SetLineWidth(2)
  h2.SetLineStyle(1)
  
  h.Draw("hist E")
  h2.Draw("hist same E")
  
  if h3:
    h3.SetLineColor(2)
    h3.SetLineWidth(2)
    h3.SetLineStyle(1)
    h3.Scale(1./nEvents, scalingOptions)
    h3.Draw("hist same")

  c.cd()
  pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
  pad2.SetTopMargin(0)
  pad2.SetBottomMargin(0.35)
  pad2.SetLeftMargin(0.15)
  pad2.SetRightMargin(0.05)
  pad2.Draw()
  pad2.cd()

  # plot ratio h/h2
  hRatio = h.Clone()
  hRatio.Divide(h2)
  hRatio.SetMarkerStyle(21)
  hRatio.SetMarkerColor(1)
  
  hRatio.GetXaxis().SetRangeUser(xRangeMin,xRangeMax)
  hRatio.GetXaxis().SetTitleSize(30)
  hRatio.GetXaxis().SetTitleFont(43)
  hRatio.GetXaxis().SetTitleOffset(4.)
  hRatio.GetXaxis().SetLabelFont(43)
  hRatio.GetXaxis().SetLabelSize(20)
  hRatio.GetXaxis().SetTitle(xAxisTitle)
  
  hRatio.GetYaxis().SetTitle(ratioYAxisTitle)
  hRatio.GetYaxis().SetTitleSize(20)
  hRatio.GetYaxis().SetTitleFont(43)
  hRatio.GetYaxis().SetTitleOffset(2.2)
  hRatio.GetYaxis().SetLabelFont(43)
  hRatio.GetYaxis().SetLabelSize(20)
  hRatio.GetYaxis().SetNdivisions(505)
  min= hRatio.GetBinContent(hRatio.GetMinimumBin())
  max= hRatio.GetBinContent(hRatio.GetMaximumBin())
  #automatic zoom-in for a very small scatter of the points
  if min>0.5 and max<1.5:
    hRatio.GetYaxis().SetRangeUser(0.5,1.5)
  elif yRatioMax>2:
    hRatio.GetYaxis().SetRangeUser(0,yRatioMax)
  else:
    hRatio.GetYaxis().SetRangeUser(2-yRatioMax,yRatioMax)
  
  hRatio.Draw("P E")
  
  # plot ratio h3/h2
  if h3:
    hRatio3 = h3.Clone()
    hRatio3.Divide(h2)
    hRatio3.SetMarkerStyle(21)
    hRatio3.SetMarkerColor(2)
    hRatio3.Draw("P E same")

  line = ROOT.TLine(xRangeMin,1,xRangeMax,1)
  line.SetLineColor(920+2)
  line.SetLineStyle(2)
  line.Draw()

  pad1.cd()

  if nEvents > 2:
    textNEvents = ROOT.TLatex()
    textNEvents.SetNDC()
    textNEvents.DrawLatex(0.55,0.6,"#it{N}_{events} = %d" % nEvents)

  leg2 = ROOT.TLegend(0.3,0.7,0.88,0.93,legendTitle)
  leg2.SetFillColor(10)
  leg2.SetBorderSize(0)
  leg2.SetFillStyle(0)
  leg2.SetTextSize(0.04)
  leg2.AddEntry(h, hLegendLabel, "l")
  if h3:
    leg2.AddEntry(h3, h3LegendLabel, "l")
  if h2:
    leg2.AddEntry(h2, h2LegendLabel, "l")
  leg2.Draw("same")
  
  c.SaveAs(outputFilename)
  c.Close()

###################################################################################################
# Plot basic histogram
###################################################################################################
def plotHist(h, outputFilename, drawOptions = "", setLogy = False, setLogz = False, text = ""):
  
  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd()
  c.cd().SetLeftMargin(0.15)
  
  if setLogy:
    c.SetLogy()
  if setLogz:
    c.SetLogz()
  h.DrawCopy(drawOptions)

  if text:
    textFit = ROOT.TLatex()
    textFit.SetTextSize(0.04)
    textFit.SetNDC()
    textFit.DrawLatex(0.6,0.8,text)

  c.SaveAs(outputFilename)
  c.Close()

#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description="Compare histograms to test the new EMCal corrections framework")
  parser.add_argument("-d", "--inputFileData", action="store",
                      type=str, metavar="inputFileData",
                      default="AnalysisResults.root",
                      help="Path of AnalysisResults.root file containing spectrum to be unfolded")
  parser.add_argument("-r", "--inputFileResponse", action="store",
                      type=str, metavar="inputFileResponse",
                      default="AnalysisResults.root",
                      help="Path of AnalysisResults.root file containing response matrix")
  parser.add_argument("-o", "--outputDir", action="store",
                      type=str, metavar="outputDir",
                      default="./unfolding_output/",
                      help="Output directory for QA plots to be written to")
  parser.add_argument("-i", "--imageFormat", action="store",
                      type=str, metavar="imageFormat",
                      default=".pdf",
                      help="Image format to save plots in, e.g. \".pdf\" or \".png\"")
                      
  # Parse the arguments
  args = parser.parse_args()
  
  print("Configuring...")
  print("inputFileData: \"{0}\"".format(args.inputFileData))
  print("inputFileResponse: \"{0}\"".format(args.inputFileResponse))
  print("ouputDir: \"{0}\"".format(args.outputDir))
  print("imageFormat: \"{0}\"".format(args.imageFormat))
  
  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFileData):
    print("File \"{0}\" does not exist! Exiting!".format(args.inputFileData))
    sys.exit(0)
  if not os.path.exists(args.inputFileResponse):
    print("File \"{0}\" does not exist! Exiting!".format(args.inputFileResponse))
    sys.exit(0)

  roounfold_inclusivejets(input_file_data = args.inputFileData, input_file_response = args.inputFileResponse, output_dir = args.outputDir, file_format = args.imageFormat)
