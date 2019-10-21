#! /usr/bin/env python

import sys
import os
import argparse
import itertools
from array import *
import numpy
import ROOT
import yaml

# Analysis utilities
#sys.path.append('/mnt/pyjetty')
#print(sys.path)
from pyjetty.alice_analysis.analysis.base import analysis_utils
from pyjetty.alice_analysis.analysis.base import analysis_base

# Load the RooUnfold library
#ROOT.gSystem.Load("$ALIBUILD_WORK_DIR/slc7_x86-64/RooUnfold/latest/lib/libRooUnfold.so")
ROOT.gSystem.Load('/Users/jamesmulligan/RooUnfold/build/libRooUnfold.dylib')

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

# Suppress a lot of standard output
ROOT.gErrorIgnoreLevel = ROOT.kWarning

################################################################
class roounfold_sd(analysis_base.analysis_base):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, observable='', input_file_data='', input_file_response='', config_file='', output_dir='', file_format='', rebin_response=False, truncation=False, binning=False, power_law_offset=0., **kwargs):
    
    super(roounfold_sd, self).__init__(input_file_data, input_file_response, config_file, output_dir, file_format, rebin_response, power_law_offset, **kwargs)
    
    self.fData = ROOT.TFile(self.input_file_data, 'READ')
    self.fResponse = ROOT.TFile(self.input_file_response, 'READ')
    self.observable = observable
    self.truncation = truncation
    self.binning = binning
    
    self.initialize_config()
  
    self.get_responses(self.rebin_response)
    
    # Create output files to store results
    for jetR in self.jetR_list:
      for sd_setting in self.sd_settings:
        
        zcut = sd_setting[0]
        beta = sd_setting[1]
        sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)
        
        fResult_name = os.path.join(self.output_dir, 'fResult_R{}_{}.root'.format(jetR, sd_label))
        setattr(self, 'fResult_name_R{}_{}'.format(jetR, sd_label), fResult_name)
        fResult = ROOT.TFile(fResult_name, 'RECREATE')
        fResult.Close()
    
    print(self)
  
  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def roounfold_sd(self):
    
    for jetR in self.jetR_list:
      for sd_setting in self.sd_settings:
        
        zcut = sd_setting[0]
        beta = sd_setting[1]
        sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)
      
        self.unfoldSingleOutputList(jetR, sd_label, zcut, beta)

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
    
    # Call base class initialization
    analysis_base.analysis_base.initialize_config(self)
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
      
      # Write tree output (default is to write only histograms)
      self.write_tree_output = config['write_tree_output']
      
      # Retrieve list of SD grooming settings
      self.jetR_list = config['jetR']
      self.sd_config_dict = config['SoftDrop']
      self.sd_config_list = list(self.sd_config_dict.keys())
      self.sd_settings = [[self.sd_config_dict[name]['zcut'], self.sd_config_dict[name]['beta']] for name in self.sd_config_list]
      
      if self.observable == 'theta_g':
        self.xtitle = '#theta_{g}'
        self.ytitle = '#frac{1}{#it{N}_{jets, inc}} #frac{d#it{N}}{d#theta_{g}}'
      if self.observable == 'zg':
        self.xtitle = '#it{z}_{g}'
        self.ytitle = '#frac{1}{#it{N}_{jets, inc}} #frac{d#it{N}}{d#it{z}_{g}}'
    
      # Retrieve histogram binnings for each SD setting
      for i, sd_setting in enumerate(self.sd_settings):
        
        zcut = sd_setting[0]
        beta = sd_setting[1]
        sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)
        config_name = self.sd_config_list[i]
        
        pt_det_bins_name = 'pt_bins_det'
        if self.truncation:
          pt_det_bins_name += '_sys_truncation'
        pt_bins_det = (self.sd_config_dict[config_name][pt_det_bins_name])
        pt_bins_truth = (self.sd_config_dict[config_name]['pt_bins_truth'])
        n_pt_bins_det = len(pt_bins_det) - 1
        setattr(self, 'n_pt_bins_det_{}'.format(sd_label), n_pt_bins_det)
        n_pt_bins_truth = len(pt_bins_truth) - 1
        setattr(self, 'n_pt_bins_truth_{}'.format(sd_label), n_pt_bins_truth)
        det_pt_bin_array = array('d',pt_bins_det)
        setattr(self, 'det_pt_bin_array_{}'.format(sd_label), det_pt_bin_array)
        truth_pt_bin_array = array('d',pt_bins_truth)
        setattr(self, 'truth_pt_bin_array_{}'.format(sd_label), truth_pt_bin_array)
          
        if self.observable == 'theta_g':
          det_bins_name = 'rg_bins_det'
          if self.binning:
            det_bins_name += '_sys_binning'
          rg_bins_det = (self.sd_config_dict[config_name][det_bins_name])
          rg_bins_truth = (self.sd_config_dict[config_name]['rg_bins_truth'])
          n_rg_bins_det = len(rg_bins_det) - 1
          setattr(self, 'n_bins_det_{}'.format(sd_label), n_rg_bins_det)
          n_rg_bins_truth = len(rg_bins_truth) - 1
          setattr(self, 'n_bins_truth_{}'.format(sd_label), n_rg_bins_truth)
          det_rg_bin_array = array('d',rg_bins_det)
          setattr(self, 'det_bin_array_{}'.format(sd_label), det_rg_bin_array)
          truth_rg_bin_array = array('d',rg_bins_truth)
          setattr(self, 'truth_bin_array_{}'.format(sd_label), truth_rg_bin_array)
        if self.observable == 'zg':
          det_bins_name = 'zg_bins_det'
          if self.binning:
            det_bins_name += '_sys_binning'
          zg_bins_det = (self.sd_config_dict[config_name][det_bins_name])
          zg_bins_truth = (self.sd_config_dict[config_name]['zg_bins_truth'])
          n_zg_bins_det = len(zg_bins_det) - 1
          setattr(self, 'n_bins_det_{}'.format(sd_label), n_zg_bins_det)
          n_zg_bins_truth = len(zg_bins_truth) - 1
          setattr(self, 'n_bins_truth_{}'.format(sd_label), n_zg_bins_truth)
          det_zg_bin_array = array('d',zg_bins_det)
          setattr(self, 'det_bin_array_{}'.format(sd_label), det_zg_bin_array)
          truth_zg_bin_array = array('d',zg_bins_truth)
          setattr(self, 'truth_bin_array_{}'.format(sd_label), truth_zg_bin_array)
        
        for jetR in self.jetR_list:

          if self.observable == 'theta_g':
            name_thn = 'hResponse_JetPt_ThetaG_R{}_{}Scaled'.format(jetR, sd_label)
            name_thn_rebinned = 'hResponse_JetPt_ThetaG_R{}_{}_rebinned'.format(jetR, sd_label)
            name_data = 'hThetaG_JetPt_R{}_{}'.format(jetR, sd_label)
            name_data_rebinned = 'hThetaG_JetPt_R{}_{}_rebinned'.format(jetR, sd_label)
          if self.observable == 'zg':
            name_thn = 'hResponse_JetPt_zg_R{}_{}Scaled'.format(jetR, sd_label)
            name_thn_rebinned = 'hResponse_JetPt_zg_R{}_{}_rebinned'.format(jetR, sd_label)
            name_data = 'hZg_JetPt_R{}_{}'.format(jetR, sd_label)
            name_data_rebinned = 'hZg_JetPt_R{}_{}_rebinned'.format(jetR, sd_label)
          name_roounfold = 'roounfold_response_R{}_{}'.format(jetR, sd_label)
          setattr(self, 'name_thn_R{}_{}'.format(jetR, sd_label), name_thn)
          setattr(self, 'name_thn_rebinned_R{}_{}'.format(jetR, sd_label), name_thn_rebinned)
          setattr(self, 'name_data_R{}_{}'.format(jetR, sd_label), name_data)
          setattr(self, 'name_data_rebinned_R{}_{}'.format(jetR, sd_label), name_data_rebinned)
          setattr(self, 'name_roounfold_R{}_{}'.format(jetR, sd_label), name_roounfold)
            
      self.min_pt_reported = 20
      self.max_pt_reported = 80
      self.regularizationParamName = 'n_iter'
      self.errorType = ROOT.RooUnfold.kCovToy

  #---------------------------------------------------------------
  # Get responses, either from file or manually rebin
  #---------------------------------------------------------------
  def get_responses(self, rebin_response=False):
    
    response_file_name = os.path.join(self.output_dir, 'response.root')
    if rebin_response:
      f = ROOT.TFile(response_file_name, 'RECREATE')
      f.Close()
        
    # Rebin response matrix, and create RooUnfoldResponse object
    # THn response matrix is: (pt-det, pt-true, obs-det, obs-true)
    for jetR in self.jetR_list:
      for sd_setting in self.sd_settings:
        
        zcut = sd_setting[0]
        beta = sd_setting[1]
        sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)
        
        name_thn = getattr(self, 'name_thn_R{}_{}'.format(jetR, sd_label))
        name_thn_rebinned = getattr(self, 'name_thn_rebinned_R{}_{}'.format(jetR, sd_label))
        name_data = getattr(self, 'name_data_R{}_{}'.format(jetR, sd_label))
        name_roounfold = getattr(self, 'name_roounfold_R{}_{}'.format(jetR, sd_label))
        
        # Retrieve desired binnings
        n_pt_bins_det = getattr(self, 'n_pt_bins_det_{}'.format(sd_label))
        det_pt_bin_array = getattr(self, 'det_pt_bin_array_{}'.format(sd_label))
        n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(sd_label))
        truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(sd_label))

        n_bins_det = getattr(self, 'n_bins_det_{}'.format(sd_label))
        det_bin_array = getattr(self, 'det_bin_array_{}'.format(sd_label))
        n_bins_truth = getattr(self, 'n_bins_truth_{}'.format(sd_label))
        truth_bin_array = getattr(self, 'truth_bin_array_{}'.format(sd_label))
        
        # Rebin if requested, and write to file
        thn = self.fResponse.Get(name_thn)
        thn.SetName(name_thn)
        setattr(self, name_thn, thn)
        if rebin_response:
          # Create rebinned THn and RooUnfoldResponse with these binnings, and write to file
          
          if self.write_tree_output:
            tree_file_name = '/Users/jamesmulligan/alidock/theta_g/rganalysis_embed_PbPb/output_alpha_0_dRmax_0.25_SDzcut_0.2_emb.root'
            
            self.utils.construct_response_from_ntuple(response_file_name, name_thn_rebinned, name_roounfold, jetR, beta, n_pt_bins_det, det_pt_bin_array, n_bins_det, det_bin_array, n_pt_bins_truth, truth_pt_bin_array, n_bins_truth, truth_bin_array, self.power_law_offset)

          else:
            self.utils.rebin_response(response_file_name, thn, name_thn_rebinned, name_roounfold, jetR, sd_label, n_pt_bins_det, det_pt_bin_array, n_bins_det, det_bin_array, n_pt_bins_truth, truth_pt_bin_array, n_bins_truth, truth_bin_array, observable, self.power_law_offset)
          
        # Also re-bin the data histogram
        hData = self.fData.Get(name_data)
        h = self.utils.rebin_data(hData, name_data, n_pt_bins_det, det_pt_bin_array, n_bins_det, det_bin_array)
        h.SetDirectory(0)
        name = '{}_{}'.format(name_data, 'rebinned')
        setattr(self, name, h)

        # Retrieve responses from file
        f = ROOT.TFile(response_file_name, 'READ')
        thn_rebinned = f.Get(name_thn_rebinned)
        #thn_rebinned.SetDirectory(0)
        roounfold_response = f.Get(name_roounfold)
        roounfold_response.UseOverflow(False)

        setattr(self, name_thn_rebinned, thn_rebinned)
        setattr(self, name_roounfold, roounfold_response)
        f.Close()

  ###################################################################################################
  # Unfold jet spectrum from a single output list
  ###################################################################################################
  def unfoldSingleOutputList(self, jetR, sd_label, zcut, beta):

    print('jetR = {}, {}'.format(jetR, sd_label))
    
    # Get data jet spectrum
    name = getattr(self, 'name_data_rebinned_R{}_{}'.format(jetR, sd_label))
    hData_PerBin = getattr(self, name)
    hData_PerBin.GetXaxis().SetRangeUser(0., 100.)
    outputFilename = os.path.join(self.output_dir, 'hData_R{}_{}{}'.format(self.utils.remove_periods(jetR), sd_label, self.file_format))
    self.utils.plot_hist(hData_PerBin, outputFilename, 'colz', False, True)
    
    hDataProj = hData_PerBin.ProjectionY()
    hDataProj.SetName('{}_{}'.format(hDataProj.GetName(), 'proj'))
    outputFilename = os.path.join(self.output_dir, 'hData_proj_R{}_{}{}'.format(self.utils.remove_periods(jetR), sd_label, self.file_format))
    self.utils.plot_hist(hDataProj, outputFilename)
    
    # Plot various slices of the response matrix (from the THn)
    self.plot_RM_slices(jetR, sd_label)
    
    # Plot the kinematic efficiency from the response THn
    self.plot_kinematic_efficiency(jetR, sd_label, zcut, beta)
    # Can either pass this to the RooUnfoldResponse object, or else can apply it afterward
    # (Need to think...)
    
    # Get MC-det and MC-truth 2D projections for unfolding closure test
    name = 'hMC_Det_R{}_{}'.format(jetR, sd_label)
    hMC_Det = self.get_MCdet2D(jetR, sd_label)
    setattr(self, name, hMC_Det)
    
    name = 'hMC_Truth_R{}_{}'.format(jetR, sd_label)
    hMC_Truth = self.get_MCtruth2D(jetR, sd_label)
    setattr(self, name, hMC_Truth)
    
    # Compute SD tagging rate
    self.compute_tagging_rate(jetR, sd_label)

    # Unfold spectrum
    if hData_PerBin and hMC_Det and hMC_Truth:
      
      self.unfoldJetSpectrum(jetR, sd_label, zcut, beta)

  #################################################################################################
  # Compute SD tagging rate, based on MC correction
  #################################################################################################
  def compute_tagging_rate(self, jetR, sd_label):
    
    # Get truth binnings, since we need to rebin the det-level histograms
    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(sd_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(sd_label))
    n_bins_truth = getattr(self, 'n_bins_truth_{}'.format(sd_label))
    truth_bin_array = getattr(self, 'truth_bin_array_{}'.format(sd_label))

    # Create a histogram to store the tagging fractions
    name = 'hTaggingFractions_R{}_{}'.format(jetR, sd_label)
    h = ROOT.TH1F(name, name, n_pt_bins_truth, truth_pt_bin_array)
    
    # Get relevant histograms
    name = getattr(self, 'name_data_rebinned_R{}_{}'.format(jetR, sd_label))
    hData2D = getattr(self, name).Clone()
    hData2D.SetName('{}_tagging'.format(hData2D.GetName()))
    hData2D_rebinned = self.utils.rebin_data(hData2D, name, n_pt_bins_truth, truth_pt_bin_array, n_bins_truth, truth_bin_array)
    
    name = 'hMC_Det_R{}_{}'.format(jetR, sd_label)
    hMC_Det = getattr(self, name).Clone()
    hMC_Det.SetName('{}_tagging'.format(hMC_Det.GetName()))
    hMC_Det_rebinned = self.utils.rebin_data(hMC_Det, name, n_pt_bins_truth, truth_pt_bin_array, n_bins_truth, truth_bin_array)
    
    name = 'hMC_Truth_R{}_{}'.format(jetR, sd_label)
    hMC_Truth = getattr(self, name).Clone()
    hMC_Truth.SetName('{}_tagging'.format(hMC_Truth.GetName()))
    
    # Compute the tagging fraction for each pt bin
    for bin in range(1, n_pt_bins_truth-3):
      min_pt_truth = truth_pt_bin_array[bin]
      max_pt_truth = truth_pt_bin_array[bin+1]
    
      fraction_tagged = self.get_tagging_rate(jetR, min_pt_truth, max_pt_truth, hData2D_rebinned, hMC_Det_rebinned, hMC_Truth)

      x = (min_pt_truth + max_pt_truth)/2.
      bin = h.FindBin(x)
      h.SetBinContent(bin, fraction_tagged)

    fResult_name = getattr(self, 'fResult_name_R{}_{}'.format(jetR, sd_label))
    fResult = ROOT.TFile(fResult_name, 'UPDATE')
    h.Write()
    fResult.Close()

  #################################################################################################
  # Compute SD tagging rate, based on MC correction
  #################################################################################################
  def get_tagging_rate(self, jetR, min_pt_truth, max_pt_truth, hData2D, hMC_Det2D, hMC_Truth2D):

    hData2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    hData = hData2D.ProjectionY()
    n_jets_inclusive = hData.Integral(0, hData.GetNbinsX()+1)
    n_jets_tagged = hData.Integral(1, hData.GetNbinsX())
    fraction_tagged_data =  n_jets_tagged/n_jets_inclusive
    #print('fraction_tagged_data: {}'.format(fraction_tagged_data))

    hMC_Det2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    hMC_Det = hMC_Det2D.ProjectionY()
    n_jets_inclusive = hMC_Det.Integral(0, hMC_Det.GetNbinsX()+1)
    n_jets_tagged = hMC_Det.Integral(1, hMC_Det.GetNbinsX())
    fraction_tagged_mc_det =  n_jets_tagged/n_jets_inclusive
    #print('fraction_tagged_mc_det: {}'.format(fraction_tagged_mc_det))
    
    hMC_Truth2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    hMC_Truth = hMC_Truth2D.ProjectionY()
    n_jets_inclusive = hMC_Truth.Integral(0, hMC_Truth.GetNbinsX()+1)
    n_jets_tagged = hMC_Truth.Integral(1, hMC_Truth.GetNbinsX())
    fraction_tagged_mc_truth =  n_jets_tagged/n_jets_inclusive
    #print('fraction_tagged_mc_truth: {}'.format(fraction_tagged_mc_truth))

    fraction_tagged = fraction_tagged_data * fraction_tagged_mc_truth / fraction_tagged_mc_det
    #print('fraction_tagged: {}'.format(fraction_tagged))

    return fraction_tagged

  #################################################################################################
  # Unfold jet spectrum
  #################################################################################################
  def unfoldJetSpectrum(self, jetR, sd_label, zcut, beta):
    
    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()
    
    # Create canvas to superpose iterations on one plot, to examine convergence
    c1 = ROOT.TCanvas('c1_{}_{}'.format(jetR, sd_label),'c1_{}_{}: histos'.format(jetR, sd_label),600,450)
    c1.cd()
    c1.SetLogy()
    ROOT.gPad.SetLeftMargin(0.15)
    leg = ROOT.TLegend(0.6,0.5,0.88,0.83,"{} Unfolding".format(type))
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    
    # Loop over values of regularization parameter
    response = getattr(self, 'roounfold_response_R{}_{}'.format(jetR, sd_label))

    name_data = getattr(self, 'name_data_rebinned_R{}_{}'.format(jetR, sd_label))
    hData = getattr(self, name_data)

    fResult_name = getattr(self, 'fResult_name_R{}_{}'.format(jetR, sd_label))
    fResult = ROOT.TFile(fResult_name, 'UPDATE')
    
    reg_param_final = self.utils.get_reg_param(self.sd_settings, self.sd_config_list, self.sd_config_dict, sd_label, self.observable, jetR)
    
    for i in range(1, reg_param_final + 3):
      
      # Set up the Bayesian unfolding object
      unfoldBayes = ROOT.RooUnfoldBayes(response, hData, i)
      #unfoldBayes.SetNToys(1000)
      
      # Perform the unfolding
      print('Unfolding with {} = {}'.format(self.regularizationParamName, i))
      hUnfolded = unfoldBayes.Hreco(self.errorType) # Produces the truth distribution, with errors, PerBin (will scale by bin width below, after refolding checks)
      
      # Save unfolded solution as class member
      name = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, sd_label, i)
      hUnfolded.SetName(name)
      setattr(self, name, hUnfolded)
      hUnfolded.SetDirectory(0)

      # Correct by kinematic efficiency
      hKinematicEfficiency = getattr(self, 'hKinematicEfficiency_R{}_{}'.format(jetR, sd_label))
      hUnfolded.Divide(hKinematicEfficiency)
        
      # Write result to file
      hUnfolded.Write()
      
      # Plot Pearson correlation coeffs for each k, to get a measure of the correlation between the bins
      covarianceMatrix = unfoldBayes.Ereco(self.errorType) # Get the covariance matrix
      #plotCorrelationCoefficients(covarianceMatrix, i, output_dir, file_format)
      
    fResult.Close()
    print('Done unfolding')
      
    self.plot_unfolded_observable(jetR, sd_label, zcut, beta)

    self.plot_unfolded_pt(jetR, sd_label, zcut, beta)

    #--------------------------------------------------------------
    
    # Smear MC truth spectrum according to the error bars on the measured spectrum, for closure test
    
    name_data = getattr(self, 'name_data_rebinned_R{}_{}'.format(jetR, sd_label))
    hData = getattr(self, name_data)
    hMC_Det = getattr(self, 'hMC_Det_R{}_{}'.format(jetR, sd_label))

    measuredErrors = self.getMeasuredErrors(hData)
    self.smearSpectrum(hMC_Det, measuredErrors)

    # Loop over values of regularization parameter to do unfolding checks
    for i in range(1, reg_param_final + 3):
      
      # Apply RM to unfolded result, and check that I obtain measured spectrum (simple technical check)
      self.plot_refolding_test(i, jetR, sd_label, zcut, beta)

      # Unfold the smeared det-level result with response, and compare to truth-level MC.
      self.plot_closure_test(i, jetR, sd_label, zcut, beta)

  #################################################################################################
  # Plot various slices of the response matrix (from the THn)
  #################################################################################################
  def plotRegParamSystematic(self):
    
    #Plot the spectra comparing only k=+1 and k-1 to the main result
    xRangeMin = min_pt_reported
    xRangeMax = max_pt_reported
    yAxisTitle = "#frac{d#sigma}{dp_{T}} [mb (GeV/c)^{-1}]"
    ratioYAxisTitle = "Ratio to k={}".format(reg_param)
    outputFilename = os.path.join(output_dir, "hJetSpectraUnfoldedRatio" + file_format)
    legendTitle = "{} Unfolding".format(type)
    hLegendLabel = "k = {}".format(reg_param-1)
    h2LegendLabel = "k = {}".format(reg_param)
    h3LegendLabel = "k = {}".format(reg_param+1)
    if hLowerkResult and hMainResult and hHigherkResult:
      # To get sensible error bars, assume main result has no errors
      for bin in range(1, hMainResult.GetNbinsX() + 1):
        hMainResult.SetBinError(bin, 0)
      plotSpectra(hLowerkResult, hMainResult, hHigherkResult, 1., xRangeMin, xRangeMax, yAxisTitle, ratioYAxisTitle, outputFilename, "", legendTitle, hLegendLabel, h2LegendLabel, h3LegendLabel)

  #################################################################################################
  # Plot various slices of the response matrix (from the THn)
  #################################################################################################
  def plot_unfolded_observable(self, jetR, sd_label, zcut, beta):
    
    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(sd_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(sd_label))

    for bin in range(1, n_pt_bins_truth-3):
      min_pt_truth = truth_pt_bin_array[bin]
      max_pt_truth = truth_pt_bin_array[bin+1]
      
      self.plot_observable(jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth)
      self.plot_observable(jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth, option = 'ratio')
      self.plot_observable(jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth, option = 'uncertainties')

  #################################################################################################
  # Plot various slices of the response matrix (from the THn)
  #################################################################################################
  def plot_observable(self, jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth, option = ''):

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()
    
    name = 'cResult{}_R{}_{}_{}-{}'.format(option, jetR, sd_label, min_pt_truth, max_pt_truth)
    c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()
    
    c.cd()
    myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    myPad.SetLeftMargin(0.2)
    myPad.SetTopMargin(0.07)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.13)
    myPad.Draw()
    myPad.cd()
    
    n_bins_truth = getattr(self, 'n_bins_truth_{}'.format(sd_label))
    truth_bin_array = getattr(self, 'truth_bin_array_{}'.format(sd_label))
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_bins_truth, truth_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle(self.xtitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle(self.ytitle)
    myBlankHisto.SetMaximum(4)
    myBlankHisto.SetMinimum(0.)
    if option == 'ratio':
      myBlankHisto.SetMaximum(1.2)
      myBlankHisto.SetMinimum(0.9)
      myBlankHisto.SetYTitle('#frac{n_{iter}}{(n_{iter}-1)}')
    if option == 'uncertainties':
      myBlankHisto.SetMaximum(40)
      myBlankHisto.SetMinimum(0.)
      myBlankHisto.SetYTitle('statistical uncertainty (%)')
    myBlankHisto.Draw("E")
    
    leg = ROOT.TLegend(0.75,0.65,0.88,0.92)
    self.utils.setup_legend(leg,0.04)

    reg_param_final = self.utils.get_reg_param(self.sd_settings, self.sd_config_list, self.sd_config_dict, sd_label, self.observable, jetR)
    
    for i in range(1, reg_param_final + 3):
      
      h = self.get_unfolded_result(jetR, sd_label, i, min_pt_truth, max_pt_truth, option)
      
      if option == 'ratio':
        if i > 1:
          h_previous = self.get_unfolded_result(jetR, sd_label, i-1, min_pt_truth, max_pt_truth, option)
          h.Divide(h_previous)
        else:
          continue
    
      if option == 'uncertainties':
        h = self.get_unfolded_result_uncertainties(jetR, sd_label, i, min_pt_truth, max_pt_truth, option)
      
      if i == 1:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(20)
      elif i == 2:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(21)
      elif i == 3:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(22)
      elif i == 4:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(23)
      elif i == 5:
        h.SetMarkerSize(2)
        h.SetMarkerStyle(33)
      elif i == 6:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(34)
      elif i == 7:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(35)
      elif i == 8:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(36)
      else:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(19)
      h.SetMarkerColor(600-6+i)
      h.SetLineStyle(1)
      h.SetLineWidth(2)
      h.SetLineColor(600-6+i)

      # Doesn't work for some reason...
      #shift = 0.5*h.GetBinWidth(1)
      #h.GetXaxis().SetLimits(truth_bin_array[0]+shift,truth_bin_array[-1]+shift)

      if option == 'ratio':
        h.DrawCopy('P hist X0 same')
      else:
        h.DrawCopy('PE X0 same')
      
      label = '{} = {}'.format(self.regularizationParamName, i)
      leg.AddEntry(h, label, 'Pe')

    leg.Draw()
    
    if option == 'ratio':
      line = ROOT.TLine(truth_bin_array[0], 1, truth_bin_array[-1], 1)
      line.SetLineColor(920+2)
      line.SetLineStyle(2)
      line.SetLineWidth(4)
      line.Draw()

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth)
    text_latex.DrawLatex(0.25, 0.85, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'R = ' + str(jetR) + '   z_{cut} = ' + str(zcut) + '   #beta = ' + str(beta)
    text_latex.DrawLatex(0.25, 0.75, text)
    
    outputFilename = os.path.join(self.output_dir, 'hUnfolded{}_{}_R{}_{}_{}-{}{}'.format(option, self.observable, self.utils.remove_periods(jetR), sd_label, int(min_pt_truth), int(max_pt_truth), self.file_format))
    c.SaveAs(outputFilename)
    c.Close()

  #################################################################################################
  # Get unfolded result in 1D, for fixed slice of pt
  #################################################################################################
  def get_unfolded_result(self, jetR, sd_label, i, min_pt_truth, max_pt_truth, option):

    name = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, sd_label, i)
    h2D = getattr(self, name)
    h2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    h = h2D.ProjectionY()
    h.SetName('{}_{}'.format(h.GetName(), option))
    
    # Normalize by integral, i.e. N_jets,inclusive in this pt-bin
    n_jets_inclusive = h.Integral(0, h.GetNbinsX()+1)
    h.Scale(1./n_jets_inclusive, 'width')
    
    return h
  
  #################################################################################################
  # Get unfolded result uncertainties in 1D, for fixed slice of pt
  #################################################################################################
  def get_unfolded_result_uncertainties(self, jetR, sd_label, i, min_pt_truth, max_pt_truth, option):
  
    h = self.get_unfolded_result(jetR, sd_label, i, min_pt_truth, max_pt_truth, option)
  
    for bin in range(1, h.GetNbinsX()+1):
      content = h.GetBinContent(bin)
      uncertainty = h.GetBinError(bin)
      h.SetBinContent(bin, uncertainty/content * 100)
      h.SetBinError(bin, 0)

    return h

  #################################################################################################
  # Plot various slices of the response matrix (from the THn)
  #################################################################################################
  def plot_unfolded_pt(self, jetR, sd_label, zcut, beta):
    
    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()
    
    name = 'cResultPt_R{}_{}'.format(jetR, sd_label)
    c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()
    
    c.cd()
    myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    myPad.SetLeftMargin(0.2)
    myPad.SetTopMargin(0.07)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.13)
    myPad.Draw()
    myPad.cd()
    
    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(sd_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(sd_label))
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_pt_bins_truth, truth_pt_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle('#it{p}_{T, ch jet}')
    myBlankHisto.GetYaxis().SetTitleOffset(2.2)
    myBlankHisto.SetYTitle('#frac{dN}{d#it{p}_{T, ch jet}}')
    myBlankHisto.SetMaximum(5000)
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")
    
    leg = ROOT.TLegend(0.75,0.65,0.88,0.92)
    self.utils.setup_legend(leg,0.04)
    
    reg_param_final = self.utils.get_reg_param(self.sd_settings, self.sd_config_list, self.sd_config_dict, sd_label, self.observable, jetR)
    
    for i in range(1, reg_param_final + 3):
      
      name = 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, sd_label, i)
      h2D = getattr(self, name)
      h2D.GetXaxis().SetRangeUser(5., 120.)
      h = h2D.ProjectionX()
      
      h.Scale(1., 'width')
      
      if i == 1:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(20)
      elif i == 2:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(21)
      elif i == 3:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(22)
      elif i == 4:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(23)
      elif i == 5:
        h.SetMarkerSize(2)
        h.SetMarkerStyle(33)
      elif i == 6:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(34)
      elif i == 7:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(35)
      elif i == 8:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(36)
      else:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(19)
      h.SetMarkerColor(600-6+i)
      h.SetLineStyle(1)
      h.SetLineWidth(2)
      h.SetLineColor(600-6+i)
      
      h.DrawCopy('PE X0 same')
      
      label = '{} = {}'.format(self.regularizationParamName, i)
      leg.AddEntry(h, label, 'Pe')
    
    leg.Draw()

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'R = ' + str(jetR) + '   z_{cut} = ' + str(zcut) + '   #beta = ' + str(beta)
    text_latex.DrawLatex(0.25, 0.75, text)

    outputFilename = os.path.join(self.output_dir, 'hUnfoldedPt_{}_R{}_{}{}'.format(self.observable, self.utils.remove_periods(jetR), sd_label, self.file_format))
    c.SaveAs(outputFilename)
    c.Close()

  ###################################################################################################
  # Plot kinematic efficiency
  # The kinematic efficiency is the ratio:
  #   Numerator: 2D truth-level projection [pt-true, sd-obs-true] using no cut on det-level
  #   Denominator: 2D truth-level projection [pt-true, sd-obs-true] using [pt-det, sd-obs-det] cut on det-level
  ###################################################################################################
  def plot_kinematic_efficiency(self, jetR, sd_label, zcut, beta):
    
    # (pt-det, pt-true, obs-det, obs-true)
    name_response = getattr(self, 'name_thn_rebinned_R{}_{}'.format(jetR, sd_label))
    response = getattr(self,  name_response)
    hResponse = response.Clone()
    hResponse.SetName('hResponse_KinEff_JetPt_Obs_R{}_{}'.format(jetR, sd_label))
    
    # Denominator -- by default, under/over-flow bins are included in projection
    hDenominator = hResponse.Projection(3, 1)
    hDenominator.SetName('{}_Denominator'.format(hDenominator.GetName()))
    
    # Numerator -- cut on det-level input binning
    det_pt_bin_array = getattr(self, 'det_pt_bin_array_{}'.format(sd_label))
    min_pt_det = det_pt_bin_array[0]
    max_pt_det = det_pt_bin_array[-1]
    det_bin_array = getattr(self, 'det_bin_array_{}'.format(sd_label))
    min_obs_det = det_bin_array[0]
    max_obs_det = det_bin_array[-1]
    hResponse.GetAxis(0).SetRangeUser(min_pt_det, max_pt_det)
    hResponse.GetAxis(2).SetRangeUser(min_obs_det, max_obs_det)
    hNumerator = hResponse.Projection(3, 1)
    hNumerator.SetName('{}_Numerator'.format(hNumerator.GetName()))

    hKinematicEfficiency = hNumerator.Clone()
    hKinematicEfficiency.SetName('hKinematicEfficiency_R{}_{}'.format(jetR, sd_label))
    hKinematicEfficiency.Divide(hDenominator)
    for bin in range(0, hKinematicEfficiency.GetNcells()+1):
      hKinematicEfficiency.SetBinError(bin, 0)
    outputFilename = os.path.join(self.output_dir, 'hKinematicEfficiency2D_R{}_{}{}'.format(self.utils.remove_periods(jetR), sd_label, self.file_format))
    self.utils.plot_hist(hKinematicEfficiency, outputFilename, 'colz')
    
    # Save kinematic efficiency as class member
    setattr(self, 'hKinematicEfficiency_R{}_{}'.format(jetR, sd_label), hKinematicEfficiency)
    
    # Plot 1D kinematic efficiency
    self.plot_kinematic_efficiency_projections(hKinematicEfficiency, jetR, sd_label, zcut, beta)

  #################################################################################################
  # Unfold jet spectrum
  #################################################################################################
  def plot_kinematic_efficiency_projections(self, hKinematicEfficiency2D, jetR, sd_label, zcut, beta):
    
    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()
    
    name = 'cKinEff_R{}_{}'.format(jetR, sd_label)
    c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()
    
    c.cd()
    myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    myPad.SetLeftMargin(0.2)
    myPad.SetTopMargin(0.07)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.13)
    myPad.Draw()
    myPad.cd()
    
    n_bins_truth = getattr(self, 'n_bins_truth_{}'.format(sd_label))
    truth_bin_array = getattr(self, 'truth_bin_array_{}'.format(sd_label))
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_bins_truth, truth_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle(self.xtitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle('#varepsilon_{kin}')
    myBlankHisto.SetMaximum(2)
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")
    
    leg = ROOT.TLegend(0.6,0.65,0.72,0.92)
    self.utils.setup_legend(leg,0.04)
    
    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(sd_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(sd_label))
    for bin in range(1, n_pt_bins_truth-3):
      min_pt_truth = truth_pt_bin_array[bin]
      max_pt_truth = truth_pt_bin_array[bin+1]
      
      hKinematicEfficiency2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
      h = hKinematicEfficiency2D.ProjectionY()
      name = 'hKinematicEfficiency_R{}_{}_{}-{}'.format(jetR, sd_label, min_pt_truth, max_pt_truth)
      h.SetName(name)
      
      if bin == 1:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(20)
      elif bin == 2:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(21)
      elif bin == 3:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(22)
      elif bin == 4:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(23)
      elif bin == 5:
        h.SetMarkerSize(2)
        h.SetMarkerStyle(33)
      else:
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(19)
      h.SetMarkerColor(600-6+bin)
      h.SetLineStyle(1)
      h.SetLineWidth(2)
      h.SetLineColor(600-6+bin)
      
      h.DrawCopy('P X0 same')
      
      label = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
      leg.AddEntry(h, label, 'Pe')
    
    leg.Draw()

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'R = ' + str(jetR)
    text_latex.DrawLatex(0.3, 0.85, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'z_{cut} = ' + str(zcut) + '   #beta = ' + str(beta)
    text_latex.DrawLatex(0.3, 0.75, text)

    outputFilename = os.path.join(self.output_dir, 'hKinematicEfficiency_R{}_{}{}'.format(self.utils.remove_periods(jetR), sd_label, self.file_format))
    c.SaveAs(outputFilename)
    c.Close()

  #################################################################################################
  # Plot various slices of the response matrix (from the THn)
  #################################################################################################
  def plot_RM_slices(self, jetR, sd_label):
    
    output_dir_RM = os.path.join(self.output_dir, 'RM')
    if not os.path.isdir(output_dir_RM):
      os.makedirs(output_dir_RM)
    
    # (pt-det, pt-true, obs-det, obs-true)
    name_response = getattr(self, 'name_thn_R{}_{}'.format(jetR, sd_label))
    hResponse = getattr(self, name_response)

    # Fix pt-true, and plot the 2D sd-observable response
    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(sd_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(sd_label))
    for bin in range(1, n_pt_bins_truth-1):
      min_pt_truth = truth_pt_bin_array[bin]
      max_pt_truth = truth_pt_bin_array[bin+1]
      
      self.plot_Obs_Response(jetR, sd_label, min_pt_truth, max_pt_truth, hResponse, output_dir_RM)

  #################################################################################################
  # Plot 2D SD-observable response for a fixed range of pt-truth
  #################################################################################################
  def plot_Obs_Response(self, jetR, sd_label, min_pt_truth, max_pt_truth, hResponse, output_dir_RM):
    
    hResponse4D = hResponse.Clone()
    hResponse4D.SetName('{}_{}_{}'.format(hResponse4D.GetName(), min_pt_truth, max_pt_truth))
    
    hResponse4D.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)
    
    truth_bin_array = getattr(self, 'truth_bin_array_{}'.format(sd_label))
    hResponse4D.GetAxis(2).SetRangeUser(truth_bin_array[0], truth_bin_array[-1])
    hResponse4D.GetAxis(3).SetRangeUser(truth_bin_array[0], truth_bin_array[-1])

    hResponse_Obs = hResponse4D.Projection(3,2)
    hResponse_Obs.SetName('hResponse_Obs_R{}_{}_{}_{}'.format(jetR, sd_label, min_pt_truth, max_pt_truth))
    
    hResponse_Obs_Normalized = self.normalizeResponseMatrix(hResponse_Obs)
    
    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet}^{true} < ' + str(max_pt_truth)
    
    outputFilename = os.path.join(output_dir_RM, '{}{}'.format(hResponse_Obs.GetName(), self.file_format))
    self.utils.plot_hist(hResponse_Obs_Normalized, outputFilename, 'colz', False, True, text)

  ################################################################################################
  # Normalize response matrix
  # Normalize the truth projection to 1
  ################################################################################################
  def normalizeResponseMatrix(self, hResponseMatrix):
    
    # Make projection onto true axis (y-axis), and scale appropriately
    hTruthProjectionBefore = hResponseMatrix.ProjectionY('{}_py'.format(hResponseMatrix.GetName()),1,hResponseMatrix.GetNbinsX()) # Do exclude under and overflow bins
    hTruthProjectionBefore.SetName("hTruthProjectionBefore")
    
    # Loop through truth-level bins, and apply normalization factor to all bins.
    nBinsY = hResponseMatrix.GetNbinsY() # truth
    nBinsX = hResponseMatrix.GetNbinsX() # det
    for truthBin in range(1,nBinsY+1):
      normalizationFactor = hTruthProjectionBefore.GetBinContent(truthBin)
      if normalizationFactor > 0:
        truthBinCenter = hTruthProjectionBefore.GetXaxis().GetBinCenter(truthBin)
        
        for detBin in range(1,nBinsX+1):
          binContent = hResponseMatrix.GetBinContent(detBin, truthBin)
          hResponseMatrix.SetBinContent(detBin, truthBin, binContent/normalizationFactor)

    return hResponseMatrix

  #################################################################################################
  # Apply RM to unfolded result, and check that I obtain measured spectrum
  #################################################################################################
  def plot_refolding_test(self, i, jetR, sd_label, zcut, beta):
        
    output_dir_refolding = os.path.join(self.output_dir, 'Refolding')
    if not os.path.isdir(output_dir_refolding):
      os.makedirs(output_dir_refolding)
    
    hUnfolded = getattr(self, 'hUnfolded_{}_R{}_{}_{}'.format(self.observable, jetR, sd_label, i))
    response = getattr(self, 'roounfold_response_R{}_{}'.format(jetR, sd_label))
    
    hFoldedTruth = response.ApplyToTruth(hUnfolded) # Produces folded distribution PerBin (unfolded spectrum is also PerBin at the moment)
    
    n_pt_bins_det = getattr(self, 'n_pt_bins_det_{}'.format(sd_label))
    det_pt_bin_array = getattr(self, 'det_pt_bin_array_{}'.format(sd_label))
    for bin in range(2, n_pt_bins_det):
      min_pt_det = det_pt_bin_array[bin]
      max_pt_det = det_pt_bin_array[bin+1]
      
      self.plot_sd_refolded_slice(hFoldedTruth, i, jetR, sd_label, zcut, beta, min_pt_det, max_pt_det, output_dir_refolding)

    self.plot_pt_refolded_slice(hFoldedTruth, i, jetR, sd_label, zcut, beta, det_pt_bin_array[2], det_pt_bin_array[-1], output_dir_refolding)

  #################################################################################################
  # Apply RM to unfolded result, and check that I obtain measured spectrum (simple technical check)
  #################################################################################################
  def plot_sd_refolded_slice(self, hFoldedTruth, i, jetR, sd_label, zcut, beta, min_pt_det, max_pt_det, output_dir_refolding):
    
    hFoldedTruth.GetXaxis().SetRangeUser(min_pt_det, max_pt_det)
    hFolded_obs = hFoldedTruth.ProjectionY()
    hFolded_obs.SetName('hFolded_obs_R{}_{}_{}_{}-{}'.format(jetR, sd_label, i, min_pt_det, max_pt_det))
    
    name_data_rebinned = getattr(self, 'name_data_rebinned_R{}_{}'.format(jetR, sd_label))
    hData_PerBin = getattr(self, name_data_rebinned)
    hData_PerBin.GetXaxis().SetRangeUser(min_pt_det, max_pt_det)
    hData_obs = hData_PerBin.ProjectionY()
    hData_obs.SetName('hData_obs_R{}_{}_{}_{}-{}'.format(jetR, sd_label, i, min_pt_det, max_pt_det))
    
    legendTitle = ''
    h1LegendLabel = 'Folded truth, {} = {}'.format(self.regularizationParamName,i)
    h2LegendLabel = 'Measured pp'
    ratioYAxisTitle = 'Folded truth / Measured'
    outputFilename = os.path.join(output_dir_refolding, 'hFoldedTruth_R{}_{}_{}-{}_{}{}'.format(self.utils.remove_periods(jetR), sd_label, int(min_pt_det), int(max_pt_det), i, self.file_format))
    self.plot_sd_obs_ratio(hFolded_obs, hData_obs, None, self.ytitle, ratioYAxisTitle, int(min_pt_det), int(max_pt_det), jetR, sd_label, zcut, beta, outputFilename, 'width', legendTitle, h1LegendLabel, h2LegendLabel)

  #################################################################################################
  # Apply RM to unfolded result, and check that I obtain measured spectrum (simple technical check)
  #################################################################################################
  def plot_pt_refolded_slice(self, hFoldedTruth, i, jetR, sd_label, zcut, beta, min_pt_det, max_pt_det, output_dir_refolding):
    
    hFoldedTruth.GetXaxis().SetRangeUser(min_pt_det, max_pt_det)
    hFolded_pt = hFoldedTruth.ProjectionX()
    hFolded_pt.SetName('hFolded_pt_R{}_{}_{}_{}-{}'.format(jetR, sd_label, i, min_pt_det, max_pt_det))
    
    name_data_rebinned = getattr(self, 'name_data_rebinned_R{}_{}'.format(jetR, sd_label))
    hData_PerBin = getattr(self, name_data_rebinned)
    hData_PerBin.GetXaxis().SetRangeUser(min_pt_det, max_pt_det)
    hData_pt = hData_PerBin.ProjectionX()
    hData_pt.SetName('hData_pt_R{}_{}_{}_{}-{}'.format(jetR, sd_label, i, min_pt_det, max_pt_det))
    
    legendTitle = ''
    h1LegendLabel = 'Folded truth, {} = {}'.format(self.regularizationParamName,i)
    h2LegendLabel = 'Measured pp'
    ratioYAxisTitle = 'Folded truth / Measured'
    outputFilename = os.path.join(output_dir_refolding, 'hFoldedTruth_pt_R{}_{}_{}{}'.format(self.utils.remove_periods(jetR), sd_label, i, self.file_format))
    self.plot_sd_obs_ratio(hFolded_pt, hData_pt, None, self.ytitle, ratioYAxisTitle, 0, 0, jetR, sd_label, zcut, beta, outputFilename, 'width', legendTitle, h1LegendLabel, h2LegendLabel)

  #################################################################################################
  # Closure test
  #################################################################################################
  def plot_closure_test(self, i, jetR, sd_label, zcut, beta):
    
    output_dir_closure = os.path.join(self.output_dir, 'ClosureTest')
    if not os.path.isdir(output_dir_closure):
      os.makedirs(output_dir_closure)

    # Unfold smeared det-level spectrum with RM
    response = getattr(self, 'roounfold_response_R{}_{}'.format(jetR, sd_label))
    hMC_Det = getattr(self, 'hMC_Det_R{}_{}'.format(jetR, sd_label))
    hMC_Truth = getattr(self, 'hMC_Truth_R{}_{}'.format(jetR, sd_label))
    
    unfold2 = ROOT.RooUnfoldBayes(response, hMC_Det, i)
    hUnfolded2 = unfold2.Hreco() # Produces the truth distribution, with errors, PerBin d

    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(sd_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(sd_label))
    for bin in range(1, n_pt_bins_truth-3):
      min_pt_truth = truth_pt_bin_array[bin]
      max_pt_truth = truth_pt_bin_array[bin+1]
      
      self.plot_sd_closure_slice(hUnfolded2, hMC_Truth, i, jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth, output_dir_closure)

    # Closure test for pt dimension
    self.plot_pt_closure_slice(hUnfolded2, hMC_Truth, i, jetR, sd_label, zcut, beta, truth_pt_bin_array[1], truth_pt_bin_array[-4], output_dir_closure)

  #################################################################################################
  # Apply RM to unfolded result, and check that I obtain measured spectrum (simple technical check)
  #################################################################################################
  def plot_sd_closure_slice(self, hUnfolded, hMC_Truth, i, jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth, output_dir_closure):
    
    hUnfolded.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    hUnfolded_obs = hUnfolded.ProjectionY()
    hUnfolded_obs.SetName('hUnfolded_obs_R{}_{}_{}_{}-{}'.format(jetR, sd_label, i, min_pt_truth, max_pt_truth))
    
    hMC_Truth.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    hMCTruth_obs = hMC_Truth.ProjectionY()
    hMCTruth_obs.SetName('hMCTruth_obs_R{}_{}_{}_{}-{}'.format(jetR, sd_label, i, min_pt_truth, max_pt_truth))
    
    legendTitle = ''
    h1LegendLabel = 'Unfolded MC-det, {} = {}'.format(self.regularizationParamName,i)
    h2LegendLabel = 'MC-truth'
    ratioYAxisTitle = 'Unfolded MC det / Truth'
    outputFilename = os.path.join(output_dir_closure, 'hClosure_R{}_{}_{}-{}_{}{}'.format(self.utils.remove_periods(jetR), sd_label, int(min_pt_truth), int(max_pt_truth), i, self.file_format))
    self.plot_sd_obs_ratio(hUnfolded_obs, hMCTruth_obs, None, self.ytitle, ratioYAxisTitle, min_pt_truth, max_pt_truth, jetR, sd_label, zcut, beta, outputFilename, 'width', legendTitle, h1LegendLabel, h2LegendLabel)

  #################################################################################################
  # Apply RM to unfolded result, and check that I obtain measured spectrum (simple technical check)
  #################################################################################################
  def plot_pt_closure_slice(self, hUnfolded, hMC_Truth, i, jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth, output_dir_closure):
  
    hUnfolded.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    hUnfolded_pt = hUnfolded.ProjectionX()
    hUnfolded_pt.SetName('hUnfolded_pt_R{}_{}_{}'.format(jetR, sd_label, i))
    
    hMC_Truth.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    hMCTruth_pt = hMC_Truth.ProjectionX()
    hMCTruth_pt.SetName('hMCTruth_pt_R{}_{}_{}_'.format(jetR, sd_label, i))
    
    legendTitle = ''
    h1LegendLabel = 'Unfolded MC-det, {} = {}'.format(self.regularizationParamName,i)
    h2LegendLabel = 'MC-truth'
    ratioYAxisTitle = 'Unfolded MC det / Truth'
    outputFilename = os.path.join(output_dir_closure, 'hClosure_pt_R{}_{}_{}{}'.format(self.utils.remove_periods(jetR), sd_label, i, self.file_format))
    self.plot_sd_obs_ratio(hUnfolded_pt, hMCTruth_pt, None, self.ytitle, ratioYAxisTitle, 0, 0, jetR, sd_label, zcut, beta, outputFilename, 'width', legendTitle, h1LegendLabel, h2LegendLabel)

  #################################################################################################
  # Get errors from measured spectrum, stored as dictionary {bin:error}
  #################################################################################################
  def getMeasuredErrors(self, h):
    
    dictErrors = {}
    for binx in range(1, h.GetNbinsX()+1):
      for biny in range(1, h.GetNbinsY()+1):
        global_bin = h.GetBin(binx, biny)
        content = h.GetBinContent(global_bin)
        error = h.GetBinError(global_bin)
        if content!=0:
          dictErrors[global_bin] = error/content
        #print('Bin ({},{}) with content {} has error {}'.format(binx, biny, content, error/content))
        else:
          print('    Error: check your histogram ranges - something might be wrong')
          print('    Bin {}, content {}, error {}'.format(global_bin, content, error))
          dictErrors[global_bin] = 0
    return dictErrors

  #################################################################################################
  # Smear spectrum according to the error bars on the measured spectrum
  #################################################################################################
  def smearSpectrum(self, h, measuredErrors):
    
    # Zero errors in all bins
    for binx in range(0, h.GetNbinsX() + 1):
      for biny in range(0, h.GetNbinsY() + 1):
        global_bin = h.GetBin(binx, biny)
        h.SetBinError(global_bin, 0)
    
    # Loop through relevant bins and smear content according to errors in measured data, and set new errors
    r = ROOT.TRandom3(0)
    for bin, error in measuredErrors.items():
      content = h.GetBinContent(bin)
      errorNew = error * content
      contentNew = content + r.Gaus(0, errorNew)
      h.SetBinContent(bin, contentNew)
      h.SetBinError(bin, errorNew)
      #print('Setting bin {} to have relative error {}'.format(bin, errorNew/contentNew))

  #################################################################################################
  # Get MC-det 2D projection
  #################################################################################################
  def get_MCdet2D(self, jetR, sd_label):
    
    # (pt-det, pt-true, obs-det, obs-true)
    name_response = getattr(self, 'name_thn_rebinned_R{}_{}'.format(jetR, sd_label))
    hResponse = getattr(self, name_response)

    hResponse4D = hResponse.Clone()
    hResponse4D.SetName('{}_clone'.format(hResponse4D.GetName()))
    
    hMC_Det = hResponse4D.Projection(2,0)
    hMC_Det.SetName('hMC_Det_R{}_{}'.format(jetR, sd_label))
    return hMC_Det

  #################################################################################################
  # Get MC-det 2D projection
  #################################################################################################
  def get_MCtruth2D(self, jetR, sd_label):
    
    # (pt-det, pt-true, obs-det, obs-true)
    name_response = getattr(self, 'name_thn_rebinned_R{}_{}'.format(jetR, sd_label))
    hResponse = getattr(self, name_response)
    
    hResponse4D = hResponse.Clone()
    hResponse4D.SetName('{}_clone'.format(hResponse4D.GetName()))
    
    hMC_Truth = hResponse4D.Projection(3,1)
    hMC_Truth.SetName('hMC_Truth_R{}_{}'.format(jetR, sd_label))
    return hMC_Truth

  #################################################################################################
  # Plot spectra and ratio of h (and h3, if supplied) to h2
  #################################################################################################
  def plot_sd_obs_ratio(self, h, h2, h3, yAxisTitle, ratioYAxisTitle, min_pt_det, max_pt_det, jetR, sd_label, zcut, beta, outputFilename, scalingOptions = "", legendTitle = "",hLegendLabel = "", h2LegendLabel = "", h3LegendLabel = "", yRatioMax = 2.2):
    
    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()
    
    c = ROOT.TCanvas("c","c: pT",800,850)
    c.cd()
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)
    pad1.SetLeftMargin(0.2)
    pad1.SetRightMargin(0.05)
    pad1.SetTopMargin(0.05)
    if '_pt_' in outputFilename:
      pad1.SetLogy()
    pad1.Draw()
    pad1.cd()
    
    h.SetLineColor(1)
    h.SetLineWidth(2)
    h.SetLineStyle(1)
    
    h.Scale(1./h.Integral(), scalingOptions)
    if '_pt_' in outputFilename:
      h.GetYaxis().SetTitle('#frac{d#it{N}}{d#it{p}_{T}}')
    else:
      h.GetYaxis().SetTitle(yAxisTitle)
    
    h.GetYaxis().SetTitleSize(0.06)
    if self.observable == 'theta_g':
      if '_pt_' in outputFilename:
        h.GetYaxis().SetRangeUser(1e-4, 10.)
      else:
        h.GetYaxis().SetRangeUser(0., 3.)

    if self.observable == 'zg':
      if '_pt_' in outputFilename:
        h.GetYaxis().SetRangeUser(1e-4, 10.)
      else:
        h.GetYaxis().SetRangeUser(0., 8.)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(20)
    xAxisTitle = self.xtitle
    h.GetXaxis().SetTitle("")
    
    h2.SetLineColor(4)
    h2.SetLineWidth(2)
    h2.SetLineStyle(1)
    h2.Scale(1./h2.Integral(), scalingOptions)
    
    h.Draw("hist same E")
    h2.Draw("hist same E")
    
    if h3:
      h3.SetLineColor(2)
      h3.SetLineWidth(2)
      h3.SetLineStyle(1)
      h3.Scale(1./h3.Integral(), scalingOptions)
      h3.Draw("hist same")
    
    c.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.2)
    pad2.SetRightMargin(0.05)
    pad2.Draw()
    pad2.cd()

    # plot ratio h/h2
    hRatio = h.Clone()
    hRatio.Divide(h2)
    hRatio.SetMarkerStyle(21)
    hRatio.SetMarkerColor(1)
    
    hRatio.GetXaxis().SetTitleSize(30)
    hRatio.GetXaxis().SetTitleFont(43)
    hRatio.GetXaxis().SetTitleOffset(4.)
    hRatio.GetXaxis().SetLabelFont(43)
    hRatio.GetXaxis().SetLabelSize(20)
    if '_pt_' in outputFilename:
      hRatio.GetXaxis().SetTitle('#it{p}_{T,jet}')
    else:
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
    
    pad1.cd()
    
    leg2 = ROOT.TLegend(0.55,0.7,0.8,0.93,legendTitle)
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

    if not '_pt_' in outputFilename:
      text_latex = ROOT.TLatex()
      text_latex.SetNDC()
      text = str(min_pt_det) + ' < #it{p}_{T,ch jet} < ' + str(max_pt_det)
      text_latex.DrawLatex(0.25, 0.85, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'R = ' + str(jetR)
    text_latex.DrawLatex(0.25, 0.75, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'z_{cut} = ' + str(zcut) + '   #beta = ' + str(beta)
    text_latex.DrawLatex(0.25, 0.65, text)

    c.SaveAs(outputFilename)
    c.Close()

#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Unfold sd-observable distribution')
  parser.add_argument('-m', '--observable', action='store',
                      type=str, metavar='inputFileData',
                      default='theta_g',
                      help='Observable to be unfolded')
  parser.add_argument('-d', '--inputFileData', action='store',
                      type=str, metavar='inputFileData',
                      default='AnalysisResults.root',
                      help='Path of AnalysisResults.root file containing spectrum to be unfolded')
  parser.add_argument('-r', '--inputFileResponse', action='store',
                      type=str, metavar='inputFileResponse',
                      default='AnalysisResults.root',
                      help='Path of AnalysisResults.root file containing response matrix')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='analysis_config.yaml',
                      help='Path of config file for analysis')
  parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./unfolding_output/',
                      help='Output directory for plots to be written to')
  parser.add_argument('-i', '--imageFormat', action='store',
                      type=str, metavar='imageFormat',
                      default='.pdf',
                      help='Image format to save plots in, e.g. \".pdf\" or \".png\"')
                      
  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('inputFileData: \'{0}\''.format(args.inputFileData))
  print('inputFileResponse: \'{0}\''.format(args.inputFileResponse))
  print('configFile: \'{0}\''.format(args.configFile))
  print('ouputDir: \'{0}\''.format(args.outputDir))
  print('imageFormat: \'{0}\''.format(args.imageFormat))
  
  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFileData):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFileData))
    sys.exit(0)
  if not os.path.exists(args.inputFileResponse):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFileResponse))
    sys.exit(0)
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = roounfold_sd(observable=args.observable, input_file_data = args.inputFileData, input_file_response = args.inputFileResponse, config_file = args.configFile, output_dir = args.outputDir, file_format = args.imageFormat, rebin_response=True, truncation=False, binning=False, power_law_offset=0.)
  analysis.roounfold_sd()
