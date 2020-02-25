#!/usr/bin/env python3

"""
  Analysis utilities for jet analysis with track dataframe.
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import math

# Data analysis and plotting
import uproot
import pandas
import numpy as np
from array import *
import ROOT

# Base class
from pyjetty.alice_analysis.analysis.base import common_base

################################################################
class AnalysisUtils(common_base.CommonBase):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(AnalysisUtils, self).__init__(**kwargs)

  #---------------------------------------------------------------
  # Rebin 2D (pt, my_observable) histogram according to specified binnings
  #---------------------------------------------------------------
  def rebin_data(self, hData, name_data, n_pt_bins, pt_bin_array, n_obs_bins, obs_bin_array):
    
    # Create empty TH2 with appropriate binning
    name = "%s_rebinned" % name_data
    h = ROOT.TH2F(name, name, n_pt_bins, pt_bin_array, n_obs_bins, obs_bin_array)
    if h.GetSumw2() is 0:
      print('sumw2 not set')
    else:
      print('sumw2 set')
    
    # Loop over all bins (including under/over-flow -- needed e.g. for SD tagging rate), 
    # and fill rebinned histogram
    for bin_x in range(0, hData.GetNbinsX() + 2):
      for bin_y in range(0, hData.GetNbinsY() + 2):
        x = hData.GetXaxis().GetBinCenter(bin_x)
        y = hData.GetYaxis().GetBinCenter(bin_y)
        
        content = hData.GetBinContent(bin_x, bin_y)
        h.Fill(x, y, content)
          
    for bin in range(0, h.GetNcells()+1):
      content = h.GetBinContent(bin)
      old_uncertainty = h.GetBinError(bin)
      new_uncertainty = math.sqrt(content)
      h.SetBinError(bin, new_uncertainty)
    
    if h.GetSumw2() is 0:
      h.Sumw2()
    return h

  #---------------------------------------------------------------
  # Rebin 4D THn response according to specified binnings, saving both a rebinned
  # THn and a RooUnfoldResponse object, and write to file
  #---------------------------------------------------------------
  def rebin_response(self, response_file_name, thn, name_thn_rebinned, name_roounfold, label,
                     n_pt_bins_det, det_pt_bin_array, n_obs_bins_det, det_obs_bin_array,
                     n_pt_bins_truth, truth_pt_bin_array, n_obs_bins_truth, truth_obs_bin_array,
                     observable, power_law_offset=0.):
  
    # Create empty THn with specified binnings
    thn_rebinned = self.create_empty_thn(name_thn_rebinned, n_pt_bins_det, det_pt_bin_array, 
                                         n_obs_bins_det, det_obs_bin_array, n_pt_bins_truth,
                                         truth_pt_bin_array, n_obs_bins_truth, truth_obs_bin_array)
    
    # Create empty RooUnfoldResponse with specified binning
    hist_measured = thn_rebinned.Projection(2, 0)
    hist_measured.SetName('hist_measured_%s' % label)
    hist_truth = thn_rebinned.Projection(3, 1)
    hist_truth.SetName('hist_truth_%s' % label)
    roounfold_response = ROOT.RooUnfoldResponse(hist_measured, hist_truth, name_roounfold, 
                                                name_roounfold) # Sets up binning
    
    # Note: Using overflow bins doesn't work for 2D unfolding in RooUnfold
    #roounfold_response.UseOverflow(True)
    
    # Loop through THn and fill rebinned THn
    self.fill_new_response(response_file_name, thn, thn_rebinned, roounfold_response, 
                           observable, power_law_offset)
  
  #---------------------------------------------------------------
  # Loop through original THn, and fill new response (THn and RooUnfoldResponse)
  #---------------------------------------------------------------
  def fill_new_response(self, response_file_name, thn, thn_rebinned, roounfold_response, 
                        observable, power_law_offset=0.):
    
    # I don't find any global bin index implementation, so I manually loop through axes 
    # (including under/over-flow)...
    for bin_0 in range(0, thn.GetAxis(0).GetNbins() + 2):
      if bin_0 % 5 == 0:
        print('{} / {}'.format(bin_0, thn.GetAxis(0).GetNbins() + 2))
      pt_det = thn.GetAxis(0).GetBinCenter(bin_0)
      for bin_1 in range(0, thn.GetAxis(1).GetNbins() + 2):
        pt_true = thn.GetAxis(1).GetBinCenter(bin_1)
        for bin_2 in range(0, thn.GetAxis(2).GetNbins() + 2):
          obs_det = thn.GetAxis(2).GetBinCenter(bin_2)
          for bin_3 in range(0, thn.GetAxis(3).GetNbins() + 2):
            obs_true = thn.GetAxis(3).GetBinCenter(bin_3)
            
            x_list = (pt_det, pt_true, obs_det, obs_true)
            x = array('d', x_list)
            global_bin = thn.GetBin(x)
            content = thn.GetBinContent(global_bin)
            
            # Impose a custom prior, if desired
            if math.fabs(power_law_offset) > 1e-3 :
              #print('Scaling prior by power_law_offset={}'.format(power_law_offset))
              if pt_true > 0. and obs_true > 0.:

                scale_factor = math.pow(pt_true, power_law_offset)
                
                if observable == 'zg':
                  scale_factor *= math.pow(obs_true, power_law_offset)
                elif observable == 'theta_g':
                  scale_factor *= (1 + obs_true)

                content = content*scale_factor
          
            # THn is filled as (pt_det, pt_true, obs_det, obs_true)
            thn_rebinned.Fill(x, content)
            #print('Fill ({}, {}, {}, {}) to response'.format(pt_det, pt_true, obs_det, obs_true))
            
            # RooUnfoldResponse should be filled (pt_det, obs_det, pt_true, obs_true)
            roounfold_response.Fill(pt_det, obs_det, pt_true, obs_true, content)

    print('writing response...')
    f = ROOT.TFile(response_file_name, 'UPDATE')
    thn_rebinned.Write()
    roounfold_response.Write()
    f.Close()
    print('done')

  #---------------------------------------------------------------
  # Create an empty THn according to specified binnings
  #---------------------------------------------------------------
  def create_empty_thn(self, name, n_pt_bins_det, det_pt_bin_array, n_obs_bins_det, det_obs_bin_array, 
                       n_pt_bins_truth, truth_pt_bin_array, n_obs_bins_truth, truth_obs_bin_array):
    
    # Create THn of response
    dim = 0;
    title = []
    nbins = []
    min = []
    max = []
    bin_edges = []
      
    title.append('p_{T,det}')
    nbins.append(n_pt_bins_det)
    bin_edges.append(det_pt_bin_array)
    min.append(det_pt_bin_array[0])
    max.append(det_pt_bin_array[-1])
    dim+=1
          
    title.append('p_{T,truth}')
    nbins.append(n_pt_bins_truth)
    bin_edges.append(truth_pt_bin_array)
    min.append(truth_pt_bin_array[0])
    max.append(truth_pt_bin_array[-1])
    dim+=1
      
    title.append('#obs_{det}')
    nbins.append(n_obs_bins_det)
    bin_edges.append(det_obs_bin_array)
    min.append(det_obs_bin_array[0])
    max.append(det_obs_bin_array[-1])
    dim+=1

    title.append('#obs_{truth}')
    nbins.append(n_obs_bins_truth)
    bin_edges.append(truth_obs_bin_array)
    min.append(truth_obs_bin_array[0])
    max.append(truth_obs_bin_array[-1])
    dim+=1
      
    nbins = (nbins)
    xmin = (min)
    xmax = (max)
    nbins_array = array('i', nbins)
    xmin_array = array('d', xmin)
    xmax_array = array('d', xmax)
    h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
    for i in range(0, dim):
      h.GetAxis(i).SetTitle(title[i])
      h.SetBinEdges(i, bin_edges[i])

    return h

  #----------------------------------------------------------------------
  # Check if re-binned response file exists -- create it if not (or if forced)
  #----------------------------------------------------------------------
  def check_rebin_response(self, output_dir, force_rebin_response):
    
    rebin_response = True
    response_path = os.path.join(output_dir, 'response.root')
    if os.path.exists(response_path):
      rebin_response = False
      if force_rebin_response:
        print('Response {} exists -- force re-create...'.format(response_path))
        rebin_response = True
      else:
        print('Response {} exists -- don\'t re-create.'.format(response_path))
    else:
      print('Response {} doesn\'t exist -- create it...'.format(response_path))
    
    return rebin_response

  #----------------------------------------------------------------------
  # Add two histograms in quadrature
  #----------------------------------------------------------------------
  def add_in_quadrature(self, h1, h2):

    h_new = h1.Clone()
    h_new.SetName('{}_new'.format(h1.GetName()))

    for i in range(1, h_new.GetNbinsX()+1):
      value1 = h1.GetBinContent(i)
      value2 = h2.GetBinContent(i)
      
      new_value = math.sqrt(value1*value1 + value2*value2)
      
      h_new.SetBinContent(i, new_value)
    
    return h_new
  
  #---------------------------------------------------------------
  # Normalize a histogram by its integral
  #---------------------------------------------------------------
  def scale_by_integral(self, h):
    
    if h.GetSumw2() is 0:
      h.Sumw2()
    
    integral = h.Integral()
    if integral > 0:
      h.Scale(1./integral)
    else:
      print('Integral is 0, check for problem')

  #---------------------------------------------------------------
  # Multiply a tgraph by a histogram, point-by-point
  # Assumes the histogram has zero uncertainty.
  #---------------------------------------------------------------
  def multiply_tgraph(self, g, h):
  
    nBins = h.GetNbinsX()
    for bin in range(1, nBins+1):
      
      # Get histogram (x,y)
      h_x = h.GetBinCenter(bin)
      h_y = h.GetBinContent(bin)
      
      # Get TGraph (x,y) and errors
      g_x = ROOT.Double(0)
      g_y = ROOT.Double(0)
      g.GetPoint(bin-1, g_x, g_y)
      yErrLow = g.GetErrorYlow(bin-1)
      yErrUp  = g.GetErrorYhigh(bin-1)
      
      #print('hist x: {}'.format(h_x))
      #print('graph x: {}'.format(g_x))
      
      new_content = g_y * h_y
      new_error_low = yErrLow * h_y
      new_error_up = yErrUp * h_y
      
      g.SetPoint(bin-1, h_x, new_content)
      g.SetPointError(bin-1, 0, 0, new_error_low, new_error_up)

  #---------------------------------------------------------------
  # Divide a histogram by a tgraph, point-by-point
  #---------------------------------------------------------------
  def divide_tgraph(self, h, g, combine_errors=False):
    
    nBins = h.GetNbinsX()
    for bin in range(1, nBins+1):

      # Get histogram (x,y)
      h_x = h.GetBinCenter(bin)
      h_y = h.GetBinContent(bin)
      h_error = h.GetBinError(bin)
      
      # Get TGraph (x,y) and errors
      g_x = ROOT.Double(0)
      g_y = ROOT.Double(0)
      g.GetPoint(bin-1, g_x, g_y)
      yErrLow = g.GetErrorYlow(bin-1)
      yErrUp  = g.GetErrorYhigh(bin-1)
      
      #print('hist x: {}'.format(h_x))
      #print('graph x: {}'.format(g_x))
      
      new_content = h_y / g_y
      if combine_errors: # combine tgraph and histogram relative uncertainties in quadrature
        new_error_low = math.sqrt( pow(yErrLow/g_y,2) + pow(h_error/h_y,2) ) * new_content
        new_error_up = math.sqrt( pow(yErrUp/g_y,2) + pow(h_error/h_y,2) ) * new_content
      else: # Assume tgraph has zero uncertainty
        new_error_low = h_error / g_y
        new_error_up = h_error / g_y
  
      if combine_errors:
        g.SetPoint(bin-1, h_x, new_content)
        g.SetPointError(bin-1, 0, 0, new_error_low, new_error_up)
      else:
        h.SetBinContent(bin, new_content)
        h.SetBinError(bin, new_error_up)
  
  #---------------------------------------------------------------
  # Remove periods from a label
  #---------------------------------------------------------------
  def remove_periods(self, text):
  
    string = str(text)
    return string.replace('.', '')
  
  #---------------------------------------------------------------
  # Set legend parameters
  #---------------------------------------------------------------
  def setup_legend(self, leg, textSize):
    
    leg.SetTextFont(42);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetFillColor(0);
    leg.SetMargin(0.25);
    leg.SetTextSize(textSize);
    leg.SetEntrySeparation(0.5);
  
  #---------------------------------------------------------------
  # Configure style options
  #---------------------------------------------------------------
  def set_plotting_options(self):
    
    font = 42
    
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetFrameFillColor(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetPadColor(10)
    ROOT.gStyle.SetCanvasColor(10)
    ROOT.gStyle.SetTitleFillColor(10)
    ROOT.gStyle.SetTitleBorderSize(1)
    ROOT.gStyle.SetStatColor(10)
    ROOT.gStyle.SetStatBorderSize(1)
    ROOT.gStyle.SetLegendBorderSize(1)
    
    ROOT.gStyle.SetDrawBorder(0)
    ROOT.gStyle.SetTextFont(font)
    ROOT.gStyle.SetStatFont(font)
    ROOT.gStyle.SetStatFontSize(0.05)
    ROOT.gStyle.SetStatX(0.97)
    ROOT.gStyle.SetStatY(0.98)
    ROOT.gStyle.SetStatH(0.03)
    ROOT.gStyle.SetStatW(0.3)
    ROOT.gStyle.SetTickLength(0.02,"y")
    ROOT.gStyle.SetEndErrorSize(3)
    ROOT.gStyle.SetLabelSize(0.05,"xyz")
    ROOT.gStyle.SetLabelFont(font,"xyz")
    ROOT.gStyle.SetLabelOffset(0.01,"xyz")
    ROOT.gStyle.SetTitleFont(font,"xyz")
    ROOT.gStyle.SetTitleOffset(1.2,"xyz")
    ROOT.gStyle.SetTitleSize(0.045,"xyz")
    ROOT.gStyle.SetMarkerSize(1)
    ROOT.gStyle.SetPalette(1)
    
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)

  #---------------------------------------------------------------
  # Plot and save a 1D histogram
  #---------------------------------------------------------------
  def plot_hist(self, h, outputFilename, drawOptions = '', setLogy = False, setLogz = False, text = ''):
    
    self.set_plotting_options()
    ROOT.gROOT.ForceStyle()
    
    h.SetLineColor(1)
    h.SetLineWidth(1)
    h.SetLineStyle(1)
    
    c = ROOT.TCanvas("c","c: hist",600,450)
    c.cd()
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetRightMargin(0.15)
    if setLogy:
      c.SetLogy()
    if setLogz:
      c.SetLogz()

    h.Draw(drawOptions)
    
    if text:
      textFit = ROOT.TLatex()
      textFit.SetTextSize(0.04)
      textFit.SetNDC()
      textFit.DrawLatex(0.3,0.8,text)
    
    c.SaveAs(outputFilename)
    c.Close()
