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

from pyjetty.mputils import treereader

################################################################
class analysis_utils(common_base.common_base):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(analysis_utils, self).__init__(**kwargs)

  #---------------------------------------------------------------
  # Rebin 2D (pt, theta) histogram according to specified binnings
  #---------------------------------------------------------------
  def rebin_data(self, hData, name_data, n_pt_bins, pt_bin_array, n_rg_bins, rg_bin_array):
    
    # Create empty TH2 with appropriate binning
    name = '{}_{}'.format(name_data, 'rebinned')
    h = ROOT.TH2F(name, name, n_pt_bins, pt_bin_array, n_rg_bins, rg_bin_array)
    h.Sumw2()
    
    # Loop over all bins (including under/over-flow -- needed for SD tagging rate), and fill rebinned histogram
    for bin_x in range(0, hData.GetNbinsX() + 2):
      for bin_y in range(0, hData.GetNbinsY() + 2):
        x = hData.GetXaxis().GetBinCenter(bin_x)
        y = hData.GetYaxis().GetBinCenter(bin_y)
        
        content = hData.GetBinContent(bin_x, bin_y)
        h.Fill(x, y, content)

    return h

  #---------------------------------------------------------------
  # Rebin THn according to specified binnings, saving both a rebinned
  # THn and a RooUnfoldResponse object, and write to file
  #---------------------------------------------------------------
  def rebin_response(self, response_file_name, thn, name_thn_rebinned, name_roounfold, jetR, sd_label, n_pt_bins_det, det_pt_bin_array, n_rg_bins_det, det_rg_bin_array, n_pt_bins_truth, truth_pt_bin_array, n_rg_bins_truth, truth_rg_bin_array, power_law_offset=0.):
  
    # Create empty THn with specified binnings
    thn_rebinned = self.create_empty_thn(name_thn_rebinned, n_pt_bins_det, det_pt_bin_array, n_rg_bins_det, det_rg_bin_array, n_pt_bins_truth, truth_pt_bin_array, n_rg_bins_truth, truth_rg_bin_array)
    
    # Create empty RooUnfoldResponse with specified binning
    hist_measured = thn_rebinned.Projection(2, 0)
    hist_measured.SetName('hist_measured_R{}_{}'.format(jetR, sd_label))
    hist_truth = thn_rebinned.Projection(3, 1)
    hist_truth.SetName('hist_truth_R{}_{}'.format(jetR, sd_label))
    roounfold_response = ROOT.RooUnfoldResponse(hist_measured, hist_truth, name_roounfold, name_roounfold) # Sets up binning
    
    # Note: Using overflow bins doesn't work for 2D unfolding in RooUnfold
    #roounfold_response.UseOverflow(True)
    
    # Loop through THn and fill rebinned THn
    self.fill_new_response(response_file_name, thn, thn_rebinned, roounfold_response, power_law_offset)
  
  #---------------------------------------------------------------
  # Loop through original THn, and fill new response (THn and RooUnfoldResponse)
  #---------------------------------------------------------------
  def fill_new_response(self, response_file_name, thn, thn_rebinned, roounfold_response, power_law_offset=0.):
    
    # I don't find any global bin index implementation, so I manually loop through axes (including under/over-flow)...
    for bin_0 in range(0, thn.GetAxis(0).GetNbins() + 2):
      if bin_0 % 5 == 0:
        print('{} / {}'.format(bin_0, thn.GetAxis(0).GetNbins() + 2))
      pt_det = thn.GetAxis(0).GetBinCenter(bin_0)
      for bin_1 in range(0, thn.GetAxis(1).GetNbins() + 2):
        pt_true = thn.GetAxis(1).GetBinCenter(bin_1)
        for bin_2 in range(0, thn.GetAxis(2).GetNbins() + 2):
          theta_det = thn.GetAxis(2).GetBinCenter(bin_2)
          for bin_3 in range(0, thn.GetAxis(3).GetNbins() + 2):
            theta_true = thn.GetAxis(3).GetBinCenter(bin_3)
            
            x_list = (pt_det, pt_true, theta_det, theta_true)
            x = array('d', x_list)
            global_bin = thn.GetBin(x)
            content = thn.GetBinContent(global_bin)
            
            # Impose a custom prior, if desired
            if math.fabs(power_law_offset) > 1e-3 :
              #print('Scaling prior by power_law_offset={}'.format(power_law_offset))
              if pt_true > 0.:
                scale_factor = math.pow(pt_true, power_law_offset)
                content = content*scale_factor
          
            # THn is filled as (pt_det, pt_true, theta_det, theta_true)
            thn_rebinned.Fill(x, content)
            #print('Fill ({}, {}, {}, {}) to response'.format(pt_det, pt_true, theta_det, theta_true))
            
            # RooUnfoldResponse should be filled (pt_det, theta_det, pt_true, theta_true)
            roounfold_response.Fill(pt_det, theta_det, pt_true, theta_true, content)

    print('writing response...')
    f = ROOT.TFile(response_file_name, 'UPDATE')
    thn_rebinned.Write()
    roounfold_response.Write()
    f.Close()
    print('done')

  #---------------------------------------------------------------
  # Construct 2D (pt, theta) histogram according to specified binnings
  # Note: automatically fills under/over-flow bins (needed for SD tagging rate)
  #---------------------------------------------------------------
  def construct_data_histograms(self, tree_file_name, name_data, n_pt_bins, pt_bin_array, n_rg_bins, rg_bin_array):
    
    # Create empty TH2 with appropriate binning
    name = '{}_{}'.format(name_data, 'rebinned')
    h = ROOT.TH2F(name, name, n_pt_bins, pt_bin_array, n_rg_bins, rg_bin_array)
    h.Sumw2()
    
    # Loop through tree and fill each entry into histogram
    tr = treereader.RTreeReader(tree_name='t',
                                branches = ['j_pt', 'sd_j_dR'],
                                file_name=tree_file_name)
      
    for i in range(tr.tree.GetEntries()):
      tr.tree.GetEntry(i)
      if tr.j_pt.size() > 0:
        
        pt = tr.j_pt[0]
        theta = tr.sd_j_dR[0]
        h.Fill(pt, theta)
  
    return h

  #---------------------------------------------------------------
  # Construct THn and RooUnfoldResponse object from tree,
  # according to specified binnings, and write to file
  #---------------------------------------------------------------
  def construct_response_histograms(self, tree_file_name, response_file_name, name_thn_rebinned, name_roounfold, jetR, sd_label, n_pt_bins_det, det_pt_bin_array, n_rg_bins_det, det_rg_bin_array, n_pt_bins_truth, truth_pt_bin_array, n_rg_bins_truth, truth_rg_bin_array, power_law_offset=0.):
    
    # Create empty THn with specified binnings
    thn_rebinned = self.create_empty_thn(name_thn_rebinned, n_pt_bins_det, det_pt_bin_array, n_rg_bins_det, det_rg_bin_array, n_pt_bins_truth, truth_pt_bin_array, n_rg_bins_truth, truth_rg_bin_array)
    
    # Create empty RooUnfoldResponse with specified binning
    hist_measured = thn_rebinned.Projection(2, 0)
    hist_measured.SetName('hist_measured_R{}_{}'.format(jetR, sd_label))
    hist_truth = thn_rebinned.Projection(3, 1)
    hist_truth.SetName('hist_truth_R{}_{}'.format(jetR, sd_label))
    roounfold_response = ROOT.RooUnfoldResponse(hist_measured, hist_truth, name_roounfold, name_roounfold) # Sets up binning
    
    # Note: Using overflow bins doesn't work for 2D unfolding in RooUnfold
    #roounfold_response.UseOverflow(True)
    
    # Loop through tree and fill response objects
    self.fill_response_histograms(tree_file_name, response_file_name, thn_rebinned, roounfold_response, power_law_offset)
  
  #---------------------------------------------------------------
  # Loop through original THn, and fill new response (THn and RooUnfoldResponse)
  #---------------------------------------------------------------
  def fill_response_histograms(self, tree_file_name, response_file_name, thn_rebinned, roounfold_response, power_law_offset=0.):
    
    tr = treereader.RTreeReader(tree_name='t',
                                branches = ['j_pt', 'ej_pt', 'sd_j_dR', 'sd_ej_dR'],
                                file_name=tree_file_name)
      
    for i in range(tr.tree.GetEntries()):
      tr.tree.GetEntry(i)
      if tr.j_pt.size() > 0 and tr.ej_pt.size() > 0:
        
        pt_det = tr.ej_pt[0]
        pt_true = tr.j_pt[0]
        theta_det = tr.sd_ej_dR[0]
        theta_true = tr.sd_j_dR[0]
      
        # Impose a custom prior, if desired
        content = 1
        if math.fabs(power_law_offset) > 1e-3 :
          #print('Scaling prior by power_law_offset={}'.format(power_law_offset))
          if pt_true > 0.:
            scale_factor = math.pow(pt_true, power_law_offset)
            content = content*scale_factor
              
        # THn is filled as (pt_det, pt_true, theta_det, theta_true)
        x_list = (pt_det, pt_true, theta_det, theta_true)
        x = array('d', x_list)
        thn_rebinned.Fill(x)
        #print('Fill ({}, {}, {}, {}) to response'.format(pt_det, pt_true, theta_det, theta_true))
        
        # RooUnfoldResponse should be filled (pt_det, theta_det, pt_true, theta_true)
        roounfold_response.Fill(pt_det, theta_det, pt_true, theta_true)

      print('writing response...')
      f = ROOT.TFile(response_file_name, 'UPDATE')
      thn_rebinned.Write()
      roounfold_response.Write()
      f.Close()
      print('done')

  #---------------------------------------------------------------
  # Create an empty THn according to specified binnings
  #---------------------------------------------------------------
  def create_empty_thn(self, name, n_pt_bins_det, det_pt_bin_array, n_rg_bins_det, det_rg_bin_array, n_pt_bins_truth, truth_pt_bin_array, n_rg_bins_truth, truth_rg_bin_array):
    
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
      
    title.append('#theta_{g,det}')
    nbins.append(n_rg_bins_det)
    bin_edges.append(det_rg_bin_array)
    min.append(det_rg_bin_array[0])
    max.append(det_rg_bin_array[-1])
    dim+=1

    title.append('#theta_{g,truth}')
    nbins.append(n_rg_bins_truth)
    bin_edges.append(truth_rg_bin_array)
    min.append(truth_rg_bin_array[0])
    max.append(truth_rg_bin_array[-1])
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
  # Get regularization parameter
  #---------------------------------------------------------------
  def get_reg_param(self, sd_settings, sd_config_list, sd_config_dict, sd_label, observable, jetR):
    
    for i, sd_setting in enumerate(sd_settings):
      
      zcut = sd_setting[0]
      beta = sd_setting[1]
      label = 'zcut{}_B{}'.format(self.remove_periods(zcut), beta)
      if label == sd_label:
        
        config_name = sd_config_list[i]
        
        reg_param = sd_config_dict[config_name]['reg_param'][observable][jetR]
        #print('reg_param for {} {} jetR={}: {}'.format(sd_label, observable, jetR, reg_param))
          
        return reg_param
      
      else:
        continue

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
