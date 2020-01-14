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
from pyjetty.alice_analysis.analysis.base import analysis_utils


from pyjetty.mputils import treereader

################################################################
class tree_utils(analysis_utils.analysis_utils):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(tree_utils, self).__init__(**kwargs)

  #---------------------------------------------------------------
  # Construct 2D (pt, obs) histogram according to specified binnings
  # Note: automatically fills under/over-flow bins (needed e.g. for SD tagging rate)
  #---------------------------------------------------------------
  def construct_data_histograms(self, tree_file_name, name_data, n_pt_bins, pt_bin_array, 
                                n_obs_bins, obs_bin_array):
    
    # Create empty TH2 with appropriate binning
    name = '%s_rebinned' % name_data
    h = ROOT.TH2F(name, name, n_pt_bins, pt_bin_array, n_obs_bins, obs_bin_array)
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
  def construct_response_histograms(self, tree_file_name, response_file_name, name_thn_rebinned, 
                                    name_roounfold, label, n_pt_bins_det, det_pt_bin_array,
                                    n_obs_bins_det, det_obs_bin_array, n_pt_bins_truth,
                                    truth_pt_bin_array, n_obs_bins_truth, truth_obs_bin_array,
                                    power_law_offset=0.):
    
    # Create empty THn with specified binnings
    thn_rebinned = self.create_empty_thn(name_thn_rebinned, n_pt_bins_det, det_pt_bin_array, n_obs_bins_det, det_obs_bin_array, n_pt_bins_truth, truth_pt_bin_array, n_obs_bins_truth, truth_rg_bin_array)
    
    # Create empty RooUnfoldResponse with specified binning
    hist_measured = thn_rebinned.Projection(2, 0)
    hist_measured.SetName('hist_measured_{}'.format(label))
    hist_truth = thn_rebinned.Projection(3, 1)
    hist_truth.SetName('hist_truth_{}'.format(label))
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
