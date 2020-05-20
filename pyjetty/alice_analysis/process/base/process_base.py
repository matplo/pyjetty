#!/usr/bin/env python3

"""
  Analysis task base class.
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import time

# Data analysis and plotting
import ROOT
import yaml
import numpy as np
from array import *

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib

# Analysis utilities
from pyjetty.alice_analysis.process.base import common_base
from pyjetty.alice_analysis.process.base import process_utils
from pyjetty.alice_analysis.process.base import jet_info

################################################################
class ProcessBase(common_base.CommonBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
    super(ProcessBase, self).__init__(**kwargs)
    self.input_file = input_file
    self.config_file = config_file
    self.output_dir = output_dir
    self.debug_level = debug_level # (0 = no debug info, 1 = some debug info, 2 = all debug info)
    
    # Create output dir
    if not self.output_dir.endswith("/"):
      self.output_dir = self.output_dir + "/"
    if not os.path.exists(self.output_dir):
      #os.makedirs(self.output_dir, 775)
      os.makedirs(self.output_dir)
      
    # Create output file
    outputfilename = os.path.join(self.output_dir, 'AnalysisResults.root')
    fout = ROOT.TFile(outputfilename, 'recreate')
    fout.Close()

    # Initialize utils class
    self.utils = process_utils.ProcessUtils()
  
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
  
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
      
    if 'event_number_max' in config:
      self.event_number_max = config['event_number_max']
    else:
      self.event_number_max = sys.maxsize
      
    self.jetR_list = config['jetR']
    self.debug_level = config['debug_level']

    # Check if constituent subtractor is included, and initialize it if so
    self.do_constituent_subtraction = False
    if 'constituent_subtractor' in config:
      print('Constituent subtractor is enabled.')
      self.do_constituent_subtraction = True
      constituent_subtractor = config['constituent_subtractor']
      
      self.max_distance = constituent_subtractor['max_distance']
      self.alpha = constituent_subtractor['alpha']
      self.max_eta = constituent_subtractor['max_eta']
      self.bge_rho_grid_size = constituent_subtractor['bge_rho_grid_size']
      self.max_pt_correct = constituent_subtractor['max_pt_correct']
      self.ghost_area = constituent_subtractor['ghost_area']
    else:
      print('Constituent subtractor is disabled.')

  #---------------------------------------------------------------
  # Create thn and set as class attribute from name, dim
  #   and lists of nbins, xmin, xmax.
  #---------------------------------------------------------------
  def create_thn(self, name, title, dim, nbins, xmin, xmax):
    
    nbins_arr = (nbins)
    xmin_arr = (min)
    xmax_arr = (max)
    nbins_array = array('i', nbins)
    xmin_array = array('d', xmin)
    xmax_array = array('d', xmax)
    h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
    for i in range(0, dim):
      h.GetAxis(i).SetTitle(title[i])
    setattr(self, name, h)

  #---------------------------------------------------------------
  # Return Lund coordinates [log(1/deltaR), log(1/kt)] of a SD jet
  #---------------------------------------------------------------
  def lund_coordinates_SD(self, jet_sd):

    sd_info = fjcontrib.get_SD_jet_info(jet_sd)
    dR = sd_info.dR
    z = sd_info.z
    pt = jet_sd.pt()
    kt = z*dR*pt
    
    if dR > 1e-5:
      return [np.log(1/dR), np.log(kt)]
    else:
      return [sys.maxsize, -sys.maxsize]

  #---------------------------------------------------------------
  # Return Lund coordinates [log(1/deltaR), log(1/kt)] of a DG jet
  #---------------------------------------------------------------
  def lund_coordinates_DG(self, jet_dg):

    dR = jet_dg.Delta()
    kt = jet_dg.kt()
    
    if dR > 1e-5:
      return [np.log(1/dR), np.log(kt)]
    else:
      return [sys.maxsize, -sys.maxsize]
    
  #---------------------------------------------------------------
  # Compare two jets and store matching candidates in user_info
  #---------------------------------------------------------------
  def set_matching_candidates(self, jet1, jet2, jetR, hname, fill_jet1_matches_only=False):
  
    # Fill histogram of matching distance of all candidates
    deltaR = jet1.delta_R(jet2)
    if hname:
      getattr(self, hname.format(jetR)).Fill(jet1.pt(), deltaR)
  
    # Add a matching candidate to the list if it is within the geometrical cut
    if deltaR < self.jet_matching_distance*jetR:
      self.set_jet_info(jet1, jet2, deltaR)
      if not fill_jet1_matches_only:
        self.set_jet_info(jet2, jet1, deltaR)

  #---------------------------------------------------------------
  # Set 'jet_match' as a matching candidate in user_info of 'jet'
  #---------------------------------------------------------------
  def set_jet_info(self, jet, jet_match, deltaR):
  
    # Get/create object to store list of matching candidates
    jet_user_info = None
    if jet.has_user_info():
      jet_user_info = jet.python_info()
    else:
      jet_user_info = jet_info.JetInfo()
      
    jet_user_info.matching_candidates.append(jet_match)
    if deltaR < jet_user_info.closest_jet_deltaR:
      jet_user_info.closest_jet = jet_match
      jet_user_info.closest_jet_deltaR = deltaR
          
    jet.set_python_info(jet_user_info)

  #---------------------------------------------------------------
  # Set accepted jet matches for pp case
  #
  # hname is the name of a matching QA histogram that will be created
  # it can be anything you want, e.g. 'hJetMatchingQA_R{}'.format(jetR)
  #---------------------------------------------------------------
  def set_matches_pp(self, jet_det, hname):

    # Create pp matching QA histogram, if not already created
    if not hasattr(self, hname):
      bin_labels = ['all', 'has_matching_candidate', 'unique_match']
      nbins = len(bin_labels)
      h = ROOT.TH2F(hname, hname, nbins, 0, nbins, 30, 0., 300.)
      for i in range(1, nbins+1):
        h.GetXaxis().SetBinLabel(i,bin_labels[i-1])
      setattr(self, hname, h)
      
    h = getattr(self, hname)
    h.Fill('all', jet_det.pt(), 1)
  
    if jet_det.has_user_info():
        
      jet_info_det = jet_det.python_info()
      h.Fill('has_matching_candidate', jet_det.pt(), 1)
      
      if len(jet_info_det.matching_candidates) == 1:
        jet_truth = jet_info_det.closest_jet
          
        # Check that the match is unique
        if jet_truth.has_user_info():
          jet_info_truth = jet_truth.python_info()
          if len(jet_info_truth.matching_candidates) == 1:
                  
            # Set accepted match
            jet_info_det.match = jet_truth
            jet_det.set_python_info(jet_info_det)
            h.Fill('unique_match', jet_det.pt(), 1)

  #---------------------------------------------------------------
  # Set accepted jet matches for Pb-Pb case
  #
  # hname is the name of a matching QA histogram that will be created
  # it can be anything you want, e.g. 'hJetMatchingQA_R{}'.format(jetR)
  #---------------------------------------------------------------
  def set_matches_AA(self, jet_det_combined, jetR, hname):
  
    set_match = False
    
    # Create Pb-Pb matching QA histogram, if not already created
    if not hasattr(self, hname):
      bin_labels = ['all', 'has_matching_ppdet_candidate', 'has_matching_pptrue_candidate', 'has_matching_pptrue_unique_candidate', 'mc_fraction', 'deltaR_combined-truth', 'passed_all_cuts']
      nbins = len(bin_labels)
      h = ROOT.TH2F(hname, hname, nbins, 0, nbins,  30, 0., 300.)
      for i in range(1, nbins+1):
        h.GetXaxis().SetBinLabel(i,bin_labels[i-1])
      setattr(self, hname, h)
    
    h = getattr(self, hname)
    h.Fill('all', jet_det_combined.pt(), 1)

    # Check if combined jet has a pp-det match (uniqueness not currently enforced)
    if jet_det_combined.has_user_info():
    
      h.Fill('has_matching_ppdet_candidate', jet_det_combined.pt(), 1)

      jet_info_combined = jet_det_combined.python_info()
      jet_pp_det = jet_info_combined.closest_jet
        
      # Check if the pp-det jet has a pp-truth match
      if jet_pp_det.has_user_info():
      
        jet_info_pp_det = jet_pp_det.python_info()
        jet_pp_truth = jet_info_pp_det.closest_jet
        h.Fill('has_matching_pptrue_candidate', jet_det_combined.pt(), 1)
        set_match = True
        
        # Check that pp-det to pp-true match is unique
        if self.is_match_unique(jet_pp_det):
          h.Fill('has_matching_pptrue_unique_candidate', jet_det_combined.pt(), 1)
        else:
          set_match = False
            
        # Check matching distance between combined jet and pp-truth
        if jet_det_combined.delta_R(jet_pp_truth) < self.jet_matching_distance*jetR:
          h.Fill('deltaR_combined-truth', jet_det_combined.pt(), 1)
        else:
          set_match = False
  
        # Check if >50% of the pp-det tracks are in the Pb-Pb det jet
        mc_fraction = self.mc_fraction(jet_pp_det, jet_det_combined)
        if mc_fraction > self.mc_fraction_threshold:
          h.Fill('mc_fraction', jet_det_combined.pt(), 1)
        else:
          set_match = False
               
    # Set accepted match
    if set_match:
      jet_info_combined.match = jet_pp_truth
      jet_det_combined.set_python_info(jet_info_combined)
      
      # Set also the pp-truth match info, since we will need to access the matching pp-det jet for prong matching
      jet_info_pp_truth = jet_pp_truth.python_info()
      jet_info_pp_truth.match = jet_pp_det
      jet_pp_truth.set_python_info(jet_info_pp_truth)
      
      h.Fill('passed_all_cuts', jet_det_combined.pt(), 1)

  #---------------------------------------------------------------
  # Set accepted jet matches for truth+background (i.e. no det-level)
  #
  # hname is the name of a matching QA histogram that will be created
  # it can be anything you want, e.g. 'hJetMatchingQA_R{}'.format(jetR)
  #---------------------------------------------------------------
  def set_matches_AA_truth(self, jet_combined, jetR, hname):
  
    set_match = False
    
    # Create Pb-Pb matching QA histogram, if not already created
    if not hasattr(self, hname):
      bin_labels = ['all', 'has_matching_candidate', 'has_matching_unique_candidate', 'mc_fraction', 'deltaR_combined-truth', 'passed_all_cuts']
      nbins = len(bin_labels)
      h = ROOT.TH2F(hname, hname, nbins, 0, nbins,  30, 0., 300.)
      for i in range(1, nbins+1):
        h.GetXaxis().SetBinLabel(i,bin_labels[i-1])
      setattr(self, hname, h)
    
    h = getattr(self, hname)
    h.Fill('all', jet_combined.pt(), 1)

    # Check if combined jet has a pp-det match (uniqueness not currently enforced)
    if jet_combined.has_user_info():
    
      h.Fill('has_matching_candidate', jet_combined.pt(), 1)

      jet_info_combined = jet_combined.python_info()
      jet_truth = jet_info_combined.closest_jet
      set_match = True
        
      # Check that match is unique
      if self.is_match_unique(jet_truth):
        h.Fill('has_matching_unique_candidate', jet_combined.pt(), 1)
      else:
        set_match = False
          
      # Check matching distance between combined jet and pp-truth
      if jet_combined.delta_R(jet_truth) < self.jet_matching_distance*jetR:
        h.Fill('deltaR_combined-truth', jet_combined.pt(), 1)
      else:
        set_match = False

      # Check if >50% of the pp-truth tracks are in the combined jet
      mc_fraction = self.mc_fraction(jet_truth, jet_combined)
      if mc_fraction > self.mc_fraction_threshold:
        h.Fill('mc_fraction', jet_combined.pt(), 1)
      else:
        set_match = False
             
    # Set accepted match
    if set_match:
      jet_info_combined.match = jet_truth
      jet_combined.set_python_info(jet_info_combined)
      h.Fill('passed_all_cuts', jet_combined.pt(), 1)

  #---------------------------------------------------------------
  # Return pt-fraction of tracks in jet_pp_det that are contained in jet_det_combined
  #---------------------------------------------------------------
  def mc_fraction(self, jet_pp_det, jet_det_combined):
  
    pt_total = jet_pp_det.pt()
    
    pt_contained = 0.
    for track in jet_det_combined.constituents():
      if track.user_index() >= 0:
        pt_contained += track.pt()
            
    return pt_contained/pt_total

  #---------------------------------------------------------------
  # Return whether a jet has a unique match
  #---------------------------------------------------------------
  def is_match_unique(self, jet_pp_det):
  
    if jet_pp_det.has_user_info():

      jet_info_pp_det = jet_pp_det.python_info()
      jet_pp_truth = jet_info_pp_det.closest_jet

      if len(jet_info_pp_det.matching_candidates) == 1:

        if jet_pp_truth.has_user_info():
      
          jet_info_pp_truth = jet_pp_truth.python_info()
          if len(jet_info_pp_truth.matching_candidates) == 1:
            return True
                  
    return False

  #---------------------------------------------------------------
  # Save all histograms
  #---------------------------------------------------------------
  def save_output_objects(self):
    
    outputfilename = os.path.join(self.output_dir, 'AnalysisResults.root')
    fout = ROOT.TFile(outputfilename, 'update')
    fout.cd()
    
    for attr in dir(self):
      
      obj = getattr(self, attr)
      
      # Write all ROOT histograms and trees to file
      types = (ROOT.TH1, ROOT.THnBase, ROOT.TTree)
      if isinstance(obj, types):
        obj.Write()
  
    fout.Close()

  #---------------------------------------------------------------
  # Save all THn and TH3, and remove them as class attributes (to clear memory)
  #---------------------------------------------------------------
  def save_thn_th3_objects(self):
    
    outputfilename = os.path.join(self.output_dir, 'AnalysisResults.root')
    fout = ROOT.TFile(outputfilename, 'update')
    fout.cd()
    
    for attr in dir(self):
      
      obj = getattr(self, attr)
      
      types = (ROOT.TH3, ROOT.THnBase)
      if isinstance(obj, types):
        obj.Write()
        delattr(self, attr)

    fout.Close()
