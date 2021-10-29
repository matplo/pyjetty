#! /usr/bin/env python

""" theory_comp.py
Loads theory comparisons, preforms un/folding, makes plots
Ezra Lesser, 2020 (elesser@berkeley.edu)
"""

import sys
import os
import argparse
from array import *
import numpy as np
import math   # for exp()
import ROOT
ROOT.gSystem.Load("$HEPPY_DIR/external/roounfold/roounfold-current/lib/libRooUnfold.so")
import yaml

# For log y plots which ROOT just decides not to work for
#import matplotlib
#matplotlib.rcParams['text.usetex'] = True   # LaTeX labels
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams["yaxis.labellocation"] = 'top'
plt.rcParams["xaxis.labellocation"] = 'right'

from pyjetty.alice_analysis.analysis.user.substructure import run_analysis

# Load pyjetty ROOT utils
ROOT.gSystem.Load('libpyjetty_rutil')

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)


################################################################
# Helper functions
################################################################

#----------------------------------------------------------------------
# Extrapolate y-values for values in xlist_new given points (x,y) in xlist and ylist
# Use power=1 for linear, or power=2 for quadratic extrapolation
#----------------------------------------------------------------------
def list_interpolate(xlist, ylist, xlist_new, power=1, require_positive=False):

  if len(xlist) < (power + 1):
    raise ValueError("list_interpolate() requires at least %i points!" % (power + 1))

  ylist_new = []
  ix = 0
  for xval in xlist_new:

    while (ix + power) < len(xlist) and xlist[ix+power] <= xval:
      ix += 1

    x1 = xlist[ix]; y1 = ylist[ix]

    # Check if data point is identical
    if xval == x1:
      if require_positive and y1 < 0:
        ylist_new.append(0)
        continue
      ylist_new.append(y1)
      continue

    # Set value to 0 if out-of-range for extrapolation
    if x1 > xval or (ix + power) >= len(xlist):
      ylist_new.append(0)
      continue

    x2 = xlist[ix+1]; y2 = ylist[ix+1]

    yval = None
    if power == 1:  # linear
      yval = linear_extrapolate(x1, y1, x2, y2, xval)
    elif power == 2:  # quadratic
      x3 = xlist[ix+2]; y3 = ylist[ix+2]
      yval = quadratic_extrapolate(x1, y1, x2, y2, x3, y3, xval)
    else:
      raise ValueError("Unrecognized power", power, "/ please use either 1 or 2")

    # Require positive values
    if require_positive and yval < 0:
      ylist_new.append(0)
      continue

    ylist_new.append(yval)

  return ylist_new


#---------------------------------------------------------------
# Given two data points, find linear fit and y-value for x
#---------------------------------------------------------------
def linear_extrapolate(x1, y1, x2, y2, x):

  return (y2 - y1) / (x2 - x1) * x + (y1 - (y2 - y1) / (x2 - x1) * x1)


#---------------------------------------------------------------
# Given three data points, find quadratic fit and y-value for x
#---------------------------------------------------------------
def quadratic_extrapolate(x1, y1, x2, y2, x3, y3, x):

  a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2))
  b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3))
  c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3)

  return (a * x * x + b * x + c) / ((x1 - x2) * (x1 - x3) * (x2 - x3))


#---------------------------------------------------------------
# Set LHS of distributions to 0 if crosses to 0 at some point (prevents multiple peaks)
#---------------------------------------------------------------
def set_zero_range(yvals):

  found_nonzero_val = False

  # Step through list backwards
  for i in range(len(yvals)-1, -1, -1):
    if yvals[i] <= 0:
      if found_nonzero_val:
        for j in range(0, i+1):
          yvals[j] = 0
        break
      yvals[i] = 0
      continue
    else:
      found_nonzero_val = True
      continue

  return yvals


# Where there are single values pos/neg between two neg/pos, interpolate point
def fix_fluctuations(yvals):

  for i in range(1, len(yvals) - 1):
    if yvals[i] > 0:
      if yvals[i+1] < 0 and yvals[i-1] < 0:
        yvals[i] = (yvals[i+1] + yvals[i-1]) / 2
    else:  # yvals[i] <= 0
      if yvals[i+1] > 0 and yvals[i-1] > 0:
        yvals[i] = (yvals[i+1] + yvals[i-1]) / 2

  return yvals


################################################################
#######################  RUN ANALYSIS  #########################
################################################################
class RunAnalysisAng(run_analysis.RunAnalysis):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(RunAnalysisAng, self).__init__(config_file, **kwargs)
    
    # Initialize yaml config
    self.initialize_user_config()

    print(self)
  
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_user_config(self):
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
      
    self.figure_approval_status = config['figure_approval_status']
    self.plot_overlay_list = \
      self.obs_config_dict['common_settings']['plot_overlay_list']
    
    self.jet_matching_distance = config['jet_matching_distance']
    
    if 'constituent_subtractor' in config:
        self.is_pp = False
    else:
        self.is_pp = True
    print('is_pp: {}'.format(self.is_pp))

    # Whether or not to use the previous preliminary result in final plots
    self.use_prev_prelim = config['use_prev_prelim']

    self.histutils = ROOT.RUtil.HistUtils()

    # Grooming settings
    self.sd_zcut = config["sd_zcut"]
    self.sd_beta = config["sd_beta"]
    self.theory_grooming_settings = [{'sd': [self.sd_zcut, self.sd_beta]}]  # self.utils.grooming_settings
    self.theory_grooming_labels = [self.utils.grooming_label(gs) for gs in
                                   self.theory_grooming_settings]

    # Theory comparisons
    if 'fPythia' in config:
      self.fPythia_name = config['fPythia']

    if 'theory_dir' in config:
      self.do_theory = config['do_theory_comp']
      if self.do_theory:

        self.theory_dir = config['theory_dir']
        self.theory_alpha = config['theory_alpha']
        self.theory_pt_bins = config['theory_pt_bins']
        self.theory_pt_bins_center = [(self.theory_pt_bins[i] + self.theory_pt_bins[i+1]) / 2 for \
                                      i in range(len(self.theory_pt_bins)-1)]
        self.theory_response_files = [ROOT.TFile(f, 'READ') for f in config['response_files']]
        self.theory_response_labels = config['response_labels']
        self.theory_pt_scale_factors_filepath = os.path.join(
          self.theory_dir, config['pt_scale_factors_filename'])
        self.rebin_theory_response = config['rebin_theory_response']
        self.output_dir_theory = os.path.join(self.output_dir, self.observable, 'theory_response')
        self.Lambda = 1  # GeV -- This variable changes the NP vs P region of theory plots

        self.do_theory_F_np = True  # NP shape f'n convolution
        if self.do_theory_F_np:
          self.Omega_list = config['Omega_list']  # list of universal(?) parameters to try
          if len(self.Omega_list) <= 0:
            self.do_theory_F_np = False

        # Define observable binnings -- may want to move to config file eventually
        self.theory_obs_bins = np.concatenate((
          np.linspace(0, 0.009, 10), np.linspace(0.01, 0.1, 19), np.linspace(0.11, 0.8, 70)))
        self.theory_obs_bins_center = np.concatenate(
          (np.linspace(0.0005, 0.0095, 10), np.linspace(0.0125, 0.0975, 18),
           np.linspace(0.105, 0.795, 70)))
        self.theory_obs_bins_width = 10 * [0.001] + 18 * [0.005] + 70 * [0.01]

        # Use the old theory prediction binnings as test (ungroomed only)
        self.use_old = False
        if self.use_old:
          self.theory_grooming_settings = []
          self.theory_grooming_labels = []
          self.theory_pt_bins = list(range(10, 160, 10))
          #self.theory_obs_bins = np.linspace(0, 0.8, 81)
          #self.theory_obs_bins_center = np.linspace(0.005, 0.795, 80)
          #self.theory_obs_bins_width = 80 * [0.01]

    else:
      self.do_theory = False

    if self.do_theory:
      self.load_pt_scale_factors(self.theory_pt_scale_factors_filepath)
      print("Loading response matrix for folding theory predictions...")
      self.load_theory_response()
      print("Loading theory histograms...")
      self.load_theory_histograms()  # Loads and folds from parton -> CH level


  #---------------------------------------------------------------
  # Load 4D response matrices used for folding the theory predictions
  #---------------------------------------------------------------
  def load_theory_response(self):

    # Check to see if Roounfold file already exists
    if not os.path.exists(self.output_dir_theory):
      os.makedirs(self.output_dir_theory)
    roounfold_filename = os.path.join(self.output_dir_theory, 'fRoounfold.root')
    roounfold_exists = os.path.exists(roounfold_filename)

    # Do the grooming if desired
    label_gr = None
    gs = None; gl = None
    if len(self.theory_grooming_settings) == 1:
      gs = self.theory_grooming_settings[0]
      gl = self.theory_grooming_labels[0]
    elif len(self.theory_grooming_settings) > 1:
      raise NotImplementedError("Not implemented for more than one grooming setting.")

    for jetR in self.jetR_list:
      for alpha in self.theory_alpha:
        label = "R%s_%s" % (str(jetR).replace('.', ''), str(alpha).replace('.', ''))
        if gs:
          label_gr = label + '_' + gl

        for ri, response in enumerate(self.theory_response_files):
          # Load charged hadron level folding response matrix
          name_ch = "hResponse_JetPt_%s_ch_%sScaled" % (self.observable, label)
          thn_ch = response.Get(name_ch)
          name_ch = "hResponse_theory_ch_%s" % label
          setattr(self, '%s_%i' % (name_ch, ri), thn_ch)
          name_ch_gr = None; thn_ch_gr = None;
          if gs:
            name_ch_gr = "hResponse_JetPt_%s_ch_%sScaled" % (self.observable, label_gr)
            thn_ch_gr = response.Get(name_ch_gr)
            name_ch_gr = "hResponse_theory_ch_%s" % label_gr
            setattr(self, '%s_%i' % (name_ch_gr, ri), thn_ch_gr)

          # Load hadron-level folding response matrix (for comparison histograms)
          name_h = "hResponse_JetPt_%s_h_%sScaled" % (self.observable, label)
          thn_h = response.Get(name_h)
          name_h = "hResponse_theory_h_%s" % label
          setattr(self, '%s_%i' % (name_h, ri), thn_h)
          name_h_gr = None; thn_h_gr = None;
          if gs:
            name_h_gr = "hResponse_JetPt_%s_h_%sScaled" % (self.observable, label_gr)
            thn_h_gr = response.Get(name_h_gr)
            name_h_gr = "hResponse_theory_h_%s" % label_gr
            setattr(self, '%s_%i' % (name_h_gr, ri), thn_h_gr)

          # Load response for H --> CH with F_np-convolved predictions
          name_Fnp = None; thn_Fnp = None; name_Fnp_gr = None; thn_Fnp_gr = None;
          if self.do_theory_F_np and ri == 0:
            name_Fnp = "hResponse_JetPt_%s_Fnp_%sScaled" % (self.observable, label)
            thn_Fnp = response.Get(name_Fnp)
            name_Fnp = "hResponse_theory_Fnp_%s" % label
            setattr(self, '%s_%i' % (name_Fnp, ri), thn_Fnp)
            if gs:
              name_Fnp_gr = "hResponse_JetPt_%s_Fnp_%sScaled" % (self.observable, label_gr)
              thn_Fnp_gr = response.Get(name_Fnp_gr)
              name_Fnp_gr = "hResponse_theory_Fnp_%s" % label_gr
              setattr(self, '%s_%i' % (name_Fnp_gr, ri), thn_Fnp_gr)

          # Create Roounfold object
          name_roounfold_h = '%s_Roounfold_%i' % (name_h, ri)
          name_roounfold_ch = '%s_Roounfold_%i' % (name_ch, ri)
          name_roounfold_h_gr = None; name_roounfold_ch_gr = None;
          if gs:
            name_roounfold_h_gr = '%s_Roounfold_%i' % (name_h_gr, ri)
            name_roounfold_ch_gr = '%s_Roounfold_%i' % (name_ch_gr, ri)
          name_roounfold_Fnp = None; name_roounfold_Fnp = None;
          if self.do_theory_F_np and ri == 0:
            name_roounfold_Fnp = '%s_Roounfold_%i' % (name_Fnp, ri)
            if gs:
              name_roounfold_Fnp_gr = '%s_Roounfold_%i' % (name_Fnp_gr, ri)

          if roounfold_exists and not self.rebin_theory_response:
            fRoo = ROOT.TFile(roounfold_filename, 'READ')
            roounfold_response_ch = fRoo.Get(name_roounfold_ch)
            roounfold_response_h = fRoo.Get(name_roounfold_h)
            roounfold_response_ch_gr = None; roounfold_response_h_gr = None;
            if gs:
              roounfold_response_ch_gr = fRoo.Get(name_roounfold_ch_gr)
              roounfold_response_h_gr = fRoo.Get(name_roounfold_h_gr)
            roounfold_response_Fnp = None; roounfold_response_Fnp_gr = None;
            if self.do_theory_F_np and ri == 0:
              roounfold_response_Fnp = fRoo.Get(name_roounfold_Fnp)
              if gs:
                roounfold_response_Fnp_gr = fRoo.Get(name_roounfold_Fnp_gr)
            fRoo.Close()

          else:  # Generated theory folding matrix needs rebinning
            # Response axes: ['p_{T}^{ch jet}', 'p_{T}^{jet, parton}', 
            #                 '#lambda_{#alpha}^{ch}', '#lambda_{#alpha}^{parton}']
            # as compared to the usual
            #      ['p_{T,det}', 'p_{T,truth}', '#lambda_{#alpha,det}', '#lambda_{#alpha,truth}']
            det_pt_bin_array = array('d', self.theory_pt_bins)
            tru_pt_bin_array = det_pt_bin_array
            obs_bins = array('d', self.theory_obs_bins)
            det_obs_bin_array = array('d', obs_bins)
            tru_obs_bin_array = det_obs_bin_array
            obs_bins_gr = None; det_obs_bin_array_gr = None; tru_obs_bin_array_gr = None;
            if gs:
              if 'sd' in gs:
                # Add bin for underflow value (tagging fraction)
                obs_bins_gr = np.insert(obs_bins, 0, -0.1)
              else: 
                obs_bins_gr = obs_bins
              det_obs_bin_array_gr = array('d', obs_bins_gr)
              tru_obs_bin_array_gr = det_obs_bin_array_gr

            n_dim = 4
            self.histutils.rebin_thn(
              roounfold_filename, thn_ch, '%s_Rebinned_%i' % (name_ch, ri), name_roounfold_ch, n_dim,
              len(det_pt_bin_array)-1, det_pt_bin_array, len(det_obs_bin_array)-1, det_obs_bin_array,
              len(tru_pt_bin_array)-1, tru_pt_bin_array, len(tru_obs_bin_array)-1, tru_obs_bin_array,
              label)
            self.histutils.rebin_thn(
              roounfold_filename, thn_h, '%s_Rebinned_%i' % (name_h, ri), name_roounfold_h, n_dim,
              len(det_pt_bin_array)-1, det_pt_bin_array, len(det_obs_bin_array)-1, det_obs_bin_array,
              len(tru_pt_bin_array)-1, tru_pt_bin_array, len(tru_obs_bin_array)-1, tru_obs_bin_array,
              label)
            if gs:
              use_underflow = 'sd' in gs
              self.histutils.rebin_thn(
                roounfold_filename, thn_ch_gr, '%s_Rebinned_%i' % (name_ch_gr, ri),
                name_roounfold_ch_gr, n_dim, len(det_pt_bin_array)-1, det_pt_bin_array,
                len(det_obs_bin_array_gr)-1, det_obs_bin_array_gr, len(tru_pt_bin_array)-1,
                tru_pt_bin_array, len(tru_obs_bin_array_gr)-1, tru_obs_bin_array_gr,
                label_gr, 0, 1, use_underflow)
              self.histutils.rebin_thn(
                roounfold_filename, thn_h_gr, '%s_Rebinned_%i' % (name_h_gr, ri),
                name_roounfold_h_gr, n_dim, len(det_pt_bin_array)-1, det_pt_bin_array,
                len(det_obs_bin_array_gr)-1, det_obs_bin_array_gr, len(tru_pt_bin_array)-1,
                tru_pt_bin_array, len(tru_obs_bin_array_gr)-1, tru_obs_bin_array_gr,
                label_gr, 0, 1, use_underflow)
            if self.do_theory_F_np and ri == 0:
              self.histutils.rebin_thn(
                roounfold_filename, thn_Fnp, '%s_Rebinned_%i' % (name_Fnp, ri),
                name_roounfold_Fnp, n_dim, len(det_pt_bin_array)-1, det_pt_bin_array,
                len(det_obs_bin_array)-1, det_obs_bin_array, len(tru_pt_bin_array)-1,
                tru_pt_bin_array, len(tru_obs_bin_array)-1, tru_obs_bin_array, label)
              if gs:
                # Don't need to do the underflow here for now... #TODO
                self.histutils.rebin_thn(
                  roounfold_filename, thn_Fnp_gr, '%s_Rebinned_%i' % (name_Fnp_gr, ri),
                  name_roounfold_Fnp_gr, n_dim, len(det_pt_bin_array)-1, det_pt_bin_array,
                  len(det_obs_bin_array)-1, det_obs_bin_array, len(tru_pt_bin_array)-1,
                  tru_pt_bin_array, len(tru_obs_bin_array)-1, tru_obs_bin_array, label_gr)
              
            f_resp = ROOT.TFile(roounfold_filename, 'READ')
            roounfold_response_ch = f_resp.Get(name_roounfold_ch)
            roounfold_response_h = f_resp.Get(name_roounfold_h)
            roounfold_response_ch_gr = None; roounfold_response_h_gr = None;
            if gs:
              roounfold_response_ch_gr = f_resp.Get(name_roounfold_ch_gr)
              roounfold_response_h_gr = f_resp.Get(name_roounfold_h_gr)
            roounfold_response_Fnp = None; roounfold_response_Fnp_gr = None;
            if self.do_theory_F_np and ri == 0:
              roounfold_response_Fnp = f_resp.Get(name_roounfold_Fnp)
              if gs:
                roounfold_response_Fnp_gr = f_resp.Get(name_roounfold_Fnp_gr)
            f_resp.Close()

          setattr(self, name_roounfold_ch, roounfold_response_ch)
          setattr(self, name_roounfold_h, roounfold_response_h)
          if gs:
            setattr(self, name_roounfold_ch_gr, roounfold_response_ch_gr)
            setattr(self, name_roounfold_h_gr, roounfold_response_h_gr)
          if self.do_theory_F_np and ri == 0:
            setattr(self, name_roounfold_Fnp, roounfold_response_Fnp)
            if gs:
              setattr(self, name_roounfold_Fnp_gr, roounfold_response_Fnp_gr)


  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #---------------------------------------------------------------
  def load_theory_histograms(self):

    # Disable folding for missed tagging fraction. Useful when this is unknown or not trusted
    disable_tagging_fraction = True

    # Set central value to test distribution. Useful for testing folding resilience
    self.exp_test = False

    # Require that hard scale and jet scale are varied together. Changes theory uncertainties
    scale_req = False

    # Do the grooming if desired
    label_gr = None
    gs = None; gl = None
    if not self.use_old and len(self.theory_grooming_settings) == 1:
      gs = self.theory_grooming_settings[0]
      gl = self.theory_grooming_labels[0]
    elif len(self.theory_grooming_settings) > 1:
      raise NotImplementedError("Not implemented for more than one grooming setting.")

    pt_bins = array('d', self.theory_pt_bins)
    obs_bins = array('d', self.theory_obs_bins)
    obs_bins_center = self.theory_obs_bins_center
    obs_bins_width = self.theory_obs_bins_width

    obs_bins_Fnp = None; obs_bins_center_Fnp = None; obs_bins_width_Fnp = None
    if self.do_theory_F_np:
      # Hard code for now... since RM does not have low-ang bins
      theory_obs_bins_Fnp = np.concatenate((
        np.linspace(0, 0.0009, 10), np.linspace(0.001, 0.009, 9),
        np.linspace(0.01, 0.1, 19), np.linspace(0.11, 0.8, 70)))

      bin_width_Fnp = theory_obs_bins_Fnp[1] - theory_obs_bins_Fnp[0]
      n_bins_Fnp = round((theory_obs_bins_Fnp[-1] - theory_obs_bins_Fnp[0]) / bin_width_Fnp)
      obs_bins_width_Fnp = [bin_width_Fnp] * n_bins_Fnp
      obs_bins_Fnp = array('d', np.linspace(
        theory_obs_bins_Fnp[0], theory_obs_bins_Fnp[-1], n_bins_Fnp))
      obs_bins_center_Fnp = [(obs_bins_Fnp[i] + obs_bins_Fnp[i+1]) / 2 for \
                             i in range(len(obs_bins_Fnp)-1)]

      '''
      bin_width_Fnp = self.theory_obs_bins[1] - self.theory_obs_bins[0]
      n_bins_Fnp = round((self.theory_obs_bins[-1] - self.theory_obs_bins[0]) / bin_width_Fnp)
      obs_bins_width_Fnp = [bin_width_Fnp] * n_bins_Fnp
      obs_bins_Fnp = array('d', np.linspace(
        self.theory_obs_bins[0], self.theory_obs_bins[-1], n_bins_Fnp))
      obs_bins_center_Fnp = [(obs_bins_Fnp[i] + obs_bins_Fnp[i+1]) / 2 for \
                             i in range(len(obs_bins_Fnp)-1)]
      '''

    obs_bins_gr = None; obs_bins_Fnp_gr = None
    if gs:
      if 'sd' in gs:
        # Add extra bin for tagging fraction
        obs_bins_gr = np.insert(obs_bins, 0, -0.1)
        if self.do_theory_F_np:
          obs_bins_Fnp_gr = np.insert(obs_bins_Fnp, 0, -0.1)
      else: 
        obs_bins_gr = obs_bins
        if self.do_theory_F_np:
          obs_bins_Fnp_gr = obs_bins_Fnp

    # Create histogram for each value of R and alpha
    for jetR in self.jetR_list:
      for alpha in self.theory_alpha:   # alpha value
        label = "R%s_%s" % (str(jetR).replace('.', ''), str(alpha).replace('.', ''))
        if gs:
          label_gr = label + '_' + gl

        name_cent = "theory_cent_%s_%s_parton" % (self.observable, label)
        name_min = "theory_min_%s_%s_parton" % (self.observable, label)
        hist_min = ROOT.TH2D(name_min, name_min, len(pt_bins)-1, 
                             pt_bins, len(obs_bins)-1, obs_bins)
        name_max = "theory_max_%s_%s_parton" % (self.observable, label)
        hist_max = ROOT.TH2D(name_max, name_max, len(pt_bins)-1, 
                             pt_bins, len(obs_bins)-1, obs_bins)

        parton_hists = ( ([], [], []), ([], [], []), ([], [], []) )
        parton_hists_Fnp = None
        if self.do_theory_F_np:
          parton_hists_Fnp = ( ([], [], []), ([], [], []), ([], [], []) )

        name_cent_gr = None; name_min_gr = None; name_max_gr = None
        hist_min_gr = None; hist_max_gr = None
        parton_hists_gr = None
        if gs:
          name_cent_gr = "theory_cent_%s_%s_parton" % (self.observable, label_gr)
          name_min_gr = "theory_min_%s_%s_parton" % (self.observable, label_gr)
          hist_min_gr = ROOT.TH2D(name_min_gr, name_min_gr, len(pt_bins)-1,
                                  pt_bins, len(obs_bins_gr)-1, obs_bins_gr)
          name_max_gr = "theory_max_%s_%s_parton" % (self.observable, label_gr)
          hist_max_gr = ROOT.TH2D(name_max_gr, name_max_gr, len(pt_bins)-1,
                                  pt_bins, len(obs_bins_gr)-1, obs_bins_gr)

          parton_hists_gr = ( ([], [], []), ([], [], []), ([], [], []) )
          parton_hists_Fnp_gr = None
          if self.do_theory_F_np:
            parton_hists_Fnp_gr = ( ([], [], []), ([], [], []), ([], [], []) )

        for l in range(0, 3):
          for m in range(0, 3):
            for n in range(0, 3):

              name_hist = "theory_%i%i%i_%s_%s_parton" % (l, m, n, self.observable, label)
              hist = ROOT.TH2D(name_hist, name_hist, len(pt_bins)-1,
                               pt_bins, len(obs_bins)-1, obs_bins)
              name_hist_Fnp = None; hist_Fnp = None
              if self.do_theory_F_np:
                name_hist_Fnp = "theory_%i%i%i_%s_%s_parton_Fnp" % (l, m, n, self.observable, label)
                hist_Fnp = ROOT.TH2D(name_hist_Fnp, name_hist_Fnp, len(pt_bins)-1,
                                     pt_bins, len(obs_bins_Fnp)-1, obs_bins_Fnp)

              name_hist_gr = None; hist_gr = None
              name_hist_Fnp_gr = None; hist_Fnp_gr = None
              if gs:
                name_hist_gr = "theory_%i%i%i_%s_%s_parton" % (l, m, n, self.observable, label_gr)
                hist_gr = ROOT.TH2D(name_hist_gr, name_hist_gr, len(pt_bins)-1,
                                    pt_bins, len(obs_bins_gr)-1, obs_bins_gr)
                if self.do_theory_F_np:
                  name_hist_Fnp_gr = "theory_%i%i%i_%s_%s_parton_Fnp" % \
                                     (l, m, n, self.observable, label_gr)
                  hist_Fnp_gr = ROOT.TH2D(name_hist_Fnp_gr, name_hist_Fnp_gr, len(pt_bins)-1,
                                          pt_bins, len(obs_bins_Fnp_gr)-1, obs_bins_Fnp_gr)

              if (scale_req and m != n) or (0 in (l, m, n) and 2 in (l, m, n)):
                parton_hists[l][m].append(None)
                if self.do_theory_F_np:
                  parton_hists_Fnp[l][m].append(None)
                if gs:
                  parton_hists_gr[l][m].append(None)
                  if self.do_theory_F_np:
                    parton_hists_Fnp_gr[l][m].append(None)
                continue

              # Loop through each pT-bin
              for i, pt_min in enumerate(self.theory_pt_bins[0:-1]):
                pt_max = self.theory_pt_bins[i+1]

                # Get scale factor for this pT bin.
                # This reverses the self-normalization of 1/sigma for correct pT scaling
                #     when doing projections onto the y-axis.
                scale_f = self.pt_scale_factor_jetR(pt_min, pt_max, jetR)

                th_dir = None
                if not self.use_old:
                  th_dir = os.path.join(
                    self.theory_dir, "ungr_ALICE_R%s" % str(jetR).replace('.', ''), 
                    "alpha%s" % str(alpha).replace('.', 'p'), "pT%s_%s" % (pt_min, pt_max))
                else:
                  th_dir = os.path.join(
                    self.theory_dir, "old", "R%s" % str(jetR).replace('.', ''),
                    "pT%s_%s" % (pt_min, pt_max), "alpha%s" % str(alpha).replace('.', 'p'))

                th_dir_gr = None
                if gs:  # != None:
                  th_dir_gr = os.path.join(
                    self.theory_dir, "gr_ALICE_R%s" % str(jetR).replace('.', ''), 
                    "alpha%s" % str(alpha).replace('.', 'p'), "pT%s_%s" % (pt_min, pt_max))

                val_li = None; val_li_gr = None; val_li_Fnp = None; val_li_Fnp_gr = None
                if self.exp_test:
                  val_li = [1 * obs_bins_width[i] for i in range(len(obs_bins_center))]
                  if self.do_theory_F_np:
                    val_li_Fnp = [1 * obs_bins_width_Fnp[i] for i in range(len(obs_bins_center_Fnp))]
                  if gs:
                    val_li_gr = [0] + val_li
                    if self.do_theory_F_np:
                      val_li_Fnp_gr = [0] + val_li_Fnp
                  #np.exp(np.linspace(0, 1, 101, True))
                  #val = np.concatenate((np.full(51, 1), np.full(50, 0)))
                  #val = [0.6 - x for x in np.linspace(0, 1, 101, True)]

                else:  # Load un/groomed predictions from files
                  x_val_li = None; y_val_li = None
                  # Load theory predictions for lambda values
                  filetype = None
                  if self.use_old:
                    filetype = ".txt"
                  else:
                    filetype = ".dat"
                  with open(os.path.join(th_dir, "%i%i%i%s" % (l, m, n, filetype))) as f:
                    lines = [line for line in f.read().split('\n') if line]
                    x_val_li = [float(line.split()[0]) for line in lines]
                    y_val_li = fix_fluctuations(
                      [float(line.split()[1]) if str(line.split()[1]).lower() != 'nan' else -1 \
                       for line in lines] )

                  # Interpolate parton curve to all bins and set 0 range on LHS tail
                  # Scale by bin width (to match RM)
                  power = 1
                  val_li = [val * obs_bins_width[i] for i, val in enumerate(set_zero_range(
                    list_interpolate(x_val_li, y_val_li, obs_bins_center,
                                     power=power, require_positive=True)))]
                  if self.do_theory_F_np:
                    val_li_Fnp = [val * obs_bins_width_Fnp[i] for i, val in enumerate(set_zero_range(
                      list_interpolate(x_val_li, y_val_li, obs_bins_center_Fnp,
                                       power=power, require_positive=True)))]

                  if gs:
                    y_val_li_gr = None
                    with open(os.path.join(th_dir_gr, "%i%i%i.dat" % (l, m, n))) as f:
                      lines = [line for line in f.read().split('\n') if line]
                      y_val_li_gr = fix_fluctuations(
                        [float(line.split()[1]) if str(line.split()[1]).lower() != 'nan' else -1 \
                         for line in lines] )

                    if "sd" in gs:
                      # Copy tail of ungroomed distribution onto groomed one for lambda > z_cut
                      # NOTE: assumes x_val_li is identical in groomed & ungroomed cases
                      k = 0
                      while k < len(x_val_li) and x_val_li[k] < gs["sd"][0]:
                        k += 1
                      if y_val_li[k] <= 0:
                        y_val_li_gr[k:] = [0] * len(y_val_li_gr[k:])
                      elif y_val_li_gr[k] <= 0:
                        y_val_li_gr[k:] = y_val_li[k:]
                      else:
                        # Prefactor is to rescale tail to match correctly
                        y_val_li_gr[k:] = [(y_val_li_gr[k] / y_val_li[k]) * value for
                                           value in y_val_li[k:]]

                    # Interpolate parton curve to all bins and set 0 range on LHS tail
                    # Scale by bin width (to match RM)
                    val_li_gr = [val * obs_bins_width[i] for i, val in enumerate(set_zero_range(
                      list_interpolate(x_val_li, y_val_li_gr, obs_bins_center,
                                       power=power, require_positive=True)))]
                    if self.do_theory_F_np:
                      val_li_Fnp_gr = [val * obs_bins_width_Fnp[i] for i, val in enumerate(
                        set_zero_range(list_interpolate(x_val_li, y_val_li_gr, obs_bins_center_Fnp,
                                         power=power, require_positive=True)))]

                  # Rescale histograms by the correct integral
                  # Simultaneously scale by scale_f to get correct jet pT scaling 
                  integral_val_li = sum(val_li)
                  val_li = [val * scale_f / integral_val_li for val in val_li]
                  integral_val_li_Fnp = None
                  if self.do_theory_F_np:
                    integral_val_li_Fnp = sum(val_li_Fnp)
                    val_li_Fnp = [val * scale_f / integral_val_li_Fnp for val in val_li_Fnp]
                  missed_tagging_fraction = None
                  if gs:
                    integral_val_li_gr = sum(val_li_gr)
                    integral_val_li_Fnp_gr = None
                    if self.do_theory_F_np:
                      integral_val_li_Fnp_gr = sum(val_li_Fnp_gr)
                    if "sd" in gs:
                      # Get SD tagging fraction to save as underflow value
                      missed_tagging_fraction = 1 - integral_val_li_gr / integral_val_li
                      if disable_tagging_fraction:
                        missed_tagging_fraction = 0
                      elif missed_tagging_fraction < 0:
                        print("WARNING: missed tagging fraction %f < 0 (\\alpha = %s, R = %s)." % \
                              (missed_tagging_fraction, alpha, jetR), "Manually setting to 0.")
                        missed_tagging_fraction = 0
                      elif missed_tagging_fraction > 1:
                        print("WARNING: missed tagging fraction %f > 1 (\\alpha = %s, R = %s)." % \
                              (missed_tagging_fraction, alpha, jetR), "Manually setting to 0.")
                        missed_tagging_fraction = 0

                      if missed_tagging_fraction == 0:
                        # Disabled, so do "regular" (groomed) scaling
                        # Normalize to 1 by dividing out groomed integral
                        val_li_gr = [0] + [val * scale_f / integral_val_li_gr for val in val_li_gr]
                        if self.do_theory_F_np:
                          val_li_Fnp_gr = [0] + [val * scale_f / integral_val_li_Fnp_gr for
                                           val in val_li_Fnp_gr]

                      else:
                        # Add underflow bin and normalize by inclusive integral
                        val_li_gr = [missed_tagging_fraction * scale_f] + \
                                    [val * scale_f / integral_val_li for val in val_li_gr]
                        if self.do_theory_F_np:
                          val_li_Fnp_gr = [0] + [val * scale_f / integral_val_li_Fnp for
                                                 val in val_li_Fnp_gr]
                    else:
                      # Normalize to 1 by dividing out groomed integral
                      val_li_gr = [val * scale_f / integral_val_li_gr for val in val_li_gr]
                      if self.do_theory_F_np:
                        val_li_Fnp_gr = [val * scale_f / integral_val_li_Fnp_gr for
                                         val in val_li_Fnp_gr]

                for j, val in enumerate(val_li):
                  hist.SetBinContent(i+1, j+1, val)
                  hist.SetBinError(i+1, j+1, 0)
                  if l == m == n == 0:
                    hist_min.SetBinContent(i+1, j+1, val)
                    hist_min.SetBinError(i+1, j+1, 0)
                    hist_max.SetBinContent(i+1, j+1, val)
                    hist_max.SetBinError(i+1, j+1, 0)
                  elif float(val) < hist_min.GetBinContent(i+1, j+1):
                    hist_min.SetBinContent(i+1, j+1, val)
                  elif float(val) > hist_max.GetBinContent(i+1, j+1):
                    hist_max.SetBinContent(i+1, j+1, val)

                if self.do_theory_F_np:
                  for j, val in enumerate(val_li_Fnp):
                    hist_Fnp.SetBinContent(i+1, j+1, val)
                    hist_Fnp.SetBinError(i+1, j+1, 0)

                if gs:
                  for j, val in enumerate(val_li_gr):
                    hist_gr.SetBinContent(i+1, j+1, val)
                    hist_gr.SetBinError(i+1, j+1, 0)
                    if l == m == n == 0:
                      hist_min_gr.SetBinContent(i+1, j+1, val)
                      hist_min_gr.SetBinError(i+1, j+1, 0)
                      hist_max_gr.SetBinContent(i+1, j+1, val)
                      hist_max_gr.SetBinError(i+1, j+1, 0)
                    elif float(val) < hist_min_gr.GetBinContent(i+1, j+1):
                      hist_min_gr.SetBinContent(i+1, j+1, val)
                    elif float(val) > hist_max_gr.GetBinContent(i+1, j+1):
                      hist_max_gr.SetBinContent(i+1, j+1, val)

                  if self.do_theory_F_np:
                    for j, val in enumerate(val_li_Fnp_gr):
                      hist_Fnp_gr.SetBinContent(i+1, j+1, val)
                      hist_Fnp_gr.SetBinError(i+1, j+1, 0)

              parton_hists[l][m].append(hist)
              if self.do_theory_F_np:
                parton_hists_Fnp[l][m].append(hist_Fnp)
              if gs:
                parton_hists_gr[l][m].append(hist_gr)
                if self.do_theory_F_np:
                  parton_hists_Fnp_gr[l][m].append(hist_Fnp_gr)

        setattr(self, name_cent, parton_hists[1][1][1])
        setattr(self, name_min, hist_min)
        setattr(self, name_max, hist_max)

        if gs:
          setattr(self, name_cent_gr, parton_hists_gr[1][1][1])
          setattr(self, name_min_gr, hist_min_gr)
          setattr(self, name_max_gr, hist_max_gr)

        # Fold from parton to CH level and scale by MPI
        print("Folding theory predictions...")
        self.fold_theory(jetR, alpha, parton_hists, scale_req)
        if gs:
          print("Folding theory predictions with %s..." % gl.replace('_', ' '))
          self.fold_theory(jetR, alpha, parton_hists_gr, scale_req, gl)

        # Also do NP shape function predictions + fold from H to CH level
        if self.do_theory_F_np:
          print("Applying NP shape function...")
          self.apply_np_shape_fn(jetR, alpha, parton_hists_Fnp, scale_req)
          if gs:
            print("Applying NP shape function with %s..." % gl.replace('_', ' '))
            self.apply_np_shape_fn(jetR, alpha, parton_hists_Fnp_gr, scale_req, gl)


  #----------------------------------------------------------------------
  # Fold theoretical predictions
  #----------------------------------------------------------------------
  def fold_theory(self, jetR, alpha, parton_hists, scale_req, grooming_label=None):

    label = "R%s_%s" % (str(jetR).replace('.', ''), str(alpha).replace('.', ''))
    if grooming_label:
      label += '_' + grooming_label

    for ri, response in enumerate(self.theory_response_files):

      # Load parton-to-charged-hadron response matrix
      response_name_ch = "hResponse_theory_ch_%s_Roounfold_%i" % (label, ri)
      response_ch = getattr(self, response_name_ch)
      response_name_h = "hResponse_theory_h_%s_Roounfold_%i" % (label, ri)
      response_h = getattr(self, response_name_h)

      folded_ch_hists = ( ([], [], []), ([], [], []), ([], [], []) )
      folded_h_hists = ( ([], [], []), ([], [], []), ([], [], []) )

      for l in range(0, 3):
        for m in range(0, 3):
          for n in range(0, 3):

            if (scale_req and m != n) or (0 in (l, m, n) and 2 in (l, m, n)):
              folded_h_hists[l][m].append(None)
              folded_ch_hists[l][m].append(None)
              continue

            # Fold theory predictions
            h_folded_ch = response_ch.ApplyToTruth(parton_hists[l][m][n])
            h_folded_h = response_h.ApplyToTruth(parton_hists[l][m][n])

            name_ch = "theory_%i%i%i_%s_%s_ch_%i" % (l, m, n, self.observable, label, ri)
            name_h = "theory_%i%i%i_%s_%s_h_%i" % (l, m, n, self.observable, label, ri)

            h_folded_ch.SetNameTitle(name_ch, name_ch)
            h_folded_h.SetNameTitle(name_h, name_h)

            folded_ch_hists[l][m].append(h_folded_ch)
            folded_h_hists[l][m].append(h_folded_h)

      printstring = "Scaling theory predictions for MPI effects for %s" % \
                    self.theory_response_labels[ri]
      if grooming_label:
        printstring += " with %s..." % grooming_label.replace('_', ' ')
      else:
        printstring += "..."
      print(printstring)
      self.mpi_scale_theory(jetR, alpha, ri, response, folded_ch_hists, folded_h_hists,
                            scale_req, grooming_label)


  #----------------------------------------------------------------------
  # Fold theoretical predictions
  #----------------------------------------------------------------------
  def mpi_scale_theory(self, jetR, alpha, ri, response, folded_ch_hists, folded_h_hists,
                       scale_req, grooming_label=None):

    label = "R%s_%s" % (str(jetR).replace('.', ''), str(alpha).replace('.', ''))
    using_sd_grooming = False
    if grooming_label:
      label += '_' + grooming_label
      using_sd_grooming = "sd" in grooming_label.lower()

    # Load parton-level theory predictions
    name_cent = "theory_cent_%s_%s" % (self.observable, label)
    name_min = "theory_min_%s_%s" % (self.observable, label)
    name_max = "theory_max_%s_%s" % (self.observable, label)

    h_cent = getattr(self, name_cent+"_parton")
    h_min = getattr(self, name_min+"_parton")
    h_max = getattr(self, name_max+"_parton")

    name_mpi_off = "hAng_JetPt_ch_%sScaled" % label
    name_mpi_on = "hAng_JetPt_ch_MPIon_%sScaled" % label

    h_mpi_off = response.Get(name_mpi_off)
    h_mpi_on = response.Get(name_mpi_on)

    '''
    # Assert that the MPI histograms have the same binning as the folded histograms
    if obs_bins != h_mpi_off_y_bins or obs_bins != h_mpi_on_y_bins:
      raise ValueError("MPI scaling histograms do not match the expected binning and/or values.")
    '''

    # Set MPI histograms to have the same binning as folded histograms
    x_bins = array('d', self.theory_pt_bins)
    n_x_bins = len(x_bins) - 1
    y_bins = array('d', self.theory_obs_bins)
    if using_sd_grooming:
      y_bins = np.insert(y_bins, 0, -0.1)
    h_mpi_off = self.histutils.rebin_th2(h_mpi_off, name_mpi_off+'_Rebinned_%i' % ri, x_bins, n_x_bins,
                                         y_bins, len(y_bins)-1, using_sd_grooming)
    h_mpi_on = self.histutils.rebin_th2(h_mpi_on, name_mpi_on+'_Rebinned_%i' % ri, x_bins, n_x_bins,
                                        y_bins, len(y_bins)-1, using_sd_grooming)

    mpi_bin_edges = [h_mpi_on.GetXaxis().GetBinLowEdge(i+1)
                     for i in range(h_mpi_on.GetNbinsX()+1)]
    mpi_bin_edge_indices = [mpi_bin_edges.index(pt) for pt in self.pt_bins_reported]
    theory_pt_bin_edge_indices = [self.theory_pt_bins.index(pt) for pt in self.pt_bins_reported]
    # Loop through each reported pT-bin
    for index, i in list(enumerate(theory_pt_bin_edge_indices))[0:-1]:
      j = theory_pt_bin_edge_indices[index+1]
      pt_min = self.theory_pt_bins[i]
      pt_max = self.theory_pt_bins[j]

      # Get scale factor for this pT bin.
      # This reverses the self-normalization of 1/sigma for correct pT scaling
      #     when doing doing projections onto the y-axis.
      scale_f = self.pt_scale_factor_k(pt_min, pt_max, -5)

      # Note: bins in ROOT are 1-indexed (0 bin is underflow). Also last bin is inclusive
      h_mpi_off_proj = h_mpi_off.ProjectionY(
        '%s_py_%s_%i' % (name_mpi_off, str(i), ri),
        mpi_bin_edge_indices[index]+1, mpi_bin_edge_indices[index+1])
      h_mpi_on_proj = h_mpi_on.ProjectionY(
        '%s_py_%s_%i' % (name_mpi_on, str(i), ri),
        mpi_bin_edge_indices[index]+1, mpi_bin_edge_indices[index+1])

      # Create ratio plot for scaling
      h_mpi_ratio = h_mpi_on_proj.Clone()
      title = 'h_theory_mpi_scaling_%s_%s_%i' % (label, str(i), ri)
      h_mpi_ratio.SetNameTitle(title, title)
      h_mpi_ratio.Divide(h_mpi_off_proj)
      h_mpi_ratio.SetDirectory(0)
      setattr(self, 'h_mpi_ratio_%s_PtBin%i-%i_%i' % (label, pt_min, pt_max, ri), h_mpi_ratio)

      pt_label = '_PtBin'+str(pt_min)+'-'+str(pt_max)

      h_cent_p_bin = h_cent.ProjectionY("%s_parton%s" % (name_cent, pt_label), i+1, j)
      h_min_p_bin = h_min.ProjectionY("%s_parton%s" % (name_min, pt_label), i+1, j)
      h_max_p_bin = h_max.ProjectionY("%s_parton%s" % (name_max, pt_label), i+1, j)

      # Keep normalization by dividing by the scale factor
      #h_cent_bin.Scale(1/(j-i))
      #h_min_bin.Scale(1/(j-i))
      #h_max_bin.Scale(1/(j-i))
      h_cent_p_bin.Scale(1/scale_f)
      h_min_p_bin.Scale(1/scale_f)
      h_max_p_bin.Scale(1/scale_f)

      # Undo the bin width scaling and set correct normalization
      norm_factor = h_cent_p_bin.Integral()
      h_cent_p_bin.Scale(1/norm_factor, "width")
      h_min_p_bin.Scale(1/norm_factor, "width")
      h_max_p_bin.Scale(1/norm_factor, "width")

      # Set new names and titles to prevent issues with saving histograms
      h_cent_p_bin.SetNameTitle("%s_parton%s" % (name_cent, pt_label),
                                "%s_parton%s" % (name_cent, pt_label))
      h_min_p_bin.SetNameTitle("%s_parton%s" % (name_min, pt_label),
                               "%s_parton%s" % (name_min, pt_label))
      h_max_p_bin.SetNameTitle("%s_parton%s" % (name_max, pt_label),
                               "%s_parton%s" % (name_max, pt_label))

      # Initialize h and ch histograms
      h_cent_ch_bin = None; h_min_ch_bin = None; h_max_ch_bin = None;
      h_cent_h_bin = None; h_min_h_bin = None; h_max_h_bin = None;

      # Initialize h and ch tagging fraction for folded distributions
      tagging_frac_cent_h = None; tagging_frac_min_h = None; tagging_frac_max_h = None;
      tagging_frac_cent_ch = None; tagging_frac_min_ch = None; tagging_frac_max_ch = None;

      # Scale the theory bin and save for plotting
      for l in range(0, 3):
        for m in range(0, 3):
          for n in range(0, 3):

            if (scale_req and m != n) or (0 in (l, m, n) and 2 in (l, m, n)):
              continue

            h_folded_h_bin = folded_h_hists[l][m][n].ProjectionY(
              name_cent + ("_h%s_%i" % (pt_label, ri)), i+1, j)
            h_folded_ch_bin = folded_ch_hists[l][m][n].ProjectionY(
              name_cent + ("_ch%s_%i" % (pt_label, ri)), i+1, j)

            #h_folded_h_bin.Scale(1/(j-i))
            #h_folded_ch_bin.Scale(1/(j-i))

            # Undo the bin width scaling and scale all by integral
            h_folded_h_bin.Scale(1/h_folded_h_bin.Integral(), "width")
            h_folded_ch_bin.Scale(1/h_folded_ch_bin.Integral(), "width")

            # Save ratio plots for seeing the change at each level
            if l == m == n == 1:
              h_cent_ratio_ch_bin = h_folded_ch_bin.Clone()
              title = 'h_cent_ratio_ch_%s%s_%i' % (label, pt_label, ri)
              h_cent_ratio_ch_bin.SetNameTitle(title, title)
              h_cent_ratio_ch_bin.Divide(h_folded_h_bin)
              h_cent_ratio_ch_bin.SetDirectory(0)
              setattr(self, title, h_cent_ratio_ch_bin)

              h_cent_ratio_h_bin = h_folded_h_bin.Clone()
              title = 'h_cent_ratio_h_%s%s_%i' % (label, pt_label, ri)
              h_cent_ratio_h_bin.SetNameTitle(title, title)
              h_cent_ratio_h_bin.Divide(h_cent_p_bin)
              h_cent_ratio_h_bin.SetDirectory(0)
              setattr(self, title, h_cent_ratio_h_bin)

            # Apply the MPI scaling to the charged distribution
            h_folded_ch_bin.Multiply(h_mpi_ratio)
            h_folded_ch_bin.Scale(1/h_folded_ch_bin.Integral("width"))

            if l == m == n == 1:
              h_cent_ch_bin = h_folded_ch_bin.Clone()
              h_cent_h_bin = h_folded_h_bin.Clone()
              h_cent_ch_bin.SetNameTitle(
                name_cent + ("_ch%s_%i" % (pt_label, ri)), name_cent + ("_ch%s_%i" % (pt_label, ri)))
              h_cent_h_bin.SetNameTitle(
                name_cent + ("_h%s_%i" % (pt_label, ri)), name_cent + ("_h%s_%i" % (pt_label, ri)))

            # Save / update the min/max/cent histograms
            if l == m == n == 0:  # Initialize the min/max histograms
              h_min_ch_bin = h_folded_ch_bin.Clone()
              h_max_ch_bin = h_folded_ch_bin.Clone()
              h_min_ch_bin.SetNameTitle(
                name_min + ("_ch%s_%i" % (pt_label, ri)), name_min + ("_ch%s_%i" % (pt_label, ri)))
              h_max_ch_bin.SetNameTitle(
                name_max + ("_ch%s_%i" % (pt_label, ri)), name_max + ("_ch%s_%i" % (pt_label, ri)))

              h_min_h_bin = h_folded_h_bin.Clone()
              h_max_h_bin = h_folded_h_bin.Clone()
              h_min_h_bin.SetNameTitle(
                name_min + ("_h%s_%i" % (pt_label, ri)), name_min + ("_h%s_%i" % (pt_label, ri)))
              h_max_h_bin.SetNameTitle(
                name_max + ("_h%s_%i" % (pt_label, ri)), name_max + ("_h%s_%i" % (pt_label, ri)))
            else:  # Update the min/max histograms
              for b in range(1, h_folded_h_bin.GetNbinsX()+1):

                h_val = h_folded_h_bin.GetBinContent(b)
                if h_val < h_min_h_bin.GetBinContent(b):
                  h_min_h_bin.SetBinContent(b, h_val)
                elif h_val > h_max_h_bin.GetBinContent(b):
                  h_max_h_bin.SetBinContent(b, h_val)

                ch_val = h_folded_ch_bin.GetBinContent(b)
                if ch_val < h_min_ch_bin.GetBinContent(b):
                  h_min_ch_bin.SetBinContent(b, ch_val)
                elif ch_val > h_max_ch_bin.GetBinContent(b):
                  h_max_ch_bin.SetBinContent(b, ch_val)


      # Steal ownership from ROOT
      h_cent_p_bin.SetDirectory(0)
      h_min_p_bin.SetDirectory(0)
      h_max_p_bin.SetDirectory(0)

      h_cent_h_bin.SetDirectory(0)
      h_min_h_bin.SetDirectory(0)
      h_max_h_bin.SetDirectory(0)

      h_cent_ch_bin.SetDirectory(0)
      h_min_ch_bin.SetDirectory(0)
      h_max_ch_bin.SetDirectory(0)

      # Save as attributes for later access (plotting)
      setattr(self, "%s_parton%s" % (name_cent, pt_label), h_cent_p_bin)
      setattr(self, "%s_parton%s" % (name_min, pt_label), h_min_p_bin)
      setattr(self, "%s_parton%s" % (name_max, pt_label), h_max_p_bin)

      setattr(self, "%s_h%s_%i" % (name_cent, pt_label, ri), h_cent_h_bin)
      setattr(self, "%s_h%s_%i" % (name_min, pt_label, ri), h_min_h_bin)
      setattr(self, "%s_h%s_%i" % (name_max, pt_label, ri), h_max_h_bin)

      setattr(self, "%s_ch%s_%i" % (name_cent, pt_label, ri), h_cent_ch_bin)
      setattr(self, "%s_ch%s_%i" % (name_min, pt_label, ri), h_min_ch_bin)
      setattr(self, "%s_ch%s_%i" % (name_max, pt_label, ri), h_max_ch_bin)


  #---------------------------------------------------------------
  # Apply NP corrections via shape function (includes MPI & hadronization)
  #---------------------------------------------------------------
  def apply_np_shape_fn(self, jetR, alpha, parton_hists, scale_req, gl=None):

    label = "R%s_%s" % (str(jetR).replace('.', ''), str(alpha).replace('.', ''))
    grooming = False
    if gl:
      label += '_' + gl
      grooming = True

    # Hard code for now... since RM does not have low-ang bins
    theory_obs_bins_Fnp = np.concatenate((
      np.linspace(0, 0.0009, 10), np.linspace(0.001, 0.009, 9),
      np.linspace(0.01, 0.1, 19), np.linspace(0.11, 0.8, 70)))

    bin_width_Fnp = theory_obs_bins_Fnp[1] - theory_obs_bins_Fnp[0]
    n_bins_Fnp = round((theory_obs_bins_Fnp[-1] - theory_obs_bins_Fnp[0]) / bin_width_Fnp)
    obs_bins_width_Fnp = [bin_width_Fnp] * n_bins_Fnp
    obs_bins_Fnp = array('d', np.linspace(
      theory_obs_bins_Fnp[0], theory_obs_bins_Fnp[-1], n_bins_Fnp))
    obs_bins_center_Fnp = [(obs_bins_Fnp[i] + obs_bins_Fnp[i+1]) / 2 for \
                           i in range(len(obs_bins_Fnp)-1)]

    '''
    bin_width_Fnp = self.theory_obs_bins[1] - self.theory_obs_bins[0]
    n_bins_Fnp = round((self.theory_obs_bins[-1] - self.theory_obs_bins[0]) / bin_width_Fnp)
    obs_bins_width_Fnp = [bin_width_Fnp] * n_bins_Fnp
    obs_bins_Fnp = array('d', np.linspace(
      self.theory_obs_bins[0], self.theory_obs_bins[-1], n_bins_Fnp))
    obs_bins_center_Fnp = [(obs_bins_Fnp[i] + obs_bins_Fnp[i+1]) / 2 for \
                           i in range(len(obs_bins_Fnp)-1)]
    '''

    for Omega in self.Omega_list:
      print("  ... Omega = %s" % str(Omega))

      # Simultaneously fold using H --> CH response (with MPI on)
      #for ri in range(len(self.theory_response_files)):
      # just use the first ri for now...
      ri = 0
      # Load hadron-to-charged-hadron response matrix
      response_name = "hResponse_theory_Fnp_%s_Roounfold_%i" % (label, ri)
      response = getattr(self, response_name)
      h_cent_ch = None; h_min_ch = None; h_max_ch = None;
      h_min_ch_name = "theory_min_%s_%s_ch_Fnp_Omega%s_%i" % \
                      (self.observable, label, str(Omega), ri)
      h_max_ch_name = "theory_max_%s_%s_ch_Fnp_Omega%s_%i" % \
                      (self.observable, label, str(Omega), ri)
      h_cent_ch_name = "theory_cent_%s_%s_ch_Fnp_Omega%s_%i" % \
                       (self.observable, label, str(Omega), ri)

      # Also save without ch folding for comparison plots
      h_cent_h = None; h_min_h = None; h_max_h = None;
      h_min_h_name = "theory_min_%s_%s_h_Fnp_Omega%s_%i" % \
                     (self.observable, label, str(Omega), ri)
      h_max_h_name = "theory_max_%s_%s_h_Fnp_Omega%s_%i" % \
                     (self.observable, label, str(Omega), ri)
      h_cent_h_name = "theory_cent_%s_%s_h_Fnp_Omega%s_%i" % \
                      (self.observable, label, str(Omega), ri)

      for l in range(0, 3):
        for m in range(0, 3):
          for n in range(0, 3):

            if (scale_req and m != n) or (0 in (l, m, n) and 2 in (l, m, n)):
              continue

            # For some reason there are some binning-dependant things happening, so
            # rebin everything to have the smallest of bins everywhere
            parton_hist = parton_hists[l][m][n].Clone()
            parton_hist.SetNameTitle(parton_hist.GetName() + "_finebin",
                                     parton_hist.GetName() + "_finebin")

            # Apply shape function
            name_h = "theory_%i%i%i_%s_%s_h_Fnp_Omega%s_%i" % \
                     (l, m, n, self.observable, label, str(Omega), ri)
            pTs = [self.pt_avg_jetR(self.theory_pt_bins[i], self.theory_pt_bins[i+1], jetR) for
                   i in range(0, len(self.theory_pt_bins) - 1)]
            h_np = self.histutils.convolve_F_np(
              Omega, jetR, alpha, array('d', obs_bins_Fnp),
              len(obs_bins_center_Fnp), array('d', obs_bins_center_Fnp),
              array('d', obs_bins_width_Fnp),
              array('d', self.theory_pt_bins), len(self.theory_pt_bins_center),
              array('d', pTs), parton_hists[l][m][n], name_h, grooming)

            # Rebin to correct binning
            # note: currently not using underflow in groomed RM here
            h_np_rebinned = self.histutils.rebin_th2(
              h_np, h_np.GetName()+"_rebinned", array('d', self.theory_pt_bins),
              len(self.theory_pt_bins)-1, array('d', self.theory_obs_bins), 
              len(self.theory_obs_bins)-1)

            # Fold hadron-level predictions to CH level
            h_folded_Fnp = response.ApplyToTruth(h_np_rebinned)
            name_ch = "theory_%i%i%i_%s_%s_ch_Fnp_Omega%s_%i" % \
                      (l, m, n, self.observable, label, str(Omega), ri)
            h_folded_Fnp.SetNameTitle(name_ch, name_ch)

            # Update min/max/cent histograms
            if l == m == n == 1: # set central variation
              h_cent_ch = h_folded_Fnp.Clone()
              h_cent_ch.SetNameTitle(h_cent_ch_name, h_cent_ch_name)
              h_cent_h = h_np_rebinned.Clone()
              h_cent_h.SetNameTitle(h_cent_h_name, h_cent_h_name)

            if l == m == n == 0:  # initialize min/max histograms
              h_min_ch = h_folded_Fnp.Clone()
              h_min_ch.SetNameTitle(h_min_ch_name, h_min_ch_name)
              h_min_h = h_np_rebinned.Clone()
              h_min_h.SetNameTitle(h_min_h_name, h_min_h_name)
              h_max_ch = h_folded_Fnp.Clone()
              h_max_ch.SetNameTitle(h_max_ch_name, h_max_ch_name)
              h_max_h = h_np_rebinned.Clone()
              h_max_h.SetNameTitle(h_max_h_name, h_max_h_name)

            else:  # loop through pT & obs bins and update each min/max value
              for i in range(h_folded_Fnp.GetNbinsX()):
                for j in range(h_folded_Fnp.GetNbinsY()):
                  val_ch = h_folded_Fnp.GetBinContent(i+1, j+1)
                  if float(val_ch) < h_min_ch.GetBinContent(i+1, j+1):
                    h_min_ch.SetBinContent(i+1, j+1, val_ch)
                  if float(val_ch) > h_max_ch.GetBinContent(i+1, j+1):
                    h_max_ch.SetBinContent(i+1, j+1, val_ch)
                  val_h = h_np_rebinned.GetBinContent(i+1, j+1)
                  if float(val_h) < h_min_h.GetBinContent(i+1, j+1):
                    h_min_h.SetBinContent(i+1, j+1, val_h)
                  if float(val_h) > h_max_h.GetBinContent(i+1, j+1):
                    h_max_h.SetBinContent(i+1, j+1, val_h)

      setattr(self, h_cent_ch_name, h_cent_ch)
      setattr(self, h_min_ch_name, h_min_ch)
      setattr(self, h_max_ch_name, h_max_ch)

      setattr(self, h_cent_h_name, h_cent_h)
      setattr(self, h_min_h_name, h_min_h)
      setattr(self, h_max_h_name, h_max_h)


  #---------------------------------------------------------------
  # Loads pT scale factors from q/g fraction theory predictions
  #---------------------------------------------------------------
  def load_pt_scale_factors(self, filepath):

    # Open file and save pT distribution
    pt_li = None; val_li = None;
    with open(filepath) as f:
      lines = [line.split() for line in f.read().split('\n') if (line and line[0] != '#')]
      pt_li = [int(float(line[0])) for line in lines]
      val_li = [float(line[1]) + float(line[2]) for line in lines]

    # Split list based on the number of measured jet R
    n_jetR = len(self.jetR_list)
    n_entries = len(val_li)
    for i, jetR in enumerate(self.jetR_list):
      val_li_jetR = val_li[i*n_entries:(i+1)*n_entries]
      pt_li_jetR = pt_li[i*n_entries:(i+1)*n_entries]
      setattr(self, "pt_scale_factor_R%s" % jetR, (pt_li_jetR, val_li_jetR))


  #---------------------------------------------------------------
  # Returns number proportional to shape of inclusive jet pT
  #     distribution theory prediction (val = jetR)
  #---------------------------------------------------------------
  def pt_scale_factor_jetR(self, ptmin, ptmax, jetR):

    pt_li, val_li = getattr(self, "pt_scale_factor_R%s" % jetR)

    # Fit a log function between the two endpoints and approx avg integral for bin
    start_i = pt_li.index(ptmin)
    end_i = pt_li.index(ptmax)
    # y = a * x^k
    k = np.log(val_li[start_i] / val_li[end_i]) / np.log(ptmin / ptmax)
    a = val_li[start_i] / ptmin**k
    return self.pt_scale_factor_k(ptmin, ptmax, k, a)


  #---------------------------------------------------------------
  # Returns approximate avg pT for bin given min and max pT of that bin
  #---------------------------------------------------------------
  def pt_avg_jetR(self, ptmin, ptmax, jetR):

    pt_li, val_li = getattr(self, "pt_scale_factor_R%s" % jetR)

    # Fit a log function between the two endpoints and approx avg integral for bin
    start_i = pt_li.index(ptmin)
    end_i = pt_li.index(ptmax)
    # y = a * x^k
    k = np.log(val_li[start_i] / val_li[end_i]) / np.log(ptmin / ptmax)
    a = val_li[start_i] / ptmin**k
    y = self.pt_scale_factor_k(ptmin, ptmax, k, a) / (ptmax - ptmin)
    return (y / a) ** (1 / k)


  #---------------------------------------------------------------
  # Returns number proportional to the integral of power law pTjet^k
  #---------------------------------------------------------------
  def pt_scale_factor_k(self, ptmin, ptmax, k, a=1e9):
    if k == -1:
      return a * np.log(ptmax / ptmin)
    return a * (ptmax**(k + 1) - ptmin**(k + 1)) / (k + 1)


  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #---------------------------------------------------------------
  def plot_single_result(self, jetR, obs_label, obs_setting, grooming_setting):
    #print('Plotting each individual result...')

    # Plot final result for each 1D substructure distribution (with PYTHIA)
    self.plot_final_result(jetR, obs_label, obs_setting, grooming_setting)


  #---------------------------------------------------------------
  # This function is called once after all subconfigurations have been looped over, for each R
  #---------------------------------------------------------------
  def plot_all_results(self, jetR):

    print('Plotting overlay of all results...')
    
    for i_config, overlay_list in enumerate(self.plot_overlay_list):
    
      if len(overlay_list) > 1:
      
        self.plot_final_result_overlay(i_config, jetR, overlay_list)


  #----------------------------------------------------------------------
  # This function is called once after all subconfigurations and jetR have been looped over
  #----------------------------------------------------------------------
  def plot_performance(self):
    print('Plotting performance plots...')
    print("ERROR: Implement plot_performance() in python script")
    return


  #----------------------------------------------------------------------
  def plot_final_result(self, jetR, obs_label, obs_setting, grooming_setting):
    print('Plot final results for {}: R = {}, {}'.format(self.observable, jetR, obs_label))

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    if self.do_theory and float(obs_label.split('_')[0]) in self.theory_alpha and \
       ( (self.use_old and not grooming_setting) or not self.use_old ):
      # Compare parton-level theory to parton-level event generators
      print("Plotting parton-level theory comparisons for", obs_label)
      self.plot_parton_comp(jetR, obs_label, obs_setting, grooming_setting)

    if self.do_theory_F_np and float(obs_label.split('_')[0]) in self.theory_alpha and \
       ( (self.use_old and not grooming_setting) or not self.use_old ):
      # Compare parton-level theory to parton-level event generators
      print("Plotting F_NP-convolved theory comparisons for", obs_label)
      self.plot_Fnp_comp(jetR, obs_label, obs_setting, grooming_setting)

    # Loop through pt slices, and plot final result for each 1D theta_g distribution
    for i in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[i]
      max_pt_truth = self.pt_bins_reported[i+1]
      maxbin = self.obs_max_bins(obs_label)[i]

      self.plot_observable(jetR, obs_label, obs_setting, grooming_setting,
                           min_pt_truth, max_pt_truth, maxbin, plot_MC=True)

      if self.do_theory and float(obs_label.split('_')[0]) in self.theory_alpha and \
         ( (self.use_old and not grooming_setting) or not self.use_old ):
        self.plot_observable(jetR, obs_label, obs_setting, grooming_setting,
                             min_pt_truth, max_pt_truth, maxbin, plot_MC=False, plot_theory=True)
        if self.do_theory_F_np:
          self.plot_observable(jetR, obs_label, obs_setting, grooming_setting,
                               min_pt_truth, max_pt_truth, maxbin, plot_MC=False,
                               plot_theory=True, plot_theory_Fnp=True)
        self.plot_theory_ratios(jetR, obs_label, obs_setting, grooming_setting,
                                min_pt_truth, max_pt_truth, maxbin)
        self.plot_theory_response(jetR, obs_label, obs_setting, grooming_setting,
                                  min_pt_truth, max_pt_truth, maxbin)

      if min_pt_truth == 40 and (jetR == 0.2 or jetR == 0.4):
        # Only want to compare to girth with \alpha=1
        if obs_label == '1':
          self.plot_obs_comp(jetR, obs_label, obs_setting, grooming_setting,
                             min_pt_truth, max_pt_truth, maxbin)

  #----------------------------------------------------------------------
  def plot_observable(self, jetR, obs_label, obs_setting, grooming_setting,
                      min_pt_truth, max_pt_truth, maxbin, plot_MC=False,
                      plot_theory=False, plot_theory_Fnp=False):

    self.set_logy = False
    #if grooming_setting:
    #  self.set_logy = True

    # For theory plots, whether or not to show original parton-level predictions
    show_parton_theory = True
    show_everything_else = True  # set 'False' to show *only* parton-level theory
    rebin_folded = False  # combine every 2 bins to reduce statistical fluctuations

    # For theory plots, whether or not to show the NP / P region
    show_np_region = True
    lambda_np_cutoff = None   # initialize variable for setting later

    # For theory plots, whether or not to show the folded uncertainty bands
    show_folded_uncertainty = True

    name = 'cResult_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt_truth, max_pt_truth)
    c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()
    
    c.cd()
    myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    myPad.SetLeftMargin(0.2)
    myPad.SetTopMargin(0.07)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.13)
    if self.set_logy:
      myPad.SetLogy()
    myPad.Draw()
    myPad.cd()
    
    xtitle = getattr(self, 'xtitle')
    ytitle = getattr(self, 'ytitle')
    color = 1   # black for data
    
    # Get histograms
    name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label,
                                              min_pt_truth, max_pt_truth)
    if grooming_setting:
      fraction_tagged = getattr(self, 'tagging_fraction_R{}_{}_{}-{}'.format(
        jetR, obs_label, min_pt_truth, max_pt_truth))
      #fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
      # maxbin+1 in grooming case to account for extra tagging bin
    if grooming_setting and maxbin:
      h = self.truncate_hist(getattr(self, name), maxbin+1, name+'_trunc')
    else:
      h = self.truncate_hist(getattr(self, name), maxbin, name+'_trunc')
    h.SetMarkerSize(1.5)
    h.SetMarkerStyle(20)
    h.SetMarkerColor(color)
    h.SetLineStyle(1)
    h.SetLineWidth(2)
    h.SetLineColor(color)
    
    h_sys = getattr(self, 'hResult_{}_systotal_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth))
    h_sys.SetLineColor(0)
    h_sys.SetFillColor(color)
    h_sys.SetFillColorAlpha(color, 0.3)
    h_sys.SetFillStyle(1001)
    h_sys.SetLineWidth(0)

    n_obs_bins_truth = self.n_bins_truth(obs_label)
    truth_bin_array = self.truth_bin_array(obs_label)
    if maxbin:
      truth_bin_array = truth_bin_array[0:maxbin+1]
      n_obs_bins_truth = len(truth_bin_array)-1
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_obs_bins_truth, truth_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle(xtitle)
    myBlankHisto.GetXaxis().SetTitleOffset(1.02)
    myBlankHisto.GetXaxis().SetTitleSize(0.055)
    myBlankHisto.SetYTitle(ytitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.1)
    myBlankHisto.GetYaxis().SetTitleSize(0.055)
    ymin = 1e-4 if self.set_logy else 0
    myBlankHisto.SetMinimum(ymin)
    if not plot_theory or not show_parton_theory:
      maxval = max(2.3*h.GetBinContent(int(0.4*h.GetNbinsX())), 1.7*h.GetMaximum())
      myBlankHisto.SetMaximum(maxval)
      myBlankHisto.Draw("E")

    if plot_theory:
      label = "R%s_%s" % (str(jetR).replace('.', ''), str(obs_label).replace('.', ''))

      # Also plot the theory predictions pre-folding for comparison...
      name_cent = "theory_cent_%s_%s_parton_PtBin%i-%i" % \
                  (self.observable, label, min_pt_truth, max_pt_truth)
      name_min = "theory_min_%s_%s_parton_PtBin%i-%i" % \
                 (self.observable, label, min_pt_truth, max_pt_truth)
      name_max = "theory_max_%s_%s_parton_PtBin%i-%i" % \
                 (self.observable, label, min_pt_truth, max_pt_truth)

      hcent_p = getattr(self, name_cent)
      hmin_p = getattr(self, name_min)
      hmax_p = getattr(self, name_max)

      if show_parton_theory:
        maxval = 1.7*max([hmax_p.GetBinContent(i) for
                          i in range(1, round(len(self.theory_obs_bins)/2))])
        myBlankHisto.SetMaximum(maxval)
        myBlankHisto.Draw("E")

      n = hcent_p.GetNbinsX()
      x = array('d', [hcent_p.GetXaxis().GetBinCenter(i) for i in range(1, hcent_p.GetNbinsX()+1)])
      y = array('d', [hcent_p.GetBinContent(i) for i in range(1, hcent_p.GetNbinsX()+1)])
      xerrup = array('d', [(x[i+1] - x[i]) / 2. for i in range(n-1)] + [0])
      xerrdn = array('d', [0] + [(x[i+1] - x[i]) / 2. for i in range(n-1)])
      yerrup = array('d', [hmax_p.GetBinContent(i)-y[i-1] for i in range(1, hmax_p.GetNbinsX()+1)])
      yerrdn = array('d', [y[i-1]-hmin_p.GetBinContent(i) for i in range(1, hmin_p.GetNbinsX()+1)])

      color = self.ColorArray[4]

      if show_np_region:
        # P vs NP cutoff point: lambda_alpha ~ Lambda / (pT * R) -- use avg value of pT for the bin.
        # Formula assumes that jet pT xsec falls like pT^(-5.5)
        formula_pt = (4.5/3.5)*(min_pt_truth**-3.5 - max_pt_truth**-3.5) / \
                     (min_pt_truth**-4.5 - max_pt_truth**-4.5)
        lambda_np_cutoff = min(self.theory_obs_bins,
                               key=lambda x:abs(x - self.Lambda / (formula_pt * jetR)))
        if lambda_np_cutoff > 0.8:
          lambda_np_cutoff = 0.8
        index_np_cutoff = list(self.theory_obs_bins).index(lambda_np_cutoff)

        if show_parton_theory:
          #+1 to include lambda in NP
          h_parton_theory_np = ROOT.TGraphAsymmErrors(
            index_np_cutoff+1, x[:index_np_cutoff+1], y[:index_np_cutoff+1], xerrdn[:index_np_cutoff+1],
            xerrup[:index_np_cutoff+1], yerrdn[:index_np_cutoff+1], yerrup[:index_np_cutoff+1])
          h_parton_theory_p = ROOT.TGraphAsymmErrors(
            n-index_np_cutoff, x[index_np_cutoff:], y[index_np_cutoff:], xerrdn[index_np_cutoff:],
            xerrup[index_np_cutoff:], yerrdn[index_np_cutoff:], yerrup[index_np_cutoff:])

          h_parton_theory_np.SetFillColorAlpha(color, 0.5)
          h_parton_theory_np.SetFillStyle(3002)
          h_parton_theory_np.SetLineStyle(5)
          h_parton_theory_np.SetLineColor(color)
          h_parton_theory_np.SetLineWidth(3)
          h_parton_theory_np.Draw('L 3 same')

          h_parton_theory_p.SetFillColorAlpha(color, 0.25)
          h_parton_theory_p.SetLineColor(color)
          h_parton_theory_p.SetLineWidth(3)
          h_parton_theory_p.Draw('L 3 same')

      elif show_parton_theory:   # don't show NP / P region on the plot
        h_parton_theory = ROOT.TGraphAsymmErrors(n, x, y, xerrdn, xerrup, yerrdn, yerrup)
        h_parton_theory.SetFillColorAlpha(color, 0.25)
        h_parton_theory.SetLineColor(color)
        h_parton_theory.SetLineWidth(3)
        h_parton_theory.Draw('L 3 same')

        # Dotted lines for error bars
        #hmin_p.SetLineColor(color)
        #hmin_p.SetLineWidth(1)
        #hmin_p.SetLineStyle(2)
        #hmin_p.Draw('L hist same')

        #hmax_p.SetLineColor(color)
        #hmax_p.SetLineWidth(1)
        #hmax_p.SetLineStyle(2)
        #hmax_p.Draw('L hist same')

      hcent_list = []   # List of items for legend
      hasym_list = []   # List to prevent hists from being garbage collected by python
      expected_list = []
      for ri in range(len(self.theory_response_files)):
        hcent = None; hmin = None; hmax = None;

        if plot_theory_Fnp:  # Have to have an extra loop for Omega
          if ri != 0:
            continue

          for oi, Omega in enumerate(self.Omega_list):
            color = self.ColorArray[5 + oi]

            # For testing if passed in flat distribution
            if self.exp_test:
              expected = ROOT.TH1D("expected Fnp", "expected Fnp", len(self.theory_obs_bins) - 1, 
                                   array('d', self.theory_obs_bins))
              Omega_a = Omega / (float(obs_label.split('_')[0]) - 1)
              for i in range(1, len(self.theory_obs_bins)):
                ob = self.theory_obs_bins_center[i-1]
                x1 = 2 / Omega_a * jetR * 65 * ob
                val = Omega_a * (1 - np.exp(-1*x1) * (1 + x1))
                expected.SetBinContent(i, val)
                expected.SetBinError(i, 0)
              expected.Scale(1/expected.Integral("width"))
              expected.SetLineColor(color)
              expected.SetLineStyle(2)
              expected.SetLineWidth(3)
              expected.Draw('L same')
              expected_list.append(expected)

            # Load 2D convolved/folded histograms
            name_cent = "theory_cent_%s_%s_ch_Fnp_Omega%s_%i" % \
                        (self.observable, label, str(Omega), ri)
            name_min = "theory_min_%s_%s_ch_Fnp_Omega%s_%i" % \
                       (self.observable, label, str(Omega), ri)
            name_max = "theory_max_%s_%s_ch_Fnp_Omega%s_%i" % \
                       (self.observable, label, str(Omega), ri)

            hcent_2D = getattr(self, name_cent)
            hmin_2D = getattr(self, name_min)
            hmax_2D = getattr(self, name_max)

            # Project onto correct pT bin with new name
            name_cent = "theory_cent_%s_%s_ch_Fnp_Omega%s_PtBin%i-%i_%i" % \
                        (self.observable, label, str(Omega), min_pt_truth, max_pt_truth, ri)
            name_min = "theory_min_%s_%s_ch_Fnp_Omega%s_PtBin%i-%i_%i" % \
                       (self.observable, label, str(Omega), min_pt_truth, max_pt_truth, ri)
            name_max = "theory_max_%s_%s_ch_Fnp_Omega%s_PtBin%i-%i_%i" % \
                       (self.observable, label, str(Omega), min_pt_truth, max_pt_truth, ri)
            
            # Have to do +1 for ROOT 1-indexing
            min_pt_i = self.theory_pt_bins.index(min_pt_truth)+1
            max_pt_i = self.theory_pt_bins.index(max_pt_truth)+1
            hcent = hcent_2D.ProjectionY(name_cent, min_pt_i, max_pt_i)
            hmin = hmin_2D.ProjectionY(name_min, min_pt_i, max_pt_i)
            hmax = hmax_2D.ProjectionY(name_max, min_pt_i, max_pt_i)
            hmin.Scale(1/hcent.Integral(), "width")
            hmax.Scale(1/hcent.Integral(), "width")
            hcent.Scale(1/hcent.Integral(), "width")

            if rebin_folded:
              hcent = hcent.Rebin(); hcent.Scale(1/2)
              hmin = hmin.Rebin(); hmin.Scale(1/2)
              hmax = hmax.Rebin(); hmax.Scale(1/2)

            # Write folded theory result to ROOT file
            output_dir = getattr(self, 'output_dir_final_results')
            final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
            fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
            if plot_theory:
              hcent.Write()
              hmin.Write()
              hmax.Write()
            fFinalResults.Close()

            if show_folded_uncertainty:
              n = hcent.GetNbinsX()
              x = array('d', [hcent.GetXaxis().GetBinCenter(i) for
                              i in range(1, hcent.GetNbinsX()+1)])
              y = array('d', [hcent.GetBinContent(i) for
                              i in range(1, hcent.GetNbinsX()+1)])
              xerrup = array('d', [(x[i+1] - x[i]) / 2. for i in range(n-1)] + [0])
              xerrdn = array('d', [0] + [(x[i+1] - x[i]) / 2. for i in range(n-1)])
              yerrup = array('d', [hmax.GetBinContent(i)-y[i-1] for
                                   i in range(1, hmax.GetNbinsX()+1)])
              yerrdn = array('d', [y[i-1]-hmin.GetBinContent(i) for
                                   i in range(1, hmin.GetNbinsX()+1)])
              h_ch_theory = ROOT.TGraphAsymmErrors(n, x, y, xerrdn, xerrup, yerrdn, yerrup)

              h_ch_theory.SetFillColorAlpha(color, 0.25)
              h_ch_theory.SetLineColor(color)
              h_ch_theory.SetLineWidth(3)
              if show_everything_else:
                h_ch_theory.Draw('L 3 same')
              hasym_list.append(h_ch_theory)

            if show_everything_else and not show_folded_uncertainty:
              hcent.SetLineColor(color)
              hcent.SetLineWidth(3)
              hcent.Draw('L hist same')
              hcent_list.append(hcent)  # Save for legend

        else:
          # Get and plot the folded & MPI-scaled theory predictions
          name_cent = "theory_cent_%s_%s_ch_PtBin%i-%i_%i" % \
                      (self.observable, label, min_pt_truth, max_pt_truth, ri)
          name_min = "theory_min_%s_%s_ch_PtBin%i-%i_%i" % \
                     (self.observable, label, min_pt_truth, max_pt_truth, ri)
          name_max = "theory_max_%s_%s_ch_PtBin%i-%i_%i" % \
                     (self.observable, label, min_pt_truth, max_pt_truth, ri)

          hcent = getattr(self, name_cent)
          hmin = getattr(self, name_min)
          hmax = getattr(self, name_max)
          color = self.ColorArray[5 + ri]

          if rebin_folded:
            hcent = hcent.Rebin(); hcent.Scale(1/2)
            hmin = hmin.Rebin(); hmin.Scale(1/2)
            hmax = hmax.Rebin(); hmax.Scale(1/2)

          # Write folded theory result to ROOT file
          output_dir = getattr(self, 'output_dir_final_results')
          final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
          fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
          if plot_theory:
            hcent.Write()
            hmin.Write()
            hmax.Write()
          fFinalResults.Close()

          if show_folded_uncertainty:
            n = hcent.GetNbinsX()
            x = array('d', [hcent.GetXaxis().GetBinCenter(i) for i in range(1, hcent.GetNbinsX()+1)])
            y = array('d', [hcent.GetBinContent(i) for i in range(1, hcent.GetNbinsX()+1)])
            xerrup = array('d', [(x[i+1] - x[i]) / 2. for i in range(n-1)] + [0])
            xerrdn = array('d', [0] + [(x[i+1] - x[i]) / 2. for i in range(n-1)])
            yerrup = array('d', [hmax.GetBinContent(i)-y[i-1] for i in range(1, hmax.GetNbinsX()+1)])
            yerrdn = array('d', [y[i-1]-hmin.GetBinContent(i) for i in range(1, hmin.GetNbinsX()+1)])
            h_ch_theory = ROOT.TGraphAsymmErrors(n, x, y, xerrdn, xerrup, yerrdn, yerrup)

            h_ch_theory.SetFillColorAlpha(color, 0.25)
            h_ch_theory.SetLineColor(color)
            h_ch_theory.SetLineWidth(3)
            if show_everything_else:
              h_ch_theory.Draw('L 3 same')
            hasym_list.append(h_ch_theory)

          if show_everything_else and not show_folded_uncertainty:
            hcent.SetLineColor(color)
            hcent.SetLineWidth(3)
            hcent.Draw('L hist same')
            hcent_list.append(hcent)  # Save for legend

    plot_pythia = False; plot_herwig = False;
    if plot_MC:

      hPythia = None; fraction_tagged_pythia = None;
      if grooming_setting:
        hPythia, fraction_tagged_pythia = self.MC_prediction(
          jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin+1, 'Pythia')
      else:
        hPythia, fraction_tagged_pythia = self.MC_prediction(
          jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin, 'Pythia')

      hHerwig = None; fraction_tagged_herwig = None;
      if self.do_theory and not self.use_prev_prelim and \
         'fastsim_generator1' in self.systematics_list:
        # Load Herwig comparison as well
        if grooming_setting:
          hHerwig, fraction_tagged_herwig = self.MC_prediction(
            jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin+1, 'Herwig')
        else:
          hHerwig, fraction_tagged_herwig = self.MC_prediction(
            jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin, 'Herwig')

      if hPythia:
        hPythia.SetFillStyle(0)
        hPythia.SetMarkerSize(1.5)
        hPythia.SetMarkerStyle(21)
        hPythia.SetMarkerColor(600-6)
        hPythia.SetLineColor(600-6)
        hPythia.SetLineWidth(1)
        hPythia.Draw('E2 same')
        plot_pythia = True
      else:
        print('No PYTHIA prediction for %s %s' % (self.observable, obs_label))

      if hHerwig:
        hHerwig.SetFillStyle(0)
        hHerwig.SetMarkerSize(1.5)
        hHerwig.SetMarkerStyle(34)
        hHerwig.SetMarkerColor(600+3)
        hHerwig.SetLineColor(600+3)
        hHerwig.SetLineWidth(1)
        hHerwig.Draw('E2 same')
        plot_herwig = True

    # Vertical line for perturbative / non-perturbative region
    if plot_theory and show_everything_else:
      # Can provide different max values between legend and bin text
      if lambda_np_cutoff < 0.58 * truth_bin_array[-1]:
        m = 0.53 * maxval
      else:
        m = 0.53 * maxval
      line = ROOT.TLine(lambda_np_cutoff + xerrup[index_np_cutoff], 0,
                        lambda_np_cutoff + xerrup[index_np_cutoff], m)
      line.SetLineColor(self.ColorArray[4])
      line.SetLineStyle(2)
      line.SetLineWidth(2)
      line.Draw()

    if show_everything_else:
      h_sys.DrawCopy('E2 same')
      h.DrawCopy('PE X0 same')
  
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text_xval = 0.6 if grooming_setting else 0.63
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(text_xval, 0.87, text)
    
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(text_xval, 0.8, text)

    text = "anti-#it{k}_{T} jets,   #it{R} = %s" % str(jetR)
    text_latex.DrawLatex(text_xval, 0.73, text)

    if plot_theory and show_parton_theory and not show_everything_else:
      text = str(min_pt_truth) + ' < #it{p}_{T,jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    else:
      text = str(min_pt_truth) + ' < #it{p}_{T,jet}^{ch} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.DrawLatex(text_xval, 0.66, text)

    text = '| #it{#eta}_{jet}| < %s' % str(0.9 - jetR)
    subobs_label = self.utils.formatted_subobs_label(self.observable)
    if subobs_label:
      text += ',   %s = %s' % (subobs_label, obs_setting)
    delta = 0.07
    text_latex.DrawLatex(text_xval, 0.66-delta, text)
    
    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting)#.replace("#beta}", "#beta}_{SD}")
      text_latex.DrawLatex(text_xval, 0.66-2*delta, text)
      
      text_latex.SetTextSize(0.04)
      text = '#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged
      if plot_pythia:
        text += (', #it{f}_{tagged}^{pythia} = %3.3f' % fraction_tagged_pythia)
      if plot_herwig:
        text += (', #it{f}_{tagged}^{herwig} = %3.3f' % fraction_tagged_herwig)
      text_latex.DrawLatex(text_xval, 0.66-3*delta, text)

    if plot_theory and show_parton_theory and not show_everything_else:
      myLegend = ROOT.TLegend(0.21, 0.79, 0.45, 0.91)
    else:
      miny = 0.57
      if plot_pythia:
        if plot_herwig:
          miny = 0.72  #TODO
        else:
          miny = 0.72
      myLegend = ROOT.TLegend(0.23, miny, text_xval-0.02, 0.92)
    self.utils.setup_legend(myLegend, 0.035)
    if show_everything_else:
      myLegend.AddEntry(h, 'ALICE pp', 'pe')
      myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    if plot_pythia:
      myLegend.AddEntry(hPythia, 'PYTHIA8 Monash2013', 'pe')
    if plot_herwig:
      myLegend.AddEntry(hHerwig, 'Herwig7 Default', 'pe')
    if plot_theory:
      if show_parton_theory:
        if show_np_region:
          myLegend.AddEntry(h_parton_theory_np, 'NLL\' (non-perturbative)', 'lf')
          myLegend.AddEntry(h_parton_theory_p, 'NLL\' (perturbative)', 'lf')
        else:
          myLegend.AddEntry(hcent_p, 'NLL\' (Parton)', 'lf')
      if show_everything_else:
        if show_folded_uncertainty:
          if plot_theory_Fnp:
            for oi, Omega in enumerate(self.Omega_list):
              myLegend.AddEntry(hasym_list[oi],
                                'NLL\' #otimes F_{NP}^{#Omega=%s} #otimes %s' %
                                (str(Omega), self.theory_response_labels[0]), 'lf')
          else:
            for ri, lab in enumerate(self.theory_response_labels):
              myLegend.AddEntry(hasym_list[ri], 'NLL\' #otimes '+lab, 'lf')
        else:
          if plot_theory_Fnp:
            for oi, Omega in enumerate(self.Omega_list):
              myLegend.AddEntry(hcent_list[oi],
                                'NLL\' #otimes F_{NP}^{#Omega=%s} #otimes %s' %
                                (str(Omega), self.theory_response_labels[0]), 'lf')
          else:
            for ri, lab in enumerate(self.theory_response_labels):
              myLegend.AddEntry(hcent_list[ri], 'NLL\' #otimes '+lab, 'lf')
        myLegend.AddEntry(line, '#it{#lambda}_{#it{#alpha}}^{NP region} #leq' + \
                          '#Lambda / (#it{p}_{T,jet}^{ch} #it{R})', 'lf')
    myLegend.Draw()

    name = 'hUnfolded_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), obs_label,
                                             int(min_pt_truth), int(max_pt_truth), self.file_format)
    if plot_MC:
      name = 'hUnfolded_R{}_{}_{}-{}_MC{}'.format(
        self.utils.remove_periods(jetR), obs_label, int(min_pt_truth),
        int(max_pt_truth), self.file_format)

    if plot_theory:
      if plot_theory_Fnp:
        name = 'hUnfolded_R{}_{}_{}-{}_FnpTheory{}'.format(
          self.utils.remove_periods(jetR), obs_label, int(min_pt_truth),
          int(max_pt_truth), self.file_format)
      else:
        name = 'hUnfolded_R{}_{}_{}-{}_Theory{}'.format(
          self.utils.remove_periods(jetR), obs_label, int(min_pt_truth),
          int(max_pt_truth), self.file_format)

    output_dir = getattr(self, 'output_dir_final_results')
    output_dir_single = output_dir + '/single_results'
    if not os.path.exists(output_dir_single):
      os.mkdir(output_dir_single)
    outputFilename = os.path.join(output_dir_single, name)
    c.SaveAs(outputFilename)
    c.Close()

    # Write result to ROOT file
    final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
    fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
    h.Write()
    h_sys.Write()
    if plot_pythia:
      hPythia.Write()
    if plot_herwig:
      hHerwig.Write()
    fFinalResults.Close()

  
  #----------------------------------------------------------------------
  def plot_theory_ratios(self, jetR, obs_label, obs_setting, grooming_setting,
                         min_pt_truth, max_pt_truth, maxbin):
    
    # Get histograms
    name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label,
                                              min_pt_truth, max_pt_truth)
    if grooming_setting and maxbin:
      h = self.truncate_hist(getattr(self, name), maxbin+1, name+'_trunc')
    else:
      h = self.truncate_hist(getattr(self, name), maxbin, name+'_trunc')
    
    n_obs_bins_truth = self.n_bins_truth(obs_label)
    truth_bin_array = self.truth_bin_array(obs_label)
    if maxbin:
      truth_bin_array = truth_bin_array[0:maxbin+1]
      n_obs_bins_truth = len(truth_bin_array)-1

    label = "R%s_%s" % (str(jetR).replace('.', ''), str(obs_label).replace('.', ''))

    for ri in range(len(self.theory_response_files)):
      
      name = 'cResult_R{}_{}_{}-{}_{}'.format(jetR, obs_label, min_pt_truth, max_pt_truth, ri)
      c = ROOT.TCanvas(name, name, 600, 450)
      c.Draw()

      c.cd()
      myPad = ROOT.TPad('myPad_%i' % ri, 'The pad %i' % ri, 0, 0, 1, 1)
      myPad.SetLeftMargin(0.2)
      myPad.SetTopMargin(0.07)
      myPad.SetRightMargin(0.04)
      myPad.SetBottomMargin(0.13)
      myPad.Draw()
      myPad.cd()

      xtitle = getattr(self, 'xtitle')
      ytitle = getattr(self, 'ytitle')
      color = 600-6

      h_ratio_h = getattr(self, 'h_cent_ratio_h_%s_PtBin%i-%i_%i' % \
                          (label, min_pt_truth, max_pt_truth, ri))
      h_ratio_h.SetMarkerSize(0.5)
      h_ratio_h.SetMarkerStyle(20)
      h_ratio_h.SetMarkerColor(self.ColorArray[3])
      h_ratio_h.SetLineStyle(1)
      h_ratio_h.SetLineWidth(2)
      h_ratio_h.SetLineColor(self.ColorArray[3])

      h_ratio_ch = getattr(self, 'h_cent_ratio_ch_%s_PtBin%i-%i_%i' % \
                           (label, min_pt_truth, max_pt_truth, ri))
      h_ratio_ch.SetMarkerSize(0.5)
      h_ratio_ch.SetMarkerStyle(20)
      h_ratio_ch.SetMarkerColor(self.ColorArray[5])
      h_ratio_ch.SetLineStyle(1)
      h_ratio_ch.SetLineWidth(2)
      h_ratio_ch.SetLineColor(self.ColorArray[5])

      h_ratio_mpi = getattr(self, "h_mpi_ratio_%s_PtBin%i-%i_%i" % \
                            (label, min_pt_truth, max_pt_truth, ri))
      h_ratio_mpi.SetMarkerSize(0.5)
      h_ratio_mpi.SetMarkerStyle(20)
      h_ratio_mpi.SetMarkerColor(self.ColorArray[4])
      #h_ratio_mpi.SetLineStyle(1)
      #h_ratio_mpi.SetLineWidth(2)
      #h_ratio_mpi.SetLineColor(self.ColorArray[4])

      h_ratio_total = h_ratio_h.Clone()
      h_ratio_total.SetNameTitle('h_theory_total_ratio_scaling_%i' % ri,
                                 'h_theory_total_ratio_scaling_%i' % ri)
      h_ratio_total.Multiply(h_ratio_ch)
      h_ratio_total.Multiply(h_ratio_mpi)
      h_ratio_total.SetMarkerSize(1)
      h_ratio_total.SetMarkerStyle(20)
      h_ratio_total.SetMarkerColor(1)
      h_ratio_total.SetLineStyle(1)
      h_ratio_total.SetLineWidth(3)
      h_ratio_total.SetLineColor(1)

      myBlankHisto = ROOT.TH1F('myBlankHisto_%i' % ri,'Blank Histogram %i' % ri,
                               n_obs_bins_truth, truth_bin_array)
      myBlankHisto.SetNdivisions(505)
      myBlankHisto.SetXTitle(xtitle)
      myBlankHisto.GetYaxis().SetTitleOffset(1.5)
      myBlankHisto.SetYTitle(ytitle)
      val = 2.5*max([h_ratio_total.GetBinContent(i) for i in range(5, h_ratio_total.GetNbinsX())])
      myBlankHisto.SetMaximum(val)
      myBlankHisto.SetMinimum(0)
      myBlankHisto.Draw("E")

      h_ratio_h.Draw('hist E same')
      h_ratio_ch.Draw('hist E same')
      h_ratio_mpi.Draw('P E same')
      h_ratio_total.Draw('hist same')

      # Draw a horizontal line at 1
      line = ROOT.TLine(0, 1, truth_bin_array[-1], 1)
      line.SetLineColor(920+2)
      line.SetLineStyle(2)
      line.SetLineWidth(2)
      line.Draw()

      text_latex = ROOT.TLatex()
      text_latex.SetNDC()
      text = 'ALICE {}'.format(self.figure_approval_status)
      text_latex.DrawLatex(0.6, 0.87, text)

      text = 'pp #sqrt{#it{s}} = 5.02 TeV'
      text_latex.SetTextSize(0.045)
      text_latex.DrawLatex(0.6, 0.8, text)

      text = str(min_pt_truth) + ' < #it{p}_{T,jet}^{ch} < ' + str(max_pt_truth) + ' GeV/#it{c}'
      text_latex.DrawLatex(0.6, 0.73, text)

      text = '#it{R} = ' + str(jetR) + ',   | #it{#eta}_{jet}| < %s' % str(0.9 - jetR)
      text_latex.DrawLatex(0.6, 0.66, text)

      subobs_label = self.utils.formatted_subobs_label(self.observable)
      delta = 0.
      if subobs_label:
        text = '%s = %s' % (subobs_label, obs_setting)
        text += ',   %s' % self.theory_response_labels[ri]
        text_latex.DrawLatex(0.6, 0.59, text)
        delta = 0.07

      myLegend = ROOT.TLegend(0.27, 0.7, 0.55, 0.9)
      self.utils.setup_legend(myLegend, 0.035)
      myLegend.AddEntry(h_ratio_total, 'Total ratio', 'pe')
      myLegend.AddEntry(h_ratio_h, 'p to h ratio', 'pe')
      myLegend.AddEntry(h_ratio_ch, 'h to ch ratio', 'pe')
      myLegend.AddEntry(h_ratio_mpi, 'MPI scaling ratio', 'pe')
      myLegend.Draw()

      name = 'hTheoryRatio_R{}_{}_{}-{}_{}{}'.format(
        self.utils.remove_periods(jetR), obs_label, int(min_pt_truth), int(max_pt_truth),
        self.theory_response_labels[ri], self.file_format)

      output_dir = self.output_dir_theory
      if not os.path.exists(output_dir):
        os.mkdir(output_dir)
      outputFilename = os.path.join(output_dir, name)
      c.SaveAs(outputFilename)
      c.Close()


  #----------------------------------------------------------------------
  def plot_theory_response(self, jetR, obs_label, obs_setting, grooming_setting,
                           min_pt_truth, max_pt_truth, maxbin):

    label = "R%s_%s" % (str(jetR).replace('.', ''), str(obs_label).replace('.', ''))
    outf = ROOT.TFile(os.path.join(self.output_dir_theory, 'fTheoryResponseProj.root'), 'UPDATE')

    for ri in range(len(self.theory_response_files)):
      # Get histograms
      thn_ch = getattr(self, "hResponse_theory_ch_%s_%i" % (label, ri))
      thn_h = getattr(self, "hResponse_theory_h_%s_%i" % (label, ri))

      # Make projections in pT bins at (charged-)/hadron level
      thn_ch.GetAxis(0).SetRangeUser(int(min_pt_truth), int(max_pt_truth))
      thn_h.GetAxis(0).SetRangeUser(int(min_pt_truth), int(max_pt_truth))

      #print(thn_ch.GetBinContent(array('i', [3, 3, 20, 20])))

      hTheoryProjection_ch = thn_ch.Projection(2, 3)
      hTheoryProjection_h = thn_h.Projection(2, 3)

      #print(hTheoryProjection_ch.GetBinContent(20, 20))
      #exit()

      name_ch = "hResponse_theory_ch_%s_PtBin%i-%i_%s" % \
                (label, min_pt_truth, max_pt_truth, self.theory_response_labels[ri])
      name_h = "hResponse_theory_h_%s_PtBin%i-%i_%s" % \
               (label, min_pt_truth, max_pt_truth, self.theory_response_labels[ri])

      hTheoryProjection_ch.SetNameTitle(name_ch, name_ch)
      hTheoryProjection_h.SetNameTitle(name_h, name_h)

      # Save the histograms
      output_dir = self.output_dir_theory
      if not os.path.exists(output_dir):
        os.mkdir(output_dir)

      text_h = str(min_pt_truth) + ' < #it{p}_{T, h jet} < ' + str(max_pt_truth)
      text_ch = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth)
      self.utils.plot_hist(hTheoryProjection_h, os.path.join(self.output_dir_theory, name_h+'.pdf'),
                           'colz', False, True, text_h)
      self.utils.plot_hist(hTheoryProjection_ch, os.path.join(self.output_dir_theory, name_ch+'.pdf'),
                           'colz', False, True, text_ch)

      hTheoryProjection_h.Write()
      hTheoryProjection_ch.Write()

      # Reset axes zoom in case histograms used later
      thn_ch.GetAxis(0).UnZoom()
      thn_h.GetAxis(0).UnZoom()

    outf.Close()


  #----------------------------------------------------------------------
  def plot_parton_comp(self, jetR, obs_label, obs_setting, grooming_setting):

    label = "R%s_%s" % (str(jetR).replace('.', ''), str(obs_label).replace('.', ''))

    # Directory to save the histograms
    output_dir = os.path.join(self.output_dir_theory, "parton_comp")
    if not os.path.exists(output_dir):
      os.mkdir(output_dir)

    hcent_p = getattr(self, "theory_cent_%s_%s_parton" % (self.observable, label))
    hmin_p = getattr(self, "theory_min_%s_%s_parton" % (self.observable, label))
    hmax_p = getattr(self, "theory_max_%s_%s_parton" % (self.observable, label))

    n_obs_bins = len(self.theory_obs_bins) - 1
    obs_edges = self.theory_obs_bins

    # Make projections in pT bins at parton level
    for i, min_pt in list(enumerate(self.theory_pt_bins))[:-1]:
      max_pt = self.theory_pt_bins[i+1]

      # Get the theory prediction
      hcent_p.GetXaxis().SetRangeUser(int(min_pt), int(max_pt))
      hmin_p.GetXaxis().SetRangeUser(int(min_pt), int(max_pt))
      hmax_p.GetXaxis().SetRangeUser(int(min_pt), int(max_pt))

      hcent_proj_p = hcent_p.ProjectionY()
      hmin_proj_p = hmin_p.ProjectionY()
      hmax_proj_p = hmax_p.ProjectionY()

      # Fix normalization
      hmin_proj_p.Scale(1/hcent_proj_p.Integral(), "width")
      hmax_proj_p.Scale(1/hcent_proj_p.Integral(), "width")
      hcent_proj_p.Scale(1/hcent_proj_p.Integral(), "width")

      # Initialize canvas & pad for plotting
      name = 'cTheoryComp_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt, max_pt)
      c = ROOT.TCanvas(name, name, 600, 450)
      c.Draw()
      c.cd()

      name = 'pTheoryComp_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt, max_pt)
      myPad = ROOT.TPad(name, name, 0, 0, 1, 1)
      myPad.SetLeftMargin(0.2)
      myPad.SetTopMargin(0.07)
      myPad.SetRightMargin(0.04)
      myPad.SetBottomMargin(0.13)
      #if grooming_setting:
      #  myPad.SetLogy()
      myPad.Draw()
      myPad.cd()

      # Find the max bin: the last bin where the angularity is 0
      maxbin = hcent_proj_p.GetNbinsX()-1   # initialize
      while hmax_proj_p.GetBinContent(maxbin) < 1e-3 and maxbin > 1:
        maxbin -= 1

      # Use blank histogram to initialize this range
      bin_array = array('d', obs_edges[0:maxbin+1])
      name = 'hTheoryComp_R{}_{}_{}-{}_Blank'.format(jetR, obs_label, min_pt, max_pt)
      myBlankHisto = ROOT.TH1F(name, name, maxbin, bin_array)
      myBlankHisto.SetNdivisions(505)
      myBlankHisto.SetXTitle(self.xtitle)
      myBlankHisto.GetXaxis().SetTitleOffset(1.02)
      myBlankHisto.GetXaxis().SetTitleSize(0.055)
      myBlankHisto.SetYTitle(self.ytitle)
      myBlankHisto.GetYaxis().SetTitleOffset(1.1)
      myBlankHisto.GetYaxis().SetTitleSize(0.055)
      myBlankHisto.SetMinimum(0.)
      myBlankHisto.SetMaximum(1.7*hmax_proj_p.GetMaximum())
      myBlankHisto.Draw()

      x = array('d', [round(hcent_proj_p.GetXaxis().GetBinCenter(i), 5) for i in range(1, maxbin+1)])
      y = array('d', [hcent_proj_p.GetBinContent(i) for i in range(1, maxbin+1)])
      xerrup = array('d', [(x[i+1] - x[i]) / 2. for i in range(maxbin-1)] + [0])
      xerrdn = array('d', [0] + [(x[i+1] - x[i]) / 2. for i in range(maxbin-1)])
      yerrup = array('d', [hmax_proj_p.GetBinContent(i)-y[i-1] for i in range(1, maxbin+1)])
      yerrdn = array('d', [y[i-1]-hmin_proj_p.GetBinContent(i) for i in range(1, maxbin+1)])
      h_theory = ROOT.TGraphAsymmErrors(maxbin, x, y, xerrdn, xerrup, yerrdn, yerrup)
      color = self.ColorArray[4]
      h_theory.SetFillColorAlpha(color, 0.25)
      h_theory.SetLineColor(color)
      h_theory.SetLineWidth(3)
      h_theory.Draw('L 3 same')

      h_resp_list = []
      for ri in range(len(self.theory_response_files)):
        # Get event generator histogram
        thn = getattr(self, "hResponse_theory_ch_%s_%i" % (label, ri))

        # Get the response matrix prediction
        thn.GetAxis(1).SetRangeUser(int(min_pt), int(max_pt))
        h_response_projection = thn.Projection(3)
        name = "hTheoryComp_p_%s_PtBin%i-%i_%s" % \
                  (label, min_pt, max_pt, self.theory_response_labels[ri])
        h_response_projection.SetNameTitle(name, name)
        h_response_projection.SetDirectory(0)

        # Rescale by integral for correct normalization
        h_response_projection.Scale(1/h_response_projection.Integral(), "width")

        color = self.ColorArray[4+1+ri]
        h_response_projection.SetLineColor(color)
        h_response_projection.SetLineWidth(3)
        h_response_projection.Draw('L hist same')
        h_resp_list.append(h_response_projection)

        # Reset thn range in case used later
        thn.GetAxis(1).UnZoom()

      text_latex = ROOT.TLatex()
      text_latex.SetNDC()
      text_xval = 0.61
      text = 'ALICE {}'.format(self.figure_approval_status)
      text_latex.DrawLatex(text_xval, 0.87, text)

      text = 'pp #sqrt{#it{s}} = 5.02 TeV'
      text_latex.SetTextSize(0.045)
      text_latex.DrawLatex(text_xval, 0.8, text)

      text = "anti-#it{k}_{T} jets,   #it{R} = %s" % str(jetR)
      text_latex.DrawLatex(text_xval, 0.73, text)

      text = str(min_pt) + ' < #it{p}_{T,jet}^{parton} < ' + str(max_pt) + ' GeV/#it{c}'
      text_latex.DrawLatex(text_xval, 0.66, text)

      text = '| #it{#eta}_{jet}| < %s' % str(0.9 - jetR)
      subobs_label = self.utils.formatted_subobs_label(self.observable)
      if subobs_label:
        text += ',   %s = %s' % (subobs_label, obs_setting)
      delta = 0.07
      text_latex.DrawLatex(text_xval, 0.66-delta, text)

      myLegend = ROOT.TLegend(0.27, 0.7, 0.55, 0.9)
      self.utils.setup_legend(myLegend, 0.035)
      myLegend.AddEntry(h_theory, "NLL'", 'lf')
      for rl, l in enumerate(self.theory_response_labels):
        myLegend.AddEntry(h_resp_list[rl], l, 'lf')
      myLegend.Draw()

      name = 'hTheoryRatio_R{}_{}_{}-{}{}'.format(
        self.utils.remove_periods(jetR), obs_label,
        int(min_pt), int(max_pt), self.file_format)
      outputFilename = os.path.join(output_dir, name)
      c.SaveAs(outputFilename)
      c.Close()

    # Reset parton-level range in case used later
    hcent_p.GetXaxis().UnZoom()
    hmin_p.GetXaxis().UnZoom()
    hmax_p.GetXaxis().UnZoom()


  #----------------------------------------------------------------------
  # Plot comparison of F_np-convolved theory to RM projection at hadron+MPI level
  #----------------------------------------------------------------------
  def plot_Fnp_comp(self, jetR, obs_label, obs_setting, grooming_setting):

    label = "R%s_%s" % (str(jetR).replace('.', ''), str(obs_label).replace('.', ''))

    # Directory to save the histograms
    output_dir = os.path.join(self.output_dir_theory, "Fnp_comp")
    if not os.path.exists(output_dir):
      os.mkdir(output_dir)

    n_obs_bins = len(self.theory_obs_bins) - 1
    obs_edges = self.theory_obs_bins

    # red, blue, green, babyblue, purple
    ColorArray = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  #["#db0000", "#0000bf", "#00bd00", "#00b6bd", "#8a008a"]

    # Make projections in pT bins at parton level
    for i, min_pt in list(enumerate(self.theory_pt_bins))[:-1]:
      max_pt = self.theory_pt_bins[i+1]

      '''
      # Initialize canvas & pad for plotting
      name = 'cTheoryComp_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt, max_pt)
      c = ROOT.TCanvas(name, name, 600, 450)
      c.Draw()
      c.cd()

      name = 'pTheoryComp_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt, max_pt)
      myPad = ROOT.TPad(name, name, 0, 0, 1, 1)
      myPad.SetLeftMargin(0.2)
      myPad.SetTopMargin(0.07)
      myPad.SetRightMargin(0.04)
      myPad.SetBottomMargin(0.13)
      if grooming_setting:
        myPad.SetLogy()
      myPad.Draw()
      myPad.cd()
      '''

      x_min = 5e-4 if grooming_setting else 0
      y_min = 1e-3 if grooming_setting else 0
      axMain = plt.gca(); axLine = None
      if grooming_setting:
        plt.yscale("log")
        plt.ylim(y_min, 1e8)

        axMain.set_xscale('linear')
        axMain.set_xlim((0.09, 1))
        axMain.spines['left'].set_visible(False)
        axMain.yaxis.set_ticks_position('right')
        axMain.yaxis.set_visible(False)

        divider = make_axes_locatable(axMain)
        axLin = divider.append_axes("left", size=2.0, pad=0, sharey=axMain)
        axLin.set_xscale('log')
        axLin.set_xlim((x_min, 0.09))
        axLin.spines['right'].set_visible(False)
        axLin.yaxis.set_ticks_position('left')
        plt.setp(axLin.get_xticklabels(), visible=True)
      else:
        plt.xscale("linear")
        plt.yscale("linear")

      h_resp_list = []
      #for ri in range(len(self.theory_response_files)):
      # Only using PYTHIA8 for now
      ri = 0
      # Get event generator histogram
      thn = getattr(self, "hResponse_theory_Fnp_%s_%i" % (label, ri))

      # Get the response matrix prediction
      thn.GetAxis(1).SetRangeUser(int(min_pt), int(max_pt))
      h_response_projection = thn.Projection(3)
      name = "hTheoryComp_Fnp_%s_PtBin%i-%i_%s" % \
             (label, min_pt, max_pt, self.theory_response_labels[ri])
      h_response_projection.SetNameTitle(name, name)
      h_response_projection.SetDirectory(0)

      # Rescale by integral for correct normalization
      h_response_projection.Scale(1/h_response_projection.Integral(), "width")

      '''
      color = self.ColorArray[4+len(self.Omega_list)+ri]
      h_response_projection.SetLineColor(color)
      h_response_projection.SetLineWidth(3)
      h_resp_list.append(h_response_projection)
      '''
      nbinsx = h_response_projection.GetNbinsX()
      if grooming_setting:
        axMain.plot([h_response_projection.GetBinCenter(i) for i in range(1, nbinsx+1)],
                    [h_response_projection.GetBinContent(i) for i in range(1, nbinsx+1)],
                    color=ColorArray[len(self.Omega_list)+ri], label="PYTHIA8")
        axLin.plot([h_response_projection.GetBinCenter(i) for i in range(1, nbinsx+1)],
                   [h_response_projection.GetBinContent(i) for i in range(1, nbinsx+1)],
                   color=ColorArray[len(self.Omega_list)+ri], label="PYTHIA8")
      else:
        plt.plot([h_response_projection.GetBinCenter(i) for i in range(1, nbinsx+1)],
                 [h_response_projection.GetBinContent(i) for i in range(1, nbinsx+1)],
                 color=ColorArray[len(self.Omega_list)+ri], label="PYTHIA8")

      y_maximum = 0
      maxbin = None
      h_theory_list = []
      for oi, Omega in enumerate(self.Omega_list):
        hcent_h = getattr(self, "theory_cent_%s_%s_h_Fnp_Omega%s_%i" % \
                          (self.observable, label, str(Omega), ri))
        hmin_h = getattr(self, "theory_min_%s_%s_h_Fnp_Omega%s_%i" % \
                         (self.observable, label, str(Omega), ri))
        hmax_h = getattr(self, "theory_max_%s_%s_h_Fnp_Omega%s_%i" % \
                         (self.observable, label, str(Omega), ri))

        # Get the theory prediction
        hcent_h.GetXaxis().SetRangeUser(int(min_pt), int(max_pt))
        hmin_h.GetXaxis().SetRangeUser(int(min_pt), int(max_pt))
        hmax_h.GetXaxis().SetRangeUser(int(min_pt), int(max_pt))

        hcent_proj_h = hcent_h.ProjectionY()
        hmin_proj_h = hmin_h.ProjectionY()
        hmax_proj_h = hmax_h.ProjectionY()

        # Fix normalization
        hmin_proj_h.Scale(1/hcent_proj_h.Integral(), "width")
        hmax_proj_h.Scale(1/hcent_proj_h.Integral(), "width")
        hcent_proj_h.Scale(1/hcent_proj_h.Integral(), "width")

        if hmax_proj_h.GetMaximum() > y_maximum:
          y_maximum = hmax_proj_h.GetMaximum()

        # Find the max bin: the last bin where the angularity is 0
        if not maxbin:
          maxbin = hcent_proj_h.GetNbinsX()-1   # initialize
          while hmax_proj_h.GetBinContent(maxbin) < 1e-3 and maxbin > 1:
            maxbin -= 1

        x = array('d', [round(hcent_proj_h.GetXaxis().GetBinCenter(i), 5) for
                        i in range(1, maxbin+1)])
        y = array('d', [hcent_proj_h.GetBinContent(i) for i in range(1, maxbin+1)])
        if grooming_setting:
          axMain.plot(x, y, 'k', color=ColorArray[oi])
          axLin.plot(x, y, 'k', color=ColorArray[oi])
        else:
          plt.plot(x, y, 'k', color=ColorArray[oi])

        xerrup = array('d', [(x[i+1] - x[i]) / 2. for i in range(maxbin-1)] + [0])
        xerrdn = array('d', [0] + [(x[i+1] - x[i]) / 2. for i in range(maxbin-1)])
        # ROOT version
        #yerrup = array('d', [hmax_proj_h.GetBinContent(i)-y[i-1] for i in range(1, maxbin+1)])
        #yerrdn = array('d', [y[i-1]-hmin_proj_h.GetBinContent(i) for i in range(1, maxbin+1)])
        # PyPlot version
        yerrup = [hmax_proj_h.GetBinContent(i) for i in range(1, maxbin+1)]
        yerrdn = [hmin_proj_h.GetBinContent(i) for i in range(1, maxbin+1)]
        if grooming_setting:
          axMain.fill_between(x, yerrdn, yerrup, alpha=0.25, facecolor=ColorArray[oi],
                           label=(r"NLL$^{\prime}$ $\otimes$ F$_\mathrm{NP}^{\Omega=%.1f}$" % Omega))
          axLin.fill_between(x, yerrdn, yerrup, alpha=0.25, facecolor=ColorArray[oi],
                           label=(r"NLL$^{\prime}$ $\otimes$ F$_\mathrm{NP}^{\Omega=%.1f}$" % Omega))
        else:
          plt.fill_between(x, yerrdn, yerrup, alpha=0.25, facecolor=ColorArray[oi],
                           label=(r"NLL$^{\prime}$ $\otimes$ F$_\mathrm{NP}^{\Omega=%.1f}$" % Omega))

        '''
        h_theory = ROOT.TGraphAsymmErrors(maxbin, x, y, xerrdn, xerrup, yerrdn, yerrup)
        color = self.ColorArray[4+oi]
        h_theory.SetFillColorAlpha(color, 0.25)
        h_theory.SetLineColor(color)
        h_theory.SetLineWidth(3)
        h_theory_list.append(h_theory)
        '''

      plt.legend(loc="upper left")

      y_max = 5e2 * y_maximum if grooming_setting else 1.7 * y_maximum
      if y_max <= 0:
        print(y_maximum, y_max)
      ax2 = None
      if grooming_setting:
        axMain.set_xlim((0.09, obs_edges[maxbin]))
        axMain.set_ylim((y_min, y_max))
        axLin.set_ylim((y_min, y_max))
      else:
        plt.xlim(x_min, obs_edges[maxbin])
        plt.ylim(y_min, y_max)
      plt.ylabel(r"$%s$" % self.ytitle.replace("#", "\\").replace("\it", "\\mathit").replace(
        '{d', '{\\mathrm{d}').replace('jet', '\\mathrm{jet}'), fontsize='20')
      plt.xlabel(r"$%s$" % self.xtitle.replace("#", "\\").replace("\it", "\\mathit"), fontsize='16')
      plt.subplots_adjust(left=(0.15 if not grooming_setting else 0.2),
                          bottom=0.15, right=0.98, top=0.98)

      ''' ROOT version
      # Use blank histogram to initialize this range
      bin_array = array('d', obs_edges[0:maxbin+1])
      name = 'hTheoryComp_R{}_{}_{}-{}_Blank'.format(jetR, obs_label, min_pt, max_pt)
      myBlankHisto = ROOT.TH1F(name, name, maxbin, bin_array)
      myBlankHisto.SetNdivisions(505)
      myBlankHisto.SetXTitle(self.xtitle)
      myBlankHisto.GetXaxis().SetTitleOffset(1.02)
      myBlankHisto.GetXaxis().SetTitleSize(0.055)
      myBlankHisto.SetYTitle(self.ytitle)
      myBlankHisto.GetYaxis().SetTitleOffset(1.1)
      myBlankHisto.GetYaxis().SetTitleSize(0.055)
      myBlankHisto.SetMinimum(y_min)
      myBlankHisto.SetMaximum(y_max)
      myBlankHisto.Draw()

      for h_theory in h_theory_list:
        h_theory.SetMinimum(y_min)
        h_theory.Draw('L 3 same')
      for h_response_projection in h_resp_list:
        h_response_projection.SetMinimum(y_min)
        h_response_projection.Draw('L hist same')
      '''

      # Reset thn range in case used later
      thn.GetAxis(1).UnZoom()

      # matplotlib version
      text = '\n'.join((
        "ALICE {}".format(self.figure_approval_status),
        r"pp $\sqrt{s} = 5.02$ TeV",
        r"anti-${k}_{\mathrm{T}}$ jets,   $R = %s$" % str(jetR),
        r"$%i < p_\mathrm{T,jet}^\mathrm{full} < %i$ GeV/$c$" % (min_pt, max_pt),
        r"$| \eta_\mathrm{jet}| < %s$" % str(0.9 - jetR)))
      subobs_label = self.utils.formatted_subobs_label(self.observable)
      if subobs_label:
        text += r',   $%s = %s$' % (
          subobs_label.replace("#", "\\").replace("\it", "\\mathit"), obs_setting)
      plt.text(0.52*obs_edges[maxbin], 0.65*y_max if
               not grooming_setting else 3*math.log(y_max/y_min),
               text, fontsize='14', fontname='sans-serif')


      ''' ROOT version
      text_latex = ROOT.TLatex()
      text_latex.SetNDC()
      text_xval = 0.61
      text = 'ALICE {}'.format(self.figure_approval_status)
      text_latex.DrawLatex(text_xval, 0.87, text)

      text = 'pp #sqrt{#it{s}} = 5.02 TeV'
      text_latex.SetTextSize(0.045)
      text_latex.DrawLatex(text_xval, 0.8, text)

      text = "anti-#it{k}_{T} jets,   #it{R} = %s" % str(jetR)
      text_latex.DrawLatex(text_xval, 0.73, text)

      text = str(min_pt) + ' < #it{p}_{T,jet}^{full} < ' + str(max_pt) + ' GeV/#it{c}'
      text_latex.DrawLatex(text_xval, 0.66, text)

      text = '| #it{#eta}_{jet}| < %s' % str(0.9 - jetR)
      subobs_label = self.utils.formatted_subobs_label(self.observable)
      if subobs_label:
        text += ',   %s = %s' % (subobs_label, obs_setting)
      delta = 0.07
      text_latex.DrawLatex(text_xval, 0.66-delta, text)

      myLegend = ROOT.TLegend(0.27, 0.7, 0.55, 0.9)
      self.utils.setup_legend(myLegend, 0.035)
      for oi, Omega in enumerate(self.Omega_list):
        myLegend.AddEntry(h_theory_list[oi],
                          "NLL' #otimes #it{F}_{NP}^{#Omega=%s}" % str(Omega), 'lf')
      #for rl, l in enumerate(self.theory_response_labels):
      rl = 0; l = self.theory_response_labels[rl]
      myLegend.AddEntry(h_resp_list[rl], l, 'lf')
      myLegend.Draw()
      '''

      if grooming_setting:
        ax2 = axLin.twinx()
        ax2.spines['left'].set_visible(False)
        ax2.tick_params(axis='y',which='both',labelright='off')

      name = 'hTheoryRatio_Fnp_R{}_{}_{}-{}{}'.format(
        self.utils.remove_periods(jetR), obs_label,
        int(min_pt), int(max_pt), self.file_format)
      outputFilename = os.path.join(output_dir, name)
      #c.SaveAs(outputFilename)
      #c.Close()
      plt.savefig(outputFilename)
      plt.close()

    # Reset hadron-level ranges in case used later
    for oi, Omega in enumerate(self.Omega_list):
      ri = 0
      hcent_h = getattr(self, "theory_cent_%s_%s_h_Fnp_Omega%s_%i" % \
                        (self.observable, label, str(Omega), ri))
      hmin_h = getattr(self, "theory_min_%s_%s_h_Fnp_Omega%s_%i" % \
                       (self.observable, label, str(Omega), ri))
      hmax_h = getattr(self, "theory_max_%s_%s_h_Fnp_Omega%s_%i" % \
                       (self.observable, label, str(Omega), ri))
      hcent_h.GetXaxis().UnZoom()
      hmin_h.GetXaxis().UnZoom()
      hmax_h.GetXaxis().UnZoom()


  #----------------------------------------------------------------------
  # Get unfolded data from the previous preliminary result (40-60 GeV/c)
  def get_h_prelim(self, jetR):
    
    if jetR == 0.2:
      xedges = [0., 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.12]
      yvals = [2.892837, 11.9108, 17.3579, 17.65965, 15.25709, 13.00818, 9.1359, 2.471203]
      staterror = [0.09723696, 0.2879163, 0.3482209, 0.3487025, 0.3212577,
                   0.2975396, 0.2503627, 0.06595427]
      syserrorl = [0.1654749, 0.4376057, 0.536369, 0.1987916, 0.3712922,
                   0.3223265, 0.3906225, 0.1588837]
      syserrorh = [0.1673992, 0.4359767, 0.5354239, 0.200042, 0.3804927,
                   0.3368305, 0.3948841, 0.1588432]
    else:  #jetR == 0.4:
      xedges = [0, 0.02, 0.03, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.18, 0.22]
      yvals = [0.6078014, 2.815131, 5.960223, 9.770085, 9.899658, 8.603309,
               6.119539, 4.60788, 2.300467, 0.7015587]
      staterror = [0.0430097, 0.1299623, 0.1900576, 0.1710785, 0.1712931, 0.1581808,
                   0.1320032, 0.1172254, 0.059633, 0.03359087]
      syserrorl = [0.1025423, 0.1138034, 0.04146903, 0.096956, 0.06155705, 0.06077894,
                   0.08091901, 0.06775198, 0.03015912, 0.04115888]
      syserrorh = [0.1024212, 0.1204349, 0.07618093, 0.1491075, 0.07706482, 0.03761498,
                   0.1431532, 0.1033103, 0.02661073, 0.0411509]

    # Scale values by R due to observable definition
    xedges_scaled = array('d', [round(val / jetR, 2) for val in xedges])
    setattr(self, "xedges_prev_prelim_%s" % jetR, xedges_scaled)
    scale_factor = [(xedges[i+1]-xedges[i]) / (xedges_scaled[i+1]-xedges_scaled[i])
                    for i in range(0, len(yvals))]
    yvals_scaled = [yvals[i] * scale_factor[i] for i in range(0, len(yvals))]
    staterror_scaled = [staterror[i] * scale_factor[i] for i in range(0, len(staterror))]
    syserrorl_scaled = [syserrorl[i] * scale_factor[i] for i in range(0, len(syserrorl))]
    syserrorh_scaled = [syserrorh[i] * scale_factor[i] for i in range(0, len(syserrorh))]

    # Create histograms
    name = "ALI−PREL−339374"
    hCompStat = ROOT.TH1D(name, name, len(yvals), xedges_scaled)
    hCompSys = ROOT.TH1D(name+'_sys', name+'_sys', len(yvals), xedges_scaled)
    for i in range(1, len(xedges), 1):
      hCompStat.SetBinContent(i, yvals_scaled[i-1])
      hCompSys.SetBinContent(i, yvals_scaled[i-1])
      hCompStat.SetBinError(i, staterror_scaled[i-1])
      hCompSys.SetBinError(i, max(syserrorl_scaled[i-1], syserrorh_scaled[i-1]))

    return hCompStat, hCompSys

  #----------------------------------------------------------------------
  def plot_obs_comp(self, jetR, obs_label, obs_setting, grooming_setting,
                    min_pt_truth, max_pt_truth, maxbin):
    
    # Scale both distributions by integrals
    scale_by_int = False

    name = 'cResult_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt_truth, max_pt_truth)
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
    
    xtitle = getattr(self, 'xtitle')
    ytitle = getattr(self, 'ytitle')
    color = 600-6
    
    # Get histograms
    name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label,
                                              min_pt_truth, max_pt_truth)
    if grooming_setting and show_everything_else:
      fraction_tagged = getattr(self, 'tagging_fraction_R{}_{}_{}-{}'.format(
        jetR, obs_label, min_pt_truth, max_pt_truth))
      #fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
    h = self.truncate_hist(getattr(self, name), maxbin, name+'_trunc')
    h.SetMarkerSize(1.5)
    h.SetMarkerStyle(20)
    h.SetMarkerColor(color)
    h.SetLineStyle(1)
    h.SetLineWidth(2)
    h.SetLineColor(color)
    if scale_by_int:
      h.Scale(1/h.Integral())
    
    h_sys = getattr(self, 'hResult_{}_systotal_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth))
    h_sys.SetLineColor(0)
    h_sys.SetFillColor(color)
    h_sys.SetFillColorAlpha(color, 0.3)
    h_sys.SetFillStyle(1001)
    h_sys.SetLineWidth(0)
    if scale_by_int:
      h_sys.Scale(1/h_sys.Integral())
    
    n_obs_bins_truth = self.n_bins_truth(obs_label)
    truth_bin_array = self.truth_bin_array(obs_label)
    if maxbin:
      truth_bin_array = truth_bin_array[0:maxbin+1]
      if jetR == 0.2:
        truth_bin_array[-1] = 0.6
      n_obs_bins_truth = len(truth_bin_array)-1
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_obs_bins_truth, truth_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle(xtitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle(ytitle)
    if scale_by_int:
      myBlankHisto.SetYTitle(ytitle + ' / #int(' + ytitle + ')')
    myBlankHisto.SetMaximum(2*h.GetMaximum())
    if self.observable == 'subjet_z' or self.observable == 'jet_axis':
      myBlankHisto.SetMaximum(1.5*h.GetMaximum())
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")

    hCompStat, hCompSys = self.get_h_prelim(jetR)

    #formatting
    hCompStat.SetFillStyle(0)
    hCompStat.SetMarkerSize(1.5)
    hCompStat.SetMarkerStyle(21)
    hCompStat.SetMarkerColor(1)
    hCompStat.SetLineColor(1)
    hCompStat.SetLineWidth(1)
    hCompSys.SetLineColor(0)
    hCompSys.SetFillColor(1)
    hCompSys.SetFillColorAlpha(1, 0.3)
    hCompSys.SetFillStyle(1001)
    hCompSys.SetLineWidth(0)
    if scale_by_int:
      hCompStat.Scale(1/hCompStat.Integral())
      hCompSys.Scale(1/hCompSys.Integral())

    hCompSys.Draw('E2 same')
    hCompStat.Draw('PE X0 same')
    h_sys.DrawCopy('E2 same')
    h.DrawCopy('PE X0 same')
  
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(0.57, 0.87, text)
    
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.8, text)

    text = str(min_pt_truth) + ' < #it{p}_{T,jet}^{ch} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.DrawLatex(0.57, 0.73, text)

    text = '#it{R} = ' + str(jetR) + ',   | #it{#eta}_{jet}| < %s' % str(0.9 - jetR)
    text_latex.DrawLatex(0.57, 0.66, text)
    
    subobs_label = self.utils.formatted_subobs_label(self.observable)
    delta = 0.
    if subobs_label:
      text = '{} = {}'.format(subobs_label, obs_setting)
      text_latex.DrawLatex(0.57, 0.59, text)
      delta = 0.07
    
    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting)#.replace("#beta}", "#beta}_{SD}")
      text_latex.DrawLatex(0.57, 0.59-delta, text)
      
      text_latex.SetTextSize(0.04)
      text = '#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged
      text_latex.DrawLatex(0.57, 0.52-delta, text)

    myLegend = ROOT.TLegend(0.25, 0.7, 0.45, 0.85)
    self.utils.setup_legend(myLegend,0.035)
    myLegend.AddEntry(h, 'This measurement', 'pe')
    myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    myLegend.AddEntry(hCompStat, 'ALI-PREL-339374', 'pe')
    myLegend.Draw()

    name = 'hUnfoldedComp_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR),
                                                 obs_label, int(min_pt_truth),
                                                 int(max_pt_truth), self.file_format)
    output_dir = getattr(self, 'output_dir_final_results')
    output_dir_single = output_dir + '/single_results'
    if not os.path.exists(output_dir_single):
      os.mkdir(output_dir_single)
    outputFilename = os.path.join(output_dir_single, name)
    c.SaveAs(outputFilename)
    c.Close()

    # Write result to ROOT file
    final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
    fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
    h.Write()
    h_sys.Write()
    hCompStat.Write()
    hCompSys.Write()
    fFinalResults.Close()

  #----------------------------------------------------------------------
  def MC_prediction(self, jetR, obs_setting, obs_label, min_pt_truth,
                    max_pt_truth, maxbin, MC='Pythia', overlay=False):
  
    if MC.lower() == 'pythia':
      hMC = self.get_pythia_from_response(jetR, obs_label, min_pt_truth,
                                          max_pt_truth, maxbin, overlay)
    elif MC.lower() == 'herwig':
      hMC = self.get_herwig_from_response(jetR, obs_label, min_pt_truth,
                                          max_pt_truth, maxbin, overlay)
    else:
      raise NotImplementedError("MC must be Pythia or Herwig.")

    n_jets_inclusive = hMC.Integral(0, hMC.GetNbinsX()+1)
    n_jets_tagged = hMC.Integral(hMC.FindBin(
      self.truth_bin_array(obs_label)[0]), hMC.GetNbinsX())

    fraction_tagged_MC =  n_jets_tagged/n_jets_inclusive
    hMC.Scale(1./n_jets_inclusive, 'width')
      
    return [hMC, fraction_tagged_MC]

  #----------------------------------------------------------------------
  def get_pythia_from_response(self, jetR, obs_label, min_pt_truth, max_pt_truth,
                               maxbin, overlay=False):

    output_dir = getattr(self, 'output_dir_main')

    prev_prelim = False
    if self.use_prev_prelim and overlay and (jetR == 0.2 or jetR == 0.4) \
       and min_pt_truth == 40 and obs_label == '1':
      prev_prelim = True
      # Need to rebin response for the binning used by previous preliminary result
      filepath = os.path.join(output_dir, 'response_prev_prelim.root')

      if not os.path.exists(filepath):
        # Create rebinned THn with these binnings, and write to file
        print("Rebinning response matrix for previous preliminary masurement...")
        name_thn = self.utils.name_thn(self.observable, jetR, obs_label)
        name_thn_rebinned = self.utils.name_thn_rebinned(self.observable, jetR, obs_label)
        name_roounfold = 'roounfold_response_R{}_{}'.format(jetR, obs_label)
        thn = ROOT.TFile(self.main_response, 'READ').Get(name_thn)
        thn.SetName(name_thn)
        label = 'R{}_{}'.format(jetR, obs_label)
        pt_bins_truth = array('d', [5, 20, 40, 60, 80, 100, 150, 200])
        pt_bins_det = array('d', [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 120, 150])
        obs_bins = getattr(self, "xedges_prev_prelim_%s" % jetR)
        self.utils.rebin_response(
          filepath, thn, name_thn_rebinned, name_roounfold, label, len(pt_bins_det)-1,
          pt_bins_det, len(obs_bins)-1, obs_bins, len(pt_bins_truth)-1, pt_bins_truth,
          len(obs_bins)-1, obs_bins, self.observable, do_roounfoldresponse=False)
    else:
      filepath = os.path.join(output_dir, 'response.root')
    f = ROOT.TFile(filepath, 'READ')

    thn_name = 'hResponse_JetPt_{}_R{}_{}_rebinned'.format(self.observable, jetR, obs_label)
    thn = f.Get(thn_name)
    thn.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)

    name = 'hPythia_{}_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    if prev_prelim:
      h = thn.Projection(3)
    else:
      h = self.truncate_hist(thn.Projection(3), maxbin, name)
    h.SetDirectory(0)

    return h

  #----------------------------------------------------------------------
  def get_herwig_from_response(self, jetR, obs_label, min_pt_truth, max_pt_truth,
                               maxbin, overlay=False):

    filepath = os.path.join(self.output_dir_fastsim_generator1, 'response.root')
    f = ROOT.TFile(filepath, 'READ')

    thn_name = 'hResponse_JetPt_{}_R{}_{}_rebinned'.format(self.observable, jetR, obs_label)
    thn = f.Get(thn_name)
    thn.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)

    name = 'hHerwig_{}_R{}_{}_{}-{}'.format(
      self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    h = self.truncate_hist(thn.Projection(3), maxbin, name)
    h.SetDirectory(0)

    return h

  #----------------------------------------------------------------------
  def plot_final_result_overlay(self, i_config, jetR, overlay_list):
    print('Plotting overlay of', overlay_list)

    # Plot overlay of different subconfigs, for fixed pt bin
    for i in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[i]
      max_pt_truth = self.pt_bins_reported[i+1]
      maxbins = [self.obs_max_bins(obs_label)[i] for obs_label in self.obs_labels]

      # Plot PYTHIA
      self.plot_observable_overlay_subconfigs(
        i_config, jetR, overlay_list, min_pt_truth,
        max_pt_truth, maxbins, plot_MC=True, MC='PYTHIA', plot_ratio=True)

      if self.do_theory and not self.use_prev_prelim and \
         'fastsim_generator1' in self.systematics_list:
        self.plot_observable_overlay_subconfigs(
          i_config, jetR, overlay_list, min_pt_truth,
          max_pt_truth, maxbins, plot_MC=True, MC='Herwig', plot_ratio=True)


  #----------------------------------------------------------------------
  def plot_observable_overlay_subconfigs(self, i_config, jetR, overlay_list, min_pt_truth,
                                         max_pt_truth, maxbins, plot_MC=False,
                                         MC='PYTHIA', plot_nll=False, plot_ratio=False):

    # Flag to plot ratio all on the same scale, 0 to 2.2
    plot_ratio_same_scale = True

    name = 'cResult_overlay_R{}_allpt_{}-{}'.format(jetR, min_pt_truth, max_pt_truth)
    if plot_ratio:
      c = ROOT.TCanvas(name, name, 600, 650)
    else:
      c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()
    
    c.cd()
    if plot_ratio:
      pad1 = ROOT.TPad('myPad', 'The pad',0,0.3,1,1)
    else:
      pad1 = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    pad1.SetLeftMargin(0.2)
    pad1.SetTopMargin(0.07)
    pad1.SetRightMargin(0.04)
    pad1.SetBottomMargin(0.13)
    pad1.SetTicks(0,1)
    if plot_ratio:
      pad1.SetBottomMargin(0.)
    # set log y axis if all configs are SD
    setlogy = False
    for name in overlay_list:
      if "SD" not in name:
        break
      elif name == overlay_list[-1]:
        setlogy = True
        pad1.SetLogy()
    pad1.Draw()
    pad1.cd()

    myLegend = ROOT.TLegend(0.22, 0.66, 0.55, 0.91)
    self.utils.setup_legend(myLegend, 0.045)
    myLegend2 = ROOT.TLegend(0.43, 0.786, 0.65, 0.91)
    self.utils.setup_legend(myLegend2, 0.045)
    
    name = 'hmain_{}_R{}_{{}}_{}-{}'.format(self.observable, jetR, min_pt_truth, max_pt_truth)
    ymax, ymin = self.get_max_min(name, overlay_list, maxbins)

    h_list = []
    text_list = []

    for i, subconfig_name in enumerate(self.obs_subconfig_list):
    
      if subconfig_name not in overlay_list:
        continue

      obs_setting = self.obs_settings[i]
      grooming_setting = self.grooming_settings[i]
      obs_label = self.obs_labels[i]
      maxbin = maxbins[i]
      
      if subconfig_name == overlay_list[0]:
        marker = 20
        marker_pythia = 24
        color = 1
      elif subconfig_name == overlay_list[1]:
        marker = 21
        marker_pythia = 25
        color = 600-6
      elif subconfig_name == overlay_list[2]:
        marker = 33
        marker_pythia = 27
        color = 632-4
      else:  # subconfig_name == overlay_list[3]:
        marker = 34
        marker_pythia = 28
        color = 416-2

      name = 'hmain_{}_R{}_{}_{}-{}'.format(
        self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
      if grooming_setting:
        fraction_tagged = getattr(self, 'tagging_fraction_R{}_{}_{}-{}'.format(
          jetR, obs_label, min_pt_truth, max_pt_truth))

      if self.use_prev_prelim and (jetR == 0.2 or jetR == 0.4) \
         and min_pt_truth == 40 and obs_label == '1':
        # Use previous preliminary result
        h, h_sys = self.get_h_prelim(jetR)
        # Move error bars to different histogram
        h_sys.SetNameTitle("hSysPrelim%s", "hSysPrelim%s")
        for i in range(1, h.GetNbinsX()+1, 1):
          h.SetBinError(i, 0.)
      else:
        if grooming_setting and maxbin:
          h = self.truncate_hist(getattr(self, name), maxbin+1, name+'_trunc')
        else:
          h = self.truncate_hist(getattr(self, name), maxbin, name+'_trunc')
        h_sys = getattr(self, 'hResult_{}_systotal_R{}_{}_{}-{}'.format(
          self.observable, jetR, obs_label, min_pt_truth, max_pt_truth))

      h.SetDirectory(0)
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(marker)
      h.SetMarkerColor(color)
      h.SetLineStyle(1)
      h.SetLineWidth(2)
      h.SetLineColor(color)
      
      h_sys.SetLineColor(0)
      h_sys.SetFillColor(color)
      h_sys.SetFillColorAlpha(color, 0.3)
      h_sys.SetFillStyle(1001)
      h_sys.SetLineWidth(0)
      
      if subconfig_name == overlay_list[0]:

        pad1.cd()
        xtitle = getattr(self, 'xtitle')
        ytitle = getattr(self, 'ytitle')
        xmin = self.obs_config_dict[subconfig_name]['obs_bins_truth'][0]
        xmax = self.obs_config_dict[subconfig_name]['obs_bins_truth'][-1]
        if maxbin:
          if self.use_prev_prelim and (jetR == 0.2 or jetR == 0.4) \
             and min_pt_truth == 40 and obs_label == '1':
            xmax = getattr(self, "xedges_prev_prelim_%s" % jetR)[-1]
          else:
            xmax = self.obs_config_dict[subconfig_name]['obs_bins_truth'][maxbin]
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, xmin, xmax)
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(xtitle)
        myBlankHisto.GetYaxis().SetTitleOffset(1.3)
        myBlankHisto.SetYTitle(ytitle)
        if jetR == 0.2:
          if min_pt_truth == 20:
            myBlankHisto.SetMaximum(1.1*ymax)
          elif min_pt_truth == 40:
            myBlankHisto.SetMaximum(1.2*ymax)
          elif min_pt_truth == 60:
            myBlankHisto.SetMaximum(1.15*ymax)
          else:
            myBlankHisto.SetMaximum(1.3*ymax)
        elif jetR == 0.4:
          if min_pt_truth == 20:
            myBlankHisto.SetMaximum(1.5*ymax)
          elif min_pt_truth == 40:
            myBlankHisto.SetMaximum(1.35*ymax)
          elif min_pt_truth == 60:
            myBlankHisto.SetMaximum(1.15*ymax)
          else:
            myBlankHisto.SetMaximum(1.5*ymax)
        else:
          myBlankHisto.SetMaximum(1.5*ymax)
        if setlogy:
          myBlankHisto.SetMaximum(5*ymax)
        myBlankHisto.SetMinimum(0.)
        if plot_ratio:
          myBlankHisto.SetMinimum(2e-4) # Don't draw 0 on top panel
          if setlogy:  # the minimum value matters
            myBlankHisto.SetMinimum(2*ymin/3)
          myBlankHisto.GetYaxis().SetTitleSize(0.065)
          myBlankHisto.GetYaxis().SetTitleOffset(1.1)
          myBlankHisto.GetYaxis().SetLabelSize(0.06)
        myBlankHisto.Draw('E')
        
        # Plot ratio
        if plot_ratio:
          
          c.cd()
          pad2 = ROOT.TPad("pad2", "pad2", 0, 0.02, 1, 0.3)
          pad2.SetTopMargin(0)
          pad2.SetBottomMargin(0.5)
          pad2.SetLeftMargin(0.2)
          pad2.SetRightMargin(0.04)
          pad2.SetTicks(1,1)
          pad2.Draw()
          pad2.cd()
          
          myBlankHisto2 = myBlankHisto.Clone("myBlankHisto_C")
          myBlankHisto2.SetYTitle("#frac{Data}{%s}" % MC)
          myBlankHisto2.SetXTitle(xtitle)
          myBlankHisto2.GetXaxis().SetTitleSize(30)
          myBlankHisto2.GetXaxis().SetTitleFont(43)
          myBlankHisto2.GetXaxis().SetTitleOffset(4.)
          myBlankHisto2.GetXaxis().SetLabelFont(43)
          myBlankHisto2.GetXaxis().SetLabelSize(25)
          myBlankHisto2.GetYaxis().SetTitleSize(20)
          myBlankHisto2.GetYaxis().SetTitleFont(43)
          myBlankHisto2.GetYaxis().SetTitleOffset(2.2)
          myBlankHisto2.GetYaxis().SetLabelFont(43)
          myBlankHisto2.GetYaxis().SetLabelSize(25)
          myBlankHisto2.GetYaxis().SetNdivisions(505)
          if plot_ratio_same_scale:
            if jetR == 0.2:
              myBlankHisto2.GetYaxis().SetRangeUser(0.5, 1.75)
            elif jetR == 0.4:
              myBlankHisto2.GetYaxis().SetRangeUser(0.5, 1.9)
            else: 
              myBlankHisto2.GetYaxis().SetRangeUser(0, 2.2)
          elif jetR == 0.2:
            if min_pt_truth == 20:
              myBlankHisto2.GetYaxis().SetRangeUser(0.6, 1.75)
            elif min_pt_truth == 40:
              myBlankHisto2.GetYaxis().SetRangeUser(0.78, 1.299)
            elif min_pt_truth == 60:
              myBlankHisto2.GetYaxis().SetRangeUser(0.55, 1.499)
            else: 
              myBlankHisto2.GetYaxis().SetRangeUser(0.5, 1.99)
          elif jetR == 0.4:
            if min_pt_truth == 20:
              myBlankHisto2.GetYaxis().SetRangeUser(0.81, 1.72)
            elif min_pt_truth == 40:
              myBlankHisto2.GetYaxis().SetRangeUser(0.7, 2.1)
            elif min_pt_truth == 60:
              myBlankHisto2.GetYaxis().SetRangeUser(0.75, 1.55)
            else: 
              myBlankHisto2.GetYaxis().SetRangeUser(0.5, 1.99)
          else: 
            myBlankHisto2.GetYaxis().SetRangeUser(0.5, 1.99)
          myBlankHisto2.Draw()
        
          line = ROOT.TLine(0,1,xmax,1)
          line.SetLineColor(920+2)
          line.SetLineStyle(2)
          line.Draw()
      
      hMC = None; fraction_tagged_MC = None;
      if plot_MC:
        if MC.lower() == "pythia":
          if grooming_setting and maxbin:
            hMC, fraction_tagged_MC = self.MC_prediction(
              jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin+1,
              MC='Pythia', overlay=True)
          else:
            hMC, fraction_tagged_MC = self.MC_prediction(
              jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin,
              MC='Pythia', overlay=True)

        elif MC.lower() == "herwig":
          if grooming_setting:
            hMC, fraction_tagged_MC = self.MC_prediction(
              jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin+1,
              'Herwig', overlay=True)
          else:
            hMC, fraction_tagged_MC = self.MC_prediction(
              jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin,
              'Herwig', overlay=True)

        else:
          raise NotImplementedError("MC must be either Pythia or Herwig.")

        plot_errors = False
        if plot_errors:
          hMC.SetMarkerSize(0)
          hMC.SetMarkerStyle(0)
          hMC.SetMarkerColor(color)
          hMC.SetFillColor(color)
        else:
          hMC.SetLineColor(color)
          hMC.SetLineColorAlpha(color, 0.5)
          hMC.SetLineWidth(4)

      if plot_ratio:
        hRatioSys = h_sys.Clone()
        hRatioSys.SetName('{}_Ratio'.format(h_sys.GetName()))
        if plot_MC:
          hRatioSys.Divide(hMC)
        hRatioSys.SetLineColor(0)
        hRatioSys.SetFillColor(color)
        hRatioSys.SetFillColorAlpha(color, 0.3)
        hRatioSys.SetFillStyle(1001)
        hRatioSys.SetLineWidth(0)
        hRatioSys.SetMaximum(1.99)
          
        hRatioStat = h.Clone()
        hRatioStat.SetName('{}_Ratio'.format(h.GetName()))
        if plot_MC:
          hRatioStat.Divide(hMC)
        hRatioStat.SetMarkerSize(1.5)
        hRatioStat.SetMarkerStyle(marker)
        hRatioStat.SetMarkerColor(color)
        hRatioStat.SetLineStyle(1)
        hRatioStat.SetLineWidth(2)
        hRatioStat.SetLineColor(color)
        hRatioStat.SetMaximum(1.99)

      pad1.cd()
      if plot_MC:
        plot_errors = False
        if plot_errors:
          hMC.DrawCopy('E3 same')
        else:
          hMC.DrawCopy('L hist same')

      h_sys.DrawCopy('E2 same')
      h.DrawCopy('PE X0 same')
      
      if plot_ratio:
        pad2.cd()
        if plot_MC:
          hRatioSys.DrawCopy('E2 same')
          hRatioStat.DrawCopy('PE X0 same')

      subobs_label = self.utils.formatted_subobs_label(self.observable)
      text = ''
      if subobs_label:
        text += '%s = %s' % (subobs_label, obs_setting)
      text_list.append(text)
      h_list.append(h)
        
    pad1.cd()
    for i, h, text in zip(range(len(h_list)), h_list, text_list):
      if i < 2:
        myLegend.AddEntry(h, text, 'pe')
      else:
        myLegend2.AddEntry(h, text, 'pe')
    myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    if plot_MC:
      if MC.lower() == "pythia":
        myLegend.AddEntry(hMC, 'PYTHIA8 Monash2013', 'l')
      elif MC.lower() == "herwig":
        myLegend.AddEntry(hMC, 'Herwig7 Default', 'l')

    text_xval = 0.63
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(text_xval, 0.87, text)
    
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(text_xval, 0.81, text)

    text = 'Charged jets   anti-#it{k}_{T}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(text_xval, 0.75, text)
    
    text = '#it{R} = ' + str(jetR) + ',   | #it{#eta}_{jet}| < %s' % str(0.9 - jetR)
    text_latex.DrawLatex(text_xval, 0.69, text)
    
    text = str(min_pt_truth) + ' < #it{p}_{T,jet}^{ch} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(text_xval, 0.63, text)

    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting)#.replace("#beta}", "#beta}_{SD}")
      text_latex.DrawLatex(text_xval, 0.57, text)

    myLegend.Draw()
    if len(h_list) > 2:
      myLegend2.Draw()

    name = 'h_{}_R{}_{}-{}_{}{}'.format(self.observable, 
                                        self.utils.remove_periods(jetR), int(min_pt_truth), 
                                        int(max_pt_truth), i_config, self.file_format)
    if plot_MC:
      name = 'h_{}_R{}_{}-{}_{}_{}{}'.format(self.observable, self.utils.remove_periods(jetR),
                                                 int(min_pt_truth), int(max_pt_truth), MC,
                                                 i_config, self.file_format)


    output_dir = getattr(self, 'output_dir_final_results')
    if not os.path.exists(os.path.join(output_dir, 'all_results')):
      os.mkdir(os.path.join(output_dir, 'all_results'))
    outputFilename = os.path.join(output_dir, 'all_results', name)
    c.SaveAs(outputFilename)

    # Write result to ROOT file
    final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
    fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
    c.Write()

    c.Close()

  #----------------------------------------------------------------------
  # Return maximum & minimum y-values of unfolded results in a subconfig list
  def get_max_min(self, name, overlay_list, maxbins):

    total_min = 1e10
    total_max = -1e10

    for i, subconfig_name in enumerate(self.obs_subconfig_list):
    
      if subconfig_name not in overlay_list:
        continue

      obs_setting = self.obs_settings[i]
      grooming_setting = self.grooming_settings[i]
      obs_label = self.utils.obs_label(obs_setting, grooming_setting)
      maxbin = maxbins[i]
      
      h = getattr(self, name.format(obs_label))
      if 'SD' in obs_label:
        content = [ h.GetBinContent(j) for j in range(2, maxbin+2) ]
      else:
        content = [ h.GetBinContent(j) for j in range(1, maxbin+1) ]

      min_val = min(content)
      if min_val < total_min:
        total_min = min_val
      max_val = max(content)
      if max_val > total_max:
        total_max = max_val

        
    return (total_max, total_min)


#----------------------------------------------------------------------
if __name__ == '__main__':

  # Define arguments
  parser = argparse.ArgumentParser(description='Jet substructure analysis')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='analysis_config.yaml',
                      help='Path of config file for analysis')

  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('configFile: \'{0}\''.format(args.configFile))
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = RunAnalysisAng(config_file = args.configFile)
  analysis.run_analysis()
