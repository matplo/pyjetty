#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of MC track information
  and do jet-finding, and save response histograms.

  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse

# Data analysis and plotting
import numpy as np
import ROOT
import yaml
from array import *

# Fastjet via python (from external library heppy)
import fjext

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_io_emb
from pyjetty.alice_analysis.process.base import jet_info
from pyjetty.alice_analysis.process.user.substructure import process_mc_base
from pyjetty.alice_analysis.process.base import thermal_generator
from pyjetty.mputils import CEventSubtractor

################################################################
class ProcessMC_ang(process_mc_base.ProcessMCBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):

    # Initialize base class
    super(ProcessMC_ang, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)

    self.pt_bins = array('d', list(range(5, 100, 5)) + list(range(100, 210, 10)))
    self.obs_bins_ang = np.concatenate((np.linspace(0, 0.009, 10), np.linspace(0.01, 0.1, 19),
                                        np.linspace(0.11, 0.8, 70)))
    self.obs_bins_mass = array('d', list(range(0, 61, 1)))

    # Override default behavior to create delta-observable histograms in Pb-Pb case
    self.fill_delta_obs = True

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects_R(self, jetR):

      for observable in self.observable_list:
        # Should only be two: observable == "ang" or "mass"
        if observable != "ang" and observable != "mass":
          raise ValueError("Observable %s is not implemented in this script" % observable)

        obs_name = self.obs_names[observable]

        # Loop over subobservable (alpha value)
        for i in range(len(self.obs_settings[observable])):

          obs_setting = self.obs_settings[observable][i]
          grooming_setting = self.obs_grooming_settings[observable][i]
          obs_label = self.utils.obs_label(obs_setting, grooming_setting)

          # Create RM histograms
          if self.is_pp:
            self.create_histograms(observable, jetR, obs_label)

          else:

            max_distance = self.max_distance if isinstance(self.max_distance, list) else \
                           self.max_distance[jetR]
            for R_max in max_distance:
              self.create_histograms(observable, jetR, obs_label, R_max)

          if self.thermal_model:
            max_distance = self.max_distance if isinstance(self.max_distance, list) else \
                           self.max_distance[jetR]
            for R_max in max_distance:
              name = 'h_%s_JetPt_R%s_%s_Rmax%s' % (observable, jetR, obs_label, R_max) if \
                len(obs_label) else ('h_%s_JetPt_R%s_Rmax%s' % (observable, jetR, R_max))
              h = ROOT.TH2F(name, name, 200, 0, 200, 100, 0, 1)
              h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
              h.GetYaxis().SetTitle(obs_name)
              setattr(self, name, h)

          name = ('h_%s_JetPt_Truth_R%s_%s' % (observable, jetR, obs_label)) if \
            len(obs_label) else ('h_%s_JetPt_Truth_R%s' % (observable, jetR))
          h = ROOT.TH2F(name, name, 40, 0, 200, 100, 0, 1)
          h.GetXaxis().SetTitle('#it{p}_{T,truth}^{ch jet}')
          h.GetYaxis().SetTitle(obs_name + '^{truth}')
          setattr(self, name, h)

  #---------------------------------------------------------------
  # Create angularity response histograms
  #---------------------------------------------------------------
  def create_histograms(self, observable, jetR, obs_label, R_max = None):

    if R_max:
      suffix = '_Rmax' + str(R_max)
    else:
      suffix = ''

    # LaTeX formatted observable name
    obs_name = self.obs_names[observable]

    # Retrieve binnings from memory
    pt_bins = self.pt_bins
    obs_bins = getattr(self, "obs_bins_" + observable)

    # Create THn of response for ang
    if self.fill_RM_histograms:
      dim = 4;
      title = ['#it{p}_{T,det}^{ch jet}', '#it{p}_{T,truth}^{ch jet}',
               obs_name + '^{det}',       obs_name + '^{truth}']
      nbins  = [len(pt_bins)-1, len(pt_bins)-1, len(obs_bins)-1, len(obs_bins)-1]
      min_li = [pt_bins[0],     pt_bins[0],     obs_bins[0],     obs_bins[0]    ]
      max_li = [pt_bins[-1],    pt_bins[-1],    obs_bins[-1],    obs_bins[-1]   ]

      name = ('hResponse_JetPt_%s_R%s_%s%s' % (observable, jetR, obs_label, suffix)) if \
        len(obs_label) else ('hResponse_JetPt_%s_R%s%s' % (observable, jetR, suffix))
      nbins = (nbins)
      xmin = (min_li)
      xmax = (max_li)
      nbins_array = array('i', nbins)
      xmin_array = array('d', xmin)
      xmax_array = array('d', xmax)
      h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
      for i in range(0, dim):
        h.GetAxis(i).SetTitle(title[i])
        if i == 0 or i == 1:
          h.SetBinEdges(i, pt_bins)
        else:  # i == 2 or i == 3
          h.SetBinEdges(i, obs_bins)
      h.Sumw2()  # enables calculation of errors
      setattr(self, name, h)

    name = 'hResidual_JetPt_%s_R%s_%s%s' % (observable, jetR, obs_label, suffix) if \
      len(obs_label) else ('hResidual_JetPt_%s_R%s%s' % (observable, jetR, suffix))
    h = ROOT.TH3F(name, name, 20, 0, 200, 50, 0., 0.5, 200, -2., 2.)
    h.GetXaxis().SetTitle('#it{p}_{T,truth}^{ch jet}')
    h.GetYaxis().SetTitle(obs_name + '^{truth}')
    h.GetZaxis().SetTitle('#frac{%s^{det}-%s^{truth}}{%s^{truth}}' % \
      (obs_name, obs_name, obs_name))
    setattr(self, name, h)

    if not self.is_pp and self.fill_delta_obs:
      # Delta-observable histograms for studying background subtraction effects
      name = 'hDeltaObs_%s_emb_R%s_%s%s' % (observable, jetR, obs_label, suffix)
      h = None
      if observable == "ang":
        h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1, 1)
      elif observable == "mass":
        h = ROOT.TH2F(name, name, 300, 0, 300, 400, -100, 100)
      setattr(self, name, h)

  #---------------------------------------------------------------
  # Calculate the observable given a jet
  #---------------------------------------------------------------
  def calculate_observable(self, observable, jet, jet_groomed_lund,
      jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed):

    if observable == "ang":

      return fjext.lambda_beta_kappa(jet, jet_groomed_lund.pair(), obs_setting, 1, jetR) \
             if grooming_setting else fjext.lambda_beta_kappa(jet, obs_setting, 1, jetR)

    elif observable == "mass":

      return jet_groomed_lund.pair().m() if grooming_setting else jet.m()

    # Should not be any other observable
    raise ValueError("Observable %s not implemented" % observable)


  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  # Fill 2D histogram of (pt, obs)
  #---------------------------------------------------------------
  def fill_observable_histograms(self, observable, hname, jet, jet_groomed_lund,
      jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed):

    obs = self.calculate_observable(observable, jet, jet_groomed_lund, jetR,
        obs_setting, grooming_setting, obs_label, jet_pt_ungroomed)

    # Fill histograms
    name = hname.format(observable, obs_label).replace("__", "_")
    name = name[:-1] if name[-1] == "_" else name
    getattr(self, name).Fill(jet_pt_ungroomed, obs)


  #---------------------------------------------------------------
  # Fill matched jet histograms
  #---------------------------------------------------------------
  def fill_matched_jet_histograms(self, observable, jet_det, jet_det_groomed_lund,
      jet_truth, jet_truth_groomed_lund, jet_pp_det, jetR, obs_setting,
      grooming_setting, obs_label, jet_pt_det_ungroomed, jet_pt_truth_ungroomed,
      R_max, suffix, **kwargs):

    obs_det = self.calculate_observable(observable, jet_det, jet_det_groomed_lund,
        jetR, obs_setting, grooming_setting, obs_label, jet_pt_det_ungroomed)

    obs_tru = self.calculate_observable(observable, jet_truth, jet_truth_groomed_lund,
        jetR, obs_setting, grooming_setting, obs_label, jet_pt_truth_ungroomed)

    # Fill histograms
    self.fill_response(observable, jetR, jet_det.pt(), jet_truth.pt(),
                       obs_det, obs_tru, obs_label, R_max)


##################################################################
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Process MC')
  parser.add_argument('-f', '--inputFile', action='store',
                      type=str, metavar='inputFile',
                      default='AnalysisResults.root',
                      help='Path of ROOT file containing TTrees')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='config/analysis_config.yaml',
                      help="Path of config file for analysis")
  parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./TestOutput',
                      help='Output directory for output to be written to')

  # Parse the arguments
  args = parser.parse_args()

  print('Configuring...')
  print('inputFile: \'{0}\''.format(args.inputFile))
  print('configFile: \'{0}\''.format(args.configFile))
  print('ouputDir: \'{0}\"'.format(args.outputDir))

  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
    sys.exit(0)

  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = ProcessMC_ang(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_mc()
