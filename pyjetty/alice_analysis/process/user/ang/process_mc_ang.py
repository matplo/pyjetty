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

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects_R(self, jetR):

      for observable in self.observable_list:
        # Should only be one: observable == "ang"
        if observable != "ang":
          raise ValueError("Observable %s is not implemented in this script" % observable)

        # Loop over subobservable (alpha value)
        for i in range(len(self.obs_settings[observable])):

          obs_setting = self.obs_settings[observable][i]
          grooming_setting = self.obs_grooming_settings[observable][i]
          obs_label = self.utils.obs_label(obs_setting, grooming_setting)

          # Create RM histograms
          if self.is_pp:
            self.create_ang_histograms(observable, jetR, obs_label)
          else:
            for R_max in self.max_distance:
              self.create_ang_histograms(observable, jetR, obs_label, R_max)
              if R_max == self.main_R_max:
                self.create_ang_histograms(observable, jetR, obs_label, str(R_max)+'_matched')

          if self.thermal_model:
            for R_max in self.max_distance:
              name = 'h_%s_JetPt_R%s_%s_Rmax%s' % (observable, jetR, obs_label, R_max)
              h = ROOT.TH2F(name, name, 200, 0, 200, 100, 0, 1)
              h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
              h.GetYaxis().SetTitle('#it{#lambda}_{#it{#alpha}}')
              setattr(self, name, h)

          name = 'h_%s_JetPt_Truth_R%s_%s' % (observable, jetR, obs_label)
          h = ROOT.TH2F(name, name, 40, 0, 200, 100, 0, 1)
          h.GetXaxis().SetTitle('#it{p}_{T,truth}^{ch jet}')
          h.GetYaxis().SetTitle('#it{#lambda}_{#it{#alpha},truth}')
          setattr(self, name, h)

  #---------------------------------------------------------------
  # Create angularity response histograms
  #---------------------------------------------------------------
  def create_ang_histograms(self, observable, jetR, obs_label, R_max = None):

    if R_max:
      suffix = '_Rmax' + str(R_max)
    else:
      suffix = ''

    # Create THn of response for ang
    if self.fill_RM_histograms:
      dim = 4;
      title = ['#it{p}_{T,det}^{ch jet}', '#it{p}_{T,truth}^{ch jet}',
               '#it{#lambda}_{#it{#alpha},det}', '#it{#lambda}_{#it{#alpha},truth}']
      pt_bins = array('d', list(range(5, 100, 5)) + list(range(100, 210, 10)))
      obs_bins = np.concatenate((np.linspace(0, 0.009, 10), np.linspace(0.01, 0.1, 19),
                                np.linspace(0.11, 0.8, 70)))
      nbins  = [len(pt_bins)-1, len(pt_bins)-1, len(obs_bins)-1, len(obs_bins)-1]
      min_li = [pt_bins[0],     pt_bins[0],     obs_bins[0],     obs_bins[0]    ]
      max_li = [pt_bins[-1],    pt_bins[-1],    obs_bins[-1],    obs_bins[-1]   ]

      name = 'hResponse_JetPt_%s_R%s_%s%s' % (observable, jetR, obs_label, suffix)
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

    name = 'hResidual_JetPt_%s_R%s_%s%s' % (observable, jetR, obs_label, suffix)
    h = ROOT.TH3F(name, name, 20, 0, 200, 50, 0., 0.5, 200, -2., 2.)
    h.GetXaxis().SetTitle('#it{p}_{T,truth}^{ch jet}')
    h.GetYaxis().SetTitle('#it{#lambda}_{#it{#alpha},truth}')
    h.GetZaxis().SetTitle(
      '#frac{#it{#lambda}_{#it{#alpha},det}-#it{#lambda}_{#it{#alpha},truth}}{#it{#lambda}_{#it{#alpha},truth}}')
    setattr(self, name, h)


  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  # Fill 2D histogram of (pt, obs)
  #---------------------------------------------------------------
  def fill_observable_histograms(self, hname, jet, jet_groomed_lund, jetR, obs_setting,
                                 grooming_setting, obs_label, jet_pt_ungroomed):

    # Calculate angularity and fill MPI scaling histogram
    ang = fjext.lambda_beta_kappa(jet, jet_groomed_lund.pair(), obs_setting, 1, jetR) \
          if grooming_setting else fjext.lambda_beta_kappa(jet, obs_setting, 1, jetR)

    # Fill histograms
    getattr(self, hname.format('ang', obs_label)).Fill(jet_pt_ungroomed, ang)


  #---------------------------------------------------------------
  # Fill matched jet histograms
  #---------------------------------------------------------------
  def fill_matched_jet_histograms(self, jet_det, jet_det_groomed_lund, jet_truth,
                                  jet_truth_groomed_lund, jet_pp_det, jetR,
                                  obs_setting, grooming_setting, obs_label,
                                  jet_pt_det_ungroomed, jet_pt_truth_ungroomed, R_max, suffix, **kwargs):

    # Calculate angularity and fill MPI scaling histogram
    ang_det = fjext.lambda_beta_kappa(jet_det, jet_det_groomed_lund.pair(), obs_setting, 1, jetR) \
              if grooming_setting else fjext.lambda_beta_kappa(jet_det, obs_setting, 1, jetR)
    ang_truth = fjext.lambda_beta_kappa(jet_truth, jet_truth_groomed_lund.pair(), obs_setting, 1, jetR) \
                if grooming_setting else fjext.lambda_beta_kappa(jet_truth, obs_setting, 1, jetR)

    # Fill histograms
    self.fill_response('ang', jetR, jet_pt_det_ungroomed, jet_pt_truth_ungroomed, ang_det, ang_truth, obs_label, R_max)


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
