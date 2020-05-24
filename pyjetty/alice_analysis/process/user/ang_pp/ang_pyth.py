#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of pythia parton-hadron-matched jets
  and build RMs for folding theory predictions
  
  Written by Ezra Lesser (elesser@berkeley.edu)
  Based on code by James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse
import math
import time
from tqdm import tqdm

# Data analysis and plotting
import uproot
import pandas
import numpy as np
import ROOT
import yaml
from array import *

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io_pyth, process_utils, process_base
from pyjetty.alice_analysis.process.user.ang_pp.helpers import lambda_beta_kappa, pT_bin

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class process_ang_data(process_base.ProcessBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_dir='', config_file='', output_dir='', debug_level=0, **kwargs):
    super(process_ang_data, self).__init__(input_dir, config_file, output_dir, debug_level, **kwargs)
    self.initialize_config()
    self.input_dir = input_dir

  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def process_ang_data(self):
    
    start_time = time.time()

    # Initialize configuration and histograms
    self.initializeHistograms()
    print(self)

    # Use IO helper class to convert ROOT TTree into a SeriesGroupBy object
    # of fastjet particles per event
    print('--- {} seconds ---'.format(time.time() - start_time))
    for jetR in self.jetR_list:
      inf_MPIon = self.input_dir + "PythiaResults_R%s.root" % str(jetR)
      inf_MPIoff = self.input_dir + "PythiaResults_R%s_MPIoff.root" % str(jetR)
      io = process_io_pyth.ProcessIO(input_file_MPIon=inf_MPIon, input_file_MPIoff=inf_MPIoff,
                                     betas=self.beta_list)
      self.df_ang_jets = io.load_data()
      print('--- {} seconds ---'.format(time.time() - start_time))

      # Fill histograms
      print('Fill response histograms for R=%s...' % jetR)
      self.analyzeJets(jetR)

    # Save ROOT histograms
    print('Save histograms...')
    process_base.ProcessBase.save_output_objects(self)

    print('--- {} seconds ---'.format(time.time() - start_time))

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
    
    # Call base class initialization
    process_base.ProcessBase.initialize_config(self)
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)

    # Set configuration for analysis
    self.jetR_list = config["jetR"]
    self.beta_list = config["betas"]

    self.n_pt_bins = config["n_pt_bins"]
    self.pt_limits = config["pt_limits"]
    self.pTbins = np.arange(self.pt_limits[0], self.pt_limits[1] + 1, 
                            (self.pt_limits[1] - self.pt_limits[0]) / self.n_pt_bins)
    self.n_lambda_bins = config["n_lambda_bins"]
    self.lambda_limits = config["lambda_limits"]
    self.n_rap_bins = config["n_rap_bins"]
    self.rap_limits = config["rap_limits"]

    # Configs for each jetR / beta
    self.config_dict = config["ang"]
  
  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initializeHistograms(self):
  
    for jetR in self.jetR_list:

      name = 'hJetPt_ch_R%s' % str(jetR).replace('.', '')
      h = ROOT.TH1F(name, name+';p_{T}^{jet, ch};#frac{dN}{dp_{T}^{jet, ch};', 300, 0, 300)
      setattr(self, name, h)

      name = 'hJetPt_p_R%s' % str(jetR).replace('.', '')
      h = ROOT.TH1F(name, name+';p_{T}^{jet, parton};#frac{dN}{dp_{T}^{jet, parton};', 300, 0, 300)
      setattr(self, name, h)

      name = 'hJetPtRes_R%s' % str(jetR).replace('.', '')
      h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
      h.GetXaxis().SetTitle('p_{T}^{jet, parton}')
      h.GetYaxis().SetTitle('#frac{p_{T}^{jet, parton}-p_{T}^{jet, ch}}{p_{T}^{jet, parton}}')
      setattr(self, name, h)

      name = 'hResponse_JetPt_R%s' % str(jetR).replace('.', '')
      h = ROOT.TH2F(name, name, 200, 0, 200, 200, 0, 200)
      h.GetXaxis().SetTitle('p_{T}^{jet, parton}')
      h.GetYaxis().SetTitle('p_{T}^{jet, ch}')
      setattr(self, name, h)

      for beta in self.beta_list:

        label = ("R%s_%s" % (str(jetR), str(beta))).replace('.', '')

        name = 'hResponse_ang_p_%s' % label
        h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
        h.GetXaxis().SetTitle('#lambda_{#beta=%s}^{parton}' % beta)
        h.GetYaxis().SetTitle('#lambda_{#beta=%s}^{ch}' % beta)
        setattr(self, name, h)

        name = 'hAng_JetPt_ch_%s' % label
        h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
        h.GetXaxis().SetTitle('p_{T}^{jet, ch}')
        h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{#beta=%s}^{ch}}' % str(beta))
        setattr(self, name, h)

        name = 'hAng_JetPt_p_%s' % label
        h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
        h.GetXaxis().SetTitle('p_{T}^{jet, parton}')
        h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{#beta=%s}^{parton}}' % str(beta))
        setattr(self, name, h)

        name = "hAngResidual_JetPt_%s" % label
        h = ROOT.TH2F(name, name, 300, 0, 300, 200, -2., 2.)
        h.GetXaxis().SetTitle('p_{T}^{jet, parton}')
        h.GetYaxis().SetTitle('#frac{#lambda_{#beta}^{jet, parton}-#lambda_{#beta}' + \
                              '^{jet, ch}}{#lambda_{#beta}^{jet, parton}}')
        setattr(self, name, h)

        # Create THn of response
        dim = 4
        title = ['p_{T}^{jet, parton}', 'p_{T}^{jet, ch}',
                 '#lambda_{#beta}^{parton}', '#lambda_{#beta}^{ch}']
        nbins = [10, 20, 100, 100]
        min_li = [0.,   0.,   0.,  0.]
        max_li = [100., 200., 1.0, 1.0]
        
        name = 'hResponse_JetPt_ang_%s' % label
        nbins = (nbins)
        xmin = (min_li)
        xmax = (max_li)
        nbins_array = array('i', nbins)
        xmin_array = array('d', xmin)
        xmax_array = array('d', xmax)
        h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
        for i in range(0, dim):
          h.GetAxis(i).SetTitle(title[i])
        setattr(self, name, h)

  #---------------------------------------------------------------
  # Main function to loop through and analyze events
  #---------------------------------------------------------------
  def analyzeJets(self, jetR):

    for index, ang_jets in tqdm(self.df_ang_jets.iterrows()):

      if abs(ang_jets["ch_eta"]) > 0.9:
        continue

      self.fillJetHistograms(ang_jets, jetR)

      for beta in self.beta_list:
        
        self.fillRMs(ang_jets, jetR, beta)

  #---------------------------------------------------------------
  # Fill jet histograms for each value of R
  #---------------------------------------------------------------
  def fillJetHistograms(self, ang_jets, jetR):

    label = "R%s" % str(jetR).replace('.', '')

    # Fill plots
    getattr(self, "hJetPt_ch_%s" % label).Fill(ang_jets["ch_pt"])
    getattr(self, "hJetPt_p_%s" % label).Fill(ang_jets["p_pt"])

    jet_pt_res = (ang_jets["p_pt"] - ang_jets["ch_pt"]) / ang_jets["p_pt"]
    getattr(self, "hJetPtRes_%s" % label).Fill(ang_jets["p_pt"], jet_pt_res)

    getattr(self, 'hResponse_JetPt_%s' % label).Fill(ang_jets["p_pt"], ang_jets["ch_pt"])


  #---------------------------------------------------------------
  # Fill jet response matrices for each value of R
  #---------------------------------------------------------------
  def fillRMs(self, ang_jets, jetR, beta):

    label = ("R%s_%s" % (str(jetR), str(beta))).replace('.', '')

    b = str(beta).replace('.', '')
    l_p = ang_jets["l_p_%s" % b]
    l_ch = ang_jets["l_ch_%s" % b]

    getattr(self, 'hResponse_ang_p_%s' % label).Fill(l_p, l_ch)
    getattr(self, 'hAng_JetPt_ch_%s' % label).Fill(ang_jets["ch_pt"], l_ch)
    getattr(self, 'hAng_JetPt_p_%s' % label).Fill(ang_jets["p_pt"], l_p)

    if l_p != 0:
      res = (l_p - l_ch) / l_p
      getattr(self, "hAngResidual_JetPt_%s" % label).Fill(ang_jets["p_pt"], res)

    x = array('d', [ang_jets["p_pt"], ang_jets["ch_pt"], l_p, l_ch])
    getattr(self, 'hResponse_JetPt_ang_%s' % label).Fill(x)

##################################################################
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Plot analysis histograms')
  parser.add_argument('-i', '--inputDir', action='store',
                      type=str, metavar='inputDir',
                      default='./',
                      help='Path of ROOT file containing TTrees')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='config/angularity.yaml',
                      help="Path of config file for jetscape analysis")
  parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./TestOutput',
                      help='Output directory for QA plots to be written to')
  
  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('inputDir: \'{0}\''.format(args.inputDir))
  print('configFile: \'{0}\''.format(args.configFile))
  print('ouputDir: \'{0}\"'.format(args.outputDir))
  print('----------------------------------------------------------------')
  
  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputDir):
    print('Directory \"{0}\" does not exist! Exiting!'.format(args.inputDir))
    sys.exit(0)
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = process_ang_data(input_dir=args.inputDir, config_file=args.configFile, 
                              output_dir=args.outputDir)
  analysis.process_ang_data()
