#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of MC track information
  and do jet-finding, and save response histograms.
  
  Most code adapted from rg analysis by James Mulligan (james.mulligan@berkeley.edu)
  Ezra Lesser (elesser@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse
import math
import time

# Data analysis and plotting
import uproot
import pandas
import numpy as np
from array import *
import ROOT
import yaml

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io, process_utils, jet_info, process_base
from pyjetty.alice_analysis.process.user.ang_pp.helpers import lambda_beta_kappa, pT_bin

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class process_ang_mc(process_base.ProcessBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
    super(process_ang_mc, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)

  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def process_ang_mc(self):
    
    start_time = time.time()
    
    # Initialize configuration
    self.initialize_config()
    
    # ------------------------------------------------------------------------
    
    # Use IO helper class to convert detector-level ROOT TTree into
    # a SeriesGroupBy object of fastjet particles per event
    print('--- {} seconds ---'.format(time.time() - start_time))
    io_det = process_io.ProcessIO(input_file=self.input_file, tree_dir="PWGHF_TreeCreator",
                                   track_tree_name="tree_Particle", event_tree_name="tree_event_char")
    df_fjparticles_det = io_det.load_data(self.reject_tracks_fraction)
    self.nEvents_det = len(df_fjparticles_det.index)
    self.nTracks_det = len(io_det.track_df.index)
    print('--- {} seconds ---'.format(time.time() - start_time))
    
    # ------------------------------------------------------------------------

    # Use IO helper class to convert truth-level ROOT TTree into
    # a SeriesGroupBy object of fastjet particles per event
    io_truth = process_io.ProcessIO(input_file=self.input_file, tree_dir="PWGHF_TreeCreator",
                                     track_tree_name="tree_Particle_gen", 
                                     event_tree_name="tree_event_char")
    df_fjparticles_truth = io_truth.load_data()
    self.nEvents_truth = len(df_fjparticles_truth.index)
    self.nTracks_truth = len(io_truth.track_df.index)
    print('--- {} seconds ---'.format(time.time() - start_time))
    
    # ------------------------------------------------------------------------

    # Now merge the two SeriesGroupBy to create a groupby df with [ev_id, run_number, fj_1, fj_2]
    # (Need a structure such that we can iterate event-by-event through both fj_1, fj_2 simultaneously)
    print('Merge det-level and truth-level into a single dataframe grouped by event...')
    self.df_fjparticles = pandas.concat([df_fjparticles_det, df_fjparticles_truth], axis=1)
    self.df_fjparticles.columns = ['fj_particles_det', 'fj_particles_truth']
    print('--- {} seconds ---'.format(time.time() - start_time))

    # ------------------------------------------------------------------------

    # Initialize histograms
    self.initializeHistograms()
    
    # Create constituent subtractor, if configured
    if self.do_constituent_subtraction:
      self.constituent_subtractor = \
        CEventSubtractor(max_distance=self.max_distance, alpha=self.alpha, max_eta=self.max_eta,
                         bge_rho_grid_size=self.bge_rho_grid_size, max_pt_correct=self.max_pt_correct,
                         ghost_area=self.ghost_area, distance_type=fjcontrib.ConstituentSubtractor.deltaR)
    
    print(self)
  
    # Find jets and fill histograms
    print('Find jets...')
    self.analyzeEvents()
    
    # Plot histograms
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
    
    self.jet_matching_distance = config['jet_matching_distance']
    self.reject_tracks_fraction = config['reject_tracks_fraction']

    # Set configuration for analysis
    self.jetR_list = config['jetR']
    self.beta_list = config['betas']

    self.n_pt_bins = config["n_pt_bins"]
    self.pt_limits = config["pt_limits"]
    self.pTbins = np.arange(self.pt_limits[0], self.pt_limits[1] + 1, 
                            (self.pt_limits[1] - self.pt_limits[0]) / self.n_pt_bins)
    self.n_lambda_bins = config['n_lambda_bins']
    self.lambda_limits = config['lambda_limits']

    self.sd_zcut = config["sd_zcut"]
    self.sd_beta = config["sd_beta"]

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initializeHistograms(self):

    self.hNevents = ROOT.TH1F('hNevents', 'hNevents', 2, -0.5, 1.5)
    self.hNevents.Fill(1, self.nEvents_det)

    self.hTrackEtaPhi = ROOT.TH2F('hTrackEtaPhi', 'hTrackEtaPhi_det', 200, -1., 1., 628, 0., 6.28)
    self.hTrackPt = ROOT.TH1F('hTrackPt', 'hTrackPt_det', 300, 0., 300.)

    for jetR in self.jetR_list:

      name = 'hJetPt_Truth_R%s' % str(jetR).replace('.', '')
      h = ROOT.TH1F(name, name+';p_{T,tru,ch jet};#frac{dN}{dp_{T,tru,ch jet};', 300, 0, 300)
      setattr(self, name, h)

      name = 'hJES_R%s' % str(jetR).replace('.', '')
      h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
      h.GetXaxis().SetTitle('p_{T,truth}')
      h.GetYaxis().SetTitle('#frac{p_{T,det,ch jet}-p_{T,tru,ch jet}}{p_{T,tru,ch jet}}')
      setattr(self, name, h)

      name = 'hDeltaR_All_R%s' % str(jetR).replace('.', '')
      h = ROOT.TH2F(name, name, 100, 0, 100, 100, 0., 2.)
      setattr(self, name, h)
      '''
      name = 'hZ_Truth_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 1.)
      setattr(self, name, h)
      
      name = 'hZ_Det_R{}'.format(jetR)
      h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 1.)
      setattr(self, name, h)
      '''
      name = 'hResponse_JetPt_R%s' % str(jetR).replace('.', '')
      h = ROOT.TH2F(name, name, 200, 0, 200, 200, 0, 200)
      h.GetXaxis().SetTitle('p_{T,det,ch jet}')
      h.GetYaxis().SetTitle('p_{T,tru,ch jet}')
      setattr(self, name, h)

      name = 'hJetPt_N_R%s' % str(jetR).replace('.', '')
      h = ROOT.TH2F(name, name, 200, 0, 200, 100, 0, 100)
      h.GetXaxis().SetTitle('p_{T,tru,ch jet}')
      h.GetYaxis().SetTitle('N_{constit}')
      setattr(self, name, h)

      for beta in self.beta_list:

        label = "R%s_%s" % (str(jetR), str(beta))

        name = 'hResponse_ang_%s' % label
        h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
        h.GetXaxis().SetTitle('#lambda_{%s,det}' % beta)
        h.GetYaxis().SetTitle('#lambda_{%s,tru}' % beta)
        setattr(self, name, h)

        name = 'hResponse_ang_%s_SD' % label
        h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
        h.GetXaxis().SetTitle('#lambda_{%s,det,SD}' % beta)
        h.GetYaxis().SetTitle('#lambda_{%s,tru,SD}' % beta)
        setattr(self, name, h)

        name = 'hAng_JetPt_det_%s' % label
        h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
        h.GetXaxis().SetTitle('p_{T,det}')
        h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{det,%s}}' % str(beta))
        setattr(self, name, h)

        name = 'hAng_JetPt_tru_%s' % label
        h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
        h.GetXaxis().SetTitle('p_{T,tru}')
        h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{tru,%s}}' % str(beta))
        setattr(self, name, h)

        '''
        for i, pTmin in list(enumerate(self.pTbins))[0:-1]:

          # Angularities, \lambda_{\beta}^{\kappa}
          pTmax = self.pTbins[i+1]
          name = ("hLambda_pT%i-%i_R%s_B%s_mcdet" % (pTmin, pTmax, jetR, beta)).replace('.', '')
          h = ROOT.TH1F(name, name, self.n_lambda_bins, 0, 1.0)
          h.GetXaxis().SetTitle('#lambda_{%s}' % beta)
          h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{%s}}' % beta)
          setattr(self, name, h)

          # Angularities with soft drop
          name = ("hLambda_pT%i-%i_R%s_B%s_mcdet_SD" % (pTmin, pTmax, jetR, beta)).replace('.', '')
          h = ROOT.TH1F(name, name, self.n_lambda_bins, 0, 1.0)
          h.GetXaxis().SetTitle('#lambda_{%s}' % beta)
          h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{%s}}' % beta)
          setattr(self, name, h)
        '''

        name = "hAngResidual_JetPt_%s" % label
        h = ROOT.TH2F(name, name, 300, 0, 300, 200, -2., 2.)
        h.GetXaxis().SetTitle('p_{T,truth}')
        h.GetYaxis().SetTitle('#frac{#lambda_{det}-#lambda_{truth}}{#lambda_{truth}}')
        setattr(self, name, h)

        # Create THn of response
        dim = 4;
        title = ['p_{T,det}', 'p_{T,truth}', '#lambda_{#beta,det}', '#lambda_{#beta,truth}']
        nbins = [20, 40, 200, 100]
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
  def analyzeEvents(self):
    
    # Fill det track histograms
    [ self.fillTrackHistograms(fj_particles_det) 
      for fj_particles_det in self.df_fjparticles['fj_particles_det'] ]
    
    fj.ClusterSequence.print_banner()
    print()
    
    for jetR in self.jetR_list:
      for beta in self.beta_list:
      
        # Set jet definition and a jet selector
        jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
        jet_selector_det = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
        jet_selector_truth_matched = fj.SelectorPtMin(5.0)
        print('jet definition is:', jet_def)
        print('jet selector for det-level is:', jet_selector_det,'\n')
        print('jet selector for truth-level matches is:', jet_selector_truth_matched,'\n')
        
        # Then can use list comprehension to iterate over the groupby and do jet-finding
        # simultaneously for fj_1 and fj_2 per event, so that I can match jets -- and fill histograms
        result = [ self.analyzeJets(fj_particles_det, fj_particles_truth, jet_def, jet_selector_det,
                                    jet_selector_truth_matched, beta)
                   for fj_particles_det, fj_particles_truth in
                   zip(self.df_fjparticles['fj_particles_det'], self.df_fjparticles['fj_particles_truth']) ]

  #---------------------------------------------------------------
  # Fill track histograms.
  #---------------------------------------------------------------
  def fillTrackHistograms(self, fj_particles_det):

    # Check that the entries exist appropriately
    # (need to check how this can happen -- but it is only a tiny fraction of events)
    if type(fj_particles_det) != fj.vectorPJ:
      return
    
    for track in fj_particles_det:
      self.hTrackEtaPhi.Fill(track.eta(), track.phi())
      self.hTrackPt.Fill(track.pt())

  #---------------------------------------------------------------
  # Analyze jets of a given event.
  # fj_particles is the list of fastjet pseudojets for a single fixed event.
  #---------------------------------------------------------------
  def analyzeJets(self, fj_particles_det, fj_particles_truth, jet_def, jet_selector_det,
                  jet_selector_truth_matched, beta):

    # Check that the entries exist appropriately
    # (need to check how this can happen -- but it is only a tiny fraction of events)
    if type(fj_particles_det) != fj.vectorPJ or type(fj_particles_truth) != fj.vectorPJ:
      print('fj_particles type mismatch -- skipping event')
      return
    
    # Perform constituent subtraction on det-level, if applicable
    if self.do_constituent_subtraction:
      fj_particles_det = self.constituent_subtractor.process_event(fj_particles_det)
      rho = self.constituent_subtractor.bge_rho.rho()

    # Do jet finding
    cs_det = fj.ClusterSequence(fj_particles_det, jet_def)
    jets_det = fj.sorted_by_pt(cs_det.inclusive_jets())
    jets_det_selected = jet_selector_det(jets_det)
    
    cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)
    jets_truth = fj.sorted_by_pt(cs_truth.inclusive_jets())
    jets_truth_selected = jet_selector_det(jets_truth)
    jets_truth_selected_matched = jet_selector_truth_matched(jets_truth)

    jetR = jet_def.R()
    
    # Fill det-level jet histograms (before matching)
    for jet_det in jets_det_selected:
      
      # Check additional acceptance criteria
      # skip event if not satisfied -- since first jet in event is highest pt
      if not self.utils.is_det_jet_accepted(jet_det):
        self.hNevents.Fill(0)
        return
      
      #self.fill_det_before_matching(jet_det, jetR)

    # Fill truth-level jet histograms (before matching)
    for jet_truth in jets_truth_selected:
      self.fill_truth_before_matching(jet_truth, jetR)

    # Loop through jets and set jet matching candidates for each jet in user_info
    hname = "hDeltaR_All_R%s" % str(jetR).replace('.', '')
    [ [ self.set_matching_candidates(jet_det, jet_truth, jetR, hname) 
        for jet_truth in jets_truth_selected_matched ] for jet_det in jets_det_selected ]
   
    # Loop through jets and set accepted matches
    hname = 'hJetMatchingQA_R%s' % str(jetR).replace('.', '')
    [ self.set_matches_pp(jet_det, hname) for jet_det in jets_det_selected ]
    
    # Loop through jets and fill matching histograms
    result = [ self.fill_matching_histograms(jet_det, jetR) for jet_det in jets_det_selected ]

    # Loop through jets and fill response if both det and truth jets are unique match
    for jetR in self.jetR_list:
      for beta in self.beta_list:

        result = [ self.fill_jet_matches(jet_det, jetR, beta) for jet_det in jets_det_selected ]

  #---------------------------------------------------------------
  # Loop through jets and fill response if both det and truth jets are unique match
  #---------------------------------------------------------------
  def fill_jet_matches(self, jet_det, jetR, beta):

    sd = fjcontrib.SoftDrop(self.sd_beta, self.sd_zcut, jetR)
    jet_def_recluster = fj.JetDefinition(fj.cambridge_algorithm, jetR)
    reclusterer = fjcontrib.Recluster(jet_def_recluster)
    sd.set_reclustering(True, reclusterer)
    if self.debug_level > 2:
      print('SoftDrop groomer is: {}'.format(sd.description()));

    # Check additional acceptance criteria
    # skip event if not satisfied -- since first jet in event is highest pt
    if not self.utils.is_det_jet_accepted(jet_det):
      return

    if jet_det.has_user_info():
      jet_truth = jet_det.python_info().match
      
      if jet_truth:
        self.fill_response_histograms(jet_det, jet_truth, sd, jetR, beta)

  '''
  #---------------------------------------------------------------
  # Loop through jets and store number of matching candidates in user_index
  # (In principle could also store matching candidate in user_info)
  #---------------------------------------------------------------
  def setNJetMatches(self, jets_det_selected, jets_truth_selected, jetR):

    # Reset user_index to 0
    for jet_det in jets_det_selected:
      jet_det.set_user_index(0)
    for jet_truth in jets_truth_selected:
      jet_truth.set_user_index(0)
    
    # Loop through jets and store number of matching candidates in user_index
    for jet_det in jets_det_selected:
      for jet_truth in jets_truth_selected:
        
        if jet_det.delta_R(jet_truth) < self.jet_matching_distance * jetR:
          
          jet_det.set_user_index(jet_det.user_index() + 1)
          jet_truth.set_user_index(jet_truth.user_index() + 1)
  '''
  #---------------------------------------------------------------
  # Fill truth jet histograms
  #---------------------------------------------------------------
  def fill_truth_before_matching(self, jet, jetR):

    getattr(self, 'hJetPt_Truth_R%s' % str(jetR).replace('.', '')).Fill(jet.pt())
    getattr(self, 'hJetPt_N_R%s' % str(jetR).replace('.', '')).Fill(jet.pt(), len(jet.constituents()))

    '''
    for constituent in jet.constituents():
      z = constituent.pt() / jet.pt()
      getattr(self, 'hZ_Truth_R{}'.format(jetR)).Fill(jet.pt(), z)
    '''


  '''
  #---------------------------------------------------------------
  # Fill det jet histograms
  #---------------------------------------------------------------
  def fill_det_before_matching(self, jet, jetR):
    
    for constituent in jet.constituents():
      z = constituent.pt() / jet.pt()
      getattr(self, 'hZ_Det_R{}'.format(jetR)).Fill(jet.pt(), z)
  '''

  #---------------------------------------------------------------
  # Loop through jets and fill matching histos
  #---------------------------------------------------------------
  def fill_matching_histograms(self, jet_det, jetR):
      
    if jet_det.has_user_info():
      jet_truth = jet_det.python_info().match
      
      if jet_truth:
        
        jet_pt_det_ungroomed = jet_det.pt()
        jet_pt_truth_ungroomed = jet_truth.pt()
        JES = (jet_pt_det_ungroomed - jet_pt_truth_ungroomed) / jet_pt_truth_ungroomed
        getattr(self, 'hJES_R%s' % str(jetR).replace('.', '')).Fill(jet_pt_truth_ungroomed, JES)

        getattr(self, 'hResponse_JetPt_R%s' % 
                str(jetR).replace('.', '')).Fill(jet_pt_det_ungroomed, jet_pt_truth_ungroomed)

  #---------------------------------------------------------------
  # Fill response histograms
  #---------------------------------------------------------------
  def fill_response_histograms(self, jet_det, jet_truth, sd, jetR, beta):

    # Ungroomed jet pT
    jet_pt_det = jet_det.pt()
    jet_pt_truth = jet_truth.pt()

    # just use kappa = 1 for now
    l_det = lambda_beta_kappa(jet_det, jetR, beta, 1) 
    l_tru = lambda_beta_kappa(jet_truth, jetR, beta, 1)

    # soft drop jet
    jet_sd_det = sd.result(jet_det)
    jet_sd_tru = sd.result(jet_truth)

    # lambda for soft drop jet
    l_sd_det = lambda_beta_kappa(jet_sd_det, jetR, beta, 1)
    l_sd_tru = lambda_beta_kappa(jet_sd_tru, jetR, beta, 1)

    label = "R%s_%s" % (str(jetR), str(beta))

    ''' Histograms per pT bin (currently unnecessary and cause clutter)
    (pTmin, pTmax) = pT_bin(jet_det.pt(), self.pTbins)
    if pTmin > -1e-3:  # pTmin will be -1 if not a valid bin
      getattr(self, ("hLambda_pT%i-%i_%s_mcdet" % (pTmin, pTmax, label)).replace('.', '')).Fill(l_det)

    (pTmin, pTmax) = pT_bin(jet_sd_det.pt(), self.pTbins)
    if pTmin > -1e-3:  # pTmin will be -1 if not a valid bin
      getattr(self, ("hLambda_pT%i-%i_%s_mcdet_SD" % (pTmin, pTmax, label)).replace('.', '')).Fill(l_sd_det)
    '''

    if l_tru != 0:
      lambda_resolution = (l_det - l_tru) / l_tru
      getattr(self, 'hAngResidual_JetPt_%s' % label).Fill(jet_pt_truth, lambda_resolution)

    # Observable plots
    getattr(self, 'hAng_JetPt_det_%s' % label).Fill(jet_pt_det, l_det)
    getattr(self, 'hAng_JetPt_tru_%s' % label).Fill(jet_pt_truth, l_tru)
    getattr(self, 'hResponse_ang_%s' % label).Fill(l_det, l_tru)
    getattr(self, 'hResponse_ang_%s_SD' % label).Fill(l_sd_det, l_sd_tru)

    x = ([jet_pt_det, jet_pt_truth, l_det, l_tru])
    x_array = array('d', x)
    getattr(self, 'hResponse_JetPt_ang_%s' % label).Fill(x_array)


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

  analysis = process_ang_mc(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_ang_mc()
