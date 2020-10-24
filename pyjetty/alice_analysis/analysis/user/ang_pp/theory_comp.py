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
import ROOT
import yaml

from pyjetty.alice_analysis.analysis.user.substructure import run_analysis
from pyjetty.alice_analysis.analysis.user.james import plotting_utils

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)


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

    # Theory comparisons
    if 'fPythia' in config:
      self.fPythia_name = config['fPythia']
    if 'theory_dir' in config:
      self.theory_dir = config['theory_dir']
      self.theory_beta = config['theory_beta']
      self.theory_pt_bins = config['theory_pt_bins']
      self.theory_response_files = [ROOT.TFile(f, 'READ') for f in config['response_files']]
      self.theory_response_labels = config['response_labels']
      self.rebin_theory_response = config['rebin_theory_response']
      self.output_dir_theory = os.path.join(self.output_dir, self.observable, 'theory_response')
      self.Lambda = 1  # GeV -- This variable changes the NP vs P region of theory plots
      # Load the RooUnfold library
      ROOT.gSystem.Load(config['roounfold_path'])
      self.do_theory = True
    else:
      self.do_theory = False

    if self.do_theory:
      print("Loading response matrix for folding theory predictions...")
      self.load_theory_response()
      print("Loading theory histograms...")
      self.load_theory_histograms()


  #---------------------------------------------------------------
  # Load 4D response matrices used for folding the theory predictions
  #---------------------------------------------------------------
  def load_theory_response(self):

    # Check to see if Roounfold file already exists
    if not os.path.exists(self.output_dir_theory):
      os.makedirs(self.output_dir_theory)
    roounfold_filename = os.path.join(self.output_dir_theory, 'fRoounfold.root')
    roounfold_exists = os.path.exists(roounfold_filename)

    for jetR in self.jetR_list:
      for beta in self.theory_beta:
        label = "R%s_%s" % (str(jetR).replace('.', ''), str(beta).replace('.', ''))

        for ri, response in enumerate(self.theory_response_files):
          # Load charged hadron level folding response matrix
          name_ch = "hResponse_JetPt_%s_ch_%sScaled" % (self.observable, label)
          thn_ch = response.Get(name_ch)
          name_ch = "hResponse_theory_ch_%s" % label
          setattr(self, '%s_%i' % (name_ch, ri), thn_ch)

          # Load hadron-level folding response matrix (for comparison histograms)
          name_h = "hResponse_JetPt_%s_h_%sScaled" % (self.observable, label)
          thn_h = response.Get(name_h)
          name_h = "hResponse_theory_h_%s" % label
          setattr(self, '%s_%i' % (name_h, ri), thn_h)

          # Create Roounfold object
          name_roounfold_h = '%s_Roounfold_%i' % (name_h, ri)
          name_roounfold_ch = '%s_Roounfold_%i' % (name_ch, ri)

          if roounfold_exists and not self.rebin_theory_response:
            fRoo = ROOT.TFile(roounfold_filename, 'READ')
            roounfold_response_ch = fRoo.Get(name_roounfold_ch)
            roounfold_response_h = fRoo.Get(name_roounfold_h)
            fRoo.Close()

          elif self.rebin_theory_response:  # Generated theory folding matrix needs rebinning
            # Response axes: ['p_{T}^{ch jet}', 'p_{T}^{jet, parton}', 
            #                 '#lambda_{#beta}^{ch}', '#lambda_{#beta}^{parton}']
            # as compared to the usual
            #      ['p_{T,det}', 'p_{T,truth}', '#lambda_{#beta,det}', '#lambda_{#beta,truth}']
            det_pt_bin_array = array('d', range(10, 110, 10))
            tru_pt_bin_array = array('d', range(10, 160, 15))
            det_obs_bin_array = array('d', np.linspace(0, 1, 101, endpoint=True))
            tru_obs_bin_array = array('d', np.linspace(-0.005, 1.005, 102, endpoint=True))
            self.utils.rebin_response(
              roounfold_filename, thn_ch, '%s_Rebinned_%i' % (name_ch, ri), name_roounfold_ch, label,
              len(det_pt_bin_array)-1, det_pt_bin_array, len(det_obs_bin_array)-1, det_obs_bin_array,
              len(tru_pt_bin_array)-1, tru_pt_bin_array, len(tru_obs_bin_array)-1, tru_obs_bin_array,
              self.observable)
            self.utils.rebin_response(
              roounfold_filename, thn_h, '%s_Rebinned_%i' % (name_h, ri), name_roounfold_h, label, 
              len(det_pt_bin_array)-1, det_pt_bin_array, len(det_obs_bin_array)-1, det_obs_bin_array,
              len(tru_pt_bin_array)-1, tru_pt_bin_array, len(tru_obs_bin_array)-1, tru_obs_bin_array,
              self.observable)

            f_resp = ROOT.TFile(self.theory_rebinned_response_file, 'READ')
            roounfold_response_ch = f_resp.Get(name_roounfold_ch)
            roounfold_response_h = f_resp.Get(name_roounfold_h)
            f_resp.Close()
          

          else:   # Theory folding matrix already has correct binning
            hist_p_jet = thn_ch.Projection(3, 1)
            hist_p_jet.SetName('hist_p_jet_%s_%i' % (label, ri))
            hist_h_jet = thn_h.Projection(2, 0)
            hist_h_jet.SetName('hist_h_jet_%s_%i' % (label, ri))
            hist_ch_jet = thn_ch.Projection(2, 0)
            hist_ch_jet.SetName('hist_ch_jet_%s_%i' % (label, ri))
            roounfold_response_ch = ROOT.RooUnfoldResponse(
              hist_ch_jet, hist_p_jet, name_roounfold_ch, name_roounfold_ch)
            roounfold_response_h = ROOT.RooUnfoldResponse(
              hist_h_jet, hist_p_jet, name_roounfold_h, name_roounfold_h)

            for bin_0 in range(1, thn_ch.GetAxis(0).GetNbins() + 1):
              if bin_0 % 5 == 0:
                print('{} / {}'.format(bin_0, thn_ch.GetAxis(0).GetNbins() + 1))
              pt_det = thn_ch.GetAxis(0).GetBinCenter(bin_0)
              for bin_1 in range(1, thn_ch.GetAxis(1).GetNbins() + 1):
                pt_true = thn_ch.GetAxis(1).GetBinCenter(bin_1)
                for bin_2 in range(0, thn_ch.GetAxis(2).GetNbins() + 1):
                  for bin_3 in range(0, thn_ch.GetAxis(3).GetNbins() + 1):
                    obs_det = thn_ch.GetAxis(2).GetBinCenter(bin_2)
                    obs_true = thn_ch.GetAxis(3).GetBinCenter(bin_3)
                    # Get content of original THn bin
                    x_list = (pt_det, pt_true, obs_det, obs_true)
                    x = array('d', x_list)
                    global_bin = thn_ch.GetBin(x)
                    content = thn_ch.GetBinContent(global_bin)
                    roounfold_response_ch.Fill(pt_det, obs_det, pt_true, obs_true, content)

                    # Assume that thn_h and thn_ch have same binning
                    obs_det = thn_h.GetAxis(2).GetBinCenter(bin_2)
                    obs_true = thn_h.GetAxis(3).GetBinCenter(bin_3)
                    # Get content of original THn bin
                    x_list = (pt_det, pt_true, obs_det, obs_true)
                    x = array('d', x_list)
                    global_bin = thn_h.GetBin(x)
                    content = thn_h.GetBinContent(global_bin)
                    roounfold_response_h.Fill(pt_det, obs_det, pt_true, obs_true, content)

            fRoo = ROOT.TFile(roounfold_filename, 'UPDATE')
            roounfold_response_ch.Write()
            roounfold_response_h.Write()
            fRoo.Close()

          setattr(self, name_roounfold_ch, roounfold_response_ch)
          setattr(self, name_roounfold_h, roounfold_response_h)


  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #---------------------------------------------------------------
  def load_theory_histograms(self):

    # Set central value to exponential distribution
    exp_test = False

    # Require that hard scale and jet scale are varied together
    scale_req = False

    # Create histogram for each value of R and beta
    for jetR in self.jetR_list:
      for beta in self.theory_beta:   # beta value
        label = "R%s_%s" % (str(jetR).replace('.', ''), str(beta).replace('.', ''))

        name_cent = "theory_cent_%s_%s_parton" % (self.observable, label)
        name_min = "theory_min_%s_%s_parton" % (self.observable, label)
        hist_min = ROOT.TH2D(name_min, name_min, len(self.theory_pt_bins) - 1, 
            self.theory_pt_bins[0], self.theory_pt_bins[-1], 101, -0.005, 1.005)
        name_max = "theory_max_%s_%s_parton" % (self.observable, label)
        hist_max = ROOT.TH2D(name_max, name_max, len(self.theory_pt_bins) - 1, 
            self.theory_pt_bins[0], self.theory_pt_bins[-1], 101, -0.005, 1.005)

        parton_hists = ( ([], [], []), ([], [], []), ([], [], []) )

        for l in range(0, 3):
          for m in range(0, 3):
            for n in range(0, 3):

              name_hist = "theory_%i%i%i_%s_%s_parton" % (l, m, n, self.observable, label)
              hist = ROOT.TH2D(name_hist, name_hist, len(self.theory_pt_bins) - 1, 
                               self.theory_pt_bins[0], self.theory_pt_bins[-1], 101, -0.005, 1.005)

              if (scale_req and m != n) or (0 in (l, m, n) and 2 in (l, m, n)):
                parton_hists[l][m].append(None)
                continue

              # Loop through each pT-bin
              for i, pt_min in enumerate(self.theory_pt_bins[0:-1]):
                pt_max = self.theory_pt_bins[i+1]
                th_dir = os.path.join(
                  self.theory_dir, "R%s" % str(jetR).replace('.', ''), 
                  "pT%s_%s" % (pt_min, pt_max), "beta%s" % str(beta).replace('.', 'p'))

                # Load theory predictions for lambda values
                with open(os.path.join(th_dir, "%i%i%i.txt" % (l, m, n))) as f:
                  lines = f.read().split('\n')
                  val_li = [line.split()[1] for line in lines]

                if exp_test:
                  val_li = [1/(x+0.4+l) for x in np.linspace(0, 0.5, 51, True)] + \
                            [0 for x in np.linspace(0.51, 1, 50, True)] 
                  #np.exp(np.linspace(0, 1, 101, True))
                  #val = np.concatenate((np.full(51, 1), np.full(50, 0)))
                  #val = [0.6 - x for x in np.linspace(0, 1, 101, True)]

                for j, val in enumerate(val_li):
                  hist.SetBinContent(i+1, j+1, float(val))
                  if l == m == n == 0:
                    hist_min.SetBinContent(i+1, j+1, float(val))
                    hist_max.SetBinContent(i+1, j+1, float(val))
                  elif float(val) < hist_min.GetBinContent(i+1, j+1):
                    hist_min.SetBinContent(i+1, j+1, float(val))
                  elif float(val) > hist_max.GetBinContent(i+1, j+1):
                    hist_max.SetBinContent(i+1, j+1, float(val))

              parton_hists[l][m].append(hist)

        setattr(self, name_cent, parton_hists[1][1][1])
        setattr(self, name_min, hist_min)
        setattr(self, name_max, hist_max)

        print("Folding theory predictions...")
        self.fold_theory(jetR, beta, parton_hists, scale_req)


  #----------------------------------------------------------------------
  # Fold theoretical predictions
  #----------------------------------------------------------------------
  def fold_theory(self, jetR, beta, parton_hists, scale_req):

    label = "R%s_%s" % (str(jetR).replace('.', ''), str(beta).replace('.', ''))

    for ri, response in enumerate(self.theory_response_files):

      # Load parton-to-charged-hadron response matrix
      response_name_ch = "hResponse_theory_ch_%s_Roounfold_%i" % (label, ri)
      response_ch = getattr(self, response_name_ch)
      response_name_h = "hResponse_theory_h_%s_Roounfold_%i" % (label, ri)
      response_h = getattr(self, response_name_h)

      folded_ch_hists = ( ([], [], []), ([], [], []), ([], [], []) )
      folded_h_hists = ( ([], [], []), ([], [], []), ([], [], []) )

      for i in range(0, 3):
        for j in range(0, 3):
          for k in range(0, 3):

            if (scale_req and j != k) or (0 in (i, j, k) and 2 in (i, j, k)):
              folded_h_hists[i][j].append(None)
              folded_ch_hists[i][j].append(None)
              continue

            # Fold theory predictions
            h_folded_ch = response_ch.ApplyToTruth(parton_hists[i][j][k])
            h_folded_h = response_h.ApplyToTruth(parton_hists[i][j][k])

            name_ch = "theory_%i%i%i_%s_%s_ch_%i" % (i, j, k, self.observable, label, ri)
            name_h = "theory_%i%i%i_%s_%s_h_%i" % (i, j, k, self.observable, label, ri)

            h_folded_ch.SetNameTitle(name_ch, name_ch)
            h_folded_h.SetNameTitle(name_h, name_h)

            folded_ch_hists[i][j].append(h_folded_ch)
            folded_h_hists[i][j].append(h_folded_h)

      print("Scaling theory predictions for MPI effects for %s..." % self.theory_response_labels[ri])
      self.mpi_scale_theory(jetR, beta, ri, response, folded_ch_hists, folded_h_hists, scale_req)


  #----------------------------------------------------------------------
  # Fold theoretical predictions
  #----------------------------------------------------------------------
  def mpi_scale_theory(self, jetR, beta, ri, response, folded_ch_hists, folded_h_hists, scale_req):

    label = "R%s_%s" % (str(jetR).replace('.', ''), str(beta).replace('.', ''))

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

    # Ensure that the scaling and theory histograms have the same binning
    y_rebin_num = h_mpi_off.GetNbinsY() / folded_ch_hists[1][1][1].GetNbinsY()
    if y_rebin_num < 1 or abs(y_rebin_num - int(y_rebin_num)) > 1e-5:
      print("ERROR: histograms for MPI scaling from response file have insufficienctly binning.")
      print("       %i versus even multiple of  %i bins required" % \
            (h_mpi_off.GetNbinsY(), folded_ch_hists[1][1][1].GetNbinsY()))
      exit(1)
    h_mpi_off.RebinY(int(y_rebin_num))
    h_mpi_on.RebinY(int(y_rebin_num))
    if h_mpi_off.GetNbinsY() != folded_ch_hists[1][1][1].GetNbinsY():
      print("ERROR: rebinning histograms for MPI scaling failed.")
      exit(1)

    mpi_bin_edges = [h_mpi_on.GetXaxis().GetBinLowEdge(i+1)
                     for i in range(h_mpi_on.GetNbinsX()+1)]
    mpi_bin_edge_indices = [mpi_bin_edges.index(pt) for pt in self.pt_bins_reported]
    theory_pt_bin_edge_indices = [self.theory_pt_bins.index(pt) for pt in self.pt_bins_reported]
    # Loop through each reported pT-bin
    for index, i in list(enumerate(theory_pt_bin_edge_indices))[0:-1]:
      j = theory_pt_bin_edge_indices[index+1]
      pt_min = self.theory_pt_bins[i]
      pt_max = self.theory_pt_bins[j]

      # Note: bins in ROOT are 1-indexed (0 bin is underflow). Also last bin is inclusive
      h_mpi_off_proj = h_mpi_off.ProjectionY(
        '%s_py_%s_%i' % (name_mpi_off, str(i), ri), mpi_bin_edge_indices[index]+1, mpi_bin_edge_indices[index+1])
      h_mpi_on_proj = h_mpi_on.ProjectionY(
        '%s_py_%s_%i' % (name_mpi_on, str(i), ri), mpi_bin_edge_indices[index]+1, mpi_bin_edge_indices[index+1])

      # Create ratio plot for scaling
      h_mpi_ratio = h_mpi_on_proj.Clone()
      title = 'h_theory_mpi_scaling_%s_%s_%i' % (label, str(i), ri)
      h_mpi_ratio.SetNameTitle(title, title)
      h_mpi_ratio.Divide(h_mpi_off_proj)
      h_mpi_ratio.SetDirectory(0)
      setattr(self, 'h_mpi_ratio_%s_PtBin%i-%i_%i' % (label, pt_min, pt_max, ri), h_mpi_ratio)

      pt_label = '_PtBin'+str(pt_min)+'-'+str(pt_max)

      h_cent_bin = h_cent.ProjectionY("%s_parton%s" % (name_cent, pt_label), i+1, j)
      h_min_bin = h_min.ProjectionY("%s_parton%s" % (name_min, pt_label), i+1, j)
      h_max_bin = h_max.ProjectionY("%s_parton%s" % (name_max, pt_label), i+1, j)

      # Keep normalization by dividing by the number of bins added in projection
      h_cent_bin.Scale(1/(j-i))
      h_min_bin.Scale(1/(j-i))
      h_max_bin.Scale(1/(j-i))

      # Set new names and titles to prevent issues with saving histograms
      h_cent_bin.SetNameTitle("%s_parton%s" % (name_cent, pt_label),
                              "%s_parton%s" % (name_cent, pt_label))
      h_min_bin.SetNameTitle("%s_parton%s" % (name_min, pt_label),
                             "%s_parton%s" % (name_min, pt_label))
      h_max_bin.SetNameTitle("%s_parton%s" % (name_max, pt_label),
                             "%s_parton%s" % (name_max, pt_label))

      # Have to shift/average bins since raw theory calculation are points
      title = 'h_cent_parton_shifted_%s%s' % (label, pt_label)
      xbins = array('d', np.linspace(0, 1, 101, endpoint=True))
      h_cent_p_bin = ROOT.TH1D(title, title, len(xbins)-1, xbins)
      for bi in range(1, h_cent_bin.GetNbinsX()):
        h_cent_p_bin.SetBinContent(
          bi, (h_cent_bin.GetBinContent(bi) + h_cent_bin.GetBinContent(bi+1)) / 2)
        h_cent_p_bin.SetBinError(bi, 0)

      # Initialize h and ch histograms
      h_cent_ch_bin = None
      h_min_ch_bin = None
      h_max_ch_bin = None

      h_cent_h_bin = None
      h_min_h_bin = None
      h_max_h_bin = None

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

            h_folded_h_bin.Scale(1/(j-i))
            h_folded_ch_bin.Scale(1/(j-i))

            # Scale all by integral of central prediction values
            # TODO: See if there's a better way to do this... shouldn't have to scale?

            h_folded_h_bin.Scale(1/h_folded_h_bin.Integral("width"))
            h_folded_ch_bin.Scale(1/h_folded_ch_bin.Integral("width"))

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
      h_cent_bin.SetDirectory(0)
      h_min_bin.SetDirectory(0)
      h_max_bin.SetDirectory(0)

      h_cent_h_bin.SetDirectory(0)
      h_min_h_bin.SetDirectory(0)
      h_max_h_bin.SetDirectory(0)

      h_cent_ch_bin.SetDirectory(0)
      h_min_ch_bin.SetDirectory(0)
      h_max_ch_bin.SetDirectory(0)

      # Save as attributes for later access (plotting)
      setattr(self, "%s_parton%s" % (name_cent, pt_label), h_cent_bin)
      setattr(self, "%s_parton%s" % (name_min, pt_label), h_min_bin)
      setattr(self, "%s_parton%s" % (name_max, pt_label), h_max_bin)

      setattr(self, "%s_h%s_%i" % (name_cent, pt_label, ri), h_cent_h_bin)
      setattr(self, "%s_h%s_%i" % (name_min, pt_label, ri), h_min_h_bin)
      setattr(self, "%s_h%s_%i" % (name_max, pt_label, ri), h_max_h_bin)

      setattr(self, "%s_ch%s_%i" % (name_cent, pt_label, ri), h_cent_ch_bin)
      setattr(self, "%s_ch%s_%i" % (name_min, pt_label, ri), h_min_ch_bin)
      setattr(self, "%s_ch%s_%i" % (name_max, pt_label, ri), h_max_ch_bin)


  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #---------------------------------------------------------------
  def plot_single_result(self, jetR, obs_label, obs_setting, grooming_setting):
    print('Plotting each individual result...')

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
    
    # Initialize performance plotting class
    if self.do_plot_performance:
      self.plotting_utils = plotting_utils.PlottingUtils(
        self.observable, self.is_pp, self.main_data, self.main_response,
        self.output_dir_performance, self.figure_approval_status)
      
    # Create output subdirectories
    self.create_output_subdir(self.output_dir_performance, 'jet')
    self.create_output_subdir(self.output_dir_performance, 'resolution')
    self.create_output_subdir(self.output_dir_performance, 'residual')
    self.create_output_subdir(self.output_dir_performance, 'residual_relative')
    self.create_output_subdir(self.output_dir_performance, 'mc_projections_det')
    self.create_output_subdir(self.output_dir_performance, 'mc_projections_truth')
    self.create_output_subdir(self.output_dir_performance, 'statistics')
    self.create_output_subdir(self.output_dir_performance, 'lund')
    if not self.is_pp:
      self.create_output_subdir(self.output_dir_performance, 'prong_matching_fraction_pt')
      self.create_output_subdir(self.output_dir_performance, 'prong_matching_fraction_ptdet')
      self.create_output_subdir(self.output_dir_performance, 'prong_matching_deltaR')
      self.create_output_subdir(self.output_dir_performance, 'prong_matching_deltaZ')
      self.create_output_subdir(self.output_dir_performance, 'prong_matching_correlation')
      self.create_output_subdir(self.output_dir_performance, 'prong_matching_N_z')
    
    # Generate performance plots
    for jetR in self.jetR_list:
  
      # Plot some subobservable-independent performance plots
      self.plotting_utils.plot_DeltaR(jetR, self.jet_matching_distance)
      self.plotting_utils.plot_JES(jetR)
      self.plotting_utils.plot_JES_proj(jetR, self.pt_bins_reported)
      self.plotting_utils.plotJER(jetR, self.utils.obs_label(self.obs_settings[0], 
                                                             self.grooming_settings[0]))
      self.plotting_utils.plot_jet_reco_efficiency(jetR, self.utils.obs_label(
        self.obs_settings[0], self.grooming_settings[0]))
      
      # Plot subobservable-dependent performance plots
      for i, _ in enumerate(self.obs_subconfig_list):

        obs_setting = self.obs_settings[i]
        grooming_setting = self.grooming_settings[i]
        obs_label = self.utils.obs_label(obs_setting, grooming_setting)
    
        self.plotting_utils.plot_obs_resolution(jetR, obs_label, self.xtitle, self.pt_bins_reported)
        self.plotting_utils.plot_obs_residual(jetR, obs_label, self.xtitle, self.pt_bins_reported)
        self.plotting_utils.plot_obs_residual(
          jetR, obs_label, self.xtitle, self.pt_bins_reported, relative=True)
        self.plotting_utils.plot_obs_projections(jetR, obs_label, obs_setting, grooming_setting,
                                                 self.xtitle, self.pt_bins_reported)
        
        if grooming_setting and self.observable != 'jet_axis':
          self.plotting_utils.plot_lund_plane(jetR, obs_label, grooming_setting)

      # Plot prong matching histograms
      if not self.is_pp:
        self.prong_match_threshold = 0.5
        min_pt = 80.
        max_pt = 100.
        prong_list = ['leading', 'subleading']
        match_list = ['leading', 'subleading', 'groomed', 'ungroomed', 'outside']
        for i, overlay_list in enumerate(self.plot_overlay_list):
          for prong in prong_list:
            for match in match_list:

              hname = 'hProngMatching_{}_{}_JetPt_R{}'.format(prong, match, jetR)
              self.plotting_utils.plot_prong_matching(
                i, jetR, hname, self.obs_subconfig_list, self.obs_settings,
                self.grooming_settings, overlay_list, self.prong_match_threshold)

              self.plotting_utils.plot_prong_matching_delta(
                i, jetR, hname, self.obs_subconfig_list, self.obs_settings,
                self.grooming_settings, overlay_list, self.prong_match_threshold,
                min_pt, max_pt, plot_deltaz=False)

              hname = 'hProngMatching_{}_{}_JetPtDet_R{}'.format(prong, match, jetR)
              self.plotting_utils.plot_prong_matching(
                i, jetR, hname, self.obs_subconfig_list, self.obs_settings,
                self.grooming_settings, overlay_list, self.prong_match_threshold)

              if 'subleading' in prong:
                hname = 'hProngMatching_{}_{}_JetPtZ_R{}'.format(prong, match, jetR)
                self.plotting_utils.plot_prong_matching_delta(
                  i, jetR, hname, self.obs_subconfig_list, self.obs_settings,
                  self.grooming_settings, overlay_list, self.prong_match_threshold,
                  min_pt, max_pt, plot_deltaz=True)

          hname = 'hProngMatching_subleading-leading_correlation_JetPtDet_R{}'.format(jetR)
          self.plotting_utils.plot_prong_matching_correlation(
            jetR, hname, self.obs_subconfig_list, self.obs_settings, self.grooming_settings,
            overlay_list, self.prong_match_threshold)

        # Plot subobservable-dependent plots
        for i, _ in enumerate(self.obs_subconfig_list):
          obs_setting = self.obs_settings[i]
          grooming_setting = self.grooming_settings[i]
          obs_label = self.utils.obs_label(obs_setting, grooming_setting)
          self.plotting_utils.plot_prong_N_vs_z(jetR, obs_label, 'tagged')
          self.plotting_utils.plot_prong_N_vs_z(jetR, obs_label, 'untagged')

  #----------------------------------------------------------------------
  def plot_final_result(self, jetR, obs_label, obs_setting, grooming_setting):
    print('Plot final results for {}: R = {}, {}'.format(self.observable, jetR, obs_label))

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    # Loop through pt slices, and plot final result for each 1D theta_g distribution
    for i in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[i]
      max_pt_truth = self.pt_bins_reported[i+1]
      maxbin = self.obs_max_bins(obs_label)[i]

      self.plot_observable(jetR, obs_label, obs_setting, grooming_setting,
                           min_pt_truth, max_pt_truth, maxbin, plot_pythia=True)

      if self.do_theory and float(obs_label) in self.theory_beta:
        self.plot_observable(jetR, obs_label, obs_setting, grooming_setting,
                             min_pt_truth, max_pt_truth, maxbin, plot_pythia=False, plot_theory=True)
        self.plot_theory_ratios(jetR, obs_label, obs_setting, grooming_setting,
                                min_pt_truth, max_pt_truth, maxbin)
        self.plot_theory_response(jetR, obs_label, obs_setting, grooming_setting,
                                  min_pt_truth, max_pt_truth, maxbin)

      if min_pt_truth == 40 and (jetR == 0.2 or jetR == 0.4):
        # Only want to compare to girth with \beta=1
        if obs_label == '1':
          self.plot_obs_comp(jetR, obs_label, obs_setting, grooming_setting,
                             min_pt_truth, max_pt_truth, maxbin)

  #----------------------------------------------------------------------
  def plot_observable(self, jetR, obs_label, obs_setting, grooming_setting,
                      min_pt_truth, max_pt_truth, maxbin, plot_pythia=False, plot_theory=False):

    # For theory plots, whether or not to show original parton-level predictions
    show_parton_theory = True

    # For theory plots, whether or not to show the NP / P region
    show_np_region = True

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
    myPad.Draw()
    myPad.cd()
    
    xtitle = getattr(self, 'xtitle')
    ytitle = getattr(self, 'ytitle')
    color = 1   # black for data
    
    # Get histograms
    name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label,
                                              min_pt_truth, max_pt_truth)
    if grooming_setting:
      fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
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
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle(ytitle)
    myBlankHisto.SetMaximum(2.5*h.GetBinContent(int(0.5*h.GetNbinsX())))
    if self.observable == 'subjet_z' or self.observable == 'jet_axis':
      myBlankHisto.SetMaximum(1.7*h.GetMaximum())
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")

    if plot_theory:
      label = "R%s_%s" % (str(jetR).replace('.', ''), str(obs_label).replace('.', ''))

      if show_parton_theory:
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

        n = hcent_p.GetNbinsX()
        x = array('d', [hcent_p.GetXaxis().GetBinCenter(i) for i in range(1, hcent_p.GetNbinsX()+1)])
        y = array('d', [hcent_p.GetBinContent(i) for i in range(1, hcent_p.GetNbinsX()+1)])
        xerrup = array('d', [(x[i+1] - x[i]) / 2. for i in range(n-1)] + [0])
        xerrdn = array('d', [0] + [(x[i+1] - x[i]) / 2. for i in range(n-1)])
        yerrup = array('d', [hmax_p.GetBinContent(i)-y[i-1] for i in range(1, hmax_p.GetNbinsX()+1)])
        yerrdn = array('d', [y[i-1]-hmin_p.GetBinContent(i) for i in range(1, hmin_p.GetNbinsX()+1)])

        color = self.ColorArray[4]
        if show_np_region:
          # P vs NP cutoff point: lambda_beta ~ Lambda / (pT * R) -- use avg value of pT for the bin.
          # Formula assumes that jet pT xsec falls like pT^(-5.5)
          formula_pt = (4.5/3.5)*(min_pt_truth**-3.5 - max_pt_truth**-3.5)/(min_pt_truth**-4.5 - max_pt_truth**-4.5)
          lambda_np_cutoff = round(self.Lambda / (formula_pt * jetR), 2)
          if lambda_np_cutoff > 1:
            lambda_np_cutoff = 1
          index_np_cutoff = [round(val, 2) for val in x].index(lambda_np_cutoff)

          #+1 to include lambda in NP
          h_parton_theory_np = ROOT.TGraphAsymmErrors(
            index_np_cutoff+1, x[:index_np_cutoff+1], y[:index_np_cutoff+1], xerrdn[:index_np_cutoff+1],
            xerrup[:index_np_cutoff+1], yerrdn[:index_np_cutoff+1], yerrup[:index_np_cutoff+1])
          h_parton_theory_p = ROOT.TGraphAsymmErrors(
            n-index_np_cutoff, x[index_np_cutoff:], y[index_np_cutoff:], xerrdn[index_np_cutoff:],
            xerrup[index_np_cutoff:], yerrdn[index_np_cutoff:], yerrup[index_np_cutoff:])

          h_parton_theory_np.SetFillColorAlpha(color, 0.5)
          h_parton_theory_p.SetFillColorAlpha(color, 0.25)
          h_parton_theory_np.SetFillStyle(3002)
          h_parton_theory_np.Draw('3 same')
          h_parton_theory_p.Draw('3 same')

          # Split central parton curve in NP and P regions
          hcent_p_np = ROOT.TH1F(
            name_cent+"_nonpert", name_cent+"_nonpert", index_np_cutoff+1, 
            array('d', [hcent_p.GetXaxis().GetBinLowEdge(i) for i in range(1, index_np_cutoff+3)]))
          hcent_p_p = ROOT.TH1F(
            name_cent+"_pert", name_cent+"_pert", n-index_np_cutoff, 
            array('d', [hcent_p.GetXaxis().GetBinLowEdge(i) for i in range(index_np_cutoff+1, n+2)]))
          for i in range(1, index_np_cutoff+2):
            hcent_p_np.SetBinContent(i, y[i-1])
          for i in range(1, n+1-index_np_cutoff):
            hcent_p_p.SetBinContent(i, y[i-1+index_np_cutoff])

          hcent_p_np.SetLineStyle(5)
          hcent_p_np.SetLineColor(color)
          hcent_p_np.SetLineWidth(3)
          hcent_p_np.Draw('L hist same')

          hcent_p_p.SetLineColor(color)
          hcent_p_p.SetLineWidth(3)
          hcent_p_p.Draw('L hist same')

        else:   # don't show NP / P region on the plot
          h_parton_theory = ROOT.TGraphAsymmErrors(n, x, y, xerrdn, xerrup, yerrdn, yerrup)
          h_parton_theory.SetFillColorAlpha(color, 0.25)
          h_parton_theory.Draw('3 same')
          hcent_p.SetLineColor(color)
          hcent_p.SetLineWidth(3)
          hcent_p.Draw('L hist same')

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
      for ri in range(len(self.theory_response_files)):
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
          h_ch_theory.Draw('3 same')
          hasym_list.append(h_ch_theory)

        hcent.SetLineColor(color)
        hcent.SetLineWidth(3)
        hcent.Draw('L hist same')
        hcent_list.append(hcent)  # Save for legend

        # Dotted lines for error bars
        #hmin.SetLineColor(color)
        #hmin.SetLineWidth(1)
        #hmin.SetLineStyle(2)
        #hmin.Draw('L hist same')

        #hmax.SetLineColor(color)
        #hmax.SetLineWidth(1)
        #hmax.SetLineStyle(2)
        #hmax.Draw('L hist same')


    if plot_pythia:
    
      hPythia, fraction_tagged_pythia = self.pythia_prediction(
        jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin)

      if hPythia:
        hPythia.SetFillStyle(0)
        hPythia.SetMarkerSize(1.5)
        hPythia.SetMarkerStyle(21)
        hPythia.SetMarkerColor(600-6)
        hPythia.SetLineColor(600-6)
        hPythia.SetLineWidth(1)
        hPythia.Draw('E2 same')
      else:
        print('No PYTHIA prediction for %s %s' % (self.observable, obs_label))
        plot_pythia = False
    
    h_sys.DrawCopy('E2 same')
    h.DrawCopy('PE X0 same')
  
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(0.58, 0.87, text)
    
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.58, 0.8, text)

    text = str(min_pt_truth) + ' < #it{p}_{T,jet}^{ch} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.DrawLatex(0.58, 0.73, text)

    text = '#it{R} = ' + str(jetR) + ',   | #it{#eta}_{jet}| < %s' % str(0.9 - jetR)
    text_latex.DrawLatex(0.58, 0.66, text)
    
    subobs_label = self.utils.formatted_subobs_label(self.observable)
    delta = 0.
    if subobs_label:
      text = '%s = %s' % (subobs_label, obs_setting)
      text_latex.DrawLatex(0.58, 0.59, text)
      delta = 0.07
    
    if grooming_setting:
      text = self.utils.formatted_grooming_label(grooming_setting)
      text_latex.DrawLatex(0.58, 0.59-delta, text)
      
      text_latex.SetTextSize(0.04)
      text = '#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged
      text_latex.DrawLatex(0.58, 0.52-delta, text)
    
      if plot_pythia:
        text_latex.SetTextSize(0.04)
        text = ('#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged) + (
          ', #it{f}_{tagged}^{pythia} = %3.3f' % fraction_tagged_pythia)
        text_latex.DrawLatex(0.57, 0.52-delta, text)

    myLegend = ROOT.TLegend(0.22, 0.65, 0.45, 0.9)
    self.utils.setup_legend(myLegend, 0.035)
    myLegend.AddEntry(h, 'ALICE pp', 'pe')
    myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    if plot_pythia:
      myLegend.AddEntry(hPythia, 'PYTHIA8 Monash2013', 'pe')
    if plot_theory:
      if show_parton_theory:
        if show_np_region:
          myLegend.AddEntry(hcent_p_np, 'NLO+NLL (non-perturbative)', 'lf')
          myLegend.AddEntry(hcent_p_p, 'NLO+NLL (perturbative)', 'lf')
        else:
          myLegend.AddEntry(hcent_p, 'NLO+NLL (Parton)', 'lf')
      for ri, lab in enumerate(self.theory_response_labels):
        myLegend.AddEntry(hcent_list[ri], 'NLO+NLL #otimes '+lab, 'lf')
    myLegend.Draw()

    name = 'hUnfolded_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), obs_label,
                                             int(min_pt_truth), int(max_pt_truth), self.file_format)
    if plot_pythia:
      name = 'hUnfolded_R{}_{}_{}-{}_Pythia{}'.format(
        self.utils.remove_periods(jetR), obs_label, int(min_pt_truth), int(max_pt_truth), self.file_format)

    if plot_theory:
      name = 'hUnfolded_R{}_{}_{}-{}_Theory{}'.format(
        self.utils.remove_periods(jetR), obs_label, int(min_pt_truth), int(max_pt_truth), self.file_format)

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
    if plot_theory:
      hcent.Write()
      hmin.Write()
      hmax.Write()
    fFinalResults.Close()

  
  #----------------------------------------------------------------------
  def plot_theory_ratios(self, jetR, obs_label, obs_setting, grooming_setting,
                         min_pt_truth, max_pt_truth, maxbin):
    
    # Get histograms
    name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label,
                                              min_pt_truth, max_pt_truth)
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
      val = abs(2.5*h_ratio_mpi.GetBinContent(int(0.5*h_ratio_mpi.GetNbinsX())))
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
        text_latex.DrawLatex(0.6, 0.59, text)
        delta = 0.07

      myLegend = ROOT.TLegend(0.27, 0.7, 0.55, 0.9)
      self.utils.setup_legend(myLegend, 0.035)
      myLegend.AddEntry(h_ratio_total, 'total ratio', 'pe')
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

    outf.Close


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
    if grooming_setting:
      fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
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
      text = self.utils.formatted_grooming_label(grooming_setting)
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
  def pythia_prediction(self, jetR, obs_setting, obs_label, min_pt_truth,
                        max_pt_truth, maxbin, overlay=False):
  
    hPythia = self.get_pythia_from_response(jetR, obs_label, min_pt_truth,
                                            max_pt_truth, maxbin, overlay)
    n_jets_inclusive = hPythia.Integral(0, hPythia.GetNbinsX()+1)
    n_jets_tagged = hPythia.Integral(hPythia.FindBin(
      self.truth_bin_array(obs_label)[0]), hPythia.GetNbinsX())

    fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
    hPythia.Scale(1./n_jets_inclusive, 'width')
      
    return [hPythia, fraction_tagged_pythia]

  #----------------------------------------------------------------------
  def get_pythia_from_response(self, jetR, obs_label, min_pt_truth, max_pt_truth,
                               maxbin, overlay=False):

    output_dir = getattr(self, 'output_dir_main')

    prev_prelim = False
    if overlay and (jetR == 0.2 or jetR == 0.4) and min_pt_truth == 40 and obs_label == '1':
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
  def plot_final_result_overlay(self, i_config, jetR, overlay_list):
    print('Plotting overlay of {}'.format(overlay_list))

    # Plot overlay of different subconfigs, for fixed pt bin
    for i in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[i]
      max_pt_truth = self.pt_bins_reported[i+1]
      maxbins = [self.obs_max_bins(obs_label)[i] for obs_label in self.obs_labels]

      # Plot PYTHIA
      self.plot_observable_overlay_subconfigs(
        i_config, jetR, overlay_list, min_pt_truth,
        max_pt_truth, maxbins, plot_pythia=True, plot_ratio = True)


  #----------------------------------------------------------------------
  def plot_observable_overlay_subconfigs(self, i_config, jetR, overlay_list, min_pt_truth,
                                         max_pt_truth, maxbins, plot_pythia=False,
                                         plot_nll=False, plot_ratio=False):

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
    pad1.Draw()
    pad1.cd()

    myLegend = ROOT.TLegend(0.56, 0.62, 0.84, 0.89)
    self.utils.setup_legend(myLegend, 0.045)
    
    name = 'hmain_{}_R{}_{{}}_{}-{}'.format(self.observable, jetR, min_pt_truth, max_pt_truth)
    ymax = self.get_maximum(name, overlay_list)

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
        fraction_tagged = getattr(self, name+'_fraction_tagged')

      if (jetR == 0.2 or jetR == 0.4) and min_pt_truth == 40 and obs_label == '1':
        # Use previous preliminary result
        h, h_sys = self.get_h_prelim(jetR)
        # Move error bars to different histogram
        h_sys.SetNameTitle("hSysPrelim%s", "hSysPrelim%s")
        for i in range(1, h.GetNbinsX()+1, 1):
          h.SetBinError(i, 0.)
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
          if (jetR == 0.2 or jetR == 0.4) and min_pt_truth == 40 and obs_label == '1':
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
        myBlankHisto.SetMinimum(0.)
        if plot_ratio:
          myBlankHisto.SetMinimum(2e-4) # Don't draw 0 on top panel
          myBlankHisto.GetYaxis().SetTitleSize(0.065)
          myBlankHisto.GetYaxis().SetTitleOffset(1.1)
          myBlankHisto.GetYaxis().SetLabelSize(0.06)
        myBlankHisto.Draw('E')
        
        # Plot ratio
        if plot_ratio:
          
          c.cd()
          pad2 = ROOT.TPad("pad2", "pad2", 0, 0.02, 1, 0.3)
          pad2.SetTopMargin(0)
          pad2.SetBottomMargin(0.4)
          pad2.SetLeftMargin(0.2)
          pad2.SetRightMargin(0.04)
          pad2.SetTicks(1,1)
          pad2.Draw()
          pad2.cd()
          
          myBlankHisto2 = myBlankHisto.Clone("myBlankHisto_C")
          myBlankHisto2.SetYTitle("#frac{Data}{PYTHIA}")
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
      
      if plot_pythia:
        hPythia, fraction_tagged_pythia = self.pythia_prediction(
          jetR, obs_setting, obs_label, min_pt_truth, max_pt_truth, maxbin, overlay=True)

        plot_errors = False
        if plot_errors:
          hPythia.SetMarkerSize(0)
          hPythia.SetMarkerStyle(0)
          hPythia.SetMarkerColor(color)
          hPythia.SetFillColor(color)
        else:
          hPythia.SetLineColor(color)
          hPythia.SetLineColorAlpha(color, 0.5)
          hPythia.SetLineWidth(4)

      if plot_ratio:
        hRatioSys = h_sys.Clone()
        hRatioSys.SetName('{}_Ratio'.format(h_sys.GetName()))
        if plot_pythia:
          hRatioSys.Divide(hPythia)
        hRatioSys.SetLineColor(0)
        hRatioSys.SetFillColor(color)
        hRatioSys.SetFillColorAlpha(color, 0.3)
        hRatioSys.SetFillStyle(1001)
        hRatioSys.SetLineWidth(0)
        hRatioSys.SetMaximum(1.99)
          
        hRatioStat = h.Clone()
        hRatioStat.SetName('{}_Ratio'.format(h.GetName()))
        if plot_pythia:
          hRatioStat.Divide(hPythia)
        hRatioStat.SetMarkerSize(1.5)
        hRatioStat.SetMarkerStyle(marker)
        hRatioStat.SetMarkerColor(color)
        hRatioStat.SetLineStyle(1)
        hRatioStat.SetLineWidth(2)
        hRatioStat.SetLineColor(color)
        hRatioStat.SetMaximum(1.99)

      pad1.cd()
      if plot_pythia:
        plot_errors = False
        if plot_errors:
          hPythia.DrawCopy('E3 same')
        else:
          hPythia.DrawCopy('L hist same')

      h_sys.DrawCopy('E2 same')
      h.DrawCopy('PE X0 same')
      
      if plot_ratio:
        pad2.cd()
        if plot_pythia:
          hRatioSys.DrawCopy('E2 same')
          hRatioStat.DrawCopy('PE X0 same')

      subobs_label = self.utils.formatted_subobs_label(self.observable)
      text = ''
      if subobs_label:
        text += '%s = %s' % (subobs_label, obs_setting)
      if grooming_setting:
        text += self.utils.formatted_grooming_label(grooming_setting)
      text_list.append(text)
      h_list.append(h)
        
    pad1.cd()
    for h, text in zip(h_list, text_list):
      myLegend.AddEntry(h, text, 'pe')
    myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    if plot_pythia:
      myLegend.AddEntry(hPythia, 'PYTHIA8 Monash2013', 'l')

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(0.25, 0.87, text)
    
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.81, text)

    text = 'Charged jets   anti-#it{k}_{T}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.75, text)
    
    text = '#it{R} = ' + str(jetR) + ',   | #it{#eta}_{jet}| < %s' % str(0.9 - jetR)
    text_latex.DrawLatex(0.25, 0.69, text)
    
    text = str(min_pt_truth) + ' < #it{p}_{T,jet}^{ch} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.63, text)

    myLegend.Draw()

    name = 'h_{}_R{}_{}-{}_{}{}'.format(self.observable, 
                                        self.utils.remove_periods(jetR), int(min_pt_truth), 
                                        int(max_pt_truth), i_config, self.file_format)
    if plot_pythia:
      name = 'h_{}_R{}_{}-{}_Pythia_{}{}'.format(self.observable, self.utils.remove_periods(jetR),
                                                 int(min_pt_truth), int(max_pt_truth),
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
  # Return maximum y-value of unfolded results in a subconfig list
  def get_maximum(self, name, overlay_list):
  
    max = 0.
    for i, subconfig_name in enumerate(self.obs_subconfig_list):
    
      if subconfig_name not in overlay_list:
        continue

      obs_setting = self.obs_settings[i]
      grooming_setting = self.grooming_settings[i]
      obs_label = self.utils.obs_label(obs_setting, grooming_setting)
      
      h = getattr(self, name.format(obs_label))
      if h.GetMaximum() > max:
        max = h.GetMaximum()
        
    return max

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
