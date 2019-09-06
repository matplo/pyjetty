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
import ROOT

# Fastjet via python (from external library fjpydev)
import fastjet as fj
import fjext

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

#---------------------------------------------------------------
def analysis_utils():

  print('These are some analysis utilities for jet analysis with track dataframe')

#---------------------------------------------------------------
# Convert ROOT TTree to pandas dataframe
# Return merged track+event dataframe from a given input file
# Returned dataframe has one row per jet constituent:
#     run_number, ev_id, ParticlePt, ParticleEta, ParticlePhi
#---------------------------------------------------------------
def load_dataframe(inputFile, tree_name):
  
  # Load event tree into dataframe, and apply event selection
  event_tree_name = 'PWGHF_TreeCreator/tree_event_char'
  event_tree = uproot.open(inputFile)[event_tree_name]
  if not event_tree:
    print('Tree {} not found in file {}'.format('tree_event_char', inputFile))
  event_df_orig = event_tree.pandas.df(['run_number', 'ev_id', 'z_vtx_reco','is_ev_rej'])
  event_df_orig.reset_index(drop=True)
  event_df = event_df_orig.query('is_ev_rej == 0')
  event_df.reset_index(drop=True)

  # Load track tree into dataframe
  track_tree_name = 'PWGHF_TreeCreator/{}'.format(tree_name)
  track_tree = uproot.open(inputFile)[track_tree_name]
  if not track_tree:
    print('Tree {} not found in file {}'.format(tree_name, inputFile))
  track_df_orig = track_tree.pandas.df()

  # Merge event info into track tree
  track_df = pandas.merge(track_df_orig, event_df, on=['run_number', 'ev_id'])
  return track_df

#---------------------------------------------------------------
# Return fastjet:PseudoJets from a given track dataframe
#---------------------------------------------------------------
def get_fjparticles(df_tracks):
  
  # Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of pt,eta,phi
  fj_particles = fjext.vectorize_pt_eta_phi(df_tracks['ParticlePt'].values, df_tracks['ParticleEta'].values, df_tracks['ParticlePhi'].values)
  
  return fj_particles

#---------------------------------------------------------------
# Plot and save a 1D histogram
#---------------------------------------------------------------
def plotHist(h, outputFilename, drawOptions = "", setLogy = False, setLogz = False):
  
  h.SetLineColor(1)
  h.SetLineWidth(1)
  h.SetLineStyle(1)
  
  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd()
  ROOT.gPad.SetLeftMargin(0.15)
  if setLogy:
    c.SetLogy()
  if setLogz:
    c.SetLogz()
  ROOT.gPad.SetLeftMargin(0.15)

  h.Draw(drawOptions)
  c.SaveAs(outputFilename)
  c.Close()

#----------------------------------------------------------------------
if __name__ == '__main__':
  
  analysis_utils()
