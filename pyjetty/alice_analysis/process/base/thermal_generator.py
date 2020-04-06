#!/usr/bin/env python3

"""
  Class to generate a thermal event.
  The pt is sampled from a Gamma distribution,
  and the eta and phi are sampled from a uniform distribution.
  
  Previously I have used the following in order to approximate the
  delta_pt of data at 5.02 TeV for full jets:
    R=0.2: N_avg  = 3500, sigma_N = 500, beta = 0.4
    R=0.4: N_avg  = 4000, sigma_N = 500, beta = 0.5
  See: https://alice-notes.web.cern.ch/system/files/notes/analysis/814

  Author: James Mulligan
"""

from __future__ import print_function

# Data analysis and plotting
import uproot
import pandas
import numpy as np
import random

import fjext

# Base class
from pyjetty.alice_analysis.process.base import common_base

################################################################
class ThermalGenerator(common_base.CommonBase):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, N_avg = 2500, sigma_N = 500, beta = 0.4, alpha = 2, eta_max = 0.9, **kwargs):
    super(ThermalGenerator, self).__init__(**kwargs)
    
    # Gamma distribution parameters: f(pt;alpba,beta) = pt^(alpha-1) exp(-pt/beta)
    # Commonly one takes b=<pT>/2
    # Note: Choosing a=1 gives exponential.
    self.alpha = alpha
    self.beta = beta
    self.N_avg = N_avg
    
    # Sample N particles from a Gaussian characterized by N_avg, sigma_N
    self.sigma_N = sigma_N
    
    # Set eta limit for uniform sampling
    self.eta_max = eta_max
    
  #---------------------------------------------------------------
  # Return a thermal event, as a SeriesGroupBy of fastjet::PseudoJet
  #---------------------------------------------------------------
  def load_event(self):
  
    # Decide how many tracks to generate
    N_tracks = int(np.random.normal(self.N_avg, self.sigma_N))
  
    # Generate tracks, and populate a dataframe
    pt_array = np.random.gamma(self.alpha, self.beta, N_tracks)
    eta_array = np.random.uniform(-self.eta_max, self.eta_max, N_tracks)
    phi_array = np.random.uniform(0., 2*np.pi, N_tracks)
    
    # Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of pt,eta,phi
    user_index_offset = int(-1e6)
    fj_particles = fjext.vectorize_pt_eta_phi(pt_array, eta_array, phi_array, user_index_offset)
    return fj_particles
