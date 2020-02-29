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

# Base class
from pyjetty.alice_analysis.process.base import common_base

################################################################
class ProcessUtils(common_base.CommonBase):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(ProcessUtils, self).__init__(**kwargs)
  
  #---------------------------------------------------------------
  # Check if det-jet passes acceptance criteria
  #---------------------------------------------------------------
  def is_det_jet_accepted(self, jet_det):

    for track in jet_det.constituents():

      if track.pt() > 100.:
        return False

    return True

  #---------------------------------------------------------------
  # Check if truth-jet passes acceptance criteria
  #---------------------------------------------------------------
  def is_truth_jet_accepted(self, jet_truth):
    
    return is_det_jet_accepted(self, jet_truth)

  #---------------------------------------------------------------
  # Compute delta-R (eta-phi) between a PseudoJet and a given eta,phi value
  #---------------------------------------------------------------
  def delta_R(self, jet, eta, phi):
  
    delta_phi = jet.phi() - phi
    delta_eta = jet.eta() - eta
    
    deltaR = np.sqrt(delta_phi*delta_phi + delta_eta*delta_eta)
    return deltaR
    
  #---------------------------------------------------------------
  # Get the leading constituent of a jet
  #---------------------------------------------------------------
  def get_leading_constituent(self, jet):
  
    leading_particle = None
    leading_particle_pt = 0.
    for particle in jet.constituents():
      if particle.pt() > leading_particle_pt:
        leading_particle = particle
      
    return leading_particle

  #---------------------------------------------------------------
  # Remove periods from a label
  #---------------------------------------------------------------
  def remove_periods(self, text):
    
    string = str(text)
    return string.replace('.', '')
