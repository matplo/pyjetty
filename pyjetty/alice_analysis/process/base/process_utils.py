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
from pyjetty.alice_analysis.process.base import common_utils

################################################################
class ProcessUtils(common_utils.CommonUtils):
  
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
  
    delta_phi = np.abs(jet.phi() - phi)
    delta_eta = jet.eta() - eta

    if delta_phi > np.pi:
      delta_phi = 2*np.pi - delta_phi
    
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
  # Perform dynamical grooming (or other custom grooming method
  # in src/fjcontrib/custom/DynamicalGroomer.hh)
  #---------------------------------------------------------------
  def dy_groom(self, dy_groomer, jet, a):
  
    if a == 'max_pt_softer':
      return dy_groomer.max_pt_softer(jet)
    elif a == 'max_z':
      return dy_groomer.max_z(jet)
    elif a == 'max_kt':
      return dy_groomer.max_kt(jet)
    elif a == 'max_kappa':
      return dy_groomer.max_kappa(jet)
    elif a == 'max_tf':
      return dy_groomer.max_tf(jet)
    elif a == 'min_tf':
      return dy_groomer.min_tf(jet)
    else:
      return dy_groomer.result(jet, a)
