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

# Fastjet via python (from external library heppy)
import fjcontrib

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
    
    return self.is_det_jet_accepted(self, jet_truth)

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
  # Return leading jet (or subjet)
  #---------------------------------------------------------------
  def leading_jet(self, jets):
  
    leading_jet = None
    for jet in jets:
    
      if not leading_jet:
        leading_jet = jet
        
      if jet.pt() > leading_jet.pt():
        leading_jet = jet

    return leading_jet

  #---------------------------------------------------------------
  # Perform grooming and return Lund declustering object
  # (cpptools/src/fjcontrib/custom/GroomerShop.hh)
  #
  # Note that the GroomerShop returns a pointer to the
  # LundDeclustering object -- this object is a class member,
  # of the GroomerShop, so it  will only stay in scope as long as
  # the GroomerShop  remains in scope.
  #---------------------------------------------------------------
  def groom(self, gshop, grooming_setting, jetR):
  
    if 'sd' in grooming_setting:
    
      zcut = grooming_setting['sd'][0]
      beta = grooming_setting['sd'][1]
      return gshop.soft_drop(beta, zcut, jetR)

    elif 'dg' in grooming_setting:
    
      if len(gshop.jet().constituents()) < 2:
        return None
        
      a = grooming_setting['dg'][0]
      
      if a == 'max_pt_softer':
        return gshop.max_pt_softer()
      elif a == 'max_z':
        return gshop.max_z()
      elif a == 'max_kt':
        return gshop.max_kt()
      elif a == 'max_kappa':
        return gshop.max_kappa()
      elif a == 'max_tf':
        return gshop.max_tf()
      elif a == 'min_tf':
        return gshop.min_tf()
      else:
        return gshop.dynamical(a)
    
    else:
      sys.exit('grooming_setting {} not recognized.'.format(grooming_setting))
    


