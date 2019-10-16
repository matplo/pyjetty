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
class process_utils(common_base.common_base):
  
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(process_utils, self).__init__(**kwargs)
  
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
    
    accept_jet = True
    
    for track in jet_truth.constituents():
      
      if track.pt() > 100.:
        accept_jet = False

    return accept_jet

  #---------------------------------------------------------------
  # Remove periods from a label
  #---------------------------------------------------------------
  def remove_periods(self, text):
    
    string = str(text)
    return string.replace('.', '')
