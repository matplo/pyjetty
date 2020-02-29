#!/usr/bin/env python3

"""
    Class to store various jet info, used primarily for jet matching.
    
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# Base class
from pyjetty.alice_analysis.process.base import common_base

################################################################
class JetInfo(common_base.CommonBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(JetInfo, self).__init__(**kwargs)

    # Store the matching candidates
    self.matching_candidates = []
    self.closest_jet = None
    self.closest_jet_deltaR = 1000.
    self.match = None
