#!/usr/bin/env python3

"""
    Class to store various jet info, used primarily for jet matching.
    
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# Base class
from pyjetty.alice_analysis.process.base import common_base

################################################################
class jet_info(common_base.common_base):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(jet_info, self).__init__(**kwargs)

    self.closest_jet_pp = None
    self.closest_jet_AA = None
    self.matching_candidates_pp = []
    self.matching_candidates_AA = []

