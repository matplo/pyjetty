#!/usr/bin/env python3

"""
  Plotting utilities for jet substructure analysis with track dataframe.

  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import math
import yaml

# Data analysis and plotting
import numpy as np
from array import *
import ROOT

# Base class
from pyjetty.alice_analysis.analysis.user.james import plotting_utils_base

################################################################
class PlottingUtils(plotting_utils_base.PlottingUtilsBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, output_dir = '.', config_file = '', R_max = None, thermal = False, groomer_studies = False, **kwargs):
    super(PlottingUtils, self).__init__(output_dir, config_file, R_max, thermal, groomer_studies, **kwargs)

    print(self)

