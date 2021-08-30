#!/usr/bin/env python

from __future__ import print_function

import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libpyjetty_TGlauberMC.dylib")

# tgmc = ROOT.TGlauberMC() # //constructor
# call something from the TGlauberMC
ROOT.TGlauberMC_f.runAndSaveNtuple(1000)
