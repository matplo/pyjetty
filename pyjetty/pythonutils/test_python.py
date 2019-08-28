#!/usr/bin/env python

import sys
print("Python version...")
print (sys.version)
print("Version info...")
print (sys.version_info)
print ("Python executable...")
print (sys.executable)
print ("Trying to import ROOT ...")
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
print ("... done.")
