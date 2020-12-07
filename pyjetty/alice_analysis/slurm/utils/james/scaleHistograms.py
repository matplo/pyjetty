# Macro to scale histograms of all Pt-hard bins, using xsec from a yaml file.
# This script expects files X/AnalysisResults.root, and will output scaled histograms
# to the same file, in a new output list with suffix "Scaled". The script will automatically loop over
# all output lists, subject to some simple criteria that covers basic use cases (can be adapted as needed).
#
# There is an option "bRemoveOutliers" to remove outliers from certain histograms. The features are
# currently hard-coded below so you will need to modify the code as needed. This feature is adapted from code of Raymond Ehlers.
#
# Author: James Mulligan (james.mulligan@berkeley.edu)
#

import ROOT
import argparse
import ctypes
import os
import sys
import yaml

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

###################################################################################
# Main function
def scaleHistograms(configFile, remove_unscaled):

  # Option to remove outliers from specified histograms
  # If the average bin content stays below the "outlierLimit" for "outlierNBinsThreshold" bins, it is removed
  bRemoveOutliers = False
  outlierLimit = 2
  outlierNBinsThreshold=4
  
  # Option to print out detailed info about scaling and outlier removal
  verbose = False
  
  # Read the cross-section, and scale histograms
  with open(configFile, 'r') as stream:
    config = yaml.safe_load(stream)

  # Find EndPtHardBin
  EndPtHardBin = 0
  for i in range(1,100):
    if i in config:
      EndPtHardBin = i
    else:
      break
  print("EndPtHardBin: {}".format(EndPtHardBin))

  # Compute average number of events per bin
  nEventsSum = 0
  for bin in range(1, EndPtHardBin+1):
    f = ROOT.TFile("{}/AnalysisResults.root".format(bin), "UPDATE")
    hNevents = f.Get("hNevents")
    nEvents = hNevents.GetBinContent(2)
    nEventsSum += nEvents
  nEventsAvg = nEventsSum*1./EndPtHardBin
  print("nEventsAvg per bin: {}".format(nEventsAvg))

  for bin in range(1, EndPtHardBin+1):
    print("ooo Scaling Pt-hard bin {} of {}".format(bin, EndPtHardBin))

    f = ROOT.TFile("{}/AnalysisResults.root".format(bin), "UPDATE")
    hCrossSection = config[bin]
    nEvents = f.Get("hNevents").GetBinContent(2)
    eventScaleFactor = nEvents / nEventsAvg
    scaleFactor = hCrossSection / eventScaleFactor
    print("ooo eventScaleFactor: {0}".format(eventScaleFactor))
    print("ooo total scaleFactor: {0}".format(scaleFactor))

    # Now, scale all the histograms
    keys = f.GetListOfKeys()
    for key in keys:
      name = key.GetName()
      if "Scaled" in name:
        continue
      if "roounfold" in name:
        continue
      obj = f.Get(name)
      if obj:
        ScaleAllHistograms(obj, scaleFactor, f, verbose, bRemoveOutliers, outlierLimit, outlierNBinsThreshold, bin-1, EndPtHardBin, name)
      else:
        print('obj not found!')
  
      obj.Write("%sScaled" % obj.GetName())
      
      if remove_unscaled:
        f.Delete('{};1'.format(obj.GetName()))  # Remove unscaled histogram

    f.Close()

###################################################################################
# Function to iterate recursively through an object to scale all TH1/TH2/THnSparse
def ScaleAllHistograms(obj, scaleFactor, f, verbose, bRemoveOutliers=False, limit=2, nBinsThreshold=4, pTHardBin=0, EndPtHardBin=20, taskName=""):
  
  # Set Sumw2 if not already done
  if obj.InheritsFrom(ROOT.THnBase.Class()):
    if obj.GetSumw2() is 0:
      obj.Sumw2()
      if verbose:
        print('Set Sumw2 on {}'.format(obj.GetName()))
  else:
    if obj.GetSumw2N() is 0:
      obj.Sumw2()
      if verbose:
        print('Set Sumw2 on {}'.format(obj.GetName()))
  
  if obj.InheritsFrom(ROOT.TProfile.Class()):
    if verbose:
      print("TProfile %s not scaled..." % obj.GetName())
  elif obj.InheritsFrom(ROOT.TH2.Class()):
    obj.Scale(scaleFactor)
    if verbose:
      print("TH2 %s was scaled..." % obj.GetName())
  elif obj.InheritsFrom(ROOT.TH1.Class()):
    if bRemoveOutliers:
      name = obj.GetName()
      #only perform outlier removal on these couple histograms
      if "Pt" in name:
        removeOutliers(pTHardBin, EndPtHardBin, obj, verbose, limit, nBinsThreshold, 1, taskName)
    obj.Scale(scaleFactor)
    if verbose:
      print("TH1 %s was scaled..." % obj.GetName())
  elif obj.InheritsFrom(ROOT.THnBase.Class()):
    obj.Scale(scaleFactor)
    if verbose:
      print("THnSparse %s was scaled..." % obj.GetName())
  else:
    if verbose:
      print("Not a histogram!")
      print(obj.GetName())
    for subobj in obj:
      ScaleAllHistograms(subobj, scaleFactor, f, verbose, bRemoveOutliers, limit, nBinsThreshold, pTHardBin, taskName)

###################################################################################
# Function to remove outliers from a TH3 (i.e. truncate the spectrum), based on projecting to the y-axis
# It truncates the 3D histogram based on when the 1D projection 4-bin moving average has been above
# "limit" for "nBinsThreshold" bins.
def removeOutliers(pTHardBin, EndPtHardBin, hist, verbose, limit=2, nBinsThreshold=4, dimension=3, taskName=""):

  #Project to the pT Truth axis
  if dimension==3:
    histToCheck = hist.ProjectionY("{}_projBefore".format(hist.GetName()))
  if dimension==2:
    histToCheck = hist.ProjectionX("{}_projBefore".format(hist.GetName()))
  if dimension==1:
    histToCheck = hist

  # Check with moving average
  foundAboveLimit = False
  cutLimitReached = False
  # The cut index is where we decided cut on that row
  cutIndex = -1
  nBinsBelowLimitAfterLimit = 0
  # nBinsThreshold= n bins that are below threshold before all bins are cut

  if verbose:
    (preMean, preMedian) = GetHistMeanAndMedian(histToCheck)
    
  for index in range(0, histToCheck.GetNcells()):
    if verbose:
      print("---------")
    avg = MovingAverage(histToCheck, index = index, numberOfCountsBelowIndex = 2, numberOfCountsAboveIndex = 2)
    if verbose:
      print("Index: {0}, Avg: {1}, BinContent: {5}, foundAboveLimit: {2}, cutIndex: {3}, cutLimitReached: {4}".format(index, avg, foundAboveLimit, cutIndex, cutLimitReached, histToCheck.GetBinContent(index)))
    if avg > limit:
      foundAboveLimit = True
        
    if not cutLimitReached:
      if foundAboveLimit and avg <= limit:
        if cutIndex == -1:
          cutIndex = index
        nBinsBelowLimitAfterLimit += 1
          
      if nBinsBelowLimitAfterLimit != 0 and avg > limit:
        # Reset
        cutIndex = -1
        nBinsBelowLimitAfterLimit = 0

      if nBinsBelowLimitAfterLimit > nBinsThreshold:
        cutLimitReached = True
        break #no need to continue the loop - we found our cut index

  # Do not perform removal here because then we miss values between the avg going below
  # the limit and crossing the nBinsThreshold
  if verbose:
    print("Hist checked: {0}, cut index: {1}".format(histToCheck.GetName(), cutIndex))
  
  # Use on both TH1 and TH2 since we don't start removing immediately, but instead only after the limit
  if cutLimitReached:
    if verbose:
      print("--> --> --> Removing outliers")
    # Check for values above which they should be removed by translating the global index
    x = ctypes.c_int(0)
    y = ctypes.c_int(0)
    z = ctypes.c_int(0)
    for index in range(0, hist.GetNcells()):
      # Get the bin x, y, z from the global bin
      hist.GetBinXYZ(index, x, y, z)
      if dimension==3:
        if y.value >= cutIndex:
          if hist.GetBinContent(index) > 1e-3:
            if verbose:
              print("Cutting for index {}. y bin {}. Cut index: {}".format(index, y, cutIndex))
            hist.SetBinContent(index, 0)
            hist.SetBinError(index, 0)
      if dimension==2:
        #for the response matrix the pT Truth is on the y-Axis
        if hist.GetName()=="hResponseMatrixEMCal":
          x.value=y.value
        if x.value >= cutIndex:
          if hist.GetBinContent(index) > 1e-3:
            if verbose:
              print("Cutting for index {}. x bin {}. Cut index: {}".format(index, x, cutIndex))
            hist.SetBinContent(index, 0)
            hist.SetBinError(index, 0)
      if dimension==1:
        if x.value >= cutIndex:
          if hist.GetBinContent(index) > 1e-3:
            if verbose:
              print("Cutting for index {}. x bin {}. Cut index: {}".format(index, x, cutIndex))
            hist.SetBinContent(index, 0)
            hist.SetBinError(index, 0)
          
  else:
    if verbose:
      print("Hist {} did not have any outliers to cut".format(hist.GetName()))

  # Check the mean and median
  # Use another temporary hist
  if dimension==3:
    histToCheckAfter = hist.ProjectionY()
  if dimension==2:
    if hist.GetName()=="hResponseMatrixEMCal":
      histToCheckAfter = hist.ProjectionY()
  else:
      histToCheckAfter = hist.ProjectionX()
  if dimension==1:
    histToCheckAfter = hist

  if verbose:
    (postMean, postMedian) = GetHistMeanAndMedian(histToCheckAfter)
    print("Pre  outliers removal mean: {}, median: {}".format(preMean, preMedian))
    print("Post outliers removal mean: {}, median: {}".format(postMean, postMedian))
  outlierFilename = "{}OutlierRemoval_{}.pdf".format(hist.GetName())
  if "Pt" in hist.GetName():
    plotOutlierPDF(histToCheck, histToCheckAfter, pTHardBin, EndPtHardBin, outlierFilename, verbose, "hist E", True)

########################################################################################################
def GetHistMeanAndMedian(hist):
  # Median
  # See: https://root-forum.cern.ch/t/median-of-histogram/7626/5
  x = ctypes.c_double(0)
  q = ctypes.c_double(0.5)
  # Apparently needed to be safe(?)
  hist.ComputeIntegral()
  hist.GetQuantiles(1, x, q)
    
  mean = hist.GetMean()
  return (mean, x.value)

########################################################################################################
def MovingAverage(hist, index, numberOfCountsBelowIndex = 0, numberOfCountsAboveIndex = 2):
  """
  # [-2, 2] includes -2, -1, 0, 1, 2
  """
  # Check inputs
  if numberOfCountsBelowIndex < 0 or numberOfCountsAboveIndex < 0:
    print("Moving average number of counts above or below must be >= 0. Please check the values!")
          
  count = 0.
  average = 0.
  for i in range(index - numberOfCountsBelowIndex, index + numberOfCountsAboveIndex + 1):
    # Avoid going over histogram limits
    if i < 0 or i >= hist.GetNcells():
      continue
    #print("Adding {}".format(hist.GetBinContent(i)))
    average += hist.GetBinContent(i)
    count += 1
    
  #if count != (numberOfCountsBelowIndex + numberOfCountsAboveIndex + 1):
  #    print("Count: {}, summed: {}".format(count, (numberOfCountsBelowIndex + numberOfCountsAboveIndex + 1)))
  #exit(0)

  return average / count

########################################################################################################
# Plot basic histogram    ##############################################################################
########################################################################################################
def plotOutlierPDF(h, hAfter, pTHardBin, EndPtHardBin, outputFilename, verbose, drawOptions = "", setLogy = False):
  
  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd()
  if setLogy:
    c.SetLogy()
  h.GetXaxis().SetRangeUser(0,250)
  h.Draw("hist")
  h.SetLineColor(616)
  hAfter.SetLineColor(820)
  hAfter.Draw("same hist")
  
  leg1 = ROOT.TLegend(0.17,0.7,0.83,0.85,"outlier removal of Bin {}".format(pTHardBin))
  leg1.SetFillColor(10)
  leg1.SetBorderSize(0)
  leg1.SetFillStyle(0)
  leg1.SetTextSize(0.04)
  leg1.AddEntry(h, "before", "l")
  leg1.AddEntry(hAfter, "after", "l")
  leg1.Draw("same")

  c.Print("{}".format(outputFilename))
  '''
  if pTHardBin == 0: #if first pt-hard bin, open a .pdf
    if verbose:
      print("Add first pT Hard bin to pdf with name: {0}".format(outputFilename))
  elif pTHardBin==EndPtHardBin-1: #otherwise add pages to the file
    if verbose:
      print("Add last pT Hard bin to pdf with name: {0}".format(outputFilename))
    c.Print("{})".format(outputFilename))
  else: #otherwise close the file
    if verbose:
      print("Add further pT Hard bin to pdf with name: {0}".format(outputFilename))
    c.Print("{}".format(outputFilename))
  '''

  c.Close()

########################################################################################################
# Get Jet radius from list analysis label                       #############################################
########################################################################################################
def getRadiusFromlistName(listName):

  radius = 0.0
  if "01" in listName:
    radius = 0.1
  if "02" in listName:
    radius = 0.2
  elif "03" in listName:
    radius = 0.3
  elif "04" in listName:
    radius = 0.4
  elif "05" in listName:
    radius = 0.5
  elif "06" in listName:
    radius = 0.6
  return radius

#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  print("Executing scaleHistograms.py...")
  print("")
  
  # Define arguments
  parser = argparse.ArgumentParser(description='Plot analysis histograms')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='analysis_config.yaml',
                      help="Path of config file for jetscape analysis")
  parser.add_argument('-r', '--remove_unscaled', action='store_true',
                      help='Remove unscaled histograms')
  
  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('configFile: \'{0}\''.format(args.configFile))
  print('----------------------------------------------------------------')
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  scaleHistograms(configFile = args.configFile, remove_unscaled = args.remove_unscaled)
