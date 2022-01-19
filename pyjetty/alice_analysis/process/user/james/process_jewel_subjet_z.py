#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of track information
  and do jet-finding, and save basic histograms.
  This specific code is to run over jewel generator data to produce histograms (at truth level) that will then be compared to
  the data. That is, this base is to run over MC, but only at truth level, without response matrices.
  Author: James Mulligan (james.mulligan@berkeley.edu)
  Based on code by: Reynier Cruz Torres (reynier@lbl.gov) 
"""

from __future__ import print_function

# General
import os
import sys
import argparse

# Data analysis and plotting
import ROOT
import numpy as np

# Fastjet via python (from external library heppy)
import fastjet as fj

# Load pyjetty ROOT utils
ROOT.gSystem.Load('libpyjetty_rutil')

# Base class
from pyjetty.alice_analysis.process.user.substructure import process_jewel_generated_base

################################################################
class Process_CurvesFromJewelTracks_subjet_z(process_jewel_generated_base.CurvesFromJewelTracks):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', **kwargs):
  
    # Initialize base class
    super(Process_CurvesFromJewelTracks_subjet_z, self).__init__(input_file, config_file, output_dir, **kwargs)
    
    self.observable = self.observable_list[0] #define as first element in each 'config' in the input file

    # Define subjet finders (from first observable defined)
    self.subjet_def = {}
    for subjetR in self.obs_settings[self.observable_list[0]]:
      self.subjet_def[subjetR] = fj.JetDefinition(fj.antikt_algorithm, subjetR)

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects(self,label=''):
    
    for jetR in self.jetR_list:
      for r in self.obs_settings[self.observable]:

        # Histogram name based on recoil subtraction method
        if not self.thermal_subtraction_method:
          name = f'h_{self.observable}_JetPt_R{jetR}_{r}{label}'
          h = ROOT.TH2F(name, name+';p_{T,ch jet};#it{z}_{#it{r}}', 200, 0, 200, 100, 0, 1.)
          setattr(self, name, h) 
        elif 'gridsub' in self.thermal_subtraction_method.lower():
          for gridsize in self.gridsizes:
            name = f'h_{self.observable}_JetPt_R{jetR}_{r}_gridsub_{gridsize}{label}'
            h = ROOT.TH2F(name, name+';p_{T,ch jet};#it{z}_{#it{r}}', 200, 0, 200, 100, 0, 1.)
            setattr(self, name, h) 
        elif '4momsub' in self.thermal_subtraction_method.lower():
          name = f'h_{self.observable}_JetPt_R{jetR}_{r}_4momsub{label}'
          h = ROOT.TH2F(name, name+';p_{T,ch jet};#it{z}_{#it{r}}', 200, 0, 200, 100, 0, 1.)
          setattr(self, name, h) 

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  #---------------------------------------------------------------
  def fill_jet_histograms(self, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
                          obs_label, jet_pt_ungroomed, suffix=None, label=''):

    # Get hist names
    if not self.thermal_subtraction_method:
      name = f'h_{self.observable}_JetPt_R{jetR}_{obs_label}{label}'
    elif 'gridsub' in self.thermal_subtraction_method.lower():
      name = f'h_{self.observable}_JetPt_R{jetR}_{obs_label}_gridsub_{suffix}{label}'
    elif '4momsub' in self.thermal_subtraction_method.lower():
      name = f'h_{self.observable}_JetPt_R{jetR}_{obs_label}_4momsub{label}'

    # Convert df of pt/eta/phi of thermal particles to list of fastjet particles, for convenience
    df_thermal = []
    if len(self.event_bck_df) != 0:
      df_thermal = self.get_fjparticles(self.event_bck_df)

    # Correct the jet pt by subtracting thermal particles within R of jet axis
    holes_in_jet = []
    negative_pt = 0.
    for thermal_particle in df_thermal:
      if jet.delta_R(thermal_particle) < jetR:
        if np.random.uniform() >= self.thermal_rejection_fraction:
          negative_pt += thermal_particle.pt()
          holes_in_jet.append(thermal_particle)
    jet_pt = jet.pt() - negative_pt   # corrected pt: shower+recoil-holes

    # For a given jet, find all inclusive subjets of a given subjet radius
    cs_subjet = fj.ClusterSequence(jet.constituents(), self.subjet_def[obs_setting])
    subjets = fj.sorted_by_pt(cs_subjet.inclusive_jets())

    # Get leading subjet (accounts for holes)
    _, leading_subjet_pt = self.leading_jet(subjets, holes_in_jet, obs_setting)
    z = leading_subjet_pt / jet_pt
    
    # If z=1, it will be default be placed in overflow bin -- prevent this
    if np.isclose(z, 1.):
      z = 0.999

    getattr(self, name).Fill(jet_pt, z)

  #---------------------------------------------------------------
  # Return leading jet (or subjet)
  #---------------------------------------------------------------
  def leading_jet(self, jets, fj_hadrons_negative, jetR):

    leading_jet = None
    leading_jet_pt = 0.
    for jet in jets:
                      
      # Get the corrected jet pt by subtracting the negative recoils within R
      jet_pt = jet.pt()
      for temp_hadron in fj_hadrons_negative:
        if jet.delta_R(temp_hadron) < jetR:
          jet_pt -= temp_hadron.pt()
              
      if not leading_jet:
        leading_jet = jet
        leading_jet_pt = jet_pt
      
      if jet_pt > leading_jet_pt:
        leading_jet = jet
        leading_jet_pt = jet_pt

    return leading_jet, leading_jet_pt

##################################################################
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Process Generator For Theory Comparison')
  parser.add_argument('-f', '--inputFile', action='store',
                      type=str, metavar='inputFile',
                      default='AnalysisResults.root',
                      help='Path of ROOT file containing TTrees')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='config/analysis_config.yaml',
                      help="Path of config file for analysis")
  parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./TestOutput',
                      help='Output directory for output to be written to')

  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('inputFile: \'{0}\''.format(args.inputFile))
  print('configFile: \'{0}\''.format(args.configFile))
  print('ouputDir: \'{0}\"'.format(args.outputDir))
  print('----------------------------------------------------------------')
  
  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
    sys.exit(0)
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = Process_CurvesFromJewelTracks_subjet_z(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_gen()