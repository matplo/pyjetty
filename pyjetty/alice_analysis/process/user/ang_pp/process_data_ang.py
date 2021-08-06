#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of track information
  and do jet-finding, and save basic histograms.

  Author: Ezra Lesser (elesser@berkeley.edu) with much code borrowed 
          from original script by James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import argparse

# Data analysis and plotting
import ROOT

# Fastjet via python (from external library heppy)
import fjext

# Base class
from pyjetty.alice_analysis.process.user.substructure import process_data_base

################################################################
class ProcessData_ang(process_data_base.ProcessDataBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):

    # Initialize base class
    super(ProcessData_ang, self).__init__(
      input_file, config_file, output_dir, debug_level, **kwargs)

    # Dictionary containing the optimal R_max values for each jet R of interest
    self.R_max = { "0.2" : "0.1", "0.4" : "0.25" }

    self.pt_bins = array('d', list(range(5, 305, 5)))
    self.obs_bins = np.concatenate((np.linspace(0, 0.009, 10), np.linspace(0.01, 0.1, 19),
                                    np.linspace(0.11, 0.8, 70)))

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects(self):

    for jetR in self.jetR_list:

      for observable in self.observable_list:
        # Should only be one: observable == "ang"
        if observable != "ang":
          raise ValueError("Observable %s is not implemented in this script" % observable)

        for grooming_setting in self.obs_grooming_settings[observable]:

          if grooming_setting:
            grooming_label = self.utils.grooming_label(grooming_setting)
            if self.is_pp:
              name = 'h_%s_JetPt_R%s_%s' % (observable, jetR, grooming_label)
              h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                            len(self.obs_bins)-1, self.obs_bins)
              h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
              h.GetYaxis().SetTitle('#it{#lambda}_{g,ch}^{#it{#alpha}}')
              setattr(self, name, h)

            else:
              # Pb-Pb: have several R_max for contituent subtraction
              max_distance = self.max_distance if isinstance(self.max_distance, list) \
                             else self.max_distance[jetR]
              for R_max in max_distance:
                name = 'h_%s_JetPt_R%s_%s_Rmax%s' % (
                  observable, jetR, grooming_label, self.R_max[jetR])
                h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                        len(self.obs_bins)-1, self.obs_bins)
                h.GetXaxis().SetTitle('#it{p}_{T,ch jet}')
                h.GetYaxis().SetTitle('#it{#lambda}_{g,ch}^{#it{#alpha}}')
                setattr(self, name, h)

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  #---------------------------------------------------------------
  def fill_jet_histograms(self, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
                          obs_label, jet_pt_ungroomed, suffix):

    # Calculate angularity and fill MPI scaling histogram
    ang = fjext.lambda_beta_kappa(jet, jet_groomed_lund.pair(), obs_setting, 1, jetR) \
          if grooming_setting else fjext.lambda_beta_kappa(jet, obs_setting, 1, jetR)

    # Fill histograms
    getattr(self, "h_ang_JetPt_R%s_%s%s" % (jetR, obs_label, suffix)).Fill(
      jet_pt_ungroomed, ang)


##################################################################
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Process data')
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
    raise ValueError("File \"%s\" does not exist" % args.inputFile)

  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    raise ValueError("File \"%s\" does not exist" % args.configFile)

  analysis = ProcessData_ang(
    input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
  analysis.process_data()
