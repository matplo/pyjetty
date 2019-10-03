#! /usr/bin/env python

import argparse
import os
import yaml
import ROOT

import roounfold_rg

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

#----------------------------------------------------------------------
def run_rg_analysis(config_file):

  # Set which actions to perform
  do_unfolding = True
  do_systematics = True
  do_plot_final_result = True

  # List of possible systematic variations
  # Regularization parameter is automatically included
  kMain, kUnfoldingRange1, kUnfoldingRange2, kPrior1, kPrior2, kTrackEff = range(0, 6)

  # Set which systematics should be performed
  systematics_list = [kMain, kTrackEff]

  # Load yaml config
  with open(config_file, 'r') as stream:
    config = yaml.safe_load(stream)

  output_dir = config['output_dir']
  file_format = config['file_format']

  # Create output dir
  if not output_dir.endswith("/"):
    output_dir = output_dir + "/"
  if not os.path.exists(output_dir):
    os.makedirs(output_dir)

  # Load paths to processing output, to be unfolded
  main_data = config['main_data']
  main_response = config['main_response']

  if kTrackEff in systematics_list:
    trkeff_data = config['trkeff_data']
    trkeff_response = config['trkeff_response']

  #----------------------------------------------------------------------
  if do_unfolding:
    print('Perform unfolding for all systematic variations...')

    # Main result
    output_dir_main = os.path.join(output_dir, 'main')
    if not os.path.isdir(output_dir_main):
      os.makedirs(output_dir_main)

    analysis_main = roounfold_rg.roounfold_rg(main_data, main_response, config_file, output_dir_main, file_format)
    analysis_main.roounfold_rg()

    # Tracking efficiency variation
    if kTrackEff in systematics_list:
      output_dir_trkeff = os.path.join(output_dir, 'trkeff')
      if not os.path.isdir(output_dir_trkeff):
        os.makedirs(output_dir_trkeff)

      analysis_trkeff = roounfold_rg.roounfold_rg(trkeff_data, trkeff_response, config_file, output_dir_trkeff, file_format)
      analysis_trkeff.roounfold_rg()

  #----------------------------------------------------------------------
  if do_systematics:
    print('Compute systematics...')

  #----------------------------------------------------------------------
  if do_plot_final_result:
    print('Plot final results...')


#----------------------------------------------------------------------
if __name__ == '__main__':

  # Define arguments
  parser = argparse.ArgumentParser(description='Unfold theta_g distribution')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='analysis_config.yaml',
                      help='Path of config file for analysis')

  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('configFile: \'{0}\''.format(args.configFile))
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  run_rg_analysis(config_file = args.configFile)
