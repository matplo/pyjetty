#!/usr/bin/env python3

'''
Script to download Pb-Pb train output from AliEn.

python download_data.py -p LHC18q

On hiccup:
  - ssh to hiccupds
  - start a screen session
  - enter alidock, then `alienv enter AliRoot/latest`, then get token
  - python download_data.py -p LHC18q
  
Note that if token expires or otherwise crashes, the script will automatically detect where to start copying again

'''

import argparse
import os
import sys
import yaml
import subprocess
import multiprocessing as mp

#---------------------------------------------------------------------------
def download_data(config_file):

    # Initialize config
    with open(config_file, 'r') as stream:
      config = yaml.safe_load(stream)

    period = config['period']
    parent_dir = config['parent_dir']
    year = config['year']
    train_name = config['train_name']
    train_PWG = config['train_PWG']
    train_number = config['train_number']
    runlist = config['runlist']
    output_dir = config['output_dir']
    
    if 'pt_hat_bins' in config:
        pt_hat_bins = config['pt_hat_bins']
    else:
        pt_hat_bins = None

    # Create output dir and cd into it
    output_dir = os.path.join(output_dir, period)
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    os.chdir(output_dir)
    print('output dir: {}'.format(outputdir))
    
    # Loop through runs, and start a download for each run in parallel
    for run in runlist:
        p = mp.Process(target=download_run, args=(year, period, run, train_PWG, train_name, train_number, pt_hat_bins))
        p.start()

#---------------------------------------------------------------------------
def download_run(parent_dir, year, period, run, train_PWG, train_name, train_number, pt_hat_bins):

    if parent_dir == 'data':
    
        train_output_dir = '/alice/{}/{}/{}/{}/{}/{}/{}'.format(parent_dir, year, period, run, train_PWG, train_name, train_number)
        
        download(train_output_dir, run)
        
    elif parent_dir == 'sim':
    
        for pt_hat_bin in pt_hat_bins:
        
            train_output_dir = '/alice/{}/{}/{}/{}/{}/{}/{}/{}'.format(parent_dir, year, period, pt_hat_bin, run, train_PWG, train_name, train_number)
            
            download(train_output_dir, run, pt_hat_bin)

#---------------------------------------------------------------------------
def download(train_output_dir, run_path, pt_hat_bin=None):

    print('train_output_dir: {}'.format(train_output_dir))
    
    if pt_hat_bin:
        run_path = '{}/{}'.format(pt_hat_bin, run)
    else:
        run_path = run

    # Construct list of subdirectories (i.e. list of files to download)
    temp_filelist_name = 'subdirs_temp_{}.txt'.format(run)
    cmd = 'alien_ls {} > {}'.format(train_output_dir, temp_filelist_name)
    os.system(cmd)
    with open(temp_filelist_name) as f:
        subdirs_all = f.read().splitlines()
        subdirs = [ x for x in subdirs_all if x.isdigit() ]
    os.remove(temp_filelist_name)

    # Remove any empty directories
    if os.path.exists(run_path):
        cmd = 'find {} -empty -type d -delete'.format(run_path)
        os.system(cmd)

    # Copy the files
    for subdir in subdirs:
        
        # Skip any directory that already exists
        subdir_path = '{}/{}'.format(run_path, subdir)
        if not os.path.exists(subdir_path):
            os.makedirs(subdir_path)
            print('downloading: {}'.format(subdir_path))
        else:
            continue

        logfile_name = "log_{}.txt".format(run)
        with open('log_{}.txt'.format(run), "a") as logfile:
            cmd = 'alien_cp alien://{}/{}/AnalysisResults.root {}'.format(train_output_dir, subdir, subdir_path)
            print(cmd, file=logfile)
            subprocess.run(cmd, check=False, shell=True, stdout=logfile, stderr=logfile)

#----------------------------------------------------------------------
if __name__ == '__main__':

    # Define arguments
    parser = argparse.ArgumentParser(description='Download train output')
    parser.add_argument('-c', '--configFile', action='store',
                        type=str, metavar='configFile',
                        default='config.yaml',
                        help='Path of config file')

    # Parse the arguments
    args = parser.parse_args()

    print('Configuring...')
    print('configFile: \'{0}\''.format(args.configFile))

    # If invalid configFile is given, exit
    if not os.path.exists(args.configFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    download_data(config_file = args.configFile)
