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
import subprocess
import multiprocessing as mp

#---------------------------------------------------------------------------
def download_data(period):

    # Set train info
    year = '2018'
    train_name = 'HF_TreeCreator'
    train_PWG = 'pass3/PWGHF'
    filename = 'AnalysisResults.root'
    train_number = '550_20201223-0456'
    
    if period == 'LHC18q': # pass3
        train_number += '_child_1'
        runlist = ['000296414', '000296510', '000296379', '000296377', '000296309', '000296433', '000296068', '000296133', '000296423', '000296065', '000296550', '000295588', '000295586', '000296270']
    elif period == 'LHC18r': # pass3
        train_number += '_child_2'
        runlist = ['000296894', '000297446', '000297544', '000296899', '000297479', '000297442', '000297415', '000296934', '000297590', '000297380', '000297123', '000296694', '000296903', '000297218']
    else:
        sys.exit('Error: period {} unrecognized'.format(period))

    # Create output dir and cd into it
    output_dir = '/mnt/rstorage/alice/data/LHC18qr/550'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    os.chdir(output_dir)
    print('output dir: {}'.format(output_dir))

    # Loop through runs, and start a download for each run in parallel
    for run in runlist:
        p = mp.Process(target=download_run, args=(year, period, run, train_PWG, train_name, train_number))
        p.start()

#---------------------------------------------------------------------------
def download_run(year, period, run, train_PWG, train_name, train_number):

    train_output_dir = '/alice/data/{}/{}/{}/{}/{}/{}'.format(year, period, run, train_PWG, train_name, train_number)
    print('train_output_dir: {}'.format(train_output_dir))

    # Construct list of subdirectories (i.e. list of files to download)
    temp_filelist_name = 'subdirs_temp_{}.txt'.format(run)
    cmd = 'alien_ls {} > {}'.format(train_output_dir, temp_filelist_name)
    os.system(cmd)
    with open(temp_filelist_name) as f:
        subdirs_all = f.read().splitlines()
        subdirs = [ x for x in subdirs_all if x.isdigit() ]
    os.remove(temp_filelist_name)

    # Remove any empty directories
    if os.path.exists(run):
        cmd = 'find {} -empty -type d -delete'.format(run)
        os.system(cmd)

    # Copy the files
    for subdir in subdirs:
        
        # Skip any directory that already exists
        subdir_path = '{}/{}'.format(run, subdir)
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
  parser.add_argument('-p', '--period', action='store',
                      type=str, metavar='period',
                      default='LHC18q',
                      help='Run period')

  # Parse the arguments
  args = parser.parse_args()

  download_data(period = args.period)
