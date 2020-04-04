#!/usr/bin/env python3

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt

#----------------------------------------------------------------------
def plot_files_sizes(files):

    # Get list of files
    list_of_files = []
    with open(files) as f:
        list_of_files = [fn.strip() for fn in f.readlines()]

    # Loop through file list, and histogram the size
    file_sizes = [os.path.getsize(file)/1e6 for file in list_of_files]
    
    # Plot distribution
    plt.hist(file_sizes, bins=50)
    plt.yscale('log')
    plt.ylabel('N files')
    plt.xlabel('file size (MB)')
    plt.savefig('h_file_sizes.pdf')
    plt.close()

#----------------------------------------------------------------------
if __name__ == '__main__':

  # Define arguments
  parser = argparse.ArgumentParser(description='Plot distribution of file sizes of a list of file paths')
  parser.add_argument('-f', '--files', action='store',
                      type=str, metavar='files',
                      default='files.txt',
                      help='Path of file paths')

  # Parse the arguments
  args = parser.parse_args()
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.files):
    print('File \"{0}\" does not exist! Exiting!'.format(args.files))
    sys.exit(0)

  plot_files_sizes(files = args.files)
