#!/usr/bin/env python3

'''
Creates files.txt file for use in processing script. 
REQUIRES that eff_smear has been run and fs data exists.

Written by Ezra Lesser (elesser@berkeley.edu), Spring 2020
'''

from __future__ import division, print_function

import os
import argparse


#####################################################################
## Parse input arguments
#####################################################################

parser = argparse.ArgumentParser(description="Creates files.txt, list of all *.root")
parser.add_argument("-i", "--inputDir", action="store", type=str, metavar="inputDir",
                    default="./TestInput", help="Master directory where ROOT files are stored")
parser.add_argument("-o", "--outputDir", action="store", type=str, metavar="outputDir",
                    default="./TestOutput", help="Output directory for files.txt")
args = parser.parse_args()

print('Configuring...')
print('inputDir: "{0}"'.format(args.inputDir))
print('ouputDir: "{0}"'.format(args.outputDir))

# If invalid inputDir is given, exit
if not os.path.exists(args.inputDir):
    print('File "{0}" does not exist! Exiting!'.format(args.inputDir))
    sys.exit(0)

data_dir = args.inputDir
if data_dir[-1] != '/':
    data_dir += '/'
out_dir = args.outputDir
if out_dir[-1] != '/':
    out_dir += '/'

#####################################################################
## Find all file paths and save to files.txt
#####################################################################

filepath = out_dir + "files.txt"
with open(filepath, 'w') as txtf:
    for dirpath, dirnames, filenames in os.walk(data_dir):
        for filename in [f for f in filenames if f.endswith(".root")]:
            txtf.write(os.path.join(dirpath, filename) + '\n')
print("ROOT file list has been saved to %s." % filepath)
