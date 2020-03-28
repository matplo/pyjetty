#!/usr/bin/env python3

'''
This script loops through a file list and checks whether an event tree and track tree exist.
It writes out a list of files that pass: files_passed.txt
And a list of files that fail: files_failed.txt

To run:
  python test_ROOT_trees.py -f files.txt

You may want to then remove the failed files by running 'cat files_failed.txt | xargs rm' or similar.
And update the file list by replaces files.txt with files_passed.txt
'''

import os
import sys
import argparse

import ROOT
import uproot

#-----------------------------------------------------------------
def test_ROOT_trees(file_list, n, index, output_dir):

    tree_dir = 'PWGHF_TreeCreator'
    event_tree_name = 'tree_event_char'
    track_tree_name='tree_Particle'
    
    if len(tree_dir) and tree_dir[-1] != '/':
      tree_dir += '/'

    # Loop through files in the file list, and write problematic file to files_to_remove.txt
    outfile_passed = os.path.join(output_dir, 'files_passed_{}.txt'.format(index))
    outfile_failed = os.path.join(output_dir, 'files_failed_{}.txt'.format(index))
    with open(file_list,  'r') as f, open(outfile_passed, 'w') as f_passed, open(outfile_failed, 'w') as f_failed:
        lines = f.read().splitlines()

        print('printing out all files that are missing event tree or track tree...')
        print('')

        nfiles = len(lines)
        nfiles_per_process = int(nfiles / n + 1)
        i_min = (index-1) * nfiles_per_process
        i_max = i_min + nfiles_per_process
        print('nfiles: {}'.format(nfiles))
        print('nfiles_per_process: {}'.format(nfiles_per_process))
        print('i_min: {}'.format(i_min))
        print('i_max: {}'.format(i_max))
        for i, file_name in enumerate(lines):

            if i < i_min or i >= i_max:
                continue
                
            if i % 10 == 0:
                print('{} / {}'.format(i, nfiles))

            # Check if file exists
            if not os.path.exists(file_name):
                print('file does not exist: {}'.format(file_name))
                f_failed.write('{}\n'.format(file_name))
                continue

            # Check if file is a zombie (i.e. exists but empty)
            file = ROOT.TFile(file_name)
            if file.IsZombie():
                print('file is empty: {}'.format(file_name))
                f_failed.write('{}\n'.format(file_name))
                continue
            file.Close()
                
            # Check for event tree
            file_dir = uproot.open(file_name)
            event_tree_path = tree_dir + event_tree_name
            event_tree = file_dir[event_tree_path]
            if not event_tree:
                print(file_name)
                f_failed.write('{}\n'.format(file_name))
                continue
            
            # Check for track tree
            track_tree_path = tree_dir + track_tree_name
            track_tree = file_dir[track_tree_path]
            if not track_tree:
                print(file_name)
                f_failed.write('{}\n'.format(file_name))
            else:
                f_passed.write('{}\n'.format(file_name))

##################################################################
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Test ROOT trees in file list')
  parser.add_argument('-f', '--file_list', action='store',
                      type=str, metavar='file_list',
                      default='files.txt',
                      help='Path of file list')
  parser.add_argument('-n', '--n', action='store',
                      type=int, metavar='n',
                      default='1',
                      help='total number of parallel processes')
  parser.add_argument('-i', '--index', action='store',
                      type=int, metavar='index',
                      default='1',
                      help='index of current process')
  parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./TestOutput',
                      help='Output directory for output to be written to')

 
  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('file list to check: \'{0}\''.format(args.file_list))
  print('Write output file lists to: {}'.format(args.outputDir))
  print('process {} of {}'.format(args.index, args.n))

  # If invalid file list is given, exit
  if not os.path.exists(args.file_list):
    print('File \"{0}\" does not exist! Exiting!'.format(args.file_list))
    sys.exit(0)
  
  test_ROOT_trees(file_list = args.file_list, n = args.n, index = args.index, output_dir = args.outputDir)
