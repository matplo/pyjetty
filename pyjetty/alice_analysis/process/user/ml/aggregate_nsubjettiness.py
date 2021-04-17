# This script aggregates the hdf5 nsubjettiness processing output
# into a single file.

import os
import argparse
import h5py
import numpy as np
from numba import jit

def aggregate(filelist, output_dir, include_four_vector):

    # Loop through file list
    if include_four_vector:
        X_total = None
    X_Nsub_total = None
    y = None
    N_list = None
    beta_list = None
    with open(filelist) as f:
        files = [line.rstrip() for line in f]
        n_files = len(files)

    for i,filename in enumerate(files):

        with h5py.File(filename,'r') as hdf:

            if include_four_vector:
                X = hdf['X'][:]
            X_Nsub = hdf['X_Nsub'][:]
            y = hdf['y'][:]

            if not N_list:
                if include_four_vector:
                    X_total = X
                X_Nsub_total = X_Nsub
                y_total = y
                N_list = list(hdf['N_list'][:])
                beta_list = list(hdf['beta_list'][:])
            else:
                if include_four_vector:                
                    X_total = np.concatenate((X_total, X),axis=0)
                X_Nsub_total = np.concatenate((X_Nsub_total, X_Nsub),axis=0)
                y_total = np.concatenate((y_total, y),axis=0)

        if i%100 == 0:
            if include_four_vector:
                print(f'{i}/{n_files} -- {X_total.shape} -- {X_Nsub_total.shape} -- {y_total.shape[0]}')
            else:
                print(f'{i}/{n_files} -- {X_Nsub_total.shape} -- {y_total.shape[0]}')

    #  Shuffle data set # check again!
    idx = np.random.permutation(len(y_total))
    if include_four_vector:
        X_shuffled = X_total[idx] 
    X_Nsub_shuffled, y_shuffled = X_Nsub_total[idx], y_total[idx]
    print(f'shuffled: {X_Nsub_shuffled.shape} -- {y_shuffled.shape}')
    # Write jet arrays to file
    with h5py.File(os.path.join(output_dir, 'nsubjettiness.h5'), 'w') as hf:
        hf.create_dataset('y', data=y_shuffled)
        if include_four_vector:
            hf.create_dataset('X', data=X_shuffled)
        hf.create_dataset('X_Nsub', data=X_Nsub_shuffled)
        hf.create_dataset('N_list', data=np.array(N_list))
        hf.create_dataset('beta_list', data=np.array(beta_list))
    
    print('done.')

##################################################################
if __name__ == '__main__':

    # Define arguments
    parser = argparse.ArgumentParser(description='Aggregate pp AA')
    parser.add_argument('-o', '--outputDir', action='store',
                        type=str, metavar='outputDir',
                        default='./TestOutput',
                        help='Output directory for output to be written to')
    parser.add_argument('--include_four_vector',action='store_true')

    # Parse the arguments
    args = parser.parse_args()

    print('Configuring...')
    print('ouputDir: \'{0}\"'.format(args.outputDir))
    print(f'aggregate four-vectors: {args.include_four_vector}')

    # If invalid configFile is given, exit
    fileList = os.path.join(args.outputDir, 'files.txt')
    if not os.path.exists(fileList):
        print('File \"{0}\" does not exist! Exiting!'.format(fileList))
        sys.exit(0)

    aggregate(filelist=fileList, output_dir=args.outputDir, include_four_vector=args.include_four_vector)
