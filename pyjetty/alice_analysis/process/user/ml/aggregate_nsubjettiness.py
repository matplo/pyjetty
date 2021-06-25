# This script aggregates the hdf5 nsubjettiness processing output
# into a single file.

import os
import sys
import yaml
import argparse
import h5py
import numpy as np

def aggregate(config_file, filelist, output_dir, include_four_vector):

    # Read config file
    with open(config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    jetR_list = config['jetR']
    max_distance_list = config['constituent_subtractor']['max_distance']
    event_types = ['hard', 'combined']

    # List of arrays to aggregate
    observables = ['X_Nsub', 'y', 'pt', 'delta_pt']

    # We have separate training data for:
    # - the hard-event and the combined-event
    # - different jet R
    # - different constituent subtraction R_max
    
    # Create a flat dict of empty objects for each combination
    # We use keys equal to the keys in the input file
    output = {}
    for event_type in event_types:
        for jetR in jetR_list:
            for R_max in max_distance_list:
                for observable in observables:
                    output_key = f'{observable}_{event_type}_R{jetR}_Rmax{R_max}'
                    output[output_key] = np.array([])
    N_list = None
    beta_list = None

    # Loop through file list
    with open(filelist) as f:
        files = [line.rstrip() for line in f]
        n_files = len(files)
    for i,filename in enumerate(files):
    
        if i%100 == 0:
            print(f'{i}/{n_files}')

        with h5py.File(filename,'r') as hdf:
        
            if not N_list:
                N_list = list(hdf['N_list'][:])
                beta_list = list(hdf['beta_list'][:])
        
            for event_type in event_types:
                for jetR in jetR_list:
                    for R_max in max_distance_list:
                        for observable in observables:

                            # Skip R_max for hard event
                            if event_type == 'hard':
                                if not np.isclose(R_max, max_distance_list[0]):
                                    continue

                            # Get arrays from file
                            output_key = f'{observable}_{event_type}_R{jetR}_Rmax{R_max}'
                            if observable == 'X':
                                if include_four_vector:
                                    output_i = hdf[output_key][:]
                            else:
                                output_i = hdf[output_key][:]

                            # Concatenate to master array
                            output_aggregated = output[output_key]
                            if output_aggregated.any():
                                if observable == 'X':
                                    if include_four_vector:
                                        output[output_key] = np.concatenate((output_aggregated, output_i),axis=0)
                                else:
                                    output[output_key] = np.concatenate((output_aggregated, output_i),axis=0)
                            else:
                                if observable == 'X':
                                    if include_four_vector:
                                        output[output_key] = output_i
                                else:
                                    output[output_key] = output_i

                            if i%100 == 0 and observable == 'X_Nsub':
                                print(f'{output_key} -- {output_aggregated.shape}')

    #  Shuffle data set
    for event_type in event_types:
        for jetR in jetR_list:
            for R_max in max_distance_list:
                
                output_key_X = f'X_{event_type}_R{jetR}_Rmax{R_max}'
                output_aggregated_X = output[output_key_X]
                output_key_X_Nsub = f'X_Nsub_{event_type}_R{jetR}_Rmax{R_max}'
                output_aggregated_X_Nsub = output[output_key_X_Nsub]
                output_key_y = f'y_{event_type}_R{jetR}_Rmax{R_max}'
                output_aggregated_y = output[output_key_y]

                idx = np.random.permutation(len(output_aggregated_y))
                output[output_key_X] = output_aggregated_X[idx]
                output[output_key_X_Nsub] = output_aggregated_X_Nsub[idx]
                output[output_key_y] = output_aggregated_y[idx]
                print(f'shuffled {output_key_X_Nsub}: {output[output_key_X_Nsub].shape}')
                print(f'shuffled {output_key_y}: {output[output_key_y].shape}')
        
    # Write jet arrays to file
    with h5py.File(os.path.join(output_dir, 'nsubjettiness.h5'), 'w') as hf:
    
        hf.create_dataset('N_list', data=np.array(N_list))
        hf.create_dataset('beta_list', data=np.array(beta_list))
        
        for event_type in event_types:
            for jetR in jetR_list:
                for R_max in max_distance_list:
                    for observable in observables:
                    
                        # Skip R_max for hard event
                        if event_type == 'hard':
                            if not np.isclose(R_max, max_distance_list[0]):
                                continue

                        output_key = f'{observable}_{event_type}_R{jetR}_Rmax{R_max}'
                        output_aggregated = output[output_key]
                        if observable == 'X':
                            if include_four_vector:
                                hf.create_dataset(output_key, data=output_aggregated)
                        else:
                            hf.create_dataset(output_key, data=output_aggregated)
    
    print('done.')

##################################################################
if __name__ == '__main__':

    # Define arguments
    parser = argparse.ArgumentParser(description='Aggregate pp AA')
    parser.add_argument('-c', '--configFile', action='store',
                        type=str, metavar='configFile',
                        default='../../../config/ml/ppAA.yaml',
                        help='Path of config file for analysis')
    parser.add_argument('-o', '--outputDir', action='store',
                        type=str, metavar='outputDir',
                        default='./TestOutput',
                        help='Output directory for output to be written to')
    parser.add_argument('--include_four_vector',action='store_true')

    # Parse the arguments
    args = parser.parse_args()

    print('Configuring...')
    print('configFile: \'{0}\''.format(args.configFile))
    print('ouputDir: \'{0}\"'.format(args.outputDir))
    print(f'aggregate four-vectors: {args.include_four_vector}')

    # If invalid configFile is given, exit
    if not os.path.exists(args.configFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    # If invalid outputdir is given, exit
    fileList = os.path.join(args.outputDir, 'files.txt')
    if not os.path.exists(fileList):
        print('File \"{0}\" does not exist! Exiting!'.format(fileList))
        sys.exit(0)

    aggregate(config_file=args.configFile, filelist=fileList, output_dir=args.outputDir, include_four_vector=args.include_four_vector)
