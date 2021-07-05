# This script aggregates the hdf5 nsubjettiness processing output
# into a single file.

import os
import sys
import yaml
import argparse
import h5py
import numpy as np

def aggregate(config_file, filelist, output_dir):

    # List of arrays to aggregate
    observables = ['X_four_vectors', 'X_Nsub', 'y', 'jet_pt', 'delta_pt', 'matched_pt', 'matched_deltaR', 'jet_angularity', 'jet_mass', 'jet_theta_g', 'jet_subjet_z', 'hadron_z', 'multiplicity_0000', 'multiplicity_0150', 'multiplicity_0500', 'multiplicity_1000']

    # Read config file
    with open(config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    jetR_list = config['jetR']
    jet_pt_bins = config['jet_pt_bins']
    max_distance_list = config['constituent_subtractor']['max_distance']
    event_types = ['hard', 'combined', 'combined_matched']

    # We have separate training data for:
    # - the hard-event and the combined-event
    # - different jet R
    # - different constituent subtraction R_max
    
    # Create a flat dict of empty objects for each combination
    # We use keys equal to the keys in the input file
    output = {}
    output['delta_pt_random_cone'] = np.array([])
    for event_type in event_types:
        for jetR in jetR_list:
            for jet_pt_bin in jet_pt_bins:
                for R_max in max_distance_list:
                    for observable in observables:
                        if accepted_setting(observables, event_type, R_max):
                            output_key = f'{observable}_{event_type}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                            output[output_key] = np.array([])
    N_list = None
    beta_list = None

    # Loop through file list
    with open(filelist) as f:
        files = [line.rstrip() for line in f]
        n_files = len(files)
    for i,filename in enumerate(files):
        
        #if i > 300:
        #    break

        if i%100 == 0:
            print(f'{i}/{n_files}')

        with h5py.File(filename,'r') as hdf:
        
            if not N_list:
                N_list = list(hdf['N_list'][:])
                beta_list = list(hdf['beta_list'][:])

            delta_pt_random_cone_i = hdf['delta_pt_random_cone'][:]
            if output['delta_pt_random_cone'].any():
                output['delta_pt_random_cone'] = np.concatenate((output['delta_pt_random_cone'], delta_pt_random_cone_i),axis=0)
            else:
                output['delta_pt_random_cone'] = delta_pt_random_cone_i
    
            for event_type in event_types:
                for jetR in jetR_list:
                    for jet_pt_bin in jet_pt_bins:
                        for R_max in max_distance_list:
                            for observable in observables:
                                if accepted_setting(observables, event_type, R_max):
                            
                                    # Get arrays from file
                                    output_key = f'{observable}_{event_type}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                                    output_i = hdf[output_key][:]
                                    if output_i.shape[0] > 0:

                                        # Concatenate to master array
                                        output_aggregated = output[output_key]
                                        if output_aggregated.any():
                                            output[output_key] = np.concatenate((output_aggregated, output_i),axis=0)
                                        else:
                                            output[output_key] = output_i

                                        if i%100 == 0:
                                            if observable in ['X_Nsub']:
                                                print(f'{output_key} -- {output_aggregated.shape}')

    #  Shuffle data set
    for event_type in event_types:
        for jetR in jetR_list:
            for jet_pt_bin in jet_pt_bins:
                for R_max in max_distance_list:
                    if accepted_setting(observables, event_type, R_max):
                    
                        output_key_y = f'y_{event_type}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                        idx = np.random.permutation(len(output[output_key_y]))
                    
                        for observable in observables:

                            output_key = f'{observable}_{event_type}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                            output_aggregated = output[output_key]
                            if output_aggregated.shape[0] == idx.shape[0]:
                                output[output_key] = output_aggregated[idx]
                                print(f'shuffled {output_key}: {output[output_key].shape}')

                                # Remove nan
                                where_are_nan = np.isnan(output[output_key])
                                output[output_key][where_are_nan] = 0.
                            elif output_aggregated.shape[0] > 0:
                                print(f'MISMATCH of shape for {output_key}: {output_aggregated.shape} vs. {idx.shape}')

    # Write jet arrays to file
    if 'X_four_vectors' in observables:
        output_filename = 'nsubjettiness_with_four_vectors.h5'
    else:
        output_filename = 'nsubjettiness_without_four_vectors.h5' 
        
    with h5py.File(os.path.join(output_dir, output_filename), 'w') as hf:
    
        hf.create_dataset('N_list', data=np.array(N_list))
        hf.create_dataset('beta_list', data=np.array(beta_list))
        hf.create_dataset('delta_pt_random_cone', data=np.array(output['delta_pt_random_cone']))
        
        for event_type in event_types:
            for jetR in jetR_list:
                for jet_pt_bin in jet_pt_bins:
                    for R_max in max_distance_list:
                        for observable in observables:
                            if accepted_setting(observables, event_type, R_max):

                                output_key = f'{observable}_{event_type}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                                output_aggregated = output[output_key]
                                hf.create_dataset(output_key, data=output_aggregated)
    
    print('done.')

def accepted_setting(observables, event_type, R_max):
    
    if 'hard' in event_type and np.isclose(R_max, 0.):
        return True

    if 'combined' in event_type:
        if 'X_four_vectors' in observables:
            if not np.isclose(R_max, 0.):
                return True
        else:
            return True

    return False

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

    # Parse the arguments
    args = parser.parse_args()

    print('Configuring...')
    print('configFile: \'{0}\''.format(args.configFile))
    print('ouputDir: \'{0}\"'.format(args.outputDir))

    # If invalid configFile is given, exit
    if not os.path.exists(args.configFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    # If invalid outputdir is given, exit
    fileList = os.path.join(args.outputDir, 'files.txt')
    if not os.path.exists(fileList):
        print('File \"{0}\" does not exist! Exiting!'.format(fileList))
        sys.exit(0)

    aggregate(config_file=args.configFile, filelist=fileList, output_dir=args.outputDir)
