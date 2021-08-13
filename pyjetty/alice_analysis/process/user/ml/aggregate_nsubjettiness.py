# This script aggregates the hdf5 nsubjettiness processing output
# into a single file.

import os
import sys
import yaml
import argparse
import h5py
import numpy as np

def aggregate(config_file, filelist, output_dir):

    n_files_max = 10000

    # List of arrays to aggregate
    observables = ['X_four_vectors', 'X_Nsub', 'y', 'jet_pt', 'delta_pt', 'matched_pt', 'matched_deltaR', 'jet_angularity', 'jet_mass', 'jet_theta_g', 'jet_subjet_z', 'hadron_z', 'multiplicity_0000', 'multiplicity_0150', 'multiplicity_0500', 'multiplicity_1000']

    # Read config file
    with open(config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    jetR_list = config['jetR']
    jet_pt_bins = config['jet_pt_bins']
    max_distance_list = config['constituent_subtractor']['max_distance']
    event_types = ['hard', 'combined_matched']

    # Create a list of keys to loop over
    # We have separate training data for:
    # - the hard-event and the combined-event
    # - different jet R
    # - different constituent subtraction R_max
    output_keys = []
    for event_type in event_types:
        for jetR in jetR_list:
            for jet_pt_bin in jet_pt_bins:
                for R_max in max_distance_list:
                    for observable in observables:
                        if accepted_setting(observable, event_type, R_max):
                            output_key = f'{observable}_{event_type}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                            output_keys.append(output_key)

    # We will create a virtual dataset for each combination
    # See: https://docs.h5py.org/en/stable/vds.html
    
    # First, we need to find the total shape for each observable set
    shapes, total_shapes = determine_shapes(output_keys, filelist, n_files_max)
    print('Determined shapes.')

    # Now, create the virtual dataset
    # We use keys equal to the keys in the input file
    if 'X_four_vectors' in observables:
        output_filename_unshuffled = 'nsubjettiness_with_four_vectors_unshuffled.h5'
    else:
        output_filename_unshuffled = 'nsubjettiness_without_four_vectors_unshuffled.h5'
    with h5py.File(os.path.join(output_dir, output_filename_unshuffled), 'w') as hf:
        
        for output_key in output_keys:
            print(f'Creating virtual dataset for {output_key}')
            print(f'  Total shape: {total_shapes[output_key]}')
            layout = h5py.VirtualLayout(shape=total_shapes[output_key], dtype=np.float64)

            # Loop through file list
            layout_index = 0
            with open(filelist) as f:
                files = [line.rstrip() for line in f]
                n_files = len(files)
            for i,filename in enumerate(files):
                if not accept_file(i, n_files, n_files_max, log=False):
                    break

                # Create virtual source
                #print(f'  Creating virtual source for file {i} with shape {shapes[output_key][i]}')
                source = h5py.VirtualSource(filename, output_key, shapes[output_key][i], dtype=np.float64)

                # Insert the source into the layout
                new_layout_index = layout_index + shapes[output_key][i][0]
                if len(total_shapes[output_key]) == 1:
                    layout[layout_index:new_layout_index] = source
                elif len(total_shapes[output_key]) == 2:
                    layout[layout_index:new_layout_index,:] = source
                else:
                    layout[layout_index:new_layout_index,:,:] = source
                
                layout_index = new_layout_index
                #print(f'new layout_index: {layout_index}')
                
            # Add virtual dataset to output file
            hf.create_virtual_dataset(output_key, layout)

    print('Virtual dataset created.')
    print()
    
    # Write N_list, beta_list
    if 'X_four_vectors' in observables:
        output_filename_shuffled = 'nsubjettiness_with_four_vectors_shuffled.h5'
    else:
        output_filename_shuffled = 'nsubjettiness_without_four_vectors_shuffled.h5'
    with h5py.File(filename, 'r') as hf:
        N_list = list(hf['N_list'][:])
        beta_list = list(hf['beta_list'][:])
    with h5py.File(os.path.join(output_dir, output_filename_shuffled), 'w') as hf_shuffled:
        hf_shuffled.create_dataset('N_list', data=np.array(N_list))
        hf_shuffled.create_dataset('beta_list', data=np.array(beta_list))
    
    #  Shuffle data set -- create randomized indices for each set of observables
    idx = {}
    N_list = None
    beta_list = None
    for output_key in output_keys:
        print(f'Shuffling dataset for {output_key}')
        
        suffix = output_key.split('_', 1)[1]
        if suffix not in idx:
            idx[suffix] = np.array([])
        
        with h5py.File(os.path.join(output_dir, output_filename_unshuffled), 'r') as hf_unshuffled:
            output = hf_unshuffled[output_key][:]
            
            if not np.any(idx[suffix]):
                idx[suffix] = np.random.permutation(len(output))

            if output.shape[0] == idx[suffix].shape[0]:
                output_shuffled = output[idx[suffix]]
                print(f'shuffled {output_key}: {output_shuffled.shape}')

                # Remove nan
                where_are_nan = np.isnan(output)
                output[where_are_nan] = 0.
            elif output.shape[0] > 0:
                print(f'MISMATCH of shape for {output_key}: {output.shape} vs. {idx[suffix].shape}')
                
            with h5py.File(os.path.join(output_dir, output_filename_shuffled), 'w') as hf_shuffled:
                hf_shuffled.create_dataset(output_key, data=output)

    print('done.')

# Determine shapes of lists in all files
def determine_shapes(output_keys, filelist, n_files_max):

    shapes = {}
    shapes['delta_pt_random_cone'] = []
    for output_key in output_keys:
        shapes[output_key] = []

    with open(filelist) as f:
        files = [line.rstrip() for line in f]
        n_files = len(files)
    for i,filename in enumerate(files):
        with h5py.File(filename,'r') as hdf:
            if not accept_file(i, n_files, n_files_max):
                break

            delta_pt_random_cone_i = hdf['delta_pt_random_cone'][:]
            shapes['delta_pt_random_cone'].append(hdf['delta_pt_random_cone'][:].shape)
            for output_key in output_keys:
                shapes[output_key].append(hdf[output_key][:].shape)

    total_shapes = {}
    for key,val in shapes.items():
        total_shape = np.sum([shape[0] for shape in val])
        if len(val[0]) == 1:
            total_shapes[key] = (total_shape,)
        elif len(val[0]) == 2:
            total_shapes[key] = (total_shape, 146)
        else:
            total_shapes[key] = (total_shape, 800, 4)

    return shapes,total_shapes

def accept_file(i, n_files, n_files_max, log=True):
    if log:
        if i%10 == 0:
            print(f'{i}/{n_files}')
    if i > n_files_max:
        return False
    return True

def accepted_setting(observable, event_type, R_max):
    
    if 'hard' in event_type and np.isclose(R_max, 0.):
        if 'delta_pt' not in observable and 'matched' not in observable:
            return True

    if 'combined' in event_type:
        if 'X_four_vectors' in observable:
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
