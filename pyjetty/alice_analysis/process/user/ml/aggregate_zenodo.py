'''
Script to create a virtual dataset from a list of hdf5 files.
The virtual dataset can then be accessed as you would a normal single hdf5 file.

Author: James Mulligan <james.mulligan@berkeley.edu>
'''

import os
import h5py
import numpy as np

def aggregate(output_dir):

    # Create files with list of paths of the individual hdf5 files
    # You can use e.g.: find <file_dir> -name "*.h5" >"files.txt"
    filelist = './files.txt'

    # Location of aggragated file
    output_dir = '.'

    # List of arrays to aggregate
    observables = ['X_four_vectors', 'X_Nsub', 'y']
    
    jetR_list = [0.4]
    jet_pt_bins = [[100, 125]]
    max_distance_list = [0, 0.25]
    event_types = ['hard', 'combined_matched']
    K_max = 30

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
                        if accepted_setting(event_type, R_max):
                            output_key = f'{observable}_{event_type}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                            output_keys.append(output_key)

    # We will create a virtual dataset for each combination
    # See: https://docs.h5py.org/en/stable/vds.html
    
    # First, we need to find the total shape for each observable set
    shapes, total_shapes, N_list, beta_list = determine_shapes(output_keys, filelist, K_max)
    print('Determined shapes.')

    # Now, create the virtual dataset
    # We use keys equal to the keys in the input file
    output_filename_unshuffled = 'nsubjettiness_with_four_vectors_unshuffled.h5'
    with h5py.File(os.path.join(output_dir, output_filename_unshuffled), 'w') as hf:
        
        for output_key in output_keys:
            print(f'Creating virtual dataset for {output_key}')
            print(f'  Total shape: {total_shapes[output_key]}')
            layout = h5py.VirtualLayout(shape=total_shapes[output_key], dtype=np.float64)

            # Loop through file list
            layout_index = 0
            with open(filelist) as f:
                files = [line.rstrip() for line in f]
            for i,filename in enumerate(files):

                # Create virtual source
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
                
            # Add virtual dataset to output file
            hf.create_virtual_dataset(output_key, layout)

        # Write N_list, beta_list
        hf.create_dataset('N_list', data=np.array(N_list))
        hf.create_dataset('beta_list', data=np.array(beta_list))

    print('Virtual dataset created.')
    print()
    
# Determine shapes of lists in all files
def determine_shapes(output_keys, filelist, K_max):

    shapes = {}
    for output_key in output_keys:
        shapes[output_key] = []

    with open(filelist) as f:
        files = [line.rstrip() for line in f]
        n_files = len(files)
    for i,filename in enumerate(files):
        with h5py.File(filename,'r') as hdf:
            if i%10 == 0:
                print(f'{i}/{n_files}')

            for output_key in output_keys:
                shapes[output_key].append(hdf[output_key][:].shape)

            if i==0:
                N_list = list(hdf['N_list'][:])
                beta_list = list(hdf['beta_list'][:])

    total_shapes = {}
    for key,val in shapes.items():
        total_shape = np.sum([shape[0] for shape in val])
        if len(val[0]) == 1:
            total_shapes[key] = (total_shape,)
        elif len(val[0]) == 2:
            total_shapes[key] = (total_shape, 3*K_max-4)
        else:
            total_shapes[key] = (total_shape, 800, 4)

    return shapes, total_shapes, N_list, beta_list

def accepted_setting(event_type, R_max):
    
    if 'hard' in event_type and np.isclose(R_max, 0.):
        return True

    if 'combined' in event_type:
        if not np.isclose(R_max, 0.):
            return True

    return False

##################################################################
if __name__ == '__main__':
    aggregate()
