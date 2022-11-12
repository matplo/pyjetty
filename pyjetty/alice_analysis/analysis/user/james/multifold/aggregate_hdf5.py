# This script aggregates the hdf5 processing output into a single file.

import os
import sys
import yaml
import argparse
import h5py
import numpy as np
import shutil

from pyjetty.alice_analysis.process.base import process_utils
utils = process_utils.ProcessUtils()
 
#---------------------------------------------------------------
def aggregate(config_file, filename_list, output_dir):
    '''
    Main

    We merge all observables present in the processing output here.
    At analysis time, one can then select a subset of observables to analyze.

    Note: We only consider a single jet radius for now.
    '''

    n_files_max = 100
    n_files_total = len(filename_list)

    # Read config file
    with open(config_file, 'r') as stream:
        config = yaml.safe_load(stream)
    #jetR_list = config['jetR']

    # First, get the (nested) keys that we want to merge
    print('Loading sample result to get keys to merge:')
    sample_result = utils.read_data(filename_list[0])
    jetR = list(sample_result.keys())[0]
    observables = list(sample_result[jetR].keys())
    print(f'Merging: ')
    print(f'  jetR: {jetR}')
    print(f'  observables: {observables}')
    print()

    # Do a preprocessing stage of merging, to condense down to ~100 files
    # This avoids bumping into the ulimit, and also is probably more efficient
    filename_list_stage0 =  merge_stage0(filename_list, observables, jetR, output_dir)

    # Now will create a virtual dataset for each observable
    # See: https://docs.h5py.org/en/stable/vds.html
    
    # First, we need to find the total shape for each observable set
    print('Determining shapes...')
    shapes, total_shape = determine_shapes(jetR, observables, filename_list_stage0, n_files_max, n_files_total)
    print(f'Done determining shapes: n_jets = {total_shape}')
    print()

    # Now, create the virtual dataset
    # We use keys equal to the keys in the input file
    output_filename = 'processing_output_merged.h5'
    with h5py.File(os.path.join(output_dir, output_filename), 'w') as hf:
        
        for observable in observables:

            print(f'Creating virtual dataset for {observable} (total shape = {total_shape})')
            layout = h5py.VirtualLayout(shape=total_shape, dtype=np.float64)

            # Loop through file list
            layout_index = 0
            for i,filename in enumerate(filename_list_stage0):
                if not accept_file(i, n_files_max, n_files_total, log=False):
                    break

                # Create virtual source
                #print(f'  Creating virtual source for file {i} with shape {shapes[i]}')
                source = h5py.VirtualSource(filename, observable, shapes[i], dtype=np.float64)

                # Insert the source into the layout
                new_layout_index = layout_index + shapes[i][0]
                layout[layout_index:new_layout_index] = source
                
                layout_index = new_layout_index
                #print(f'new layout_index: {layout_index}')
                
            # Add virtual dataset to output file
            hf.create_virtual_dataset(observable, layout)

    print('Virtual dataset created.')
    print()

    # Check that we can read in the virtual dataset successfully.
    # Note: This is mainly to protect against the subtle issue with  
    #       HDF5 virtual datasets that the files are not closed and one can therefore  
    #       hit the ulimit which results in some empty data without any error.
    print('Checking that we can read the virtual dataset...')
    with h5py.File(os.path.join(output_dir, output_filename), 'r') as hf:
        result_total = hf['jet_pt'][:]
        print(f'  n_jets: {result_total.shape[0]}')
        print(f'  mean pt: {np.mean(result_total)}')
        print()
        index=0
        for i,filename in enumerate(filename_list_stage0):
            result = result_total[index:index+shapes[i][0]]
            if np.mean(result) < 5.:
                print(f'mean={np.mean(result)}')
                print(f'index={i}')
                print(f'{filename}')
                print(f'{result.shape}')
                print(f'{result}')
                sys.exit('ERROR -- check the max number of files open (ulimit -Hn)')
            index += shapes[i][0]
    print('Done!')
    print()

#---------------------------------------------------------------
def merge_stage0(filename_list, observables, jetR, output_dir, n_chunks=100):
    '''
    Merge the processing output into a smaller number of files.
    This avoids bumping into the ulimit, and also is probably more efficient.
    '''
    print('Performing Stage0 merge...')

    # Create directory (and remove it if it already exists)
    outputdir_stage0 = os.path.join(output_dir, 'Stage0')
    if os.path.exists(outputdir_stage0):
        shutil.rmtree(outputdir_stage0)
    os.makedirs(outputdir_stage0)        

    # Split the filelist up into chunks
    n_files = len(filename_list)
    n_files_per_chunk = int(n_files/n_chunks)+1
    filename_chunk_list = [filename_list[x:x+n_files_per_chunk] for x in range(0,n_files,n_files_per_chunk)]    

    # Loop through chunks, and for each chunk concatenate the arrays for each observable
    for i_chunk,filename_chunk in enumerate(filename_chunk_list):
        output_dict = {}
        for filename in filename_chunk:
            with h5py.File(filename, 'r') as hf:
                for observable in observables:
                    if observable in output_dict:
                        output_dict[observable] = np.concatenate((output_dict[observable], hf[f'{jetR}/{observable}'][:]))
                    else:
                        output_dict[observable] = hf[f'{jetR}/{observable}'][:]

        utils.write_data(output_dict, outputdir_stage0, filename = f'AnalysisResults{i_chunk}.h5')

    # Return the new filelist
    filename_list_stage0 = [os.path.join(outputdir_stage0, f'AnalysisResults{i}.h5') for i in range(len(filename_chunk_list))]
    return filename_list_stage0

#---------------------------------------------------------------
def determine_shapes(jetR, observables, filename_list, n_files_max, n_files_total):
    '''
    Determine the shape of the data (i.e. number of jets) in all files.

    Note: we have already verified that each observable has the same n_jets,
    so we can just grab the length of any observable.

    :returns:
      - shapes - list of shape of each file [shape_file0, shape_file1, ...]
      - total_shapes - shape of total number of jets
    '''

    # First, get a list of n_jets for each file
    shapes = []
    for i,filename in enumerate(filename_list):
        with h5py.File(filename,'r') as hdf:
            if not accept_file(i, n_files_max, n_files_total):
                break
            shapes.append(hdf[observables[0]][:].shape)

    # Also return the total n_jets over all files
    total_shape = (np.sum([shape[0] for shape in shapes]),)

    return shapes, total_shape

#---------------------------------------------------------------
def accept_file(i, n_files_max, n_files_total, log=True):
    '''
    Function called for each file while aggregating.
    '''
    if log:
        if i%10 == 0:
            print(f'{i}/{n_files_total}')
    if i > n_files_max:
        return False
    return True

##################################################################
if __name__ == '__main__':

    # Define arguments
    parser = argparse.ArgumentParser(description='Aggregate data')
    parser.add_argument('-c', '--configFile', action='store',
                        type=str, metavar='configFile',
                        default='../../../../config/multifold/pp/process.yaml',
                        help='Path of config file for analysis')
    parser.add_argument('-o', '--outputDir', action='store',
                        type=str, metavar='outputDir',
                        default='/rstorage/alice/AnalysisResults/james/XXXXXX/',
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
    fileList = os.path.join(args.outputDir, 'processing_output_files.txt')
    if not os.path.exists(fileList):
        print('File \"{0}\" does not exist! Exiting!'.format(fileList))
        sys.exit(0)

    # Create list of filenames
    with open(fileList) as f:
        filename_list = [line.rstrip() for line in f.readlines()]

    aggregate(config_file=args.configFile, filename_list=filename_list, output_dir=args.outputDir)
