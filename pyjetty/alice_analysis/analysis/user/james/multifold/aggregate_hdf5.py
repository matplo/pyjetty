# This script aggregates the hdf5 processing output into a single file.

import os
import sys
import yaml
import argparse
import h5py
import numpy as np
import shutil
from collections import defaultdict

sys.path.append('../../../../../..')
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs
utils = analysis_utils_obs.AnalysisUtils_Obs()

#---------------------------------------------------------------
def aggregate(config_file, filename_list, is_mc, pt_hat_cross_section_dict, output_dir):
    '''
    Main

    We merge all observables present in the processing output here.
    At analysis time, one can then select a subset of observables to analyze.

    Note: We only consider a single jet radius for now.
    '''

    n_files_max = 100000                       # Option to only aggregate n_files_max files
    n_files_total = len(filename_list)

    # Read config file and construct list of observable keys
    # in order to only aggregate the observables we will analyze
    with open(config_file, 'r') as stream:
        config = yaml.safe_load(stream)
        analysis_observable_subconfigs = config['analysis_observables']
        observables = list(analysis_observable_subconfigs.keys())
        observable_keys = []
        for observable in observables:
            obs_subconfig_list = analysis_observable_subconfigs[observable]
            obs_config_dict = dict((key,config[observable][key]) for key in obs_subconfig_list)
            obs_settings = utils.obs_settings(observable, obs_config_dict, obs_subconfig_list)
            obs_grooming_settings = utils.grooming_settings(obs_config_dict)    
            obs_labels = [utils.obs_label(obs_settings[i], obs_grooming_settings[i]) for i in range(len(obs_subconfig_list))]
            obs_keys = [f'{observable}_{obs_label}' if obs_label else observable for obs_label in obs_labels]
            observable_keys += obs_keys
        jetR = str(config['jetR'][0])

    if is_mc:
        levels = ['det_matched', 'truth_matched']
        #jetR = list(sample_result.keys())[0]
        #print('Loading sample result to get keys to merge:')
        #sample_result = utils.read_data(filename_list[0])
        #levels = [level for level in sample_result[jetR].keys() if len(list(sample_result[jetR][level].values()))>0]
    else:
        levels = ['data']
    print(f'Merging: ')
    print(f'  jetR: {jetR}')
    print(f'  levels: {levels}')
    print(f'  observables: {observable_keys}')
    print()
    if n_files_max < n_files_total:
        print(f'Only aggregating {n_files_max} out of {n_files_total} files')
        print()

    # Do a preprocessing stage of merging
    # This avoids bumping into the ulimit, and also is probably more efficient
    # For data: Condense down to ~100 files
    # For MC: Condense each pt-hat bin to 1 file
    filename_list_stage0 =  merge_stage0(filename_list, levels, observable_keys, jetR, 
                                         is_mc, pt_hat_cross_section_dict, 
                                         n_files_max, output_dir)

    # Now will create a virtual dataset for each observable
    # See: https://docs.h5py.org/en/stable/vds.html

    # First, we need to find the total shape for each observable set
    print('Determining shapes...')
    shapes, total_shape = determine_shapes(levels, observable_keys, filename_list_stage0)
    print(f'Done determining shapes:')
    for level in levels:
        print(f'  n_jets = {total_shape[level]} ({level})')
    print()

    # Now, create the virtual dataset
    # We use output keys equal to the keys in the input file
    output_filename = 'processing_output_merged.h5'
    with h5py.File(os.path.join(output_dir, output_filename), 'w') as hf:
        
        if is_mc:
            observable_keys += ['pt_hat_scale_factors']
        
        for level in levels:
            for observable in observable_keys:

                print(f'Creating virtual dataset for {level} {observable} (total shape = {total_shape[level]})')
                layout = h5py.VirtualLayout(shape=total_shape[level], dtype=np.float64)

                # Loop through file list
                layout_index = 0
                for i,filename in enumerate(filename_list_stage0):

                    # Create virtual source
                    #print(f'  Creating virtual source for file {i} with shape {shapes[i]}')
                    source = h5py.VirtualSource(filename, f'{level}/{observable}', shapes[level][i], dtype=np.float64)

                    # Insert the source into the layout
                    new_layout_index = layout_index + shapes[level][i][0]
                    layout[layout_index:new_layout_index] = source
                    
                    layout_index = new_layout_index
                    #print(f'new layout_index: {layout_index}')
                    
                # Add virtual dataset to output file
                hf.create_virtual_dataset(f'{level}/{observable}', layout)

    print('Virtual dataset created.')
    print()

    # Check that we can read in the virtual dataset successfully.
    # Note: This is mainly to protect against the subtle issue with  
    #       HDF5 virtual datasets that the files are not closed and one can therefore  
    #       hit the ulimit which results in some empty data without any error.
    print('Checking that we can read the virtual dataset with standard file opening...')
    with h5py.File(os.path.join(output_dir, output_filename), 'r') as hf:
        for level in levels:
            if 'jet_pt' in hf[level].keys():
                result_total = hf[level]['jet_pt'][:]
                print(f'  {level}')
                print(f'  n_jets: {result_total.shape[0]}')
                print(f'  mean pt: {np.mean(result_total)}')
                if is_mc:
                    scale_factor = hf[level]['pt_hat_scale_factors'][:]
                    print(f'  pt-hat scale factors: {scale_factor}')
                print()
                index=0
                for i,filename in enumerate(filename_list_stage0):
                    result = result_total[index:index+shapes[level][i][0]]
                    if np.mean(result) < 5.:
                        print(f'mean={np.mean(result)}')
                        print(f'index={i}')
                        print(f'{filename}')
                        print(level)
                        print(f'{result.shape}')
                        print(f'{result}')
                        sys.exit('ERROR -- check the max number of files open (ulimit -Hn)')
                    index += shapes[level][i][0]

    print('Checking that we can read the virtual dataset with silx...')
    utils.read_data(os.path.join(output_dir, output_filename))

    print('Done!')
    print()

#---------------------------------------------------------------
def merge_stage0(filename_list, levels, observable_keys, jetR, is_mc, pt_hat_cross_section_dict, n_files_max, output_dir, n_chunks=100):
    '''
    Merge the processing output into a smaller number of files.
    This avoids bumping into the ulimit, and also is probably more efficient.

    For data: Condense down to n_chunks files
    For MC: Condense each pt-hat bin to 1 file
    '''
    print('Performing Stage0 merge...')

    # Create directory (and remove it if it already exists)
    outputdir_stage0 = os.path.join(output_dir, 'Stage0')
    if os.path.exists(outputdir_stage0):
        shutil.rmtree(outputdir_stage0)
    os.makedirs(outputdir_stage0)    

    if is_mc:

        # First, create a dict that separates the filelist by pt-hat bin
        file_dict = defaultdict(list)
        for filename in filename_list:
            #pt_hat_bin = int(filename.split('/')[9])
            pt_hat_bin = int(filename.split('/')[-4])
            if pt_hat_bin not in range(1,21):
                raise ValueError(f'Unexpected pt-ht bin {pt_hat_bin} -- check parsing of path')
            file_dict[pt_hat_bin].append(filename)

        # Check that we found 20 pt hat bins
        pt_hat_bins = sorted(list(file_dict.keys()))
        n_pt_hat_bins = len(pt_hat_bins)
        print(f'Found {n_pt_hat_bins} pt-hat bins: {pt_hat_bins}')
        print()
        if n_pt_hat_bins != 20:
            raise KeyError(f'We only found {n_pt_hat_bins}!')

        # Get the total number of events
        print('Computing total number of events...')
        n_events_total = 0
        for filename in filename_list[:n_files_max]:
            with h5py.File(filename, 'r') as hf:
                n_events_total += hf['n_events'][()]
        n_events_avg = n_events_total / n_pt_hat_bins
        print(f'n_events_total: {n_events_total}')
        print(f'n_events_avg per pt-hat bin: {n_events_avg}')
        print()

        # Then concatenate and write the arrays for each pt-hat bin
        for pt_hat_bin in pt_hat_bins:
            print(f'  pt-hat bin {pt_hat_bin}: {len(file_dict[pt_hat_bin])} files')
            n_events_pt_hat_bin = 0
            output_dict = {}            
            for level in levels:
                output_dict[level] = {}

            for filename in file_dict[pt_hat_bin][:int(n_files_max/n_pt_hat_bins)]:
                with h5py.File(filename, 'r') as hf:
                    for level in levels:
                        n_events_pt_hat_bin += hf['n_events'][()]
                        for observable in observable_keys:
                            if observable in output_dict[level].keys():
                                output_dict[level][observable] = np.concatenate((output_dict[level][observable], hf[f'{jetR}/{level}/{observable}'][:]))
                            else:
                                output_dict[level][observable] = hf[f'{jetR}/{level}/{observable}'][:]

            # Now that we have computed n_events_pt_hat_bin, write pt_hat_scale_factor = (pt-hat cross-section)/(event_scale_factor)
            observable = 'pt_hat_scale_factors'
            pt_hat_scale_factors = pt_hat_cross_section_dict[pt_hat_bin] / (n_events_pt_hat_bin / n_events_avg)
            for filename in file_dict[pt_hat_bin][:int(n_files_max/n_pt_hat_bins)]:
                with h5py.File(filename, 'r') as hf:
                    for level in levels:
                        n_jets = next(iter(hf[jetR][level].values())).shape[0]
                        if observable in output_dict[level].keys():
                            output_dict[level][observable] = np.concatenate((output_dict[level][observable], np.repeat(pt_hat_scale_factors, n_jets)))
                        else:
                            output_dict[level][observable] = np.repeat(pt_hat_scale_factors, n_jets)
            
            utils.write_data(output_dict, outputdir_stage0, filename = f'MergedResults{pt_hat_bin}.h5')

            #print('Checking that we can read the intermediate dataset...')
            #utils.read_data(os.path.join(outputdir_stage0, f'MergedResults{pt_hat_bin}.h5'))

        filename_list_stage0 = [os.path.join(outputdir_stage0, f'MergedResults{i}.h5') for i in range(1, n_pt_hat_bins+1)]

    else:

        # Split the filelist up into chunks
        filename_list = filename_list[:n_files_max]
        n_files = len(filename_list)
        n_files_per_chunk = int(n_files/n_chunks)+1
        filename_chunk_list = [filename_list[x:x+n_files_per_chunk] for x in range(0,n_files,n_files_per_chunk)]    

        # Loop through chunks, and for each chunk concatenate the arrays for each observable
        for i_chunk,filename_chunk in enumerate(filename_chunk_list):
            output_dict = {}
            for level in levels:
                output_dict[level] = {}

            for filename in filename_chunk:
                with h5py.File(filename, 'r') as hf:
                    for level in levels:
                        for observable in observable_keys:
                            if observable in output_dict[level].keys():
                                output_dict[level][observable] = np.concatenate((output_dict[level][observable], hf[f'{jetR}/{level}/{observable}'][:]))
                            else:
                                output_dict[level][observable] = hf[f'{jetR}/{level}/{observable}'][:]

            utils.write_data(output_dict, outputdir_stage0, filename = f'MergedResults{i_chunk}.h5')

        filename_list_stage0 = [os.path.join(outputdir_stage0, f'MergedResults{i}.h5') for i in range(len(filename_chunk_list))]

    # Return the new filelist
    return filename_list_stage0

#---------------------------------------------------------------
def determine_shapes(levels, observable_keys, filename_list):
    '''
    Determine the shape of the data (i.e. number of jets) in all files.

    Note: we have already verified that each observable has the same n_jets,
    so we can just grab the length of any observable.

    :returns:
      - shapes - list of shape of each file [shape_file0, shape_file1, ...]
      - total_shapes - shape of total number of jets
    '''

    # First, get a list of n_jets for each file
    shapes = defaultdict(list)
    for filename in filename_list:
        with h5py.File(filename,'r') as hdf:
            for level in levels:
                shapes[level].append(hdf[level][observable_keys[0]][:].shape)

    # Also return the total n_jets over all files
    total_shape = {}
    for level in levels:
        total_shape[level] = (np.sum([shape[0] for shape in shapes[level]]),)

    return shapes, total_shape

#---------------------------------------------------------------
def check_if_mc(file):
    '''
    Function to check whether the production is data or MC
    '''
    if 'LHC18b8' in file:
        pt_hat_cross_section_filename = '/global/cfs/cdirs/alice/jdmull/multifold/scale_factors/LHC18b8/scaleFactors.yaml'
        #pt_hat_cross_section_filename = '/rstorage/alice/data/LHC18b8/scaleFactors.yaml'
    elif 'LHC18b8' in file:
        pt_hat_cross_section_filename = '/rstorage/alice/data/LHC20g4/scaleFactors.yaml'
    else:
        pt_hat_cross_section_filename = None

    if pt_hat_cross_section_filename:
        is_mc = True
        with open(pt_hat_cross_section_filename, 'r') as stream:
            pt_hat_cross_section_dict = yaml.safe_load(stream)
    else:
        is_mc = False
        pt_hat_cross_section_dict = None

    print(f'is_mc: {is_mc}')
    if is_mc:
        print(f'Using pt-hat scale factors from {pt_hat_cross_section_filename}')
    print()

    return is_mc, pt_hat_cross_section_dict

##################################################################
if __name__ == '__main__':

    # Define arguments
    parser = argparse.ArgumentParser(description='Aggregate data')
    parser.add_argument('-c', '--configFile', action='store',
                         type=str, metavar='configFile',
                         default='../../../../config/multifold/pp/analysis_R04.yaml',
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

    # Check if data or MC, and get the pt-hat scale factors if MC
    is_mc, pt_hat_cross_section_dict = check_if_mc(filename_list[0])

    aggregate(config_file=args.configFile,
              filename_list=filename_list, 
              is_mc=is_mc, 
              pt_hat_cross_section_dict=pt_hat_cross_section_dict, 
              output_dir=args.outputDir)