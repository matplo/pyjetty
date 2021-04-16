# This file takes the data generated from Yue Shi (many h5 files) and reformats it such that it can be used in the ML code
# The format at the end is similar to the quark/gluon test data set
# The data set at the end containes labeled pp and AA jets
# There will be an equal number of pp and AA jets and the data set is shuffled already
# The final format is (jets, particles, 4-vectors), where the 4-vector format is (px,py,pz,E)

# Yue Shi's files contain 5988 Pythia h5 files and 2176 Jewel h5 files
# Note: The number of jets in each file is different and there are generally more jets in Pythia files

import os
import argparse
import h5py
import numpy as np

def skim(input_file, output_dir):

    # ---------------------------------------------------------
    # First: read in Pythia data files and do some reformatting
    # ---------------------------------------------------------
    with h5py.File(input_file, 'r') as hdf:
            
        # jetsubs_cons are the 4-vectors of particles inside the jets
        dataset0 = hdf.get('jetsubs_cons')

        # Change format: jet, particle, 4-vector
        dataset = np.transpose(np.array(dataset0),(0,2,1)) 
            
        # Reorder 4-vector from (px,py,pz,E) to (E,px,py,pz)
        dataset[:,:,[0,1,2,3]] = dataset[:,:,[3,0,1,2]]
            
        # Replace nan with zero
        where_are_nan = np.isnan(dataset)
        dataset[where_are_nan] = 0.
            
        # Zero pad such that they all have the same shape (.., 800, 4)
        data_final = np.pad(dataset,((0,0),(0,800-dataset.shape[1]),(0,0)),'constant',constant_values=(0))
        print(f'final data shape: {data_final.shape}')

    # Create labels: Pythia 0, Jewel 1
    if 'pythia8' in input_file:
        y = np.zeros(data_final.shape[0])
    elif 'jewel' in input_file:
        y = np.ones(data_final.shape[0])
    print(f'label shape: {y.shape[0]}, value: {y[0]}')

    # -----------------------------------------------------------------------------
    # Write labeled, shuffled and equal sample size data to a single h5 file
    # -----------------------------------------------------------------------------

    # Shuffle data set # check again!
    idx = np.random.permutation(len(y))
    data_ppAA_jets = data_final[idx]

    # Write result to a new hdf5 file
    input_filename = input_file.rsplit('/', 1)[-1]
    output_filename = os.path.join(output_dir, f'skim_{input_filename}')
    print(f'writing to new file: {output_filename}')
    with h5py.File(output_filename,'w') as hdf:
        hdf.create_dataset('data',data = data_ppAA_jets)
        hdf.create_dataset('labels',data = y)
    
    print('done.')

##################################################################
if __name__ == '__main__':

    # Define arguments
    parser = argparse.ArgumentParser(description='Skim pp AA')
    parser.add_argument('-f', '--inputFile', action='store',
                        type=str, metavar='inputFile',
                        default='./blah.h5',
                        help='Path of input file for skimming')
    parser.add_argument('-o', '--outputDir', action='store',
                        type=str, metavar='outputDir',
                        default='./TestOutput',
                        help='Output directory for output to be written to')

    # Parse the arguments
    args = parser.parse_args()

    print('Configuring...')
    print('inputFile: \'{0}\''.format(args.inputFile))
    print('ouputDir: \'{0}\"'.format(args.outputDir))

    # If invalid inputFile is given, exit
    if not os.path.exists(args.inputFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
        sys.exit(0)

    skim(input_file=args.inputFile, output_dir=args.outputDir)

