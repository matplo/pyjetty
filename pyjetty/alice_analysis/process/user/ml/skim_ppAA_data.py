'''
This script takes a hdf5 file generated from Yue Shi (containing all particles in event)
and writes out dataframes of particle four-vectors per event,
for both the hard process and the background
'''

import os
import argparse
import h5py
import numpy as np
import pandas as pd
import pickle

# ---------------------------------------------------------
# Main
# ---------------------------------------------------------
def skim(input_file, output_dir):

    # Read in Pythia/Jewel data files and do some reformatting
    with h5py.File(input_file, 'r') as hdf:
        
        # Get four-vectors from file
        particles_hard = np.array(hdf.get('event_hard_p'))
        particles_background = np.array(hdf.get('event_ue_p'))
        status = np.array(hdf.get('event_hard_status_hepevt'))
        
        # Write 4-vectors into a dataframe
        df_particles_hard = construct_dataframe(particles_hard, is_hard=True, status=status)
        df_particles_background = construct_dataframe(particles_background)
        print(f'df_particles_hard: {df_particles_hard.describe()}')
        print(f'df_particles_background: {df_particles_background.describe()}')

    # Write result to a new hdf5 file
    input_filename = input_file.rsplit('/', 1)[-1].replace('.h5', '')
    output_filename = os.path.join(output_dir, f'skim_{input_filename}.pkl')
    print(f'writing to new file: {output_filename}')
    with open(output_filename, 'wb') as f:
        pickle.dump(df_particles_hard, f)
        pickle.dump(df_particles_background, f)
    print('done.')
 
# ---------------------------------------------------------
# Construct dataframe of particles from ndarray and status
# ---------------------------------------------------------
def construct_dataframe(particles, is_hard=False, status=None):

    # particles has format: (event, px/py/pz/E, particle index)
    #print(particles)
    
    # Change format: event, particle, 4-vector
    dataset = np.transpose(np.array(particles),(0,2,1))
    #print(dataset)
    
    # Reorder 4-vector from (px,py,pz,E) to (E,px,py,pz)
    dataset[:,:,[0,1,2,3]] = dataset[:,:,[3,0,1,2]]
    #print(dataset.shape)
    
    # Replace nan with zero
    where_are_nan = np.isnan(dataset)
    #print(f'where_are_nan: {where_are_nan}')
    #print(dataset)
    
    # Translate 3D numpy array into a dataframe
    # Define a unique index for each jet
    columns = ['E', 'px', 'py', 'pz']
    df_particles = pd.DataFrame(dataset.reshape(-1, 4), columns=columns)
    df_particles.index = np.repeat(np.arange(dataset.shape[0]), dataset.shape[1]) + 1
    df_particles.index.name = 'event_id'
    #print(df_particles)
    
    # For the hard particles, we need to select the final-state particles
    if is_hard:
        #print(status)
        #print(status.flatten())
        df_particles['status'] = status.flatten()
        df_particles = df_particles.query('status == 1')
        df_particles = df_particles.drop(columns=['status'])
        #print(df_particles)
        
    return df_particles

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
        
    # If output dir does not exist, create it
    if not os.path.exists(args.outputDir):
        os.makedirs(args.outputDir)

    skim(input_file=args.inputFile, output_dir=args.outputDir)

