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
        #print(f'hdf5 keys: {list(hdf.keys())}')
        event_key = 'event_hard_p'
        event_shape = hdf[event_key].shape
        print(f'event_hard_p shape: {event_shape}')

        status_key = 'event_hard_status_hepevt'
        status_shape = hdf[status_key].shape
        print(f'event_hard_status_hepevt shape: {status_shape}')

        particles_hard = np.empty(event_shape)
        status = np.empty(status_shape)

        n_events = hdf[event_key].shape[0]
        for i in range(n_events):
            if i%10 == 0:
                print(f'event {i}/{n_events}')

            particles_hard[i] = np.array(hdf[event_key][i,:,:])
            status[i] = np.array(hdf[status_key][i,:])

    # Write 4-vectors into a dataframe
    df_particles_hard = construct_dataframe(particles_hard, is_hard=True, status=status)
    #df_particles_background = construct_dataframe(particles_background)
    print(f'df_particles_hard: {df_particles_hard.describe()}')
    #print(f'df_particles_background: {df_particles_background.describe()}')

    # Write result to a new hdf5 file
    input_filename = input_file.rsplit('/', 1)[-1].replace('.h5', '')
    output_filename = os.path.join(output_dir, f'skim_{input_filename}.pkl')
    print(f'writing to new file: {output_filename}')
    with open(output_filename, 'wb') as f:
        pickle.dump(df_particles_hard, f)
        #pickle.dump(df_particles_background, f)
    print('done.')
 
# ---------------------------------------------------------
# Construct dataframe of particles from ndarray and status
# ---------------------------------------------------------
def construct_dataframe(particles, is_hard=False, status=None):

    # particles has format: (event, px/py/pz/E, particle index)
    #print(f'particles: {particles}')

    # Change format: event, particle, 4-vector
    dataset = np.transpose(np.array(particles),(0,2,1))
    #print(f'dataset: {dataset}')
    
    # Reorder 4-vector from (px,py,pz,E) to (E,px,py,pz)
    dataset[:,:,[0,1,2,3]] = dataset[:,:,[3,0,1,2]]
    #print(f'dataset: {dataset}')
    #print(f'dataset shape: {dataset}')
    
    # Translate 3D numpy array into a dataframe
    # Define a unique index for each jet
    columns = ['E', 'px', 'py', 'pz']
    df_particles = pd.DataFrame(dataset.reshape(-1, 4), columns=columns)
    df_particles.index = np.repeat(np.arange(dataset.shape[0]), dataset.shape[1]) + 1
    df_particles.index.name = 'event_id'
    #print(f'df_particles: {df_particles}')

    # For the hard particles, we need to select the final-state particles
    if is_hard:
        #print(f'status: {status}')
        #print(f'status.flatten: {status.flatten()}')
        df_particles['status'] = status.flatten()
        df_particles = df_particles.query('status == 1')
        df_particles = df_particles.drop(columns=['status'])
        #print(f'df_particles: {df_particles}')
        
    # Drop entries with nan
    df_particles = df_particles[df_particles['E'].notna()]
    #print(f'df_particles: {df_particles}')
        
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

