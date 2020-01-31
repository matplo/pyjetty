#!/usr/bin/env python3

'''
Loads truth-level MC track information and applies effciency cuts and pT smearing
to emulate a fast simulation, and saves resulting ROOT files.

Written by Ezra Lesser (elesser@berkeley.edu), Spring 2020
Some code taken from processing script by James Mulligan
'''

from __future__ import print_function, division

import os
import argparse
import time
import pandas as pd
#import fastjet as fj
#import fjext

from pyjetty.alice_analysis.process.base import process_io
from pyjetty.cstoy import alice_efficiency


#################################################################################
class eff_smear:

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, inputFile='', outputDir=''):
        self.input_file = inputFile
        self.output_dir = outputDir

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def eff_smear(self):

        start_time = time.time()

        # ------------------------------------------------------------------------

        # Initialize dataframes from data
        self.init_df()
        print('--- {} seconds ---'.format(time.time() - start_time))

        # ------------------------------------------------------------------------

        # Applying efficiency cuts and pT smearing
        self.apply_eff_smear(self.df_fjparticles)
        print('--- {} seconds ---'.format(time.time() - start_time))

    #---------------------------------------------------------------
    # Initialize dataframe from input data
    #---------------------------------------------------------------
    def init_df(self):
        # Use IO helper class to convert truth-level ROOT TTree into
        # a SeriesGroupBy object of fastjet particles per event
        io_truth = process_io.process_io(input_file=self.input_file, track_tree_name='tree_Particle_gen')
        self.df_fjparticles = io_truth.load_dataframe()
        #self.df_fjparticles = io_truth.load_data(group_by_evid=False)  # Get particle info in fj format (slow)
        self.nTracks_truth = len(self.df_fjparticles)
        

    #---------------------------------------------------------------
    # Apply fast simulation procedure
    #---------------------------------------------------------------
    def apply_eff_smear(self, df):
        print(df, len(df))
        # Apply efficiency cut for fastjet particles
        eff_smearer = alice_efficiency.AliceChargedParticleEfficiency()
        df = df[df["ParticlePt"].map(lambda x: eff_smearer.pass_eff_cut(x))]
        '''
        for fj_particles_truth, fj_particles_fs in \
            zip(self.df_fjparticles['fj_particles_det'], self.df_fjparticles['fj_particles_truth']):

            # Check that the entries exist appropriately
            # (need to check how this can happen -- but it is only a tiny fraction of events)
            if type(fj_particles_det) != fj.vectorPJ or type(fj_particles_truth) != fj.vectorPJ:
                print('fj_particles type mismatch -- skipping event')
                return

        '''
        print(df)
        print(len(df), type(df))


##################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Fast Simulation Generator")
    parser.add_argument("-i", "--inputFile", action="store", type=str, metavar="inputFile",
                        default="AnalysisResults.root", help="Path of ROOT file containing MC TTrees")
    parser.add_argument("-o", "--outputDir", action="store", type=str, metavar="outputDir",
                        default="./TestOutput", help="Output path for fast sim ROOT TTree")
    args = parser.parse_args()

    print('Configuring...')
    print('inputFile: \'{0}\''.format(args.inputFile))
    print('ouputDir: \'{0}\"'.format(args.outputDir))

    # If invalid inputFile is given, exit
    if not os.path.exists(args.inputFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
        sys.exit(0)

    processor = eff_smear(inputFile=args.inputFile, outputDir=args.outputDir)
    processor.eff_smear()
