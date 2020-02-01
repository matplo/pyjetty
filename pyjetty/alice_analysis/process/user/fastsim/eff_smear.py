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
import numpy as np
#import fastjet as fj
#import fjext

from pyjetty.alice_analysis.process.base import process_io
from pyjetty.cstoy import alice_efficiency


# Approximated fit for MC calculated (LHC18b8) pT resolution
def sigma_pt(pt):
    if pt < 1:
        return pt * (-0.035 * pt + 0.04)
    elif pt < 60:
        return pt * (0.00085 * pt + 0.00415)
    # else: approximately valid for at least pt < 90 GeV
    return pt * (0.0015 * pt - 0.035)


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

        # Apply eta cut at the end of the TPC
        self.df_fjparticles = self.apply_eta_cut(self.df_fjparticles)
        print('--- {} seconds ---'.format(time.time() - start_time))

        # Apply efficiency cut
        self.df_fjparticles = self.apply_eff_cut(self.df_fjparticles)
        print('--- {} seconds ---'.format(time.time() - start_time))

        # Apply pT smearing
        self.df_fjparticles = self.apply_pt_smear(self.df_fjparticles)
        print('--- {} seconds ---'.format(time.time() - start_time))

        # ------------------------------------------------------------------------

        print(self.df_fjparticles)


    #---------------------------------------------------------------
    # Initialize dataframe from input data
    #---------------------------------------------------------------
    def init_df(self):
        # Use IO helper class to convert truth-level ROOT TTree into
        # a SeriesGroupBy object of fastjet particles per event
        io_truth = process_io.process_io(input_file=self.input_file, 
                                         track_tree_name='tree_Particle_gen')
        self.df_fjparticles = io_truth.load_dataframe()
        #self.df_fjparticles = io_truth.load_data(group_by_evid=False)  # Get particle info in fj format (slow)
        self.nTracks_truth = len(self.df_fjparticles)
        print("DataFrame loaded from data.")
        print(self.df_fjparticles)
        
    #---------------------------------------------------------------
    # Apply eta cuts
    #---------------------------------------------------------------
    def apply_eta_cut(self, df):
        df = df[df["ParticleEta"].map(abs) < 0.9]
        print("%i out of %i total truth tracks deleted after eta cut." % \
              (self.nTracks_truth - len(df), self.nTracks_truth))
        return df

    #---------------------------------------------------------------
    # Apply efficiency cuts
    #---------------------------------------------------------------
    def apply_eff_cut(self, df):

        # Apply efficiency cut for fastjet particles
        eff_smearer = alice_efficiency.AliceChargedParticleEfficiency()
        df = df[df["ParticlePt"].map(lambda x: eff_smearer.pass_eff_cut(x))]
        print("%i out of %i total truth tracks deleted after efficiency cut." % \
              (self.nTracks_truth - len(df), self.nTracks_truth))
        return df

    #---------------------------------------------------------------
    # Apply pt smearing
    #---------------------------------------------------------------
    def apply_pt_smear(self, df):
        smeared_pt = [ np.random.normal(pt, sigma_pt(pt)) for pt in df["ParticlePt"] ]
        df = pd.DataFrame({"run_number": df["run_number"], "ev_id": df["ev_id"], 
                           "ParticlePt": smeared_pt, "ParticleEta": df["ParticleEta"], 
                           "ParticlePhi": df["ParticlePhi"], "z_vtx_reco": df["z_vtx_reco"],
                           "is_ev_rej": df["is_ev_rej"]})
        print("pT has been smeared for all tracks.")
        return df

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
