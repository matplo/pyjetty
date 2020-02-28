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

        # Build truth-level histogram of track pT multiplicity
        print("Building truth-level track pT histogram...")
        self.hist_list.append( ("truth_pt", self.build_pt_hist()) )
        print('--- {} seconds ---'.format(time.time() - start_time))

        # Apply eta cut at the end of the TPC
        self.df_fjparticles = self.apply_eta_cut(self.df_fjparticles)
        print('--- {} seconds ---'.format(time.time() - start_time))

        # Apply efficiency cut
        self.df_fjparticles = self.apply_eff_cut(self.df_fjparticles)
        print('--- {} seconds ---'.format(time.time() - start_time))

        # Build truth-level histogram of track pT multiplicity after efficiency cuts
        print("Building truth-level track pT histogram after efficiency cuts...")
        self.hist_list.append( ("truth_pt_eff_cuts", self.build_pt_hist()) )
        print('--- {} seconds ---'.format(time.time() - start_time))

        # Apply pT smearing
        self.df_fjparticles = self.apply_pt_smear(self.df_fjparticles)
        print('--- {} seconds ---'.format(time.time() - start_time))

        # Build truth-level histogram of track pT multiplicity
        print("Building detector-level track pT histogram...")
        self.hist_list.append( ("fastsim_pt", self.build_pt_hist()) )
        print('--- {} seconds ---'.format(time.time() - start_time))

        # ------------------------------------------------------------------------

        # Write data to file
        print(self.df_fjparticles)
        print("Writing fast simulation to ROOT TTree...")
        self.io.save_dataframe("AnalysisResultsFastSim.root", self.df_fjparticles,
                               df_true=True, histograms=self.hist_list)
        print('--- {} seconds ---'.format(time.time() - start_time))


    #---------------------------------------------------------------
    # Initialize dataframe from input data
    #---------------------------------------------------------------
    def init_df(self):
        # Use IO helper class to convert truth-level ROOT TTree into
        # a SeriesGroupBy object of fastjet particles per event
        self.io = process_io.process_io(input_file=self.input_file, output_dir=self.output_dir,
                                        tree_dir='PWGHF_TreeCreator',
                                        track_tree_name='tree_Particle_gen')
        self.df_fjparticles = self.io.load_dataframe()
        self.nTracks_truth = len(self.df_fjparticles)
        print("DataFrame loaded from data.")
        print(self.df_fjparticles)

        # Initialize a list of histograms to be written to file
        self.hist_list = []

    #---------------------------------------------------------------
    # Apply eta cuts
    #---------------------------------------------------------------
    def apply_eta_cut(self, df):
        df = df[df["ParticleEta"].map(abs) < 0.9]
        print("%i out of %i total truth tracks deleted after eta cut." % \
              (self.nTracks_truth - len(df), self.nTracks_truth))
        return df

    #---------------------------------------------------------------
    # Build histogram of pT values and return it
    #---------------------------------------------------------------
    def build_pt_hist(self):
        bins = np.concatenate((np.arange(0, 5, 0.1), np.arange(5, 20, 1), np.arange(20, 52, 2)))
        return np.histogram(self.df_fjparticles["ParticlePt"], bins=bins)

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
        true_pt = df["ParticlePt"]
        smeared_pt = [ np.random.normal(pt, sigma_pt(pt)) for pt in true_pt ]
        df = pd.DataFrame({"run_number": df["run_number"], "ev_id": df["ev_id"], 
                           "ParticlePt": smeared_pt, "ParticleEta": df["ParticleEta"], 
                           "ParticlePhi": df["ParticlePhi"], "z_vtx_reco": df["z_vtx_reco"],
                           "is_ev_rej": df["is_ev_rej"]})
        print("pT has been smeared for all tracks.")

        # Create histogram to verify pt smearing distribution
        pt_dif = [ (pt_smear - pt_true) / pt_true for pt_smear, pt_true in zip(smeared_pt, true_pt) ]
        pt_bins = np.concatenate((np.arange(0, 1, 0.1), np.arange(1, 10, .5), np.arange(10, 20, 1),
                                  np.arange(20, 50, 2), np.arange(50, 95, 5)))
        dif_bins = np.arange(0, 0.5, .001)
        pt_smearing_dists = np.histogram2d(true_pt, pt_dif, bins=[pt_bins, dif_bins])
        self.hist_list.append( ("pt_smearing", pt_smearing_dists) )

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
