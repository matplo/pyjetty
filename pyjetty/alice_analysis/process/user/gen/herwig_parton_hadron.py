#!/usr/bin/env python

'''
For generating particle TTrees at final-state parton & hadron level
using pre-created Herwig 7 Monte Carlo simulation output files
Ezra Lesser (elesser@berkeley.edu)
'''

from __future__ import print_function

import ROOT

import copy
import argparse
import os

from pyjetty.mputils import *

from pyjetty.alice_analysis.process.user.gen import parton_hadron_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)
# Automatically set Sumw2 when creating new histograms
#ROOT.TH1.SetDefaultSumw2()

################################################################
class HerwigPartonHadron(parton_hadron_base.PartonHadronBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, output_dir='', tree_output_fname='', no_ALICE_eta_cut=False,
                 args=None, **kwargs):

        super(HerwigPartonHadron, self).__init__(
            output_dir, tree_output_fname, no_ALICE_eta_cut, args, **kwargs)

        self.initialize_config(args)


    #---------------------------------------------------------------
    # Initialize config file into class members
    #---------------------------------------------------------------
    def initialize_config(self, args):

        self.herwig_file = args.input_file
        self.herwig_file_MPI = args.input_file_mpi

        self.run_number = args.run_number

        # Initialize event counter
        self.N_events_MPIon = 0
        self.N_events_MPIoff = 0

        # Initialize variables for final cross sections from event generator
        self.xsec = None
        self.xsec_MPI = None


    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def HerwigPartonHadron(self, args):

        # Create ROOT TTree file for storing raw PYTHIA particle information
        outf_path = os.path.join(self.output_dir, args.tree_output_fname)
        outf = ROOT.TFile(outf_path, 'recreate')
        outf.cd()

        self.parse_events()
        self.parse_events(MPIon=True)

        ####################################################################
        # Save output trees in ROOT file
        self.fill_write_trees()
        self.save_xsec_N(self.xsec, self.N_events_MPIoff,
                         self.xsec_MPI, self.N_events_MPIon)

        outf.Write()
        outf.Close()


    #---------------------------------------------------------------
    # Read events from output, find jets, and fill histograms
    #---------------------------------------------------------------
    def parse_events(self, MPIon=False):

        if MPIon:
            nev_name = "N_events_MPIon"
            infile = self.herwig_file_MPI
        else:
            nev_name = "N_events_MPIoff"
            infile = self.herwig_file

        print("Reading events from %s..." % infile)

        with open(infile, 'r') as f:
            ev_num = 0

            # Flags to assist with keeping track of place within file
            reading_ev = False
            parton = False
            parton_final = False
            parton_finished = False
            hadron = False
            hadron_final = False

            partons_px = []
            partons_py = []
            partons_pz = []
            partons_e = []
            #partons_q = []
            hadrons_px = []
            hadrons_py = []
            hadrons_pz = []
            hadrons_e = []
            #hadrons_q = []
            h_is_charged = []

            for line in f:

                # Waiting to start reading event
                if not reading_ev:
                    if "Event number" in line:
                        reading_ev = True
                        ev_num = int(line.split()[2])
                        if not ev_num % 1000:
                            print("Event number", ev_num, end="\r")
                        setattr(self, nev_name, getattr(self, nev_name)+1)
                    elif "Total integrated xsec:" in line:
                        if MPIon:
                            self.xsec_MPI = float(line.split()[3])
                        else:
                            self.xsec = float(line.split()[3])
                    continue
                
                # Reading event
                # First step is to read the parton info
                elif not parton:
                    if "ShowerHandler" in line:
                        parton = True
                    continue
                elif not parton_final:
                    if "final" in line:
                        parton_final = True
                    continue
                
                # Get showered partons
                elif not parton_finished:
                    # Read parton information
                    vals = line.split()
                    if line[0] == '-':
                        parton_finished = True
                    elif len(vals) == 5 and line[2] == ' ':
                        partons_px.append(vals[0])
                        partons_py.append(vals[1])
                        partons_pz.append(vals[2])
                        partons_e.append(vals[3])
                        #partons_q.append(vals[4])
                    continue

                # Get final hadrons
                elif not hadron:
                    if "DecayHandler" in line:
                        hadron = True
                    continue
                elif not hadron_final:
                    if "final" in line:
                        hadron_final = True
                    continue

                # Check if event is over
                elif line[0] == '-':
                    # Finished reading hadron info
                    reading_ev = False
                    parton = False
                    parton_final = False
                    parton_finished = False
                    hadron = False
                    hadron_final = False
                    
                    # Get correct structure for finding jets
                    partons = fjext.vectorize_px_py_pz_e(
                        partons_px, partons_py, partons_pz, partons_e)

                    hadrons = fjext.vectorize_px_py_pz_e(
                        hadrons_px, hadrons_py, hadrons_pz, hadrons_e)

                    self.fill_branches(
                        partons, hadrons, h_is_charged, ev_num, self.run_number, MPIon)

                    partons_px = []
                    partons_py = []
                    partons_pz = []
                    partons_e = []
                    #partons_q = []
                    hadrons_px = []
                    hadrons_py = []
                    hadrons_pz = []
                    hadrons_e = []
                    #hadrons_q = []
                    h_is_charged = []

                    continue

                elif line[2].isnumeric():
                    # Save the hadron name for charge descrimination
                    i = 1
                    while line[i].isnumeric():
                        i += 1
                    hadron_type = line[i:].split()[0]
                    continue

                elif line[2] == ' ':  # and len(line.split()) == 5:
                    # Reading hadron information
                    vals = line.split()
                    hadrons_px.append(vals[0])
                    hadrons_py.append(vals[1])
                    hadrons_pz.append(vals[2])
                    hadrons_e.append(vals[3])
                    #hadrons_q.append(vals[4])

                    h_is_charged.append('+' in hadron_type or '-' in hadron_type)


################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Herwig7 debug level-1 output parser',
                                     prog=os.path.basename(__file__))
    parser.add_argument('-i', '--input-file', action='store', type=str, default='LHC.log',
                        help='Input .log file from Herwig7 analysis')
    parser.add_argument('-m', '--input-file-mpi', action='store', type=str, default=None,
                        help='Input .log file with MPI on from Herwig7 analysis')
    parser.add_argument('-r', '--run-number', action='store', type=int, default=1111,
                        help='Run number for these TTrees (use the Herwig7 seed value)')
    parser.add_argument('-o', '--output-dir', action='store', type=str, default='./', 
                        help='Output directory for generated ROOT file(s)')
    parser.add_argument('--tree-output-fname', default="AnalysisResults.root", type=str,
                        help="Filename for the (unscaled) generated particle ROOT TTree")
    parser.add_argument('--no-ALICE-eta-cut', action="store_true",
                        help="Disable track cut for eta<0.9 at both parton and hadron level") 
    args = parser.parse_args()

    process = HerwigPartonHadron(
        output_dir=args.output_dir, tree_output_fname=args.tree_output_fname,
        no_ALICE_eta_cut=args.no_ALICE_eta_cut, args=args)
    process.HerwigPartonHadron(args)
