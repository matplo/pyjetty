#!/usr/bin/env python3

'''
For generating particle TTrees at final-state parton & hadron level
using the PYTHIA 8 Monte Carlo simulation
Ezra Lesser (elesser@berkeley.edu)
'''

from __future__ import print_function

import ROOT

import copy
import argparse
import os

from pyjetty.mputils import *

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

from pyjetty.alice_analysis.process.user.gen import parton_hadron_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)
# Automatically set Sumw2 when creating new histograms
#ROOT.TH1.SetDefaultSumw2()


################################################################
class PythiaPartonHadron(parton_hadron_base.PartonHadronBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, output_dir='', tree_output_fname='', no_ALICE_eta_cut=False,
                 args=None, **kwargs):

        super(PythiaPartonHadron, self).__init__(
            output_dir, tree_output_fname, no_ALICE_eta_cut, args, **kwargs)

        self.nev = args.nev
        self.user_seed = args.user_seed


    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def PythiaPartonHadron(self, args):

        # Create ROOT TTree file for storing raw PYTHIA particle information
        outf_path = os.path.join(self.output_dir, args.tree_output_fname)
        outf = ROOT.TFile(outf_path, 'recreate')
        outf.cd()

        pinfo('user seed for pythia', self.user_seed)
        # mycfg = ['PhaseSpace:pThatMin = 100']
        mycfg = ['Random:setSeed=on', 'Random:seed={}'.format(self.user_seed)]
        mycfg.append('HadronLevel:all=off')

        ####################################################################
        # PYTHIA instance with MPI off
        setattr(args, "py_noMPI", True)
        pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
        self.calculate_events(pythia)
        pythia.stat()
        print()

        ####################################################################
        # PYTHIA instance with MPI on
        setattr(args, "py_noMPI", False)
        pythia_MPI = pyconf.create_and_init_pythia_from_args(args, mycfg)
        self.calculate_events(pythia_MPI, MPIon=True)
        pythia_MPI.stat()
        print()

        ####################################################################
        # Print out some final statistics
        print("N total final MPI-off events:", int(self.N_events_MPIoff), "with",
              int(pythia.info.nAccepted() - self.N_events_MPIoff),
              "events rejected at hadronization step",

              "\nN total final MPI-on events:", int(self.N_events_MPIon), "with",
              int(pythia_MPI.info.nAccepted() - self.N_events_MPIon),
              "events rejected at hadronization step")

        ####################################################################
        # Save output trees in ROOT file
        self.fill_write_trees()
        self.save_xsec_N(pythia.info.sigmaGen(), self.N_events_MPIoff,
                         pythia_MPI.info.sigmaGen(), self.N_events_MPIon)

        outf.Write()
        outf.Close()


    #---------------------------------------------------------------
    # Calculate events and pass information on to jet finding
    #---------------------------------------------------------------
    def calculate_events(self, pythia, MPIon=False):
        
        iev = 0  # Event loop count

        while iev < self.nev:
            if not pythia.next():
                continue

            # Creates std::vector<fastjet::PseudoJet> of final-state partons
            parts_pythia_p = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)
            
            hstatus = pythia.forceHadronLevel()
            if not hstatus:
                #pwarning('forceHadronLevel false event', iev)
                continue

            # Creates std::vector<fastjet::PseudoJet> of final-state hadrons
            parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)

            # h_is_charged is a std::vector<bool>
            h_is_charged = pythiafjext.is_charged(pythia, [pythiafjext.kFinal], 0, True)

            self.fill_branches(parts_pythia_p, parts_pythia_h, h_is_charged, iev, self.user_seed, MPIon)

            # Some "accepted" events don't survive hadronization step -- keep track here
            iev += 1

        setattr(self, "N_events_MPI%s" % ("on" if MPIon else "off"), iev)


################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly',
                                     prog=os.path.basename(__file__))
    pyconf.add_standard_pythia_args(parser)
    # Could use --py-seed
    parser.add_argument('--user-seed', help='PYTHIA starting seed', default=1111, type=int)
    parser.add_argument('-o', '--output-dir', action='store', type=str, default='./', 
                        help='Output directory for generated ROOT file(s)')
    parser.add_argument('--tree-output-fname', default="AnalysisResults.root", type=str,
                        help="Filename for the generated particle ROOT TTree")
    parser.add_argument('--no-ALICE-eta-cut', action="store_true",
                        help="Disable track cut for eta<0.9 at both parton and hadron level") 
    args = parser.parse_args()

    # Use PYTHIA seed for event generation
    if args.user_seed < 0:
        args.user_seed = 1111

    # Have at least 1 event
    if args.nev < 1:
        args.nev = 1

    if args.py_noMPI:
        print("\033[91m%s\033[00m" % "WARNING: py-noMPI flag ignored for this program")
        time.sleep(3)
        print()

    process = PythiaPartonHadron(
        output_dir=args.output_dir, tree_output_fname=args.tree_output_fname,
        no_ALICE_eta_cut=args.no_ALICE_eta_cut, args=args)
    process.PythiaPartonHadron(args)
