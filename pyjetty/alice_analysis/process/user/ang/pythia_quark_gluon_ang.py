#!/usr/bin/env python
'''
Script for looking at the quark vs gluon dependence of substructure observables
Author: Ezra Lesser (elesser@berkeley.edu)
'''

from __future__ import print_function

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext

import ROOT

import tqdm
import yaml
import copy
import argparse
import os
import array
import numpy as np

from pyjetty.mputils import *

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

from pyjetty.alice_analysis.process.base import process_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)
# Automatically set Sumw2 when creating new histograms
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

################################################################
class PythiaQuarkGluon(process_base.ProcessBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, args=None, **kwargs):

        super(PythiaQuarkGluon, self).__init__(
            input_file, config_file, output_dir, debug_level, **kwargs)

        # Call base class initialization
        process_base.ProcessBase.initialize_config(self)

        # Read config file
        with open(self.config_file, 'r') as stream:
            config = yaml.safe_load(stream)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.jetR_list = config["jetR"]

        self.user_seed = args.user_seed
        self.nev = args.nev

        self.noMPI = (bool)(1-args.MPIon)
        self.noISR = (bool)(1-args.ISRon)

        # PDG ID values for quarks and gluons
        self.quark_pdg_ids = [1, 2, 3, 4, 5, 6, 7, 8]
        self.charm_pdg_ids = [4]
        self.gluon_pdg_ids = [9, 21]

        # hadron level - ALICE tracking restriction
        self.max_eta_hadron = 0.9

        self.min_leading_track_pT = config["min_leading_track_pT"] if "min_leading_track_pT" in config else None

        self.pt_bins = array.array('d', list(range(5, 100, 5)) + list(range(100, 210, 10)))
        self.obs_bins_ang = np.concatenate((np.linspace(0, 0.009, 10), np.linspace(0.01, 0.1, 19),
                                            np.linspace(0.11, 0.8, 70)))
        self.obs_bins_mass = np.concatenate(
          (np.linspace(0, 0.9, 10), np.linspace(1, 9.8, 45), np.linspace(10, 14.5, 10),
           np.linspace(15, 19, 5), np.linspace(20, 60, 9)))

        self.observable_list = config['process_observables']
        self.obs_settings = {}
        self.obs_grooming_settings = {}
        self.obs_names = {}
        for observable in self.observable_list:

          obs_config_dict = config[observable]
          obs_config_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]

          obs_subconfig_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
          self.obs_settings[observable] = self.utils.obs_settings(observable, obs_config_dict, obs_subconfig_list)
          self.obs_grooming_settings[observable] = self.utils.grooming_settings(obs_config_dict)

          self.obs_names[observable] = obs_config_dict["common_settings"]["xtitle"]

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def pythia_quark_gluon(self, args):

        # Create ROOT TTree file for storing raw PYTHIA particle information
        outf_path = os.path.join(self.output_dir, args.tree_output_fname)
        outf = ROOT.TFile(outf_path, 'recreate')
        outf.cd()

        # Initialize response histograms
        self.initialize_hist()

        pinfo('user seed for pythia', self.user_seed)
        mycfg = ['Random:setSeed=on', 'Random:seed={}'.format(self.user_seed)]
        mycfg.append('HadronLevel:all=off')

        # print the banner first
        fj.ClusterSequence.print_banner()
        print()

        # -------------------------------
        # Setting MPIs and ISRs
        print('Will run no MPI:',self.noMPI)
        print('Will run no ISR:',self.noISR)
        setattr(args, "py_noMPI", self.noMPI)
        setattr(args, "py_noISR", self.noISR)
        # -------------------------------

        pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)

        self.init_jet_tools()
        self.calculate_events(pythia)
        pythia.stat()
        print()

        self.scale_print_final_info(pythia)

        outf.Write()
        outf.Close()

        self.save_output_objects()

    #---------------------------------------------------------------
    # Initialize histograms
    #---------------------------------------------------------------
    def initialize_hist(self):

        self.hNevents = ROOT.TH1I("hNevents", 'Number accepted events (unscaled)', 2, -0.5, 1.5)

        for jetR in self.jetR_list:

            # Store a list of all the histograms just so that we can rescale them later
            hist_list_name = "hist_list_R%s" % str(jetR).replace('.', '')
            setattr(self, hist_list_name, [])

            R_label = str(jetR).replace('.', '') + 'Scaled'

            for observable in self.observable_list:
                # Should only be two: observable == "ang" or "mass"
                if observable != "ang" and observable != "mass":
                    raise ValueError("Observable %s is not implemented in this script" % observable)

                obs_name = self.obs_names[observable]
                obs_bins = getattr(self, "obs_bins_" + observable)
                # Use more finely binned pT bins for TH2s than for the RMs
                pt_bins = array.array('d', list(range(0, 201, 1)))

                # Loop over subobservable (alpha value)
                for i in range(len(self.obs_settings[observable])):

                    obs_setting = self.obs_settings[observable][i]
                    grooming_setting = self.obs_grooming_settings[observable][i]
                    obs_label = self.utils.obs_label(obs_setting, grooming_setting)

                    for parton_type in ["quark", "charm", "gluon"]:

                        name = ('h_%s_JetPt_%s_R%s_%s' % (observable, parton_type, jetR, obs_label)) if \
                            len(obs_label) else ('h_%s_JetPt_%s_R%s' % (observable, parton_type, jetR))
                        h = ROOT.TH2F(name, name, len(pt_bins)-1, pt_bins, len(obs_bins)-1, obs_bins)
                        h.GetXaxis().SetTitle('#it{p}_{T,%s}^{ch jet}' % (parton_type[0] + "-init"))
                        h.GetYaxis().SetTitle(obs_name + '^{%s}' % (parton_type[0] + "-init"))
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

    #---------------------------------------------------------------
    # Initiate jet defs, selectors, and sd (if required)
    #---------------------------------------------------------------
    def init_jet_tools(self):

        for jetR in self.jetR_list:
            jetR_str = str(jetR).replace('.', '')

            # set up our jet definition and a jet selector
            jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
            setattr(self, "jet_def_R%s" % jetR_str, jet_def)
            print(jet_def)

        pwarning('max eta for particles after hadronization set to', self.max_eta_hadron)
        parts_selector_h = fj.SelectorAbsEtaMax(self.max_eta_hadron)

        for jetR in self.jetR_list:
            jetR_str = str(jetR).replace('.', '')

            jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsEtaMax(self.max_eta_hadron - jetR)
            #jet_selector = fj.SelectorPtMin(0.) & fj.SelectorAbsEtaMax(self.max_eta_hadron - jetR)
            setattr(self, "jet_selector_R%s" % jetR_str, jet_selector)

            count1 = 0  # Number of partonic parents which match to >1 ch-jets
            setattr(self, "count1_R%s" % jetR_str, count1)
            count2 = 0  # Number of partonic parents which match to zero ch-jets
            setattr(self, "count2_R%s" % jetR_str, count2)

    #---------------------------------------------------------------
    # Calculate events and pass information on to jet finding
    #---------------------------------------------------------------
    def calculate_events(self, pythia):

        iev = 0  # Event loop count

        while iev < self.nev:
            if not pythia.next():
                continue

            self.parents = []
            self.event = pythia.event
            #print(self.event)
            fs_parton_5 = fj.PseudoJet(pythia.event[5].px(), pythia.event[5].py(), pythia.event[5].pz(), pythia.event[5].e())
            fs_parton_6 = fj.PseudoJet(pythia.event[6].px(), pythia.event[6].py(), pythia.event[6].pz(), pythia.event[6].e())
            self.parents = [fs_parton_5, fs_parton_6] # parent partons in dijet

            # Save PDG code of the parent partons
            self.parent_ids = [pythia.event[5].id(), pythia.event[6].id()]

            # parton level
            #parts_pythia_p = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)

            hstatus = pythia.forceHadronLevel()
            if not hstatus:
                continue

            # full-hadron level
            #parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)

            # charged-hadron level
            parts_pythia_hch = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, True)

            # Some "accepted" events don't survive hadronization step -- keep track here
            self.hNevents.Fill(0)
            self.find_jets_fill_histograms(parts_pythia_hch, iev)

            iev += 1

    #---------------------------------------------------------------
    # Find primordial parent
    #---------------------------------------------------------------
    def primordial_parent(self,p):
        parent1 = parent2 = -10
        while p > 6:
            parent1 = self.event[p].mother1()
            parent2 = self.event[p].mother2()
            if parent1 != parent2:
                p = max(parent1,parent2)
            else:
                p = parent1
        return p

    #---------------------------------------------------------------
    # Find jets, do matching between levels, and fill histograms
    #---------------------------------------------------------------
    def find_jets_fill_histograms(self, parts_pythia_hch, iev):
        # Loop over jet radii
        for jetR in self.jetR_list:

            jetR_str = str(jetR).replace('.', '')
            jet_selector = getattr(self, "jet_selector_R%s" % jetR_str)
            jet_def = getattr(self, "jet_def_R%s" % jetR_str)

            count1 = getattr(self, "count1_R%s" % jetR_str)
            count2 = getattr(self, "count2_R%s" % jetR_str)

            # Get the jets at different levels
            #jets_p  = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_p  ))) # parton level
            #jets_h  = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h  ))) # full hadron level
            jets_ch = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_hch))) # charged hadron level

            R_label = str(jetR).replace('.', '') + 'Scaled'

            # Find the charged jet closest to the axis of the original parton
            # Require that the match is within some small angle, and that it is unique
            jet_matching_distance = 0.6  # Match jets with deltaR < jet_matching_distance*jetR
            self.parent0match, self.parent1match = None, None
            for i_jch, jch in enumerate(jets_ch):
                # Do constituent pT cut
                if self.min_leading_track_pT and not \
                    self.utils.is_truth_jet_accepted(jch, self.min_leading_track_pT):
                    continue
                for i_parent, parent in enumerate(self.parents):
                    parentmatch_name = "parent%imatch" % i_parent
                    if jch.delta_R(parent) < jet_matching_distance * jetR:
                        match = getattr(self, parentmatch_name)
                        if not match:
                            setattr(self, parentmatch_name, jch)
                        else:  # Already found a match
                            # Set flag value so that we know to ignore this one
                            setattr(self, parentmatch_name, 0)

            # If we have matches, fill histograms
            for i_parent, parent in enumerate(self.parents):
                jet = getattr(self, "parent%imatch" % i_parent)
                if not jet:
                  if jet == 0:
                      # More than one match -- take note and continue
                      count1 += 1
                      continue
                  else:  # jet == None
                      # No matches -- take note and continue
                      count2 += 1
                      continue

                # One unique match
                # Identify the histograms which need to be filled
                parton_id = self.parent_ids[i_parent]
                parton_types = []
                if parton_id in self.quark_pdg_ids:
                    parton_types += ["quark"]
                    if parton_id in self.charm_pdg_ids:
                        parton_types += ["charm"]
                elif parton_id in self.gluon_pdg_ids:
                    parton_types += ["gluon"]

                # If parent parton not identified, skip for now
                if not len(parton_types):
                    continue

                # Fill histograms
                for observable in self.observable_list:
                    for i in range(len(self.obs_settings[observable])):

                        obs_setting = self.obs_settings[observable][i]
                        grooming_setting = self.obs_grooming_settings[observable][i]
                        obs_label = self.utils.obs_label(obs_setting, grooming_setting)

                        # Groom jet, if applicable
                        jet_groomed_lund = None
                        if grooming_setting:
                            gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
                            jet_groomed_lund = self.utils.groom(gshop, grooming_setting, jetR)
                            if not jet_groomed_lund:
                                continue

                        obs = self.calculate_observable(
                            observable, jet, jet_groomed_lund, jetR, obs_setting,
                            grooming_setting, obs_label, jet.pt())

                        for parton_type in parton_types:
                            getattr(self, ('h_%s_JetPt_%s_R%s_%s' % (observable, parton_type, jetR, obs_label)) if \
                                len(obs_label) else ('h_%s_JetPt_%s_R%s' % (observable, parton_type, jetR))).Fill(
                                jet.pt(), obs)

            setattr(self, "count1_R%s" % jetR_str, count1)
            setattr(self, "count2_R%s" % jetR_str, count2)

    #---------------------------------------------------------------
    # Calculate the observable given a jet
    #---------------------------------------------------------------
    def calculate_observable(self, observable, jet, jet_groomed_lund,
        jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed):

        if observable == "ang":

            return fjext.lambda_beta_kappa(jet, jet_groomed_lund.pair(), obs_setting, 1, jetR) \
                   if grooming_setting else fjext.lambda_beta_kappa(jet, obs_setting, 1, jetR)

        elif observable == "mass":

            if grooming_setting:
                j_groomed = jet_groomed_lund.pair()
                if not j_groomed.has_constituents():
                    # Untagged jet -- record underflow value
                    return -1
                else:
                    return j_groomed.m()

            return jet.m()

        # Should not be any other observable
        raise ValueError("Observable %s not implemented" % observable)

    #---------------------------------------------------------------
    # Initiate scaling of all histograms and print final simulation info
    #---------------------------------------------------------------
    def scale_print_final_info(self, pythia):
        # Scale all jet histograms by the appropriate factor from generated cross section and the number of accepted events
        scale_f = pythia.info.sigmaGen() / self.hNevents.GetBinContent(1)

        for jetR in self.jetR_list:
            hist_list_name = "hist_list_R%s" % str(jetR).replace('.', '')
            for h in getattr(self, hist_list_name):
                h.Scale(scale_f)

        print("N total final events:", int(self.hNevents.GetBinContent(1)), "with",
              int(pythia.info.nAccepted() - self.hNevents.GetBinContent(1)),
              "events rejected at hadronization step")
        self.hNevents.SetBinError(1, 0)

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
                        help="Filename for the (unscaled) generated particle ROOT TTree")
    parser.add_argument('--MPIon', action='store', type=int, default=1,
                        help="MPI on or off")
    parser.add_argument('--ISRon', action='store', type=int, default=1,
                        help="ISR on or off")
    parser.add_argument('-c', '--config_file', action='store', type=str, default='config/angularity.yaml',
                        help="Path of config file for observable configurations")

    args = parser.parse_args()

    # If invalid configFile is given, exit
    if not os.path.exists(args.config_file):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    # Use PYTHIA seed for event generation
    if args.user_seed < 0:
        args.user_seed = 1111

    # Have at least 1 event
    if args.nev < 1:
        args.nev = 1

    process = PythiaQuarkGluon(config_file=args.config_file, output_dir=args.output_dir, args=args)
    process.pythia_quark_gluon(args)
