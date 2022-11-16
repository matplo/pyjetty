#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of track information
  and do jet-finding, and save output observables as a numpy array in an hdf5 file.
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse
import time

# Data analysis and plotting
import ROOT
import yaml
import numpy as np
from collections import defaultdict

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.mputils import CEventSubtractor

# Base class
from pyjetty.alice_analysis.process.base import process_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class ProcessDataMultifold(process_base.ProcessBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):

        # Initialize base class
        super(ProcessDataMultifold, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)

        # Initialize configuration
        self.initialize_config()

        # User-specific initialization
        self.initialize_user_config()

        # Create output dict
        # It will have the format: output_dict[R][observable_label] = [obs_jet1, obs_jet2, ...]
        # Every observable must get filled for every jet
        # We will check that the size is the same for all observables
        self.output_dict = {}
        for R in self.jetR_list:
            self.output_dict[R] = {}
            self.output_dict[R]['data'] = defaultdict(list)

    #---------------------------------------------------------------
    # Initialize config file into class members
    #---------------------------------------------------------------
    def initialize_config(self):
    
        # Call base class initialization
        process_base.ProcessBase.initialize_config(self)
    
        # Read config file
        with open(self.config_file, 'r') as stream:
            config = yaml.safe_load(stream)

        self.is_pp = not self.do_constituent_subtraction
    
        # Create dictionaries to store grooming settings and observable settings for each observable
        # Each dictionary entry stores a list of subconfiguration parameters
        #   The observable list stores the observable setting, e.g. subjetR
        #   The grooming list stores a list of grooming settings {'sd': [zcut, beta]} or {'dg': [a]}
        self.observable_list = config['process_observables']
        self.obs_settings = {}
        self.obs_grooming_settings = {}
        for observable in self.observable_list:
            obs_config_dict = config[observable]            
            obs_subconfig_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
            self.obs_settings[observable] = self.utils.obs_settings(observable, obs_config_dict, obs_subconfig_list)
            self.obs_grooming_settings[observable] = self.utils.grooming_settings(obs_config_dict)
            
        # Construct set of unique grooming settings
        self.grooming_settings = []
        lists_grooming = [self.obs_grooming_settings[obs] for obs in self.observable_list]
        for observable in lists_grooming:
            for setting in observable:
                if setting not in self.grooming_settings and setting != None:
                    self.grooming_settings.append(setting)

    #---------------------------------------------------------------
    # Initialize config file into class members
    #---------------------------------------------------------------
    def initialize_user_config(self):
      
        for jetR in self.jetR_list:
            for observable in self.observable_list:

                # Anti-kt subjets -- Define subjet finders
                if observable == 'leading_subjet_z':
                    self.subjet_def = self.recursive_defaultdict()
                    for subjetR in self.obs_settings[observable]:
                        self.subjet_def[jetR][subjetR] = fj.JetDefinition(fj.antikt_algorithm, subjetR)

                if observable == 'jet_axis':
                    self.reclusterer_wta = self.recursive_defaultdict()
                    for jetR in self.jetR_list:
                        jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
                        jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
                        self.reclusterer_wta[jetR] =  fjcontrib.Recluster(jet_def_wta)

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def process_data(self):
    
        self.start_time = time.time()

        # Use IO helper class to convert ROOT TTree into a SeriesGroupBy object of fastjet particles per event
        print('--- {} seconds ---'.format(time.time() - self.start_time))
        io = process_io.ProcessIO(input_file=self.input_file, track_tree_name='tree_Particle',
                                is_pp=self.is_pp, use_ev_id_ext=True)
        self.df_fjparticles = io.load_data(m=self.m)
        self.nEvents = len(self.df_fjparticles.index)
        self.nTracks = len(io.track_df.index)
        print('--- {} seconds ---'.format(time.time() - self.start_time))

        # Create constituent subtractor, if configured
        if not self.is_pp:
            self.constituent_subtractor = [CEventSubtractor(max_distance=R_max, alpha=self.alpha, max_eta=self.max_eta, bge_rho_grid_size=self.bge_rho_grid_size, max_pt_correct=self.max_pt_correct, ghost_area=self.ghost_area, distance_type=fjcontrib.ConstituentSubtractor.deltaR) for R_max in self.max_distance]
        
        print(self)

        # Find jets and fill histograms
        print('Analyze events...')
        self.analyze_events()
    
        # Save output
        n_jets = {len(next(iter(self.output_dict[self.jetR_list[0]]['data'].values())))}
        print(f'Save output arrays (n_jets={n_jets})...')
        self.utils.write_data(self.output_dict, self.output_dir, filename = 'AnalysisResults.h5')

        print('--- {} seconds ---'.format(time.time() - self.start_time))

    #---------------------------------------------------------------
    # Main function to loop through and analyze events
    #---------------------------------------------------------------
    def analyze_events(self):

        print('Find jets...')
        fj.ClusterSequence.print_banner()
        print()
        self.event_number = 0
    
        # Use list comprehension to do jet-finding and fill histograms
        [self.analyze_event(fj_particles) for fj_particles in self.df_fjparticles]
        
        print('--- {} seconds ---'.format(time.time() - self.start_time))
    
    #---------------------------------------------------------------
    # Analyze jets of a given event.
    # fj_particles is the list of fastjet pseudojets for a single fixed event.
    #---------------------------------------------------------------
    def analyze_event(self, fj_particles):
  
        self.event_number += 1
        if self.event_number > self.event_number_max:
            return
        if self.debug_level > 1:
            print('-------------------------------------------------')
            print('event {}'.format(self.event_number))
    
        if len(fj_particles) > 1:
            if np.abs(fj_particles[0].pt() - fj_particles[1].pt()) <  1e-10:
                print('WARNING: Duplicate particles may be present')
                print([p.user_index() for p in fj_particles])
                print([p.pt() for p in fj_particles])
        
        # Perform constituent subtraction for each R_max (do this once, for all jetR)
        if not self.is_pp:
            fj_particles_subtracted = [self.constituent_subtractor[i].process_event(fj_particles) for i, R_max in enumerate(self.max_distance)]
    
        # Loop through jetR, and process event for each R
        for jetR in self.jetR_list:

            # Set jet definition and a jet selector
            jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
            jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
            if self.debug_level > 2:
                print('jet definition is:', jet_def)
                print('jet selector is:', jet_selector,'\n')
        
            # Analyze
            if self.is_pp:
            
                # Do jet finding
                cs = fj.ClusterSequence(fj_particles, jet_def)
                jets = fj.sorted_by_pt(cs.inclusive_jets())
                jets_selected = jet_selector(jets)
            
                self.analyze_jets(jets_selected, jetR)
        
            else:
            
                for i, R_max in enumerate(self.max_distance):
                            
                    if self.debug_level > 1:
                        print('R_max: {}'.format(R_max))
                        
                    # Keep track of whether to fill R_max-independent histograms
                    self.fill_Rmax_indep_hists = (i == 0)
                    
                    # Perform constituent subtraction
                    rho = self.constituent_subtractor[i].bge_rho.rho()
                    if self.fill_R_indep_hists and self.fill_Rmax_indep_hists:
                        getattr(self, 'hRho').Fill(rho)
                    
                    # Do jet finding (re-do each time, to make sure matching info gets reset)
                    cs = fj.ClusterSequence(fj_particles_subtracted[i], jet_def)
                    jets = fj.sorted_by_pt(cs.inclusive_jets())
                    jets_selected = jet_selector(jets)
                    
                    self.analyze_jets(jets_selected, jetR, R_max = R_max)
                
    #---------------------------------------------------------------
    # Analyze jets of a given event.
    #---------------------------------------------------------------
    def analyze_jets(self, jets_selected, jetR, R_max = None):

        # Set suffix for filling histograms
        if R_max:
            suffix = '_Rmax{}'.format(R_max)
        else:
            suffix = ''
        
        # Loop through jets and call user function to fill histos
        [self.analyze_accepted_jet(jet, jetR, suffix) for jet in jets_selected]
  
    #---------------------------------------------------------------
    # Fill histograms
    #---------------------------------------------------------------
    def analyze_accepted_jet(self, jet, jetR, suffix):

        # Check additional acceptance criteria
        if not self.utils.is_det_jet_accepted(jet):
            return
          
        # Fill base histograms
        jet_pt_ungroomed = jet.pt()

        # Append the jet pt, eta
        self.output_dict[jetR]['data'][f'jet_pt{suffix}'].append(jet_pt_ungroomed)
        self.output_dict[jetR]['data'][f'jet_eta{suffix}'].append(jet.eta())

        # Loop through each jet subconfiguration (i.e. subobservable / grooming setting)
        # Note that the subconfigurations are defined by the first observable, if multiple are defined
        for observable in self.observable_list:
            for i in range(len(self.obs_settings[observable])):
        
                obs_setting = self.obs_settings[observable][i]
                grooming_setting = self.obs_grooming_settings[observable][i]
        
                # Groom jet, if applicable
                if grooming_setting:
                    gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
                    jet_groomed_lund = self.utils.groom(gshop, grooming_setting, jetR)
                    if not jet_groomed_lund:
                        continue
                else:
                    jet_groomed_lund = None

                # Call user function to fill histograms
                self.fill_jet_observables(jet, jet_groomed_lund, jetR, observable, obs_setting, 
                                          grooming_setting, jet_pt_ungroomed, suffix)

        # Check that the output size is correct, i.e. that we have filled each observable for every jet
        expected_size = len(next(iter(self.output_dict[jetR]['data'].values())))
        for jetR in self.jetR_list:
            for observable_label,result in self.output_dict[jetR]['data'].items():
                if len(result) != expected_size:
                    raise ValueError(f'Observable {observable_label} does not have expected length of {expected_size}! {self.output_dict}')

    #---------------------------------------------------------------
    # This function is called once for each jet subconfiguration
    #---------------------------------------------------------------
    def fill_jet_observables(self, jet, jet_groomed_lund, jetR, observable, obs_setting, 
                             grooming_setting, jet_pt_ungroomed, suffix):
    
        observable_label = f'{observable}_{self.utils.obs_label(obs_setting, grooming_setting)}{suffix}'
        if self.debug_level > 1:
            print(f'filling {observable_label}...')

        if observable == 'theta_g':
            result = jet_groomed_lund.Delta() / jetR
                        
        elif observable == 'zg':
            result = jet_groomed_lund.z()

        elif observable == 'leading_subjet_z':
            cs_subjet = fj.ClusterSequence(jet.constituents(), self.subjet_def[jetR][obs_setting])
            subjets = fj.sorted_by_pt(cs_subjet.inclusive_jets())
            leading_subjet = self.utils.leading_jet(subjets)
            result = leading_subjet.pt() / jet_pt_ungroomed
            if np.isclose(result, 1.):
                result = 0.999
        
        elif observable == 'jet_axis':
            if obs_setting == 'Standard_SD':
                jet_groomed = jet_groomed_lund.pair()
                result = jet.delta_R(jet_groomed)
            if obs_setting == 'Standard_WTA':
                jet_wta = self.reclusterer_wta[jetR].result(jet)
                result = jet.delta_R(jet_wta)
            if obs_setting == 'WTA_SD':
                jet_wta = self.reclusterer_wta[jetR].result(jet)
                jet_groomed = jet_groomed_lund.pair()
                result = jet_groomed.delta_R(jet_wta)

        elif observable == 'ang':
            kappa = 1
            if grooming_setting:
                jet_groomed = jet_groomed_lund.pair()
                result = fjext.lambda_beta_kappa(jet, jet_groomed, obs_setting, kappa, jetR)
            else:
                result = fjext.lambda_beta_kappa(jet, obs_setting, kappa, jetR)

        else:
            sys.exit(f'Observable {observable} not found!')

        # Append to the output
        self.output_dict[jetR]['data'][observable_label].append(result)

##################################################################
if __name__ == '__main__':
    # Define arguments
    parser = argparse.ArgumentParser(description='Process data')
    parser.add_argument('-f', '--inputFile', action='store',
                        type=str, metavar='inputFile',
                        default='AnalysisResults.root',
                        help='Path of ROOT file containing TTrees')
    parser.add_argument('-c', '--configFile', action='store',
                        type=str, metavar='configFile',
                        default='config/analysis_config.yaml',
                        help="Path of config file for analysis")
    parser.add_argument('-o', '--outputDir', action='store',
                        type=str, metavar='outputDir',
                        default='./TestOutput',
                        help='Output directory for output to be written to')
    
    # Parse the arguments
    args = parser.parse_args()
    
    print('Configuring...')
    print('inputFile: \'{0}\''.format(args.inputFile))
    print('configFile: \'{0}\''.format(args.configFile))
    print('ouputDir: \'{0}\"'.format(args.outputDir))
    print('----------------------------------------------------------------')
  
    # If invalid inputFile is given, exit
    if not os.path.exists(args.inputFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
        sys.exit(0)
    
    # If invalid configFile is given, exit
    if not os.path.exists(args.configFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    analysis = ProcessDataMultifold(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
    analysis.process_data()