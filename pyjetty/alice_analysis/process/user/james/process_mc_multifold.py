#!/usr/bin/env python3

"""
  Analysis class to read a ROOT TTree of MC track information
  and do jet-finding, and save output observables as a numpy array in an hdf5 file.
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import time
import os
import sys
import argparse

# Data analysis and plotting
import pandas
import numpy as np
import ROOT
import yaml
from array import *
from collections import defaultdict

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_io_emb
from pyjetty.alice_analysis.process.base import thermal_generator
from pyjetty.mputils import CEventSubtractor

# Base class
from pyjetty.alice_analysis.process.base import process_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class ProcessMCMultifold(process_base.ProcessBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):

        # Initialize base class
        super(ProcessMCMultifold, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
            
        # Initialize configuration
        self.initialize_config()

        # User-specific initialization
        self.initialize_user_config()

        # Create output dict
        # It will have the format: output_dict[R][reconstruction_level][observable_label] = [obs_jet1, obs_jet2, ...]
        # Every observable must get filled for every jet
        # We will check that the size is the same for all observables
        self.output_dict = {}
        for R in self.jetR_list:
            self.output_dict[R] = {}
            self.output_dict[R]['det_matched'] = defaultdict(list)
            self.output_dict[R]['truth_matched'] = defaultdict(list)
            self.output_dict[R]['det_before_matching'] = defaultdict(list)
            self.output_dict[R]['truth_before_matching'] = defaultdict(list)

    #---------------------------------------------------------------
    # Initialize config file into class members
    #---------------------------------------------------------------
    def initialize_config(self):

        # Call base class initialization
        process_base.ProcessBase.initialize_config(self)

        # Read config file
        with open(self.config_file, 'r') as stream:
            config = yaml.safe_load(stream)

        self.fast_simulation = config['fast_simulation']
        self.jet_matching_distance = config['jet_matching_distance']
        self.reject_tracks_fraction = config['reject_tracks_fraction']
        if 'mc_fraction_threshold' in config:
            self.mc_fraction_threshold = config['mc_fraction_threshold']
        
        if self.do_constituent_subtraction:
            self.is_pp = False
            self.emb_file_list = config['emb_file_list']
            self.main_R_max = config['constituent_subtractor']['main_R_max']
        else:
            self.is_pp = True
            
        if 'thermal_model' in config:
            self.thermal_model = True
            beta = config['thermal_model']['beta']
            N_avg = config['thermal_model']['N_avg']
            sigma_N = config['thermal_model']['sigma_N']
            self.thermal_generator = thermal_generator.ThermalGenerator(N_avg, sigma_N, beta)
        else:
            self.thermal_model = False
         
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
                        self.reclusterer_wta[jetR] = fjcontrib.Recluster(jet_def_wta)

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def process_mc(self):
    
        self.start_time = time.time()
    
        # ------------------------------------------------------------------------
    
        # Use IO helper class to convert detector-level ROOT TTree into
        # a SeriesGroupBy object of fastjet particles per event
        print(f'--- {time.time()-self.start_time} seconds ---')
        if self.fast_simulation:
            tree_dir = ''
        else:
            tree_dir = 'PWGHF_TreeCreator'
        io_det = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
                                      track_tree_name='tree_Particle', use_ev_id_ext=False)
        df_fjparticles_det = io_det.load_data(m=self.m, reject_tracks_fraction=self.reject_tracks_fraction)
        self.nEvents_det = len(df_fjparticles_det.index)
        self.nTracks_det = len(io_det.track_df.index)
        print(f'--- {time.time()-self.start_time} seconds ---')

        # ------------------------------------------------------------------------

        # Use IO helper class to convert truth-level ROOT TTree into
        # a SeriesGroupBy object of fastjet particles per event
        io_truth = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
                                        track_tree_name='tree_Particle_gen', use_ev_id_ext=False)
        df_fjparticles_truth = io_truth.load_data(m=self.m)
        self.nEvents_truth = len(df_fjparticles_truth.index)
        self.nTracks_truth = len(io_truth.track_df.index)
        print(f'--- {time.time()-self.start_time} seconds ---')

        # ------------------------------------------------------------------------

        # Now merge the two SeriesGroupBy to create a groupby df with [ev_id, run_number, fj_1, fj_2]
        # (Need a structure such that we can iterate event-by-event through both fj_1, fj_2 simultaneously)
        print('Merge det-level and truth-level into a single dataframe grouped by event...')
        self.df_fjparticles = pandas.concat([df_fjparticles_det, df_fjparticles_truth], axis=1)
        self.df_fjparticles.columns = ['fj_particles_det', 'fj_particles_truth']
        print(f'--- {time.time()-self.start_time} seconds ---')

        # ------------------------------------------------------------------------
    
        # Set up the Pb-Pb embedding object
        if not self.is_pp and not self.thermal_model:
            self.process_io_emb = process_io_emb.ProcessIO_Emb(self.emb_file_list, track_tree_name='tree_Particle', m=self.m)
        
        # ------------------------------------------------------------------------

        # Create constituent subtractor, if configured
        if self.do_constituent_subtraction:
            self.constituent_subtractor = [CEventSubtractor(max_distance=R_max, alpha=self.alpha, max_eta=self.max_eta, bge_rho_grid_size=self.bge_rho_grid_size, max_pt_correct=self.max_pt_correct, ghost_area=self.ghost_area, distance_type=fjcontrib.ConstituentSubtractor.deltaR) for R_max in self.max_distance]
        
        print(self)
        
        # Find jets and fill histograms
        print('Analyze events...')
        self.analyze_events()
        
        # Save output
        n_jets = len(next(iter(self.output_dict[self.jetR_list[0]]['det_matched'].values())))
        print(f'Save output arrays (n_jets_matched={n_jets} from n_events={self.event_number})...')
        self.output_dict['n_events'] = self.event_number
        self.utils.write_data(self.output_dict, self.output_dir, filename = 'AnalysisResults.h5')
        
        print(f'--- {time.time()-self.start_time} seconds ---')

    #---------------------------------------------------------------
    # Main function to loop through and analyze events
    #---------------------------------------------------------------
    def analyze_events(self):
    
        fj.ClusterSequence.print_banner()
        print()
            
        self.event_number = 0
        
        # Then can use list comprehension to iterate over the groupby and do jet-finding
        # simultaneously for fj_1 and fj_2 per event, so that I can match jets -- and fill histograms
        [self.analyze_event(fj_particles_det, fj_particles_truth) for fj_particles_det, fj_particles_truth in zip(self.df_fjparticles['fj_particles_det'], self.df_fjparticles['fj_particles_truth'])]
      
    #---------------------------------------------------------------
    # Analyze jets of a given event.
    # fj_particles is the list of fastjet pseudojets for a single fixed event.
    #---------------------------------------------------------------
    def analyze_event(self, fj_particles_det, fj_particles_truth):
  
        self.event_number += 1
        if self.event_number > self.event_number_max:
            return
        if self.debug_level > 1:
            print('-------------------------------------------------')
            print(f'event {self.event_number}')
      
        # Check that the entries exist appropriately
        # (need to check how this can happen -- but it is only a tiny fraction of events)
        if type(fj_particles_det) != fj.vectorPJ or type(fj_particles_truth) != fj.vectorPJ:
            print('fj_particles type mismatch -- skipping event')
            return
      
        if len(fj_particles_truth) > 1:
            if np.abs(fj_particles_truth[0].pt() - fj_particles_truth[1].pt()) <  1e-10:
                print('WARNING: Duplicate particles may be present')
                print([p.user_index() for p in fj_particles_truth])
                print([p.pt() for p in fj_particles_truth])
                
        # If Pb-Pb, construct embedded event (do this once, for all jetR)
        if not self.is_pp:
        
            # If thermal model, generate a thermal event and add it to the det-level particle list
            if self.thermal_model:
                fj_particles_combined_beforeCS = self.thermal_generator.load_event()
                
                # Form the combined det-level event
                # The pp-det tracks are each stored with a unique user_index >= 0
                #   (same index in fj_particles_combined and fj_particles_det -- which will be used in prong-matching)
                # The thermal tracks are each stored with a unique user_index < 0
                [fj_particles_combined_beforeCS.push_back(p) for p in fj_particles_det]

            # Main case: Get Pb-Pb event and embed it into the det-level particle list
            else:
                fj_particles_combined_beforeCS = self.process_io_emb.load_event()
              
                # Form the combined det-level event
                # The pp-det tracks are each stored with a unique user_index >= 0
                #   (same index in fj_particles_combined and fj_particles_det -- which will be used in prong-matching)
                # The Pb-Pb tracks are each stored with a unique user_index < 0
                [fj_particles_combined_beforeCS.push_back(p) for p in fj_particles_det]
         
            # Perform constituent subtraction for each R_max
            fj_particles_combined = [self.constituent_subtractor[i].process_event(fj_particles_combined_beforeCS) for i, R_max in enumerate(self.max_distance)]
        
            if self.debug_level > 3:
                print([p.user_index() for p in fj_particles_truth])
                print([p.pt() for p in fj_particles_truth])
                print([p.user_index() for p in fj_particles_det])
                print([p.pt() for p in fj_particles_det])
                print([p.user_index() for p in fj_particles_combined_beforeCS])
                print([p.pt() for p in fj_particles_combined_beforeCS])

        # Loop through jetR, and process event for each R
        for jetR in self.jetR_list:
    
            # Set jet definition and a jet selector
            jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
            jet_selector_det = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
            jet_selector_truth_matched = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9)
            if self.debug_level > 2:
                print('')
                print('jet definition is:', jet_def)
                print('jet selector for det-level is:', jet_selector_det)
                print('jet selector for truth-level matches is:', jet_selector_truth_matched)
      
            # Analyze
            if self.is_pp:
            
                # Find pp det and truth jets
                cs_det = fj.ClusterSequence(fj_particles_det, jet_def)
                jets_det_pp = fj.sorted_by_pt(cs_det.inclusive_jets())
                jets_det_pp_selected = jet_selector_det(jets_det_pp)
                
                cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)
                jets_truth = fj.sorted_by_pt(cs_truth.inclusive_jets())
                jets_truth_selected = jet_selector_det(jets_truth)
                jets_truth_selected_matched = jet_selector_truth_matched(jets_truth)
            
                self.analyze_jets(jets_det_pp_selected, jets_truth_selected, jets_truth_selected_matched, jetR)
        
            else:
      
                for i, R_max in enumerate(self.max_distance):
            
                    if self.debug_level > 1:
                        print('')
                        print(f'R_max: {R_max}')
                        print(f'Total number of combined particles: {len([p.pt() for p in fj_particles_combined_beforeCS])}')
                        print(f'After constituent subtraction {i}: {len([p.pt() for p in fj_particles_combined[i]])}')
                        
                    # Keep track of whether to fill R_max-independent histograms
                    self.fill_Rmax_indep_hists = (i == 0)

                    # Do jet finding (re-do each time, to make sure matching info gets reset)
                    cs_det = fj.ClusterSequence(fj_particles_det, jet_def)
                    jets_det_pp = fj.sorted_by_pt(cs_det.inclusive_jets())
                    jets_det_pp_selected = jet_selector_det(jets_det_pp)
                    
                    cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)
                    jets_truth = fj.sorted_by_pt(cs_truth.inclusive_jets())
                    jets_truth_selected = jet_selector_det(jets_truth)
                    jets_truth_selected_matched = jet_selector_truth_matched(jets_truth)
                    
                    cs_combined = fj.ClusterSequence(fj_particles_combined[i], jet_def)
                    jets_combined = fj.sorted_by_pt(cs_combined.inclusive_jets())
                    jets_combined_selected = jet_selector_det(jets_combined)

                    self.analyze_jets(jets_combined_selected, jets_truth_selected, jets_truth_selected_matched, jetR,
                                        jets_det_pp_selected = jets_det_pp_selected, R_max = R_max)

    #---------------------------------------------------------------
    # Analyze jets of a given event.
    #---------------------------------------------------------------
    def analyze_jets(self, jets_det_selected, jets_truth_selected, jets_truth_selected_matched, jetR,
                     jets_det_pp_selected = None, R_max = None):

        if self.debug_level > 1:
            print(f'Number of det-level jets: {len(jets_det_selected)}')

        # Set suffix for filling histograms
        if R_max:
            suffix = f'_Rmax{R_max}'
        else:
            suffix = ''
    
        # Fill det-level jet histograms (before matching)
        for jet_det in jets_det_selected:
        
            # Check additional acceptance criteria
            # skip event if not satisfied -- since first jet in event is highest pt
            if not self.utils.is_det_jet_accepted(jet_det):
                if self.debug_level > 1:
                    print('event rejected due to jet acceptance')
                return
        
            # Fill det-level observables (before matching)
            if self.thermal_model:
                self.fill_jet_observables_before_matching(jet_det, jetR, suffix, key='det_before_matching')
  
        # Fill truth-level jet observables (before matching)
        for jet_truth in jets_truth_selected:
        
            if self.is_pp or self.fill_Rmax_indep_hists:
                self.fill_jet_observables_before_matching(jet_truth, jetR, suffix, key='truth_before_matching')
        
        if self.debug_level > 0:
            print(f'Finished filling observables before matching')

        # Loop through jets and set jet matching candidates for each jet in user_info
        if self.is_pp:
            [[self.set_matching_candidates(jet_det, jet_truth, jetR, None) for jet_truth in jets_truth_selected_matched] for jet_det in jets_det_selected]
        else:
            # First fill the combined-to-pp matches, then the pp-to-pp matches
            [[self.set_matching_candidates(jet_det_combined, jet_det_pp, jetR, None, fill_jet1_matches_only=True) for jet_det_pp in jets_det_pp_selected] for jet_det_combined in jets_det_selected]
            [[self.set_matching_candidates(jet_det_pp, jet_truth, jetR, None) for jet_truth in jets_truth_selected_matched] for jet_det_pp in jets_det_pp_selected]
        
        # Loop through jets and set accepted matches
        hname = f'hJetMatchingQA_R{jetR}{suffix}'
        if self.is_pp:
            [self.set_matches_pp(jet_det, hname) for jet_det in jets_det_selected]
        else:
            [self.set_matches_AA(jet_det_combined, jetR, hname) for jet_det_combined in jets_det_selected]
          
        # Loop through jets and fill response histograms if both det and truth jets are unique match
        [self.fill_matched_jet_observables(jet_det, jetR, suffix) for jet_det in jets_det_selected]

        if self.debug_level > 0:
            print(f'Finished filling matched jet observables')
  
    #---------------------------------------------------------------
    # This function is called once for each jet
    #---------------------------------------------------------------
    def fill_jet_observables_before_matching(self, jet, jetR, suffix, key=''):

        self.output_dict[jetR]['truth_before_matching'][f'jet_pt{suffix}'].append(jet.pt())

        # Loop through each jet subconfiguration (i.e. subobservable / grooming setting)
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
                self.fill_observable(jet, jet_groomed_lund, jetR, observable, obs_setting,
                                     grooming_setting, jet.pt(), suffix, key)

        # Check that the output size is correct, i.e. that we have filled each observable for every jet
        expected_size = len(next(iter(self.output_dict[jetR][key].values())))
        for jetR in self.jetR_list:
            for observable_label,result in self.output_dict[jetR][key].items():
                if len(result) != expected_size:
                    raise ValueError(f'Observable {observable_label} {key} does not have expected length of {expected_size}! {self.output_dict[jetR][key]}')

    #---------------------------------------------------------------
    # This function is called once for each jet subconfiguration
    #---------------------------------------------------------------
    def fill_observable(self, jet, jet_groomed_lund, jetR, observable, obs_setting,
                        grooming_setting, jet_pt_ungroomed, suffix, key):

        observable_label = f'{observable}_{self.utils.obs_label(obs_setting, grooming_setting)}{suffix}'
        if self.debug_level > 1:
            print(f'filling {key} {observable_label}...')

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

        if self.debug_level > 1:
            print(result)

        # Append to the output
        self.output_dict[jetR][key][observable_label].append(result)

    #---------------------------------------------------------------
    # Loop through jets and call user function to fill matched
    # histos if both det and truth jets are unique match.
    #---------------------------------------------------------------
    def fill_matched_jet_observables(self, jet_det, jetR, suffix):
  
        # Get matched truth jet
        if jet_det.has_user_info():
            jet_truth = jet_det.python_info().match
      
            if jet_truth:
        
                jet_pt_det_ungroomed = jet_det.pt()
                jet_pt_truth_ungroomed = jet_truth.pt()
                self.output_dict[jetR]['det_matched'][f'jet_pt{suffix}'].append(jet_pt_det_ungroomed)
                self.output_dict[jetR]['truth_matched'][f'jet_pt{suffix}'].append(jet_pt_truth_ungroomed)
        
                # If Pb-Pb case, we need to keep jet_det, jet_truth, jet_pp_det
                jet_pp_det = None
                if not self.is_pp:
        
                    # Get pp-det jet
                    jet_pp_det = jet_truth.python_info().match
                
                # Loop through each jet subconfiguration (i.e. subobservable / grooming setting)
                for observable in self.observable_list:
                    for i in range(len(self.obs_settings[observable])):
            
                        obs_setting = self.obs_settings[observable][i]
                        grooming_setting = self.obs_grooming_settings[observable][i]
                        obs_label = self.utils.obs_label(obs_setting, grooming_setting)
                        
                        if self.debug_level > 3:
                            print(f'obs_label: {obs_label}')
          
                        # Groom jets, if applicable
                        if grooming_setting:
                                    
                            # Groom det jet
                            gshop_det = fjcontrib.GroomerShop(jet_det, jetR, self.reclustering_algorithm)
                            jet_det_groomed_lund = self.utils.groom(gshop_det, grooming_setting, jetR)
                            if not jet_det_groomed_lund:
                                continue

                            # Groom truth jet
                            gshop_truth = fjcontrib.GroomerShop(jet_truth, jetR, self.reclustering_algorithm)
                            jet_truth_groomed_lund = self.utils.groom(gshop_truth, grooming_setting, jetR)
                            if not jet_truth_groomed_lund:
                                continue
                        
                        else:
          
                            jet_det_groomed_lund = None
                            jet_truth_groomed_lund = None
              
                        # Call user function to fill histos
                        self.fill_matched_jet_observable(jet_det, jet_det_groomed_lund, jet_truth,
                                                         jet_truth_groomed_lund, jet_pp_det, jetR,
                                                         observable, obs_setting, grooming_setting,
                                                         jet_pt_det_ungroomed, jet_pt_truth_ungroomed,
                                                         suffix)

    #---------------------------------------------------------------
    # Fill matched jet histograms
    #---------------------------------------------------------------
    def fill_matched_jet_observable(self, jet_det, jet_det_groomed_lund, jet_truth,
                                    jet_truth_groomed_lund, jet_pp_det, jetR,
                                    observable, obs_setting, grooming_setting,
                                    jet_pt_det_ungroomed, jet_pt_truth_ungroomed, suffix):
        
        observable_label = f'{observable}_{self.utils.obs_label(obs_setting, grooming_setting)}{suffix}'
        if self.debug_level > 1:
            print(f'filling {observable_label}...')

        if observable == 'theta_g':
            result_det = jet_det_groomed_lund.Delta() / jetR
            result_truth = jet_truth_groomed_lund.Delta() / jetR

        elif observable == 'zg':
            result_det = jet_det_groomed_lund.z()
            result_truth = jet_truth_groomed_lund.z()

        elif observable == 'leading_subjet_z':

            # Find all subjets
            subjetR = obs_setting
            cs_subjet_det = fj.ClusterSequence(jet_det.constituents(), self.subjet_def[jetR][subjetR])
            subjets_det = fj.sorted_by_pt(cs_subjet_det.inclusive_jets())

            cs_subjet_truth = fj.ClusterSequence(jet_truth.constituents(), self.subjet_def[jetR][subjetR])
            subjets_truth = fj.sorted_by_pt(cs_subjet_truth.inclusive_jets())
    
            if not self.is_pp:
                cs_subjet_det_pp = fj.ClusterSequence(jet_pp_det.constituents(), self.subjet_def[jetR][subjetR])
                subjets_det_pp = fj.sorted_by_pt(cs_subjet_det_pp.inclusive_jets())

            # Loop through subjets and set subjet matching candidates for each subjet in user_info
            if self.is_pp:
                [[self.set_matching_candidates(subjet_det, subjet_truth, subjetR, None) for subjet_truth in subjets_truth] for subjet_det in subjets_det]
            else:
                # First fill the combined-to-pp matches, then the pp-to-pp matches
                [[self.set_matching_candidates(subjet_det_combined, subjet_det_pp, subjetR, None, fill_jet1_matches_only=True) for subjet_det_pp in subjets_det_pp] for subjet_det_combined in subjets_det]
                [[self.set_matching_candidates(subjet_det_pp, subjet_truth, subjetR, None) for subjet_truth in subjets_truth] for subjet_det_pp in subjets_det_pp]
      
            # Loop through subjets and set accepted matches
            if self.is_pp:
                [self.set_matches_pp(subjet_det, 'hSubjetMatchingQA_R{}_{}'.format(jetR, subjetR)) for subjet_det in subjets_det]
            else:
                [self.set_matches_AA(subjet_det_combined, subjetR, 'hSubjetMatchingQA_R{}_{}'.format(jetR, subjetR)) for subjet_det_combined in subjets_det]

            # Get leading subjet and fill histograms      
            leading_subjet_det = self.utils.leading_jet(subjets_det)
            leading_subjet_truth = self.utils.leading_jet(subjets_truth)
        
            # Note that we don't want to check whether they are geometrical matches
            # We rather want to correct the measured leading subjet to the true leading subjet
            if leading_subjet_det and leading_subjet_truth:
          
                z_leading_det = leading_subjet_det.pt() / jet_pt_det_ungroomed
                z_leading_truth = leading_subjet_truth.pt() / jet_pt_truth_ungroomed
                
                # If z=1, it will be default be placed in overflow bin -- prevent this
                if np.isclose(z_leading_det, 1.):
                    z_leading_det = 0.999
                if np.isclose(z_leading_truth, 1.):
                    z_leading_truth = 0.999
          
                result_det = z_leading_det
                result_truth = z_leading_truth
            else:
                raise ValueError('Leading subjet not found...')

        elif observable == 'jet_axis':
            if obs_setting == 'Standard_SD':

                jet_groomed_det = jet_det_groomed_lund.pair()
                result_det = jet_det.delta_R(jet_groomed_det)

                jet_groomed_truth = jet_truth_groomed_lund.pair()
                result_truth = jet_truth.delta_R(jet_groomed_truth)

            if obs_setting == 'Standard_WTA':

                jet_wta_det = self.reclusterer_wta[jetR].result(jet_det)
                result_det = jet_det.delta_R(jet_wta_det)

                jet_wta_truth = self.reclusterer_wta[jetR].result(jet_truth)
                result_truth = jet_truth.delta_R(jet_wta_truth)

            if obs_setting == 'WTA_SD':

                jet_wta_det = self.reclusterer_wta[jetR].result(jet_det)
                jet_groomed_det = jet_det_groomed_lund.pair()
                result_det = jet_groomed_det.delta_R(jet_wta_det)

                jet_wta_truth = self.reclusterer_wta[jetR].result(jet_truth)
                jet_groomed_truth = jet_truth_groomed_lund.pair()
                result_truth = jet_groomed_truth.delta_R(jet_wta_truth)

        elif observable == 'ang':

            kappa = 1
            if grooming_setting:
                jet_groomed_det = jet_det_groomed_lund.pair()
                result_det = fjext.lambda_beta_kappa(jet_det, jet_groomed_det, obs_setting, kappa, jetR)
                jet_groomed_truth = jet_truth_groomed_lund.pair()
                result_truth = fjext.lambda_beta_kappa(jet_truth, jet_groomed_truth, obs_setting, kappa, jetR)
            else:
                result_det = fjext.lambda_beta_kappa(jet_det, obs_setting, kappa, jetR)
                result_truth = fjext.lambda_beta_kappa(jet_truth, obs_setting, kappa, jetR)

        else:
            sys.exit(f'Observable {observable} not found!')

        # Append to the output
        self.output_dict[jetR]['det_matched'][observable_label].append(result_det)
        self.output_dict[jetR]['truth_matched'][observable_label].append(result_truth)

##################################################################
if __name__ == '__main__':
    # Define arguments
    parser = argparse.ArgumentParser(description='Process MC')
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

    # If invalid inputFile is given, exit
    if not os.path.exists(args.inputFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
        sys.exit(0)
  
    # If invalid configFile is given, exit
    if not os.path.exists(args.configFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    analysis = ProcessMCMultifold(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
    analysis.process_mc()
