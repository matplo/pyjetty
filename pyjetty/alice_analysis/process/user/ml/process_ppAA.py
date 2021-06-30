#!/usr/bin/env python3

"""
Class to read pp vs AA data set, do jet finding, and compute N-subjettiness
"""

import os
import sys
import argparse
import yaml
import h5py
import pickle

# Data analysis and plotting
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext

# Energy flow package
import energyflow

# Analysis utilities
from pyjetty.mputils import CEventSubtractor
from pyjetty.alice_analysis.process.base import thermal_generator
from pyjetty.alice_analysis.process.base import process_base
from pyjetty.alice_analysis.process.base import process_utils

# Base class
from pyjetty.alice_analysis.process.base import common_base

################################################################
class ProcessppAA(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, config_file='', input_file='', output_dir='', **kwargs):
        super(common_base.CommonBase, self).__init__(**kwargs)
        
        self.config_file = config_file
        self.input_file = input_file
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        
        # Initialize config file
        self.initialize_config()
        
        # Initialize utils class
        self.utils = process_utils.ProcessUtils()
        
        # Load dataframe of particle four-vectors for all particles in the event
        # (separate dataframes for hard process and background)
        print()
        print('Loading particle dataframes')
        with open(self.input_file, 'rb') as f:
            df_particles_hard = pickle.load(f)
            if not self.thermal_model:
                df_particles_background = pickle.load(f)
        print('Done.')
        print()
        
        # Apply eta cut
        #print(f'df_particles_hard (before filtering pseudorapidity): {df_particles_hard.describe()}')
        df_particles_hard = self.apply_eta_cut(df_particles_hard)
        print(f'df_particles_hard (after filtering pseudorapidity): {df_particles_hard.describe()}')
        
        if self.thermal_model:
            # Construct dummy combined event -- we will add thermal particles in event loop
            df_particles_combined = df_particles_hard
        else:
            #print(f'df_particles_background (before filtering pseudorapidity): {df_particles_background.describe()}')
            df_particles_background = self.apply_eta_cut(df_particles_background)
            print(f'df_particles_background (after filtering pseudorapidity): {df_particles_background.describe()}')
        
            # Construct a dataframe of the combined hard+background particles
            df_particles_combined = pd.concat([df_particles_hard, df_particles_background])

        # Next, we will transform these into fastjet::PseudoJet objects.
        # This allows us to do jet finding and use the fastjet contrib to compute Nsubjettiness
        
        # (i) Group the particle dataframe by event id
        #     df_particles_grouped is a DataFrameGroupBy object with one particle dataframe per event
        df_particles_hard_grouped = df_particles_hard.groupby('event_id')
        df_particles_combined_grouped = df_particles_combined.groupby('event_id')
        
        # (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet::PseudoJets
        print('Converting particle dataframes to fastjet::PseudoJets...')
        df_fjparticles_hard = df_particles_hard_grouped.apply(self.get_fjparticles)
        df_fjparticles_combined = df_particles_combined_grouped.apply(self.get_fjparticles)
        self.df_fjparticles = pd.concat([df_fjparticles_hard, df_fjparticles_combined], axis=1)
        self.df_fjparticles.columns = ['fj_particles_hard', 'fj_particles_combined']
        print('Done.')
        print()
        
        # Create list of N-subjettiness observables: number of axes and beta values
        self.N_list = []
        self.beta_list = []
        for i in range(self.K-2):
            self.N_list += [i+1] * 3
            self.beta_list += [0.5,1,2]
        self.N_list += [self.K-1] * 2  
        self.beta_list += [1,2]
        
        # Construct dictionary to store all jet quantities of interest
        self.jet_variables = {'hard': {}, 'combined': {}}
        self.four_vectors = {'hard': {}, 'combined': {}}
        self.jet_qa_variables = {'hard': {}, 'combined': {}}
        self.delta_pt_random_cone = []
        for label in self.jet_variables.keys():
            for jetR in self.jetR_list:
                self.jet_variables[label][f'R{jetR}'] = {}
                self.four_vectors[label][f'R{jetR}'] = {}
                self.jet_qa_variables[label][f'R{jetR}'] = {}
                for jet_pt_bin in self.jet_pt_bins:
                    self.jet_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'] = {}
                    self.four_vectors[label][f'R{jetR}'][f'pt{jet_pt_bin}'] = {}
                    self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'] = {}
                    for R_max in self.max_distance:
                        
                        if label=='combined' or np.isclose(R_max, 0.):
                            self.jet_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'] = {}
                            self.four_vectors[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'] = {}
                            self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'] = {}
                            for i,N in enumerate(self.N_list):
                                beta = self.beta_list[i]
                                self.jet_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'][f'n_subjettiness_N{N}_beta{beta}'] = []
                            self.four_vectors[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['jet_constituent_four_vectors'] = []
                            
                            # Some QA
                            self.qa_observables = ['delta_pt', 'jet_pt', 'jet_angularity', 'jet_mass', 'jet_theta_g', 'jet_subjet_z', 'hadron_z']
                            for qa_observable in self.qa_observables:
                                self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'][qa_observable] = []

        # Create constituent subtractors
        self.constituent_subtractor = [CEventSubtractor(max_distance=R_max, alpha=self.alpha, max_eta=self.eta_max, bge_rho_grid_size=self.bge_rho_grid_size, max_pt_correct=self.max_pt_correct, ghost_area=self.ghost_area, distance_type=fjcontrib.ConstituentSubtractor.deltaR) for R_max in self.max_distance]
                
        print(self)
        print()

    #---------------------------------------------------------------
    # Initialize config file into class members
    #---------------------------------------------------------------
    def initialize_config(self):
    
        # Read config file
        with open(self.config_file, 'r') as stream:
          config = yaml.safe_load(stream)
          
        self.K = max(config['K'])
        
        self.jetR_list = config['jetR']
        self.jet_pt_bins = config['jet_pt_bins']
        self.eta_max = config['eta_max']
        self.jet_matching_distance = config['jet_matching_distance']
        
        # Initialize thermal model
        if 'thermal_model' in config:
          self.thermal_model = True
          self.thermal_generator = thermal_generator.ThermalGenerator(N_avg=config['thermal_model']['N_avg'],
                                                                      sigma_N=config['thermal_model']['sigma_N'],
                                                                      beta=config['thermal_model']['beta'],
                                                                      eta_max=self.eta_max)
        else:
            self.thermal_model = False
                    
        # Initialize constituent subtractor
        constituent_subtractor = config['constituent_subtractor']
        self.max_distance = constituent_subtractor['max_distance']
        self.alpha = constituent_subtractor['alpha']
        self.bge_rho_grid_size = constituent_subtractor['bge_rho_grid_size']
        self.max_pt_correct = constituent_subtractor['max_pt_correct']
        self.ghost_area = constituent_subtractor['ghost_area']

    #---------------------------------------------------------------
    # Compute pseudorapidity from four-vector
    #---------------------------------------------------------------
    def apply_eta_cut(self, df):
   
        df['eta'] = df.apply(self.eta, axis=1)
        df = df[np.abs(df['eta']) < self.eta_max].drop(['eta'],axis=1)
        return df

    #---------------------------------------------------------------
    # Compute pseudorapidity from four-vector
    #---------------------------------------------------------------
    def eta(self, df):
            
        p = np.sqrt(df['px']*df['px'] + df['py']*df['py'] + df['pz']*df['pz'])
        numerator = p + df['pz']
        denominator = p - df['pz']
        
        if not np.isclose(numerator, 0.) and not np.isclose(denominator, 0.):
            eta = 0.5*np.log( (p + df['pz']) / (p - df['pz']) )
        else:
            eta = 1000
            
        if np.abs(eta) < 1. and df['E'] > 1000.:
            print(df)
            
        return eta

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def process_ppAA(self):
    
        # Loop over events and do jet finding
        # Fill each of the jet_variables into a list
        fj.ClusterSequence.print_banner()
        print('Finding jets and computing N-subjettiness...')
        result = [self.analyze_event(fj_particles_hard, fj_particles_combined) for fj_particles_hard, fj_particles_combined in zip(self.df_fjparticles['fj_particles_hard'], self.df_fjparticles['fj_particles_combined'])]
        
        # Transform the dictionary of lists into a dictionary of numpy arrays
        self.jet_variables_numpy = self.transform_to_numpy(self.jet_variables)
        self.four_vectors_numpy = self.transform_to_numpy(self.four_vectors)
        self.jet_qa_variables_numpy = self.transform_to_numpy(self.jet_qa_variables)
        print()
        
        # Reformat output for ML algorithms (array with 1 array per jet which contain all N-subjettiness values)
        X_Nsub = {'hard': {}, 'combined': {}}
        for label in X_Nsub.keys():
            for jetR in self.jetR_list:
                X_Nsub[label][f'R{jetR}'] = {}
                for jet_pt_bin in self.jet_pt_bins:
                    X_Nsub[label][f'R{jetR}'][f'pt{jet_pt_bin}'] = {}
                    for R_max in self.max_distance:
                        if label=='combined' or np.isclose(R_max, 0.):
                            X_Nsub_reformatted = np.array([list(self.jet_variables_numpy[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'].values())])[0].T
                            X_Nsub[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'] = X_Nsub_reformatted
            
        # Write jet arrays to file
        with h5py.File(os.path.join(self.output_dir, 'nsubjettiness.h5'), 'w') as hf:
            for label in X_Nsub.keys():
                for jetR in self.jetR_list:
                    for jet_pt_bin in self.jet_pt_bins:
                        for R_max in self.max_distance:
                            if label=='combined' or np.isclose(R_max, 0.):
                            
                                suffix = f'_{label}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                        
                                # Write Nsubjettiness
                                hf.create_dataset(f'X_Nsub{suffix}', data=X_Nsub[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'])
                                
                                # Write labels: Pythia 0, Jewel 1
                                if 'jewel' in self.input_file:
                                    y = np.ones(X_Nsub[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'].shape[0])
                                elif 'pythia' in self.input_file:
                                    y = np.zeros(X_Nsub[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'].shape[0])
                                hf.create_dataset(f'y{suffix}', data=y)
                                
                                # Write four-vectors
                                X = self.four_vectors_numpy[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['jet_constituent_four_vectors']
                                hf.create_dataset(f'X_four_vectors{suffix}', data=X)
                                print(label)
                                print(f'R{jetR}')
                                print(f'pt{jet_pt_bin}')
                                print(f'Rmax{R_max}')
                                print(X.shape)
                                
                                for qa_observable in self.qa_observables:
                                    hf.create_dataset(f'{qa_observable}{suffix}', data=self.jet_qa_variables_numpy[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'][qa_observable])
                                    data=self.jet_qa_variables_numpy[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'][qa_observable]
  
            hf.create_dataset('N_list', data=self.N_list)
            hf.create_dataset('beta_list', data=self.beta_list)
            hf.create_dataset('delta_pt_random_cone', data=self.delta_pt_random_cone)
                        
    #---------------------------------------------------------------
    # Process an event (in this case, just a single jet per event)
    #---------------------------------------------------------------
    def analyze_event(self, fj_particles_hard, fj_particles_combined_beforeCS):
    
        # Check that the entries exist appropriately
        if type(fj_particles_hard) != fj.vectorPJ:
            print('fj_particles type mismatch -- skipping event')
            return
        if not self.thermal_model and type(fj_particles_combined_beforeCS) != fj.vectorPJ:
            print('fj_particles type mismatch -- skipping event')
            return
            
        # If thermal model, generate a thermal event and add it to the combined particle list
        if self.thermal_model:
          fj_particles_background = self.thermal_generator.load_event()
          
          # Form the combined event
          # The hard event tracks are each stored with a unique user_index >= 0
          # The thermal tracks are each stored with a unique user_index < 0
          [fj_particles_combined_beforeCS.push_back(p) for p in fj_particles_background]
          
        # Compute delta-pt by random cone method
        self.delta_pt_RC(fj_particles_combined_beforeCS)
            
        # Perform constituent subtraction for each R_max
        fj_particles_combined = []
        for i, R_max in enumerate(self.max_distance):
            if R_max == 0:
                fj_particles_combined.append(fj_particles_combined_beforeCS)
            else:
                fj_particles_combined.append(self.constituent_subtractor[i].process_event(fj_particles_combined_beforeCS))
        
        # Loop through jetR, and process event for each R
        for jetR in self.jetR_list:
             
            for jet_pt_bin in self.jet_pt_bins:
                min_jet_pt = jet_pt_bin[0]
                max_jet_pt = jet_pt_bin[1]
                
                # Set jet definition and a jet selector
                jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
                jet_selector = fj.SelectorPtMin(min_jet_pt) & fj.SelectorPtMax(max_jet_pt) & fj.SelectorAbsRapMax(self.eta_max - jetR)
            
                for i, R_max in enumerate(self.max_distance):
                    #print()
                    #print(R_max)
                    #print('Total number of combined particles: {}'.format(len([p.pt() for p in fj_particles_combined_beforeCS])))
                    #print('After constituent subtraction {}: {}'.format(i, len([p.pt() for p in fj_particles_combined[i]])))

                    # Do jet finding (re-do each time, to make sure matching info gets reset)
                    cs_combined = fj.ClusterSequence(fj_particles_combined[i], jet_def)
                    jets_combined = fj.sorted_by_pt(cs_combined.inclusive_jets())
                    jets_combined_selected = jet_selector(jets_combined)
                    cs_hard = fj.ClusterSequence(fj_particles_hard, jet_def)
                    jets_hard = fj.sorted_by_pt(cs_hard.inclusive_jets())
                    jets_hard_selected = jet_selector(jets_hard)
                    
                    self.analyze_jets(jets_combined_selected, jets_hard_selected, jetR, jet_pt_bin, R_max = R_max)

    #---------------------------------------------------------------
    # Compute delta-pt by random cone method
    #---------------------------------------------------------------
    def delta_pt_RC(self, fj_particles_combined_beforeCS):
    
        event_pt = 0.
        cone_pt = 0.
        R_cone = 0.4
        eta_random = np.random.uniform(-self.eta_max+R_cone, self.eta_max-R_cone)
        phi_random = np.random.uniform(0, 2*np.pi)
        for particle in fj_particles_combined_beforeCS:
            event_pt += particle.pt()
            delta_R = self.utils.delta_R(particle, eta_random, phi_random)
            if delta_R < R_cone:
                cone_pt += particle.pt()
        rho = event_pt / (2*self.eta_max *2*np.pi)
        delta_pt = cone_pt - rho*np.pi*R_cone*R_cone
        
        self.delta_pt_random_cone.append(delta_pt)

    #---------------------------------------------------------------
    # Analyze jets of a given event.
    #---------------------------------------------------------------
    def analyze_jets(self, jets_combined_selected, jets_hard_selected, jetR, jet_pt_bin, R_max = None):

        # Fill combined jet info
        for jet_combined in jets_combined_selected:
            self.fill_nsubjettiness(jet_combined, jetR, jet_pt_bin, R_max, 'combined')
            self.fill_four_vectors(jet_combined, jetR, jet_pt_bin, R_max, 'combined')
            
        # Fill hard jet info
        if np.isclose(R_max, 0.):
            for jet_hard in jets_hard_selected:
                self.fill_nsubjettiness(jet_hard, jetR, jet_pt_bin, R_max, 'hard')
                self.fill_four_vectors(jet_hard, jetR, jet_pt_bin, R_max, 'hard')
                        
        #----------------------------------
        # Match jets

        # Loop through jets and set jet matching candidates for each jet in user_info
        for jet_combined in jets_combined_selected:
            for jet_hard in jets_hard_selected:
                
                # Add a matching candidate to the list if it is within the geometrical cut
                deltaR = jet_combined.delta_R(jet_hard)
                if deltaR < self.jet_matching_distance*jetR:
                    process_base.ProcessBase.set_jet_info(None, jet_combined, jet_hard, deltaR)
                    process_base.ProcessBase.set_jet_info(None, jet_hard, jet_combined, deltaR)
        
        # Loop through jets and set accepted matches
        for jet_combined in jets_combined_selected:
            process_base.ProcessBase.set_matches_pp(None, jet_combined, None)
              
        # Loop through jets and fill response histograms if both det and truth jets are unique match
        result = [self.fill_jet_matches(jet_combined, jetR, jet_pt_bin, R_max, 'combined') for jet_combined in jets_combined_selected]
                
    #---------------------------------------------------------------
    # Analyze jets of a given event.
    #---------------------------------------------------------------
    def fill_jet_matches(self, jet_combined, jetR, jet_pt_bin, R_max, label):
    
        # Get matched pp jet
        if jet_combined.has_user_info():
          jet_pp = jet_combined.python_info().match
          
          if jet_pp:
            jet_pt_combined = jet_combined.pt()
            jet_pt_pp = jet_pp.pt()
            
            delta_pt = (jet_pt_combined - jet_pt_pp)
            self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['delta_pt'].append(delta_pt)
    
    #---------------------------------------------------------------
    # Analyze jets of a given event.
    #---------------------------------------------------------------
    def fill_nsubjettiness(self, jet, jetR, jet_pt_bin, R_max = None, label = ''):

        # Compute N-subjettiness
        axis_definition = fjcontrib.KT_Axes()
        for i,N in enumerate(self.N_list):
            beta = self.beta_list[i]
        
            measure_definition = fjcontrib.UnnormalizedMeasure(beta)
            n_subjettiness_calculator = fjcontrib.Nsubjettiness(N, axis_definition, measure_definition)
            n_subjettiness = n_subjettiness_calculator.result(jet)/jet.pt()
            self.jet_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'][f'n_subjettiness_N{N}_beta{beta}'].append(n_subjettiness)
            
        # Fill some jet QA
        self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['jet_pt'].append(jet.pt())
        
        # angularity
        alpha = 1
        kappa = 1
        angularity = fjext.lambda_beta_kappa(jet, alpha, kappa, jetR)
        self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['jet_angularity'].append(angularity)
        
        # mass
        self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['jet_mass'].append(jet.m())
        
        # theta_g
        beta = 0
        zcut = 0.2
        gshop = fjcontrib.GroomerShop(jet, jetR, fj.cambridge_algorithm)
        jet_groomed_lund = gshop.soft_drop(beta, zcut, jetR)
        theta_g = jet_groomed_lund.Delta() / jetR
        self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['jet_theta_g'].append(theta_g)
        
        # subjet z
        subjetR = 0.1
        subjet_def = fj.JetDefinition(fj.antikt_algorithm, subjetR)
        cs_subjet = fj.ClusterSequence(jet.constituents(), subjet_def)
        subjets = fj.sorted_by_pt(cs_subjet.inclusive_jets())
        leading_subjet = self.utils.leading_jet(subjets)
        z_leading = leading_subjet.pt() / jet.pt()
        self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['jet_subjet_z'].append(z_leading)
        
        # leading hadron z
        leading_particle = self.utils.leading_jet(jet.constituents())
        z_leading = leading_particle.pt() / jet.pt()
        self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['hadron_z'].append(z_leading)

    #---------------------------------------------------------------
    # Write four-vectors of jet constituents
    #---------------------------------------------------------------
    def fill_four_vectors(self, jet, jetR, jet_pt_bin, R_max = None, label = ''):

        # Loop through jet constituents
        particle_list = []
        for i,particle in enumerate(jet.constituents()):
            particle_list.append(np.array([particle.E(), particle.px(), particle.py(), particle.pz()]))
        
        # Zero pad such that all jets have the same number of four-vectors
        n_max = 800
        if len(particle_list) > n_max:
            sys.exit('ERROR: particle list has {len(particle_list)} entries before zero-padding')
        particle_list += [np.array([0,0,0,0])]*(n_max-len(particle_list))
        
        # Append list of four-vectors to output
        self.four_vectors[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['jet_constituent_four_vectors'].append(particle_list)
          
    #---------------------------------------------------------------
    # Transform dictionary of lists into a dictionary of numpy arrays
    #---------------------------------------------------------------
    def transform_to_numpy(self, jet_variables_list):

        jet_variables_numpy = {'hard': {}, 'combined': {}}

        for label in jet_variables_numpy.keys():
            for jetR in self.jetR_list:
                jet_variables_numpy[label][f'R{jetR}'] = {}
                for jet_pt_bin in self.jet_pt_bins:
                    jet_variables_numpy[label][f'R{jetR}'][f'pt{jet_pt_bin}'] = {}
               
                    for R_max in self.max_distance:
                        if label=='combined' or np.isclose(R_max, 0.):

                            jet_variables_numpy[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'] = {}
                  
                            for key,val in jet_variables_list[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'].items():
                                jet_variables_numpy[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'][key] = np.array(val)
                            
        return jet_variables_numpy
            
    #---------------------------------------------------------------
    # Construct fastjet::PseudoJets from four-vectors
    #---------------------------------------------------------------
    def get_fjparticles(self, df_particles_grouped):
                                                 
        user_index_offset = 0
        return fjext.vectorize_px_py_pz_e(df_particles_grouped['px'].values,
                                          df_particles_grouped['py'].values,
                                          df_particles_grouped['pz'].values,
                                          df_particles_grouped['E'].values,
                                          user_index_offset)        

##################################################################
if __name__ == '__main__':

    # Define arguments
    parser = argparse.ArgumentParser(description='Process pp AA')
    parser.add_argument('-c', '--configFile', action='store',
                        type=str, metavar='configFile',
                        default='./config/ml/ppAA.yaml',
                        help='Path of config file for analysis')
    parser.add_argument('-f', '--inputFile', action='store',
                        type=str, metavar='inputFile',
                default='./skim_blah.pkl',
                        help='Path of input file for analysis')
    parser.add_argument('-o', '--outputDir', action='store',
                        type=str, metavar='outputDir',
                        default='./TestOutput',
                        help='Output directory for output to be written to')

    # Parse the arguments
    args = parser.parse_args()

    print('Configuring...')
    print('configFile: \'{0}\''.format(args.configFile))
    print('inputFile: \'{0}\''.format(args.inputFile))
    print('ouputDir: \'{0}\"'.format(args.outputDir))

    # If invalid configFile is given, exit
    if not os.path.exists(args.configFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    # If invalid inputFile is given, exit
    if not os.path.exists(args.inputFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
        sys.exit(0)

    analysis = ProcessppAA(config_file=args.configFile, input_file=args.inputFile, output_dir=args.outputDir)
    analysis.process_ppAA()
