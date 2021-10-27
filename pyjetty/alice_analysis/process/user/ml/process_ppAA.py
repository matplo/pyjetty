#!/usr/bin/env python3

"""
Class to read pp vs AA data set, do jet finding, and compute N-subjettiness
"""

import os
import sys
import argparse
import yaml
import h5py
import time

# Data analysis and plotting
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext
import fjtools

# Energy flow package
try:
    import energyflow
except ImportError:
    pass

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_io_emb
from pyjetty.alice_analysis.process.base import process_base
from pyjetty.alice_analysis.process.base import process_utils
from pyjetty.alice_analysis.process.base import thermal_generator
from pyjetty.mputils import CEventSubtractor

# Base class
from pyjetty.alice_analysis.process.base import common_base

################################################################
class ProcessppAA(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, config_file='', input_file='', output_dir='', test=False, **kwargs):
        super(common_base.CommonBase, self).__init__(**kwargs)
       
        self.start_time = time.time()
        
        self.test = test
        self.config_file = config_file
        self.input_file = input_file
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        
        # Initialize config file
        self.initialize_config()
        
        # Initialize utils class
        self.utils = process_utils.ProcessUtils()

        #---------------------------------------------------------------
        # Pb-Pb ALICE data
        if self.PbPb_data:

            # Use IO helper class to convert ROOT TTree into a SeriesGroupBy object of fastjet particles per event
            io = process_io.ProcessIO(input_file=self.input_file, track_tree_name='tree_Particle',
                                    is_pp=False, use_ev_id_ext=True)
            self.df_fjparticles = io.load_data(min_pt=self.min_pt)
            print('Done.')

        #---------------------------------------------------------------
        # pp ALICE data -- embed into Pb-Pb background
        elif self.pp_data:

            # Use IO helper class to convert pp ROOT TTree into a SeriesGroupBy object of fastjet particles per event
            io_pp = process_io.ProcessIO(input_file=self.input_file, track_tree_name='tree_Particle', 
                                           is_pp=True, use_ev_id_ext=True)
            self.df_fjparticles = io_pp.load_data(min_pt=self.min_pt, reject_tracks_fraction=self.reject_tracks_fraction)
            self.df_fjparticles.columns = ['fj_particles_combined']
            
            # Set up the Pb-Pb embedding object
            self.process_io_emb = process_io_emb.ProcessIO_Emb(self.emb_file_list, track_tree_name='tree_Particle')

        #---------------------------------------------------------------
        # MC -- here we want to have two dataframes, one for the hard event, and one for the combined event
        else:

            #---------------------------------------------------------------
            # If test mode, load quark-gluon dataset -- otherwise, load ROOT TTree of generated events
            if self.test:
                # https://energyflow.network/docs/datasets/#quark-and-gluon-jets
                # X : a three-dimensional numpy array of jets:
                #     list of jets with list of particles for each jet, with (pt,y,phi,pid) values for each particle
                # y : a numpy array of quark/gluon jet labels (quark=1 and gluon=0).
                # The jets are padded with zero-particles in order to make a contiguous array.
                print()
                print('Loading qg dataset:')
                self.X, self.y = energyflow.datasets.qg_jets.load(self.n_total)
                print('(n_jets, n_particles per jet, n_variables): {}'.format(self.X.shape))
                print()

                # Next, we will transform these into fastjet::PseudoJet objects.
                # This allows us to use the fastjet contrib to compute Nsubjetiness, and in general it
                # will be needed in order to perform jet finding on an event (in data or MC).

                # Translate 3D numpy array (100,000 x 556 particles x 4 vars) into a dataframe
                # Define a unique index for each jet
                columns = ['pt', 'y', 'phi', 'pid']
                df_particles = pd.DataFrame(self.X.reshape(-1, 4), columns=columns)
                df_particles.index = np.repeat(np.arange(self.X.shape[0]), self.X.shape[1]) + 1
                df_particles.index.name = 'jet_id'

                # Construct dummy combined event -- we will add thermal particles in event loop
                #df_particles_copy = df_particles.copy(deep=True)
                
                # (i) Group the particle dataframe by jet id
                #     df_particles_grouped is a DataFrameGroupBy object with one particle dataframe per jet
                df_fjparticles_grouped_hard = df_particles.groupby('jet_id')
                #df_fjparticles_grouped_combined = df_particles_copy.groupby('jet_id')
                
                # (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet::PseudoJets
                # NOTE: for now we neglect the mass -- and assume y=eta
                # TO DO: Add y to https://github.com/matplo/heppy/blob/master/cpptools/src/fjext/fjtools.cxx
                # TO DO: Add mass vector using pdg
                print('Converting particle dataframe to fastjet::PseudoJets...')
                df_fjparticles_hard = df_fjparticles_grouped_hard.apply(self.get_fjparticles)
                df_fjparticles_combined = df_fjparticles_hard
                print('Done.')
                print()

            #---------------------------------------------------------------
            else:
                # Load dataframe of particle four-vectors for all particles in the event
                # (separate dataframes for hard process and background)
                # Use IO helper class to convert truth-level ROOT TTree into
                # a SeriesGroupBy object of fastjet particles per event
                tree_dir = 'PWGHF_TreeCreator'
                skip_event_tree = 'pyquen' in self.input_file or 'pythia_gen_qorg' in self.input_file
                if skip_event_tree:
                    tree_dir = ''
                io_hard = process_io.ProcessIO(input_file=self.input_file, tree_dir=tree_dir,
                                                track_tree_name='tree_Particle_gen', use_ev_id_ext=False,
                                                skip_event_tree=skip_event_tree)
                df_fjparticles_hard = io_hard.load_data(min_pt=self.min_pt)
                self.nEvents_truth = len(df_fjparticles_hard.index)
                self.nTracks_truth = len(io_hard.track_df.index)
                print('--- {} seconds ---'.format(time.time() - self.start_time))

                # Construct dummy combined event -- we will add thermal particles in event loop
                # (Note: deep copy of df_fjparticles does not copy underlying objects)
                df_fjparticles_combined = io_hard.load_data(min_pt=self.min_pt)

            # Merge hard event and background dataframes
            if self.thermal_model:
                self.df_fjparticles = pd.concat([df_fjparticles_hard, df_fjparticles_combined], axis=1)
                self.df_fjparticles.columns = ['fj_particles_hard', 'fj_particles_combined']
                print('Done.')
                print()
            else:
                print('Only thermal model is supported at the moment.')
        
        # Create list of N-subjettiness observables: number of axes and beta values
        self.N_list = []
        self.beta_list = []
        for i in range(self.K-2):
            self.N_list += [i+1] * 3
            self.beta_list += [0.5,1,2]
        self.N_list += [self.K-1] * 2  
        self.beta_list += [1,2]
        
        # Construct dictionary to store all jet quantities of interest
        self.jet_variables = {'hard': {}, 'combined': {}, 'combined_matched': {}}
        self.four_vectors = {'hard': {}, 'combined': {}, 'combined_matched': {}}
        self.jet_qa_variables = {'hard': {}, 'combined': {}, 'combined_matched': {}}
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
                        
                        if 'combined' in label or np.isclose(R_max, 0.):
                            self.jet_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'] = {}
                            self.four_vectors[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'] = {}
                            self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'] = {}
                            for i,N in enumerate(self.N_list):
                                beta = self.beta_list[i]
                                self.jet_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'][f'n_subjettiness_N{N}_beta{beta}'] = []
                            self.four_vectors[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['jet_constituent_four_vectors'] = []
                            
                            # Some QA
                            self.qa_observables = ['delta_pt', 'matched_pt', 'matched_deltaR', 'jet_pt', 'jet_angularity', 'jet_mass', 'jet_theta_g', 'jet_subjet_z', 'hadron_z', 'multiplicity_0000', 'multiplicity_0150', 'multiplicity_0500', 'multiplicity_1000']
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
          
        self.K = config['K_max']
        
        self.jetR_list = config['jetR']
        self.jet_pt_bins = config['jet_pt_bins']
        self.photon_jet = config['photon_jet']
        self.min_pt = config['min_pt']
        self.leading_track_pt = config['leading_track_pt']
        self.eta_max = config['eta_max']
        self.jet_matching_distance = config['jet_matching_distance']
        self.n_total = config['n_train'] + config['n_val'] + config['n_test']
        self.event_index = 0
        
        # Initialize thermal model
        if 'thermal_model' in config:
          self.thermal_model = True
          self.thermal_generator = thermal_generator.ThermalGenerator(N_avg=config['thermal_model']['N_avg'],
                                                                      sigma_N=config['thermal_model']['sigma_N'],
                                                                      beta=config['thermal_model']['beta'],
                                                                      eta_max=self.eta_max)
        else:
            self.thermal_model = False

        # Determine input data type
        self.PbPb_data = 'LHC18qr' in self.input_file
        self.pp_data = 'LHC17pq' in self.input_file
        self.is_data = self.pp_data or self.PbPb_data
        self.is_mc = not self.is_data
            
        if self.pp_data:
            self.emb_file_list = config['emb_file_list']
            self.reject_tracks_fraction = config['reject_tracks_fraction']
                    
        # Initialize constituent subtractor
        constituent_subtractor = config['constituent_subtractor']
        self.max_distance = constituent_subtractor['max_distance']
        self.alpha = constituent_subtractor['alpha']
        self.bge_rho_grid_size = constituent_subtractor['bge_rho_grid_size']
        self.max_pt_correct = constituent_subtractor['max_pt_correct']
        self.ghost_area = constituent_subtractor['ghost_area']

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def process_ppAA(self):
    
        # Loop over events and do jet finding
        # Fill each of the jet_variables into a list
        fj.ClusterSequence.print_banner()
        print('Finding jets and computing N-subjettiness...')
        if self.is_data:
            result = [self.analyze_event(None, fj_particles_combined) for fj_particles_combined in self.df_fjparticles]        
        else:
            result = [self.analyze_event(fj_particles_hard, fj_particles_combined) for fj_particles_hard, fj_particles_combined in zip(self.df_fjparticles['fj_particles_hard'], self.df_fjparticles['fj_particles_combined'])]        
        
        # Transform the dictionary of lists into a dictionary of numpy arrays
        self.jet_variables_numpy = self.transform_to_numpy(self.jet_variables)
        self.four_vectors_numpy = self.transform_to_numpy(self.four_vectors)
        self.jet_qa_variables_numpy = self.transform_to_numpy(self.jet_qa_variables)
        print()
        
        # Reformat output for ML algorithms (array with 1 array per jet which contain all N-subjettiness values)
        if self.test:
            X_Nsub = {'hard': {}}
        else:
            X_Nsub = {'hard': {}, 'combined': {}, 'combined_matched': {}}
        for label in X_Nsub.keys():
            for jetR in self.jetR_list:
                X_Nsub[label][f'R{jetR}'] = {}
                for jet_pt_bin in self.jet_pt_bins:
                    X_Nsub[label][f'R{jetR}'][f'pt{jet_pt_bin}'] = {}
                    for R_max in self.max_distance:
                        if 'combined' in label or np.isclose(R_max, 0.):
                            X_Nsub_reformatted = np.array([list(self.jet_variables_numpy[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'].values())])[0].T
                            X_Nsub[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'] = X_Nsub_reformatted
            
        # Write jet arrays to file
        with h5py.File(os.path.join(self.output_dir, 'nsubjettiness.h5'), 'w') as hf:
            for label in X_Nsub.keys():
                for jetR in self.jetR_list:
                    for jet_pt_bin in self.jet_pt_bins:
                        for R_max in self.max_distance:
                            if 'combined' in label or np.isclose(R_max, 0.):
                            
                                suffix = f'_{label}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                        
                                # Write Nsubjettiness
                                hf.create_dataset(f'X_Nsub{suffix}', data=X_Nsub[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'])

                                # Check whether any training entries are empty
                                [print(f'WARNING: input entry {i} is empty') for i,x in enumerate(X_Nsub[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']) if not x.any()]
  
                                # Write labels: Pythia 0, Jewel 1
                                if self.PbPb_data:
                                    self.y = np.ones(X_Nsub[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'].shape[0])
                                elif self.pp_data:
                                    self.y = np.zeros(X_Nsub[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'].shape[0])
                                else:
                                    if not self.test:
                                        if 'jewel_PbPb' in self.input_file or 'pyquen' in self.input_file or 'qorg_glue' in self.input_file:
                                            self.y = np.ones(X_Nsub[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'].shape[0])
                                        else:
                                            self.y = np.zeros(X_Nsub[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'].shape[0])
                                hf.create_dataset(f'y{suffix}', data=self.y)
                                
                                # Write four-vectors
                                X = self.four_vectors_numpy[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['jet_constituent_four_vectors']
                                hf.create_dataset(f'X_four_vectors{suffix}', data=X)
                                print(label)
                                print(f'R{jetR}')
                                print(f'pt{jet_pt_bin}')
                                print(f'Rmax{R_max}')
                                print(X.shape)

                                # If test, compare four-vectors after jet-finding to those before jet-finding
                                if self.test:
                                    print(f'four-vectors before jet finding: {self.X}')
                                    print(f'four-vectors after jet finding: {X}')
                                
                                for qa_observable in self.qa_observables:
                                    hf.create_dataset(f'{qa_observable}{suffix}', data=self.jet_qa_variables_numpy[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'][qa_observable])

                                # Make some QA plots
                                self.output_dir_i = os.path.join(self.output_dir, f'{label}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}')
                                if not os.path.exists(self.output_dir_i):
                                    os.makedirs(self.output_dir_i)
                                self.plot_QA(label, jetR, jet_pt_bin, R_max)
  
            hf.create_dataset('N_list', data=self.N_list)
            hf.create_dataset('beta_list', data=self.beta_list)
            hf.create_dataset('delta_pt_random_cone', data=self.delta_pt_random_cone)
                            
    #---------------------------------------------------------------
    # Process an event
    #---------------------------------------------------------------
    def analyze_event(self, fj_particles_hard, fj_particles_combined_beforeCS):
    
        # Check that the entries exist appropriately
        if fj_particles_hard and type(fj_particles_hard) != fj.vectorPJ:
            print('fj_particles type mismatch -- skipping event')
            return
        if not self.thermal_model and type(fj_particles_combined_beforeCS) != fj.vectorPJ:
            print('fj_particles type mismatch -- skipping event')
            return
        
        # Embed if necessary
        if self.pp_data or self.thermal_model and not self.test:

            # If pp data, get Pb-Pb event and embed it into pp event
            # The pp tracks are each stored with a unique user_index >= 0
            # The Pb-Pb tracks are each stored with a unique user_index < 0
            if self.pp_data:
                fj_particles_background = [p for p in self.process_io_emb.load_event() if p.pt() > self.min_pt]

            # If thermal model, generate a thermal event and add it to the combined particle list
            elif self.thermal_model:
                fj_particles_background = [p for p in self.thermal_generator.load_event() if p.pt() > self.min_pt]
          
            # Form the combined event
            # The hard event tracks are each stored with a unique user_index >= 0
            # The thermal tracks are each stored with a unique user_index < 0
            [fj_particles_combined_beforeCS.push_back(p) for p in fj_particles_background]

            # Compute delta-pt by random cone method
            self.delta_pt_RC(fj_particles_background)

        # Perform constituent subtraction for each R_max
        if not self.test:

            fj_particles_combined = []
            for i, R_max in enumerate(self.max_distance):
                if R_max == 0:
                    fj_particles_combined.append(fj_particles_combined_beforeCS)
                else:
                    fj_particles_combined.append(self.constituent_subtractor[i].process_event(fj_particles_combined_beforeCS))

        # Loop through jetR, and process event for each R
        for jetR in self.jetR_list:
             
            for jet_pt_bin in self.jet_pt_bins:
                
                # Set jet definition and a jet selector
                # For the hard jets, they should satisfy the pt interval
                # For the combined jets, they can go outside, since they will be matched to hard jets
                jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
                if self.photon_jet:
                    min_photon_E = jet_pt_bin[0]
                    min_jet_pt = jet_pt_bin[1]

                    # Check for photon
                    found_photon = False
                    for particle in fj_particles_hard:
                        if particle.pid() == 22 and particle.pt() > min_photon_E:
                            found_photon = True
                            break

                    # Define jets
                    jet_selector_hard = fj.SelectorPtMin(min_jet_pt) & fj.SelectorAbsRapMax(self.eta_max - jetR)
                    jet_selector_combined = fj.SelectorPtMin(min_jet_pt/5.) & fj.SelectorAbsRapMax(self.eta_max - jetR)

                else:
                    min_jet_pt = jet_pt_bin[0]
                    max_jet_pt = jet_pt_bin[1]

                    jet_selector_hard = fj.SelectorPtMin(min_jet_pt) & fj.SelectorPtMax(max_jet_pt) & fj.SelectorAbsRapMax(self.eta_max - jetR)
                    if self.is_data:
                        jet_selector_combined = fj.SelectorPtMin(min_jet_pt) & fj.SelectorAbsRapMax(self.eta_max - jetR)
                    elif self.is_mc:
                        jet_selector_combined = fj.SelectorPtMin(min_jet_pt/5.) & fj.SelectorAbsRapMax(self.eta_max - jetR)
            
                for i, R_max in enumerate(self.max_distance):
                    #print()
                    #print(R_max)
                    #print('Total number of combined particles: {}'.format(len([p.pt() for p in fj_particles_combined_beforeCS])))
                    #print('After constituent subtraction {}: {}'.format(i, len([p.pt() for p in fj_particles_combined[i]])))

                    # Do jet finding (re-do each time, to make sure matching info gets reset)
                    if self.test:
                        # Cluster each jet with R=infinity
                        jet_def = fj.JetDefinition(fj.antikt_algorithm, fj.JetDefinition.max_allowable_R)
                        cs = fj.ClusterSequence(fj_particles_hard, jet_def)
                        jets_hard_selected = [fj.sorted_by_pt(cs.inclusive_jets())[0]]
                        jets_combined_selected = []
                    else:

                        cs_combined = fj.ClusterSequence(fj_particles_combined[i], jet_def)
                        jets_combined = fj.sorted_by_pt(cs_combined.inclusive_jets())
                        jets_combined_selected = jet_selector_combined(jets_combined)

                        if self.is_mc:
                            cs_hard = fj.ClusterSequence(fj_particles_hard, jet_def)
                            jets_hard = fj.sorted_by_pt(cs_hard.inclusive_jets())
                            jets_hard_selected = jet_selector_hard(jets_hard)
                        else:
                            jets_hard_selected = None
        
                    self.analyze_jets(jets_combined_selected, jets_hard_selected, jetR, jet_pt_bin, R_max = R_max)

        self.event_index += 1
        if self.event_index%100 == 0:
            print(f'event: {self.event_index}')

    #---------------------------------------------------------------
    # Compute delta-pt by random cone method
    #---------------------------------------------------------------
    def delta_pt_RC(self, fj_particles_background):
    
        event_pt = 0.
        cone_pt = 0.
        R_cone = 0.4
        eta_random = np.random.uniform(-self.eta_max+R_cone, self.eta_max-R_cone)
        phi_random = np.random.uniform(0, 2*np.pi)
        for particle in fj_particles_background:
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

        # If photon-jet, look for leading jet that passes criteria
        if self.photon_jet:

            # Find photon
            min_photon_E = jet_pt_bin[0]
            photon = None
            for particle in fj_particles_hard:
                if particle.pid() == 22 and particle.pt() > min_photon_E:
                    photon = particle
            
            # Check that jet is in opposite hemisphere, otherwise remove it
            jets_hard_selected = [jet for jet in jets_hard_selected if np.abs(photon.Delta_R(jet)) > np.pi/2]

        # Fill hard jet info
        if self.is_mc:
            if np.isclose(R_max, 0.):
                for jet_hard in jets_hard_selected:
                    self.fill_nsubjettiness(jet_hard, jetR, jet_pt_bin, R_max, 'hard')
                    self.fill_four_vectors(jet_hard, jetR, jet_pt_bin, R_max, 'hard')

        if not self.test:

            # Fill combined jet info, if inside pt interval
            for jet_combined in jets_combined_selected:

                # In data, how to select jets? In particular in pp data -- want to require that combined jet contains pp jet
                # For now: leading track requirement (after constituent subtraction)
                if self.is_data:
                    accept_jet = False
                    for p in jet_combined.constituents():
                        if self.pp_data:
                            if p.pt() > self.leading_track_pt and p.user_index() >= 0:
                                accept_jet = True  
                        elif self.PbPb_data:
                            if p.pt() > self.leading_track_pt:
                                accept_jet = True  
                    if not accept_jet:
                        continue 

                if jet_pt_bin[0] < jet_combined.pt() < jet_pt_bin[1]:
                    self.fill_nsubjettiness(jet_combined, jetR, jet_pt_bin, R_max, 'combined')
                    self.fill_four_vectors(jet_combined, jetR, jet_pt_bin, R_max, 'combined')
                        
            #----------------------------------
            # Match jets
            if self.is_mc:
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
                    
                # Loop through jets and fill output if both det and truth jets are unique match
                result = [self.fill_jet_matches(jet_combined, jetR, jet_pt_bin, R_max, 'combined_matched') for jet_combined in jets_combined_selected]
                
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
            
            # Accept jet if hard jet is within pt bin
            if jet_pt_bin[0] < jet_pt_pp < jet_pt_bin[1]:
                
                delta_pt = (jet_pt_combined - jet_pt_pp)
                self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['delta_pt'].append(delta_pt)
                
                self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['matched_deltaR'].append(jet_combined.delta_R(jet_pp))
                matched_pt = fjtools.matched_pt(jet_combined, jet_pp)
                self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['matched_pt'].append(matched_pt)

                self.fill_nsubjettiness(jet_combined, jetR, jet_pt_bin, R_max, label)
                self.fill_four_vectors(jet_combined, jetR, jet_pt_bin, R_max, label)
    
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
        
        # multiplicity
        n_constituents = len(jet.constituents())
        self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['multiplicity_0000'].append(n_constituents)
        multiplicity_0150 = 0
        multiplicity_0500 = 0
        multiplicity_1000 = 0
        for constituent in jet.constituents():
            if constituent.pt() > 0.15:
                multiplicity_0150 += 1
            if constituent.pt() > 0.5:
                multiplicity_0500 += 1
            if constituent.pt() > 1.:
                multiplicity_1000 += 1
        self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['multiplicity_0150'].append(multiplicity_0150)
        self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['multiplicity_0500'].append(multiplicity_0500)
        self.jet_qa_variables[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}']['multiplicity_1000'].append(multiplicity_1000)
        
    #---------------------------------------------------------------
    # Write four-vectors of jet constituents
    #---------------------------------------------------------------
    def fill_four_vectors(self, jet, jetR, jet_pt_bin, R_max = None, label = ''):

        # Loop through jet constituents
        particle_list = []
        for i,particle in enumerate(jet.constituents()):
            #particle_list.append(np.array([particle.E(), particle.px(), particle.py(), particle.pz()]))
            particle_list.append(np.array([particle.perp(), particle.rapidity(), particle.phi_02pi(), 0.]))
        
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

        if self.test:
            jet_variables_numpy = {'hard': {}}
        else:
            jet_variables_numpy = {'hard': {}, 'combined': {}, 'combined_matched': {}}

        for label in jet_variables_numpy.keys():
            for jetR in self.jetR_list:
                jet_variables_numpy[label][f'R{jetR}'] = {}
                for jet_pt_bin in self.jet_pt_bins:
                    jet_variables_numpy[label][f'R{jetR}'][f'pt{jet_pt_bin}'] = {}
               
                    for R_max in self.max_distance:
                        if 'combined' in label or np.isclose(R_max, 0.):

                            jet_variables_numpy[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'] = {}
                  
                            for key,val in jet_variables_list[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'].items():
                                jet_variables_numpy[label][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'][key] = np.array(val)
                            
        return jet_variables_numpy

    #---------------------------------------------------------------
    # Plot QA
    #---------------------------------------------------------------
    def plot_QA(self, event_type, jetR, jet_pt_bin, R_max):
    
        for qa_observable in self.qa_observables:
            
            qa_result = self.jet_qa_variables_numpy[event_type][f'R{jetR}'][f'pt{jet_pt_bin}'][f'Rmax{R_max}'][qa_observable]
            qa_observable_shape = qa_result.shape
            if qa_observable_shape[0] == 0:
                continue

            # Plot distributions
            plt.xlabel(rf'{qa_observable}', fontsize=14)
            max = np.amax(qa_result)*1.2
            bins = np.linspace(0, max, 20)
            plt.hist(qa_result,
                     bins,
                     histtype='step',
                     density=True,
                     label = 'pp or AA',
                     linewidth=2,
                     linestyle='-',
                     alpha=0.5)
            plt.legend(loc='best', fontsize=14, frameon=False)
            
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir_i, f'{qa_observable}.pdf'))
            plt.close()   
                                        
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
    # Transform particles to fastjet::PseudoJets
    #---------------------------------------------------------------
    def get_fjparticles(self, df_particles_grouped):

        user_index_offset = 0
        return fjext.vectorize_pt_eta_phi(df_particles_grouped['pt'].values,
                                          df_particles_grouped['y'].values,
                                          df_particles_grouped['phi'].values,
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
                        default='./jewel.root',
                        help='Path of input file for analysis')
    parser.add_argument('-o', '--outputDir', action='store',
                        type=str, metavar='outputDir',
                        default='./TestOutput',
                        help='Output directory for output to be written to')
    parser.add_argument('--test', 
                        help='test mode using qg dataset', 
                        action='store_true', default=False)

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
    if not os.path.exists(args.inputFile) and not args.test:
        print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
        sys.exit(0)

    analysis = ProcessppAA(config_file=args.configFile, input_file=args.inputFile, output_dir=args.outputDir, test=args.test)
    analysis.process_ppAA()
