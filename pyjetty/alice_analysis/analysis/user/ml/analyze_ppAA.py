#!/usr/bin/env python3

"""
Example class to read quark-gluon dataset
"""

import os
import sys
import argparse
import yaml
import h5py
import pickle
import subprocess
from numba import jit, prange
import functools
import shutil

# Data analysis and plotting
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_context('paper', rc={'font.size':18,'axes.titlesize':18,'axes.labelsize':18})

# Energy flow package
import energyflow
import energyflow.archs

# sklearn
import sklearn
import sklearn.linear_model
import sklearn.ensemble
import sklearn.model_selection
import sklearn.pipeline

# Tensorflow and Keras
import tensorflow as tf
from tensorflow import keras
import keras_tuner

# Base class
from pyjetty.alice_analysis.process.base import common_base

#--------------------------------------------------------------- 
# Create a copy of four-vectors with a min-pt cut
#---------------------------------------------------------------        
@jit(nopython=True) 
def filter_four_vectors(X_particles, min_pt=0.):

    n_jets = X_particles.shape[0]
    n_particles = 800
    
    for i in prange(n_jets):
        jet = X_particles[i]

        new_jet_index = 0
        new_jet = np.zeros(jet.shape)

        for j in prange(n_particles):
            if jet[j][0] > min_pt:
                 new_jet[new_jet_index] = jet[j]
                 new_jet_index += 1
        
        X_particles[i] = new_jet.copy()

    return X_particles

################################################################
class AnalyzePPAA(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, config_file='', output_dir='', old_pp=False, old_AA=False, **kwargs):
        super(common_base.CommonBase, self).__init__(**kwargs)
        
        self.config_file = config_file
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        
        # Settings for tests against old input data
        self.old_pp = old_pp
        self.old_AA = old_AA
        if old_pp or old_AA:
            self.old_input_file = '/home/james/pyjetty/pyjetty/alice_analysis/analysis/user/ml/551797/nsubjettiness_with_four_vectors.h5'

        # Initialize config file
        self.initialize_config()
            
        self.filename = 'nsubjettiness_with_four_vectors_unshuffled.h5'
        with h5py.File(os.path.join(self.output_dir, self.filename), 'r') as hf:
            self.N_list = hf['N_list'][:]
            self.beta_list = hf['beta_list'][:]
            self.delta_pt_random_cone = hf['delta_pt_random_cone'][:]
        
        self.qa_observables = ['matched_pt', 'matched_deltaR', 'jet_pt', 'jet_angularity', 'thrust', 'LHA', 'pTD', 'jet_mass', 'jet_theta_g', 'zg', 'jet_subjet_z', 'hadron_z', 'multiplicity_0000', 'multiplicity_0150', 'multiplicity_0500', 'multiplicity_1000']

        # Remove keras-tuner folder, if it exists
        if os.path.exists('keras_tuner'):
            shutil.rmtree('keras_tuner')

        print(self)
        print()
        
    #---------------------------------------------------------------
    # Initialize config file into class members
    #---------------------------------------------------------------
    def initialize_config(self):
    
        # Read config file
        with open(self.config_file, 'r') as stream:
          config = yaml.safe_load(stream)
          
        self.jetR_list = config['jetR']
        self.jet_pt_bins = config['jet_pt_bins']
        self.max_distance_list = config['constituent_subtractor']['max_distance']
        self.event_types = ['hard', 'combined_matched']
          
        self.n_train = config['n_train']
        self.n_val = config['n_val']
        self.n_test = config['n_test']
        self.n_total = self.n_train + self.n_val + self.n_test
        self.test_frac = 1. * self.n_test / self.n_total
        self.val_frac = 1. * self.n_val / (self.n_train + self.n_val)
        
        self.K_list = config['K']
        self.K_max = config['K_max']
        self.dmax = config['dmax']
        self.efp_measure = config['efp_measure']
        self.efp_beta = config['efp_beta']

        self.random_state = None  # seed for shuffling data (set to an int to have reproducible results)

        self.constituent_subtraction_study = config['constituent_subtraction_study']

        self.pp_label = config['pp_label']
        self.AA_label = config['AA_label']
        
        # Initialize model-specific settings
        self.config = config
        self.models = config['models']
        self.model_settings = {}
        for model in self.models:
            self.model_settings[model] = {}
            
            if 'dnn' in model:
                self.model_settings[model]['loss'] = config[model]['loss']
                self.model_settings[model]['learning_rate'] = config[model]['learning_rate']
                self.model_settings[model]['epochs'] = config[model]['epochs']
                self.model_settings[model]['batch_size'] = config[model]['batch_size']
                self.model_settings[model]['metrics'] = config[model]['metrics']
            
            if 'linear' in model:
                self.model_settings[model]['sgd_loss'] = config[model]['sgd_loss']
                self.model_settings[model]['sgd_penalty'] = config[model]['sgd_penalty']
                self.model_settings[model]['sgd_alpha'] = [float(x) for x in config[model]['sgd_alpha']]
                self.model_settings[model]['sgd_max_iter'] = config[model]['sgd_max_iter']
                self.model_settings[model]['sgd_tol'] = [float(x) for x in config[model]['sgd_tol']]
                self.model_settings[model]['sgd_learning_rate'] = config[model]['sgd_learning_rate']
                self.model_settings[model]['sgd_early_stopping'] = config[model]['sgd_early_stopping']
                self.model_settings[model]['n_iter'] = config[model]['n_iter']
                self.model_settings[model]['cv'] = config[model]['cv']
                self.model_settings[model]['lda_tol'] = [float(x) for x in config[model]['lda_tol']]

            if 'lasso' in model:
                self.model_settings[model]['alpha'] = config[model]['alpha']
                self.model_settings[model]['max_iter'] = config[model]['max_iter']
                self.model_settings[model]['tol'] = float(config[model]['tol'])
                self.model_settings[model]['n_iter'] = config[model]['n_iter']
                self.model_settings[model]['cv'] = config[model]['cv']
                if 'nsub' in model:
                    self.K_lasso = config[model]['K_lasso']
                if 'efp' in model:
                    self.d_lasso = config[model]['d_lasso']

            if model == 'pfn':
                self.model_settings[model]['Phi_sizes'] = tuple(config[model]['Phi_sizes'])
                self.model_settings[model]['F_sizes'] = tuple(config[model]['F_sizes'])
                self.model_settings[model]['epochs'] = config[model]['epochs']
                self.model_settings[model]['batch_size'] = config[model]['batch_size']
                self.model_settings[model]['use_pids'] = config[model]['use_pids']
                self.min_pt = config[model]['min_pt']
                
            if model == 'efn':
                self.model_settings[model]['Phi_sizes'] = tuple(config[model]['Phi_sizes'])
                self.model_settings[model]['F_sizes'] = tuple(config[model]['F_sizes'])
                self.model_settings[model]['epochs'] = config[model]['epochs']
                self.model_settings[model]['batch_size'] = config[model]['batch_size']
                self.model_settings[model]['learning_rate'] = config[model]['learning_rate']

        if self.constituent_subtraction_study:

            self.model_settings['pfn'] = {}
            self.model_settings['pfn']['Phi_sizes'] = tuple(config['pfn']['Phi_sizes'])
            self.model_settings['pfn']['F_sizes'] = tuple(config['pfn']['F_sizes'])
            self.model_settings['pfn']['epochs'] = config['pfn']['epochs']
            self.model_settings['pfn']['batch_size'] = config['pfn']['batch_size']
            self.model_settings['pfn']['use_pids'] = config['pfn']['use_pids']
            self.min_pt = config['pfn']['min_pt']

            self.model_settings['nsub_dnn'] = {}
            self.model_settings['nsub_dnn']['loss'] = config['nsub_dnn']['loss']
            self.model_settings['nsub_dnn']['learning_rate'] = config['nsub_dnn']['learning_rate']
            self.model_settings['nsub_dnn']['epochs'] = config['nsub_dnn']['epochs']
            self.model_settings['nsub_dnn']['batch_size'] = config['nsub_dnn']['batch_size']
            self.model_settings['nsub_dnn']['metrics'] = config['nsub_dnn']['metrics']

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def analyze_pp_aa(self):

        # Loop through combinations of event type, jetR, and R_max
        self.AUC = {}
        for event_type in self.event_types:
            for jetR in self.jetR_list:
                for jet_pt_bin in self.jet_pt_bins:
                    for R_max in self.max_distance_list:

                        # Skip if no models are selected
                        if not self.models:
                            continue
                    
                        # For hard event, skip constituent subtraction variations
                        if event_type=='hard' and not np.isclose(R_max, 0.):
                            continue
                            
                        # If four-vectors are included, R_max=0 is skipped for combined event,
                        # since due to size/time constraints, we skip merging four-vectors for Rmax=0
                        if 'with_four_vectors' in self.filename and 'combined' in event_type and np.isclose(R_max,0):
                            continue
                    
                        # Clear variables
                        self.y = None
                        self.y_train = None
                        self.y_test = None
                        self.X_particles = None
                        self.X_Nsub = None
                        self.pt = None
                        self.delta_pt = None
   
                        # Create output dir
                        self.output_dir_i = os.path.join(self.output_dir, f'{event_type}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}')
                        if not os.path.exists(self.output_dir_i):
                            os.makedirs(self.output_dir_i)
                    
                        # Read input variables
                        key_suffix = f'_{event_type}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                        with h5py.File(os.path.join(self.output_dir, self.filename), 'r') as hf:

                            # First, get the full input arrays
                            self.y_total = hf[f'y{key_suffix}'][:3*self.n_total]
                            X_particles_total = hf[f'X_four_vectors{key_suffix}'][:3*self.n_total]
                            X_Nsub_total = hf[f'X_Nsub{key_suffix}'][:3*self.n_total]

                            # Check whether any training entries are empty
                            [print(f'WARNING: input entry {i} is empty') for i,x in enumerate(X_Nsub_total) if not x.any()]

                            # Determine total number of jets
                            total_jets = int(self.y_total.size)
                            total_jets_AA = int(np.sum(self.y_total))
                            total_jets_pp = total_jets - total_jets_AA 
                            print(f'Total number of jets available: {total_jets_pp} (pp), {total_jets_AA} (AA)')

                            # If there is an imbalance, remove excess jets
                            if total_jets_pp / total_jets_AA > 1.:
                                indices_to_remove = np.where( np.isclose(self.y_total,0) )[0][total_jets_AA:]
                            elif total_jets_pp / total_jets_AA < 1.:
                                indices_to_remove = np.where( np.isclose(self.y_total,1) )[0][total_jets_pp:]
                            y_balanced = np.delete(self.y_total, indices_to_remove)
                            X_particles_balanced = np.delete(X_particles_total, indices_to_remove, axis=0)
                            X_Nsub_balanced = np.delete(X_Nsub_total, indices_to_remove, axis=0)
                            total_jets = int(y_balanced.size)
                            total_jets_AA = int(np.sum(y_balanced))
                            total_jets_pp = total_jets - total_jets_AA 
                            print(f'Total number of jets available after balancing: {total_jets_pp} (pp), {total_jets_AA} (AA)')

                            # Shuffle dataset 
                            idx = np.random.permutation(len(y_balanced))
                            if y_balanced.shape[0] == idx.shape[0]:
                                y_shuffled = y_balanced[idx]
                                X_particles_shuffled = X_particles_balanced[idx]
                                X_Nsub_shuffled = X_Nsub_balanced[idx]
                            else:
                                print(f'MISMATCH of shape: {y_shuffled.shape} vs. {idx.shape}')

                            # Truncate the input arrays to the requested size
                            self.y = y_shuffled[:self.n_total]
                            self.X_particles = X_particles_shuffled[:self.n_total]
                            self.X_Nsub = X_Nsub_shuffled[:self.n_total]
                            print(f'y_shuffled sum: {np.sum(self.y)}')
                            print(f'y_shuffled shape: {self.y.shape}')

                            # If testing against old dataset, load and swap the pp or AA values 
                            # (does not support balancing the number of jets at the moment)
                            if self.old_pp or self.old_AA:
                                
                                # Load old dataset, and get appropriate entries
                                with h5py.File(os.path.join(self.output_dir, self.old_input_file), 'r') as hf_old:
                                    y_old = hf_old[f'y{key_suffix}'][:self.n_total]
                                    sum = np.sum(y_old)
                                    print(f'y_old sum: {sum}')
                                    print(f'y_old shape: {y_old.shape}')

                                    # Filter desired entries
                                    filter_array = np.zeros(y_old.shape)
                                    if self.old_pp:
                                        filter_array = np.add(filter_array, 1-y_old)
                                    if self.old_AA:
                                        filter_array = np.add(filter_array, y_old)
                                    print(f'original y_old: {y_old}')
                                    filter_array = filter_array.astype(dtype=bool)
                                    print(f'filter array: {filter_array}')
                                    X_particles_old = hf_old[f'X_four_vectors{key_suffix}'][:self.n_total][filter_array]
                                    X_Nsub_old = hf_old[f'X_Nsub{key_suffix}'][:self.n_total][filter_array]

                                    # Set labels for these compressed entries
                                    if self.old_pp and self.old_AA:
                                        y_old_compressed = y_old
                                    elif self.old_pp:
                                        y_old_compressed = np.zeros(X_particles_old.shape[0])
                                    elif self.old_AA:
                                        y_old_compressed = np.ones(X_particles_old.shape[0])
                                    sum = np.sum(y_old_compressed)
                                    print(f'y_old_compressed sum: {sum}')
                                    print(f'y_old_compressed shape: {y_old_compressed.shape}')
                                    print(f'y_old_compressed: {y_old_compressed}')

                                    # Remove corresponding entries from new input
                                    filter_array = np.ones(y_unshuffled.shape, dtype=bool)
                                    if self.old_pp:
                                        filter_array = np.add(filter_array, y_unshuffled-1) 
                                    if self.old_AA:
                                        filter_array = np.add(filter_array, -1*y_unshuffled) 
                                    filter_array = filter_array.astype(dtype=bool)
                                    sum = np.sum(filter_array)
                                    print(f'filter_array sum: {sum}')
                                    print(f'filter_array shape: {filter_array.shape}')
                                    print(f'filter_array: {filter_array}')
                                    X_particles_unshuffled = np.concatenate((X_particles_unshuffled[filter_array], 
                                                                             X_particles_old))[:self.n_total]
                                    X_Nsub_unshuffled = np.concatenate((X_Nsub_unshuffled[filter_array],
                                                                        X_Nsub_old))[:self.n_total]
                                    y_unshuffled = np.concatenate((y_unshuffled[filter_array],
                                                                   y_old_compressed))[:self.n_total]
                                    sum = np.sum(y_unshuffled)
                                    print(f'y_new sum: {sum}')
                                    print(f'y_new shape: {y_unshuffled.shape}')
                                    print(f'y_new: {y_unshuffled}')

                                    # Shuffle dataset 
                                    idx = np.random.permutation(len(y_unshuffled))
                                    if y_unshuffled.shape[0] == idx.shape[0]:
                                        self.y = y_unshuffled[idx]
                                        self.X_particles = X_particles_unshuffled[idx]
                                        self.X_Nsub = X_Nsub_unshuffled[idx]
                                    else:
                                        print(f'MISMATCH of shape: {y_unshuffled.shape} vs. {idx.shape}')
                            
                            # Create a second set of four-vectors in which a min-pt cut is applied -- the labels can stay the same
                            if 'pfn' in self.models:
                                self.X_particles_min_pt = filter_four_vectors(np.copy(self.X_particles), min_pt=self.min_pt)

                            # Also get some QA info
                            self.qa_results = {}
                            for qa_observable in self.qa_observables:
                                qa_result = hf[f'{qa_observable}{key_suffix}'][:self.n_total]
                                if qa_result.shape[0] == 0:
                                    continue
                                self.qa_results[qa_observable] = qa_result
                                self.qa_results[qa_observable][np.isnan(self.qa_results[qa_observable])] = -1.
                            delta_pt_result = hf[f'delta_pt{key_suffix}'][:self.n_total]
                            if delta_pt_result.shape[0] != 0:
                                self.delta_pt = delta_pt_result
                                self.delta_pt[np.isnan(self.delta_pt)] = 0.

                        # Define formatted labels for features
                        self.feature_labels = []
                        for i,N in enumerate(self.N_list):
                            beta = self.beta_list[i]
                            self.feature_labels.append(r'$\tau_{}^{{{}}}$'.format(N,beta))

                        # Split into training and test sets
                        # We will split into validatation sets (for tuning hyperparameters) separately for each model
                        X_Nsub_train, X_Nsub_test, self.y_train, self.y_test = sklearn.model_selection.train_test_split(self.X_Nsub, self.y, test_size=self.test_frac)
                        test_jets = int(self.y_test.size)
                        test_jets_AA = int(np.sum(self.y_test))
                        test_jets_pp = test_jets - test_jets_AA
                        print(f'Total number of test jets: {test_jets_pp} (pp), {test_jets_AA} (AA)')

                        # Construct training/test sets for each K
                        self.training_data = {}
                        for K in self.K_list:
                            n = 3*K-4
                            self.training_data[K] = {}
                            self.training_data[K]['X_Nsub_train'] = X_Nsub_train[:,:n]
                            self.training_data[K]['X_Nsub_test'] = X_Nsub_test[:,:n]
                            self.training_data[K]['N_list'] = self.N_list[:n]
                            self.training_data[K]['beta_list'] = self.beta_list[:n]
                            self.training_data[K]['feature_labels'] = self.feature_labels[:n]

                        # Set up dict to store roc curves
                        self.roc_curve_dict = {}
                        self.roc_curve_dict_lasso = {}
                        self.N_terms_lasso = {}
                        self.observable_lasso = {}
                        for model in self.models:
                            self.roc_curve_dict[model] = {}

                        # Plot the input data
                        jet_pt_bin_rounded = [int(pt) for pt in jet_pt_bin]
                        self.plot_QA(event_type, jetR, jet_pt_bin_rounded, R_max)

                        # Plot first few K (before and after scaling)
                        if 'neural_network' in self.models and K < 5:
                            self.plot_nsubjettiness_distributions(K, self.training_data[K]['X_Nsub_train'], self.y_train, self.training_data[K]['feature_labels'], 'before_scaling')
                            self.plot_nsubjettiness_distributions(K, sklearn.preprocessing.scale(self.training_data[K]['X_Nsub_train']), self.y_train, self.training_data[K]['feature_labels'], 'after_scaling')

                        # Compute EFPs
                        if 'efp_dnn' in self.models or 'efp_linear' in self.models or 'efp_lasso' in self.models:

                            print()
                            print(f'Calculating d <= {self.dmax} EFPs for {self.n_total} jets... ')
                            
                            # Specify parameters of EFPs
                            # TODO: check beta dependence !!
                            efpset = energyflow.EFPSet(('d<=', self.dmax), measure=self.efp_measure, beta=self.efp_beta)

                            # Load labels and data, four vectors. Format: (pT,y,phi,m=0). Note: no PID yet which would be 5th entry... check later!
                            # To make sure, don't need any further preprocessing like for EFNs?
                            X_EFP = self.X_particles
                            Y_EFP = self.y #Note not "to_categorical" here... 
                            # Switch here to Jesse's quark/gluon data set.
                            #X_EFP, self.Y_EFP = energyflow.datasets.qg_jets.load(self.n_train + self.n_val + self.n_test)
                            
                            # Convert to list of np.arrays of jets in format (pT,y,phi,mass or PID) -> dim: (# jets, # particles in jets, #4)
                            # and remove zero entries
                            masked_X_EFP = [x[x[:,0] > 0] for x in X_EFP]
                            
                            # Now compute EFPs
                            X_EFP = efpset.batch_compute(masked_X_EFP)

                            # Record which EFPs correspond to which indices
                            # Note: graph images are available here: https://github.com/pkomiske/EnergyFlow/tree/images/graphs
                            self.graphs = efpset.graphs()[1:]
                            for i,efp in enumerate(self.graphs):
                                print(f'  efp {i} -- edges: {efp}')

                            # Preprocess, plot, and store the EFPs for each d
                            self.X_EFP_train = {}
                            self.X_EFP_test = {}
                            self.Y_EFP_train = {}
                            self.Y_EFP_test = {}
                            for d in range(1, self.dmax+1):

                                # Select EFPs with degree <= d
                                X_EFP_d = X_EFP[:,efpset.sel(('d<=', d))]

                                # Remove the 0th EFP (=1)
                                X_EFP_d = X_EFP_d[:,1:]
                                print(f'There are {X_EFP_d.shape[1]} terms for d<={d} (connected + disconnected, and excluding d=0)')

                                # Plot EFPs
                                if d == 2:
                                    self.plot_efp_distributions(d, X_EFP_d, suffix='before_scaling')
                                    self.plot_efp_distributions(d, sklearn.preprocessing.scale(X_EFP_d.astype(np.float128)), suffix='after_scaling')

                                # Do train/val/test split (Note: separate val_set generated in DNN training.)
                                (X_EFP_train_d, X_EFP_val, 
                                 self.X_EFP_test[d], self.Y_EFP_train[d], 
                                 Y_EFP_val, self.Y_EFP_test[d]) = energyflow.utils.data_split(X_EFP_d, Y_EFP, val=self.n_val, test=self.n_test)
                                
                                # Preprocessing: zero mean unit variance
                                self.X_EFP_train[d] = sklearn.preprocessing.scale(X_EFP_train_d.astype(np.float128))


                            print('Done.') 

                            # Plot a few single observables

                            # EFP Lasso for paper -- run with d = 4
                            #observable = '[(0, 1)] + 3.54* [(0, 1), (0, 2)] + 1.72 * [(0, 1), (0, 2), (0, 3), (0, 4)] -3.82 * [(0, 1), (0, 1), (2, 3), (2, 3)]'
                            if self.dmax == 4:
                                observable = rf'$\mathcal{{O}}^{{\mathrm{{ML}}}}_{{\mathrm{{EFP}}}}$ (4 terms)' 
                                ylabel = rf'$\frac{{1}}{{\sigma}} \frac{{d\sigma}}{{ d \mathcal{{O}}^{{\mathrm{{ML}}}}_{{\mathrm{{EFP}}}} }}$'
                                X = self.X_EFP_train[self.dmax][:,0] + 3.54*self.X_EFP_train[self.dmax][:,4] + 1.72*self.X_EFP_train[self.dmax][:,17] - 3.82*self.X_EFP_train[self.dmax][:,23]
                                y = self.Y_EFP_train[self.dmax]
                                self.plot_observable(X, y, xlabel=observable, ylabel=ylabel, filename='EFP_0_4_17_23.pdf')
                            
                            # Nsub for paper
                            #observable = r'$\tau_{10}^{1})^{0.152} (\tau_{11}^{1})^{0.335} (\tau_{14}^{1})^{1.382} (\tau_{14}^{2})^{2.13}$'
                            observable = rf'$\mathcal{{O}}^{{\mathrm{{ML}}}}_{{N-\mathrm{{sub}}}}$ (4 terms)' 
                            ylabel = rf'$\frac{{1}}{{\sigma}} \frac{{d\sigma}}{{ d \mathcal{{O}}^{{\mathrm{{ML}}}}_{{N-\mathrm{{sub}}}} }}$'
                            #print(self.training_data[self.K_max]['feature_labels'])
                            #print(self.training_data[self.K_max]['X_Nsub_train'].shape)
                            X_train = self.training_data[self.K_max]['X_Nsub_train']
                            eps = 0.
                            term1 = np.multiply( np.power(X_train[:,28]+eps, 0.152), np.power(X_train[:,31]+eps, 0.335))
                            term2 = np.multiply( np.power(X_train[:,40]+eps, 1.382), np.power(X_train[:,41]+eps, 2.13))
                            tau = np.multiply(term1, term2)
                            print(np.mean(term1))
                            print(np.mean(term2))
                            print(np.mean(tau))
                            np.log(tau)
                            y = self.y_train
                            self.plot_observable(tau, y, xlabel=observable, ylabel=ylabel, filename='tau_10_11_14_14.pdf', logx=True, logy=True)

                        # Train models
                        self.train_models(event_type, jetR, jet_pt_bin, R_max)
        
        # If enabled, train models with input in cone before and after constituent subtraction
        # (need input from "combined" event type, and labels from "hard" event type, so we do this separately from above loop)
        if self.constituent_subtraction_study:
            self.perform_constituent_subtraction_study()

        # Run plotting script
        print()
        print('Run plotting script...')
        cmd = f'python plot_ppAA.py -c {self.config_file} -o {self.output_dir}'
        subprocess.run(cmd, check=True, shell=True)

    #---------------------------------------------------------------
    # Train models
    #---------------------------------------------------------------
    def train_models(self, event_type, jetR, jet_pt_bin, R_max):

        # Train ML models
        self.key_suffix = f'_{event_type}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
        for model in self.models:
            print()
        
            # Dict to store AUC
            self.AUC[f'{model}{self.key_suffix}'] = []
        
            model_settings = self.model_settings[model]

            # Nsubjettiness  
            for K in self.K_list:
                if model == 'nsub_linear':
                    self.fit_nsub_linear(model, model_settings, K)
                if model == 'nsub_dnn':
                    self.fit_nsub_dnn(model, model_settings, K)
            if model == 'nsub_lasso':
                self.fit_nsubjettiness_lasso(model, model_settings, self.K_lasso)
                
            # EFPs
            for d in range(1, self.dmax+1):
                if model == 'efp_linear':
                    self.fit_efp_linear(model, model_settings, d)
                if model == 'efp_dnn':
                    self.fit_efp_dnn(model, model_settings, d)
            if model == 'efp_lasso':
                self.fit_efp_lasso(model, model_settings, self.d_lasso)

            # Deep sets
            if model == 'pfn':
                self.fit_pfn(model, model_settings, self.y, self.X_particles)

                model_label = 'pfn_min_pt'
                self.AUC[f'{model_label}{self.key_suffix}'] = []
                self.fit_pfn(model_label, model_settings, self.y, self.X_particles_min_pt)
                
            if model == 'efn':
                self.fit_efn(model, model_settings)

        # Plot traditional observables
        for observable in self.qa_observables:
            if 'matched' not in observable:
                self.roc_curve_dict_lasso[observable] = sklearn.metrics.roc_curve(self.y_total[:self.n_total], -self.qa_results[observable])

        # Save ROC curves to file
        if 'nsub_dnn' in self.models or 'efp_dnn' in self.models or 'nsub_linear' in self.models or 'efp_linear' in self.models or 'pfn' in self.models or 'efn' in self.models:
            output_filename = os.path.join(self.output_dir_i, f'ROC{self.key_suffix}.pkl')
            with open(output_filename, 'wb') as f:
                pickle.dump(self.roc_curve_dict, f)
                pickle.dump(self.AUC, f)

        # Separate lasso from others, so that we can re-run it quickly
        if 'nsub_lasso' in self.models or 'efp_lasso' in self.models:
            output_filename = os.path.join(self.output_dir_i, f'ROC{self.key_suffix}_lasso.pkl')
            with open(output_filename, 'wb') as f_lasso:
                pickle.dump(self.roc_curve_dict_lasso, f_lasso)
                pickle.dump(self.N_terms_lasso, f_lasso)
                pickle.dump(self.observable_lasso, f_lasso)

    #---------------------------------------------------------------
    # Fit linear model for Nsubjettiness
    #---------------------------------------------------------------
    def fit_nsub_linear(self, model, model_settings, K):

        X_train = self.training_data[K]['X_Nsub_train']
        X_test = self.training_data[K]['X_Nsub_test']
        y_train = self.y_train
        y_test = self.y_test
        self.fit_linear_model(X_train, y_train, X_test, y_test, model, model_settings, dim_label='K', dim=K, type='LDA_search')

    #---------------------------------------------------------------
    # Fit linear model for EFPs
    #---------------------------------------------------------------
    def fit_efp_linear(self, model, model_settings, d):

        X_train = self.X_EFP_train[d]
        X_test = self.X_EFP_test[d]
        y_train = self.Y_EFP_train[d]
        y_test = self.Y_EFP_test[d]
        self.fit_linear_model(X_train, y_train, X_test, y_test, model, model_settings, dim_label='d', dim=d, type='LDA_search')

    #---------------------------------------------------------------
    # Fit Lasso for Nsubjettiness
    #---------------------------------------------------------------
    def fit_nsubjettiness_lasso(self, model, model_settings, K):

        X_train = self.training_data[K]['X_Nsub_train']
        X_test = self.training_data[K]['X_Nsub_test']
        y_train = self.y_train
        y_test = self.y_test
        self.fit_lasso(X_train, y_train, X_test, y_test, model, model_settings, dim_label='K', dim=K, observable_type='product')

    #---------------------------------------------------------------
    # Fit Lasso for EFPs
    #---------------------------------------------------------------
    def fit_efp_lasso(self, model, model_settings, d):

        X_train = self.X_EFP_train[self.d_lasso]
        X_test = self.X_EFP_test[self.d_lasso]
        y_train = self.Y_EFP_train[self.d_lasso]
        y_test = self.Y_EFP_test[self.d_lasso]
        self.fit_lasso(X_train, y_train, X_test, y_test, model, model_settings, dim_label='d', dim=d, observable_type='sum')

    #---------------------------------------------------------------
    # Fit Dense Neural Network for Nsubjettiness
    #---------------------------------------------------------------
    def fit_nsub_dnn(self, model, model_settings, K):

        # Preprocessing: zero mean unit variance
        X_Nsub_train = sklearn.preprocessing.scale(self.training_data[K]['X_Nsub_train'])
        X_Nsub_test = sklearn.preprocessing.scale(self.training_data[K]['X_Nsub_test'])

        self.fit_dnn(X_Nsub_train, self.y_train, X_Nsub_test, self.y_test, model, model_settings, dim_label='K', dim=K)

    #---------------------------------------------------------------
    # Fit Dense Neural Network for EFPs
    #---------------------------------------------------------------
    def fit_efp_dnn(self, model, model_settings, d):

        X_train = self.X_EFP_train[d]
        X_test = self.X_EFP_test[d]
        y_train = self.Y_EFP_train[d]
        y_test = self.Y_EFP_test[d]
        self.fit_dnn(X_train, y_train, X_test, y_test, model, model_settings, dim_label='d', dim=d)

    #---------------------------------------------------------------
    # Fit ML model -- SGDClassifier or LinearDiscriminant
    #   - SGDClassifier: Linear model (SVM by default, w/o kernel) with SGD training
    #   - For best performance, data should have zero mean and unit variance
    #---------------------------------------------------------------
    def fit_linear_model(self, X_train, y_train, X_test, y_test, model, model_settings, dim_label='', dim=None, type='SGD'):
        print(f'Training {model} ({type}), {dim_label}={dim}...')
        
        if type == 'SGD':
        
            # Define model
            clf = sklearn.linear_model.SGDClassifier(loss=model_settings['sgd_loss'],
                                                        max_iter=model_settings['sgd_max_iter'],
                                                        learning_rate=model_settings['sgd_learning_rate'],
                                                        early_stopping=model_settings['sgd_early_stopping'],
                                                        random_state=self.random_state)

            # Optimize hyperparameters with random search, using cross-validation to determine best set
            # Here we just search over discrete values, although can also easily specify a distribution
            param_distributions = {'penalty': model_settings['sgd_penalty'],
                                'alpha': model_settings['sgd_alpha'],
                                'tol': model_settings['sgd_tol']}

            randomized_search = sklearn.model_selection.RandomizedSearchCV(clf, param_distributions,
                                                                           n_iter=model_settings['n_iter'],
                                                                           cv=model_settings['cv'],
                                                                           random_state=self.random_state)
            search_result = randomized_search.fit(X_train, y_train)
            final_model = search_result.best_estimator_
            result_info = search_result.cv_results_
            print(f'Best params: {search_result.best_params_}')

            # Get predictions for the test set
            #y_predict_train = final_model.predict(X_train)
            #y_predict_test = final_model.predict(X_test)
            
            y_predict_train = sklearn.model_selection.cross_val_predict(clf, X_train, y_train, cv=3, method="decision_function")
            
            # Compare AUC on train set and test set
            AUC_train = sklearn.metrics.roc_auc_score(y_train, y_predict_train)
            print(f'AUC = {AUC_train} (cross-val train set)')
            print()

            # Compute ROC curve: the roc_curve() function expects labels and scores
            self.roc_curve_dict[model][dim] = sklearn.metrics.roc_curve(y_train, y_predict_train)
        
            # Check number of threhsolds used for ROC curve
            # print('thresholds: {}'.format(self.roc_curve_dict[model][K][2]))
            
            # Plot confusion matrix
            #self.plot_confusion_matrix(self.y_train, y_predict_train, f'{model}_K{K}')

        elif type == 'LDA':

            # energyflow implementation
            clf = energyflow.archs.LinearClassifier(linclass_type='lda')
            history = clf.fit(X_train, y_train)
            preds_EFP = clf.predict(X_test)        
            auc_EFP = sklearn.metrics.roc_auc_score(y_test,preds_EFP[:,1])
            print(f'  AUC = {auc_EFP} (test set)')
            self.AUC[f'{model}{self.key_suffix}'].append(auc_EFP)
            self.roc_curve_dict[model][dim] = sklearn.metrics.roc_curve(y_test, preds_EFP[:,1])

        elif type == 'LDA_search':

            # Define model
            clf = sklearn.discriminant_analysis.LinearDiscriminantAnalysis()

            # Optimize hyperparameters
            param_distributions = {'tol': model_settings['lda_tol']}

            randomized_search = sklearn.model_selection.GridSearchCV(clf, param_distributions)
            search_result = randomized_search.fit(X_train, y_train)
            final_model = search_result.best_estimator_
            result_info = search_result.cv_results_
            print(f'Best params: {search_result.best_params_}')

            y_predict_train = sklearn.model_selection.cross_val_predict(clf, X_train, y_train, cv=3, method="decision_function")
            
            # Compare AUC on train set and test set
            AUC_train = sklearn.metrics.roc_auc_score(y_train, y_predict_train)
            print(f'AUC = {AUC_train} (cross-val train set)')
            print()

            # Compute ROC curve: the roc_curve() function expects labels and scores
            self.roc_curve_dict[model][dim] = sklearn.metrics.roc_curve(y_train, y_predict_train)

    #---------------------------------------------------------------
    # Fit Lasso
    #
    #   The parameter alpha multiplies to L1 term
    #   If convergence error: can increase max_iter and/or tol, and/or set normalize=True
    # 
    # observable_type: ['product', 'sum']
    #   - if sum observable, we assume preprocessing is done beforehand
    #   - if product observable, can uncomment preprocessing below after taking log (but off by default)
    #---------------------------------------------------------------
    def fit_lasso(self, X_train, y_train, X_test, y_test, model, model_settings, dim_label='', dim=None, observable_type='product'):
        print()
        print(f'Training {model} ({observable_type} observable), {dim_label}={dim}...')
        
        # First, copy the test training labels, which we will need for ROC curve
        # This is needed because for product observable we don't want to take the log
        y_test_roc = y_test.copy()
        y_train_roc = y_train.copy()

        # If product observable, take the logarithm of the data and labels, such that the product observable 
        # becomes a sum and the exponents in the product observable become the regression weights
        if observable_type == 'product':

            offset = 1.e-4
            X_train = np.log(X_train + offset)
            X_test = np.log(X_test + offset)

            eps = .01
            y_train = np.log(eps + (1. - 2. * eps) * y_train)
            y_test = np.log(eps + (1. - 2. * eps) * y_test)

            # Preprocessing: zero mean unit variance
            #X_train = sklearn.preprocessing.scale(X_train)
            #X_test = sklearn.preprocessing.scale(X_test)

        # Loop through values of regularization parameter
        self.roc_curve_dict_lasso[model] = {}
        self.N_terms_lasso[model] = {}
        self.observable_lasso[model] = {}

        for alpha in model_settings['alpha']:
            self.fit_lasso_single_alpha(alpha, X_train, y_train, X_test, y_test, y_test_roc, y_train_roc, model, model_settings, 
                                        dim_label=dim_label, dim=dim, observable_type=observable_type)

    #---------------------------------------------------------------
    # Fit Lasso for a single alpha value
    #---------------------------------------------------------------
    def fit_lasso_single_alpha(self, alpha, X_train, y_train, X_test, y_test, y_test_roc, y_train_roc, model, model_settings, 
                                dim_label='', dim=None, observable_type='product'):
        print()
        print(f'Fitting lasso regression with alpha = {alpha}')
    
        lasso_clf = sklearn.linear_model.Lasso(alpha=alpha, max_iter=model_settings['max_iter'],
                                                tol=model_settings['tol'])
                                                
        plot_learning_curve = False
        if plot_learning_curve:
            # Split into validation set
            X_train, X_val, y_train, y_val = sklearn.model_selection.train_test_split(X_train, y_train, test_size=0.2)
            train_errors = []
            validation_errors = []
            train_sizes = np.linspace(0, len(X_train)/2, 50)[1:]
            print('Compute Lasso learning curve...')
            for train_size in train_sizes:
                train_size = int(train_size)
                lasso_clf.fit(X_train[:train_size], y_train[:train_size])
                y_predict_train = lasso_clf.predict(X_train[:train_size])
                y_predict_val = lasso_clf.predict(X_val)
                train_errors.append(sklearn.metrics.mean_squared_error(y_predict_train, y_train[:train_size]))
                validation_errors.append(sklearn.metrics.mean_squared_error(y_predict_val, y_val))
        else:
            # Cross-validation
            lasso_clf.fit(X_train, y_train)
            scores = sklearn.model_selection.cross_val_score(lasso_clf, X_train, y_train,
                                                                        scoring='neg_mean_squared_error',
                                                                        cv=model_settings['cv'])
            print(f'cross-validation scores: {scores}')
            y_predict_train = lasso_clf.predict(X_train)
            rmse = sklearn.metrics.mean_squared_error(y_train, y_predict_train)
            print(f'training rmse: {rmse}')

        # Compute AUC on test set
        y_predict_test = lasso_clf.predict(X_test)
        auc_test = sklearn.metrics.roc_auc_score(y_test_roc, y_predict_test)
        rmse_test = sklearn.metrics.mean_squared_error(y_test, y_predict_test)
        print(f'AUC = {auc_test} (test set)')
        print(f'test rmse: {rmse_test}')
        
        # ROC curve
        self.roc_curve_dict_lasso[model][alpha] = sklearn.metrics.roc_curve(y_test_roc, y_predict_test)
        
        if plot_learning_curve:
            plt.axis([0, train_sizes[-1], 0, 10])
            plt.xlabel('training size', fontsize=16)
            plt.ylabel('MSE', fontsize=16)
            plt.plot(train_sizes, train_errors, linewidth=2,
                        linestyle='solid', alpha=0.9, color=sns.xkcd_rgb['dark sky blue'], label='train')
            plt.plot(train_sizes, validation_errors, linewidth=2,
                        linestyle='solid', alpha=0.9, color=sns.xkcd_rgb['watermelon'], label='val')
            plt.axline((0, rmse_test), (len(X_train), rmse_test), linewidth=4, label='test',
                        linestyle='dotted', alpha=0.9, color=sns.xkcd_rgb['medium green'])
            plt.legend(loc='best', fontsize=12)
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir_i, f'Lasso_learning_curve_a{alpha}.pdf'))
            plt.close()

        # Print out observable
        observable = ''
        n_terms = 0
        coeffs = lasso_clf.coef_
        nonzero_coeffs = coeffs[np.absolute(coeffs)>1e-10]
        mean_coeff = np.mean(np.absolute(nonzero_coeffs))
        if observable_type == 'product' and np.mean(nonzero_coeffs) < 0.:
            mean_coeff *= -1
        elif observable_type == 'sum' and np.mean( np.dot(X_test, coeffs) ) < 0:
            mean_coeff *= -1
        print(f'mean_coeff: {mean_coeff}')
        coeffs = np.divide(coeffs, mean_coeff)
        for i,_ in enumerate(coeffs):
            coeff = np.round(coeffs[i], 3)
            if not np.isclose(coeff, 0., atol=1e-10):
                n_terms += 1

                if observable_type == 'product':
                    if 'nsub' in model:
                        N,beta = self.N_beta_from_index(i)
                        observable += rf'(\tau_{{{N}}}^{{{beta}}})^{{{coeff}}} '
                    elif 'efp' in model:
                        observable += rf'({self.graphs[i]})^{{{coeff}}} '

                elif observable_type == 'sum':
                    if 'efp' in model:
                        if n_terms > 0:
                            observable += ' + '
                        observable += f'{coeff} * {self.graphs[i]}'

        print(f'Observable: {observable}')

        self.N_terms_lasso[model][alpha] = n_terms
        self.observable_lasso[model][alpha] = observable

        # Plot observable
        if observable_type == 'product':
            designed_observable = np.exp( np.dot(X_train, coeffs) )
            y = y_train_roc
        elif observable_type == 'sum':
            designed_observable = np.dot(X_test, coeffs) # For EFPs use X_test since is not preprocessed
            y = y_test_roc
        xlabel = rf'$\mathcal{{O}} = {observable}$'
        ylabel = rf'$\frac{{1}}{{\sigma}} \frac{{d\sigma}}{{ d \mathcal{{O}} }}$'
        if 'nsub' in model:
            if n_terms < 5:
                xfontsize=12
            else:
                xfontsize=8
        elif 'efp' in model:
            if n_terms < 2:
                xfontsize=12
            else:
                xfontsize=6
        logy = 'nsub' in model

        self.plot_observable(designed_observable, y, xlabel=xlabel, ylabel=ylabel, filename=f'{model}_{alpha}.pdf', logy=logy)

        if 'nsub' in model:
            # Plot tau_(N-1)_1 directly
            observable = rf'$\tau_{{{19}}}^{{{1}}}$'
            index = 3*self.K_lasso - 6
            tau_max_1 = np.exp(X_train[:,index])
            self.plot_observable(tau_max_1, y, xlabel=observable, filename='tau_max_1.pdf')

            # Plot 1/tau_(N-1)_1
            observable = rf'$1/\tau_{{{19}}}^{{{1}}}$'
            index = 3*self.K_lasso - 6
            tau_max_1_inv = np.reciprocal(np.exp(X_train[:,index]) + 0.0001)
            self.plot_observable(tau_max_1_inv, y, xlabel=observable, filename='tau_max_1_inv.pdf')

    #---------------------------------------------------------------
    # Train DNN, using hyperparameter optimization with keras tuner
    #---------------------------------------------------------------
    def fit_dnn(self, X_train, Y_train, X_test, Y_test, model, model_settings, dim_label='', dim=None):
        print()
        print(f'Training {model}, {dim_label}={dim}...')

        tuner = keras_tuner.Hyperband(functools.partial(self.dnn_builder, input_shape=[X_train.shape[1]], model_settings=model_settings),
                                        objective='val_accuracy',
                                        max_epochs=10,
                                        factor=3,
                                        directory='keras_tuner',
                                        project_name=f'{model}{dim}')

        tuner.search(X_train, Y_train, 
                        batch_size=model_settings['batch_size'],
                        epochs=model_settings['epochs'], 
                        validation_split=self.val_frac)
        
        # Get the optimal hyperparameters
        best_hps=tuner.get_best_hyperparameters(num_trials=1)[0]
        units1 = best_hps.get('units1')
        units2 = best_hps.get('units2')
        units3 = best_hps.get('units3')
        learning_rate = best_hps.get('learning_rate')
        print()
        print(f'Best hyperparameters:')
        print(f'   units: ({units1}, {units2}, {units3})')
        print(f'   learning_rate: {learning_rate}')
        print()

        # Retrain the model with best number of epochs
        hypermodel = tuner.hypermodel.build(best_hps)
        history = hypermodel.fit(X_train, Y_train, epochs=model_settings['epochs'], validation_split=self.val_frac)

        # Plot metrics as a function of epochs
        self.plot_NN_epochs(model_settings['epochs'], history, model, dim_label=dim_label, dim=dim) 

        # Get predictions for test data set
        preds_DNN = hypermodel.predict(X_test).reshape(-1)
        
        # Get AUC
        auc_DNN = sklearn.metrics.roc_auc_score(Y_test, preds_DNN)
        print(f'  AUC = {auc_DNN} (test set)')
        
        # Store AUC
        self.AUC[f'{model}{self.key_suffix}'].append(auc_DNN)
        
        # Get & store ROC curve
        self.roc_curve_dict[model][dim] = sklearn.metrics.roc_curve(Y_test, preds_DNN)

    #---------------------------------------------------------------
    # Construct model for hyperparameter tuning with keras tuner
    #---------------------------------------------------------------
    def dnn_builder(self, hp, input_shape, model_settings):

        model = keras.models.Sequential()
        model.add(keras.layers.Flatten(input_shape=input_shape))

        # Tune size of first dense layer
        hp_units1 = hp.Int('units1', min_value=32, max_value=512, step=32)
        hp_units2 = hp.Int('units2', min_value=32, max_value=512, step=32)
        hp_units3 = hp.Int('units3', min_value=32, max_value=512, step=32)
        model.add(keras.layers.Dense(units=hp_units1, activation='relu'))
        model.add(keras.layers.Dense(units=hp_units2, activation='relu'))
        model.add(keras.layers.Dense(units=hp_units3, activation='relu'))
        model.add(keras.layers.Dense(1,activation='sigmoid'))  # softmax? # Last layer has to be 1 or 2 for binary classification?

        # Print DNN summary
        model.summary()

        # Tune the learning rate for the optimizer
        hp_learning_rate = hp.Choice('learning_rate', values=model_settings['learning_rate']) # if error, change name to lr or learning_rate

        model.compile(optimizer=keras.optimizers.Adam(learning_rate=hp_learning_rate),  # For Stochastic gradient descent use: SGD
                      loss=model_settings['loss'],
                      metrics=model_settings['metrics'])

        return model

    #---------------------------------------------------------------
    # Fit ML model -- Deep Set/Particle Flow Networks
    #---------------------------------------------------------------
    def fit_pfn(self, model, model_settings, y, X_particles):
    
        # Convert labels to categorical
        Y_PFN = energyflow.utils.to_categorical(y, num_classes=2)
                        
        # (pT,y,phi,m=0)
        X_PFN = X_particles
        
        # Preprocess by centering jets and normalizing pts
        for x_PFN in X_PFN:
            mask = x_PFN[:,0] > 0
            yphi_avg = np.average(x_PFN[mask,1:3], weights=x_PFN[mask,0], axis=0)
            x_PFN[mask,1:3] -= yphi_avg
            x_PFN[mask,0] /= x_PFN[:,0].sum()
        
        # Handle particle id channel !! Note: If changed to pT,y,phi,m the 4th component is not PID but m .. fix later
        #if model_settings['use_pids']:
        #    self.my_remap_pids(X_PFN)
        #else:
        X_PFN = X_PFN[:,:,:3]
        
        # Check shape
        if y.shape[0] != X_PFN.shape[0]:
            print(f'Number of labels {y.shape} does not match number of jets {X_PFN.shape} ! ')

        # Split data into train, val and test sets
        (X_PFN_train, X_PFN_val, X_PFN_test, Y_PFN_train, Y_PFN_val, Y_PFN_test) = energyflow.utils.data_split(X_PFN, Y_PFN,
                                                                                             val=self.n_val, test=self.n_test)
        # Build architecture
        pfn = energyflow.archs.PFN(input_dim=X_PFN.shape[-1],
                                   Phi_sizes=model_settings['Phi_sizes'],
                                   F_sizes=model_settings['F_sizes'])

        # Train model
        history = pfn.fit(X_PFN_train,
                          Y_PFN_train,
                          epochs=model_settings['epochs'],
                          batch_size=model_settings['batch_size'],
                          validation_data=(X_PFN_val, Y_PFN_val),
                          verbose=1)
                          
        # Plot metrics are a function of epochs
        self.plot_NN_epochs(model_settings['epochs'], history, model)
        
        # Get predictions on test data
        preds_PFN = pfn.predict(X_PFN_test, batch_size=1000)

        # Get AUC and ROC curve + make plot
        auc_PFN = sklearn.metrics.roc_auc_score(Y_PFN_test[:,1], preds_PFN[:,1])
        print('Particle Flow Networks/Deep Sets: AUC = {} (test set)'.format(auc_PFN))
        self.AUC[f'{model}{self.key_suffix}'].append(auc_PFN)
        
        self.roc_curve_dict[model] = sklearn.metrics.roc_curve(Y_PFN_test[:,1], preds_PFN[:,1])
        
    #---------------------------------------------------------------
    # Fit ML model -- (IRC safe) Energy Flow Networks
    #---------------------------------------------------------------
    def fit_efn(self, model, model_settings):
    
        # Convert labels to categorical
        Y_EFN = energyflow.utils.to_categorical(self.y, num_classes=2)
                        
        # (pT,y,phi,m=0)
        X_EFN = self.X_particles
        
        # Can switch here to quark vs gluon data set
        #X_EFN, y_EFN = energyflow.datasets.qg_jets.load(self.n_train + self.n_val + self.n_test)
        #Y_EFN = energyflow.utils.to_categorical(y_EFN, num_classes=2)
        #print('(n_jets, n_particles per jet, n_variables): {}'.format(X_EFN.shape))

        # For now just use the first 30 entries of the 4-vectors per jet (instead of 800)
        #np.set_printoptions(threshold=sys.maxsize)        
        #X_EFN = X_EFN[:,:30]        
        
        # Preprocess data set by centering jets and normalizing pts
        # Note: this step is somewhat different for pp/AA compared to the quark/gluon data set
        for x_EFN in X_EFN:
            mask = x_EFN[:,0] > 0
            
            # Compute y,phi averages
            yphi_avg = np.average(x_EFN[mask,1:3], weights=x_EFN[mask,0], axis=0)

            # Adjust phi range: Initially it is [0,2Pi], now allow for negative values and >2Pi 
            # so there are no gaps for a given jet.
            # Mask particles that are far away from the average phi & cross the 2Pi<->0 boundary
            mask_phi_1 = ((x_EFN[:,2] - yphi_avg[1] >  np.pi) & (x_EFN[:,2] != 0.))
            mask_phi_2 = ((x_EFN[:,2] - yphi_avg[1] < -np.pi) & (x_EFN[:,2] != 0.))
            
            x_EFN[mask_phi_1,2] -= 2*np.pi
            x_EFN[mask_phi_2,2] += 2*np.pi            
            
            # Now recompute y,phi averages after adjusting the phi range
            yphi_avg1 = np.average(x_EFN[mask,1:3], weights=x_EFN[mask,0], axis=0)            
            
            # And center jets in the y,phi plane
            x_EFN[mask,1:3] -= yphi_avg1

            # Normalize transverse momenta p_Ti -> z_i
            x_EFN[mask,0] /= x_EFN[:,0].sum()
            
            # Set particle four-vectors to zero if the z value is below a certain threshold.
            mask2 = x_EFN[:,0]<0.00001
            x_EFN[mask2,:]=0
        
        # Do not use PID for EFNs
        X_EFN = X_EFN[:,:,:3]
        
        # Make 800 four-vector array smaller, e.g. only 150. Ok w/o background
        X_EFN = X_EFN[:,:150]
        
        # Check shape
        if self.y.shape[0] != X_EFN.shape[0]:
            print(f'Number of labels {self.y.shape} does not match number of jets {X_EFN.shape} ! ')
            
        # Split data into train, val and test sets 
        # and separate momentum fraction z and angles (y,phi)
        (z_EFN_train, z_EFN_val, z_EFN_test, 
         p_EFN_train, p_EFN_val, p_EFN_test,
         Y_EFN_train, Y_EFN_val, Y_EFN_test) = energyflow.utils.data_split(X_EFN[:,:,0], X_EFN[:,:,1:], Y_EFN, 
                                                                           val=self.n_val, test=self.n_test)
        
        # Build architecture
        opt = keras.optimizers.Adam(learning_rate=model_settings['learning_rate']) # if error, change name to learning_rate
        efn = energyflow.archs.EFN(input_dim=2,
                                   Phi_sizes=model_settings['Phi_sizes'],
                                   F_sizes=model_settings['F_sizes'],
                                   optimizer=opt)
        
        # Train model
        history = efn.fit([z_EFN_train,p_EFN_train],
                          Y_EFN_train,
                          epochs=model_settings['epochs'],
                          batch_size=model_settings['batch_size'],
                          validation_data=([z_EFN_val,p_EFN_val], Y_EFN_val),
                          verbose=1)
                          
        # Plot metrics are a function of epochs
        self.plot_NN_epochs(model_settings['epochs'], history, model)
        
        # Get predictions on test data
        preds_EFN = efn.predict([z_EFN_test,p_EFN_test], batch_size=1000)     

        # Get AUC and ROC curve + make plot
        auc_EFN = sklearn.metrics.roc_auc_score(Y_EFN_test[:,1], preds_EFN[:,1])
        print('(IRC safe) Energy Flow Networks: AUC = {} (test set)'.format(auc_EFN))
        self.AUC[f'{model}{self.key_suffix}'].append(auc_EFN)
        
        self.roc_curve_dict[model] = sklearn.metrics.roc_curve(Y_EFN_test[:,1], preds_EFN[:,1])
        
    #---------------------------------------------------------------
    # Perform constituent subtraction study with PFN
    #---------------------------------------------------------------
    def perform_constituent_subtraction_study(self):
        print()
        print('Performing constituent subtraction study...')
        print()        

        # Create output dir
        jetR = self.jetR_list[0]
        jet_pt_bin = self.jet_pt_bins[0]
        R_max = self.max_distance_list[-1]
        self.output_dir_i = os.path.join(self.output_dir, f'constituent_subtraction_R{jetR}_pt{jet_pt_bin}_Rmax0')
        if not os.path.exists(self.output_dir_i):
            os.makedirs(self.output_dir_i)

        # Get the labels from the hard event, and the four-vectors from the combined event
        key_suffix_hard = f'_hard_R{jetR}_pt{jet_pt_bin}_Rmax0'
        key_suffix_combined = f'_combined_matched_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
        self.key_suffix = key_suffix_hard
        with h5py.File(os.path.join(self.output_dir, self.filename), 'r') as hf:

            # First, get the full input arrays
            y_total = hf[f'y{key_suffix_hard}'][:3*self.n_total]
            if self.constituent_subtraction_study['pfn']:
                X_particles_hard_total = hf[f'cone_four_vectors_hard{key_suffix_combined}'][:3*self.n_total]
                X_particles_beforeCS_total = hf[f'cone_four_vectors_beforeCS{key_suffix_combined}'][:3*self.n_total]
                X_particles_afterCS_total = hf[f'cone_four_vectors_afterCS{key_suffix_combined}'][:3*self.n_total]
            if self.constituent_subtraction_study['nsub']:
                X_nsub_hard_total = hf[f'X_Nsub_cone_hard{key_suffix_combined}'][:3*self.n_total]
                X_nsub_beforeCS_total = hf[f'X_Nsub_cone_beforeCS{key_suffix_combined}'][:3*self.n_total]
                X_nsub_afterCS_total = hf[f'X_Nsub_cone_afterCS{key_suffix_combined}'][:3*self.n_total]

            # Determine total number of jets
            total_jets = int(y_total.size)
            total_jets_AA = int(np.sum(y_total))
            total_jets_pp = total_jets - total_jets_AA 
            print(f'Total number of jets available: {total_jets_pp} (pp), {total_jets_AA} (AA)')

            # If there is an imbalance, remove excess jets
            if total_jets_pp / total_jets_AA > 1.:
                indices_to_remove = np.where( np.isclose(y_total,0) )[0][total_jets_AA:]
            elif total_jets_pp / total_jets_AA < 1.:
                indices_to_remove = np.where( np.isclose(y_total,1) )[0][total_jets_pp:]
            y_balanced = np.delete(y_total, indices_to_remove)
            if self.constituent_subtraction_study['pfn']:
                X_particles_hard_balanced = np.delete(X_particles_hard_total, indices_to_remove, axis=0)
                X_particles_beforeCS_balanced = np.delete(X_particles_beforeCS_total, indices_to_remove, axis=0)
                X_particles_afterCS_balanced = np.delete(X_particles_afterCS_total, indices_to_remove, axis=0)
            if self.constituent_subtraction_study['nsub']:
                X_nsub_hard_balanced = np.delete(X_nsub_hard_total, indices_to_remove, axis=0)
                X_nsub_beforeCS_balanced = np.delete(X_nsub_beforeCS_total, indices_to_remove, axis=0)
                X_nsub_afterCS_balanced = np.delete(X_nsub_afterCS_total, indices_to_remove, axis=0)

            total_jets = int(y_balanced.size)
            total_jets_AA = int(np.sum(y_balanced))
            total_jets_pp = total_jets - total_jets_AA 
            print(f'Total number of jets available after balancing: {total_jets_pp} (pp), {total_jets_AA} (AA)')

            # Shuffle dataset 
            idx = np.random.permutation(len(y_balanced))
            if y_balanced.shape[0] == idx.shape[0]:
                y_shuffled = y_balanced[idx]
                if self.constituent_subtraction_study['pfn']:
                    X_particles_hard_shuffled = X_particles_hard_balanced[idx]
                    X_particles_beforeCS_shuffled = X_particles_beforeCS_balanced[idx]
                    X_particles_afterCS_shuffled = X_particles_afterCS_balanced[idx]
                if self.constituent_subtraction_study['nsub']:
                    X_nsub_hard_shuffled = X_nsub_hard_balanced[idx]
                    X_nsub_beforeCS_shuffled = X_nsub_beforeCS_balanced[idx]
                    X_nsub_afterCS_shuffled = X_nsub_afterCS_balanced[idx]
            else:
                print(f'MISMATCH of shape: {y_balanced.shape} vs. {idx.shape}')

            # Truncate the input arrays to the requested size
            y = y_shuffled[:self.n_total]
            if self.constituent_subtraction_study['pfn']:
                X_particles_hard = X_particles_hard_shuffled[:self.n_total]
                X_particles_beforeCS = X_particles_beforeCS_shuffled[:self.n_total]
                X_particles_afterCS = X_particles_afterCS_shuffled[:self.n_total]
            if self.constituent_subtraction_study['nsub']:
                X_nsub_hard = X_nsub_hard_shuffled[:self.n_total]
                X_nsub_beforeCS = X_nsub_beforeCS_shuffled[:self.n_total]
                X_nsub_afterCS = X_nsub_afterCS_shuffled[:self.n_total]
            print(f'y_shuffled sum: {np.sum(y)}')
            print(f'y_shuffled shape: {y.shape}')

            # Create a second set of four-vectors in which a min-pt cut is applied -- the labels can stay the same
            if self.constituent_subtraction_study['pfn']:
                X_particles_hard_min_pt = filter_four_vectors(np.copy(X_particles_hard), min_pt=self.min_pt)
                X_particles_beforeCS_min_pt = filter_four_vectors(np.copy(X_particles_beforeCS), min_pt=self.min_pt)
                X_particles_afterCS_min_pt = filter_four_vectors(np.copy(X_particles_afterCS), min_pt=self.min_pt)

            # Split N-subjettiness input into training and test sets
            if self.constituent_subtraction_study['nsub']:
                X_nsub_hard_train, X_nsub_hard_test, y_train, y_test = sklearn.model_selection.train_test_split(X_nsub_hard, y, test_size=self.test_frac)
                X_nsub_beforeCS_train, X_nsub_beforeCS_test, _,_ = sklearn.model_selection.train_test_split(X_nsub_beforeCS, y, test_size=self.test_frac)
                X_nsub_afterCS_train, X_nsub_afterCS_test, _,_ = sklearn.model_selection.train_test_split(X_nsub_afterCS, y, test_size=self.test_frac)

        # Plot the number of nonzero particles for each of the 3 sets of four-vectors
        if self.constituent_subtraction_study['pfn']:
            self.plot_constituent_subtraction_multiplicity(X_particles_hard, X_particles_beforeCS, X_particles_afterCS, 'pfn')
            self.plot_constituent_subtraction_multiplicity(X_particles_hard_min_pt, X_particles_beforeCS_min_pt, X_particles_afterCS_min_pt, 'pfn_min_pt')

        # Plot first few K    
        if self.constituent_subtraction_study['nsub']:

            # Define formatted labels for features
            feature_labels = []
            for i,N in enumerate(self.N_list):
                if i < 5:
                    beta = self.beta_list[i]
                    feature_labels.append(r'$\tau_{}^{{{}}}$'.format(N,beta))

            self.plot_nsubjettiness_distributions(3, X_nsub_hard_train[:,:5], y_train, feature_labels, 'hard')
            self.plot_nsubjettiness_distributions(3, X_nsub_beforeCS_train[:,:5], y_train, feature_labels, 'beforeCS')
            self.plot_nsubjettiness_distributions(3, X_nsub_afterCS_train[:,:5], y_train, feature_labels, 'afterCS')

        # ------------------------------------------

        # Set up dict to store roc curves
        self.roc_curve_dict = {}
        if self.constituent_subtraction_study['pfn']:
            models = ['pfn_hard', 'pfn_beforeCS', 'pfn_afterCS']
            for model in models:
                self.roc_curve_dict[model] = {}
                self.AUC[f'{model}{self.key_suffix}'] = []
        if self.constituent_subtraction_study['pfn_min_pt']:
            models = ['pfn_hard_min_pt', 'pfn_beforeCS_min_pt', 'pfn_afterCS_min_pt']
            for model in models:
                self.roc_curve_dict[model] = {}
                self.AUC[f'{model}{self.key_suffix}'] = []
        if self.constituent_subtraction_study['nsub']:
            models = ['nsub_hard', 'nsub_beforeCS', 'nsub_afterCS']
            for model in models:
                self.roc_curve_dict[model] = {}
                self.AUC[f'{model}{self.key_suffix}'] = []

        # Train models
        if self.constituent_subtraction_study['pfn']:
            model_settings = self.model_settings['pfn']
            self.fit_pfn('pfn_hard', model_settings, y, X_particles_hard)
            self.fit_pfn('pfn_beforeCS', model_settings, y, X_particles_beforeCS)
            self.fit_pfn('pfn_afterCS', model_settings, y, X_particles_afterCS)

        if self.constituent_subtraction_study['pfn_min_pt']:
            model_settings = self.model_settings['pfn']
            self.fit_pfn('pfn_hard_min_pt', model_settings, y, X_particles_hard_min_pt)
            self.fit_pfn('pfn_beforeCS_min_pt', model_settings, y, X_particles_beforeCS_min_pt)
            self.fit_pfn('pfn_afterCS_min_pt', model_settings, y, X_particles_afterCS_min_pt)

        if self.constituent_subtraction_study['nsub']:
            model_settings = self.model_settings['nsub_dnn']
            self.fit_dnn(X_nsub_hard_train, y_train, X_nsub_hard_test, y_test, 'nsub_hard', model_settings, dim_label='K', dim=self.K_max)
            self.fit_dnn(X_nsub_beforeCS_train, y_train, X_nsub_beforeCS_test, y_test, 'nsub_beforeCS', model_settings, dim_label='K', dim=self.K_max)
            self.fit_dnn(X_nsub_afterCS_train, y_train, X_nsub_afterCS_test, y_test, 'nsub_afterCS', model_settings, dim_label='K', dim=self.K_max)

        # Save ROC curves to file
        output_filename = os.path.join(self.output_dir_i, f'ROC_constituent_subtraction.pkl')
        with open(output_filename, 'wb') as f:
            pickle.dump(self.roc_curve_dict, f)

    #---------------------------------------------------------------
    # Plot multiplicity before and after constituent subtraction
    #---------------------------------------------------------------
    def plot_constituent_subtraction_multiplicity(self, X_particles_hard, X_particles_beforeCS, X_particles_afterCS, label):

        n_particles_hard = self.jet_multiplicity(X_particles_hard)
        n_particles_beforeCS = self.jet_multiplicity(X_particles_beforeCS)
        n_particles_afterCS = self.jet_multiplicity(X_particles_afterCS)

        max = np.amax(n_particles_beforeCS)*1.2
        bins = np.linspace(0, max, 50)
        plt.hist(n_particles_hard,
                 bins,
                 histtype='step',
                 density=True,
                 label = 'Jet',
                 linewidth=2,
                 linestyle='-',
                 alpha=0.5)
        plt.hist(n_particles_beforeCS,
                 bins,
                 histtype='step',
                 density=True,
                 label = 'Jet + Background (before subtraction)',
                 linewidth=2,
                 linestyle='-',
                 alpha=0.5)
        plt.hist(n_particles_afterCS,
                 bins,
                 histtype='step',
                 density=True,
                 label = 'Jet + Background (after subtraction)',
                 linewidth=2,
                 linestyle='-',
                 alpha=0.5)
        xtitle = rf'$N_{{\mathrm{{constituents}}}}$'
        plt.xlabel(xtitle, fontsize=14)
        plt.ylabel(rf'$\frac{{ dN }}{{ d{{N_{{ \mathrm{{constituents}} }} }} }}$', fontsize=16)
        legend = plt.legend(loc='best', fontsize=10, frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir_i, f'constituent_subtraction_multiplicity_{label}.pdf'))
        plt.close()

    #---------------------------------------------------------------
    # Get multiplicity of 3d array of four-vectors
    #---------------------------------------------------------------
    def jet_multiplicity(self, X_particles): 

        n_particles_hard = []
        for jet in X_particles[:int(X_particles.shape[0]/100.)]:
            n_particles = 0
            for four_vector in jet:
                if not np.isclose(np.sum(np.absolute(four_vector)), 0.):
                    n_particles += 1
            n_particles_hard.append(n_particles)
        return n_particles_hard

    #---------------------------------------------------------------
    # Return N,beta from N-subjettiness index
    #---------------------------------------------------------------
    def N_beta_from_index(self, i):
               
        N = 1 + int((i+1)/3)
        if i%3 == 0:
            beta = 1
        elif i%3 == 1:
            beta = 2
        elif i%3 == 2:
            beta = 0.5
        
        return N,beta

    #--------------------------------------------------------------- 
    # My own remap PID routine (similar to remap_pids from energyflow)
    #---------------------------------------------------------------         
    def my_remap_pids(self,events, pid_i=3, error_on_unknown=True):
        # PDGid to small float dictionary
        PID2FLOAT_MAP = {0: 0.0, 22: 1.4,
                         211: .1, -211: .2,
                         321: .3, -321: .4,
                         130: .5,
                         2112: .6, -2112: .7,
                         2212: .8, -2212: .9,
                         11: 1.0, -11: 1.1,
                         13: 1.2, -13: 1.3}
        
        """Remaps PDG id numbers to small floats for use in a neural network.
        `events` are modified in place and nothing is returned.
    
        **Arguments**
    
        - **events** : _numpy.ndarray_
            - The events as an array of arrays of particles.
        - **pid_i** : _int_
            - The column index corresponding to pid information in an event.
        - **error_on_unknown** : _bool_
            - Controls whether a `KeyError` is raised if an unknown PDG ID is
            encountered. If `False`, unknown PDG IDs will map to zero.
        """
    
        if events.ndim == 3:
            pids = events[:,:,pid_i].astype(int).reshape((events.shape[0]*events.shape[1]))
            if error_on_unknown:
                events[:,:,pid_i] = np.asarray([PID2FLOAT_MAP[pid]
                                                for pid in pids]).reshape(events.shape[:2])
            else:
                events[:,:,pid_i] = np.asarray([PID2FLOAT_MAP.get(pid, 0)
                                                for pid in pids]).reshape(events.shape[:2])
        else:
            if error_on_unknown:
                for event in events:
                    event[:,pid_i] = np.asarray([PID2FLOAT_MAP[pid]
                                                 for pid in event[:,pid_i].astype(int)])
            else:
                for event in events:
                    event[:,pid_i] = np.asarray([PID2FLOAT_MAP.get(pid, 0)
                                                 for pid in event[:,pid_i].astype(int)])        

    #---------------------------------------------------------------
    # Plot NN metrics are a function of epochs
    #---------------------------------------------------------------
    def plot_NN_epochs(self, n_epochs, history, label, dim_label='', dim=None):
    
        epoch_list = range(1, n_epochs+1)
        loss = history.history['loss']
        val_loss = history.history['val_loss']

        if 'acc' in history.history:
            acc = history.history['acc']
            val_acc = history.history['val_acc']
        else:
            acc = history.history['accuracy']
            val_acc = history.history['val_accuracy']
        
        plt.axis([0, n_epochs, 0, 1])
        plt.xlabel('epochs', fontsize=16)
        plt.plot(epoch_list, loss, linewidth=2,
                 linestyle='solid', alpha=0.9, color=sns.xkcd_rgb['dark sky blue'],
                 label='loss')
        plt.plot(epoch_list, val_loss, linewidth=2,
                 linestyle='solid', alpha=0.9, color=sns.xkcd_rgb['faded purple'],
                 label='val_loss')
        plt.plot(epoch_list, acc, linewidth=2,
                 linestyle='dotted', alpha=0.9, color=sns.xkcd_rgb['watermelon'],
                 label='acc')
        plt.plot(epoch_list, val_acc, linewidth=2,
                 linestyle='dotted', alpha=0.9, color=sns.xkcd_rgb['medium green'],
                 label='val_acc')
        
        plt.legend(loc='best', fontsize=12)
        plt.tight_layout()
        if dim:
            plt.savefig(os.path.join(self.output_dir_i, f'DNN_epoch_{label}_{dim_label}{dim}.pdf'))
        else:
            plt.savefig(os.path.join(self.output_dir_i, f'PFN_epoch_{label}.pdf'))
        plt.close()

    #---------------------------------------------------------------
    # Plot confusion matrix
    # Note: not normalized to relative error
    #---------------------------------------------------------------
    def plot_confusion_matrix(self, y_train, y_predict_train, label):
    
        confusion_matrix = sklearn.metrics.confusion_matrix(y_train, y_predict_train)
        sns.heatmap(confusion_matrix)
        plt.savefig(os.path.join(self.output_dir_i, f'confusion_matrix_{label}.pdf'))
        plt.close()
        
    #---------------------------------------------------------------
    # Plot QA
    #---------------------------------------------------------------
    def plot_QA(self, event_type, jetR, jet_pt_bin, R_max):
    
        for qa_observable in self.qa_observables:
            if qa_observable not in self.qa_results:
                continue
            qa_observable_shape = self.qa_results[qa_observable].shape
            if qa_observable_shape[0] == 0:
                continue
            
            if  self.y_total[:self.n_total].shape[0] != qa_observable_shape[0]:
                sys.exit(f'ERROR: {qa_observable}: {qa_observable_shape}, y shape: {self.y_total[self.n_total].shape}')
               
            jewel_indices = self.y_total[:self.n_total]
            pythia_indices = 1 - self.y_total[:self.n_total]
            result_jewel = self.qa_results[qa_observable][jewel_indices.astype(bool)]
            result_pythia = self.qa_results[qa_observable][pythia_indices.astype(bool)]

            # Set some labels
            if qa_observable == 'jet_theta_g':
                xlabel = rf'$\theta_{{\rm{{g}} }}$'
                ylabel = rf'$\frac{{1}}{{\sigma}} \frac{{d\sigma}}{{d\theta_{{\rm{{g}} }} }}$'
                bins = np.linspace(0, 1., 50)
            elif qa_observable == 'thrust':
                xlabel = rf'$\lambda_2$'
                ylabel = rf'$\frac{{1}}{{\sigma}} \frac{{d\sigma}}{{d\lambda_2 }} }}$'
                bins = np.linspace(0, 0.3, 50)
            elif qa_observable == 'zg':
                xlabel = rf'$z_{{\rm{{g}} }}$'
                ylabel = rf'$\frac{{1}}{{\sigma}} \frac{{d\sigma}}{{dz_{{\rm{{g}} }} }}$'
                bins = np.linspace(0.2, 0.5, 15)
            else:
                ylabel = ''
                xlabel = rf'{qa_observable}'
                bins = np.linspace(0, np.amax(result_pythia), 20)
            plt.xlabel(xlabel, fontsize=14)
            plt.ylabel(ylabel, fontsize=16)

            if qa_observable == 'jet_pt':
                stat='count'
            else:
                stat='density'

            # Construct dataframes for histplot
            df_jewel = pd.DataFrame(result_jewel, columns=[xlabel])
            df_pythia = pd.DataFrame(result_pythia, columns=[xlabel])
            
            # Add label columns to each df to differentiate them for plotting
            df_jewel['generator'] = np.repeat(self.AA_label, result_jewel.shape[0])
            df_pythia['generator'] = np.repeat(self.pp_label, result_pythia.shape[0])
            df = df_jewel.append(df_pythia, ignore_index=True)

            # Histplot
            h = sns.histplot(df, x=xlabel, hue='generator', stat=stat, bins=bins, element='step', common_norm=False)
            h.legend_.set_title(None)
            plt.setp(h.get_legend().get_texts(), fontsize='14') # for legend text

            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir_i, f'{qa_observable}.pdf'))
            plt.close()

        # delta pt
        if event_type == 'combined_matched':
            bins = np.linspace(-50, 50, 100)
            mean = np.round(np.mean(self.delta_pt),2)
            delta_pt_centered = np.subtract(self.delta_pt, mean)
            sigma = np.round(np.std(delta_pt_centered),2)
            plt.hist(delta_pt_centered,
                     bins,
                     histtype='stepfilled',
                     label = rf'$\mathrm{{mean}} = {mean},\;\sigma = {sigma}$',
                     linewidth=2,
                     linestyle='-',
                     alpha=0.5)
            plt.title(rf'{event_type} event: $R = {jetR}, p_T = {jet_pt_bin}, R_{{max}} = {R_max}$', fontsize=14)
            plt.xlabel(r'$\delta p_{T}$', fontsize=14)
            plt.yscale('log')
            legend = plt.legend(loc='best', fontsize=14, frameon=False)
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir_i, f'delta_pt.pdf'))
            plt.close()
            
        # delta pt -- random cone
        delta_pt = False
        if delta_pt:
            if 'combined' in event_type:
                bins = np.linspace(-50, 50, 100)
                mean = np.round(np.mean(self.delta_pt_random_cone),2)
                sigma = np.round(np.std(self.delta_pt_random_cone),2)
                plt.hist(self.delta_pt_random_cone,
                        bins,
                        histtype='stepfilled',
                        label = rf'$\mathrm{{mean}} = {mean},\;\sigma = {sigma}$',
                        linewidth=2,
                        linestyle='-',
                        alpha=0.5)
                plt.title(rf'{event_type} event: $R = {jetR}, p_T = {jet_pt_bin}, R_{{max}} = {R_max}$', fontsize=14)
                plt.xlabel(r'$\delta p_{T}$', fontsize=14)
                plt.yscale('log')
                legend = plt.legend(loc='best', fontsize=14, frameon=False)
                plt.tight_layout()
                plt.savefig(os.path.join(self.output_dir_i, f'delta_pt_random_cone.pdf'))
                plt.close()

    #---------------------------------------------------------------
    # Plot JEWEL vs. PYTHIA
    #---------------------------------------------------------------
    def plot_observable(self, X, y_train, xlabel='', ylabel='', filename='', xfontsize=12, yfontsize=16, logx=False, logy=False):
            
        jewel_indices = y_train
        pythia_indices = 1 - y_train

        observable_jewel = X[jewel_indices.astype(bool)]
        observable_pythia = X[pythia_indices.astype(bool)]

        df_jewel = pd.DataFrame(observable_jewel, columns=[xlabel])
        df_pythia = pd.DataFrame(observable_pythia, columns=[xlabel])

        df_jewel['generator'] = np.repeat(self.AA_label, observable_jewel.shape[0])
        df_pythia['generator'] = np.repeat(self.pp_label, observable_pythia.shape[0])
        df = df_jewel.append(df_pythia, ignore_index=True)

        if filename == 'tau_10_11_14_14.pdf':
            #bins = 10 ** np.linspace(np.log10(1.e-16), np.log10(1.e-10), 50)
            bins = np.linspace(0., 1.e-12, 50)
            print(bins)
            print(df.describe())
        else:
            bins = np.linspace(np.amin(X), np.amax(X), 50)
        h = sns.histplot(df, x=xlabel, hue='generator', stat='density', bins=bins, element='step', common_norm=False, log_scale=[False,logy])
        if h.legend_:
            #h.legend_.set_bbox_to_anchor((0.85, 0.85))
            h.legend_.set_title(None)
            plt.setp(h.get_legend().get_texts(), fontsize='14') # for legend text

        plt.xlabel(xlabel, fontsize=xfontsize)
        plt.ylabel(ylabel, fontsize=yfontsize)

        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir_i, f'{filename}'))
        plt.close()

    #---------------------------------------------------------------
    # Plot nsubjettiness
    #---------------------------------------------------------------
    def plot_nsubjettiness_distributions(self, K, X_train, y_train, feature_labels, suffix=''):
            
        print(f'Plotting input nsubjettiness data, K={K}...')
        print()
        
        # Separate PYTHIA/JEWEL
        jewel_indices = y_train
        pythia_indices = 1 - y_train
        n_plot = int(self.n_train) # Plot a subset to save time/space
        X_Nsub_jewel = X_train[jewel_indices.astype(bool)][:n_plot]
        X_Nsub_pythia = X_train[pythia_indices.astype(bool)][:n_plot]

        # Construct dataframes for scatter matrix plotting
        df_jewel = pd.DataFrame(X_Nsub_jewel, columns=feature_labels)
        df_pythia = pd.DataFrame(X_Nsub_pythia, columns=feature_labels)
        
        # Add label columns to each df to differentiate them for plotting
        df_jewel['generator'] = np.repeat(self.AA_label, X_Nsub_jewel.shape[0])
        df_pythia['generator'] = np.repeat(self.pp_label, X_Nsub_pythia.shape[0])
        df = df_jewel.append(df_pythia, ignore_index=True)

        # Plot scatter matrix
        g = sns.pairplot(df, corner=True, hue='generator', plot_kws={'alpha':0.1})
        #g.legend.set_bbox_to_anchor((0.75, 0.75))
        #plt.savefig(os.path.join(self.output_dir_i, f'training_data_scatter_matrix_K{K}.png'), dpi=50)
        plt.savefig(os.path.join(self.output_dir_i, f'training_data_scatter_matrix_K{K}_{suffix}.pdf'))
        plt.close()
        
        # Plot correlation matrix
        df.drop(columns=['generator'])
        corr_matrix = df.corr()
        sns.heatmap(corr_matrix)
        plt.savefig(os.path.join(self.output_dir_i, f'training_data_correlation_matrix_K{K}_{suffix}.pdf'))
        plt.close()

            
    #---------------------------------------------------------------
    # Plot EFPs
    #---------------------------------------------------------------
    def plot_efp_distributions(self, d, X_EFP_d, suffix='before_scaling'):
            
        print(f'Plotting input EFP data {suffix}, d={d}...')

        # Separate PYTHIA/JEWEL
        jewel_indices = self.y
        pythia_indices = 1 - self.y
        X_jewel = X_EFP_d[jewel_indices.astype(bool)]
        X_pythia = X_EFP_d[pythia_indices.astype(bool)]

        # Get labels
        graphs = [str(x) for x in self.graphs[:4]]

        # Construct dataframes for scatter matrix plotting
        df_jewel = pd.DataFrame(X_jewel, columns=graphs)
        df_pythia = pd.DataFrame(X_pythia, columns=graphs)
        
        # Add label columns to each df to differentiate them for plotting
        df_jewel['generator'] = np.repeat(self.AA_label, X_jewel.shape[0])
        df_pythia['generator'] = np.repeat(self.pp_label, X_pythia.shape[0])
        df = df_jewel.append(df_pythia, ignore_index=True)

        # Plot scatter matrix
        g = sns.pairplot(df, corner=True, hue='generator', plot_kws={'alpha':0.1})
        #g.legend.set_bbox_to_anchor((0.75, 0.75))
        #plt.savefig(os.path.join(self.output_dir_i, f'training_data_scatter_matrix_K{K}.png'), dpi=50)
        plt.savefig(os.path.join(self.output_dir_i, f'training_data_scatter_matrix_d{d}_{suffix}.pdf'))
        plt.close()
        
        # Plot correlation matrix
        df.drop(columns=['generator'])
        corr_matrix = df.corr()
        sns.heatmap(corr_matrix)
        plt.savefig(os.path.join(self.output_dir_i, f'training_data_correlation_matrix_d{d}_{suffix}.pdf'))
        plt.close()
            
##################################################################
if __name__ == '__main__':

    # Define arguments
    parser = argparse.ArgumentParser(description='Process qg')
    parser.add_argument('-c', '--configFile', action='store',
                        type=str, metavar='configFile',
                        default='../../../config/ml/ppAA.yaml',
                        help='Path of config file for analysis')
    parser.add_argument('-o', '--outputDir', action='store',
                        type=str, metavar='outputDir',
                        default='./TestOutput',
                        help='Output directory for output to be written to')
    parser.add_argument('--old_pp', 
                         help='use old pp aggregated input data', 
                         action='store_true', default=False)
    parser.add_argument('--old_AA',
                         help='use old AA aggregated input data',
                         action='store_true', default=False)

    # Parse the arguments
    args = parser.parse_args()

    print('Configuring...')
    print('configFile: \'{0}\''.format(args.configFile))
    print('ouputDir: \'{0}\"'.format(args.outputDir))

    # If invalid configFile is given, exit
    if not os.path.exists(args.configFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    analysis = AnalyzePPAA(config_file=args.configFile, output_dir=args.outputDir, old_pp=args.old_pp, old_AA=args.old_AA)
    analysis.analyze_pp_aa()
