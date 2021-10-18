#!/usr/bin/env python3

"""
Example class to read quark-gluon dataset
"""

import os
import sys
import argparse
import yaml
import h5py
import re

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

# Base class
from pyjetty.alice_analysis.process.base import common_base

################################################################
class AnalyzePPAA(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, input_file='', config_file='', output_dir='', old_pp=False, old_AA=False, debug_level=0, **kwargs):
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
        
        alpha_list = self.config['lasso']['alpha']
        self.colors = {'PFN': sns.xkcd_rgb['faded purple'],
                       'EFN': sns.xkcd_rgb['faded red'],
                       'DNN (K = 5)': sns.xkcd_rgb['watermelon'],
                       'DNN (K = 10)': sns.xkcd_rgb['light brown'],
                       'DNN (K = 20)': sns.xkcd_rgb['medium green'],
                       'DNN (K = 30)': sns.xkcd_rgb['dark sky blue'],
                       'DNN (K = 40)': sns.xkcd_rgb['faded purple'],
                       'Lasso': sns.xkcd_rgb['watermelon'],
                       rf'Lasso $(\alpha = {alpha_list[0]})$': sns.xkcd_rgb['watermelon'],
                       rf'Lasso $(\alpha = {alpha_list[1]})$': sns.xkcd_rgb['light brown'],
                       rf'Lasso $(\alpha = {alpha_list[2]})$': sns.xkcd_rgb['faded purple'],
                       'jet_mass': sns.xkcd_rgb['faded red'],
                       'jet_angularity': sns.xkcd_rgb['medium green'],
                       'jet_theta_g': sns.xkcd_rgb['medium brown']
                      }
        self.linestyles = {'PFN': 'solid',
                           'EFN': 'solid',
                           'DNN': 'solid',
                           'jet_mass': 'dotted',
                           'jet_angularity': 'dotted',
                           'jet_theta_g': 'dotted',
                           'multiplicity': 'dotted',
                           'Lasso': 'dotted'
                          }
                       
        self.filename = 'nsubjettiness_with_four_vectors_unshuffled.h5'
        with h5py.File(os.path.join(self.output_dir, self.filename), 'r') as hf:
            self.N_list = hf['N_list'][:]
            self.beta_list = hf['beta_list'][:]
            self.delta_pt_random_cone = hf['delta_pt_random_cone'][:]
        
        self.qa_observables = ['matched_pt', 'matched_deltaR', 'jet_pt', 'jet_angularity', 'jet_mass', 'jet_theta_g', 'jet_subjet_z', 'hadron_z', 'multiplicity_0000', 'multiplicity_0150', 'multiplicity_0500', 'multiplicity_1000']

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
        self.K_ROC_list = config['K_ROC']
        
        # Initialize model-specific settings
        self.config = config
        self.models = config['models']
        self.model_settings = {}
        for model in self.models:
            self.model_settings[model] = {}
            
            if model == 'linear':
                self.model_settings[model]['loss'] = config[model]['loss']
                self.model_settings[model]['penalty'] = config[model]['penalty']
                self.model_settings[model]['alpha'] = [float(x) for x in config[model]['alpha']]
                self.model_settings[model]['max_iter'] = config[model]['max_iter']
                self.model_settings[model]['tol'] = [float(x) for x in config[model]['tol']]
                self.model_settings[model]['learning_rate'] = config[model]['learning_rate']
                self.model_settings[model]['early_stopping'] = config[model]['early_stopping']
                self.model_settings[model]['n_iter'] = config[model]['n_iter']
                self.model_settings[model]['cv'] = config[model]['cv']
                self.model_settings[model]['random_state'] = config[model]['random_state']

            if model == 'random_forest':
                self.model_settings[model]['random_state'] = config[model]['random_state']

            if model == 'neural_network':
                self.model_settings[model]['loss'] = config[model]['loss']
                self.model_settings[model]['learning_rate'] = config[model]['learning_rate']
                self.model_settings[model]['epochs'] = config[model]['epochs']
                self.model_settings[model]['batch_size'] = config[model]['batch_size']
                self.model_settings[model]['metrics'] = config[model]['metrics']
                self.model_settings[model]['random_state'] = config[model]['random_state']
            
            if model == 'pfn':
                self.model_settings[model]['Phi_sizes'] = tuple(config[model]['Phi_sizes'])
                self.model_settings[model]['F_sizes'] = tuple(config[model]['F_sizes'])
                self.model_settings[model]['epochs'] = config[model]['epochs']
                self.model_settings[model]['batch_size'] = config[model]['batch_size']
                self.model_settings[model]['use_pids'] = config[model]['use_pids']
                
            if model == 'efn':
                self.model_settings[model]['Phi_sizes'] = tuple(config[model]['Phi_sizes'])
                self.model_settings[model]['F_sizes'] = tuple(config[model]['F_sizes'])
                self.model_settings[model]['epochs'] = config[model]['epochs']
                self.model_settings[model]['batch_size'] = config[model]['batch_size']
                
            if model == 'lasso':
                self.K_lasso = config[model]['K_lasso']
                self.model_settings[model]['alpha'] = config[model]['alpha']
                self.model_settings[model]['max_iter'] = config[model]['max_iter']
                self.model_settings[model]['tol'] = float(config[model]['tol'])
                self.model_settings[model]['n_iter'] = config[model]['n_iter']
                self.model_settings[model]['cv'] = config[model]['cv']
                self.model_settings[model]['random_state'] = config[model]['random_state']

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
                            self.y_total = hf[f'y{key_suffix}'][:400000]
                            X_particles_total = hf[f'X_four_vectors{key_suffix}'][:400000]
                            X_Nsub_total = hf[f'X_Nsub{key_suffix}'][:400000]

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
                        if 'linear' in self.models:
                            self.roc_curve_dict['SGDClassifier'] = {}
                        if 'random_forest' in self.models:
                            self.roc_curve_dict['RandomForest'] = {}
                        if 'neural_network' in self.models:
                            self.roc_curve_dict['DNN'] = {}
                        if 'lasso' in self.models:
                            self.roc_curve_dict['Lasso'] = {}

                        # Plot the input data
                        jet_pt_bin_rounded = [int(pt) for pt in jet_pt_bin]
                        self.plot_QA(event_type, jetR, jet_pt_bin_rounded, R_max)
                        
                        for K in self.K_list:
                            if K <= 4:
                                self.plot_training_data(K)
                           
                        # Train models
                        self.train_models(event_type, jetR, jet_pt_bin, R_max)
          
        # Plot AUC vs. K for hard vs. combined
        for jetR in self.jetR_list:
            for jet_pt_bin in self.jet_pt_bins:
                for R_max in self.max_distance_list:
                
                    if np.isclose(R_max, 0.):
                        continue
                
                    if 'pfn' in self.models and 'neural_network' in self.models:
                        auc_list = {}
                        suffix = f'_hard_R{jetR}_pt{jet_pt_bin}_Rmax0'
                        auc_list[f'neural_network{suffix}'] = self.AUC[f'neural_network{suffix}']
                        auc_list[f'pfn{suffix}'] = self.AUC[f'pfn{suffix}']
                        suffix = f'_combined_matched_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                        auc_list[f'neural_network{suffix}'] = self.AUC[f'neural_network{suffix}']
                        auc_list[f'pfn{suffix}'] = self.AUC[f'pfn{suffix}']
                        outputdir = os.path.join(self.output_dir, f'combined_matched_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}')
                        self.plot_AUC_convergence(outputdir, auc_list, event_type, jetR, jet_pt_bin, R_max, suffix='3')
    
    #---------------------------------------------------------------
    # Train models
    #---------------------------------------------------------------
    def train_models(self, event_type, jetR, jet_pt_bin, R_max):

        # Train ML models
        self.key_suffix = f'_{event_type}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
        for model in self.models:
        
            # Dict to store AUC
            self.AUC[f'{model}{self.key_suffix}'] = []
        
            model_settings = self.model_settings[model]
            
            for K in self.K_list:
            
                if model == 'linear':
                    self.fit_linear_model(K, model_settings)
                if model == 'random_forest':
                    self.fit_random_forest(K, model_settings)
                if model == 'neural_network':
                    self.fit_neural_network(K, model_settings)
                
            if model == 'lasso':
                self.fit_lasso(model_settings)
            
            if model == 'pfn':
                self.fit_pfn(model_settings)
                
            if model == 'efn':
                self.fit_efn(model_settings)
                
        # Plot traditional observables
        self.roc_curve_dict['jet_mass'] = sklearn.metrics.roc_curve(self.y_total[:self.n_total], -self.qa_results['jet_mass'])
        self.roc_curve_dict['jet_angularity'] = sklearn.metrics.roc_curve(self.y_total[:self.n_total], -self.qa_results['jet_angularity'])
        self.roc_curve_dict['jet_theta_g'] = sklearn.metrics.roc_curve(self.y_total[:self.n_total], -self.qa_results['jet_theta_g'])
        #self.roc_curve_dict['multiplicity_0000'] = sklearn.metrics.roc_curve(self.y, -self.qa_results['multiplicity_0000'])

        # Plot several verisons of AUC vs. K curve
        if 'neural_network' in self.models:
            auc_list = {}
            auc_list[f'neural_network{self.key_suffix}'] = self.AUC[f'neural_network{self.key_suffix}']
            self.plot_AUC_convergence(self.output_dir_i, auc_list, event_type, jetR, jet_pt_bin, R_max, suffix='1')
        if 'pfn' in self.models and 'neural_network' in self.models:
            auc_list = {}
            auc_list[f'neural_network{self.key_suffix}'] = self.AUC[f'neural_network{self.key_suffix}']
            auc_list[f'pfn{self.key_suffix}'] = self.AUC[f'pfn{self.key_suffix}']
            self.plot_AUC_convergence(self.output_dir_i, auc_list, event_type, jetR, jet_pt_bin, R_max, suffix='2')
        
        # Plot several versions of ROC curves and significance improvement
        if 'pfn' in self.models and 'neural_network' in self.models and 'efn' in self.models:
            roc_list = {}
            roc_list['PFN'] = self.roc_curve_dict['PFN']
            roc_list['EFN'] = self.roc_curve_dict['EFN']
            for K in self.K_ROC_list:
                roc_list[f'DNN (K = {K})'] = self.roc_curve_dict['DNN'][K]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='1')
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='1')
        
        if 'pfn' in self.models and 'neural_network' in self.models:
            roc_list = {}
            roc_list['PFN'] = self.roc_curve_dict['PFN']
            for K in self.K_ROC_list:
                roc_list[f'DNN (K = {K})'] = self.roc_curve_dict['DNN'][K]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='2')
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='2')
        
        if 'neural_network' in self.models:
            roc_list = {}
            for K in self.K_ROC_list:
                roc_list[f'DNN (K = {K})'] = self.roc_curve_dict['DNN'][K]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='3')
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='3')
        
        if 'lasso' in self.models:
            roc_list = {}
            roc_list['jet_mass'] = self.roc_curve_dict['jet_mass']
            roc_list['jet_angularity'] = self.roc_curve_dict['jet_angularity']
            roc_list['jet_theta_g'] = self.roc_curve_dict['jet_theta_g']
            for alpha in self.config['lasso']['alpha']:
                roc_list[rf'Lasso $(\alpha = {alpha})$'] = self.roc_curve_dict['Lasso'][self.K_lasso][alpha]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='4')
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='4')
        
        if 'neural_network' in self.models and 'lasso' in self.models:
            roc_list = {}
            roc_list[f'DNN (K = {self.K_lasso})'] = self.roc_curve_dict['DNN'][K]
            roc_list['jet_mass'] = self.roc_curve_dict['jet_mass']
            roc_list['jet_angularity'] = self.roc_curve_dict['jet_angularity']
            roc_list['jet_theta_g'] = self.roc_curve_dict['jet_theta_g']
            for alpha in self.config['lasso']['alpha']:
                roc_list[rf'Lasso $(\alpha = {alpha})$'] = self.roc_curve_dict['Lasso'][self.K_lasso][alpha]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='5')
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='5')
        
    #---------------------------------------------------------------
    # Fit ML model -- 1. SGDClassifier
    #   - Linear model (SVM by default, w/o kernel) with SGD training
    #   - For best performance, data should have zero mean and unit variance
    #---------------------------------------------------------------
    def fit_linear_model(self, K, model_settings):
        print(f'Training SGDClassifier, K={K}...')
        
        # Define model
        sgd_clf = sklearn.linear_model.SGDClassifier(loss=model_settings['loss'],
                                                     max_iter=model_settings['max_iter'],
                                                     learning_rate=model_settings['learning_rate'],
                                                     early_stopping=model_settings['early_stopping'],
                                                     random_state=model_settings['random_state'])

        # Optimize hyperparameters with random search, using cross-validation to determine best set
        # Here we just search over discrete values, although can also easily specify a distribution
        param_distributions = {'penalty': model_settings['penalty'],
                               'alpha': model_settings['alpha'],
                               'tol': model_settings['tol']}
        randomized_search = sklearn.model_selection.RandomizedSearchCV(sgd_clf, param_distributions,
                                                                       n_iter=model_settings['n_iter'],
                                                                       cv=model_settings['cv'],
                                                                       random_state=model_settings['random_state'])
        search_result = randomized_search.fit(self.training_data[K]['X_Nsub_train'], self.y_train)
        final_model = search_result.best_estimator_
        result_info = search_result.cv_results_
        print(f'Best params: {search_result.best_params_}')

        # Get predictions for the test set
        #y_predict_train = final_model.predict(self.training_data[K]['X_Nsub_train'])
        #y_predict_test = final_model.predict(self.training_data[K]['X_Nsub_test'])
        
        y_predict_train = sklearn.model_selection.cross_val_predict(sgd_clf, self.training_data[K]['X_Nsub_train'], self.y_train, cv=3,method="decision_function")
        
        # Compare AUC on train set and test set
        AUC_train = sklearn.metrics.roc_auc_score(self.y_train, y_predict_train)
        print(f'AUC = {AUC_train} (train set)')
        AUC_test = sklearn.metrics.roc_auc_score(self.y_train, y_predict_train)
        print(f'AUC = {AUC_test} (test set)')
        print()

        # Compute ROC curve: the roc_curve() function expects labels and scores
        self.roc_curve_dict['SGDClassifier'][K] = sklearn.metrics.roc_curve(self.y_train, y_predict_train)
    
        # Check number of threhsolds used for ROC curve
        # print('thresholds: {}'.format(self.roc_curve_dict['SGDClassifier'][K][2]))
        
        # Plot confusion matrix
        #self.plot_confusion_matrix(self.y_train, y_predict_train, f'linear_K{K}')
 
    #---------------------------------------------------------------
    # Fit ML model -- 2. Random Forest Classifier
    #---------------------------------------------------------------
    def fit_random_forest(self, K, model_settings):
        print(f'Training Random Forest Classifier, K={K}...')
        
        forest_clf = sklearn.ensemble.RandomForestClassifier(random_state=model_settings['random_state'])
        y_Nsub_probas_forest = sklearn.model_selection.cross_val_predict(forest_clf, self.training_data[K]['X_Nsub_train'], self.y_train, cv=3,method="predict_proba")
        
        # The output here are class probabilities. We us use the positive class's probability for the ROC curve
        y_Nsub_scores_forest = y_Nsub_probas_forest[:,1]
        
        print(y_Nsub_scores_forest)
        
        # Compute AUC & ROC curve
        Nsub_auc_RFC = sklearn.metrics.roc_auc_score(self.y_train,y_Nsub_scores_forest)
        print(f'AUC = {Nsub_auc_RFC} (cross validation)')
        print()

        self.roc_curve_dict['RandomForest'][K] = sklearn.metrics.roc_curve(self.y_train,y_Nsub_scores_forest)

    #---------------------------------------------------------------
    # Fit ML model -- 3. Dense Neural network with Keras
    #---------------------------------------------------------------
    def fit_neural_network(self, K, model_settings):
        print(f'Training Dense Neural Network, K={K}...')
     
        # input_shape expects shape of an instance (not including batch size)
        DNN = keras.models.Sequential()
        DNN.add(keras.layers.Flatten(input_shape=[self.training_data[K]['X_Nsub_train'].shape[1]]))
        DNN.add(keras.layers.Dense(300,activation='relu'))
        DNN.add(keras.layers.Dense(300,activation='relu'))
        DNN.add(keras.layers.Dense(100,activation='relu'))
        DNN.add(keras.layers.Dense(1,activation='sigmoid')) # softmax? # Last layer has to be 1 or 2 for binary classification?

        # Print DNN summary
        DNN.summary()
        
        # Compile DNN
        opt = keras.optimizers.Adam(lr=model_settings['learning_rate']) # if error, change name to learning_rate
        DNN.compile(loss=model_settings['loss'],
                    optimizer=opt,                       # For Stochastic gradient descent use: "sgd"
                    metrics=model_settings['metrics'])

        # Train DNN - use validation_split to split into validation set
        history = DNN.fit(self.training_data[K]['X_Nsub_train'],
                          self.y_train,
                          batch_size=model_settings['batch_size'],
                          epochs=model_settings['epochs'],
                          validation_split=self.val_frac)
                          
        # Plot metrics are a function of epochs
        if K in self.K_ROC_list:
            self.plot_NN_epochs(model_settings['epochs'], history, 'DNN', K)
        
        # Get predictions for test data set
        y_Nsub_test_preds_DNN = DNN.predict(self.training_data[K]['X_Nsub_test']).reshape(-1)
        
        # Get AUC
        Nsub_auc_DNN = sklearn.metrics.roc_auc_score(self.y_test, y_Nsub_test_preds_DNN)
        print(f'AUC = {Nsub_auc_DNN} (test set)')
        self.AUC[f'neural_network{self.key_suffix}'].append(Nsub_auc_DNN)
        
        # Get ROC curve results
        self.roc_curve_dict['DNN'][K] = sklearn.metrics.roc_curve(self.y_test, y_Nsub_test_preds_DNN)

    #---------------------------------------------------------------
    # Fit ML model -- 4. Deep Set/Particle Flow Networks
    #---------------------------------------------------------------
    def fit_pfn(self, model_settings):
    
        # Convert labels to categorical
        Y_PFN = energyflow.utils.to_categorical(self.y, num_classes=2)
                        
        # (pT,y,phi,m=0)
        X_PFN = self.X_particles
        
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
        if self.y.shape[0] != X_PFN.shape[0]:
            print(f'Number of labels {self.y.shape} does not match number of jets {X_PFN.shape} ! ')

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
        self.plot_NN_epochs(model_settings['epochs'], history, 'PFN')
        
        # Get predictions on test data
        preds_PFN = pfn.predict(X_PFN_test, batch_size=1000)

        # Get AUC and ROC curve + make plot
        auc_PFN = sklearn.metrics.roc_auc_score(Y_PFN_test[:,1], preds_PFN[:,1])
        print('Particle Flow Networks/Deep Sets: AUC = {} (test set)'.format(auc_PFN))
        self.AUC[f'pfn{self.key_suffix}'].append(auc_PFN)
        
        self.roc_curve_dict['PFN'] = sklearn.metrics.roc_curve(Y_PFN_test[:,1], preds_PFN[:,1])
        
    #---------------------------------------------------------------
    # Fit ML model -- 5. (IRC safe) Energy Flow Networks
    #---------------------------------------------------------------
    def fit_efn(self, model_settings):
    
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
        opt = keras.optimizers.Adam(learning_rate=0.01)
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
        self.plot_NN_epochs(model_settings['epochs'], history, 'EFN')
        
        # Get predictions on test data
        preds_EFN = efn.predict([z_EFN_test,p_EFN_test], batch_size=1000)     

        # Get AUC and ROC curve + make plot
        auc_EFN = sklearn.metrics.roc_auc_score(Y_EFN_test[:,1], preds_EFN[:,1])
        print('(IRC safe) Energy Flow Networks: AUC = {} (test set)'.format(auc_EFN))
        self.AUC[f'efn{self.key_suffix}'].append(auc_EFN)
        
        self.roc_curve_dict['EFN'] = sklearn.metrics.roc_curve(Y_EFN_test[:,1], preds_EFN[:,1])
        
    #---------------------------------------------------------------
    # Fit ML model -- 6. Lasso regression
    #   The parameter alpha multiplies to L1 term
    #   If convergence error: can increase max_iter and/or tol, and/or set normalize=True
    #---------------------------------------------------------------
    def fit_lasso(self, model_settings):
        print(f'Training Lasso regression...')
        
        # Take the logarithm of the data and labels, such that the product observable becomes a sum
        # and the exponents in the product observable become the regression weights
        X_train_lasso = np.log(self.training_data[self.K_lasso]['X_Nsub_train']+0.0001)
        X_test_lasso = np.log(self.training_data[self.K_lasso]['X_Nsub_test']+0.0001)
        
        # Not taking log of test labels to make ROC curve later which requires integers for the test data
        eps = .01
        y_train_lasso = np.log(eps + (1. - 2. * eps) * self.y_train)
        y_test_lasso = np.log(eps + (1. - 2. * eps) * self.y_test)
        y_test_lasso_roc = self.y_test
        
        # Loop through values of regularization parameter
        self.roc_curve_dict['Lasso'][self.K_lasso] = {}
        self.N_terms_lasso = {}
        for alpha in model_settings['alpha']:
            print(f'Fitting lasso regression with alpha = {alpha}')
        
            lasso_clf = sklearn.linear_model.Lasso(alpha=alpha, max_iter=model_settings['max_iter'],
                                                   tol=model_settings['tol'])
                                                   
            plot_learning_curve = False
            if plot_learning_curve:
                # Split into validation set
                X_train, X_val, y_train, y_val = sklearn.model_selection.train_test_split(X_train_lasso, y_train_lasso, test_size=0.2)
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
                lasso_clf.fit(X_train_lasso, y_train_lasso)
                scores = sklearn.model_selection.cross_val_score(lasso_clf, X_train_lasso, y_train_lasso,
                                                                            scoring='neg_mean_squared_error',
                                                                            cv=model_settings['cv'])
                print(f'cross-validation scores: {scores}')
                y_predict_train = lasso_clf.predict(X_train_lasso)
                rmse = sklearn.metrics.mean_squared_error(y_train_lasso, y_predict_train)
                print(f'training rmse: {rmse}')
            
            # Compute AUC on test set
            y_predict_test = lasso_clf.predict(X_test_lasso)
            auc_lasso_test = sklearn.metrics.roc_auc_score(y_test_lasso_roc, y_predict_test)
            rmse_lasso_test = sklearn.metrics.mean_squared_error(y_test_lasso, y_predict_test)
            print(f'AUC = {auc_lasso_test} (test set)')
            print(f'test rmse: {rmse_lasso_test}')
            
            # ROC curve
            self.roc_curve_dict['Lasso'][self.K_lasso][alpha] = sklearn.metrics.roc_curve(y_test_lasso_roc, y_predict_test)
            
            if plot_learning_curve:
                plt.axis([0, train_sizes[-1], 0, 10])
                plt.xlabel('training size', fontsize=16)
                plt.ylabel('MSE', fontsize=16)
                plt.plot(train_sizes, train_errors, linewidth=2,
                         linestyle='solid', alpha=0.9, color=sns.xkcd_rgb['dark sky blue'], label='train')
                plt.plot(train_sizes, validation_errors, linewidth=2,
                         linestyle='solid', alpha=0.9, color=sns.xkcd_rgb['watermelon'], label='val')
                plt.axline((0, rmse_lasso_test), (len(X_train), rmse_lasso_test), linewidth=4, label='test',
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
            coeffs = np.divide(coeffs, mean_coeff)
            for i,_ in enumerate(coeffs):
                coeff = np.round(coeffs[i], 3)
                if not np.isclose(coeff, 0., atol=1e-10):
                    N,beta = self.N_beta_from_index(i)
                    observable += rf'(\tau_{{{N}}}^{{{beta}}})^{{{coeff}}} '
                    n_terms += 1
            print(f'Observable: {observable}')
            self.N_terms_lasso[alpha] = n_terms

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
    def plot_NN_epochs(self, n_epochs, history, label, K=None):
    
        epoch_list = range(1, n_epochs+1)
        loss = history.history['loss']
        acc = history.history['acc']
        val_loss = history.history['val_loss']
        val_acc = history.history['val_acc']
        
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
        if K:
            plt.savefig(os.path.join(self.output_dir_i, f'DNN_epoch_{label}_K{K}.pdf'))
        else:
            plt.savefig(os.path.join(self.output_dir_i, f'PFN_epoch_{label}.pdf'))
        plt.close()

    #---------------------------------------------------------------
    # Plot AUC as a function of K
    #---------------------------------------------------------------
    def plot_AUC_convergence(self, output_dir, auc_list, event_type, jetR, jet_pt_bin, R_max, suffix):
    
        plt.axis([0, self.K_list[-1], 0, 1])
        plt.title(rf'{event_type} event: $R = {jetR}, p_T = {jet_pt_bin}, R_{{max}} = {R_max}$', fontsize=14)
        plt.xlabel('K', fontsize=16)
        plt.ylabel('AUC', fontsize=16)

        for label,value in auc_list.items():
        
            if 'hard' in label:
                color = sns.xkcd_rgb['dark sky blue']
                label_suffix = ' (no background)'
            if 'combined' in label:
                color = sns.xkcd_rgb['medium green']
                label_suffix = ' (thermal background)'

            if 'pfn' in label:
                AUC_PFN = value[0]
                label = f'PFN{label_suffix}'
                plt.axline((0, AUC_PFN), (1, AUC_PFN), linewidth=4, label=label,
                           linestyle='solid', alpha=0.5, color=color)
            elif 'efn' in label:
                AUC_EFN = value[0]
                label = f'EFN{label_suffix}'
                plt.axline((0, AUC_EFN), (1, AUC_EFN), linewidth=4, label=label,
                           linestyle='solid', alpha=0.5, color=color)
            elif 'neural_network' in label:
                label = f'DNN{label_suffix}'
                plt.plot(self.K_list, value, linewidth=2,
                         linestyle='solid', alpha=0.9, color=color,
                         label=label)
                    
        plt.legend(loc='lower right', fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'AUC_convergence{suffix}.pdf'))
        plt.close()

    #--------------------------------------------------------------- 
    # Plot ROC curves
    #--------------------------------------------------------------- 
    def plot_roc_curves(self, roc_list, event_type, jetR, jet_pt_bin, R_max, suffix):
    
        plt.plot([0, 1], [0, 1], 'k--') # dashed diagonal
        plt.axis([0, 1, 0, 1])
        plt.title(rf'{event_type} event: $R = {jetR}, p_T = {jet_pt_bin}, R_{{max}} = {R_max}$', fontsize=14)
        plt.xlabel('False Positive Rate', fontsize=16)
        plt.ylabel('True Positive Rate', fontsize=16)
        plt.grid(True)
    
        for label,value in roc_list.items():
            if label in ['PFN', 'EFN', 'jet_mass', 'jet_angularity', 'jet_theta_g'] or 'multiplicity' in label:
                linewidth = 4
                alpha = 0.5
                linestyle = self.linestyles[label]
                color=self.colors[label]
                if label == 'jet_mass':
                    label = r'$m_{\mathrm{jet}}$'
                if label == 'jet_angularity':
                    label = r'$\lambda_1$'
                if label == 'jet_theta_g':
                    label = r'$\theta_{\mathrm{g}}$'
            elif 'DNN' in label:
                linewidth = 2
                alpha = 0.9
                linestyle = 'solid'
                color=self.colors[label]
            elif 'Lasso' in label:
                linewidth = 2
                alpha = 0.9
                linestyle = 'solid'
                color=self.colors[label]
                reg_param = float(re.search('= (.*)\)', label).group(1))
                n_terms = self.N_terms_lasso[reg_param]
                label += f', {n_terms} terms'
                
            FPR = value[0]
            TPR = value[1]
            plt.plot(FPR, TPR, linewidth=linewidth, label=label,
                     linestyle=linestyle, alpha=alpha, color=color)
                    
        plt.legend(loc='lower right', fontsize=10)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir_i, f'ROC_{suffix}.pdf'))
        plt.close()
 
    #--------------------------------------------------------------- 
    # Plot Significance improvement
    #--------------------------------------------------------------- 
    def plot_significance_improvement(self, roc_list, event_type, jetR, jet_pt_bin, R_max, suffix):
    
        plt.axis([0, 1, 0, 3])
        plt.title(rf'{event_type} event: $R = {jetR}, p_T = {jet_pt_bin}, R_{{max}} = {R_max}$', fontsize=14)
        plt.xlabel('True Positive Rate', fontsize=16)
        plt.ylabel('Significance improvement', fontsize=16)
        plt.grid(True)
            
        for label,value in roc_list.items():
            if label in ['PFN', 'EFN', 'jet_mass', 'jet_angularity', 'jet_theta_g'] or 'multiplicity' in label:
                linewidth = 4
                alpha = 0.5
                linestyle = self.linestyles[label]
            elif 'DNN' in label:
                linewidth = 2
                alpha = 0.9
                linestyle = 'solid'
            elif 'Lasso' in label:
                linewidth = 4
                alpha = 0.5
                linestyle = 'solid'
                
            FPR = value[0]
            TPR = value[1]
            plt.plot(TPR, TPR/np.sqrt(FPR+0.001), linewidth=linewidth, label=label,
                     linestyle=linestyle, alpha=alpha, color=self.colors[label])
         
        plt.legend(loc='best', fontsize=10)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir_i, f'Significance_improvement_{suffix}.pdf'))
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

            # Plot distributions
            fig, (ax1, ax2) = plt.subplots(nrows=2)
            plt.xlabel(rf'{qa_observable}', fontsize=14)
            max = np.amax(result_pythia)*1.2
            bins = np.linspace(0, max, 20)
            n_jewel,_,_ = ax1.hist(result_jewel,
                                   bins,
                                   histtype='step',
                                   density=True,
                                   label = 'JEWEL',
                                   linewidth=2,
                                   linestyle='-',
                                   alpha=0.5)
            n_pythia,_,_ = ax1.hist(result_pythia,
                                    bins,
                                    histtype='step',
                                    density=True,
                                    label = 'PYTHIA',
                                    linewidth=2,
                                    linestyle='-',
                                    alpha=0.5)
            legend = ax1.legend(loc='best', fontsize=14, frameon=False)

            # Plot ratio
            ratio = n_jewel/n_pythia
            ratio[np.isnan(ratio)] = 0
            x = (bins[1:] + bins[:-1]) / 2
            ax2.plot(x, ratio)
            
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
    # Main processing function
    #---------------------------------------------------------------
    def plot_training_data(self, K):
            
        # Print fraction of AA jets
        if K == min(self.K_list):
            print('Fraction of AA jets: {}'.format(np.sum(self.y_train)/self.y_train.size))
            print()
        
        print(f'Plotting input data, K={K}...')
        print()
        
        # Separate PYTHIA/JEWEL
        jewel_indices = self.y_train
        pythia_indices = 1 - self.y_train
        n_plot = int(self.n_train/10) # Plot a subset to save time/space
        X_Nsub_jewel = self.training_data[K]['X_Nsub_train'][jewel_indices.astype(bool)][:n_plot]
        X_Nsub_pythia = self.training_data[K]['X_Nsub_train'][pythia_indices.astype(bool)][:n_plot]

        # Construct dataframes for scatter matrix plotting
        df_jewel = pd.DataFrame(X_Nsub_jewel, columns=self.training_data[K]['feature_labels'])
        df_pythia = pd.DataFrame(X_Nsub_pythia, columns=self.training_data[K]['feature_labels'])
        
        # Add label columns to each df to differentiate them for plotting
        df_jewel['generator'] = np.repeat('JEWEL', X_Nsub_jewel.shape[0])
        df_pythia['generator'] = np.repeat('PYTHIA', X_Nsub_pythia.shape[0])
                
        # Plot scatter matrix
        df = df_jewel.append(df_pythia)
        g = sns.pairplot(df, corner=True, hue='generator')
        #g.legend.set_bbox_to_anchor((0.75, 0.75))
        #plt.savefig(os.path.join(self.output_dir_i, f'training_data_scatter_matrix_K{K}.png'), dpi=50)
        plt.savefig(os.path.join(self.output_dir_i, f'training_data_scatter_matrix_K{K}.pdf'))
        plt.close()
        
        # Plot correlation matrix
        df.drop(columns=['generator'])
        corr_matrix = df.corr()
        sns.heatmap(corr_matrix)
        plt.savefig(os.path.join(self.output_dir_i, f'training_data_correlation_matrix_K{K}.pdf'))
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
