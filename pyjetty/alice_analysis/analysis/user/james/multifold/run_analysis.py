#! /usr/bin/env python

'''
This class steers analysis of a jet substructure analyses using multifold.

Author: James Mulligan (james.mulligan@berkeley.edu)
'''

import sys
import os
import argparse
import yaml
import numpy as np

try:
    import ROOT
except ImportError:
    pass

import seaborn as sns
sns.set_context('paper', rc={'font.size':18, 'axes.titlesize':18, 'axes.labelsize':18, 'text.usetex':True})

from sklearn.model_selection import train_test_split

import tensorflow as tf
from tensorflow import keras

# Suppress some annoying warnings
np.finfo(np.dtype("float32")) 
np.finfo(np.dtype("float64"))

sys.path.append('../../../../../..')
from pyjetty.alice_analysis.analysis.base import common_base
from pyjetty.alice_analysis.analysis.user.james.multifold import multifold_utils
from pyjetty.alice_analysis.analysis.user.james.multifold import multifold_plot_utils

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class RunAnalysis(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, config_file='', **kwargs):
        super(RunAnalysis, self).__init__(**kwargs)
        self.config_file = config_file

        # Initialize utils class
        self.utils = multifold_utils.MultiFoldUtils()

        # Initialize utils class
        self.plot_utils = multifold_plot_utils.MultiFoldPlotUtils()

        # Initialize yaml config
        self.initialize_config()
    
        print(self)
        print()

    #---------------------------------------------------------------
    # Initialize config file into class members
    #---------------------------------------------------------------
    def initialize_config(self):

        # Read config file
        with open(self.config_file, 'r') as stream:
            config = yaml.safe_load(stream)

        self.jetR = config['jetR'][0]
        if 'constituent_subtractor' in config:
            self.is_pp = False
            self.R_max = config['constituent_subtractor']['main_R_max']
        else:
            self.is_pp = True
            self.R_max = None

        # Construct a dict of information for each observable:
        #     self.observable_info[observable]['obs_settings'] = [obs_setting1, obs_setting2, ...]
        #                                     ['obs_grooming_settings'] = [grooming_setting1, grooming_setting1, ...]
        #                                     ['obs_labels'] = [obs_label1, obs_label2, ...]
        #                                     ['obs_keys'] = [obs_key1, obs_key2, ...]
        #                                     ['pt_bins_truth'] = [pt_bins_truth1, pt_bins_truth2, ...]
        #                                     ['pt_bins_det'] = [pt_bins_det1, pt_bins_det2, ...]
        #                                     ['obs_bins_truth'] = [obs_bins_truth1, obs_bins_truth2, ...]
        #                                     ['obs_bins_det'] = [obs_bins_det1, obs_bins_det2, ...]
        #                                     ['xtitle'] = xtitle
        #                                     ['ytitle'] = ytitle
        #                                     ['pt_bins_reported'] = pt_bins_reported
        # Each list corresponds to the list of observable subconfigs, and each key with a single entry
        #   corresponds to common settings for all subconfigs.
        analysis_observable_subconfigs = config['analysis_observables']
        self.observables = list(analysis_observable_subconfigs.keys())

        self.observable_info = self.recursive_defaultdict()
        for observable in self.observables:
            
            obs_subconfig_list = analysis_observable_subconfigs[observable]
            obs_config_dict = dict((key,config[observable][key]) for key in obs_subconfig_list)
            self.observable_info[observable]['n_subobservables'] = len(obs_subconfig_list)
            self.observable_info[observable]['obs_settings'] = self.utils.obs_settings(observable, obs_config_dict, obs_subconfig_list)
            self.observable_info[observable]['obs_grooming_settings'] = self.utils.grooming_settings(obs_config_dict)    
            self.observable_info[observable]['obs_labels'] = [self.utils.obs_label(self.observable_info[observable]['obs_settings'][i], 
                                                                                   self.observable_info[observable]['obs_grooming_settings'][i])
                                                              for i in range(len(obs_subconfig_list))]
            self.observable_info[observable]['obs_keys'] = [f'{observable}_{obs_label}' if obs_label else observable for obs_label in self.observable_info[observable]['obs_labels']]
            self.observable_info[observable]['xtitle'] = config[observable]['common_settings']['xtitle']
            self.observable_info[observable]['ytitle'] = config[observable]['common_settings']['ytitle']

            # Load HEPData binnings and TGraphs
            # These consist of nested lists: 
            #   self.observable_info[observable]['hepdata_tgraph'] = [ [config1_pt1, config1_pt2, ...], [config2_pt1, config2_pt2, ...], ... ]
            self.observable_info[observable]['obs_bins'] = [self.utils.bins_from_hepdata(obs_config_dict[obs_subconfig_list[i]], self.jetR) 
                                                                                   for i in range(len(obs_subconfig_list))]
            self.observable_info[observable]['hepdata_tgraph'] = [self.utils.tgraph_from_hepdata(obs_config_dict[obs_subconfig_list[i]], self.jetR) 
                                                                                           for i in range(len(obs_subconfig_list))]
                                                                  
            # For jet_pt, store the binnings we will report
            if 'pt_bins_reported' in config[observable]['common_settings']:
                pt_bins = config[observable]['common_settings']['pt_bins_reported']
                pt_bin_pairs = [[pt_bins[i],pt_bins[i+1]] for i in range(len(pt_bins)-1)]
                pt_bin_pairs_nested = [[pt_bin_pair] for pt_bin_pair in pt_bin_pairs]
                self.observable_info[observable]['pt_bins_reported'] = pt_bins
                self.observable_info[observable]['pt_bins_reported_pairs'] = pt_bin_pairs
                self.observable_info[observable]['pt_bins_reported_pairs_nested'] = pt_bin_pairs_nested

        # Directory to write analysis output to
        self.output_dir = os.path.expandvars(config['output_dir'])

        # Set which analysis steps to perform
        self.do_plot_observables = config['do_plot_observables']
        self.do_unfolding = config['do_unfolding']
        self.do_plot_unfolding_results = config['do_plot_unfolding_results']
        self.do_systematics = config['do_systematics']
        self.do_plot_final_result = config['do_plot_final_result']

        # List of systematic variations to perform
        self.systematics_list = config['systematics_list']

        self.n_iterations = config['n_iterations']

        # Load paths to processing output, to be unfolded
        self.main_data = os.path.expandvars(config['main_data'])
        self.response = {}
        self.response['main'] = os.path.expandvars(config['main_response'])
        if 'trkeff' in self.systematics_list:
            self.response['trkeff'] = os.path.expandvars(config['trkeff_response'])
        if 'fastsim_generator0' in self.systematics_list:
            self.response['fastsim_list'] = [os.path.expandvars(x) for x in config['fastsim_response']]

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def run_analysis(self):

        # Load response for each systematic, and plot/unfold
        for systematic in self.systematics_list:

            if self.do_plot_observables or self.do_unfolding:
                self.load_data(systematic)
                
            if self.do_plot_observables:
                self.plot_observables(systematic)

            if self.do_unfolding:
                self.perform_unfolding(systematic)

            if self.do_plot_unfolding_results:
                self.plot_unfolding_results(systematic)

        # Compute systematic uncertainties
        if self.do_systematics:
            self.compute_systematics()

        # Plot final results
        if self.do_plot_final_result:
            self.plot_results()

    #---------------------------------------------------------------
    # Load aggregated results for a given systematic for both data and MC
    # The results will be stored in a dict of 1D numpy arrays:
    #   - self.results[data_type][obs_key] = np.array([...])
    #       where data_type is: 'data', 'mc_det_matched', 'mc_truth_matched'
    #       and obs_key is self.observable_info[observable]['obs_key'][i] (plus 'pt_hat_scale_factors' for MC)
    #---------------------------------------------------------------
    def load_data(self, systematic):
        print()
        print('Loading data...')
        results = {}
        results['data'] = self.utils.read_data(self.main_data)['data']
        print()
        print(f'Loading MC (systematic={systematic})...')
        results_mc = self.utils.read_data(self.response[systematic])
        results['mc_det_matched'] = results_mc['det_matched']
        results['mc_truth_matched'] = results_mc['truth_matched']
        print()

        self.n_jets = {}
        for data_type in results.keys():
            self.n_jets[data_type] = results[data_type][self.observable_info[self.observables[0]]['obs_keys'][0]].shape[0]
        n_jets_total_data = self.n_jets['data'] 
        n_jets_total_mc = self.n_jets['mc_truth_matched']
        print(f'The following observables are available to analyze: (n_jets={n_jets_total_data} in data, n_jets={n_jets_total_mc} in mc)')
        [print(obs_key) for obs_key in results['data'].keys()]
        print()

        # Shuffle all of the arrays 
        #   - Data/MC shuffled separately (same shuffling for matched det/truth MC)
        #   - Same shuffling for all observables
        print('Shuffling dataset...')
        idx_data = np.random.permutation(self.n_jets['data'])
        idx_mc = np.random.permutation(self.n_jets['mc_truth_matched'])
        idx = {'data': idx_data,
               'mc_det_matched': idx_mc,
               'mc_truth_matched': idx_mc}

        self.results = self.recursive_defaultdict()
        for data_type in results.keys():
            for observable in self.observables:
                for i in range(self.observable_info[observable]['n_subobservables']):
                    obs_key = self.observable_info[observable]['obs_keys'][i]
                    print(f'    {obs_key}')

                    # Check number of jets is as expected for each observable
                    if results[data_type][obs_key].shape[0] != self.n_jets[data_type] or idx[data_type].shape[0] != self.n_jets[data_type]:
                        shape = results[data_type][obs_key].shape
                        raise ValueError(f'{data_type} {obs_key} in data has unexpected number of jets: {shape}')

                    # Shuffle
                    self.results[data_type][obs_key] = results[data_type][obs_key][idx[data_type]]
            if 'mc' in data_type:
                self.results[data_type]['pt_hat_scale_factors'] = results[data_type]['pt_hat_scale_factors'][idx[data_type]]
        print('Done.')
        print()
        print('We will analyze the following observables:')
        [print(obs_key) for obs_key in self.results['data'].keys() if obs_key != 'pt_hat_scale_factors']
        print()

    #---------------------------------------------------------------
    # Plot observables before unfolding
    #---------------------------------------------------------------
    def plot_observables(self, systematic):
        output_dir = self.utils.create_output_subdir(self.output_dir, f'{systematic}/plot_observables')

        # First, apply pt cuts to a copy of results
        # TODO: bin both mc_det_matched and mc_truth_matched by e.g. truth-level pt
        print('Binning results into pt bins to make some plots...')
        print()
        pt_bin_pairs = self.observable_info['jet_pt']['pt_bins_reported_pairs']
        pt_bin_pairs_nested = self.observable_info['jet_pt']['pt_bins_reported_pairs_nested']
        variables = [['jet_pt']] * len(pt_bin_pairs_nested)
        mask_data_type_dict = {'sim_det': 'sim_truth'}
        results_dict = self.utils.apply_cut(self.results, self.n_jets, variables, pt_bin_pairs_nested,
                                            mask_data_type_dict=mask_data_type_dict, n_jets_max=1000000)

        # Print number of jets passing each cut
        for data_type in results_dict.keys():
            n_jets = [results_dict[data_type]['jet_pt'][f'jet_pt{x[0]}-{x[1]}'].shape[0] for x in pt_bin_pairs]
            print(data_type)
            [print(f'  pt={pt_bin_pairs[i]} has n_jets={n_jets[i]}') for i in range(len(n_jets))]

        # Plot 1D projections of each observable
        print()
        print(f'Plotting 1D projections...')
        self.plot_utils.plot_1D_projections_before_unfolding(results_dict, pt_bin_pairs, self.observables, self.observable_info, self.jetR, output_dir)

        # Plot pairplot between each pair of observables
        print()
        print(f'Plotting pairplot...')
        self.plot_utils.plot_pairplot(results_dict, pt_bin_pairs, self.observables, self.observable_info, output_dir)
        print('Done!')
        print()

    #---------------------------------------------------------------
    # Perform unfolding for a single systematic
    #
    # The basic idea of OmniFold is that since an ML classifier approximates a likelihood 
    # ratio between the two classes, we can use the classifier output thresholds to
    # reweight one class into the other.
    #
    # Throughout, one can manipulate weighted distributions by simply propagating the 
    # weights with the unweighted distributions, and applying them during appropriate 
    # operations, namely when computing the loss or when filling histograms.
    #
    # We will follow notation of the OmniFold paper, where:
    #    - omega_n: MC-det weights for iteration n
    #    - nu_n: MC-truth weights for iteration n -- such that P_unfolded(n) = P(MC-truth, weights=nu_n)
    # We interchangeably refer to 'omega' as 'weights_det', and 'nu' and 'weights_truth'.
    #
    # I personally find the 'push'/'pull' notation a bit confusing, so I neglect those labels and just assume
    #   it is understood that there is an event-by-event association of MC-truth and MC-det which share
    #   the same weight.
    #
    # OmniFold code based in part on:
    #  - https://github.com/ftoralesacosta/AI4EIC_Omnfold
    #---------------------------------------------------------------  
    def perform_unfolding(self, systematic):
        output_dir = self.utils.create_output_subdir(self.output_dir, f'{systematic}/unfolding_results')
        print('Perform unfolding...')
        print()

        # Construct ndarray of observables for each data_type
        ndarray_dict, cols_dict, n_jets_dict = self.utils.convert_results_to_ndarray(self.results, self.observables, self.observable_info)

        # Option to balance samples between data and MC
        balance_samples = True
        if balance_samples:
            n_jets = min(n_jets_dict['data'], n_jets_dict['mc_det_matched'])
            print(f'n_jets: {n_jets}')
            for data_type in self.results.keys():
                ndarray_dict[data_type] = ndarray_dict[data_type][:n_jets,:]
                n_jets_dict[data_type] = n_jets

        # Perform MultiFold
        print('  Performing MultiFold...')
        n_observables = ndarray_dict['data'].shape[1]

        # Initialize training data
        # These are fixed independent of iteration -- the weights will vary with iteration
        multifold_result_dict = self.recursive_defaultdict()
        multifold_result_dict['nature_det'] = ndarray_dict['data']
        multifold_result_dict['sim_det'] = ndarray_dict['mc_det_matched']
        multifold_result_dict['sim_truth'] = ndarray_dict['mc_truth_matched']
        multifold_result_dict['columns'] = cols_dict[list(self.results.keys())[0]]
        multifold_result_dict['pt_hat_scale_factors'] = self.results['mc_truth_matched']['pt_hat_scale_factors'][:n_jets]

        # Initialize labels
        labels_dict = {}
        labels_dict['nature_det'] = np.ones(n_jets_dict['data'])
        labels_dict['sim_det'] = np.zeros(n_jets_dict['mc_det_matched'])
        labels_dict['sim_truth'] = np.zeros(n_jets_dict['mc_truth_matched'])
        labels_dict['nature_truth'] = np.ones(n_jets_dict['mc_truth_matched'])

        # Define ML model
        inputs = keras.layers.Input((n_observables, ))
        hidden_layer_1 = keras.layers.Dense(50, activation='relu')(inputs)
        hidden_layer_2 = keras.layers.Dense(50, activation='relu')(hidden_layer_1)
        hidden_layer_3 = keras.layers.Dense(50, activation='relu')(hidden_layer_2)
        outputs = keras.layers.Dense(1, activation='sigmoid')(hidden_layer_3)
        model = keras.models.Model(inputs=inputs, outputs=outputs)

        # Iterate
        for n in range(1, self.n_iterations+1):
            print(f'    Iteration {n}')

            #------------------
            # Step 1: 
            #   We wish to compute: omega_{n} = nu_{n-1} L[(1,nature_det),(nu_{n-1},sim_det)]
            #      where L[(),()] is the likelihood ratio between the two (weighted) event ensembles,
            #      which we approximate by an ML classifier.
            #
            #   To compute L, we therefore train a classifier to discriminate 'nature_det' from 'sim_det',
            #      where 'sim_det' is reweighted by 'weights_truth', which are the truth weights from the previous iteration
            #
            #   We set the initial truth weights nu_0 to the pt-hat-scale-factors.
            #   Ultimately, these will then be implicitly included in our final truth weights.
            if n == 1:
                multifold_result_dict[f'weights_truth_iteration{n-1}'] = multifold_result_dict['pt_hat_scale_factors']

            training_data = np.concatenate((multifold_result_dict['nature_det'], multifold_result_dict['sim_det']))
            labels = np.concatenate((labels_dict['nature_det'],labels_dict['sim_det']))
            weights = np.concatenate((multifold_result_dict[f'weights_truth_iteration{n-1}'], np.ones(n_jets_dict['data'])))
            #print(f'    training_data: {training_data.shape} -- {training_data}')
            #print(f'    labels: {labels.shape} -- {labels}')
            #print(f'    weights: {weights.shape} -- {weights}')

            test_size = 0.2
            X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(training_data, labels, weights, test_size=test_size)

            # Stack the weights with the labels, in order to apply the weights when computing the loss
            y_train = np.stack((y_train, w_train), axis=1)
            y_test = np.stack((y_test, w_test), axis=1)   

            model.compile(loss=self.weighted_binary_crossentropy,
                          optimizer='Adam',
                          metrics=['accuracy'])

            n_epochs = 2
            batch_size = 10000
            model.fit(X_train,
                      y_train,
                      epochs=n_epochs,
                      batch_size=batch_size,
                      validation_data=(X_test, y_test))

            # Compute omega_{n} ('weights_det') from the previous iteration's truth weights and the classifier:
            #    omega_{n} = nu_{n-1} * L[(),()]
            # To do this, get the likelihood ratio from the classifier for a sample of simulated events.
            # The classifier returns a probability f of being true --  we want to convert this 
            #   to a ratio P(true)/P(false) = P(true)/(1 - P(true))
            f = model.predict(multifold_result_dict['sim_det'], batch_size=batch_size)
            eps = 1.e-6
            weights_update = f / (1. - f + eps)
            weights_update = np.squeeze(np.nan_to_num(weights_update))
            multifold_result_dict[f'weights_det_iteration{n}'] = multifold_result_dict[f'weights_truth_iteration{n-1}'] * weights_update

            #------------------
            # Step 2: 
            #   We wish to compute: nu_{n} = nu_{n-1} L[(omega_n,sim_truth),(nu_{n-1},sim_truth)]
            #      where L[(),()] is the likelihood ratio between the two (weighted) event ensembles,
            #      which we approximate by an ML classifier.
            #
            #   To compute L, we therefore train a classifier to discriminate the two 'sim_truth' samples,
            #       one weighted by 'weights_det' from Step 1, and the other by 'weights_truth' from the previous iteration.
            training_data = np.concatenate((multifold_result_dict['sim_truth'], multifold_result_dict['sim_truth']))
            labels = np.concatenate((labels_dict['sim_truth'], labels_dict['nature_truth']))
            weights = np.concatenate((multifold_result_dict[f'weights_det_iteration{n}'], multifold_result_dict[f'weights_truth_iteration{n-1}']))
            #print(f'    training_data: {training_data.shape} -- {training_data}')
            #print(f'    labels: {labels.shape} -- {labels}')
            #print(f'    weights: {weights.shape} -- {weights}')

            X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(training_data, labels, weights, test_size=test_size)

            # Stack the weights with the labels, in order to apply the weights when computing the loss
            y_train = np.stack((y_train, w_train), axis=1)
            y_test = np.stack((y_test, w_test), axis=1)   

            model.compile(loss=self.weighted_binary_crossentropy,
                          optimizer='Adam',
                          metrics=['accuracy'])

            batch_size = 2000
            model.fit(X_train,
                      y_train,
                      epochs=n_epochs,
                      batch_size=batch_size,
                      validation_data=(X_test, y_test))

            # Compute nu_{n} ('weights_truth') from the previous iteration's truth weights and the classifier:
            #    nu_{n} = nu_{n-1} * L[(),()]
            # To do this, get the likelihood ratio from the classifier for a sample of simulated events.
            # The classifier returns a probability f of being true --  we want to convert this 
            #   to a ratio P(true)/P(false) = P(true)/(1 - P(true))
            # In general, this step causes the new generator to drift away from the data from step 1
            #
            # Note: we could alternately have trained the classifier using nu_0 rather than nu_{n-1} for
            #       the truth generator weights, in which case the coeffs below would be obtained simply
            #       as the classifier output rather without multiplying by 'weights_truth_iteration{n-1}',
            #       i.e. determining the coefficient in one go rather than iteratively updating the coeff.
            f = model.predict(multifold_result_dict['sim_truth'], batch_size=batch_size)
            weights_update = f / (1. - f + eps)
            weights_update = np.squeeze(np.nan_to_num(weights_update))
            multifold_result_dict[f'weights_truth_iteration{n}'] = multifold_result_dict[f'weights_truth_iteration{n-1}'] * weights_update

        # Save results to file
        print()
        print('Saving unfolding result...')
        self.utils.write_data(multifold_result_dict, output_dir, filename = f'unfolding_results.h5')
        print()

    #---------------------------------------------------------------
    # Binary crossentropy for classifying two samples with weights
    # Weights are "hidden" by stacking both the labels and weights in y
    #---------------------------------------------------------------
    def weighted_binary_crossentropy(self, y_true, y_pred):

        weights = tf.gather(y_true, [1], axis=1) # event weights
        y_true = tf.gather(y_true, [0], axis=1) # actual y_true for loss
        
        # Clip the prediction value to prevent NaN's and Inf's
        epsilon = keras.backend.epsilon()
        y_pred = keras.backend.clip(y_pred, epsilon, 1. - epsilon)

        t_loss = -weights * ((y_true) * keras.backend.log(y_pred) +
                            (1 - y_true) * keras.backend.log(1 - y_pred))
        
        return keras.backend.mean(t_loss)

    #---------------------------------------------------------------
    # Plot unfolding results
    #---------------------------------------------------------------
    def plot_unfolding_results(self, systematic):
        print('Plotting unfolding results...')
        output_dir = self.utils.create_output_subdir(self.output_dir, f'{systematic}/unfolding_results')

        unfolding_results = self.utils.read_data(os.path.join(output_dir, 'unfolding_results.h5'))
        print()

        # Convert the results to a dict of 1D numpy arrays, in order to perform cuts etc
        results_dict, n_jets = self.utils.convert_unfolding_results_to_dict(unfolding_results)

        # Apply pt cuts to a copy of results
        print('Binning results into pt bins to make some plots...')
        print()
        pt_bin_pairs_nested = self.observable_info['jet_pt']['pt_bins_reported_pairs_nested']
        variables = [['jet_pt']] * len(pt_bin_pairs_nested)
        results_dict = self.utils.apply_cut(results_dict, n_jets, variables, pt_bin_pairs_nested)

        # Plot 1D projections of each observable
        print(f'Plotting 1D projections...')
        pt_bin_pairs = self.observable_info['jet_pt']['pt_bins_reported_pairs']
        self.plot_utils.plot_unfolding_results(results_dict, pt_bin_pairs, self.observables, self.observable_info, self.jetR, output_dir)

    #---------------------------------------------------------------
    # Plot final results
    #---------------------------------------------------------------
    def plot_results(self):
        print('Plotting final results...')
        output_dir = self.utils.create_output_subdir(self.output_dir, 'final_results')

    #---------------------------------------------------------------
    # This function is called once for each subconfiguration
    #---------------------------------------------------------------
    def plot_single_result(self, jetR, obs_label, obs_setting, grooming_setting):    
        print(f'Plot final results for {self.observable}: R = {jetR}, {obs_label}')

#----------------------------------------------------------------------
if __name__ == '__main__':

  # Define arguments
  parser = argparse.ArgumentParser(description='Unfold theta_g distribution')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='../../../../config/multifold/pp/analysis.yaml',
                      help='Path of config file for analysis')

  # Parse the arguments
  args = parser.parse_args()

  print('Configuring...')
  print('configFile: \'{0}\''.format(args.configFile))

  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = RunAnalysis(config_file = args.configFile)
  analysis.run_analysis()