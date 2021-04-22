#!/usr/bin/env python3

"""
Example class to read quark-gluon dataset
"""

import os
import argparse
import yaml
import h5py

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
    def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
        super(common_base.CommonBase, self).__init__(**kwargs)
        
        self.config_file = config_file
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        
        # Initialize config file
        self.initialize_config()
        
        print(self)
        print()
        
        # Read input variables
        with h5py.File(os.path.join(self.output_dir, 'nsubjettiness.h5'), 'r') as hf:
            self.y = hf['y'][:]
            if 'X' in hf:
                self.X_particles = hf['X'][:]
            self.X_Nsub = hf['X_Nsub'][:]
            self.N_list = hf['N_list'][:]
            self.beta_list = hf['beta_list'][:]

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
                    
    #---------------------------------------------------------------
    # Initialize config file into class members
    #---------------------------------------------------------------
    def initialize_config(self):
    
        # Read config file
        with open(self.config_file, 'r') as stream:
          config = yaml.safe_load(stream)
          
        self.n_train = config['n_train']
        self.n_val = config['n_val']
        self.n_test = config['n_test']
        self.n_total = self.n_train + self.n_val + self.n_test
        self.test_frac = 1. * self.n_test / self.n_total
        
        self.K_list = config['K']
        
        # Initialize model-specific settings
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
                self.model_settings[model]['metrics'] = config[model]['metrics']
                self.model_settings[model]['random_state'] = config[model]['random_state']
            
            if model == 'pfn':
                self.model_settings[model]['Phi_sizes'] = tuple(config[model]['Phi_sizes'])
                self.model_settings[model]['F_sizes'] = tuple(config[model]['F_sizes'])
                self.model_settings[model]['epochs'] = config[model]['epochs']
                self.model_settings[model]['batch_size'] = config[model]['batch_size']
                self.model_settings[model]['use_pids'] = config[model]['use_pids']
                
            if model == 'lasso':
                self.model_settings[model]['alpha'] = config[model]['alpha']

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def analyze_pp_aa(self):
    
        # Plot the input data
        [self.plot_training_data(K) for K in self.K_list]

        # Train ML models
        for model in self.models:
        
            model_settings = self.model_settings[model]
            
            for K in self.K_list:
            
                if model == 'linear':
                    self.fit_linear_model(K, model_settings)
                if model == 'random_forest':
                    self.fit_random_forest(K, model_settings)
                if model == 'neural_network':
                    self.fit_neural_network(K, model_settings)
                if model == 'pfn':
                    self.fit_pfn(model_settings)
                if model == 'lasso':
                    self.fit_lasso(K, model_settings)
                
        # Plot ROC curve and significance improvement
        self.plot_roc_curves()
        self.plot_significance_improvement()
        
    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def plot_training_data(self, K):
            
        # Print fraction of q/g jets
        if K == min(self.K_list):
            print('Fraction of AA jets: {}'.format(np.sum(self.y_train)/self.y_train.size))
            print()
        
        print(f'Plotting input data, K={K}...')
        print()

        # Plot scatter matrix
        df = pd.DataFrame(self.training_data[K]['X_Nsub_train'], columns=self.training_data[K]['feature_labels'])
        sns.pairplot(df, kind='hist', corner=True)
        plt.savefig(os.path.join(self.output_dir, f'training_data_scatter_matrix_K{K}.png'), dpi=200)
        plt.close()
        
        # Plot correlation matrix
        corr_matrix = df.corr()
        sns.heatmap(corr_matrix)
        plt.savefig(os.path.join(self.output_dir, f'training_data_correlation_matrix_K{K}.pdf'))
        plt.close()
       
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

        # Train DNN - need validation set here (use test set for now)
        DNN.fit(self.training_data[K]['X_Nsub_train'], self.y_train, epochs=model_settings['epochs'],
                validation_data=(self.training_data[K]['X_Nsub_test'], self.y_test))
        
        # Get predictions for validation data set
        y_Nsub_test_preds_DNN = DNN.predict(self.training_data[K]['X_Nsub_test']).reshape(-1)
        
        # Get AUC
        Nsub_auc_DNN = sklearn.metrics.roc_auc_score(self.y_test, y_Nsub_test_preds_DNN)
        print(f'AUC = {Nsub_auc_DNN} (validation set)')
        
        # Get ROC curve results
        self.roc_curve_dict['DNN'][K] = sklearn.metrics.roc_curve(self.y_test, y_Nsub_test_preds_DNN)
    
    #---------------------------------------------------------------
    # Fit ML model -- 4. Deep Set/Particle Flow Networks
    #---------------------------------------------------------------
    def fit_pfn(self, model_settings):

        # convert labels to categorical
        Y_PFN = energyflow.utils.to_categorical(self.y, num_classes=2)
        
        # preprocess by centering jets and normalizing pts
        X_PFN = self.X_particles
        for x_PFN in X_PFN:
            mask = x_PFN[:,0] > 0
            yphi_avg = np.average(x_PFN[mask,1:3], weights=x_PFN[mask,0], axis=0)
            x_PFN[mask,1:3] -= yphi_avg
            x_PFN[mask,0] /= x_PFN[:,0].sum()
            
        # Convert (E,px,py,pz) to (pT,y,phi,m)
        X_PFN = energyflow.ptyphims_from_p4s(self.X_particles)[:,:,:] # Note: 4th entry is m, not PID .. could rewrite routine
        
        # handle particle id channel !! Note: If changed to pT,y,phi,m the 4th component is not PID but m .. fix later
        #if model_settings['use_pids']:
        #    self.my_remap_pids(X_PFN)
        #else:
        #    X_PFN = X_PFN[:,:,:3]

        # Split data into train, val and test sets
        (X_PFN_train, X_PFN_val, X_PFN_test,Y_PFN_train, Y_PFN_val, Y_PFN_test) = energyflow.utils.data_split(X_PFN, Y_PFN,
                                                                                             val=self.n_val, test=self.n_test)
        # build architecture
        pfn = energyflow.archs.PFN(input_dim=X_PFN.shape[-1],
                                   Phi_sizes=model_settings['Phi_sizes'],
                                   F_sizes=model_settings['F_sizes'])

        # train model
        pfn.fit(X_PFN_train, Y_PFN_train,
                epochs=model_settings['epochs'],
                batch_size=model_settings['batch_size'],
                validation_data=(X_PFN_val, Y_PFN_val),
                verbose=1)
        
        # get predictions on test data
        preds_PFN = pfn.predict(X_PFN_test, batch_size=1000)

        # Get AUC and ROC curve + make plot
        auc_PFN = sklearn.metrics.roc_auc_score(Y_PFN_test[:,1], preds_PFN[:,1])
        print('Particle Flow Networks/Deep Sets: AUC = {} (test set)'.format(auc_PFN))
        
        self.roc_curve_dict['PFN'] = sklearn.metrics.roc_curve(Y_PFN_test[:,1], preds_PFN[:,1])
        
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        
        # Now we compare the PFN ROC curve to single observables

        # 1. Jet mass (Note: This takes in (pt,y,phi) and converts it to 4-vectors and computes jet mass)
        #             (Note: X_PFN_train is centered and normalized .. should be ok)
        masses = np.asarray([energyflow.ms_from_p4s(energyflow.p4s_from_ptyphims(x).sum(axis=0)) for x in X_PFN_train])
        self.roc_curve_dict['Jet_mass'] = sklearn.metrics.roc_curve(Y_PFN_train[:,1], -masses)
        
        # 2. Multiplicity (Is this a useful observable for pp vs AA?)
        mults = np.asarray([np.count_nonzero(x[:,0]) for x in X_PFN_train])
        self.roc_curve_dict['Multiplicity'] = sklearn.metrics.roc_curve(Y_PFN_train[:,1], -mults)

        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        
    #---------------------------------------------------------------
    # Fit ML model -- 5. Lasso regression
    #---------------------------------------------------------------
    def fit_lasso(self, K, model_settings):
        print(f'Training Lasso regression...')
        
        # Take the logarithm of the data and labels
        # Not taking log of test labels to make ROC curve later which requires integers for the test data
        X_train_lasso = np.log(self.training_data[K]['X_Nsub_train'])
        X_test_lasso = np.log(self.training_data[K]['X_Nsub_test'])
        
        eps = .01
        y_train_lasso = np.log(eps + (1. - 2. * eps) * self.y_train)
        y_test_lasso = self.y_test
        
        # Use Lasso regression
        # The parameter alpha multiplies to L1 term
        lasso_clf = sklearn.linear_model.Lasso(alpha = model_settings['alpha'])
        
        # Fit model
        lasso_clf.fit(X_train_lasso, y_train_lasso)
        
        # Print coefficients of log(tau_N^beta) or exponents of the product observable
        print('Coefficients of lasso regression: c_N^beta * log(tau_n^beta):\n {}'.format(lasso_clf.coef_))

        # Make predictions for test set
        preds_lasso = lasso_clf.predict(X_test_lasso)
        
        # Compute AUC
        auc_lasso = sklearn.metrics.roc_auc_score(y_test_lasso, preds_lasso)
        print(f'AUC = {auc_lasso} (test set)')
        
        # ROC curve
        self.roc_curve_dict['Lasso'] = sklearn.metrics.roc_curve(y_test_lasso, preds_lasso)

        # Seems we would have to scan different alpha and K and see if some observable has a nice interpretation??

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
    # Plot ROC curves
    #--------------------------------------------------------------- 
    def plot_roc_curves(self):
    
        plt.plot([0, 1], [0, 1], 'k--') # dashed diagonal
        plt.axis([0, 1, 0, 1])
        plt.xlabel('False Positive Rate', fontsize=16)
        plt.ylabel('True Positive Rate', fontsize=16)
        plt.grid(True)
    
        for label,value in self.roc_curve_dict.items():
            
            if label in ['PFN', 'Jet_mass', 'Multiplicity','Lasso']:
                FPR = value[0]
                TPR = value[1]
                plt.plot(FPR, TPR, linewidth=2, label=label)
            else:
                for K in self.K_list:
                    FPR = value[K][0]
                    TPR = value[K][1]
                    plt.plot(FPR, TPR, linewidth=2, label=f'{label}, K={K}')
                    
        plt.legend(loc='lower right')
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'ROC.pdf'))
        plt.close()
        
    #--------------------------------------------------------------- 
    # Plot Significance improvement
    #--------------------------------------------------------------- 
    def plot_significance_improvement(self):
    
        plt.axis([0, 1, 0, 3])
        plt.xlabel('True Positive Rate', fontsize=16)
        plt.ylabel('Significance improvement', fontsize=16)
        plt.grid(True)
    
        for label,value in self.roc_curve_dict.items():
            
            if label in ['PFN', 'Jet_mass', 'Multiplicity','Lasso']:
                FPR = value[0]
                TPR = value[1]
                plt.plot(TPR, TPR/np.sqrt(FPR+0.001), linewidth=2, label=label)
            else:
                for K in self.K_list:
                    FPR = value[K][0]
                    TPR = value[K][1]
                    plt.plot(TPR, TPR/np.sqrt(FPR+0.001), linewidth=2, label=f'{label}, K={K}')
                    
        plt.legend(loc='lower right')
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'Significance_improvement.pdf'))
        plt.close()
        
    #---------------------------------------------------------------
    # Plot confusion matrix
    # Note: not normalized to relative error
    #---------------------------------------------------------------
    def plot_confusion_matrix(self, y_train, y_predict_train, label):
    
        confusion_matrix = sklearn.metrics.confusion_matrix(y_train, y_predict_train)
        sns.heatmap(confusion_matrix)
        plt.savefig(os.path.join(self.output_dir, f'confusion_matrix_{label}.pdf'))
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

    # Parse the arguments
    args = parser.parse_args()

    print('Configuring...')
    print('configFile: \'{0}\''.format(args.configFile))
    print('ouputDir: \'{0}\"'.format(args.outputDir))

    # If invalid configFile is given, exit
    if not os.path.exists(args.configFile):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    analysis = AnalyzePPAA(config_file=args.configFile, output_dir=args.outputDir)
    analysis.analyze_pp_aa()
