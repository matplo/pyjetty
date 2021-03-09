#!/usr/bin/env python3

"""
Example class to read quark-gluon dataset
"""

import os

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
from energyflow.utils import data_split, standardize, to_categorical, remap_pids
from energyflow.archs import PFN

# sklearn
import sklearn
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import SGDClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import RandomForestClassifier

# Tensorflow and Keras
import tensorflow as tf
from tensorflow import keras

# Base class
from pyjetty.alice_analysis.process.base import common_base

################################################################
class ProcessQG(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
        super(common_base.CommonBase, self).__init__(**kwargs)
        
        self.output_dir = '.'
        
        # Load labeled data -- quark/gluon data set from energyflow data sets is up to 200.000 (?)
        self.train = 150000
        self.val = 20000
        self.test = 30000
        #self.train = 1500#0
        #self.val = 200#0
        #self.test = 300#0
        
        # https://energyflow.network/docs/datasets/#quark-and-gluon-jets
        # X : a three-dimensional numpy array of jets:
        #     list of jets with list of particles for each jet, with (pt,y,phi,pid) values for each particle
        # y : a numpy array of quark/gluon jet labels (quark=1 and gluon=0).
        # The jets are padded with zero-particles in order to make a contiguous array.
        print()
        print('Loading qg dataset:')
        self.X, self.y = energyflow.datasets.qg_jets.load(self.train + self.val + self.test)
        print('(n_jets, n_particles per jet, n_variables): {}'.format(self.X.shape))
        print()

        # Next, we will transform these into fastjet::PseudoJet objects.
        # This allows us use the fastjet contrib to compute Nsubjetiness, and in general it
        # will be needed in order to perform jet finding on an event (in data or MC).

        # Translate 3D numpy array (100,000 x 556 particles x 4 vars) into a dataframe
        # Define a unique index for each jet
        columns = ['pt', 'y', 'phi', 'pid']
        df_particles = pd.DataFrame(self.X.reshape(-1, 4), columns=columns)
        df_particles.index = np.repeat(np.arange(self.X.shape[0]), self.X.shape[1]) + 1
        df_particles.index.name = 'jet_id'
        
        # (i) Group the particle dataframe by jet id
        #     df_particles_grouped is a DataFrameGroupBy object with one particle dataframe per jet
        df_particles_grouped = df_particles.groupby('jet_id')
        
        # (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet::PseudoJets
        # NOTE: for now we neglect the mass -- and assume y=eta
        # TO DO: Add y to https://github.com/matplo/heppy/blob/master/cpptools/src/fjext/fjtools.cxx
        # TO DO: Add mass vector using pdg
        print('Converting particle dataframe to fastjet::PseudoJets...')
        self.df_fjparticles = df_particles_grouped.apply(self.get_fjparticles)
        print('Done.')
        print()
        
        # Define Nsubjettiness observables to compute
        # The K-body phase space is (3K-4)-dimensional
        # Choose K and create list of N-subjettiness observables: number of axes and beta values
        
        self.K = 24 # Note: N-subjettiness plot only works up to K=6 (1810.05165 used K=24)
        self.N_list = []
        self.beta_list = []

        for i in range(self.K-2):
            self.N_list += [i+1] * 3
            self.beta_list += [0.5,1,2]

        self.N_list += [self.K-1] * 2  
        self.beta_list += [1,2]
        
        # Construct dictionary to store all jet quantities of interest
        self.jet_variables = {}
        for i,N in enumerate(self.N_list):
            beta = self.beta_list[i]
            self.jet_variables['n_subjettiness_N{}_beta{}'.format(N,beta)] = []
        
        print(self)
        print()

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def process_qg(self):
    
        # Loop over jets and compute quantities of interest
        # Fill each of the jet_variables into a list
        fj.ClusterSequence.print_banner()
        print('Finding jets and computing N-subjettiness...')
        result = [self.analyze_event(fj_particles) for fj_particles in self.df_fjparticles]
        
        # Transform the dictionary of lists into a dictionary of numpy arrays
        self.jet_variables_numpy = self.transform_to_numpy(self.jet_variables)
        n_subjettiness = self.jet_variables_numpy['n_subjettiness_N{}_beta{}'.format(self.N_list[0], self.beta_list[0])]
        print('Done. Number of clustered jets: {}'.format(len(n_subjettiness)))
        print()
        
        # Plot jet quantities
        if self.K <= 6:
            self.plot_nsubjettiness()

        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        
        # Reformat output for ML algorithms (array with 1 array per jet which contain all N-subjettiness values)
        X_Nsub = np.array([list(self.jet_variables_numpy.values())])[0].T
        
        # Split data into train and test sets (extend to include a validation set as well. See: data_split?)
        test_frac = 0.2
        (X_Nsub_train, X_Nsub_test, y_Nsub_train, y_Nsub_test) = data_split(X_Nsub, self.y, val=0, test=test_frac)

        # Fit ML model -- 1. SGDClassifier
        print('Training SGDClassifier')
        sgd_clf = SGDClassifier(max_iter=1000, tol=1e-3, random_state=42)
        sgd_clf.fit(X_Nsub_train, y_Nsub_train)
        
        # Use cross validation predict (here split training set) and compute the confusion matrix from the predictions
        y_Nsub_crossval_SGD = cross_val_predict(sgd_clf, X_Nsub_train, y_Nsub_train, cv=3,method="decision_function")
        #confusion_SGD = confusion_matrix(y_Nsub_train, y_Nsub_crossval_SGD)
        #print('Confusion matrix for SGD Classifier (test set): \n {}'.format(confusion_SGD))        

        # Get predictions for the test set .. actually don't need this when using cross_val_predict (?)
        # preds_Nsub_SGD = sgd_clf.predict(X_Nsub_test)
        
        # Get AUC from training process
        Nsub_auc_SGD = roc_auc_score(y_Nsub_train,y_Nsub_crossval_SGD)
        print('SGDClassifier: AUC = {} (cross validation)'.format(Nsub_auc_SGD))
        
        # Compute ROC curve: the roc_curve() function expects labels and scores
        Nsub_fpr_SGD, Nsub_tpr_SGD, thresholds = roc_curve(y_Nsub_train,y_Nsub_crossval_SGD)
        
        # Plot ROC curve for SGDClassifier
        # self.plot_roc_curve(Nsub_fpr_SGD,Nsub_tpr_SGD)
        
        # Check number of threhsolds used for ROC curve
        # n_thresholds = len(np.unique(y_Nsub_scores_SGD)) + 1

        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        
        # Fit ML model -- 2. Random Forest Classifier
        forest_clf = RandomForestClassifier(random_state=42)
        y_Nsub_probas_forest = cross_val_predict(forest_clf, X_Nsub_train, y_Nsub_train, cv=3,method="predict_proba")
        
        # The output here are class probabilities. We us use the positive class's probability for the ROC curve
        y_Nsub_scores_forest = y_Nsub_probas_forest[:,1]
        
        print(y_Nsub_scores_forest)
        
        # Compute AUC & ROC curve
        Nsub_auc_RFC = roc_auc_score(y_Nsub_train,y_Nsub_scores_forest)
        print('Random Forest Classifier: AUC = {} (cross validation)'.format(Nsub_auc_RFC))
        Nsub_fpr_forest, Nsub_tpr_forest, thresholds_forest = roc_curve(y_Nsub_train,y_Nsub_scores_forest)
        
        # Plot ROC curve
        self.plot_roc_curve(Nsub_fpr_SGD,Nsub_tpr_SGD,Nsub_fpr_forest, Nsub_tpr_forest,"SGD_Nsub","RF_Nsub")

        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        
        # Fit ML model -- 3. Dense Neural network with Keras
        # input_shape expects shape of an instance (not including batch size)
        DNN = keras.models.Sequential()
        DNN.add(keras.layers.Flatten(input_shape=[X_Nsub_train.shape[1]]))
        DNN.add(keras.layers.Dense(300,activation='relu'))
        DNN.add(keras.layers.Dense(300,activation='relu'))
        DNN.add(keras.layers.Dense(100,activation='relu'))
        DNN.add(keras.layers.Dense(1,activation='sigmoid')) # softmax? # Last layer has to be 1 or 2 for binary classification?

        # Print DNN summary
        DNN.summary()
        
        # Compile DNN
        opt = keras.optimizers.Adam(learning_rate=0.001) # lr = 0.001 (cf 1810.05165)
        DNN.compile(loss="binary_crossentropy",          # Loss function - use categorical_crossentropy instead ?
                    optimizer=opt,                       # For Stochastic gradient descent use: "sgd"
                    metrics=["accuracy"])                # Measure accuracy during training

        # Train DNN - need validation set here (use test set for now)
        DNN.fit(X_Nsub_train,y_Nsub_train, epochs=39, validation_data=(X_Nsub_test,y_Nsub_test))
        
        # Get predictions for validation data set
        y_Nsub_test_preds_DNN = DNN.predict(X_Nsub_test).reshape(-1)
        
        # Get AUC
        Nsub_auc_DNN = roc_auc_score(y_Nsub_test,y_Nsub_test_preds_DNN)
        print('Dense Neural Network: AUC = {} (validation set)'.format(Nsub_auc_DNN))        
        
        # Get ROC curve results
        Nsub_fpr_DNN, Nsub_tpr_DNN, thresholds = roc_curve(y_Nsub_test,y_Nsub_test_preds_DNN)
        
        # Plot ROC curve
        self.plot_roc_curve(Nsub_fpr_SGD,Nsub_tpr_SGD,Nsub_fpr_DNN, Nsub_tpr_DNN,"SGD_Nsub","DNN_Nsub")
        
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        
        # Fit ML model -- 4. Deep Set/Particle Flow Networks

        # network architecture parameters
        Phi_sizes, F_sizes = (100, 100, 256), (100, 100, 100)
        
        # network training parameters
        num_epoch = 34
        batch_size = 500
        
        # Use PID information
        use_pids = True
        
        # convert labels to categorical
        Y_PFN = to_categorical(self.y, num_classes=2)
        
        # preprocess by centering jets and normalizing pts
        X_PFN = self.X
        for x_PFN in X_PFN:
            mask = x_PFN[:,0] > 0
            yphi_avg = np.average(x_PFN[mask,1:3], weights=x_PFN[mask,0], axis=0)
            x_PFN[mask,1:3] -= yphi_avg
            x_PFN[mask,0] /= x_PFN[:,0].sum()
        
        # handle particle id channel [?? ... remap_pids is not working] 
        #if use_pids:
        #    remap_pids(X_PFN, pid_i=3)
        #else:
        X_PFN = X_PFN[:,:,:3]

        # Split data into train, val and test sets
        (X_PFN_train, X_PFN_val, X_PFN_test,Y_PFN_train, Y_PFN_val, Y_PFN_test) = data_split(X_PFN, Y_PFN, 
                                                                                             val=self.val, test=self.test)
        # build architecture
        pfn = PFN(input_dim=X_PFN.shape[-1], Phi_sizes=Phi_sizes, F_sizes=F_sizes)

        # train model
        pfn.fit(X_PFN_train, Y_PFN_train,
                epochs=num_epoch,
                batch_size=batch_size,
                validation_data=(X_PFN_val, Y_PFN_val),
                verbose=1)
        
        # get predictions on test data
        preds_PFN = pfn.predict(X_PFN_test, batch_size=1000)

        # Get AUC and ROC curve + make plot
        auc_PFN = roc_auc_score(Y_PFN_test[:,1], preds_PFN[:,1])
        print('Particle Flow Networks/Deep Sets: AUC = {} (test set)'.format(auc_PFN))   
        
        fpr_PFN, tpr_PFN, threshs = roc_curve(Y_PFN_test[:,1], preds_PFN[:,1])
        self.plot_roc_curve(Nsub_fpr_SGD,Nsub_tpr_SGD,fpr_PFN, tpr_PFN,"SGD_Nsub","PFN_woPID")
        
    #---------------------------------------------------------------
    # Process an event (in this case, just a single jet per event)
    #---------------------------------------------------------------
    def analyze_event(self, fj_particles):
        
        # Cluster each jet with R=infinity
        jetR = fj.JetDefinition.max_allowable_R
        jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
        cs = fj.ClusterSequence(fj_particles, jet_def)
        jet = fj.sorted_by_pt(cs.inclusive_jets())[0]

        # Compute N-subjettiness
        axis_definition = fjcontrib.KT_Axes()
        for i,N in enumerate(self.N_list):
            beta = self.beta_list[i]
        
            measure_definition = fjcontrib.UnnormalizedMeasure(beta)
            n_subjettiness_calculator = fjcontrib.Nsubjettiness(N, axis_definition, measure_definition)
            n_subjettiness = n_subjettiness_calculator.result(jet)/jet.pt()
            self.jet_variables['n_subjettiness_N{}_beta{}'.format(N, beta)].append(n_subjettiness)
        
        # Compute four-vector...
        # ...
            
    #---------------------------------------------------------------
    # Plot N-subjettiness
    #---------------------------------------------------------------
    def plot_nsubjettiness(self):
    
        linestyles = ['-', '--', ':', '-.', '-']
    
        bins = np.linspace(0, 0.7, 100)
        for i,N in enumerate(self.N_list):
            beta = self.beta_list[i]
            
            plt.hist(self.jet_variables_numpy['n_subjettiness_N{}_beta{}'.format(N,beta)],
                     bins,
                     histtype='stepfilled',
                     label = r'$N={}, \beta={}$'.format(N, beta),
                     linewidth=2,
                     linestyle=linestyles[N-1],
                     alpha=0.5)
                     
        plt.xlabel(r'$\tau_{N}^{\beta}$', fontsize=14)
        plt.yscale('log')
        legend = plt.legend(loc='best', fontsize=10, frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'Nsubjettiness.pdf'))
        plt.close()


    #--------------------------------------------------------------- 
    # Plot ROC curve                                                 
    #--------------------------------------------------------------- 
    def plot_roc_curve(self, fpr1, tpr1, fpr2, tpr2, label1=None, label2=None):
        plt.plot(fpr1, tpr1, "b:", label=label1)
        plt.plot(fpr2, tpr2, linewidth=2, label=label2)
        plt.plot([0, 1], [0, 1], 'k--') # dashed diagonal
        plt.axis([0, 1, 0, 1])                                   
        plt.xlabel('False Positive Rate', fontsize=16)
        plt.ylabel('True Positive Rate', fontsize=16) 
        plt.grid(True)    
        plt.legend(loc="lower right")        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'ROC_{}_{}.pdf'.format(label1,label2)))
        plt.close()
          
    #---------------------------------------------------------------
    # Transform dictionary of lists into a dictionary of numpy arrays
    #---------------------------------------------------------------
    def transform_to_numpy(self, jet_variables_list):

        jet_variables_numpy = {}
        for key,val in jet_variables_list.items():
            jet_variables_numpy[key] = np.array(val)
                    
        return jet_variables_numpy
            
    #---------------------------------------------------------------
    # Cluster jets
    #---------------------------------------------------------------
    def get_fjparticles(self, df_particles_grouped):
                                                 
        user_index_offset = 0
        return fjext.vectorize_pt_eta_phi(df_particles_grouped['pt'].values,
                                          df_particles_grouped['y'].values,
                                          df_particles_grouped['phi'].values,
                                          user_index_offset)

##################################################################
if __name__ == '__main__':

    analysis = ProcessQG()
    analysis.process_qg()
