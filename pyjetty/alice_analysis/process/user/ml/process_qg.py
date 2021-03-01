#!/usr/bin/env python3

"""
Example class to read quark-gluon dataset
"""

# Data analysis and plotting
import pandas as pd
import numpy as np

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext

# Energy flow package
import energyflow

# Base class
from pyjetty.alice_analysis.process.base import common_base

################################################################
class ProcessQG(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
        super(common_base.CommonBase, self).__init__(**kwargs)
        
        # Load labeled data
        self.n_train = 75000
        self.n_validation = 10000
        self.n_test = 15000
        
        # https://energyflow.network/docs/datasets/#quark-and-gluon-jets
        # X : a three-dimensional numpy array of jets:
        #     list of jets with list of particles for each jet, with (pt,y,phi,pid) values for each particle
        # y : a numpy array of quark/gluon jet labels (quark=1 and gluon=0).
        # The jets are padded with zero-particles in order to make a contiguous array.
        print()
        print('Loading qg dataset:')
        X, self.y = energyflow.datasets.qg_jets.load(self.n_train + self.n_validation + self.n_test)
        print('(n_jets, n_particles per jet, n_variables): {}'.format(X.shape))
        print()

        # Next, we will transform these into fastjet::PseudoJet objects.
        # This allows us use the fastjet contrib to compute Nsubjetiness, and in general it
        # will be needed in order to perform jet finding on an event (in data or MC).

        # Translate 3D numpy array (100,000 x 556 particles x 4 vars) into a dataframe
        # Define a unique index for each jet
        columns = ['pt', 'y', 'phi', 'pid']
        df_particles = pd.DataFrame(X.reshape(-1, 4), columns=columns)
        df_particles.index = np.repeat(np.arange(X.shape[0]), X.shape[1]) + 1
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
        
        print(self)
        print()

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def process_qg(self):
    
        # Cluster each jet with R=infinity
        fj.ClusterSequence.print_banner()
        print('Finding jets...')
        jets = [self.find_jets(fj_particles[:1]) for fj_particles in self.df_fjparticles]
        print('Done. Number of clustered jets: {}'.format(len(jets)))
        print()
   
        # Compute N-subjetiness into numpy array
        N = 1
        beta = 1.
        # axis_definition = fjcontrib.Nsubjettiness.NjettinessDefinition.KT_Axes
        axis_definition = fjcontrib.KT_Axes()
        # measure_definition = fjcontrib.Nsubjettiness.NjettinessDefinition.UnnormalizedMeasure(beta)
        measure_definition = fjcontrib.UnnormalizedMeasure(beta)
        n_subjetiness_calculator = fjcontrib.Nsubjettiness(N, axis_definition, measure_definition)
        n_subjetiness = np.array([n_subjetiness_calculator.result(j) for jet in jets for j in jet])
        print(n_subjetiness)
        
        # Compute four-vector numpy array
        # ...
        
        # Fit ML model
        # ...
        
    #---------------------------------------------------------------
    # Cluster jets
    #---------------------------------------------------------------
    def find_jets(self, fj_particles):
        
        jetR = fj.JetDefinition.max_allowable_R
        jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
        jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(1.0)
        cs = fj.ClusterSequence(fj_particles, jet_def)
        jets = fj.sorted_by_pt(cs.inclusive_jets())
        jets_selected = jet_selector(jets)
        return jets_selected
        
    #---------------------------------------------------------------
    # Cluster jets
    #---------------------------------------------------------------
    def get_fjparticles(self, df_particles_grouped):
                                                 
        return fjext.vectorize_pt_eta_phi(df_particles_grouped['pt'].values,
                                          df_particles_grouped['y'].values,
                                          df_particles_grouped['phi'].values,
                                          user_index_offset=0)

##################################################################
if __name__ == '__main__':

    analysis = ProcessQG()
    analysis.process_qg()
