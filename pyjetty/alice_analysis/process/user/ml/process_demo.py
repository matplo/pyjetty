#!/usr/bin/env python3

"""
Example class to read a ROOT TTree of track information

Embed PYTHIA MC into Pb-Pb data
"""

# Data analysis and plotting
import pandas
import numpy as np

# Fastjet via python (from external library heppy)
import fastjet as fj

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_io_emb
from pyjetty.alice_analysis.process.base import common_base

################################################################
class ProcessDemo(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
        super(common_base.CommonBase, self).__init__(**kwargs)
        
        # Some example files on hiccup, from the latest datasets
        
        # PYTHIA anchored to PbPb period
        self.input_file_Pythia = '/rstorage/alice/data/LHC20g4/568/LHC20g4/13/296191/0003/AnalysisResults.root'
        
        # List of 0-10% Pb-Pb data files
        self.emb_file_list = '/rstorage/alice/data/LHC18qr/570/files.txt'

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def process_demo(self):
    
        # Use IO helper class to convert detector-level ROOT TTree into
        # a SeriesGroupBy object of fastjet particles per event
        tree_dir = 'PWGHF_TreeCreator'
        event_tree_name='tree_event_char'
        io_det = process_io.ProcessIO(input_file=self.input_file_Pythia, tree_dir=tree_dir,
                                      track_tree_name='tree_Particle', event_tree_name=event_tree_name,
                                      is_pp=True, min_cent=0., max_cent=10., use_ev_id_ext=False)
        df_fjparticles_det = io_det.load_data(reject_tracks_fraction=0.)
        self.nEvents_det = len(df_fjparticles_det.index)
        self.nTracks_det = len(io_det.track_df.index)
        print('nEvents in detector-level PYTHIA: {}'.format(self.nEvents_det))
        print('nTracks in detector-level PYTHIA: {}'.format(self.nTracks_det))
        
        # ------------------------------------------------------------------------

        # Use IO helper class to convert truth-level ROOT TTree into
        # a SeriesGroupBy object of fastjet particles per event
        io_truth = process_io.ProcessIO(input_file=self.input_file_Pythia, tree_dir=tree_dir,
                                        track_tree_name='tree_Particle_gen', event_tree_name=event_tree_name,
                                        is_pp=True, use_ev_id_ext=False)
        df_fjparticles_truth = io_truth.load_data()
        self.nEvents_truth = len(df_fjparticles_truth.index)
        self.nTracks_truth = len(io_truth.track_df.index)
        print('nEvents in truth-level PYTHIA: {}'.format(self.nEvents_truth))
        print('nTracks in truth-level PYTHIA: {}'.format(self.nTracks_truth))

        # ------------------------------------------------------------------------

        # Now merge the two SeriesGroupBy to create a groupby df with [ev_id, run_number, fj_1, fj_2]
        # (Need a structure such that we can iterate event-by-event through both fj_1, fj_2 simultaneously)
        print('Merge det-level and truth-level into a single dataframe grouped by event...')
        self.df_fjparticles = pandas.concat([df_fjparticles_det, df_fjparticles_truth], axis=1)
        self.df_fjparticles.columns = ['fj_particles_det', 'fj_particles_truth']
        print(self.df_fjparticles)
        # ------------------------------------------------------------------------
        
        # Set up the Pb-Pb embedding object
        self.process_io_emb = process_io_emb.ProcessIO_Emb(self.emb_file_list, track_tree_name='tree_Particle',
                                                           min_cent=0., max_cent=10.)
                       
        # ------------------------------------------------------------------------
                
        self.analyze_events()

    #---------------------------------------------------------------
    # Main function to loop through and analyze events
    #---------------------------------------------------------------
    def analyze_events(self):
        
        self.n_event = 0

        # Then can use list comprehension to iterate over the groupby and do jet-finding
        # simultaneously for fj_1 and fj_2 per event, so that I can match jets
        result = [self.analyze_event(fj_particles_det, fj_particles_truth) for fj_particles_det, fj_particles_truth in zip(self.df_fjparticles['fj_particles_det'], self.df_fjparticles['fj_particles_truth'])]
        
    #---------------------------------------------------------------
    # Analyze jets of a given event.
    # fj_particles is the list of fastjet pseudojets for a single fixed event.
    #---------------------------------------------------------------
    def analyze_event(self, fj_particles_det, fj_particles_truth):
          
        if self.n_event !=0:
            return

        # Check that the entries exist appropriately
        if type(fj_particles_det) != fj.vectorPJ or type(fj_particles_truth) != fj.vectorPJ:
          print('fj_particles type mismatch -- skipping event')
          return
          
        # If Pb-Pb, construct embedded event
        print('Load Pb-Pb event...')
        fj_particles_combined = self.process_io_emb.load_event()
            
        # Form the combined det-level event -- std::vector of fj.PseudoJet tracks
        # The pp-det tracks are each stored with a unique user_index >= 0
        #   (same index in fj_particles_combined and fj_particles_det)
        # The Pb-Pb tracks are each stored with a unique user_index < 0
        print('embed PYTHIA detector-level event into Pb-Pb event')
        [fj_particles_combined.push_back(p) for p in fj_particles_det]

        print('user_index of PYTHIA truth event:')
        print([p.user_index() for p in fj_particles_truth])

        print('user_index of combined detector-level event:')
        print([p.user_index() for p in fj_particles_combined])

        self.n_event +=1

##################################################################
if __name__ == '__main__':

    analysis = ProcessDemo()
    analysis.process_demo()
