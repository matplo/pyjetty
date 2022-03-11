#!/usr/bin/env python3

"""
Class to loop over tracks generated from JEWEL simulations, do jet finding,
and potentially do recoil subtraction.
author: Rey Cruz-Torres (reynier@lbl.gov)
several functions based on code by James Mulligan
"""

import uproot
import pandas as pd
import numpy as np
import time
import math
import ROOT
import yaml
import os
from array import *

# Fastjet via python (from external library fjpydev)
import fastjet as fj
import fjext
import fjcontrib

from pyjetty.alice_analysis.process.base import process_utils
from pyjetty.mputils import CEventSubtractor

# Load pyjetty ROOT utils
ROOT.gSystem.Load('libpyjetty_rutil')

################################################################
################################################################
################################################################
class CurvesFromJewelTracks():
    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, input_file='', config_file='', output_dir='', **kwargs):
        self.config_file = config_file
        self.input_file = input_file
        self.output_dir = output_dir
        self.track_tree = 'tree_Particle_gen'
        self.dirname = 'PWGHF_TreeCreator'

        self.histutils = ROOT.RUtil.HistUtils()

        # Initialize utils class
        self.utils = process_utils.ProcessUtils()

        # Read config file
        with open(self.config_file, 'r') as stream:
            config = yaml.safe_load(stream)

        self.track_df = None

        self.jetR_list = config['jetR']
        self.reclustering_algorithm = fj.cambridge_algorithm

        self.recoils_off = False
        if 'recoils_off' in config:
            self.recoils_off = config['recoils_off']

        #----------------------------------------------
        # Loading recoil subtraction parameters
        self.thermal_subtraction_method = None
        self.gridsizes = None

        if 'thermal_subtraction_method' in config:
            self.thermal_subtraction_method = config['thermal_subtraction_method'].lower()

        if not self.thermal_subtraction_method:
            print('Will not do recoil subtraction')
        elif 'gridsub' in self.thermal_subtraction_method:
            if 'gridsizes' in config:
                self.gridsizes = config['gridsizes']
            else:
                print('User requested gridsub subtraction method, but no gridsize was provided. Bailing out!')
                exit()

        self.grid_dict = {}
        self.cell_phi = {}
        self.cell_eta = {}
        self.populated_cells = []
        self.populated_cells_w_constit = []
        self.total_thermal_momentum = 0
        self.unsubtracted_thermal_momentum = 0

        self.run_diagnostics = None
        if 'run_diagnostics' in config:
            self.run_diagnostics = config['run_diagnostics']

        #----------------------------------------------
        # If specified in the config file, randomly reject this fraction of thermals
        self.thermal_rejection_fraction = 0.0
        if 'thermal_rejection_fraction' in config:
            self.thermal_rejection_fraction = config['thermal_rejection_fraction']

        #----------------------------------------------
        # If constituent subtractor is present, initialize it
        self.constituent_subtractor = None
        if 'constituent_subtractor' in config:

            print('Constituent subtractor is enabled.')
            constituent_subtractor = config['constituent_subtractor']

            max_distance = constituent_subtractor['max_distance']
            alpha = constituent_subtractor['alpha']
            max_eta = constituent_subtractor['max_eta']
            bge_rho_grid_size = constituent_subtractor['bge_rho_grid_size']
            max_pt_correct = constituent_subtractor['max_pt_correct']
            ghost_area = constituent_subtractor['ghost_area']

            self.constituent_subtractor = CEventSubtractor(max_distance=max_distance, 
                                                           alpha=alpha, 
                                                           max_eta=max_eta, 
                                                           bge_rho_grid_size=bge_rho_grid_size, 
                                                           max_pt_correct=max_pt_correct, 
                                                           ghost_area=ghost_area, 
                                                           distance_type=fjcontrib.ConstituentSubtractor.deltaR)

        else:
            print('Constituent subtractor is disabled.')

        #----------------------------------------------
        # Create dictionaries to store grooming settings and observable settings for each observable
        # Each dictionary entry stores a list of subconfiguration parameters
        # The observable list stores the observable setting, e.g. subjetR
        # The grooming list stores a list of grooming settings {'sd': [zcut, beta]} or {'dg': [a]}
        self.observable_list = config['process_observables']
        self.obs_settings = {}
        self.obs_grooming_settings = {}
        self.obs_names = {}
        for observable in self.observable_list:

            obs_config_dict = config[observable]
            obs_config_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]

            if "common_settings" in list(obs_config_dict.keys()) and \
              "xtitle" in list(obs_config_dict["common_settings"].keys()):

              self.obs_names[observable] = obs_config_dict["common_settings"]["xtitle"]

            obs_subconfig_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
            self.obs_settings[observable] = self.utils.obs_settings(observable, obs_config_dict, obs_subconfig_list)
            self.obs_grooming_settings[observable] = self.utils.grooming_settings(obs_config_dict)

            self.initialize_histos()

    #---------------------------------------------------------------
    # Initialize some histograms
    #---------------------------------------------------------------
    def initialize_histos(self):
        # Initialize base histograms
        self.hNevents = ROOT.TH1F('hNevents', 'hNevents', 2, -0.5, 1.5)

        if not self.thermal_subtraction_method:
            print('no histograms will be initialized for this setting in initialize_histos')

        elif 'gridsub' in self.thermal_subtraction_method:
            for jetR in self.jetR_list:
                for gridsize in self.gridsizes:
                    name = 'h_thermal_fraction_not_subtracted_v_pT_R{}_gridsize{}'.format(jetR,gridsize)
                    h = ROOT.TH2F(name,name,100,0,200,100,0,1.01)
                    setattr(self,name,h)

        elif '4momsub' in self.thermal_subtraction_method or \
             'negative_recombiner' in self.thermal_subtraction_method:

            print('no histograms will be initialized for this setting in initialize_histos')

        else:
            raise NotImplementedError(
                "Subtraction method %s not implemented" % self.thermal_subtraction_method)

    #---------------------------------------------------------------
    # Main function
    #---------------------------------------------------------------
    def process_gen(self):
        self.start_time = time.time()
        print('**********************')
        print('Entering main function')
        print('**********************')
        print('--- {} seconds ---'.format(time.time() - self.start_time))

        print('Loading the data')
        self.load_data(self.input_file)
        self.hNevents.Fill(1, self.nEvents)
        print('--- {} seconds ---'.format(time.time() - self.start_time))

        print('Initializing histograms')
        self.initialize_user_output_objects()
        if self.run_diagnostics:
            self.initialize_user_output_objects('_diagnostics')
        print('--- {} seconds ---'.format(time.time() - self.start_time))

        print('Doing the analysis and filling histograms')
        self.analyze_events()
        print('--- {} seconds ---'.format(time.time() - self.start_time))

        print('Saving histograms to output file')
        self.save_output_objects()
        print('--- {} seconds ---'.format(time.time() - self.start_time))

    #---------------------------------------------------------------
    # Load data from root files
    #---------------------------------------------------------------
    def load_data(self,input_fname):
        print('Loading data from file:',input_fname)

        # Load trees from input root file
        track_tree = uproot.open(input_fname)[self.dirname+'/'+self.track_tree]

        # Store the trees in pandas dataframes
        self.track_df = track_tree.arrays(library='pd')

        # If we are using the GridSub method, then we won't be
        # using the dummy particles
        if not self.thermal_subtraction_method:
            track_criteria = 'ParticlePt > 1e-5 and ParticleEta < 0.9 and ParticleEta > -0.9'
            if not self.recoils_off:
                track_criteria += ' and Status != 3'
            # Not used by default, unless constituent subtraction is activated or
            #     user accesses thermal particles in user function
            bck_criteria = 'ParticlePt > 1e-5 and ParticleEta < 0.9 and ParticleEta > -0.9 and Status == 3'
        elif 'gridsub' in self.thermal_subtraction_method:
            track_criteria = 'ParticlePt > 1e-5 and ParticleEta < 0.9 and ParticleEta > -0.9 and Status != 3'
            bck_criteria = 'ParticlePt > 1e-5 and ParticleEta < 0.9 and ParticleEta > -0.9 and Status == 3'
        elif '4momsub' in self.thermal_subtraction_method:
            track_criteria = 'ParticleEta < 0.9 and ParticleEta > -0.9 and Status != 3'
            bck_criteria = 'ParticleEta < 0.9 and ParticleEta > -0.9 and Status == 3'
        elif 'negative_recombiner' in self.thermal_subtraction_method:
            track_criteria = 'ParticlePt > 1e-5 and ParticleEta < 0.9 and ParticleEta > -0.9 and Status != 3'
            bck_criteria = 'ParticlePt > 1e-5 and ParticleEta < 0.9 and ParticleEta > -0.9 and Status == 3'

        self.jet_df = self.track_df.query(track_criteria)
        self.bck_df = None if self.recoils_off else self.track_df.query(bck_criteria)
        self.ev_idx = []

        print('Transforming track dataframe into SeriesGroupBy object of fastjet particles per event.')
        self.df_fjparticles = self.group_fjparticles(self.jet_df)
        self.nEvents = len(self.df_fjparticles.index)

    #---------------------------------------------------------------
    # Transform the track dataframe into a SeriesGroupBy object of fastjet particles per event.
    # This function was adapted from alice_analysis/process/base/process_io.py
    #---------------------------------------------------------------
    def group_fjparticles(self,df):
        print("Transform the track dataframe into a series object of fastjet particles per event...")
        # (i) Group dataframe by event  track_df_grouped is a DataFrameGroupBy object with one track dataframe per event
        df_grouped = None
        df_grouped = df.groupby(['run_number', 'ev_id'])

        # (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet particles
        df_fjparticles = None
        df_fjparticles = df_grouped.apply(self.get_fjparticles)

        return df_fjparticles

    #---------------------------------------------------------------
    # Return fastjet:PseudoJets from a given track dataframe
    #---------------------------------------------------------------
    def get_fjparticles(self, df_tracks):
        m=0.1396
        user_index_offset = 0
        m_array = np.full((df_tracks['ParticlePt'].values.size), m)

        # Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of pt,eta,phi
        fj_particles = fjext.vectorize_pt_eta_phi_m(
          df_tracks['ParticlePt'].values, df_tracks['ParticleEta'].values,
          df_tracks['ParticlePhi'].values, m_array, user_index_offset)

        self.ev_idx.append(df_tracks['ev_id'].values[0])

        return fj_particles

    #---------------------------------------------------------------
    # Main function to loop through and analyze events
    #---------------------------------------------------------------
    def analyze_events(self):

        print('Finding jets...')
        fj.ClusterSequence.print_banner()
        print()
        self.event_number = 0

        # Do jet-finding and fill histograms
        if 'gridsub' in self.thermal_subtraction_method:
            print('Using GridSub method to remove recoils')
            for gridsize in self.gridsizes:
                print('Now doing the analysis for a gridsize =',gridsize)
                self.populated_cells = []
                self.populated_cells_w_constit = []
                self.grid_dict[gridsize] = self.create_grids(gridsize)
                self.event_number = 0
                for fj_particles in self.df_fjparticles:
                    self.analyze_event(fj_particles,gridsize)
                    if self.run_diagnostics:
                       self.event_number = 0
                       self.analyze_event(fj_particles, gridsize, True)
                del self.grid_dict[gridsize]
                print('--- {} seconds ---'.format(time.time() - self.start_time))

        else:
            if not self.thermal_subtraction_method:
                print('No recoil removal method selected')
            else:
                print('Using %s method to remove recoils' % self.thermal_subtraction_method)
            for fj_particles in self.df_fjparticles:
                self.analyze_event(fj_particles)
            print('--- {} seconds ---'.format(time.time() - self.start_time))

    #---------------------------------------------------------------
    # Analyze jets of a given event.
    # fj_particles is the list of fastjet pseudojets for a single fixed event.
    #---------------------------------------------------------------
    def analyze_event(self, fj_particles, gridsize=None, diagnostic=False):

        if self.event_number > self.nEvents:
            return

        self.event_bck_df = None
        self.event_bck_df = self.bck_df[self.bck_df['ev_id']==self.ev_idx[self.event_number]]
        #print(self.event_number)

        self.event_number += 1

        # Convert df of pt/eta/phi of thermal particles to list of
        #     fastjet particles, for convenience
        thermal_particles = []
        thermal_particles_selected = []
        if len(self.event_bck_df):
            thermal_particles = self.get_fjparticles(self.event_bck_df)
            # Drop specified fraction of thermal particles
            # Loop manually since the wrapped functions are a bit funky
            for i, p in enumerate(thermal_particles):
                if np.random.uniform() >= self.thermal_rejection_fraction:
                    thermal_particles_selected.append(p)
            #print(f'n_thermals before: {len(thermal_particles)}')
            #print(f'n_thermals after: {len(thermal_particles_selected)}')

        if len(fj_particles) > 1:
            if np.abs(fj_particles[0].eta() - fj_particles[1].eta()) < 1e-10:
                if fj_particles[0].pt() < 1e-5:
                    print('WARNING: Duplicate DUMMY particles may be present in event',
                          self.event_number)
                else:
                    print('WARNING: Duplicate JET particles may be present in event',
                          self.event_number)
                    #[print(p.eta(),p.pt()) for p in fj_particles]

        # If constituent subtraction is enabled, perform subtraction on the event
        if self.constituent_subtractor:

            # Determine rho from thermal particles
            self.constituent_subtractor.bge_rho.set_particles(thermal_particles_selected)

            # Perform subtraction over full event (jet+recoil)
            fj_particles = self.constituent_subtractor.subtractor.subtract_event(fj_particles)

            #rho = self.constituent_subtractor.bge_rho.rho()
            #print(f'rho: {rho}')
            #print()

        # Loop through jetR, and process event for each R
        for jetR in self.jetR_list:
            # Set jet definition and a jet selector
            jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
            jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)

            particles = fj_particles
            # For negative pT treatment, add thermals and negative recombiner
            if "negative_recombiner" in self.thermal_subtraction_method:
                for part in thermal_particles_selected:
                    part.set_user_index(-1)
                    particles.push_back(part)
                recombiner = fjext.NegativeEnergyRecombiner(-1)
                jet_def.set_recombiner(recombiner)

            # Do jet finding
            cs = fj.ClusterSequence(particles, jet_def)
            jets = fj.sorted_by_pt(cs.inclusive_jets())
            jets_selected = jet_selector(jets)

            # If no jets were selected, move on to next event
            if len(jets_selected) == 0:
                continue

            #----------------------------------------------
            fj_subtracted_constituents = None
            subtracted_jets = []

            if not self.thermal_subtraction_method:
                subtracted_jets = jets_selected

            elif "negative_recombiner" in self.thermal_subtraction_method:
                subtracted_jets = jets_selected

            elif '4momsub' in self.thermal_subtraction_method:
                '''
                ------------------------------------------------
                Thermal subtraction using 4MomSub method
                ------------------------------------------------
                1. Cluster the initial jet collection from the final state particles (including dummies).
                2. Compile a list of the thermal momenta (particles in the HepMC event record with status code 3).
                For each jet, get the list of thermal momenta that have DeltaR < 1x10-5 with one of the jet constituents, i.e a dummy particle.
                4. Sum up the four-momenta of the matched thermal momenta. This constitutes the background.
                5. For each jet subtract the background four-momentum from the jet's four momentum, this provides the corrected jet collection.
                6. Calculate jet observables from corrected jet four-momenta.
                '''
                for jet in jets_selected:
                    fj_subtracted_constituents = None
                    fj_subtracted_constituents = self.subtract_thermal_4momsub(jet,jetR)
                    jet_def_largerR = fj.JetDefinition(fj.antikt_algorithm, 2.*jetR)
                    cs_largerR = fj.ClusterSequence(fj_subtracted_constituents, jet_def_largerR)
                    subtracted_jets = fj.sorted_by_pt(cs_largerR.inclusive_jets())
                    if len(subtracted_jets) > 1:
                        print('WARNING: Got more than one subtracted jet out of one input jet')

            elif 'gridsub' in self.thermal_subtraction_method:
                '''
                ------------------------------------------------
                Thermal subtraction using GridSub1 method
                ------------------------------------------------
                1. Cluster the initial jet collection from the final state particles.
                2. Compile a list of the thermal momenta (particles in the HepMC event record with status code 3).
                3. Define the grid resolution and place grid over jets.
                4. Inside each cell sum the jet constituents' four-momenta and subtract the thermal
                four-momenta that fall into the cell (note: no matching is required, thermal four momenta
                with distance DeltaR < R from the jet axis are considered), providing a single four momentum for each cell.
                5. In case a cell contains more thermal momentum than jet constituents, the cells is set
                to have zero four-momentum. This is deemed to be the case when the (scalar) pT of
                the thermal component is larger than the pT of the particle component.
                6. Re-cluster the jets with the cell four-momenta as input to get the final, subtracted jets.
                7. Calculate jet observables from re-clustered jets.
                '''
                fj_subtracted_constituents = self.subtract_thermal_gridsub1(jets_selected,jetR,gridsize,diagnostic)
                if not fj_subtracted_constituents:
                    continue
                cs = fj.ClusterSequence(fj_subtracted_constituents, jet_def)
                subtracted_jets = fj.sorted_by_pt(cs.inclusive_jets())

            else:
                raise NotImplementedError("Thermal subtraction method %s not recognized" % \
                                          self.thermal_subtraction_method)

            #----------------------------------------------

            for i, jet in enumerate(subtracted_jets):
                self.analyze_accepted_jet(jet, jetR, gridsize, diagnostic)

    #---------------------------------------------------------------
    # Fill histograms
    #---------------------------------------------------------------
    def analyze_accepted_jet(self, jet, jetR, gridsize, diagnostic=False):

        # Fill base histograms
        jet_pt_ungroomed = jet.pt()

        # Loop through each jet subconfiguration (i.e. subobservable / grooming setting)
        # Note that the subconfigurations are defined by the first observable, if multiple are defined
        for observable in self.observable_list:
          for i in range(len(self.obs_settings[observable])):
            obs_setting = self.obs_settings[observable][i]
            grooming_setting = self.obs_grooming_settings[observable][i]
            obs_label = self.utils.obs_label(obs_setting, grooming_setting)

            # Groom jet, if applicable
            jet_def = fj.JetDefinition(self.reclustering_algorithm, jetR)
            if grooming_setting:
              # For negative_recombiner case, we set the negative recombiner
              if "negative_recombiner" in self.thermal_subtraction_method:
                recombiner = fjext.NegativeEnergyRecombiner(-1)
                jet_def.set_recombiner(recombiner)
              gshop = fjcontrib.GroomerShop(jet, jet_def)
              jet_groomed_lund = self.utils.groom(gshop, grooming_setting, jetR)
              if not jet_groomed_lund:
                continue

            else:
              jet_groomed_lund = None

            # Call user function to fill histograms
            if not diagnostic:
              self.fill_jet_histograms(
                observable, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
                obs_label, jet_pt_ungroomed, gridsize)
            else:
              self.fill_jet_histograms(
                observable, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
                obs_label, jet_pt_ungroomed, gridsize, suffix='_diagnostics')

    #---------------------------------------------------------------
    # GridSub1 subtraction method
    # (jet clustering before discretization)
    #---------------------------------------------------------------
    def subtract_thermal_gridsub1(self,jets,jetR,gridsize,diagnostic=False):
        m = 0.1396
        #-----------------------------------
        # Clean up what's leftover from the previous event
        self.reset_grids(gridsize) # empty grid from previous event
        self.populated_cells = [] # empty this list, which was populated in the previous event
        self.populated_cells_w_constit = [] # empty this list, which was populated in the previous event
        unsubtracted_thermal_pT_fraction = -1000
        #-----------------------------------

        # Loop over all jets in the event
        for jet in jets:
            jet_pt  = jet.pt()
            jet_eta = jet.eta()
            jet_phi = jet.phi()

            if jet_phi > math.pi:
                jet_phi -= 2. * math.pi

            # Add jet particles to grid
            for constit in jet.constituents():
                self.populate_grid_with_constituents(constit, gridsize)

            # Subtract thermal particles from grid
            if not diagnostic:
                for th in range(0, len(self.event_bck_df)):
                    self.subtract_thermals_from_grid(th, gridsize, jet_eta, jet_phi, jetR)

                if self.total_thermal_momentum > 0:
                    unsubtracted_thermal_pT_fraction = \
                        self.unsubtracted_thermal_momentum / self.total_thermal_momentum
                    name = 'h_thermal_fraction_not_subtracted_v_pT_R{}_gridsize{}'.format(
                        jetR, gridsize)
                    getattr(self, name).Fill(jet_pt, unsubtracted_thermal_pT_fraction)

        # Create dataframe to store the surviving constituents
        jet_df = pd.DataFrame(columns = ['ParticlePt', 'ParticleEta', 'ParticlePhi','m'])

        # Loop over grid
        #for k1 in self.grid_dict[gridsize].keys():
        #    for k2 in self.grid_dict[gridsize][k1].keys():
        for pair in self.populated_cells:
            grid_cell = self.grid_dict[gridsize][pair[0]][pair[1]]

            # If a cell has negative pT, set the 4 momentum to 0
            if grid_cell.Pt()<0:
                grid_cell.SetPtEtaPhiM(0,0,0,0)

            pt  = grid_cell.Pt()
            eta = grid_cell.Eta()
            phi = grid_cell.Phi()

            # Append to dataframe all grid entries with positive pT
            if pt > 0:
                jet_df = jet_df.append({'ParticlePt':grid_cell.Pt(),
                'ParticleEta':grid_cell.Eta(),
                'ParticlePhi':grid_cell.Phi(),
                'm':m},
                ignore_index = True)

        # If no entry was stored in grid at this point, return None
        if len(jet_df) == 0:
            print('This should never happen! The end of the world must be close')
            return None

        # Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of pt,eta,phi
        fj_particles = fjext.vectorize_pt_eta_phi_m(
          jet_df['ParticlePt'].values, jet_df['ParticleEta'].values,
          jet_df['ParticlePhi'].values, jet_df['m'].values, 0)

        return fj_particles

    #---------------------------------------------------------------
    # Add jet constituents to grid
    #---------------------------------------------------------------
    def populate_grid_with_constituents(self,constit,gridsize):
        m = 0.1396
        c_eta_idx = self.histutils.find_cell(constit.eta(),array('d',self.cell_eta[gridsize]),len(self.cell_phi[gridsize]),False)
        c_phi_idx = self.histutils.find_cell(constit.phi(),array('d',self.cell_phi[gridsize]),len(self.cell_phi[gridsize]),True)
        
        pair = [c_eta_idx,c_phi_idx]
        if pair not in self.populated_cells:
            self.populated_cells.append(pair)
            self.populated_cells_w_constit.append(pair)

        # Store the information about this constituent in a Lorentz vector and add to grid
        const_mom = ROOT.TLorentzVector()
        const_mom.SetPtEtaPhiM(constit.pt(),constit.eta(),constit.phi(),m)
        self.grid_dict[gridsize][c_eta_idx][c_phi_idx] += const_mom

    #---------------------------------------------------------------
    # Subtract thermal constituents from grid
    #---------------------------------------------------------------
    def subtract_thermals_from_grid(self,th,gridsize,jet_eta,jet_phi,jetR):
        m = 0.1396
        self.total_thermal_momentum = 0
        self.unsubtracted_thermal_momentum = 0

        if np.random.uniform() >= self.thermal_rejection_fraction:

            th_eta = self.event_bck_df['ParticleEta'].iloc[th]
            th_phi = self.event_bck_df['ParticlePhi'].iloc[th]
            th_pt  = self.event_bck_df['ParticlePt' ].iloc[th]

            c_eta_idx = self.histutils.find_cell(th_eta,array('d',self.cell_eta[gridsize]),len(self.cell_phi[gridsize]),False)
            c_phi_idx = self.histutils.find_cell(th_phi,array('d',self.cell_phi[gridsize]),len(self.cell_phi[gridsize]),True )
            pair = [c_eta_idx,c_phi_idx]

            stored_pT = self.grid_dict[gridsize][c_eta_idx][c_phi_idx].Pt()
            # If angular distance between thermal and jet axis < jet R, subtract thermal 4 momentum from cell
            if math.sqrt( pow(jet_eta-th_eta,2) + pow(jet_phi-th_phi,2)) < jetR:
                temp_background = ROOT.TLorentzVector()
                temp_background.SetPtEtaPhiM(th_pt,th_eta,th_phi,m)

                if stored_pT - temp_background.Pt() > 0:
                    self.grid_dict[gridsize][c_eta_idx][c_phi_idx] -= temp_background
                else:
                    self.grid_dict[gridsize][c_eta_idx][c_phi_idx] = ROOT.TLorentzVector(0,0,0,0)

                self.total_thermal_momentum += temp_background.Pt()
                if pair not in self.populated_cells_w_constit:
                    self.unsubtracted_thermal_momentum += temp_background.Pt()

    #---------------------------------------------------------------
    # Create grids
    #---------------------------------------------------------------
    def create_grids(self,gridsize=0.05):
        print('-------------------------------------------')
        print('Creating grid with size: ',gridsize)
        # phi binning [-pi,pi]
        min_phi = -3.141593
        max_phi = 3.141593
        steps_phi = (int)((max_phi-min_phi)/gridsize)
        self.cell_phi[gridsize]=np.linspace(min_phi,max_phi,steps_phi)

        # eta binning [-0.9,0.9] <- ALICE tracking acceptance
        min_eta = -0.9
        max_eta = 0.9
        steps_eta = (int)((max_eta-min_eta)/gridsize)
        self.cell_eta[gridsize]=np.linspace(min_eta,max_eta,steps_eta)

        jet_dict = {} # Grid for the jet particles 
        for eta in self.cell_eta[gridsize]:
            dict_jet_phi = {}
            jet_dict[eta] = dict_jet_phi
            for phi in self.cell_phi[gridsize]:
                dict_jet_phi[phi] = None
                dict_jet_phi[phi] = ROOT.TLorentzVector(0,0,0,0)

        print('Grid edges in eta:\n',self.cell_eta[gridsize])
        print('Grid edges in phi:\n',self.cell_phi[gridsize])
        print('-------------------------------------------')
        return jet_dict

    #---------------------------------------------------------------
    # Reset grids content
    #---------------------------------------------------------------
    def reset_grids(self,gridsize):
        #for k1 in self.grid_dict[gridsize].keys():
        #    for k2 in self.grid_dict[gridsize][k1].keys():
        #        self.grid_dict[gridsize][k1][k2] = None
        #        self.grid_dict[gridsize][k1][k2] = ROOT.TLorentzVector(0,0,0,0)
        for pair in self.populated_cells:
            self.grid_dict[gridsize][pair[0]][pair[1]] = ROOT.TLorentzVector(0,0,0,0)

    #---------------------------------------------------------------
    # Find index of jet constituent in cell
    # (Deprecated. Moved to c++ in RUtils)
    #---------------------------------------------------------------
    '''
    def find_cell(self,val,cell,phi):
        if phi and val > math.pi:
            #print('got input for phi =',val)
            val = val - 2*math.pi
            #print('changed it to',val)
        for e in range(0,len(cell)):
            if val >= cell[e] and val < cell[e+1]:
                return cell[e]
            elif val==cell[len(cell)-1]:
                return cell[len(cell)-2]
        print('value is out of bounds')
        print(val,cell)
        exit()
    '''
    #---------------------------------------------------------------
    # Save all histograms
    #---------------------------------------------------------------
    def save_output_objects(self,option='recreate'):

      outputfilename = os.path.join(self.output_dir, 'AnalysisResults.root')
      fout = ROOT.TFile(outputfilename, option)
      fout.cd()

      for attr in dir(self):

        obj = getattr(self, attr)

        # Write all ROOT histograms and trees to file
        types = (ROOT.TH1, ROOT.THnBase, ROOT.TTree)
        if isinstance(obj, types):
          obj.Write()

      fout.Close()

    #---------------------------------------------------------------
    # This function is called once for each jet subconfiguration
    # You must implement this
    #---------------------------------------------------------------
    def fill_jet_histograms(
      self, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting, obs_label,
      jet_pt_ungroomed, suffix=None):

      raise NotImplementedError('You must implement fill_jet_histograms()!')

    #---------------------------------------------------------------
    # This function is called once
    # You must implement this
    #---------------------------------------------------------------
    def initialize_user_output_objects(self):
        raise NotImplementedError('You must implement initialize_user_output_objects()!')

    #---------------------------------------------------------------
    # 4MomSub subtraction method
    #---------------------------------------------------------------
    def subtract_thermal_4momsub(self,jet,jetR,diagnostic=False):
        m = 0.1396
        jet_particles = 0.
        dummies = 0.
        background = ROOT.TLorentzVector()

        for constit in jet.constituents():
            if constit.pt() < 1e-5:
                for th in range(0, len(self.event_bck_df)):
                    th_eta = self.event_bck_df['ParticleEta'].iloc[th]
                    th_phi = self.event_bck_df['ParticlePhi'].iloc[th]
                    th_pt  = self.event_bck_df['ParticlePt' ].iloc[th]

                    if math.sqrt( pow(constit.eta()-th_eta,2) + pow(constit.phi()-th_phi,2)) < 1e-5:
                        temp_background = ROOT.TLorentzVector()
                        temp_background.SetPtEtaPhiM(th_pt,th_eta,th_phi,m)
                        background += temp_background
                        dummies += 1.
            else:
                jet_particles += 1.

        if jet_particles == 0:
            print('Clustered jet only contains dummies. Bailing out!')
            exit()
        
        jet_df = pd.DataFrame(columns = ['ParticlePt', 'ParticleEta', 'ParticlePhi','m'])
        for constit in jet.constituents():
            if constit.pt() > 1e-5 and dummies > 0:
                jet_df = jet_df.append({'ParticlePt':constit.pt()-background.Pt()/((float)(jet_particles)),
                'ParticleEta':constit.eta(),
                'ParticlePhi':constit.phi(),
                'm':m},
                ignore_index = True)
            elif constit.pt() > 1e-5:
                jet_df = jet_df.append({'ParticlePt':constit.pt(),
                'ParticleEta':constit.eta(),
                'ParticlePhi':constit.phi(),
                'm':m},
                ignore_index = True)

        if len(jet_df) == 0:
            return None

        # Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of pt,eta,phi
        fj_particles = fjext.vectorize_pt_eta_phi_m(
          jet_df['ParticlePt'].values, jet_df['ParticleEta'].values,
          jet_df['ParticlePhi'].values, jet_df['m'].values, 0)

        return fj_particles
