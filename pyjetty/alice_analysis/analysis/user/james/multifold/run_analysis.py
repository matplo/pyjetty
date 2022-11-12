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
from array import *
import ROOT
from numba import jit
import h5py

# Suppress some annoying warnings
np.finfo(np.dtype("float32")) 
np.finfo(np.dtype("float64"))

from pyjetty.alice_analysis.analysis.base import common_base
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

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
class RunAnalysis(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, config_file='', **kwargs):
        super(RunAnalysis, self).__init__(**kwargs)
        self.config_file = config_file

        # Initialize utils class
        self.utils = analysis_utils_obs.AnalysisUtils_Obs()

        # Initialize yaml config
        self.initialize_config()

        # Create output dirs
        self.create_output_dirs()
        
        self.ColorArray = [ROOT.kBlue-4, ROOT.kAzure+7, ROOT.kCyan-2, ROOT.kViolet-8,
                        ROOT.kBlue-6, ROOT.kGreen+3, ROOT.kPink-4, ROOT.kRed-4,
                        ROOT.kOrange-3]
        self.MarkerArray = [20, 21, 22, 23, 33, 34, 24, 25, 26, 32]

        print(self)

        # Load aggregated results
        # The results are stored with keys: self.observable_info[observable]['obs_key'][i]
        # i.e. f'{observable}_{self.utils.obs_label(obs_setting, grooming_setting)}'
        self.results = self.utils.read_data(self.main_data)

        self.n_jets = next(iter(self.results.values())).shape
        print(f'Analyzing the following observables: (n_jets={self.n_jets})')
        for observable in self.observables:
            for i in range(self.observable_info[observable]['n_subobservables']):
                obs_key = self.observable_info[observable]['obs_keys'][i]
                print(f'  {obs_key}   (mean={np.mean(self.results[obs_key])})')
                if self.results[obs_key].shape != self.n_jets:
                    raise ValueError(f'{obs_key} has unexpected number of jets: {self.results[obs_key].shape}')
        print()

    #---------------------------------------------------------------
    # Initialize config file into class members
    #---------------------------------------------------------------
    def initialize_config(self):

        # Read config file
        with open(self.config_file, 'r') as stream:
            config = yaml.safe_load(stream)

        self.jetR_list = config['jetR']
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
            self.observable_info[observable]['obs_bins_truth'] = [obs_config_dict[obs_subconfig_list[i]]['obs_bins_truth']
                                                                  for i in range(len(obs_subconfig_list))]
            self.observable_info[observable]['obs_bins_det'] = [obs_config_dict[obs_subconfig_list[i]]['obs_bins_det']
                                                                for i in range(len(obs_subconfig_list))]
            self.observable_info[observable]['xtitle'] = config[observable]['common_settings']['xtitle']
            self.observable_info[observable]['ytitle'] = config[observable]['common_settings']['ytitle']

            # For jet_pt, store the binnings we will report
            if 'pt_bins_reported' in config[observable]['common_settings']:
                self.observable_info[observable]['pt_bins_reported'] = config[observable]['common_settings']['pt_bins_reported']

        # Directory to write analysis output to
        self.output_dir = config['output_dir']

        # Set which analysis steps to perform
        self.do_plot_observables = config['do_plot_observables']
        self.do_unfolding = config['do_unfolding']
        self.do_systematics = config['do_systematics']
        self.do_plot_final_result = config['do_plot_final_result']

        # List of systematic variations to perform
        self.systematics_list = config['systematics_list']

        # Load paths to processing output, to be unfolded
        self.main_data = config['main_data']
        self.main_response = config['main_response']

        if 'trkeff' in self.systematics_list:
            self.trkeff_response = config['trkeff_response']
        if 'fastsim_generator0' in self.systematics_list:
            self.fastsim_response_list = config['fastsim_response']

        # Figure info
        self.file_format = config['file_format']
        self.figure_approval_status = config['figure_approval_status']

    #---------------------------------------------------------------
    # Create a set of output directories for a given observable
    #---------------------------------------------------------------
    def create_output_dirs(self):

        if self.do_plot_observables:
            self.create_output_subdir(self.output_dir, 'plot_observables')

        if self.do_unfolding:
            for systematic in self.systematics_list:
                self.create_output_subdir(self.output_dir, systematic)

        if self.do_plot_final_result:
            self.create_output_subdir(self.output_dir, 'final_results')

    #---------------------------------------------------------------
    # Create a single output subdirectory
    #---------------------------------------------------------------
    def create_output_subdir(self, output_dir, name):

        output_subdir = os.path.join(output_dir, name)
        setattr(self, f'output_dir_{name}', output_subdir)
        if not os.path.isdir(output_subdir):
            os.makedirs(output_subdir)

        return output_subdir

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def run_analysis(self):

        if self.do_plot_observables:
            self.plot_observables()

        if self.do_unfolding:
            self.perform_unfolding()

        if self.do_systematics:
            self.compute_systematics()

        if self.do_plot_final_result:
            self.plot_results()

    #---------------------------------------------------------------
    # Plot observables before unfolding
    #---------------------------------------------------------------
    def plot_observables(self):

        # First, apply pt cuts to a copy of results
        pt_bins = self.observable_info['jet_pt']['pt_bins_reported']
        pt_bin_pairs = [[[pt_bins[i],pt_bins[i+1]]] for i in range(len(pt_bins)-1)]
        variables = [['jet_pt']] * len(pt_bin_pairs)
        results_dict = self.apply_cut(variables, pt_bin_pairs)

        for key,val in results_dict.items():
            n_jets = val['jet_pt'].shape[0]
            print(f'{key} has n_jets={n_jets}')

        print()
        print(f'Plotting...')
        for observable in self.observables:
            for i in range(self.observable_info[observable]['n_subobservables']):
                obs_key = self.observable_info[observable]['obs_keys'][i]
                self.results[obs_key]

        # Plot 1D projections of each observable
        self.plot_1D_projections()

    #---------------------------------------------------------------
    # Mask all results according to specified conditions
    #   variables = [[variable1], [variable2, variable3], ...]
    #   cuts = [ [[variable1_min,variable1_max]], [[variable2_min,variable2_max], [variable3_min,variable3_max]], ...]
    #
    # Each entry in the variables/cuts list corresponds to a set of cuts that will be simultanously appled
    #
    # A dictionary will be returned, containing the different cut combinations specified:
    #   result_dict[f'{variable1}{variable1_min}-{variable1_max}] = result1
    #   result_dict[f'{variable2}{variable2_min}-{variable2_max}_{variable3}{variable3_min}-{variable3_max}] = result2
    #
    # The reason for this is that we want to support both:
    #   - Returning multiple different cuts on the same variable (e.g. pt=[20,40,60,80])
    #   - Cutting on multiple variables simultaneously
    #---------------------------------------------------------------
    def apply_cut(self, variables_list, cuts_list):

        if len(variables_list) != len(cuts_list):
            raise ValueError(f'variables_list has different length than cuts_list! {variables_list} vs. {cuts_list}')

        results_dict = self.recursive_defaultdict()
        
        # Loop over all cut combinations
        for variables, cuts in zip(variables_list, cuts_list):

            # Loop through all cuts in a given combination and construct mask
            total_mask = np.zeros(self.n_jets, dtype=bool) # For simplicity: cut_mask[i]=True means we will discard index i
            cut_key = ''
            for i in range(len(variables)):
                variable = variables[i]
                cut_min, cut_max = cuts[i]
                cut_mask = (self.results[variable] < cut_min) | (self.results[variable] > cut_max)
                total_mask = (total_mask | cut_mask)

                if i>0:
                    cut_key += '_'
                cut_key += f'{variable}{cut_min}-{cut_max}'

            # Now that we have the complete mask for this combination, apply it to all arrays
            for obs_key,result in self.results.copy().items():
                results_dict[cut_key][obs_key] = result[~total_mask]
                    
        return results_dict

    #---------------------------------------------------------------
    # This function is called once for each subconfiguration
    #---------------------------------------------------------------
    def plot_1D_projections(self):   

        self.create_output_subdir(self.output_dir_plot_observables, '1D')

        # Normalize by integral, i.e. N_jets,inclusive in this pt-bin
    
        # First scale by bin width -- then normalize by integral
        # (where integral weights by bin width)
        #h.Scale(1., 'width')
        #
        #if grooming_setting and 'sd' in grooming_setting:
        #    # If SD, the untagged jets are in the first bin
        #    n_jets_inclusive = h.Integral(minbin, maxbin, 'width')
        #    n_jets_tagged = h.Integral(minbin+1, maxbin, 'width')
        #    f_tagging = n_jets_tagged/n_jets_inclusive
        #else:
        #    n_jets_inclusive = h.Integral(minbin, maxbin, 'width')
        # 
        #h.Scale(1./n_jets_inclusive)

        #self.utils.formatted_grooming_label(grooming_setting)
        #self.utils.formatted_subobs_label(self.observable)

    #---------------------------------------------------------------
    # This function is called once for each subconfiguration
    #---------------------------------------------------------------
    def plot_results(self):

        # Loop through jet radii
        for jetR in self.jetR_list:

            # Loop through subconfigurations
            for i, subconfig in enumerate(self.obs_subconfig_list):

                obs_setting = self.obs_settings[i]
                grooming_setting = self.grooming_settings[i]
                obs_label = self.obs_labels[i]

                # Plot result for each individial subconfiguration
                self.plot_single_result(jetR, obs_label, obs_setting, grooming_setting)
                
    #---------------------------------------------------------------
    # This function is called once for each subconfiguration
    #---------------------------------------------------------------
    def plot_single_result(self, jetR, obs_label, obs_setting, grooming_setting):    
        print(f'Plot final results for {self.observable}: R = {jetR}, {obs_label}')

        self.utils.set_plotting_options()
        ROOT.gROOT.ForceStyle()
        
        # Loop through pt slices, and plot final result for each 1D theta_g distribution
        #for bin in range(0, len(self.pt_bins_reported) - 1):
        #    min_pt_truth = self.pt_bins_reported[bin]
        #    max_pt_truth = self.pt_bins_reported[bin+1]
        #
        #    self.plot_observable(jetR, obs_label, obs_setting, grooming_setting, min_pt_truth, max_pt_truth, plot_pythia=True)

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