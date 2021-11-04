#!/usr/bin/env python3

"""
Example class to read quark-gluon dataset
"""

import os
import sys
import argparse
import yaml
import re
import pickle

# Data analysis and plotting
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_context('paper', rc={'font.size':18,'axes.titlesize':18,'axes.labelsize':18})

# Base class
from pyjetty.alice_analysis.process.base import common_base

################################################################
class PlotPPAA(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, config_file='', output_dir='', **kwargs):
        super(common_base.CommonBase, self).__init__(**kwargs)
        
        self.config_file = config_file
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        
        # Initialize config file
        self.initialize_config()
        
        alpha_list = self.config['lasso']['alpha']
        self.colors = {'PFN': sns.xkcd_rgb['faded purple'],
                       'EFN': sns.xkcd_rgb['faded red'],
                       'PFN_hard': sns.xkcd_rgb['faded purple'],
                       'EFN_hard': sns.xkcd_rgb['faded red'],
                       'PFN_background': sns.xkcd_rgb['dark sky blue'],
                       'EFN_background': sns.xkcd_rgb['medium green'],
                       rf'DNN (K = {self.K_list[0]})': sns.xkcd_rgb['watermelon'],
                       rf'DNN (K = {self.K_list[1]})': sns.xkcd_rgb['light brown'],
                       rf'DNN (K = {self.K_list[2]})': sns.xkcd_rgb['medium green'],
                       rf'DNN (K = {self.K_list[3]})': sns.xkcd_rgb['dark sky blue'],
                       rf'DNN (K = {self.K_list[4]})': sns.xkcd_rgb['faded purple'],
                       'Lasso': sns.xkcd_rgb['watermelon'],
                       rf'Lasso $(\alpha = {alpha_list[0]})$': sns.xkcd_rgb['watermelon'],
                       rf'Lasso $(\alpha = {alpha_list[1]})$': sns.xkcd_rgb['light brown'],
                       rf'Lasso $(\alpha = {alpha_list[2]})$': sns.xkcd_rgb['faded purple'],
                       'jet_mass': sns.xkcd_rgb['faded red'],
                       'LHA': sns.xkcd_rgb['faded red'],
                       'thrust': sns.xkcd_rgb['faded red'],
                       'hadron_z': sns.xkcd_rgb['faded red'],
                       'zg': sns.xkcd_rgb['faded red'],
                       'pTD': sns.xkcd_rgb['faded red'],
                       'multiplicity_0150': sns.xkcd_rgb['faded red'],
                       'jet_angularity': sns.xkcd_rgb['medium green'],
                       'jet_theta_g': sns.xkcd_rgb['medium brown'],
                       'PFN_hard': sns.xkcd_rgb['faded purple'],
                       'PFN_beforeCS': sns.xkcd_rgb['faded red'],
                       'PFN_afterCS': sns.xkcd_rgb['medium green'],
                       'EFP (d = 3)': sns.xkcd_rgb['watermelon'],
                       'EFP (d = 4)': sns.xkcd_rgb['light brown'],
                       'EFP (d = 5)': sns.xkcd_rgb['medium green'],
                       'EFP (d = 6)': sns.xkcd_rgb['dark sky blue'],
                       'EFP (d = 7)': sns.xkcd_rgb['faded purple']
                      }
        self.linestyles = {'PFN': 'solid',
                           'EFN': 'solid',
                           'DNN': 'solid',
                           'jet_mass': 'dotted',
                           'jet_angularity': 'dotted',
                           'jet_theta_g': 'dotted',
                           'LHA': 'dotted',
                           'thrust': 'dotted',
                           'hadron_z': 'dotted',
                           'zg': 'dotted',
                           'pTD': 'dotted',
                           'multiplicity_0150': 'dotted',
                           'Lasso': 'dotted',
                           'PFN_hard': 'solid',
                           'PFN_beforeCS': 'solid',
                           'PFN_afterCS': 'solid',
                           'EFN_hard': 'solid',
                           'PFN_background': 'solid',
                           'EFN_background': 'solid'
                          }
        
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
        self.K_list = config['K']
        self.models = config['models']
        self.K_lasso = config['lasso']['K_lasso']
        self.dmax = config['efp']['dmax']
        self.constituent_subtraction_study = config['constituent_subtraction_study']
        
        # Initialize model-specific settings
        self.config = config

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def plot_pp_aa(self):
    
        # Loop through combinations of event type, jetR, and R_max
        for event_type in self.event_types:
            for jetR in self.jetR_list:
                for jet_pt_bin in self.jet_pt_bins:
                    for R_max in self.max_distance_list:
                    
                        # For hard event, skip constituent subtraction variations
                        if event_type=='hard' and not np.isclose(R_max, 0.):
                            continue
                            
                        # If four-vectors are included, R_max=0 is skipped for combined event,
                        # since due to size/time constraints, we skip merging four-vectors for Rmax=0
                        if 'combined' in event_type and np.isclose(R_max,0):
                            continue

                        # Create output dir
                        self.output_dir_i = os.path.join(self.output_dir, f'{event_type}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}')
                        if not os.path.exists(self.output_dir_i):
                            os.makedirs(self.output_dir_i)

                        # Train models
                        self.plot_models(event_type, jetR, jet_pt_bin, R_max)

        # Plots for hard vs. combined
        for jetR in self.jetR_list:
            for jet_pt_bin in self.jet_pt_bins:
                for R_max in self.max_distance_list:
                
                    if np.isclose(R_max, 0.):
                        continue

                    output_dir_hard = os.path.join(self.output_dir, f'hard_R{jetR}_pt{jet_pt_bin}_Rmax0')
                    output_dir_combined = os.path.join(self.output_dir, f'combined_matched_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}')

                
                    # Plot AUC vs. K
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

                    # Plot PFN/EFN
                    if 'pfn' in self.models and 'efn' in self.models:

                        suffix_Rmax0 = f'_hard_R{jetR}_pt{jet_pt_bin}_Rmax0'
                        roc_filename = os.path.join(output_dir_hard, f'ROC{suffix_Rmax0}.pkl')
                        with open(roc_filename, 'rb') as f_Rmax0:
                            roc_curve_dict_Rmax0 = pickle.load(f_Rmax0)

                        suffix_Rmax025 = f'_combined_matched_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                        roc_filename = os.path.join(output_dir_combined, f'ROC{suffix_Rmax025}.pkl')
                        with open(roc_filename, 'rb') as f_Rmax025:
                            roc_curve_dict_Rmax025 = pickle.load(f_Rmax025)

                        roc_list = {}
                        roc_list[f'PFN_hard'] = roc_curve_dict_Rmax0['PFN']
                        roc_list[f'EFN_hard'] = roc_curve_dict_Rmax0['EFN']
                        roc_list[f'PFN_background'] = roc_curve_dict_Rmax025['PFN']
                        roc_list[f'EFN_background'] = roc_curve_dict_Rmax025['EFN']
                        
                        outputdir = os.path.join(self.output_dir, f'combined_matched_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}')
                        self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='6')

                    # Plot PFN only
                    if 'pfn' in self.models:

                        suffix_Rmax0 = f'_hard_R{jetR}_pt{jet_pt_bin}_Rmax0'
                        roc_filename = os.path.join(output_dir_hard, f'ROC{suffix_Rmax0}.pkl')
                        with open(roc_filename, 'rb') as f_Rmax0:
                            roc_curve_dict_Rmax0 = pickle.load(f_Rmax0)

                        suffix_Rmax025 = f'_combined_matched_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                        roc_filename = os.path.join(output_dir_combined, f'ROC{suffix_Rmax025}.pkl')
                        with open(roc_filename, 'rb') as f_Rmax025:
                            roc_curve_dict_Rmax025 = pickle.load(f_Rmax025)

                        roc_list = {}
                        roc_list[f'PFN_hard'] = roc_curve_dict_Rmax0['PFN']
                        roc_list[f'PFN_background'] = roc_curve_dict_Rmax025['PFN']
                        
                        outputdir = os.path.join(self.output_dir, f'combined_matched_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}')
                        self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='9')

        # Plot ROC curves of constituent subtraction study
        if self.constituent_subtraction_study:

            jetR = self.jetR_list[0]
            jet_pt_bin = self.jet_pt_bins[0]
            R_max = self.max_distance_list[-1]
            self.output_dir_i = os.path.join(self.output_dir, f'hard_R{jetR}_pt{jet_pt_bin}_Rmax0')
                
            self.key_suffix = f'_hard_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
            roc_filename = os.path.join(self.output_dir_i, f'ROC_constituent_subtraction.pkl')
            with open(roc_filename, 'rb') as f:
                roc_curve_dict = pickle.load(f)

            roc_list = {}
            for key in roc_curve_dict.keys():
                roc_list[key] = roc_curve_dict[key]
                self.plot_roc_curves(roc_list, 'hard', jetR, jet_pt_bin, R_max, suffix='7')
            
    #---------------------------------------------------------------
    # Plot models
    #---------------------------------------------------------------
    def plot_models(self, event_type, jetR, jet_pt_bin, R_max):

        # Load ML results from file
        self.key_suffix = f'_{event_type}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
        roc_filename = os.path.join(self.output_dir_i, f'ROC{self.key_suffix}.pkl')
        with open(roc_filename, 'rb') as f:
            self.roc_curve_dict = pickle.load(f)
            self.AUC = pickle.load(f)
            if 'lasso' in self.models:
                self.N_terms_lasso = pickle.load(f)

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
            for K in self.K_list:
                roc_list[f'DNN (K = {K})'] = self.roc_curve_dict['DNN'][K]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='1')
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='1')
        
        if 'pfn' in self.models and 'neural_network' in self.models:
            roc_list = {}
            roc_list['PFN'] = self.roc_curve_dict['PFN']
            for K in self.K_list:
                roc_list[f'DNN (K = {K})'] = self.roc_curve_dict['DNN'][K]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='2')
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='2')
        
        if 'neural_network' in self.models:
            roc_list = {}
            for K in self.K_list:
                roc_list[f'DNN (K = {K})'] = self.roc_curve_dict['DNN'][K]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='3')
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='3')
        
        if 'lasso' in self.models:
            roc_list = {}
            #roc_list['jet_mass'] = self.roc_curve_dict['jet_mass']
            roc_list['jet_angularity'] = self.roc_curve_dict['jet_angularity']
            roc_list['jet_theta_g'] = self.roc_curve_dict['jet_theta_g']
            #roc_list['hadron_z'] = self.roc_curve_dict['hadron_z']
            for alpha in self.config['lasso']['alpha']:
                roc_list[rf'Lasso $(\alpha = {alpha})$'] = self.roc_curve_dict['Lasso'][self.K_lasso][alpha]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='4')
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='4')
        
        if 'neural_network' in self.models and 'lasso' in self.models:
            roc_list = {}
            roc_list[f'DNN (K = {self.K_lasso})'] = self.roc_curve_dict['DNN'][K]
            roc_list['thrust'] = self.roc_curve_dict['thrust']
            roc_list['jet_angularity'] = self.roc_curve_dict['jet_angularity']
            roc_list['jet_theta_g'] = self.roc_curve_dict['jet_theta_g']
            for alpha in self.config['lasso']['alpha']:
                roc_list[rf'Lasso $(\alpha = {alpha})$'] = self.roc_curve_dict['Lasso'][self.K_lasso][alpha]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='5')
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='5')

        if 'efp' in self.models:
             roc_list = {}
             for d in range(3, self.dmax+1):
                 roc_list[f'EFP (d = {d})'] = self.roc_curve_dict['efp'][d]
             self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='8')

        # Plot several versions of ROC curves and significance improvement
        if 'pfn' in self.models and 'efn' in self.models:
            roc_list = {}
            roc_list['PFN'] = self.roc_curve_dict['PFN']
            roc_list['EFN'] = self.roc_curve_dict['EFN']
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max, suffix='10')

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
            index=0
            if label in ['PFN', 'EFN', 'jet_mass', 'jet_angularity', 'LHA', 'thrust', 'pTD', 'hadron_z', 'zg', 'jet_theta_g'] or 'multiplicity' in label or 'PFN' in label or 'EFN' in label:
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
            elif 'DNN' in label or 'EFP' in label:
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
            else:
                linewidth = 2
                linestyle = 'solid'
                alpha = 0.9
                color = sns.xkcd_rgb['almost black']
  
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
            
##################################################################
if __name__ == '__main__':

    # Define arguments
    parser = argparse.ArgumentParser(description='Plot ROC curves')
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

    analysis = PlotPPAA(config_file=args.configFile, output_dir=args.outputDir)
    analysis.plot_pp_aa()