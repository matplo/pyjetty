#!/usr/bin/env python3

"""
Plot pp vs. AA classification performance
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
import matplotlib as mpl
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

        # Suffix for plot outputfile names
        self.roc_plot_index = 0
        self.significance_plot_index = 0
        self.auc_plot_index = 0

        self.plot_title = False
                
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

        self.models = config['models']
        self.K_list = config['K']
        self.K_max = config['K_max']
        self.dmax = config['dmax']
        self.efp_measure = config['efp_measure']
        self.efp_beta = config['efp_beta']

        self.constituent_subtraction_study = config['constituent_subtraction_study']
        
        # Initialize model-specific settings
        self.nsubjettiness_alpha_list = config['nsub_lasso']['alpha']
        self.nsubjettiness_alpha_plot_list = config['nsub_lasso']['alpha_plot']
        self.K_lasso_nsubjettiness = config['nsub_lasso']['K_lasso']

        self.efp_alpha_list = config['efp_lasso']['alpha']
        self.d_lasso_efp = config['efp_lasso']['d_lasso']

        self.alpha_lasso_external = 0.0007

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def plot_pp_aa(self):
    
        # Loop through combinations of event type, jetR, and R_max
        for event_type in self.event_types:
            for jetR in self.jetR_list:
                for jet_pt_bin in self.jet_pt_bins:
                    for R_max in self.max_distance_list:

                        if not self.models:
                            continue
                    
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

                        # Load ML results from file
                        self.key_suffix = f'_{event_type}_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                        roc_filename = os.path.join(self.output_dir_i, f'ROC{self.key_suffix}.pkl')
                        with open(roc_filename, 'rb') as f:
                            self.roc_curve_dict = pickle.load(f)
                            self.AUC = pickle.load(f)

                        lasso_filename = os.path.join(self.output_dir_i, f'ROC{self.key_suffix}_lasso.pkl')
                        if os.path.exists(lasso_filename):
                            with open(lasso_filename, 'rb') as f_lasso:
                                self.roc_curve_dict_lasso = pickle.load(f_lasso)
                                self.N_terms_lasso = pickle.load(f_lasso)
                                self.observable_lasso = pickle.load(f_lasso)

                        # Load external lasso file
                        if 'efp_lasso' in self.models:
                            external_lasso_file = '/rstorage/ml/egml/nsubjettiness/807291/felix_lasso_paper/ROC_hard_R0.4_pt[100.0, 125.0]_Rmax0_lasso_687.pkl'
                            if os.path.exists(external_lasso_file):
                                with open(external_lasso_file, 'rb') as f_external_lasso:
                                    self.roc_curve_dict_lasso_external = pickle.load(f_external_lasso)
                                    self.N_terms_lasso_external = pickle.load(f_external_lasso)
                                    self.observable_lasso_external = pickle.load(f_external_lasso)
                                self.N_terms_lasso['efp_lasso'][self.alpha_lasso_external] = self.N_terms_lasso_external[self.alpha_lasso_external]
                                self.observable_lasso['efp_lasso'][self.alpha_lasso_external] = self.observable_lasso_external[self.alpha_lasso_external]
                                print(self.observable_lasso_external[self.alpha_lasso_external])

                        # Plot models for a single setting
                        self.plot_models(event_type, jetR, jet_pt_bin, R_max)

        # Plots that compare hard vs. combined results
        for jetR in self.jetR_list:
            for jet_pt_bin in self.jet_pt_bins:
                for R_max in self.max_distance_list:
                
                    if not self.models:
                        continue

                    if np.isclose(R_max, 0.):
                        continue

                    output_dir_hard = os.path.join(self.output_dir, f'hard_R{jetR}_pt{jet_pt_bin}_Rmax0')
                    output_dir_combined = os.path.join(self.output_dir, f'combined_matched_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}')

                    # Plot PFN
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
                        roc_list[f'PFN_hard'] = roc_curve_dict_Rmax0['pfn']
                        roc_list[f'PFN_background'] = roc_curve_dict_Rmax025['pfn']
                        roc_list[f'PFN_hard_min_pt'] = roc_curve_dict_Rmax0['pfn_min_pt']
                        
                        outputdir = os.path.join(self.output_dir, f'combined_matched_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}')
                        self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max)

                    # Plot AUC vs. K
                    if 'pfn' in self.models and 'nsub_dnn' in self.models:
                        auc_list = {}
                        suffix = f'_hard_R{jetR}_pt{jet_pt_bin}_Rmax0'
                        auc_list[f'nsub_dnn{suffix}'] = self.AUC[f'nsub_dnn{suffix}']
                        auc_list[f'pfn{suffix}'] = self.AUC[f'pfn{suffix}']
                        suffix = f'_combined_matched_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
                        auc_list[f'nsub_dnn{suffix}'] = self.AUC[f'nsub_dnn{suffix}']
                        auc_list[f'pfn{suffix}'] = self.AUC[f'pfn{suffix}']
                        outputdir = os.path.join(self.output_dir, f'combined_matched_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}')
                        self.plot_AUC_convergence(outputdir, auc_list, event_type, jetR, jet_pt_bin, R_max)

        # Plot ROC curves of constituent subtraction study
        if self.constituent_subtraction_study['pfn'] or self.constituent_subtraction_study['nsub']:

            jetR = self.jetR_list[0]
            jet_pt_bin = self.jet_pt_bins[0]
            R_max = self.max_distance_list[-1]
            self.output_dir_i = os.path.join(self.output_dir, f'constituent_subtraction_R{jetR}_pt{jet_pt_bin}_Rmax0')
                
            self.key_suffix = f'_hard_R{jetR}_pt{jet_pt_bin}_Rmax{R_max}'
            roc_filename = os.path.join(self.output_dir_i, f'ROC_constituent_subtraction.pkl')
            with open(roc_filename, 'rb') as f:
                roc_curve_dict = pickle.load(f)

            if self.constituent_subtraction_study['pfn']:

                # Plot PFN
                roc_list = {}
                roc_list['pfn_hard'] = roc_curve_dict['pfn_hard']
                roc_list['pfn_beforeCS'] = roc_curve_dict['pfn_beforeCS']
                roc_list['pfn_afterCS'] = roc_curve_dict['pfn_afterCS']
                self.plot_roc_curves(roc_list, 'CS', jetR, jet_pt_bin, R_max)

                # Plot PFN w/min_pt cut
                roc_list = {}
                roc_list['pfn_hard_min_pt'] = roc_curve_dict['pfn_hard_min_pt']
                roc_list['pfn_beforeCS_min_pt'] = roc_curve_dict['pfn_beforeCS_min_pt']
                roc_list['pfn_afterCS_min_pt'] = roc_curve_dict['pfn_afterCS_min_pt']
                self.plot_roc_curves(roc_list, 'CS', jetR, jet_pt_bin, R_max)

            if self.constituent_subtraction_study['nsub']:

                # Plot Nsubjettiness
                roc_list = {}
                roc_list['nsub_hard'] = roc_curve_dict['nsub_hard'][self.K_max]
                roc_list['nsub_beforeCS'] = roc_curve_dict['nsub_beforeCS'][self.K_max]
                roc_list['nsub_afterCS'] = roc_curve_dict['nsub_afterCS'][self.K_max]
                self.plot_roc_curves(roc_list, 'CS', jetR, jet_pt_bin, R_max)
                
    #---------------------------------------------------------------
    # Plot several versions of ROC curves and significance improvement
    #---------------------------------------------------------------
    def plot_models(self, event_type, jetR, jet_pt_bin, R_max):

        if 'pfn' in self.models and 'efn' in self.models:
            roc_list = {}
            roc_list['PFN'] = self.roc_curve_dict['pfn']
            roc_list['EFN'] = self.roc_curve_dict['efn']
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max)

        if 'pfn' in self.models and 'nsub_dnn' in self.models:
            roc_list = {}
            roc_list['PFN'] = self.roc_curve_dict['pfn']
            for K in self.K_list:
                roc_list[f'Nsub (M = {K}), DNN'] = self.roc_curve_dict['nsub_dnn'][K]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max)
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max)
        
        if 'efp_linear' in self.models:
             roc_list = {}
             for d in range(3, self.dmax+1):
                 roc_list[f'EFP (d = {d}), Linear'] = self.roc_curve_dict['efp_linear'][d]
             self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max)

        if 'nsub_lasso' in self.models and 'efp_lasso' in self.models:
            roc_list = {}
            roc_list['jet_angularity'] = self.roc_curve_dict_lasso['jet_angularity']
            roc_list['thrust'] = self.roc_curve_dict_lasso['thrust']
            roc_list['jet_theta_g'] = self.roc_curve_dict_lasso['jet_theta_g']
            roc_list['zg'] = self.roc_curve_dict_lasso['zg']
            for alpha in self.nsubjettiness_alpha_plot_list:
                roc_list[rf'Lasso $(\alpha = {alpha})$, Nsub'] = self.roc_curve_dict_lasso['nsub_lasso'][alpha]
            roc_list[rf'Lasso $(\alpha = {self.alpha_lasso_external})$, EFP'] = self.roc_curve_dict_lasso_external['efp'][self.alpha_lasso_external]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max)

        if 'nsub_lasso' in self.models and 'efp_lasso' in self.models:
            roc_list = {}
            roc_list['jet_angularity'] = self.roc_curve_dict_lasso['jet_angularity']
            roc_list['thrust'] = self.roc_curve_dict_lasso['thrust']
            roc_list['jet_theta_g'] = self.roc_curve_dict_lasso['jet_theta_g']
            roc_list['zg'] = self.roc_curve_dict_lasso['zg']
            for alpha in self.nsubjettiness_alpha_plot_list:
                roc_list[rf'Lasso $(\alpha = {alpha})$, Nsub'] = self.roc_curve_dict_lasso['nsub_lasso'][alpha]
            alpha_efp = self.efp_alpha_list[1]
            roc_list[rf'Lasso $(\alpha = {alpha_efp})$, EFP'] = self.roc_curve_dict_lasso['efp_lasso'][alpha_efp]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max)
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max)

        if 'nsub_dnn' in self.models and 'nsub_lasso' in self.models and 'efp_dnn' in self.models and 'efp_lasso' in self.models:
            roc_list = {}
            roc_list[f'Nsub (M = {self.K_lasso_nsubjettiness}), DNN'] = self.roc_curve_dict['nsub_dnn'][self.K_lasso_nsubjettiness]
            roc_list[f'EFP (d = {self.d_lasso_efp})'] = self.roc_curve_dict['efp_linear'][self.d_lasso_efp]
            roc_list['thrust'] = self.roc_curve_dict_lasso['thrust']
            roc_list['jet_angularity'] = self.roc_curve_dict_lasso['jet_angularity']
            roc_list['jet_theta_g'] = self.roc_curve_dict_lasso['jet_theta_g']
            roc_list['zg'] = self.roc_curve_dict_lasso['zg']
            alpha_nsubjettiness = self.nsubjettiness_alpha_list[0]
            roc_list[rf'Lasso $(\alpha = {alpha_nsubjettiness})$, Nsub'] = self.roc_curve_dict_lasso['nsub_lasso'][alpha_nsubjettiness]
            alpha_efp = self.efp_alpha_list[0]
            roc_list[rf'Lasso $(\alpha = {alpha_efp})$, EFP'] = self.roc_curve_dict_lasso['efp_lasso'][alpha_efp]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max)
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max)

        if 'efp_dnn' in self.models:
             roc_list = {}
             for d in range(3, self.dmax+1):
                 roc_list[f'EFP (d = {d}), DNN'] = self.roc_curve_dict['efp_dnn'][d]
             self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max)

        if 'nsub_dnn' in self.models:
            roc_list = {}
            for K in self.K_list:
                roc_list[f'Nsub (M = {K}), DNN'] = self.roc_curve_dict['nsub_dnn'][K]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max)
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max)

        if 'nsub_linear' in self.models:
            roc_list = {}
            for K in self.K_list:
                roc_list[f'Nsub (M {K}), Linear'] = self.roc_curve_dict['nsub_linear'][K]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max)
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max)

        if 'nsub_dnn' in self.models and 'nsub_lasso' in self.models:
            roc_list = {}
            roc_list[f'Nsub (M = {self.K_lasso_nsubjettiness}), DNN'] = self.roc_curve_dict['nsub_dnn'][self.K_lasso_nsubjettiness]
            roc_list['thrust'] = self.roc_curve_dict_lasso['thrust']
            roc_list['jet_angularity'] = self.roc_curve_dict_lasso['jet_angularity']
            roc_list['jet_theta_g'] = self.roc_curve_dict_lasso['jet_theta_g']
            roc_list['zg'] = self.roc_curve_dict_lasso['zg']
            for alpha in self.nsubjettiness_alpha_list:
                roc_list[rf'Lasso $(\alpha = {alpha})$, Nsub'] = self.roc_curve_dict_lasso['nsub_lasso'][alpha]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max)
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max)

        if 'efp_linear' in self.models and 'efp_lasso' in self.models:
            roc_list = {}
            roc_list[f'EFP (d = {self.d_lasso_efp}), Linear'] = self.roc_curve_dict['efp_linear'][self.d_lasso_efp]
            roc_list['thrust'] = self.roc_curve_dict_lasso['thrust']
            roc_list['jet_angularity'] = self.roc_curve_dict_lasso['jet_angularity']
            roc_list['jet_theta_g'] = self.roc_curve_dict_lasso['jet_theta_g']
            roc_list['zg'] = self.roc_curve_dict_lasso['zg']
            for alpha in self.efp_alpha_list:
                roc_list[rf'Lasso $(\alpha = {alpha})$, EFP'] = self.roc_curve_dict_lasso['efp_lasso'][alpha]
            self.plot_roc_curves(roc_list, event_type, jetR, jet_pt_bin, R_max)
            self.plot_significance_improvement(roc_list, event_type, jetR, jet_pt_bin, R_max)

    #--------------------------------------------------------------- 
    # Plot ROC curves
    #--------------------------------------------------------------- 
    def plot_roc_curves(self, roc_list, event_type, jetR, jet_pt_bin, R_max):
    
        plt.plot([0, 1], [0, 1], 'k--') # dashed diagonal
        plt.axis([0, 1, 0, 1])
        plt.title('JEWEL vs. PYTHIA8     ' + r'$100<p_{\mathrm{T,jet}}<125\;\mathrm{GeV}$', fontsize=14)
        if self.plot_title:
            plt.title(rf'{event_type} event: $R = {jetR}, p_T = {jet_pt_bin}, R_{{max}} = {R_max}$', fontsize=14)
        plt.xlabel('False Positive Rate (AA)', fontsize=16)
        plt.ylabel('True Positive Rate (AA)', fontsize=16)
        #plt.ylabel(rf'$\varepsilon_{{\rm{{true positive}} }}^{{\rm{{AA}} }}$', fontsize=16)
        plt.grid(True)
    
        for label,value in roc_list.items():
            index=0
            if label in ['PFN', 'EFN', 'jet_mass', 'jet_angularity', 'LHA', 'thrust', 'pTD', 'hadron_z', 'zg', 'jet_theta_g'] or 'multiplicity' in label or 'PFN' in label or 'EFN' in label or 'pfn' in label:
                linewidth = 4
                alpha = 0.5
                linestyle = self.linestyle(label)
                color=self.color(label)
                legend_fontsize = 12
                if label == 'jet_mass':
                    label = r'$m_{\mathrm{jet}}$'
                if label == 'jet_angularity':
                    label = r'$\lambda_1$ (girth)'
                    linewidth = 2
                if label == 'thrust':
                    label = r'$\lambda_2$ (thrust)'
                    linewidth = 2
                if label == 'jet_theta_g':
                    label = r'$\theta_{\mathrm{g}}$'
                    linewidth = 2
                if label == 'zg':
                    label = r'$z_{\mathrm{g}}$'
                    linewidth = 2
                if label == 'PFN':
                    label = 'Particle Flow Network'
                if label == 'EFN':
                    label = 'Energy Flow Network'
                if label == 'PFN_hard':
                    label = 'Jet'
                if label == 'PFN_hard_min_pt':
                    label = rf'Jet ($p_{{ \rm{{T}} }}^{{ \rm{{particle}} }} > 1$ GeV)'
                if label == 'PFN_background':
                    label = 'Jet + Background'
                if label == 'pfn_hard':
                    label = 'Jet'
                if label == 'pfn_beforeCS':
                    label = 'Jet + Background (before subtraction)'                
                if label == 'pfn_afterCS':
                    label = 'Jet + Background (after subtraction)'                
            elif 'Lasso' in label:
                linewidth = 4
                alpha = 0.5
                linestyle = 'solid'
                color=self.color(label)
                legend_fontsize = 10
                reg_param = float(re.search('= (.*)\)', label).group(1))
                if 'Nsub' in label:
                    n_terms = self.N_terms_lasso['nsub_lasso'][reg_param]
                    print(self.observable_lasso['nsub_lasso'][reg_param])
                elif 'EFP' in label:
                    n_terms = self.N_terms_lasso['efp_lasso'][reg_param]

                if 'Nsub' in label:
                    label = rf'$\prod_{{N,\beta}} \left( \tau_N^{{\beta}} \right) ^{{c_{{N,\beta}} }}$'
                elif 'EFP' in label:
                    label = rf'$\sum_{{G}} c_{{G}} \rm{{EFP}}_{{G}}$'
                label += f', {n_terms} terms'

            elif 'DNN' in label or 'EFP' in label or 'nsub' in label:
                linewidth = 2
                alpha = 0.9
                linestyle = 'solid'
                color=self.color(label)
                legend_fontsize = 12
            else:
                linewidth = 2
                linestyle = 'solid'
                alpha = 0.9
                color = sns.xkcd_rgb['almost black']
                legend_fontsize = 12
  
            FPR = value[0]
            TPR = value[1]
            plt.plot(FPR, TPR, linewidth=linewidth, label=label,
                     linestyle=linestyle, alpha=alpha, color=color)
                    
        plt.legend(loc='lower right', fontsize=legend_fontsize)

        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir_i, f'ROC_{self.roc_plot_index}.pdf'))
        plt.close()

        self.roc_plot_index += 1
 
    #--------------------------------------------------------------- 
    # Plot Significance improvement
    #--------------------------------------------------------------- 
    def plot_significance_improvement(self, roc_list, event_type, jetR, jet_pt_bin, R_max):
    
        plt.axis([0, 1, 0, 3])
        if self.plot_title:
            plt.title(rf'{event_type} event: $R = {jetR}, p_T = {jet_pt_bin}, R_{{max}} = {R_max}$', fontsize=14)
        plt.xlabel('True Positive Rate', fontsize=16)
        plt.ylabel('Significance improvement', fontsize=16)
        plt.grid(True)
            
        for label,value in roc_list.items():
            if label in ['PFN', 'EFN', 'jet_mass', 'jet_angularity', 'jet_theta_g'] or 'multiplicity' in label:
                linewidth = 4
                alpha = 0.5
                linestyle = self.linestyle(label)
            elif 'DNN' in label:
                linewidth = 2
                alpha = 0.9
                linestyle = 'solid'
            elif 'Lasso' in label:
                linewidth = 4
                alpha = 0.5
                linestyle = 'solid'
            else:
                linewidth = 2
                linestyle = 'solid'
                alpha = 0.9
                color = sns.xkcd_rgb['almost black']
                
            FPR = value[0]
            TPR = value[1]
            plt.plot(TPR, TPR/np.sqrt(FPR+0.001), linewidth=linewidth, label=label,
                     linestyle=linestyle, alpha=alpha, color=self.color(label))
         
        plt.legend(loc='best', fontsize=10)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir_i, f'Significance_improvement_{self.significance_plot_index}.pdf'))
        plt.close()

        self.significance_plot_index += 1

    #---------------------------------------------------------------
    # Plot AUC as a function of K
    #---------------------------------------------------------------
    def plot_AUC_convergence(self, output_dir, auc_list, event_type, jetR, jet_pt_bin, R_max):
    
        plt.axis([0, self.K_list[-1], 0, 1])
        if self.plot_title:
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
        plt.savefig(os.path.join(output_dir, f'AUC_convergence{self.auc_plot_index}.pdf'))
        plt.close()

        self.auc_plot_index += 1

    #---------------------------------------------------------------
    # Get color for a given label
    #---------------------------------------------------------------
    def color(self, label):

        color = None
        if label in ['jet_angularity', 'PFN', 'PFN_hard', 'pfn_hard', 'pfn_hard_min_pt', rf'Nsub (M = {self.K_list[4]}), DNN', 'EFP (d = 7), Linear']:
            color = sns.xkcd_rgb['faded purple'] 
        elif label in ['EFN', 'EFN_hard', 'jet_mass', 'LHA', 'hadron_z', 'thrust', 'pTD', 'multiplicity_0150', rf'Lasso $(\alpha = {self.alpha_lasso_external})$, EFP']:
            color = sns.xkcd_rgb['faded red']    
        elif label in ['PFN_background', 'pfn_afterCS', 'pfn_afterCS_min_pt', rf'Nsub (M = {self.K_list[3]}), DNN', 'EFP (d = 6), Linear', rf'Lasso $(\alpha = {self.nsubjettiness_alpha_plot_list[1]})$, Nsub', 'zg']:
            color = sns.xkcd_rgb['dark sky blue']    
        elif label in ['EFN_background', 'PFN_hard_min_pt', rf'Nsub (M = {self.K_list[2]}), DNN', 'EFP (d = 5), Linear', rf'Lasso $(\alpha = {self.nsubjettiness_alpha_plot_list[0]})$, Nsub']:
            color = sns.xkcd_rgb['medium green']  
        elif label in [rf'Nsub (M = {self.K_list[0]}), DNN', 'EFP (d = 3), Linear', rf'Lasso $(\alpha = {self.efp_alpha_list[1]})$, EFP']:
            color = sns.xkcd_rgb['watermelon'] 
        elif label in [rf'Nsub (M = {self.K_list[1]}), DNN', 'EFP (d = 4), Linear', rf'Lasso $(\alpha = {self.efp_alpha_list[0]})$, EFP']:
            color = sns.xkcd_rgb['light brown'] 
        elif label in ['jet_theta_g']:
            color = sns.xkcd_rgb['medium brown']
        else:
            color = sns.xkcd_rgb['almost black']
        #  'pfn_beforeCS', 'pfn_beforeCS_min_pt', 

        return color

    #---------------------------------------------------------------
    # Get linestyle for a given label
    #---------------------------------------------------------------
    def linestyle(self, label):
 
        linestyle = None
        if 'PFN' in label or 'EFN' in label or 'DNN' in label or 'pfn' in label:
            linestyle = 'solid'
        else:
            linestyle = 'dotted'

        return linestyle
            
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