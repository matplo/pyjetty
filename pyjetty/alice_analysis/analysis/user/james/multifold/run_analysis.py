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
import pandas as pd
import array
import ROOT
import ctypes

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_context('paper', rc={'font.size':18, 'axes.titlesize':18, 'axes.labelsize':18, 'text.usetex':True})

import sklearn.preprocessing

# Suppress some annoying warnings
np.finfo(np.dtype("float32")) 
np.finfo(np.dtype("float64"))

from pyjetty.alice_analysis.analysis.base import common_base
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs

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
        self.utils = analysis_utils_obs.AnalysisUtils_Obs()

        # Initialize yaml config
        self.initialize_config()
        
        self.colors = [ROOT.kBlue-4, ROOT.kGreen-2, ROOT.kRed-4, 
                       ROOT.kAzure+7, ROOT.kCyan-2, ROOT.kViolet-8,
                       ROOT.kBlue-6, ROOT.kPink-4, ROOT.kOrange-3]
        self.markers = [20, 21, 22, 23, 33, 34, 24, 25, 26, 32]
        self.fillstyles = [1001, 3144, 1001, 3144]

        self.data_color = ROOT.kGray+3
        self.data_marker = 21
        self.theory_marker = 20
        self.marker_size = 1.5
        self.alpha = 0.6
        self.line_width = 2
        self.line_style = 1

        print(self)

        # Load aggregated results for both data and MC
        # The results will be stored with keys: 'data', 'mc_det_matched', 'mc_truth_matched'
        # The arrays are stored with keys specified in self.observable_info: self.observable_info[observable]['obs_key'][i]
        #   i.e. f'{observable}_{self.utils.obs_label(obs_setting, grooming_setting)}'
        print()
        print('Loading data...')
        results = {}
        results['data'] = self.utils.read_data(self.main_data)['data']
        print()
        print('Loading MC...')
        results_mc = self.utils.read_data(self.main_response)
        results['mc_det_matched'] = results_mc['det_matched']
        results['mc_truth_matched'] = results_mc['truth_matched']
        print()

        self.n_jets = {}
        for data_type in results.keys():
            self.n_jets[data_type] = next(iter(results[data_type].values())).shape[0]
        n_jets_total_data = self.n_jets['data'] 
        n_jets_total_mc = self.n_jets['mc_truth_matched']
        print(f'Analyzing the following observables: (n_jets={n_jets_total_data} in data, n_jets={n_jets_total_mc} in mc)')
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
        for data_type in results:
            for observable in self.observables:
                for i in range(self.observable_info[observable]['n_subobservables']):
                    obs_key = self.observable_info[observable]['obs_keys'][i]

                    # Check number of jets is as expected for each observable
                    if results[data_type][obs_key].shape[0] != self.n_jets[data_type] or idx[data_type].shape[0] != self.n_jets[data_type]:
                        shape = results[data_type][obs_key].shape
                        raise ValueError(f'{data_type} {obs_key} in data has unexpected number of jets: {shape}')

                    # Shuffle
                    self.results[data_type][obs_key] = results[data_type][obs_key][idx[data_type]]
                    if 'mc' in data_type:
                        self.results[data_type]['total_scale_factor'] = results[data_type]['total_scale_factor'][idx[data_type]]
        print('Done.')

        # Make some labels for plotting
        self.data_type_label = {'data': 'Uncorrected Data',
                                'mc_det_matched': 'MC det-level (matched)',
                                'mc_truth_matched': 'MC truth-level (matched)',
                                }

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
            print(observable)

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
            self.observable_info[observable]['obs_bins'] = [self.bins_from_hepdata(obs_config_dict[obs_subconfig_list[i]]) 
                                                                                   for i in range(len(obs_subconfig_list))]
            self.observable_info[observable]['hepdata_tgraph'] = [self.tgraph_from_hepdata(obs_config_dict[obs_subconfig_list[i]]) 
                                                                                           for i in range(len(obs_subconfig_list))]
                                                                  
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

        # Create output dir
        self.create_output_subdir(self.output_dir, 'plot_observables')

        # First, apply pt cuts to a copy of results
        # TODO: bin both mc_det_matched and mc_truth_matched by e.g. truth-level pt
        print('Binning results into pt bins to make some plots...')
        print()
        pt_bins = self.observable_info['jet_pt']['pt_bins_reported']
        pt_bin_pairs = [[pt_bins[i],pt_bins[i+1]] for i in range(len(pt_bins)-1)]
        pt_bin_pairs_nested = [[pt_bin_pair] for pt_bin_pair in pt_bin_pairs]
        variables = [['jet_pt']] * len(pt_bin_pairs_nested)
        results_dict = self.apply_cut(variables, pt_bin_pairs_nested, n_jets_max=1000000)

        # Print number of jets passing each cut
        for data_type in results_dict.keys():
            n_jets = [results_dict[data_type]['jet_pt'][f'jet_pt{x[0]}-{x[1]}'].shape[0] for x in pt_bin_pairs]
            print(data_type)
            [print(f'  pt={pt_bin_pairs[i]} has n_jets={n_jets[i]}') for i in range(len(n_jets))]

        # Plot 1D projections of each observable
        print()
        print(f'Plotting 1D projections...')
        self.plot_1D_projections(results_dict, pt_bin_pairs)

        # Plot pairplot between each pair of observables
        print()
        print(f'Plotting pairplot...')
        self.plot_pairplot(results_dict, pt_bin_pairs)
        print('Done!')

    #---------------------------------------------------------------
    # Mask all results according to specified conditions
    #   variables = [[variable1], [variable2, variable3], ...]
    #   cuts = [ [[variable1_min,variable1_max]], [[variable2_min,variable2_max], [variable3_min,variable3_max]], ...]
    #
    # Each entry in the variables/cuts list corresponds to a set of cuts that will be simultanously appled
    #
    # A dictionary will be returned, containing the different cut combinations specified, e.g.:
    #   result_dict['data'][f'{variable1}{variable1_min}-{variable1_max}] = result1
    #   result_dict['data'][f'{variable2}{variable2_min}-{variable2_max}_{variable3}{variable3_min}-{variable3_max}] = result2
    # The dictionary will include both data and MC and det/truth levels.
    #
    # The reason for this is that we want to support both:
    #   - Returning multiple different cuts on the same variable (e.g. pt=[20,40,60,80])
    #   - Cutting on multiple variables simultaneously
    #---------------------------------------------------------------
    def apply_cut(self, variables_list, cuts_list, n_jets_max=100000000):

        if len(variables_list) != len(cuts_list):
            raise ValueError(f'variables_list has different length than cuts_list! {variables_list} vs. {cuts_list}')

        results_dict = self.recursive_defaultdict()
        
        # Loop through both data and MC
        for data_type in self.results.keys():

            # Loop over all cut combinations
            for variables, cuts in zip(variables_list, cuts_list):

                # Loop through all cuts in a given combination and construct mask
                total_mask = np.zeros(self.n_jets[data_type], dtype=bool) # For simplicity: cut_mask[i]=True means we will discard index i
                cut_key = ''
                for i in range(len(variables)):
                    variable = variables[i]
                    cut_min, cut_max = cuts[i]
                    cut_mask = (self.results[data_type][variable] < cut_min) | (self.results[data_type][variable] > cut_max)
                    total_mask = (total_mask | cut_mask)

                    if i>0:
                        cut_key += '_'
                    cut_key += f'{variable}{cut_min}-{cut_max}'

                # Now that we have the complete mask for this combination, apply it to all arrays
                for obs_key,result in self.results[data_type].copy().items():
                    results_dict[data_type][obs_key][cut_key] = result[~total_mask][:n_jets_max]
                if 'mc' in data_type:
                    results_dict[data_type]['total_scale_factor'][cut_key] = self.results[data_type]['total_scale_factor'].copy()[~total_mask][:n_jets_max]
                    
        return results_dict

    #---------------------------------------------------------------
    # Plot a 1D histogram from numpy array of result
    #---------------------------------------------------------------
    def plot_1D_projections(self, result_dict, pt_bin_pairs):  
        
        self.create_output_subdir(self.output_dir_plot_observables, '1D_projections')

        # Loop over observables
        for observable in self.observables:

            # Loop over subobservables
            for i in range(self.observable_info[observable]['n_subobservables']):
                obs_key = self.observable_info[observable]['obs_keys'][i]

                # Loop over pt bins
                for j,pt_bin_pair in enumerate(pt_bin_pairs):
                    pt_min = pt_bin_pair[0]
                    pt_max = pt_bin_pair[1]
                    pt_cut_key = f'jet_pt{pt_min}-{pt_max}'
                    self.plot_1D_projection(result=result_dict, obs_key=obs_key, pt_cut_key=pt_cut_key, 
                                            observable=observable, subobservable_index=i, jet_pt_index=j,
                                            pt_min=pt_min, pt_max=pt_max)

    #---------------------------------------------------------------
    # Plot a 1D histogram from numpy array of result
    #---------------------------------------------------------------
    def plot_1D_projection(self, result=None, obs_key='', pt_cut_key='', observable='', subobservable_index=0, jet_pt_index=0, pt_min=0, pt_max=0):  

        self.utils.set_plotting_options()
        ROOT.gROOT.ForceStyle()

        # Get bins -- only continue if we find bins in hepdata
        bins_list = self.observable_info[observable]['obs_bins'][subobservable_index]
        if not bins_list:
            return
        bins = bins_list[jet_pt_index]
        if not bins.any():
            return

        # Get labels
        obs_key = self.observable_info[observable]['obs_keys'][subobservable_index]
        obs_setting = self.observable_info[observable]['obs_settings'][subobservable_index]
        grooming_setting = self.observable_info[observable]['obs_grooming_settings'][subobservable_index]
        subobs_label = self.utils.formatted_subobs_label(observable)
        if grooming_setting:
            grooming_label = self.utils.formatted_grooming_label(grooming_setting)
        else:
            grooming_label = ''
        xtitle = self.observable_info[observable]['xtitle']
        ytitle = self.observable_info[observable]['ytitle']

        # Construct dict of hists for data and MC
        h_dict = {}
        for data_type in self.results.keys():

            # Create hist   
            hname = f'h_{data_type}_{obs_key}_pt{pt_min}-{pt_max}'
            #h = ROOT.TH1F(hname, hname, bins.size-1, bins)
            h = ROOT.TH1F(hname, hname, bins.size-1, array.array('d', bins))
            h.Sumw2()

            # Loop through array and fill weighted hist in order to apply pt-hat sum weights for MC
            n_jets = result[data_type][obs_key][pt_cut_key].shape[0]
            if 'mc' in data_type:
                weights = result[data_type]['total_scale_factor'][pt_cut_key]
            else:
                weights = np.ones(n_jets)

            for i in range(n_jets):
                h.Fill(result[data_type][obs_key][pt_cut_key][i], weights[i])

            # Normalize by integral, i.e. N_jets,inclusive in this pt-bin
            # First scale by bin width -- then normalize by integral
            # (where integral weights by bin width)
            h.Scale(1., 'width')
            nbins = h.GetNbinsX()
            if 'SD' in grooming_label:
                # If SD, the untagged jets are in the first bin
                min_bin = 0
                n_jets_inclusive = h.Integral(min_bin, nbins, 'width')
                n_jets_tagged = h.Integral(min_bin+1, nbins, 'width')
                f_tagging = n_jets_tagged/n_jets_inclusive
            else:
                min_bin = 1
                n_jets_inclusive = h.Integral(min_bin, nbins+1, 'width')
                f_tagging = None
            h.Scale(1./n_jets_inclusive)

            h_dict[data_type] = h

        # Generate and save plot
        self.plot_distribution_and_ratio(h_dict, obs_key=obs_key, observable=observable, 
                                         subobservable_index=subobservable_index,
                                         jet_pt_index=jet_pt_index, pt_min=pt_min, pt_max=pt_max, 
                                         subobs_label=subobs_label, obs_setting=obs_setting, 
                                         grooming_label=grooming_label, f_tagging=f_tagging,
                                         xtitle=xtitle, ytitle=ytitle)

    #-------------------------------------------------------------------------------------------
    # Plot distributions in upper panel, and ratio in lower panel
    #-------------------------------------------------------------------------------------------
    def plot_distribution_and_ratio(self, h_dict, obs_key='', observable='', subobservable_index=0, jet_pt_index=0, 
                                    pt_min=0, pt_max=0, subobs_label='',
                                    obs_setting='', grooming_label='', f_tagging=0,
                                    xtitle='', ytitle='', logy = False):

        c_name = f'c_{obs_key}_pt{pt_min}-{pt_max}'
        c = ROOT.TCanvas(c_name, c_name, 600, 650)
        c.Draw()
        c.cd()

        # Distribution
        pad2_dy = 0.45
        pad1 = ROOT.TPad('myPad', 'The pad',0,pad2_dy,1,1)
        pad1.SetLeftMargin(0.2)
        pad1.SetTopMargin(0.08)
        pad1.SetRightMargin(0.04)
        pad1.SetBottomMargin(0.)
        pad1.SetTicks(0,1)
        pad1.Draw()
        if logy:
            pad1.SetLogy()
        pad1.cd()

        legend = ROOT.TLegend(0.25,0.6,0.4,0.88)
        self.utils.setup_legend(legend, 0.05, sep=-0.1)

        legend_ratio = ROOT.TLegend(0.25,0.25,0.45,0.4)
        self.utils.setup_legend(legend_ratio, 0.06, sep=-0.1)

        bins = np.array(next(iter(h_dict.values())).GetXaxis().GetXbins())
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, bins[0], bins[-1])
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(xtitle)
        myBlankHisto.SetYTitle(ytitle)
        myBlankHisto.SetMaximum(1.7*next(iter(h_dict.values())).GetMaximum())
        myBlankHisto.SetMinimum(1e-3)
        myBlankHisto.GetYaxis().SetTitleSize(0.08)
        myBlankHisto.GetYaxis().SetTitleOffset(1.1)
        myBlankHisto.GetYaxis().SetLabelSize(0.06)
        myBlankHisto.Draw('E')

        # Ratio
        c.cd()
        pad2 = ROOT.TPad('pad2', 'pad2', 0, 0.02, 1, pad2_dy)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.21)
        pad2.SetLeftMargin(0.2)
        pad2.SetRightMargin(0.04)
        pad2.SetTicks(0,1)
        pad2.Draw()
        pad2.cd()

        myBlankHisto2 = myBlankHisto.Clone('myBlankHisto_C')
        myBlankHisto2.SetYTitle('Ratio')
        myBlankHisto2.SetXTitle(xtitle)
        myBlankHisto2.GetXaxis().SetTitleSize(26)
        myBlankHisto2.GetXaxis().SetTitleFont(43)
        myBlankHisto2.GetXaxis().SetTitleOffset(2.3)
        myBlankHisto2.GetXaxis().SetLabelFont(43)
        myBlankHisto2.GetXaxis().SetLabelSize(22)
        myBlankHisto2.GetYaxis().SetTitleSize(28)
        myBlankHisto2.GetYaxis().SetTitleFont(43)
        myBlankHisto2.GetYaxis().SetTitleOffset(2.)
        myBlankHisto2.GetYaxis().SetLabelFont(43)
        myBlankHisto2.GetYaxis().SetLabelSize(20)
        myBlankHisto2.GetYaxis().SetNdivisions(505)
        y_ratio_min = 0.
        y_ratio_max = 2.1
        myBlankHisto2.GetYaxis().SetRangeUser(y_ratio_min, y_ratio_max)
        myBlankHisto2.Draw('')

        # Draw distribution
        pad1.cd()

        # Draw HEPData
        g = self.observable_info[observable]['hepdata_tgraph'][subobservable_index][jet_pt_index]
        if g:
            g.SetMarkerColorAlpha(self.data_color, 1)
            g.SetMarkerSize(self.marker_size)
            g.SetMarkerStyle(self.data_marker)
            g.SetLineColorAlpha(self.data_color, 1)
            g.SetLineWidth(self.line_width)
            g.Draw('PE Z same')
            legend.AddEntry(g, 'Published Data', 'PE')

        # Draw our data, MC distributions
        for i,data_type in enumerate(h_dict.keys()):
            h = h_dict[data_type]
            if 'mc' in data_type:
                h.SetFillColorAlpha(self.colors[i], self.alpha)
                h.SetFillStyle(self.fillstyles[0])
                h.SetLineColorAlpha(self.colors[i], self.alpha)
                h.SetLineStyle(self.line_style)
                h.SetMarkerSize(0.)
                h.SetMarkerStyle(0)
                h.SetLineWidth(0)
                h.Draw('E3 same') # X0
                legend.AddEntry(h, self.data_type_label[data_type], 'f')
            else:
                h.SetMarkerColorAlpha(self.colors[i], self.alpha)
                h.SetMarkerSize(self.marker_size)
                h.SetMarkerStyle(self.data_marker)
                h.SetLineColorAlpha(self.colors[i], self.alpha)
                h.SetLineStyle(self.line_style)
                h.SetLineWidth(self.line_width)
                h.Draw('PE same') # X0
                legend.AddEntry(h, self.data_type_label[data_type], 'PE')

        legend.Draw()

        # Draw ratio
        pad2.cd()

        # Plot published data uncertainties at ratio=1
        g_ratio = self.divide_tgraph_by_tgraph(g, g)
        g_ratio.SetMarkerColorAlpha(self.data_color, 1)
        g_ratio.SetMarkerSize(self.marker_size)
        g_ratio.SetMarkerStyle(self.data_marker)
        g_ratio.SetLineColorAlpha(self.data_color, 1)
        g_ratio.SetLineWidth(self.line_width)
        g_ratio.Draw('PE Z same')

        # Plot ratios to published data
        for i,data_type in enumerate(h_dict.keys()):
            h = h_dict[data_type]
            g_ratio = self.divide_histogram_by_tgraph(h, g, include_tgraph_uncertainties=False)
            if g_ratio:
                if 'mc' in data_type:
                    g_ratio.SetFillColorAlpha(self.colors[i], self.alpha)
                    g_ratio.SetFillStyle(self.fillstyles[0])
                    g_ratio.SetLineColorAlpha(self.colors[i], self.alpha)
                    g_ratio.SetMarkerSize(0.)
                    g_ratio.SetMarkerStyle(0)
                    g_ratio.SetLineWidth(0)
                    g_ratio.Draw('E3 same')
                    legend_ratio.AddEntry(g_ratio, f'{self.data_type_label[data_type]} / Published', 'f')
                    legend_ratio.Draw()
                else:
                    g_ratio.SetMarkerColorAlpha(self.colors[i], self.alpha)
                    g_ratio.SetMarkerSize(self.marker_size)
                    g_ratio.SetMarkerStyle(self.data_marker)
                    g_ratio.SetLineColorAlpha(self.colors[i], self.alpha)
                    g_ratio.SetLineWidth(self.line_width)
                    g_ratio.Draw('PE Z same')
                    legend_ratio.AddEntry(g_ratio, f'{self.data_type_label[data_type]} / Published', 'PE')
                    legend_ratio.Draw()

        line = ROOT.TLine(bins[0], 1, bins[-1], 1)
        line.SetLineColor(920+2)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()

        pad1.cd()
        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        text_latex.SetTextSize(0.06);

        x = 0.65
        y = 0.85
        text = f'ALICE {self.figure_approval_status}'
        text_latex.DrawLatex(x, y, text)

        text = 'pp #sqrt{#it{s}} = 5.02 TeV'
        text_latex.SetTextSize(0.045)
        text_latex.DrawLatex(x, y-0.07, text)

        text = str(pt_min) + ' < #it{p}_{T, ch jet} < ' + str(pt_max) + ' GeV/#it{c}'
        text_latex.DrawLatex(x, y-0.14, text)

        text = '#it{R} = ' + str(self.jetR) + '   | #eta_{jet}| < 0.5'
        text_latex.DrawLatex(x, y-0.22, text)
        
        delta = 0.
        if subobs_label:
            text = f'{subobs_label} = {obs_setting}'
            text_latex.DrawLatex(x, y-0.29, text)
            delta = 0.07
        
        if grooming_label:
            text = grooming_label
            text_latex.DrawLatex(x, y-0.29-delta, text)
            
        if f_tagging:
            text = f'#it{{f}}_{{tagged}}^{{data}} ={f_tagging:.2f}'
            text_latex.DrawLatex(x, y-0.37-delta, text)

        c.SaveAs(os.path.join(self.output_dir_1D_projections, f'{obs_key}_pt{pt_min}-{pt_max}{self.file_format}'))
        c.Close()

    #---------------------------------------------------------------
    # Plot pairplot for all observables
    # We will overlay data, mc_det_matched, mc_truth_matched
    #---------------------------------------------------------------
    def plot_pairplot(self, result_dict, pt_bin_pairs):

        self.create_output_subdir(self.output_dir_plot_observables, f'2D_projections')

        # For each pt bin, construct dataframe containing all observables (except pt)

        # Plot only a subset of values for speed and image size
        n_values = 1000

        # Loop over pt bins
        dataframes_dict = self.recursive_defaultdict()
        for pt_bin_pair in pt_bin_pairs:
            pt_min = pt_bin_pair[0]
            pt_max = pt_bin_pair[1]
            pt_cut_key = f'jet_pt{pt_min}-{pt_max}'

            # Construct a separate dataframe for each data_type, then we will append them
            for data_type in result_dict.keys():
                df = pd.DataFrame()

                # Loop over data_type and observables
                for observable in self.observables:

                    # Loop over subobservables
                    for i in range(self.observable_info[observable]['n_subobservables']):
                        obs_key = self.observable_info[observable]['obs_keys'][i]
                        obs_key_new = self.obs_key_formatted(observable, i)
                        if obs_key != 'jet_pt':
                            df[obs_key_new] = result_dict[data_type][obs_key][pt_cut_key][:n_values]

                # For simplicity, remove jets with SD untagged entries
                for col in df.columns: 
                    df = df.drop(df.index[df[col] < 0])

                # Center and scale all columns of the df
                dataframes_dict[pt_cut_key][data_type] = pd.DataFrame(sklearn.preprocessing.scale(df,with_mean=True, with_std=True),
                                                                      columns = df.columns)

                # Also store the data_type
                n_jets = dataframes_dict[pt_cut_key][data_type].iloc[:,0].shape[0]
                dataframes_dict[pt_cut_key][data_type]['data_type'] = np.repeat(data_type, n_jets)

            # Append the data_types into a single dataframe
            df_total = pd.concat([dataframes_dict[pt_cut_key][data_type] for data_type in result_dict.keys()], ignore_index=True)

            # Plot scatter matrix
            print(f'  Plotting {pt_cut_key}...')
            g = sns.pairplot(df_total, corner=True, hue='data_type', plot_kws={'alpha':0.1})
            plt.tight_layout()
            g.legend.set_bbox_to_anchor((0.75, 0.75))
            #for lh in g._legend.legendHandles:   # Trying to set legend fontsize...
            #    print(dir(lh)) 
            plt.savefig(os.path.join(self.output_dir_2D_projections, f'pairplot_{pt_cut_key}.pdf'))
            plt.close()

            # Plot correlation matrix
            for data_type in result_dict.keys():

                dataframes_dict[pt_cut_key][data_type].drop(columns=['data_type'])
                corr_matrix = dataframes_dict[pt_cut_key][data_type].corr(method='pearson', numeric_only=True)
                sns.heatmap(corr_matrix)
                plt.tight_layout()
                plt.savefig(os.path.join(self.output_dir_2D_projections, f'correlation_matrix_{pt_cut_key}_{data_type}.pdf'))
                plt.close()

    #---------------------------------------------------------------
    # Return formatted observable key
    #---------------------------------------------------------------
    def obs_key_formatted(self, observable, i_subobservable, tlatex=False):

        s = self.observable_info[observable]['xtitle']

        obs_setting = self.observable_info[observable]['obs_settings'][i_subobservable]
        subobs_label_list = self.utils.formatted_subobs_label(observable)
        if subobs_label_list:
            subobs_label = subobs_label_list
            s += f' {subobs_label} = {obs_setting}'

        grooming_setting = self.observable_info[observable]['obs_grooming_settings'][i_subobservable]
        if grooming_setting:
            s += f' {self.utils.formatted_grooming_label(grooming_setting)}'

        # Convert from tlatex to latex
        if not tlatex:
            s = f'${s}$'
            s = s.replace('#it','')
            s = s.replace('} {','},\;{')
            s = s.replace('#','\\')
            s = s.replace('SD',',\;SD')
            s = s.replace(', {\\beta} = 0', '')
            s = s.replace('{\Delta R}','')
            s = s.replace('Standard_WTA','\mathrm{Standard-WTA}')
            s = s.replace('{\\lambda}_{{\\alpha}},\;{\\alpha} = ','\lambda_')

        return s

    #---------------------------------------------------------------
    # Perform unfolding
    #---------------------------------------------------------------
    def perform_unfolding(self):
        print('Perform unfolding...')

        # Create output dirs
        for systematic in self.systematics_list:
            self.create_output_subdir(self.output_dir, systematic)

    #---------------------------------------------------------------
    # Plot final results
    #---------------------------------------------------------------
    def plot_results(self):
        print('Plotting final results...')

        # Create output dir
        self.create_output_subdir(self.output_dir, 'final_results')

    #---------------------------------------------------------------
    # This function is called once for each subconfiguration
    #---------------------------------------------------------------
    def plot_single_result(self, jetR, obs_label, obs_setting, grooming_setting):    
        print(f'Plot final results for {self.observable}: R = {jetR}, {obs_label}')

    # ---------------------------------------------------------------
    # Get bin array (for each pt bin) from hepdata file
    # ---------------------------------------------------------------
    def bins_from_hepdata(self, block):

        if 'hepdata' in block:

            # Open the HEPData file
            hepdata_filename = block['hepdata']['file']        
            f = ROOT.TFile(hepdata_filename, 'READ')

            # The list of tables contains one entry per pt bin
            tables = block['hepdata']['tables'][self.jetR]
            h_name = block['hepdata']['hname']
            bins_list = []
            for table in tables:

                # Get the histogram, and return the bins
                if table:
                    dir = f.Get(table)
                    h = dir.Get(h_name)
                    bins = np.array(h.GetXaxis().GetXbins())

                    # For Soft Drop observables, we need to exclude the "untagged" bin so that it will become underflow
                    if 'SoftDrop' in block:
                        bins = bins[1:]

                else:
                    bins = np.array([])

                bins_list.append(bins)

            f.Close()

        else: 
            bins_list = None

        return bins_list

    # ---------------------------------------------------------------
    # Get tgraph (for each pt bin) from hepdata file
    # ---------------------------------------------------------------
    def tgraph_from_hepdata(self, block):

        if 'hepdata' in block:

            # Open the HEPData file
            hepdata_filename = block['hepdata']['file']        
            f = ROOT.TFile(hepdata_filename, 'READ')

            # The list of tables contains one entry per pt bin
            tables = block['hepdata']['tables'][self.jetR]
            g_name = block['hepdata']['gname']
            tgraph_list = []
            for table in tables:

                # Get the TGraph
                if table:
                    dir = f.Get(table)
                    g = dir.Get(g_name)
                else:
                    g = None

                tgraph_list.append(g)

            f.Close()

        else: 
            tgraph_list = None

        return tgraph_list

    #---------------------------------------------------------------
    # Divide a histogram by a tgraph, point-by-point
    #---------------------------------------------------------------
    def divide_histogram_by_tgraph(self, h, g, include_tgraph_uncertainties=True):

        # Truncate tgraph to range of histogram bins
        g_truncated = self.truncate_tgraph(g, h)
        if not g_truncated:
            return None

        # Clone tgraph, in order to return a new one
        g_new = g_truncated.Clone(f'{g_truncated.GetName()}_divided')

        nBins = h.GetNbinsX()
        for bin in range(1, nBins+1):

            # Get histogram (x,y)
            h_x = h.GetBinCenter(bin)
            h_y = h.GetBinContent(bin)
            h_error = h.GetBinError(bin)

            # Get TGraph (x,y) and errors
            gx, gy, yErrLow, yErrUp = self.get_gx_gy(g_truncated, bin-1)

            #print(f'h_x: {h_x}')
            #print(f'gx: {gx}')
            #print(f'h_y: {h_y}')
            #print(f'gy: {gy}')

            if not np.isclose(h_x, gx):
                print(f'WARNING: hist x: {h_x}, graph x: {gx} -- will not plot ratio')
                return None

            new_content = h_y / gy

            # Combine tgraph and histogram relative uncertainties in quadrature
            if gy > 0. and h_y > 0.:
                if include_tgraph_uncertainties:
                    new_error_low = np.sqrt( pow(yErrLow/gy,2) + pow(h_error/h_y,2) ) * new_content
                    new_error_up = np.sqrt( pow(yErrUp/gy,2) + pow(h_error/h_y,2) ) * new_content
                else:
                    new_error_low = h_error/h_y * new_content
                    new_error_up = h_error/h_y * new_content
            else:
                new_error_low = (yErrLow/gy) * new_content
                new_error_up = (yErrUp/gy) * new_content

            g_new.SetPoint(bin-1, h_x, new_content)
            g_new.SetPointError(bin-1, 0, 0, new_error_low, new_error_up)

        return g_new

    #---------------------------------------------------------------
    # Divide a tgraph by a tgraph, point-by-point: g1/g2
    # NOTE: Ignore uncertainties on denominator
    #---------------------------------------------------------------
    def divide_tgraph_by_tgraph(self, g1, g2):

        # Clone tgraph, in order to return a new one
        g_new = g1.Clone(f'{g1.GetName()}_divided')

        if g1.GetN() != g2.GetN():
            sys.exit(f'ERROR: TGraph {g1.GetName()} has {g1.GetN()} points, but {g2.GetName()} has {g2.GetN()} points')

        for i in range(0, g1.GetN()):

            # Get TGraph (x,y) and errors
            g1_x = ctypes.c_double(0)
            g1_y = ctypes.c_double(0)
            g1.GetPoint(i, g1_x, g1_y)
            y1ErrLow = g1.GetErrorYlow(i)
            y1ErrUp  = g1.GetErrorYhigh(i)
            g1x = g1_x.value
            g1y = g1_y.value

            g2_x = ctypes.c_double(0)
            g2_y = ctypes.c_double(0)
            g2.GetPoint(i, g2_x, g2_y)
            g2x = g2_x.value
            g2y = g2_y.value

            if not np.isclose(g1x, g2x):
                sys.exit(f'ERROR: TGraph {g1.GetName()} point {i} at {g1x}, but {g2.GetName()} at {g2x}')

            new_content = g1y / g2y
            new_error_low = y1ErrLow/g1y * new_content
            new_error_up = y1ErrUp/g1y * new_content

            g_new.SetPoint(i, g1x, new_content)
            g_new.SetPointError(i, 0, 0, new_error_low, new_error_up)
        return g_new

    #---------------------------------------------------------------
    # Truncate data tgraph to histogram binning range
    #---------------------------------------------------------------
    def truncate_tgraph(self, g, h):

        #print('truncate_tgraph')
        #print(h.GetName())
        #print(np.array(h.GetXaxis().GetXbins()))

        # Create new TGraph with number of points equal to number of histogram bins
        nBins = h.GetNbinsX()
        g_new = ROOT.TGraphAsymmErrors(nBins)
        g_new.SetName(f'{g.GetName()}_truncated')

        h_offset = 0
        for bin in range(1, nBins+1):

            # Get histogram (x)
            h_x = h.GetBinCenter(bin)

            # Get TGraph (x,y) and errors
            gx, gy, yErrLow, yErrUp = self.get_gx_gy(g, bin-1)

            #print(f'h_x: {h_x}')
            #print(f'gx: {gx}')
            #print(f'gy: {gy}')

            # If traph is offset from center of the bin, center it
            xErrLow = g.GetErrorXlow(bin-1)
            xErrUp = g.GetErrorXhigh(bin-1)
            if xErrLow > 0 and xErrUp > 0:
                x_min = gx - xErrLow
                x_max = gx + xErrUp
                x_center = (x_min + x_max)/2.
                if h_x > x_min and h_x < x_max:
                    if not np.isclose(gx, x_center):
                        gx = x_center

            # If tgraph starts below hist (e.g. when hist has min cut), try to get next tgraph point
            g_offset = 0
            while gx+1e-8 < h_x and g_offset < g.GetN()+1:
                g_offset += 1
                gx, gy, yErrLow, yErrUp = self.get_gx_gy(g, bin-1+g_offset)
            #print(f'new gx: {gx}')

            # If tgraph started above hist (see below) and we exhausted the tgraph points, skip
            if h_offset > 0 and np.isclose(gx, 0):
                continue

            # If tgraph starts above hist, try to get next hist bin
            h_offset = 0
            while gx-1e-8 > h_x and h_offset < nBins+1:
                h_offset += 1
                h_x = h.GetBinCenter(bin+h_offset)
                #print(f'h_x: {h_x}')
                #print(f'gx: {gx}')

            if not np.isclose(h_x, gx):
                print(f'WARNING: hist x: {h_x}, graph x: {gx}')
                return None

            g_new.SetPoint(bin-1, gx, gy)
            g_new.SetPointError(bin-1, 0, 0, yErrLow, yErrUp)
            #print()

        return g_new

    #---------------------------------------------------------------
    # Get points from tgraph by index
    #---------------------------------------------------------------
    def get_gx_gy(self, g, index):

        g_x = ctypes.c_double(0)
        g_y = ctypes.c_double(0)
        g.GetPoint(index, g_x, g_y)
        yErrLow = g.GetErrorYlow(index)
        yErrUp  = g.GetErrorYhigh(index)

        gx = g_x.value
        gy = g_y.value

        return gx, gy, yErrLow, yErrUp

    #---------------------------------------------------------------
    # Create a single output subdirectory
    #---------------------------------------------------------------
    def create_output_subdir(self, output_dir, name):

        output_subdir = os.path.join(output_dir, name)
        setattr(self, f'output_dir_{name}', output_subdir)
        if not os.path.isdir(output_subdir):
            os.makedirs(output_subdir)

        return output_subdir

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