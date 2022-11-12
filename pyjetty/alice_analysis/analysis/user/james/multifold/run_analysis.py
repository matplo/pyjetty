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
import array
import ROOT

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
        
        self.colors = [ROOT.kBlue-4, ROOT.kAzure+7, ROOT.kCyan-2, ROOT.kViolet-8,
                       ROOT.kBlue-6, ROOT.kGreen+3, ROOT.kPink-4, ROOT.kRed-4, ROOT.kOrange-3]
        self.markers = [20, 21, 22, 23, 33, 34, 24, 25, 26, 32]
        self.fillstyles = [1001, 3144, 1001, 3144]

        self.data_color = ROOT.kGray+3
        self.data_marker = 21
        self.theory_marker = 20
        self.marker_size = 1.5
        self.alpha = 0.7
        self.line_width = 2
        self.line_style = 1

        print(self)

        # Load aggregated results
        # The results are stored with keys: self.observable_info[observable]['obs_key'][i]
        # i.e. f'{observable}_{self.utils.obs_label(obs_setting, grooming_setting)}'
        self.results = self.utils.read_data(self.main_data)

        self.n_jets_total = next(iter(self.results.values())).shape[0]
        print(f'Analyzing the following observables: (n_jets={self.n_jets_total})')
        for observable in self.observables:
            for i in range(self.observable_info[observable]['n_subobservables']):
                obs_key = self.observable_info[observable]['obs_keys'][i]
                print(f'  {obs_key}')
                if self.results[obs_key].shape[0] != self.n_jets_total:
                    raise ValueError(f'{obs_key} has unexpected number of jets: {self.results[obs_key].shape}')
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
            self.observable_info[observable]['obs_bins_truth'] = [np.array(obs_config_dict[obs_subconfig_list[i]]['obs_bins_truth'])
                                                                  for i in range(len(obs_subconfig_list))]
            self.observable_info[observable]['obs_bins_det'] = [np.array(obs_config_dict[obs_subconfig_list[i]]['obs_bins_det'])
                                                                for i in range(len(obs_subconfig_list))]
            self.observable_info[observable]['xtitle'] = config[observable]['common_settings']['xtitle']
            self.observable_info[observable]['ytitle'] = config[observable]['common_settings']['ytitle']
            #self.observable_info[observable]['y_min'] = [obs_config_dict[obs_subconfig_list[i]]['y_min'] for i in range(len(obs_subconfig_list))]
            #self.observable_info[observable]['y_max'] = [obs_config_dict[obs_subconfig_list[i]]['y_max'] for i in range(len(obs_subconfig_list))]

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

        # Create output dir
        self.create_output_subdir(self.output_dir, 'plot_observables')

        # First, apply pt cuts to a copy of results
        pt_bins = self.observable_info['jet_pt']['pt_bins_reported']
        pt_bin_pairs = [[pt_bins[i],pt_bins[i+1]] for i in range(len(pt_bins)-1)]
        pt_bin_pairs_nested = [[pt_bin_pair] for pt_bin_pair in pt_bin_pairs]
        variables = [['jet_pt']] * len(pt_bin_pairs_nested)
        results_dict = self.apply_cut(variables, pt_bin_pairs_nested)

        # Print number of jets passing each cut
        n_jets = [results_dict['jet_pt'][f'jet_pt{x[0]}-{x[1]}'].shape[0] for x in pt_bin_pairs]
        [print(f'pt={pt_bin_pairs[i]} has n_jets={n_jets[i]}') for i in range(len(n_jets))]

        # Plot 1D projections of each observable
        print()
        print(f'Plotting 1D projections...')
        self.plot_1D_projections(results_dict, pt_bin_pairs)

        # Plot pairplot between each pair of observables
        print()
        print(f'Plotting pairplot...')
        self.plot_pairplot(results_dict, pt_bin_pairs)

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
            total_mask = np.zeros(self.n_jets_total, dtype=bool) # For simplicity: cut_mask[i]=True means we will discard index i
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
                results_dict[obs_key][cut_key] = result[~total_mask]
                    
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
                for pt_bin_pair in pt_bin_pairs:
                    pt_min = pt_bin_pair[0]
                    pt_max = pt_bin_pair[1]
                    pt_cut_key = f'jet_pt{pt_min}-{pt_max}'
                    self.plot_1D_projection(result=result_dict[obs_key][pt_cut_key], observable=observable, subobservable_index=i, 
                                            pt_min=pt_min, pt_max=pt_max)

    #---------------------------------------------------------------
    # Plot a 1D histogram from numpy array of result
    #---------------------------------------------------------------
    def plot_1D_projection(self, result=None, observable='', subobservable_index=0, pt_min=0, pt_max=0):  

        self.utils.set_plotting_options()
        ROOT.gROOT.ForceStyle()

        # Get bins and labels
        bins = self.observable_info[observable]['obs_bins_det'][subobservable_index]
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
        #y_min = self.observable_info[observable]['y_min'][subobservable_index]
        #y_max = self.observable_info[observable]['y_max'][subobservable_index]

        # Create hist
        
        hname = f'h_{obs_key}_pt{pt_min}-{pt_max}'
        #h = ROOT.TH1F(hname, hname, bins.size-1, bins)
        h = ROOT.TH1F(hname, hname, bins.size-1, array.array('d', bins))
        h.Sumw2()

        # Loop through array and fill weighted hist
        n_jets = result.shape[0]
        weights = np.ones(n_jets)
        for i in range(n_jets):
            h.Fill(result[i], weights[i])

        # TODO: Scale by pt-hat sum weights?

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

        # Plot HEPData
        # Load them into class members at initialization?
        # Can use plot_utils.tgraph_from_hepdata(), divide_histogram_by_tgraph() from jetscape plot script

        # Generate and save plot
        self.plot_distribution_and_ratio(h, pt_min=pt_min, pt_max=pt_max, 
                                         subobs_label=subobs_label, obs_setting=obs_setting, 
                                         grooming_label=grooming_label, f_tagging=f_tagging,
                                         xtitle=xtitle, ytitle=ytitle)

    #-------------------------------------------------------------------------------------------
    # Plot distributions in upper panel, and ratio in lower panel
    #-------------------------------------------------------------------------------------------
    def plot_distribution_and_ratio(self, h, pt_min=0, pt_max=0, subobs_label='',
                                    obs_setting='', grooming_label='', f_tagging=0,
                                    xtitle='', ytitle='', logy = False):
        c = ROOT.TCanvas(f'c_{h.GetName()}', f'c_{h.GetName()}', 600, 650)
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

        legend = ROOT.TLegend(0.25,0.8,0.4,0.9)
        self.utils.setup_legend(legend, 0.05, sep=-0.1)

        legend_ratio = ROOT.TLegend(0.25,0.8,0.45,1.07)
        self.utils.setup_legend(legend_ratio, 0.07, sep=-0.1)

        bins = np.array(h.GetXaxis().GetXbins())
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, bins[0], bins[-1])
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(xtitle)
        myBlankHisto.SetYTitle(ytitle)
        myBlankHisto.SetMaximum(1.7*h.GetMaximum())
        myBlankHisto.SetMinimum(0.)
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
        y_ratio_max = 2.
        myBlankHisto2.GetYaxis().SetRangeUser(y_ratio_min, y_ratio_max)
        myBlankHisto2.Draw('')

        # Draw distribution
        pad1.cd()
        h.SetMarkerSize(self.marker_size)
        h.SetMarkerStyle(self.data_marker)
        h.SetMarkerColor(self.data_color)
        h.SetLineStyle(self.line_style)
        h.SetLineWidth(self.line_width)
        h.SetLineColor(self.data_color)
        h.Draw('PE Z same')
        legend.AddEntry(h, 'Uncorrected Data', 'PE')

        # Draw HEPData
        #self.output_dict[f'jetscape_distribution_{label}'] =  self.observable_settings[f'jetscape_distribution']
        #if self.observable_settings[f'jetscape_distribution'].GetNbinsX() > 1:
        #    self.observable_settings[f'jetscape_distribution'].SetFillColor(self.jetscape_color[0])
        #    self.observable_settings[f'jetscape_distribution'].SetFillColorAlpha(self.jetscape_color[0], self.jetscape_alpha[0])
        #    self.observable_settings[f'jetscape_distribution'].SetFillStyle(1001)
        #    self.observable_settings[f'jetscape_distribution'].SetMarkerSize(0.)
        #    self.observable_settings[f'jetscape_distribution'].SetMarkerStyle(0)
        #    self.observable_settings[f'jetscape_distribution'].SetLineWidth(0)
        #    self.observable_settings[f'jetscape_distribution'].DrawCopy('E3 same')
        #elif self.observable_settings[f'jetscape_distribution'].GetNbinsX() == 1:
        #    self.observable_settings[f'jetscape_distribution'].SetMarkerSize(self.marker_size)
        #    self.observable_settings[f'jetscape_distribution'].SetMarkerStyle(self.data_marker+1)
        #    self.observable_settings[f'jetscape_distribution'].SetMarkerColor(self.jetscape_color[0])
        #    self.observable_settings[f'jetscape_distribution'].SetLineStyle(self.line_style)
        #    self.observable_settings[f'jetscape_distribution'].SetLineWidth(self.line_width)
        #    self.observable_settings[f'jetscape_distribution'].SetLineColor(self.jetscape_color[0])
        #    self.observable_settings[f'jetscape_distribution'].DrawCopy('PE same')
        #legend.AddEntry(self.observable_settings[f'jetscape_distribution'], 'JETSCAPE', 'f')

        legend.Draw()

        # Draw ratio
        pad2.cd()
        #if self.observable_settings['data_distribution']:
        #    data_ratio = self.utils.divide_tgraph_by_tgraph(self.observable_settings['data_distribution'],
        #                                                         self.observable_settings['data_distribution'])
        #    data_ratio.Draw('PE Z same')
        #    self.output_dict[f'data_ratio_{label}'] = data_ratio
        #if self.observable_settings['ratio']:
        #    self.output_dict[f'ratio_{label}'] = self.observable_settings['ratio']
        #    if self.observable_settings['ratio'].GetN() > 1:
        #        self.observable_settings['ratio'].SetFillColor(self.jetscape_color[0])
        #        self.observable_settings['ratio'].SetFillColorAlpha(self.jetscape_color[0], self.jetscape_alpha[0])
        #        self.observable_settings['ratio'].SetFillStyle(1001)
        #        self.observable_settings['ratio'].SetMarkerSize(0.)
        #        self.observable_settings['ratio'].SetMarkerStyle(0)
        #        self.observable_settings['ratio'].SetLineWidth(0)
        #        self.observable_settings['ratio'].Draw('E3 same')
        #    elif self.observable_settings['ratio'].GetN() == 1:
        #        self.observable_settings['ratio'].SetMarkerSize(self.marker_size)
        #        self.observable_settings['ratio'].SetMarkerStyle(self.data_marker+1)
        #        self.observable_settings['ratio'].SetMarkerColor(self.jetscape_color[0])
        #        self.observable_settings['ratio'].SetLineStyle(self.line_style)
        #        self.observable_settings['ratio'].SetLineWidth(self.line_width)
        #        self.observable_settings['ratio'].SetLineColor(self.jetscape_color[0])
        #        self.observable_settings['ratio'].Draw('PE same')

        #if self.observable_settings['ratio']:
        #    legend_ratio.AddEntry(self.observable_settings['ratio'], 'JETSCAPE/Data', 'f')
        #if self.observable_settings['data_distribution']:
        #    legend_ratio.AddEntry(data_ratio, 'Data uncertainties', 'PE')
        #legend_ratio.Draw()

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
            text = '{} = {}'.format(subobs_label, obs_setting)
            text_latex.DrawLatex(x, y-0.29, text)
            delta = 0.07
        
        if grooming_label:
            text = grooming_label
            text_latex.DrawLatex(x, y-0.29-delta, text)
            
        if f_tagging:
            text = f'#it{{f}}_{{tagged}}^{{data}} ={f_tagging:.2f}'
            text_latex.DrawLatex(x, y-0.37-delta, text)

        c.SaveAs(os.path.join(self.output_dir_1D_projections, f'{h.GetName()}{self.file_format}'))
        c.Close()

    #---------------------------------------------------------------
    # Plot pairplot for all observables
    #---------------------------------------------------------------
    def plot_pairplot(self, result_dict, pt_bin_pairs):
        
        # Build a dataframe of all observables
        df = ''  

        # Separate PYTHIA/JEWEL
        #jewel_indices = y_train
        #pythia_indices = 1 - y_train
        #n_plot = int(self.n_train) # Plot a subset to save time/space
        #X_Nsub_jewel = X_train[jewel_indices.astype(bool)][:n_plot]
        #X_Nsub_pythia = X_train[pythia_indices.astype(bool)][:n_plot]

        ## Construct dataframes for scatter matrix plotting
        #df_jewel = pd.DataFrame(X_Nsub_jewel, columns=feature_labels)
        #df_pythia = pd.DataFrame(X_Nsub_pythia, columns=feature_labels)
        
        ## Add label columns to each df to differentiate them for plotting
        #df_jewel['generator'] = np.repeat(self.AA_label, X_Nsub_jewel.shape[0])
        #df_pythia['generator'] = np.repeat(self.pp_label, X_Nsub_pythia.shape[0])
        #df = df_jewel.append(df_pythia, ignore_index=True)

        ## Plot scatter matrix
        #g = sns.pairplot(df, corner=True, hue='generator', plot_kws={'alpha':0.1})
        ##g.legend.set_bbox_to_anchor((0.75, 0.75))
        ##plt.savefig(os.path.join(self.output_dir_i, f'training_data_scatter_matrix_K{K}.png'), dpi=50)
        #plt.savefig(os.path.join(self.output_dir_i, f'training_data_scatter_matrix_K{K}_{suffix}.pdf'))
        #plt.close()
        
        ## Plot correlation matrix
        #df.drop(columns=['generator'])
        #corr_matrix = df.corr()
        #sns.heatmap(corr_matrix)
        #plt.savefig(os.path.join(self.output_dir_i, f'training_data_correlation_matrix_K{K}_{suffix}.pdf'))
        #plt.close()


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

        name = 'cResult_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt_truth, max_pt_truth)
        c = ROOT.TCanvas(name, name, 600, 450)
        c.Draw()
        
        c.cd()
        myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
        myPad.SetLeftMargin(0.2)
        myPad.SetTopMargin(0.07)
        myPad.SetRightMargin(0.04)
        myPad.SetBottomMargin(0.13)
        myPad.Draw()
        myPad.cd()
        
        xtitle = getattr(self, 'xtitle')
        ytitle = getattr(self, 'ytitle')
        color = 600-6
        
        # Get histograms
        name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
        h = getattr(self, name)
        h.SetName(name)
        h.SetMarkerSize(1.5)
        h.SetMarkerStyle(20)
        h.SetMarkerColor(color)
        h.SetLineStyle(1)
        h.SetLineWidth(2)
        h.SetLineColor(color)
        
        name = 'hResult_{}_systotal_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
        h_sys = getattr(self, name)
        h_sys.SetName(name)
        h_sys.SetLineColor(0)
        h_sys.SetFillColor(color)
        h_sys.SetFillColorAlpha(color, 0.3)
        h_sys.SetFillStyle(1001)
        h_sys.SetLineWidth(0)
        
        n_obs_bins_truth = self.n_bins_truth(obs_label)
        truth_bin_array = self.truth_bin_array(obs_label)
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_obs_bins_truth, truth_bin_array)
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(xtitle)
        myBlankHisto.GetYaxis().SetTitleOffset(1.5)
        myBlankHisto.SetYTitle(ytitle)
        myBlankHisto.SetMaximum(1.7*h.GetMaximum())
        myBlankHisto.SetMinimum(0.)
        myBlankHisto.Draw("E")

        if plot_pythia:
        
            hPythia, fraction_tagged_pythia = self.pythia_prediction(jetR, obs_setting, grooming_setting, obs_label, min_pt_truth, max_pt_truth)
            if hPythia:
                hPythia.SetFillStyle(0)
                hPythia.SetMarkerSize(1.5)
                hPythia.SetMarkerStyle(21)
                hPythia.SetMarkerColor(1)
                hPythia.SetLineColor(1)
                hPythia.SetLineWidth(1)
                hPythia.Draw('E2 same')
            else:
                print('No PYTHIA prediction for {} {}'.format(self.observable, obs_label))
                plot_pythia = False
            
        h_sys.DrawCopy('E2 same')
        h.DrawCopy('PE X0 same')
    
        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        text = 'ALICE {}'.format(self.figure_approval_status)
        text_latex.DrawLatex(0.57, 0.87, text)
        
        text = 'pp #sqrt{#it{s}} = 5.02 TeV'
        text_latex.SetTextSize(0.045)
        text_latex.DrawLatex(0.57, 0.8, text)

        text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
        text_latex.DrawLatex(0.57, 0.73, text)

        text = '#it{R} = ' + str(jetR) + '   | #eta_{jet}| < 0.5'
        text_latex.DrawLatex(0.57, 0.66, text)
        
        subobs_label = self.utils.formatted_subobs_label(self.observable)
        delta = 0.
        if subobs_label:
            text = '{} = {}'.format(subobs_label, obs_setting)
            text_latex.DrawLatex(0.57, 0.59, text)
            delta = 0.07
        
        if grooming_setting:
            text = self.utils.formatted_grooming_label(grooming_setting)
            text_latex.DrawLatex(0.57, 0.59-delta, text)
        
        if 'sd' in grooming_setting:
            fraction_tagged = getattr(self, 'tagging_fraction_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt_truth, max_pt_truth))
            text_latex.SetTextSize(0.04)
            text = '#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged
            text_latex.DrawLatex(0.57, 0.52-delta, text)
        
            if plot_pythia:
                text_latex.SetTextSize(0.04)
                text = ('#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged) + (', #it{f}_{tagged}^{pythia} = %3.3f' % fraction_tagged_pythia)
                text_latex.DrawLatex(0.57, 0.52-delta, text)

        myLegend = ROOT.TLegend(0.25,0.7,0.45,0.85)
        self.utils.setup_legend(myLegend,0.035)
        myLegend.AddEntry(h, 'ALICE pp', 'pe')
        myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
        if plot_pythia:
            myLegend.AddEntry(hPythia, 'PYTHIA8 Monash2013', 'pe')
        myLegend.Draw()
            
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