#! /usr/bin/env python

'''
Plotting class.

Author: James Mulligan (james.mulligan@berkeley.edu)
'''

import os
import numpy as np
import pandas as pd

import sklearn.preprocessing

try:
    import ROOT
except ImportError:
    pass

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_context('paper', rc={'font.size':18, 'axes.titlesize':18, 'axes.labelsize':18, 'text.usetex':True})

# Suppress some annoying warnings
np.finfo(np.dtype("float32")) 
np.finfo(np.dtype("float64"))

from pyjetty.alice_analysis.analysis.base import common_base
from pyjetty.alice_analysis.analysis.user.james.multifold import multifold_utils

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class MultiFoldPlotUtils(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, **kwargs):
        super(MultiFoldPlotUtils, self).__init__(**kwargs)

        # Initialize utils class
        self.utils = multifold_utils.MultiFoldUtils()
        self.utils.set_plotting_options()
        ROOT.gROOT.ForceStyle()
        
        self.colors = [ROOT.kBlue-4, ROOT.kGreen-2, ROOT.kRed-4, 
                       ROOT.kAzure+7, ROOT.kCyan-2, ROOT.kViolet-8,
                       ROOT.kBlue-6, ROOT.kPink-4, ROOT.kOrange-3]
        self.markers = [20, 22, 23, 33, 34, 24, 25, 26, 32]
        self.fillstyles = [1001, 3144, 1001, 3144]

        self.data_color = ROOT.kGray+3
        self.data_marker = 21
        self.marker_size = 1.5
        self.alpha = 0.6
        self.line_width = 2
        self.line_style = 1

        self.file_format = '.pdf'
        self.figure_approval_status = 'WIP'

        print()
        print(self)
        print()

    #---------------------------------------------------------------
    # Plot a 1D histogram from numpy array of result
    #---------------------------------------------------------------
    def plot_1D_projections_before_unfolding(self, result_dict, pt_bin_pairs, observables, observable_info, jetR, output_dir):  
        output_dir = self.utils.create_output_subdir(output_dir, '1D_projections')
        
        # Loop over observables
        for observable in observables:

            # Loop over subobservables
            for i in range(observable_info[observable]['n_subobservables']):
                obs_key = observable_info[observable]['obs_keys'][i]

                # Loop over pt bins
                for j,pt_bin_pair in enumerate(pt_bin_pairs):
                    pt_min = pt_bin_pair[0]
                    pt_max = pt_bin_pair[1]
                    pt_cut_key = f'jet_pt{pt_min}-{pt_max}'

                    # Set which results to plot -- data,mc_det_matched,mc_truth_matched 
                    #                              (and published hepdata will be plotted automatically)
                    h_keys = [[data_type,obs_key] for data_type in result_dict.keys()]

                    data_type_label = {'data': 'Uncorrected Data',
                                       'mc_det_matched': 'MC det-level (matched)',
                                       'mc_truth_matched': 'MC truth-level (matched)',
                                      }
                    h_legend_labels = [data_type_label[data_type] for data_type in result_dict.keys()]
                    h_plotmarkers = ['mc' not in data_type for data_type in result_dict.keys()]
                    weights = [result_dict[data_type]['pt_hat_scale_factors'][pt_cut_key] if 'mc' in data_type else np.array([]) for data_type in result_dict.keys()]

                    self.plot_1D_projection(result_dict, h_keys, h_legend_labels, h_plotmarkers, weights,
                                            observable_info=observable_info, obs_key=obs_key, cut_key=pt_cut_key, 
                                            observable=observable, subobservable_index=i, jet_pt_index=j,
                                            pt_min=pt_min, pt_max=pt_max, jetR=jetR, output_dir=output_dir)

    #---------------------------------------------------------------
    # Plot pairplot for all observables
    # We will overlay data, mc_det_matched, mc_truth_matched
    #---------------------------------------------------------------
    def plot_pairplot(self, result_dict, pt_bin_pairs, observables, observable_info, output_dir):
        output_dir = self.utils.create_output_subdir(output_dir, '2D_projections')

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
                for observable in observables:

                    # Loop over subobservables
                    for i in range(observable_info[observable]['n_subobservables']):
                        obs_key = observable_info[observable]['obs_keys'][i]
                        obs_key_new = self.utils.obs_key_formatted(observable_info, observable, i)
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
            plt.savefig(os.path.join(output_dir, f'pairplot_{pt_cut_key}.pdf'))
            plt.close()

            # Plot correlation matrix
            for data_type in result_dict.keys():

                dataframes_dict[pt_cut_key][data_type].drop(columns=['data_type'])
                corr_matrix = dataframes_dict[pt_cut_key][data_type].corr(method='pearson', numeric_only=True)
                sns.heatmap(corr_matrix)
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, f'correlation_matrix_{pt_cut_key}_{data_type}.pdf'))
                plt.close()

    #---------------------------------------------------------------
    # Plot unfolding results as a function of iterations
    #
    # unfolding_results contains:
    #   'nature_det' -- experimental data 
    #       - [obs_key]
    #   'sim_det' -- det-level simulation
    #       - [obs_key]
    #       - [weights_det_iteration{i}] -- weights that map sim_det to nature_det
    #       - [pt_hat_scale_factors]
    #   'sim_truth' -- truth-level simulation
    #       - [obs_key]
    #       - [weights_truth_iteration{i}] -- weights that map sim_truth to nature_truth
    #       - [pt_hat_scale_factors]
    #---------------------------------------------------------------
    def plot_unfolding_results(self, unfolding_results, pt_bin_pairs, observables, observable_info, jetR, output_dir):
        output_dir_truth = self.utils.create_output_subdir(output_dir, 'Truth')
        output_dir_det = self.utils.create_output_subdir(output_dir, 'Det')

        iterations = [key[-1] for key in unfolding_results['sim_truth'].keys() if 'weights_truth' in key]

        # Loop over observables
        for observable in observables:

            # Loop over subobservables
            for i in range(observable_info[observable]['n_subobservables']):
                obs_key = observable_info[observable]['obs_keys'][i]

                # Loop over pt bins
                for j,pt_bin_pair in enumerate(pt_bin_pairs):
                    pt_min = pt_bin_pair[0]
                    pt_max = pt_bin_pair[1]
                    pt_cut_key = f'jet_pt{pt_min}-{pt_max}'

                    # TODO: combine unfolding weights and pt-hat-scale-factors?
                    #       or combine them already during the unfolding?

                    # ---------------------          
                    # Set which results to plot -- unfolded result
                    #  (published hepdata will be plotted automatically)
                    h_keys, h_legend_labels, h_plotmarkers, weights = [], [], [], []
                    # Plot unfolded iterations: sim_truth weighted by weights_truth 
                    for iteration in iterations:
                        data_type = 'sim_truth'
                        label = f'Unfolded, iteration {iteration}'
                        h_keys.append([data_type,obs_key])
                        h_legend_labels.append(label)
                        h_plotmarkers.append(False)
                        weights.append(unfolding_results[data_type][f'weights_truth_iteration{iteration}'][pt_cut_key])
                    # Plot original data
                    data_type = 'nature_det'
                    h_keys.append([data_type,obs_key])
                    h_legend_labels.append('Uncorrected Data')
                    h_plotmarkers.append(True)
                    weights.append(np.array([]))

                    self.plot_1D_projection(unfolding_results, h_keys, h_legend_labels, h_plotmarkers, weights,
                                            observable_info=observable_info, obs_key=obs_key, cut_key=pt_cut_key, 
                                            observable=observable, subobservable_index=i, jet_pt_index=j,
                                            pt_min=pt_min, pt_max=pt_max, jetR=jetR, output_dir=output_dir_truth, prefix='Truth_')

                    # ---------------------          
                    # Set which results to plot -- reweighted det-level 
                    #  (published hepdata will be plotted automatically)
                    h_keys, h_legend_labels, h_plotmarkers, weights = [], [], [], []
                    # Plot reweighted det-level iterations: sim_det weighted by weights_det
                    for iteration in iterations:
                        if iteration == '0':
                            continue
                        data_type = 'sim_det'
                        label = f'MC det-level, iteration {iteration}'
                        h_keys.append([data_type,obs_key])
                        h_legend_labels.append(label)
                        h_plotmarkers.append(False)
                        weights.append(unfolding_results[data_type][f'weights_det_iteration{iteration}'][pt_cut_key])
                    # Plot original data
                    data_type = 'nature_det'
                    h_keys.append([data_type,obs_key])
                    h_legend_labels.append('Uncorrected Data')
                    h_plotmarkers.append(True)
                    weights.append(np.array([]))

                    self.plot_1D_projection(unfolding_results, h_keys, h_legend_labels, h_plotmarkers, weights,
                                            observable_info=observable_info, obs_key=obs_key, cut_key=pt_cut_key, 
                                            observable=observable, subobservable_index=i, jet_pt_index=j,
                                            pt_min=pt_min, pt_max=pt_max, jetR=jetR, output_dir=output_dir_det, prefix='Det_')

    #---------------------------------------------------------------
    # Plot a 1D histogram for a given obs_key and pt_cut from a dict of numpy arrays
    #
    # We want to support flexibility to e.g.:
    #  (1) Plot MC-truth/MC-det/data before unfolding
    #  (2) Plot Reweighted det-level vs. data after unfolding, for multiple iterations
    # This function should therefore be flexible in the inputs that it takes. 
    #
    # We take the following inputs:
    #   - result: dict of numpy arrays
    #       result[data_type][obs_key][cut_key] = 1darray
    #   - h_keys: (ordered) list of histogram keys to plot
    #       h_keys = [[data_type1, obs_key1], [data_type1, obs_key2], ...]
    #   - h_weights: (ordered) list of 1darray's of weights
    #   - h_legend_labels: (ordered) list of histogram legend labels
    #   - h_plotmarkers: (ordered) list of whether to plot markers or fill, e.g. [1, 1, 0, 0, ...]
    #---------------------------------------------------------------
    def plot_1D_projection(self, result, h_key_pairs, h_legend_labels, h_plotmarkers, weights,
                           observable_info=None, obs_key='', cut_key='', observable='', subobservable_index=0, 
                           jet_pt_index=0, pt_min=0, pt_max=0, jetR=0, output_dir='', prefix=''):  

        # Get bins -- only continue if we find bins in hepdata
        bins_list = observable_info[observable]['obs_bins'][subobservable_index]
        if not bins_list:
            return
        bins = bins_list[jet_pt_index]
        if not bins.any():
            return

        # Get labels
        obs_key = observable_info[observable]['obs_keys'][subobservable_index]
        obs_setting = observable_info[observable]['obs_settings'][subobservable_index]
        grooming_setting = observable_info[observable]['obs_grooming_settings'][subobservable_index]
        subobs_label = self.utils.formatted_subobs_label(observable)
        if grooming_setting:
            grooming_label = self.utils.formatted_grooming_label(grooming_setting)
        else:
            grooming_label = ''
        xtitle = observable_info[observable]['xtitle']
        ytitle = observable_info[observable]['ytitle']

        # Construct dict of hists that we will plot
        h_dict = {}
        h_keys = []
        for i_hist,keys in enumerate(h_key_pairs):
            data_type, obs_key = keys
            
            # Create hist   
            hname = f'h_{h_legend_labels[i_hist]}_{obs_key}_pt{pt_min}-{pt_max}{prefix}'
            h = ROOT.TH1F(hname, hname, bins.size-1, bins)
            h.Sumw2()

            # Loop through array and fill weighted hist in order to apply pt-hat sum weights for MC            
            n_jets = result[data_type][obs_key][cut_key].shape[0]
            if weights[i_hist].any():
                for i_jet in range(n_jets):
                    h.Fill(result[data_type][obs_key][cut_key][i_jet], weights[i_hist][i_jet])
            else:
                for i_jet in range(n_jets):
                    h.Fill(result[data_type][obs_key][cut_key][i_jet])

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

            h_dict[f'{data_type}_{h_legend_labels[i_hist]}'] = h
            h_keys.append(f'{data_type}_{h_legend_labels[i_hist]}')

        # Generate and save plot
        self.plot_distribution_and_ratio(h_dict, h_keys, h_legend_labels, h_plotmarkers,
                                         hepdata_tgraph=observable_info[observable]['hepdata_tgraph'][subobservable_index][jet_pt_index],
                                         obs_key=obs_key, pt_min=pt_min, pt_max=pt_max, jetR=jetR,
                                         subobs_label=subobs_label, obs_setting=obs_setting, 
                                         grooming_label=grooming_label, f_tagging=f_tagging,
                                         xtitle=xtitle, ytitle=ytitle, 
                                         output_dir=output_dir, prefix=prefix)

    #-------------------------------------------------------------------------------------------
    # Plot a set of histograms
    #   - In upper panel: Plot histograms as well as hepdata tgraph
    #   - In low panel: ratio of all histograms to hepdata tgraph
    #
    # The main parameters specify the data to be plotted
    #   - h_dict: dictionary containing histograms to be plotted
    #   - h_keys: (ordered) list of histogram keys to plot
    #   - h_legend_labels: (ordered) list of histogram legend labels
    #   - h_plotmarkers: (ordered) list of whether to plot markers or fill, e.g. [1, 1, 0, 0, ...]
    #   - hepdata_tgraph (optional): tgraph of published data
    # 
    # The rest of the parameters are for book-keeping/plotting the observable, kinematics, etc.
    #  - obs_key: label specifying the observable+subsobservable+grooming_settings
    #  - pt_min,pt_max
    #  - jetR
    #  - obs_setting, grooming_label
    #  - f_tagging
    #  - ...
    #-------------------------------------------------------------------------------------------
    def plot_distribution_and_ratio(self, h_dict, h_keys, h_legend_labels, h_plotmarkers,
                                    hepdata_tgraph = None,
                                    obs_key='', pt_min=0, pt_max=0, jetR=0, 
                                    subobs_label='', obs_setting='', grooming_label='', f_tagging=0,
                                    xtitle='', ytitle='', logy = False,
                                    output_dir='', prefix=''):

        if len(h_keys) != len(h_legend_labels) or len(h_keys) != len(h_plotmarkers):
            raise ValueError(f'Length of inputs is not the same! h_keys={h_keys}, h_legend_labels={h_legend_labels}, h_plotmarkers={h_plotmarkers}')

        c_name = f'c_{obs_key}_pt{pt_min}-{pt_max}{prefix}'
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
        self.utils.setup_legend(legend_ratio, 0.05, sep=-0.1)

        bins = np.array(next(iter(h_dict.values())).GetXaxis().GetXbins())
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, bins[0], bins[-1])
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(xtitle)
        myBlankHisto.SetYTitle(ytitle)
        myBlankHisto.SetMaximum(2.*next(iter(h_dict.values())).GetMaximum())
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
        if hepdata_tgraph:
            hepdata_tgraph.SetMarkerColorAlpha(self.data_color, 1)
            hepdata_tgraph.SetMarkerSize(self.marker_size)
            hepdata_tgraph.SetMarkerStyle(self.data_marker)
            hepdata_tgraph.SetLineColorAlpha(self.data_color, 1)
            hepdata_tgraph.SetLineWidth(self.line_width)
            hepdata_tgraph.Draw('PE Z same')
            legend.AddEntry(hepdata_tgraph, 'Published Data', 'PE')

        # Draw our data, MC distributions
        i_marker = 0
        for i,h_key in enumerate(h_keys):
            h = h_dict[h_key]
            if h_plotmarkers[i]:
                h.SetMarkerColorAlpha(self.colors[i], self.alpha)
                h.SetMarkerSize(self.marker_size)
                h.SetMarkerStyle(self.markers[i_marker])
                h.SetLineColorAlpha(self.colors[i], self.alpha)
                h.SetLineStyle(self.line_style)
                h.SetLineWidth(self.line_width)
                h.Draw('PE same') # X0
                legend.AddEntry(h, h_legend_labels[i], 'PE')
                i_marker += 1
            else:
                h.SetFillColorAlpha(self.colors[i], self.alpha)
                h.SetFillStyle(self.fillstyles[0])
                h.SetLineColorAlpha(self.colors[i], self.alpha)
                h.SetLineStyle(self.line_style)
                h.SetMarkerSize(0.)
                h.SetMarkerStyle(0)
                h.SetLineWidth(0)
                h.Draw('E3 same') # X0
                legend.AddEntry(h, h_legend_labels[i], 'f')

        legend.Draw()

        # Draw ratio
        pad2.cd()

        # Plot published data uncertainties at ratio=1
        g_ratio = self.utils.divide_tgraph_by_tgraph(hepdata_tgraph, hepdata_tgraph)
        g_ratio.SetMarkerColorAlpha(self.data_color, 1)
        g_ratio.SetMarkerSize(self.marker_size)
        g_ratio.SetMarkerStyle(self.data_marker)
        g_ratio.SetLineColorAlpha(self.data_color, 1)
        g_ratio.SetLineWidth(self.line_width)
        g_ratio.Draw('PE Z same')

        # Plot ratios to published data
        i_marker = 0
        for i,h_key in enumerate(h_keys):
            h = h_dict[h_key]
            g_ratio = self.utils.divide_histogram_by_tgraph(h, hepdata_tgraph, include_tgraph_uncertainties=False)
            if g_ratio:
                if h_plotmarkers[i]:
                    g_ratio.SetMarkerColorAlpha(self.colors[i], self.alpha)
                    g_ratio.SetMarkerSize(self.marker_size)
                    g_ratio.SetMarkerStyle(self.markers[i_marker])
                    g_ratio.SetLineColorAlpha(self.colors[i], self.alpha)
                    g_ratio.SetLineWidth(self.line_width)
                    g_ratio.Draw('PE Z same')
                    legend_ratio.AddEntry(g_ratio, f'{h_legend_labels[i]} / Published', 'PE')
                    legend_ratio.Draw()
                    i_marker += 1
                else:
                    g_ratio.SetFillColorAlpha(self.colors[i], self.alpha)
                    g_ratio.SetFillStyle(self.fillstyles[0])
                    g_ratio.SetLineColorAlpha(self.colors[i], self.alpha)
                    g_ratio.SetMarkerSize(0.)
                    g_ratio.SetMarkerStyle(0)
                    g_ratio.SetLineWidth(0)
                    g_ratio.Draw('E3 same')
                    legend_ratio.AddEntry(g_ratio, f'{h_legend_labels[i]} / Published', 'f')
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

        text = '#it{R} = ' + str(jetR) + '   | #eta_{jet}| < 0.5'
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

        c.SaveAs(os.path.join(output_dir, f'{prefix}{obs_key}_pt{pt_min}-{pt_max}{self.file_format}'))
        c.Close()