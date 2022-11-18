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
        self.markers = [20, 21, 22, 23, 33, 34, 24, 25, 26, 32]
        self.fillstyles = [1001, 3144, 1001, 3144]

        self.data_color = ROOT.kGray+3
        self.data_marker = 21
        self.theory_marker = 20
        self.marker_size = 1.5
        self.alpha = 0.6
        self.line_width = 2
        self.line_style = 1

        # Make some labels for plotting
        self.data_type_label = {'data': 'Uncorrected Data',
                                'mc_det_matched': 'MC det-level (matched)',
                                'mc_truth_matched': 'MC truth-level (matched)',
                               }

        self.file_format = '.pdf'
        self.figure_approval_status = 'WIP'

        print(self)
        print()

    #---------------------------------------------------------------
    # Plot a 1D histogram from numpy array of result
    #---------------------------------------------------------------
    def plot_1D_projections(self, result_dict, pt_bin_pairs, observables, observable_info, jetR, output_dir):  
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
                    self.plot_1D_projection(result=result_dict, observable_info=observable_info, obs_key=obs_key, pt_cut_key=pt_cut_key, 
                                            observable=observable, subobservable_index=i, jet_pt_index=j,
                                            pt_min=pt_min, pt_max=pt_max, jetR=jetR, output_dir=output_dir)

    #---------------------------------------------------------------
    # Plot a 1D histogram from numpy array of result
    #---------------------------------------------------------------
    def plot_1D_projection(self, result=None, observable_info=None, obs_key='', pt_cut_key='', observable='', subobservable_index=0, 
                           jet_pt_index=0, pt_min=0, pt_max=0, jetR=0, output_dir=''):  

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

        # Construct dict of hists for data and MC
        h_dict = {}
        for data_type in result.keys():

            # Create hist   
            hname = f'h_{data_type}_{obs_key}_pt{pt_min}-{pt_max}'
            h = ROOT.TH1F(hname, hname, bins.size-1, bins)
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
        self.plot_distribution_and_ratio(h_dict, observable_info=observable_info, 
                                         obs_key=obs_key, observable=observable, 
                                         subobservable_index=subobservable_index,
                                         jet_pt_index=jet_pt_index, pt_min=pt_min, pt_max=pt_max, jetR=jetR,
                                         subobs_label=subobs_label, obs_setting=obs_setting, 
                                         grooming_label=grooming_label, f_tagging=f_tagging,
                                         xtitle=xtitle, ytitle=ytitle, output_dir=output_dir)

    #-------------------------------------------------------------------------------------------
    # Plot distributions in upper panel, and ratio in lower panel
    #-------------------------------------------------------------------------------------------
    def plot_distribution_and_ratio(self, h_dict, observable_info={}, obs_key='', observable='', subobservable_index=0, jet_pt_index=0, 
                                    pt_min=0, pt_max=0, jetR=0, subobs_label='',
                                    obs_setting='', grooming_label='', f_tagging=0,
                                    xtitle='', ytitle='', logy = False,
                                    output_dir=''):

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
        g = observable_info[observable]['hepdata_tgraph'][subobservable_index][jet_pt_index]
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
        g_ratio = self.utils.divide_tgraph_by_tgraph(g, g)
        g_ratio.SetMarkerColorAlpha(self.data_color, 1)
        g_ratio.SetMarkerSize(self.marker_size)
        g_ratio.SetMarkerStyle(self.data_marker)
        g_ratio.SetLineColorAlpha(self.data_color, 1)
        g_ratio.SetLineWidth(self.line_width)
        g_ratio.Draw('PE Z same')

        # Plot ratios to published data
        for i,data_type in enumerate(h_dict.keys()):
            h = h_dict[data_type]
            g_ratio = self.utils.divide_histogram_by_tgraph(h, g, include_tgraph_uncertainties=False)
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

        c.SaveAs(os.path.join(output_dir, f'{obs_key}_pt{pt_min}-{pt_max}{self.file_format}'))
        c.Close()

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