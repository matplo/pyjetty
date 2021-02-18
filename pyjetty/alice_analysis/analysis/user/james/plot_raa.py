#! /usr/bin/env python

import sys
import os
import argparse
from array import *
import numpy as np
import ROOT
import yaml

from pyjetty.alice_analysis.process.base import common_base
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class PlotRAA(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, **kwargs):
        super(PlotRAA, self).__init__(**kwargs)

        self.utils = analysis_utils_obs.AnalysisUtils_Obs()
            
        self.initialize()

    #---------------------------------------------------------------
    # Initialize config file into class members
    #---------------------------------------------------------------
    def initialize(self):

        self.output_dir = './PRL_results'
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        # Load config file
        config_file = './PRL_results/plot.yaml'
        with open(config_file, 'r') as stream:
            config = yaml.safe_load(stream)

        self.config_james = config['james']
        self.config_james_results = [self.config_james[name] for name in self.config_james if 'result' in name ]

        self.config_laura = config['laura']
        self.config_laura_results = [self.config_james[name] for name in self.config_laura if 'result' in name ]

        self.plot_data = True
        self.plot_theory = True

        self.colors = [600-6, 632-4]
        self.ratio_color = ROOT.kGray+3
        self.markers = [20, 21, 33]
        self.theory_colors = [ROOT.kViolet-8, ROOT.kAzure-4, ROOT.kTeal-8, ROOT.kOrange+6, ROOT.kRed-7, ROOT.kPink+1, ROOT.kCyan-2, ROOT.kBlue-10]

        ROOT.gStyle.SetLineStyleString(11,'30 12');
        self.line_style = [1, 1, 1, 1, 1, 1, 11, 1, 1]
        self.line_width = [4, 4, 4, 4, 4, 4, 6, 4, 4]
         
        self.formatted_grooming_label = 'Soft Drop #it{{z}}_{{cut}}={}, #it{{#beta}}=0'
        self.figure_approval_status = ''
        
        self.debug_level = 0

    #---------------------------------------------------------------
    # This function is called once for each subconfiguration
    #----------------------------------------------------------------------
    def plot_raa(self):

        for result in self.config_james_results:
            self.init_result_james(result)
            self.plot_result()

        #for result in self.config_laura_results:
        #    self.init_result_laura(result)
        #    self.plot_result()

    #---------------------------------------------------------------
    # Initialize
    #---------------------------------------------------------------
    def init_result_james(self, result):
    
        # Init from config file
        self.observable = result['observable']
        self.jetR = result['jetR']
        self.zcut = result['zcut']
        self.obs_label = result['obs_label']
        self.min_pt = result['min_pt']
        self.max_pt = result['max_pt']
        self.file_pp_name = result['file_pp']
        self.file_AA_name = result['file_AA']
        self.set_xy_titles()

        # Get hists from ROOT file
        self.file_pp = ROOT.TFile(self.file_pp_name, 'READ')
        self.file_AA = ROOT.TFile(self.file_AA_name, 'READ')
        
        self.main_result_name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, self.jetR, self.obs_label, self.min_pt, self.max_pt)
        self.sys_total_name = 'hResult_{}_systotal_R{}_{}_{}-{}'.format(self.observable, self.jetR, self.obs_label, self.min_pt, self.max_pt)
        self.h_tagging_fraction_name = 'h_tagging_fraction_R{}_{}'.format(self.jetR, self.obs_label)

        self.h_main_pp = self.file_pp.Get(self.main_result_name)
        self.h_sys_pp = self.file_pp.Get(self.sys_total_name)
        self.h_tagging_fraction_pp = self.file_pp.Get(self.h_tagging_fraction_name)
        self.h_main_AA = self.file_AA.Get(self.main_result_name)
        self.h_sys_AA = self.file_AA.Get(self.sys_total_name)
        self.h_tagging_fraction_AA = self.file_AA.Get(self.h_tagging_fraction_name)
    
        if self.debug_level > 0:
            integral_inclusive = self.h_main_pp.Integral(1, self.h_main_pp.GetNbinsX()+1, 'width')
            integral_reported = self.h_main_pp.Integral(2, self.h_main_pp.GetNbinsX()+1, 'width')
            tagging_fraction = integral_reported / integral_inclusive
            print('tagging fraction pp re-computed: {}'.format(tagging_fraction))
            
            integral_inclusive = self.h_main_AA.Integral(1, self.h_main_pp.GetNbinsX()+1, 'width')
            integral_reported = self.h_main_AA.Integral(2, self.h_main_pp.GetNbinsX()+1, 'width')
            tagging_fraction = integral_reported / integral_inclusive
            print('tagging fraction AA re-computed: {}'.format(tagging_fraction))
    
        pt = (self.min_pt + self.max_pt)/2
        bin = self.h_tagging_fraction_pp.GetXaxis().FindBin(pt)
        self.f_tagging_pp = self.h_tagging_fraction_pp.GetBinContent(bin)
        pt = (self.min_pt + self.max_pt)/2
        bin = self.h_tagging_fraction_AA.GetXaxis().FindBin(pt)
        self.f_tagging_AA = self.h_tagging_fraction_AA.GetBinContent(bin)

        self.ymax = self.h_main_pp.GetMaximum()
    
        if self.plot_theory:
            if self.plot_data:
                name = 'hRAA_{}_R{}_Theory.pdf'.format(self.observable, self.utils.remove_periods(self.jetR))
            else:
                name = 'hRAA_{}_R{}_TheoryOnly.pdf'.format(self.observable, self.utils.remove_periods(self.jetR))
        else:
            name = 'hRAA_{}_{}_R{}.pdf'.format(self.observable, self.obs_label, self.utils.remove_periods(self.jetR))
        self.output_filename = os.path.join(self.output_dir, name)
        
        # Load theory predictions
        if self.plot_theory:
            self.init_theory_james()
    
    #---------------------------------------------------------------
    # This function is called once for each subconfiguration
    #----------------------------------------------------------------------
    def plot_result(self):
        print('Plot final results for {}: R = {}, {}'.format(self.observable, self.jetR, self.obs_label))

        self.utils.set_plotting_options()
        ROOT.gROOT.ForceStyle()

        c = ROOT.TCanvas('c', 'c', 600, 650)
        c.Draw()
        
        c.cd()
        pad1 = ROOT.TPad('myPad', 'The pad',0,0.45,1,1)
        pad1.SetLeftMargin(0.2)
        pad1.SetTopMargin(0.08)
        pad1.SetRightMargin(0.04)
        pad1.SetBottomMargin(0.)
        pad1.Draw()
        pad1.cd()

        myLegend = ROOT.TLegend(0.25,0.65,0.45,0.85)
        self.utils.setup_legend(myLegend,0.055)
        
        self.h_main_pp.SetMarkerSize(1.5)
        self.h_main_pp.SetMarkerStyle(self.markers[0])
        self.h_main_pp.SetMarkerColor(self.colors[0])
        self.h_main_pp.SetLineStyle(1)
        self.h_main_pp.SetLineWidth(2)
        self.h_main_pp.SetLineColor(self.colors[0])
          
        self.h_sys_pp.SetLineColor(0)
        self.h_sys_pp.SetFillColor(self.colors[0])
        self.h_sys_pp.SetFillColorAlpha(self.colors[0], 0.3)
        self.h_sys_pp.SetFillStyle(1001)
        self.h_sys_pp.SetLineWidth(0)
        
        self.h_main_AA.SetMarkerSize(1.5)
        self.h_main_AA.SetMarkerStyle(self.markers[1])
        self.h_main_AA.SetMarkerColor(self.colors[1])
        self.h_main_AA.SetLineStyle(1)
        self.h_main_AA.SetLineWidth(2)
        self.h_main_AA.SetLineColor(self.colors[1])
          
        self.h_sys_AA.SetLineColor(0)
        self.h_sys_AA.SetFillColor(self.colors[1])
        self.h_sys_AA.SetFillColorAlpha(self.colors[1], 0.3)
        self.h_sys_AA.SetFillStyle(1001)
        self.h_sys_AA.SetLineWidth(0)
          
        xmin = self.h_sys_pp.GetBinLowEdge(2)
        xmax = self.h_sys_pp.GetBinLowEdge(self.h_sys_pp.GetNbinsX()+1)
            
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, xmin, xmax)
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(self.xtitle)
        myBlankHisto.SetYTitle(self.ytitle)
        if self.observable == 'theta_g':
            if self.jetR == 0.2:
                ymax = 1.9*self.ymax
            else:
                ymax = 2.5*self.ymax
        else:
            ymax = 3.*self.ymax
        myBlankHisto.SetMaximum(ymax)
        myBlankHisto.SetMinimum(2e-4) # Don't draw 0 on top panel
        myBlankHisto.GetYaxis().SetTitleSize(0.065)
        myBlankHisto.GetYaxis().SetTitleSize(0.08)
        myBlankHisto.GetYaxis().SetTitleOffset(1.15)
        myBlankHisto.GetYaxis().SetLabelSize(0.06)
        myBlankHisto.Draw('E')
        
        if self.observable == 'theta_g':
            rg_axis_tf1 = ROOT.TF1('rg_axis_tf1', 'x', 0, self.jetR*xmax-0.01)
            rg_axis = ROOT.TGaxis(xmin, ymax, xmax, ymax, 'rg_axis_tf1', 505, '- S')
            rg_axis.SetTitle('#it{R}_{g}')
            rg_axis.SetTitleSize(25)
            rg_axis.SetTitleFont(43)
            rg_axis.SetTitleOffset(0.62)
            rg_axis.SetLabelFont(43)
            rg_axis.SetLabelSize(25)
            rg_axis.SetTickSize(0.025)
            rg_axis.SetLabelOffset(0.016)
            rg_axis.Draw()
     
        c.cd()
        pad2 = ROOT.TPad('pad2', 'pad2', 0, 0.02, 1, 0.45)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.21)
        pad2.SetLeftMargin(0.2)
        pad2.SetRightMargin(0.04)
        pad2.SetTicks(0,1)
        pad2.Draw()
        pad2.cd()
              
        myBlankHisto2 = myBlankHisto.Clone('myBlankHisto_C')
        myBlankHisto2.SetYTitle('#frac{Pb#font[122]{-}Pb}{pp}')
        myBlankHisto2.SetXTitle(self.xtitle)
        myBlankHisto2.GetXaxis().SetTitleSize(30)
        myBlankHisto2.GetXaxis().SetTitleFont(43)
        myBlankHisto2.GetXaxis().SetTitleOffset(1.95)
        myBlankHisto2.GetXaxis().SetLabelFont(43)
        myBlankHisto2.GetXaxis().SetLabelSize(25)
        myBlankHisto2.GetYaxis().SetTitleSize(25)
        myBlankHisto2.GetYaxis().SetTitleFont(43)
        myBlankHisto2.GetYaxis().SetTitleOffset(2.)
        myBlankHisto2.GetYaxis().SetLabelFont(43)
        myBlankHisto2.GetYaxis().SetLabelSize(25)
        myBlankHisto2.GetYaxis().SetNdivisions(505)
        
        if self.observable == 'theta_g':
        
            if self.jetR == 0.2:
          
                myBlankHisto2.GetYaxis().SetRangeUser(0., 2.7)
                
                ratio_legend = ROOT.TLegend(0.25,0.72,0.43,0.975)
                self.utils.setup_legend(ratio_legend,0.052)
                ratio_legend2 = ROOT.TLegend(0.54,0.72,0.7,0.87)
                self.utils.setup_legend(ratio_legend2,0.05)
                
            else:
          
                myBlankHisto2.GetYaxis().SetRangeUser(0., 2.7)
                
                ratio_legend = ROOT.TLegend(0.22,0.8,0.35,0.97)
                self.utils.setup_legend(ratio_legend,0.05)
                ratio_legend2 = ROOT.TLegend(0.57,0.8,0.69,0.97)
                self.utils.setup_legend(ratio_legend2,0.05)
            
        else:
        
            if self.jetR == 0.2:
                myBlankHisto2.GetYaxis().SetRangeUser(0.61, 1.59)
            else:
                myBlankHisto2.GetYaxis().SetRangeUser(0.41, 1.79)

            ratio_legend = ROOT.TLegend(0.23,0.75,0.4,0.97)
            self.utils.setup_legend(ratio_legend,0.054)
            ratio_legend2 = ROOT.TLegend(0.65,0.75,0.8,0.97)
            self.utils.setup_legend(ratio_legend2,0.054)

        myBlankHisto2.Draw('')
          
        h_ratio = self.h_main_AA.Clone()
        h_ratio.Divide(self.h_main_pp)
        h_ratio.SetMarkerColor(self.ratio_color)
        h_ratio.SetLineColor(self.ratio_color)

        h_ratio_sys = self.h_sys_AA.Clone()
        h_ratio_sys.Divide(self.h_sys_pp)
        h_ratio_sys.SetFillColor(self.ratio_color)
        h_ratio_sys.SetFillColorAlpha(self.ratio_color, 0.3)

        pad1.cd()
        if self.plot_data:
            self.h_sys_pp.DrawCopy('E2 same')
            self.h_sys_AA.DrawCopy('E2 same')
            self.h_main_pp.DrawCopy('PE X0 same')
            self.h_main_AA.DrawCopy('PE X0 same')
          
        pad2.cd()

        if self.plot_theory:

            for i, g in enumerate(self.prediction_g_list):
                  
                label = self.label_list[i]
                sublabel = self.sublabel_list[i]

                color = self.theory_colors[i]
                g.SetLineColor(color)
                g.SetFillColor(color)
                if type(g) in [ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                    g.SetLineColor(0)
                    #g.Draw("3 same")
                    if self.observable == 'theta_g':
                        if self.jetR == 0.4:
                            if 'Pablos' in label or 'JETSCAPE' in label:
                                ratio_legend.AddEntry(g, '{}{}'.format(label, sublabel), 'F')
                            else:
                                ratio_legend2.AddEntry(g, '{}{}'.format(label, sublabel), 'F')
                        elif self.jetR == 0.2:
                            if 'Yuan' in label:
                                ratio_legend2.AddEntry(g, '{}{}'.format(label, sublabel), 'F')
                            else:
                                ratio_legend.AddEntry(g, '{}{}'.format(label, sublabel), 'F')
                    elif self.observable == 'zg':
                        if 'Pablos' in label:
                            ratio_legend2.AddEntry(g, '{}{}'.format(label, sublabel), 'F')
                        else:
                            ratio_legend.AddEntry(g, '{}{}'.format(label, sublabel), 'F')
                elif type(g) == ROOT.TGraph:
                    g.SetLineStyle(self.line_style[i])
                    g.SetLineWidth(4)
                    #g.Draw("same")
                    if self.observable == 'theta_g':
                        if self.jetR == 0.4:
                            if 'Pablos' in label or 'JETSCAPE' in label:
                                ratio_legend.AddEntry(g, '{}{}'.format(label, sublabel), 'L')
                            else:
                                ratio_legend2.AddEntry(g, '{}{}'.format(label, sublabel), 'L')
                        if self.jetR == 0.2:
                            if 'Yuan' in label:
                                ratio_legend2.AddEntry(g, '{}{}'.format(label, sublabel), 'L')
                            else:
                                ratio_legend.AddEntry(g, '{}{}'.format(label, sublabel), 'L')
                    elif self.observable == 'zg':
                        if 'Pablos' in label:
                            ratio_legend2.AddEntry(g, '{}{}'.format(label, sublabel), 'L')
                        else:
                            ratio_legend.AddEntry(g, '{}{}'.format(label, sublabel), 'L')

        # Draw curves in specified order
        if self.plot_theory:
            for i in self.draw_order:
                g = self.prediction_g_list[i]
                if type(g) in [ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                    g.Draw("3 same")
                elif type(g) == ROOT.TGraph:
                    g.Draw("same")
          
        ratio_legend.Draw()
        ratio_legend2.Draw()

        line = ROOT.TLine(xmin,1,xmax,1)
        line.SetLineColor(920+2)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
         
        if self.plot_data:
            h_ratio_sys.DrawCopy('E2 same')
            h_ratio.DrawCopy('PE X0 same')
       
        pad1.cd()
        myLegend.AddEntry(self.h_main_pp, 'pp', 'PE')
        myLegend.AddEntry(self.h_main_AA, 'Pb#font[122]{-}Pb 0#font[122]{-}10%', 'PE')
        myLegend.AddEntry(self.h_sys_pp, 'Sys. uncertainty', 'f')
        
        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        
        text_latex.SetTextSize(0.055)
        text = '#it{{f}}_{{tagged}}^{{pp}} = {:.2f}, #it{{f}}_{{tagged}}^{{AA}} = {:.2f}'.format(self.f_tagging_pp, self.f_tagging_AA)
        text_latex.DrawLatex(0.56, 0.36, text)

        text_latex.SetTextSize(0.065)
        text = 'ALICE {}'.format(self.figure_approval_status)
        text_latex.DrawLatex(0.56, 0.83, text)
        
        text = '#sqrt{#it{s_{#it{NN}}}} = 5.02 TeV'
        text_latex.DrawLatex(0.56, 0.75, text)

        text = 'Charged jets   anti-#it{k}_{T}'
        text_latex.DrawLatex(0.56, 0.67, text)
        
        text = '#it{R} = ' + str(self.jetR) + '   | #it{{#eta}}_{{jet}}| < {}'.format(0.9-self.jetR)
        text_latex.DrawLatex(0.56, 0.6, text)
        
        text = str(self.min_pt) + ' < #it{p}_{T, ch jet} < ' + str(self.max_pt) + ' GeV/#it{c}'
        text_latex.DrawLatex(0.56, 0.53, text)
        
        text = self.formatted_grooming_label.format(self.zcut)
        text_latex.DrawLatex(0.56, 0.44, text)
        
        myLegend.Draw()

        c.SaveAs(self.output_filename)
        c.Close()
        
    #---------------------------------------------------------------
    # Set axis title for given observable
    #----------------------------------------------------------------------
    def set_xy_titles(self):
    
        if self.observable == 'rg':
            self.xtitle = '#it{R}_{g}'
        if self.observable == 'theta_g':
            self.xtitle = '#it{#theta}_{g}'
        if self.observable == 'zg':
            self.xtitle = '#it{z}_{g}'
            
        self.ytitle = '#frac{{1}}{{#it{{#sigma}}_{{jet, inc}}}} #frac{{d#it{{#sigma}}}}{{d{}}}'.format(self.xtitle)
        
    #---------------------------------------------------------------
    # Initialize theory predictions
    #---------------------------------------------------------------
    def init_theory_james(self):
    
        # Get theory predictions
        if self.plot_theory:
            if self.observable == 'theta_g':
                config_file = 'theory_predictions_rg.yaml'
            if self.observable == 'zg':
                config_file = 'theory_predictions_zg.yaml'
        
        # Get binning, for re-binning
        bin_edges = np.array(self.h_main_pp.GetXaxis().GetXbins())[1:]
        self.bin_array = array('d', bin_edges)

        with open(config_file, 'r') as stream:
            config = yaml.safe_load(stream)
              
        self.label_list = []
        self.observable_name_list = []
        self.sublabel_list = []
        self.prediction_g_list = []
      
        self.prediction_types = [name for name in list(config.keys())]
        for type in self.prediction_types:
      
            theory = config[type]
            configs = [name for name in list(theory.keys()) if 'config' in name ]
            for prediction in configs:
                theory_prediction = theory[prediction]
                config_jetR = theory_prediction['jetR']
                if config_jetR == self.jetR:
                
                    plot_list = theory['plot_list']
                
                    if type == 'lbnl':
                    
                        x = np.array(theory_prediction['x'])
                        y_pp = np.array(theory_prediction['y_pp'])

                        n = len(x)
                        xerr = np.zeros(n)

                        if 'f_q' in theory_prediction:
                            f_q = theory_prediction['f_q']
                            f_q_min = theory_prediction['f_q_min']
                            f_q_max = theory_prediction['f_q_max']
                            y_AA_quark = np.array(theory_prediction['y_AA_quark'])
                            y_AA_gluon = np.array(theory_prediction['y_AA_gluon'])

                            y_AA = f_q*y_AA_quark + (1-f_q)*y_AA_gluon
                            y_AA_min = f_q_min*y_AA_quark + (1-f_q_min)*y_AA_gluon
                            y_AA_max = f_q_max*y_AA_quark + (1-f_q_max)*y_AA_gluon

                            ratio = np.divide(y_AA, y_pp)
                            ratio_lower = np.divide(y_AA_min, y_pp)
                            ratio_upper = np.divide(y_AA_max, y_pp)
                            g = ROOT.TGraphAsymmErrors(n, x, ratio, xerr, xerr, ratio-ratio_lower, ratio_upper-ratio)

                        else:
                            y_AA = np.array(theory_prediction['y_AA'])
                            ratio = np.divide(y_AA, y_pp)
                            g = ROOT.TGraph(n, x, ratio)

                    elif type == 'caucal':
                                
                        x = np.array(theory_prediction['x'])
                        ratio = np.array(theory_prediction['ratio'])
                        ratio_neg_unc_tot = np.array(theory_prediction['ratio_neg_unc_tot'])
                        ratio_pos_unc_tot = np.array(theory_prediction['ratio_pos_unc_tot'])

                        n = len(x)
                        xerr = np.zeros(n)
                        g = ROOT.TGraphAsymmErrors(n, x, ratio, xerr, xerr, ratio_neg_unc_tot, ratio_pos_unc_tot)
                    
                    elif type == 'guangyou':
                    
                        x = np.array(theory_prediction['x'])
                        ratio = np.array(theory_prediction['ratio'])

                        n = len(x)
                        g = ROOT.TGraph(n, x, ratio)
                    
                    elif type in ['hybrid_model', 'yang_ting']:
                    
                        x = np.array(theory_prediction['x'])
                        ratio_lower = np.array(theory_prediction['ratio_lower'])
                        ratio_upper = np.array(theory_prediction['ratio_upper'])
                        ratio = (ratio_lower + ratio_upper) / 2.

                        n = len(x)
                        xerr = np.zeros(n)
                        g = ROOT.TGraphAsymmErrors(n, x, ratio, xerr, xerr, ratio-ratio_lower, ratio_upper-ratio)
                    
                    elif type == 'jetscape':

                        plot_ratio_directly = True
                        if plot_ratio_directly:
                            x = np.array(theory_prediction['x_ratio'])
                            ratio = np.array(theory_prediction['ratio'])
                            ratio_err = np.array(theory_prediction['ratio_err'])
                        else:
                            # Get distributions
                            xbins = np.array(theory_prediction['xbins'])
                            y_pp = np.array(theory_prediction['y_pp'])
                            y_pp_err = np.array(theory_prediction['y_pp_err'])
                            y_AA = np.array(theory_prediction['y_AA'])
                            y_AA_err = np.array(theory_prediction['y_AA_err'])
                            
                            # Rebin distributions
                            x, y_pp, y_pp_err = self.rebin_arrays(xbins, y_pp, y_pp_err)
                            _, y_AA, y_AA_err = self.rebin_arrays(xbins, y_AA, y_AA_err)

                            # Form ratio and propagate uncertainty
                            y_pp_err_fraction = np.divide(y_pp_err, y_pp)
                            y_AA_err_fraction = np.divide(y_AA_err, y_AA)
                            
                            ratio = np.divide(y_AA, y_pp)
                            ratio_err_fraction = np.sqrt(np.square(y_AA_err_fraction) + np.square(y_pp_err_fraction))
                            ratio_err = np.multiply(ratio, ratio_err_fraction)
              
                        n = len(x)
                        xerr = np.zeros(n)
                        g = ROOT.TGraphErrors(n, x, ratio, xerr, ratio_err)
    
                    if prediction in plot_list:
                        self.prediction_g_list.append(g)
                        self.label_list.append(theory['label'])
                        self.sublabel_list.append(theory_prediction['sublabel'])
                        self.observable_name_list.append(theory['observable'])

        # Set draw order (slightly hacky-- just move selected curves to start or end of draw list)
        self.draw_order = list(range(0, len(self.prediction_g_list)))
        for i, g in enumerate(self.prediction_g_list):

            label = self.label_list[i]
            sublabel = self.sublabel_list[i]
            if 'Pablos' in label:
                index = self.draw_order.index(i)
                self.draw_order.insert(0, self.draw_order.pop(index))

            if 'JETSCAPE' in label or 'quark' in sublabel or 'med' in sublabel:
                index = self.draw_order.index(i)
                self.draw_order.append(self.draw_order.pop(index))
          
        if self.observable == 'zg':
            for i, g in enumerate(self.prediction_g_list):

                label = self.label_list[i]
                if 'JETSCAPE' in label:
                    index = self.draw_order.index(i)
                    self.draw_order.insert(0, self.draw_order.pop(index))

    #---------------------------------------------------------------
    # Rebin numpy arrays (xbins,y) representing a histogram,
    # where x represents the bin edges, and y the bin content
    #
    # Return (x_rebinned, y_rebinned) where x_rebinned is the bin centers
    #
    # I don't find any easy way to do this...so I
    #        construct TH1, rebin, and then extract array)
    #----------------------------------------------------------------------
    def rebin_arrays(self, x, y, yerr):

        x_array = array('d', x)
        h = ROOT.TH1F('h', 'h', len(x)-1, x_array)
        for i, content in enumerate(y):
            h.SetBinContent(i+1, y[i])
            h.SetBinError(i+1, yerr[i])
        h_rebinned = h.Rebin(len(self.bin_array)-1, 'h_rebinned', self.bin_array)
        x_rebinned = []
        y_rebinned = []
        y_err_rebinned = []
        for i in range(len(self.bin_array)-1):
            x_rebinned.append(h_rebinned.GetBinCenter(i+1))
            y_rebinned.append(h_rebinned.GetBinContent(i+1))
            y_err_rebinned.append(h_rebinned.GetBinError(i+1))
          
        return (np.array(x_rebinned), np.array(y_rebinned), np.array(y_err_rebinned))

#----------------------------------------------------------------------
if __name__ == '__main__':

  analysis = PlotRAA()
  analysis.plot_raa()
