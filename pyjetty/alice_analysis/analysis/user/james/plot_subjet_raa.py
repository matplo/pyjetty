#! /usr/bin/env python

import sys
import os
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

        self.output_dir = './subjet_z'
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        # Load config file
        config_file = './subjet_z/plot.yaml'
        with open(config_file, 'r') as stream:
            config = yaml.safe_load(stream)
        self.config_results = [config[name] for name in config if 'result' in name ]

        self.plot_data = True
        self.plot_theory = True

        self.colors = [600-6, 632-4]
        self.ratio_color = ROOT.kGray+3
        self.markers = [20, 21]
        ROOT.gStyle.SetLineStyleString(11,'30 12')
        
        self.theory_colors = [ROOT.kViolet-8, ROOT.kAzure-4, ROOT.kTeal-8, ROOT.kOrange+6, ROOT.kRed-7, ROOT.kPink+1, ROOT.kCyan-2, ROOT.kBlue-10]
        
        self.line_style = [1, 1, 1, 1, 1, 1, 11, 1, 1]
        self.line_width = [4, 4, 4, 4, 4, 4, 6, 4, 4]
         
        self.figure_approval_status = 'Preliminary'
        
        self.xtitle = '#it{z}_{r}'
        self.ytitle = '#frac{{1}}{{#it{{#sigma}}_{{#it{{z}}_{{r}} > 0.7}}}} #frac{{d#it{{#sigma}}}}{{d{}}}'.format(self.xtitle)

        self.debug_level = 0
        
    #---------------------------------------------------------------
    # This function is called once for each subconfiguration
    #----------------------------------------------------------------------
    def plot_raa(self):

        for result in self.config_results:
            self.centrality = [0, 10]
            self.init_result(result)
            self.plot_result()

    #---------------------------------------------------------------
    # Initialize
    #---------------------------------------------------------------
    def init_result(self, result):
    
        # Init from config file
        self.observable = result['observable']
        self.jetR = result['jetR']
        self.obs_label = result['obs_label']
        self.min_pt = result['min_pt']
        self.max_pt = result['max_pt']
        self.file_pp_name = result['file_pp']
        self.file_AA_name = result['file_AA']
        self.theory_colors = self.theory_colors
        self.line_style = self.line_style
        self.line_width = self.line_width
        
        # Get binning to plot
        self.bins = np.array(result['bins'])
        self.n_bins = len(self.bins) - 1

        # Get hists from ROOT file
        self.file_pp = ROOT.TFile(self.file_pp_name, 'READ')
        self.file_AA = ROOT.TFile(self.file_AA_name, 'READ')
        
        self.main_result_name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, self.jetR, self.obs_label, self.min_pt, self.max_pt)
        self.sys_total_name = 'hResult_{}_systotal_R{}_{}_{}-{}'.format(self.observable, self.jetR, self.obs_label, self.min_pt, self.max_pt)

        h_main_AA = self.file_AA.Get(self.main_result_name)
        h_sys_AA = self.file_AA.Get(self.sys_total_name)
        #print(bins)
        #print([h_main_AA.GetXaxis().GetXbins().GetArray()[i] for i in range(n_bins+1)])
        self.h_main_AA = h_main_AA.Rebin(self.n_bins, f'{h_main_AA.GetName()}_rebinned', self.bins)
        self.h_sys_AA = h_sys_AA.Rebin(self.n_bins, f'{h_sys_AA.GetName()}_rebinned', self.bins)

        h_main_pp = self.file_pp.Get(self.main_result_name)
        h_sys_pp = self.file_pp.Get(self.sys_total_name)
        self.h_main_pp = h_main_pp.Rebin(self.n_bins, f'{h_main_pp.GetName()}_rebinned', self.bins)
        self.h_sys_pp = h_sys_pp.Rebin(self.n_bins, f'{h_sys_pp.GetName()}_rebinned', self.bins)
        
        # Normalize to the integral over the reported range
        # Note that histograms are already scaled for bin width in run_analysis.get_obs_distribution()
        integral_AA = self.h_main_AA.Integral(1, self.h_main_AA.GetNbinsX(), 'width')
        self.h_main_AA.Scale(1./integral_AA)
        self.h_sys_AA.Scale(1./integral_AA)
        print(f'integral AA: {integral_AA}')
    
        integral_pp = self.h_main_pp.Integral(1, self.h_main_pp.GetNbinsX(), 'width')
        self.h_main_pp.Scale(1./integral_pp)
        self.h_sys_pp.Scale(1./integral_pp)
        print(f'integral pp: {integral_pp}')
        
        self.h_ratio = self.h_main_AA.Clone()
        self.h_ratio.Divide(self.h_main_pp)
        self.h_ratio_sys = self.h_sys_AA.Clone()
        self.h_ratio_sys.Divide(self.h_sys_pp)
    
        self.ymax = self.h_main_pp.GetMaximum()
        
        # Load theory predictions
        if self.plot_theory:
            self.init_theory(result)
        
    #---------------------------------------------------------------
    # This function is called once for each subconfiguration
    #
    # Required histograms:
    #   - self.h_main_pp
    #   - self.h_sys_pp
    #   - self.h_main_AA
    #   - self.h_sys_AA
    #   - self.h_ratio
    #   - self.h_ratio_sys
    #
    # And the theory tgraphs should be placed in:
    #   - self.label_list
    #   - self.prediction_g_list
    #
    # As well as the tagging rates:
    #   - self.f_tagging_pp (number)
    #   - self.f_tagging_AA (number)
    #----------------------------------------------------------------------
    def plot_result(self):
        print('Plot final results for {}: R = {}, {}'.format(self.observable, self.jetR, self.obs_label))
        
        observable_label = self.observable

        if self.plot_theory:
            if self.plot_data:
                name = 'h_{}_{}-{}_R{}_r{}Theory.pdf'.format(observable_label, self.centrality[0], self.centrality[1], self.utils.remove_periods(self.jetR), self.utils.remove_periods(self.obs_label))
            else:
                name = 'h_{}_{}-{}_R{}_r{}TheoryOnly.pdf'.format(observable_label, self.centrality[0], self.centrality[1], self.utils.remove_periods(self.jetR), self.utils.remove_periods(self.obs_label))
        else:
            name = 'h_{}_{}-{}_R{}_r{}.pdf'.format(observable_label, self.centrality[0], self.centrality[1], self.utils.remove_periods(self.jetR), self.utils.remove_periods(self.obs_label))
        self.output_filename = os.path.join(self.output_dir, name)

        self.utils.set_plotting_options()
        ROOT.gROOT.ForceStyle()

        c = ROOT.TCanvas('c', 'c', 600, 650)
        c.Draw()
        
        c.cd()
        pad2_dy = 0.45
        pad1 = ROOT.TPad('myPad', 'The pad',0,pad2_dy,1,1)
        pad1.SetLeftMargin(0.2)
        pad1.SetTopMargin(0.08)
        pad1.SetRightMargin(0.04)
        pad1.SetBottomMargin(0.)
        pad1.Draw()
        pad1.cd()

        myLegend = ROOT.TLegend(0.65,0.65,0.8,0.85)
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
          
        #xmin = self.h_main_pp.GetBinLowEdge(1)
        #xmax = self.h_main_pp.GetBinLowEdge(self.h_main_pp.GetNbinsX()+1)
            
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, self.bins[0], self.bins[-1])
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(self.xtitle)
        myBlankHisto.SetYTitle(self.ytitle)
        ymax = 2.*self.ymax
        myBlankHisto.SetMaximum(ymax)
        myBlankHisto.SetMinimum(2e-4) # Don't draw 0 on top panel
        myBlankHisto.GetYaxis().SetTitleSize(0.065)
        myBlankHisto.GetYaxis().SetTitleSize(0.08)
        myBlankHisto.GetYaxis().SetTitleOffset(1.15)
        myBlankHisto.GetYaxis().SetLabelSize(0.06)
        myBlankHisto.Draw('E')

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
        
        myBlankHisto2.GetYaxis().SetRangeUser(0., 1.99)

        ratio_legend = ROOT.TLegend(0.43,0.75,0.6,0.97)
        self.utils.setup_legend(ratio_legend,0.07)

        myBlankHisto2.Draw('')
          
        self.h_ratio.SetMarkerSize(1.5)
        self.h_ratio.SetMarkerStyle(self.markers[1])
        self.h_ratio.SetMarkerColor(self.ratio_color)
        self.h_ratio.SetLineColor(self.ratio_color)

        self.h_ratio_sys.SetFillColor(self.ratio_color)
        self.h_ratio_sys.SetFillColorAlpha(self.ratio_color, 0.3)

        pad1.cd()
        if self.plot_data:
            self.h_sys_pp.Draw('E2 same')
            self.h_sys_AA.Draw('E2 same')
            self.h_main_pp.DrawCopy('PE X0 same')
            self.h_main_AA.DrawCopy('PE X0 same')
          
        pad2.cd()

        if self.plot_theory:

            for i, prediction in enumerate(self.prediction_list):
                  
                label = self.label_list[i]
                sublabel = self.sublabel_list[i]
                color = self.theory_colors[i]
                
                print(type(prediction))
                if type(prediction) in [ROOT.TH1F]:

                    prediction.SetLineColor(0)
                    prediction.SetMarkerSize(0)
                    prediction.SetMarkerStyle(0)
                    prediction.SetFillColor(color)
                    
                    prediction.Draw('E3 same')
                    ratio_legend.AddEntry(prediction, f'{label}{sublabel}', 'F')

                if type(prediction) in [ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                
                    prediction.SetFillColor(color)
                    prediction.SetLineColor(0)
                    
                    prediction.Draw("3 same")
                    ratio_legend.AddEntry(prediction, f'{label}{sublabel}', 'F')
          
        ratio_legend.Draw()
        line = ROOT.TLine(self.bins[0], 1, self.bins[-1], 1)
        line.SetLineColor(920+2)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
         
        if self.plot_data:
            self.h_ratio_sys.Draw('E2 same')
            self.h_ratio.Draw('PE X0 same')
       
        pad1.cd()
        myLegend.AddEntry(self.h_main_pp, 'pp', 'PE')
        myLegend.AddEntry(self.h_main_AA, 'Pb#font[122]{{-}}Pb {}#font[122]{{-}}{}%'.format(self.centrality[0], self.centrality[1]), 'PE')
        myLegend.AddEntry(self.h_sys_pp, 'Sys. uncertainty', 'f')
        
        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        
        x = 0.25
        text_latex.SetTextSize(0.065)
        text = 'ALICE {}'.format(self.figure_approval_status)
        text_latex.DrawLatex(x, 0.83, text)
        
        text = '#sqrt{#it{s_{#it{NN}}}} = 5.02 TeV'
        text_latex.DrawLatex(x, 0.75, text)

        text = 'Charged jets   anti-#it{k}_{T}'
        text_latex.DrawLatex(x, 0.67, text)
        
        text = '#it{R} = ' + str(self.jetR) + '   #it{{#eta}}_{{jet}}| < {}'.format(0.9-self.jetR)
        text_latex.DrawLatex(x, 0.6, text)
        
        text = str(self.min_pt) + ' < #it{p}_{T, ch jet} < ' + str(self.max_pt) + ' GeV/#it{c}'
        text_latex.DrawLatex(x, 0.51, text)
        
        text = 'Leading subjets' + '   #it{r} = ' + str(self.obs_label)
        text_latex.DrawLatex(x, 0.42, text)
        
        myLegend.Draw()

        c.SaveAs(self.output_filename)
        c.Close()

    #---------------------------------------------------------------
    # Initialize theory predictions
    #---------------------------------------------------------------
    def init_theory(self, result):
    
        self.label_list = []
        self.sublabel_list = []
        self.prediction_list = []
    
        #----------------------------------------------------------
        # JETSCAPE
        f_pp = ROOT.TFile(result['jetscape_pp'], 'READ')
        f_AA = ROOT.TFile(result['jetscape_AA'], 'READ')
        
        #h_name = f'h_chjet_subjetz_alice_R{self.jetR}_r{self.obs_label}_pt0.0Scaled'
        h_name = f'h_chjet_subjetz_alice_R{self.jetR}_r{self.obs_label}Scaled'
        h_pp = f_pp.Get(h_name)
        h_pp.SetDirectory(0)
        f_pp.Close()
        
        h_AA = f_AA.Get(h_name)
        h_AA.SetDirectory(0)
        f_AA.Close()

        # Rebin to data binning
        bin_edges = np.array(self.h_main_pp.GetXaxis().GetXbins())[0:]
        h_pp = h_pp.Rebin(bin_edges.size-1, f'{h_pp.GetName()}_rebinned', bin_edges)
        h_AA = h_AA.Rebin(bin_edges.size-1, f'{h_AA.GetName()}_rebinned', bin_edges)
        
        # Normalization
        h_pp.Scale(1., 'width')
        h_AA.Scale(1., 'width')
        h_pp.Scale(1./h_pp.Integral(1, h_pp.GetNbinsX()))
        h_AA.Scale(1./h_AA.Integral(1, h_AA.GetNbinsX()))

        # Form ratio
        hRatio = h_AA.Clone()
        hRatio.SetName(f'hRatio_{h_AA.GetName()}')
        hRatio.Divide(h_pp)
      
        # Add to theory list
        self.prediction_list.append(hRatio)
        self.label_list.append('JETSCAPE')
        self.sublabel_list.append('')
      
        #----------------------------------------------------------
        

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
