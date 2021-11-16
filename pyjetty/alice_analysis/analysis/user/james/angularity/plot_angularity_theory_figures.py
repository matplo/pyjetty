"""
  macro for plotting multi-paneled angularity figures
  """

# General
import os
import sys
import math
import yaml
import argparse

# Data analysis and plotting
import ROOT
import numpy as np
from array import *

# Base class
from pyjetty.alice_analysis.analysis.base import common_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class PlotAngularityFigures(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, output_dir='', **kwargs):
        super(PlotAngularityFigures, self).__init__(**kwargs)

        self.output_dir = output_dir
        self.file_format = '.pdf'

        #------------------------------------------------------

        #self.base_dir = '/home/james/plot-angularity/'
        self.base_dir = '/Users/jamesmulligan/Analysis_Angularity/plot-angularity'
        self.file = 'ang/final_results/fFinalResults.root'
        self.beta_list = [1.5, 2, 3]
        self.R_list = [0.2, 0.4]
        self.min_pt = 60
        self.max_pt = 80
        self.folding_labels = ['PYTHIA8', 'P8time', 'Herwig7']

        self.xmin = -0.01
        self.ymax = 18.99
        self.ymin_ratio = 0.1
        self.ymax_ratio = 15
        self.scale_factor_groomed_beta3 = 0.04
        self.scale_factor_groomed_beta2 = 0.2
        self.scale_factor_groomed_beta15 = 0.5
        self.scale_factor_ungroomed_beta2 = 0.5
        self.scale_factor_ungroomed_beta3 = 0.1

        self.xtitle =  '#it{#lambda}_{#it{#beta}}'
        self.ytitle = '#frac{d#it{#sigma}}{d#it{#lambda}_{#it{#beta}}} #/#int_{#it{#lambda}_{#it{#beta}}^{NP}}^{1} #frac{d#it{#sigma}}{d#it{#lambda}_{#it{#beta}}} d#it{#lambda}_{#it{#beta}}'
        self.xtitle_groomed =  '#it{#lambda}_{#it{#beta},g}'
        self.ytitle_groomed = '#frac{1}{#it{#sigma}_{jet}} #frac{d#it{#sigma}}{d#it{#lambda}_{#it{#beta},g}}'

        self.left_offset = 0.4
        self.bottom_offset = 0.15
        self.ratio_height = 0.25
        self.top_bottom_scale_factor = (1 - self.ratio_height) / (1 - self.bottom_offset - self.ratio_height)

        #------------------------------------------------------

        self.marker_data = 21
        self.markers = [20, 34, 33, 22, 23]
        self.marker_size = 3
        self.marker_size_ratio = 2
        self.alpha = 0.7
        self.color_data = 1
        self.colors = [ROOT.kRed-7, ROOT.kTeal-8, ROOT.kViolet-8, ROOT.kBlue-9]
        
        #------------------------------------------------------

        # P vs NP cutoff point: lambda_beta ~ Lambda / (pT * R) -- use avg value of pT for the bin.
        # Formula assumes that jet pT xsec falls like pT^(-5.5)
        formula_pt = (4.5/3.5)*(self.min_pt**-3.5 - self.max_pt**-3.5)/(self.min_pt**-4.5 - self.max_pt**-4.5)
        self.Lambda = 1.
        self.lambda_np = {}
        for R in self.R_list:
            self.lambda_np[R] = round(self.Lambda / (formula_pt * R), 2)

        #------------------------------------------------------
        # Store paths to all final results in a dictionary
        self.predictions = {}
        
        subdir_path = os.path.join(self.base_dir, 'AngR02_ptbin3')
        filename = os.path.join(subdir_path, 'ang/final_results/fFinalResults.root')
        self.predictions['0.2'] = filename

        subdir_path = os.path.join(self.base_dir, 'AngR04_ptbin3')
        filename = os.path.join(subdir_path, 'ang/final_results/fFinalResults.root')
        self.predictions['0.4'] = filename

        print(self)

    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    def plot_results(self):

        self.setOptions()
        ROOT.gROOT.ForceStyle()

        self.plot_multipanel(R=0.2, groomed=False)
        self.plot_multipanel(R=0.4, groomed=False)

    #-------------------------------------------------------------------------------------------
    def plot_multipanel(self, R=1, groomed=False):

        # Create multi-panel canvas
        cname = 'c{}{}'.format(R,groomed)
        c = ROOT.TCanvas(cname,cname,2400,1400)
        c.SetRightMargin(0.05);
        c.SetLeftMargin(self.left_offset);
        c.SetTopMargin(0.03);
        c.SetBottomMargin(self.bottom_offset/2);
        c.cd()
        c.Divide(3, 1, 0.01, 0.)

        # Keep histograms in memory, otherwise there can be problems with double deletes (i.e. ROOT then python deletes)
        self.plot_list = []
        self.g_theory_dict = {'1.5': [], '2': [], '3': []}
        self.h_ratio_dict = {'1.5': [], '2': [], '3': []}
        self.g_ratio_dict = {'1.5': [], '2': [], '3': []}

        # Plot each pt bin in its own pad
        self.plot_angularity(c, pad=1, R=R, beta = '1.5', groomed=groomed)
        self.plot_angularity(c, pad=2, R=R, beta = '2',   groomed=groomed)
        self.plot_angularity(c, pad=3, R=R, beta = '3',   groomed=groomed)

        if groomed:
            outfilename = 'hJetAngularity_Theory_{}_SD{}'.format(R, self.file_format)
        else:
            outfilename = 'hJetAngularity_Theory_{}{}'.format(R, self.file_format)
        output_filename = os.path.join(self.output_dir, outfilename)
        c.SaveAs(output_filename)

    #-------------------------------------------------------------------------------------------
    # Get beta histograms from file, and call plot_beta_overlay to draw them
    #-------------------------------------------------------------------------------------------
    def plot_angularity(self, c, pad=0, R=1, beta = '', groomed=False):
        
        filename = self.predictions[str(R)]
        f = ROOT.TFile(filename, 'READ')

        self.h = None
        self.h_sys = None
        self.blank_histo_list = []

        grooming_label = ''
        if groomed:
            grooming_label = '_SD_zcut02_B0'
            
        # Get data hist
        h_name ='hmain_ang_R{}_{}{}_{}-{}_trunc'.format(R, beta, grooming_label, self.min_pt, self.max_pt)
        h_sys_name =  'hResult_ang_systotal_R{}_{}{}_n3_{}-{}'.format(R, beta, grooming_label, self.min_pt, self.max_pt)
        self.h = f.Get(h_name)
        self.h_sys = f.Get(h_sys_name)
        self.h.SetDirectory(0)
        self.h_sys.SetDirectory(0)
        
        # Normalize such that integral in perturbative region is 1
        n = self.h.GetNbinsX()
        min_bin = self.h.FindBin(self.lambda_np[R] - 0.001) + 1
        integral_perturbative = self.h.Integral(min_bin, n+1, 'width')
        self.h.Scale(1./integral_perturbative)
        self.h_sys.Scale(1./integral_perturbative)
        
        self.scale_label = self.scale_histogram_for_visualization(self.h_sys, beta, groomed)
        self.scale_histogram_for_visualization(self.h, beta, groomed)
        
        # Get folded theory predictions
        for i,folding_label in enumerate(self.folding_labels):

            beta_label = beta
            if beta == '1.5':
                beta_label = 15
            name_cent = 'theory_cent_ang_R{}_{}_ch_PtBin60-80_{}'.format(self.remove_periods(str(R)), beta_label, i)
            name_min = 'theory_min_ang_R{}_{}_ch_PtBin60-80_{}'.format(self.remove_periods(str(R)), beta_label, i)
            name_max = 'theory_max_ang_R{}_{}_ch_PtBin60-80_{}'.format(self.remove_periods(str(R)), beta_label, i)
            
            h_theory_cent = f.Get(name_cent)
            h_theory_min = f.Get(name_min)
            h_theory_max = f.Get(name_max)
            h_theory_cent.SetDirectory(0)
            h_theory_min.SetDirectory(0)
            h_theory_max.SetDirectory(0)
            
            # Rebin to data binning
            h_theory_cent_rebinned = self.rebin_theory(h_theory_cent, self.h)
            h_theory_min_rebinned = self.rebin_theory(h_theory_min, self.h)
            h_theory_max_rebinned = self.rebin_theory(h_theory_max, self.h)
            
            # Normalize such that integral in perturbative region is 1
            min_bin =  h_theory_cent_rebinned.FindBin(self.lambda_np[R]-0.001) + 1
            integral_perturbative_theory =  h_theory_cent_rebinned.Integral(min_bin, n+1, 'width')
            integral_total =  h_theory_cent_rebinned.Integral(1, n+1, 'width')

            h_theory_cent_rebinned.Scale(1./integral_perturbative_theory)
            h_theory_min_rebinned.Scale(1./integral_perturbative_theory)
            h_theory_max_rebinned.Scale(1./integral_perturbative_theory)
                        
            # Scale additionally for visualization
            self.scale_histogram_for_visualization(h_theory_cent_rebinned, beta, groomed)
            self.scale_histogram_for_visualization(h_theory_min_rebinned, beta, groomed)
            self.scale_histogram_for_visualization(h_theory_max_rebinned, beta, groomed)
            
            x = np.array([h_theory_cent_rebinned.GetXaxis().GetBinCenter(i) for i in range(1, n+1)])
            y = np.array([h_theory_cent_rebinned.GetBinContent(i) for i in range(1, n+1)])
            xerrup = xerrdn = np.array([0. for i in range(n)])
            yerrup = np.array([h_theory_max_rebinned.GetBinContent(i)-y[i-1] for i in range(1, n+1)])
            yerrdn = np.array([y[i-1]-h_theory_min_rebinned.GetBinContent(i) for i in range(1, n+1)])
            g_theory = ROOT.TGraphAsymmErrors(n, x, y, xerrdn, xerrup, yerrdn, yerrup)
            g_theory.SetName('g_theory_{}_{}'.format(i, beta))
            self.g_theory_dict[beta].append(g_theory)

            # Construct ratios in self.h_ratio_dict, self.g_ratio_dict
            self.construct_ratio(self.h, self.h_sys, h_theory_cent_rebinned, h_theory_min_rebinned, h_theory_max_rebinned, n, x, xerrup, xerrdn, y, yerrup, yerrdn, R, beta, i, pad)
            
        f.Close()

        # Plot overlay of beta values
        self.plot_beta_overlay(c, pad, R, beta, groomed)

        # Keep histograms in memory
        self.plot_list.append(self.h)
        self.plot_list.append(self.h_sys)
        self.plot_list.append(self.g_theory_dict)
        self.plot_list.append(self.h_ratio_dict)
        self.plot_list.append(self.g_ratio_dict)
        self.plot_list.append(self.blank_histo_list)

    #-------------------------------------------------------------------------------------------
    # Draw beta histograms in given pad
    #-------------------------------------------------------------------------------------------
    def plot_beta_overlay(self, c, pad, R, beta, groomed):

        # Create canvas
        c.cd(pad)

        # Set pad to plot distributions
        pad1 = ROOT.TPad("pad1_{}".format(R), "pad1{}".format(R),0, self.bottom_offset + self.ratio_height,1,1)
        self.plot_list.append(pad1)
        if pad in [1]:
            pad1.SetLeftMargin(self.left_offset)
        else:
            pad1.SetLeftMargin(0.)
        pad1.SetRightMargin(0.)
        pad1.SetTopMargin(0.0)
        pad1.SetBottomMargin(0.)
        pad1.SetTicks(0,1)
        pad1.Draw()
        pad1.cd()

        # Draw blank histos
        blankname = 'myBlankHisto_{}_{}'.format(pad, R)
        xmax = self.h.GetXaxis().GetBinUpEdge(self.h.GetXaxis().GetNbins())
        myBlankHisto = ROOT.TH1F(blankname,blankname, 1, self.xmin, xmax)
        myBlankHisto.SetNdivisions(505)
        if groomed:
            myBlankHisto.SetYTitle(self.ytitle_groomed)
        else:
            myBlankHisto.SetYTitle(self.ytitle)    
        myBlankHisto.GetYaxis().SetTitleSize(0.1)
        myBlankHisto.GetYaxis().SetTitleOffset(1.4)
        myBlankHisto.GetYaxis().SetLabelSize(0.06)
        myBlankHisto.SetMinimum(0.001)
        myBlankHisto.SetMaximum(self.ymax)
        myBlankHisto.Draw()
        self.blank_histo_list.append(myBlankHisto)

        scale_factor = 1.
        if pad in [1]:
            shift = 0.0
            shift2 = 0.0
        else:
            shift = -0.08
            shift2 = -0.08 - shift
        
        leg = ROOT.TLegend(0.2+shift,0.96-0.072*scale_factor*4,0.6,0.96)
        self.setupLegend(leg,0.07)
        self.plot_list.append(leg)

        leg3 = ROOT.TLegend(0.2+shift,0.95-0.072*scale_factor*2,0.55,0.95)
        if pad == 3:
            self.setupLegend(leg3,0.07)
            self.plot_list.append(leg3)
            
        line_lambda_np = ROOT.TLine(self.lambda_np[R], 0, self.lambda_np[R], self.ymax*0.6)
        line_lambda_np.SetLineColor(self.colors[-1])
        line_lambda_np.SetLineStyle(2)
        line_lambda_np.SetLineWidth(2)
        self.plot_list.append(line_lambda_np)

        # Draw data
        self.h.SetMarkerColor(self.color_data)
        self.h.SetLineColor(self.color_data)
        self.h.SetLineWidth(2)
        self.h.SetMarkerStyle(self.marker_data)
        self.h.SetMarkerSize(self.marker_size)

        self.h_sys.SetLineColor(0)
        self.h_sys.SetMarkerSize(0)
        self.h_sys.SetMarkerColor(0)
        self.h_sys.SetFillColor(self.color_data)
        self.h_sys.SetFillColorAlpha(self.color_data, 0.3)
        self.h_sys.SetFillStyle(1001)
        self.h_sys.SetLineWidth(0)
        
        leg.AddEntry(self.h,'ALICE','PE')
        
        # Draw theory
        for i,folding_label in enumerate(self.folding_labels):
        
            self.g_theory_dict[beta][i].SetFillColorAlpha(self.colors[i], 0.25)
            self.g_theory_dict[beta][i].SetLineColor(self.colors[i])
            self.g_theory_dict[beta][i].SetLineWidth(3)
            self.g_theory_dict[beta][i].Draw('L 3 same')
            
            leg.AddEntry(self.g_theory_dict[beta][i],'NLO+NLL #otimes {}'.format(folding_label),'lf')
            
        self.h_sys.Draw('E2 same')
        line_lambda_np.Draw()
        self.h.Draw('PE X0 same')

        if pad == 2:
            leg.Draw('same')
        if pad == 3:
            leg3.Draw('same')
        
        if pad == 3:
            leg3.AddEntry(line_lambda_np, '#it{#lambda}_{#it{#beta}}^{NP} #leq #Lambda / (#it{p}_{T,jet}^{ch} #it{R})', 'lf')
            leg3.AddEntry(self.h_sys, 'Sys. uncertainty', 'f')

        # # # # # # # # # # # # # # # # # # # # # # # #
        # text
        # # # # # # # # # # # # # # # # # # # # # # # #
        ymax = 0.93
        dy = 0.07
        x = 0.45 + shift + shift2
        size = 0.06
    
        if pad == 1:
            system0 = ROOT.TLatex(x,ymax,'ALICE')
            system0.SetNDC()
            system0.SetTextSize(size*scale_factor)
            system0.Draw()

            system1 = ROOT.TLatex(x,ymax-dy,'pp  #sqrt{#it{s}} = 5.02 TeV')
            system1.SetNDC()
            system1.SetTextSize(size*scale_factor)
            system1.Draw()

            system2 = ROOT.TLatex(x,ymax-2*dy,'charged jets   anti-#it{k}_{T}')
            system2.SetNDC()
            system2.SetTextSize(size*scale_factor)
            system2.Draw()

            system3 = ROOT.TLatex(x,ymax-3*dy, '#it{{R}} = {}    | #it{{#eta}}_{{jet}}| < {}'.format(R, 0.9-R))
            system3.SetNDC()
            system3.SetTextSize(size*scale_factor)
            system3.Draw()
            
            system4 = ROOT.TLatex(x,ymax-4.*dy-0.02, str(self.min_pt) + ' < #it{p}_{T,jet}^{ch} < ' + str(self.max_pt) + ' GeV/#it{c}')
            system4.SetNDC()
            system4.SetTextSize(size*scale_factor)
            system4.Draw()

            self.plot_list.append(system0)
            self.plot_list.append(system1)
            self.plot_list.append(system2)
            self.plot_list.append(system3)
            self.plot_list.append(system4)
            
        system5 = ROOT.TLatex(x+0.2,ymax-6.*dy, '#it{{#beta}} = {}{}'.format(beta,self.scale_label))
        system5.SetNDC()
        if pad in [1]:
            beta_size = size
        else:
            beta_size = 1.3*size
        system5.SetTextSize(beta_size)
        system5.Draw()
        self.plot_list.append(system5)

        if groomed and pad == 1:
            system5 = ROOT.TLatex(x, ymax-4.2*dy, 'Soft Drop #it{z}_{cut} = 0.2 #beta_{SD} = 0')
            system5.SetNDC()
            system5.SetTextSize(size*scale_factor)
            system5.Draw()
            self.plot_list.append(system5)

        # Set pad for ratio
        c.cd(pad)
        pad2 = ROOT.TPad("pad2_{}".format(R), "pad2{}".format(R),0,0,1,self.bottom_offset+self.ratio_height)
        self.plot_list.append(pad2)
        if pad in [1]:
            pad2.SetLeftMargin(self.left_offset)
        else:
            pad2.SetLeftMargin(0.)
        pad2.SetRightMargin(0.)
        pad2.SetTopMargin(0.)
        pad2.SetBottomMargin(self.bottom_offset/(self.bottom_offset+self.ratio_height))
        pad2.SetTicks(1,1)
        pad2.SetLogy()
        pad2.Draw()
        pad2.cd()

        # Draw blank histos
        blankname = 'myBlankHisto2_{}_{}'.format(pad, R)
        myBlankHisto2 = ROOT.TH1F(blankname,blankname, 1, self.xmin, xmax-0.001)
        myBlankHisto2.SetNdivisions(505, "y")
        myBlankHisto2.SetNdivisions(505, "x")
        if groomed:
            myBlankHisto2.SetXTitle(self.xtitle_groomed)
        else:
            myBlankHisto2.SetXTitle(self.xtitle)
        myBlankHisto2.SetYTitle('#frac{Data}{Theory}')
        myBlankHisto2.GetXaxis().SetTitleSize(0.15)
        myBlankHisto2.GetXaxis().SetTitleOffset(1.)
        myBlankHisto2.GetXaxis().SetLabelSize(0.1)
        myBlankHisto2.GetYaxis().SetTitleSize(0.12)
        myBlankHisto2.GetYaxis().SetTitleOffset(1.)
        myBlankHisto2.GetYaxis().SetLabelSize(0.1)
        myBlankHisto2.SetMinimum(self.ymin_ratio)
        myBlankHisto2.SetMaximum(self.ymax_ratio)
        myBlankHisto2.Draw()
        self.blank_histo_list.append(myBlankHisto2)

        line = ROOT.TLine(self.xmin,1,xmax,1)
        line.SetLineColor(1)
        line.SetLineStyle(2)
        line.Draw('same')
        self.plot_list.append(line)

        # Draw ratio
        for i,folding_label in enumerate(self.folding_labels):
        
            # Draw tgraph with sys uncertainty
            self.g_ratio_dict[beta][i].SetFillColorAlpha(self.colors[i], 0.25)
            self.g_ratio_dict[beta][i].SetLineColor(self.colors[i])
            self.g_ratio_dict[beta][i].SetLineWidth(3)
            self.g_ratio_dict[beta][i].Draw('3 same')
        
            # Draw th1 with stat uncertainty
            self.h_ratio_dict[beta][i].SetMarkerColorAlpha(self.colors[i], self.alpha)
            self.h_ratio_dict[beta][i].SetLineColorAlpha(self.colors[i], self.alpha)
            self.h_ratio_dict[beta][i].SetFillColorAlpha(self.colors[i], self.alpha)
            self.h_ratio_dict[beta][i].SetLineColor(self.colors[i])
            self.h_ratio_dict[beta][i].SetLineWidth(2)
            self.h_ratio_dict[beta][i].SetMarkerStyle(self.markers[i])
            self.h_ratio_dict[beta][i].SetMarkerSize(self.marker_size)
            self.h_ratio_dict[beta][i].Draw('PE same')
            
        line_lambda_np_ratio = ROOT.TLine(self.lambda_np[R], self.ymin_ratio, self.lambda_np[R], self.ymax_ratio)
        line_lambda_np_ratio.SetLineColor(self.colors[-1])
        line_lambda_np_ratio.SetLineStyle(2)
        line_lambda_np_ratio.SetLineWidth(2)
        line_lambda_np_ratio.Draw()
        self.plot_list.append(line_lambda_np_ratio)

    #-------------------------------------------------------------------------------------------
    # Construct ratio data/theory as TGraph
    # Fills:
    #   - self.h_ratio_list with histogram of ratio with stat uncertainties
    #   - self.g_theory_dict with tgraph of ratio with sys uncertainties
    #-------------------------------------------------------------------------------------------
    def construct_ratio(self, h, h_sys, h_theory_cent, h_theory_min, h_theory_max, n, x, xerrup, xerrdn, y, yerrup, yerrdn, R, beta, i, pad):

        # Construct central value
        h_ratio = h.Clone()
        h_ratio.SetName('{}_{}_{}_{}'.format(h_ratio.GetName(), R, beta, pad))
        h_ratio.SetDirectory(0)
        h_ratio.Divide(h_theory_cent)
        self.plot_list.append(h_ratio)
        self.h_ratio_dict[beta].append(h_ratio)
        y_ratio = np.array([h_ratio.GetBinContent(i) for i in range(1, n+1)])
        
        # Construct systematic uncertainties: combine data and theory uncertainties
        
        # Get relative systematic from data
        y_data = np.array([h.GetBinContent(i) for i in range(1, n+1)])
        y_sys_data = np.array([h_sys.GetBinError(i) for i in range(1, n+1)])
        y_sys_data_relative = np.divide(y_sys_data, y_data)
        
        # Get relative systematics from theory
        yerr_up_relative = np.divide(yerrup, y)
        yerr_dn_relative = np.divide(yerrdn, y)
        
        # Combine systematics in quadrature
        y_sys_total_up_relative = np.sqrt( np.square(y_sys_data_relative) + np.square(yerr_up_relative))
        y_sys_total_dn_relative = np.sqrt( np.square(y_sys_data_relative) + np.square(yerr_dn_relative))
        
        y_sys_total_up = np.multiply(y_sys_total_up_relative, y_ratio)
        y_sys_total_dn = np.multiply(y_sys_total_dn_relative, y_ratio)
        
        # Note: invert direction of asymmetric uncertainty
        g_ratio = ROOT.TGraphAsymmErrors(n, x, y_ratio, xerrdn, xerrup, y_sys_total_up, y_sys_total_dn)
        g_ratio.SetName('g_ratio_{}_{}'.format(i, beta))
        self.g_ratio_dict[beta].append(g_ratio)
        
    #-------------------------------------------------------------------------------------------
    # Rebin theory histogram according to data binning
    # Set statistical uncertainties to 0 (since we will neglect them)
    #-------------------------------------------------------------------------------------------
    def rebin_theory(self, h_theory, h):
    
        n = h.GetNbinsX()
        bins = []
        for bin in range(1, n+2):
            bins.append(h.GetBinLowEdge(bin))
        xbins = array('d', bins)
        
        h_theory_rebinned = h_theory.Rebin(n, '{}_rebinned'.format(h_theory.GetName()), xbins)
        h_theory_rebinned.SetDirectory(0)
        h_theory_rebinned.Scale(1., 'width')
        
        for bin in range(1, n+2):
            h_theory_rebinned.SetBinError(bin, 0)
    
        return h_theory_rebinned

    #-------------------------------------------------------------------------------------------
    # Scale vertical amplitude of histogram, for visualization
    #-------------------------------------------------------------------------------------------
    def scale_histogram_for_visualization(self, h, beta, groomed):

        scale_factor = 1.
        if groomed:
            if beta == '1.5':
                scale_factor = self.scale_factor_groomed_beta15
            if beta == '2':
                scale_factor = self.scale_factor_groomed_beta2
            if beta == '3':
                scale_factor = self.scale_factor_groomed_beta3       
        else:
            if beta == '2':
                scale_factor = self.scale_factor_ungroomed_beta2
            if beta == '3':
                scale_factor = self.scale_factor_ungroomed_beta3

        h.Scale(scale_factor)

        if math.isclose(scale_factor, 1.):
            plot_label = ''
        else:
            plot_label = ' (#times{})'.format(scale_factor)
        
        return plot_label

    #-------------------------------------------------------------------------------------------
    # Set legend parameters
    #-------------------------------------------------------------------------------------------
    def setupLegend(self, leg, textSize):

        leg.SetTextFont(42);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetFillColor(0);
        leg.SetMargin(0.25);
        leg.SetTextSize(textSize);
        leg.SetEntrySeparation(0.5);

    #---------------------------------------------------------------
    # Remove periods from a label
    #---------------------------------------------------------------
    def remove_periods(self, text):

        string = str(text)
        return string.replace('.', '')

    #-------------------------------------------------------------------------------------------
    def setOptions(self):

        font = 42

        ROOT.gStyle.SetFrameBorderMode(0)
        ROOT.gStyle.SetFrameFillColor(0)
        ROOT.gStyle.SetCanvasBorderMode(0)
        ROOT.gStyle.SetPadBorderMode(0)
        ROOT.gStyle.SetPadColor(10)
        ROOT.gStyle.SetCanvasColor(10)
        ROOT.gStyle.SetTitleFillColor(10)
        ROOT.gStyle.SetTitleBorderSize(1)
        ROOT.gStyle.SetStatColor(10)
        ROOT.gStyle.SetStatBorderSize(1)
        ROOT.gStyle.SetLegendBorderSize(1)

        ROOT.gStyle.SetDrawBorder(0)
        ROOT.gStyle.SetTextFont(font)
        ROOT.gStyle.SetStatFont(font)
        ROOT.gStyle.SetStatFontSize(0.05)
        ROOT.gStyle.SetStatX(0.97)
        ROOT.gStyle.SetStatY(0.98)
        ROOT.gStyle.SetStatH(0.03)
        ROOT.gStyle.SetStatW(0.3)
        ROOT.gStyle.SetTickLength(0.02,"y")
        ROOT.gStyle.SetEndErrorSize(3)
        ROOT.gStyle.SetLabelSize(0.05,"xyz")
        ROOT.gStyle.SetLabelFont(font,"xyz")
        ROOT.gStyle.SetLabelOffset(0.01,"xyz")
        ROOT.gStyle.SetTitleFont(font,"xyz")
        ROOT.gStyle.SetTitleOffset(1.2,"xyz")
        ROOT.gStyle.SetTitleSize(0.045,"xyz")
        ROOT.gStyle.SetMarkerSize(1)
        ROOT.gStyle.SetPalette(1)

        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
if __name__ == '__main__':
    print('Executing plot_angularity.py...')
    print('')

    # Define arguments
    parser = argparse.ArgumentParser(description='Plot angularity')
    parser.add_argument(
        '-o',
        '--outputDir',
        action='store',
        type=str,
        metavar='outputDir',
        default='.',
        help='Output directory for output to be written to'
    )

    # Parse the arguments
    args = parser.parse_args()

    analysis = PlotAngularityFigures(output_dir=args.outputDir)
    analysis.plot_results()
