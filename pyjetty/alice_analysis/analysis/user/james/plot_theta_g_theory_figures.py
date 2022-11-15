"""
macro for plotting multi-paneled angularity theory figures
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

# Prevent segmentation fault from C code (doesn't seem to help?)
sys.settrace

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class PlotAngularityFigures(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, input_dir='', output_dir='', **kwargs):
        super(PlotAngularityFigures, self).__init__(**kwargs)

        self.output_dir = output_dir
        self.file_format = '.pdf'

        #------------------------------------------------------

        self.base_dir = '/home/james/pyjetty/pyjetty/alice_analysis/analysis/user/james/ezra_theta_g/'
        self.data_file = 'fFinalResults_theta_g_data.root'
        self.theory_file = 'folded_scet_calculations.root'
        self.R_list = [0.4]
        self.pt_list = [60, 80]
        self.folding_labels = ['PYTHIA8', 'Herwig7']

        # Note: logx doesn't work well, since lowest bins are different, and dominate some plots
        self.logx = False 
        # Note: logy is also not great, since the right-hand tails drop to different values
        self.logy = False 

        # Grooming settings
        self.sd_zcut = [0.1]
        self.sd_beta = [0, 1, 2]
        #self.theory_grooming_settings = [{'sd': [self.sd_zcut, self.sd_beta]}]
        # self.utils.grooming_settings
        #self.theory_grooming_labels = [self.utils.grooming_label(gs) for\
        #                               gs in self.theory_grooming_settings]

        # Don't use theory predictions below this value (to disable, set to 0)
        self.thetag_cut = 0.

        self.xmin = -0.03
        self.ymax = 13.5
        self.ymin_ratio = 0.4
        self.ymax_ratio = 2.5

        self.xtitle = '#it{#theta}_{g}'
        self.ytitle = '#frac{1}{#it{#sigma}_{#it{#theta}_{g}>#it{#theta}_{g}^{NP}}} ' + \
                      '#frac{d#it{#sigma}}{d#it{#theta}_{g}}'

        # y title with integral in denominator
        #self.ytitle = '#frac{d#it{#sigma}}{d#it{#theta}_{g}} ' + \
        #              '#/#int_{#it{#theta}_{g}^{NP}}^{1} ' + \
        #              '#frac{d#it{#sigma}}{d#it{#theta}_{g}} d#it{#theta}_{g}'

        self.left_offset = 0.4
        self.bottom_offset = 0.15
        self.ratio_height = 0.25
        self.top_bottom_scale_factor = (1 - self.ratio_height) / \
                                       (1 - self.bottom_offset - self.ratio_height)

        #------------------------------------------------------

        self.marker_data = 21
        self.markers = [20, 34, 33, 22, 23]
        self.marker_size = 3
        self.marker_size_ratio = 2
        self.alpha = 0.7
        self.color_data = 1
        self.colors = [ROOT.kRed-7, ROOT.kTeal-8, ROOT.kViolet-8, ROOT.kBlue-9]

        #------------------------------------------------------

        # Soft scale
	# Peturbative region defined by zcut * pT * theta_g^{1 + beta} * R = Lambda
        self.theta_np = {}
        self.Lambda = 1

        for i, min_pt in list(enumerate(self.pt_list))[:-1]:
            max_pt = self.pt_list[i+1]

            # P vs NP cutoff point ~ Lambda / (pT * R) -- use avg value of pT for the bin.
            # Formula assumes that jet pT xsec falls like pT^(-5.5)
            formula_pt = (4.5/3.5)*(min_pt**-3.5 - max_pt**-3.5) / \
                         (min_pt**-4.5 - max_pt**-4.5)
            # Also scale by ~20% to account for shift in full --> charged pT spectrum
            formula_pt *= 1.2

            self.theta_np[min_pt] = { R : { zcut : { beta : (self.Lambda / (zcut * formula_pt * R))**(1 / (1 + beta)) \
                                                     for beta in self.sd_beta }
                                            for zcut in self.sd_zcut }
                                      for R in self.R_list }

            #for R in self.R_list:
            #    self.theta_np[min_pt][R] = {}
            #    for zcut in self.sd_zcut:
            #        self.theta_np[min_pt][R][zcut] = {}
            #        for beta in self.sd_beta:
            #            self.theta_np[min_pt][R][zcut][beta] = (
            #                self.Lambda / (zcut * formula_pt * R) ) ** (1 / (1 + beta)) 


        #------------------------------------------------------
        # Store paths to all final results in a dictionary
        self.predictions = {}

        for i, min_pt in list(enumerate(self.pt_list))[:-1]:
            self.predictions[min_pt] = {}
            for R in self.R_list:
                self.predictions[min_pt][str(R)] = os.path.join(self.base_dir, self.theory_file)

        print(self)

    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    def plot_results(self):

        self.setOptions()
        ROOT.gROOT.ForceStyle()

        for R in self.R_list:
            for i, min_pt in list(enumerate(self.pt_list))[:-1]:
                max_pt = self.pt_list[i+1]
                for zcut in self.sd_zcut:
                    self.plot_multipanel(R, min_pt, max_pt, zcut)

    #-------------------------------------------------------------------------------------------
    def plot_multipanel(self, R, min_pt, max_pt, zcut):

        # Create multi-panel canvas
        cname = "c_" + str(R) + '_PtBin' + str(min_pt) + '-' + str(max_pt) + '_zcut' + str(zcut)
        c = ROOT.TCanvas(cname, cname, 2400, 1400)
        c.SetRightMargin(0.05);
        c.SetLeftMargin(self.left_offset);
        c.SetTopMargin(0.03);
        c.SetBottomMargin(self.bottom_offset / 2);
        c.cd()
        c.Divide(3, 1, 0.01, 0.)

        # Keep histograms in memory, otherwise there can be problems
        #   with double deletes (i.e. ROOT then python deletes)
        self.plot_list = []
        self.g_theory_dict = { beta : [] for beta in self.sd_beta }
        self.h_ratio_dict = { beta : [] for beta in self.sd_beta }
        self.g_ratio_dict = { beta : [] for beta in self.sd_beta }

        # Plot each pt bin in its own pad
        for i, beta in enumerate(self.sd_beta):
            self.plot_theta_g(
                c, pad=i+1, R=R, min_pt=min_pt, max_pt=max_pt, zcut=zcut, beta=beta)

        outfilename = "hJet_theta_g_Theory_%s_PtBin%i-%i_zcut%s" % \
            (str(R), min_pt, max_pt, str(zcut).replace('.', ''))
        outfilename += self.file_format
        output_filename = os.path.join(self.output_dir, outfilename)
        c.SaveAs(output_filename)

    #-------------------------------------------------------------------------------------------
    # Get alpha histograms from file, and call plot_alpha_overlay to draw them
    #-------------------------------------------------------------------------------------------
    def plot_theta_g(self, c, pad=0, R=1, min_pt=60, max_pt=80, zcut=0.2, beta=0):
        
        filename = os.path.join(self.base_dir, self.data_file)
        f_data = ROOT.TFile(filename, 'READ')

        self.h = None
        self.h_sys = None
        self.blank_histo_list = []

        grooming_label = 'SD_zcut%s_B%s' % (str(zcut).replace('.',''), beta)

        # Get data hist
        h_name ='hmain_theta_g_R{}_{}_{}-{}'.format(
            R, grooming_label, min_pt, max_pt)
        h_sys_name =  'hResult_theta_g_systotal_R{}_{}_{}-{}'.format(
            R, grooming_label, min_pt, max_pt)
        self.h = f_data.Get(h_name)
        self.h_sys = f_data.Get(h_sys_name)
        self.h.SetDirectory(0)
        self.h_sys.SetDirectory(0)

        # Remove negative bin edges (which has the tagging fraction)
        self.h = self.remove_negative_bin_edges(self.h)
        self.h_sys = self.remove_negative_bin_edges(self.h_sys)

        # Normalize such that integral in perturbative region is 1
        n = self.h.GetNbinsX()
        min_bin = self.h.FindBin(self.theta_np[min_pt][R][zcut][beta]) + 1
        #if self.h.GetBinCenter(min_bin) <= self.lambda_np_groomed[min_pt][R][alpha]:
        #    min_bin += 1
        if min_bin > n:
            min_bin = n
        integral_perturbative = self.h.Integral(min_bin, n, 'width')
        self.h.Scale(1./integral_perturbative)
        self.h_sys.Scale(1./integral_perturbative)

        #self.ymax = 2 * self.h.GetMaximum() * (1 + (pad == 3) * 0.2)

        self.scale_label = self.scale_histogram_for_visualization(
            self.h_sys, R, min_pt, zcut, beta)
        self.scale_histogram_for_visualization(self.h, R, min_pt, zcut, beta)

        filename = os.path.join(self.base_dir, self.theory_file)
        f_theory = ROOT.TFile(filename, 'READ')

        # Get folded theory predictions
        for i, folding_label in enumerate(self.folding_labels):

            name_cent = 'h1_folded_theta_g_ch_MPIon_R{}_{}_{}_pT_{}_{}'.format(
                str(R).replace('.', ''), grooming_label, folding_label, min_pt, max_pt)
            name_min = 'h1_min_folded_theta_g_ch_MPIon_R{}_{}_{}_pT_{}_{}'.format(
                str(R).replace('.', ''), grooming_label, folding_label, min_pt, max_pt)
            name_max = 'h1_max_folded_theta_g_ch_MPIon_R{}_{}_{}_pT_{}_{}'.format(
                str(R).replace('.', ''), grooming_label, folding_label, min_pt, max_pt)

            h_theory_cent = f_theory.Get(name_cent)
            h_theory_min = f_theory.Get(name_min)
            h_theory_max = f_theory.Get(name_max)
            h_theory_cent.SetDirectory(0)
            h_theory_min.SetDirectory(0)
            h_theory_max.SetDirectory(0)
            
            # Rebin to data binning
            h_theory_cent_rebinned = self.rebin_theory(h_theory_cent, self.h)
            h_theory_min_rebinned = self.rebin_theory(h_theory_min, self.h)
            h_theory_max_rebinned = self.rebin_theory(h_theory_max, self.h)

            h_theory_cent_rebinned = self.remove_negative_bin_edges(h_theory_cent_rebinned)
            h_theory_min_rebinned = self.remove_negative_bin_edges(h_theory_min_rebinned)
            h_theory_max_rebinned = self.remove_negative_bin_edges(h_theory_max_rebinned)

            # Normalize such that integral in perturbative region is 1
            min_bin = h_theory_cent_rebinned.FindBin(self.theta_np[min_pt][R][zcut][beta]) + 1
            #if self.h.GetBinCenter(min_bin) <= self.theta_np[min_pt][R][zcut][beta]:
            #    min_bin += 1
            if min_bin > n:
                # At least normalize by one bin (the last bin)
                min_bin = n
            integral_perturbative_theory =  h_theory_cent_rebinned.Integral(min_bin, n, 'width')
            integral_total =  h_theory_cent_rebinned.Integral(1, n, 'width')

            h_theory_cent_rebinned.Scale(1./integral_perturbative_theory)
            h_theory_min_rebinned.Scale(1./integral_perturbative_theory)
            h_theory_max_rebinned.Scale(1./integral_perturbative_theory)

            # Scale additionally for visualization
            self.scale_histogram_for_visualization(h_theory_cent_rebinned, R, min_pt, zcut, beta)
            self.scale_histogram_for_visualization(h_theory_min_rebinned, R, min_pt, zcut, beta)
            self.scale_histogram_for_visualization(h_theory_max_rebinned, R, min_pt, zcut, beta)

            n_theory_bins = h_theory_cent_rebinned.GetNbinsX()
            x = np.array([h_theory_cent_rebinned.GetXaxis().GetBinCenter(i) \
                          for i in range(1, n_theory_bins+1)])
            y = np.array([h_theory_cent_rebinned.GetBinContent(i) \
                          for i in range(1, n_theory_bins+1)])
            xerrup = xerrdn = np.array([0. for i in range(n_theory_bins)])
            yerrup = np.array([h_theory_max_rebinned.GetBinContent(i)-y[i-1] \
                               for i in range(1, n_theory_bins+1)])
            yerrdn = np.array([y[i-1]-h_theory_min_rebinned.GetBinContent(i) \
                               for i in range(1, n_theory_bins+1)])
            g_theory = ROOT.TGraphAsymmErrors(n_theory_bins, x, y, xerrdn, xerrup, yerrdn, yerrup)
            g_theory_name = 'g_theory_%i_%s' % (i, grooming_label)
            g_theory.SetNameTitle(g_theory_name, g_theory_name)
            self.g_theory_dict[beta].append(g_theory)

            # Construct ratios in self.h_ratio_dict, self.g_ratio_dict
            self.construct_ratio(self.h, self.h_sys, h_theory_cent_rebinned, h_theory_min_rebinned,
                                 h_theory_max_rebinned, n_theory_bins, x, xerrup, xerrdn, y, yerrup, yerrdn,
                                 R, zcut, beta, i, pad)

        f_data.Close()
        f_theory.Close()

        # Plot overlay of alpha values
        self.plot_beta_overlay(c, pad, R, zcut, beta, min_pt, max_pt)

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
    def plot_beta_overlay(self, c, pad, R, zcut, beta, min_pt, max_pt):

        self.logy = False

        # Create canvas
        c.cd(pad)

        grooming_label = "zcut%s_B%s" % (str(zcut).replace('.',''), beta)

        # Set pad to plot distributions
        setattr(self, "pad1_{}_PtBin{}-{}_{}_{}".format(R, min_pt, max_pt, pad, grooming_label), 
                ROOT.TPad("pad1_{}_PtBin{}-{}_{}_{}".format(R, min_pt, max_pt, pad, grooming_label),
                          "pad1_{}_PtBin{}-{}_{}_{}".format(R, min_pt, max_pt, pad, grooming_label),
                          0, self.bottom_offset + self.ratio_height,1,1))
        pad1 = getattr(self, "pad1_{}_PtBin{}-{}_{}_{}".format(R, min_pt, max_pt, pad, grooming_label))
        self.plot_list.append(pad1)
        if pad in [1]:
            pad1.SetLeftMargin(self.left_offset)
        else:
            pad1.SetLeftMargin(0.)
        pad1.SetRightMargin(0.)
        pad1.SetTopMargin(0.0)
        pad1.SetBottomMargin(0.)
        pad1.SetTicks(2,2)
        if self.logy:
            pad1.SetLogy()
        if self.logx:
            pad1.SetLogx()
        pad1.Draw()
        pad1.cd()

        # Draw blank histos
        blankname = 'myBlankHisto_{}_PtBin{}-{}_{}'.format(pad, min_pt, max_pt, R)
        xmax = self.h.GetXaxis().GetBinUpEdge(self.h.GetXaxis().GetNbins())
        myBlankHisto = ROOT.TH1F(blankname,blankname, 1, self.xmin, xmax)
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetYTitle(self.ytitle)
        myBlankHisto.GetYaxis().SetTitleSize(0.09)
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
            shift2 = -0.15 - shift

        # Legend for the center pad (2)
        leg = ROOT.TLegend(0.2+shift, 0.7, 0.6, 0.96)
        size = 0.08
        self.setupLegend(leg, size)

        self.plot_list.append(leg)

        leg3 = ROOT.TLegend(0.2+shift/2, 0.72, 0.55, 0.96)
        if pad == 2:
            self.setupLegend(leg3, size)
            self.plot_list.append(leg3)

        line_np = ROOT.TLine(self.theta_np[min_pt][R][zcut][beta], 0, 
                             self.theta_np[min_pt][R][zcut][beta], self.ymax*0.6)
        line_np.SetLineColor(self.colors[-1])
        line_np.SetLineStyle(2)
        line_np.SetLineWidth(2)
        self.plot_list.append(line_np)

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
        
        leg.AddEntry(self.h, 'Data','PE')
        
        # Draw theory
        for i, folding_label in enumerate(self.folding_labels):

            self.g_theory_dict[beta][i].SetFillColorAlpha(self.colors[i], 0.25)
            self.g_theory_dict[beta][i].SetLineColor(self.colors[i])
            self.g_theory_dict[beta][i].SetLineWidth(3)
            self.g_theory_dict[beta][i].Draw('L 3 same')

            leg.AddEntry(self.g_theory_dict[beta][i], 'NLL\' #otimes '+folding_label, 'lf')

        self.h_sys.Draw('E2 same')
        line_np.Draw()
        self.h.Draw('PE X0 same')

        if pad == 3:
            leg.Draw('same')
        if pad == 2:
            leg3.Draw('same')
        
        if pad == 2:
            leg3.AddEntry(line_np, "#it{#theta}_{g}^{NP} #leq " + \
                          "(#frac{#Lambda}{#it{z}_{cut} #it{p}_{T} #it{R}})^" + \
                          "{#frac{1}{1+#it{#beta}}}", 'lf')
            leg3.AddEntry(self.h_sys, 'Syst. uncertainty', 'f')

        # Reset for ratio plot
        self.logy = True

        # # # # # # # # # # # # # # # # # # # # # # # #
        # text
        # # # # # # # # # # # # # # # # # # # # # # # #
        ymax = 0.93
        dy = 0.07
        x = 0.45 + shift + shift2
    
        if pad == 1:
            system0 = ROOT.TLatex(x,ymax,'ALICE')
            system0.SetNDC()
            system0.SetTextSize(size / 1.2)#*scale_factor)
            system0.Draw()

            system1 = ROOT.TLatex(x,ymax-dy,'pp  #sqrt{#it{s}} = 5.02 TeV')
            system1.SetNDC()
            system1.SetTextSize(size / 1.2)#*scale_factor)
            system1.Draw()

            system2 = ROOT.TLatex(x,ymax-2*dy,'charged jets   anti-#it{k}_{T}')
            system2.SetNDC()
            system2.SetTextSize(size / 1.2)#*scale_factor)
            system2.Draw()

            system3 = ROOT.TLatex(x,ymax-3*dy,
                                  '#it{{R}} = {}    |#it{{#eta}}_{{jet}}| < {}'.format(R, 0.9-R))
            system3.SetNDC()
            system3.SetTextSize(size / 1.2)#*scale_factor)
            system3.Draw()
            
            system4 = ROOT.TLatex(x,ymax-4.*dy-0.02,
                                  str(min_pt) + ' < #it{p}_{T}^{ch jet} < ' + \
                                  str(max_pt) + ' GeV/#it{c}')
            system4.SetNDC()
            system4.SetTextSize(size / 1.2)#*scale_factor)
            system4.Draw()

            self.plot_list.append(system0)
            self.plot_list.append(system1)
            self.plot_list.append(system2)
            self.plot_list.append(system3)
            self.plot_list.append(system4)

        #system5y = ymax-6.*dy if pad == 1 else ymax-4.6*dy
        system5y = ymax-6.*dy
        system5x = x+0.35 if pad == 1 else x+0.4 if pad==2 else x+0.45
        system5 = ROOT.TLatex(system5x, system5y,
                              '#it{{#beta}} = {}{}'.format(beta, self.scale_label))
        system5.SetNDC()
        if pad in [1]:
            beta_size = size / 1.3
        else:
            beta_size = size
        system5.SetTextSize(beta_size)
        system5.Draw()
        self.plot_list.append(system5)

        if pad == 2:
            system5 = ROOT.TLatex(0.2+(shift/2), ymax-3.8*dy,
                                  'Soft drop: #it{z}_{cut} = %s' % str(zcut))
            system5.SetNDC()
            system5.SetTextSize(size)
            system5.Draw()
            self.plot_list.append(system5)

        # Set pad for ratio
        c.cd(pad)
        pad2 = ROOT.TPad("pad2_{}".format(R), "pad2{}".format(R),
                         0,0,1,self.bottom_offset+self.ratio_height)
        self.plot_list.append(pad2)
        if pad in [1]:
            pad2.SetLeftMargin(self.left_offset)
        else:
            pad2.SetLeftMargin(0.)
        pad2.SetRightMargin(0.)
        pad2.SetTopMargin(0.)
        pad2.SetBottomMargin(self.bottom_offset/(self.bottom_offset+self.ratio_height))
        pad2.SetTicks(1,2)
        if self.logy:
            pad2.SetLogy()
        if self.logx:
            pad2.SetLogx()
        pad2.Draw()
        pad2.cd()

        # Draw blank histos
        blankname = 'myBlankHisto2_{}_PtBin{}-{}_{}'.format(pad, min_pt, max_pt, R)
        myBlankHisto2 = ROOT.TH1F(blankname,blankname, 1, self.xmin, xmax-0.001)
        myBlankHisto2.SetMinimum(self.ymin_ratio)
        myBlankHisto2.SetMaximum(self.ymax_ratio)
        myBlankHisto2.SetNdivisions(510, "y")
        myBlankHisto2.SetNdivisions(505, "x")
        myBlankHisto2.SetXTitle(self.xtitle)
        myBlankHisto2.SetYTitle('#frac{Data}{Theory}')
        myBlankHisto2.GetXaxis().SetTitleSize(0.15)
        myBlankHisto2.GetXaxis().SetTitleOffset(0.7)
        myBlankHisto2.GetXaxis().SetLabelSize(0.1)
        myBlankHisto2.GetYaxis().SetTitleSize(0.12)
        myBlankHisto2.GetYaxis().SetTitleOffset(1.)
        myBlankHisto2.GetYaxis().SetLabelSize(0.1)
        myBlankHisto2.Draw()
        self.blank_histo_list.append(myBlankHisto2)

        line = ROOT.TLine(self.xmin,1,xmax,1)
        line.SetLineColor(1)
        line.SetLineStyle(2)
        line.Draw('same')
        self.plot_list.append(line)

        if pad in [1]:
            # Add y-axis numerical label text because ROOT is incapable
            ymin_text = ROOT.TLatex(0.327, 0.42, '0.5')
            ymin_text.SetNDC()
            ymax_text = ROOT.TLatex(0.367, 0.89, '2')
            ymax_text.SetNDC()

            ymin_text.SetTextSize(0.1)
            ymin_text.Draw()
            self.plot_list.append(ymin_text)

            ymax_text.SetTextSize(0.1)
            ymax_text.Draw()
            self.plot_list.append(ymax_text)

        # Draw ratio
        for i, folding_label in enumerate(self.folding_labels):
        
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

        line_np_ratio = ROOT.TLine(
            self.theta_np[min_pt][R][zcut][beta], self.ymin_ratio,
            self.theta_np[min_pt][R][zcut][beta], self.ymax_ratio)
        line_np_ratio.SetLineColor(self.colors[-1])
        line_np_ratio.SetLineStyle(2)
        line_np_ratio.SetLineWidth(2)
        line_np_ratio.Draw()
        self.plot_list.append(line_np_ratio)

    #-------------------------------------------------------------------------------------------
    # Construct ratio data/theory as TGraph
    # Fills:
    #   - self.h_ratio_list with histogram of ratio with stat uncertainties
    #   - self.g_theory_dict with tgraph of ratio with sys uncertainties
    #-------------------------------------------------------------------------------------------
    def construct_ratio(self, h, h_sys, h_theory_cent, h_theory_min, h_theory_max, n, x,
                        xerrup, xerrdn, y, yerrup, yerrdn, R, zcut, beta, i, pad):

        grooming_label = "zcut%s_B%s" % (str(zcut).replace('.', ''), str(beta))

        # Construct central value
        h_trim = self.trim_data(h, h_theory_cent)
        h_ratio = h_trim.Clone()
        h_ratio.SetName('{}_{}_{}_{}'.format(h_ratio.GetName(), R, grooming_label, pad))
        h_ratio.SetDirectory(0)
        h_ratio.Divide(h_theory_cent)
        self.plot_list.append(h_ratio)
        self.h_ratio_dict[beta].append(h_ratio)
        y_ratio = np.array([h_ratio.GetBinContent(i) for i in range(1, h_trim.GetNbinsX()+1)])

        # Construct systematic uncertainties: combine data and theory uncertainties

        # Get relative systematic from data
        y_data = np.array([h_trim.GetBinContent(i) for i in range(1, h_trim.GetNbinsX()+1)])
        h_sys_trim = self.trim_data(h_sys, h_theory_cent)
        y_sys_data = np.array([h_sys_trim.GetBinError(i) for i in range(1, h_trim.GetNbinsX()+1)])
        y_sys_data_relative = np.divide(y_sys_data, y_data)

        # Get relative systematics from theory
        yerr_up_relative = np.divide(yerrup, y)
        yerr_dn_relative = np.divide(yerrdn, y)

        # Trim according to theta_g cut
        #if self.thetag_cut != 0:
        #    yerr_up_relative = [yerr_up_relative[i] for i in range(len(yerr_up_relative)) \
        #                        if len(yerr_up_relative)-len(y_data) <= i]
        #    yerr_dn_relative = [yerr_dn_relative[i] for i in range(len(yerr_dn_relative)) \
        #                        if len(yerr_dn_relative)-len(y_data) <= i]

        # Combine systematics in quadrature
        y_sys_total_up_relative = np.sqrt( np.square(y_sys_data_relative) + \
                                           np.square(yerr_up_relative))
        y_sys_total_dn_relative = np.sqrt( np.square(y_sys_data_relative) + \
                                           np.square(yerr_dn_relative))

        y_sys_total_up = np.multiply(y_sys_total_up_relative, y_ratio)
        y_sys_total_dn = np.multiply(y_sys_total_dn_relative, y_ratio)

        # Note: invert direction of asymmetric uncertainty
        g_ratio = ROOT.TGraphAsymmErrors(n, x, y_ratio, xerrdn, xerrup,
                                         y_sys_total_up, y_sys_total_dn)
        g_ratio.SetName('g_ratio_{}_{}'.format(i, grooming_label))
        self.g_ratio_dict[beta].append(g_ratio)

    #-------------------------------------------------------------------------------------------
    # Rebin theory histogram according to data binning
    # Set statistical uncertainties to 0 (since we will neglect them)
    #-------------------------------------------------------------------------------------------
    def rebin_theory(self, h_theory, h):

        xbins = array('d', [h.GetBinLowEdge(bi) for bi in range(1, h.GetNbinsX()+2) \
           if h.GetBinLowEdge(bi) >= self.thetag_cut])
        n = len(xbins) - 1

        for bi in range(1, h_theory.GetNbinsX()+2):
            # Undo scaling by bin width before rebinning
            h_theory.SetBinContent(bi, h_theory.GetBinContent(bi) * h_theory.GetBinWidth(bi))

        h_theory_rebinned = h_theory.Rebin(n, h_theory.GetName()+'_rebinned', xbins)
        h_theory_rebinned.SetDirectory(0)
        h_theory_rebinned.Scale(1., 'width')
        
        for bi in range(1, n+2):
            h_theory_rebinned.SetBinError(bi, 0)
    
        return h_theory_rebinned

    #-------------------------------------------------------------------------------------------
    # Rebin data histogram according to theory binning
    # Ensures that theta_g cut does not affect ratio plots
    #-------------------------------------------------------------------------------------------
    def trim_data(self, h, h_theory):

        xbins = array('d', [h_theory.GetBinLowEdge(bi) \
                      for bi in range(1, h_theory.GetNbinsX()+2)])
        n = len(xbins) - 1
        h_trimmed = h.Rebin(n, h.GetName()+'_trimmed', xbins)
        return h_trimmed

    #-------------------------------------------------------------------------------------------
    # Remove bins from histogram that have negative edges for plotting purposes
    #-------------------------------------------------------------------------------------------
    def remove_negative_bin_edges(self, h):

        n = h.GetNbinsX()
        xbins = array('d', [h.GetBinLowEdge(bi) for bi in range(1, n+2)])

        while xbins[0] < 0:
            xbins = xbins[1:]
            n -= 1

        h_rebinned = h.Rebin(n, h.GetName()+'_negrm', xbins)
        h_rebinned.SetDirectory(0)

        return h_rebinned

    #-------------------------------------------------------------------------------------------
    # Scale vertical amplitude of histogram, for visualization
    #-------------------------------------------------------------------------------------------
    def scale_histogram_for_visualization(self, h, R, min_pt, zcut, beta):

        scale_factor = 1.
        ''' Can implement formatting for theta_g later
        if groomed:
            if R == 0.2:
                if min_pt == 20:
                    if alpha == '3':
                        scale_factor = 0.1
                #elif min_pt == 40:
            #elif R == 0.4:
        else:
            if R == 0.2:
                if min_pt == 20:
                    if alpha == '2':
                        scale_factor = 0.25
                    elif alpha == '3':
                        scale_factor = 0.02
                elif min_pt == 40:
                    if alpha == '2':
                        scale_factor = 0.45
                    elif alpha == '3':
                        scale_factor = 0.06
                elif min_pt == 60:
                    if alpha == '2':
                        scale_factor = 0.3
                    elif alpha == '3':
                        scale_factor = 0.12
                elif min_pt == 80:
                    if alpha == '2':
                        scale_factor = 0.42
                    if alpha == '3':
                        scale_factor = 0.035
            elif R == 0.4:
                if min_pt == 20:
                    if alpha == '2':
                        scale_factor = 0.45
                    elif alpha == '3':
                        scale_factor = 0.2
                elif min_pt == 40:
                    if alpha == '2':
                        scale_factor = 0.45
                    elif alpha == '3':
                        scale_factor = 0.25
                elif min_pt == 60:
                    if alpha == '2':
                        scale_factor = 0.65
                    elif alpha == '3':
                        scale_factor = 0.27
                elif min_pt == 80:
                    if alpha == '2':
                        scale_factor = 0.33
                    elif alpha == '3':
                        scale_factor = 0.15
        '''

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
        '-i',
        '--inputDir',
        action='store',
        type=str,
        metavar='inputDir',
        default='.',
        help='input directory containing ROOT file'
    )
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

    analysis = PlotAngularityFigures(input_dir=args.inputDir, output_dir=args.outputDir)
    analysis.plot_results()
