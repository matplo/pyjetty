"""
macro for plotting multi-paneled data-pQCD figures
"""

# General
import os
import sys
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
class PlotSubjetZFigures(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, output_dir='', **kwargs):
        super(PlotSubjetZFigures, self).__init__(**kwargs)

        self.output_dir = output_dir
        self.file_format = '.pdf'

        #------------------------------------------------------

        self.observable = 'inclusive_subjet_z'
        self.base_dir = f'/home/james/pyjetty/pyjetty/alice_analysis/analysis/user/james/subjet_z/{self.observable}'
        self.data_file = 'fFinalResults.root'
        self.theory_file = 'folded_scet_R04_pT_80_120.root'
        self.R_list = [0.4]
        self.r_list = [0.1, 0.2]
        self.pt_list = [80, 120]
        self.folding_labels = ['PYTHIA8', 'Herwig7']
        self.folding_plot_labels = ['PYTHIA8', 'HERWIG7']
        self.z_leading_min = 0.7001

        self.xmin = -0.03
        self.xmax = 1.03
        if self.observable ==  'inclusive_subjet_z':
            self.logy = True 
            self.ymin = 0.5
            self.ymax = 2000
            self.ymin_ratio = 0.1
            self.ymax_ratio = 20
        if self.observable ==  'leading_subjet_z':
            self.logy = True 
            self.ymin = 0.16
            self.ymax = 1400
            self.ymin_ratio = 0.3
            self.ymax_ratio = 12.9

        self.xtitle = '#it{z}_{#it{r}}'
        self.ytitle = '#frac{1}{ #it{#sigma}_{ 0.7 < #it{z}_{#it{r}} < #it{z}_{#it{r}}^{NP} } } ' + \
                      '#frac{d#it{#sigma}}{d#it{z}_{#it{r}}}'

        self.left_offset = 0.26
        self.bottom_offset = 0.15
        self.ratio_height = 0.3
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
        # Peturbative region defined by Lambda=(1-z)pt*r
        self.zr_np = {}
        self.Lambda = 1
        for i, min_pt in list(enumerate(self.pt_list))[:-1]:
            max_pt = self.pt_list[i+1]

            # P vs NP cutoff point: zr ~ 1-Lambda/(pt*r) -- use avg value of pT for the bin.
            # Formula assumes that jet pT xsec falls like pT^(-5.5)
            formula_pt = (4.5/3.5)*(min_pt**-3.5 - max_pt**-3.5) / \
                         (min_pt**-4.5 - max_pt**-4.5)
            # Also scale by ~20% to account for shift in full --> charged pT spectrum
            formula_pt *= 1.2

            self.zr_np[min_pt] = { R : { r : (1 - self.Lambda / (formula_pt * r)) \
                                        for r in self.r_list }
                                        for R in self.R_list }

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
                self.plot_multipanel(R, min_pt, max_pt)

    #-------------------------------------------------------------------------------------------
    def plot_multipanel(self, R, min_pt, max_pt):

        # Create multi-panel canvas
        cname = "c_" + str(R) + '_PtBin' + str(min_pt) + '-' + str(max_pt)
        c = ROOT.TCanvas(cname, cname, 2000, 1400)
        c.SetRightMargin(0.1);
        c.SetLeftMargin(self.left_offset);
        c.SetTopMargin(0.03);
        c.SetBottomMargin(self.bottom_offset / 2);
        c.cd()
        c.Divide(2, 1, 0.01, 0.)

        # Keep histograms in memory, otherwise there can be problems
        #   with double deletes (i.e. ROOT then python deletes)
        self.plot_list = []
        self.g_theory_dict = { r : [] for r in self.r_list }
        self.h_ratio_dict = { r : [] for r in self.r_list }
        self.h_ratio_sys_dict = { r : [] for r in self.r_list }
        self.g_ratio_dict = { r : [] for r in self.r_list }

        # Plot each pt bin in its own pad
        for i, r in enumerate(self.r_list):
            self.plot_subjet_z(
                c, pad=i+1, R=R, min_pt=min_pt, max_pt=max_pt, r=r)

        outfilename = f"hJet_{self.observable}_Theory_{R}_PtBin{min_pt}-{max_pt}{self.file_format}"
        output_filename = os.path.join(self.output_dir, outfilename)
        c.SaveAs(output_filename)

    #-------------------------------------------------------------------------------------------
    # Get alpha histograms from file, and call plot_alpha_overlay to draw them
    #-------------------------------------------------------------------------------------------
    def plot_subjet_z(self, c, pad=0, R=1, min_pt=60, max_pt=80, r=0.):
        
        filename = os.path.join(self.base_dir, self.data_file)
        f_data = ROOT.TFile(filename, 'READ')

        self.h = None
        self.h_sys = None
        self.blank_histo_list = []

        # Get data hist
        h_name =f'hmain_{self.observable}_R{R}_{r}_{min_pt}-{max_pt}'
        h_sys_name = f'hResult_{self.observable}_systotal_R{R}_{r}_{min_pt}-{max_pt}'
        self.h = f_data.Get(h_name)
        self.h_sys = f_data.Get(h_sys_name)
        self.h.SetDirectory(0)
        self.h_sys.SetDirectory(0)

        # Normalize such that integral in perturbative region is 1
        n = self.h.GetNbinsX()
        min_bin = self.h.FindBin(self.z_leading_min)             
        max_bin = self.h.FindBin(self.zr_np[min_pt][R][r]) - 1 # Take bin below zr_NP
        integral_total = self.h.Integral(1, n, 'width')
        integral_perturbative = self.h.Integral(min_bin, max_bin, 'width')
        print(f'integral_total, leading_subjet_z, r={r}: {integral_total}')
        print(f'integral_perturbative, leading_subjet_z, r={r}: {integral_perturbative}')
        self.h.Scale(1./integral_perturbative)
        self.h_sys.Scale(1./integral_perturbative)

        # Note: we don't compute z_loss or <N_subjets>, since we do not include 
        # the NP region in the normalization

        filename = os.path.join(self.base_dir, self.theory_file)
        f_theory = ROOT.TFile(filename, 'READ')

        # Get folded theory predictions
        for i, folding_label in enumerate(self.folding_labels):

            name_cent = 'h1_folded_{}_ch_MPIon_R{}_{}_{}_pT_{}_{}'.format(
                self.observable, str(R).replace('.', ''), r, folding_label, min_pt, max_pt)
            name_min = 'h1_min_folded_{}_ch_MPIon_R{}_{}_{}_pT_{}_{}'.format(
                self.observable, str(R).replace('.', ''), r, folding_label, min_pt, max_pt)
            name_max = 'h1_max_folded_{}_ch_MPIon_R{}_{}_{}_pT_{}_{}'.format(
                self.observable, str(R).replace('.', ''), r, folding_label, min_pt, max_pt)

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

            # Scale to the leading subjet perturbative region 0.5<zr<zr_NP
            n = self.h.GetNbinsX()
            min_bin = h_theory_cent_rebinned.FindBin(self.z_leading_min)
            max_bin = h_theory_cent_rebinned.FindBin(self.zr_np[min_pt][R][r]) - 1 # Take bin below zr_NP
            integral_perturbative_theory = h_theory_cent_rebinned.Integral(min_bin, max_bin, 'width')
            print(f'integral_perturbative, leading_subjet_z, r={r}: {integral_perturbative_theory}')

            h_theory_cent_rebinned.Scale(1./integral_perturbative_theory)
            h_theory_min_rebinned.Scale(1./integral_perturbative_theory)
            h_theory_max_rebinned.Scale(1./integral_perturbative_theory)

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
            g_theory_name = 'g_theory_%i_%s' % (i, r)
            g_theory.SetNameTitle(g_theory_name, g_theory_name)
            self.g_theory_dict[r].append(g_theory)

            # Construct ratios in self.h_ratio_dict, self.h_ratio_sys_dict, self.g_ratio_dict
            self.construct_ratio(self.h, self.h_sys, h_theory_cent_rebinned,
                                 n_theory_bins, x, xerrup, xerrdn, y, yerrup, yerrdn, R, r, pad)

        f_data.Close()
        f_theory.Close()

        # Plot overlay of alpha values
        self.plot_beta_overlay(c, pad, R,r, min_pt, max_pt)

        # Keep histograms in memory
        self.plot_list.append(self.h)
        self.plot_list.append(self.h_sys)
        self.plot_list.append(self.g_theory_dict)
        self.plot_list.append(self.h_ratio_dict)
        self.plot_list.append(self.h_ratio_sys_dict)
        self.plot_list.append(self.g_ratio_dict)
        self.plot_list.append(self.blank_histo_list)

    #-------------------------------------------------------------------------------------------
    # Draw beta histograms in given pad
    #-------------------------------------------------------------------------------------------
    def plot_beta_overlay(self, c, pad, R, r, min_pt, max_pt):

        # Create canvas
        c.cd(pad)

        grooming_label = f'{r}'

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
        pad1.SetTicks(0,1)
        if self.logy:
            pad1.SetLogy()
        pad1.Draw()
        pad1.cd()

        # Draw blank histos
        blankname = 'myBlankHisto_{}_PtBin{}-{}_{}'.format(pad, min_pt, max_pt, R)
        myBlankHisto = ROOT.TH1F(blankname,blankname, 1, self.xmin, self.xmax)
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetYTitle(self.ytitle)
        myBlankHisto.GetYaxis().SetTitleSize(0.09)
        myBlankHisto.GetYaxis().SetTitleOffset(1.2)
        myBlankHisto.GetYaxis().SetLabelSize(0.08)
        myBlankHisto.SetMinimum(self.ymin)
        myBlankHisto.SetMaximum(self.ymax)

        myBlankHisto.Draw()
        self.blank_histo_list.append(myBlankHisto)

        if pad in [1]:
            shift = 0.0
            shift2 = 0.0
        else:
            if self.observable == 'leading_subjet_z':
                shift = 0.0
                shift2 = 0.0
            elif self.observable == 'inclusive_subjet_z':
                shift = -0.08
                shift2 = -0.15 - shift

        # Legend for the center pad (2)
        if self.observable == 'leading_subjet_z':
            x_legend = 0.1
        elif self.observable == 'inclusive_subjet_z':
            x_legend = 0.23
        leg = ROOT.TLegend(x_legend+shift, 0.55, 0.6, 0.96)
        size = 0.07
        self.setupLegend(leg, size, sep=0.)
        self.plot_list.append(leg)

        line_np = ROOT.TLine(self.zr_np[min_pt][R][r], 0, 
                             self.zr_np[min_pt][R][r], self.ymax)
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
        
        leg.AddEntry(self.h,'ALICE','PE')

        # Draw theory
        for i, folding_label in enumerate(self.folding_plot_labels):

            self.g_theory_dict[r][i].SetFillColorAlpha(self.colors[i], 0.25)
            self.g_theory_dict[r][i].SetLineColor(self.colors[i])
            self.g_theory_dict[r][i].SetLineWidth(3)
            self.g_theory_dict[r][i].Draw('L 3 same')

            leg.AddEntry(self.g_theory_dict[r][i], 'NLL\' #otimes '+folding_label, 'lf')

        self.h_sys.Draw('E2 same')
        line_np.Draw()
        self.h.Draw('PE X0 same')

        if pad == 2:
            leg.AddEntry(self.h_sys, 'Syst. uncertainties', 'f')
            leg.AddEntry(line_np, "#it{z}_{#it{r}}^{NP} = " + \
                          "1 #font[122]{-} #frac{#Lambda}{ #it{p}_{T} #it{z}_{#it{r}} }", 'lf')
            leg.Draw('same')

        # # # # # # # # # # # # # # # # # # # # # # # #
        # text
        # # # # # # # # # # # # # # # # # # # # # # # #
        ymax = 0.92
        dy = 0.09
        if self.observable == 'leading_subjet_z':
            x = 0.31 + shift + shift2
        elif self.observable == 'inclusive_subjet_z':
            x = 0.36 + shift + shift2
    
        if pad == 1:

            system0 = ROOT.TLatex(x,ymax,'ALICE')
            system0.SetNDC()
            system0.SetTextSize(size)
            system0.Draw()

            system1 = ROOT.TLatex(x,ymax-dy,'pp  #sqrt{#it{s}} = 5.02 TeV')
            system1.SetNDC()
            system1.SetTextSize(size)
            system1.Draw()

            system2 = ROOT.TLatex(x,ymax-2*dy,'Charged-particle anti-#it{k}_{T} jets')
            system2.SetNDC()
            system2.SetTextSize(size)
            system2.Draw()

            system3 = ROOT.TLatex(x,ymax-3*dy,
                                  '#it{{R}} = {}    |#it{{#eta}}_{{jet}}| < {}'.format(R, 0.9-R))
            system3.SetNDC()
            system3.SetTextSize(size)
            system3.Draw()
            
            system4 = ROOT.TLatex(x,ymax-4.*dy-0.02,
                                  str(min_pt) + ' < #it{p}_{T}^{ch jet} < ' + \
                                  str(max_pt) + ' GeV/#it{c}')
            system4.SetNDC()
            system4.SetTextSize(size)
            system4.Draw()

            self.plot_list.append(system0)
            self.plot_list.append(system1)
            self.plot_list.append(system2)
            self.plot_list.append(system3)
            self.plot_list.append(system4)

            if self.observable == 'leading_subjet_z':
                system5y = ymax-5.3*dy
                system5x = x - shift - shift2
                system5 = ROOT.TLatex(system5x, system5y, 'Leading anti-#it{k}_{T} subjets')
            elif self.observable == 'inclusive_subjet_z':
                system5y = ymax-5.3*dy
                system5x = x - shift - shift2
                system5 = ROOT.TLatex(system5x, system5y, 'Inclusive anti-#it{k}_{T} subjets')
            system5.SetNDC()
            system5.SetTextSize(size)
            system5.Draw()
            self.plot_list.append(system5)

            if self.observable == 'inclusive_subjet_z':
                system5x += 0.2

        else:
            if self.observable == 'leading_subjet_z':
                system5x = x_legend
            elif self.observable == 'inclusive_subjet_z':
                system5x = x_legend + 0.2
                    
        system6y = ymax-7.*dy
        system6 = ROOT.TLatex(system5x, system6y, f'#it{{r}} = {r}')
        system6.SetNDC()
        system6.SetTextSize(size)
        system6.Draw()
        self.plot_list.append(system6)

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
        pad2.Draw()
        pad2.cd()

        # Draw blank histos
        blankname = 'myBlankHisto2_{}_PtBin{}-{}_{}'.format(pad, min_pt, max_pt, R)
        myBlankHisto2 = ROOT.TH1F(blankname,blankname, 1, self.xmin, self.xmax)
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

        line = ROOT.TLine(self.xmin, 1, self.xmax, 1)
        line.SetLineColor(1)
        line.SetLineStyle(2)
        line.Draw('same')
        self.plot_list.append(line)

        # Draw ratio
        for i, folding_label in enumerate(self.folding_labels):
        
            # Draw tgraph with sys uncertainty
            self.g_ratio_dict[r][i].SetFillColorAlpha(self.colors[i], 0.25)
            self.g_ratio_dict[r][i].SetLineColor(self.colors[i])
            self.g_ratio_dict[r][i].SetLineWidth(3)
            self.g_ratio_dict[r][i].Draw('3 same')
        
            # Draw th1 with sys uncertainty
            self.h_ratio_sys_dict[r][i].SetLineColor(0)
            self.h_ratio_sys_dict[r][i].SetMarkerSize(0)
            self.h_ratio_sys_dict[r][i].SetMarkerColor(0)
            self.h_ratio_sys_dict[r][i].SetFillColor(self.color_data)
            self.h_ratio_sys_dict[r][i].SetFillColorAlpha(self.color_data, 0.3)
            self.h_ratio_sys_dict[r][i].SetFillStyle(1001)
            self.h_ratio_sys_dict[r][i].SetLineWidth(0)
            self.h_ratio_sys_dict[r][i].Draw('E2 same')

            # Draw th1 with stat uncertainty
            self.h_ratio_dict[r][i].SetMarkerColor(self.colors[i])
            self.h_ratio_dict[r][i].SetLineColor(self.colors[i])
            self.h_ratio_dict[r][i].SetFillColor(self.colors[i])
            self.h_ratio_dict[r][i].SetLineWidth(1)
            self.h_ratio_dict[r][i].SetMarkerStyle(self.marker_data)
            self.h_ratio_dict[r][i].SetMarkerSize(self.marker_size)
            self.h_ratio_dict[r][i].Draw('PE X0 same')

        line_np_ratio = ROOT.TLine(
            self.zr_np[min_pt][R][r], self.ymin_ratio,
            self.zr_np[min_pt][R][r], self.ymax_ratio)
        line_np_ratio.SetLineColor(self.colors[-1])
        line_np_ratio.SetLineStyle(2)
        line_np_ratio.SetLineWidth(2)
        line_np_ratio.Draw()
        self.plot_list.append(line_np_ratio)

    #-------------------------------------------------------------------------------------------
    # Construct ratio data/theory as TGraph
    #  Plot data and theory uncertainties separately
    # Fills:
    #   - self.h_ratio_dict with histogram of ratio with stat uncertainties
    #   - self.h_ratio_sys_dict with histogram of ratio with sys uncertainties
    #   - self.g_theory_dict with tgraph of ratio with sys uncertainties
    #-------------------------------------------------------------------------------------------
    def construct_ratio(self, h, h_sys, h_theory_cent, n, x,
                        xerrup, xerrdn, y, yerrup, yerrdn, R, r, pad):

        grooming_label = f'{r}'

        # Construct ratio central value (and data statistical uncertainties)
        h_trim = self.trim_data(h, h_theory_cent)
        h_ratio = h_trim.Clone()
        h_ratio.SetName('{}_{}_{}_{}'.format(h_ratio.GetName(), R, grooming_label, pad))
        h_ratio.SetDirectory(0)
        h_ratio.Divide(h_theory_cent)
        self.plot_list.append(h_ratio)
        self.h_ratio_dict[r].append(h_ratio)

        # Construct ratio with data systematic uncertainties
        h_sys_trim = self.trim_data(h_sys, h_theory_cent)
        h_sys_ratio = h_sys_trim.Clone()
        h_sys_ratio.SetName('{}_{}_{}_{}'.format(h_sys_ratio.GetName(), R, grooming_label, pad))
        h_sys_ratio.SetDirectory(0)
        h_sys_ratio.Divide(h_theory_cent)
        self.plot_list.append(h_sys_ratio)
        self.h_ratio_sys_dict[r].append(h_sys_ratio)

        # Get relative systematics from theory, and plot at ratio=1
        yerr_up_relative = np.divide(yerrup, y)
        yerr_dn_relative = np.divide(yerrdn, y)
        xerrup = xerrdn = np.array([0. for i in range(n)])
        y_ratio = np.ones(n)
        g_ratio = ROOT.TGraphAsymmErrors(n, x, y_ratio, xerrdn, xerrup, yerr_dn_relative, yerr_up_relative)
        g_ratio.SetName(f'g_ratio_{r}')
        self.g_ratio_dict[r].append(g_ratio)

    #-------------------------------------------------------------------------------------------
    # Rebin theory histogram according to data binning
    # Set statistical uncertainties to 0 (since we will neglect them)
    #-------------------------------------------------------------------------------------------
    def rebin_theory(self, h_theory, h):

        xbins = array('d', [h.GetBinLowEdge(bi) for bi in range(1, h.GetNbinsX()+2)])
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
    # Ensures that NP cut does not affect ratio plots
    #-------------------------------------------------------------------------------------------
    def trim_data(self, h, h_theory):

        xbins = array('d', [h_theory.GetBinLowEdge(bi) \
                      for bi in range(1, h_theory.GetNbinsX()+2)])
        n = len(xbins) - 1
        h_trimmed = h.Rebin(n, h.GetName()+'_trimmed', xbins)
        return h_trimmed

    #-------------------------------------------------------------------------------------------
    # Set legend parameters
    #-------------------------------------------------------------------------------------------
    def setupLegend(self, leg, textSize, sep=0.5):

        leg.SetTextFont(42);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetFillColor(0);
        leg.SetMargin(0.25);
        leg.SetTextSize(textSize);
        leg.SetEntrySeparation(sep);

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
        ROOT.gStyle.SetTickLength(0.04,"x")
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

    analysis = PlotSubjetZFigures(output_dir=args.outputDir)
    analysis.plot_results()