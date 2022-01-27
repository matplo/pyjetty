"""
  macro for plotting multi-paneled angularity figures
  """

# General
import os
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
class PlotDynamicalGroomingFigures(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, output_dir='', **kwargs):
        super(PlotDynamicalGroomingFigures, self).__init__(**kwargs)

        self.output_dir = output_dir
        self.file_format = '.pdf'

        #------------------------------------------------------
        # Store paths to all final results in a dictionary
        self.results = {}        
        base_dir = '/Users/jamesmulligan/Analysis_theta_g_pp/roounfold_output/pp/343450'
        self.results['zg'] = os.path.join(base_dir, 'zg/final_results/fFinalResults.root')
        self.results['theta_g'] = os.path.join(base_dir, 'theta_g/final_results/fFinalResults.root')

        #------------------------------------------------------

        self.observable_list = ['zg', 'theta_g']
        self.R = 0.4
        self.min_pt = 60
        self.max_pt = 80

        #------------------------------------------------------
        # Get theory predictions
        self.predictions = {'zg': {}, 'theta_g': {}}        
        with open('./alba.yaml', 'r') as stream:
            self.theory_config = yaml.safe_load(stream)

        #------------------------------------------------------

        self.left_offset = 0.2
        self.bottom_offset = 0.15
        self.ratio_height = 0.25
        self.top_bottom_scale_factor = (1 - self.ratio_height) / (1 - self.bottom_offset - self.ratio_height)

        self.marker_data = 21
        self.marker_size = 3
        self.marker_size_ratio = 2
        self.alpha = 0.7
        self.color_data = 1
        self.color_theory = ROOT.kBlue-9

        print(self)

    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    def plot_results(self):
        self.setOptions()
        ROOT.gROOT.ForceStyle()

        for observable in self.observable_list:
            print(f'Plotting {observable}...')

            self.initialize_observable(observable)
            self.plot_multipanel(observable)

    #-------------------------------------------------------------------------------------------
    def initialize_observable(self, observable):

        if observable == 'zg':
            self.xmin = -0.01
            self.xmax = 0.5
            self.ymax = 8.99
            self.ymin_ratio = 1e-3
            self.ymax_ratio = 1.99
            self.scale_factor_a2 = 1
            self.xtitle =  '#it{z}_{g}'
            self.ytitle = '#frac{1}{#sigma_{jet}} #frac{d#it{#sigma}}{d#it{z}_{g}}'
        elif observable == 'theta_g':
            self.xmin = -0.015
            self.xmax = 0.999
            self.ymax = 3.99
            self.ymin_ratio = 1e-3
            self.ymax_ratio = 1.99
            self.scale_factor_a2 = 1
            self.xtitle =  '#it{#theta}_{g}'
            self.ytitle = '#frac{1}{#sigma_{jet}} #frac{d#it{#sigma}}{d#it{#theta}_{g}}'

    #-------------------------------------------------------------------------------------------
    def plot_multipanel(self, observable):

        # Create multi-panel canvas
        cname = f'c{observable}'
        c = ROOT.TCanvas(cname,cname,2400,1400)
        c.SetRightMargin(0.05);
        c.SetLeftMargin(self.left_offset);
        c.SetTopMargin(0.03);
        c.SetBottomMargin(self.bottom_offset/2);
        c.cd()
        c.Divide(2, 1, 0.01, 0.)

        # Keep histograms in memory, otherwise there can be problems with double deletes (i.e. ROOT then python deletes)
        self.plot_list = []
        self.g_theory_dict = {}
        self.h_ratio_dict = {}
        self.h_ratio_sys_dict = {}
        self.g_ratio_dict = {}

        # Plot each pt bin in its own pad
        self.plot_observable(c, pad=1, observable=observable, a='1')
        self.plot_observable(c, pad=2, observable=observable, a='2')

        outfilename = f'h_{observable}_DyG_Theory_{self.file_format}'
        output_filename = os.path.join(self.output_dir, outfilename)
        c.SaveAs(output_filename)

    #-------------------------------------------------------------------------------------------
    # Get data histograms from file, and call plot_overlay to draw them
    #-------------------------------------------------------------------------------------------
    def plot_observable(self, c, pad=0, observable=None, a=''):
        print(f'  Plotting a = {a}...')
        
        # Reset some things
        self.h = None
        self.h_sys = None
        self.blank_histo_list = []

        # Get data hist
        filename = self.results[observable]
        f = ROOT.TFile(filename, 'READ')
        grooming_label = f'_DG_a{a}0'
        h_name = f'hmain_{observable}_R{self.R}{grooming_label}_{self.min_pt}-{self.max_pt}'
        h_sys_name =  f'hResult_{observable}_systotal_R{self.R}{grooming_label}_{self.min_pt}-{self.max_pt}'
        self.h = f.Get(h_name)
        self.h_sys = f.Get(h_sys_name)
        self.h.SetDirectory(0)
        self.h_sys.SetDirectory(0)
        f.Close()
        
        # Normalize such that integral is 1
        n = self.h.GetNbinsX()
        integral = self.h.Integral(1, n+1, 'width')
        self.h.Scale(1./integral)
        self.h_sys.Scale(1./integral)
        
        # Get theory predictions and create tgraph -- already normalized
        self.predictions[observable][a] = self.theory_config[f'{observable}_a{a}']

        bins = np.array(self.predictions[observable][a]['bins'])
        x = (bins[1:] + bins[:-1]) / 2
        central = np.array(self.predictions[observable][a]['central'])
        lower = np.array(self.predictions[observable][a]['lower'])
        upper = np.array(self.predictions[observable][a]['upper'])
        n = len(x)
        xerr = np.zeros(n)
        g_theory = ROOT.TGraphAsymmErrors(n, x, central, xerr, xerr, central-lower, upper-central)
        g_theory.SetName(f'g_theory_{observable}_{a}')
        self.g_theory_dict[a] = g_theory

        # Construct ratio data/theory (divide histogram by tgraph)
        self.construct_ratio(self.h, self.h_sys, g_theory, x, central, lower, upper, a)

        # Plot distributions in given pad
        # It will draw:
        #   - self.h
        #   - self.h_sys
        #   - self.g_theory_dict[a]
        #   - self.g_ratio_dict[a] -- ratio w/sys uncertainty
        #   - self.h_ratio_dict[a] -- ratio w/stat uncertainty
        self.plot_overlay(c, pad, a, observable)

        # Keep histograms in memory
        self.plot_list.append(self.h)
        self.plot_list.append(self.h_sys)
        self.plot_list.append(self.g_theory_dict)
        self.plot_list.append(self.h_ratio_dict)
        self.plot_list.append(self.h_ratio_sys_dict)
        self.plot_list.append(self.g_ratio_dict)
        self.plot_list.append(self.blank_histo_list)

    #-------------------------------------------------------------------------------------------
    # Construct ratio data/theory as TGraph (see alternately: analysis/base/analysis_utils.divide_tgraph(self, h, g, combine_errors=False))
    # Fills:
    #   - self.h_ratio_dict with histogram of ratio with stat uncertainties
    #   - self.h_ratio_sys_dict with histogram of ratio with sys uncertainties
    #   - self.g_theory_dict with tgraph of ratio with sys uncertainties
    #-------------------------------------------------------------------------------------------
    def construct_ratio(self, h, h_sys, g_theory, x, central, lower, upper, a):

        # Construct central value (only include stat uncertainty from data)
        h_ratio = h.Clone()
        h_ratio.SetName(f'{h_ratio.GetName()}_{a}')
        h_ratio.SetDirectory(0)
        h_ratio_sys = h.Clone()
        h_ratio_sys.SetName(f'{h_ratio_sys.GetName()}_{a}')
        h_ratio_sys.SetDirectory(0)
        n = h.GetNbinsX()
        for bin in range(1, n+1):
            # Get histogram (x,y)
            h_x = h.GetBinCenter(bin)
            h_y = h.GetBinContent(bin)
            h_error = h.GetBinError(bin)
            h_sys_error = h_sys.GetBinError(bin)

            # Get TGraph (x,y) and errors
            g_x = ROOT.Double(0)
            g_y = ROOT.Double(0)
            g_theory.GetPoint(bin-1, g_x, g_y)
            
            if not np.isclose(h_x, g_x):
                print('hist x: {}'.format(h_x))
                print('graph x: {}'.format(g_x))
                print()
            new_content = h_y / g_y
            h_ratio.SetBinContent(bin, new_content)
            h_ratio.SetBinError(bin, h_error/h_y * new_content)
            h_ratio_sys.SetBinContent(bin, new_content)
            h_ratio_sys.SetBinError(bin, h_sys_error/h_y * new_content)
        self.h_ratio_dict[a] = h_ratio
        self.h_ratio_sys_dict[a] = h_ratio_sys

        # Construct theory systematic uncertainties on ratio, and plot at ratio=1

        # Get relative systematics from theory
        yerr_up_relative = np.divide(upper-central, central)
        yerr_dn_relative = np.divide(central-lower, central)

        xerrup = xerrdn = np.array([0. for i in range(n)])
        y_ratio = np.ones(n)
        g_ratio = ROOT.TGraphAsymmErrors(n, x, y_ratio, xerrdn, xerrup, yerr_dn_relative, yerr_up_relative)
        g_ratio.SetName(f'g_ratio_{a}')
        self.g_ratio_dict[a] = g_ratio

    #-------------------------------------------------------------------------------------------
    # Draw beta histograms in given pad
    #-------------------------------------------------------------------------------------------
    def plot_overlay(self, c, pad, a, observable):

        # Create canvas
        c.cd(pad)

        # Set pad to plot distributions
        pad1 = ROOT.TPad("pad1", "pad1", 0, self.bottom_offset + self.ratio_height,1,1)
        self.plot_list.append(pad1)
        if pad in [1]:
            pad1.SetLeftMargin(self.left_offset)
        else:
            pad1.SetLeftMargin(0.)
        pad1.SetRightMargin(0.)
        pad1.SetTopMargin(0.0)
        pad1.SetBottomMargin(0.)
        pad1.SetTicks(1,1)
        pad1.Draw()
        pad1.cd()

        # Draw blank histos
        blankname = f'myBlankHisto_{pad}'
        xmax = self.h.GetXaxis().GetBinUpEdge(self.h.GetXaxis().GetNbins())
        myBlankHisto = ROOT.TH1F(blankname,blankname, 1, self.xmin, self.xmax)
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetYTitle(self.ytitle)    
        myBlankHisto.GetYaxis().SetTitleSize(0.1)
        myBlankHisto.GetYaxis().SetTitleOffset(0.8)
        myBlankHisto.GetYaxis().SetLabelSize(0.07)
        myBlankHisto.SetMinimum(0.001)
        myBlankHisto.SetMaximum(self.ymax)
        myBlankHisto.Draw()
        self.blank_histo_list.append(myBlankHisto)

        scale_factor = 1.3
        if pad in [1]:
            shift = 0.0
            shift2 = 0.0
        else:
            shift = -0.08
            shift2 = -0.08 - shift
        
        leg = ROOT.TLegend(0.3+shift,0.96-0.072*scale_factor*4,0.6,0.96)
        self.setupLegend(leg,0.07)
        self.plot_list.append(leg)

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
        self.g_theory_dict[a].SetFillColorAlpha(self.color_theory, 0.25)
        self.g_theory_dict[a].SetLineColor(self.color_theory)
        self.g_theory_dict[a].SetLineWidth(3)
        self.g_theory_dict[a].Draw('L 3 same')
        leg.AddEntry(self.g_theory_dict[a],'LO+N^{2}DL+NP','lf')        
        leg.AddEntry(self.h_sys, 'Sys. uncertainty', 'f')

        self.h_sys.Draw('E2 same')
        self.h.Draw('PE X0 same')

        if pad == 2:
            leg.Draw('same')
        
        # # # # # # # # # # # # # # # # # # # # # # # #
        # text
        # # # # # # # # # # # # # # # # # # # # # # # #
        ymax = 0.9
        dy = 0.08
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

            system3 = ROOT.TLatex(x,ymax-3*dy, '#it{{R}} = {}    | #it{{#eta}}_{{jet}}| < {}'.format(self.R, 0.9-self.R))
            system3.SetNDC()
            system3.SetTextSize(size*scale_factor)
            system3.Draw()
            
            system4 = ROOT.TLatex(x,ymax-4.*dy-0.025, str(self.min_pt) + ' < #it{p}_{T,jet}^{ch} < ' + str(self.max_pt) + ' GeV/#it{c}')
            system4.SetNDC()
            system4.SetTextSize(size*scale_factor)
            system4.Draw()

            system5 = ROOT.TLatex(x, ymax-5*dy-0.06, 'dynamical grooming')
            system5.SetNDC()
            system5.SetTextSize(size*scale_factor)
            system5.Draw()

            self.plot_list.append(system0)
            self.plot_list.append(system1)
            self.plot_list.append(system2)
            self.plot_list.append(system3)
            self.plot_list.append(system4)
            self.plot_list.append(system5)
            
        xlabel = x+0.2
        if pad in [1]:
            beta_size = 1.3*size
        else:
            beta_size = 1.3*size
            if observable == 'theta_g':
                xlabel -= 0.4
        system6 = ROOT.TLatex(xlabel,ymax-7.*dy-0.0, '#it{{a}} = {}'.format(a))
        system6.SetNDC()
        system6.SetTextSize(beta_size)
        system6.Draw()
        self.plot_list.append(system6)

        # Set pad for ratio
        c.cd(pad)
        pad2 = ROOT.TPad(f"pad2_{self.R}", f"pad2{self.R}",0,0,1,self.bottom_offset+self.ratio_height)
        self.plot_list.append(pad2)
        if pad in [1]:
            pad2.SetLeftMargin(self.left_offset)
        else:
            pad2.SetLeftMargin(0.)
        pad2.SetRightMargin(0.)
        pad2.SetTopMargin(0.)
        pad2.SetBottomMargin(self.bottom_offset/(self.bottom_offset+self.ratio_height))
        pad2.SetTicks(1,1)
        #pad2.SetLogy()
        pad2.Draw()
        pad2.cd()

        # Draw blank histos
        blankname = f'myBlankHisto2_{pad}'
        myBlankHisto2 = ROOT.TH1F(blankname,blankname, 1, self.xmin, self.xmax)
        myBlankHisto2.SetNdivisions(505, "y")
        myBlankHisto2.SetNdivisions(505, "x")
        myBlankHisto2.SetXTitle(self.xtitle)
        myBlankHisto2.SetYTitle('#frac{Data}{Theory}')
        myBlankHisto2.GetXaxis().SetTitleSize(0.15)
        myBlankHisto2.GetXaxis().SetTitleOffset(1.)
        myBlankHisto2.GetXaxis().SetLabelSize(0.1)
        myBlankHisto2.GetYaxis().SetTitleSize(0.12)
        myBlankHisto2.GetYaxis().SetTitleOffset(0.7)
        myBlankHisto2.GetYaxis().SetLabelSize(0.1)
        myBlankHisto2.SetMinimum(self.ymin_ratio)
        myBlankHisto2.SetMaximum(self.ymax_ratio)
        myBlankHisto2.Draw()
        self.blank_histo_list.append(myBlankHisto2)

        line = ROOT.TLine(self.xmin,1,self.xmax,1)
        line.SetLineColor(1)
        line.SetLineStyle(2)
        line.Draw('same')
        self.plot_list.append(line)

        # Draw ratio
        
        # Draw tgraph with sys uncertainty
        self.g_ratio_dict[a].SetFillColorAlpha(self.color_theory, 0.25)
        self.g_ratio_dict[a].SetLineColor(self.color_theory)
        self.g_ratio_dict[a].SetLineWidth(3)
        self.g_ratio_dict[a].Draw('3 same')

        # Draw th1 with sys uncertainty
        self.h_ratio_sys_dict[a].SetLineColor(0)
        self.h_ratio_sys_dict[a].SetMarkerSize(0)
        self.h_ratio_sys_dict[a].SetMarkerColor(0)
        self.h_ratio_sys_dict[a].SetFillColor(self.color_data)
        self.h_ratio_sys_dict[a].SetFillColorAlpha(self.color_data, 0.3)
        self.h_ratio_sys_dict[a].SetFillStyle(1001)
        self.h_ratio_sys_dict[a].SetLineWidth(0)
        self.h_ratio_sys_dict[a].Draw('E2 same')

        # Draw th1 with stat uncertainty
        self.h_ratio_dict[a].SetMarkerColor(self.color_data)
        self.h_ratio_dict[a].SetLineColor(self.color_data)
        self.h_ratio_dict[a].SetFillColor(self.color_theory)
        self.h_ratio_dict[a].SetLineWidth(2)
        self.h_ratio_dict[a].SetMarkerStyle(self.marker_data)
        self.h_ratio_dict[a].SetMarkerSize(self.marker_size)
        self.h_ratio_dict[a].Draw('PE X0 same')

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
    print('Executing plot_dynamical_grooming_theory_figures.py...')
    print('')

    # Define arguments
    parser = argparse.ArgumentParser(description='Plot dynamical grooming')
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

    analysis = PlotDynamicalGroomingFigures(output_dir=args.outputDir)
    analysis.plot_results()
