"""
  macro for plotting multi-paneled angularity figures
  """

# General
import os
import sys
import yaml
import argparse

# Data analysis and plotting
import ROOT
import numpy as np

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

        self.base_dir = '/home/james/plot-angularity/'
        self.file = 'ang/final_results/fFinalResults.root'
        self.beta_list = [1, 1.5, 2, 3]

        self.xmin = -0.01
        self.xmax = 0.65
        self.xtitle =  '#it{#lambda}_{#it{#beta}}'
        self.ytitle = '#frac{1}{#it{#sigma}_{jet}} #frac{d#it{#sigma}}{d#it{#lambda}_{#it{#beta}}}'

        self.left_offset = 0.2

        #------------------------------------------------------

        self.markers = [21, 20, 34, 33, 22, 23]
        self.alpha = 1.
        self.colors = [ROOT.kViolet-8, ROOT.kBlue-7, ROOT.kRed-7, ROOT.kTeal-8]

        #------------------------------------------------------
        # Store paths to all final results in a dictionary
        self.predictions = {'0.2': [], '0.4': []}
        subdirs = [x for _,x,_ in os.walk(self.base_dir)][0]
        for subdir in subdirs:
            subdir_path = os.path.join(self.base_dir, subdir)
            filename = os.path.join(subdir_path, 'ang/final_results/fFinalResults.root')
            if 'R02' in subdir:
                self.predictions['0.2'].append(filename)
            if 'R04' in subdir:
                self.predictions['0.4'].append(filename)
            continue

        print(self)

    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    def plot_results(self):

        self.setOptions()
        ROOT.gROOT.ForceStyle()

        self.plot_multipanel(R=0.2, groomed=False)
        #self.plot_multipanel(R=0.2, groomed=True)

    #-------------------------------------------------------------------------------------------
    def plot_multipanel(self, R=1, groomed=False):

        # Create multi-panel canvas
        cname = 'c'
        c = ROOT.TCanvas(cname,cname,1800,1600)
        c.SetRightMargin(0.05);
        c.SetLeftMargin(self.left_offset);
        c.SetTopMargin(0.05);
        c.SetBottomMargin(0.);
        c.cd()
        c.Divide(2, 2, 0.01, 0.)

        # Keep histograms in memory, otherwise there can be problems with double deletes (i.e. ROOT then python deletes)
        self.plot_list = []

        # Plot each pt bin in its own pad
        self.plot_angularity(c, pad=1, R=R, ptbin=1, pt_label='20-40', groomed=groomed)
        self.plot_angularity(c, pad=2, R=R, ptbin=2, pt_label='40-60', groomed=groomed)
        self.plot_angularity(c, pad=3, R=R, ptbin=3, pt_label='60-80', groomed=groomed)
        self.plot_angularity(c, pad=4, R=R, ptbin=4, pt_label='80-100', groomed=groomed)

        output_filename = os.path.join(self.output_dir, 'hJetAngularity_{}{}'.format(R, self.file_format))
        c.SaveAs(output_filename)

    #-------------------------------------------------------------------------------------------
    # Get beta histograms from file, and call plot_beta_overlay to draw them
    #-------------------------------------------------------------------------------------------
    def plot_angularity(self, c, pad=0, R=1, ptbin=0, pt_label='', groomed=False):
        
        filename = self.predictions[str(R)][ptbin-1]
        f = ROOT.TFile(filename, 'READ')

        h = None
        h_sys = None
        h_pythia = None
        self.h_list = []
        self.h_sys_list = []
        self.h_pythia_list = []
        self.h_ratio_list = []
        self.blank_histo_list = []

        for i,beta in enumerate(self.beta_list):

            h_name ='hmain_ang_R{}_{}_{}_trunc'.format(R, beta, pt_label)
            h_sys_name =  'hResult_ang_systotal_R{}_{}_n3_{}'.format(R, beta, pt_label)
            h_pythia_name = 'hPythia_ang_R{}_{}_{}'.format(R, beta, pt_label)

            h = f.Get(h_name)
            h_sys = f.Get(h_sys_name)
            h_pythia = f.Get(h_pythia_name)
            h.SetDirectory(0)
            h_sys.SetDirectory(0)
            h_pythia.SetDirectory(0)

            self.h_list.append(h)
            self.h_sys_list.append(h_sys)
            self.h_pythia_list.append(h_pythia)

            # Plot the ratio
            output_filename = os.path.join(self.output_dir, 'hRatio_R{}_pt{}_{}{}'.format(R, ptbin, beta, self.file_format))
            h_ratio = self.plot_ratio(h_pythia, h, output_filename, self.xtitle, self.ytitle, R, pt_label, beta)
            h_ratio.SetName('{}_{}_{}_{}_{}'.format(h_ratio.GetName(), R, ptbin, beta, pad))
            self.plot_list.append(h_ratio)
            self.h_ratio_list.append(h_ratio)
            
        f.Close()

        # Plot overlay of beta values
        self.plot_beta_overlay(c, pad, R)

        # Keep histograms in memory
        self.plot_list.append(self.h_list)
        self.plot_list.append(self.h_sys_list)
        self.plot_list.append(self.h_pythia_list)
        self.plot_list.append(self.h_ratio_list)
        self.plot_list.append(self.blank_histo_list)

    #-------------------------------------------------------------------------------------------
    # Draw beta histograms in given pad
    #-------------------------------------------------------------------------------------------
    def plot_beta_overlay(self, c, pad, R):

        # Create canvas
        c.cd(pad)

        # Set pad and histo arrangement
        myPad = ROOT.TPad("myPad_{}".format(R), "The pad{}".format(R),0,0,1,1)
        self.plot_list.append(myPad)
        if pad in [1,3]:
            myPad.SetLeftMargin(self.left_offset)
        else:
            myPad.SetLeftMargin(0.)
        myPad.SetRightMargin(0.)
        myPad.SetTopMargin(0.01)
        myPad.SetBottomMargin(0.15)
        myPad.SetTicks(0,1)
        #myPad.SetLogy()
        myPad.Draw()
        myPad.cd()

        # Draw blank histos
        blankname = 'myBlankHisto_{}_{}'.format(pad, R)
        myBlankHisto = ROOT.TH1F(blankname,blankname, 1, self.xmin, self.xmax)
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(self.xtitle)
        myBlankHisto.GetXaxis().SetTitleSize(0.06)
        myBlankHisto.GetXaxis().SetTitleOffset(1.1)
        myBlankHisto.SetYTitle(self.ytitle)
        myBlankHisto.GetYaxis().SetTitleSize(0.06)
        myBlankHisto.GetYaxis().SetTitleOffset(1.4)
        myBlankHisto.SetMinimum(0.)
        if pad in [1,2]:
            myBlankHisto.SetMaximum(19.99)
        else:
            myBlankHisto.SetMaximum(34.99)
        myBlankHisto.Draw()
        self.blank_histo_list.append(myBlankHisto)

        scale_factor = 1.
        if pad in [1,3]:
            shift = 0.0
            shift2 = 0.0
        else:
            shift = -0.08
            shift2 = -0.08 - shift
        
        leg = ROOT.TLegend(0.7+shift,0.93-0.042*scale_factor*4,0.85,0.96)
        self.setupLegend(leg,0.03)
        self.plot_list.append(leg)

        # Draw data
        for i,beta in enumerate(self.beta_list):

            self.h_list[i].SetMarkerColorAlpha(self.colors[i], self.alpha)
            self.h_list[i].SetLineColorAlpha(self.colors[i], self.alpha)
            self.h_list[i].SetLineWidth(2)
            self.h_list[i].SetMarkerStyle(self.markers[i])
            self.h_list[i].Draw('PE same')

            self.h_sys_list[i].SetLineColor(0)
            self.h_sys_list[i].SetFillColor(self.colors[i])
            self.h_sys_list[i].SetFillColorAlpha(self.colors[i], 0.3)
            self.h_sys_list[i].SetFillStyle(1001)
            self.h_sys_list[i].SetLineWidth(0)
            self.h_sys_list[i].Draw('E2 same')

            self.h_pythia_list[i].SetLineColor(self.colors[i])
            self.h_pythia_list[i].SetLineColorAlpha(self.colors[i], 0.5)
            self.h_pythia_list[i].SetLineWidth(4)
            self.h_pythia_list[i].Draw('L hist same')

            leg.AddEntry(self.h_list[i],'#beta = {}'.format(beta),'P')

        for i,beta in enumerate(self.beta_list):
            self.h_list[i].Draw('PE same')

        leg.Draw('same')

        #line = ROOT.TLine(self.h_list[i].GetXaxis().GetXmin(),1,self.h_list[i].GetXaxis().GetXmax(),1)
        #line.SetLineColor(1)
        #line.SetLineStyle(2)
        #line.Draw('same')
        #self.plot_list.append(line)

        # # # # # # # # # # # # # # # # # # # # # # # #
        # text
        # # # # # # # # # # # # # # # # # # # # # # # #
        if pad == 1:
            ymax = 0.93
            dy = 0.05
            x = 0.35 + shift + shift2
            system0 = ROOT.TLatex(x,ymax,'#bf{ALICE}')
            system0.SetNDC()
            system0.SetTextSize(0.04*scale_factor)
            system0.Draw()

            system1 = ROOT.TLatex(x,ymax-dy,'pp  #sqrt{#it{s}} = 5.02 TeV')
            system1.SetNDC()
            system1.SetTextSize(0.04*scale_factor)
            system1.Draw()

            system2 = ROOT.TLatex(x,ymax-2*dy,'Charged jets   anti-#it{k}_{T}')
            system2.SetNDC()
            system2.SetTextSize(0.04*scale_factor)
            system2.Draw()

            system3 = ROOT.TLatex(x,ymax-3*dy, '#it{{R}} = {}    | #it{{#eta}}_{{jet}}| = {}'.format(R, 0.9-R))
            system3.SetNDC()
            system3.SetTextSize(0.04*scale_factor)
            system3.Draw()

            self.plot_list.append(system0)
            self.plot_list.append(system1)
            self.plot_list.append(system2)
            self.plot_list.append(system3)

    #-------------------------------------------------------------------------------------------
    # Plot ratio h1/h2
    def plot_ratio(self, h_pp, h_AA, outputFilename, xtitle, ytitle, R, pt_label, beta):

        # Create canvas
        cname = 'cratio_{}_{}_{}'.format(R, pt_label, beta)
        c = ROOT.TCanvas(cname,cname,800,850)
        ROOT.SetOwnership(c, False) # For some reason this is necessary to avoid a segfault...some bug in ROOT or pyroot
                                    # Supposedly fixed in https://github.com/root-project/root/pull/3787
        c.cd()
        pad1 = ROOT.TPad('pad1', 'pad1', 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.SetLeftMargin(0.2)
        pad1.SetRightMargin(0.05)
        pad1.SetTopMargin(0.05)
        pad1.SetLogy()
        pad1.Draw()
        pad1.cd()

        # Set pad and histo arrangement
        myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
        myPad.SetLeftMargin(0.22)
        myPad.SetTopMargin(0.04)
        myPad.SetRightMargin(0.04)
        myPad.SetBottomMargin(0.15)

        # Set spectra styles
        h_AA.SetMarkerSize(3)
        h_AA.SetMarkerStyle(33)
        h_AA.SetMarkerColor(600-6)
        h_AA.SetLineStyle(1)
        h_AA.SetLineWidth(2)
        h_AA.SetLineColor(600-6)

        h_pp.SetMarkerSize(2)
        h_pp.SetMarkerStyle(21)
        h_pp.SetMarkerColor(1)
        h_pp.SetLineStyle(1)
        h_pp.SetLineWidth(2)
        h_pp.SetLineColor(1)

        # Draw spectra
        h_pp.SetXTitle(xtitle)
        h_pp.GetYaxis().SetTitleOffset(2.2)
        h_pp.SetYTitle(ytitle)
        h_pp.SetMaximum(h_pp.GetMaximum()*10.)
        h_pp.SetMinimum(h_pp.GetMinimum()/10.)

        h_pp.Draw('PE X0 same')
        h_AA.Draw('PE X0 same')

        # # # # # # # # # # # # # # # # # # # # # # # #
        # Add legends and text
        # # # # # # # # # # # # # # # # # # # # # # # #
        c.cd()
        pad2 = ROOT.TPad('pad2', 'pad2', 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.35)
        pad2.SetLeftMargin(0.2)
        pad2.SetRightMargin(0.05)
        pad2.Draw()
        pad2.cd()

        # plot ratio
        hRatio = h_AA.Clone()
        hRatio.SetName('hRatio_{}'.format(cname))
        hRatio.Divide(h_pp)
        hRatio.SetMarkerStyle(21)
        hRatio.SetMarkerSize(2)

        hRatio.GetXaxis().SetTitleSize(30)
        hRatio.GetXaxis().SetTitleFont(43)
        hRatio.GetXaxis().SetTitleOffset(4.)
        hRatio.GetXaxis().SetLabelFont(43)
        hRatio.GetXaxis().SetLabelSize(20)
        hRatio.GetXaxis().SetTitle(xtitle)

        hRatio.GetYaxis().SetTitle('#it{R}_{AA}')
        hRatio.GetYaxis().SetTitleSize(20)
        hRatio.GetYaxis().SetTitleFont(43)
        hRatio.GetYaxis().SetTitleOffset(2.2)
        hRatio.GetYaxis().SetLabelFont(43)
        hRatio.GetYaxis().SetLabelSize(20)
        hRatio.GetYaxis().SetNdivisions(505)

        hRatio.SetMinimum(0.)
        hRatio.SetMaximum(1.49)
        hRatio.Draw('P E')

        line = ROOT.TLine(hRatio.GetXaxis().GetXmin(), 1, hRatio.GetXaxis().GetXmax(), 1)
        line.SetLineColor(920+2)
        line.SetLineStyle(2)
        line.SetLineWidth(4)
        line.Draw()

        #c.SaveAs(outputFilename)

        return hRatio.Clone('{}_{}'.format(hRatio.GetName(), 'new'))

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
