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
        self.xtitle_groomed =  '#it{#lambda}_{#it{#beta},g}'
        self.ytitle_groomed = '#frac{1}{#it{#sigma}_{jet}} #frac{d#it{#sigma}}{d#it{#lambda}_{#it{#beta},g}}'
        self.logy = False

        self.left_offset = 0.2
        self.bottom_offset = 0.15
        self.ratio_height = 0.25
        self.top_bottom_scale_factor = (1 - self.ratio_height) / (1 - self.bottom_offset - self.ratio_height)

        #------------------------------------------------------

        self.markers = [21, 20, 34, 33, 22, 23]
        self.marker_size = 2.5
        self.marker_size_ratio = 2
        self.alpha = 0.7
        self.colors = [ROOT.kBlue+3, ROOT.kBlue-2, ROOT.kBlue-6, ROOT.kBlue-9]
        #ROOT.kViolet-8, ROOT.kBlue-4, ROOT.kRed-7, ROOT.kTeal-8]

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
        self.plot_multipanel(R=0.4, groomed=False)
        self.plot_multipanel(R=0.2, groomed=True)
        self.plot_multipanel(R=0.4, groomed=True)

    #-------------------------------------------------------------------------------------------
    def plot_multipanel(self, R=1, groomed=False):

        # Create multi-panel canvas
        cname = 'c{}{}'.format(R,groomed)
        c = ROOT.TCanvas(cname,cname,1800,1600)
        c.SetRightMargin(0.05);
        c.SetLeftMargin(self.left_offset);
        c.SetTopMargin(0.05);
        c.SetBottomMargin(self.bottom_offset/2);
        c.cd()
        c.Divide(2, 2, 0.01, 0.)

        # Keep histograms in memory, otherwise there can be problems with double deletes (i.e. ROOT then python deletes)
        self.plot_list = []

        # Plot each pt bin in its own pad
        self.plot_angularity(c, pad=1, R=R, ptbin=1, minpt=20, maxpt=40, groomed=groomed)
        self.plot_angularity(c, pad=2, R=R, ptbin=2, minpt=40, maxpt=60, groomed=groomed)
        self.plot_angularity(c, pad=3, R=R, ptbin=3, minpt=60, maxpt=80, groomed=groomed)
        self.plot_angularity(c, pad=4, R=R, ptbin=4, minpt=80, maxpt=100, groomed=groomed)

        if groomed:
            outfilename = 'hJetAngularity_{}_SD{}'.format(R, self.file_format)
        else:
            outfilename = 'hJetAngularity_{}{}'.format(R, self.file_format)
        output_filename = os.path.join(self.output_dir, outfilename)
        c.SaveAs(output_filename)

    #-------------------------------------------------------------------------------------------
    # Get beta histograms from file, and call plot_beta_overlay to draw them
    #-------------------------------------------------------------------------------------------
    def plot_angularity(self, c, pad=0, R=1, ptbin=0, minpt=0, maxpt=0, groomed=False):
        
        filename = self.predictions[str(R)][ptbin-1]
        f = ROOT.TFile(filename, 'READ')

        h = None
        h_sys = None
        h_pythia = None
        self.h_list = []
        self.h_sys_list = []
        self.h_pythia_list = []
        self.h_ratio_list = []
        self.h_ratio_sys_list = []
        self.blank_histo_list = []

        grooming_label = ''
        if groomed:
            grooming_label = '_SD_zcut02_B0'

        for i,beta in enumerate(self.beta_list):



            h_name ='hmain_ang_R{}_{}{}_{}-{}_trunc'.format(R, beta, grooming_label, minpt, maxpt)
            h_sys_name =  'hResult_ang_systotal_R{}_{}{}_n3_{}-{}'.format(R, beta, grooming_label, minpt, maxpt)
            h_pythia_name = 'hPythia_ang_R{}_{}{}_{}-{}'.format(R, beta, grooming_label, minpt, maxpt)

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
            h_ratio = h.Clone()
            h_ratio.SetName('{}_{}_{}_{}_{}'.format(h_ratio.GetName(), R, ptbin, beta, pad))
            h_ratio.SetDirectory(0)
            h_ratio.Divide(h_pythia)
            self.plot_list.append(h_ratio)
            self.h_ratio_list.append(h_ratio)

            h_ratio_sys = h_sys.Clone()
            h_ratio_sys.SetName('sys{}_{}_{}_{}_{}'.format(h_ratio.GetName(), R, ptbin, beta, pad))
            h_ratio_sys.SetDirectory(0)
            h_ratio_sys.Divide(h_pythia)
            self.plot_list.append(h_ratio_sys)
            self.h_ratio_sys_list.append(h_ratio_sys)

        f.Close()

        # Plot overlay of beta values
        self.plot_beta_overlay(c, pad, R, minpt, maxpt, groomed)

        # Keep histograms in memory
        self.plot_list.append(self.h_list)
        self.plot_list.append(self.h_sys_list)
        self.plot_list.append(self.h_pythia_list)
        self.plot_list.append(self.h_ratio_list)
        self.plot_list.append(self.h_ratio_sys_list)
        self.plot_list.append(self.blank_histo_list)

    #-------------------------------------------------------------------------------------------
    # Draw beta histograms in given pad
    #-------------------------------------------------------------------------------------------
    def plot_beta_overlay(self, c, pad, R, minpt, maxpt, groomed):

        # Create canvas
        c.cd(pad)

        # Set pad to plot distributions
        if pad in [3,4]:
            pad1 = ROOT.TPad("pad1_{}".format(R), "pad1{}".format(R),0, self.bottom_offset + self.ratio_height,1,1)
        else:
            pad1 = ROOT.TPad("pad1_{}".format(R), "pad1{}".format(R),0, self.ratio_height,1,1)
        self.plot_list.append(pad1)
        if pad in [1,3]:
            pad1.SetLeftMargin(self.left_offset)
        else:
            pad1.SetLeftMargin(0.)
        pad1.SetRightMargin(0.)
        pad1.SetTopMargin(0.0)
        pad1.SetBottomMargin(0.)
        pad1.SetTicks(0,1)
        if self.logy:
            pad1.SetLogy()
            #pad1.SetLogx()
        pad1.Draw()
        pad1.cd()

        # Draw blank histos
        blankname = 'myBlankHisto_{}_{}'.format(pad, R)
        myBlankHisto = ROOT.TH1F(blankname,blankname, 1, self.xmin, self.xmax)
        myBlankHisto.SetNdivisions(505)
        if groomed:
            myBlankHisto.SetYTitle(self.ytitle_groomed)
        else:
            myBlankHisto.SetYTitle(self.ytitle)    
        if pad in [1,2]:
            myBlankHisto.GetYaxis().SetTitleSize(0.1)
            myBlankHisto.GetYaxis().SetTitleOffset(0.7)
            myBlankHisto.GetYaxis().SetLabelSize(0.06)
        else:
            myBlankHisto.GetYaxis().SetTitleSize(0.115)
            myBlankHisto.GetYaxis().SetTitleOffset(0.6)
            myBlankHisto.GetYaxis().SetLabelSize(0.07)
        if self.logy:
            myBlankHisto.SetMinimum(0.01)
            myBlankHisto.SetMaximum(1000)
        else:
            myBlankHisto.SetMinimum(0.001)
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
        
        leg = ROOT.TLegend(0.75+shift,0.9-0.072*scale_factor*4,0.9,0.96)
        if pad in [1,2]:
            self.setupLegend(leg,0.07)
        else:
            self.setupLegend(leg,0.08)
        self.plot_list.append(leg)

        leg2 = ROOT.TLegend(0.18+shift,0.9-0.072*scale_factor*2,0.35,0.9)
        if pad == 2:
            self.setupLegend(leg2,0.07)
            self.plot_list.append(leg2)

        # Draw data
        for i,beta in enumerate(self.beta_list):

            self.h_list[i].SetMarkerColorAlpha(self.colors[i], self.alpha)
            self.h_list[i].SetLineColorAlpha(self.colors[i], self.alpha)
            self.h_list[i].SetLineWidth(2)
            self.h_list[i].SetMarkerStyle(self.markers[i])
            self.h_list[i].SetMarkerSize(self.marker_size)
            self.h_list[i].Draw('PE same')

            self.h_sys_list[i].SetLineColor(0)
            self.h_sys_list[i].SetMarkerSize(0)
            self.h_sys_list[i].SetMarkerColor(0)
            self.h_sys_list[i].SetFillColor(self.colors[i])
            self.h_sys_list[i].SetFillColorAlpha(self.colors[i], 0.3)
            self.h_sys_list[i].SetFillStyle(1001)
            self.h_sys_list[i].SetLineWidth(0)
            self.h_sys_list[i].Draw('E2 same')

            self.h_pythia_list[i].SetLineColor(self.colors[i])
            self.h_pythia_list[i].SetLineColorAlpha(self.colors[i], 0.5)
            self.h_pythia_list[i].SetLineWidth(3)
            self.h_pythia_list[i].Draw('L hist same')

            leg.AddEntry(self.h_list[i],'#beta = {}'.format(beta),'P')
            
            if i == 2:
                leg2.AddEntry(self.h_sys_list[i], 'Sys. uncertainty', 'f')
                leg2.AddEntry(self.h_pythia_list[i], 'PYTHIA8 Monash 2013', 'l')

        for i,beta in enumerate(self.beta_list):
            self.h_list[i].Draw('PE same')

        leg.Draw('same')
        if pad == 2:
            leg2.Draw('same')

        # # # # # # # # # # # # # # # # # # # # # # # #
        # text
        # # # # # # # # # # # # # # # # # # # # # # # #
        ymax = 0.9
        dy = 0.09
        x = 0.3 + shift + shift2
        if pad in [1,2]:
            size = 0.08
        else:
            size = 0.09
    
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

            system3 = ROOT.TLatex(x,ymax-3*dy, '#it{{R}} = {}    | #it{{#eta}}_{{jet}}| = {}'.format(R, 0.9-R))
            system3.SetNDC()
            system3.SetTextSize(size*scale_factor)
            system3.Draw()

            self.plot_list.append(system0)
            self.plot_list.append(system1)
            self.plot_list.append(system2)
            self.plot_list.append(system3)

        system4 = ROOT.TLatex(x+0.2,ymax-6.*dy, str(minpt) + ' < #it{p}_{T,jet}^{ch} < ' + str(maxpt) + ' GeV/#it{c}')
        system4.SetNDC()
        system4.SetTextSize(size*scale_factor)
        system4.Draw()
        self.plot_list.append(system4)

        if groomed and pad == 1:
            system5 = ROOT.TLatex(x, ymax-4.2*dy, 'Soft Drop #it{z}_{cut} = 0.2 #beta_{SD} = 0')
            system5.SetNDC()
            system5.SetTextSize(size*scale_factor)
            system5.Draw()
            self.plot_list.append(system5)

        # Set pad for ratio
        c.cd(pad)
        if pad in [3,4]:
            pad2 = ROOT.TPad("pad2_{}".format(R), "pad2{}".format(R),0,0,1,self.bottom_offset+self.ratio_height)
        else:
            pad2 = ROOT.TPad("pad2_{}".format(R), "pad2{}".format(R),0,0,1,self.ratio_height)
        self.plot_list.append(pad2)
        if pad in [1,3]:
            pad2.SetLeftMargin(self.left_offset)
        else:
            pad2.SetLeftMargin(0.)
        pad2.SetRightMargin(0.)
        pad2.SetTopMargin(0.)
        if pad in [3,4]:
            pad2.SetBottomMargin(self.bottom_offset/(self.bottom_offset+self.ratio_height))
        else:
            pad2.SetBottomMargin(0.)
        pad2.SetTicks(1,1)
        pad2.Draw()
        pad2.cd()

        # Draw blank histos
        blankname = 'myBlankHisto2_{}_{}'.format(pad, R)
        myBlankHisto2 = ROOT.TH1F(blankname,blankname, 1, self.xmin, self.xmax)
        myBlankHisto2.SetNdivisions(505, "y")
        if groomed:
            myBlankHisto2.SetXTitle(self.xtitle_groomed)
        else:
            myBlankHisto2.SetXTitle(self.xtitle)
        myBlankHisto2.SetYTitle('#frac{Data}{MC}')        
        if pad in [1,2]:
            myBlankHisto2.GetYaxis().SetTitleSize(0.2)
            myBlankHisto2.GetYaxis().SetTitleOffset(0.3)
            myBlankHisto2.GetYaxis().SetLabelSize(0.2)
            myBlankHisto2.SetMinimum(0.01)
        else:
            myBlankHisto2.GetXaxis().SetTitleSize(0.2)
            myBlankHisto2.GetXaxis().SetTitleOffset(0.8)
            myBlankHisto2.GetXaxis().SetLabelSize(0.12)      
            myBlankHisto2.GetYaxis().SetTitleSize(0.12)
            myBlankHisto2.GetYaxis().SetTitleOffset(0.5)
            myBlankHisto2.GetYaxis().SetLabelSize(0.12)
            myBlankHisto2.SetMinimum(0.)
        myBlankHisto2.SetMaximum(1.99)
        myBlankHisto2.Draw()
        self.blank_histo_list.append(myBlankHisto2)

        line = ROOT.TLine(self.xmin,1,self.xmax,1)
        line.SetLineColor(1)
        line.SetLineStyle(2)
        line.Draw('same')
        self.plot_list.append(line)

        # Draw ratio
        for i,beta in enumerate(self.beta_list):

            self.h_ratio_sys_list[i].SetLineColor(0)
            self.h_ratio_sys_list[i].SetFillColor(self.colors[i])
            self.h_ratio_sys_list[i].SetFillColorAlpha(self.colors[i], 0.3)
            self.h_ratio_sys_list[i].SetFillStyle(1001)
            self.h_ratio_sys_list[i].SetLineWidth(0)
            self.h_ratio_sys_list[i].Draw('E2 same')

            self.h_ratio_list[i].SetMarkerColorAlpha(self.colors[i], self.alpha)
            self.h_ratio_list[i].SetLineColorAlpha(self.colors[i], self.alpha)
            self.h_ratio_list[i].SetLineWidth(2)
            self.h_ratio_list[i].SetMarkerStyle(self.markers[i])
            self.h_ratio_list[i].SetMarkerSize(self.marker_size)
            self.h_ratio_list[i].Draw('PE same')

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
