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
        
        self.jetR = 0.2

        self.base_dir = '/Users/jamesmulligan/Analysis_Angularity/plot-angularity/v2'
        self.file = 'fFinalResults.root'
        self.beta_list = [1, 1.5, 2, 3]

        self.logx = False # Note: logx doesn't work well, since lowest bins are different, and dominate some plots
        self.logy = False # Note: logy is also not great, since the right-hand tails drop to different values

        self.xmin = -0.01
        if self.logx:
            self.xmin = 0.001
            
        self.scale_factor_groomed_beta3_lowpt = 0.03
        self.scale_factor_groomed_beta2_lowpt = 0.2
        self.scale_factor_groomed_beta15_lowpt = 0.5
        self.scale_factor_ungroomed_beta3_lowpt = 0.5
        self.scale_factor_ungroomed_beta2_lowpt = 1.
        self.scale_factor_ungroomed_beta15_lowpt = 1.
        
        self.scale_factor_groomed_beta3_highpt = 0.02
        self.scale_factor_groomed_beta2_highpt = 0.1
        self.scale_factor_groomed_beta15_highpt = 0.25
        self.scale_factor_ungroomed_beta3_highpt = 0.3
        self.scale_factor_ungroomed_beta2_highpt = 0.5
        self.scale_factor_ungroomed_beta15_highpt = 0.7

        self.xtitle =  '#it{#lambda}_{#it{#alpha}}'
        self.ytitle = '#frac{1}{#it{#sigma}} #frac{d#it{#sigma}}{d#it{#lambda}_{#it{#alpha}}}'
        self.xtitle_groomed =  '#it{#lambda}_{#it{#alpha},g}'
        self.ytitle_groomed = '#frac{1}{#it{#sigma}_{inc}} #frac{d#it{#sigma}}{d#it{#lambda}_{#it{#alpha},g}}'

        self.left_offset = 0.2
        self.bottom_offset = 0.15
        self.ratio_height = 0.2

        #------------------------------------------------------

        self.markers = [21, 20, 34, 33, 22, 23]
        self.marker_size = 2.5
        self.marker_size_ratio = 2
        self.alpha = 0.7
        self.colors = [ROOT.kViolet-8, ROOT.kBlue-9, ROOT.kRed-7, ROOT.kTeal-8]

        #------------------------------------------------------
        # Store paths to all final results in a dictionary
        self.predictions = {'0.2': [], '0.4': []}
        subdirs = [x for _,x,_ in os.walk(self.base_dir)][0]
        for subdir in subdirs:
            subdir_path = os.path.join(self.base_dir, subdir)
            filename = os.path.join(subdir_path, self.file)
            if 'R02' in subdir:
                self.predictions['0.2'].append(filename)
            if 'R04' in subdir:
                self.predictions['0.4'].append(filename)
            continue
        self.predictions['0.2'].sort()
        self.predictions['0.4'].sort()

        print(self)

    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------
    def plot_results(self):

        self.setOptions()
        ROOT.gROOT.ForceStyle()

        if self.jetR == 0.2:
            self.plot_multipanel(R=0.2, groomed=False)
            self.plot_multipanel(R=0.2, groomed=True)
        elif self.jetR == 0.4:
            self.plot_multipanel(R=0.4, groomed=False)
            self.plot_multipanel(R=0.4, groomed=True)

    #-------------------------------------------------------------------------------------------
    def plot_multipanel(self, R=1, groomed=False):

        # Create multi-panel canvas
        cname = 'c{}{}'.format(R,groomed)
        c = ROOT.TCanvas(cname,cname,1800,2200)
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
        self.h_herwig_list = []
        self.h_pythia_ratio_list = []
        self.h_pythia_ratio_sys_list = []
        self.h_herwig_ratio_list = []
        self.h_herwig_ratio_sys_list = []
        self.blank_histo_list = []

        grooming_label = ''
        if groomed:
            grooming_label = '_SD_zcut02_B0'

        for i,beta in enumerate(self.beta_list):

            h_name ='hmain_ang_R{}_{}{}_{}-{}_trunc'.format(R, beta, grooming_label, minpt, maxpt)
            h_sys_name =  'hResult_ang_systotal_R{}_{}{}_n3_{}-{}'.format(R, beta, grooming_label, minpt, maxpt)
            h_pythia_name = 'hPythia_ang_R{}_{}{}_{}-{}'.format(R, beta, grooming_label, minpt, maxpt)
            h_herwig_name = 'hHerwig_ang_R{}_{}{}_{}-{}'.format(R, beta, grooming_label, minpt, maxpt)

            h = f.Get(h_name)
            h_sys = f.Get(h_sys_name)
            h_pythia = f.Get(h_pythia_name)
            h_herwig = f.Get(h_herwig_name)
            h.SetDirectory(0)
            h_sys.SetDirectory(0)
            h_pythia.SetDirectory(0)
            h_herwig.SetDirectory(0)

            self.h_list.append(h)
            self.h_sys_list.append(h_sys)
            self.h_pythia_list.append(h_pythia)
            self.h_herwig_list.append(h_herwig)

            # Plot the pythia ratio
            h_pythia_ratio = h.Clone()
            h_pythia_ratio.SetName('{}_{}_{}_{}_{}'.format(h_pythia_ratio.GetName(), R, ptbin, beta, pad))
            h_pythia_ratio.SetDirectory(0)
            h_pythia_ratio.Divide(h_pythia)
            self.plot_list.append(h_pythia_ratio)
            self.h_pythia_ratio_list.append(h_pythia_ratio)

            h_pythia_ratio_sys = h_sys.Clone()
            h_pythia_ratio_sys.SetName('sys{}_{}_{}_{}_{}'.format(h_pythia_ratio.GetName(), R, ptbin, beta, pad))
            h_pythia_ratio_sys.SetDirectory(0)
            h_pythia_ratio_sys.Divide(h_pythia)
            self.plot_list.append(h_pythia_ratio_sys)
            self.h_pythia_ratio_sys_list.append(h_pythia_ratio_sys)
            
            # Plot the herwig ratio
            h_herwig_ratio = h.Clone()
            h_herwig_ratio.SetName('{}_{}_{}_{}_{}'.format(h_herwig_ratio.GetName(), R, ptbin, beta, pad))
            h_herwig_ratio.SetDirectory(0)
            h_herwig_ratio.Divide(h_herwig)
            self.plot_list.append(h_herwig_ratio)
            self.h_herwig_ratio_list.append(h_herwig_ratio)

            h_herwig_ratio_sys = h_sys.Clone()
            h_herwig_ratio_sys.SetName('sys{}_{}_{}_{}_{}'.format(h_herwig_ratio.GetName(), R, ptbin, beta, pad))
            h_herwig_ratio_sys.SetDirectory(0)
            h_herwig_ratio_sys.Divide(h_herwig)
            self.plot_list.append(h_herwig_ratio_sys)
            self.h_herwig_ratio_sys_list.append(h_herwig_ratio_sys)

        f.Close()

        # Plot overlay of beta values
        self.plot_beta_overlay(c, pad, R, minpt, maxpt, groomed)

        # Keep histograms in memory
        self.plot_list.append(self.h_list)
        self.plot_list.append(self.h_sys_list)
        self.plot_list.append(self.h_pythia_list)
        self.plot_list.append(self.h_pythia_ratio_list)
        self.plot_list.append(self.h_pythia_ratio_sys_list)
        self.plot_list.append(self.h_herwig_list)
        self.plot_list.append(self.h_herwig_ratio_list)
        self.plot_list.append(self.h_herwig_ratio_sys_list)
        self.plot_list.append(self.blank_histo_list)

    #-------------------------------------------------------------------------------------------
    # Draw beta histograms in given pad
    #-------------------------------------------------------------------------------------------
    def plot_beta_overlay(self, c, pad, R, minpt, maxpt, groomed):

        if groomed:
            self.logy = True
            if pad in [2,4]:
                if R == 0.4:
                    self.xmax = 0.5999
                else:
                    self.xmax = 0.5499
            else:
                if R == 0.4:
                    self.xmax = 0.6999
                else:
                    self.xmax = 0.6499
        else:
            if pad in [2,4]:
                if R == 0.4:
                    self.xmax = 0.5999
                else:
                    self.xmax = 0.5499
            else:
                self.xmax = 0.6499


        # Create canvas
        c.cd(pad)

        # Set pad to plot distributions
        if pad in [3,4]:
            pad1 = ROOT.TPad("pad1_{}".format(R), "pad1{}".format(R),0, self.bottom_offset + 2*self.ratio_height,1,1)
        else:
            pad1 = ROOT.TPad("pad1_{}".format(R), "pad1{}".format(R),0, 2*self.ratio_height,1,1)
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
        if self.logx:
            pad1.SetLogx()
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
            if groomed:
                myBlankHisto.GetYaxis().SetTitleSize(1.05*0.1)
                myBlankHisto.GetYaxis().SetTitleOffset(0.78)
                myBlankHisto.GetYaxis().SetLabelSize(1.*0.06)
            else:
                myBlankHisto.GetYaxis().SetTitleSize(1.1*0.1)
                myBlankHisto.GetYaxis().SetTitleOffset(0.7)
                myBlankHisto.GetYaxis().SetLabelSize(1.1*0.06)
        else:
            if groomed:
                myBlankHisto.GetYaxis().SetTitleSize(1.1*0.115)
                myBlankHisto.GetYaxis().SetTitleOffset(0.65)
                myBlankHisto.GetYaxis().SetLabelSize(1.1*0.07)
            else:
                myBlankHisto.GetYaxis().SetTitleSize(1.1*0.115)
                myBlankHisto.GetYaxis().SetTitleOffset(0.6)
                myBlankHisto.GetYaxis().SetLabelSize(1.1*0.07)
        if self.logy:
            if pad in [1,2]:
                if R == 0.2:
                    myBlankHisto.SetMinimum(0.0011)
                    myBlankHisto.SetMaximum(5e7)
                elif R == 0.4:
                    myBlankHisto.SetMinimum(0.0002)
                    myBlankHisto.SetMaximum(5e7)
            else:
                if R == 0.2:
                    myBlankHisto.SetMinimum(2e-4)
                    myBlankHisto.SetMaximum(2e2)
                elif R == 0.4:
                    myBlankHisto.SetMinimum(4e-3)
                    myBlankHisto.SetMaximum(2e1)
        else:
            myBlankHisto.SetMinimum(0.001)
            if pad in [1,2]:
                if groomed:
                    myBlankHisto.SetMaximum(12.99)
                else:
                    myBlankHisto.SetMaximum(14.99)
            else:
                if groomed:
                    myBlankHisto.SetMaximum(12.99)
                else:
                    myBlankHisto.SetMaximum(9.99)
        myBlankHisto.Draw()
        self.blank_histo_list.append(myBlankHisto)

        scale_factor = 1.
        if pad in [1,3]:
            shift = 0.0
            shift2 = 0.0
        else:
            shift = -0.08
            shift2 = -0.08 - shift
        
        leg = ROOT.TLegend(0.68+shift,0.9-0.072*scale_factor*4,0.83,0.96)
        if pad in [1,2]:
            self.setupLegend(leg,0.07)
        else:
            self.setupLegend(leg,0.08)
        self.plot_list.append(leg)

        leg2 = ROOT.TLegend(0.13+shift,0.9-0.072*scale_factor*3,0.3,0.9)
        if pad == 2:
            self.setupLegend(leg2,0.06)
            self.plot_list.append(leg2)

        # Draw data
        for i,beta in enumerate(self.beta_list):

            self.h_list[i].SetMarkerColorAlpha(self.colors[i], 1)
            self.h_list[i].SetLineColorAlpha(self.colors[i], self.alpha)
            self.h_list[i].SetLineWidth(2)
            self.h_list[i].SetMarkerStyle(self.markers[i])
            self.h_list[i].SetMarkerSize(self.marker_size)
            scale_label = self.scale_histogram_for_visualization(self.h_list[i], i, minpt, groomed)
            self.h_list[i].Draw('PE X0 same')

            self.h_sys_list[i].SetLineColor(0)
            self.h_sys_list[i].SetMarkerSize(0)
            self.h_sys_list[i].SetMarkerColor(0)
            self.h_sys_list[i].SetFillColor(self.colors[i])
            self.h_sys_list[i].SetFillColorAlpha(self.colors[i], 0.3)
            self.h_sys_list[i].SetFillStyle(1001)
            self.h_sys_list[i].SetLineWidth(0)
            scale_label = self.scale_histogram_for_visualization(self.h_sys_list[i], i, minpt, groomed)
            self.h_sys_list[i].Draw('E2 same')

            self.h_pythia_list[i].SetLineColor(self.colors[i])
            self.h_pythia_list[i].SetLineColorAlpha(self.colors[i], 0.7)
            self.h_pythia_list[i].SetLineWidth(3)
            scale_label = self.scale_histogram_for_visualization(self.h_pythia_list[i], i, minpt, groomed)
            self.h_pythia_list[i].Draw('L hist same')
            
            self.h_herwig_list[i].SetLineColor(self.colors[i])
            self.h_herwig_list[i].SetLineColorAlpha(self.colors[i], 0.7)
            self.h_herwig_list[i].SetLineWidth(3)
            self.h_herwig_list[i].SetLineStyle(2)
            scale_label = self.scale_histogram_for_visualization(self.h_herwig_list[i], i, minpt, groomed)
            self.h_herwig_list[i].Draw('L hist same')

            leg.AddEntry(self.h_list[i],'#it{{#alpha}} = {}{}'.format(beta, scale_label),'P')
            
            if i == 2:
                leg2.AddEntry(self.h_sys_list[i], 'Syst. uncertainty', 'f')
                leg2.AddEntry(self.h_pythia_list[i], 'PYTHIA8 Monash 2013', 'l')
                leg2.AddEntry(self.h_herwig_list[i], 'Herwig7', 'l')

        for i,beta in enumerate(self.beta_list):
            self.h_list[i].Draw('PE X0 same')

        leg.Draw('same')
        if pad == 2:
            leg2.Draw('same')

        # # # # # # # # # # # # # # # # # # # # # # # #
        # text
        # # # # # # # # # # # # # # # # # # # # # # # #
        ymax = 0.9
        dy = 0.09
        x = 0.24 + shift + shift2
        if pad in [1,2]:
            size = 0.08
        else:
            size = 0.093
    
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

            system3 = ROOT.TLatex(x,ymax-3*dy, '#it{{R}} = {}    |#it{{#eta}}_{{jet}}| < {}'.format(R, 0.9-R))
            system3.SetNDC()
            system3.SetTextSize(size*scale_factor)
            system3.Draw()

            self.plot_list.append(system0)
            self.plot_list.append(system1)
            self.plot_list.append(system2)
            self.plot_list.append(system3)

        if groomed:
            if pad in [1,2]:
                x_pt = x+0.2
                y_pt = ymax-5.5*dy
            else:
                x_pt = x+0.22
                y_pt = ymax-9.*dy
        else:
            x_pt = x+0.2
            y_pt = ymax-4.5*dy
        system4 = ROOT.TLatex(x_pt, y_pt, str(minpt) + ' < #it{p}_{T}^{ch jet} < ' + str(maxpt) + ' GeV/#it{c}')
        system4.SetNDC()
        system4.SetTextSize(size*scale_factor)
        system4.Draw()
        self.plot_list.append(system4)

        if groomed and pad == 1:
            system5 = ROOT.TLatex(x, ymax-4.2*dy, 'Soft Drop #it{z}_{cut} = 0.2 #beta = 0')
            system5.SetNDC()
            system5.SetTextSize(size*scale_factor)
            system5.Draw()
            self.plot_list.append(system5)

        # Set pad for ratio
        c.cd(pad)
        if pad in [3,4]:
            pad2 = ROOT.TPad("pad2_{}".format(R), "pad2{}".format(R),0,self.bottom_offset+self.ratio_height,1,self.bottom_offset+2*self.ratio_height)
            pad3 = ROOT.TPad("pad3_{}".format(R), "pad3{}".format(R),0,0,1,self.bottom_offset+self.ratio_height)
        else:
            pad2 = ROOT.TPad("pad2_{}".format(R), "pad2{}".format(R),0,self.ratio_height,1,2*self.ratio_height)
            pad3 = ROOT.TPad("pad3_{}".format(R), "pad3{}".format(R),0,0,1,self.ratio_height)
        self.plot_list.append(pad2)
        self.plot_list.append(pad3)
        if pad in [1,3]:
            pad2.SetLeftMargin(self.left_offset)
            pad3.SetLeftMargin(self.left_offset)
        else:
            pad2.SetLeftMargin(0.)
            pad3.SetLeftMargin(0.)
        pad2.SetRightMargin(0.)
        pad2.SetTopMargin(0.)
        pad3.SetRightMargin(0.)
        pad3.SetTopMargin(0.)
        if pad in [3,4]:
            pad2.SetBottomMargin(0.)
            pad3.SetBottomMargin(self.bottom_offset/(self.bottom_offset+self.ratio_height))
        else:
            pad2.SetBottomMargin(0.)
            pad3.SetBottomMargin(0.)
        if self.logx:
            pad2.SetLogx()
            pad3.SetLogx()
        pad2.SetTicks(1,1)
        pad2.Draw()
        pad3.SetTicks(1,1)
        pad3.Draw()

        # Draw blank histos
        pad2.cd()
        blankname = 'myBlankHisto2_{}_{}'.format(pad, R)
        myBlankHisto2 = ROOT.TH1F(blankname,blankname, 1, self.xmin, self.xmax)
        myBlankHisto2.SetNdivisions(505, "y")
        if groomed:
            myBlankHisto2.SetXTitle(self.xtitle_groomed)
        else:
            myBlankHisto2.SetXTitle(self.xtitle)
        myBlankHisto2.SetYTitle('#frac{Data}{PYTHIA8}')
        if pad in [1,2]:
            myBlankHisto2.GetYaxis().SetTitleSize(0.2)
            myBlankHisto2.GetYaxis().SetTitleOffset(0.4)
            myBlankHisto2.GetYaxis().SetLabelSize(0.2)
        else:
            myBlankHisto2.GetYaxis().SetTitleSize(0.2)
            myBlankHisto2.GetYaxis().SetTitleOffset(0.4)
            myBlankHisto2.GetYaxis().SetLabelSize(0.2)
        myBlankHisto2.SetMinimum(0.01)
        myBlankHisto2.SetMaximum(1.99)
        myBlankHisto2.Draw()
        self.blank_histo_list.append(myBlankHisto2)
        
        pad3.cd()
        blankname = 'myBlankHisto3_{}_{}'.format(pad, R)
        myBlankHisto3 = myBlankHisto2.Clone()
        myBlankHisto3.SetNdivisions(505, "y")
        if groomed:
            myBlankHisto3.SetXTitle(self.xtitle_groomed)
        else:
            myBlankHisto3.SetXTitle(self.xtitle)
        myBlankHisto3.SetYTitle('#frac{Data}{Herwig7}')
        if pad in [1,2]:
            myBlankHisto3.GetYaxis().SetTitleSize(0.2)
            myBlankHisto3.GetYaxis().SetTitleOffset(0.4)
            myBlankHisto3.GetYaxis().SetLabelSize(0.2)
            myBlankHisto3.SetMinimum(0.01)
        else:
            myBlankHisto3.GetXaxis().SetTitleSize(0.2)
            myBlankHisto3.GetXaxis().SetTitleOffset(0.8)
            myBlankHisto3.GetXaxis().SetLabelSize(0.12)
            myBlankHisto3.GetYaxis().SetTitleSize(0.11)
            myBlankHisto3.GetYaxis().SetTitleOffset(0.7)
            myBlankHisto3.GetYaxis().SetLabelSize(0.11)
            myBlankHisto3.SetMinimum(0.)
        myBlankHisto3.Draw()
        self.blank_histo_list.append(myBlankHisto3)

        line = ROOT.TLine(self.xmin,1,self.xmax,1)
        line.SetLineColor(1)
        line.SetLineStyle(2)
        pad2.cd()
        line.Draw('same')
        pad3.cd()
        line.Draw('same')
        self.plot_list.append(line)

        # Draw ratio
        for i,beta in enumerate(self.beta_list):
        
            pad2.cd()
            self.h_pythia_ratio_sys_list[i].SetLineColor(0)
            self.h_pythia_ratio_sys_list[i].SetFillColor(self.colors[i])
            self.h_pythia_ratio_sys_list[i].SetFillColorAlpha(self.colors[i], 0.3)
            self.h_pythia_ratio_sys_list[i].SetFillStyle(1001)
            self.h_pythia_ratio_sys_list[i].SetLineWidth(0)
            self.h_pythia_ratio_sys_list[i].Draw('E2 same')

            self.h_pythia_ratio_list[i].SetMarkerColorAlpha(self.colors[i], 1)
            self.h_pythia_ratio_list[i].SetLineColorAlpha(self.colors[i], self.alpha)
            self.h_pythia_ratio_list[i].SetLineWidth(2)
            self.h_pythia_ratio_list[i].SetMarkerStyle(self.markers[i])
            self.h_pythia_ratio_list[i].SetMarkerSize(self.marker_size)
            self.h_pythia_ratio_list[i].Draw('PE X0 same')
            
            pad3.cd()
            self.h_herwig_ratio_sys_list[i].SetLineColor(0)
            self.h_herwig_ratio_sys_list[i].SetFillColor(self.colors[i])
            self.h_herwig_ratio_sys_list[i].SetFillColorAlpha(self.colors[i], 0.3)
            self.h_herwig_ratio_sys_list[i].SetFillStyle(1001)
            self.h_herwig_ratio_sys_list[i].SetLineWidth(0)
            self.h_herwig_ratio_sys_list[i].Draw('E2 same')

            self.h_herwig_ratio_list[i].SetMarkerColorAlpha(self.colors[i], 1)
            self.h_herwig_ratio_list[i].SetLineColorAlpha(self.colors[i], self.alpha)
            self.h_herwig_ratio_list[i].SetLineWidth(2)
            self.h_herwig_ratio_list[i].SetMarkerStyle(self.markers[i])
            self.h_herwig_ratio_list[i].SetMarkerSize(self.marker_size)
            self.h_herwig_ratio_list[i].Draw('PE X0 same')

    # Scale vertical amplitude of histogram, for visualization
    #-------------------------------------------------------------------------------------------
    def scale_histogram_for_visualization(self, h, i, minpt, groomed):

        scale_factor = 1.
        if groomed:
            if i == 1:
                if minpt > 50.:
                    scale_factor = self.scale_factor_groomed_beta15_highpt
                else:
                    scale_factor = self.scale_factor_groomed_beta15_lowpt
            if i == 2:
                if minpt > 50.:
                    scale_factor = self.scale_factor_groomed_beta2_highpt
                else:
                    scale_factor = self.scale_factor_groomed_beta2_lowpt
            if i == 3:
                if minpt > 50.:
                    scale_factor = self.scale_factor_groomed_beta3_highpt
                else:
                    scale_factor = self.scale_factor_groomed_beta3_lowpt
        else:
            if i == 1:
                if minpt > 50.:
                    scale_factor = self.scale_factor_ungroomed_beta15_highpt
                else:
                    scale_factor = self.scale_factor_ungroomed_beta15_lowpt
            if i == 2:
                if minpt > 50.:
                    scale_factor = self.scale_factor_ungroomed_beta2_highpt
                else:
                    scale_factor = self.scale_factor_ungroomed_beta2_lowpt
            if i == 3:
                if minpt > 50.:
                    scale_factor = self.scale_factor_ungroomed_beta3_highpt
                else:
                    scale_factor = self.scale_factor_ungroomed_beta3_lowpt

        h.Scale(scale_factor)

        if math.isclose(scale_factor, 1.):
            plot_label = ''
        else:
            plot_label = ' (#times{})'.format(scale_factor)
        
        return plot_label

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
