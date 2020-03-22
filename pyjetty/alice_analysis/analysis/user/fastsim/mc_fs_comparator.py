#!/usr/bin/env python3

'''
For comparing fast simulation results to the full simulation.
Ezra Lesser (elesser@berkeley.edu)
'''

from __future__ import print_function

import os
import sys
import argparse
import yaml
from array import *
import ROOT

from pyjetty.alice_analysis.analysis.base import analysis_base, analysis_utils


################################################################
class mc_fs_comparator(analysis_base.AnalysisBase):

    ################################################################
    # Analysis class constructor
    def __init__(self, inputFileMC, inputFileFS, configFile, outputDir, imageFormat, **kwargs):
        super(mc_fs_comparator, self).__init__(inputFileFS, inputFileMC, configFile, 
                                               outputDir, imageFormat, **kwargs)

        # Load configuration and open files
        with open(args.configFile, 'r') as stream:
            self.config = yaml.safe_load(stream)

        self.jetR_list = self.config['jetR']
        self.beta_list = self.config['betas']

        self.f_mc = ROOT.TFile(args.inputFileMC, 'READ')
        self.f_fs = ROOT.TFile(args.inputFileFS, 'READ')

        if not args.outputDir[-1] == '/':
            args.outputDir += '/'
        f_out_name = args.outputDir + "AnalysisResultsFStoMC.root"
        self.f_out = ROOT.TFile(f_out_name, 'RECREATE')

    ################################################################
    # Create comparison plots
    def mc_fs_comparator(self):

        for jetR in self.jetR_list:
            for beta in self.beta_list:
                label = ("R%s_B%s" % (str(jetR), str(beta))).replace('.', '')

                # Get relevant info from configuration file
                pt_bins = self.config["configs"][label]["pt_bins_det"]
                ang_bins = self.config["configs"][label]["lambda_bins_det"]

                # Rebin histograms for comparions
                name = "hLambda_JetpT_det_%sScaled" % label
                h_mc = self.f_mc.Get(name)
                print("self.utils", self.utils)
                h_mc_rebinned = self.utils.rebin_data(h_mc, h_mc.GetName()+'_mc', len(pt_bins)-1, 
                                                      array('d', pt_bins), len(ang_bins)-1,
                                                      array('d', ang_bins))
                h_fs = self.f_fs.Get(name)
                h_fs_rebinned = self.utils.rebin_data(h_fs, h_fs.GetName()+'_fs', len(pt_bins)-1,
                                                      array('d', pt_bins), len(ang_bins)-1,
                                                      array('d', ang_bins))
                h_tru = self.f_mc.Get("hLambda_JetpT_tru_%sScaled" % label)
                h_tru_rebinned = self.utils.rebin_data(h_tru, h_tru.GetName()+'_fs', len(pt_bins)-1,
                                                       array('d', pt_bins), len(ang_bins)-1,
                                                       array('d', ang_bins))

                for i in range(0, len(pt_bins)-1):
                    title = "FStoMC_%s_pT_%i-%i" % (label, pt_bins[i], pt_bins[i+1])

                    # Canvas to accompany all plots
                    c = ROOT.TCanvas('c'+title, 'c'+title, 900, 1100)
                    c.cd()

                    # "Pad" to contain the top histograms
                    topPad = ROOT.TPad('topPad', 'topPad', 0, 0.3, 1, 1.0)
                    topPad.SetLeftMargin(0.15)
                    topPad.SetTopMargin(0.1)
                    topPad.SetRightMargin(0.05)
                    topPad.SetBottomMargin(0)
                    topPad.Draw()
                    topPad.cd()

                    # Draw histograms to top pad
                    h_tru_projy = h_tru_rebinned.ProjectionY(title+'_tru', i, i+1)
                    h_tru_projy.SetMarkerSize(1.5)
                    h_tru_projy.SetMarkerStyle(21)
                    h_tru_projy.SetMarkerColor(1)
                    h_tru_projy.SetStats(0)
                    h_tru_projy.Scale(1. / h_tru_projy.GetEntries(), "width")
                    h_tru_projy.SetTitle("Angularity distribution for full and fast MC")
                    h_tru_projy.SetYTitle('#frac{1}{N} #frac{dN}{d#lambda}')
                    h_tru_projy.Draw("hist same E")
                    h_mc_projy = h_mc_rebinned.ProjectionY(title+'_mc', i, i+1)
                    h_mc_projy.SetMarkerSize(1.5)
                    h_mc_projy.SetMarkerStyle(20)
                    h_mc_projy.SetMarkerColor(2)
                    h_mc_projy.SetStats(0)
                    h_mc_projy.Scale(1. / h_mc_projy.GetEntries(), "width")
                    h_mc_projy.Draw("hist same E")
                    h_fs_projy = h_fs_rebinned.ProjectionY(title+'_fs', i, i+1)
                    h_fs_projy.SetMarkerSize(1.5)
                    h_fs_projy.SetMarkerStyle(21)
                    h_fs_projy.SetMarkerColor(4)
                    h_fs_projy.SetStats(0)
                    h_fs_projy.Scale(1. / h_fs_projy.GetEntries(), "width")
                    h_fs_projy.Draw("hist same E")

                    # Create and format the legend
                    leg = ROOT.TLegend(0.75, 0.8, 0.95, 0.9)
                    leg.SetHeader(title, 'C')
                    leg.SetTextSize(0.)
                    leg.AddEntry(h_tru_projy, "Full MC, truth", "lep")
                    leg.AddEntry(h_mc_projy, "Full MC, det level", "lep")
                    leg.AddEntry(h_fs_projy, "Fast simulation", "lep")
                    leg.Draw()

                    # Text overlaid on the plot
                    text_latex = ROOT.TLatex()
                    text_latex.SetNDC()
                    text = "%i < p_{T} < %s GeV/c" % (pt_bins[i], pt_bins[i+1])
                    text_latex.DrawLatex(0.2, 0.85, text)
                    text = "R = " + str(jetR)
                    text_latex.DrawLatex(0.2, 0.75, text)

                    # Create bottom pad for ratio plot
                    c.cd()
                    botPad = ROOT.TPad("botPad", "botPad", 0, 0.05, 1, 0.3)
                    botPad.SetTopMargin(0)
                    botPad.SetBottomMargin(0.25)
                    botPad.SetLeftMargin(0.15)
                    botPad.SetRightMargin(0.05)
                    botPad.Draw()
                    botPad.cd()

                    # Create and format the ratio plot
                    h_ratio = h_fs_projy.Clone()
                    h_ratio.Divide(h_mc_projy)
                    h_ratio.SetTitle("")#'h'+title)
                    h_ratio.GetXaxis().SetTitleSize(30)
                    h_ratio.GetXaxis().SetTitleFont(43)
                    h_ratio.GetXaxis().SetTitleOffset(4.)
                    h_ratio.GetXaxis().SetLabelFont(43)
                    h_ratio.GetXaxis().SetLabelSize(20)
                    h_ratio.SetXTitle('#lambda_{#beta=%s}' % beta)
                    h_ratio.GetYaxis().SetTitleSize(20)
                    h_ratio.GetYaxis().SetTitleFont(43)
                    h_ratio.GetYaxis().SetTitleOffset(2.2)
                    h_ratio.GetYaxis().SetLabelFont(43)
                    h_ratio.GetYaxis().SetLabelSize(20)
                    h_ratio.SetYTitle("fast sim / full MC det")
                    h_ratio.SetMarkerColor(51)
                    h_ratio.SetStats(0)
                    h_ratio.SetMinimum(0.8)
                    h_ratio.SetMaximum(1.2)
                    h_ratio.Draw('E')

                    # Line at 1 for visual reference
                    line = ROOT.TLine(0, 1, ang_bins[-1], 1)
                    line.Draw()

                    input("Press Enter to continue...")

        # Close files after reading / writing
        self.close_files()

    ################################################################
    def close_files(self):
        self.f_mc.Close()
        self.f_fs.Close()
        self.f_out.Close()


################################################################
if __name__ == '__main__':
    # Parse and check arguments
    parser = argparse.ArgumentParser(description='Compare full and fast simulations')
    parser.add_argument('-m', '--inputFileMC', action='store',
                        type=str, metavar='inputFileMC', default='AnalysisResults.root',
                        help='Path of ROOT file containing full simulation TTrees')
    parser.add_argument('-f', '--inputFileFS', action='store',
                        type=str, metavar='inputFileFS', default='AnalysisResults.root',
                        help='Path of ROOT file containing fast simulation TTrees')
    parser.add_argument('-c', '--configFile', action='store', type=str, metavar='configFile',
                        default='config/analysis_config.yaml',
                        help="Path of config file for analysis")
    parser.add_argument('-o', '--outputDir', action='store',
                        type=str, metavar='outputDir', default='./TestOutput',
                        help='Output directory for output to be written to')
    parser.add_argument('-i', '--imageFormat', action='store',
                        type=str, metavar='imageFormat', default='.png',
                        help='Image format to save plots in, e.g. \".pdf\" or \".png\"')

    # Parse the arguments
    args = parser.parse_args()

    print('Configuring...')
    print('inputFileMC: "{0}"'.format(args.inputFileMC))
    print('inputFileFS: "{0}"'.format(args.inputFileFS))
    print('configFile: "{0}"'.format(args.configFile))
    print('ouputDir: "{0}"'.format(args.outputDir))

    # If invalid inputFile is given, exit
    if not os.path.exists(args.inputFileMC):
        print('File \"{0}\" does not exist! Exiting!'.format(args.inputFileMC))
        sys.exit(0)
    elif not os.path.exists(args.inputFileFS):
        print('File \"{0}\" does not exist! Exiting!'.format(args.inputFileFS))
        sys.exit(0)

    # Initiate class and start analysis
    analysis = mc_fs_comparator(args.inputFileMC, args.inputFileFS, args.configFile,
                                args.outputDir, args.imageFormat) 
    analysis.mc_fs_comparator()
