#! /usr/bin/env python

import sys
import os
import argparse
from array import *
import numpy
import ROOT
import yaml

from pyjetty.alice_analysis.analysis.user.substructure import run_analysis

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class RunAnalysisJames(run_analysis.RunAnalysis):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(RunAnalysisJames, self).__init__(config_file, **kwargs)

    # Initialize yaml config
    self.initialize_user_config()

    print(self)

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_user_config(self):

    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)

    # Theory comparisons
    if 'fPythia' in config:
      self.fPythia_name = config['fPythia']
    if 'fNLL' in config:
      self.fNLL = config['fNLL']
    if 'fNPcorrection_numerator' in config and 'fNPcorrection_denominator' in config:
      self.fNPcorrection_numerator = config['fNPcorrection_numerator']
      self.fNPcorrection_denominator = config['fNPcorrection_denominator']

  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #---------------------------------------------------------------
  def plot_single_result(self, jetR, obs_label, obs_setting, sd_setting):

    self.plot_final_result(jetR, obs_label, obs_setting, sd_setting)

  #---------------------------------------------------------------
  # This function is called once after all subconfigurations have been looped over, for each R
  #---------------------------------------------------------------
  def plot_all_results(self, jetR, obs_label, obs_setting, sd_setting):
    #TODO
    pass

  #----------------------------------------------------------------------
  # This function is called once after all subconfigurations and jetR have been looped over
  #----------------------------------------------------------------------
  def plot_performance(self):

    print('Plotting performance plots...')

  #----------------------------------------------------------------------
  def plot_final_result(self, jetR, obs_label, obs_setting, sd_setting):
    print('Plot final results for {}: R = {}, {} ...'.format(self.observable, jetR, obs_label))

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    # Loop through pt slices, and compute systematics for each 1D theta_g distribution
    for bin in range(0, len(self.pt_bins_reported) - 1):
      min_pt_truth = self.pt_bins_reported[bin]
      max_pt_truth = self.pt_bins_reported[bin+1]

      #self.get_NPcorrection(self.observable, jetR, obs_label, obs_setting, min_pt_truth, max_pt_truth)

      self.plot_observable(jetR, obs_label, obs_setting, sd_setting, min_pt_truth, max_pt_truth, plot_pythia=False)
      #self.plot_observable(jetR, obs_label, obs_setting, sd_setting, min_pt_truth, max_pt_truth, plot_pythia=True)

  #----------------------------------------------------------------------
  def get_NPcorrection(self, observable, jetR, obs_label, obs_setting, min_pt_truth, max_pt_truth):

    fNumerator = ROOT.TFile(self.fNPcorrection_numerator, 'READ')
    fDenominator = ROOT.TFile(self.fNPcorrection_denominator, 'READ')

    n_bins_truth = self.n_bins_truth(obs_label)
    truth_bin_array = self.truth_bin_array(obs_label)

    hname = 'histogram_h_{}_B{}_{}-{}'.format(observable, obs_setting[1], int(min_pt_truth), int(max_pt_truth))

    hNumerator = fNumerator.Get(hname)
    hNumerator.SetDirectory(0)
    n_jets_inclusive = hNumerator.Integral(0, hNumerator.GetNbinsX()+1)
    n_jets_tagged = hNumerator.Integral(hNumerator.FindBin(truth_bin_array[0]), hNumerator.GetNbinsX())
    fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
    hNumerator.Scale(1./n_jets_inclusive, 'width')

    hDenominator = fDenominator.Get(hname)
    hDenominator.SetDirectory(0)
    n_jets_inclusive = hDenominator.Integral(0, hDenominator.GetNbinsX()+1)
    n_jets_tagged = hDenominator.Integral(hDenominator.FindBin(truth_bin_array[0]), hDenominator.GetNbinsX())
    fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
    hDenominator.Scale(1./n_jets_inclusive, 'width')

    hNumerator.Divide(hDenominator)
    hNumerator.SetName('hNPcorrection_{}_{}_{}-{}'.format(observable, obs_label, min_pt_truth, max_pt_truth))
    setattr(self, 'hNPcorrection_{}_{}_{}-{}'.format(observable, obs_label, min_pt_truth, max_pt_truth), hNumerator)

    #output_dir = getattr(self, 'output_dir_final_results')
    #outputFilename = os.path.join(output_dir, 'hNPcorrection_{}_{}_{}-{}.pdf'.format(observable, obs_label, min_pt_truth, max_pt_truth))
    #self.utils.plot_hist(hNumerator, outputFilename)
  
  #----------------------------------------------------------------------
  def plot_observable(self, jetR, obs_label, obs_setting, sd_setting, min_pt_truth, max_pt_truth, plot_pythia=False):

    name = 'cResult_R{}_{}_{}-{}'.format(jetR, obs_label, min_pt_truth, max_pt_truth)
    c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()

    c.cd()
    myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    myPad.SetLeftMargin(0.2)
    myPad.SetTopMargin(0.07)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.13)
    myPad.Draw()
    myPad.cd()

    ymax = getattr(self, 'ymax')
    xtitle = getattr(self, 'xtitle')
    ytitle = getattr(self, 'ytitle')

    n_obs_bins_truth = self.n_bins_truth(obs_label)
    truth_bin_array = self.truth_bin_array(obs_label)
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_obs_bins_truth, truth_bin_array)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle( xtitle )
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle(ytitle)
    myBlankHisto.SetMaximum(ymax)
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")

    if plot_pythia:

      hPythia = self.get_pythia_from_response(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
      n_jets_inclusive = hPythia.Integral(0, hPythia.GetNbinsX()+1)
      n_jets_tagged = hPythia.Integral(hPythia.FindBin(truth_bin_array[0]), hPythia.GetNbinsX())
      fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
      hPythia.Scale(1./n_jets_inclusive, 'width')
      hPythia.SetFillStyle(0)
      hPythia.SetMarkerSize(1.5)
      hPythia.SetMarkerStyle(21)
      hPythia.SetMarkerColor(2)
      hPythia.SetLineColor(2)
      hPythia.SetLineWidth(1)
      hPythia.Draw('E2 same')

    color = 600-6
    h_sys = getattr(self, 'hResult_{}_systotal_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth))
    h_sys.SetLineColor(0)
    h_sys.SetFillColor(color)
    h_sys.SetFillColorAlpha(color, 0.3)
    h_sys.SetFillStyle(1001)
    h_sys.SetLineWidth(0)
    h_sys.DrawCopy('E2 same')

    name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, jetR, obs_label, min_pt_truth, max_pt_truth)
    if sd_setting:
      fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
    h = getattr(self, name)
    h.SetMarkerSize(1.5)
    h.SetMarkerStyle(20)
    h.SetMarkerColor(color)
    h.SetLineStyle(1)
    h.SetLineWidth(2)
    h.SetLineColor(color)
    h.DrawCopy('PE X0 same')

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE Preliminary'
    text_latex.DrawLatex(0.57, 0.87, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.8, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.73, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = '#it{R} = ' + str(jetR) + '   | #eta_{jet}| < 0.5'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.66, text)

    subobs_label = self.utils.formatted_subobs_label(self.observable)
    if subobs_label:
      text = '{} = {}'.format(subobs_label, obs_setting)
      text_latex.DrawLatex(0.57, 0.59, text)

    if sd_setting:
      text = self.utils.formatted_sd_label(sd_setting)
      text_latex.DrawLatex(0.57, 0.52, text)

      text_latex = ROOT.TLatex()
      text_latex.SetNDC()
      text_latex.SetTextSize(0.04)
      text = '#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged
      text_latex.DrawLatex(0.57, 0.45, text)

      if plot_pythia:
        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        text_latex.SetTextSize(0.04)
        text = ('#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged) + (', #it{f}_{tagged}^{pythia} = %3.3f' % fraction_tagged_pythia)
        text_latex.DrawLatex(0.57, 0.38, text)

    myLegend = ROOT.TLegend(0.27,0.7,0.5,0.85)
    self.utils.setup_legend(myLegend,0.035)
    myLegend.AddEntry(h, 'ALICE pp', 'pe')
    myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    if plot_pythia:
      myLegend.AddEntry(hPythia, 'PYTHIA Monash2013', 'pe')
    myLegend.Draw()

    name = 'hUnfolded_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), obs_label, int(min_pt_truth), int(max_pt_truth), self.file_format)
    if plot_pythia:
      name = 'hUnfolded_R{}_{}_{}-{}_Pythia{}'.format(self.utils.remove_periods(jetR), obs_label, int(min_pt_truth), int(max_pt_truth), self.file_format)
    output_dir = getattr(self, 'output_dir_final_results')
    outputFilename = os.path.join(output_dir, name)
    c.SaveAs(outputFilename)
    c.Close()

    # Write result to ROOT file
    if not plot_pythia:
      final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
      fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
      h.Write()
      fFinalResults.Close()

  #----------------------------------------------------------------------
  def get_pythia_from_response(self, observable, jetR, obs_label, min_pt_truth, max_pt_truth):

    output_dir = getattr(self, 'output_dir_main')
    file = os.path.join(output_dir, 'response.root')
    f = ROOT.TFile(file, 'READ')

    if observable == 'theta_g':
      thn_name = 'hResponse_JetPt_ThetaG_R{}_{}_rebinned'.format(jetR, obs_label)
    if observable == 'zg':
      thn_name = 'hResponse_JetPt_zg_R{}_{}_rebinned'.format(jetR, obs_label)

    thn = f.Get(thn_name)
    thn.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)

    h = thn.Projection(3)
    h.SetName('hPythia_{}_R{}_{}_{}-{}'.format(observable, jetR, obs_label, min_pt_truth, max_pt_truth))
    h.SetDirectory(0)

    return h

#----------------------------------------------------------------------
if __name__ == '__main__':

  # Define arguments
  parser = argparse.ArgumentParser(description='Jet substructure analysis')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='analysis_config.yaml',
                      help='Path of config file for analysis')

  # Parse the arguments
  args = parser.parse_args()

  print('Configuring...')
  print('configFile: \'{0}\''.format(args.configFile))

  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = RunAnalysisJames(config_file = args.configFile)
  analysis.run_analysis()
