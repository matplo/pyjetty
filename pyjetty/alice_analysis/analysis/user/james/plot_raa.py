#! /usr/bin/env python

import sys
import os
import argparse
from array import *
import numpy
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
  def __init__(self, config_file='', **kwargs):
    super(PlotRAA, self).__init__(**kwargs)
    
    self.utils = analysis_utils_obs.AnalysisUtils_Obs()
    
    self.initialize()

    print(self)
  
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize(self):
    
    # Read config file
    #with open(self.config_file, 'r') as stream:
    #  config = yaml.safe_load(stream)
      
    #self.figure_approval_status = config['figure_approval_status']
    self.figure_approval_status = 'Work in progress'
    
    self.output_dir = '/Users/jamesmulligan/Analysis_theta_g/roounfold_output/'

    self.observables = ['theta_g']
    self.obs_label = 'SD_zcut02_B0'
    self.formatted_grooming_label = 'SD #it{z}_{cut}=0.2, #beta=0'
    self.jetR = 0.2
    if self.jetR == 0.4:
      self.filedir_pp_name = '/Users/jamesmulligan/Analysis_theta_g/roounfold_output/pp_ref/67940-04/'
      self.filedir_AA_name = '/Users/jamesmulligan/Analysis_theta_g/roounfold_output/PbPb/78488-04/'
      self.min_pt = 80
      self.max_pt = 100
    else:
      self.filedir_pp_name = '/Users/jamesmulligan/Analysis_theta_g/roounfold_output/pp_ref/67940-02/'
      self.filedir_AA_name = '/Users/jamesmulligan/Analysis_theta_g/roounfold_output/PbPb/78488-02/'
      self.min_pt = 60
      self.max_pt = 80
    
    if not os.path.exists(self.output_dir):
      os.makedirs(self.output_dir)

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def init_observable(self, observable):
  
    self.file_pp_name = os.path.join(self.filedir_pp_name, '{}/final_results/fFinalResults.root'.format(observable))
    self.file_AA_name = os.path.join(self.filedir_AA_name, '{}/final_results/fFinalResults.root'.format(observable))

    if observable == 'theta_g':
      self.xtitle = '#theta_{g}'
    if observable == 'zg':
      self.xtitle = '#it{z}_{g}'
    self.ytitle = '#frac{{1}}{{#sigma_{{jet, inc}}}} #frac{{d#sigma}}{{d{}}}'.format(self.xtitle)
 
    self.main_result_name = 'hmain_{}_R{}_{}_{}-{}'.format(observable, self.jetR, self.obs_label, self.min_pt, self.max_pt)
    self.sys_total_name = 'hResult_{}_systotal_R{}_{}_{}-{}'.format(observable, self.jetR, self.obs_label, self.min_pt, self.max_pt)
    self.h_tagging_fraction_name = 'h_tagging_fraction_R{}_{}'.format(self.jetR, self.obs_label)
 
    self.file_pp = ROOT.TFile(self.file_pp_name, 'READ')
    self.file_AA = ROOT.TFile(self.file_AA_name, 'READ')

    self.h_main_pp = self.file_pp.Get(self.main_result_name)
    self.h_sys_pp = self.file_pp.Get(self.sys_total_name)
    self.h_tagging_fraction_pp = self.file_pp.Get(self.h_tagging_fraction_name)
    self.h_main_AA = self.file_AA.Get(self.main_result_name)
    self.h_sys_AA = self.file_AA.Get(self.sys_total_name)
    self.h_tagging_fraction_AA = self.file_AA.Get(self.h_tagging_fraction_name)
    
    pt = (self.min_pt + self.max_pt)/2
    bin = self.h_tagging_fraction_pp.GetXaxis().FindBin(pt)
    self.f_tagging_pp = self.h_tagging_fraction_pp.GetBinContent(bin)
    pt = (self.min_pt + self.max_pt)/2
    bin = self.h_tagging_fraction_AA.GetXaxis().FindBin(pt)
    self.f_tagging_AA = self.h_tagging_fraction_AA.GetBinContent(bin)

    self.ymax = self.h_main_pp.GetMaximum()
    self.colors = [600-6, 632-4, 416-2]
    self.markers = [20, 21, 33]
    
    name = 'hRAA_{}_R{}.pdf'.format(observable, self.utils.remove_periods(self.jetR))
    self.output_filename = os.path.join(self.output_dir, name)
      
  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #----------------------------------------------------------------------
  def plot_raa(self):
  
    for observable in self.observables:
  
      self.init_observable(observable)
      self.plot_final_result(observable)
  
  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #----------------------------------------------------------------------
  def plot_final_result(self, observable):
    print('Plot final results for {}: R = {}, {}'.format(observable, self.jetR, self.obs_label))

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    c = ROOT.TCanvas('c', 'c', 600, 650)
    c.Draw()
    
    c.cd()
    pad1 = ROOT.TPad('myPad', 'The pad',0,0.3,1,1)
    pad1.SetLeftMargin(0.2)
    pad1.SetTopMargin(0.07)
    pad1.SetRightMargin(0.04)
    pad1.SetBottomMargin(0.)
    pad1.Draw()
    pad1.cd()

    myLegend = ROOT.TLegend(0.66,0.65,0.8,0.85)
    self.utils.setup_legend(myLegend,0.04)
    
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
      
    xmin = self.h_sys_pp.GetBinLowEdge(2)
    xmax = self.h_sys_pp.GetBinLowEdge(self.h_sys_pp.GetNbinsX()+1)
        
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, xmin, xmax)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle(self.xtitle)
    myBlankHisto.SetYTitle(self.ytitle)
    if observable == 'theta_g':
      ymax = 3*self.ymax
    else:
      ymax = 2.3*self.ymax
    myBlankHisto.SetMaximum(ymax)
    myBlankHisto.SetMinimum(2e-4) # Don't draw 0 on top panel
    myBlankHisto.GetYaxis().SetTitleSize(0.065)
    myBlankHisto.GetYaxis().SetTitleOffset(1.4)
    myBlankHisto.GetYaxis().SetLabelSize(0.06)
    myBlankHisto.Draw('E')
    
    if observable == 'theta_g':
      rg_axis_tf1 = ROOT.TF1('rg_axis_tf1', 'x', 0, self.jetR-0.01)
      rg_axis = ROOT.TGaxis(xmin, ymax, xmax, ymax, 'rg_axis_tf1', 505, '- S')
      rg_axis.SetTitle('#it{R}_{g}')
      rg_axis.SetTitleSize(25)
      rg_axis.SetTitleFont(43)
      rg_axis.SetTitleOffset(0.6)
      rg_axis.SetLabelFont(43)
      rg_axis.SetLabelSize(25)
      rg_axis.SetTickSize(0.015)
      rg_axis.SetLabelOffset(0.015)
      rg_axis.Draw()
 
    c.cd()
    pad2 = ROOT.TPad('pad2', 'pad2', 0, 0.02, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.4)
    pad2.SetLeftMargin(0.2)
    pad2.SetRightMargin(0.04)
    pad2.Draw()
    pad2.cd()
          
    myBlankHisto2 = myBlankHisto.Clone('myBlankHisto_C')
    myBlankHisto2.SetYTitle('#frac{Pb-Pb}{pp}')
    myBlankHisto2.SetXTitle(self.xtitle)
    myBlankHisto2.GetXaxis().SetTitleSize(30)
    myBlankHisto2.GetXaxis().SetTitleFont(43)
    myBlankHisto2.GetXaxis().SetTitleOffset(4.)
    myBlankHisto2.GetXaxis().SetLabelFont(43)
    myBlankHisto2.GetXaxis().SetLabelSize(25)
    myBlankHisto2.GetYaxis().SetTitleSize(20)
    myBlankHisto2.GetYaxis().SetTitleFont(43)
    myBlankHisto2.GetYaxis().SetTitleOffset(2.2)
    myBlankHisto2.GetYaxis().SetLabelFont(43)
    myBlankHisto2.GetYaxis().SetLabelSize(25)
    myBlankHisto2.GetYaxis().SetNdivisions(505)
    if observable == 'theta_g':
      myBlankHisto2.GetYaxis().SetRangeUser(0., 2.49)
    else:
      myBlankHisto2.GetYaxis().SetRangeUser(0., 1.99)

    myBlankHisto2.Draw()
  
    line = ROOT.TLine(xmin,1,xmax,1)
    line.SetLineColor(920+2)
    line.SetLineStyle(2)
    line.Draw()
      
    h_ratio = self.h_main_AA.Clone()
    h_ratio.Divide(self.h_main_pp)
      
    h_ratio_sys = self.h_sys_AA.Clone()
    h_ratio_sys.Divide(self.h_sys_pp)

    pad1.cd()
    self.h_sys_pp.DrawCopy('E2 same')
    self.h_sys_AA.DrawCopy('E2 same')
    self.h_main_pp.DrawCopy('PE X0 same')
    self.h_main_AA.DrawCopy('PE X0 same')
      
    pad2.cd()
    h_ratio_sys.DrawCopy('E2 same')
    h_ratio.DrawCopy('PE X0 same')
          
    pad1.cd()
    myLegend.AddEntry(self.h_main_pp, 'pp', 'PE')
    myLegend.AddEntry(self.h_main_AA, 'Pb-Pb 0-10%', 'PE')
    myLegend.AddEntry(self.h_sys_pp, 'Sys. uncertainty', 'f')
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    
    text_latex.SetTextSize(0.04)
    text = '#it{{f}}_{{tagged}}^{{pp}} = {:.2f}, #it{{f}}_{{tagged}}^{{AA}} = {:.2f}'.format(self.f_tagging_pp, self.f_tagging_AA)
    text_latex.DrawLatex(0.57, 0.57, text)

    text_latex.SetTextSize(0.045)
    text = 'ALICE {}'.format(self.figure_approval_status)
    text_latex.DrawLatex(0.25, 0.87, text)
    
    text = '#sqrt{#it{s_{#it{NN}}}} = 5.02 TeV'
    text_latex.DrawLatex(0.25, 0.81, text)

    text = 'Charged jets   anti-#it{k}_{T}'
    text_latex.DrawLatex(0.25, 0.75, text)
    
    text = '#it{R} = ' + str(self.jetR) + '   | #eta_{jet}| < 0.5'
    text_latex.DrawLatex(0.25, 0.69, text)
    
    text = str(self.min_pt) + ' < #it{p}_{T, ch jet} < ' + str(self.max_pt) + ' GeV/#it{c}'
    text_latex.DrawLatex(0.25, 0.63, text)
    
    text = self.formatted_grooming_label
    text_latex.DrawLatex(0.25, 0.57, text)
    
    myLegend.Draw()

    c.SaveAs(self.output_filename)
    c.Close()
    
    # Write result to ROOT file
    #final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
    #fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
    #h.Write()
    #h_sys.Write()
    #hPythia.Write()
    #fFinalResults.Close()

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
  
  analysis = PlotRAA(config_file = args.configFile)
  analysis.plot_raa()
