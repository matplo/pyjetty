#! /usr/bin/env python

import sys
import os
from array import *
import numpy
import ROOT
import yaml

from pyjetty.alice_analysis.process.base import common_base
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class PlotLaura(common_base.CommonBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(PlotLaura, self).__init__(**kwargs)
    
    self.utils = analysis_utils_obs.AnalysisUtils_Obs()
    
    self.initialize()

    print(self)
  
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize(self):
    
    self.output_dir = '/Users/jamesmulligan/Analysis_theta_g/roounfold_output/'
 
    self.observable = 'zg'
    self.obs_label = 'SD_zcut02_B0'
    self.formatted_grooming_label = 'SD #it{z}_{cut}=0.2, #beta=0'
    self.jetR = 0.4
    self.min_pt = 60
    self.max_pt = 80
      
    # James
    self.filedir_james_name = '/Users/jamesmulligan/Analysis_theta_g/roounfold_output/pp_ref/67940_dg/'
    self.file_james_name = os.path.join(self.filedir_james_name, '{}/final_results/fFinalResults.root'.format(self.observable))
   
    self.main_result_name_james = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, self.jetR, self.obs_label, self.min_pt, self.max_pt)
    self.sys_total_name_james = 'hResult_{}_systotal_R{}_{}_{}-{}'.format(self.observable, self.jetR, self.obs_label, self.min_pt, self.max_pt)
    self.h_tagging_fraction_name = 'h_tagging_fraction_R{}_{}'.format(self.jetR, self.obs_label)

    # Laura
    self.filedir_laura_name = '/Users/jamesmulligan/Analysis_theta_g/roounfold_output/laura'
    if self.observable == 'zg':
      self.file_laura_name = os.path.join(self.filedir_laura_name, 'PPzg_LandL.root')
    elif self.observable == 'theta_g':
      self.file_laura_name = os.path.join(self.filedir_laura_name, 'PPrg_LandL.root')
      
    self.main_result_name_laura = 'h1_pp'
    self.sys_total_name_laura = 'h1_pp_err'
    
    if not os.path.exists(self.output_dir):
      os.makedirs(self.output_dir)
      
    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def init_observable(self):
  
    if self.observable == 'theta_g':
      self.xtitle = '#theta_{g}'
    if self.observable == 'zg':
      self.xtitle = '#it{z}_{g}'
    self.ytitle = '#frac{{1}}{{#sigma_{{jet, inc}}}} #frac{{d#sigma}}{{d{}}}'.format(self.xtitle)
 
    self.file_james = ROOT.TFile(self.file_james_name, 'READ')
    self.file_laura = ROOT.TFile(self.file_laura_name, 'READ')

    self.h_main_james = self.file_james.Get(self.main_result_name_james)
    self.h_sys_james = self.file_james.Get(self.sys_total_name_james)
    self.h_tagging_fraction_james = self.file_james.Get(self.h_tagging_fraction_name)

    self.h_main_laura = self.file_laura.Get(self.main_result_name_laura)
    self.h_sys_laura = self.file_laura.Get(self.sys_total_name_laura)
    
    # Translate Laura's Rg to theta_g
    if self.observable == 'theta_g':
      rg_bins = self.h_main_laura.GetXaxis().GetXbins()
      rg_bins = [bin/self.jetR for bin in rg_bins]
      h_new = ROOT.TH1F('hnew', 'hnew', len(rg_bins), array('d',(rg_bins)))
      for i in range(1, len(rg_bins)+1):
        rg = rg_bins[i-1]
        content = self.h_main_laura.GetBinContent(i)
        uncertainty =  self.h_main_laura.GetBinError(i)
        h_new.SetBinContent(i, content)
        h_new.SetBinError(i, uncertainty)
      h_new.Scale(self.jetR)
      self.h_main_laura = h_new

    pt = (self.min_pt + self.max_pt)/2
    #bin = self.h_tagging_fraction_james.GetXaxis().FindBin(pt)
    #self.f_tagging_james = self.h_tagging_fraction_james.GetBinContent(bin)

    self.ymax = self.h_main_james.GetMaximum()
    self.colors = [600-6, 632-4, 416-2]
    self.markers = [20, 21, 33]
    
    name = 'h_crosscheck_pp_{}.pdf'.format(self.observable)
    self.output_filename = os.path.join(self.output_dir, name)
      
  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #----------------------------------------------------------------------
  def plot_crosscheck(self):
  
    self.init_observable()
    self.plot_final_result()
  
  #---------------------------------------------------------------
  # This function is called once for each subconfiguration
  #----------------------------------------------------------------------
  def plot_final_result(self):
    print('Plot final results for {}: R = {}, {}'.format(self.observable, self.jetR, self.obs_label))

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
    
    self.h_main_james.SetMarkerSize(1.5)
    self.h_main_james.SetMarkerStyle(self.markers[0])
    self.h_main_james.SetMarkerColor(self.colors[0])
    self.h_main_james.SetLineStyle(1)
    self.h_main_james.SetLineWidth(2)
    self.h_main_james.SetLineColor(self.colors[0])
      
    self.h_sys_james.SetLineColor(0)
    self.h_sys_james.SetFillColor(self.colors[0])
    self.h_sys_james.SetFillColorAlpha(self.colors[0], 0.3)
    self.h_sys_james.SetFillStyle(1001)
    self.h_sys_james.SetLineWidth(0)
    
    self.h_main_laura.SetMarkerSize(1.5)
    self.h_main_laura.SetMarkerStyle(self.markers[1])
    self.h_main_laura.SetMarkerColor(self.colors[1])
    self.h_main_laura.SetLineStyle(1)
    self.h_main_laura.SetLineWidth(2)
    self.h_main_laura.SetLineColor(self.colors[1])
      
    self.h_sys_laura.SetLineColor(0)
    self.h_sys_laura.SetFillColor(self.colors[1])
    self.h_sys_laura.SetFillColorAlpha(self.colors[1], 0.3)
    self.h_sys_laura.SetFillStyle(1001)
    self.h_sys_laura.SetLineWidth(0)
 
    xmin = self.h_sys_james.GetBinLowEdge(2)
    xmax = self.h_sys_james.GetBinLowEdge(self.h_sys_james.GetNbinsX()+1)
        
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, xmin, xmax)
    myBlankHisto.SetNdivisions(505)
    myBlankHisto.SetXTitle(self.xtitle)
    myBlankHisto.SetYTitle(self.ytitle)
    if self.observable == 'theta_g':
      ymax = 3*self.ymax
    else:
      ymax = 2.3*self.ymax
    myBlankHisto.SetMaximum(ymax)
    myBlankHisto.SetMinimum(2e-4) # Don't draw 0 on top panel
    myBlankHisto.GetYaxis().SetTitleSize(0.065)
    myBlankHisto.GetYaxis().SetTitleOffset(1.4)
    myBlankHisto.GetYaxis().SetLabelSize(0.06)
    myBlankHisto.Draw('E')
    
    if self.observable == 'theta_g':
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
    myBlankHisto2.SetYTitle('')
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
    if self.observable == 'theta_g':
      myBlankHisto2.GetYaxis().SetRangeUser(0., 2.49)
    else:
      myBlankHisto2.GetYaxis().SetRangeUser(0., 1.99)

    myBlankHisto2.Draw()
  
    line = ROOT.TLine(xmin,1,xmax,1)
    line.SetLineColor(920+2)
    line.SetLineStyle(2)
    line.Draw()
      
    h_ratio = self.h_main_laura.Clone()
    h_ratio.Divide(self.h_main_james)
   
    #h_ratio_sys = self.h_sys_laura.Clone()
    #h_ratio_sys.Divide(self.h_sys_james)

    pad1.cd()
    self.h_sys_james.DrawCopy('E2 same')
    #self.h_sys_laura.DrawCopy('E2 same')
    self.h_main_james.DrawCopy('PE X0 same')
    self.h_main_laura.DrawCopy('PE X0 same')
   
    pad1.cd()
    myLegend.AddEntry(self.h_main_james, 'james', 'PE')
    myLegend.AddEntry(self.h_main_laura, 'laura', 'PE')
    myLegend.AddEntry(self.h_sys_james, 'Sys. uncertainty', 'f')
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    
    text_latex.SetTextSize(0.04)
    #text = '#it{{f}}_{{tagged}}^{{pp}} = {:.2f}, #it{{f}}_{{tagged}}^{{AA}} = {:.2f}'.format(self.f_tagging_pp, self.f_tagging_AA)
    #text_latex.DrawLatex(0.57, 0.57, text)

    text_latex.SetTextSize(0.045)
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

  analysis = PlotLaura()
  analysis.plot_crosscheck()
