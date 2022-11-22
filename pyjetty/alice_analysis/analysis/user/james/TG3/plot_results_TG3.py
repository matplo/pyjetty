"""
macro for constructing TG3 plots

The script has to support a few things:
 - Plot either ratio alone, or distribution+ratio
 - Plot either data from hepdata, or custom data file
 - Plot a variety of theory predictions, either for ratio or distribution
 
To do this, we write a separate plotting function for each observable
We include generic functions for:
  - Plotting ratio
  - Plotting distribution+ratio
  - Plotting distribution
  - Initialize common config items

Author: James Mulligan (james.mulligan@berkeley.edu)
"""

# General
import os
import sys
import yaml
import argparse

# Data analysis and plotting
import ROOT
import numpy as np
import root_numpy

# Base class
sys.path.append('.')
from pyjetty.alice_analysis.analysis.base import common_base
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class PlotResults(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file='', output_dir='', **kwargs):
        super(PlotResults, self).__init__(**kwargs)
        self.config_file = config_file
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        
        self.utils = analysis_utils_obs.AnalysisUtils_Obs()
        self.utils.set_plotting_options()
        ROOT.gROOT.ForceStyle()
            
        self.data_color = ROOT.kGray+3
        self.data_markers = [25, 21]
        self.alpha = 0.7
        self.marker_size = 1.5
        self.line_width = 2
        self.line_style = 1
        self.theory_colors = [ROOT.kViolet-8, ROOT.kAzure-4, ROOT.kTeal-8, ROOT.kRed-7, ROOT.kPink+1, ROOT.kBlue-10]
                            # ROOT.kBlue-10, ROOT.kOrange+6, ROOT.kCyan-2
        
        print(self)
        
    #-------------------------------------------------------------------------------------------
    #
    #-------------------------------------------------------------------------------------------
    def plot_results(self):
        
        # Jet and hadron RAA -- data only
        self.plot_combined_raa()
        
        # Charged particle RAA
        self.plot_hadron_raa()

        # Inclusive full jet RAA
        self.plot_jet_raa()
        
        # Charged jet g
        self.plot_girth()
        
        # Charged jet mass
        self.plot_mass()
        
        # Charged Soft Drop zg
        #self.plot_zg()
        
        # Charged Soft Drop theta_g
        #self.plot_tg()
        
        # h-jet
        self.plot_hjet_IAA()
        self.plot_hjet_IAA_ratio()
        self.plot_hjet_dphi()
        
    #-------------------------------------------------------------------------------------------
    def plot_combined_raa(self):
    
        # Initialize data and theory info
        self.observable = 'combined_raa'
        self.init_result()
             
        # Plot ratio
        self.x_min = 5.
        self.x_max = 1e3
        self.y_ratio_min = 0.
        self.y_ratio_max = 2.
        self.ytitle = '#it{R}_{AA}'
        self.plot_ratio()
    
    #-------------------------------------------------------------------------------------------
    def plot_hadron_raa(self):
    
        # Initialize data and theory info
        self.observable = 'hadron_raa'
        self.init_result()
             
        # Plot ratio
        self.x_min = 5.
        self.x_max = 50.
        self.y_ratio_min = 0.
        self.y_ratio_max = 1.
        self.ytitle = '#it{R}_{AA}'
        self.plot_ratio()

    #-------------------------------------------------------------------------------------------
    def plot_jet_raa(self):

        # Initialize data and theory info
        self.observable = 'jet_raa'
        self.init_result()
             
        # Plot ratio
        self.x_min = 40.
        self.x_max = 140.
        self.y_ratio_min = 0.
        self.y_ratio_max = 1.4
        self.ytitle = '#it{R}_{AA}'
        self.plot_ratio()

    #-------------------------------------------------------------------------------------------
    def plot_girth(self):
    
        # Initialize data and theory info
        self.observable = 'girth'
        self.init_result(self_normalize=True)
             
        # Plot distribution + ratio
        self.y_max = 37
        self.y_ratio_min = 0.
        self.y_ratio_max = 2.49
        pp_label = 'PYTHIA6'
        self.ytitle = f'#frac{{1}}{{#it{{#sigma}}}} #frac{{d#it{{#sigma}}}}{{ d{self.xtitle} }}'
        self.plot_distribution_and_ratio(pp_label=pp_label)

    #-------------------------------------------------------------------------------------------
    def plot_mass(self):
    
        # Initialize data and theory info
        self.observable = 'mass'
        self.init_result(self_normalize=True)
             
        # Plot distribution
        self.x_min = 0.
        self.x_max = 22.
        self.y_ratio_min = -0.01
        self.y_ratio_max = 0.34
        pp_label = 'PYTHIA6'
        self.ytitle = f'#frac{{1}}{{#it{{#sigma}}}} #frac{{d#it{{#sigma}}}}{{ d{self.xtitle} }} (GeV^{{-1}}#it{{c}}^{{2}})'
        self.xtitle = f'{self.xtitle} (GeV/#it{{c}}^{{2}})'
        self.plot_distribution(pp_label=pp_label)
        #self.plot_distribution_and_ratio(pp_label=pp_label, plot_ratio=False)

    #-------------------------------------------------------------------------------------------
    def plot_zg(self):
        
        # Initialize data and theory info
        self.observable = 'zg'
        self.init_result()
             
        # Plot distribution + ratio
        self.y_max = 10.
        self.y_ratio_min = 0.6
        self.y_ratio_max = 1.59
        self.ytitle = f'#frac{{1}}{{#it{{#sigma}}_{{jet, inc}}}} #frac{{d#it{{#sigma}}}}{{ d{self.xtitle} }}'
        self.plot_distribution_and_ratio()

    #-------------------------------------------------------------------------------------------
    def plot_tg(self):
        
        # Initialize data and theory info
        self.observable = 'theta_g'
        self.init_result()
             
        # Plot distribution + ratio
        self.y_max = 4.1
        self.y_ratio_min = 0.3
        self.y_ratio_max = 1.89
        self.ytitle = f'#frac{{1}}{{#it{{#sigma}}_{{jet, inc}}}} #frac{{d#it{{#sigma}}}}{{ d{self.xtitle} }}'
        self.plot_distribution_and_ratio()

    #-------------------------------------------------------------------------------------------
    def plot_hjet_IAA(self):

        # Initialize data and theory info
        self.observable = 'hjet_IAA'
        self.init_result()
             
        # Plot ratio
        self.x_min = 20.
        self.x_max = 100.
        self.y_ratio_min = 0.
        self.y_ratio_max = 2.1
        self.ytitle = '#DeltaI_{AA} = #Delta_{recoil}^{PbPb} / #Delta_{recoil}^{PYTHIA}'
        self.plot_ratio()
        
    #-------------------------------------------------------------------------------------------
    def plot_hjet_IAA_ratio(self):

        # Initialize data and theory info
        self.observable = 'hjet_IAA_ratio'
        self.init_result()
             
        # Plot ratio
        self.x_min = 20.
        self.x_max = 100.
        self.y_ratio_min = 0.
        self.y_ratio_max = 2.1
        self.ytitle = '#Delta_{recoil}^{#it{R}=0.2} / #Delta_{recoil}^{#it{R}=0.5}'
        self.plot_ratio()
        
    #-------------------------------------------------------------------------------------------
    def plot_hjet_dphi(self):

        # Initialize data and theory info
        self.observable = 'hjet_dphi'
        self.init_result()
             
        # Plot ratio
        self.x_min = 1.5
        self.x_max = 3.2
        self.y_ratio_min = -0.005
        self.y_ratio_max = 0.1
        self.ytitle = '#Delta_{recoil}' # '#Phi(#Delta#it{#varphi})'
        self.plot_distribution()

    #---------------------------------------------------------------
    # Initialize a given observable
    #---------------------------------------------------------------
    def init_result(self, self_normalize=False):
    
        # Initialize an empty dict containing relevant info
        self.observable_settings = {}
        
        # Store theoretical predictions and labels in list
        self.observable_settings['prediction_distribution_list'] = []
        self.observable_settings['prediction_distribution_labels'] = []
        self.observable_settings['prediction_ratio_list'] = []
        self.observable_settings['prediction_ratio_labels'] = []
    
        # Load config file
        with open(self.config_file, 'r') as stream:
            config = yaml.safe_load(stream)
            result = config[self.observable]
            
        # Common settings
        self.xtitle = result['xtitle']
        self.sqrts = result['sqrts']
        self.centrality = result['centrality']
        self.jetR = result['jetR']
        if 'pt' in result:
            self.pt = result['pt']
        if 'eta_R' in result:
            self.eta_R = result['eta_R']
            self.eta_cut = np.round(self.eta_R - self.jetR, decimals=1)
        if 'eta_cut' in result:
            self.eta_cut = result['eta_cut']
        if 'obs_label' in result:
            self.obs_label = result['obs_label']
        if 'zcut' in result:
            self.zcut = result['zcut']
        if 'theory_config' in result:
            self.theory_config = result['theory_config']
  
        # Data -- distribution
        if self.observable in ['combined_raa']:
            file = ROOT.TFile(result['file'], 'READ')
            self.observable_settings['g_hadron_raa_ALICE'] = file.Get(result['hadron_raa_ALICE_name'])
            self.observable_settings['g_hadron_raa_sys_ALICE'] = file.Get(result['hadron_raa_sys_ALICE_name'])
            self.observable_settings['g_hadron_raa_CMS'] = file.Get(result['hadron_raa_CMS_name'])
            self.observable_settings['g_hadron_raa_sys_CMS'] = file.Get(result['hadron_raa_sys_CMS_name'])
            self.observable_settings['g_jet_raa_ALICE'] = file.Get(result['jet_raa_ALICE_name'])
            self.observable_settings['g_jet_raa_sys_corr_ALICE'] = file.Get(result['jet_raa_corr_ALICE_name'])
            self.observable_settings['g_jet_raa_sys_shape_ALICE'] = file.Get(result['jet_raa_shape_ALICE_name'])
            self.observable_settings['g_jet_raa_ATLAS'] = file.Get(result['jet_raa_ATLAS_name'])
            self.observable_settings['g_jet_raa_sys_ATLAS'] = file.Get(result['jet_raa_sys_ATLAS_name'])
        elif self.observable in ['hadron_raa']:
            file = ROOT.TFile(result['file'], 'READ')
            self.observable_settings['g_hadron_raa_ALICE'] = file.Get(result['hadron_raa_ALICE_name'])
            self.observable_settings['g_hadron_raa_sys_ALICE'] = file.Get(result['hadron_raa_sys_ALICE_name'])
        elif self.observable == 'jet_raa':
            f = ROOT.TFile(result['hepdata'], 'READ')
            dir = f.Get('Table 30')
            h = dir.Get('Hist1D_y1')
            y = np.asarray(root_numpy.hist2array(h), dtype=np.float64)
            stat = np.asarray(root_numpy.hist2array(dir.Get('Hist1D_y1_e1')), dtype=np.float64)
            sys_corr = np.asarray(root_numpy.hist2array(dir.Get('Hist1D_y1_e2')), dtype=np.float64)
            sys_shape = np.asarray(root_numpy.hist2array(dir.Get('Hist1D_y1_e3')), dtype=np.float64)
            
            bins = np.array(h.GetXaxis().GetXbins())
            x = (bins[1:] + bins[:-1]) / 2
            x_err = np.ediff1d(bins)/2
            n = len(x)
            f.Close()

            self.observable_settings['g_jet_raa_ALICE'] = ROOT.TGraphErrors(n, x, y, np.zeros(n), stat)
            self.observable_settings['g_jet_raa_sys_corr_ALICE']  = ROOT.TGraphErrors(n, x, y, x_err, sys_corr)
            self.observable_settings['g_jet_raa_sys_shape_ALICE']  = ROOT.TGraphErrors(n, x, y, x_err, sys_shape)
        elif self.observable in ['girth', 'mass']:
            file = ROOT.TFile(result['file'], 'READ')
            self.observable_settings['g_pp'] = file.Get(result['name_pp'])
            self.observable_settings['h_AA_stat'] = file.Get(result['name_AA_stat'])
            self.observable_settings['g_AA_sys'] = file.Get(result['name_AA_sys'])
            self.observable_settings['h_AA_stat'].SetDirectory(0)
            self.bins = np.array(self.observable_settings['h_AA_stat'].GetXaxis().GetXbins())
        elif self.observable in ['zg', 'theta_g']:
            self.file_pp_name = result['file_pp']
            self.file_AA_name = result['file_AA']
            self.file_pp = ROOT.TFile(self.file_pp_name, 'READ')
            self.file_AA = ROOT.TFile(self.file_AA_name, 'READ')
            main_result_name = f'hmain_{self.observable}_R{self.jetR}_{self.obs_label}_{self.pt[0]}-{self.pt[1]}'
            sys_total_name = f'hResult_{self.observable}_systotal_R{self.jetR}_{self.obs_label}_{self.pt[0]}-{self.pt[1]}'
            self.observable_settings['h_pp'] = self.file_pp.Get(main_result_name)
            self.observable_settings['h_pp_sys'] = self.file_pp.Get(sys_total_name)
            self.observable_settings['h_AA_stat'] = self.file_AA.Get(main_result_name)
            self.observable_settings['h_AA_sys'] = self.file_AA.Get(sys_total_name)
            self.observable_settings['h_pp'].SetDirectory(0)
            self.observable_settings['h_pp_sys'].SetDirectory(0)
            self.observable_settings['h_AA_stat'].SetDirectory(0)
            self.observable_settings['h_AA_sys'].SetDirectory(0)
            self.bins = np.array(self.observable_settings['h_AA_stat'].GetXaxis().GetXbins())[1:]

            h_tagging_fraction_name = f'h_tagging_fraction_R{self.jetR}_{self.obs_label}'
            self.h_tagging_fraction_pp = self.file_pp.Get(h_tagging_fraction_name)
            self.h_tagging_fraction_AA = self.file_AA.Get(h_tagging_fraction_name)
            bin = self.h_tagging_fraction_pp.GetXaxis().FindBin((self.pt[0] + self.pt[1])/2)
            self.f_tagging_pp = self.h_tagging_fraction_pp.GetBinContent(bin)
            bin = self.h_tagging_fraction_AA.GetXaxis().FindBin((self.pt[0] + self.pt[1])/2)
            self.f_tagging_AA = self.h_tagging_fraction_AA.GetBinContent(bin)
        elif self.observable == 'hjet_dphi':
            f = ROOT.TFile(result['hepdata'], 'READ')
            dir = f.Get(result['dir'])
            self.observable_settings['g_pp'] = dir.Get('Graph1D_y1')
            
            bins = np.array(result['x_bins'])
            y = np.array(result['y_PbPb'])
            y_stat = np.array(result['y_stat_PbPb'])
            x = (bins[1:] + bins[:-1]) / 2
            x_err = (bins[1:] - bins[:-1]) / 2
            n = len(x)
            self.observable_settings['g_AA'] = ROOT.TGraphErrors(n, x, y, x_err, y_stat)
            
            integral_pp_data = dir.Get('Hist1D_y1').Integral()
            integral_PbPb_data = np.sum(y)

        # Data -- ratio
        if self.observable == 'girth':
            #self.observable_settings['hepdata_ratio'] = result['hepdata_ratio']
            #f = ROOT.TFile(self.inclusive_chjet_observables['g_alice']['hepdata'], 'READ')
            #dir = f.Get('Table 11')
            #h_data = dir.Get('Graph1D_y1')
            #h_data_list.append([h_data, '0-10%'])
            #f.Close()
            self.observable_settings['h_ratio_stat'] = file.Get(result['name_ratio_stat'])
            self.observable_settings['g_ratio_sys'] = file.Get(result['name_ratio_sys'])
            self.observable_settings['h_ratio_stat'].SetDirectory(0)
        elif self.observable in ['zg', 'theta_g']:
            self.observable_settings['h_ratio_stat'] = self.observable_settings['h_AA_stat'].Clone()
            self.observable_settings['h_ratio_stat'].Divide(self.observable_settings['h_pp'])
            self.observable_settings['h_ratio_sys'] = self.observable_settings['h_AA_sys'].Clone()
            self.observable_settings['h_ratio_sys'].Divide(self.observable_settings['h_pp_sys'])
        elif 'hjet_IAA' in self.observable:
            f = ROOT.TFile(result['hepdata'], 'READ')
            dir = f.Get(result['dir'])
            h = dir.Get('Hist1D_y1')
            y = np.asarray(root_numpy.hist2array(h), dtype=np.float64)
            stat = np.asarray(root_numpy.hist2array(dir.Get('Hist1D_y1_e1')), dtype=np.float64)
            sys_shape_plus = np.asarray(root_numpy.hist2array(dir.Get('Hist1D_y1_e2plus')), dtype=np.float64)
            sys_shape_minus = np.absolute(np.asarray(root_numpy.hist2array(dir.Get('Hist1D_y1_e2minus')), dtype=np.float64))
            sys_corr_plus = np.asarray(root_numpy.hist2array(dir.Get('Hist1D_y1_e3plus')), dtype=np.float64)
            sys_corr_minus = np.absolute(np.asarray(root_numpy.hist2array(dir.Get('Hist1D_y1_e3minus')), dtype=np.float64))

            bins = np.array(h.GetXaxis().GetXbins())
            x = (bins[1:] + bins[:-1]) / 2
            x_err = np.ediff1d(bins)/2
            n = len(x)
            f.Close()

            self.observable_settings['g_jet_hjet_IAA_ALICE'] = ROOT.TGraphErrors(n, x, y, np.zeros(n), stat)
            self.observable_settings['g_jet_hjet_IAA_sys_corr_ALICE']  = ROOT.TGraphAsymmErrors(n, x, y, x_err, x_err, sys_corr_minus, sys_corr_plus)
            self.observable_settings['g_jet_hjet_IAA_sys_shape_ALICE']  = ROOT.TGraphAsymmErrors(n, x, y, x_err, x_err, sys_shape_minus, sys_shape_plus)
        
        # Theory -- JETSCAPE
        if 'file_jetscape_pp' in result:
            
            f_jetscape_pp = ROOT.TFile(result['file_jetscape_pp'], 'READ')
            f_jetscape_AA = ROOT.TFile(result['file_jetscape_AA'], 'READ')
            if 'hname_jetscape' in result:
                hname = result['hname_jetscape']
                h_jetscape_pp = f_jetscape_pp.Get(hname)
                h_jetscape_pp.SetDirectory(0)
                f_jetscape_pp.Close()
                
                h_jetscape_AA = f_jetscape_AA.Get(hname)
                h_jetscape_AA.SetDirectory(0)
                #if self.observable == 'jet_raa':
                #    h_jetscape_AA.Scale(0.5) # Since we hadd [0,5] and [5,10]
                f_jetscape_AA.Close()

            ## For hadrons, impose a 1 GeV minimum, and subtract the recoil hadrons
            if self.observable == 'hadron_raa':
                # Impose 1 Gev minimum
                h_jetscape_pp_xbins = np.array(h_jetscape_pp.GetXaxis().GetXbins())
                h_jetscape_pp_xbins = h_jetscape_pp_xbins[(h_jetscape_pp_xbins>=8.)]
                h_jetscape_pp = h_jetscape_pp.Rebin(h_jetscape_pp_xbins.size-1, f'{hname}_pp_rebinned', h_jetscape_pp_xbins)

                h_jetscape_AA_xbins = np.array(h_jetscape_AA.GetXaxis().GetXbins())
                h_jetscape_AA_xbins = h_jetscape_AA_xbins[(h_jetscape_AA_xbins>=8.)]
                h_jetscape_AA = h_jetscape_AA.Rebin(h_jetscape_AA_xbins.size-1, f'{hname}_AA_rebinned', h_jetscape_AA_xbins)

                # Subtract holes
                f_jetscape_AA = ROOT.TFile(result['file_jetscape_AA'], 'READ')
                h_recoil_name = 'h_hadron_pt_alice_holesScaled'
                h_recoil = f_jetscape_AA.Get(h_recoil_name)
                h_recoil.SetDirectory(0)
                h_recoil_rebinned = h_recoil.Rebin(h_jetscape_AA_xbins.size-1, f'{h_recoil_name}_AA', h_jetscape_AA_xbins)
                h_jetscape_AA.Add(h_recoil_rebinned, -1)
                f_jetscape_AA.Close()
            
            # Construct IAA
            if 'hjet' in self.observable:
           
                if self.observable == 'hjet_IAA_ratio':
                    h_pp_ntrigger = f_jetscape_pp.Get(result['hname_jetscape_ntrigger_R05'])
                    h_pp_ntrigger.SetDirectory(0)
                    h_pp_high = f_jetscape_pp.Get(result['hname_jetscape_high_R05'])
                    h_pp_high.SetDirectory(0)
                    h_pp_low = f_jetscape_pp.Get(result['hname_jetscape_low_R05'])
                    h_pp_low.SetDirectory(0)
                    f_jetscape_pp.Close()
                    
                    h_AA_ntrigger = f_jetscape_AA.Get(result['hname_jetscape_ntrigger_R02'])
                    h_AA_ntrigger.SetDirectory(0)
                    h_AA_high = f_jetscape_AA.Get(result['hname_jetscape_high_R02'])
                    h_AA_high.SetDirectory(0)
                    h_AA_low = f_jetscape_AA.Get(result['hname_jetscape_low_R02'])
                    h_AA_low.SetDirectory(0)
                    f_jetscape_AA.Close()
                    
                else:
                    h_pp_ntrigger = f_jetscape_pp.Get(result['hname_jetscape_ntrigger'])
                    h_pp_ntrigger.SetDirectory(0)
                    h_pp_high = f_jetscape_pp.Get(result['hname_jetscape_high'])
                    h_pp_high.SetDirectory(0)
                    h_pp_low = f_jetscape_pp.Get(result['hname_jetscape_low'])
                    h_pp_low.SetDirectory(0)
                    f_jetscape_pp.Close()
                    
                    h_AA_ntrigger = f_jetscape_AA.Get(result['hname_jetscape_ntrigger'])
                    h_AA_ntrigger.SetDirectory(0)
                    h_AA_high = f_jetscape_AA.Get(result['hname_jetscape_high'])
                    h_AA_high.SetDirectory(0)
                    h_AA_low = f_jetscape_AA.Get(result['hname_jetscape_low'])
                    h_AA_low.SetDirectory(0)
                    f_jetscape_AA.Close()

                # Delta recoil
                n_trig_high_pp = h_pp_ntrigger.GetBinContent(h_pp_ntrigger.FindBin(30.))
                n_trig_low_pp = h_pp_ntrigger.GetBinContent(h_pp_ntrigger.FindBin(8.5))
                n_trig_high_AA = h_AA_ntrigger.GetBinContent(h_AA_ntrigger.FindBin(30.))
                n_trig_low_AA = h_AA_ntrigger.GetBinContent(h_AA_ntrigger.FindBin(8.5))

                h_pp_high.Scale(1./n_trig_high_pp, 'width')
                h_pp_low.Scale(1./n_trig_low_pp, 'width')
                h_AA_high.Scale(1./n_trig_high_AA, 'width')
                h_AA_low.Scale(1./n_trig_low_AA, 'width')
                
                h_jetscape_pp = h_pp_high.Clone('h_delta_recoil_pp')
                h_jetscape_AA = h_AA_high.Clone('h_delta_recoil_AA')
                if self.observable == 'hjet_IAA_ratio':
                    h_jetscape_pp.Add(h_pp_low, -1*result['c_ref_R05'])
                    h_jetscape_AA.Add(h_AA_low, -1*result['c_ref_R02'])
                else:
                    h_jetscape_pp.Add(h_pp_low, -1)
                    h_jetscape_AA.Add(h_AA_low, -1*result['c_ref'])
            
            # Normalization
            h_jetscape_pp.Scale(1., 'width')
            h_jetscape_AA.Scale(1., 'width')
            if self.observable == 'hjet_IAA_ratio':
                # Eta acceptance -- (0.9-R)*2
                h_jetscape_pp.Scale(1./(2*(0.9-0.5))) # R=0.5
                h_jetscape_AA.Scale(1./(2*(0.9-0.2))) # R=0.2
            if self_normalize:
                if self.observable in ['zg', 'theta_g']:
                    min_bin = 0
                else:
                    min_bin = 1
                h_jetscape_AA.Scale(1./h_jetscape_AA.Integral(min_bin, h_jetscape_AA.GetNbinsX()))
                h_jetscape_pp.Scale(1./h_jetscape_pp.Integral(min_bin, h_jetscape_pp.GetNbinsX()))
                    
            # Form ratio
            h_ratio_jetscape = h_jetscape_AA
            h_ratio_jetscape.Divide(h_jetscape_pp)
            self.observable_settings['prediction_ratio_list'].append(h_ratio_jetscape)
            self.observable_settings['prediction_ratio_labels'].append('JETSCAPE')

            # Load PYTHIA theory prediction for IAA R=0.2/0.5 ratio
            if self.observable == 'hjet_IAA_ratio':

                bins = np.array(h_ratio_jetscape.GetXaxis().GetXbins())
                x = (bins[1:] + bins[:-1]) / 2
                n = len(x)
                y_pythia = np.array(result['y_pythia'])
                y_err_pythia = np.array(result['y_err_pythia'])
                g_pythia = ROOT.TGraphErrors(n, x, y_pythia, np.zeros(n), y_err_pythia)
                self.observable_settings['prediction_ratio_list'].append(g_pythia)
                self.observable_settings['prediction_ratio_labels'].append('PYTHIA')
            
            if self.observable in ['mass']:
                self.observable_settings['prediction_distribution_list'].append(h_jetscape_AA)
                self.observable_settings['prediction_distribution_labels'].append('JETSCAPE')
                
        # Theory -- folded JETSCAPE
        if self.observable == 'hjet_dphi':
            file_pp = ROOT.TFile(result['file_pp'], 'READ')
            file_AA = ROOT.TFile(result['file_AA'], 'READ')
            h_pp = file_pp.Get(result['hname'])
            h_AA = file_AA.Get(result['hname'])
            h_pp.SetDirectory(0)
            h_AA.SetDirectory(0)
            
            # Also add versions that are scaled to integral of data
            integral_pp_jetscape = h_pp.Integral(1, h_pp.GetNbinsX()+1)
            integral_PbPb_jetscape = h_AA.Integral(1, h_AA.GetNbinsX()+1)
            h_pp_scaled = h_pp.Clone(f'{h_pp.GetName()}_scaled')
            h_AA_scaled = h_AA.Clone(f'{h_AA.GetName()}_scaled')
            h_pp_scaled.SetDirectory(0)
            h_AA_scaled.SetDirectory(0)
            h_pp_scaled.Scale(integral_pp_data/integral_pp_jetscape)
            h_AA_scaled.Scale(integral_PbPb_data/integral_PbPb_jetscape)
            
            self.observable_settings['prediction_distribution_list'].append(h_pp)
            self.observable_settings['prediction_distribution_labels'].append('JETSCAPE, pp')
            #self.observable_settings['prediction_distribution_list'].append(h_pp_scaled)
            #self.observable_settings['prediction_distribution_labels'].append('JETSCAPE, pp (scaled)')
            self.observable_settings['prediction_distribution_list'].append(h_AA)
            self.observable_settings['prediction_distribution_labels'].append('JETSCAPE, Pb-Pb')
            self.observable_settings['prediction_distribution_list'].append(h_AA_scaled)
            self.observable_settings['prediction_distribution_labels'].append('JETSCAPE, Pb-Pb (scaled)')

        # Theory -- Hybrid
        if self.observable == 'girth':
            self.observable_settings['prediction_ratio_list'].append(file.Get(result['name_ratio_hybrid_wake0_lres0']))
            self.observable_settings['prediction_ratio_labels'].append('Hybrid Model, #it{L}_{res} = 0, wake off')
            
            self.observable_settings['prediction_ratio_list'].append(file.Get(result['name_ratio_hybrid_wake1_lres0']))
            self.observable_settings['prediction_ratio_labels'].append('Hybrid Model, #it{L}_{res} = 0, wake on')

            #self.observable_settings['prediction_ratio_list'].append(file.Get(result['name_ratio_hybrid_wake0_lres2']))
            #self.observable_settings['prediction_ratio_labels'].append('Hybrid Model, #it{L}_{res} = 2#pi/T, wake off')

            self.observable_settings['prediction_ratio_list'].append(file.Get(result['name_ratio_hybrid_wake0_lresinf']))
            self.observable_settings['prediction_ratio_labels'].append('Hybrid Model, #it{L}_{res} = #infty, wake off')
            file.Close()
            
        if self.observable == 'mass':
            self.observable_settings['prediction_distribution_list'].append(file.Get(result['name_hybrid_wake0_lres0']))
            self.observable_settings['prediction_distribution_labels'].append('Hybrid Model, #it{L}_{res} = 0, wake off')
            
            self.observable_settings['prediction_distribution_list'].append(file.Get(result['name_hybrid_wake1_lres0']))
            self.observable_settings['prediction_distribution_labels'].append('Hybrid Model, #it{L}_{res} = 0, wake on')
            
        if self.observable == 'jet_raa':
            hybrid_x = np.array(result['hybrid_x'])
            n = len(hybrid_x)
            hybrid_lres0_upper = np.array(result['hybrid_lres0_upper'])
            hybrid_lres0_lower = np.array(result['hybrid_lres0_lower'])
            g_hybrid_lres0 = ROOT.TGraphAsymmErrors(n, hybrid_x, hybrid_lres0_upper, np.zeros(n), np.zeros(n), hybrid_lres0_upper-hybrid_lres0_lower, np.zeros(n))
            self.observable_settings['prediction_ratio_list'].append(g_hybrid_lres0)
            self.observable_settings['prediction_ratio_labels'].append('Hybrid model, #it{L}_{res}=0')
            
            hybrid_lres2piT_upper = np.array(result['hybrid_lres2piT_upper'])
            hybrid_lres2piT_lower = np.array(result['hybrid_lres2piT_lower'])
            g_hybrid_lres2piT = ROOT.TGraphAsymmErrors(n, hybrid_x, hybrid_lres2piT_upper, np.zeros(n), np.zeros(n), hybrid_lres2piT_upper-hybrid_lres2piT_lower, np.zeros(n))
            self.observable_settings['prediction_ratio_list'].append(g_hybrid_lres2piT)
            self.observable_settings['prediction_ratio_labels'].append('Hybrid model, #it{L}_{res}=2#pi/#it{T}')

        # Theory -- JEWEL
        if self.observable == 'mass':
            self.observable_settings['prediction_distribution_list'].append(file.Get(result['name_jewel_rec']))
            self.observable_settings['prediction_distribution_labels'].append('JEWEL, recoils on')
            
            self.observable_settings['prediction_distribution_list'].append(file.Get(result['name_jewel_norec']))
            self.observable_settings['prediction_distribution_labels'].append('JEWEL, recoils off')
            
        if self.observable == 'hadron_raa':
            if 'jewel_x' in result:
                jewel_x = np.array(result['jewel_x'])
                jewel_RAA = np.array(result['jewel_RAA'])
                n = len(jewel_x)
                g = ROOT.TGraph(n, jewel_x, jewel_RAA)
                self.observable_settings['prediction_ratio_list'].append(g)
                self.observable_settings['prediction_ratio_labels'].append('JEWEL, recoils on?')
            
        if self.observable == 'jet_raa':
            jewel_x = np.array(result['jewel_x'])
            n = len(jewel_x)
            jewel_recoils_on_RAA = np.array(result['jewel_recoils_on_RAA'])
            g_recoils_on = ROOT.TGraph(n, jewel_x, jewel_recoils_on_RAA)
            self.observable_settings['prediction_ratio_list'].append(g_recoils_on)
            self.observable_settings['prediction_ratio_labels'].append('JEWEL, recoils on')
            
            jewel_recoils_off_RAA = np.array(result['jewel_recoils_off_RAA'])
            g_recoils_off = ROOT.TGraph(n, jewel_x, jewel_recoils_off_RAA)
            self.observable_settings['prediction_ratio_list'].append(g_recoils_off)
            self.observable_settings['prediction_ratio_labels'].append('JEWEL, recoils off')
            
        # Theory -- zg/tg
        if self.observable in ['zg', 'theta_g']:
            self.bin_array = np.array(self.observable_settings['h_pp'].GetXaxis().GetXbins())[1:]

            with open(self.theory_config, 'r') as stream:
                config = yaml.safe_load(stream)
                  
            self.label_list = []
            self.sublabel_list = []
            self.prediction_g_list = []
          
            self.prediction_types = [name for name in list(config.keys())]
            for type in self.prediction_types:
          
                theory = config[type]
                configs = [name for name in list(theory.keys()) if 'config' in name ]
                for prediction in configs:
                    theory_prediction = theory[prediction]
                    config_jetR = theory_prediction['jetR']
                    if config_jetR == self.jetR:
                    
                        plot_list = theory['plot_list']
                        
                        if type in ['yang_ting', 'guangyou']:
                            continue
                    
                        if type == 'lbnl':
                        
                            if not 'med #it{q}/#it{g}' in theory_prediction['sublabel']:
                                continue
                        
                            x = np.array(theory_prediction['x'])
                            y_pp = np.array(theory_prediction['y_pp'])

                            n = len(x)
                            xerr = np.zeros(n)

                            if 'f_q' in theory_prediction:
                                f_q = theory_prediction['f_q']
                                f_q_min = theory_prediction['f_q_min']
                                f_q_max = theory_prediction['f_q_max']
                                y_AA_quark = np.array(theory_prediction['y_AA_quark'])
                                y_AA_gluon = np.array(theory_prediction['y_AA_gluon'])

                                y_AA = f_q*y_AA_quark + (1-f_q)*y_AA_gluon
                                y_AA_min = f_q_min*y_AA_quark + (1-f_q_min)*y_AA_gluon
                                y_AA_max = f_q_max*y_AA_quark + (1-f_q_max)*y_AA_gluon

                                ratio = np.divide(y_AA, y_pp)
                                ratio_lower = np.divide(y_AA_min, y_pp)
                                ratio_upper = np.divide(y_AA_max, y_pp)
                                g = ROOT.TGraphAsymmErrors(n, x, ratio, xerr, xerr, ratio-ratio_lower, ratio_upper-ratio)

                            else:
                                y_AA = np.array(theory_prediction['y_AA'])
                                ratio = np.divide(y_AA, y_pp)
                                g = ROOT.TGraph(n, x, ratio)

                        elif type == 'caucal':
                                    
                            x = np.array(theory_prediction['x'])
                            ratio = np.array(theory_prediction['ratio'])
                            ratio_neg_unc_tot = np.array(theory_prediction['ratio_neg_unc_tot'])
                            ratio_pos_unc_tot = np.array(theory_prediction['ratio_pos_unc_tot'])

                            n = len(x)
                            xerr = np.zeros(n)
                            g = ROOT.TGraphAsymmErrors(n, x, ratio, xerr, xerr, ratio_neg_unc_tot, ratio_pos_unc_tot)
                        
                        elif type == 'guangyou':
                        
                            x = np.array(theory_prediction['x'])
                            ratio = np.array(theory_prediction['ratio'])

                            n = len(x)
                            g = ROOT.TGraph(n, x, ratio)
                        
                        elif type in ['hybrid_model', 'yang_ting']:
                        
                            x = np.array(theory_prediction['x'])
                            ratio_lower = np.array(theory_prediction['ratio_lower'])
                            ratio_upper = np.array(theory_prediction['ratio_upper'])
                            ratio = (ratio_lower + ratio_upper) / 2.

                            n = len(x)
                            xerr = np.zeros(n)
                            g = ROOT.TGraphAsymmErrors(n, x, ratio, xerr, xerr, ratio-ratio_lower, ratio_upper-ratio)
                        
                        elif type == 'jetscape':

                            bins = np.array(theory_prediction['xbins'])
                            x = (bins[1:] + bins[:-1]) / 2
                            ratio = np.array(theory_prediction['ratio'])
                            ratio_err = np.array(theory_prediction['ratio_err'])
                                
                            n = len(x)
                            xerr = np.zeros(n)
                            g = ROOT.TGraphErrors(n, x, ratio, xerr, ratio_err)
        
                        if prediction in plot_list:
                            self.observable_settings['prediction_ratio_list'].append(g)
                            label = theory['label']
                            sublabel = theory_prediction['sublabel']
                            self.observable_settings['prediction_ratio_labels'].append(f'{label}{sublabel}')

    #-------------------------------------------------------------------------------------------
    # Plot ratio in single panel
    #-------------------------------------------------------------------------------------------
    def plot_ratio(self):
                  
        cname = 'c'
        c = ROOT.TCanvas(cname, cname, 600, 450)
        c.SetRightMargin(0.05)
        c.SetLeftMargin(0.15)
        c.SetTopMargin(0.05)
        c.SetBottomMargin(0.17)
        c.SetTicks(0,1)
        if self.observable == 'combined_raa':
            c.SetLogx()
        c.cd()
        
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, self.x_min, self.x_max)
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(self.xtitle)
        myBlankHisto.SetYTitle(self.ytitle)
        myBlankHisto.SetMaximum(self.y_ratio_max)
        myBlankHisto.GetXaxis().SetTitleSize(0.06)
        myBlankHisto.GetXaxis().SetTitleOffset(1.2)
        myBlankHisto.GetXaxis().SetLabelSize(0.05)
        myBlankHisto.GetYaxis().SetTitleSize(0.07)
        myBlankHisto.GetYaxis().SetTitleOffset(1.)
        myBlankHisto.GetYaxis().SetLabelSize(0.05)
        myBlankHisto.Draw('E')
        
        if self.observable == 'combined_raa':
            y_leg = 0.62
        elif self.observable in ['hadron_raa', 'jet_raa']:
            y_leg = 0.67
        elif self.observable in ['hjet_IAA', 'hjet_IAA_ratio']:
            y_leg = 0.58
        leg = ROOT.TLegend(0.2,y_leg,0.4,y_leg+0.11)
        self.utils.setup_legend(leg, 0.04, sep=-0.1)
        
        if self.observable == 'jet_raa':
            leg2 = ROOT.TLegend(0.58,y_leg-0.09,0.78,y_leg+0.11)
        elif self.observable in ['hjet_IAA', 'hjet_IAA_ratio']:
            leg2 = ROOT.TLegend(0.53,y_leg,0.73,y_leg+0.11)
        else:
            leg2 = ROOT.TLegend(0.6,y_leg,0.8,y_leg+0.11)
        self.utils.setup_legend(leg2, 0.04, sep=-0.1)
        
        # Draw theory predictions
        for i, prediction in enumerate(self.observable_settings['prediction_ratio_list']):
        
            label = self.observable_settings['prediction_ratio_labels'][i]
            color = self.theory_colors[i]
            if self.observable == 'hjet_IAA_ratio' and label == 'PYTHIA':
                color = self.theory_colors[i+1]
            
            prediction.SetFillColor(color)
            prediction.SetFillColorAlpha(color, self.alpha)
            prediction.SetLineColor(color)
            prediction.SetLineColorAlpha(color, self.alpha)
            prediction.SetLineWidth(8)
            prediction.SetMarkerSize(0)
            prediction.SetMarkerStyle(0)
        
            if type(prediction) in [ROOT.TH1F]:
            
                prediction.Draw('E3 same')
                leg2.AddEntry(prediction, label, 'L')

            elif type(prediction) in [ROOT.TGraph, ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                 
                if type(prediction) in [ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                    prediction.Draw('3 same')
                elif type(prediction) == ROOT.TGraph:
                    prediction.Draw('same')

                leg2.AddEntry(prediction, label, 'L')
           
        # Draw experimental data
        if self.observable == 'combined_raa':

            for name in self.observable_settings.keys():
                if 'hadron_raa' in name:
                    if 'ALICE' in name:
                        color = self.theory_colors[0]
                    elif 'CMS' in name:
                        color = self.theory_colors[1]
                if 'jet_raa' in name:
                    if 'ALICE' in name:
                        color = self.theory_colors[2]
                    elif 'ATLAS' in name:
                        color = self.theory_colors[3]
                
                if 'raa_sys' in name:
                    self.observable_settings[name].SetFillColor(color)
                    self.observable_settings[name].SetFillColorAlpha(color, 0.3)
                    self.observable_settings[name].SetMarkerStyle(0)
                    if 'corr' in name:
                        self.observable_settings[name].SetFillStyle(0)
                        self.observable_settings[name].SetLineWidth(1)
                        self.observable_settings[name].SetLineColor(color)
                    else:
                        self.observable_settings[name].SetFillStyle(1001)
                        self.observable_settings[name].SetLineWidth(0)
                        self.observable_settings[name].SetLineColor(0)
                elif 'raa' in name:
                    self.observable_settings[name].SetMarkerColor(color)
                    self.observable_settings[name].SetLineColor(color)
                    self.observable_settings[name].SetLineWidth(1)
                    self.observable_settings[name].SetMarkerSize(1.3)
                    # Need to set xerrors to 0...
                    for i in range(self.observable_settings[name].GetN()):
                        yerr = self.observable_settings[name].GetErrorY(i)
                        if type(self.observable_settings[name]) == ROOT.TGraphAsymmErrors:
                            self.observable_settings[name].SetPointError(i, 0, 0, yerr, yerr)
                        elif type(self.observable_settings[name]) == ROOT.TGraphErrors:
                            self.observable_settings[name].SetPointError(i, 0, yerr)
            
            self.observable_settings['g_hadron_raa_CMS'].SetMarkerStyle(25)
            self.observable_settings['g_hadron_raa_CMS'].Draw('PE Z same')
            self.observable_settings['g_hadron_raa_sys_CMS'].Draw('E2 same')
            
            self.observable_settings['g_jet_raa_ATLAS'].SetMarkerStyle(24)
            self.observable_settings['g_jet_raa_ATLAS'].Draw('PE Z same')
            self.observable_settings['g_jet_raa_sys_ATLAS'].Draw('E2 same')
            
            self.observable_settings['g_hadron_raa_ALICE'].SetMarkerStyle(21)
            self.observable_settings['g_hadron_raa_ALICE'].Draw('PE Z same')
            self.observable_settings['g_hadron_raa_sys_ALICE'].Draw('E2 same')
            
            self.observable_settings['g_jet_raa_ALICE'].SetMarkerStyle(20)
            self.observable_settings['g_jet_raa_ALICE'].Draw('PE Z same')
            self.observable_settings['g_jet_raa_sys_corr_ALICE'].Draw('E2 same')
            self.observable_settings['g_jet_raa_sys_shape_ALICE'].Draw('E2 same')
            
            leg.AddEntry(self.observable_settings['g_hadron_raa_ALICE'], 'ALICE','PE')
            leg.AddEntry(self.observable_settings['g_hadron_raa_CMS'], 'CMS','PE')
            leg2.AddEntry(self.observable_settings['g_jet_raa_ALICE'], 'ALICE','PE')
            leg2.AddEntry(self.observable_settings['g_jet_raa_ATLAS'], 'ATLAS','PE')
            
        elif self.observable == 'hadron_raa':

            self.observable_settings['g_hadron_raa_sys_ALICE'].SetFillColor(self.data_color)
            self.observable_settings['g_hadron_raa_sys_ALICE'].SetFillColorAlpha(self.data_color, 0.3)
            self.observable_settings['g_hadron_raa_sys_ALICE'].SetMarkerStyle(0)
            self.observable_settings['g_hadron_raa_sys_ALICE'].SetFillStyle(1001)
            self.observable_settings['g_hadron_raa_sys_ALICE'].SetLineWidth(0)
            self.observable_settings['g_hadron_raa_sys_ALICE'].SetLineColor(0)
            
            self.observable_settings['g_hadron_raa_ALICE'].SetMarkerColor(self.data_color)
            self.observable_settings['g_hadron_raa_ALICE'].SetLineColor(self.data_color)
            self.observable_settings['g_hadron_raa_ALICE'].SetLineWidth(1)
            self.observable_settings['g_hadron_raa_ALICE'].SetMarkerSize(1.3)
            # Need to set xerrors to 0...
            for i in range(self.observable_settings['g_hadron_raa_ALICE'].GetN()):
                yerr = self.observable_settings['g_hadron_raa_ALICE'].GetErrorY(i)
                if type(self.observable_settings['g_hadron_raa_ALICE']) == ROOT.TGraphAsymmErrors:
                    self.observable_settings['g_hadron_raa_ALICE'].SetPointError(i, 0, 0, yerr, yerr)

            self.observable_settings['g_hadron_raa_ALICE'].SetMarkerStyle(21)
            self.observable_settings['g_hadron_raa_ALICE'].Draw('PE Z same')
            self.observable_settings['g_hadron_raa_sys_ALICE'].Draw('E2 same')

            leg.AddEntry(self.observable_settings['g_hadron_raa_ALICE'], '0-10%','PE')
            leg.AddEntry(self.observable_settings['g_hadron_raa_sys_ALICE'], 'Sys. uncertainty','f')
            
        elif self.observable == 'jet_raa':
                    
            self.observable_settings['g_jet_raa_sys_shape_ALICE'].SetFillColor(self.data_color)
            self.observable_settings['g_jet_raa_sys_shape_ALICE'].SetFillColorAlpha(self.data_color, 0.3)
            self.observable_settings['g_jet_raa_sys_shape_ALICE'].SetMarkerStyle(0)
            self.observable_settings['g_jet_raa_sys_shape_ALICE'].SetFillStyle(1001)
            self.observable_settings['g_jet_raa_sys_shape_ALICE'].SetLineWidth(0)
            self.observable_settings['g_jet_raa_sys_shape_ALICE'].SetLineColor(0)
            
            self.observable_settings['g_jet_raa_sys_corr_ALICE'].SetFillColor(self.data_color)
            self.observable_settings['g_jet_raa_sys_corr_ALICE'].SetFillColorAlpha(self.data_color, 0.3)
            self.observable_settings['g_jet_raa_sys_corr_ALICE'].SetMarkerStyle(0)
            self.observable_settings['g_jet_raa_sys_corr_ALICE'].SetFillStyle(0)
            self.observable_settings['g_jet_raa_sys_corr_ALICE'].SetLineWidth(1)
            self.observable_settings['g_jet_raa_sys_corr_ALICE'].SetLineColor(self.data_color)
            
            self.observable_settings['g_jet_raa_ALICE'].SetMarkerColor(self.data_color)
            self.observable_settings['g_jet_raa_ALICE'].SetLineColor(self.data_color)
            self.observable_settings['g_jet_raa_ALICE'].SetLineWidth(1)
            self.observable_settings['g_jet_raa_ALICE'].SetMarkerSize(1.3)
            
            self.observable_settings['g_jet_raa_ALICE'].SetMarkerStyle(21)
            self.observable_settings['g_jet_raa_ALICE'].Draw('PE Z same')
            self.observable_settings['g_jet_raa_sys_corr_ALICE'].Draw('E2 same')
            self.observable_settings['g_jet_raa_sys_shape_ALICE'].Draw('E2 same')

            leg.AddEntry(self.observable_settings['g_jet_raa_ALICE'], '0-10%','PE')
            leg.AddEntry(self.observable_settings['g_jet_raa_sys_corr_ALICE'], 'Corr. uncertainty','f')
            leg.AddEntry(self.observable_settings['g_jet_raa_sys_shape_ALICE'], 'Shape uncertainty','f')
            
        elif self.observable in ['hjet_IAA', 'hjet_IAA_ratio']:
        
            self.observable_settings['g_jet_hjet_IAA_sys_shape_ALICE'].SetFillColor(self.data_color)
            self.observable_settings['g_jet_hjet_IAA_sys_shape_ALICE'].SetFillColorAlpha(self.data_color, 0.3)
            self.observable_settings['g_jet_hjet_IAA_sys_shape_ALICE'].SetMarkerStyle(0)
            self.observable_settings['g_jet_hjet_IAA_sys_shape_ALICE'].SetFillStyle(1001)
            self.observable_settings['g_jet_hjet_IAA_sys_shape_ALICE'].SetLineWidth(0)
            self.observable_settings['g_jet_hjet_IAA_sys_shape_ALICE'].SetLineColor(0)
            
            self.observable_settings['g_jet_hjet_IAA_sys_corr_ALICE'].SetFillColor(self.data_color)
            self.observable_settings['g_jet_hjet_IAA_sys_corr_ALICE'].SetFillColorAlpha(self.data_color, 0.3)
            self.observable_settings['g_jet_hjet_IAA_sys_corr_ALICE'].SetMarkerStyle(0)
            self.observable_settings['g_jet_hjet_IAA_sys_corr_ALICE'].SetFillStyle(0)
            self.observable_settings['g_jet_hjet_IAA_sys_corr_ALICE'].SetLineWidth(1)
            self.observable_settings['g_jet_hjet_IAA_sys_corr_ALICE'].SetLineColor(self.data_color)
            
            self.observable_settings['g_jet_hjet_IAA_ALICE'].SetMarkerColor(self.data_color)
            self.observable_settings['g_jet_hjet_IAA_ALICE'].SetLineColor(self.data_color)
            self.observable_settings['g_jet_hjet_IAA_ALICE'].SetLineWidth(1)
            self.observable_settings['g_jet_hjet_IAA_ALICE'].SetMarkerSize(1.3)
            
            self.observable_settings['g_jet_hjet_IAA_ALICE'].SetMarkerStyle(21)
            self.observable_settings['g_jet_hjet_IAA_ALICE'].Draw('PE Z same')
            self.observable_settings['g_jet_hjet_IAA_sys_corr_ALICE'].Draw('E2 same')
            self.observable_settings['g_jet_hjet_IAA_sys_shape_ALICE'].Draw('E2 same')

            leg.AddEntry(self.observable_settings['g_jet_hjet_IAA_ALICE'], '0-10%','PE')
            leg.AddEntry(self.observable_settings['g_jet_hjet_IAA_sys_corr_ALICE'], 'Corr. uncertainty','f')
            leg.AddEntry(self.observable_settings['g_jet_hjet_IAA_sys_shape_ALICE'], 'Shape uncertainty','f')
            
        leg.Draw('same')
        leg2.Draw('same')

        if self.observable in ['combined_raa', 'hadron_raa', 'hjet_IAA', 'hjet_IAA_ratio']:
            line = ROOT.TLine(myBlankHisto.GetXaxis().GetXmin(),1,myBlankHisto.GetXaxis().GetXmax(),1)
            line.SetLineColor(1)
            line.SetLineStyle(2)
            line.Draw('same')

        # # # # # # # # # # # # # # # # # # # # # # # #
        # text
        # # # # # # # # # # # # # # # # # # # # # # # #
        x = 0.2
        dx = 0.4
        y = 0.87
        dy = 0.06
        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        text_latex.SetTextSize(0.05)
        text = f'#bf{{ALICE}} #sqrt{{#it{{s_{{#it{{NN}}}}}}}} = {self.sqrts/1000.} TeV'
        text_latex.DrawLatex(x, y, text)

        if self.observable == 'combined_raa':
            text = f'Charged particles'
            text_latex.DrawLatex(x, y-dy, text)
            text = '0-5%'
            text_latex.DrawLatex(x, y-2*dy, text)
            text = f'Jets, anti-#it{{k}}_{{T}}, #it{{R}} = {self.jetR}'
            text_latex.DrawLatex(x+dx, y-dy, text)
            text = '0-10%'
            text_latex.DrawLatex(x+dx, y-2*dy, text)
        
        elif self.observable == 'hadron_raa':
            text = f'Charged particles  |#eta| < {self.eta_cut}'
            text_latex.DrawLatex(x,y-dy, text)
        elif self.observable == 'jet_raa':
            text = f'#it{{R}} = {self.jetR}, anti-#it{{k}}_{{T}}, |#it{{#eta}}_{{jet}}| < {self.eta_cut}'
            text_latex.DrawLatex(x, y-dy, text)
        elif self.observable in ['hjet_IAA']:
            text = f'#it{{R}} = {self.jetR}, anti-#it{{k}}_{{T}}, |#it{{#eta}}_{{jet}}| < {self.eta_cut}, #it{{#pi}} - #Delta#it{{#varphi}} < 0.6'
            text_latex.DrawLatex(x, y-dy, text)
            text = 'TT{20,50} - TT{8,9}'
            text_latex.DrawLatex(x, y-2.2*dy, text)
        elif self.observable in ['hjet_IAA_ratio']:
            text = f'anti-#it{{k}}_{{T}}, |#it{{#eta}}_{{jet}}| < 0.9 - #it{{R}}, #it{{#pi}} - #Delta#it{{#varphi}} < 0.6'
            text_latex.DrawLatex(x, y-dy, text)
            text = 'TT{20,50} - TT{8,9}'
            text_latex.DrawLatex(x, y-2.2*dy, text)

        output_filename = os.path.join(self.output_dir, f'h_{self.observable}.pdf')
        c.SaveAs(output_filename)
        c.Close()

    #-------------------------------------------------------------------------------------------
    # Plot distribution
    #-------------------------------------------------------------------------------------------
    def plot_distribution(self, pp_label=None, logy = False):
    
        cname = 'c'
        c = ROOT.TCanvas(cname, cname, 600, 450)
        c.SetRightMargin(0.05)
        if self.observable == 'mass':
            c.SetLeftMargin(0.17)
        else:
            c.SetLeftMargin(0.15)
        c.SetTopMargin(0.05)
        c.SetBottomMargin(0.17)
        c.SetTicks(0,1)
        if self.observable == 'combined_raa':
            c.SetLogx()
        c.cd()
        
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, self.x_min, self.x_max)
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(self.xtitle)
        myBlankHisto.SetYTitle(self.ytitle)
        myBlankHisto.SetMaximum(self.y_ratio_max)
        myBlankHisto.SetMinimum(self.y_ratio_min)
        myBlankHisto.GetXaxis().SetTitleSize(0.06)
        myBlankHisto.GetXaxis().SetTitleOffset(1.2)
        myBlankHisto.GetXaxis().SetLabelSize(0.05)
        myBlankHisto.GetYaxis().SetTitleSize(0.07)
        if self.observable == 'mass':
            myBlankHisto.GetYaxis().SetTitleOffset(1.1)
        else:
            myBlankHisto.GetYaxis().SetTitleOffset(1.)
        myBlankHisto.GetYaxis().SetLabelSize(0.05)
        myBlankHisto.SetLineWidth(0)
        myBlankHisto.Draw()
        
        if self.observable in ['hjet_dphi']:
            y_leg = 0.52
            leg = ROOT.TLegend(0.2,y_leg,0.4,y_leg+0.11)
        if self.observable in ['mass']:
            y_leg = 0.73
            leg = ROOT.TLegend(0.65,y_leg,0.8,y_leg+0.13)
        self.utils.setup_legend(leg, 0.04, sep=-0.1)
        
        if self.observable == 'hjet_dphi':
            y_leg = 0.32
            leg2 = ROOT.TLegend(0.2,y_leg,0.5,y_leg+0.16)
        if self.observable == 'mass':
            y_leg = 0.47
            leg2 = ROOT.TLegend(0.51,y_leg,0.7,y_leg+0.2)
        self.utils.setup_legend(leg2, 0.04, sep=-0.1)
        
        # Draw theory predictions
        for i, prediction in enumerate(self.observable_settings['prediction_distribution_list']):

            label = self.observable_settings['prediction_distribution_labels'][i]
            print(label)
            color = self.theory_colors[i]
           
            prediction.SetFillColor(color)
            prediction.SetFillColorAlpha(color, self.alpha)
            prediction.SetLineColor(color)
            prediction.SetLineColorAlpha(color, self.alpha)
            prediction.SetLineWidth(8)
            prediction.SetMarkerSize(0)
            prediction.SetMarkerStyle(0)
           
            if self.observable == 'hjet_dphi':
                prediction.SetLineColor(0)
                prediction.SetFillColorAlpha(color, 1.)
                if 'pp' in label:
                    prediction.SetFillColor(ROOT.kViolet-8)
                else:
                    prediction.SetFillColor(ROOT.kRed-7)
                prediction.SetLineWidth(8)
                if 'scaled' in label:
                    prediction.SetFillStyle(3001)
              
            if type(prediction) in [ROOT.TH1D, ROOT.TH1F]:
                prediction.Draw('E3 same')
                if self.observable == 'hjet_dphi':
                    leg2.AddEntry(prediction, label, 'F')
                else:
                    leg2.AddEntry(prediction, label, 'L')

            elif type(prediction) in [ROOT.TGraph, ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                 
                if type(prediction) in [ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                    prediction.Draw('3 same')
                elif type(prediction) == ROOT.TGraph:
                    prediction.Draw('same')

                leg2.AddEntry(prediction, label, 'L')
                
        # Draw experimental data
        if self.observable == 'hjet_dphi':
        
            self.observable_settings['g_pp'].SetMarkerSize(self.marker_size)
            self.observable_settings['g_pp'].SetMarkerStyle(self.data_markers[0])
            self.observable_settings['g_pp'].SetMarkerColor(self.data_color)
            self.observable_settings['g_pp'].SetLineStyle(self.line_style)
            self.observable_settings['g_pp'].SetLineWidth(self.line_width)
            self.observable_settings['g_pp'].SetLineColor(self.data_color)
        
            self.observable_settings['g_AA'].SetMarkerSize(self.marker_size)
            self.observable_settings['g_AA'].SetMarkerStyle(self.data_markers[1])
            self.observable_settings['g_AA'].SetMarkerColor(self.data_color)
            self.observable_settings['g_AA'].SetLineStyle(self.line_style)
            self.observable_settings['g_AA'].SetLineWidth(self.line_width)
            self.observable_settings['g_AA'].SetLineColor(self.data_color)

            leg.AddEntry(self.observable_settings['g_AA'], 'Pb-Pb 0-10%','PE')
            leg.AddEntry(self.observable_settings['g_pp'], 'PYTHIA + Pb-Pb','PE')

            self.observable_settings['g_pp'].Draw('PE Z same')
            self.observable_settings['g_AA'].Draw('PE Z same')

        else:

            # pp distribution
            if 'h_pp_sys' in self.observable_settings:
                self.observable_settings['h_pp_sys'].SetLineColor(0)
                self.observable_settings['h_pp_sys'].SetFillColor(self.data_color)
                #self.observable_settings['h_pp_sys'].SetFillColorAlpha(self.data_color, 0.3)
                self.observable_settings['h_pp_sys'].SetFillStyle(3004)
                self.observable_settings['h_pp_sys'].SetLineWidth(0)
                self.observable_settings['h_pp_sys'].SetMarkerStyle(0)
            
            if 'g_pp' in self.observable_settings:
                self.observable_settings['g_pp'].SetMarkerSize(self.marker_size)
                self.observable_settings['g_pp'].SetMarkerStyle(self.data_markers[0])
                self.observable_settings['g_pp'].SetMarkerColor(self.data_color)
                self.observable_settings['g_pp'].SetLineStyle(self.line_style)
                self.observable_settings['g_pp'].SetLineWidth(self.line_width)
                self.observable_settings['g_pp'].SetLineColor(self.data_color)
            elif 'h_pp' in self.observable_settings:
                self.observable_settings['h_pp'].SetMarkerSize(self.marker_size)
                self.observable_settings['h_pp'].SetMarkerStyle(self.data_markers[0])
                self.observable_settings['h_pp'].SetMarkerColor(self.data_color)
                self.observable_settings['h_pp'].SetLineStyle(self.line_style)
                self.observable_settings['h_pp'].SetLineWidth(self.line_width)
                self.observable_settings['h_pp'].SetLineColor(self.data_color)
            
            # AA distribution
            self.observable_settings['h_AA_stat'].SetMarkerSize(self.marker_size)
            self.observable_settings['h_AA_stat'].SetMarkerStyle(self.data_markers[1])
            self.observable_settings['h_AA_stat'].SetMarkerColor(self.data_color)
            self.observable_settings['h_AA_stat'].SetLineStyle(self.line_style)
            self.observable_settings['h_AA_stat'].SetLineWidth(self.line_width)
            self.observable_settings['h_AA_stat'].SetLineColor(self.data_color)
            
            if 'g_AA_sys' in self.observable_settings:
                self.observable_settings['g_AA_sys'].SetLineColor(0)
                self.observable_settings['g_AA_sys'].SetFillColor(self.data_color)
                self.observable_settings['g_AA_sys'].SetFillColorAlpha(self.data_color, 0.3)
                self.observable_settings['g_AA_sys'].SetFillStyle(1001)
                self.observable_settings['g_AA_sys'].SetLineWidth(0)
            if 'h_AA_sys' in self.observable_settings:
                self.observable_settings['h_AA_sys'].SetLineColor(0)
                self.observable_settings['h_AA_sys'].SetFillColor(self.data_color)
                self.observable_settings['h_AA_sys'].SetFillColorAlpha(self.data_color, 0.3)
                self.observable_settings['h_AA_sys'].SetFillStyle(1001)
                self.observable_settings['h_AA_sys'].SetLineWidth(0)

            # Add legend
            if not pp_label:
                pp_label = 'pp'
            if 'g_pp' in self.observable_settings:
                leg.AddEntry(self.observable_settings['g_pp'], pp_label, 'P')
            else:
                leg.AddEntry(self.observable_settings['h_pp'], pp_label, 'PE')
            leg.AddEntry(self.observable_settings['h_AA_stat'], 'Pb#font[122]{{-}}Pb {}#font[122]{{-}}{}%'.format(int(self.centrality[0]), int(self.centrality[1])), 'PE')
            if 'g_AA_sys' in self.observable_settings:
                leg.AddEntry(self.observable_settings['g_AA_sys'], 'Sys. uncertainty', 'f')
            elif 'h_AA_sys' in self.observable_settings:
                leg.AddEntry(self.observable_settings['h_AA_sys'], 'Sys. uncertainty', 'f')

            # Draw distribution
            if 'h_pp_sys' in self.observable_settings:
                self.observable_settings['h_pp_sys'].Draw('E2 same')
            if 'g_pp' in self.observable_settings:
                self.observable_settings['g_pp'].Draw('P same')
            elif 'h_pp' in self.observable_settings:
                self.observable_settings['h_pp'].Draw('PE X0 same')
                
            if 'g_AA_sys' in self.observable_settings:
                self.observable_settings['g_AA_sys'].Draw('E2 same')
            elif 'h_AA_sys' in self.observable_settings:
                self.observable_settings['h_AA_sys'].Draw('E2 same')
            self.observable_settings['h_AA_stat'].DrawCopy('PE X0 same')


        leg.Draw('same')
        leg2.Draw('same')
        
        # # # # # # # # # # # # # # # # # # # # # # # #
        # text
        # # # # # # # # # # # # # # # # # # # # # # # #
        x = 0.2
        dx = 0.4
        if self.observable == 'mass':
            y = 0.89
        else:
            y = 0.87
        dy = 0.06
        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        text_latex.SetTextSize(0.05)
        text = f'#bf{{ALICE}} #sqrt{{#it{{s_{{#it{{NN}}}}}}}} = {self.sqrts/1000.} TeV'
        text_latex.DrawLatex(x, y, text)

        if self.observable in ['hjet_dphi']:
            text = f'#it{{R}} = {self.jetR}, anti-#it{{k}}_{{T}}, |#it{{#eta}}_{{jet}}| < {self.eta_cut}, #it{{#pi}} - #Delta#it{{#varphi}} < 0.6'
            text_latex.DrawLatex(x, y-dy, text)
            text = '40 < #it{p}_{T,jet}^{ch} < 60 GeV/#it{c}'
            text_latex.DrawLatex(x, y-2.2*dy, text)
            text = 'TT{20,50} - TT{8,9}'
            text_latex.DrawLatex(x, y-3.4*dy, text)

        if self.observable in ['mass']:
            text = 'Charged-particle jets'
            text_latex.DrawLatex(x, y-dy, text)
            text = f'#it{{R}} = {self.jetR}, anti-#it{{k}}_{{T}}, |#it{{#eta}}_{{jet}}| < {0.9-self.jetR}'
            text_latex.DrawLatex(x, y-2.*dy, text)
            pt_label = f'{self.pt[0]} < #it{{p}}_{{T, ch jet}} < {self.pt[1]} GeV/#it{{c}}'
            text_latex.DrawLatex(x, y-3.3*dy, pt_label)

        output_filename = os.path.join(self.output_dir, f'h_{self.observable}.pdf')
        c.SaveAs(output_filename)
        c.Close()
        
    #-------------------------------------------------------------------------------------------
    # Plot distributions in upper panel, and ratio in lower panel
    #-------------------------------------------------------------------------------------------
    def plot_distribution_and_ratio(self, plot_ratio=True, pp_label=None, logy = False):
    
        # Output filename
        output_filename = os.path.join(self.output_dir, f'h_{self.observable}.pdf')
    
        c = ROOT.TCanvas('c', 'c', 600, 650)
        c.Draw()
               
        c.cd()
        pad2_dy = 0.45
        pad1 = ROOT.TPad('myPad', 'The pad',0,pad2_dy,1,1)
        pad1.SetLeftMargin(0.2)
        pad1.SetTopMargin(0.08)
        pad1.SetRightMargin(0.04)
        pad1.SetBottomMargin(0.)
        pad1.SetTicks(0,1)
        pad1.Draw()
        if logy:
            pad1.SetLogy()
        pad1.cd()
        
        if self.observable in ['theta_g']:
            myLegend = ROOT.TLegend(0.24,0.65,0.4,0.88)
        else:
            myLegend = ROOT.TLegend(0.68,0.65,0.85,0.88)
        self.utils.setup_legend(myLegend, 0.055, sep=-0.1)
        
        # pp distribution
        if 'h_pp_sys' in self.observable_settings:
            self.observable_settings['h_pp_sys'].SetLineColor(0)
            self.observable_settings['h_pp_sys'].SetFillColor(self.data_color)
            #self.observable_settings['h_pp_sys'].SetFillColorAlpha(self.data_color, 0.3)
            self.observable_settings['h_pp_sys'].SetFillStyle(3004)
            self.observable_settings['h_pp_sys'].SetLineWidth(0)
            self.observable_settings['h_pp_sys'].SetMarkerStyle(0)
        
        if 'g_pp' in self.observable_settings:
            self.observable_settings['g_pp'].SetMarkerSize(self.marker_size)
            self.observable_settings['g_pp'].SetMarkerStyle(self.data_markers[0])
            self.observable_settings['g_pp'].SetMarkerColor(self.data_color)
            self.observable_settings['g_pp'].SetLineStyle(self.line_style)
            self.observable_settings['g_pp'].SetLineWidth(self.line_width)
            self.observable_settings['g_pp'].SetLineColor(self.data_color)
        elif 'h_pp' in self.observable_settings:
            self.observable_settings['h_pp'].SetMarkerSize(self.marker_size)
            self.observable_settings['h_pp'].SetMarkerStyle(self.data_markers[0])
            self.observable_settings['h_pp'].SetMarkerColor(self.data_color)
            self.observable_settings['h_pp'].SetLineStyle(self.line_style)
            self.observable_settings['h_pp'].SetLineWidth(self.line_width)
            self.observable_settings['h_pp'].SetLineColor(self.data_color)
        
        # AA distribution
        self.observable_settings['h_AA_stat'].SetMarkerSize(self.marker_size)
        self.observable_settings['h_AA_stat'].SetMarkerStyle(self.data_markers[1])
        self.observable_settings['h_AA_stat'].SetMarkerColor(self.data_color)
        self.observable_settings['h_AA_stat'].SetLineStyle(self.line_style)
        self.observable_settings['h_AA_stat'].SetLineWidth(self.line_width)
        self.observable_settings['h_AA_stat'].SetLineColor(self.data_color)
          
        if 'g_AA_sys' in self.observable_settings:
            self.observable_settings['g_AA_sys'].SetLineColor(0)
            self.observable_settings['g_AA_sys'].SetFillColor(self.data_color)
            self.observable_settings['g_AA_sys'].SetFillColorAlpha(self.data_color, 0.3)
            self.observable_settings['g_AA_sys'].SetFillStyle(1001)
            self.observable_settings['g_AA_sys'].SetLineWidth(0)
        if 'h_AA_sys' in self.observable_settings:
            self.observable_settings['h_AA_sys'].SetLineColor(0)
            self.observable_settings['h_AA_sys'].SetFillColor(self.data_color)
            self.observable_settings['h_AA_sys'].SetFillColorAlpha(self.data_color, 0.3)
            self.observable_settings['h_AA_sys'].SetFillStyle(1001)
            self.observable_settings['h_AA_sys'].SetLineWidth(0)
             
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, self.bins[0], self.bins[-1])
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(self.xtitle)
        myBlankHisto.SetYTitle(self.ytitle)
        myBlankHisto.SetMaximum(self.y_max)
        myBlankHisto.SetMinimum(2e-4) # Don't draw 0 on top panel
        myBlankHisto.GetYaxis().SetTitleSize(0.09)
        myBlankHisto.GetYaxis().SetTitleOffset(1.)
        myBlankHisto.GetYaxis().SetLabelSize(0.07)
        myBlankHisto.Draw('E')

        # Ratio
        c.cd()
        pad2 = ROOT.TPad('pad2', 'pad2', 0, 0.02, 1, pad2_dy)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.21)
        pad2.SetLeftMargin(0.2)
        pad2.SetRightMargin(0.04)
        pad2.SetTicks(0,1)
        pad2.Draw()
        pad2.cd()
              
        myBlankHisto2 = myBlankHisto.Clone('myBlankHisto_C')
        myBlankHisto2.SetYTitle('#frac{Pb#font[122]{-}Pb}{pp}')
        myBlankHisto2.SetXTitle(self.xtitle)
        myBlankHisto2.GetXaxis().SetTitleSize(30)
        myBlankHisto2.GetXaxis().SetTitleFont(43)
        myBlankHisto2.GetXaxis().SetTitleOffset(1.95)
        myBlankHisto2.GetXaxis().SetLabelFont(43)
        myBlankHisto2.GetXaxis().SetLabelSize(25)
        myBlankHisto2.GetYaxis().SetTitleSize(28)
        myBlankHisto2.GetYaxis().SetTitleFont(43)
        myBlankHisto2.GetYaxis().SetTitleOffset(2.)
        myBlankHisto2.GetYaxis().SetLabelFont(43)
        myBlankHisto2.GetYaxis().SetLabelSize(25)
        myBlankHisto2.GetYaxis().SetNdivisions(505)
        
        myBlankHisto2.GetYaxis().SetRangeUser(self.y_ratio_min, self.y_ratio_max)

        if self.observable == 'girth':
            ratio_legend = ROOT.TLegend(0.45,0.65,0.6,0.95)
        elif self.observable == 'mass':
            ratio_legend = ROOT.TLegend(0.25,0.65,0.4,0.95)
        if self.observable == 'zg':
            ratio_legend = ROOT.TLegend(0.29,0.75,0.4,0.95)
        if self.observable == 'theta_g':
            ratio_legend = ROOT.TLegend(0.34,0.72,0.52,0.96)
        self.utils.setup_legend(ratio_legend, 0.07, sep=0.2)
        
        if self.observable == 'zg':
            ratio_legend2 = ROOT.TLegend(0.5,0.72,0.7,0.975)
            self.utils.setup_legend(ratio_legend2, 0.07, sep=0.2)
        elif self.observable == 'theta_g':
            ratio_legend2 = ROOT.TLegend(0.6,0.72,0.8,0.96)
            self.utils.setup_legend(ratio_legend2, 0.07, sep=0.2)

        myBlankHisto2.Draw('')
        
        if plot_ratio:
            self.observable_settings['h_ratio_stat'].SetMarkerSize(1.5)
            self.observable_settings['h_ratio_stat'].SetMarkerStyle(self.data_markers[1])
            self.observable_settings['h_ratio_stat'].SetMarkerColor(self.data_color)
            self.observable_settings['h_ratio_stat'].SetLineColor(self.data_color)
            self.observable_settings['h_ratio_stat'].SetLineWidth(self.line_width)

            if 'g_ratio_sys' in self.observable_settings:
                self.observable_settings['g_ratio_sys'].SetFillColor(self.data_color)
                self.observable_settings['g_ratio_sys'].SetFillColorAlpha(self.data_color, 0.3)
            elif 'h_ratio_sys' in self.observable_settings:
                self.observable_settings['h_ratio_sys'].SetFillColor(self.data_color)
                self.observable_settings['h_ratio_sys'].SetFillColorAlpha(self.data_color, 0.3)

        pad1.cd()
        
        # Draw theory predictions for distribution
        for i, prediction in enumerate(self.observable_settings['prediction_distribution_list']):
        
            label = self.observable_settings['prediction_distribution_labels'][i]
            color = self.theory_colors[i]
            
            prediction.SetFillColor(color)
            prediction.SetFillColorAlpha(color, self.alpha)
            prediction.SetLineColor(color)
            prediction.SetLineWidth(8)
            prediction.SetMarkerSize(0)
            prediction.SetMarkerStyle(0)
        
            if type(prediction) in [ROOT.TH1F]:
            
                prediction.Draw('E3 same')
                ratio_legend.AddEntry(prediction, label, 'L')

            elif type(prediction) in [ROOT.TGraph, ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                 
                if type(prediction) in [ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                    prediction.Draw('3 same')
                elif type(prediction) == ROOT.TGraph:
                    prediction.Draw('same')

                ratio_legend.AddEntry(prediction, label, 'L')
           
        # Draw distribution
        if 'h_pp_sys' in self.observable_settings:
            self.observable_settings['h_pp_sys'].Draw('E2 same')
        if 'g_pp' in self.observable_settings:
            self.observable_settings['g_pp'].Draw('P same')
        elif 'h_pp' in self.observable_settings:
            self.observable_settings['h_pp'].Draw('PE X0 same')
            
        if 'g_AA_sys' in self.observable_settings:
            self.observable_settings['g_AA_sys'].Draw('E2 same')
        elif 'h_AA_sys' in self.observable_settings:
            self.observable_settings['h_AA_sys'].Draw('E2 same')
        self.observable_settings['h_AA_stat'].DrawCopy('PE X0 same')
        
        pad2.cd()

        # Draw theory predictions for ratio
        for i, prediction in enumerate(self.observable_settings['prediction_ratio_list']):
        
            label = self.observable_settings['prediction_ratio_labels'][i]
            color = self.theory_colors[i]
            
            prediction.SetFillColor(color)
            prediction.SetFillColorAlpha(color, self.alpha)
            prediction.SetLineColor(color)
            prediction.SetLineWidth(8)
            prediction.SetMarkerSize(0)
            prediction.SetMarkerStyle(0)
        
            if type(prediction) in [ROOT.TH1F]:
            
                prediction.Draw('E3 same')
                if plot_ratio:
                    ratio_legend.AddEntry(prediction, label, 'L')

            elif type(prediction) in [ROOT.TGraph, ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                 
                if type(prediction) in [ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                    prediction.Draw('3 same')
                elif type(prediction) == ROOT.TGraph:
                    prediction.Draw('same')

                if self.observable in ['zg', 'theta_g'] and 'Pablos' in label:
                    ratio_legend2.AddEntry(prediction, label, 'L')
                else:
                    ratio_legend.AddEntry(prediction, label, 'L')
                    
        # Re-draw some curves
        if self.observable in ['girth']:
            self.observable_settings['prediction_ratio_list'][0].Draw('E3 same')
          
        ratio_legend.Draw()
        if self.observable in ['zg', 'theta_g']:
            ratio_legend2.Draw()
        line = ROOT.TLine(self.bins[0], 1, self.bins[-1], 1)
        line.SetLineColor(920+2)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
        
        if plot_ratio:
            if 'g_ratio_sys' in self.observable_settings:
                self.observable_settings['g_ratio_sys'].Draw('E2 same')
            if 'h_ratio_sys' in self.observable_settings:
                self.observable_settings['h_ratio_sys'].Draw('E2 same')
            self.observable_settings['h_ratio_stat'].Draw('PE X0 same')
       
        pad1.cd()
        if not pp_label:
            pp_label = 'pp'
        if 'g_pp' in self.observable_settings:
            myLegend.AddEntry(self.observable_settings['g_pp'], pp_label, 'P')
        else:
            myLegend.AddEntry(self.observable_settings['h_pp'], pp_label, 'PE')
        myLegend.AddEntry(self.observable_settings['h_AA_stat'], 'Pb#font[122]{{-}}Pb {}#font[122]{{-}}{}%'.format(int(self.centrality[0]), int(self.centrality[1])), 'PE')
        if 'g_AA_sys' in self.observable_settings:
            myLegend.AddEntry(self.observable_settings['g_AA_sys'], 'Sys. uncertainty', 'f')
        elif 'h_AA_sys' in self.observable_settings:
            myLegend.AddEntry(self.observable_settings['h_AA_sys'], 'Sys. uncertainty', 'f')
       
        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        
        if self.observable in ['theta_g']:
            x = 0.52
        else:
            x = 0.25
        text_latex.SetTextSize(0.065)
        text = f'#bf{{ALICE}} #sqrt{{#it{{s_{{#it{{NN}}}}}}}} = {self.sqrts/1000.} TeV'
        text_latex.DrawLatex(x, 0.83, text)

        text = 'Charged-particle jets'
        text_latex.DrawLatex(x, 0.75, text)
        
        text = f'#it{{R}} = {self.jetR}, anti-#it{{k}}_{{T}}, |#it{{#eta}}_{{jet}}| < {0.9-self.jetR}'
        text_latex.DrawLatex(x, 0.67, text)
        
        # pt bin label
        pt_label = None
        if self.observable in ['girth', 'mass', 'zg', 'theta_g']:
            pt_label = f'{self.pt[0]} < #it{{p}}_{{T, ch jet}} < {self.pt[1]} GeV/#it{{c}}'
            text_latex.DrawLatex(x, 0.57, pt_label)
            
        if self.observable in ['zg', 'theta_g']:
            grooming_label = f'Soft Drop #it{{z}}_{{cut}}={self.zcut}, #it{{#beta}}=0'
            text_latex.DrawLatex(x, 0.47, grooming_label)
        
        myLegend.Draw()

        c.SaveAs(output_filename)
        c.Close()
 
    #---------------------------------------------------------------
    # Rebin numpy arrays (xbins,y) representing a histogram,
    # where x represents the bin edges, and y the bin content
    #
    # Return (x_rebinned, y_rebinned) where x_rebinned is the bin centers
    #
    # I don't find any easy way to do this...so I
    #        construct TH1, rebin, and then extract array)
    #----------------------------------------------------------------------
    def rebin_arrays(self, x, y, yerr):

        x_array = array('d', x)
        h = ROOT.TH1F('h', 'h', len(x)-1, x_array)
        for i, content in enumerate(y):
            h.SetBinContent(i+1, y[i])
            h.SetBinError(i+1, yerr[i])
        h_rebinned = h.Rebin(len(self.bin_array)-1, 'h_rebinned', self.bin_array)
        x_rebinned = []
        y_rebinned = []
        y_err_rebinned = []
        for i in range(len(self.bin_array)-1):
            x_rebinned.append(h_rebinned.GetBinCenter(i+1))
            y_rebinned.append(h_rebinned.GetBinContent(i+1))
            y_err_rebinned.append(h_rebinned.GetBinError(i+1))
          
        return (np.array(x_rebinned), np.array(y_rebinned), np.array(y_err_rebinned))
  
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
if __name__ == '__main__':
    print('Executing plot_results_TG3.py...')
    print('')

    # Define arguments
    parser = argparse.ArgumentParser(description='Plot TG3 observables')
    parser.add_argument(
        '-c',
        '--configFile',
        action='store',
        type=str,
        metavar='configFile',
        default='TG3.yaml',
        help='Config file'
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

    analysis = PlotResults(config_file=args.configFile, output_dir=args.outputDir)
    analysis.plot_results()
