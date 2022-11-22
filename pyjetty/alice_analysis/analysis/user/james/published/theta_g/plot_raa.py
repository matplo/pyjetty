#! /usr/bin/env python

import sys
import os
from array import *
import numpy as np
import ROOT
import yaml
import hepdata_lib

from pyjetty.alice_analysis.process.base import common_base
from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class PlotRAA(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, **kwargs):
        super(PlotRAA, self).__init__(**kwargs)

        self.utils = analysis_utils_obs.AnalysisUtils_Obs()
            
        self.initialize()

    #---------------------------------------------------------------
    # Initialize config file into class members
    #---------------------------------------------------------------
    def initialize(self):

        self.output_dir = '.'

        # Load config file
        config_file = './plot.yaml'
        with open(config_file, 'r') as stream:
            config = yaml.safe_load(stream)

        self.config_james = config['james']
        self.config_james_results = [self.config_james[name] for name in self.config_james if 'result' in name ]

        self.config_laura = config['laura']
        self.config_laura_results = [self.config_laura[name] for name in self.config_laura if 'result' in name ]

        self.plot_data = True
        self.plot_theory = True
        self.plot_axis_rg = False

        self.colors = [600-6, 632-4]
        self.ratio_color = ROOT.kGray+3
        self.markers = [20, 21]
        ROOT.gStyle.SetLineStyleString(11,'30 12');
        
        self.theory_colors = {}
        self.theory_colors['JETSCAPE'] = ROOT.kViolet-8
        self.theory_colors['JEWEL, recoils off'] = ROOT.kCyan-2
        self.theory_colors['JEWEL, recoils on'] = ROOT.kGreen-6
        self.theory_colors['Pablos, #it{L}_{res} = 0'] = ROOT.kTeal-8
        self.theory_colors['Pablos, #it{L}_{res} = 2/#pi#it{T}'] = ROOT.kOrange+6
        self.theory_colors['Pablos, #it{L}_{res} = #infty'] = ROOT.kRed-7
        self.theory_colors['Caucal'] = ROOT.kAzure-4
        self.theory_colors['Qin'] = ROOT.kCyan-2
        self.theory_colors['Chien'] = ROOT.kBlue-10
        self.theory_colors['Yuan, #hat{q}L = 5 GeV^{2}'] = ROOT.kPink+1
        self.theory_colors['Yuan, med #it{q}/#it{g}'] = ROOT.kBlue-10
        self.theory_colors['Yuan, quark'] = ROOT.kGray

        self.line_style = {}
        self.line_style['Qin'] = 11
        self.line_style['Yuan, #hat{q}L = 5 GeV^{2}'] = 1
        self.line_style['Yuan, med #it{q}/#it{g}'] = 11
        self.line_style['Yuan, quark'] = 11

        self.formatted_grooming_label = 'Soft Drop #it{{z}}_{{cut}}={}, #it{{#beta}}=0'
        self.figure_approval_status = ''
        
        self.debug_level = 0

    #---------------------------------------------------------------
    # This function is called once for each subconfiguration
    #----------------------------------------------------------------------
    def plot_raa(self):
    
        # Create HEPData submission
        self.hepdata_dir = os.path.join(self.output_dir, 'hepdata')
        self.hepdata_submission = hepdata_lib.Submission()

        for result in self.config_james_results:
            self.centrality = [0, 10]
            self.init_result_james(result)
            self.plot_result()

        for result in self.config_laura_results:
            self.centrality = [30, 50]
            self.init_result_laura(result)
            self.plot_result()
            
        # Write submission files
        self.hepdata_submission.create_files(self.hepdata_dir)

    #---------------------------------------------------------------
    # Initialize
    #---------------------------------------------------------------
    def init_result_james(self, result):
    
        # Init from config file
        self.observable = result['observable']
        self.jetR = result['jetR']
        self.zcut = result['zcut']
        self.obs_label = result['obs_label']
        self.min_pt = result['min_pt']
        self.max_pt = result['max_pt']
        self.file_pp_name = result['file_pp']
        self.file_AA_name = result['file_AA']
        self.set_xy_titles()

        # Get hists from ROOT file
        self.file_pp = ROOT.TFile(self.file_pp_name, 'READ')
        self.file_AA = ROOT.TFile(self.file_AA_name, 'READ')
        
        self.main_result_name = 'hmain_{}_R{}_{}_{}-{}'.format(self.observable, self.jetR, self.obs_label, self.min_pt, self.max_pt)
        self.sys_total_name = 'hResult_{}_systotal_R{}_{}_{}-{}'.format(self.observable, self.jetR, self.obs_label, self.min_pt, self.max_pt)
        self.h_tagging_fraction_name = 'h_tagging_fraction_R{}_{}'.format(self.jetR, self.obs_label)

        self.h_main_pp = self.file_pp.Get(self.main_result_name)
        self.h_sys_pp = self.file_pp.Get(self.sys_total_name)
        self.h_tagging_fraction_pp = self.file_pp.Get(self.h_tagging_fraction_name)
        self.h_main_AA = self.file_AA.Get(self.main_result_name)
        self.h_sys_AA = self.file_AA.Get(self.sys_total_name)
        self.h_tagging_fraction_AA = self.file_AA.Get(self.h_tagging_fraction_name)

        # Add JEWEL systematic, if applicable
        if 'jewel_systematic' in result:
            jewel_filename = result['jewel_systematic']
            self.file_jewel = ROOT.TFile(jewel_filename, 'READ')
            if self.observable == 'zg':
                sys_jewel_name = f'hSystematic_{self.observable}_generator_R{self.jetR}_{self.obs_label}_n6_{self.min_pt}-{self.max_pt}'
            else:
                sys_jewel_name = f'hSystematic_{self.observable}_generator_R{self.jetR}_{self.obs_label}_n3_{self.min_pt}-{self.max_pt}'
            self.h_sys_jewel = self.file_jewel.Get(sys_jewel_name)
            for bin in range(1, self.h_sys_AA.GetNbinsX()+1):
                if not np.isclose(self.h_sys_AA.GetBinCenter(bin), self.h_sys_jewel.GetBinCenter(bin)):
                    sys.exit('ERROR: JEWEL systematic and data do not align')
                content = self.h_sys_AA.GetBinContent(bin)
                jewel_error_percent = self.h_sys_jewel.GetBinContent(bin)
                jewel_error = jewel_error_percent/100. * content
                old_systematic = self.h_sys_AA.GetBinError(bin)
                new_systematic = np.sqrt(jewel_error*jewel_error + old_systematic*old_systematic)
                self.h_sys_AA.SetBinError(bin, new_systematic)
                #print(self.h_sys_AA.GetBinCenter(bin))
                #print(f'content: {content}')
                #print(f'jewel_error_percent: {jewel_error_percent}')
                #print(f'jewel_error: {jewel_error}')
                #print(f'old_systematic: {old_systematic}')
                #print(f'new_systematic: {new_systematic}')
                #print()
                self.plot_systematic_uncertainties()

        # Option: Translate theta_g to R_g
        if self.observable == 'theta_g' and self.plot_axis_rg:
            self.h_main_pp = self.translate_rg_theta_g(self.h_main_pp, 'rg')
            self.h_sys_pp = self.translate_rg_theta_g(self.h_sys_pp, 'rg')
            self.h_main_AA = self.translate_rg_theta_g(self.h_main_AA, 'rg')
            self.h_sys_AA = self.translate_rg_theta_g(self.h_sys_AA, 'rg')
        
        self.h_ratio = self.h_main_AA.Clone()
        self.h_ratio.Divide(self.h_main_pp)
        self.h_ratio_sys = self.h_sys_AA.Clone()
        self.h_ratio_sys.Divide(self.h_sys_pp)
    
        if self.debug_level > 0:
            integral_inclusive = self.h_main_pp.Integral(1, self.h_main_pp.GetNbinsX()+1, 'width')
            integral_reported = self.h_main_pp.Integral(2, self.h_main_pp.GetNbinsX()+1, 'width')
            tagging_fraction = integral_reported / integral_inclusive
            print('tagging fraction pp re-computed: {}'.format(tagging_fraction))
            
            integral_inclusive = self.h_main_AA.Integral(1, self.h_main_pp.GetNbinsX()+1, 'width')
            integral_reported = self.h_main_AA.Integral(2, self.h_main_pp.GetNbinsX()+1, 'width')
            tagging_fraction = integral_reported / integral_inclusive
            print('tagging fraction AA re-computed: {}'.format(tagging_fraction))
    
        pt = (self.min_pt + self.max_pt)/2
        bin = self.h_tagging_fraction_pp.GetXaxis().FindBin(pt)
        self.f_tagging_pp = self.h_tagging_fraction_pp.GetBinContent(bin)
        pt = (self.min_pt + self.max_pt)/2
        bin = self.h_tagging_fraction_AA.GetXaxis().FindBin(pt)
        self.f_tagging_AA = self.h_tagging_fraction_AA.GetBinContent(bin)

        self.ymax = self.h_main_pp.GetMaximum()
        
        # Load theory predictions
        if self.plot_theory:
            self.init_theory_james()
            
        # Add HEPData tables
        self.add_hepdata_table(is_pp=True)
        self.add_hepdata_table(is_pp=False)
        self.add_hepdata_table(is_ratio=True)
        
    #---------------------------------------------------------------
    # Initialize
    #---------------------------------------------------------------
    def init_result_laura(self, result):
    
        # Init from config file
        self.observable = result['observable']
        self.jetR = result['jetR']
        self.zcut = result['zcut']
        self.min_pt = result['min_pt']
        self.max_pt = result['max_pt']
        if 'xmax' in result:
            self.rg_xmax = result['xmax']
        self.file_name = result['file']
        self.set_xy_titles()
        
        h_PbPb_name = self.config_laura['h_PbPb']
        h_pp_name = self.config_laura['h_pp']
        g_PbPb_sys_name = self.config_laura['g_PbPb_sys']
        g_pp_sys_name = self.config_laura['g_pp_sys']
        h_ratio_name = self.config_laura['h_ratio']
        g_ratio_sys_name = self.config_laura['g_ratio_sys']

        # Get hists from ROOT file
        self.file = ROOT.TFile(self.file_name, 'READ')
        
        self.h_main_pp = self.file.Get(h_pp_name)
        self.h_sys_pp = self.file.Get(g_pp_sys_name)
        self.h_main_AA = self.file.Get(h_PbPb_name)
        self.h_sys_AA = self.file.Get(g_PbPb_sys_name)
        self.h_ratio = self.file.Get(h_ratio_name)
        self.h_ratio_sys = self.file.Get(g_ratio_sys_name)
    
        # Option: Translate R_g to theta_g
        if self.observable == 'rg' and not self.plot_axis_rg:
            self.h_main_pp = self.translate_rg_theta_g(self.h_main_pp, 'theta_g', tgraph=False)
            self.h_sys_pp = self.translate_rg_theta_g(self.h_sys_pp, 'theta_g', tgraph=True)
            self.h_main_AA = self.translate_rg_theta_g(self.h_main_AA, 'theta_g', tgraph=False)
            self.h_sys_AA = self.translate_rg_theta_g(self.h_sys_AA, 'theta_g', tgraph=True)
            self.h_ratio = self.translate_rg_theta_g(self.h_ratio, 'theta_g', tgraph=False, ratio=True)
            self.h_ratio_sys = self.translate_rg_theta_g(self.h_ratio_sys, 'theta_g', tgraph=True, ratio=True)

        if self.debug_level > -1:
            integral_inclusive = self.h_main_pp.Integral(1, self.h_main_pp.GetNbinsX()+1, 'width')
            integral_reported = self.h_main_pp.Integral(2, self.h_main_pp.GetNbinsX()+1, 'width')
            self.f_tagging_pp = integral_reported / integral_inclusive
            print('tagging fraction pp re-computed: {}'.format(self.f_tagging_pp))
            
            integral_inclusive = self.h_main_AA.Integral(1, self.h_main_pp.GetNbinsX()+1, 'width')
            integral_reported = self.h_main_AA.Integral(2, self.h_main_pp.GetNbinsX()+1, 'width')
            self.f_tagging_AA = integral_reported / integral_inclusive
            print('tagging fraction AA re-computed: {}'.format(self.f_tagging_AA))
    
        self.ymax = self.h_main_pp.GetMaximum()
        if 'ymax' in result:
            self.ymax = result['ymax']
        
        # Load theory predictions
        if self.plot_theory:
        
            self.label_list = []
            self.sublabel_list = []
            self.prediction_g_list = []
            
            filename = self.config_laura['file_theory']
            self.file_theory = ROOT.TFile(filename, 'READ')
            
            # JETSCAPE
            bins = np.array(result['jetscape']['xbins'])
            x = (bins[1:] + bins[:-1]) / 2
            ratio = np.array(result['jetscape']['ratio'])
            ratio_err = np.array(result['jetscape']['ratio_err'])
            n = len(x)
            xerr = np.zeros(n)
            g = ROOT.TGraphErrors(n, x, ratio, xerr, ratio_err)
            self.prediction_g_list.append(g)
            self.label_list.append('JETSCAPE')
            self.sublabel_list.append('')
            
            theory_list = result['theory']
            for theory in theory_list:
                key, label = list(theory.items())[0]
                self.label_list.append(label)
                self.sublabel_list.append('')
                
                g = self.file_theory.Get(key)
                
                # Option: Translate R_g to theta_g
                if self.observable == 'rg' and not self.plot_axis_rg:
                    g = self.translate_rg_theta_g(g, 'theta_g', tgraph=True, ratio=True)
            
                self.prediction_g_list.append(g)
            
            self.set_draw_order()
            
        # Add HEPData tables
        self.add_hepdata_table(is_pp=True)
        self.add_hepdata_table(is_pp=False)
        self.add_hepdata_table(is_ratio=True)

    #----------------------------------------------------------------------
    def add_hepdata_table(self, is_pp=False, is_ratio=False):
  
        result_index = self.set_hepdata_table_index()
        if is_ratio:
            h = self.h_ratio
            h_sys = self.h_ratio_sys
            index = 3*result_index
        elif is_pp:
            h = self.h_main_pp
            h_sys = self.h_sys_pp
            index = 3*result_index-2
        else:
            h = self.h_main_AA
            h_sys = self.h_sys_AA
            index = 3*result_index-1
        table = hepdata_lib.Table(f'Table {index}')

        if is_ratio:
            table.keywords["reactions"] = ['P P --> jet+X, Pb Pb --> jet+X']
        elif is_pp:
            table.keywords["reactions"] = ['P P --> jet+X']
        else:
            table.keywords["reactions"] = ['Pb Pb --> jet+X']
        table.keywords["cmenergies"] = ['5020']
        
        # Set observable-specific info in table
        x_label, y_label = self.set_hepdata_table_descriptors(table, result_index, is_pp, is_ratio)
     
        # Define variables
        h = hepdata_lib.root_utils.get_hist_1d_points(h)
     
        n_significant_digits = 3 # Note: use 2 significant digits for uncertainties
        
        x = hepdata_lib.Variable(x_label, is_independent=True, is_binned=True, units='')
        x.digits = n_significant_digits
        x.values = h['x_edges']
     
        y = hepdata_lib.Variable(y_label, is_independent=False, is_binned=False, units='')
        y.digits = n_significant_digits
        y.values = h['y']
        if is_ratio:
            y.add_qualifier('RE', 'P P --> jet+X, Pb Pb --> jet+X')
            y.add_qualifier('CENTRALITY', str(self.centrality))
        elif is_pp:
            y.add_qualifier('RE', 'P P --> jet+X')
        else:
            y.add_qualifier('RE', 'Pb Pb --> jet+X')
            y.add_qualifier('CENTRALITY', str(self.centrality))
        y.add_qualifier('SQRT(S)', 5.02, 'TeV')
        y.add_qualifier('ETARAP', '|0.9-R|')
        y.add_qualifier('jet radius', self.jetR)
        y.add_qualifier('jet method', 'Anti-$k_{T}$')
 
        # Define uncertainties
        stat = hepdata_lib.Uncertainty('stat', is_symmetric=True)
        stat.values = [float('{:.2g}'.format(dy)) for dy in h['dy']]
     
        # Add tables to submission
        table.add_variable(x)
        table.add_variable(y)
        y.add_uncertainty(stat)
     
        # Add unfolding systematic
        if self.centrality == [0,10]:
            h_sys = hepdata_lib.root_utils.get_hist_1d_points(h_sys)
            sys = hepdata_lib.Uncertainty('sys', is_symmetric=True)
            sys.values = [float('{:.2g}'.format(dy)) for dy in h_sys['dy']]
        elif self.centrality == [30,50]:
            h_sys = hepdata_lib.root_utils.get_graph_points(h_sys)
            sys = hepdata_lib.Uncertainty('sys', is_symmetric=False)
            sys.values = [(float('{:.2g}'.format(dy[0])), float('{:.2g}'.format(dy[1]))) for dy in h_sys['dy']]
        y.add_uncertainty(sys)
     
        # Add table to the submission
        self.hepdata_submission.add_table(table)

    #----------------------------------------------------------------------
    def set_hepdata_table_index(self):
    
        if self.observable == 'zg':
            if np.isclose(self.jetR, 0.2):
                return 1
            elif np.isclose(self.jetR, 0.4):
                if self.centrality == [0,10]:
                    return 5
                elif self.centrality == [30,50]:
                    return 2
        elif self.observable in ['theta_g', 'rg']:
            if np.isclose(self.jetR, 0.2):
                return 3
            elif np.isclose(self.jetR, 0.4):
                if np.isclose(self.zcut, 0.2):
                    return 4
                elif np.isclose(self.zcut, 0.4):
                    return 6

    #----------------------------------------------------------------------
    def set_hepdata_table_descriptors(self, table, result_index, is_pp=False, is_ratio=False):
        
        if self.observable == 'zg':
            table.description = r'Groomed jet momentum splitting fraction $z_{{\mathrm{g}}}$'
            x_label = r'$z_{{\mathrm{g}}}$'
            if is_ratio:
                y_label = 'PbPb/pp'
            else:
                y_label = r'$\frac{1}{\sigma_{inc}} \frac{d\sigma}{dz_{{\mathrm{g}}}}$'
            if result_index == 1:
                table.location = 'Figure 2 (left)'
            elif result_index == 2:
                table.location = 'Figure 2 (right)'
            elif result_index == 5:
                table.location = 'Figure 2 (right) in Ref. [Public Note]'
            
        elif self.observable in ['theta_g', 'rg']:
            table.description = r'Groomed jet radius (scaled) $\theta_{{\mathrm{g}}}$'
            x_label = r'$\theta_{{\mathrm{g}}}$'
            if is_ratio:
                y_label = 'PbPb/pp'
            else:
                y_label = r'$\frac{1}{\sigma_{inc}} \frac{d\sigma}{d\theta_{{\mathrm{g}}}}$'
            if result_index == 3:
                table.location = 'Figure 3 (left)'
            elif result_index == 4:
                table.location = 'Figure 3 (right)'
            elif result_index == 6:
                table.location = 'Figure 5 (right) in Ref. [Public Note]'
                
        if is_ratio:
            table.description += r' $-$ ratio of Pb-Pb to pp collisions.'
        elif is_pp:
            table.description += ' in pp collisions.'
        else:
            table.description += ' in Pb-Pb collisions.'
                
        table.description += '\n'
        table.description += r'${}<p_{{\mathrm{{T}}}}^{{\mathrm{{ch jet}}}}<{} \;\mathrm{{GeV}}/c$, Soft Drop $z_{{\mathrm{{cut}}}}={}, \beta=0$.'.format(self.min_pt, self.max_pt, self.zcut)
        table.description += '\n\nNote: The first bin corresponds to the Soft Drop untagged fraction.'
        
        if is_pp:
            table.description += '\n\n'
            table.description += r'For the "trkeff" and "generator" systematic uncertainty sources, the signed systematic uncertainty breakdowns ($\pm$ vs. $\mp$), denote correlation across bins (both within this table, and across tables for a given centrality). For the remaining sources ("unfolding") no correlation information is specified ($\pm$ is always used).'
        
        return x_label, y_label

    #---------------------------------------------------------------
    # This function is called once for each subconfiguration
    #
    # Required histograms:
    #   - self.h_main_pp
    #   - self.h_sys_pp (tgraph for laura)
    #   - self.h_main_AA
    #   - self.h_sys_AA (tgraph for laura)
    #   - self.h_ratio
    #   - self.h_ratio_sys (tgraph for laura)
    #
    # And the theory tgraphs should be placed in:
    #   - self.label_list
    #   - self.prediction_g_list
    #
    # As well as the tagging rates:
    #   - self.f_tagging_pp (number)
    #   - self.f_tagging_AA (number)
    #----------------------------------------------------------------------
    def plot_result(self):
        print('Plot final results for {}: R = {}, {}'.format(self.observable, self.jetR, self.obs_label))
        
        observable_label = self.observable
        if self.observable in ['rg', 'theta_g']:
            if self.plot_axis_rg:
                observable_label = 'rg'
            else:
                observable_label = 'theta_g'
                
        if self.plot_theory:
            if self.plot_data:
                name = 'h_{}_{}-{}_R{}_zcut{}_Theory.pdf'.format(observable_label, self.centrality[0], self.centrality[1], self.utils.remove_periods(self.jetR), self.utils.remove_periods(self.zcut))
            else:
                name = 'h_{}_{}-{}_R{}_zcut{}_TheoryOnly.pdf'.format(observable_label, self.centrality[0], self.centrality[1], self.utils.remove_periods(self.jetR), self.utils.remove_periods(self.zcut))
        else:
            name = 'h_{}_{}-{}_R{}_zcut{}.pdf'.format(observable_label, self.centrality[0], self.centrality[1], self.utils.remove_periods(self.jetR), self.utils.remove_periods(self.zcut))
        self.output_filename = os.path.join(self.output_dir, name)

        self.utils.set_plotting_options()
        ROOT.gROOT.ForceStyle()

        c = ROOT.TCanvas('c', 'c', 600, 650)
        c.Draw()
        
        c.cd()
        pad2_dy = 0.45
        pad1 = ROOT.TPad('myPad', 'The pad',0,pad2_dy,1,1)
        pad1.SetLeftMargin(0.2)
        pad1.SetTopMargin(0.08)
        pad1.SetRightMargin(0.04)
        pad1.SetBottomMargin(0.)
        if self.observable in ['rg', 'theta_g']:
            pad1.SetTicks(0,1)
        else:
            pad1.SetTicks(1,1)
        pad1.Draw()
        pad1.cd()

        myLegend = ROOT.TLegend(0.24,0.65,0.44,0.85)
        self.utils.setup_legend(myLegend,0.06, sep=0.)
        
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
          
        xmin = self.h_main_pp.GetBinLowEdge(2)
        xmax = self.h_main_pp.GetBinLowEdge(self.h_main_pp.GetNbinsX()+1)
        if self.centrality[0] == 30 and self.observable == 'rg':
            if self.plot_axis_rg:
                xmax = self.rg_xmax
            else:
                xmax = self.rg_xmax/self.jetR
            
        myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, xmin, xmax)
        myBlankHisto.SetNdivisions(505)
        myBlankHisto.SetXTitle(self.xtitle)
        myBlankHisto.SetYTitle(self.ytitle)
        if self.observable in ['rg', 'theta_g']:
            ymax = 1.9*self.ymax
        else:
            ymax = 3.*self.ymax
        myBlankHisto.SetMaximum(ymax)
        myBlankHisto.SetMinimum(2e-4) # Don't draw 0 on top panel
        myBlankHisto.GetYaxis().SetTitleSize(0.065)
        myBlankHisto.GetYaxis().SetTitleSize(0.08)
        myBlankHisto.GetYaxis().SetTitleOffset(1.15)
        myBlankHisto.GetYaxis().SetLabelSize(0.075)
        myBlankHisto.Draw('E')
        
        # Draw top x-axis
        if self.observable in ['rg', 'theta_g']:
            if self.plot_axis_rg: # Plot theta_g axis on top
                top_axis_tf1 = ROOT.TF1('top_axis_tf1', 'x', 0, xmax/self.jetR-0.001)
                top_axis = ROOT.TGaxis(xmin, ymax, xmax, ymax, 'top_axis_tf1', 505, '- S')
                top_axis.SetTitle('#it{#theta}_{g}')
            else: # Plot rg axis on top
                top_axis_tf1 = ROOT.TF1('top_axis_tf1', 'x', 0, self.jetR*xmax-0.001)
                top_axis = ROOT.TGaxis(xmin, ymax, xmax, ymax, 'top_axis_tf1', 505, '- S')
                top_axis.SetTitle('#it{R}_{g}')
            top_axis.SetTitleSize(24)
            top_axis.SetTitleFont(43)
            top_axis.SetTitleOffset(0.69)
            top_axis.SetLabelFont(43)
            top_axis.SetLabelSize(25)
            top_axis.SetTickSize(0.025)
            top_axis.SetLabelOffset(0.016)
            top_axis.Draw()

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
        myBlankHisto2.GetYaxis().SetTitleSize(25)
        myBlankHisto2.GetYaxis().SetTitleFont(43)
        myBlankHisto2.GetYaxis().SetTitleOffset(2.)
        myBlankHisto2.GetYaxis().SetLabelFont(43)
        myBlankHisto2.GetYaxis().SetLabelSize(25)
        myBlankHisto2.GetYaxis().SetNdivisions(505)
        
        if self.observable in ['rg', 'theta_g']:
        
            if self.centrality[0] == 0:
          
                myBlankHisto2.GetYaxis().SetRangeUser(0.3, 2.49)
                
                ratio_legend = ROOT.TLegend(0.305,0.68,0.42,0.985)
                self.utils.setup_legend(ratio_legend, 0.06, sep=-0.2)
                ratio_legend2 = ROOT.TLegend(0.55,0.68,0.69,0.985)
                self.utils.setup_legend(ratio_legend2, 0.06, sep=-0.2)
                
            elif self.centrality[0] == 30:
          
                if np.isclose(self.zcut, 0.2):
                    myBlankHisto2.GetYaxis().SetRangeUser(0.4, 1.79)
                else:
                    myBlankHisto2.GetYaxis().SetRangeUser(0.4, 1.89)
                
                ratio_legend = ROOT.TLegend(0.54,0.65,0.7,0.92)
                self.utils.setup_legend(ratio_legend, 0.06, sep=-0.2)
                ratio_legend2 = ROOT.TLegend(0.54,0.65,0.7,0.92)
                self.utils.setup_legend(ratio_legend2, 0.06, sep=-0.2)
            
        else:
        
            if self.centrality[0] == 0:

                if self.jetR == 0.2:
                    myBlankHisto2.GetYaxis().SetRangeUser(0.61, 1.59)
                else:
                    myBlankHisto2.GetYaxis().SetRangeUser(0.41, 1.79)
                    
            elif self.centrality[0] == 30:
                    myBlankHisto2.GetYaxis().SetRangeUser(0.61, 1.59)

            if self.centrality[0] == 0 and self.jetR == 0.4:
                ratio_legend = ROOT.TLegend(0.23,0.68,0.38,0.96)
                self.utils.setup_legend(ratio_legend, 0.07, sep=-0.4)
                ratio_legend2 = ROOT.TLegend(0.52,0.82,0.67,0.96)
                self.utils.setup_legend(ratio_legend2, 0.07, sep=-0.2)
            elif self.centrality[0] == 30:
                ratio_legend = ROOT.TLegend(0.3,0.74,0.45,0.91)
                self.utils.setup_legend(ratio_legend, 0.07, sep=-0.4)
                ratio_legend2 = ROOT.TLegend(0.5,0.72,0.65,0.95)
                self.utils.setup_legend(ratio_legend2, 0.07, sep=-0.4)
            else:
                ratio_legend = ROOT.TLegend(0.3,0.68,0.45,0.96)
                self.utils.setup_legend(ratio_legend, 0.07, sep=-0.4)
                ratio_legend2 = ROOT.TLegend(0.6,0.73,0.75,0.96)
                self.utils.setup_legend(ratio_legend2, 0.07, sep=-0.2)

        myBlankHisto2.Draw('')
          
        self.h_ratio.SetMarkerSize(1.5)
        self.h_ratio.SetMarkerStyle(self.markers[1])
        self.h_ratio.SetMarkerColor(self.ratio_color)
        self.h_ratio.SetLineColor(self.ratio_color)

        self.h_ratio_sys.SetFillColor(self.ratio_color)
        self.h_ratio_sys.SetFillColorAlpha(self.ratio_color, 0.3)

        pad1.cd()
        if self.plot_data:
            self.h_sys_pp.Draw('E2 same')
            self.h_sys_AA.Draw('E2 same')
            self.h_main_pp.DrawCopy('PE X0 same')
            self.h_main_AA.DrawCopy('PE X0 same')
          
        pad2.cd()

        if self.plot_theory:

            for i, g in enumerate(self.prediction_g_list):
                  
                label = self.label_list[i]
                sublabel = self.sublabel_list[i]

                color = self.theory_colors[label+sublabel]
                g.SetLineColor(color)
                g.SetFillColor(color)
                g.SetFillColorAlpha(color, 0.85)
                if type(g) in [ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                    g.SetLineColor(0)
                    if self.observable in ['rg', 'theta_g']:
                        if self.jetR == 0.4:
                            if 'Pablos' in label or 'JETSCAPE' in label or 'Caucal' in label:
                                ratio_legend.AddEntry(g, '{}{}'.format(label, sublabel), 'F')
                            else:
                                ratio_legend2.AddEntry(g, '{}{}'.format(label, sublabel), 'F')
                        elif self.jetR == 0.2:
                            if 'JEWEL' in label:
                                ratio_legend2.AddEntry(g, '{}{}'.format(label, sublabel), 'F')
                                if 'recoils off' in sublabel:
                                    g.SetFillStyle(3144)
                                if 'recoils on' in sublabel:
                                    g.SetFillStyle(3144)
                            else:
                                ratio_legend.AddEntry(g, '{}{}'.format(label, sublabel), 'F')
                    elif self.observable == 'zg':
                        if 'Pablos' in label or 'JETSCAPE' in label or 'Caucal' in label:
                            ratio_legend.AddEntry(g, '{}{}'.format(label, sublabel), 'F')
                        else:
                            if 'JEWEL' in label:
                                if 'recoils off' in sublabel:
                                    g.SetFillStyle(3144)
                                if 'recoils on' in sublabel:
                                    g.SetFillStyle(3144)
                            ratio_legend2.AddEntry(g, '{}{}'.format(label, sublabel), 'F')
                elif type(g) == ROOT.TGraph:
                    g.SetLineStyle(self.line_style[label+sublabel])
                    g.SetLineWidth(4)
                    if self.observable in ['rg', 'theta_g']:
                        if self.jetR == 0.4:
                            if 'Pablos' in label or 'JETSCAPE' in label:
                                ratio_legend.AddEntry(g, '{}{}'.format(label, sublabel), 'L')
                            else:
                                ratio_legend2.AddEntry(g, '{}{}'.format(label, sublabel), 'L')
                        elif self.jetR == 0.2:
                            if 'JEWEL' in label:
                                ratio_legend2.AddEntry(g, '{}{}'.format(label, sublabel), 'F')
                                if 'recoils off' in sublabel:
                                    g.SetFillStyle(3144)
                                if 'recoils on' in sublabel:
                                    g.SetFillStyle(3144)
                            if 'Yuan' in label:
                                ratio_legend2.AddEntry(g, '{}{}'.format(label, sublabel), 'L')
                            else:
                                ratio_legend.AddEntry(g, '{}{}'.format(label, sublabel), 'L')
                    elif self.observable == 'zg':
                        if 'Pablos' in label:
                            ratio_legend.AddEntry(g, '{}{}'.format(label, sublabel), 'L')
                        else:
                            ratio_legend2.AddEntry(g, '{}{}'.format(label, sublabel), 'L')

        # Draw curves in specified order
        if self.plot_theory:
            for i in self.draw_order:
                g = self.prediction_g_list[i]
                if type(g) in [ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                    g.Draw("3 same")
                elif type(g) == ROOT.TGraph:
                    g.Draw("same")
          
        ratio_legend.Draw()
        ratio_legend2.Draw()

        line = ROOT.TLine(xmin,1,xmax,1)
        line.SetLineColor(920+2)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
         
        if self.plot_data:
            self.h_ratio_sys.Draw('E2 same')
            self.h_ratio.Draw('PE X0 same')
       
        pad1.cd()
        myLegend.AddEntry(self.h_main_pp, 'pp', 'PE')
        myLegend.AddEntry(self.h_main_AA, 'Pb#font[122]{{-}}Pb {}#font[122]{{-}}{}%'.format(self.centrality[0], self.centrality[1]), 'PE')
        myLegend.AddEntry(self.h_sys_pp, 'Sys. uncertainty', 'f')
        
        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        
        text_latex.SetTextSize(0.055)
        text = '#it{{f}}_{{tagged}}^{{pp}} = {:.2f}, #it{{f}}_{{tagged}}^{{AA}} = {:.2f}'.format(self.f_tagging_pp, self.f_tagging_AA)
        text_latex.DrawLatex(0.55, 0.36, text)

        text_latex.SetTextSize(0.065)
        text = f'ALICE {self.figure_approval_status} #sqrt{{#it{{s_{{#it{{NN}}}}}}}} = 5.02 TeV'
        text_latex.DrawLatex(0.55, 0.83, text)
        
        text = 'Charged-particle jets'#, anti-#it{k}_{T}'
        text_latex.DrawLatex(0.55, 0.75, text)
        
        text = '#it{R} = ' + str(self.jetR) + ',   | #it{{#eta}}_{{jet}}| < {}'.format(0.9-self.jetR)
        text_latex.DrawLatex(0.55, 0.67, text)
        
        text = str(self.min_pt) + ' < #it{p}_{T, ch jet} < ' + str(self.max_pt) + ' GeV/#it{c}'
        text_latex.DrawLatex(0.55, 0.58, text)
        
        text = self.formatted_grooming_label.format(self.zcut)
        text_latex.DrawLatex(0.55, 0.48, text)
        
        myLegend.Draw()

        c.SaveAs(self.output_filename)
        c.Close()

    #---------------------------------------------------------------
    # Set draw order (slightly hacky-- just move selected curves to start or end of draw list)
    #----------------------------------------------------------------------
    def set_draw_order(self):
        
        self.draw_order = list(range(0, len(self.prediction_g_list)))
        for i, g in enumerate(self.prediction_g_list):

            label = self.label_list[i]
            sublabel = self.sublabel_list[i]

            if 'Chien' in label:
                index = self.draw_order.index(i)
                self.draw_order.append(self.draw_order.pop(index))

            if 'JEWEL' in label:
                index = self.draw_order.index(i)
                self.draw_order.insert(0, self.draw_order.pop(index))

            if 'Pablos' in label:
                index = self.draw_order.index(i)
                self.draw_order.insert(0, self.draw_order.pop(index))

            if 'JETSCAPE' in label:
                index = self.draw_order.index(i)
                self.draw_order.append(self.draw_order.pop(index))
          
    #---------------------------------------------------------------
    # Translate R_g <--> theta_g
    # Options: 'rg': theta_g --> r_g
    #          'theta_g': r_g --> theta_g
    #----------------------------------------------------------------------
    def translate_rg_theta_g(self, h, option='', tgraph=False, ratio=False):
    
        if option == 'rg':
            scale_factor = self.jetR
        elif option == 'theta_g':
            scale_factor = 1/self.jetR

        if tgraph:
        
            x = ROOT.Double()
            y = ROOT.Double()
            for i in range(h.GetN()):
                h.GetPoint(i,x,y)
                
                if type(h) in [ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:

                    if type(h) in [ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                        x_low = h.GetErrorXlow(i)
                        x_high = h.GetErrorXhigh(i)
                        y_low = h.GetErrorYlow(i)
                        y_high = h.GetErrorYhigh(i)
                        
                        x_low_new = x_low*scale_factor
                        x_high_new = x_high*scale_factor
                    
                if ratio:
                    y_new = y
                    if type(h) in [ROOT.TGraphErrors, ROOT.TGraphAsymmErrors]:
                        y_low_new = y_low
                        y_high_new = y_high
                else:
                    y_new = y/scale_factor
                    y_low_new = y_low/scale_factor
                    y_high_new = y_high/scale_factor
                            
                h.SetPoint(i,x*scale_factor,y_new)
                if type(h) in [ROOT.TGraphAsymmErrors]:
                    h.SetPointError(i, x_low_new, x_high_new, y_low_new, y_high_new)
                elif type(h) in [ROOT.TGraphErrors]:
                    h.SetPointError(i, x_low_new, y_low_new)
                
            return h
            
        else:
            original_bins = h.GetXaxis().GetXbins()
            new_bins = [bin*scale_factor for bin in original_bins]
            
            h_new = ROOT.TH1F('hnew{}'.format(h.GetName()), 'hnew{}'.format(h.GetName()),
                              len(new_bins)-1, array('d',(new_bins)))
            
            for i in range(1, len(original_bins)+1):
                content = h.GetBinContent(i)
                uncertainty =  h.GetBinError(i)
                h_new.SetBinContent(i, content)
                h_new.SetBinError(i, uncertainty)
    
            if not ratio:
                h_new.Scale(1./scale_factor)
        
            return h_new

    #---------------------------------------------------------------
    # Set axis title for given observable
    #----------------------------------------------------------------------
    def set_xy_titles(self):
    
        if self.observable in ['rg', 'theta_g']:
            if self.plot_axis_rg:
                self.xtitle = '#it{R}_{g}'
            else:
                self.xtitle = '#it{#theta}_{g}'
        if self.observable == 'zg':
            self.xtitle = '#it{z}_{g}'
        
        self.ytitle = '#frac{{1}}{{#it{{#sigma}}_{{jet, inc}}}} #frac{{d#it{{#sigma}}}}{{d{}}}'.format(self.xtitle)
        
    #---------------------------------------------------------------
    # Initialize theory predictions
    #---------------------------------------------------------------
    def init_theory_james(self):
    
        # Get theory predictions
        if self.plot_theory:
            if self.observable == 'theta_g':
                config_file = 'theory_predictions_rg.yaml'
            if self.observable == 'zg':
                config_file = 'theory_predictions_zg.yaml'
        
        # Get binning, for re-binning
        bin_edges = np.array(self.h_main_pp.GetXaxis().GetXbins())[1:]
        self.bin_array = array('d', bin_edges)

        with open(config_file, 'r') as stream:
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
                
                    if type == 'lbnl':
                    
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
                            #g = ROOT.TGraphAsymmErrors(n, x, ratio, xerr, xerr, ratio-ratio_lower, ratio_upper-ratio)
                            g = ROOT.TGraph(n, x, ratio)

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

                    elif type == 'jewel':

                        pp_filename = theory_prediction['pp_filename']
                        AA_filename = theory_prediction['AA_filename']
                        f_pp = ROOT.TFile(pp_filename, 'READ')
                        f_AA = ROOT.TFile(AA_filename, 'READ')

                        sublabel = theory_prediction['sublabel']

                        h_name_pp = f'h_{self.observable}_JetPt_R{self.jetR}Scaled'
                        if ', recoils on (GridSub' in sublabel:
                            gridsize = theory_prediction['gridsize']
                            h_name_AA = f'h_{self.observable}_JetPt_R{self.jetR}_gridsub_{gridsize}Scaled'
                        else:
                            h_name_AA = f'h_{self.observable}_JetPt_R{self.jetR}Scaled'
                        h_pp_2D = f_pp.Get(h_name_pp)
                        h_pp_2D.SetDirectory(0)
                        f_pp.Close()
                        h_AA_2D = f_AA.Get(h_name_AA)
                        h_AA_2D.SetDirectory(0)
                        f_AA.Close()

                        # Select pt bin and project onto z_r axis
                        h_pp_2D.GetXaxis().SetRangeUser(self.min_pt, self.max_pt)
                        h_pp = h_pp_2D.ProjectionY()
                        h_pp.SetName(f'{h_pp.GetName()}_pp_{sublabel}')        

                        h_AA_2D.GetXaxis().SetRangeUser(self.min_pt, self.max_pt)
                        h_AA = h_AA_2D.ProjectionY()  
                        h_AA.SetName(f'{h_AA.GetName()}_AA_{sublabel}')        

                        # Rebin to data binning
                        h_pp = h_pp.Rebin(bin_edges.size-1, f'{h_pp.GetName()}_rebinned', bin_edges)
                        h_AA = h_AA.Rebin(bin_edges.size-1, f'{h_AA.GetName()}_rebinned', bin_edges)
                        
                        # Normalization
                        h_pp.Scale(1., 'width')
                        h_AA.Scale(1., 'width')
                        h_pp.Scale(1./h_pp.Integral(1, h_pp.GetNbinsX()))
                        h_AA.Scale(1./h_AA.Integral(1, h_pp.GetNbinsX()))

                        # Form ratio
                        hRatio = h_AA.Clone()
                        hRatio.SetName(f'hRatio_{h_AA.GetName()}')
                        hRatio.Divide(h_pp)

                        # Construct TGraph
                        g = ROOT.TGraphErrors(hRatio)

                    # Option: Translate R_g to theta_g
                    if self.observable == 'theta_g' and self.plot_axis_rg:
                        g = self.translate_rg_theta_g(g, 'rg', tgraph=True, ratio=True)
    
                    if prediction in plot_list:
                        self.prediction_g_list.append(g)
                        self.label_list.append(theory['label'])
                        self.sublabel_list.append(theory_prediction['sublabel'])

        self.set_draw_order()

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

    #----------------------------------------------------------------------
    def plot_systematic_uncertainties(self):

        name = 'cSys'
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

        leg = ROOT.TLegend(0.67,0.75,0.8,0.92)
        self.utils.setup_legend(leg,0.04)

        h_systematic_total = self.h_sys_AA.Clone(f'{self.h_sys_AA.GetName()}_systotal')
        for bin in range(1, self.h_sys_jewel.GetNbinsX()+1):
            content = self.h_sys_AA.GetBinContent(bin)
            uncertainty = self.h_sys_AA.GetBinError(bin)
            h_systematic_total.SetBinContent(bin, np.abs(uncertainty/content*100.))
            h_systematic_total.SetBinError(bin, 0)
        h_systematic_total.SetTitle('')
        h_systematic_total.SetXTitle( getattr(self, 'xtitle') )
        h_systematic_total.GetXaxis().SetRange(2, h_systematic_total.GetNbinsX())
        h_systematic_total.GetYaxis().SetTitleOffset(1.5)
        h_systematic_total.SetYTitle('Systematic uncertainty (%)')
        h_systematic_total.SetMaximum(2.7*h_systematic_total.GetMaximum(50))
        h_systematic_total.SetMinimum(1e-3)
        h_systematic_total.SetLineColor(ROOT.kTeal-8)
        h_systematic_total.SetLineStyle(1)
        h_systematic_total.SetLineWidth(2)
        h_systematic_total.Draw("hist X0")
        leg.AddEntry(h_systematic_total, 'Total', 'l')

        self.h_sys_jewel.SetMarkerSize(0)
        self.h_sys_jewel.SetLineColor(ROOT.kViolet-8)
        self.h_sys_jewel.SetLineStyle(1)
        self.h_sys_jewel.SetLineWidth(2)
        self.h_sys_jewel.DrawCopy('X0 same')
        self.h_sys_jewel.Scale(-1)
        self.h_sys_jewel.DrawCopy('X0 same')
        leg.AddEntry(self.h_sys_jewel, 'jet quenching model', 'l')

        leg.Draw()

        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        text = str(self.min_pt) + ' < #it{p}_{T, ch jet} < ' + str(self.max_pt)
        text_latex.DrawLatex(0.3, 0.85, text)

        text_latex = ROOT.TLatex()
        text_latex.SetNDC()
        text = '#it{R} = ' + str(self.jetR)
        text_latex.DrawLatex(0.3, 0.78, text)

        text_latex.DrawLatex(0.3, 0.64, self.formatted_grooming_label.format(self.zcut))

        outputFilename = os.path.join(self.output_dir, f'hSystematics_{self.observable}_R{self.utils.remove_periods(self.jetR)}_{self.obs_label}_{int(self.min_pt)}-{int(self.max_pt)}.pdf')
        c.SaveAs(outputFilename)
        c.Close()

#----------------------------------------------------------------------
if __name__ == '__main__':

  analysis = PlotRAA()
  analysis.plot_raa()
