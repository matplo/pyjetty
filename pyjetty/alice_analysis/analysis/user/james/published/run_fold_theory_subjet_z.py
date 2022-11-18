#! /usr/bin/env python

"""
Code to do theory folding in order to compare to the measured distributions
The class 'TheoryFolding' below inherits from the 'TheoryFolding' class in:
pyjetty/alice_analysis/analysis/user/substructure/run_fold_theory.py

Based on code from Rey Cruz-Torres.
"""

import sys
import os
import argparse
from array import *
import numpy as np
import ROOT
ROOT.gSystem.Load("$HEPPY_DIR/external/roounfold/roounfold-current/lib/libRooUnfold.so")

from pyjetty.alice_analysis.analysis.user.substructure import run_fold_theory

# Load pyjetty ROOT utils
ROOT.gSystem.Load('libpyjetty_rutil')

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
################################################################
################################################################
class TheoryFolding(run_fold_theory.TheoryFolding):

    #----------------------------------------------------------------------
    # Load theory curves
    #----------------------------------------------------------------------
    def load_theory_curves(self):

        self.theory_scale_vars = {}

        # open theory file, load data
        th_file = os.path.join(self.theory_dir, 'subjet_ALICE.npy') 
        print(f'reading from files in: {th_file}')
        SCET_data=np.load(th_file,allow_pickle=True).item(0)
        # print(SCET_data.keys())

        # Loop over each jet R specified in the config file
        for jetR in self.jetR_list:
            scale_var = []

            # Loop through subconfigurations to fold (e.g. in the jet-axis analysis there Standard_WTA, Standard_SD_1, ...)
            for i,conf in enumerate(self.obs_subconfig_list):
                obs_setting = self.obs_settings[i]            # labels such as 'Standard_WTA'
                grooming_setting = self.grooming_settings[i]  # grooming parameters
                label = obs_setting

                pt_bins = array('d', self.theory_pt_bins)

                # The commented-out section below gives the input-level observable distribution the same
                # binning as the output-level binning, but here I will default to the RM binning, so that
                # the input level has as fine binning as possible
                #if self.obs_bin_option == 'nominal_ana_binning':
                #   obs_bins = array('d',self.subobs_specific_binning[jetR][conf])
                #elif self.obs_bin_option == 'theory_obs_binning':
                #   obs_bins = array('d', self.theory_obs_bins)   # bins which we want to have in the result
                #else:
                #   obs_bins = array('d',getattr(self,'binning_R%s_'%((str)(jetR).replace('.',''))+obs_setting))

                obs_bins = array('d',getattr(self,'binning_R%s_'%((str)(jetR).replace('.',''))+str(obs_setting)))
                obs_width = np.subtract(obs_bins[1:],obs_bins[:-1])
                #print(f'obs_bins: {obs_bins}')

                # -----------------------------------------------------        
                # Create histograms where theory curves will be stored
                th_hists_no_scaling = []          # Basically a copy of the theory calculations, but binned
                th_hists = []                     # Histograms that will actually be used in the folding
                hist_names = []

                # -----------------------------------------------------
                # Get appropriate subobservable from SCET calculation
                SCET_data_obs_setting = SCET_data['r=%f'%obs_setting]
                #print(SCET_data_obs_setting.keys())

                # loop over pT bins
                for p, pt in enumerate(pt_bins[:-1]):

                    # Get scale factor for this pT bin.
                    # This reverses the self-normalization of 1/sigma for correct pT scaling when doing projections onto the y-axis.
                    scale_f = self.pt_scale_factor_jetR(pt,pt_bins[p+1],jetR)

                    # Get appropriate pt bin from SCET calculation
                    SCET_data_obs_setting_pt_bin = SCET_data_obs_setting['pT=%f'%pt]
                    #print(SCET_data_obs_setting_pt_bin.keys())
                    
                    if 'leading' in self.observable:
                        SCET_data_variations = SCET_data_obs_setting_pt_bin['leading']
                    elif 'inclusive' in self.observable:
                        SCET_data_variations = SCET_data_obs_setting_pt_bin['inclusive']
                    # print('scales, zbins=',SCET_data_variations.shape)

                    # Get z_r values
                    # Remove every other value in the array, since they are given as [0., 0.001, 0.001, 0.002, 0.002, ..., 0.999, 0.999, 1.0]
                    # Then take center value as "bin center"
                    zr_edges = SCET_data_obs_setting_pt_bin['z'][::20]
                    zr_edges = np.append(zr_edges, SCET_data_obs_setting_pt_bin['z'][-1])
                    #print(zr_edges)
                    zr_values = (zr_edges[1:] + zr_edges[:-1]) / 2
                    #print('zr_values:',zr_values.shape)
                    #print([zr_values[i] for i in range(zr_values.size)])

                    # loop over scale variations and fill histograms
                    n_scale_variations = SCET_data_variations.shape[0]
                    for sv in range(0,n_scale_variations):

                        # Get cross-section
                        # Remove every other value in the array, since they are given as [0., 0.001, 0.001, 0.002, 0.002, ..., 0.999, 0.999, 1.0]
                        y_val_n = SCET_data_variations[sv][1::2]
                        # Then average every 10 elements together, in order to match the binning of our RM [0, 0.01, 0.02, ..., 0.99, 1.0]
                        y_val_n =np.mean(y_val_n.reshape(-1, 10), axis=1)
                        #print(y_val_n)
                        #print([y_val_n[i] for i in range(y_val_n.size)])

                        # Interpolate the given values and return the value at the requested bin center
                        y_val_bin_ctr = self.interpolate_values_linear(zr_values,y_val_n,obs_bins)

                        if p==0:
                            hist_name = 'h2_input_%s_R%s_obs_pT_%s' % ( self.observable , (str)(jetR).replace('.','') , obs_setting )
                            hist_name += '_sv%i' % (sv)

                            hist_name_no_scaling = hist_name + '_no_scaling'

                            th_hist                     = ROOT.TH2D(hist_name                    ,';p_{T}^{jet};%s'%(self.observable), len(pt_bins)-1, pt_bins, len(obs_bins)-1, obs_bins)
                            th_hist_no_scaling          = ROOT.TH2D(hist_name_no_scaling         ,';p_{T}^{jet};%s'%(self.observable), len(pt_bins)-1, pt_bins, len(obs_bins)-1, obs_bins)

                            th_hists.append(th_hist)
                            hist_names.append(hist_name)
                            th_hists_no_scaling.append(th_hist_no_scaling)

                        # Save content into histogram before any scaling has been applied (to compare to the theory curves and make sure everything went fine)
                        for ob in range(0,len(obs_bins)-1):
                            th_hists_no_scaling[sv].SetBinContent(p+1,ob+1,y_val_bin_ctr[ob])

                        # Multiply by bin width and scale with pT-dependent factor
                        y_val_bin_ctr = np.multiply(y_val_bin_ctr,obs_width)
                        integral_y_val_bin_ctr = sum(y_val_bin_ctr)
                        y_val_bin_ctr = [ val * scale_f / integral_y_val_bin_ctr for val in y_val_bin_ctr ]

                        # Save scaled content into the histograms
                        for ob in range(0,len(obs_bins)-1):
                            th_hists[sv].SetBinContent(p+1,ob+1,y_val_bin_ctr[ob])

                # ------------------------------------------------------------------------------------------------------------  
                new_obs_lab = obs_setting
                for n_pt in range(0,len(self.final_pt_bins)-1):
                    histo_list = []
                    for sv in range(0,n_scale_variations):
                        projection_name = 'h1_input_%s_R%s_%s_sv%i_pT_%i_%i' % ( self.observable,(str)(jetR).replace('.',''),obs_setting,sv,(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]))

                        # Determine the bin number that corresponds to the pT edges given           
                        min_bin, max_bin = self.bin_position( self.theory_pt_bins , self.final_pt_bins[n_pt],self.final_pt_bins[n_pt+1] )

                        h1_input_hist = th_hists[sv].ProjectionY(projection_name,min_bin,max_bin)
                        h1_input_hist.SetTitle(projection_name)
                        h1_input_hist.SetDirectory(0)
                    
                        # Undo the bin width scaling and set correct normalization
                        norm_factor = h1_input_hist.Integral()
                        if norm_factor == 0: norm_factor = 1
                        h1_input_hist.Scale(1./norm_factor, "width")
                    
                        for b in range(0,h1_input_hist.GetNbinsX()):
                            h1_input_hist.SetBinError(b+1,0)

                        histo_list.append(h1_input_hist)

                    # Create envelope histograms
                    hist_min, hist_max = self.min_max( histo_list )

                    # Rename some objects
                    name_h_cent = 'h1_input_%s_R%s_%s_pT_%i_%i'     % ( self.observable,(str)(jetR).replace('.',''),new_obs_lab,(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]))
                    name_h_min  = 'h1_min_input_%s_R%s_%s_pT_%i_%i' % ( self.observable,(str)(jetR).replace('.',''),new_obs_lab,(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]))
                    name_h_max  = 'h1_max_input_%s_R%s_%s_pT_%i_%i' % ( self.observable,(str)(jetR).replace('.',''),new_obs_lab,(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]))

                    h_central = histo_list[0]
                    h_central.SetName(name_h_cent)
                    hist_min .SetName(name_h_min )
                    hist_max .SetName(name_h_max )

                    # Create a graph out of these histograms
                    graph_cent = self.histo_to_graph(h_central,hist_min,hist_max)
                    graph_min  = ROOT.TGraph(hist_min)
                    graph_max  = ROOT.TGraph(hist_max)
                    graph_frac = self.fractional_error(h_central,hist_min,hist_max)

                    graph_cent.SetName('g_input_%s_R%s_%s_pT_%i_%i'     % ( self.observable,(str)(jetR).replace('.',''),new_obs_lab,(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
                    graph_min .SetName('g_min_input_%s_R%s_%s_pT_%i_%i' % ( self.observable,(str)(jetR).replace('.',''),new_obs_lab,(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
                    graph_max .SetName('g_max_input_%s_R%s_%s_pT_%i_%i' % ( self.observable,(str)(jetR).replace('.',''),new_obs_lab,(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
                    graph_frac.SetName('g_frac_input_%s_R%s_%s_pT_%i_%i'% ( self.observable,(str)(jetR).replace('.',''),new_obs_lab,(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))

                    xtit = self.obs_label
                    ytit = '#frac{1}{#sigma} #frac{d#sigma}{d'+xtit+'}'
                    tit = 'input (hadron-level, no MPI) %i < #it{p}_{T}^{jet} < %i GeV/#it{c}'%((int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]))

                    self.pretty_1D_object(graph_cent,2,2,1,tit, xtit, ytit, True)
                    self.pretty_1D_object(graph_min ,1,1,2,tit, xtit, ytit)
                    self.pretty_1D_object(graph_max ,1,1,2,tit, xtit, ytit)
                    self.pretty_1D_object(graph_frac,2,2,1,tit, xtit, ytit, True)

                    outpdfname = os.path.join(self.output_dir, 'control_plots' , 'processed_plots' )
                    if not os.path.exists(outpdfname):
                        os.makedirs(outpdfname)
                    outpdfname_1 = os.path.join(outpdfname, 'theory_%s_pT_%i_%i_GeVc_input.pdf'%(self.create_label(jetR,obs_setting,grooming_setting),(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])) )
                    self.plot_processed_functions( graph_cent, graph_min, graph_max, outpdfname_1)

                    # loop over response files (e.g. Pythia, Herwig, ...)
                    for ri,_ in enumerate(self.theory_response_files):
                        for lev in self.response_levels:
                            outpdfname_2 = os.path.join(outpdfname, 'comp_gen_input_theory_%s_pT_%i_%i_GeVc_'%(self.create_label(jetR,obs_setting,grooming_setting),(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])) )
                            outpdfname_2 += lev[0]+"_"+lev[1]+"_MPI"+lev[2]+"_"+self.theory_response_labels[ri]+".pdf" 
                            self.plot_comparison_SCET_gen_input( graph_cent, jetR , obs_setting , grooming_setting, lev[0], lev[1], lev[2], self.theory_response_labels[ri], self.final_pt_bins[n_pt], self.final_pt_bins[n_pt+1], outpdfname_2)

                    self.outfile.cd()
                    h_central .Write()
                    hist_min  .Write()
                    hist_max  .Write()
                    graph_cent.Write()
                    graph_min .Write()
                    graph_max .Write()
                    graph_frac.Write()

                # -----------------------------------------------------
                # Setting the filled histograms as attributes
                self.outfile.cd()
                for sv in range(0,n_scale_variations):
                    setattr(self,hist_names[sv],th_hists[sv])
                
                outpdfname = os.path.join(self.output_dir, 'control_plots' , 'input' )
                if not os.path.exists(outpdfname):
                    os.makedirs(outpdfname)
                outpdfname = os.path.join(outpdfname, 'theory_input_%s.pdf'%(label) )
                self.plot_input_theory( th_hists_no_scaling , th_hists , outpdfname )

                scale_var.append(n_scale_variations)
            self.theory_scale_vars[jetR] = scale_var

    #----------------------------------------------------------------------
    # Plot input theory curves both as histograms and curves
    #----------------------------------------------------------------------
    def plot_input_theory( self , h_list_no_scaling , h_list , outpdfname ):

        for i in range(0,len(h_list_no_scaling)):

            c1 = ROOT.TCanvas('c1','c1',1000,800)
            c1.Divide(2,2)
        
            for j in range(0,4):
                c1.cd(j+1).SetLogz()
                c1.cd(j+1).SetLeftMargin(0.20)
                if j > 1:
                    c1.cd(j+1).SetTheta(50)
                    c1.cd(j+1).SetPhi(220)
                else:
                    c1.cd(j+1).SetBottomMargin(0.20)
                    c1.cd(j+1).SetRightMargin(0.24)

            self.pretty_TH2D(h_list_no_scaling[i],'input theory curves, scale variation %i'%(i),'#it{p}_{T}^{jet} [GeV/#it{c}]',self.obs_label,'#frac{1}{#sigma} #frac{d#sigma}{d'+self.obs_label+'}')
            self.pretty_TH2D(h_list           [i],'scaled input, scale variation %i'%(i)       ,'#it{p}_{T}^{jet} [GeV/#it{c}]',self.obs_label,'~ #frac{d#sigma}{d'+self.obs_label+'}')
            
            h_list_no_scaling[i].GetXaxis().SetTitleOffset(1.6)
            h_list           [i].GetXaxis().SetTitleOffset(1.6)
            h_list_no_scaling[i].GetYaxis().SetTitleOffset(1.5)
            h_list           [i].GetYaxis().SetTitleOffset(1.5)
            h_list_no_scaling[i].GetZaxis().SetTitleOffset(1.4)
            h_list           [i].GetZaxis().SetTitleOffset(1.4)

            c1.cd(1)
            h_list_no_scaling[i].Draw('COLZ')
            
            c1.cd(2)
            h_list[i].Draw('COLZ')
        
            c1.cd(3)
            h_list_no_scaling[i].Draw('LEGO1')

            c1.cd(4)
            h_list[i].Draw('LEGO1')
        
            c1.Draw()   

            if len(h_list_no_scaling)==1:
                c1.Print(outpdfname)
            else: 
                if i == 0: c1.Print(outpdfname+'(')
                elif i == len(h_list_no_scaling)-1:  c1.Print(outpdfname+')')
                else:  c1.Print(outpdfname)

            del c1
    
#----------------------------------------------------------------------
if __name__ == '__main__':

  # Define arguments
  parser = argparse.ArgumentParser(description='Folding theory predictions')
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

  analysis = TheoryFolding(config_file = args.configFile)
  analysis.run_theory_folding()