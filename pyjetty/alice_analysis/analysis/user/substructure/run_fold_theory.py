#! /usr/bin/env python
"""
run_fold_theory.py
Code to fold theory curves from a given to a desired 'level'
Adapted from Ezra Lesser's code by Rey Cruz-Torres (reynier@lbl.gov)
"""

import sys
import os
import argparse
from array import *
import numpy as np
import ROOT
ROOT.gSystem.Load("$HEPPY_DIR/external/roounfold/roounfold-current/lib/libRooUnfold.so")
import yaml

from pyjetty.alice_analysis.analysis.user.substructure import analysis_utils_obs

# Load pyjetty ROOT utils
ROOT.gSystem.Load('libpyjetty_rutil')

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 2000

################################################################
################################################################
################################################################
class TheoryFolding():
  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(TheoryFolding, self).__init__(**kwargs)

    # Initialize utils class
    self.utils = analysis_utils_obs.AnalysisUtils_Obs()
    self.config_file = config_file

    # Initialize yaml config
    self.initialize_user_config()

    ROOT.gStyle.SetOptStat(0)
  
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_user_config(self):
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
      
    self.histutils = ROOT.RUtil.HistUtils()

    if 'theory_dir' in config:
      self.jetR_list = config['jetR']
      self.observable = config['th_fold_observable']
      self.obs_config_dict = config[self.observable]

      # If the user specifies certain subconfigurations to fold via the th_subconfigs parameter,
      # only fold those. Otherwise, assume we want to ufold all subconfigs in the config file
      if 'th_subconfigs' in config:
        self.obs_subconfig_list = config['th_subconfigs']
      else:
        self.obs_subconfig_list = [name for name in list(self.obs_config_dict.keys()) if 'config' in name ]

      self.obs_settings = self.utils.obs_settings(self.observable, self.obs_config_dict, self.obs_subconfig_list) 
      self.grooming_settings = [ self.get_grooming_setting(self.obs_config_dict[cf]) for cf in self.obs_subconfig_list]
      #self.obs_labels = [self.utils.obs_label(self.obs_settings[i], self.grooming_settings[i]) for i in range(len(self.obs_subconfig_list))]

      # this corresponds to a label in the RM name which indicates whether we are going from parton to charged-hadron level, or
      # from hadron to charged-hadron level, ...
      self.response_levels = config['response_levels']

      self.theory_dir = config['theory_dir']

      subobs_binning_in_config = True
      self.subobs_specific_binning = {}
      for R,jetR in enumerate(self.jetR_list):
        temp_bin_dic = {}
        for conf in self.obs_subconfig_list:
          if 'obs_bins_truth' in config[self.observable][conf]:
            list1 = config[self.observable][conf]['obs_bins_truth']
            if 'obs_max_reported' in config[self.observable][conf]:
              temp_bin_dic[conf] = list1[:list1.index(config[self.observable][conf]['obs_max_reported'][R])+1]
            else:
              temp_bin_dic[conf] = list1
          else:
           subobs_binning_in_config = False
        self.subobs_specific_binning[jetR] = temp_bin_dic

      # binning
      self.theory_obs_bins = None
      use_nominal_ana_binning = False
      if 'use_nominal_ana_binning' in config:
        use_nominal_ana_binning = config['use_nominal_ana_binning']

      #---------------------------------
      # Specify which binning to use for the observable
      self.obs_bin_option = None
      if 'theory_obs_bins' in config and not use_nominal_ana_binning:
        self.theory_obs_bins = config['theory_obs_bins']
        print('Will use binning defined in variable theory_obs_bins')
        self.obs_bin_option = 'theory_obs_binning'
      elif not 'theory_obs_bins' in config and use_nominal_ana_binning and subobs_binning_in_config:
        print('Will use binning used in the analysis')
        self.obs_bin_option = 'nominal_ana_binning'
      elif 'theory_obs_bins' in config and use_nominal_ana_binning:
        print('User requested both binning used in the analysis and binning stored in theory_obs_bins. Conflict! Will default to the latter')
        self.obs_bin_option = 'theory_obs_binning'
        self.theory_obs_bins = config['theory_obs_bins']
      elif use_nominal_ana_binning and not subobs_binning_in_config:
        print('User requested observable binning used in the analysis, but one or more set of bins were not found in config file. Defaulting to RM binning.')
        self.obs_bin_option = 'response_binning'
      else:
        print('No observable binning was specified, so will take whichever binning the RM comes with.')
        self.obs_bin_option = 'response_binning'
      #---------------------------------

      self.theory_pt_bins  = config['theory_pt_bins' ] # pT binning of theory calculations
      self.final_pt_bins   = config['final_pt_bins'  ] # pT binning wanted for the final curves

      # response matrices for the folding, and labels describing them
      self.theory_response_fname = config['response_files']
      self.theory_response_files = [ROOT.TFile(f, 'READ') for f in config['response_files']]
      self.theory_response_labels = config['response_labels']

      # scale factors needed to scale distributions
      self.theory_pt_scale_factors_filepath = os.path.join(self.theory_dir, config['pt_scale_factors_path'])

      self.use_tagging_fraction = False
      if 'use_tagging_fraction' in config:
        self.use_tagging_fraction = config['use_tagging_fraction']

      #self.output_dir = config['output_dir']
      self.output_dir = config['theory_dir']
      self.output_dir_theory = os.path.join(self.output_dir, self.observable, 'theory_response') 

      self.set_observable_label()
    else:
      print('Missing several parameters in the config file!')
      exit()

  #---------------------------------------------------------------
  # Set Observable Label
  def set_observable_label(self):
    obs = ''
    if self.observable == 'jet_axis':
      obs = '#it{#DeltaR}_{axis}'
    elif self.observable == 'ang':
      obs = '#it{#lambda}_{#it{#beta}}^{#it{#kappa}=1}'
    else:
      obs = 'observable'
    
    self.obs_label = obs

  #---------------------------------------------------------------

  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def run_theory_folding(self):
      # Creating a root file to store results
      outfilename = 'folded_scet_pT_{}_{}.root'.format(self.final_pt_bins[0],self.final_pt_bins[-1])
      if len(self.jetR_list) == 1:
        outfilename = 'folded_scet_R{}_pT_{}_{}.root'.format(((str)(self.jetR_list[0])).replace('.',''),self.final_pt_bins[0],self.final_pt_bins[-1])
      outfilename = os.path.join( self.theory_dir , outfilename )
      self.outfile = ROOT.TFile(outfilename,'recreate')
      print('===========================================================')
      print('Loading pT scale factors...')
      self.load_pt_scale_factors(self.theory_pt_scale_factors_filepath)
      print('===========================================================')
      print('Loading response matrix for folding theory predictions...')
      self.load_theory_response()
      print('===========================================================')
      print('Loading theory curves...')
      self.load_theory_curves()
      print('===========================================================')
      print("Folding theory histograms...")
      self.fold_theory()
      print('===========================================================')
      print("Applying some final scalings...")
      self.final_processing()
      # ------------
      # Closing the root file with all results from this code
      self.outfile.Close()
      print('***********************************************************************************************\nDone!')
      print('output produced by this code can be found in:',self.output_dir)
      print('***********************************************************************************************')

  #---------------------------------------------------------------
  # Loads pT scale factors from q/g fraction theory predictions
  #---------------------------------------------------------------
  def load_pt_scale_factors(self, filepath):

    for i, jetR in enumerate(self.jetR_list):
      full_path = os.path.join(filepath,'qg_fractions-ALICE-R%s.txt' % ((str)(jetR).replace('.','')) )
      # Open file and save pT distribution
      pt_li = None; val_li = None;
      with open(full_path) as f:
        lines = [line.split() for line in f.read().split('\n') if (line and line[0] != '#')]
        pt_li = [int(float(line[0])) for line in lines]
        val_li = [float(line[1]) + float(line[2]) for line in lines]

      n_entries = len(val_li)
      val_li_jetR = val_li[i*n_entries:(i+1)*n_entries]
      pt_li_jetR = pt_li[i*n_entries:(i+1)*n_entries]
      setattr(self, "pt_scale_factor_R%s" % jetR, (pt_li_jetR, val_li_jetR))

  #---------------------------------------------------------------
  # Load theory calculations
  #---------------------------------------------------------------
  def load_theory_curves(self):
    # The user needs to implement this function
    raise NotImplementedError('You must implement initialize_user_output_objects()!')

    # The theory curves are given as 1/σ dσ/d(obs). Before doing the folding,
    # we need to convert these quantities to dσ/d(obs). To do this, we scale
    # the theory curves by scale factors (loaded in the previous function)

  #---------------------------------------------------------------
  # Linear interpolation
  #---------------------------------------------------------------
  def interpolate_values_linear(self,x_val,y_val_n,new_bin_edges):
    # given the input bin edges, determine the bin center by taking the average
    bin_ctr = np.add(new_bin_edges[1:],new_bin_edges[:-1])/2.
    # do the interpolation
    return np.interp(bin_ctr,x_val,y_val_n,left=0,right=0,period=None)

  #---------------------------------------------------------------
  # Load 4D response matrices used for folding the theory predictions
  #---------------------------------------------------------------
  def load_theory_response(self):

    # Check to see if Roounfold file already exists
    if not os.path.exists(self.output_dir_theory):
      os.makedirs(self.output_dir_theory)
    roounfold_filename = os.path.join(self.output_dir_theory, 'fRoounfold.root')

    # Loop over jet R
    for jetR in self.jetR_list:
      # Loop through subconfigurations to fold
      # (e.g. Standard_WTA, Standard_SD_1, ... in the jet-axis analysis, or
      # beta = 1.5, 2, 3, ... in the angularities analysis)
      for i,conf in enumerate(self.obs_subconfig_list):

        obs_setting = self.obs_settings[i]           # labels such as 'Standard_WTA'
        grooming_setting = self.grooming_settings[i] # grooming parameters
        label = self.create_label( jetR , obs_setting , grooming_setting ) 
        
        # loop over response files (e.g. Pythia, Herwig, ...)
        for ri, response in enumerate(self.theory_response_files):

          for li, lev in enumerate(self.response_levels):
            
            # Load response matrix 
            name_RM = "hResponse_JetPt_"+self.observable+"_"+lev[0]+"_"+lev[1]+"_MPI"+lev[2]+"_"+label
            thn = response.Get(name_RM)
            if thn == None:
              print('Could not find RM:',name_RM,'in',self.theory_response_fname[ri])
              exit()
            thn.SetName(name_RM+'_{}'.format(self.theory_response_labels[ri]))
          
            print('Loading',self.theory_response_labels[ri],'RM:',name_RM)

            # Create Roounfold object
            name_roounfold_obj = '%s_Roounfold_%s' % (name_RM, self.theory_response_labels[ri])
            name_roounfold_thn = '%s_Rebinned_%s'  % (name_RM, self.theory_response_labels[ri])
          
            '''
            Response axes:
            ['p_{T}^{final}', 'p_{T}^{initial}', 'obs^{final}', 'obs^{initial}']
            e.g. ['p_{T}^{ch}', 'p_{T}^{h}', 'obs^{ch}', 'obs^{h}']
            
            as compared to the usual
            ['p_{T}^{initial}', 'p_{T}^{final}', 'obs^{initial}', 'obs^{final}']
            e.g. ['p_{T}^{det}', 'p_{T}^{truth}', 'obs^{det}', 'obs_{truth}']
            '''

            # Regardless of what the user asks for the binning,
            # default to the RM binning at input level, so that
            # the input level has as fine binning as possible
            RM_binning = self.return_histo_binning_1D(thn.Projection(3))
            RM_binning = [l for l in RM_binning if l >= 0]
            setattr(self,'binning_R%s_'%((str)(jetR).replace('.',''))+str(obs_setting), RM_binning)

            # If no binning was specified by the user, take the RM binning
            if self.obs_bin_option == 'nominal_ana_binning':
              binning = self.subobs_specific_binning[jetR][conf]
            else:
              if ri == 0 and li == 0:
                if self.obs_bin_option == 'response_binning':
                  binning = RM_binning
                elif self.obs_bin_option == 'theory_obs_binning':
                  binning = self.theory_obs_bins

            det_pt_bin_array = array('d', self.theory_pt_bins)
            tru_pt_bin_array = det_pt_bin_array
            det_obs_bin_array = array('d', binning   ) # Binning we are folding to
            tru_obs_bin_array = array('d', RM_binning) # Binning we are folding from (Use the RM binning at input level)
          
            if grooming_setting and self.use_tagging_fraction:
              # Add bin for underflow value (tagging fraction)
              det_obs_bin_array = np.insert(det_obs_bin_array, 0, -0.001)
              tru_obs_bin_array = det_obs_bin_array
          
            n_dim = 4
            self.histutils.rebin_thn( roounfold_filename, thn, name_roounfold_thn , name_roounfold_obj, n_dim,
              len(det_pt_bin_array )-1, det_pt_bin_array ,
              len(det_obs_bin_array)-1, det_obs_bin_array,
              len(tru_pt_bin_array )-1, tru_pt_bin_array ,
              len(tru_obs_bin_array)-1, tru_obs_bin_array,
              label,0,1) #, grooming_setting!=None )
          
            f_resp = ROOT.TFile(roounfold_filename, 'READ')
            roounfold_response = f_resp.Get(name_roounfold_obj)
            roounfold_thn      = f_resp.Get(name_roounfold_thn )
            f_resp.Close() 
          
            setattr(self, name_roounfold_obj, roounfold_response)
            setattr(self, name_roounfold_thn, roounfold_thn     ) 

            # Save the response matrix to pdf file
            outpdfname = os.path.join(self.output_dir, 'control_plots', 'RM_and_MPI' )
            if not os.path.exists(outpdfname):
              os.makedirs(outpdfname)
            outpdfname = os.path.join(outpdfname, 'RM_slices_%s.pdf'%(label) )
            self.plot_RM_slice_histograms( roounfold_thn, outpdfname)

            #self.outfile.cd()
            #roounfold_thn.Write()

  #----------------------------------------------------------------------
  # Extract binning from a 1D histogram
  #----------------------------------------------------------------------
  def return_histo_binning_1D(self,h1):
    nBins = h1.GetNbinsX()
    binning = []
    for b in range(0,nBins):
      binning.append(h1.GetBinLowEdge(b+1))
      if b == nBins-1:
        binning.append(h1.GetBinLowEdge(b+1)+h1.GetBinWidth(b+1))
    return array('d',binning)

  #----------------------------------------------------------------------
  # Fold theoretical predictions
  #----------------------------------------------------------------------
  def fold_theory(self):

    # Loop over jet R
    for jetR in self.jetR_list:
     # Loop through subconfigurations to fold
     # (e.g. Standard_WTA, Standard_SD_1, ... in the jet-axis analysis, or
     # beta = 1.5, 2, 3, ... in the angularities analysis)
     for i, subconfig in enumerate(self.obs_subconfig_list):

       obs_setting = self.obs_settings[i]
       grooming_setting = self.grooming_settings[i]
       label = self.create_label( jetR , obs_setting , grooming_setting ) 

       # Retrieve theory histograms to be folded
       th_hists = []

       for sv in range(0,self.theory_scale_vars[jetR][i]):
         hist_name = 'h2_input_%s_R%s_obs_pT_%s' % ( self.observable , (str)(jetR).replace('.','') , obs_setting )
         if grooming_setting:
           hist_name += '_'
           hist_name += self.utils.grooming_label(grooming_setting)
         hist_name += '_sv%i' % (sv)
 
         th_hist = getattr(self,hist_name)
         th_hists.append(th_hist) 

       # loop over response files (e.g. Pythia, Herwig, ...)
       for ri, response in enumerate(self.theory_response_files):

         for lev in self.response_levels:
           name_RM = "hResponse_JetPt_"+self.observable+"_"+lev[0]+"_"+lev[1]+"_MPI"+lev[2]+"_"+label
           name_roounfold_obj = '%s_Roounfold_%s' % (name_RM, self.theory_response_labels[ri])
           response = getattr(self,name_roounfold_obj)

           for sv in range(0,self.theory_scale_vars[jetR][i]):
             h_folded_ch = response.ApplyToTruth(th_hists[sv])
             folded_hist_name = 'h2_folded_'+self.observable+'_'+lev[1]+"_MPI"+lev[2]+'_R%s_obs_pT_%s_%s_sv%i' % ((str)(jetR).replace('.','') , obs_setting , self.theory_response_labels[ri] , sv )
             h_folded_ch.SetNameTitle(folded_hist_name,folded_hist_name)

             setattr(self, folded_hist_name, h_folded_ch)

  #----------------------------------------------------------------------
  # Undoing some scalings that had been introduced before
  #----------------------------------------------------------------------
  def final_processing(self):

    # Loop over jet R
    for jetR in self.jetR_list:
      # Loop through subconfigurations to fold
      # (e.g. Standard_WTA, Standard_SD_1, ... in the jet-axis analysis, or
      # beta = 1.5, 2, 3, ... in the angularities analysis)
      for i, subconfig in enumerate(self.obs_subconfig_list):

        obs_setting = self.obs_settings[i]
        grooming_setting = self.grooming_settings[i]
        label = self.create_label( jetR , obs_setting , grooming_setting ) 

        pt_bins = array('d', self.theory_pt_bins)

        # loop over response files (e.g. Pythia, Herwig, ...)
        for ri, response in enumerate(self.theory_response_files):
         
          # ----------------------------------------------------------------------------------------
          # Preparing MPI correction 
          # Grab the two histograms that will be used for the MPI correction
          name_mpi_off = 'h_JetPt_'+self.observable+'_ch_MPIoff_'+label
          name_mpi_on = 'h_JetPt_'+self.observable+'_ch_MPIon_'+label
 
          h2_mpi_off = response.Get(name_mpi_off)
          h2_mpi_off.SetName(h2_mpi_off.GetName()+'_{}'.format(self.theory_response_labels[ri]))
          h2_mpi_on = response.Get(name_mpi_on)
          h2_mpi_on.SetName(h2_mpi_on.GetName()+'_{}'.format(self.theory_response_labels[ri]))

          # Gotta make sure the histograms we will use for the correction have the proper binning
          if self.obs_bin_option == 'nominal_ana_binning':
            y_bins = array('d', self.subobs_specific_binning[jetR][subconfig])
          elif self.obs_bin_option == 'theory_obs_binning':
            y_bins = array('d', self.theory_obs_bins)
          else:
            y_bins = array('d',getattr(self,'binning_R%s_'%((str)(jetR).replace('.',''))+str(obs_setting)))
  
          if grooming_setting and self.use_tagging_fraction:
            y_bins = np.insert(y_bins, 0, -0.001)
 
          xtit_MPIoff = h2_mpi_off.GetXaxis().GetTitle()
          ytit_MPIoff = h2_mpi_off.GetYaxis().GetTitle()
          
          h2_mpi_off = self.histutils.rebin_th2(h2_mpi_off, name_mpi_off+'_Rebinned_%s' % self.theory_response_labels[ri], pt_bins, len(pt_bins)-1, y_bins, len(y_bins)-1, grooming_setting!=None )
          h2_mpi_on  = self.histutils.rebin_th2(h2_mpi_on , name_mpi_on +'_Rebinned_%s' % self.theory_response_labels[ri], pt_bins, len(pt_bins)-1, y_bins, len(y_bins)-1, grooming_setting!=None )
          
          h2_mpi_off.GetXaxis().SetTitle(xtit_MPIoff)
          h2_mpi_off.GetYaxis().SetTitle(ytit_MPIoff)
          h2_mpi_on .GetXaxis().SetTitle(xtit_MPIoff) # The titles should be the same in both histograms
          h2_mpi_on .GetYaxis().SetTitle(ytit_MPIoff) # The titles should be the same in both histograms
 
          h2_mpi_ratio = h2_mpi_on.Clone()
          title = 'h_mpi_on_over_off_'+self.observable+'_JetPt_ch_'+label
          h2_mpi_ratio.SetNameTitle(title,title)
          h2_mpi_ratio.Divide(h2_mpi_off)
          h2_mpi_ratio.SetDirectory(0)
 
          # ----------------------------------------------------------------------------------------
          # Loop over response levels
          for lev in self.response_levels:
            # Loop over scale variations
            for sv in range(0,self.theory_scale_vars[jetR][i]):
               
              folded_hist_name = 'h2_folded_'+self.observable+'_'+lev[1]+"_MPI"+lev[2]+'_R%s_obs_pT_%s_%s_sv%i' % ((str)(jetR).replace('.','') , obs_setting , self.theory_response_labels[ri] , sv )
              h2_folded_hist = getattr(self,folded_hist_name)
 
              if lev[2]=='off':
                yesMPI_hist_name = 'h2_folded_'+self.observable+"_"+lev[1]+"_MPIon"+'_R%s_obs_pT_%s_%s_sv%i' % ((str)(jetR).replace('.','') , obs_setting , self.theory_response_labels[ri] , sv )
                h2_folded_hist.SetNameTitle(yesMPI_hist_name,yesMPI_hist_name)
         
                # Copy that won't have MPI corrections
                h2_folded_hist_noMPI = h2_folded_hist.Clone() 
                noMPI_hist_name = 'h2_folded_'+self.observable+"_"+lev[1]+"_MPI"+lev[2]+'_R%s_obs_pT_%s_%s_sv%i' % ((str)(jetR).replace('.','') , obs_setting , self.theory_response_labels[ri] , sv )
                h2_folded_hist_noMPI.SetNameTitle(noMPI_hist_name,noMPI_hist_name)
  
                h2_folded_hist.Multiply(h2_mpi_ratio)
           
              self.outfile.cd()
              h2_folded_hist.Write()
              if lev[2]=='off':
                h2_folded_hist_noMPI.Write()
           
              # If desired binning is different from what was used for the folding, need to take that into account before changing the pT normalization
              for n_pt in range(0,len(self.final_pt_bins)-1):
                # Get the bins that correspond to the pT edges given
                min_bin, max_bin = self.bin_position( self.theory_pt_bins, self.final_pt_bins[n_pt], self.final_pt_bins[n_pt+1] )
           
                projection_name = 'h1_folded_%s_%s_MPIon_R%s_%s_%s_sv%i_pT_%i_%i' % ( self.observable,lev[1],(str)(jetR).replace('.',''),obs_setting,self.theory_response_labels[ri],sv,(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]))            
                h1_folded_hist = h2_folded_hist.ProjectionY(projection_name,min_bin,max_bin)
                h1_folded_hist.SetTitle(projection_name)
                h1_folded_hist.SetDirectory(0)
           
                # Undo the bin width scaling and set correct normalization
                norm_factor = h1_folded_hist.Integral()
                if norm_factor == 0: norm_factor = 1
                h1_folded_hist.Scale(1./norm_factor, "width")
           
                for b in range(0,h1_folded_hist.GetNbinsX()):
                  h1_folded_hist.SetBinError(b+1,0)
 
                setattr(self,projection_name, h1_folded_hist)          
 
                if lev[2]=='off':
                  # Now doing the same, for the histograms with no MPI corrections
                  projection_name_noMPI = 'h1_folded_%s_%s_MPIoff_R%s_%s_%s_sv%i_pT_%i_%i' % ( self.observable,lev[1],(str)(jetR).replace('.',''),obs_setting,self.theory_response_labels[ri],sv,(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]))
                  h1_folded_hist_noMPI = h2_folded_hist_noMPI.ProjectionY(projection_name_noMPI,min_bin,max_bin)
                  h1_folded_hist_noMPI.SetTitle(projection_name_noMPI)
                  h1_folded_hist_noMPI.SetDirectory(0)
           
                  # Undo the bin width scaling and set correct normalization
                  norm_factor_noMPI = h1_folded_hist_noMPI.Integral()
                  if norm_factor_noMPI == 0: norm_factor_noMPI = 1
                  h1_folded_hist_noMPI.Scale(1./norm_factor_noMPI, "width")
           
                  for b in range(0,h1_folded_hist_noMPI.GetNbinsX()):
                    h1_folded_hist_noMPI.SetBinError(b+1,0)
          
                  setattr(self,projection_name_noMPI,h1_folded_hist_noMPI)
           
            new_obs_lab = obs_setting
            if grooming_setting:
              new_obs_lab += '_'
              new_obs_lab += self.utils.grooming_label(grooming_setting)
           
            # Do the loop backwards and find min and max histograms
            for n_pt in range(0,len(self.final_pt_bins)-1):
              histo_list = []
              histo_list_noMPI = []
              for sv in range(0,self.theory_scale_vars[jetR][i]):
                projection_name = 'h1_folded_%s_%s_MPIon_R%s_%s_%s_sv%i_pT_%i_%i' % ( self.observable,lev[1],(str)(jetR).replace('.',''),obs_setting,self.theory_response_labels[ri],sv,(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]))   
                histo_list.append(getattr(self,projection_name)) 
 
                if lev[2]=='off':          
                  projection_name_noMPI = 'h1_folded_%s_%s_MPIoff_R%s_%s_%s_sv%i_pT_%i_%i' % ( self.observable,lev[1],(str)(jetR).replace('.',''),obs_setting,self.theory_response_labels[ri],sv,(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]))
                  histo_list_noMPI.append(getattr(self,projection_name_noMPI))
           
              hist_min, hist_max = self.min_max( histo_list )
           
              # Create a graph out of these histograms
              name_central = 'h1_folded_%s_%s_MPIon_R%s_%s_%s_sv0_pT_%i_%i' % ( self.observable,lev[1],(str)(jetR).replace('.',''),obs_setting,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]))
              h_central = getattr(self,name_central)
              graph_cent = self.histo_to_graph(h_central,hist_min,hist_max)
              graph_min = ROOT.TGraph(hist_min)
              graph_max = ROOT.TGraph(hist_max)
              graph_frac = self.fractional_error(h_central,hist_min,hist_max)
           
              h_central.SetName('h1_folded_%s_%s_MPIon_R%s_%s_%s_pT_%i_%i'     % ( self.observable,lev[1],(str)(jetR).replace('.',''),new_obs_lab,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
              hist_min .SetName('h1_min_folded_%s_%s_MPIon_R%s_%s_%s_pT_%i_%i' % ( self.observable,lev[1],(str)(jetR).replace('.',''),new_obs_lab,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
              hist_max .SetName('h1_max_folded_%s_%s_MPIon_R%s_%s_%s_pT_%i_%i' % ( self.observable,lev[1],(str)(jetR).replace('.',''),new_obs_lab,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
           
              graph_cent.SetName('g_folded_%s_%s_MPIon_R%s_%s_%s_pT_%i_%i'     % ( self.observable,lev[1],(str)(jetR).replace('.',''),new_obs_lab,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
              graph_min .SetName('g_min_folded_%s_%s_MPIon_R%s_%s_%s_pT_%i_%i' % ( self.observable,lev[1],(str)(jetR).replace('.',''),new_obs_lab,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
              graph_max .SetName('g_max_folded_%s_%s_MPIon_R%s_%s_%s_pT_%i_%i' % ( self.observable,lev[1],(str)(jetR).replace('.',''),new_obs_lab,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
              graph_frac.SetName('g_frac_folded_%s_%s_MPIon_R%s_%s_%s_pT_%i_%i'% ( self.observable,lev[1],(str)(jetR).replace('.',''),new_obs_lab,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
           
              # Now doing the same, for the histograms with no MPI corrections
              if lev[2]=='off':
                hist_min_noMPI, hist_max_noMPI  = self.min_max( histo_list_noMPI )
                name_central_noMPI = 'h1_folded_%s_%s_MPIoff_R%s_%s_%s_sv0_pT_%i_%i' % ( self.observable,lev[1],(str)(jetR).replace('.',''),obs_setting,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]))
                h_central_noMPI = getattr(self,name_central_noMPI)
                graph_cent_noMPI = self.histo_to_graph(h_central_noMPI,hist_min_noMPI, hist_max_noMPI)
                graph_min_noMPI = ROOT.TGraph(hist_min_noMPI)
                graph_max_noMPI = ROOT.TGraph(hist_max_noMPI)
                graph_frac_noMPI = self.fractional_error(h_central_noMPI,hist_min_noMPI, hist_max_noMPI)
           
                h_central_noMPI.SetName('h1_folded_%s_%s_MPIoff_R%s_%s_%s_pT_%i_%i'     % ( self.observable,lev[1],(str)(jetR).replace('.',''),new_obs_lab,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
                hist_min_noMPI .SetName('h1_min_folded_%s_%s_MPIoff_R%s_%s_%s_pT_%i_%i' % ( self.observable,lev[1],(str)(jetR).replace('.',''),new_obs_lab,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
                hist_max_noMPI .SetName('h1_max_folded_%s_%s_MPIoff_R%s_%s_%s_pT_%i_%i' % ( self.observable,lev[1],(str)(jetR).replace('.',''),new_obs_lab,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
           
                graph_cent_noMPI.SetName('g_folded_%s_%s_MPIoff_R%s_%s_%s_pT_%i_%i'     % ( self.observable,lev[1],(str)(jetR).replace('.',''),new_obs_lab,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
                graph_min_noMPI .SetName('g_min_folded_%s_%s_MPIoff_R%s_%s_%s_pT_%i_%i' % ( self.observable,lev[1],(str)(jetR).replace('.',''),new_obs_lab,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
                graph_max_noMPI .SetName('g_max_folded_%s_%s_MPIoff_R%s_%s_%s_pT_%i_%i' % ( self.observable,lev[1],(str)(jetR).replace('.',''),new_obs_lab,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
                graph_frac_noMPI.SetName('g_frac_folded_%s_%s_MPIoff_R%s_%s_%s_pT_%i_%i'% ( self.observable,lev[1],(str)(jetR).replace('.',''),new_obs_lab,self.theory_response_labels[ri],(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1])))
           
              xtit = self.obs_label
              ytit = '#frac{1}{#sigma} #frac{d#sigma}{d'+xtit+'}'
              tit_noMPI = 'output (charged-level, MPI off) %i < #it{p}_{T}^{jet} < %i GeV/#it{c}'%((int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]))
              tit = 'output (charged-level, MPI on) %i < #it{p}_{T}^{jet} < %i GeV/#it{c}'%((int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]))
              
              if lev[2]=='off':
                self.pretty_1D_object(graph_cent_noMPI, 8,2,1,tit_noMPI, xtit, ytit, True)
                self.pretty_1D_object(graph_min_noMPI , 1,1,2,tit_noMPI, xtit, ytit)
                self.pretty_1D_object(graph_max_noMPI , 1,1,2,tit_noMPI, xtit, ytit)
                self.pretty_1D_object(graph_frac_noMPI, 8,2,1,tit_noMPI, xtit,'frac', True)
           
              self.pretty_1D_object(graph_cent      ,62,2,1,tit      , xtit, ytit, True)
              self.pretty_1D_object(graph_min       , 1,1,2,tit      , xtit, ytit)
              self.pretty_1D_object(graph_max       , 1,1,2,tit      , xtit, ytit)
              self.pretty_1D_object(graph_frac      ,62,2,1,tit      , xtit,'frac', True)
           
              outpdfname = os.path.join(self.output_dir, 'control_plots' , 'processed_plots' )
              if not os.path.exists(outpdfname):
                os.makedirs(outpdfname)
 
              if lev[2]=='off':
                outpdfname_noMPI = os.path.join(outpdfname, 'theory_%s_pT_%i_%i_GeVc_output_noMPI_%s.pdf'%(label,(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]),self.theory_response_labels[ri]) )
                self.plot_processed_functions( graph_cent_noMPI, graph_min_noMPI, graph_max_noMPI, outpdfname_noMPI )
 
              outpdfname       = os.path.join(outpdfname, 'theory_%s_pT_%i_%i_GeVc_output_%s.pdf'      %(label,(int)(self.final_pt_bins[n_pt]),(int)(self.final_pt_bins[n_pt+1]),self.theory_response_labels[ri]) )
              self.plot_processed_functions( graph_cent      , graph_min      , graph_max      , outpdfname       )
           
              # Saving results to root file
              self.outfile.cd()
 
              if lev[2]=='off':          
                h_central_noMPI .Write()
                hist_min_noMPI  .Write()
                hist_max_noMPI  .Write()
                graph_cent_noMPI.Write()
                graph_min_noMPI .Write()
                graph_max_noMPI .Write()
                graph_frac_noMPI.Write()
           
              h_central .Write()
              hist_min  .Write()
              hist_max  .Write()
              graph_cent.Write()
              graph_min .Write()
              graph_max .Write()
              graph_frac.Write()
           
            #---------------------------
            # Create some pdf with plots
            outpdfname = os.path.join(self.output_dir, 'control_plots', 'RM_and_MPI' )
            if not os.path.exists(outpdfname):
              os.makedirs(outpdfname)
            outpdfname = os.path.join(outpdfname, 'mpi_corr_%s.pdf'%(label) )
            self.plot_MPI_correction_histograms( h2_mpi_on, h2_mpi_off, h2_mpi_ratio, outpdfname)

  #---------------------------------------------------------------
  # Given a pair of bin-edge values, return their index
  #---------------------------------------------------------------
  def bin_position( self , list_pt_th , min_p , max_p ):
    min_b = list_pt_th.index(min_p)+1
    max_b = list_pt_th.index(max_p)
    return min_b, max_b

  #---------------------------------------------------------------
  # Put together a label commonly used by several functions in code
  #---------------------------------------------------------------
  def create_label(self,jetR,obs_setting,grooming_setting):

    label = "R%s_" % (str(jetR))
    label += self.subobs_label(obs_setting)

    if grooming_setting:
      label += '_'
      label += self.utils.grooming_label(grooming_setting)

    label += 'Scaled'

    return label

  #---------------------------------------------------------------
  # Returns number proportional to shape of inclusive jet pT
  #     distribution theory prediction (val = jetR)
  #---------------------------------------------------------------
  def pt_scale_factor_jetR(self, ptmin, ptmax, jetR):

    pt_li, val_li = getattr(self, "pt_scale_factor_R%s" % jetR)

    # Fit a log function between the two endpoints and approx avg integral for bin
    start_i = pt_li.index(ptmin)
    end_i = pt_li.index(ptmax)
    # y = a * x^k
    k = np.log(val_li[start_i] / val_li[end_i]) / np.log(ptmin / ptmax)
    a = val_li[start_i] / ptmin**k

    return self.pt_scale_factor_k(ptmin, ptmax, k, a)

  #---------------------------------------------------------------
  # Returns number proportional to the integral of power law pTjet^k
  #---------------------------------------------------------------
  def pt_scale_factor_k(self, ptmin, ptmax, k, a=1e9):
    # The scale factor is defined as the area under the curve (i.e. the integral of the pT spectrum)
    # which corresponds to int_ptmin^ptmax (a * x^k dx) = a*x^(k+1)/(k+1) |_ptmin^ptmax = a*(ptmax^(k+1)-ptmin^(k + 1))/(k+1)
    # and this quantifies the relative number of jets in that pT range
    if k == -1:
      return a * np.log(ptmax / ptmin)
    return a * (ptmax**(k + 1) - ptmin**(k + 1)) / (k + 1)

  #----------------------------------------------------------------------
  # Extrapolate y-values for values in xlist_new given points (x,y) in xlist and ylist
  # Use power=1 for linear, or power=2 for quadratic extrapolation
  #----------------------------------------------------------------------
  def list_interpolate(self,xlist, ylist, xlist_new, power=1, require_positive=False):

    if len(xlist) < (power + 1):
      raise ValueError("list_interpolate() requires at least %i points!" % (power + 1))

    ylist_new = []
    ix = 0
    for xval in xlist_new:

      while (ix + power) < len(xlist) and xlist[ix+power] <= xval:
        ix += 1

      x1 = xlist[ix]; y1 = ylist[ix]

      # Check if data point is identical
      if xval == x1:
        if require_positive and y1 < 0:
          ylist_new.append(0)
          continue
        ylist_new.append(y1)
        continue

      # Set value to 0 if out-of-range for extrapolation
      if x1 > xval or (ix + power) >= len(xlist):
        ylist_new.append(0)
        continue

      x2 = xlist[ix+1]; y2 = ylist[ix+1]

      yval = None
      if power == 1:  # linear
        yval = self.linear_extrapolate(x1, y1, x2, y2, xval)
      else:
        raise ValueError("Unrecognized power", power, "/ please use either 1 or 2")

      # Require positive values
      if require_positive and yval < 0:
        ylist_new.append(0)
        continue

      ylist_new.append(yval)

    return ylist_new

  #---------------------------------------------------------------
  # Given two data points, find linear fit and y-value for x
  #---------------------------------------------------------------
  def linear_extrapolate(self,x1, y1, x2, y2, x):
    return (y2 - y1) / (x2 - x1) * x + (y1 - (y2 - y1) / (x2 - x1) * x1)

  #---------------------------------------------------------------
  # Set LHS of distributions to 0 if crosses to 0 at some point (prevents multiple peaks)
  #---------------------------------------------------------------
  def set_zero_range(self,yvals):

    found_nonzero_val = False

    # Step through list backwards
    for i in range(len(yvals)-1, -1, -1):
      if yvals[i] <= 0:
        if found_nonzero_val:
          for j in range(0, i+1):
            yvals[j] = 0
          break
        yvals[i] = 0
        continue
      else:
        found_nonzero_val = True
        continue

    return yvals

  #---------------------------------------------------------------
  # Where there are single values pos/neg between two neg/pos,
  # interpolate point
  #---------------------------------------------------------------
  def fix_fluctuations(self,yvals):
 
    for i in range(1, len(yvals) - 1):
      if yvals[i] > 0:
        if yvals[i+1] < 0 and yvals[i-1] < 0:
          yvals[i] = (yvals[i+1] + yvals[i-1]) / 2
      else:  # yvals[i] <= 0
        if yvals[i+1] > 0 and yvals[i-1] > 0:
          yvals[i] = (yvals[i+1] + yvals[i-1]) / 2

    return yvals

  #---------------------------------------------------------------
  # Given a list of histograms, return two histograms with the min
  # and max histograms
  #---------------------------------------------------------------
  def min_max( self , histo_list ):
    nHist = len(histo_list)
    nBins = histo_list[0].GetNbinsX()
    
    hist_min = histo_list[0].Clone()
    hist_min.SetName(hist_min.GetName().replace('_sv0','')+'_min')

    hist_max = histo_list[0].Clone()
    hist_max.SetName(hist_max.GetName().replace('_sv0','')+'_max')
    
    for b in range(0,nBins):
        bin_content = []
        for h in range(0,nHist):
            bin_content.append(histo_list[h].GetBinContent(b+1))
        min_cont = min(bin_content)
        max_cont = max(bin_content)
    
        hist_min.SetBinContent(b+1,min_cont)
        hist_max.SetBinContent(b+1,max_cont)
    
    return hist_min, hist_max

  #---------------------------------------------------------------
  # Given a central, min, and max histograms, return a graph
  #---------------------------------------------------------------
  def histo_to_graph(self,hc,hmin,hmax):
    tlx = []
    tly = []
    tl_min = []
    tl_max = []
    
    nBins = hc.GetNbinsX()
    listofzeros = [0] * nBins
 
    for b in range(0,nBins):
      tlx.append(hc.GetXaxis().GetBinCenter(b+1))
      cent = hc.GetBinContent(b+1)
      tly.append(cent)
      tl_min.append(cent-hmin.GetBinContent(b+1))
      tl_max.append(hmax.GetBinContent(b+1)-cent)
    
    lx = array('d',tlx)
    ly = array('d',tly)
    l_min = array('d',tl_min)
    l_max = array('d',tl_max)
    l_0 = array('d',listofzeros) 

    graph = ROOT.TGraphAsymmErrors(nBins,lx,ly,l_0,l_0,l_min,l_max)
    
    return graph 

  #----------------------------------------------------------------------
  def subobs_label( self , subobs ):
    label = ''
    if type(subobs)==str:
      label = subobs
    elif type(subobs)==int:
      label = '%i' % subobs
    elif type(subobs)==float:
      #label = str(subobs).replace('.','')
      label = str(subobs)
    else:
      'This option has not been created yet. Bailing out!'
      exit()
    return label

  #----------------------------------------------------------------------
  # Return Grooming Parameters
  #----------------------------------------------------------------------
  def get_grooming_setting( self , subconf ):
    if 'SoftDrop' in subconf:
      zcut = subconf['SoftDrop']['zcut']
      beta = subconf['SoftDrop']['beta']
      return {'sd':[zcut,beta]}
    else:
      return None

  #----------------------------------------------------------------------
  # Plot processed functions
  #----------------------------------------------------------------------
  def plot_processed_functions( self, graphs, graph_min, graph_max, outpdfname): 
    c1 = ROOT.TCanvas('c1','c1',900,600)
    c1.SetLeftMargin(0.22)
    c1.SetRightMargin(0.03)
    c1.SetBottomMargin(0.15)
    graphs.Draw('ALE3')
    graph_min.Draw('sameL')
    graph_max.Draw('sameL')
    c1.Draw()
    c1.Print(outpdfname)
    del c1

  #----------------------------------------------------------------------
  # Plot the MPI-correction histograms
  #----------------------------------------------------------------------
  def plot_MPI_correction_histograms( self, h2MPIon, h2MPIoff, h2MPIcorr, outpdfname):
    c1 = ROOT.TCanvas('c1','c1',1000,300)
    c1.Divide(3,1)

    h2MPIon  .GetXaxis().SetRangeUser(self.final_pt_bins[0],self.final_pt_bins[-1]) # Only plot in the pT range that will be considered ultimately
    h2MPIoff .GetXaxis().SetRangeUser(self.final_pt_bins[0],self.final_pt_bins[-1])
    h2MPIcorr.GetXaxis().SetRangeUser(self.final_pt_bins[0],self.final_pt_bins[-1])

    for i in range(0,3):
      if i < 2:
        c1.cd(i+1).SetLogz()
      c1.cd(i+1).SetBottomMargin(0.19)
      c1.cd(i+1).SetLeftMargin(0.18) 
      c1.cd(i+1).SetRightMargin(0.15)

    c1.cd(1)
    self.pretty_TH2D(h2MPIon, 'MPI on')
    h2MPIon.Draw('COLZ')
    
    c1.cd(2)
    self.pretty_TH2D(h2MPIoff, 'MPI off')
    h2MPIoff.Draw('COLZ')
    
    c1.cd(3)
    self.pretty_TH2D(h2MPIcorr, '(MPI on)/(MPI off)')
    h2MPIcorr.SetMaximum(3) # Hopefully most of the correction is small
    h2MPIcorr.Draw('COLZ')
    
    c1.Draw()
    c1.Print(outpdfname)
    del c1

  #----------------------------------------------------------------------
  # Plot the RM slices
  #----------------------------------------------------------------------
  def plot_RM_slice_histograms( self, hRM , outpdfname):
    c1 = ROOT.TCanvas('c1','c1',1000,900)
    c1.Divide(2,2)

    hRM.GetAxis(0).SetRangeUser(self.final_pt_bins[0],self.final_pt_bins[-1]) # Only plot in the pT range that will be considered ultimately
    hRM.GetAxis(1).SetRangeUser(self.final_pt_bins[0],self.final_pt_bins[-1]) # Only plot in the pT range that will be considered ultimately

    h2_had_obs_pT  = hRM.Projection(3,1)
    h2_chr_obs_pT  = hRM.Projection(2,0)
    h2_obs_chr_had = hRM.Projection(3,2)
    h2_pT_chr_had  = hRM.Projection(1,0)

    self.pretty_TH2D( h2_had_obs_pT  ,h2_had_obs_pT .GetYaxis().GetTitle()+' vs. '+h2_had_obs_pT .GetXaxis().GetTitle())
    self.pretty_TH2D( h2_chr_obs_pT  ,h2_chr_obs_pT .GetYaxis().GetTitle()+' vs. '+h2_chr_obs_pT .GetXaxis().GetTitle())
    self.pretty_TH2D( h2_obs_chr_had ,h2_obs_chr_had.GetYaxis().GetTitle()+' vs. '+h2_obs_chr_had.GetXaxis().GetTitle())
    self.pretty_TH2D( h2_pT_chr_had  ,h2_pT_chr_had .GetYaxis().GetTitle()+' vs. '+h2_pT_chr_had .GetXaxis().GetTitle())

    h2_had_obs_pT .SetMinimum(1e-7)
    h2_chr_obs_pT .SetMinimum(1e-7)
    h2_obs_chr_had.SetMinimum(1e-7)
    h2_pT_chr_had .SetMinimum(1e-7)

    for i in range(0,4):
      c1.cd(i+1).SetLogz()
      if i < 2:
        c1.cd(i+1).SetTheta(50)
        c1.cd(i+1).SetPhi(220)
        c1.cd(i+1).SetLeftMargin(0.12)
      else:
        c1.cd(i+1).SetBottomMargin(0.18)
        c1.cd(i+1).SetLeftMargin(0.17)
        c1.cd(i+1).SetRightMargin(0.16)

    c1.cd(1)
    h2_had_obs_pT.GetXaxis().SetTitleOffset(1.5)
    h2_had_obs_pT.GetYaxis().SetTitleOffset(1.5)
    h2_had_obs_pT.Draw('LEGO0')

    c1.cd(2)
    h2_had_obs_pT.GetXaxis().SetTitleOffset(1.6)
    h2_had_obs_pT.GetYaxis().SetTitleOffset(1.6)
    h2_chr_obs_pT.Draw('LEGO0')

    c1.cd(3)
    h2_obs_chr_had.Draw('COLZ')

    c1.cd(4)
    h2_pT_chr_had.Draw('COLZ')

    c1.Draw()
    c1.Print(outpdfname)
    del c1

  #----------------------------------------------------------------------
  # Edit 2D histogram
  #----------------------------------------------------------------------
  def pretty_TH2D(self, h2 ,tit='', xtit='',ytit='',ztit=''):
    h2.GetXaxis().SetLabelSize(0.06)
    h2.GetXaxis().SetTitleSize(0.06)
    h2.GetXaxis().CenterTitle()
    h2.GetXaxis().SetNdivisions(107)
    if xtit!='':
        h2.GetXaxis().SetTitle(xtit)
    h2.GetXaxis().SetTitleOffset(1.4)
    
    h2.GetYaxis().SetLabelSize(0.06)
    h2.GetYaxis().SetTitleSize(0.06)
    h2.GetYaxis().CenterTitle()
    h2.GetYaxis().SetNdivisions(107)
    if ytit!='':
        h2.GetYaxis().SetTitle(ytit)
    h2.GetYaxis().SetTitleOffset(1.4)
    
    h2.GetZaxis().SetLabelSize(0.06)
    h2.GetZaxis().SetTitleSize(0.06)
    h2.GetZaxis().CenterTitle()
    h2.GetZaxis().SetNdivisions(106)
    if ztit!='':
      h2.GetZaxis().SetTitle(ztit)
    h2.GetZaxis().SetTitleOffset(1.3) 
        
    h2.SetTitle(tit)

  #----------------------------------------------------------------------
  # Edit 1D histogram and graphs
  #----------------------------------------------------------------------
  def pretty_1D_object(self,h1,color=1,lwidth=1,lstyle=1,tit='',xtit='',ytit='',Fill=False):
    h1.SetLineColor(color)
    h1.SetMarkerColor(color)
    if Fill:
        h1.SetFillColorAlpha(color,0.2)
    h1.SetLineStyle(lstyle)
    h1.SetLineWidth(lwidth)
    
    h1.GetXaxis().SetLabelSize(0.07)
    h1.GetXaxis().SetTitleSize(0.07)
    h1.GetXaxis().CenterTitle()
    h1.GetXaxis().SetNdivisions(107)
    if xtit != '':
        h1.GetXaxis().SetTitle(xtit)
    h1.GetXaxis().SetTitleOffset(1)
    
    h1.GetYaxis().SetLabelSize(0.08)
    h1.GetYaxis().SetTitleSize(0.08)
    h1.GetYaxis().CenterTitle()
    h1.GetYaxis().SetNdivisions(107)
    if ytit != '':
        h1.GetYaxis().SetTitle(ytit)
    h1.GetYaxis().SetTitleOffset(1.2)
    
    if tit != '':    
        h1.SetTitle(tit)

  #----------------------------------------------------------------------
  # Produce graph with fractional error
  #----------------------------------------------------------------------
  def fractional_error(self, h_cen, h_min, h_max):
    bin_ctr = []
    min_val = []
    max_val = []
    npts = 0
    for b in range(0,h_cen.GetNbinsX()):
      npts += 1
      bin_ctr.append(h_cen.GetBinCenter (b+1))
      cent_val = h_cen.GetBinContent(b+1)
      if cent_val==0.:
        cent_val = 1.
      min_val.append((h_cen.GetBinContent(b+1)-h_min.GetBinContent(b+1))/cent_val)
      max_val.append((h_max.GetBinContent(b+1)-h_cen.GetBinContent(b+1))/cent_val)
    listofones = [1] * npts
    listofzeros = [0] * npts

    bin_ctr = array('d',bin_ctr)
    listofones = array('d',listofones)
    listofzeros = array('d',listofzeros)
    min_val = array('d',min_val)
    max_val = array('d',max_val)

    graph = ROOT.TGraphAsymmErrors(npts,bin_ctr,listofones,listofzeros,listofzeros,min_val,max_val)
    return graph

  #----------------------------------------------------------------------
  # Plot a comparison of the input SCET calculation to generator projections
  #----------------------------------------------------------------------
  def plot_comparison_SCET_gen_input( self, g_scet, jetR , obs_setting , grooming_setting, lev0, lev1, lev2, gen, pTmin, pTmax, outpdfname):

    label = self.create_label( jetR , obs_setting , grooming_setting )
    name_RM = "hResponse_JetPt_"+self.observable+"_"+lev0+"_"+lev1+"_MPI"+lev2+"_"+label
    hRM_name = '%s_Rebinned_%s'  % (name_RM, gen)

    hRM = getattr(self,hRM_name)
    hRM.GetAxis(1).SetRangeUser(pTmin,pTmax)
    h1_gen = hRM.Projection(3)
    h1_gen.SetName(hRM.GetName()+'proj_obs_pT_%i_%i'%(pTmin,pTmax))
    for b in range(0,h1_gen.GetNbinsX()):
      h1_gen.SetBinContent(b+1,(h1_gen.GetBinContent(b+1))/(h1_gen.GetBinWidth(b+1)))
    h1_gen.Scale(1./h1_gen.Integral('width'))
    h1_gen.SetLineColor(1)
    h1_gen.SetMinimum(0)

    xtit = self.obs_label
    ytit = '#frac{1}{#sigma} #frac{d#sigma}{d'+xtit+'}'
    tit = str(obs_setting)
    if grooming_setting:
      if 'sd' in grooming_setting: 
        tit+=' SD (z_{cut}=%.1f,#beta=%i)'%(grooming_setting['sd'][0],grooming_setting['sd'][1])
    tit += ', input %i < #it{p}_{T}^{jet} < %i GeV/#it{c}'%((int)(pTmin),(int)(pTmax))
    self.pretty_1D_object(h1_gen,1,2,1,tit, xtit, ytit)

    c1 = ROOT.TCanvas('c1','c1',900,600)
    c1.SetLeftMargin(0.22)
    c1.SetRightMargin(0.03)
    c1.SetBottomMargin(0.15)
    h1_gen.Draw()
    h1_gen.Draw('samehist')
    
    g_scet.Draw('sameLE3')

    leg = ROOT.TLegend(0.7,0.6,0.85,0.85)
    leg.SetLineColor(0)
    leg.AddEntry(g_scet,'SCET')
    leg.AddEntry(h1_gen,gen+', '+lev0+', MPI '+lev2)
    leg.Draw('same')

    c1.Draw()
    c1.Print(outpdfname)
    del h1_gen
    del leg
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