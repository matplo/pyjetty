#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import ROOT

import tqdm
import yaml
import copy
import argparse
import os

from pyjetty.mputils import *

from pyjetty.alice_analysis.process.base import process_base
from pyjetty.alice_analysis.process.user.ang_pp.helpers import lambda_alpha_kappa_i

from array import array
import numpy as np

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)
# Automatically set Sumw2 when creating new histograms
ROOT.TH1.SetDefaultSumw2()

################################################################
class herwig_parton_hadron(process_base.ProcessBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, input_file='', config_file='', output_dir='',
                 debug_level=0, args=None, **kwargs):

        super(herwig_parton_hadron, self).__init__(
            input_file, config_file, output_dir, debug_level, **kwargs)

        self.initialize_config(args)


    #---------------------------------------------------------------
    # Initialize config file into class members
    #---------------------------------------------------------------
    def initialize_config(self, args):
    
        # Call base class initialization
        process_base.ProcessBase.initialize_config(self)

        # Read config file
        with open(self.config_file, 'r') as stream:
            config = yaml.safe_load(stream)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.herwig_file = args.input_file
        self.herwig_file_MPI = args.input_file_mpi

        # Defaults to None if not in use
        self.level = args.no_match_level

        self.jetR_list = config["jetR"]
        self.alpha_list = config["alphas"]

        # SoftDrop parameters
        self.use_SD = True   # Change this to use SD
        self.sd_beta = config["sd_beta"]
        self.sd_zcut = config["sd_zcut"]
        self.grooming_settings = [{'sd': [self.sd_zcut, self.sd_beta]}]  # self.utils.grooming_settings
        self.grooming_labels = [self.utils.grooming_label(gs) for gs in self.grooming_settings]

        self.n_pt_bins = config["n_pt_bins"]
        self.pt_limits = config["pt_limits"]
        self.n_lambda_bins = config['n_lambda_bins']
        self.lambda_limits = config['lambda_limits']

        # Manually added binnings for RM and scaling histograms
        self.pt_bins = array('d', list(range(10, 50, 5)) + list(range(50, 210, 10)))
        self.obs_bins = np.concatenate((np.linspace(0, 0.0009, 10), np.linspace(0.001, 0.009, 9),
                                        np.linspace(0.01, 0.1, 19), np.linspace(0.11, 1., 90)))

        # hadron level - ALICE tracking restriction
        self.max_eta_hadron = 0.9

        # Whether or not to rescale final jet histograms based on sigma/N
        self.no_scale = args.no_scale

        # Whether or not to save particle info in raw tree structure
        self.no_tree = args.no_tree

        # Initialize variables for final cross sections from event generator
        self.xsec = None
        self.xsec_MPI = None


    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def herwig_parton_hadron(self, args):

        # Create ROOT TTree file for storing raw PYTHIA particle information
        outf_path = os.path.join(self.output_dir, args.tree_output_fname)
        outf = ROOT.TFile(outf_path, 'recreate')
        outf.cd()

        # Initialize response histograms
        self.initialize_hist()

        # Print the banner first
        fj.ClusterSequence.print_banner()
        print()

        self.init_jet_tools()
        self.parse_events()
        if self.herwig_file_MPI: 
            self.parse_events(MPIon=True)

        if not self.no_tree:
            for jetR in self.jetR_list:
                getattr(self, "tw_R%s" % str(jetR).replace('.', '')).fill_tree()

        # Scale histograms
        self.scale_print_final_info()

        outf.Write()
        outf.Close()

        self.save_output_objects()


    #---------------------------------------------------------------
    # Initialize histograms
    #---------------------------------------------------------------
    def initialize_hist(self):

        self.hNevents = ROOT.TH1I("hNevents", 'Number accepted events (unscaled)', 2, -0.5, 1.5)
        self.hNeventsMPI = ROOT.TH1I("hNeventsMPI", 'Number accepted events (unscaled)', 2, -0.5, 1.5)

        for jetR in self.jetR_list:

            # Store a list of all the histograms just so that we can rescale them later
            hist_list_name = "hist_list_R%s" % str(jetR).replace('.', '')
            setattr(self, hist_list_name, [])
            hist_list_name_MPIon = "hist_list_MPIon_R%s" % str(jetR).replace('.', '')
            setattr(self, hist_list_name_MPIon, [])

            R_label = str(jetR).replace('.', '') + 'Scaled'

            if self.level in [None, 'ch']:
                name = 'hJetPt_ch_R%s' % R_label
                h = ROOT.TH1F(name, name+';p_{T}^{ch jet};#frac{dN}{dp_{T}^{ch jet}};', 300, 0, 300)
                h.Sumw2()  # enables calculation of errors
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hNconstit_Pt_ch_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{ch jet}')
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

            if self.level in [None, 'h']:
                name = 'hJetPt_h_R%s' % R_label
                h = ROOT.TH1F(name, name+';p_{T}^{jet, h};#frac{dN}{dp_{T}^{jet, h}};', 300, 0, 300)
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hNconstit_Pt_h_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{h jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{h jet}')
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

            if self.level in [None, 'p']:
                name = 'hJetPt_p_R%s' % R_label
                h = ROOT.TH1F(name, name+';p_{T}^{jet, parton};#frac{dN}{dp_{T}^{jet, parton}};',
                              300, 0, 300)
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hNconstit_Pt_p_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{p jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{p jet}')
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

            if self.level == None:
                name = 'hJetPtRes_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
                h.GetXaxis().SetTitle('#it{p}_{T}^{parton jet}')
                h.GetYaxis().SetTitle(
                    '#frac{#it{p}_{T}^{parton jet}-#it{p}_{T}^{ch jet}}{#it{p}_{T}^{parton jet}}')
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hResponse_JetPt_R%s' % R_label
                h = ROOT.TH2F(name, name, 200, 0, 200, 200, 0, 200)
                h.GetXaxis().SetTitle('#it{p}_{T}^{parton jet}')
                h.GetYaxis().SetTitle('#it{p}_{T}^{ch jet}')
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                '''
                # Jet multiplicity for matched jets with a cut at ch-jet level
                name = 'hNconstit_Pt_ch_PtBinCH60-80_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{ch jet}')
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hNconstit_Pt_h_PtBinCH60-80_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{h jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{h jet}')
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hNconstit_Pt_p_PtBinCH60-80_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{parton jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{parton jet}')
                h.Sumw2()
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)
                '''

            for alpha in self.alpha_list:

                label = ("R%s_%s" % (str(jetR), str(alpha))).replace('.', '')

                if self.level in [None, 'ch']:
                    name = 'hAng_JetPt_ch_%sScaled' % label
                    h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                  len(self.obs_bins)-1, self.obs_bins)
                    h.GetXaxis().SetTitle('p_{T}^{ch jet}')
                    h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{#alpha=%s}^{ch}}' % str(alpha))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hAng_JetPt_ch_MPIon_%sScaled' % label
                    h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                  len(self.obs_bins)-1, self.obs_bins)
                    h.GetXaxis().SetTitle('p_{T}^{ch jet}')
                    h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{#alpha=%s}^{ch}}' % str(alpha))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name_MPIon).append(h)

                    if self.use_SD:
                        # SoftDrop groomed jet histograms for MPI scaling
                        for gl in self.grooming_labels:
                            name = 'hAng_JetPt_ch_%s_%sScaled' % (label, gl)
                            h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                          len(self.obs_bins)-1, self.obs_bins)
                            h.GetXaxis().SetTitle('p_{T}^{ch jet}')
                            h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{#alpha=%s}^{ch}}' % str(alpha))
                            h.Sumw2()
                            setattr(self, name, h)
                            getattr(self, hist_list_name).append(h)

                            name = 'hAng_JetPt_ch_MPIon_%s_%sScaled' % (label, gl)
                            h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                          len(self.obs_bins)-1, self.obs_bins)
                            h.GetXaxis().SetTitle('p_{T}^{ch jet}')
                            h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{#alpha=%s}^{ch}}' % str(alpha))
                            h.Sumw2()
                            setattr(self, name, h)
                            getattr(self, hist_list_name_MPIon).append(h)

                if self.level in [None, 'h']:
                    name = 'hAng_JetPt_h_%sScaled' % label
                    h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                  len(self.obs_bins)-1, self.obs_bins)
                    h.GetXaxis().SetTitle('p_{T}^{jet, h}')
                    h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{#alpha=%s}^{h}}' % str(alpha))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    if self.use_SD:
                        for gl in self.grooming_labels:
                            name = 'hAng_JetPt_h_%s_%sScaled' % (label, gl)
                            h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                          len(self.obs_bins)-1, self.obs_bins)
                            h.GetXaxis().SetTitle('p_{T}^{jet, h}')
                            h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{#alpha=%s}^{h}}' % str(alpha))
                            h.Sumw2()
                            setattr(self, name, h)
                            getattr(self, hist_list_name).append(h)
                            

                if self.level in [None, 'p']:
                    name = 'hAng_JetPt_p_%sScaled' % label
                    h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                  len(self.obs_bins)-1, self.obs_bins)
                    h.GetXaxis().SetTitle('p_{T}^{jet, parton}')
                    h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{#alpha=%s}^{parton}}' % str(alpha))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    if self.use_SD:
                        for gl in self.grooming_labels:
                            name = 'hAng_JetPt_p_%s_%sScaled' % (label, gl)
                            h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                          len(self.obs_bins)-1, self.obs_bins)
                            h.GetXaxis().SetTitle('p_{T}^{jet, parton}')
                            h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{#alpha=%s}^{parton}}' % str(alpha))
                            h.Sumw2()
                            setattr(self, name, h)
                            getattr(self, hist_list_name).append(h)

                if self.level == None:
                    name = 'hResponse_ang_%sScaled' % label
                    h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
                    h.GetXaxis().SetTitle('#lambda_{#alpha=%s}^{parton}' % alpha)
                    h.GetYaxis().SetTitle('#lambda_{#alpha=%s}^{ch}' % alpha)
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    if self.use_SD:
                        for gl in self.grooming_labels:
                            name = 'hResponse_ang_%s_%sScaled' % (label, gl)
                            h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
                            h.GetXaxis().SetTitle('#lambda_{#alpha=%s}^{parton}' % alpha)
                            h.GetYaxis().SetTitle('#lambda_{#alpha=%s}^{ch}' % alpha)
                            h.Sumw2()
                            setattr(self, name, h)
                            getattr(self, hist_list_name).append(h)

                    '''
                    name = 'hResponse_ang_PtBinCH20-40_%sScaled' % label
                    h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
                    h.GetXaxis().SetTitle('#lambda_{#alpha=%s}^{parton}' % alpha)
                    h.GetYaxis().SetTitle('#lambda_{#alpha=%s}^{ch}' % alpha)
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hResponse_ang_PtBinCH40-60_%sScaled' % label
                    h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
                    h.GetXaxis().SetTitle('#lambda_{#alpha=%s}^{parton}' % alpha)
                    h.GetYaxis().SetTitle('#lambda_{#alpha=%s}^{ch}' % alpha)
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hResponse_ang_PtBinCH60-80_%sScaled' % label
                    h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
                    h.GetXaxis().SetTitle('#lambda_{#alpha=%s}^{parton}' % alpha)
                    h.GetYaxis().SetTitle('#lambda_{#alpha=%s}^{ch}' % alpha)
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    # Phase space plots integrated over all pT bins
                    name = 'hPhaseSpace_DeltaR_Pt_ch_%sScaled' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  150, 0, 1.5)
                    h.GetXaxis().SetTitle('(p_{T, i})_{ch jet}')
                    h.GetYaxis().SetTitle('(#Delta R_{i})_{ch jet} / R')
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_DeltaR_ch_%sScaled' % label
                    h = ROOT.TH2F(name, name, 150, 0, 1.5,
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(#Delta R_{i})_{ch jet} / R')
                    h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{ch jet}' % str(alpha))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_Pt_ch_%sScaled' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(p_{T, i})_{ch jet}')
                    h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{ch jet}' % str(alpha))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_DeltaR_Pt_p_%sScaled' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  150, 0, 1.5)
                    h.GetXaxis().SetTitle('(p_{T, i})_{parton jet}')
                    h.GetYaxis().SetTitle('(#Delta R_{i})_{parton jet} / R')
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_DeltaR_p_%sScaled' % label
                    h = ROOT.TH2F(name, name, 150, 0, 1.5,
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(#Delta R_{i})_{parton jet} / R')
                    h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{parton jet}' % str(alpha))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_Pt_p_%sScaled' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(p_{T, i})_{parton jet}')
                    h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{parton jet}' % str(alpha))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    # Phase space plots binned in ch jet pT
                    name = 'hPhaseSpace_DeltaR_Pt_ch_PtBinCH60-80_%sScaled' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  150, 0, 1.5)
                    h.GetXaxis().SetTitle('(p_{T, i})_{ch jet}')
                    h.GetYaxis().SetTitle('(#Delta R_{i})_{ch jet} / R')
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_DeltaR_Pt_p_PtBinCH60-80_%sScaled' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  150, 0, 1.5)
                    h.GetXaxis().SetTitle('(p_{T, i})_{parton jet}')
                    h.GetYaxis().SetTitle('(#Delta R_{i})_{parton jet} / R')
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_DeltaR_ch_PtBinCH60-80_%sScaled' % label
                    h = ROOT.TH2F(name, name, 150, 0, 1.5,
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(#Delta R_{i})_{ch jet} / R')
                    h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{ch jet}' % str(alpha))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_DeltaR_p_PtBinCH60-80_%sScaled' % label
                    h = ROOT.TH2F(name, name, 150, 0, 1.5,
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(#Delta R_{i})_{parton jet} / R')
                    h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{parton jet}' % str(alpha))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_Pt_ch_PtBinCH60-80_%sScaled' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(p_{T, i})_{ch jet}')
                    h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{ch jet}' % str(alpha))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_Pt_p_PtBinCH60-80_%sScaled' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(p_{T, i})_{parton jet}')
                    h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{parton jet}' % str(alpha))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    # Annulus plots for amount of lambda contained within some r < R
                    self.annulus_plots_num_r = 150
                    self.annulus_plots_max_x = 1.5
                    low_bound = self.annulus_plots_max_x / self.annulus_plots_num_r / 2.
                    up_bound = self.annulus_plots_max_x + low_bound

                    name = 'hAnnulus_ang_ch_%sScaled' % label
                    h = ROOT.TH2F(name, name, self.annulus_plots_num_r, low_bound, up_bound,
                                  100, 0, 1.)
                    h.GetXaxis().SetTitle('(#it{r} / #it{R})_{ch jet}')
                    h.GetYaxis().SetTitle(
                        ('(#frac{#lambda_{#alpha=%s}(#it{r})}' + \
                         '{#lambda_{#alpha=%s}(#it{R})})_{ch jet}') % (str(alpha), str(alpha)))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hAnnulus_ang_ch_PtBinCH60-80_%sScaled' % label
                    h = ROOT.TH2F(name, name, self.annulus_plots_num_r, low_bound, up_bound,
                                  100, 0, 1.)
                    h.GetXaxis().SetTitle('(#it{r} / #it{R})_{ch jet}')
                    h.GetYaxis().SetTitle(
                        ('(#frac{#lambda_{#alpha=%s}(#it{r})}' + \
                         '{#lambda_{#alpha=%s}(#it{R})})_{ch jet}') % (str(alpha), str(alpha)))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hAnnulus_ang_p_%sScaled' % label
                    h = ROOT.TH2F(name, name, self.annulus_plots_num_r, low_bound, up_bound,
                                  100, 0, 1.)
                    h.GetXaxis().SetTitle('(#it{r} / #it{R})_{parton jet}')
                    h.GetYaxis().SetTitle(
                        ('(#frac{#lambda_{#alpha=%s}(#it{r})}' + \
                         '{#lambda_{#alpha=%s}(#it{R})})_{parton jet}') % (str(alpha), str(alpha)))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hAnnulus_ang_p_PtBinCH60-80_%sScaled' % label
                    h = ROOT.TH2F(name, name, self.annulus_plots_num_r, low_bound, up_bound,
                                  100, 0, 1.)
                    h.GetXaxis().SetTitle('(#it{r} / #it{R})_{parton jet}')
                    h.GetYaxis().SetTitle(
                        ('(#frac{#lambda_{#alpha=%s}(#it{r})}' + \
                         '{#lambda_{#alpha=%s}(#it{R})})_{parton jet}') % (str(alpha), str(alpha)))
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)
                    '''

                    name = "hAngResidual_JetPt_%sScaled" % label
                    h = ROOT.TH2F(name, name, 300, 0, 300, 200, -3., 1.)
                    h.GetXaxis().SetTitle('p_{T}^{jet, parton}')
                    h.GetYaxis().SetTitle('#frac{#lambda_{#alpha}^{jet, parton}-#lambda_{#alpha}' + \
                                          '^{ch jet}}{#lambda_{#alpha}^{jet, parton}}')
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = "hAngDiff_JetPt_%sScaled" % label
                    h = ROOT.TH2F(name, name, 300, 0, 300, 200, -2., 2.)
                    h.GetXaxis().SetTitle('#it{p}_{T}^{jet, ch}')
                    h.GetYaxis().SetTitle('#it{#lambda}_{#it{#alpha}}^{jet, parton}-' + \
                                          '#it{#lambda}_{#it{#alpha}}^{jet, ch}')
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    # Create THn of response
                    dim = 4
                    title = ['p_{T}^{ch jet}', 'p_{T}^{parton jet}', 
                             '#lambda_{#alpha}^{ch}', '#lambda_{#alpha}^{parton}']
                    nbins  = [len(self.pt_bins)-1,  len(self.pt_bins)-1,
                              len(self.obs_bins)-1, len(self.obs_bins)-1]
                    min_li = [self.pt_bins[0],      self.pt_bins[0],
                              self.obs_bins[0],     self.obs_bins[0]    ]
                    max_li = [self.pt_bins[-1],     self.pt_bins[-1],
                              self.obs_bins[-1],    self.obs_bins[-1]   ]

                    name = 'hResponse_JetPt_ang_ch_%sScaled' % label
                    nbins = (nbins)
                    xmin = (min_li)
                    xmax = (max_li)
                    nbins_array = array('i', nbins)
                    xmin_array = array('d', xmin)
                    xmax_array = array('d', xmax)
                    h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
                    for i in range(0, dim):
                        h.GetAxis(i).SetTitle(title[i])
                        if i == 0 or i == 1:
                            h.SetBinEdges(i, self.pt_bins)
                        else:  # i == 2 or i == 3
                            h.SetBinEdges(i, self.obs_bins)
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    if self.use_SD:
                        # SoftDrop groomed jet response matrices
                        for gl in self.grooming_labels:
                            name = 'hResponse_JetPt_ang_ch_%s_%sScaled' % (label, gl)
                            h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
                            for i in range(0, dim):
                                h.GetAxis(i).SetTitle(title[i])
                                if i == 0 or i == 1:
                                    h.SetBinEdges(i, self.pt_bins)
                                else:  # i == 2 or i == 3
                                    h.SetBinEdges(i, self.obs_bins)
                            h.Sumw2()
                            setattr(self, name, h)
                            getattr(self, hist_list_name).append(h)

                    # Another set of THn for full hadron folding
                    title = ['p_{T}^{h jet}', 'p_{T}^{parton jet}', 
                             '#lambda_{#alpha}^{h}', '#lambda_{#alpha}^{parton}']

                    name = 'hResponse_JetPt_ang_h_%sScaled' % label
                    h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
                    for i in range(0, dim):
                        h.GetAxis(i).SetTitle(title[i])
                        if i == 0 or i == 1:
                            h.SetBinEdges(i, self.pt_bins)
                        else:  # i == 2 or i == 3
                            h.SetBinEdges(i, self.obs_bins)
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    if self.use_SD:
                        # SoftDrop groomed jet response matrices
                        for gl in self.grooming_labels:
                            name = 'hResponse_JetPt_ang_h_%s_%sScaled' % (label, gl)
                            h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
                            for i in range(0, dim):
                                h.GetAxis(i).SetTitle(title[i])
                                if i == 0 or i == 1:
                                    h.SetBinEdges(i, self.pt_bins)
                                else:  # i == 2 or i == 3
                                    h.SetBinEdges(i, self.obs_bins)
                            h.Sumw2()
                            setattr(self, name, h)
                            getattr(self, hist_list_name).append(h)

                    # Finally, a set of THn for folding H --> CH (with MPI on)
                    title = ['p_{T}^{ch jet}', 'p_{T}^{h jet}', 
                             '#lambda_{#alpha}^{ch}', '#lambda_{#alpha}^{h}']

                    name = 'hResponse_JetPt_ang_Fnp_%sScaled' % label
                    h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
                    for i in range(0, dim):
                        h.GetAxis(i).SetTitle(title[i])
                        if i == 0 or i == 1:
                            h.SetBinEdges(i, self.pt_bins)
                        else:  # i == 2 or i == 3
                            h.SetBinEdges(i, self.obs_bins)
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, hist_list_name_MPIon).append(h)

                    if self.use_SD:
                        # SoftDrop groomed jet response matrices
                        for gl in self.grooming_labels:
                            name = 'hResponse_JetPt_ang_Fnp_%s_%sScaled' % (label, gl)
                            h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
                            for i in range(0, dim):
                                h.GetAxis(i).SetTitle(title[i])
                                if i == 0 or i == 1:
                                    h.SetBinEdges(i, self.pt_bins)
                                else:  # i == 2 or i == 3
                                    h.SetBinEdges(i, self.obs_bins)
                            h.Sumw2()
                            setattr(self, name, h)
                            getattr(self, hist_list_name_MPIon).append(h)


    #---------------------------------------------------------------
    # Initiate jet defs, selectors, and sd (if required)
    #---------------------------------------------------------------
    def init_jet_tools(self):
        
        for jetR in self.jetR_list:
            jetR_str = str(jetR).replace('.', '')
            
            if not self.no_tree:
                # Initialize tree writer
                name = 'particle_unscaled_R%s' % jetR_str
                t = ROOT.TTree(name, name)
                setattr(self, "t_R%s" % jetR_str, t)
                tw = RTreeWriter(tree=t)
                setattr(self, "tw_R%s" % jetR_str, tw)
            
            # set up our jet definition and a jet selector
            jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
            setattr(self, "jet_def_R%s" % jetR_str, jet_def)
            print(jet_def)

        pwarning('max eta for particles after hadronization set to', self.max_eta_hadron)
        parts_selector_h = fj.SelectorAbsEtaMax(self.max_eta_hadron)

        for jetR in self.jetR_list:
            jetR_str = str(jetR).replace('.', '')
            
            jet_selector = fj.SelectorPtMin(5.0) & \
                           fj.SelectorAbsEtaMax(self.max_eta_hadron - jetR)
            setattr(self, "jet_selector_R%s" % jetR_str, jet_selector)

            #max_eta_parton = self.max_eta_hadron + 2. * jetR
            #setattr(self, "max_eta_parton_R%s" % jetR_str, max_eta_parton)
            #pwarning("Max eta for partons with jet R =", jetR, "set to", max_eta_parton)
            #parts_selector_p = fj.SelectorAbsEtaMax(max_eta_parton)
            #setattr(self, "parts_selector_p_R%s" % jetR_str, parts_selector_p)

            count1 = 0  # Number of jets rejected from ch-h matching
            setattr(self, "count1_R%s" % jetR_str, count1)
            count2 = 0  # Number of jets rejected from h-p matching
            setattr(self, "count2_R%s" % jetR_str, count2)

    #---------------------------------------------------------------
    # Read events from output, find jets, and fill histograms
    #---------------------------------------------------------------
    def parse_events(self, MPIon=False):

        if MPIon:
            hNevents = self.hNeventsMPI
            infile = self.herwig_file_MPI
        else:
            hNevents = self.hNevents
            infile = self.herwig_file

        print("Reading events from %s..." % infile)

        with open(infile, 'r') as f:
            ev_num = 0

            # Flags to assist with keeping track of place within file
            reading_ev = False
            parton = False
            parton_final = False
            parton_finished = False
            hadron = False
            hadron_final = False

            partons_px = []
            partons_py = []
            partons_pz = []
            partons_e = []
            #partons_q = []
            hadrons_px = []
            hadrons_py = []
            hadrons_pz = []
            hadrons_e = []
            #hadrons_q = []
            ch_hadrons_px = []
            ch_hadrons_py = []
            ch_hadrons_pz = []
            ch_hadrons_e = []
            #ch_hadrons_q = []

            for line in f:

                # Waiting to start reading event
                if not reading_ev:
                    if "Event number" in line:
                        reading_ev = True
                        ev_num = int(line.split()[2])
                        if not ev_num % 1000:
                            print("Event number", ev_num, end="\r")
                        hNevents.Fill(0)
                    elif "Total integrated xsec:" in line:
                        if MPIon:
                            self.xsec_MPI = float(line.split()[3])
                        else:
                            self.xsec = float(line.split()[3])
                    continue
                
                # Reading event
                # First step is to read the parton info
                elif not parton:
                    if "ShowerHandler" in line:
                        parton = True
                    continue
                elif not parton_final:
                    if "final" in line:
                        parton_final = True
                    continue
                
                # Get showered partons
                elif not MPIon and not parton_finished:
                    # Read parton information
                    vals = line.split()
                    if line[0] == '-':
                        parton_finished = True
                    elif len(vals) == 5 and line[2] == ' ':
                        partons_px.append(vals[0])
                        partons_py.append(vals[1])
                        partons_pz.append(vals[2])
                        partons_e.append(vals[3])
                        #partons_q.append(vals[4])
                    continue

                # Get final hadrons
                elif not hadron:
                    if "DecayHandler" in line:
                        hadron = True
                    continue
                elif not hadron_final:
                    if "final" in line:
                        hadron_final = True
                    continue

                # Check if event is over
                elif line[0] == '-':
                    # Finished reading hadron info
                    reading_ev = False
                    parton = False
                    parton_final = False
                    parton_finished = False
                    hadron = False
                    hadron_final = False
                    
                    # Get correct structure for finding jets
                    partons = None
                    hadrons = None
                    if not MPIon:
                        partons = fjext.vectorize_px_py_pz_e(
                            partons_px, partons_py, partons_pz, partons_e)

                        hadrons = fjext.vectorize_px_py_pz_e(
                            hadrons_px, hadrons_py, hadrons_pz, hadrons_e)

                    ch_hadrons = fjext.vectorize_px_py_pz_e(
                        ch_hadrons_px, ch_hadrons_py, ch_hadrons_pz, ch_hadrons_e)

                    self.find_jets_fill_hist(partons, hadrons, ch_hadrons, ev_num, MPIon)

                    partons_px = []
                    partons_py = []
                    partons_pz = []
                    partons_e = []
                    #partons_q = []
                    hadrons_px = []
                    hadrons_py = []
                    hadrons_pz = []
                    hadrons_e = []
                    #hadrons_q = []
                    ch_hadrons_px = []
                    ch_hadrons_py = []
                    ch_hadrons_pz = []
                    ch_hadrons_e = []
                    #ch_hadrons_q = []

                    continue

                elif line[2].isnumeric():
                    # Save the hadron name for charge descrimination
                    i = 1
                    while line[i].isnumeric():
                        i += 1
                    hadron_type = line[i:].split()[0]
                    continue

                elif line[2] == ' ':  # and len(line.split()) == 5:
                    # Reading hadron information
                    vals = line.split()
                    if not MPIon:
                        hadrons_px.append(vals[0])
                        hadrons_py.append(vals[1])
                        hadrons_pz.append(vals[2])
                        hadrons_e.append(vals[3])
                        #hadrons_q.append(vals[4])

                    if '+' in hadron_type or '-' in hadron_type:
                        ch_hadrons_px.append(vals[0])
                        ch_hadrons_py.append(vals[1])
                        ch_hadrons_pz.append(vals[2])
                        ch_hadrons_e.append(vals[3])
                        #ch_hadrons_q.append(vals[4])

        return partons, hadrons, ch_hadrons

    #---------------------------------------------------------------
    # Read events from output, find jets, and fill histograms
    #---------------------------------------------------------------
    def find_jets_fill_hist(self, partons, hadrons,
                            ch_hadrons, iev, MPIon=False):

        for jetR in self.jetR_list:
            #print("Filling jet histograms for R = %s..." % str(jetR))
            
            jetR_str = str(jetR).replace('.', '')
            jet_selector = getattr(self, "jet_selector_R%s" % jetR_str)
            jet_def = getattr(self, "jet_def_R%s" % jetR_str)
            t = None; tw = None;
            if not self.no_tree:
                t = getattr(self, "t_R%s" % jetR_str)
                tw = getattr(self, "tw_R%s" % jetR_str)
            count1 = getattr(self, "count1_R%s" % jetR_str)
            count2 = getattr(self, "count2_R%s" % jetR_str)

            #if not (iev+1) % 1000:
            #    print("Event number %s" % str(iev+1), end='\r')

            jets_ch = fj.sorted_by_pt(jet_selector(jet_def(ch_hadrons)))
            jets_p = fj.sorted_by_pt(jet_selector(jet_def(partons)))
            jets_h = fj.sorted_by_pt(jet_selector(jet_def(hadrons)))

            if MPIon:
                for jet in jets_ch:
                    self.fill_MPI_histograms(jetR, jet)

            if self.level and not MPIon:  # Only save info at one level w/o matching
                if not self.no_tree:
                    jets = locals()["jets_%s" % self.level]
                    for jet in jets:
                        self.fill_unmatched_jet_tree(tw, jetR, iev, jet)
                continue

            for i,jchh in enumerate(jets_ch):

                # match hadron (full) jet
                drhh_list = []
                for j, jh in enumerate(jets_h):
                    drhh = jchh.delta_R(jh)
                    if drhh < jetR / 2.:
                        drhh_list.append((j,jh))
                if len(drhh_list) != 1:
                    count1 += 1
                else:  # Require unique match
                    j, jh = drhh_list[0]

                    # match parton level jet
                    dr_list = []
                    for k, jp in enumerate(jets_p):
                        dr = jh.delta_R(jp)
                        if dr < jetR / 2.:
                            dr_list.append((k, jp))
                    if len(dr_list) != 1:
                        count2 += 1
                    else:
                        k, jp = dr_list[0]

                        if self.debug_level > 0:
                            pwarning('event', iev)
                            pinfo('matched jets: ch.h:', jchh.pt(), 'h:', jh.pt(),
                                  'p:', jp.pt(), 'dr:', dr)

                        if not MPIon:
                            self.fill_jet_histograms(jetR, jp, jh, jchh)
                            if not self.no_tree:
                                self.fill_matched_jet_tree(tw, jetR, iev, jp, jh, jchh)
                        else:
                            self.fill_jet_histograms_MPI(jetR, jp, jh, jchh)

                #print("  |-> SD jet params z={0:10.3f} dR={1:10.3f} mu={2:10.3f}".format(
                #    sd_info.z, sd_info.dR, sd_info.mu))

            if MPIon:
                setattr(self, "count1_R%s_MPIon" % jetR_str, count1)
                setattr(self, "count2_R%s_MPIon" % jetR_str, count2)
            else:
                setattr(self, "count1_R%s" % jetR_str, count1)
                setattr(self, "count2_R%s" % jetR_str, count2)


    #---------------------------------------------------------------
    # Fill jet tree with (unscaled/raw) matched parton/hadron tracks
    #---------------------------------------------------------------
    def fill_matched_jet_tree(self, tw, jetR, iev, jp, jh, jchh):

        tw.fill_branch('iev', iev)
        tw.fill_branch('ch', jchh)
        tw.fill_branch('h', jh)
        tw.fill_branch('p', jp)

        kappa = 1
        for alpha in self.alpha_list:
            label = str(alpha).replace('.', '')
            tw.fill_branch("l_ch_%s" % label,
                           fjext.lambda_beta_kappa(jchh, alpha, kappa, jetR))
            tw.fill_branch("l_h_%s" % label,
                           fjext.lambda_beta_kappa(jh, alpha, kappa, jetR))
            tw.fill_branch("l_p_%s" % label,
                           fjext.lambda_beta_kappa(jp, alpha, kappa, jetR))

            # Save SoftDrop variables as well if desired
            if self.use_SD:
                for i, gs in enumerate(self.grooming_settings):
                    gl = self.grooming_labels[i]

                    # SoftDrop jets
                    gshop_chh = fjcontrib.GroomerShop(jchh, jetR, self.reclustering_algorithm)
                    jet_sd_chh = self.utils.groom(gshop_chh, gs, jetR).pair()
                    gshop_h = fjcontrib.GroomerShop(jh, jetR, self.reclustering_algorithm)
                    jet_sd_h = self.utils.groom(gshop_h, gs, jetR).pair()
                    gshop_p = fjcontrib.GroomerShop(jp, jetR, self.reclustering_algorithm)
                    jet_sd_p = self.utils.groom(gshop_p, gs, jetR).pair()

                    tw.fill_branch("l_ch_%s_%s" % (label, gl), fjext.lambda_beta_kappa(
                        jchh, jet_sd_chh, alpha, kappa, jetR))
                    tw.fill_branch("l_h_%s_%s" % (label, gl), fjext.lambda_beta_kappa(
                        jh, jet_sd_h, alpha, kappa, jetR))
                    tw.fill_branch("l_p_%s_%s" % (label, gl), fjext.lambda_beta_kappa(
                        jp, jet_sd_p, alpha, kappa, jetR))


    #---------------------------------------------------------------
    # Fill jet tree with (unscaled/raw) unmatched parton/hadron tracks
    #---------------------------------------------------------------
    def fill_unmatched_jet_tree(self, tw, jetR, iev, jet):

        tw.fill_branch('iev', iev)
        tw.fill_branch(self.level, jet)

        kappa = 1
        for alpha in self.alpha_list:
            label = str(alpha).replace('.', '')
            tw.fill_branch('l_%s_%s' % (self.level, label),
                           fjext.lambda_beta_kappa(jet, alpha, kappa, jetR))

            if self.use_SD:
                for i, gs in enumerate(self.grooming_settings):
                    gl = self.grooming_labels[i]

                    # SoftDrop jets
                    gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
                    jet_sd = self.utils.groom(gshop, gs, jetR).pair()

                    tw.fill_branch("l_ch_%s_%s" % (label, gl), fjext.lambda_beta_kappa(
                        jet, jet_sd, alpha, kappa, jetR))

    
    #---------------------------------------------------------------
    # Fill jet histograms for MPI-on PYTHIA run-through
    #---------------------------------------------------------------
    def fill_MPI_histograms(self, jetR, jet):

        for alpha in self.alpha_list:
            label = ("R%s_%s" % (str(jetR), str(alpha))).replace('.', '')
            h = getattr(self, 'hAng_JetPt_ch_MPIon_%sScaled' % label)

            kappa = 1
            h.Fill(jet.pt(), fjext.lambda_beta_kappa(jet, alpha, kappa, jetR))

            if self.use_SD:
                for i, gs in enumerate(self.grooming_settings):
                    gl = self.grooming_labels[i]
                    gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
                    jet_sd = self.utils.groom(gshop, gs, jetR).pair()

                    getattr(self, 'hAng_JetPt_ch_MPIon_%s_%sScaled' % (label, gl)).Fill(
                        jet.pt(), fjext.lambda_beta_kappa(jet, jet_sd, alpha, kappa, jetR))


    #---------------------------------------------------------------
    # Fill jet histograms
    #---------------------------------------------------------------
    def fill_jet_histograms(self, jetR, jp, jh, jch):

        R_label = str(jetR).replace('.', '') + 'Scaled'

        # Fill jet histograms which are not dependant on angualrity
        if self.level in [None, 'ch']:
            getattr(self, 'hJetPt_ch_R%s' % R_label).Fill(jch.pt())
            getattr(self, 'hNconstit_Pt_ch_R%s' % R_label).Fill(jch.pt(), len(jch.constituents()))
        if self.level in [None, 'h']:
            getattr(self, 'hJetPt_h_R%s' % R_label).Fill(jh.pt())
            getattr(self, 'hNconstit_Pt_h_R%s' % R_label).Fill(jh.pt(), len(jh.constituents()))
        if self.level in [None, 'p']:
            getattr(self, 'hJetPt_p_R%s' % R_label).Fill(jp.pt())
            getattr(self, 'hNconstit_Pt_p_R%s' % R_label).Fill(jp.pt(), len(jp.constituents()))

        if self.level == None:
            if jp.pt():  # prevent divide by 0
                getattr(self, 'hJetPtRes_R%s' % R_label).Fill(jp.pt(), (jp.pt() - jch.pt()) / jp.pt())
            getattr(self, 'hResponse_JetPt_R%s' % R_label).Fill(jp.pt(), jch.pt())

            '''
            if 60 <= jch.pt() < 80:
                getattr(self, 'hNconstit_Pt_ch_PtBinCH60-80_R%s' % R_label).Fill(
                    jch.pt(), len(jch.constituents()))
                getattr(self, 'hNconstit_Pt_h_PtBinCH60-80_R%s' % R_label).Fill(
                    jh.pt(), len(jh.constituents()))
                getattr(self, 'hNconstit_Pt_p_PtBinCH60-80_R%s' % R_label).Fill(
                    jp.pt(), len(jp.constituents()))
            '''

        # Fill angularity histograms and response matrices
        for alpha in self.alpha_list:
            self.fill_RMs(jetR, alpha, jp, jh, jch)


    #---------------------------------------------------------------
    # Fill jet histograms
    #---------------------------------------------------------------
    def fill_RMs(self, jetR, alpha, jp, jh, jch):

        # Calculate angularities
        kappa = 1
        lp = fjext.lambda_beta_kappa(jp, alpha, kappa, jetR)
        lh = fjext.lambda_beta_kappa(jh, alpha, kappa, jetR)
        lch = fjext.lambda_beta_kappa(jch, alpha, kappa, jetR)

        label = ("R%s_%s" % (str(jetR), str(alpha))).replace('.', '')

        if self.level in [None, 'ch']:
            getattr(self, 'hAng_JetPt_ch_%sScaled' % label).Fill(jch.pt(), lch)
            if self.use_SD:
                for i, gs in enumerate(self.grooming_settings):
                    gl = self.grooming_labels[i]
                    gshop = fjcontrib.GroomerShop(jch, jetR, self.reclustering_algorithm)
                    jch_sd = self.utils.groom(gshop, gs, jetR).pair()
                    getattr(self, 'hAng_JetPt_ch_%s_%sScaled' % (label, gl)).Fill(
                        jch.pt(), fjext.lambda_beta_kappa(jch, jch_sd, alpha, kappa, jetR))

        if self.level in [None, 'h']:
            getattr(self, 'hAng_JetPt_h_%sScaled' % label).Fill(jh.pt(), lh)
            if self.use_SD:
                for i, gs in enumerate(self.grooming_settings):
                    gl = self.grooming_labels[i]
                    gshop = fjcontrib.GroomerShop(jh, jetR, self.reclustering_algorithm)
                    jh_sd = self.utils.groom(gshop, gs, jetR).pair()
                    getattr(self, 'hAng_JetPt_h_%s_%sScaled' % (label, gl)).Fill(
                        jh.pt(), fjext.lambda_beta_kappa(jh, jh_sd, alpha, kappa, jetR))

        if self.level in [None, 'p']:
            getattr(self, 'hAng_JetPt_p_%sScaled' % label).Fill(jp.pt(), lp)
            if self.use_SD:
                for i, gs in enumerate(self.grooming_settings):
                    gl = self.grooming_labels[i]
                    gshop = fjcontrib.GroomerShop(jp, jetR, self.reclustering_algorithm)
                    jp_sd = self.utils.groom(gshop, gs, jetR).pair()
                    getattr(self, 'hAng_JetPt_p_%s_%sScaled' % (label, gl)).Fill(
                        jp.pt(), fjext.lambda_beta_kappa(jp, jp_sd, alpha, kappa, jetR))

        if self.level == None:
            getattr(self, 'hResponse_ang_%sScaled' % label).Fill(lp, lch)
            if self.use_SD:
                for i, gs in enumerate(self.grooming_settings):
                    gl = self.grooming_labels[i]
                    gshop_p = fjcontrib.GroomerShop(jp, jetR, self.reclustering_algorithm)
                    jp_sd = self.utils.groom(gshop_p, gs, jetR).pair()
                    gshop_ch = fjcontrib.GroomerShop(jp, jetR, self.reclustering_algorithm)
                    jch_sd = self.utils.groom(gshop_ch, gs, jetR).pair()
                    getattr(self, 'hResponse_ang_%s_%sScaled' % (label, gl)).Fill(
                        fjext.lambda_beta_kappa(jp, jp_sd, alpha, kappa, jetR), \
                        fjext.lambda_beta_kappa(jch, jch_sd, alpha, kappa, jetR))

            '''
            # Lambda at p-vs-ch-level for various bins in ch jet pT 
            if 20 <= jch.pt() < 40:
                getattr(self, 'hResponse_ang_PtBinCH20-40_%sScaled' % label).Fill(lp, lch)
            elif 40 <= jch.pt() < 60:
                getattr(self, 'hResponse_ang_PtBinCH40-60_%sScaled' % label).Fill(lp, lch)
            elif 60 <= jch.pt() < 80:
                getattr(self, 'hResponse_ang_PtBinCH60-80_%sScaled' % label).Fill(lp, lch)

            # Phase space plots and annulus histograms, including those binned in ch jet pT
            num_r = self.annulus_plots_num_r
            ang_per_r_ch = [0] * num_r
            for particle in jch.constituents():
                deltaR = particle.delta_R(jch)
                getattr(self, 'hPhaseSpace_DeltaR_Pt_ch_%sScaled' % label).Fill(
                    particle.pt(), deltaR / jetR)

                lambda_i = lambda_beta_kappa_i(particle, jch, jetR, alpha, 1)
                getattr(self, 'hPhaseSpace_ang_DeltaR_ch_%sScaled' % label).Fill(deltaR / jetR, lambda_i)
                getattr(self, 'hPhaseSpace_ang_Pt_ch_%sScaled' % label).Fill(particle.pt(), lambda_i)

                if 60 <= jch.pt() < 80:
                    getattr(self, 'hPhaseSpace_DeltaR_Pt_ch_PtBinCH60-80_%sScaled' % label).Fill(
                        particle.pt(), deltaR / jetR)
                    getattr(self, 'hPhaseSpace_ang_DeltaR_ch_PtBinCH60-80_%sScaled' % label).Fill(
                        deltaR / jetR, lambda_i)
                    getattr(self, 'hPhaseSpace_ang_Pt_ch_PtBinCH60-80_%sScaled' % label).Fill(
                        particle.pt(), lambda_i)

                ang_per_r_ch = [ang_per_r_ch[i] + lambda_i * 
                                (deltaR <= ((i+1) * jetR * self.annulus_plots_max_x / num_r))
                                for i in range(0, num_r, 1)]

            ang_per_r_p = [0] * num_r
            for particle in jp.constituents():
                deltaR = particle.delta_R(jp)
                getattr(self, 'hPhaseSpace_DeltaR_Pt_p_%sScaled' % label).Fill(
                    particle.pt(), deltaR / jetR)

                lambda_i = lambda_beta_kappa_i(particle, jp, jetR, alpha, 1)
                getattr(self, 'hPhaseSpace_ang_DeltaR_p_%sScaled' % label).Fill(deltaR / jetR, lambda_i)
                getattr(self, 'hPhaseSpace_ang_Pt_p_%sScaled' % label).Fill(particle.pt(), lambda_i)

                if 60 <= jch.pt() < 80:
                    getattr(self, 'hPhaseSpace_DeltaR_Pt_p_PtBinCH60-80_%sScaled' % label).Fill(
                        particle.pt(), deltaR / jetR)
                    getattr(self, 'hPhaseSpace_ang_DeltaR_p_PtBinCH60-80_%sScaled' % label).Fill(
                        deltaR / jetR, lambda_i)
                    getattr(self, 'hPhaseSpace_ang_Pt_p_PtBinCH60-80_%sScaled' % label).Fill(
                        particle.pt(), lambda_i)

                ang_per_r_p = [ang_per_r_p[i] + lambda_i *
                               (deltaR <= ((i+1) * jetR * self.annulus_plots_max_x / num_r))
                             for i in range(0, num_r, 1)]

            for i in range(0, num_r, 1):
                getattr(self, 'hAnnulus_ang_p_%sScaled' % label).Fill(
                    (i+1) * self.annulus_plots_max_x / num_r, ang_per_r_p[i] / (lp + 1e-11))
                getattr(self, 'hAnnulus_ang_ch_%sScaled' % label).Fill(
                    (i+1) * self.annulus_plots_max_x / num_r, ang_per_r_ch[i] / (lch + 1e-11))
                if 60 <= jch.pt() < 80:
                    getattr(self, 'hAnnulus_ang_p_PtBinCH60-80_%sScaled' % label).Fill(
                        (i+1) * self.annulus_plots_max_x / num_r, ang_per_r_p[i] / (lp + 1e-11))
                    getattr(self, 'hAnnulus_ang_ch_PtBinCH60-80_%sScaled' % label).Fill(
                        (i+1) * self.annulus_plots_max_x / num_r, ang_per_r_ch[i] / (lch + 1e-11))
            '''

            # Residual plots (with and without divisor in y-axis)
            getattr(self, "hAngDiff_JetPt_%sScaled" % label).Fill(jch.pt(), lp - lch)
            if lp:  # prevent divide by 0
                getattr(self, "hAngResidual_JetPt_%sScaled" % label).Fill(jp.pt(), (lp - lch) / lp)

            # 4D response matrices for "forward folding" to ch level
            x = ([jch.pt(), jp.pt(), lch, lp])
            x_array = array('d', x)
            getattr(self, 'hResponse_JetPt_ang_ch_%sScaled' % label).Fill(x_array)

            x = ([jh.pt(), jp.pt(), lh, lp])
            x_array = array('d', x)
            getattr(self, 'hResponse_JetPt_ang_h_%sScaled' % label).Fill(x_array)

            if self.use_SD:
                for i, gs in enumerate(self.grooming_settings):
                    gl = self.grooming_labels[i]

                    # SoftDrop jet angularities
                    gshop_ch = fjcontrib.GroomerShop(jch, jetR, self.reclustering_algorithm)
                    jet_sd_ch = self.utils.groom(gshop_ch, gs, jetR).pair()
                    lch_sd = fjext.lambda_beta_kappa(jch, jet_sd_ch, alpha, kappa, jetR)
                    gshop_h = fjcontrib.GroomerShop(jh, jetR, self.reclustering_algorithm)
                    jet_sd_h = self.utils.groom(gshop_h, gs, jetR).pair()
                    lh_sd = fjext.lambda_beta_kappa(jh, jet_sd_h, alpha, kappa, jetR)
                    gshop_p = fjcontrib.GroomerShop(jp, jetR, self.reclustering_algorithm)
                    jet_sd_p = self.utils.groom(gshop_p, gs, jetR).pair()
                    lp_sd = fjext.lambda_beta_kappa(jp, jet_sd_p, alpha, kappa, jetR)

                    x = ([jch.pt(), jp.pt(), lch_sd, lp_sd])
                    x_array = array('d', x)
                    getattr(self, 'hResponse_JetPt_ang_ch_%s_%sScaled' % (label, gl)).Fill(x_array)

                    x = ([jh.pt(), jp.pt(), lh, lp])
                    x_array = array('d', x)
                    getattr(self, 'hResponse_JetPt_ang_h_%s_%sScaled' % (label, gl)).Fill(x_array)


    #---------------------------------------------------------------
    # Fill jet histograms for MPI (which are just the H-->CH RMs)
    #---------------------------------------------------------------
    def fill_jet_histograms_MPI(self, jetR, jp, jh, jch):

        for alpha in self.alpha_list:

            # Calculate angularities
            kappa = 1
            lh = fjext.lambda_beta_kappa(jh, alpha, kappa, jetR)
            lch = fjext.lambda_beta_kappa(jch, alpha, kappa, jetR)

            label = ("R%s_%s" % (str(jetR), str(alpha))).replace('.', '')

            # 4D response matrices for "forward folding" from h to ch level
            x = ([jch.pt(), jh.pt(), lch, lh])
            x_array = array('d', x)
            getattr(self, 'hResponse_JetPt_ang_Fnp_%sScaled' % label).Fill(x_array)

            if self.use_SD:
                for i, gs in enumerate(self.grooming_settings):
                    gl = self.grooming_labels[i]

                    # SoftDrop jet angularities
                    gshop_ch = fjcontrib.GroomerShop(jch, jetR, self.reclustering_algorithm)
                    jet_sd_ch = self.utils.groom(gshop_ch, gs, jetR).pair()
                    lch_sd = fjext.lambda_beta_kappa(jch, jet_sd_ch, alpha, kappa, jetR)
                    gshop_h = fjcontrib.GroomerShop(jh, jetR, self.reclustering_algorithm)
                    jet_sd_h = self.utils.groom(gshop_h, gs, jetR).pair()
                    lh_sd = fjext.lambda_beta_kappa(jh, jet_sd_h, alpha, kappa, jetR)

                    x = ([jch.pt(), jh.pt(), lch_sd, lh_sd])
                    x_array = array('d', x)
                    getattr(self, 'hResponse_JetPt_ang_Fnp_%s_%sScaled' % (label, gl)).Fill(x_array)


    #---------------------------------------------------------------
    # Initiate scaling of all histograms and print final simulation info
    #---------------------------------------------------------------
    def scale_print_final_info(self):

        # Scale all jet histograms by the appropriate factor from generated cross section
        # and the number of accepted events
        if not self.no_scale:
            scale_f = self.xsec / self.hNevents.GetBinContent(1)
            print("Weight MPIoff histograms by (cross section)/(N events) =", scale_f)

            MPI_scale_f = None
            if self.herwig_file_MPI:
                MPI_scale_f = self.xsec_MPI / self.hNeventsMPI.GetBinContent(1)
                print("Weight MPIon histograms by (cross section)/(N events) =", MPI_scale_f)

            self.scale_jet_histograms(scale_f, MPI_scale_f)
        print()

        self.hNevents.SetBinError(1, 0)
        if self.herwig_file_MPI:
            self.hNeventsMPI.SetBinError(1, 0)

        for jetR in self.jetR_list:
            jetR_str = str(jetR).replace('.', '')
            count1 = getattr(self, "count1_R%s" % jetR_str)
            count2 = getattr(self, "count2_R%s" % jetR_str)
            print(("For R=%s:  %i jets cut at first match criteria; " + \
                  "%i jets cut at second match criteria.") % 
                  (str(jetR), count1, count2))
        print()


    #---------------------------------------------------------------
    # Scale all jet histograms by sigma/N
    #---------------------------------------------------------------
    def scale_jet_histograms(self, scale_f, MPI_scale_f):

        for jetR in self.jetR_list:
            hist_list_name = "hist_list_R%s" % str(jetR).replace('.', '')
            for h in getattr(self, hist_list_name):
                h.Scale(scale_f)

            hist_list_MPIon_name = "hist_list_MPIon_R%s" % str(jetR).replace('.', '')
            for h in getattr(self, hist_list_MPIon_name):
                h.Scale(MPI_scale_f)


################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Herwig7 debug level-1 output parser',
                                     prog=os.path.basename(__file__))
    parser.add_argument('-i', '--input-file', action='store', type=str, default='LHC.log',
                        help='Input .log file from Herwig7 analysis')
    parser.add_argument('-m', '--input-file-mpi', action='store', type=str, default=None,
                        help='Input .log file with MPI on from Herwig7 analysis')
    parser.add_argument('-o', '--output-dir', action='store', type=str, default='./', 
                        help='Output directory for generated ROOT file(s)')
    parser.add_argument('--tree-output-fname', default="AnalysisResults.root", type=str,
                        help="Filename for the (unscaled) generated particle ROOT TTree")
    parser.add_argument('--no-tree', default=False, action='store_true', 
                        help="Do not save tree of particle information, only create histograms")
    parser.add_argument('--no-match-level', help="Save simulation for only one level with " + \
                        "no matching. Options: 'p', 'h', 'ch'", default=None, type=str)
    parser.add_argument('--no-scale', help="Turn off rescaling all histograms by cross section / N",
                        action='store_true', default=False)
    parser.add_argument('-c', '--config-file', action='store', type=str,
                        default='config/angularity.yaml',
                        help="Path of config file for observable configurations")
    args = parser.parse_args()

    if args.no_match_level not in [None, 'p', 'h', 'ch']:
        print("ERROR: Unrecognized type %s. Please use 'p', 'h', or 'ch'" % args.type_only)
        exit(1)

    # If invalid configFile is given, exit
    if not os.path.exists(args.config_file):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    process = herwig_parton_hadron(
        config_file=args.config_file, output_dir=args.output_dir, args=args)
    process.herwig_parton_hadron(args)
