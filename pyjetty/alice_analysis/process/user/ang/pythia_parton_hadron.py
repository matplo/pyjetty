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

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

from pyjetty.alice_analysis.process.base import process_base
from pyjetty.alice_analysis.process.user.ang.helpers import lambda_alpha_kappa_i

from array import array
import numpy as np

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)
# Automatically set Sumw2 when creating new histograms
ROOT.TH1.SetDefaultSumw2()

################################################################
class pythia_parton_hadron(process_base.ProcessBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, input_file='', config_file='', output_dir='',
                 debug_level=0, args=None, **kwargs):

        super(pythia_parton_hadron, self).__init__(
            input_file, config_file, output_dir, debug_level, **kwargs)

        self.initialize_config(args)

    #---------------------------------------------------------------
    # Main processing function
    #---------------------------------------------------------------
    def pythia_parton_hadron(self, args):

        # Create ROOT TTree file for storing raw PYTHIA particle information
        outf_path = os.path.join(self.output_dir, args.tree_output_fname)
        outf = ROOT.TFile(outf_path, 'recreate')
        outf.cd()

        # Initialize response histograms
        self.initialize_hist()

        pinfo('user seed for pythia', self.user_seed)
        # mycfg = ['PhaseSpace:pThatMin = 100']
        mycfg = ['Random:setSeed=on', 'Random:seed={}'.format(self.user_seed)]
        mycfg.append('HadronLevel:all=off')

        # PYTHIA instance with MPI off
        setattr(args, "py_noMPI", True)
        pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)

        # print the banner first
        fj.ClusterSequence.print_banner()
        print()

        self.init_jet_tools()
        self.calculate_events(pythia)
        pythia.stat()
        print()

        # PYTHIA instance with MPI on
        setattr(args, "py_noMPI", False)
        pythia_MPI = pyconf.create_and_init_pythia_from_args(args, mycfg)
        self.calculate_events(pythia_MPI, MPIon=True)
        print()

        if not self.no_tree:
            for jetR in self.jetR_list:
                getattr(self, "tw_R%s" % str(jetR).replace('.', '')).fill_tree()

        self.scale_print_final_info(pythia, pythia_MPI)

        outf.Write()
        outf.Close()

        self.save_output_objects()

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

        # Defaults to None if not in use
        self.level = args.no_match_level

        self.jetR_list = config["jetR"]

        self.user_seed = args.user_seed
        self.nev = args.nev

        self.observable_list = config["process_observables"]
        self.obs_settings = {}
        self.obs_grooming_settings = {}
        self.obs_names = {}
        for obs in self.observable_list:

            obs_config_dict = config[obs]
            obs_config_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]

            obs_subconfig_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
            self.obs_settings[obs] = self.utils.obs_settings(
                obs, obs_config_dict, obs_subconfig_list)
            self.obs_grooming_settings[obs] = self.utils.grooming_settings(obs_config_dict)

            self.obs_names[obs] = obs_config_dict["common_settings"]["xtitle"]

        # Construct set of unique grooming settings
        self.grooming_settings = []
        lists_grooming = [self.obs_grooming_settings[obs] for obs in self.observable_list]
        for observable in lists_grooming:
            for setting in observable:
                if setting not in self.grooming_settings and setting != None:
                    self.grooming_settings.append(setting)

        # Manually added binnings for RM and scaling histograms
        self.pt_bins = array('d', list(range(5, 210, 5)))
        self.obs_bins_ang = np.concatenate((np.linspace(0, 0.0009, 10), np.linspace(0.001, 0.009, 9),
                                        np.linspace(0.01, 0.1, 19), np.linspace(0.11, 1., 90)))
        self.obs_bins_mass = np.concatenate(
            (np.linspace(0, 0.9, 10), np.linspace(1, 9.8, 45), np.linspace(10, 14.5, 10),
            np.linspace(15, 19, 5), np.linspace(20, 60, 9)))

        # hadron level - ALICE tracking restriction
        self.max_eta_hadron = 0.9

        # Whether or not to rescale final jet histograms based on sigma/N
        self.no_scale = args.no_scale

        # Whether or not to save particle info in raw tree structure
        self.no_tree = args.no_tree

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

            R_label = "R" + str(jetR).replace('.', '') + 'Scaled'

            for MPI in ["", "MPIon_"]:
                R_label = MPI + R_label
                list_name = hist_list_name_MPIon if MPI else hist_list_name
                if self.level in [None, 'ch']:
                    name = 'hJetPt_ch_%s' % R_label
                    h = ROOT.TH1F(name, name+';p_{T}^{ch jet};#frac{dN}{dp_{T}^{ch jet}};', 300, 0, 300)
                    h.Sumw2()  # enables calculation of errors
                    setattr(self, name, h)
                    getattr(self, list_name).append(h)

                    name = 'hNconstit_Pt_ch_%s' % R_label
                    h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                    h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
                    h.GetYaxis().SetTitle('#it{N}_{constit}^{ch jet}')
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, list_name).append(h)

                if self.level in [None, 'h']:
                    name = 'hJetPt_h_%s' % R_label
                    h = ROOT.TH1F(name, name+';p_{T}^{jet, h};#frac{dN}{dp_{T}^{jet, h}};', 300, 0, 300)
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, list_name).append(h)

                    name = 'hNconstit_Pt_h_%s' % R_label
                    h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                    h.GetXaxis().SetTitle('#it{p}_{T}^{h jet}')
                    h.GetYaxis().SetTitle('#it{N}_{constit}^{h jet}')
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, list_name).append(h)

                if self.level in [None, 'p']:
                    name = 'hJetPt_p_%s' % R_label
                    h = ROOT.TH1F(name, name+';p_{T}^{jet, parton};#frac{dN}{dp_{T}^{jet, parton}};',
                                  300, 0, 300)
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, list_name).append(h)

                    name = 'hNconstit_Pt_p_%s' % R_label
                    h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                    h.GetXaxis().SetTitle('#it{p}_{T}^{p jet}')
                    h.GetYaxis().SetTitle('#it{N}_{constit}^{p jet}')
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, list_name).append(h)

                if self.level == None:
                    name = 'hJetPtRes_%s' % R_label
                    h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
                    h.GetXaxis().SetTitle('#it{p}_{T}^{parton jet}')
                    h.GetYaxis().SetTitle(
                        '#frac{#it{p}_{T}^{parton jet}-#it{p}_{T}^{ch jet}}{#it{p}_{T}^{parton jet}}')
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, list_name).append(h)

                    name = 'hResponse_JetPt_%s' % R_label
                    h = ROOT.TH2F(name, name, 200, 0, 200, 200, 0, 200)
                    h.GetXaxis().SetTitle('#it{p}_{T}^{parton jet}')
                    h.GetYaxis().SetTitle('#it{p}_{T}^{ch jet}')
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, list_name).append(h)

                    '''
                    # Jet multiplicity for matched jets with a cut at ch-jet level
                    name = 'hNconstit_Pt_ch_PtBinCH60-80_%s' % R_label
                    h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                    h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
                    h.GetYaxis().SetTitle('#it{N}_{constit}^{ch jet}')
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, list_name).append(h)

                    name = 'hNconstit_Pt_h_PtBinCH60-80_%s' % R_label
                    h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                    h.GetXaxis().SetTitle('#it{p}_{T}^{h jet}')
                    h.GetYaxis().SetTitle('#it{N}_{constit}^{h jet}')
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, list_name).append(h)

                    name = 'hNconstit_Pt_p_PtBinCH60-80_%s' % R_label
                    h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                    h.GetXaxis().SetTitle('#it{p}_{T}^{parton jet}')
                    h.GetYaxis().SetTitle('#it{N}_{constit}^{parton jet}')
                    h.Sumw2()
                    setattr(self, name, h)
                    getattr(self, list_name).append(h)
                    '''

            for obs in self.observable_list:
                if obs != "ang" and obs != "mass":
                    raise ValueError("Observable %s not yet implemented in this script!" % obs)

                for i in range(len(self.obs_settings[obs])):

                    obs_setting = self.obs_settings[obs][i]
                    grooming_setting = self.obs_grooming_settings[obs][i]
                    obs_label = self.utils.obs_label(obs_setting, grooming_setting)
                    label = ("R%s_%s" % (jetR, obs_label)).replace('.', '')
                    obs_bins = getattr(self, "obs_bins_"+obs)

                    if self.level in [None, 'ch']:
                        name = 'h_%s_JetPt_ch_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                      len(obs_bins)-1, obs_bins)
                        h.GetXaxis().SetTitle('p_{T}^{ch jet}')
                        h.GetYaxis().SetTitle(self.obs_names[obs]+'^{ch}')
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'h_%s_JetPt_ch_MPIon_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                      len(obs_bins)-1, obs_bins)
                        h.GetXaxis().SetTitle('p_{T}^{ch jet}')
                        h.GetYaxis().SetTitle(self.obs_names[obs]+'^{ch}')
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name_MPIon).append(h)


                    if self.level in [None, 'h']:
                        name = 'h_%s_JetPt_h_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                      len(obs_bins)-1, obs_bins)
                        h.GetXaxis().SetTitle('p_{T}^{jet, h}')
                        h.GetYaxis().SetTitle(self.obs_names[obs]+'^{h}')
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                    if self.level in [None, 'p']:
                        name = 'h_%s_JetPt_p_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, len(self.pt_bins)-1, self.pt_bins,
                                      len(obs_bins)-1, obs_bins)
                        h.GetXaxis().SetTitle('p_{T}^{jet, parton}')
                        h.GetYaxis().SetTitle(self.obs_names[obs]+'^{parton}')
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                    if self.level == None:
                        name = 'hResponse_%s_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, len(obs_bins)-1, obs_bins, len(obs_bins)-1, obs_bins)
                        h.GetXaxis().SetTitle(self.obs_names[obs]+'^{parton}')
                        h.GetYaxis().SetTitle(self.obs_names[obs]+'^{ch}')
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        '''
                        name = 'hResponse_%s_PtBinCH20-40_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
                        h.GetXaxis().SetTitle(self.obs_names[obs]+'^{parton}')
                        h.GetYaxis().SetTitle(self.obs_names[obs]+'^{ch}')
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hResponse_%s_PtBinCH40-60_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
                        h.GetXaxis().SetTitle(self.obs_names[obs]+'^{parton}')
                        h.GetYaxis().SetTitle(self.obs_names[obs]+'^{ch}')
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hResponse_%s_PtBinCH60-80_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
                        h.GetXaxis().SetTitle(self.obs_names[obs]+'^{parton}')
                        h.GetYaxis().SetTitle(self.obs_names[obs]+'^{ch}')
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        # Phase space plots integrated over all pT bins
                        name = 'hPhaseSpace_DeltaR_Pt_ch_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                      150, 0, 1.5)
                        h.GetXaxis().SetTitle('(p_{T, i})_{ch jet}')
                        h.GetYaxis().SetTitle('(#Delta R_{i})_{ch jet} / R')
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hPhaseSpace_%s_DeltaR_ch_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, 150, 0, 1.5,
                                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                        h.GetXaxis().SetTitle('(#Delta R_{i})_{ch jet} / R')
                        h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{ch jet}' % str(alpha))
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hPhaseSpace_%s_Pt_ch_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                        h.GetXaxis().SetTitle('(p_{T, i})_{ch jet}')
                        h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{ch jet}' % str(alpha))
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hPhaseSpace_DeltaR_Pt_p_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                      150, 0, 1.5)
                        h.GetXaxis().SetTitle('(p_{T, i})_{parton jet}')
                        h.GetYaxis().SetTitle('(#Delta R_{i})_{parton jet} / R')
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hPhaseSpace_%s_DeltaR_p_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, 150, 0, 1.5,
                                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                        h.GetXaxis().SetTitle('(#Delta R_{i})_{parton jet} / R')
                        h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{parton jet}' % str(alpha))
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hPhaseSpace_%s_Pt_p_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                        h.GetXaxis().SetTitle('(p_{T, i})_{parton jet}')
                        h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{parton jet}' % str(alpha))
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        # Phase space plots binned in ch jet pT
                        name = 'hPhaseSpace_DeltaR_Pt_ch_PtBinCH60-80_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                      150, 0, 1.5)
                        h.GetXaxis().SetTitle('(p_{T, i})_{ch jet}')
                        h.GetYaxis().SetTitle('(#Delta R_{i})_{ch jet} / R')
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hPhaseSpace_DeltaR_Pt_p_PtBinCH60-80_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                      150, 0, 1.5)
                        h.GetXaxis().SetTitle('(p_{T, i})_{parton jet}')
                        h.GetYaxis().SetTitle('(#Delta R_{i})_{parton jet} / R')
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hPhaseSpace_%s_DeltaR_ch_PtBinCH60-80_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, 150, 0, 1.5,
                                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                        h.GetXaxis().SetTitle('(#Delta R_{i})_{ch jet} / R')
                        h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{ch jet}' % str(alpha))
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hPhaseSpace_%s_DeltaR_p_PtBinCH60-80_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, 150, 0, 1.5,
                                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                        h.GetXaxis().SetTitle('(#Delta R_{i})_{parton jet} / R')
                        h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{parton jet}' % str(alpha))
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hPhaseSpace_%s_Pt_ch_PtBinCH60-80_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                      self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                        h.GetXaxis().SetTitle('(p_{T, i})_{ch jet}')
                        h.GetYaxis().SetTitle('(#lambda_{#alpha=%s, i})_{ch jet}' % str(alpha))
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hPhaseSpace_%s_Pt_p_PtBinCH60-80_%sScaled' % (obs, label)
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

                        name = 'hAnnulus_%s_ch_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, self.annulus_plots_num_r, low_bound, up_bound,
                                      100, 0, 1.)
                        h.GetXaxis().SetTitle('(#it{r} / #it{R})_{ch jet}')
                        h.GetYaxis().SetTitle(
                            ('(#frac{#lambda_{#alpha=%s}(#it{r})}' + \
                             '{#lambda_{#alpha=%s}(#it{R})})_{ch jet}') % (str(alpha), str(alpha)))
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hAnnulus_%s_ch_PtBinCH60-80_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, self.annulus_plots_num_r, low_bound, up_bound,
                                      100, 0, 1.)
                        h.GetXaxis().SetTitle('(#it{r} / #it{R})_{ch jet}')
                        h.GetYaxis().SetTitle(
                            ('(#frac{#lambda_{#alpha=%s}(#it{r})}' + \
                             '{#lambda_{#alpha=%s}(#it{R})})_{ch jet}') % (str(alpha), str(alpha)))
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hAnnulus_%s_p_%sScaled' % (obs, label)
                        h = ROOT.TH2F(name, name, self.annulus_plots_num_r, low_bound, up_bound,
                                      100, 0, 1.)
                        h.GetXaxis().SetTitle('(#it{r} / #it{R})_{parton jet}')
                        h.GetYaxis().SetTitle(
                            ('(#frac{#lambda_{#alpha=%s}(#it{r})}' + \
                             '{#lambda_{#alpha=%s}(#it{R})})_{parton jet}') % (str(alpha), str(alpha)))
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = 'hAnnulus_%s_p_PtBinCH60-80_%sScaled' % (obs, label)
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

                        name = "h_%sResidual_JetPt_%sScaled" % (obs, label)
                        h = ROOT.TH2F(name, name, 300, 0, 300, 200, -3., 1.)
                        h.GetXaxis().SetTitle('p_{T}^{jet, parton}')
                        h.GetYaxis().SetTitle(
                            '#frac{' + self.obs_names[obs] + '^{parton}-' + \
                            self.obs_names[obs] + '^{ch}}{' + self.obs_names[obs] + '^{parton}}')
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        name = "h_%sDiff_JetPt_%sScaled" % (obs, label)
                        h = ROOT.TH2F(name, name, 300, 0, 300, 200, -2., 2.)
                        h.GetXaxis().SetTitle('#it{p}_{T}^{jet, ch}')
                        h.GetYaxis().SetTitle(self.obs_names[obs] + '^{parton}-' + \
                                              self.obs_names[obs] + '^{ch}')
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        # Create THn of response
                        dim = 4
                        title = ['p_{T}^{ch jet}', 'p_{T}^{parton jet}',
                                 self.obs_names[obs]+'^{ch}', self.obs_names[obs]+'^{parton}']
                        nbins  = [len(self.pt_bins)-1,  len(self.pt_bins)-1,
                                  len(obs_bins)-1, len(obs_bins)-1]
                        min_li = [self.pt_bins[0],      self.pt_bins[0],
                                  obs_bins[0],     obs_bins[0]    ]
                        max_li = [self.pt_bins[-1],     self.pt_bins[-1],
                                  obs_bins[-1],    obs_bins[-1]   ]

                        name = 'hResponse_JetPt_%s_ch_%sScaled' % (obs, label)
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
                                h.SetBinEdges(i, obs_bins)
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        # Another set of THn for full hadron folding
                        title = ['p_{T}^{h jet}', 'p_{T}^{parton jet}',
                                 self.obs_names[obs] + '^{h}', self.obs_names[obs] + '^{parton}']

                        name = 'hResponse_JetPt_%s_h_%sScaled' % (obs, label)
                        h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
                        for i in range(0, dim):
                            h.GetAxis(i).SetTitle(title[i])
                            if i == 0 or i == 1:
                                h.SetBinEdges(i, self.pt_bins)
                            else:  # i == 2 or i == 3
                                h.SetBinEdges(i, obs_bins)
                        h.Sumw2()
                        setattr(self, name, h)
                        getattr(self, hist_list_name).append(h)

                        # Finally, a set of THn for folding H --> CH (with MPI on)
                        title = ['p_{T}^{ch jet}', 'p_{T}^{h jet}',
                                 self.obs_names[obs] + '^{ch}', self.obs_names[obs] + '^{h}']

                        name = 'hResponse_JetPt_%s_Fnp_%sScaled' % (obs, label)
                        h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
                        for i in range(0, dim):
                            h.GetAxis(i).SetTitle(title[i])
                            if i == 0 or i == 1:
                                h.SetBinEdges(i, self.pt_bins)
                            else:  # i == 2 or i == 3
                                h.SetBinEdges(i, obs_bins)
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
    # Calculate events and pass information on to jet finding
    #---------------------------------------------------------------
    def calculate_events(self, pythia, MPIon=False):
        
        iev = 0  # Event loop count

        if MPIon:
            hNevents = self.hNeventsMPI
        else:
            hNevents = self.hNevents

        while hNevents.GetBinContent(1) < self.nev:
            if not pythia.next():
                continue

            parts_pythia_p = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)

            hstatus = pythia.forceHadronLevel()
            if not hstatus:
                #pwarning('forceHadronLevel false event', iev)
                continue
            #parts_pythia_h = pythiafjext.vectorize_select(
            #     pythia, [pythiafjext.kHadron, pythiafjext.kCharged])
            parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)

            parts_pythia_hch = pythiafjext.vectorize_select(
                pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, True)

            """ TODO: fix for multiple jet R
            parts_pythia_p_selected = parts_selector_p(parts_pythia_p)
            parts_pythia_h_selected = parts_selector_h(parts_pythia_h)
            parts_pythia_hch_selected = parts_selector_h(parts_pythia_hch)

            if self.debug_level > 1:
                pinfo('debug partons...')
                for p in parts_pythia_p_selected:
                    pyp = pythiafjext.getPythia8Particle(p)
                    print(pyp.name())
                pinfo('debug hadrons...')
                for p in parts_pythia_h_selected:
                    pyp = pythiafjext.getPythia8Particle(p)
                    print(pyp.name())
                pinfo('debug ch. hadrons...')
                for p in parts_pythia_hch_selected:
                    pyp = pythiafjext.getPythia8Particle(p)
                    print(pyp.name())
            """

            # Some "accepted" events don't survive hadronization step -- keep track here
            hNevents.Fill(0)
            self.find_jets_fill_trees(parts_pythia_p, parts_pythia_h, parts_pythia_hch, iev, MPIon)

            iev += 1


    #---------------------------------------------------------------
    # Find jets, do matching between levels, and fill histograms & trees
    #---------------------------------------------------------------
    def find_jets_fill_trees(self, parts_pythia_p, parts_pythia_h, parts_pythia_hch,
                             iev, MPIon=False):

        for jetR in self.jetR_list:
            jetR_str = str(jetR).replace('.', '')
            jet_selector = getattr(self, "jet_selector_R%s" % jetR_str)
            jet_def = getattr(self, "jet_def_R%s" % jetR_str)
            t = None; tw = None;
            if not self.no_tree:
                t = getattr(self, "t_R%s" % jetR_str)
                tw = getattr(self, "tw_R%s" % jetR_str)
            count1 = getattr(self, "count1_R%s" % jetR_str)
            count2 = getattr(self, "count2_R%s" % jetR_str)

            # parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
            jets_p = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_p)))
            jets_h = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h)))
            jets_ch = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_hch)))

            if MPIon:
                if self.level:
                    for jet in locals()["jets_" + self.level]:
                        self.fill_unmatched_histograms(jetR, jet, self.level, MPI=True)
                    continue  # Don't need to do matching
                else:
                    for jet in jets_p:
                        self.fill_unmatched_histograms(jetR, jet, 'p', MPI=True)
                    for jet in jets_h:
                        self.fill_unmatched_histograms(jetR, jet, 'h', MPI=True)
                    for jet in jets_ch:
                        self.fill_unmatched_histograms(jetR, jet, 'ch', MPI=True)

            else:  # MPI off
                if self.level:
                    for jet in locals()["jets_" + self.level]:
                        self.fill_unmatched_histograms(jetR, jet, self.level)
                        if not self.no_tree:
                            self.fill_unmatched_jet_tree(tw, jetR, iev, jet)
                    continue  # Don't need to do matching
                else:
                    for jet in jets_p:
                        self.fill_unmatched_histograms(jetR, jet, 'p')
                    for jet in jets_h:
                        self.fill_unmatched_histograms(jetR, jet, 'h')
                    for jet in jets_ch:
                        self.fill_unmatched_histograms(jetR, jet, 'ch')

            # Start the matching procedure
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
    def fill_matched_jet_tree(self, tw, jetR, iev, jp, jh, jch):

        tw.fill_branch('iev', iev)
        tw.fill_branch('ch', jch)
        tw.fill_branch('h', jh)
        tw.fill_branch('p', jp)

        self.fill_unmatched_jet_tree(tw, jetR, iev, jp, level='p', save_iev=False)
        self.fill_unmatched_jet_tree(tw, jetR, iev, jh, level='h', save_iev=False)
        self.fill_unmatched_jet_tree(tw, jetR, iev, jch, level='ch', save_iev=False)

    #---------------------------------------------------------------
    # Fill jet tree with (unscaled/raw) unmatched parton/hadron tracks
    #---------------------------------------------------------------
    def fill_unmatched_jet_tree(self, tw, jetR, iev, jet, level='ch', save_iev=True):

        if save_iev:
            tw.fill_branch('iev', iev)
            tw.fill_branch(level, jet)

        for obs in self.observable_list:
            for i in range(len(self.obs_settings[obs])):

                obs_setting = self.obs_settings[obs][i]
                grooming_setting = self.obs_grooming_settings[obs][i]
                obs_label = self.utils.obs_label(obs_setting, grooming_setting)
                label = ("R%s_%s" % (jetR, obs_label)).replace('.', '')

                jet_sd = None
                if grooming_setting:
                    gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
                    jet_sd = self.utils.groom(gshop, grooming_setting, jetR).pair()

                obs_val = None

                # Calculate angularities
                if obs == "ang":
                    alpha = obs_setting
                    kappa = 1
                    if grooming_setting:
                        obs_val = fjext.lambda_beta_kappa(jet, jet_sd, alpha, kappa, jetR)
                    else:
                        obs_val = fjext.lambda_beta_kappa(jet, alpha, kappa, jetR)

                # Jet mass histograms
                elif obs == "mass":
                    if grooming_setting:
                        # Untagged jets -- record underflow value
                        obs_val = jet_sd.m() if jet_sd.has_constituents() else -1

                    else:
                        obs_val = jet.m()

                else:
                    raise ValueError("Observable not implemented in fill_unmatched_jet_tree")

                tw.fill_branch("%s_%s_%s" % (obs, level, label), obs_val)

    #---------------------------------------------------------------
    # Fill jet histograms for PYTHIA run-through before matching
    #---------------------------------------------------------------
    def fill_unmatched_histograms(self, jetR, jet, level, MPI=False):

        # Observable-independent histograms
        R_label = str(jetR).replace('.', '') + 'Scaled'
        if MPI:
            getattr(self, 'hJetPt_%s_MPIon_R%s' % (level, R_label)).Fill(jet.pt())
            getattr(self, 'hNconstit_Pt_%s_MPIon_R%s' % (level, R_label)).Fill(jet.pt(), len(jet.constituents()))
        else:
            getattr(self, 'hJetPt_%s_R%s' % (level, R_label)).Fill(jet.pt())
            getattr(self, 'hNconstit_Pt_%s_R%s' % (level, R_label)).Fill(jet.pt(), len(jet.constituents()))

        for obs in self.observable_list:
            for i in range(len(self.obs_settings[obs])):

                obs_setting = self.obs_settings[obs][i]
                grooming_setting = self.obs_grooming_settings[obs][i]
                obs_label = self.utils.obs_label(obs_setting, grooming_setting)
                label = ("R%s_%s" % (jetR, obs_label)).replace('.', '')

                if MPI:
                    h = getattr(self, 'h_%s_JetPt_%s_MPIon_%sScaled' % (obs, level, label))
                else:
                    h = getattr(self, 'h_%s_JetPt_%s_%sScaled' % (obs, level, label))

                # Jet angularity histograms
                if obs == "ang":
                    alpha = obs_setting
                    kappa = 1

                    if grooming_setting:
                        gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
                        jet_sd = self.utils.groom(gshop, grooming_setting, jetR).pair()
                        h.Fill(jet.pt(), fjext.lambda_beta_kappa(jet, jet_sd, alpha, kappa, jetR))
                    else:
                        h.Fill(jet.pt(), fjext.lambda_beta_kappa(jet, alpha, kappa, jetR))

                # Jet mass histograms
                elif obs == "mass":
                    if grooming_setting:
                        gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
                        jet_sd = self.utils.groom(gshop, grooming_setting, jetR).pair()
                        if not jet_sd.has_constituents():
                            # Untagged jet -- record underflow value
                            h.Fill(jet.pt(), -1)
                        else:
                            h.Fill(jet.pt(), jet_sd.m())

                    else:
                        h.Fill(jet.pt(), jet.m())

                else:
                    raise ValueError("Observable not implemented in fill_unmatched_histograms")

    #---------------------------------------------------------------
    # Fill matched jet histograms
    #---------------------------------------------------------------
    def fill_jet_histograms(self, jetR, jp, jh, jch):

        R_label = str(jetR).replace('.', '') + 'Scaled'
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

        # Fill observable histograms and response matrices
        for alpha in self.alpha_list:
            self.fill_RMs(jetR, alpha, jp, jh, jch)


    #---------------------------------------------------------------
    # Fill jet response matrices
    #---------------------------------------------------------------
    def fill_RMs(self, jetR, alpha, jp, jh, jch):

        for obs in self.observable_list:
            for i in range(len(self.obs_settings[obs])):

                obs_setting = self.obs_settings[obs][i]
                grooming_setting = self.obs_grooming_settings[obs][i]
                obs_label = self.utils.obs_label(obs_setting, grooming_setting)
                label = ("R%s_%s" % (jetR, obs_label)).replace('.', '')

                jp_sd, jh_sd, jch_sd = None, None, None
                if grooming_setting:
                    gshop_p = fjcontrib.GroomerShop(jp, jetR, self.reclustering_algorithm)
                    jp_sd = self.utils.groom(gshop_p, grooming_setting, jetR).pair()
                    gshop_h = fjcontrib.GroomerShop(jh, jetR, self.reclustering_algorithm)
                    jh_sd = self.utils.groom(gshop_h, grooming_setting, jetR).pair()
                    gshop_ch = fjcontrib.GroomerShop(jch, jetR, self.reclustering_algorithm)
                    jch_sd = self.utils.groom(gshop_ch, grooming_setting, jetR).pair()

                obs_p, obs_h, obs_ch = None, None, None

                # Calculate angularities
                if obs == "ang":
                    alpha = obs_setting
                    kappa = 1
                    if grooming_setting:
                        obs_p = fjext.lambda_beta_kappa(jp, jp_sd, alpha, kappa, jetR)
                        obs_h = fjext.lambda_beta_kappa(jh, jh_sd, alpha, kappa, jetR)
                        obs_ch = fjext.lambda_beta_kappa(jch, jch_sd, alpha, kappa, jetR)
                    else:
                        obs_p = fjext.lambda_beta_kappa(jp, alpha, kappa, jetR)
                        obs_h = fjext.lambda_beta_kappa(jh, alpha, kappa, jetR)
                        obs_ch = fjext.lambda_beta_kappa(jch, alpha, kappa, jetR)

                # Jet mass histograms
                elif obs == "mass":
                    if grooming_setting:
                        # Untagged jets -- record underflow value
                        obs_p = jp_sd.m() if jp_sd.has_constituents() else -1
                        obs_h = jh_sd.m() if jh_sd.has_constituents() else -1
                        obs_ch = jch_sd.m() if jch_sd.has_constituents() else -1

                    else:
                        obs_p = jp.m()
                        obs_h = jh.m()
                        obs_ch = jch.m()

                else:
                    raise ValueError("Observable not implemented in fill_unmatched_histograms")


                for level in ['p', 'h', 'ch']:
                    if self.level in [None, level]:
                        getattr(self, 'h_%s_JetPt_%s_%sScaled' % (obs, level, label)).Fill(jch.pt(), locals()['obs_'+level])

                if self.level == None:
                    getattr(self, 'hResponse_%s_%sScaled' % (obs, label)).Fill(obs_p, obs_ch)

                    '''
                    # Lambda at p-vs-ch-level for various bins in ch jet pT
                    if 20 <= jch.pt() < 40:
                        getattr(self, 'hResponse_%s_PtBinCH20-40_%sScaled' % (obs, label)).Fill(lp, lch)
                    elif 40 <= jch.pt() < 60:
                        getattr(self, 'hResponse_%s_PtBinCH40-60_%sScaled' % (obs, label)).Fill(lp, lch)
                    elif 60 <= jch.pt() < 80:
                        getattr(self, 'hResponse_%s_PtBinCH60-80_%sScaled' % (obs, label)).Fill(lp, lch)

                    # Phase space plots and annulus histograms, including those binned in ch jet pT
                    num_r = self.annulus_plots_num_r
                    ang_per_r_ch = [0] * num_r
                    for particle in jch.constituents():
                        deltaR = particle.delta_R(jch)
                        getattr(self, 'hPhaseSpace_DeltaR_Pt_ch_%sScaled' % (obs, label)).Fill(
                            particle.pt(), deltaR / jetR)

                        lambda_i = lambda_beta_kappa_i(particle, jch, jetR, alpha, 1)
                        getattr(self, 'hPhaseSpace_%s_DeltaR_ch_%sScaled' % (obs, label)).Fill(deltaR / jetR, lambda_i)
                        getattr(self, 'hPhaseSpace_%s_Pt_ch_%sScaled' % (obs, label)).Fill(particle.pt(), lambda_i)

                        if 60 <= jch.pt() < 80:
                            getattr(self, 'hPhaseSpace_DeltaR_Pt_ch_PtBinCH60-80_%sScaled' % (obs, label)).Fill(
                                particle.pt(), deltaR / jetR)
                            getattr(self, 'hPhaseSpace_%s_DeltaR_ch_PtBinCH60-80_%sScaled' % (obs, label)).Fill(
                                deltaR / jetR, lambda_i)
                            getattr(self, 'hPhaseSpace_%s_Pt_ch_PtBinCH60-80_%sScaled' % (obs, label)).Fill(
                                particle.pt(), lambda_i)

                        ang_per_r_ch = [ang_per_r_ch[i] + lambda_i * 
                                        (deltaR <= ((i+1) * jetR * self.annulus_plots_max_x / num_r))
                                        for i in range(0, num_r, 1)]

                    ang_per_r_p = [0] * num_r
                    for particle in jp.constituents():
                        deltaR = particle.delta_R(jp)
                        getattr(self, 'hPhaseSpace_DeltaR_Pt_p_%sScaled' % (obs, label)).Fill(
                            particle.pt(), deltaR / jetR)

                        lambda_i = lambda_beta_kappa_i(particle, jp, jetR, alpha, 1)
                        getattr(self, 'hPhaseSpace_%s_DeltaR_p_%sScaled' % (obs, label)).Fill(deltaR / jetR, lambda_i)
                        getattr(self, 'hPhaseSpace_%s_Pt_p_%sScaled' % (obs, label)).Fill(particle.pt(), lambda_i)

                        if 60 <= jch.pt() < 80:
                            getattr(self, 'hPhaseSpace_DeltaR_Pt_p_PtBinCH60-80_%sScaled' % (obs, label)).Fill(
                                particle.pt(), deltaR / jetR)
                            getattr(self, 'hPhaseSpace_%s_DeltaR_p_PtBinCH60-80_%sScaled' % (obs, label)).Fill(
                                deltaR / jetR, lambda_i)
                            getattr(self, 'hPhaseSpace_%s_Pt_p_PtBinCH60-80_%sScaled' % (obs, label)).Fill(
                                particle.pt(), lambda_i)

                        ang_per_r_p = [ang_per_r_p[i] + lambda_i *
                                       (deltaR <= ((i+1) * jetR * self.annulus_plots_max_x / num_r))
                                     for i in range(0, num_r, 1)]

                    for i in range(0, num_r, 1):
                        getattr(self, 'hAnnulus_%s_p_%sScaled' % (obs, label)).Fill(
                            (i+1) * self.annulus_plots_max_x / num_r, ang_per_r_p[i] / (lp + 1e-11))
                        getattr(self, 'hAnnulus_%s_ch_%sScaled' % (obs, label)).Fill(
                            (i+1) * self.annulus_plots_max_x / num_r, ang_per_r_ch[i] / (lch + 1e-11))
                        if 60 <= jch.pt() < 80:
                            getattr(self, 'hAnnulus_%s_p_PtBinCH60-80_%sScaled' % (obs, label)).Fill(
                                (i+1) * self.annulus_plots_max_x / num_r, ang_per_r_p[i] / (lp + 1e-11))
                            getattr(self, 'hAnnulus_%s_ch_PtBinCH60-80_%sScaled' % (obs, label)).Fill(
                                (i+1) * self.annulus_plots_max_x / num_r, ang_per_r_ch[i] / (lch + 1e-11))
                    '''

                    # Residual plots (with and without divisor in y-axis)
                    getattr(self, "h_%sDiff_JetPt_%sScaled" % (obs, label)).Fill(jch.pt(), obs_p - obs_ch)
                    if obs_p:  # prevent divide by 0
                        getattr(self, "h_%sResidual_JetPt_%sScaled" % (obs, label)).Fill(jp.pt(), (obs_p - obs_ch) / obs_p)

                    # 4D response matrices for "forward folding" to ch level
                    x = ([jch.pt(), jp.pt(), obs_ch, obs_p])
                    x_array = array('d', x)
                    getattr(self, 'hResponse_JetPt_%s_ch_%sScaled' % (obs, label)).Fill(x_array)

                    x = ([jh.pt(), jp.pt(), obs_h, obs_p])
                    x_array = array('d', x)
                    getattr(self, 'hResponse_JetPt_%s_h_%sScaled' % (obs, label)).Fill(x_array)

    #---------------------------------------------------------------
    # Fill jet histograms for MPI (which are just the H-->CH RMs)
    #---------------------------------------------------------------
    def fill_jet_histograms_MPI(self, jetR, jp, jh, jch):

        for obs in self.observable_list:
            for i in range(len(self.obs_settings[obs])):

                obs_setting = self.obs_settings[obs][i]
                grooming_setting = self.obs_grooming_settings[obs][i]
                obs_label = self.utils.obs_label(obs_setting, grooming_setting)

                jp_sd, jh_sd, jch_sd = None, None, None
                if grooming_setting:
                    gshop_p = fjcontrib.GroomerShop(jp, jetR, self.reclustering_algorithm)
                    jp_sd = self.utils.groom(gshop_p, grooming_setting, jetR).pair()
                    gshop_h = fjcontrib.GroomerShop(jh, jetR, self.reclustering_algorithm)
                    jh_sd = self.utils.groom(gshop_h, grooming_setting, jetR).pair()
                    gshop_ch = fjcontrib.GroomerShop(jch, jetR, self.reclustering_algorithm)
                    jch_sd = self.utils.groom(gshop_ch, grooming_setting, jetR).pair()

                obs_p, obs_h, obs_ch = None, None, None

                # Calculate angularities
                if obs == "ang":
                    alpha = obs_setting
                    kappa = 1
                    if grooming_setting:
                        obs_p = fjext.lambda_beta_kappa(jp, jp_sd, alpha, kappa, jetR)
                        obs_h = fjext.lambda_beta_kappa(jh, jh_sd, alpha, kappa, jetR)
                        obs_ch = fjext.lambda_beta_kappa(jch, jch_sd, alpha, kappa, jetR)
                    else:
                        obs_p = fjext.lambda_beta_kappa(jp, alpha, kappa, jetR)
                        obs_h = fjext.lambda_beta_kappa(jh, alpha, kappa, jetR)
                        obs_ch = fjext.lambda_beta_kappa(jch, alpha, kappa, jetR)

                # Jet mass histograms
                elif obs == "mass":
                    if grooming_setting:
                        # Untagged jets -- record underflow value
                        obs_p = jp_sd.m() if jp_sd.has_constituents() else -1
                        obs_h = jh_sd.m() if jh_sd.has_constituents() else -1
                        obs_ch = jch_sd.m() if jch_sd.has_constituents() else -1

                    else:
                        obs_p = jp.m()
                        obs_h = jh.m()
                        obs_ch = jch.m()

                else:
                    raise ValueError("Observable not implemented in fill_unmatched_histograms")

                # 4D response matrices for "forward folding" from h to ch level
                x = ([jch.pt(), jh.pt(), obs_ch, obs_h])
                x_array = array('d', x)
                getattr(self, 'hResponse_JetPt_%s_Fnp_%sScaled' % (obs, label)).Fill(x_array)

    #---------------------------------------------------------------
    # Initiate scaling of all histograms and print final simulation info
    #---------------------------------------------------------------
    def scale_print_final_info(self, pythia, pythia_MPI):

        # Scale all jet histograms by the appropriate factor from generated cross section
        # and the number of accepted events
        if not self.no_scale:
            scale_f = pythia.info.sigmaGen() / self.hNevents.GetBinContent(1)
            print("Weight MPIoff histograms by (cross section)/(N events) =", scale_f)
            MPI_scale_f = pythia_MPI.info.sigmaGen() / self.hNeventsMPI.GetBinContent(1)
            print("Weight MPIon histograms by (cross section)/(N events) =", MPI_scale_f)
            self.scale_jet_histograms(scale_f, MPI_scale_f)
        print()

        print("N total final MPI-off events:", int(self.hNevents.GetBinContent(1)), "with",
              int(pythia.info.nAccepted() - self.hNevents.GetBinContent(1)),
              "events rejected at hadronization step")
        self.hNevents.SetBinError(1, 0)

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
    parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly',
                                     prog=os.path.basename(__file__))
    pyconf.add_standard_pythia_args(parser)
    # Could use --py-seed
    parser.add_argument('--user-seed', help='PYTHIA starting seed', default=1111, type=int)
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
    parser.add_argument('-c', '--config_file', action='store', type=str, default='config/angularity.yaml',
                        help="Path of config file for observable configurations")
    args = parser.parse_args()

    if args.no_match_level not in [None, 'p', 'h', 'ch']:
        print("ERROR: Unrecognized type %s. Please use 'p', 'h', or 'ch'" % args.type_only)
        exit(1)

    # If invalid configFile is given, exit
    if not os.path.exists(args.config_file):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    # Use PYTHIA seed for event generation
    if args.user_seed < 0:
        args.user_seed = 1111

    # Have at least 1 event
    if args.nev < 1:
        args.nev = 1

    if args.py_noMPI:
        print("\033[91m%s\033[00m" % "WARNING: py-noMPI flag ignored for this program")
        time.sleep(3)
        print()

    process = pythia_parton_hadron(config_file=args.config_file, output_dir=args.output_dir, args=args)
    process.pythia_parton_hadron(args)
