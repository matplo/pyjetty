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
import array
import numpy as np

from pyjetty.mputils import *

from pyjetty.alice_analysis.process.base import process_base
from pyjetty.alice_analysis.process.user.ang_pp.helpers \
    import lambda_beta_kappa, lambda_beta_kappa_i

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

        # Defaults to None if not in use
        self.level = args.no_match_level

        self.jetR_list = config["jetR"]
        self.beta_list = config["betas"]

        # SoftDrop parameters
        self.use_SD = False   # Change this to use SD
        self.sd_beta = config["sd_beta"]
        self.sd_zcut = config["sd_zcut"]

        self.n_pt_bins = config["n_pt_bins"]
        self.pt_limits = config["pt_limits"]
        self.pTbins = np.arange(self.pt_limits[0], self.pt_limits[1] + 1, 
                                (self.pt_limits[1] - self.pt_limits[0]) / self.n_pt_bins)
        self.n_lambda_bins = config['n_lambda_bins']
        self.lambda_limits = config['lambda_limits']

        # hadron level - ALICE tracking restriction
        self.max_eta_hadron = 0.9

        # Whether or not to rescale final jet histograms based on sigma/N
        self.no_scale = args.no_scale


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

        # print the banner first
        fj.ClusterSequence.print_banner()
        print()

        self.init_jet_tools()
        self.parse_events()

        for jetR in self.jetR_list:
            getattr(self, "tw_R%s" % str(jetR).replace('.', '')).fill_tree()

        # TODO: scale histograms
        #self.scale_print_final_info(pythia, pythia_MPI)

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

            R_label = str(jetR).replace('.', '') + 'Scaled'

            if self.level in [None, 'ch']:
                name = 'hJetPt_ch_R%s' % R_label
                h = ROOT.TH1F(name, name+';p_{T}^{ch jet};#frac{dN}{dp_{T}^{ch jet}};', 300, 0, 300)
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hNconstit_Pt_ch_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{ch jet}')
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

            if self.level in [None, 'h']:
                name = 'hJetPt_h_R%s' % R_label
                h = ROOT.TH1F(name, name+';p_{T}^{jet, h};#frac{dN}{dp_{T}^{jet, h}};', 300, 0, 300)
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hNconstit_Pt_h_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{h jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{h jet}')
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

            if self.level in [None, 'p']:
                name = 'hJetPt_p_R%s' % R_label
                h = ROOT.TH1F(name, name+';p_{T}^{jet, parton};#frac{dN}{dp_{T}^{jet, parton}};',
                              300, 0, 300)
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hNconstit_Pt_p_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{p jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{p jet}')
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

            if self.level == None:
                name = 'hJetPtRes_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 200, -1., 1.)
                h.GetXaxis().SetTitle('#it{p}_{T}^{parton jet}')
                h.GetYaxis().SetTitle(
                    '#frac{#it{p}_{T}^{parton jet}-#it{p}_{T}^{ch jet}}{#it{p}_{T}^{parton jet}}')
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hResponse_JetPt_R%s' % R_label
                h = ROOT.TH2F(name, name, 200, 0, 200, 200, 0, 200)
                h.GetXaxis().SetTitle('#it{p}_{T}^{parotn jet}')
                h.GetYaxis().SetTitle('#it{p}_{T}^{ch jet}')
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                # Jet multiplicity for matched jets with a cut at ch-jet level
                name = 'hNconstit_Pt_ch_PtBinCH60-80_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{ch jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{ch jet}')
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hNconstit_Pt_h_PtBinCH60-80_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{h jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{h jet}')
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

                name = 'hNconstit_Pt_p_PtBinCH60-80_R%s' % R_label
                h = ROOT.TH2F(name, name, 300, 0, 300, 50, 0.5, 50.5)
                h.GetXaxis().SetTitle('#it{p}_{T}^{parton jet}')
                h.GetYaxis().SetTitle('#it{N}_{constit}^{parton jet}')
                setattr(self, name, h)
                getattr(self, hist_list_name).append(h)

            for beta in self.beta_list:

                label = ("R%s_%sScaled" % (str(jetR), str(beta))).replace('.', '')

                if self.level in [None, 'ch']:
                    name = 'hAng_JetPt_ch_%s' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('p_{T}^{ch jet}')
                    h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{#beta=%s}^{ch}}' % str(beta))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hAng_JetPt_ch_MPIon_%s' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('p_{T}^{ch jet}')
                    h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{#beta=%s}^{ch}}' % str(beta))
                    setattr(self, name, h)

                if self.level in [None, 'h']:
                    name = 'hAng_JetPt_h_%s' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('p_{T}^{jet, h}')
                    h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{#beta=%s}^{h}}' % str(beta))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                if self.level in [None, 'p']:
                    name = 'hAng_JetPt_p_%s' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('p_{T}^{jet, parton}')
                    h.GetYaxis().SetTitle('#frac{dN}{d#lambda_{#beta=%s}^{parton}}' % str(beta))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                if self.level == None:
                    name = 'hResponse_ang_%s' % label
                    h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
                    h.GetXaxis().SetTitle('#lambda_{#beta=%s}^{parton}' % beta)
                    h.GetYaxis().SetTitle('#lambda_{#beta=%s}^{ch}' % beta)
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hResponse_ang_PtBinCH20-40_%s' % label
                    h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
                    h.GetXaxis().SetTitle('#lambda_{#beta=%s}^{parton}' % beta)
                    h.GetYaxis().SetTitle('#lambda_{#beta=%s}^{ch}' % beta)
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hResponse_ang_PtBinCH40-60_%s' % label
                    h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
                    h.GetXaxis().SetTitle('#lambda_{#beta=%s}^{parton}' % beta)
                    h.GetYaxis().SetTitle('#lambda_{#beta=%s}^{ch}' % beta)
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hResponse_ang_PtBinCH60-80_%s' % label
                    h = ROOT.TH2F(name, name, 100, 0, 1, 100, 0, 1)
                    h.GetXaxis().SetTitle('#lambda_{#beta=%s}^{parton}' % beta)
                    h.GetYaxis().SetTitle('#lambda_{#beta=%s}^{ch}' % beta)
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    # Phase space plots integrated over all pT bins
                    name = 'hPhaseSpace_DeltaR_Pt_ch_%s' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  150, 0, 1.5)
                    h.GetXaxis().SetTitle('(p_{T, i})_{ch jet}')
                    h.GetYaxis().SetTitle('(#Delta R_{i})_{ch jet} / R')
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_DeltaR_ch_%s' % label
                    h = ROOT.TH2F(name, name, 150, 0, 1.5,
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(#Delta R_{i})_{ch jet} / R')
                    h.GetYaxis().SetTitle('(#lambda_{#beta=%s, i})_{ch jet}' % str(beta))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_Pt_ch_%s' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(p_{T, i})_{ch jet}')
                    h.GetYaxis().SetTitle('(#lambda_{#beta=%s, i})_{ch jet}' % str(beta))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_DeltaR_Pt_p_%s' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  150, 0, 1.5)
                    h.GetXaxis().SetTitle('(p_{T, i})_{parton jet}')
                    h.GetYaxis().SetTitle('(#Delta R_{i})_{parton jet} / R')
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_DeltaR_p_%s' % label
                    h = ROOT.TH2F(name, name, 150, 0, 1.5,
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(#Delta R_{i})_{parton jet} / R')
                    h.GetYaxis().SetTitle('(#lambda_{#beta=%s, i})_{parton jet}' % str(beta))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_Pt_p_%s' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(p_{T, i})_{parton jet}')
                    h.GetYaxis().SetTitle('(#lambda_{#beta=%s, i})_{parton jet}' % str(beta))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    # Phase space plots binned in ch jet pT
                    name = 'hPhaseSpace_DeltaR_Pt_ch_PtBinCH60-80_%s' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  150, 0, 1.5)
                    h.GetXaxis().SetTitle('(p_{T, i})_{ch jet}')
                    h.GetYaxis().SetTitle('(#Delta R_{i})_{ch jet} / R')
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_DeltaR_Pt_p_PtBinCH60-80_%s' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  150, 0, 1.5)
                    h.GetXaxis().SetTitle('(p_{T, i})_{parton jet}')
                    h.GetYaxis().SetTitle('(#Delta R_{i})_{parton jet} / R')
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_DeltaR_ch_PtBinCH60-80_%s' % label
                    h = ROOT.TH2F(name, name, 150, 0, 1.5,
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(#Delta R_{i})_{ch jet} / R')
                    h.GetYaxis().SetTitle('(#lambda_{#beta=%s, i})_{ch jet}' % str(beta))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_DeltaR_p_PtBinCH60-80_%s' % label
                    h = ROOT.TH2F(name, name, 150, 0, 1.5,
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(#Delta R_{i})_{parton jet} / R')
                    h.GetYaxis().SetTitle('(#lambda_{#beta=%s, i})_{parton jet}' % str(beta))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_Pt_ch_PtBinCH60-80_%s' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(p_{T, i})_{ch jet}')
                    h.GetYaxis().SetTitle('(#lambda_{#beta=%s, i})_{ch jet}' % str(beta))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hPhaseSpace_ang_Pt_p_PtBinCH60-80_%s' % label
                    h = ROOT.TH2F(name, name, self.n_pt_bins, self.pt_limits[0], self.pt_limits[1],
                                  self.n_lambda_bins, self.lambda_limits[0], self.lambda_limits[1])
                    h.GetXaxis().SetTitle('(p_{T, i})_{parton jet}')
                    h.GetYaxis().SetTitle('(#lambda_{#beta=%s, i})_{parton jet}' % str(beta))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    # Annulus plots for amount of lambda contained within some r < R
                    self.annulus_plots_num_r = 150
                    self.annulus_plots_max_x = 1.5
                    low_bound = self.annulus_plots_max_x / self.annulus_plots_num_r / 2.
                    up_bound = self.annulus_plots_max_x + low_bound

                    name = 'hAnnulus_ang_ch_%s' % label
                    h = ROOT.TH2F(name, name, self.annulus_plots_num_r, low_bound, up_bound,
                                  100, 0, 1.)
                    h.GetXaxis().SetTitle('(#it{r} / #it{R})_{ch jet}')
                    h.GetYaxis().SetTitle(
                        ('(#frac{#lambda_{#beta=%s}(#it{r})}' + \
                         '{#lambda_{#beta=%s}(#it{R})})_{ch jet}') % (str(beta), str(beta)))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hAnnulus_ang_ch_PtBinCH60-80_%s' % label
                    h = ROOT.TH2F(name, name, self.annulus_plots_num_r, low_bound, up_bound,
                                  100, 0, 1.)
                    h.GetXaxis().SetTitle('(#it{r} / #it{R})_{ch jet}')
                    h.GetYaxis().SetTitle(
                        ('(#frac{#lambda_{#beta=%s}(#it{r})}' + \
                         '{#lambda_{#beta=%s}(#it{R})})_{ch jet}') % (str(beta), str(beta)))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hAnnulus_ang_p_%s' % label
                    h = ROOT.TH2F(name, name, self.annulus_plots_num_r, low_bound, up_bound,
                                  100, 0, 1.)
                    h.GetXaxis().SetTitle('(#it{r} / #it{R})_{parton jet}')
                    h.GetYaxis().SetTitle(
                        ('(#frac{#lambda_{#beta=%s}(#it{r})}' + \
                         '{#lambda_{#beta=%s}(#it{R})})_{parton jet}') % (str(beta), str(beta)))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = 'hAnnulus_ang_p_PtBinCH60-80_%s' % label
                    h = ROOT.TH2F(name, name, self.annulus_plots_num_r, low_bound, up_bound,
                                  100, 0, 1.)
                    h.GetXaxis().SetTitle('(#it{r} / #it{R})_{parton jet}')
                    h.GetYaxis().SetTitle(
                        ('(#frac{#lambda_{#beta=%s}(#it{r})}' + \
                         '{#lambda_{#beta=%s}(#it{R})})_{parton jet}') % (str(beta), str(beta)))
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = "hAngResidual_JetPt_%s" % label
                    h = ROOT.TH2F(name, name, 300, 0, 300, 200, -3., 1.)
                    h.GetXaxis().SetTitle('p_{T}^{jet, parton}')
                    h.GetYaxis().SetTitle('#frac{#lambda_{#beta}^{jet, parton}-#lambda_{#beta}' + \
                                          '^{ch jet}}{#lambda_{#beta}^{jet, parton}}')
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    name = "hAngDiff_JetPt_%s" % label
                    h = ROOT.TH2F(name, name, 300, 0, 300, 200, -2., 2.)
                    h.GetXaxis().SetTitle('#it{p}_{T}^{jet, ch}')
                    h.GetYaxis().SetTitle('#it{#lambda}_{#it{#beta}}^{jet, parton}-' + \
                                          '#it{#lambda}_{#it{#beta}}^{jet, ch}')
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    # Create THn of response
                    dim = 4
                    title = ['p_{T}^{ch jet}', 'p_{T}^{parton jet}', 
                             '#lambda_{#beta}^{ch}', '#lambda_{#beta}^{parton}']
                    nbins = [14, 14, 100, 101]
                    min_li = [ 10.,  10., 0.,  -0.005]
                    max_li = [150., 150., 1.0,  1.005]

                    name = 'hResponse_JetPt_ang_ch_%s' % label
                    nbins = (nbins)
                    xmin = (min_li)
                    xmax = (max_li)
                    nbins_array = array.array('i', nbins)
                    xmin_array = array.array('d', xmin)
                    xmax_array = array.array('d', xmax)
                    h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
                    for i in range(0, dim):
                        h.GetAxis(i).SetTitle(title[i])
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)

                    # Another set of THn for full hadron folding
                    title = ['p_{T}^{h jet}', 'p_{T}^{parton jet}', 
                             '#lambda_{#beta}^{h}', '#lambda_{#beta}^{parton}']

                    name = 'hResponse_JetPt_ang_h_%s' % label
                    h = ROOT.THnF(name, name, dim, nbins_array, xmin_array, xmax_array)
                    for i in range(0, dim):
                        h.GetAxis(i).SetTitle(title[i])
                    setattr(self, name, h)
                    getattr(self, hist_list_name).append(h)


    #---------------------------------------------------------------
    # Initiate jet defs, selectors, and sd (if required)
    #---------------------------------------------------------------
    def init_jet_tools(self):
        
        for jetR in self.jetR_list:
            jetR_str = str(jetR).replace('.', '')
            
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

        # SoftDrop settings
        if self.use_SD:
            args.output = args.output.replace('.root', '_sdbeta{}.root'.format(args.sd_beta))
        for jetR in self.jetR_list:
            sd = fjcontrib.SoftDrop(args.sd_beta, 0.1, jetR) if self.use_SD else None
            setattr(self, "sd_R%s" % str(jetR).replace('.', ''), sd)


    #---------------------------------------------------------------
    # Read events from output, find jets, and fill histograms
    #---------------------------------------------------------------
    def parse_events(self, MPIon=False):

        partons_per_event, hadrons_per_event, ch_hadrons_per_event = self.parse_events_from_file(MPIon)
        print()

        self.find_jets_fill_hist(partons_per_event, hadrons_per_event, ch_hadrons_per_event, MPIon)


    #---------------------------------------------------------------
    # Read events from output, find jets, and fill histograms
    #---------------------------------------------------------------
    def parse_events_from_file(self, MPIon=False):

        if MPIon:
            hNevents = self.hNeventsMPI
        else:
            hNevents = self.hNevents

        print("Reading events from %s..." % self.herwig_file)

        with open(self.herwig_file, 'r') as f:
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

            partons = []
            hadrons = []
            ch_hadrons = []

            for line in f:

                # Waiting to start reading event
                if not reading_ev:
                    if "Event number" in line:
                        reading_ev = True
                        ev_num = int(line.split()[2])
                        if not ev_num % 100:
                            print("Event number", ev_num, end="\r")
                        hNevents.Fill(0)
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
                elif not parton_finished:
                    # Read parton information
                    vals = line.split()
                    if line[0] == '-':
                        parton_finished = True
                    elif len(vals) == 5 and line[0] == ' ':
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
                    partons.append(fjext.vectorize_px_py_pz_e(
                        partons_px, partons_py, partons_pz, partons_e))

                    hadrons.append(fjext.vectorize_px_py_pz_e(
                        hadrons_px, hadrons_py, hadrons_pz, hadrons_e))

                    ch_hadrons.append(fjext.vectorize_px_py_pz_e(
                        ch_hadrons_px, ch_hadrons_py, ch_hadrons_pz, ch_hadrons_e))

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

                elif line[0].isnumeric():
                    # Save the hadron name for charge descrimination
                    i = 1
                    while line[i].isnumeric():
                        i += 1
                    hadron_type = line[i:].split()[0]
                    continue

                elif line[0] == line[1] == ' ':  # and len(vals) == 5:
                    # Reading hadron information
                    vals = line.split()
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
    def find_jets_fill_hist(self, partons_per_event, hadrons_per_event,
                            ch_hadrons_per_event, MPIon=False):

        for jetR in self.jetR_list:
            print("Filling jet histograms for R = %s..." % str(jetR))
            
            jetR_str = str(jetR).replace('.', '')
            jet_selector = getattr(self, "jet_selector_R%s" % jetR_str)
            jet_def = getattr(self, "jet_def_R%s" % jetR_str)
            sd = getattr(self, "sd_R%s" % jetR_str)
            t = getattr(self, "t_R%s" % jetR_str)
            tw = getattr(self, "tw_R%s" % jetR_str)
            count1 = getattr(self, "count1_R%s" % jetR_str)
            count2 = getattr(self, "count2_R%s" % jetR_str)

            for iev in range(len(partons_per_event)):
                if not (iev+1) % 100:
                    print("Event number %s" % str(iev+1), end='\r')
                
                jets_p = fj.sorted_by_pt(jet_selector(jet_def(partons_per_event[iev])))
                jets_h = fj.sorted_by_pt(jet_selector(jet_def(hadrons_per_event[iev])))
                jets_ch = fj.sorted_by_pt(jet_selector(jet_def(ch_hadrons_per_event[iev])))

                if MPIon:
                    for jet in jets_ch:
                        self.fill_MPI_histograms(jetR, jet, sd)
                    continue

                if self.level:  # Only save info at one level w/o matching
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

                            self.fill_matched_jet_tree(tw, jetR, iev, jp, jh, jchh, sd)
                            self.fill_jet_histograms(jetR, jp, jh, jchh, sd)

                    #print("  |-> SD jet params z={0:10.3f} dR={1:10.3f} mu={2:10.3f}".format(
                    #    sd_info.z, sd_info.dR, sd_info.mu))

            setattr(self, "count1_R%s" % jetR_str, count1)
            setattr(self, "count2_R%s" % jetR_str, count2)



    #---------------------------------------------------------------
    # Fill jet tree with (unscaled/raw) matched parton/hadron tracks
    #---------------------------------------------------------------
    def fill_matched_jet_tree(self, tw, jetR, iev, jp, jh, jchh, sd=None):

        tw.fill_branch('iev', iev)
        tw.fill_branch('ch', jchh)
        tw.fill_branch('h', jh)
        tw.fill_branch('p', jp)

        kappa = 1
        for beta in self.beta_list:
            label = str(beta).replace('.', '')
            tw.fill_branch("l_ch_%s" % label,
                           lambda_beta_kappa(jchh, jetR, beta, kappa))
            tw.fill_branch("l_h_%s" % label,
                           lambda_beta_kappa(jh, jetR, beta, kappa))
            tw.fill_branch("l_p_%s" % label,
                           lambda_beta_kappa(jp, jetR, beta, kappa))

        # Save SoftDrop variables as well if desired
        if self.use_SD:

            if sd == None:
                print("ERROR: NoneType passed as SoftDrop groomer.")
                exit(1)

            jchh_sd = sd.result(jchh)
            jchh_sd_info = fjcontrib.get_SD_jet_info(jchh_sd)
            jh_sd = sd.result(jh)
            jh_sd_info = fjcontrib.get_SD_jet_info(jh_sd)
            jp_sd = sd.result(jp)
            jp_sd_info = fjcontrib.get_SD_jet_info(jp_sd)

            tw.fill_branch('p_zg', jp_sd_info.z)
            tw.fill_branch('p_Rg', jp_sd_info.dR)
            tw.fill_branch('p_thg', jp_sd_info.dR/jetR)
            tw.fill_branch('p_mug', jp_sd_info.mu)

            tw.fill_branch('h_zg', jh_sd_info.z)
            tw.fill_branch('h_Rg', jh_sd_info.dR)
            tw.fill_branch('h_thg', jh_sd_info.dR/jetR)
            tw.fill_branch('h_mug', jh_sd_info.mu)

            tw.fill_branch('ch_zg', jchh_sd_info.z)
            tw.fill_branch('ch_Rg', jchh_sd_info.dR)
            tw.fill_branch('ch_thg', jchh_sd_info.dR/jetR)
            tw.fill_branch('ch_mug', jchh_sd_info.mu)


    #---------------------------------------------------------------
    # Fill jet tree with (unscaled/raw) unmatched parton/hadron tracks
    #---------------------------------------------------------------
    def fill_unmatched_jet_tree(self, tw, jetR, iev, jet):

        tw.fill_branch('iev', iev)
        tw.fill_branch(self.level, jet)

        kappa = 1
        for beta in self.beta_list:
            label = str(beta).replace('.', '')
            tw.fill_branch('l_%s_%s' % (self.level, label),
                           lambda_beta_kappa(jet, jetR, beta, kappa))

    
    #---------------------------------------------------------------
    # Fill jet histograms for MPI-on PYTHIA run-through
    #---------------------------------------------------------------
    def fill_MPI_histograms(self, jetR, jet, sd):

        for beta in self.beta_list:
            label = ("R%s_%sScaled" % (str(jetR), str(beta))).replace('.', '')
            h = getattr(self, 'hAng_JetPt_ch_MPIon_%s' % label)

            kappa = 1
            h.Fill(jet.pt(), lambda_beta_kappa(jet, jetR, beta, kappa))
            
            if sd != None:
                pass  # TODO


    #---------------------------------------------------------------
    # Fill jet histograms
    #---------------------------------------------------------------
    def fill_jet_histograms(self, jetR, jp, jh, jch, sd):

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
            
            if 60 <= jch.pt() < 80:
                getattr(self, 'hNconstit_Pt_ch_PtBinCH60-80_R%s' % R_label).Fill(
                    jch.pt(), len(jch.constituents()))
                getattr(self, 'hNconstit_Pt_h_PtBinCH60-80_R%s' % R_label).Fill(
                    jh.pt(), len(jh.constituents()))
                getattr(self, 'hNconstit_Pt_p_PtBinCH60-80_R%s' % R_label).Fill(
                    jp.pt(), len(jp.constituents()))

        # Fill angularity histograms and response matrices
        for beta in self.beta_list:
            self.fill_RMs(jetR, beta, jp, jh, jch, sd)


    #---------------------------------------------------------------
    # Fill jet histograms
    #---------------------------------------------------------------
    def fill_RMs(self, jetR, beta, jp, jh, jch, sd):

        # Calculate angularities
        kappa = 1
        lp = lambda_beta_kappa(jp, jetR, beta, kappa)
        lh = lambda_beta_kappa(jh, jetR, beta, kappa)
        lch = lambda_beta_kappa(jch, jetR, beta, kappa)

        label = ("R%s_%sScaled" % (str(jetR), str(beta))).replace('.', '')

        if self.level in [None, 'ch']:
            getattr(self, 'hAng_JetPt_ch_%s' % label).Fill(jch.pt(), lch)

        if self.level in [None, 'h']:
            getattr(self, 'hAng_JetPt_h_%s' % label).Fill(jh.pt(), lh)

        if self.level in [None, 'p']:
            getattr(self, 'hAng_JetPt_p_%s' % label).Fill(jp.pt(), lp)

        if self.level == None:
            getattr(self, 'hResponse_ang_%s' % label).Fill(lp, lch)

            # Lambda at p-vs-ch-level for various bins in ch jet pT 
            if 20 <= jch.pt() < 40:
                getattr(self, 'hResponse_ang_PtBinCH20-40_%s' % label).Fill(lp, lch)
            elif 40 <= jch.pt() < 60:
                getattr(self, 'hResponse_ang_PtBinCH40-60_%s' % label).Fill(lp, lch)
            elif 60 <= jch.pt() < 80:
                getattr(self, 'hResponse_ang_PtBinCH60-80_%s' % label).Fill(lp, lch)

            # Phase space plots and annulus histograms, including those binned in ch jet pT
            num_r = self.annulus_plots_num_r
            ang_per_r_ch = [0] * num_r
            for particle in jch.constituents():
                deltaR = particle.delta_R(jch)
                getattr(self, 'hPhaseSpace_DeltaR_Pt_ch_%s' % label).Fill(
                    particle.pt(), deltaR / jetR)

                lambda_i = lambda_beta_kappa_i(particle, jch, jetR, beta, 1)
                getattr(self, 'hPhaseSpace_ang_DeltaR_ch_%s' % label).Fill(deltaR / jetR, lambda_i)
                getattr(self, 'hPhaseSpace_ang_Pt_ch_%s' % label).Fill(particle.pt(), lambda_i)

                if 60 <= jch.pt() < 80:
                    getattr(self, 'hPhaseSpace_DeltaR_Pt_ch_PtBinCH60-80_%s' % label).Fill(
                        particle.pt(), deltaR / jetR)
                    getattr(self, 'hPhaseSpace_ang_DeltaR_ch_PtBinCH60-80_%s' % label).Fill(
                        deltaR / jetR, lambda_i)
                    getattr(self, 'hPhaseSpace_ang_Pt_ch_PtBinCH60-80_%s' % label).Fill(
                        particle.pt(), lambda_i)

                ang_per_r_ch = [ang_per_r_ch[i] + lambda_i * 
                                (deltaR <= ((i+1) * jetR * self.annulus_plots_max_x / num_r))
                                for i in range(0, num_r, 1)]

            ang_per_r_p = [0] * num_r
            for particle in jp.constituents():
                deltaR = particle.delta_R(jp)
                getattr(self, 'hPhaseSpace_DeltaR_Pt_p_%s' % label).Fill(
                    particle.pt(), deltaR / jetR)

                lambda_i = lambda_beta_kappa_i(particle, jp, jetR, beta, 1)
                getattr(self, 'hPhaseSpace_ang_DeltaR_p_%s' % label).Fill(deltaR / jetR, lambda_i)
                getattr(self, 'hPhaseSpace_ang_Pt_p_%s' % label).Fill(particle.pt(), lambda_i)

                if 60 <= jch.pt() < 80:
                    getattr(self, 'hPhaseSpace_DeltaR_Pt_p_PtBinCH60-80_%s' % label).Fill(
                        particle.pt(), deltaR / jetR)
                    getattr(self, 'hPhaseSpace_ang_DeltaR_p_PtBinCH60-80_%s' % label).Fill(
                        deltaR / jetR, lambda_i)
                    getattr(self, 'hPhaseSpace_ang_Pt_p_PtBinCH60-80_%s' % label).Fill(
                        particle.pt(), lambda_i)

                ang_per_r_p = [ang_per_r_p[i] + lambda_i *
                               (deltaR <= ((i+1) * jetR * self.annulus_plots_max_x / num_r))
                             for i in range(0, num_r, 1)]

            for i in range(0, num_r, 1):
                getattr(self, 'hAnnulus_ang_p_%s' % label).Fill(
                    (i+1) * self.annulus_plots_max_x / num_r, ang_per_r_p[i] / (lp + 1e-11))
                getattr(self, 'hAnnulus_ang_ch_%s' % label).Fill(
                    (i+1) * self.annulus_plots_max_x / num_r, ang_per_r_ch[i] / (lch + 1e-11))
                if 60 <= jch.pt() < 80:
                    getattr(self, 'hAnnulus_ang_p_PtBinCH60-80_%s' % label).Fill(
                        (i+1) * self.annulus_plots_max_x / num_r, ang_per_r_p[i] / (lp + 1e-11))
                    getattr(self, 'hAnnulus_ang_ch_PtBinCH60-80_%s' % label).Fill(
                        (i+1) * self.annulus_plots_max_x / num_r, ang_per_r_ch[i] / (lch + 1e-11))

            # Residual plots (with and without divisor in y-axis)
            getattr(self, "hAngDiff_JetPt_%s" % label).Fill(jch.pt(), lp - lch)
            if lp:  # prevent divide by 0
                getattr(self, "hAngResidual_JetPt_%s" % label).Fill(jp.pt(), (lp - lch) / lp)

            # 4D response matrices for "forward folding" to ch level
            x = ([jch.pt(), jp.pt(), lch, lp])
            x_array = array.array('d', x)
            getattr(self, 'hResponse_JetPt_ang_ch_%s' % label).Fill(x_array)

            x = ([jh.pt(), jp.pt(), lh, lp])
            x_array = array.array('d', x)
            getattr(self, 'hResponse_JetPt_ang_h_%s' % label).Fill(x_array)


    #---------------------------------------------------------------
    # Initiate scaling of all histograms and print final simulation info
    #---------------------------------------------------------------
    def scale_print_final_info(self, pythia, pythia_MPI):

        # Scale all jet histograms by the appropriate factor from generated cross section
        # and the number of accepted events
        if not self.no_scale:
            scale_f = pythia.info.sigmaGen() / self.hNevents.GetBinContent(1)
            print("Weight MPIoff tree by (cross section)/(N events) =", scale_f)
            MPI_scale_f = pythia_MPI.info.sigmaGen() / self.hNeventsMPI.GetBinContent(1)
            print("Weight MPIon tree by (cross section)/(N events) =", MPI_scale_f)
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

            for beta in self.beta_list:
                label = ("R%s_%sScaled" % (str(jetR), str(beta))).replace('.', '')
                getattr(self, 'hAng_JetPt_ch_MPIon_%s' % label).Scale(MPI_scale_f)


################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Herwig7 debug level-1 output parser',
                                     prog=os.path.basename(__file__))
    parser.add_argument('-i', '--input-file', action='store', type=str, default='LHC.log',
                        help='Input .log file from Herwig7 analysis')
    parser.add_argument('-o', '--output-dir', action='store', type=str, default='./', 
                        help='Output directory for generated ROOT file(s)')
    parser.add_argument('--tree-output-fname', default="AnalysisResults.root", type=str,
                        help="Filename for the (unscaled) generated particle ROOT TTree")
    parser.add_argument('--no-match-level', help="Save simulation for only one level with " + \
                        "no matching. Options: 'p', 'h', 'ch'", default=None, type=str)
    parser.add_argument('--no-scale', help="Turn off rescaling all histograms by cross section / N",
                        action='store_true', default=False)
    parser.add_argument('-c', '--config-file', action='store', type=str, default='config/angularity.yaml',
                        help="Path of config file for observable configurations")
    args = parser.parse_args()

    if args.no_match_level not in [None, 'p', 'h', 'ch']:
        print("ERROR: Unrecognized type %s. Please use 'p', 'h', or 'ch'" % args.type_only)
        exit(1)

    # If invalid configFile is given, exit
    if not os.path.exists(args.config_file):
        print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    process = herwig_parton_hadron(config_file=args.config_file, output_dir=args.output_dir, args=args)
    process.herwig_parton_hadron(args)
