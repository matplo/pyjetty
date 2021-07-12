#!/usr/bin/env python3

'''
For generating particle TTrees at final-state parton & hadron level
using the PYTHIA 8 Monte Carlo simulation
Ezra Lesser (elesser@berkeley.edu)
'''

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import ROOT

import copy
import argparse
import os

from pyjetty.mputils import *

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

from pyjetty.alice_analysis.process.base import common_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)
# Automatically set Sumw2 when creating new histograms
#ROOT.TH1.SetDefaultSumw2()


################################################################
class PartonHadronBase(common_base.CommonBase):

    #---------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------
    def __init__(self, output_dir='', tree_output_fname='', no_ALICE_eta_cut=False,
                 args=None, **kwargs):

        super(PartonHadronBase, self).__init__(**kwargs)

        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.tree_output_fname = tree_output_fname

        self.no_ALICE_eta_cut = no_ALICE_eta_cut
        self.ALICE_max_eta = 0.9

        self.init_tree_writers()


    #---------------------------------------------------------------
    # Initialize tree writers
    #---------------------------------------------------------------
    def init_tree_writers(self):

        # MPI off
        name = "tree_Particle_gen_h_MPIoff"
        t = ROOT.TTree(name, name)
        setattr(self, name, t)
        tw = RTreeWriter(tree=t)
        setattr(self, name+"_writer", tw)

        name = "tree_Particle_gen_p_MPIoff"
        t = ROOT.TTree(name, name)
        setattr(self, name, t)
        tw = RTreeWriter(tree=t)
        setattr(self, name+"_writer", tw)

        # MPI on
        name = "tree_Particle_gen_h_MPIon"
        t = ROOT.TTree(name, name)
        setattr(self, name, t)
        tw = RTreeWriter(tree=t)
        setattr(self, name+"_writer", tw)

        name = "tree_Particle_gen_p_MPIon"
        t = ROOT.TTree(name, name)
        setattr(self, name, t)
        tw = RTreeWriter(tree=t)
        setattr(self, name+"_writer", tw)


    #---------------------------------------------------------------
    # Fill all trees branches with relevant particle information
    #---------------------------------------------------------------
    def fill_branches(self, parts_p, parts_h, h_is_charged, iev, run_number, MPIon=False):

        t_p_name = "tree_Particle_gen_p_MPI%s" % ("on" if MPIon else "off")
        t_h_name = "tree_Particle_gen_h_MPI%s" % ("on" if MPIon else "off")
        tw_p = getattr(self, t_p_name+"_writer")
        tw_h = getattr(self, t_h_name+"_writer")

        for particle in parts_p:
            if not self.no_ALICE_eta_cut and abs(particle.eta()) > self.ALICE_max_eta:
                continue
            tw_p.fill_branch('run_number', run_number)
            tw_p.fill_branch('ev_id', iev)
            tw_p.fill_branch('ParticleE', particle.E())
            tw_p.fill_branch('ParticlePx', particle.px())
            tw_p.fill_branch('ParticlePy', particle.py())
            tw_p.fill_branch('ParticlePz', particle.pz())
            # Seems user_info is not working for some reason... 
            #tw_p.fill_branch('is_gluon', particle.user_info().isGluon())
            #tw_p.fill_branch('is_quark', particle.user_info().isQuark())

        for particle, is_charged in zip(parts_h, h_is_charged):
            if not self.no_ALICE_eta_cut and abs(particle.eta()) > self.ALICE_max_eta:
                continue
            tw_h.fill_branch('run_number', run_number)
            tw_h.fill_branch('ev_id', iev)
            tw_h.fill_branch('ParticleE', particle.E())
            tw_h.fill_branch('ParticlePx', particle.px())
            tw_h.fill_branch('ParticlePy', particle.py())
            tw_h.fill_branch('ParticlePz', particle.pz())
            # Seems user_info is not working for some reason... 
            #tw_h.fill_branch('is_charged', particle.user_info().isCharged())
            tw_h.fill_branch('is_charged', is_charged)


    #---------------------------------------------------------------
    # Fill trees after all branches have been filled
    #---------------------------------------------------------------
    def fill_write_trees(self):

        for MPI in ["off", "on"]:
            for level in ["p", "h"]:
                getattr(self, "tree_Particle_gen_%s_MPI%s_writer" % (level, MPI)).fill_tree()
                getattr(self, "tree_Particle_gen_%s_MPI%s" % (level, MPI)).Write()

    #---------------------------------------------------------------
    # Fill trees after all branches have been filled
    #---------------------------------------------------------------
    def save_xsec_N(self, xsec_MPIoff, Nev_MPIoff, xsec_MPIon, Nev_MPIon):

        hxsec = ROOT.TH1F('hxsec_MPIoff', 'hxsec_MPIoff', 1, -0.5, 0.5)
        hxsec.SetBinContent(1, xsec_MPIoff)
        hxsec.SetBinError(1, 0)
        hxsec.Write()

        hxsec = ROOT.TH1F('hxsec_MPIon', 'hxsec_MPIon', 1, -0.5, 0.5)
        hxsec.SetBinContent(1, xsec_MPIon)
        hxsec.SetBinError(1, 0)
        hxsec.Write()

        hNev = ROOT.TH1F('hNev_MPIoff', 'hNev_MPIoff', 1, -0.5, 0.5)
        hNev.SetBinContent(1, Nev_MPIoff)
        hNev.SetBinError(1, 0)
        hNev.Write()

        hNev = ROOT.TH1F('hNev_MPIon', 'hNev_MPIon', 1, -0.5, 0.5)
        hNev.SetBinContent(1, Nev_MPIon)
        hNev.SetBinError(1, 0)
        hNev.Write()
