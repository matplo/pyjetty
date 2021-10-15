#!/usr/bin/env python

from __future__ import print_function

import os
import tqdm

import ROOT
ROOT.gROOT.SetBatch(True)

import select_particles

from pyjetty.alice_analysis.process.base import common_base

################################################################
class HepMC2antupleBase(common_base.CommonBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input = '', output = '', as_data = False, hepmc = 2, nev = 0, gen = 'pythia', no_progress_bar = False, include_parton = False, **kwargs):
    super(HepMC2antupleBase, self).__init__(**kwargs)
    self.input = input
    self.output = output
    self.as_data = as_data
    self.hepmc = hepmc
    self.nev = nev
    self.gen = gen
    self.no_progress_bar = no_progress_bar
    self.include_parton = include_parton

  #---------------------------------------------------------------
  def init(self):

    self.outf = ROOT.TFile(self.output, 'recreate')
    self.outf.cd()
    self.tdf = ROOT.TDirectoryFile('PWGHF_TreeCreator', 'PWGHF_TreeCreator')
    self.tdf.cd()
    if self.as_data:
      self.t_p = ROOT.TNtuple('tree_Particle', 'tree_Particle', 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID')
    else:
      self.t_p = ROOT.TNtuple('tree_Particle_gen', 'tree_Particle_gen', 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID')
      if self.include_parton:
        self.t_pp = ROOT.TNtuple('tree_Particle_gen_parton', 'tree_Particle_gen_parton', 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID')
    self.t_e = ROOT.TNtuple('tree_event_char', 'tree_event_char', 'run_number:ev_id:z_vtx_reco:is_ev_rej')

    # run number will be a double - file size in MB
    self.run_number = os.path.getsize(self.input) / 1.e6
    self.ev_id = 0

    # unfortunately pyhepmc_ng does not provide the table
    # pdt = pyhepmc_ng.ParticleDataTable()
    # use ROOT instead
    self.pdg = ROOT.TDatabasePDG()
    self.particles_accepted = set([])
    if self.include_parton:
      self.partons_accepted = set([])
    
    if not self.no_progress_bar:
      if self.nev > 0:
        self.pbar = tqdm.tqdm(range(self.nev))
      else:
        self.pbar = tqdm.tqdm()
  
  #---------------------------------------------------------------
  def accept_particle(self, part, status, end_vertex, pid, pdg, gen, parton=False):

    if gen == 'pythia':
      return select_particles.accept_particle_pythia(part, status, end_vertex, pid, pdg, parton)
    elif gen == 'herwig':
      return select_particles.accept_particle_herwig(part, status, end_vertex, pid, pdg, parton)
    elif gen == 'jewel':
      return select_particles.accept_particle_jewel(part, status, end_vertex, pid, pdg, parton)
    elif gen == 'jewel_charged':
      return select_particles.accept_particle_jewel(part, status, end_vertex, pid, pdg, parton, select_charged=True)
    elif gen == 'jetscape':
      return select_particles.accept_particle_jetscape(part, pdg, parton)
    elif gen == 'martini':
      return select_particles.accept_particle_martini(part, status, end_vertex, pid, pdg, parton)
    elif gen == 'hybrid':
      return select_particles.accept_particle_hybrid(part, status, end_vertex, pid, pdg, parton)

    sys.exit('Generator type unknown: {}'.format(gen))

  #---------------------------------------------------------------
  def increment_event(self):

    self.ev_id = self.ev_id + 1
    if not self.no_progress_bar:
      self.pbar.update()
    else:
      if self.ev_id % 100 == 0:
        print('event {}'.format(self.ev_id))

  #---------------------------------------------------------------
  def finish(self):
  
    self.print_particles()
    self.outf.Write()
    self.outf.Close()
  
  #---------------------------------------------------------------
  def print_particles(self):
  
    # Print final list of particles that we accepted
    print('particles included: {}'.format(self.particles_accepted))
    reference_particles = ['Omega+', 'Xi-', 'e+', 'Sigma-', 'mu-', 'Omega-', 'antiproton', 'proton', 'mu+', 'Sigma+', 'Sigma-_bar', 'Sigma+_bar', 'K-', 'pi+', 'K+', 'pi-', 'e-', 'Xi-_bar', 'antineutron', 'neutron', 'K_L0', 'gamma']
    
    # Check that we are not missing any particles
    for particle in reference_particles:
      if particle not in self.particles_accepted:
        print('WARNING: Missing particles: {} not found in your accepted particles!'.format(particle))
        
    # Check that we do not have any extra particles
    for particle in self.particles_accepted:
      if particle not in reference_particles:
        print('WARNING: Extra particles: {} was found in your accepted particles!'.format(particle))
