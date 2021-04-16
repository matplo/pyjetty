#!/usr/bin/env python

from __future__ import print_function

import os
import argparse
import tqdm
import numpy as np
import pandas as pd

import ROOT
ROOT.gROOT.SetBatch(True)

import fastjet as fj
import fjext

from pyjetty.alice_analysis.process.base import common_base

################################################################
class Parquet2antuple(common_base.CommonBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input = '', output = '', nev = 0, no_progress_bar = False, **kwargs):
    super(Parquet2antuple, self).__init__(**kwargs)
    
    self.input = input
    self.output = output
    self.nev = nev
    self.no_progress_bar = no_progress_bar
    
    self.init()
    print(self)
   
  #---------------------------------------------------------------
  def init(self):

    self.outf = ROOT.TFile(self.output, 'recreate')
    self.outf.cd()
    self.tdf = ROOT.TDirectoryFile('PWGHF_TreeCreator', 'PWGHF_TreeCreator')
    self.tdf.cd()
    self.t_p = ROOT.TNtuple('tree_Particle_gen', 'tree_Particle_gen', 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID:status')
    self.t_e = ROOT.TNtuple('tree_event_char', 'tree_event_char', 'run_number:ev_id:z_vtx_reco:is_ev_rej:event_plane_angle')

    # run number will be a double - file size in MB
    self.run_number = os.path.getsize(self.input) / 1.e6
    self.ev_id = 0

    self.pdg = ROOT.TDatabasePDG()
    self.particles_accepted = set([])
  
  #---------------------------------------------------------------
  def main(self):
  
    # Read chunk of events into a dataframe
    # Fields: particle_ID, status, E, px, py, pz, event_plane_angle
    df_event_chunk = pd.read_parquet(self.input)
    
    self.n_event_max = df_event_chunk.shape[0]
    if not self.no_progress_bar:
      if self.nev > 0:
        self.pbar = tqdm.tqdm(range(self.nev))
      else:
        self.pbar = tqdm.tqdm(range(self.n_event_max))

    # Iterate through events
    self.analyze_event_chunk(df_event_chunk)
        
    self.finish()
    
  # ---------------------------------------------------------------
  # Analyze event chunk
  # ---------------------------------------------------------------
  def analyze_event_chunk(self, df_event_chunk):
    
    # Loop through events
    for i,event in df_event_chunk.iterrows():
  
        if self.no_progress_bar and i % 1000 == 0:
            print('event: {}'.format(i))
  
        self.fill_event(event)
        self.increment_event()
        if self.nev > 0 and self.ev_id > self.nev:
          break

  #---------------------------------------------------------------
  def fill_event(self, event):
  
    self.t_e.Fill(self.run_number, self.ev_id, 0, 0, event['event_plane_angle'])
        
    # Get arrays of particle quantities in the event
    #   We will use the status to subtract holes at analysis time
    px = event['px']
    py = event['py']
    pz = event['pz']
    e = event['E']
    pid = event['particle_ID']
    status = event['status']

    # Looop through and fill each particle to tree
    for i,pid_i in enumerate(pid):
    
      pt = np.sqrt(px[i]*px[i] + py[i]*py[i])
      eta = np.arcsinh(pz[i]/pt)
      phi = np.arctan2(py[i],px[i])
      if phi < 0:
        phi += 2*np.pi

      # (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
      if abs(pid_i) in [11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334]:
          self.particles_accepted.add(self.pdg.GetParticle(int(pid_i)).GetName())
          self.t_p.Fill(self.run_number, self.ev_id, pt, eta, phi, pid_i, status[i])
        
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
    reference_particles = ['Omega+', 'Xi-', 'e+', 'Sigma-', 'mu-', 'Omega-', 'antiproton', 'proton', 'mu+', 'Sigma+', 'Sigma-_bar', 'Sigma+_bar', 'K-', 'pi+', 'K+', 'pi-', 'e-', 'Xi-_bar']
    
    # Check that we are not missing any particles
    for particle in reference_particles:
      if particle not in self.particles_accepted:
        print('WARNING: Missing particles: {} not found in your accepted particles!'.format(particle))
        
    # Check that we do not have any extra particles
    for particle in self.particles_accepted:
      if particle not in reference_particles:
        print('WARNING: Extra particles: {} was found in your accepted particles!'.format(particle))
        
#---------------------------------------------------------------
if __name__ == '__main__':
  
  parser = argparse.ArgumentParser(description='hepmc to ALICE Ntuple format', prog=os.path.basename(__file__))
  parser.add_argument('-i', '--input', help='input file', default='', type=str, required=True)
  parser.add_argument('-o', '--output', help='output root file', default='', type=str, required=True)
  parser.add_argument('--nev', help='number of events', default=-1, type=int)
  parser.add_argument('--no-progress-bar', help='whether to print progress bar', action='store_true', default=False)
  args = parser.parse_args()
  
  converter = Parquet2antuple(input = args.input, output = args.output, nev = args.nev, no_progress_bar = args.no_progress_bar)
  converter.main()
