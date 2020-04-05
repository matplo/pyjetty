#!/usr/bin/env python

from __future__ import print_function

import os
import argparse
import tqdm

import pyhepmc_ng

import ROOT
ROOT.gROOT.SetBatch(True)

import select_particles

# jit improves execution time by 18% - tested with jetty pythia8 events
# BUT produces invalid root file
# from numba import jit
# @jit

#---------------------------------------------------------------
def fill_event(run_number, ev_id, event_hepmc, tw_e, tw_p, pdg, gen, particles_accepted):

  tw_e.Fill(run_number, ev_id, 0, 0)

  for part in event_hepmc.particles:
  
    if accept_particle(part, pdg, gen):
    
      particles_accepted.add(pdg.GetParticle(part.pid).GetName())
      tw_p.Fill(run_number, ev_id, part.momentum.pt(), part.momentum.eta(), part.momentum.phi(), part.pid)
 
#---------------------------------------------------------------
def accept_particle(part, pdg, gen):

  if gen == 'pythia':
    return select_particles.accept_particle_pythia(part, pdg)
  if gen == 'herwig':
    return select_particles.accept_particle_herwig(part, pdg)
  if gen == 'jewel':
    return select_particles.accept_particle_jewel(part, pdg)
  elif gen == 'jetscape':
    return select_particles.accept_particle_jetscape(part, pdg)
  elif gen == 'martini':
    return select_particles.accept_particle_martini(part, pdg)
  elif gen == 'hybrid':
    return select_particles.accept_particle_hybrid(part, pdg)
  else:
    sys.exit('Generator type unknown: {}'.format(gen))

#---------------------------------------------------------------
def main():
  parser = argparse.ArgumentParser(description='hepmc to ALICE Ntuple format', prog=os.path.basename(__file__))
  parser.add_argument('-i', '--input', help='input file', default='', type=str, required=True)
  parser.add_argument('-o', '--output', help='output root file', default='', type=str, required=True)
  parser.add_argument('--as-data', help='write as data - tree naming convention', action='store_true', default=False)
  parser.add_argument('--hepmc', help='what format 2 or 3', default=2, type=int)
  parser.add_argument('--nev', help='number of events', default=-1, type=int)
  parser.add_argument('-g', '--gen', help='generator type: pythia, herwig, jewel, jetscape, martini, hybrid', default='pythia', type=str, required=True)
  parser.add_argument('--no-progress-bar', help='whether to print progress bar', action='store_true', default=False)
  args = parser.parse_args()

  if args.hepmc == 3:
    input_hepmc = pyhepmc_ng.ReaderAscii(args.input)
  if args.hepmc == 2:
    input_hepmc = pyhepmc_ng.ReaderAsciiHepMC2(args.input)

  if input_hepmc.failed():
    print ("[error] unable to read from {}".format(args.input))
    sys.exit(1)

  outf = ROOT.TFile(args.output, 'recreate')
  outf.cd()
  tdf = ROOT.TDirectoryFile('PWGHF_TreeCreator', 'PWGHF_TreeCreator')
  tdf.cd()
  if args.as_data:
    t_p = ROOT.TNtuple('tree_Particle', 'tree_Particle', 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID')
  else:
    t_p = ROOT.TNtuple('tree_Particle_gen', 'tree_Particle_gen', 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID')
  t_e = ROOT.TNtuple('tree_event_char', 'tree_event_char', 'run_number:ev_id:z_vtx_reco:is_ev_rej')

  # run number will be a double - file size in MB
  run_number = os.path.getsize(args.input) / 1.e6
  ev_id = 0

  # unfortunately pyhepmc_ng does not provide the table
  # pdt = pyhepmc_ng.ParticleDataTable()
  # use ROOT instead
  pdg = ROOT.TDatabasePDG()

  event_hepmc = pyhepmc_ng.GenEvent()

  if not args.no_progress_bar:
    if args.nev > 0:
      pbar = tqdm.tqdm(range(args.nev))
    else:
      pbar = tqdm.tqdm()

  particles_accepted = set([])
  while not input_hepmc.failed():
    ev = input_hepmc.read_event(event_hepmc)
    if input_hepmc.failed():
      break

    fill_event(run_number, ev_id, event_hepmc, t_e, t_p, pdg, args.gen, particles_accepted)

    ev_id = ev_id + 1
    if not args.no_progress_bar:
      pbar.update()
    else:
      if ev_id % 100 == 0:
        print('event {}'.format(ev_id))
    if args.nev > 0 and ev_id > args.nev:
      break
      
  # Print final list of particles that we accepted
  print('particles included: {}'.format(particles_accepted))
  reference_particles = ['Omega+', 'Xi-', 'e+', 'Sigma-', 'mu-', 'Omega-', 'antiproton', 'proton', 'mu+', 'Sigma+', 'Sigma-_bar', 'Sigma+_bar', 'K-', 'pi+', 'K+', 'pi-', 'e-', 'Xi-_bar']
  
  # Check that we are not missing any particles
  for particle in reference_particles:
    if particle not in particles_accepted:
      print('WARNING: {} not found in your accepted particles!'.format(particle))
      
  # Check that we do not have any extra particles
  for particle in particles_accepted:
    if particle not in reference_particles:
      print('WARNING: {} was found in your accepted particles!'.format(particle))

  outf.Write()
  outf.Close()

#---------------------------------------------------------------
if __name__ == '__main__':
  main()
