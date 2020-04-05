#!/usr/bin/env python

from __future__ import print_function

import os
import pyhepmc_ng

import ROOT
ROOT.gROOT.SetBatch(True)

#---------------------------------------------------------------
def accept_particle_pythia(part, pdg):

  return accept_particle_status(part, pdg, status = [1])
  
#---------------------------------------------------------------
def accept_particle_herwig(part, pdg):

  return accept_particle_status(part, pdg, status = [1])
  
#---------------------------------------------------------------
def accept_particle_jewel(part, pdg):

  return accept_particle_status(part, pdg, status = [1])
  
#---------------------------------------------------------------
def accept_particle_martini(part, pdg):
  '''
  Status codes:
    parton: 1
    hadron: 201
    recoil: 301
    negative particles: 401
  '''
  
  return accept_particle_status(part, pdg, status = [201])
  
#---------------------------------------------------------------
def accept_particle_hybrid(part, pdg):
  '''
  Status codes:
    hadron: 1
    wake, positive: 6
    wake, negative: 7
  '''
  
  return accept_particle_status(part, pdg, status = [1, 6, 7])
  
#---------------------------------------------------------------
def accept_particle_status(part, pdg, status = [1]):
  
  # Check status
  #print('status: {}'.format(part.status))
  if part.status not in status:
    return False
  
  # Check that particle does not have any daughter vertex
  #print(part.end_vertex)
  if part.end_vertex:
    return False
  
  # Check PID for charged particles
  #print('pid: {} = {}'.format(part.pid, pdg.GetParticle(part.pid).GetName()))
  if pdg.GetParticle(part.pid).Charge() == 0:
    return False
    
  return True

#---------------------------------------------------------------
def accept_particle_jetscape(part, pdg):
  
  # The final-state hadrons are stored as outgoing particles in a disjoint vertex with t = 100
  parent_vertex = part.production_vertex
  vertex_time = parent_vertex.position.t
  if abs(vertex_time - 100) > 1e-3:
    return False
  
  # Check PID for charged particles, and reject neutrinos
  #print('pid: {} = {}'.format(part.pid, pdg.GetParticle(part.pid).GetName()))
  if pdg.GetParticle(part.pid).Charge() == 0:
    return False

  return True
