#!/usr/bin/env python

from __future__ import print_function

import os
import pyhepmc_ng

import ROOT
ROOT.gROOT.SetBatch(True)

#---------------------------------------------------------------
def accept_particle_pythia(part, status, end_vertex, pid, pdg, parton=False):

  if parton:
    raise NotImplementedError('Parton tree not implemented yet for PYTHIA') 

  return accept_particle_status(part, status, end_vertex, pid, pdg, status_accepted = [1])
  
#---------------------------------------------------------------
def accept_particle_herwig(part, status, end_vertex, pid, pdg, parton=False):

  if parton:
    raise NotImplementedError('Parton tree implemented but currently not working for Herwig') 
    return accept_particle_status(part, status, end_vertex, pid, pdg, parton, status_accepted = [11])

  return accept_particle_status(part, status, end_vertex, pid, pdg, status_accepted = [1])
  
#---------------------------------------------------------------
def accept_particle_jewel(part, status, end_vertex, pid, pdg, parton=False):

  if parton:
    raise NotImplementedError('Parton tree not implemented yet for JEWEL') 

  return accept_particle_status(part, status, end_vertex, pid, pdg, status_accepted = [1], select_charged=False)

#---------------------------------------------------------------
def accept_particle_jetscape(part, pdg, parton=False):

  if parton:
    raise NotImplementedError('Parton tree not implemented yet for JETSCAPE') 
  
  # The final-state hadrons are stored as outgoing particles in a disjoint vertex with t = 100
  parent_vertex = part.production_vertex
  vertex_time = parent_vertex.position.t
  if abs(vertex_time - 100) > 1e-3:
    return False
  
  # Check PID for charged particles, and reject neutrinos
  #print('pid: {} = {}'.format(part.pid, pdg.GetParticle(part.pid).GetName()))
  if pdg.GetParticle(part.pid):
    if pdg.GetParticle(part.pid).Charge() == 0:
      return False
  else:
    return False

  return True

#---------------------------------------------------------------
def accept_particle_martini(part, status, end_vertex, pid, pdg, parton=False):
  '''
  Status codes:
    parton: 1
    hadron: 201
    recoil: 301
    negative particles: 401
  '''

  if parton:
    raise NotImplementedError('Parton tree not implemented yet for MARTINI') 

  return accept_particle_status(part, status, end_vertex, pid, pdg, status_accepted = [201])
  
#---------------------------------------------------------------
def accept_particle_hybrid(part, pdg, parton):
  '''
  Status codes:
    hadron: 1
    wake, positive: 6
    wake, negative: 7
  '''

  if parton:
    raise NotImplementedError('Parton tree not implemented yet for hybrid') 
  
  return accept_particle_status(part, status, end_vertex, pid, pdg, status_accepted = [1, 6, 7])
  
#---------------------------------------------------------------
def accept_particle_status(part, status, end_vertex, pid, pdg, parton=False, status_accepted = [1], select_charged=True):
  
  # Check status
  #print('status: {} in {}'.format(status, status_accpted))
  if status not in status_accepted:
    return False
  
  # Check that particle does not have any daughter vertex
  #if not parton and end_vertex:
  #  print(part, status, pid)
  if end_vertex:
    return False
  
  # Check PID for charged particles
  if select_charged:
    if pdg.GetParticle(part.pid):
      #if parton:
      #  print(part, status)
      #  print('pid: {} = {}'.format(part.pid, pdg.GetParticle(part.pid).GetName()))
      if not parton and pdg.GetParticle(part.pid).Charge() == 0:
        return False
      else:
        return False
    
  return True
