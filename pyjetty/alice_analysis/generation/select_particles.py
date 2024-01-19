#!/usr/bin/env python

from __future__ import print_function

import os
import pyhepmc

import ROOT
ROOT.gROOT.SetBatch(True)

#---------------------------------------------------------------
def accept_particle_pythia(part, status, end_vertex, pid, pdg, parton=False, select_charged=False):

  if parton:
    raise NotImplementedError('Parton tree not implemented yet for PYTHIA') 

  return accept_particle_status(part, status, end_vertex, pid, pdg, status_accepted = [1], select_charged=select_charged)
  
#---------------------------------------------------------------
def accept_particle_herwig(part, status, end_vertex, pid, pdg, parton=False):

  if parton:
    raise NotImplementedError('Parton tree implemented but currently not working for Herwig') 
    return accept_particle_status(part, status, end_vertex, pid, pdg, parton, status_accepted = [11])

  return accept_particle_status(part, status, end_vertex, pid, pdg, status_accepted = [1])
  
#---------------------------------------------------------------
def accept_particle_jewel(part, status, end_vertex, pid, pdg, parton=False, select_charged=False):

  if parton:
    raise NotImplementedError('Parton tree not implemented yet for JEWEL')

  return accept_particle_status(part, status, end_vertex, pid, pdg, status_accepted = [1,3], 
                                select_charged=select_charged, charged_exception=[3], pt_exception=1.e-5)

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
def accept_particle_status(part, status, end_vertex, pid, pdg, parton=False, status_accepted = [1], select_charged=True, charged_exception=[], pt_exception=0.):
  
  # Check status
  #print('status: {} in {}?'.format(status, status_accepted))
  if status not in status_accepted:
    return False
  
  # Check that particle does not have any daughter vertex
  #if not parton and end_vertex:
  #  print(part, status, pid)
  if end_vertex:
    return False
  
  # Check PID for charged particles
  if select_charged:
    if pdg.GetParticle(pid):
      #if parton:
      #  print(part, status)
      #  print('pid: {} = {}'.format(part.pid, pdg.GetParticle(part.pid).GetName()))
      if not parton and pdg.GetParticle(pid).Charge() == 0:
        if status not in charged_exception and part.momentum().perp() > pt_exception: # Exceptions for JEWEL in order to keep recoils and dummies
          return False
    else:
      return False
    
  return True
