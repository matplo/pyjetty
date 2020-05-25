#!/user/bin/env python3

'''
Helper functions for angularity analysis.

Ezra Lesser (elesser@berkeley.edu)
'''

import numpy as np
from math import pi

''' # Not needed: use instead pjet1.delta_R(pjet2)
# Return \Delta{R} between two fastjet.PsuedoJet objects  
def deltaR(pjet1, pjet2):

  # Check that there's not a problem with +/- pi in phi part
  phi1 = pjet1.phi()
  if phi1 - pjet2.phi() > pi:
    phi1 -= 2*pi
  elif pjet2.phi() - phi1 > pi:
    phi1 += 2*pi

  return np.sqrt( (pjet1.eta() - pjet2.eta())**2 + (phi1 - pjet2.phi())**2 )
'''

# Return jet angularity for fastjet.PseudoJet object
def lambda_beta_kappa(jet, jetR, beta, kappa):
  return sum( [ (constit.pt() / jet.pt())**kappa * (jet.delta_R(constit) / jetR)**beta 
                for constit in jet.constituents() ] )

# Helper function for finding the correct jet pT bin
def pT_bin(jet_pT, pTbins):
  for i, pTmin in list(enumerate(pTbins))[0:-1]:
    pTmax = pTbins[i+1]
    if pTmin <= jet_pT < pTmax:
      return (pTmin, pTmax)
  return (-1, -1)
