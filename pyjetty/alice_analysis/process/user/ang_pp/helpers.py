#!/user/bin/env python3

'''
Helper functions for angularity analysis.

Ezra Lesser (elesser@berkeley.edu)
'''

import numpy as np


# Return \Delta{R} between two fastjet.PsuedoJet objects  
def deltaR(pjet1, pjet2):
  return np.sqrt( (pjet1.eta() - pjet2.eta())**2 + (pjet1.phi() - pjet2.phi())**2 )

# Return jet angularity for fastjet.PseudoJet object
def lambda_beta_kappa(jet, jetR, beta, kappa):
  return sum( [ (constit.pt() / jet.pt())**kappa * (deltaR(jet, constit) / jetR)**beta 
                for constit in jet.constituents() ] )

# Helper function for finding the correct jet pT bin
def pT_bin(self, jet_pT, pTbins):
  for i, pTmin in list(enumerate(pTbins))[0:-1]:
    pTmax = pTbins[i+1]
    if pTmin <= jet_pT < pTmax:
      return (pTmin, pTmax)
  return (None, None)
