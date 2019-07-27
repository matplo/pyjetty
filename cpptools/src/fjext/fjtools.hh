#ifndef PYJETTY_PYJETTY_HH
#define PYJETTY_PYJETTY_HH

#include <fastjet/PseudoJet.hh>
#include <vector>

namespace PyJetty
{
	std::vector<fastjet::PseudoJet> vectorize_pt_eta_phi(double *pt, int npt, double *eta, int neta, double *phi, int nphi);
}

#endif
