#ifndef PYJETTY_PYJETTY_HH
#define PYJETTY_PYJETTY_HH

#include <fastjet/PseudoJet.hh>
#include <vector>

namespace PyJetty
{
	double sum(double *seq, int n);
	std::vector<fastjet::PseudoJet> vectorize_pt_eta_phi(double *pt, int npt, double *eta, int neta, double *phi, int nphi);
	
	class TestClass
	{
	public:
		TestClass();
		~TestClass();
		int status;
	};
}

#endif