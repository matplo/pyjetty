#include "pyjetty.hh"
#include <iostream>
#include <fastjet/PseudoJet.hh>

namespace PyJetty
{
	double sum(double *seq, int n)
	{
		double _sum = 0;
		for (unsigned int i = 0; i < n; i++)
		{
			_sum += seq[i];
		}
		return _sum;
	}

	std::vector<fastjet::PseudoJet> vectorize_pt_eta_phi(double *pt, int npt, double *eta, int neta, double *phi, int nphi)
	{
		std::vector<fastjet::PseudoJet> v;
		if (npt != neta || npt != nphi) 
		{
			std::cerr << "[error] vectorize_pt_eta_phi : incompatible array sizes" << std::endl;
			return v;
		}
		for (unsigned int i = 0; i < npt; i++)
		{
		    double px = pt[i] * cos(phi[i]);
		    double py = pt[i] * sin(phi[i]);
		    double pz = pt[i] * sinh(eta[i]);
		    double e  = sqrt(px*px + py*py + pz*pz);
			fastjet::PseudoJet psj(px, py, pz, e);
			psj.set_user_index(i);
			v.push_back(psj);
		}
		return v;
	}

	TestClass::TestClass() {;}
	TestClass::~TestClass() {;}
};