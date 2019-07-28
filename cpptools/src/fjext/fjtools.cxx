#include "fjtools.hh"
#include <iostream>
#include <fastjet/PseudoJet.hh>

namespace PyJetty
{
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

	std::vector<fastjet::PseudoJet> vectorize_px_py_pz(double *px, int npx, double *py, int npy, double *pz, int npz)
	{
		std::vector<fastjet::PseudoJet> v;
		if (npx != npy || npx != npz) 
		{
			std::cerr << "[error] vectorize_px_py_pz : incompatible array sizes" << std::endl;
			return v;
		}
		for (unsigned int i = 0; i < npx; i++)
		{
		    double e  = sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
			fastjet::PseudoJet psj(px[i], py[i], pz[i], e);
			psj.set_user_index(i);
			v.push_back(psj);
		}
		return v;
	}

	std::vector<fastjet::PseudoJet> vectorize_px_py_pz_e(double *px, int npx, double *py, int npy, double *pz, int npz, double *e, int ne)
	{
		std::vector<fastjet::PseudoJet> v;
		if (npx != npy || npx != npz || npx != ne) 
		{
			std::cerr << "[error] vectorize_px_py_pz_e : incompatible array sizes" << std::endl;
			return v;
		}
		for (unsigned int i = 0; i < npx; i++)
		{
			fastjet::PseudoJet psj(px[i], py[i], pz[i], e[i]);
			psj.set_user_index(i);
			v.push_back(psj);
		}
		return v;
	}

	std::vector<fastjet::PseudoJet> vectorize_px_py_pz_m(double *px, int npx, double *py, int npy, double *pz, int npz, double *m, int nm)
	{
		std::vector<fastjet::PseudoJet> v;
		if (npx != npy || npx != npz || npx != nm) 
		{
			std::cerr << "[error] vectorize_px_py_pz_m : incompatible array sizes" << std::endl;
			return v;
		}
		for (unsigned int i = 0; i < npx; i++)
		{
		    double e  = sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i] * m[i]*m[i]);
			fastjet::PseudoJet psj(px[i], py[i], pz[i], e);
			psj.set_user_index(i);
			v.push_back(psj);
		}
		return v;
	}
};
