#include "djtools.hh"

#include <TF1.h>
#include <TMath.h>
#include <TString.h>
#include <TRandom.h>
#include <iostream>

namespace DJetPyJettyFJTools
{
	DJetMatchMaker::DJetMatchMaker()
	{
		;
	}

	void DJetMatchMaker::
	set_daughters0_pt_eta_phi_m(double *pt, int npt, double *eta, int neta, double *phi, int nphi, double *m, int nm, int user_index_offset)
	{
		daughters0.clear();
		daughters0 = vectorize_pt_eta_phi_m(pt, npt, eta, neta, phi, nphi, m, nm, user_index_offset);
	}

	void DJetMatchMaker::
	set_daughters1_pt_eta_phi_m(double *pt, int npt, double *eta, int neta, double *phi, int nphi, double *m, int nm, int user_index_offset)
	{
		daughters1.clear();
		daughters1 = vectorize_pt_eta_phi_m(pt, npt, eta, neta, phi, nphi, m, nm, user_index_offset);
	}

	void DJetMatchMaker::
	set_Ds_pt_eta_phi_m(double *pt, int npt, double *eta, int neta, double *phi, int nphi, double *m, int nm, int user_index_offset)
	{
		Ds.clear();
		Ds = vectorize_pt_eta_phi_m(pt, npt, eta, neta, phi, nphi, m, nm, user_index_offset);		
	}

	void DJetMatchMaker::
	set_ch_pt_eta_phi_m(double *pt, int npt, double *eta, int neta, double *phi, int nphi, double *m, int nm, int user_index_offset)
	{
		ch.clear();
		ch = vectorize_pt_eta_phi_m(pt, npt, eta, neta, phi, nphi, m, nm, user_index_offset);		
	}

	void DJetMatchMaker::
	set_daughters0_pt_eta_phi(double *pt, int npt, double *eta, int neta, double *phi, int nphi, int user_index_offset)
	{
		daughters0.clear();
		daughters0 = vectorize_pt_eta_phi(pt, npt, eta, neta, phi, nphi, user_index_offset);
	}

	void DJetMatchMaker::
	set_daughters1_pt_eta_phi(double *pt, int npt, double *eta, int neta, double *phi, int nphi, int user_index_offset)
	{
		daughters1.clear();
		daughters1 = vectorize_pt_eta_phi(pt, npt, eta, neta, phi, nphi, user_index_offset);
	}

	void DJetMatchMaker::
	set_Ds_pt_eta_phi(double *pt, int npt, double *eta, int neta, double *phi, int nphi, int user_index_offset)
	{
		Ds.clear();
		Ds = vectorize_pt_eta_phi(pt, npt, eta, neta, phi, nphi, user_index_offset);		
	}

	void DJetMatchMaker::
	set_ch_pt_eta_phi(double *pt, int npt, double *eta, int neta, double *phi, int nphi, int user_index_offset)
	{
		ch.clear();
		ch = vectorize_pt_eta_phi(pt, npt, eta, neta, phi, nphi, user_index_offset);		
	}

	std::vector<fastjet::PseudoJet> DJetMatchMaker::match(double r, int n)
	{
		std::vector<fastjet::PseudoJet> v;
		unsigned int nmatch = 0;
		if (n < 0)
		{
			// run for all matches in the same call
			for (auto p : ch)
			{
				for (auto &d : daughters0 )
				{
					if (d.delta_R(p) < r)
					{
						p = p * 1.e-6;
						p.set_user_index(d.user_index());
						nmatch++;
					}
				}
				for (auto &d : daughters1 )
				{
					if (d.delta_R(p) < r)
					{
						p = p * 1.e-6;
						p.set_user_index(d.user_index());
						nmatch++;
					}
				}
				v.push_back(p);
			}
		}
		else
		{
			if (n < Ds.size())
			{
				for (auto p : ch)
				{
					if (daughters0[n].delta_R(p) < r)
						{
							p = p * 1.e-6;
							p.set_user_index(daughters0[n].user_index());
							nmatch++;
						}				
					if (daughters1[n].delta_R(p) < r)
						{
							p = p * 1.e-6;
							p.set_user_index(daughters1[n].user_index());
							nmatch++;
						}				
					v.push_back(p);
				}
			}
			else
			{
				std::cerr << "[error] DJetMatchMaker::match : asking to match N beyond D0cand array size" << std::endl;
			}
		}
		std::cout << "[i] nmatch for dR=" << r << " : " << nmatch << std::endl;
		return v;
	}

	DJetMatchMaker::~DJetMatchMaker()
	{
		;
	}

	// this is copy from HEPPY (for experimental purposes pasted here)
	// https://github.com/matplo/heppy/blob/4fa9e09c20e2fc08d6e54d8b7087f36ebc595309/cpptools/src/fjext/fjtools.cxx#L20
	std::vector<fastjet::PseudoJet> vectorize_pt_eta_phi(double *pt, int npt, double *eta, int neta, double *phi, int nphi, int user_index_offset)
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
			psj.set_user_index(i + user_index_offset);
			v.push_back(psj);
		}
		return v;
	}

	// this is copy from HEPPY (for experimental purposes pasted here)
	// https://github.com/matplo/heppy/blob/4fa9e09c20e2fc08d6e54d8b7087f36ebc595309/cpptools/src/fjext/fjtools.cxx#L41
	std::vector<fastjet::PseudoJet> vectorize_pt_eta_phi_m(double *pt, int npt, double *eta, int neta, double *phi, int nphi, double *m, int nm, int user_index_offset)
	{
		std::vector<fastjet::PseudoJet> v;
		if (npt != neta || npt != nphi || npt != nm) 
		{
			std::cerr << "[error] vectorize_pt_eta_phi : incompatible array sizes" << std::endl;
			return v;
		}
		for (unsigned int i = 0; i < npt; i++)
		{
		    double px = pt[i] * cos(phi[i]);
		    double py = pt[i] * sin(phi[i]);
		    double pz = pt[i] * sinh(eta[i]);
		    double e  = sqrt(px*px + py*py + pz*pz + m[i]*m[i]);
			fastjet::PseudoJet psj(px, py, pz, e);
			psj.set_user_index(i + user_index_offset);
			v.push_back(psj);
		}
		return v;
	}

}