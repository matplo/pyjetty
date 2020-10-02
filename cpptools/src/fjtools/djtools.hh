#ifndef __PYJETTY_DJETFJTOOLS_HH
#define __PYJETTY_DJETFJTOOLS_HH

#include <vector>
#include <fastjet/PseudoJet.hh>

namespace DJetPyJettyFJTools
{
	// note - this is a copy (for allow experimental code) from HEPPY/cpptools/src/fjext/fjtools.cxx
	std::vector<fastjet::PseudoJet> vectorize_pt_eta_phi(double *pt, int npt, double *eta, int neta, double *phi, int nphi, int user_index_offset = 0);
	std::vector<fastjet::PseudoJet> vectorize_pt_eta_phi_m(double *pt, int npt, double *eta, int neta, double *phi, int nphi, double *m, int nm, int user_index_offset = 0);

	class DJetMatchMaker
	{
	public:
		DJetMatchMaker();
		void set_Ds_pt_eta_phi_m(double *pt, int npt, double *eta, int neta, double *phi, int nphi, double *m, int nm, int user_index_offset = 10000);
		void set_daughters0_pt_eta_phi_m(double *pt, int npt, double *eta, int neta, double *phi, int nphi, double *m, int nm, int user_index_offset = 20000);
		void set_daughters1_pt_eta_phi_m(double *pt, int npt, double *eta, int neta, double *phi, int nphi, double *m, int nm, int user_index_offset = 30000);
		void set_ch_pt_eta_phi_m(double *pt, int npt, double *eta, int neta, double *phi, int nphi, double *m, int nm, int user_index_offset = 0);

		void set_Ds_pt_eta_phi(double *pt, int npt, double *eta, int neta, double *phi, int nphi, int user_index_offset = 10000);
		void set_daughters0_pt_eta_phi(double *pt, int npt, double *eta, int neta, double *phi, int nphi, int user_index_offset = 20000);
		void set_daughters1_pt_eta_phi(double *pt, int npt, double *eta, int neta, double *phi, int nphi, int user_index_offset = 30000);
		void set_ch_pt_eta_phi(double *pt, int npt, double *eta, int neta, double *phi, int nphi, int user_index_offset = 0);

		std::vector<fastjet::PseudoJet> match(double r, int n = -1);
		bool is_a_djet(const fastjet::PseudoJet &j);
		std::vector<fastjet::PseudoJet> filter_D0_jets(const std::vector<fastjet::PseudoJet> &jets);
		std::vector<fastjet::PseudoJet> get_Dcand_in_jet(const fastjet::PseudoJet &j);

		~DJetMatchMaker();

		// fastjet guts here
		std::vector<fastjet::PseudoJet> ch;
		std::vector<fastjet::PseudoJet> Ds;
		std::vector<fastjet::PseudoJet> daughters0;
		std::vector<fastjet::PseudoJet> daughters1;

		int ch_user_index_offset;
		int Dcand_user_index_offset;
		int daughter0_user_index_offset;
		int daughter1_user_index_offset;
	};
}

#endif