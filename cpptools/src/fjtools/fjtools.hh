#ifndef __PYJETTY_FJTOOLS_HH
#define __PYJETTY_FJTOOLS_HH

class TF1;
class TH2F;
#include <THn.h>
#include <string>
#include <vector>
#include <fastjet/PseudoJet.hh>

#ifdef USE_ROOUNFOLD
#include <RooUnfoldResponse.h>
#endif

namespace PyJettyFJTools
{
	class BoltzmannBackground
	{
	public:
		BoltzmannBackground();
		BoltzmannBackground(double mean_pt, double min_pt = 0.15, double max_pt = 5.);
		void reset(double mean_pt, double min_pt, double max_pt);

		~BoltzmannBackground();

		double eval(double x);
		double integral();
		std::string get_formula();
		double get_mean_pt() {return fmean_pt;}
		double get_max_pt_subtract() {return fmax_pt;}
		double get_constant();

		std::vector<fastjet::PseudoJet> get_parts() {return fparts;}

		// getNrandom particles w/ indexoffset ...
		std::vector<fastjet::PseudoJet> generate(int nparts, double max_eta, int offset = 0);
		// subtract particles (according to the probability... - fixed to 1 in maxpt range)
		std::vector<fastjet::PseudoJet> subtract_recalc_from_vector(const std::vector<fastjet::PseudoJet> &v);
		std::vector<fastjet::PseudoJet> subtract(const std::vector<fastjet::PseudoJet> &v, double mean_pt = -1., int n = -1);

		//void set_mean_pt(const double &v) {mean_pt = v;}
		//void set_max_pt_subtract(const double &v) {max_pt_subtract = v;}

		std::string description();
	private:
		TF1 *funbg; //!
		double fmean_pt;
		double fmin_pt;
		double fmax_pt;
		std::vector<fastjet::PseudoJet> fparts;
	};

	// fraction of pT of jet j1 contained in j0 - constit by constit
	double matched_pt(const fastjet::PseudoJet &j0, const fastjet::PseudoJet &j1);

	// return indices of jets matched to j jet - using rapidity to calculate deltaR
	std::vector<int> matched_Ry(const fastjet::PseudoJet &j, const std::vector<fastjet::PseudoJet> &v, double Rmatch);
	// return indices of jets matched to j jet - using pseudorapidity to calculate deltaR
	std::vector<int> matched_Reta(const fastjet::PseudoJet &j, const std::vector<fastjet::PseudoJet> &v, double Rmatch);

    // Rebin 2D histogram h with name hname using axes given by x_bins and y_bins
    TH2F* rebin_th2(TH2F* h, std::string hname, int n_x_bins, double* x_bins,
                    int n_y_bins, double* y_bins, bool move_y_underflow = false);

	// Rebin N-dimensional THn to a new histogram with name name_thn_rebinned using provided axes
	// WARNING: currently requires n_dim = 4
	THnF* rebin_thn(std::string response_file_name,
                    THnF* thn,
                    const std::string & name_thn_rebinned,
                    const std::string & name_roounfold,
					const unsigned int & n_dim,
                    const int* axes_n_bins,
					double** axes_bin_arrays,
                    const std::string label="",
                    const double & prior_variation_parameter=0.,
                    int prior_option=1,
					bool move_underflow=false,
                    bool do_roounfoldresponse=true);

	// Create empty THn using provided axes
	THnF* create_empty_thn(const char* name, const unsigned int & n_dim, const char** axes_titles,
						   const int* axes_n_bins, double** axes_bin_arrays);

	// Fill empty thn_rebinned with data from thn
	void fill_rebinned_thn(std::string response_file_name, THnF* thn,
                           THnF* thn_rebinned, const unsigned int & n_dim,
                           bool do_roounfoldresponse=true,
                           RooUnfoldResponse* roounfold_response=nullptr,
						   const double & prior_variation_parameter=0.,
                           int prior_option=1,
                           bool move_underflow=false);

    // Set scaling of prior
    double prior_scale_factor_obs(double obs_true, double content,
                                  double prior_variation_parameter, int option);

};

#endif
