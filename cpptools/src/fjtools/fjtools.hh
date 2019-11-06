#ifndef __PYJETTY_FJTOOLS_HH
#define __PYJETTY_FJTOOLS_HH

class TF1;
#include <string>
#include <vector>
#include <fastjet/PseudoJet.hh>

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
		std::vector<fastjet::PseudoJet> subtract(const std::vector<fastjet::PseudoJet> &v);

		//void set_mean_pt(const double &v) {mean_pt = v;}
		//void set_max_pt_subtract(const double &v) {max_pt_subtract = v;}

	private:
		TF1 *funbg; //!
		double fmean_pt;
		double fmin_pt;
		double fmax_pt;
		std::vector<fastjet::PseudoJet> fparts;
	};
};

#endif