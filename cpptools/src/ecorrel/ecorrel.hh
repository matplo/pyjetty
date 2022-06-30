#ifndef __PYJETTY_ECORREL_HH
#define __PYJETTY_ECORREL_HH

#include <vector>
#include <fastjet/PseudoJet.hh>

namespace EnergyCorrelators
{
	class CorrelatorsContainer
	{
		public:
			CorrelatorsContainer();
			virtual ~CorrelatorsContainer();
			void addwr(const double &w, const double &r);
			std::vector<double> *weights();
			std::vector<double> *rs();
			std::vector<double> *rxw();

		private:
			std::vector<double> fr;
			std::vector<double> fw;
			std::vector<double> frxw;
	};

	CorrelatorsContainer EEC(const std::vector<fastjet::PseudoJet> &parts, const double &scale);
	CorrelatorsContainer E3C(const std::vector<fastjet::PseudoJet> &parts, const double &scale);
	CorrelatorsContainer E4C(const std::vector<fastjet::PseudoJet> &parts, const double &scale);

	std::vector<fastjet::PseudoJet> constituents_as_vector(const fastjet::PseudoJet &jet);
	
	//	class CorrelatorBuilder
	//	{
	//	public:
	//		CorrelatorBuilder();
	//		CorrelatorBuilder(std::vector<fastjet::PseudoJet> &parts, double scale);
	//		~CorrelatorBuilder();
	//
	//		std::vector<double> EEC(double scale = 0.);
	//		std::vector<double> E3C(double scale = 0.);
	//		std::vector<double> E4C(double scale = 0.);
	//
	//		void setScale(double scale);
	//		double getScale();
	//
	//		std::string description();
	//	private:
	//		std::vector<fastjet::PseudoJet> fparts;
	//		double scale;
	//	};

};

#endif
