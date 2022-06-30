#ifndef __PYJETTY_THERMALIZER__HH
#define __PYJETTY_THERMALIZER__HH

#include <TObject.h>
#include <TH2F.h>
#include <TRandom3.h>
#include <TF1.h>

#include <fastjet/PseudoJet.hh>

namespace RUtilExt
{
	class Thermalizer : public TObject
	{
	public:
		Thermalizer();
		Thermalizer(Double_t meanpt, Double_t max_multiplicity, Double_t max_delta_R = -1, Double_t maxabseta = -1);
		virtual ~Thermalizer();

		std::vector<fastjet::PseudoJet> thermalize(Double_t pt, Double_t eta, Double_t phi, Double_t mass = 0);

		std::vector<fastjet::PseudoJet> thermalize(const fastjet::PseudoJet &p)
		{
			return thermalize(p.perp(), p.eta(), p.phi(), p.m());
		}

		std::vector<fastjet::PseudoJet> thermalize_vector(const std::vector<fastjet::PseudoJet> &parts);

	private:
		Double_t fMeanPt;
		Int_t 	 fMaxMultiplicity;
		Double_t fMaxDeltaR;
		Double_t fMaxAbsEta;
		TRandom3 fRandom;
		TF1      *fBoltzmann; //!
		TF1      *fExpo; //!

	ClassDef(Thermalizer, 1)
	};
};
#endif
