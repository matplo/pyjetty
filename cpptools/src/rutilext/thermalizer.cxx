#include <thermalizer.hh>

#include <Math/Vector4D.h>
#include <TMath.h>

#include <iostream>
ClassImp(RUtilExt::Thermalizer)

namespace RUtilExt
{
	double boltzmann(double x[], double par[])
	{
		return 1. / par[0] * x[0] * TMath::Exp(-(2. / par[0]) * x[0]);
	}

	double gauss(double x[0], double par[])
	{
		return par[0] / (par[2] * TMath::Sqrt(2. * TMath::Pi())) * TMath::Exp(-1./2. * TMath::Power((x[0]-par[1]) / par[2], 2.));
	}

	Thermalizer::Thermalizer() 
		: TObject()
		, fMeanPt(0.7)
		, fMaxMultiplicity(-1)
		, fMaxDeltaR(-1)
		, fMaxAbsEta(-1)
		, fRandom(0)
		, fBoltzmann()
		, fExpo()
	{
		std::cout 	<< "[i] thermalizer:"
					<< " <pt>=" << fMeanPt 
					<< " max mult.=" << fMaxMultiplicity 
					<< " max |eta|=" << fMaxAbsEta
					<< " max dR=" << fMaxDeltaR 
					<< std::endl;
		fBoltzmann = new TF1("Thermalizer_fBoltzmann", &boltzmann, 0, fMeanPt * 10., 1);
		fBoltzmann->SetParameter(0, fMeanPt);

		fExpo = new TF1("Thermalizer_fExpo", "expo", 0, 1);
		fExpo->SetParameter(0,  1);
		fExpo->SetParameter(1, -1);
		// fGauss = new TF1("Thermalizer_fGauss", &gauss, 1, 1, 3);
		// fGauss->SetParameter(0, 1);
		// fGauss->SetParameter(1, 0);
		// fGauss->SetParameter(2, 1);
	}

	Thermalizer::Thermalizer(Double_t meanpt, Double_t max_multiplicity, Double_t max_delta_R, Double_t maxabseta)
		: TObject()
		, fMeanPt(meanpt)
		, fMaxMultiplicity(max_multiplicity)
		, fMaxDeltaR(max_delta_R)
		, fMaxAbsEta(maxabseta)
		, fRandom(0)
		, fBoltzmann()
		, fExpo(0)
	{
		std::cout 	<< "[i] thermalizer:"
					<< " <pt>=" << fMeanPt 
					<< " max mult.=" << fMaxMultiplicity 
					<< " max |eta|=" << fMaxAbsEta
					<< " max dR=" << fMaxDeltaR 
					<< std::endl;
		fBoltzmann = new TF1("Thermalizer_fBoltzmann", &boltzmann, 0, fMeanPt * 10., 1);
		fBoltzmann->SetParameter(0, fMeanPt);

		fExpo = new TF1("Thermalizer_fExpo", "expo", 0, 1);
		fExpo->SetParameter(0,  1);
		fExpo->SetParameter(1, -1);
	}

	Thermalizer::~Thermalizer()
	{
		delete fBoltzmann;
	}

	std::vector<fastjet::PseudoJet> Thermalizer::thermalize(Double_t pt, Double_t eta, Double_t phi, Double_t mass)
	{
		std::vector<fastjet::PseudoJet> outv;
		Double_t sumpt = 0;
		Int_t multiplicity = 0;
		// std::cout << "[ ] thermalize: pt=" << pt << " phi=" << phi << " eta=" << eta << " maxR=" << fMaxDeltaR << " maxmult=" << fMaxMultiplicity << std::endl;
		bool multiplicity_condition = (multiplicity < fMaxMultiplicity);
		if (fMaxMultiplicity < 0)
		{
			multiplicity_condition = true;
		}
		// while (multiplicity_condition && sumpt < pt)
		while (sumpt < pt)
		{
			Double_t _maxpt = pt - sumpt;
			if (_maxpt < 1e-4) break;
			Double_t _pt = 0;
			if (fMaxMultiplicity < 0)
			{
				_pt = fBoltzmann->GetRandom(0, _maxpt);
			}
			else
			{
				if (fMaxMultiplicity - multiplicity == 1)
				{
					_pt = pt - sumpt;
					if (_pt <= 0)
						break;
				}
				else
				{
					//_pt = fBoltzmann->GetRandom(0, _maxpt) + (pt - sumpt) / fMaxMultiplicity;
					// _pt = fBoltzmann->GetRandom(0, _maxpt / (fMaxMultiplicity - multiplicity));
					_pt = fRandom.Rndm() * _maxpt;
				}
			}

			if (_pt > _maxpt)
			{
				_pt = _maxpt;
			}

			Double_t _eta = 0;
			Double_t _phi = 0;
			Double_t _r   = 0;
			if (fMaxDeltaR > -1)
			{
				// regulate random eta and phi
				// - generate random R
				// _r = fRandom.Rndm() * fMaxDeltaR;
				fExpo->SetParameter(1, -1. * _pt);
				_r = fMaxDeltaR * fExpo->GetRandom();
				// generate random phi within R
				_phi = fRandom.Rndm() * _r * 2. - _r;
				// generate random eta within sqrt(R^2 - phi^2)
				_eta = fRandom.Rndm() * TMath::Sqrt(_r * _r - _phi * _phi);
			}
			else
			{
				// completely random event except the maxabs eta
				_eta = fRandom.Rndm() * fMaxAbsEta * 2. - fMaxAbsEta;
				_phi = fRandom.Rndm() * TMath::Pi() * 2.;
			}

			if (_pt > 1e-4)
			{
				ROOT::Math::PtEtaPhiMVector _lv(_pt, eta + _eta, phi + _phi, mass);
				fastjet::PseudoJet psj(_lv.px(), _lv.py(), _lv.pz(), _lv.e());
				outv.push_back(psj);
			}
			sumpt += _pt;
			multiplicity += 1;

			multiplicity_condition = (multiplicity < fMaxMultiplicity);
			if (fMaxMultiplicity < 0)
			{
				multiplicity_condition = true;
			}

			// std::cout << "    " << multiplicity << " pt=" << _pt << " z=" << _pt / pt << " phi=" << _phi << " eta=" << _eta << " _r=" << _r;
			// std::cout << "    sumpt=" << sumpt << std::endl;

			if (multiplicity_condition == false)
				break;
		}
		// std::cout << "   - sum pt=" << sumpt << " mult=" << multiplicity << std::endl;
		return outv;
	}

	std::vector<fastjet::PseudoJet> Thermalizer::thermalize_vector(const std::vector<fastjet::PseudoJet> &parts)
	{
		std::vector<fastjet::PseudoJet> outv;
		for (auto &p : parts)
		{
			auto _v = thermalize(p);
			outv.insert(outv.end(), _v.begin(), _v.end());
		}
		return outv;
	}

}