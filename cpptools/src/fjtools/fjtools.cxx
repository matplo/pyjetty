#include "fjtools.hh"

#include <TF1.h>
#include <TMath.h>
#include <TString.h>
#include <TRandom.h>

namespace PyJettyFJTools
{
	double boltzmann_norm(double *x, double *par)
	{
		// double fval = par[1] / x[0] * x[0] * TMath::Exp(-(2. / par[0]) * x[0]);
		// double fval = par[1] / x[0] * x[0] * TMath::Exp(-(x[0] / par[0]));
		// return fval;
		return par[1] * x[0] * TMath::Exp(-(x[0] / (par[0]/2.)));
	}

	double boltzmann(double *x, double *par)
	{
		return x[0] * TMath::Exp(-(x[0] / (par[0]/2.)));
	}

	BoltzmannBackground::BoltzmannBackground()
	: funbg(0)
	, fmean_pt(0.7)
	, fmin_pt(0.15)
	, fmax_pt(5.)
	, fparts()
	{
		reset(0.7, 0.15, 5.);
	}

	BoltzmannBackground::BoltzmannBackground(double mean_pt, double min_pt, double max_pt)
	: funbg(0)
	, fmean_pt(mean_pt)
	, fmin_pt(mean_pt)
	, fmax_pt(max_pt)
	, fparts()
	{
		reset(mean_pt, min_pt, max_pt);
	}


	void BoltzmannBackground::reset(double mean_pt, double min_pt, double max_pt)
	{
		bool _reinit = false;
		if (fmin_pt != min_pt)
		{
			_reinit = true;
			fmin_pt = min_pt;
		}
		if (fmax_pt != max_pt)
		{
			_reinit = true;
			fmax_pt = max_pt;
		}
		if (fmean_pt != mean_pt)
		{
			_reinit = true;
			fmean_pt = mean_pt;
		}
		if (_reinit)
		{
			TF1 _tmpf("_tmpf", boltzmann, fmin_pt, fmax_pt, 1);
			_tmpf.SetParameter(0, fmean_pt);
			double _int = _tmpf.Integral(fmin_pt, fmax_pt);
			if (_int == 0)
				_int = 1e-9;
			if (funbg)
			{
				delete funbg;
			}
			funbg = new TF1("BoltzmannBackground_funBG", boltzmann_norm, fmin_pt, fmax_pt, 2);
			funbg->SetParameter(0, fmean_pt);
			funbg->SetParameter(1, 1./_int);		
		}
		fparts.clear();
	}

	double BoltzmannBackground::eval(double x)
	{
		return funbg->Eval(x);
	}

	double BoltzmannBackground::integral()
	{
		return funbg->Integral(fmin_pt, fmax_pt);
	}

	std::string BoltzmannBackground::get_formula()
	{
		std::string s;
		s = TString::Format("%f * x[0] * TMath::Exp(-(x[0] / %f))", funbg->GetParameter(1), funbg->GetParameter(0)).Data();
		return s;
	}

	double BoltzmannBackground::get_constant()
	{
		return funbg->GetParameter(1);
	}

	// getNrandom particles w/ indexoffset ...
	std::vector<fastjet::PseudoJet> BoltzmannBackground::generate(int nparts, double max_eta, int offset)
	{	
		fparts.clear();
		for (int i = 0; i < nparts; i++)
		{
			double _pt  = funbg->GetRandom(fmin_pt, fmax_pt);
			double _eta = gRandom->Rndm() * max_eta * 2. - max_eta;
			double _phi = gRandom->Rndm() * TMath::Pi() * 2. - TMath::Pi();
			fastjet::PseudoJet _p;
			_p.reset_PtYPhiM (_pt, _eta, _phi, 0.0);
			_p.set_user_index(i + offset);
			fparts.push_back(_p);		
		}
		return fparts;
	}

	// subtract particles (according to the probability... - fixed to 1 in maxpt range)
	std::vector<fastjet::PseudoJet> BoltzmannBackground::subtract(const std::vector<fastjet::PseudoJet> &v)
	{
		fparts.clear();
		double _total_pt = 0;
		for (auto &p : v)
		{
			_total_pt += p.pt();
		}
		for (auto &p : v)
		{
			if (p.pt() > fmax_pt)
			{
				fparts.push_back(p);
				continue;
			}
			else
			{
				double fraction = funbg->Eval(p.pt()) * _total_pt;
				if (p.pt() - fraction > fmin_pt)
				{
					fastjet::PseudoJet _p;
					_p.reset_PtYPhiM(p.pt() - fraction, p.rap(), p.phi(), p.m());
					_p.set_user_index(p.user_index());
					fparts.push_back(_p);		
				}
			}
		}
		return fparts;
	}

	BoltzmannBackground::~BoltzmannBackground()
	{
		delete funbg;
	}
}