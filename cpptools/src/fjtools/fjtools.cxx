#include "fjtools.hh"

#include <TF1.h>
#include <TMath.h>
#include <TString.h>
#include <TRandom.h>
#include <iostream>

namespace PyJettyFJTools
{
	// fraction of pT of jet j1 contained in j0 - constit by constit
	double matched_pt(const fastjet::PseudoJet &j0, const fastjet::PseudoJet &j1)
	{
		double pt_sum = 0;
		double pt_sum_1 = 0;
		for (auto &p1 : j1.constituents())
			pt_sum_1 += p1.pt();
		for (auto &p0 : j0.constituents())
		{
			for (auto &p1 : j1.constituents())
			{
				if (p0.user_index() == p1.user_index())
				{
					pt_sum += p1.pt();
				}
			}
		}
		return (pt_sum / pt_sum_1);
	}

	// return indices of jets matched to j jet	
	std::vector<int> matched_Ry(const fastjet::PseudoJet &j, const std::vector<fastjet::PseudoJet> &v, double Rmatch)
	{
		std::vector<int> retv;
		for(std::size_t i = 0; i < v.size(); ++i)
		{
			if (j.delta_R(v[i]) < Rmatch)
			{
				retv.push_back(i);
			}
		}
		return retv;
	}

	// return indices of jets matched to j jet	
	std::vector<int> matched_Reta(const fastjet::PseudoJet &j, const std::vector<fastjet::PseudoJet> &v, double Rmatch)
	{
		std::vector<int> retv;
		for(std::size_t i = 0; i < v.size(); ++i)
		{
			double dR = TMath::Sqrt(TMath::Power(j.delta_phi_to(v[i]), 2.) + TMath::Power(j.eta() - v[i].eta(), 2.));
			if (dR < Rmatch)
			{
				retv.push_back(i);
			}
		}
		return retv;
	}

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

	std::string BoltzmannBackground::description()
	{
		std::string s;
		s = TString::Format("\n[i] BoltzmannBackground with \n    fmean_pt=%f\n    fmin_pt=%f\n    fmax_pt=%f\n", fmean_pt, fmin_pt, fmax_pt);
		return s;
	}

	BoltzmannBackground::BoltzmannBackground()
	: funbg(0)
	, fmean_pt(-1)
	, fmin_pt(-1)
	, fmax_pt(-1)
	, fparts()
	{
		reset(0.7, 0.15, 5.);
	}

	BoltzmannBackground::BoltzmannBackground(double mean_pt, double min_pt, double max_pt)
	: funbg(0)
	, fmean_pt(-1)
	, fmin_pt(-1)
	, fmax_pt(-1)
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
	std::vector<fastjet::PseudoJet> BoltzmannBackground::subtract(const std::vector<fastjet::PseudoJet> &v, double mean_pt, int n)
	{
		fparts.clear();
		double _count = n;
		double _mean_pt = mean_pt;
		double _total_pt = _count * _mean_pt;
		if (_mean_pt < 0 || _count < 0)
		{
			for (auto &p : v)
			{
				if (p.pt() <= fmax_pt)
				{
					_total_pt += p.pt();
					_count += 1.;
				}
			}
			if (_count > 0)
			{
				_mean_pt = _total_pt / _count;
			}
		}
		if (_count > 0)
			reset(_mean_pt, fmin_pt, fmax_pt);
		std::cout << "* total_pt=" << TString::Format("%f", _total_pt).Data() << std::endl;
		for (auto &p : v)
		{
			if (p.pt() > fmax_pt)
			{
				fastjet::PseudoJet _p(p.px(), p.py(), p.pz(), p.e());
				_p.set_user_index(p.user_index());
				fparts.push_back(_p);		
			}
			else
			{
				// std::cout << TString::Format("%f", fmin_pt).Data() << std::endl;
				// if (p.pt() > fmin_pt)
				// 	std::cout << TString::Format("%f %f %f", p.pt(), funbg->Eval(p.pt()), _count).Data() << std::endl;
				// double fraction = funbg->Eval(p.pt()) * funbg->GetNpx() / 2.;
				// double fraction = funbg->Eval(p.pt()) * (fmax_pt - fmean_pt) * funbg->GetNpx();
				double fraction = funbg->Eval(p.pt()) * _total_pt;
				if ((p.pt() - fraction) > fmin_pt)
				{
					fastjet::PseudoJet _p;
					_p.reset_PtYPhiM(p.pt() - fraction, p.rap(), p.phi(), p.m());
					// std::cout << TString::Format("%f - %f = %f", p.pt(), fraction, _p.pt()).Data() << std::endl;
					_p.set_user_index(p.user_index());
					fparts.push_back(_p);		
				}
			}
		}
		return fparts;
	}

	// subtract particles (according to the probability... - fixed to 1 in maxpt range)
	std::vector<fastjet::PseudoJet> BoltzmannBackground::subtract_recalc_from_vector(const std::vector<fastjet::PseudoJet> &v)
	{
		fparts.clear();
		double _total_pt = 0;
		double _count = 0;
		for (auto &p : v)
		{
			if (p.pt() <= fmax_pt)
			{
				_total_pt += p.pt();
				_count += 1.;
			}
		}
		double _mean_pt = 0;
		if (_count > 0)
		{
			_mean_pt = _total_pt / _count;
			reset(_mean_pt, fmin_pt, fmax_pt);
		}
		for (auto &p : v)
		{
			if (p.pt() > fmax_pt)
			{
				fastjet::PseudoJet _p(p.px(), p.py(), p.pz(), p.e());
				_p.set_user_index(p.user_index());
				fparts.push_back(_p);		
			}
			else
			{
				// std::cout << TString::Format("%f", fmin_pt).Data() << std::endl;
				// if (p.pt() > fmin_pt)
				// 	std::cout << TString::Format("%f %f %f", p.pt(), funbg->Eval(p.pt()), _count).Data() << std::endl;
				// double fraction = funbg->Eval(p.pt()) * funbg->GetNpx() / 2.;
				// double fraction = funbg->Eval(p.pt()) * (fmax_pt - fmean_pt) * funbg->GetNpx();
				// std::cout << TString::Format("%f", _total_pt).Data() << std::endl;
				double fraction = funbg->Eval(p.pt()) * _total_pt;
				if ((p.pt() - fraction) > fmin_pt)
				{
					fastjet::PseudoJet _p;
					_p.reset_PtYPhiM(p.pt() - fraction, p.rap(), p.phi(), p.m());
					// std::cout << TString::Format("%f - %f = %f", p.pt(), fraction, _p.pt()).Data() << std::endl;
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