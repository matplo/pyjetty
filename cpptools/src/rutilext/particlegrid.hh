#ifndef __PYJETTY_PARTICLEGRID__HH
#define __PYJETTY_PARTICLEGRID__HH

#include <TObject.h>
#include <TH2F.h>
#include <cmath>
#include <fastjet/PseudoJet.hh>

namespace RUtilExt
{
	class ParticleGrid : public TObject
	{
	public:
		ParticleGrid() : TObject(), fGrid(0)
		{;}
		ParticleGrid(Double_t fNbinsEta, Double_t fMinEta, Double_t fMaxEta, Double_t fNbinsPhi, Double_t fMinPhi, Double_t fMaxPhi)
		: TObject()
		, fGrid(0)
		{
			newGrid(fNbinsEta, fMinEta, fMaxEta, fNbinsPhi, fMinPhi, fMaxPhi);
		}

		virtual ~ParticleGrid()
		{
			delete fGrid;
		}

		void fillParticle(Double_t eta, Double_t phi, Double_t e)
		{
			fGrid->Fill(eta, phi, e);
		}

		void fillParticle(const fastjet::PseudoJet &p)
		{
			fGrid->Fill(p.eta(), p.phi(), p.e());
		}

		void fillParticles(const std::vector<fastjet::PseudoJet> &parts)
		{
			for (auto &p : parts)
			{
				fGrid->Fill(p.eta(), p.phi(), p.e());
			}
		}

		std::vector<fastjet::PseudoJet> getGridParticles()
		{
			std::vector<fastjet::PseudoJet> outv;
			for (int ieta = 1; ieta < fGrid->GetNbinsX()+1; ieta++)
			{
				for (int iphi = 1; iphi < fGrid->GetNbinsY()+1; iphi++)
				{
					if (fGrid->GetBinContent(ieta, iphi) > 0)
					{
						fastjet::PseudoJet psj;
						psj.reset_PtYPhiM(	fGrid->GetBinContent(ieta, iphi),
											fGrid->GetXaxis()->GetBinCenter(ieta),
											fGrid->GetYaxis()->GetBinCenter(iphi),
											0.);
						outv.push_back(psj);
					}
				}
			}
			return outv;
		}

		void Reset()
		{
			fGrid->Reset();
		}

	private:
		void newGrid(Double_t fNbinsEta, Double_t fMinEta, Double_t fMaxEta, Double_t fNbinsPhi, Double_t fMinPhi, Double_t fMaxPhi)
		{
			if (fGrid != 0x0)
				delete fGrid;
			fGrid = new TH2F("ParticleGrid", "ParticleGrid", fNbinsEta, fMinEta, fMaxEta, fNbinsPhi, fMinPhi, fMaxPhi);
			fGrid->SetDirectory(0);
		}
		TH2F *fGrid; //!
	ClassDef(ParticleGrid, 1)
	};
};
#endif
