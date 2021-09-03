#ifndef LOOPUTIL_HH
#define LOOPUTIL_HH

#include <TString.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <iostream>
#include <numeric>

namespace LoopUtil
{
	class TPbar
	{
	public:
	        TPbar(Long64_t n)
	                : fN(n)
	                , fNCalls(0)
	                , fSw(new TStopwatch)
	                , fTperIt()
	                , fOldPcent(-1)
	                , fOldRTime(0)
	                {
	                        fSw->Start();
	                }

	        void Update(Long64_t nup = 1)
	        {
	        	if (fN > 0)
	        	{
	                fSw->Stop();
	                fNCalls = fNCalls + nup;
	                Double_t frac = (1.*fNCalls+1)/fN;
	                Double_t pcent = frac * 100.;
	                //if (((pcent - TMath::Floor(pcent)) == 0) || (fNCalls == fN - 1) || nup > 1 || fNCalls == 1)
	                //if (pcent - fOldPcent >= 0.1 || (fNCalls == fN - 1) || nup > 1 || fNCalls == 1)
	                if (fSw->RealTime() - fOldRTime > 1 || (fNCalls == fN - 1) || nup > 1 || fNCalls == 1)
	                {
	                        Double_t etas = 1. * fSw->RealTime() / fNCalls * (fN - fNCalls);
	                        fTperIt.push_back(1. * fSw->RealTime() / fNCalls);
	                        double sum = std::accumulate(fTperIt.begin(), fTperIt.end(), 0.0);
	                        double mean_t_per_it = sum / fTperIt.size();
	                        etas = mean_t_per_it * (fN - fNCalls);
	                        Double_t fh = TMath::Floor(etas / 60. / 60.);
	                        Double_t fm = TMath::Floor(etas / 60. - fh * 60.);
	                        Double_t fs = etas - fh * 60. * 60. - fm * 60.;
	                        TString seta = TString::Format("%02d:%02d:%02d (%1.1f/s)", Int_t(fh), Int_t(fm), Int_t(fs), 1.*fNCalls/fSw->RealTime());
	                        fh = TMath::Floor(fSw->RealTime() / 60. / 60.);
	                        fm = TMath::Floor(fSw->RealTime() / 60. - fh * 60.);
	                        fs = fSw->RealTime() - fh * 60. * 60. - fm * 60.;
	                        TString sela = TString::Format("%02d:%02d:%02d", Int_t(fh), Int_t(fm), Int_t(fs));
	                        std::cout << "\r[i] event #" << fNCalls+1 << " of " << fN << " is " << TString::Format("%3.1f %% ETA:%s T:%s                       \r", pcent, seta.Data(), sela.Data()); std::cout.flush();
	                        fOldPcent = pcent;
	                        fOldRTime = fSw->RealTime();
	                }
	                fSw->Start(kFALSE);
	                if (fNCalls >= fN)
	                {
	                        std::cout << std::endl;
	                        fSw->Stop();
	                        return;
	                }
	            }
	            else
	            {
	                fSw->Stop();
	                fNCalls = fNCalls + nup;
	                if (fSw->RealTime() - fOldRTime > 1 || fNCalls == 1)
	                {
                        fTperIt.push_back(1. * fSw->RealTime() / fNCalls);
                        double sum = std::accumulate(fTperIt.begin(), fTperIt.end(), 0.0);
                        double mean_t_per_it = sum / fTperIt.size();
                        Double_t fh = TMath::Floor(fSw->RealTime() / 60. / 60.);
                        Double_t fm = TMath::Floor(fSw->RealTime() / 60. - fh * 60.);
                        Double_t fs = fSw->RealTime() - fh * 60. * 60. - fm * 60.;
                        TString sela = TString::Format("%02d:%02d:%02d (it/s %f)", Int_t(fh), Int_t(fm), Int_t(fs), 1./mean_t_per_it);
                        std::cout << "\r[i] event #" << fNCalls+1 << " " << TString::Format("T:%s                       \r", sela.Data()); std::cout.flush();
                        fOldRTime = fSw->RealTime();
	                }
	                fSw->Start(kFALSE);
	            }
	        }

	        ~TPbar()
	        {
	                delete fSw;
	        }

	        Long64_t NCalls() {return fNCalls;}

	private:
	        Long64_t fN;
	        Long64_t fNCalls;
	        TStopwatch *fSw;
	        std::vector<Double_t> fTperIt;
	        Double_t fOldPcent;
	        Double_t fOldRTime;
	};
}

#endif
