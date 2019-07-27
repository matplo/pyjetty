#include "pythiahepmc.hh"

#include <ios>

namespace HepMCTools
{
Pythia8HepMCWrapper::Pythia8HepMCWrapper()
	: fOutputName("pythia8_hepmc.dat")
	, fOutput(fOutputName.c_str(), std::ios::out)
	, fPythiaToHepMC()
	, fiEv(0)
	{
		;
	}

Pythia8HepMCWrapper::Pythia8HepMCWrapper(const char *fname)
	: fOutputName(fname)
	, fOutput(fOutputName.c_str(), std::ios::out)
	, fPythiaToHepMC()
	, fiEv(0)
	{
		;
	}

Pythia8HepMCWrapper::~Pythia8HepMCWrapper()
{
	;
}

bool Pythia8HepMCWrapper::fillEvent(Pythia8::Pythia &pythia)
{
	HepMC::GenEvent* hepmc_event = new HepMC::GenEvent();
	bool _filled = fPythiaToHepMC.fill_next_event( pythia.event, hepmc_event, fiEv, &pythia.info, &pythia.settings);
	if (_filled == false)
	{
		std::cerr << "[error] Pythia8HepMCWrapper::fillEvent false" << std::endl;
	}
	fOutput << hepmc_event;
	delete hepmc_event;
	return _filled;
}
};
