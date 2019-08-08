#include "pythiahepmc.hh"

#include <ios>

namespace HepMCTools
{
Pythia8HepMC2Wrapper::Pythia8HepMC2Wrapper()
	: fOutputName("pythia8_hepmc.dat")
	, fOutput(fOutputName.c_str(), std::ios::out)
	, fPythiaToHepMC()
	, fiEv(0)
	{
		;
	}

Pythia8HepMC2Wrapper::Pythia8HepMC2Wrapper(const char *fname)
	: fOutputName(fname)
	, fOutput(fOutputName.c_str(), std::ios::out)
	, fPythiaToHepMC()
	, fiEv(0)
	{
		;
	}

Pythia8HepMC2Wrapper::~Pythia8HepMC2Wrapper()
{
	;
}

bool Pythia8HepMC2Wrapper::fillEvent(Pythia8::Pythia &pythia)
{
	HepMC::GenEvent* hepmc_event = new HepMC::GenEvent();
	bool _filled = fPythiaToHepMC.fill_next_event( pythia.event, hepmc_event, fiEv, &pythia.info, &pythia.settings);
	if (_filled == false)
	{
		std::cerr << "[error] Pythia8HepMC2Wrapper::fillEvent false" << std::endl;
	}
	fOutput << hepmc_event;
	delete hepmc_event;
	return _filled;
}
};
