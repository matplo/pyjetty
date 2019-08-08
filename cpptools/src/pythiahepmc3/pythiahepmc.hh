#ifndef PYJETTY_PYTHIA8HEPMCWRAPPER_HH
#define PYJETTY_PYTHIA8HEPMCWRAPPER_HH

#include <Pythia8/Pythia.h>
#include "Pythia8/Pythia8ToHepMC3.h"

#if HEPMC31
#define HEPMC_ALIAS HepMC3
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAscii.h"
#else
#define HEPMC_ALIAS HepMC
#include "HepMC/GenEvent.h"
#include "HepMC3/WriterAscii.h"
#endif

#include <string>

namespace HepMCTools
{
	class Pythia8HepMC3Wrapper
	{
	public:
		Pythia8HepMC3Wrapper();
		Pythia8HepMC3Wrapper(const char *fname);
		~Pythia8HepMC3Wrapper();
		bool fillEvent(Pythia8::Pythia &pythia);
	private:
		std::string fOutputName;
		HEPMC_ALIAS::WriterAscii *fOutput;
		HEPMC_ALIAS::Pythia8ToHepMC3 fPythiaToHepMC;
		int fiEv;
	};
}

#endif
