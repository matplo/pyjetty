#ifndef PYJETTY_PYTHIA8HEPMCWRAPPER_HH
#define PYJETTY_PYTHIA8HEPMCWRAPPER_HH

#include <Pythia8/Pythia.h>
#include <Pythia8Plugins/HepMC2.h>
#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>

#include <string>

namespace HepMCTools
{
	class Pythia8HepMC2Wrapper
	{
	public:
		Pythia8HepMC2Wrapper();
		Pythia8HepMC2Wrapper(const char *fname);
		~Pythia8HepMC2Wrapper();
		bool fillEvent(Pythia8::Pythia &pythia);
	private:
		std::string fOutputName;
		HepMC::IO_GenEvent fOutput;
		HepMC::Pythia8ToHepMC fPythiaToHepMC;
		int fiEv;
	};
}

#endif
