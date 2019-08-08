#include "pythiahepmc.hh"

#include <string>
#include <iostream>

using namespace HEPMC_ALIAS;

namespace HepMCTools
{
Pythia8HepMC3Wrapper::Pythia8HepMC3Wrapper()
	: fOutputName("pythia8_hepmc.dat")
	, fOutput(0)
	, fPythiaToHepMC()
	, fiEv(0)
	{
		;
	}

Pythia8HepMC3Wrapper::Pythia8HepMC3Wrapper(const char *fname)
	: fOutputName(fname)
	, fOutput(0)
	, fPythiaToHepMC()
	, fiEv(0)
	{
		;
	}

Pythia8HepMC3Wrapper::~Pythia8HepMC3Wrapper()
{
	if (fOutput)
		delete fOutput;
}

bool Pythia8HepMC3Wrapper::fillEvent(Pythia8::Pythia &pythia)
{
	if (fiEv == 0)
	{
		shared_ptr<GenRunInfo> run = make_shared<GenRunInfo>();
		struct GenRunInfo::ToolInfo generator={std::string("Pythia8"),std::to_string(PYTHIA_VERSION).substr(0,5),std::string("Used generator")};
		run->tools().push_back(generator);
		pythia.settings.writeFile(fOutputName + ".cmnd");
		struct GenRunInfo::ToolInfo config={std::string(fOutputName + ".cmnd"),"1.0",std::string("Control cards")};
		run->tools().push_back(config);
		std::vector<std::string> names;
		for (int iWeight=0; iWeight < pythia.info.nWeights(); ++iWeight) 
		{
			std::string s=pythia.info.weightLabel(iWeight);
			if (!s.length()) s=std::to_string((long long int)iWeight);
			names.push_back(s);
		}
		if (!names.size()) 
			names.push_back("default");
		run->set_weight_names(names);
		fOutput = new WriterAscii(fOutputName, run);
	}

	bool _filled = false;
	if (fOutput)
	{
		GenEvent hepmc_event( Units::GEV, Units::MM );
		_filled = fPythiaToHepMC.fill_next_event( pythia.event, &hepmc_event, fiEv, &pythia.info, &pythia.settings);
		if (_filled == false)
		{
			std::cerr << "[error] Pythia8HepMC3Wrapper::fillEvent false" << std::endl;
		}
		fOutput->write_event(hepmc_event);
		fiEv++;
	}
	else
	{
		std::cerr << "[error] Pythia8HepMC3Wrapper::fillEvent false - no valid output 0x" << fOutput << std::endl;		
	}
	return _filled;
}
};
