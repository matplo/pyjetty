#include "test_loop.hh"

#include <iostream>

#if HEPMC31
#define HEPMC_ALIAS HepMC3
#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#else
#define HEPMC_ALIAS HepMC
#include "HepMC/GenEvent.h"
#include "HepMC/ReaderAscii.h"
#endif


// something from jetscape via James
void test_loop(const char *fname)
{
  // Read HepMC file
  std::cout << "[i] test_loop HEPMC3 on " << fname << std::endl;
  // std::string hepmcFile = outputDirBin.append(fname);
  std::string hepmcFile(fname);
  HEPMC_ALIAS ::ReaderAscii reader(hepmcFile.c_str());
  // Loop over HepMC events, and call analysis task to process them
  int nevents = 0;
  while (!reader.failed()) 
  {
    // Read event
    HEPMC_ALIAS ::GenEvent event(HEPMC_ALIAS ::Units::GEV,HEPMC_ALIAS ::Units::MM);
    reader.read_event(event);
    if (reader.failed()) 
    {
      break;
    }
    //analyzer->AnalyzeEvent(event);  
    nevents++;
    if (nevents % 10000 == 0) std::cout << " - events read " << nevents << std::endl;
  }
  reader.close();
  std::cout << "[i] read " << nevents << " events" << std::endl;
}