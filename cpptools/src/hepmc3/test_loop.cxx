#include "test_loop.hh"

#include <iostream>

#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"

// something from jetscape via James
void test_loop(const char *fname)
{
  // Read HepMC file
  std::cout << "[i] test_loop HEPMC3 on " << fname << std::endl;
  // std::string hepmcFile = outputDirBin.append(fname);
  std::string hepmcFile(fname);
  HepMC3::ReaderAscii reader(hepmcFile.c_str());
  // Loop over HepMC events, and call analysis task to process them
  int nevents = 0;
  while (!reader.failed()) 
  {
    // Read event
    HepMC3::GenEvent event(HepMC3::Units::GEV,HepMC3::Units::MM);
    reader.read_event(event);
    if (reader.failed()) 
    {
      break;
    }
    //analyzer->AnalyzeEvent(event);  
    nevents++;
  }
  reader.close();
  std::cout << "[i] read " << nevents << " events" << std::endl;
}