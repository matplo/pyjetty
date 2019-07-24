%module pyjetty

%include "std_string.i"
%include "std_vector.i"

%{
/* Includes the header in the wrapper code */
#include "pyjetty.hh"
#include "pythiahepmc.hh"

#include <fastjet/PseudoJet.hh>
#include <Pythia8/Pythia.h>
#include "pyfjtools.hh"
%}

/* Parse the header file to generate wrappers */
%include "pyjetty.hh"
%include "pythiahepmc.hh" 
%include "pyfjtools.hh"