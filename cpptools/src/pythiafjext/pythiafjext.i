%module pythiafjext
%include "std_vector.i"
%{
	#define SWIG_FILE_WITH_INIT
	#include <Pythia8/Pythia.h>
	#define SWIG
	#include <fastjet/PseudoJet.hh>
	#include "pyfjtools.hh"
%}

%include "pyfjtools.hh"

