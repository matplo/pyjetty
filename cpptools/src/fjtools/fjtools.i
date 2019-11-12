%module fjtools
%include "std_vector.i"
%{
	#define SWIG_FILE_WITH_INIT
	// #include <Pythia8/Pythia.h>
	#include <TF1.h>
 	#include <fastjet/FunctionOfPseudoJet.hh>
 	#include <fastjet/PseudoJet.hh>
	#define SWIG
	#include "fjtools.hh"
%}

%include "std_string.i"
%include "std_vector.i"
%include "typemaps.i"
%include "../numpy.i"
%init %{
	import_array();
%}
%fragment("NumPy_Fragments");

%apply (int* IN_ARRAY1, int DIM1) {(int* selection, int nsel)};
%include "fjtools.hh"
%clear (int* selection, int nsel);
