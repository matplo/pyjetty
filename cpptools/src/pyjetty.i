%module pyjetty
%{
/* Includes the header in the wrapper code */
#define SWIG_FILE_WITH_INIT
#include "pyjetty.hh"
#include "pythiahepmc.hh"
#include "pyfjtools.hh"
%}

%include "std_string.i"
%include "std_vector.i"
%include "typemaps.i"
%include "numpy.i"
%init %{
	import_array();
%}
%fragment("NumPy_Fragments");

/* Parse the header file to generate wrappers */
%apply (double* IN_ARRAY1, int DIM1) {(double* seq, int n)};
%include "pyjetty.hh"
%clear (double* seq, int n);
%include "pythiahepmc.hh" 
%include "pyfjtools.hh"