%module hepmc3ext
%{
/* Includes the header in the wrapper code */
#define SWIG_FILE_WITH_INIT
#include "test_loop.hh"
%}

%include "std_string.i"
%include "std_vector.i"
%include "typemaps.i"
%include "../numpy.i"
%init %{
	import_array();
%}
%fragment("NumPy_Fragments");

/* Parse the header file to generate wrappers */
%include "test_loop.hh"
