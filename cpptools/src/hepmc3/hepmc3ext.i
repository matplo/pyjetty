%module hepmc3ext
%{
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

%include "test_loop.hh"
