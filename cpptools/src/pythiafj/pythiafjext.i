%module pyjetty
%{
/* Includes the header in the wrapper code */
/* #define SWIG_FILE_WITH_INIT */
#include "pyfjtools.hh"
%}

%include "std_string.i"
%include "std_vector.i"
%include "typemaps.i"
/*
%include "../numpy.i"
%init %{
	import_array();
%}
%fragment("NumPy_Fragments");
*/

/* Parse the header file to generate wrappers */
%include "pyfjtools.hh"