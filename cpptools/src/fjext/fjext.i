%module pyjetty
%{
/* Includes the header in the wrapper code */
#define SWIG_FILE_WITH_INIT
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

/* Parse the header file to generate wrappers */
%apply (double* IN_ARRAY1, int DIM1) {(double* seq, int n)};
%apply (double* IN_ARRAY1, int DIM1) {(double* pt, int npt), (double* eta, int neta), (double* phi, int nphi)};
%include "fjtools.hh"
%clear (double* seq, int n);
%clear (double *pt, int npt, double *eta, int neta, double *phi, int nphi);
