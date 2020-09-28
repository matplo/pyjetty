%module fjtools
%include "std_vector.i"
%template(IntVector) std::vector<int>;
%{
	#define SWIG_FILE_WITH_INIT
	// #include <Pythia8/Pythia.h>
	#include <TF1.h>
 	#include <fastjet/FunctionOfPseudoJet.hh>
 	#include <fastjet/PseudoJet.hh>
	#define SWIG
	#include "fjtools.hh"
	#include "djtools.hh"
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

%apply (double* IN_ARRAY1, int DIM1) {(double* pt, int npt), (double* eta, int neta), (double* phi, int nphi)};
%apply (double* IN_ARRAY1, int DIM1) {(double* pt, int npt), (double* eta, int neta), (double* phi, int nphi), (double* m, int nm)};
%include "djtools.hh"
%clear (double *pt, int npt, double *eta, int neta, double *phi, int nphi);
%clear (double *pt, int npt, double *eta, int neta, double *phi, int nphi), (double* m, int nm);
