%module fjext
%{
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
%apply (double* IN_ARRAY1, int DIM1) {(double* px, int npx), (double* py, int npy), (double* pz, int npz)};
%apply (double* IN_ARRAY1, int DIM1) {(double* px, int npx), (double* py, int npy), (double* pz, int npz), (double* e, int ne)};
%apply (double* IN_ARRAY1, int DIM1) {(double* px, int npx), (double* py, int npy), (double* pz, int npz), (double* m, int nm)};
%include "fjtools.hh"
%clear (double* seq, int n);
%clear (double *pt, int npt, double *eta, int neta, double *phi, int nphi);
%clear (double* px, int npx), (double* py, int npy), (double* pz, int npz);
%clear (double* px, int npx), (double* py, int npy), (double* pz, int npz), (double* e, int ne);
%clear (double* px, int npx), (double* py, int npy), (double* pz, int npz), (double* m, int nm);
