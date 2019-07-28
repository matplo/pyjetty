%module pythiaext

%include "std_string.i"
%include "std_vector.i"

// Add necessary symbols to generated header
%{
#include "pythiahepmc.hh"
%}

// Process symbols in header

%include "pythiahepmc.hh"
