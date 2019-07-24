%module pyjetty
%{
/* Includes the header in the wrapper code */
#include "pyjetty.hh"
#include "pythiahepmc.hh"
%}

/* Parse the header file to generate wrappers */
%include "pyjetty.hh"
%include "pythiahepmc.hh" 
