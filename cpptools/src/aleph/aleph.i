%module aleph
%{
	#include "aleph.hh"
%}

%include "std_string.i"
%include "std_vector.i"

// Process symbols in header

%include "aleph.hh"

/*
%typemap(out) std::vector<std::vector<double>>& 
{
    for(int i = 0; i < $1->size(); ++i)
    {       
        int subLength = $1->data()[i].size();
        npy_intp dims[] = { subLength };
        PyObject* temp = PyArray_SimpleNewFromData(1, dims, NPY_INT, $1->data()[i].data());
        $result = SWIG_Python_AppendOutput($result, temp);
    }       
}
*/

namespace std
{
	%template(AlephParticleVector) std::vector<Aleph::Particle>;
	%template(vstring) std::vector <std::string>;
	%template(vdouble) std::vector <double>;
	%template(vvdouble) std::vector< std::vector<double> >;
}
