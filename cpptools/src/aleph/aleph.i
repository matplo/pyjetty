%module aleph
%{
	#include "aleph.hh"
%}

%include "std_string.i"
%include "std_vector.i"

// Process symbols in header

/*
%typemap(out) vector<vector<double>>& 
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
	%template(AlephParticleVector) vector<Aleph::Particle>;

    %template(IntVector) vector<int>;
    %template(DoubleVector) vector<double>;
    %template(VectorDoubleVector) vector< vector<double> >;
    %template(StringVector) vector<string>;
    %template(ConstCharVector) vector<const char*>;

	// %template(vstring) vector <string>;
	// %template(vdouble) vector <double>;
	// %template(vvdouble) vector< vector<double> >;
}

%include "aleph.hh"
