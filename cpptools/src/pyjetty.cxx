#include "pyjetty.hh"

namespace PyJetty
{
	double sum(double *seq, int n)
	{
		double _sum = 0;
		for (unsigned int i = 0; i < n; i++)
		{
			_sum += seq[i];
		}
		return _sum;
	}

	TestClass::TestClass() {;}
	TestClass::~TestClass() {;}
};