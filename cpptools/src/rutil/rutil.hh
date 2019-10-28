#ifndef __PYJETTY_RUTIL__HH
#define __PYJETTY_RUTIL__HH

#include <TObject.h>

namespace RUtil
{
	class Test : public TObject
	{
	public:
		Test() : TObject()
		{;}
		virtual ~Test() {;}
		void setMember(Double_t v) {fMember = v;}
		Double_t getMember() {return fMember;}
	private:
		Double_t fMember;

	ClassDef(Test, 1)
	};
};
#endif