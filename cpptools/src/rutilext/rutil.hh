#ifndef __PYJETTY_RUTIL__HH
#define __PYJETTY_RUTIL__HH

#include <TObject.h>
#include <TH2.h>
#include <THn.h>

#include <vector>

namespace RUtilExt
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


    typedef double (*prior_scale_func)(const double & obs_true,
                                       const double & content,
                                       const double & prior_variation_parameter);

    double prior_scale_func_0(const double & obs_true, const double & content,
                              const double & prior_variation_parameter);
    double prior_scale_func_1(const double & obs_true, const double & content,
                              const double & prior_variation_parameter);
    double prior_scale_func_2(const double & obs_true, const double & content,
                              const double & prior_variation_parameter);
    double prior_scale_func_3(const double & obs_true, const double & content,
                              const double & prior_variation_parameter);
    double prior_scale_func_4(const double & obs_true, const double & content,
                              const double & prior_variation_parameter);
    double prior_scale_func_def(const double & obs_true, const double & content,
                                const double & prior_variation_parameter);

	void delete_h(TH2* h);
	void delete_h(THn* h);

	bool* sorted_match(const int* a, const int a_len, const int* b, const int b_len);

    void fillTH1weighted(TH1 *h, const std::vector<double> &x, const std::vector<double> &w);

};
#endif
