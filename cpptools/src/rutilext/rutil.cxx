#include <rutil.hh>

#include <TFile.h>

#include <iostream>

ClassImp(RUtilExt::Test);

namespace RUtilExt
{
    typedef double (*prior_scale_func)(const double & obs_true,
                                       const double & content,
                                       const double & prior_variation_parameter);

    //---------------------------------------------------------------
    // Prior scaling functions
    //---------------------------------------------------------------
    double prior_scale_func_0(const double & obs_true, const double & content,
                              const double & prior_variation_parameter) {
        // power law
        return std::pow(obs_true, prior_variation_parameter);
    }

    double prior_scale_func_1(const double & obs_true, const double & content,
                              const double & prior_variation_parameter) {
        // linear scaling of distributions
        return std::abs(1 + prior_variation_parameter * (2 * obs_true - 1));
    }

    double prior_scale_func_2(const double & obs_true, const double & content,
                              const double & prior_variation_parameter) {
        // sharpening/smoothing the distributions
        return std::pow(content, 1 + prior_variation_parameter);
    }

    double prior_scale_func_3(const double & obs_true, const double & content,
                              const double & prior_variation_parameter) {
        return (1 + obs_true);
    }

    double prior_scale_func_4(const double & obs_true, const double & content,
                              const double & prior_variation_parameter) {
        
        // Ax+B, where A=slope, B=offset at z=0
        // For 0.7<z<1.0, dz = 0.3 --> A = 1/dz, B = 1-(1-dz/2)*A
        float dz = 0.3;
        float A = prior_variation_parameter*1./dz;
        return (A*obs_true + 1 - (1-dz/2.)*A);
    }

    double prior_scale_func_def(const double & obs_true, const double & content,
                                const double & prior_variation_parameter) {
        return obs_true;
    }


    //---------------------------------------------------------------
    // Clean up dynamic memory from a histogram on the heap
    void delete_h(TH2* h) {
        delete h;
        return;
    }

    void delete_h(THn* h) {
        delete h;
        return;
    }

	/////////////////////////////////////////////////////////////////////////////////
	// Iterate through sorted lists and return list of bools for items of a in b
	bool* sorted_match(const int* a, const int a_len, const int* b, const int b_len) {
		bool* isin = new bool[a_len]();
		if (!(b_len)) { return isin; }
		int j = 0;  // index in b
		for(int i = 0; i < a_len; i++) {
			while (b[j] < a[i]) if (++j > b_len) return isin;
			if (a[i] == b[j]) isin[i] = true;
		} // for (i)
		return isin;
	} // sorted_match

    void fillTH1weighted(TH1 *h, const std::vector<double> &x, const std::vector<double> &w)
    {
        std::cout << " ->? " << w.size() << std::endl;
        for (size_t i = 0; i < w.size(); i++)
        {
            // std::cout << w.size() << " " << i << std::endl;
            // h->Fill(x[i], w[i]);
            // double _x = x[i] - w[i];
        }
    }

}  // namespace RUtil
