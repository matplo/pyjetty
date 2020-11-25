#ifndef __PYJETTY_RUTIL__HH
#define __PYJETTY_RUTIL__HH

#include <TObject.h>
#include <TH2.h>
#include <THn.h>

//#if USE_ROOUNFOLD
#include <RooUnfoldResponse.h>
//#endif

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

    //------------------------------------------------------
    // Rebinning utilities
    class HistUtils : public TObject
    {
    public:
        HistUtils() : TObject()
        {;}
        virtual ~HistUtils() {;}
        
        // Rebin 2D histogram h with name hname using axes given by x_bins and y_bins
        TH2F* rebin_th2(TH2F & h_to_rebin, char* hname, double* x_bins, int n_x_bins,
                        double* y_bins, int n_y_bins, bool move_y_underflow = false);

        // Rebin N-dimensional THn to a new histogram with name name_thn_rebinned using provided axes
        // WARNING: currently requires n_dim = 4
        THnF* rebin_thn(std::string response_file_name,
                        THnF* thn,
                        const std::string & name_thn_rebinned,
                        const std::string & name_roounfold,
                        const unsigned int & n_dim,
                        const int* axes_n_bins,
                        double* axes_bin_array0,
                        double* axes_bin_array1,
                        double* axes_bin_array2,
                        double* axes_bin_array3,
                        const std::string label="",
                        const double & prior_variation_parameter=0.,
                        int prior_option=1,
                        bool move_underflow=false,
                        bool do_roounfoldresponse=true);

        // Create empty THn using provided axes
        THnF* create_empty_thn(const char* name, const unsigned int & n_dim, const char** axes_titles,
                               const int* axes_n_bins, std::vector<double*> axes_bin_arrays);

        // Fill empty thn_rebinned with data from thn
        void fill_rebinned_thn(std::string response_file_name, THnF* thn,
                               THnF* thn_rebinned, const unsigned int & n_dim,
                               bool do_roounfoldresponse=true,
                               RooUnfoldResponse* roounfold_response=nullptr,
                               const double & prior_variation_parameter=0.,
                               int prior_option=1,
                               bool move_underflow=false);

        // Set scaling of prior
        double prior_scale_factor_obs(double obs_true, double content,
                                      double prior_variation_parameter, int option);

    ClassDef(HistUtils, 1)
    };

};
#endif
