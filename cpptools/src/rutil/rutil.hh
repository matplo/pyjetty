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
    double prior_scale_func_def(const double & obs_true, const double & content,
                                const double & prior_variation_parameter);

	void delete_h(TH2* h);
	void delete_h(THn* h);

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
        THnF* rebin_thn(const std::string & response_file_name,
                        const THnF* thn,
                        const std::string & name_thn_rebinned,
                        const std::string & name_roounfold,
                        const int & n_dim,
                        const int & n_pt_bins_det,
                        const double* det_pt_bin_array,
                        const int & n_obs_bins_det,
                        const double* det_bin_array,
                        const int & n_pt_bins_truth,
                        const double* truth_pt_bin_array,
                        const int & n_obs_bins_truth,
                        const double* truth_bin_array,
                        const std::string & label="",
                        const double & prior_variation_parameter=0.,
                        const int & prior_option=1,
                        const bool move_underflow=false,
                        const bool do_roounfoldresponse=true);

    private:
        // Create empty THn using provided axes
        THnF* create_empty_thn(const char* name, const int & n_dim,
                               const int & n_pt_bins_det, const double* det_pt_bin_array,
                               const int & n_obs_bins_det, const double* det_bin_array,
                               const int & n_pt_bins_truth, const double* truth_pt_bin_array,
                               const int & n_obs_bins_truth, const double* truth_bin_array);

        RooUnfoldResponse* create_empty_roounfoldresponse(
            THnF* thn_rebinned, const std::string & name_roounfold, const std::string & label);

        // Fill empty thn_rebinned with data from thn
        void fill_rebinned_thn(const std::string & response_file_name, const THnF* thn,
                               THnF* thn_rebinned, const unsigned int & n_dim,
                               const prior_scale_func prior_scale_f,
                               const bool do_roounfoldresponse=true,
                               RooUnfoldResponse* roounfold_response=nullptr,
                               const double & prior_variation_parameter=0.,
                               const bool move_underflow=false);

        // Set scaling of prior
        prior_scale_func prior_scale_factor_obs(const int & option);

    ClassDef(HistUtils, 1)
    };

};
#endif
