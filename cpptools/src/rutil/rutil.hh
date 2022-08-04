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
    double prior_scale_func_4(const double & obs_true, const double & content,
                              const double & prior_variation_parameter);
    double prior_scale_func_def(const double & obs_true, const double & content,
                                const double & prior_variation_parameter);

	void delete_h(TH2* h);
	void delete_h(THn* h);

	bool* sorted_match(const int* a, const int a_len, const int* b, const int b_len);

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

        // Rebin 2D histogram h with name hname using axes given by x_bins and y_bins
        TH2D* rebin_th2(TH2D & h_to_rebin, char* hname, double* x_bins, int n_x_bins,
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
                        const bool use_miss_fake=false,
                        const bool do_roounfoldresponse=true);

        //---------------------------------------------------------------
        // Remove outliers from a TH1 via "simple" method:
        //   delete any bin contents with N counts < limit
        // Modifies histogram in-place and returns its pointer
        //---------------------------------------------------------------
        TH1* simpleRemoveOutliers(TH1* hist, bool verbose, int limit);

        //---------------------------------------------------------------
        // Remove outliers from a TH1 via pT-hat method:
        //   delete any bin contents with pT > limit
        // Modifies histogram in-place and returns its pointer
        //---------------------------------------------------------------
        TH1* pThatRemoveOutliers(TH1* hist, bool verbose, const double & limit);

        //---------------------------------------------------------------
        // Remove outliers from a THn via "simple" method:
        //   delete any bin contents with N counts < limit
        // Modifies histogram in-place and returns its pointer
        //---------------------------------------------------------------
        THn* simpleRemoveOutliersTHn(THn* hist, bool verbose, int limit, int dim);

        //---------------------------------------------------------------
        // Remove outliers from a THn via pT-hat method:
        //   delete any bin contents with pT_truth > limit
        // Modifies histogram in-place and returns its pointer
        //---------------------------------------------------------------
        THn* pThatRemoveOutliers(THn* hist, bool verbose, const double & limit, int dim, int pTdim);

        //------------------------------------------------------
        // Convolution of nonperturbative shape functions

        // Create and return 2D histogram, convolving h with shape function
        TH2D* convolve_F_np(const double & Omega, const double & R, const double & beta,
                            const double* ob_bins, const int & n_ob_bins, const double* obs,
                            const double* ob_bin_width,
                            const double* pT_bins, const int & n_pT_bins, const double* pTs,
                            const TH2D & h, const std::string & name, const bool groomed = false,
                            const double & sd_beta = 0, const double & sd_zcut = 0.2,
                            const std::string & option = "");

        double find_cell(double val, const double * cell, const int range, bool phi);

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
                               const float min_det_pt=0.,
                               const float min_truth_pt=0.,
                               const float min_det=0.,
                               const float min_truth=0.,
                               const float max_det_pt=0.,
                               const float max_truth_pt=0.,
                               const float max_det=0.,
                               const float max_truth=0.,
                               const double & prior_variation_parameter=0.,
                               const bool move_underflow=false,
                               const bool use_miss_fake=false);

        // Set scaling of prior
        prior_scale_func prior_scale_factor_obs(const int & option);

        // Recursive helper function for simpleRemoveOutliersTHn()
        void simpleRemoveOutliersTHn_recurse(
            THn* hist, int limit, int dim, int* n_bins, int* x, int dim_to_update);
        // Recursive helper function for pThatRemoveOutliersTHn()
        void pThatRemoveOutliersTHn_recurse(
            THn* hist, int limit, int dim, int pTdim, int max_bin_number,
            int* n_bins, int* x, int dim_to_update, bool verbose);

        //------------------------------------------------------
        // Convolution of nonperturbative shape functions

        // Non-perturbative parameter with factored-out beta dependence
        // Omega is Omega_{a=0} == Omega_{beta=2}  (universal[?])
        inline double Omega_beta(const double & Omega, const double & beta);

        // Shape function for convolving nonperturbative effects
        inline double F_np(const double & Omega, const double & k, const double & beta);

    ClassDef(HistUtils, 1)
    };

};
#endif
