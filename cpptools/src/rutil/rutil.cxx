#include <rutil.hh>

#include <TFile.h>

#include <iostream>

ClassImp(RUtil::Test);
ClassImp(RUtil::HistUtils);

namespace RUtil
{
    typedef double (*prior_scale_func)(const double & obs_true,
                                       const double & content,
                                       const double & prior_variation_parameter);

    //---------------------------------------------------------------
    // Rebin 2D histogram h with name hname using axes given by x_bins and y_bins
    // If move_y_underflow is set, y-underflow bins in h_to_rebin are added to (x, 1) bin
    //     [WARNING: it is NOT preserved in y-underflow bin: it is really added to (x, 1)]
    // RETURNS: pointer to a new rebinned TH2F on the heap
    //---------------------------------------------------------------
    TH2F* HistUtils::rebin_th2(TH2F & h_to_rebin, char* hname, double* x_bins, int n_x_bins,
                    double* y_bins, int n_y_bins, bool move_y_underflow /*= false*/) {

        // Initialize the new empty histogram
        std::string name = std::string(hname) + "_rebinned";
        TH2F* h = new TH2F(name.c_str(), name.c_str(), n_x_bins, x_bins, n_y_bins, y_bins);

        /*  Probably don't need this check
        // Check whether sumw2 has been previously set, just for our info
        if (h->GetSumw2() == 0) {
            std::cout << "rebin_th2() -- sumw2 has not been set" << std::endl;
        } else {
            std::cout << "rebin_th2() -- sumw2 has been set" << std::endl;
        }  // if (h->GetSumw2() == 0)
        */

        // Loop over all bins and fill rebinned histogram
        for (unsigned int bin_x = 1; bin_x <= h_to_rebin.GetNbinsX(); bin_x++) {
            for (unsigned int bin_y = 0; bin_y <= h_to_rebin.GetNbinsY(); bin_y++) {

                // If underflow bin of observable, and if use_underflow is activated,
                //   put the contents of the underflow bin into the first bin of the rebinned TH2
                if (bin_y == 0) {
                    if (move_y_underflow) {
                        h->Fill(h_to_rebin.GetXaxis()->GetBinCenter(bin_x),
                                h->GetYaxis()->GetBinCenter(1),
                                h_to_rebin.GetBinContent(bin_x, bin_y));
                    }
                    continue;
                }  // biny == 0

                h->Fill(h_to_rebin.GetXaxis()->GetBinCenter(bin_x),
                        h_to_rebin.GetYaxis()->GetBinCenter(bin_y),
                        h_to_rebin.GetBinContent(bin_x, bin_y));

            }  // for(biny)
        }  // for(binx)

        // We need to manually set the uncertainties, since sumw2 does the wrong thing in this case
        // Specifically: We fill the rebinned histo from several separate weighted fills, so sumw2
        // gives sqrt(a^2+b^2) where we simply want counting uncertainties of sqrt(a+b).
        for (unsigned int i = 0; i <= h->GetNcells(); i++) {
            h->SetBinError(i, std::sqrt(h->GetBinContent(i)));
        }

        // We want to make sure sumw2 is set after rebinning, since
        // we will scale etc. this histogram
        if (h->GetSumw2() == 0) { h->Sumw2(); }

        return h;

    }  // rebin_th2

    //---------------------------------------------------------------
    // Rebin THn according to specified binnings; return pointer to rebinned THn
    //---------------------------------------------------------------
    THnF* HistUtils::rebin_thn(const std::string & response_file_name,
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
                               const std::string & label/*=""*/,
                               const double & prior_variation_parameter/*=0.*/,
                               const int & prior_option/*=1*/,
                               const bool move_underflow/*=false*/,
                               const bool do_roounfoldresponse/*=true*/) {
 
        // ------------------------------------------------------
        // Create empty THn with specified binnings
        THnF* thn_rebinned = this->create_empty_thn(name_thn_rebinned.c_str(), n_dim,
                                                    n_pt_bins_det, det_pt_bin_array,
                                                    n_obs_bins_det, det_bin_array,
                                                    n_pt_bins_truth, truth_pt_bin_array,
                                                    n_obs_bins_truth, truth_bin_array);
        
        for (unsigned int i = 0; i < n_dim; i++) {
            thn_rebinned->GetAxis(i)->SetTitle(thn->GetAxis(i)->GetTitle());
        }
        
        // ------------------------------------------------------
        // Create RooUnfoldResponse
        RooUnfoldResponse* roounfold_response = (
            do_roounfoldresponse ? 
            this->create_empty_roounfoldresponse(thn_rebinned, name_roounfold, label) : nullptr);

        const prior_scale_func f = this->prior_scale_factor_obs(prior_option);
        
        // Loop through THn and fill rebinned THn and RooUnfoldResponse
        this->fill_rebinned_thn(response_file_name, thn, thn_rebinned, n_dim, f,
                                do_roounfoldresponse, roounfold_response,
                                prior_variation_parameter, move_underflow);

        return thn_rebinned;

    }  // rebin_thn


    //---------------------------------------------------------------
    // Create an empty THn according to specified binnings
    //---------------------------------------------------------------
    THnF* HistUtils::create_empty_thn(const char* name, const int & n_dim,
                                      const int & n_pt_bins_det, const double* det_pt_bin_array,
                                      const int & n_obs_bins_det, const double* det_bin_array,
                                      const int & n_pt_bins_truth, const double* truth_pt_bin_array,
                                      const int & n_obs_bins_truth, const double* truth_bin_array) {
 
        // Obviously only working for n_dim == 4 at the moment
        if (n_dim != 4) {
            std::cerr << "ERROR: Not Implemented: Assertion n_dim == 4 failed in "
                      << "fjtools.cxx::create_empty_thn()" << std::endl;
            std::terminate();
        }

        // Initialize array of min & max values
        int nbins[4] = { n_pt_bins_det, n_pt_bins_truth, n_obs_bins_det, n_obs_bins_truth };
        double min[4] = { det_pt_bin_array[0], truth_pt_bin_array[0], det_bin_array[0], truth_bin_array[0] };
        double max[4] = { det_pt_bin_array[n_pt_bins_det], truth_pt_bin_array[n_pt_bins_truth],
                          det_bin_array[n_obs_bins_det], truth_bin_array[n_obs_bins_truth] };

        /*  Use this in arbitrary dimension case
        double** bin_edges = new double*[n_dim];
        bin_edges[0] = det_pt_bin_array;
        bin_edges[1] = truth_pt_bin_array;
        bin_edges[2] = det_bin_array;
        bin_edges[3] = truth_bin_array;
        */

        // Create the new THn and set correct axes titles / bin edges
        THnF* h = new THnF(name, name, n_dim, nbins, min, max);
        h->SetBinEdges(0, det_pt_bin_array);
        h->SetBinEdges(1, truth_pt_bin_array);
        h->SetBinEdges(2, det_bin_array);
        h->SetBinEdges(3, truth_bin_array);

        /*  For arbitrary dimensions
        for (unsigned int i = 0; i < n_dim; i++) {
            h->SetBinEdges(i, bin_edges[i]);
        }
        */

        // Clean up memory
        //delete[] bin_edges;

        return h;

    }  // create_empty_thn


    //---------------------------------------------------------------
    // Create empty RooUnfoldResponse object and return pointer
    RooUnfoldResponse* HistUtils::create_empty_roounfoldresponse(
        THnF* thn_rebinned, const std::string & name_roounfold, const std::string & label) {

        // Create empty RooUnfoldResponse with specified binning
        TH2D* hist_measured = thn_rebinned->Projection(2, 0);
        std::string name_measured = std::string("hist_measured_") + label;
        hist_measured->SetName(name_measured.c_str());
        TH2D* hist_truth = thn_rebinned->Projection(3, 1);
        std::string name_truth = std::string("hist_truth_") + label;
        hist_truth->SetName(name_truth.c_str());
        RooUnfoldResponse* roounfold_response = new RooUnfoldResponse(
            hist_measured, hist_truth, name_roounfold.c_str(), name_roounfold.c_str());

        // Note: Using overflow bins doesn't work for 2D unfolding in RooUnfold
        //       -- we have to do it manually
        //roounfold_response.UseOverflow(True)

        return roounfold_response;

    }  // create_empty_roounfold_response


    //---------------------------------------------------------------
    // Fill thn_rebinned with data from thn
    //
    // Don't include underflow/overflow by default
    // If move_underflow = True, then fill underflow content of the observable
    // (from original THn) into first bin (of rebinned THn)
    void HistUtils::fill_rebinned_thn(
        const std::string & response_file_name, const THnF* thn, THnF* thn_rebinned, const unsigned int & n_dim,
        const prior_scale_func prior_scale_f, const bool do_roounfoldresponse/*=true*/,
        RooUnfoldResponse* roounfold_response/*=nullptr*/,
        const double & prior_variation_parameter/*=0.*/, const bool move_underflow/*=false*/) {

        // Only working for n_dim == 4 at the moment; generalizing to N dimensions
        // will require some sort of recursive implementation
        if (n_dim != 4) {
            std::cerr << "ERROR: Not Implemented: Assertion n_dim == 4 failed in "
                      << "fjtools.cxx::fill_rebinned_thn()" << std::endl;
            std::terminate();
        }

        // loop through all axes
        const unsigned int n_bins_0 = thn->GetAxis(0)->GetNbins();
        const unsigned int n_bins_1 = thn->GetAxis(1)->GetNbins();
        const unsigned int n_bins_2 = thn->GetAxis(2)->GetNbins();
        const unsigned int n_bins_3 = thn->GetAxis(3)->GetNbins();

        // I don't find any global bin index implementation, so I manually loop through axes
        int* global_bin = new int[n_dim];
        double* x = new double[n_dim];
        for (unsigned int bin_0 = 1; bin_0 < n_bins_0+1; bin_0++) {
            global_bin[0] = bin_0;
            x[0] = thn->GetAxis(0)->GetBinCenter(bin_0);

            // print helpful message while waiting
            std::cout << bin_0 << " / " << n_bins_0 << '\r' << std::flush;

            for (unsigned int bin_1 = 1; bin_1 < n_bins_1+1; bin_1++) {
                global_bin[1] = bin_1;
                x[1] = thn->GetAxis(1)->GetBinCenter(bin_1);

                for (unsigned int bin_2 = 0; bin_2 < n_bins_2+1; bin_2++) {
                    global_bin[2] = bin_2;
                    x[2] = thn->GetAxis(2)->GetBinCenter(bin_2);

                    for (unsigned int bin_3 = 0; bin_3 < n_bins_3+1; bin_3++) {
                        global_bin[3] = bin_3;
                        x[3] = thn->GetAxis(3)->GetBinCenter(bin_3);

                        int bin = thn->GetBin(global_bin);
                        double content = thn->GetBinContent(bin);
                        if (content == 0) { continue; }
                        double error = thn->GetBinError(bin);

                        // Impose a custom prior, if desired
                        if (std::abs(prior_variation_parameter) > 1e-5 && x[1] > 0 && x[3] > 0) {

                            // Scale number of counts according to variation of pt & observable prior
                            double scale_factor = std::pow(x[1], prior_variation_parameter) *
                                (*prior_scale_f)(x[3], content, prior_variation_parameter);

                            content *= scale_factor;
                            error *= scale_factor;

                        }  // scale prior

                        // If underflow bin, and if move_underflow flag is activated,
                        // put the contents of the underflow bin into first bin of thn_rebinned
                        if (bin_2 == 0 || bin_3 == 0) {
                            if (move_underflow) {
                                if (bin_2 == 0) {
                                    x[2] = thn_rebinned->GetAxis(2)->GetBinCenter(1);
                                }  // bin_2 == 0
                                if (bin_3 == 0) {
                                    x[3] = thn_rebinned->GetAxis(3)->GetBinCenter(1);
                                    std::string name(thn->GetName());
                                    if (name.find("matched") != std::string::npos) {
                                        if (bin_2 == 0) { content = 1; }
                                        else { content = 0; }
                                    }
                                }  // bin_3 == 0
                            } else { continue; }  // move_underflow
                        }  // underflow bins

                        // THn is filled as (x[0], x[1], x[2], x[3])
                        // corresponding e.g. to (pt_det, pt_true, obs_det, obs_true)
                        bin = thn_rebinned->GetBin(x);
                        double prev_content = thn_rebinned->GetBinContent(bin);
                        double prev_error2 = thn_rebinned->GetBinError2(bin);
                        thn_rebinned->SetBinContent(bin, prev_content + content);
                        thn_rebinned->SetBinError(bin, std::sqrt(prev_error2 + std::pow(error, 2)));

                        // RooUnfoldResponse should be filled as (x[0], x[2], x[1], x[3])
                        // corresponding e.g. to (pt_det, obs_det, pt_true, obs_true)
                        if (do_roounfoldresponse) {
                            roounfold_response->FillContentError(
                                x[0], x[2], x[1], x[3], content, error);
                        }

                    }  // bin_3 loop
                }  // bin_2 loop
            }  // bin_1 loop
        }  // bin_0 loop

        std::cout << "writing response..." << std::endl;
        TFile f(response_file_name.c_str(), "UPDATE");
        thn_rebinned->Write();
        if (do_roounfoldresponse) { roounfold_response->Write(); }
        f.Close();
        std::cout << "done" << std::endl;

        // clean up memory
        delete[] global_bin;
        delete[] x;

        return;

    }  // fill_rebinned_thn


    //---------------------------------------------------------------
    // Compute scale factor to vary prior of observable
    //
    // Note that prior_variation_parameter is the parameter used to scale both
    // the pt prior (assumed to be a change to the power law exponent) and the observable prior,
    // and is typically taken to be +/- 0.5.
    //---------------------------------------------------------------
    prior_scale_func HistUtils::prior_scale_factor_obs(const int & option) {

        switch(option) {
            case 0: // power law
                return prior_scale_func_0;
            case 1: // linear scaling of distributions
                return prior_scale_func_1;
            case 2: // sharpening/smoothing the distributions
                return prior_scale_func_2;
            case 3:
                return prior_scale_func_3;
            default:
                return prior_scale_func_def;
        }
    }  // prior_scale_factor_obs

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

}  // namespace RUtil
