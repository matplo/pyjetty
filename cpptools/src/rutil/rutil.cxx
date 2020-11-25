#include <rutil.hh>

#include <TFile.h>

#include <iostream>

ClassImp(RUtil::Test);
ClassImp(RUtil::HistUtils);

namespace RUtil
{

    // Rebin 2D histogram h with name hname using axes given by x_bins and y_bins
    // If move_y_underflow is set, y-underflow bins in h_to_rebin are added to (x, 1) bin
    //     [WARNING: it is NOT preserved in y-underflow bin: it is really added to (x, 1)]
    // RETURNS: pointer to a new rebinned TH2F on the heap
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
    THnF* HistUtils::rebin_thn(std::string response_file_name,
                    THnF* thn,
                    const std::string & name_thn_rebinned,
                    const std::string & name_roounfold,
                    const unsigned int & n_dim,
                    const int* axes_n_bins,
                    double* axes_bin_array0,
                    double* axes_bin_array1,
                    double* axes_bin_array2,
                    double* axes_bin_array3,
                    const std::string label/*=""*/,
                    const double & prior_variation_parameter/*=0.*/,
                    int prior_option/*=1*/,
                    bool move_underflow/*=false*/,
                    bool do_roounfoldresponse/*=true*/) {
        
        // Put axis binnings into a vector
        std::vector<double*> axes_bin_arrays;
        axes_bin_arrays.push_back(axes_bin_array0);
        axes_bin_arrays.push_back(axes_bin_array1);
        axes_bin_arrays.push_back(axes_bin_array2);
        axes_bin_arrays.push_back(axes_bin_array3);
 
        // Initialize char* array of axes titles
        const char** axes_titles = new const char*[n_dim];
        for (unsigned int i = 0; i < n_dim; i++) {
            axes_titles[i] = thn->GetAxis(i)->GetTitle();
        }

        // ------------------------------------------------------
        // Create THn
        
        // Create empty THn with specified binnings
        THnF* thn_rebinned = create_empty_thn(
            name_thn_rebinned.c_str(), n_dim, axes_titles, axes_n_bins, axes_bin_arrays);
        
        // ------------------------------------------------------
        // Create RooUnfoldResponse
        RooUnfoldResponse* roounfold_response = nullptr;
        if (do_roounfoldresponse) {

            // Create empty RooUnfoldResponse with specified binning
            TH2D* hist_measured = thn_rebinned->Projection(2, 0);
            std::string name_measured = std::string("hist_measured_") + label;
            hist_measured->SetName(name_measured.c_str());
            TH2D* hist_truth = thn_rebinned->Projection(3, 1);
            std::string name_truth = std::string("hist_truth_") + label;
            hist_truth->SetName(name_truth.c_str());
            roounfold_response = new RooUnfoldResponse(
                hist_measured, hist_truth, name_roounfold.c_str(), name_roounfold.c_str());
            
            // Note: Using overflow bins doesn't work for 2D unfolding in RooUnfold
            //       -- we have to do it manually
            //roounfold_response.UseOverflow(True)
        }
        
        // Loop through THn and fill rebinned THn and RooUnfoldResponse
        fill_rebinned_thn(response_file_name, thn, thn_rebinned, n_dim,
                          do_roounfoldresponse, roounfold_response,
                          prior_variation_parameter, prior_option, move_underflow);

        // clean up memory
        delete[] axes_titles;

        return thn_rebinned;

    }  // rebin_thn


    //---------------------------------------------------------------
    // Create an empty THn according to specified binnings
    //---------------------------------------------------------------
    THnF* HistUtils::create_empty_thn(const char* name, const unsigned int & n_dim, const char** axes_titles,
                           const int* axes_n_bins, std::vector<double*> axes_bin_arrays) {

        // Initialize array of min & max values
        double* min = new double[n_dim];
        double* max = new double[n_dim];

        for (unsigned int i = 0; i < n_dim; i++) {
            min[i] = axes_bin_arrays[i][0];
            max[i] = axes_bin_arrays[i][axes_n_bins[i]];
        }

        // Create the new THn and set correct axes titles / bin edges
        THnF* h = new THnF(name, name, n_dim, axes_n_bins, min, max);
        for (unsigned int i = 0; i < n_dim; i++) {
            h->GetAxis(i)->SetTitle(axes_titles[i]);
            h->SetBinEdges(i, axes_bin_arrays[i]);
        }

        // clean up memory
        delete[] min;
        delete[] max;

        return h;

    }  // create_empty_thn


    //---------------------------------------------------------------
    // Fill thn_rebinned with data from thn
    //
    // Don't include underflow/overflow by default
    // If move_underflow = True, then fill underflow content of the observable
    // (from original THn) into first bin (of rebinned THn)
    void HistUtils::fill_rebinned_thn(std::string response_file_name, THnF* thn, THnF* thn_rebinned,
                           const unsigned int & n_dim,
                           bool do_roounfoldresponse/*=true*/,
                           RooUnfoldResponse* roounfold_response/*=nullptr*/,
                           const double & prior_variation_parameter/*=0.*/,
                           int prior_option/*=1*/,
                           bool move_underflow/*=false*/) {
        
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
        for (unsigned int bin_0 = 1; bin_0 < n_bins_0; bin_0++) {
            global_bin[0] = bin_0;
            x[0] = thn->GetAxis(0)->GetBinCenter(bin_0);

            // print helpful message while waiting
            std::cout << x[0] << " / " << n_bins_0 << '\r' << std::flush;

            for (unsigned int bin_1 = 1; bin_1 < n_bins_1; bin_1++) {
                global_bin[1] = bin_1;
                x[1] = thn->GetAxis(1)->GetBinCenter(bin_1);

                for (unsigned int bin_2 = 0; bin_2 < n_bins_2; bin_2++) {
                    global_bin[2] = bin_2;
                    x[2] = thn->GetAxis(2)->GetBinCenter(bin_2);

                    for (unsigned int bin_3 = 0; bin_3 < n_bins_3; bin_3++) {
                        global_bin[3] = bin_3;
                        x[3] = thn->GetAxis(3)->GetBinCenter(bin_3);

                        double content = thn->GetBinContent(global_bin);

                        // Impose a custom prior, if desired
                        if (std::abs(prior_variation_parameter) > 1e-5 && x[1] > 0 && x[3] > 0) {

                            // Scale number of counts according to variation of pt prior
                            content *= std::pow(x[1], prior_variation_parameter);

                            // Scale number of counts according to variation of observable prior
                            content *= prior_scale_factor_obs(
                                x[3], content, prior_variation_parameter, prior_option);

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
                        thn_rebinned->Fill(x, content);

                        // RooUnfoldResponse should be filled as (x[0], x[2], x[1], x[3])
                        // corresponding e.g. to (pt_det, obs_det, pt_true, obs_true)
                        if (do_roounfoldresponse) {
                            roounfold_response->Fill(x[0], x[2], x[1], x[3], content);
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
    double HistUtils::prior_scale_factor_obs(double obs_true, double content,
                                  double prior_variation_parameter, int option) {

        switch(option)
        {
            case 0: // power law
                return std::pow(obs_true, prior_variation_parameter);
            case 1: // linear scaling of distributions
                return (1 + prior_variation_parameter*(2*obs_true - 1));
            case 2: // sharpening/smoothing the distributions
                return std::pow(content, 1 + prior_variation_parameter);
            case 3:
                return (1 + obs_true);
            default:
                return obs_true;
        }
    }

}
