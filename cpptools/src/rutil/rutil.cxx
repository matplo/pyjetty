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

	double HistUtils::find_cell(double val, const double * cell, const int range, bool phi){

        	if (phi && val > TMath::Pi())
            		val -= 2*TMath::Pi();

        	for(int e = 0 ; e < range ; e++){
            	if (val >= cell[e] and val < cell[e+1])
                	return cell[e];
        	}
        	std::cout << "value is out of bounds" << std::endl;
        	std::cout << val << " " << cell << std::endl;
        	exit(0);
    }
	
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
	// Same function, but overloaded for TH2D 
    //---------------------------------------------------------------
    TH2D* HistUtils::rebin_th2(TH2D & h_to_rebin, char* hname, double* x_bins, int n_x_bins,
                    double* y_bins, int n_y_bins, bool move_y_underflow /*= false*/) {

        // Initialize the new empty histogram
        std::string name = std::string(hname) + "_rebinned";
        TH2D* h = new TH2D(name.c_str(), name.c_str(), n_x_bins, x_bins, n_y_bins, y_bins);

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
                               const bool use_miss_fake/*=false*/,
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
        float min_det_pt = det_pt_bin_array[0];
        float min_truth_pt = truth_pt_bin_array[0];
        float min_det = det_bin_array[0];
        float min_truth = truth_bin_array[0];
        float max_det_pt = det_pt_bin_array[n_pt_bins_det];
        float max_truth_pt = truth_pt_bin_array[n_pt_bins_truth];
        float max_det = det_bin_array[n_obs_bins_det];
        float max_truth = truth_bin_array[n_obs_bins_truth];
        
        this->fill_rebinned_thn(response_file_name, thn, thn_rebinned, n_dim, f,
                                do_roounfoldresponse, roounfold_response,
                                min_det_pt, min_truth_pt, min_det, min_truth,
                                max_det_pt, max_truth_pt, max_det, max_truth, 
                                prior_variation_parameter, move_underflow, use_miss_fake);

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
        const float min_det_pt, const float min_truth_pt, const float min_det, const float min_truth,
        const float max_det_pt, const float max_truth_pt, const float max_det, const float max_truth,
        const double & prior_variation_parameter/*=0.*/, const bool move_underflow/*=false*/,
        const bool use_miss_fake/*=false*/) {

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
        for (unsigned int bin_0 = 0; bin_0 < n_bins_0+2; bin_0++) {
            global_bin[0] = bin_0;
            x[0] = thn->GetAxis(0)->GetBinCenter(bin_0);

            // print helpful message while waiting
            std::cout << bin_0 << " / " << n_bins_0 << '\r' << std::flush;

            for (unsigned int bin_1 = 0; bin_1 < n_bins_1+2; bin_1++) {
                global_bin[1] = bin_1;
                x[1] = thn->GetAxis(1)->GetBinCenter(bin_1);

                for (unsigned int bin_2 = 0; bin_2 < n_bins_2+2; bin_2++) {
                    global_bin[2] = bin_2;
                    x[2] = thn->GetAxis(2)->GetBinCenter(bin_2);

                    for (unsigned int bin_3 = 0; bin_3 < n_bins_3+2; bin_3++) {
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
                            }
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
                            
                            bool pt_in_det_range = x[0] > min_det_pt && x[0] < max_det_pt;
                            bool obs_in_det_range = x[2] > min_det && x[2] < max_det;
                            bool pt_in_true_range = x[1] > min_truth_pt && x[1] < max_truth_pt;
                            bool obs_in_true_range = x[3] > min_truth && x[3] < max_truth;

                            bool in_det_range = pt_in_det_range && obs_in_det_range;
                            bool in_true_range = pt_in_true_range && obs_in_true_range;
                            
                            // Fill if both det, true are in domain of RM
                            if (in_det_range and in_true_range) {
                                roounfold_response->FillContentError(x[0], x[2], x[1], x[3], content, error);
                            }
                        
                            if (use_miss_fake) {

                                // If input is not in det-range (this is our usual kinematic efficiency correction), Miss
                                if (!in_det_range && in_true_range) {
                                    roounfold_response->Miss(x[1], x[3], content);
                                }
                                // If truth-level is outside RM range (e.g. jet pt range is technically not [0,\infty]), Fake
                                // This is usually a negligible correction for us
                                else if (in_det_range && !in_true_range) {
                                    roounfold_response->Fake(x[0], x[2], content);
                                }
                            }
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
            case 4:
                return prior_scale_func_4;
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

    //---------------------------------------------------------------
    // Remove outliers from a TH1 via "simple" method:
    //   delete any bin contents with N counts < limit
    // Modifies histogram in-place and returns its pointer
    //---------------------------------------------------------------
    TH1* HistUtils::simpleRemoveOutliers(TH1* hist, bool verbose, int limit) {

        if (verbose) {
            std::cout << "Applying simple removal of outliers with counts < "
                      << limit << " for " << hist->GetName() << std::endl;
        }

        for (int i = 1; i <= hist->GetNcells(); i++) {
            if (hist->GetBinContent(i) < limit) {
                hist->SetBinContent(i, 0);
                hist->SetBinError(i, 0);
            }
        }

        return hist;
    }

    //---------------------------------------------------------------
    // Remove outliers from a TH1 via pT-hat method:
    //   delete any bin contents with pT > limit
    // Modifies histogram in-place and returns its pointer
    //---------------------------------------------------------------
    TH1* HistUtils::pThatRemoveOutliers(TH1* hist, bool verbose, const double & limit) {

        // First search for the pT bin to start the cut
        int max_cell_number = -1;
        // Assume that pT is on the x-axis
        for (int i = 1; i <= hist->GetNbinsX(); i++) {
            if (hist->GetXaxis()->GetBinLowEdge(i) > limit) {
                max_cell_number = hist->GetBin(i, 0, 0);
                break;
            }
        }

        // Apply cut if necessary
        if (max_cell_number >= 0) {
            if (verbose) {
                std::cout << "Applying pT-hat removal of outliers for pT > "
                          << limit << " for " << hist->GetName() << std::endl;
            }

            for (int i = max_cell_number; i <= hist->GetNcells(); i++) {
                if (hist->GetBinContent(i)) {
                    hist->SetBinContent(i, 0);
                    hist->SetBinError(i, 0);
                }
            }
        }

        return hist;
    }

    //---------------------------------------------------------------
    // Remove outliers from a THn via pT-hat method:
    //   delete any bin contents with pT_truth > limit
    // Modifies histogram in-place and returns its pointer
    //---------------------------------------------------------------
    THn* HistUtils::pThatRemoveOutliers(THn* hist, bool verbose, const double & limit, int dim, int pTdim) {

        // Safety check
        if (pTdim > dim) {
            std::cerr << "ERROR: cannot have pTdim = " << pTdim
                      << " for a THn of dim = " << dim << std::endl;
            throw pTdim;
        }

        // First search for the pT bin to start the cut
        int max_bin_number = -1;
        for (int i = 1; i <= hist->GetAxis(pTdim)->GetNbins(); i++) {
            if (hist->GetAxis(pTdim)->GetBinLowEdge(i) > limit) {
                max_bin_number = i;
                break;
            }
        }

        if (max_bin_number >= 0) {
            if (verbose) {
                std::cout << "Applying pT-hat removal of outliers for pT > "
                          << limit << " for " << hist->GetName() << std::endl;
            }

            int n_bins[dim] = { 0 };
            for (int d = 0; d < dim; d++) {
                n_bins[d] = hist->GetAxis(d)->GetNbins();
            }

            // Recursively iterate through dimensions to remove extra counts
            int x[dim] = { 0 };
            int dim_to_update = 0;
            pThatRemoveOutliersTHn_recurse(
                hist, limit, dim, pTdim, max_bin_number, n_bins, x, dim_to_update, verbose);

        }

        return hist;
    }

    //---------------------------------------------------------------
    void HistUtils::pThatRemoveOutliersTHn_recurse(
        THn* hist, int limit, int dim, int pTdim, int max_bin_number,
        int* n_bins, int* x, int dim_to_update, bool verbose) {

        // In the deepest recursion case, empty the bin and return
        if (dim_to_update == dim) {
            if (hist->GetBinContent(x)) {
                if (verbose) {
                   std::cout << "Erasing point " << x << " from " << hist->GetName()
                             << " with min pT-truth " << hist->GetAxis(pTdim)->GetBinLowEdge(x[pTdim])
                             << std::endl;
                }
                hist->SetBinContent(x, 0);
                hist->SetBinError(x, 0);
            }
            return;
        }

        // Otherwise, go to next dimension
        for (int i = 1; i <= n_bins[dim_to_update]; i++) {
            if (dim_to_update == pTdim && i < max_bin_number) {
                continue;
            }
            x[dim_to_update] = i;
            simpleRemoveOutliersTHn_recurse(hist, limit, dim, n_bins, x, dim_to_update+1);
        }

        return;
    }

    //---------------------------------------------------------------
    // Remove outliers from a THn via "simple" method:
    //   delete any bin contents with N counts < limit
    // Modifies histogram in-place and returns its pointer
    //---------------------------------------------------------------
    THn* HistUtils::simpleRemoveOutliersTHn(THn* hist, bool verbose, int limit, int dim) {

        if (verbose) {
            std::cout << "Applying simple removal of outliers with counts < "
                      << limit << " for " << hist->GetName() << std::endl;
        }

        int n_bins[dim] = { 0 };
        for (int d = 0; d < dim; d++) {
            n_bins[d] = hist->GetAxis(d)->GetNbins();
        }

        // Recursively iterate through dimensions to remove extra counts
        int x[dim] = { 0 };
        int dim_to_update = 0;
        simpleRemoveOutliersTHn_recurse(hist, limit, dim, n_bins, x, dim_to_update);

        return hist;
    }

    //---------------------------------------------------------------
    void HistUtils::simpleRemoveOutliersTHn_recurse(
        THn* hist, int limit, int dim, int* n_bins, int* x, int dim_to_update) {

        // In the deepest recursion case, check bin and return
        if (dim_to_update == dim) {
            if (hist->GetBinContent(x) < limit) {
                hist->SetBinContent(x, 0);
                hist->SetBinError(x, 0);
            }
            return;
        }

        // Otherwise, go to next dimension
        for (int i = 1; i <= n_bins[dim_to_update]; i++) {
            x[dim_to_update] = i;
            simpleRemoveOutliersTHn_recurse(hist, limit, dim, n_bins, x, dim_to_update+1);
        }

        return;
    }

    //---------------------------------------------------------------
    // Create and return 2D histogram, convolving h with shape function
    // ob & pT bins are identical in both old & new histograms
    // obs & pTs are arrays of the central bin values
    TH2D* HistUtils::convolve_F_np(const double & Omega, const double & R, const double & beta,
								   const double* ob_bins, const int & n_ob_bins, const double* obs,
								   const double* ob_bin_width,
								   const double* pT_bins, const int & n_pT_bins, const double* pTs,
								   const TH2D & h, const std::string & name,
								   const bool groomed/*=false*/, const double & sd_beta/*=0*/,
								   const double & sd_zcut/*=0.2*/, 
								   const std::string & option/*=""*/) {

		if (groomed && sd_beta != 0) {
			printf("ERROR: currently only implemented for sd_beta == 0\n");
			throw 1;
		}

		// Initialize NP-convolved histogram to have the same binnings as h
		TH2D* h_np = (TH2D*) h.Clone();
		h_np->SetNameTitle(name.c_str(), name.c_str());

        // Numerical integration
		// Loop over all x-bins in the final histogram
		for(int pT_i = 0; pT_i < n_pT_bins; pT_i++) {
			double pT = pTs[pT_i];
			// Loop over all y-bins in the final histogram
			for(int ob_np_i = 0; ob_np_i < n_ob_bins; ob_np_i++) {
				double ob_np = obs[ob_np_i];

				// Calculate the integral for this bin
				double integral = 0;
				int ob_i = 0;
				double ob; double k; double dk;
				do {
					ob = obs[ob_i];
					if (!groomed) { 
						k = (ob_np - ob) * pT * R;
						// Use chain rule to find dk in terms of dob and dpT
						dk = pT * R;
						if (option == "width") { dk *= ob_bin_width[ob_i]; }
					} else {  // groomed case (assumes sd_beta == 0 for simplification)
						k = pT * R * std::pow((ob_np - ob) * std::pow(sd_zcut, beta - 1), 1. / beta);
						dk = pT * R * std::pow(ob, (1 - beta) / beta) *
							std::pow(sd_zcut, beta - 1);
						if (option == "width") { dk *= ob_bin_width[ob_i]; }
					}

					integral += dk * HistUtils::F_np(Omega, k, beta) * 
						h.GetBinContent(pT_i+1, ob_i+1);

					ob_i++;
				} while (ob_i <= ob_np_i);

				if (option == "width") { h_np->SetBinContent(pT_i+1, ob_np_i+1, integral); }
				else { h_np->SetBinContent(pT_i+1, ob_np_i+1, integral * ob_bin_width[ob_np_i]); }
				h_np->SetBinError(pT_i+1, ob_np_i+1, 0);

			}  // for(ob)
		}  // for(pT)

		return h_np;

    }  // FnpUtils::convolve_F_np


    // Non-perturbative parameter with factored-out beta dependence
    // Omega is Omega_{a=0} == Omega_{beta=2}  (universal[?])
    inline double HistUtils::Omega_beta(const double & Omega, const double & beta) {
        return (Omega / (beta - 1.));
	}  // FnpUtils::Omega_beta


    // Shape function for convolving nonperturbative effects
    inline double HistUtils::F_np(const double & Omega, const double & k, const double & beta) {
        double sb = HistUtils::Omega_beta(Omega, beta);
		return (4. * k) / (sb * sb) * std::exp(-2. * k / sb);
    }  // FnpUtils::F_np


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

}  // namespace RUtil
