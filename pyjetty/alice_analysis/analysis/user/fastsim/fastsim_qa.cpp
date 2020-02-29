/* ROOT macro for creating QA plots for fast simulation
 * Ezra Lesser (elesser@berkeley.edu)
 */

// Open fast simulation file
int job_number = 35288;
std::string fastsim_dir = "/rstorage/u/alice/AnalysisResults/fastsim/";
std::string fastsim_filename = fastsim_dir + std::to_string(job_number) + "/AnalysisResultsFinal.root";
TFile* f_fastsim = new TFile(fastsim_filename.c_str(), "READ");

// Open pT histograms and create ratio
TCanvas* c1 = new TCanvas("c1", "c1");
gPad->SetLogy()
TH1D* h_tru = (TH1D*) f_fastsim->Get("truth_ptScaled");
//h_tru->Draw("hist same");
h_tru->SetBinContent(0, 1) // Prevent divide by 0 issues
h_tru->SetBinContent(1, 1)
TH1D* h_cut = (TH1D*) f_fastsim->Get("truth_pt_eff_cutsScaled");
//h_cut->Draw("hist same");

//TCanvas* c2 = new TCanvas("c2", "c2");

/*
TH1D* h_ratio = new TH1D("h_ratio", "h_ratio", h_tru->GetNbinsX(), h_tru->GetXaxis()->GetXbins());
h_ratio->Copy(h_tru);
*/

TH1D* h_ratio = (TH1D*) h_tru->Clone();
h_ratio->SetName("h_ratio");
h_ratio->SetTitle("h_ratio");

// This is what we want to do, but it's not working for some reason, so do it manually
//h_ratio->Divide(h_tru);

/*
for (int i = 0; i < h_ratio->GetNbinsX(); i++) {
    double ratio = double(h_cut->GetBinContent(i)) / double(h_tru->GetBinContent(i));
    std::cout << double(h_cut->GetBinContent(i)) << " / " << double(h_tru->GetBinContent(i)) << " = "
	      << ratio << '\n';
    h_ratio->SetBinContent(i, ratio);
    }*/

//h_ratio->DrawCopy();

TRatioPlot* rp = new TRatioPlot(h_cut, h_tru);
rp->Draw();
rp->GetLowerRefGraph()->Draw();


// Clean up memory
delete f_fastsim;
delete c1;
delete c2;
delete h_ratio;
