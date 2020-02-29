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
h_tru->Draw("hist same");
TH1D* h_cut = (TH1D*) f_fastsim->Get("truth_pt_eff_cutsScaled");
h_cut->Draw("hist same");

TCanvas* c2 = new TCanvas("c2", "c2");
TH1D* h_ratio = (TH1D*) h_cut->Clone();
h_ratio->Divide(h_tru);
h_ratio->Draw();


// Clean up memory
delete f_fastsim;
delete c1;
delete c2;
