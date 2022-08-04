// ROOT Macro for comparing ROOT THn before and after pT-hat based cut
// Written by Ezra Lesser, August 2022
#include <string>
#include "THn.h"
#include "TH2.h"

int pT_cut_comp() {

  // Open up trimmed and untrimmed ROOT files
  std::string base_dir = "/rstorage/alice/AnalysisResults/ang/877553/";
  std::string untrimmed_dir = base_dir + "Scaled_no_cut/";
  std::string trimmed_dir = base_dir + "pt_trimmed_piecewise/";
  std::string file_name = "AnalysisResults.root";

  // Jet parameters
  std::string obs = "ang";
  std::string jetR = "0.2";  // "0.4";
  std::string Rmax = "0.1";  // "0.25";
  //std::vector<std::string> alphas = {"1", "1.5", "2", "3"};
  std::vector<std::string> alphas = {"1"};

  for (int pThatBin = 1; pThatBin <= 20; pThatBin++) { // loop over pT-hat bins
    // Open ROOT files
    TFile* file_untrimmed = new TFile((untrimmed_dir + std::to_string(pThatBin) + '/' + file_name).c_str(), "r");
    TFile* file_trimmed = new TFile((trimmed_dir + std::to_string(pThatBin) + '/' + file_name).c_str(), "r");

    for (const std::string & alpha : alphas) {
      // Retrieve RM for making projections
      std::string RM_name =
        "hResponse_JetPt_" + obs + "_R" + jetR + '_' + alpha + "_Rmax" + Rmax + "Scaled";
      THn* RM_untrimmed = (THn*) file_untrimmed->Get(RM_name.c_str());
      RM_untrimmed->SetNameTitle("untrimmed", "untrimmed");
      THn* RM_trimmed = (THn*) file_trimmed->Get(RM_name.c_str());

      // Project out the x-y (pT) distributions
      TH2D* pT_untrimmed = RM_untrimmed->Projection(1, 0);
      TH2D* pT_trimmed = RM_trimmed->Projection(1, 0);

      // Create ratio trimmed/untrimmed and write to file
      TH2D* pT_ratio = (TH2D*) pT_trimmed->Clone((RM_name + "_ratio").c_str());
      pT_ratio->Divide(pT_untrimmed);

      const char canvas_name[] = "pT ratio trimmed/untrimmed";
      TCanvas c(canvas_name, canvas_name, 700, 600);
      c.Draw();

      const char pad_name[] = "the pad";
      TPad pad(pad_name, pad_name, 0, 0, 1, 1);
      pad.SetTicks(1, 1);
      pad.SetLeftMargin(0.1);
      pad.SetRightMargin(0.05);
      pad.Draw();

      pT_ratio->Draw("colz");
      gStyle->SetOptStat(0);
      c.SaveAs((RM_name + "_ratio_pThat" + std::to_string(pThatBin) + ".pdf").c_str());
      delete pT_untrimmed, pT_trimmed, pT_ratio;
    }

    file_untrimmed->Close();
    file_trimmed->Close();
    delete file_untrimmed, file_trimmed;
  }

  return 0;
}
