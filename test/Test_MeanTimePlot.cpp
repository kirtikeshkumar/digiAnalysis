/*
**	Filename : Test_Dummy.cpp
**	2024-12-03
**	username : kirtikeshkumar
*/

#include "Analysis.h"
#include "WaveForm.h"
#include "includes.hh"
#include "singleHits.h"
#include <TApplication.h>
#include <iostream>
int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);

  std::string fname =
      "/media/kirtikesh/KirtikeshSSD/Waves/CeBr_CsI_0pt5Vpp/"
      "Am_CsI_950V_Waves_CFD_Delay150ns_Fraction50_10lsb_0pt5Vpp_8IS_300PG_"
      "6000Gate_5696Holdoff_TConnector_1/FILTERED/"
      "DataF_CH2@V1730_167_Am_CsI_950V_Waves_CFD_Delay150ns_Fraction50_10lsb_"
      "0pt5Vpp_8IS_300PG_6000Gate_5696Holdoff_TConnector_1.root";

  // Read to singleHits
  digiAnalysis::Analysis an(fname, 0, 200000, 1);

  // Get the vector
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();

  TH2 *hMTPlot =
      new TH2F("MTPlot", "Energy vs MeanTime", 4096, 0, 4096, 500, 1.5, 8);
  int nentries = hitsVector.size();
  for (int i = 0; i < nentries; i++) {
    if (hitsVector[i]->GetPSD() > 0.0) {
      hMTPlot->Fill(hitsVector[i]->GetEnergy(), hitsVector[i]->GetMeanTime());
    }
  }
  TCanvas *c1 = new TCanvas("c1", "Energy vs MeanTime", 800, 600);
  hMTPlot->Draw("COLZ");
  fApp->Run();
  return 0;
}