#include "Analysis.h"
#include "Pair.h"
#include "WaveForm.h"
#include "globals.h"
#include "includes.hh"
#include "singleHits.h"
#include <RtypesCore.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TString.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <ratio>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {
#ifdef WAVES
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
      "NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp/FILTERED/"
      "SDataF_NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp.root";

  digiAnalysis::Analysis an(0, fname, 0, 00000, 20);
  std::cout << "getting the vector from an" << std::endl;

  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int nentries = hitsVector.size();
  std::cout << "got the vector from an: " << nentries << std::endl;

  bool keepGoing = true;
  std::string userInput;
  TH2 *hDecTime = new TH2F("hDecTime", "DecayTime vs Energy", 16384, 0, 16384,
                           100, 100, 300);

  TProfile *havDecTime = new TProfile("havDecTime", "Average decay time", 356,
                                      0, 16376); // for 4 keV bins

  TH2 *hEMT = new TH2F("hEMT", "hEMT", 4096, 0, 16384, 1000, 0, 4);

  digiAnalysis::WaveForm *wfPtr = nullptr;
  double meanTime;
  for (int iter = 0; iter < nentries && keepGoing; iter++) {
    if (iter % 1000 == 0)
      std::cout << iter << " : " << hitsVector[iter]->GetEnergy() << std::endl;
    meanTime = hitsVector[iter]->GetMeanTime();
    // std::cout << std::endl << iter << std::endl;
    hEMT->Fill(hitsVector[iter]->GetEnergy(), hitsVector[iter]->GetMeanTime());
    if (fabs(meanTime - 3.05) < 0.05) {
      wfPtr = hitsVector[iter]->GetWFPtr();
      wfPtr->SetSmooth(100);
      wfPtr->FitExponential(1, 1050, 3400);
      if (wfPtr->IsFit()) {
        // std::cout << iter << std::endl;
        hDecTime->Fill(hitsVector[iter]->GetEnergy(), wfPtr->GetFitPar(1) * 2);
        havDecTime->Fill(hitsVector[iter]->GetEnergy(),
                         wfPtr->GetFitPar(1) * 2);
        // std::cout << iter << std::endl;
      }

      // std::cout << iter << std::endl;
      // if (abs(hitsVector[iter]->GetEnergy() - 3514) < 3) {
      //   wfPtr->Plot();
      //   std::cout << "par val: " << wfPtr->GetFitPar(1) * 2
      //             << " meanTime: " << hitsVector[iter]->GetMeanTime()
      //             << std::endl;
      //   std::cout << "Do you want to see the next waveform? (y/n): ";
      //   std::getline(std::cin, userInput);
      //   if (userInput != "y" && userInput != "Y") {
      //     keepGoing = false;
      //   }
      // }
    }
  }

  TCanvas *c1 = new TCanvas("c1", "DecayTime vs Energy", 800, 600);
  hDecTime->Draw("COLZ");
  TCanvas *c2 = new TCanvas("c2", "Average DecayTime vs Energy", 800, 600);
  havDecTime->Draw();
  TCanvas *c3 = new TCanvas("c3", "Average DecayTime vs MeanTime", 800, 600);
  hEMT->Draw("COLZ");
  fApp->Run();
#else
  std::cout << "Build with WAVES set to ON" << std::endl;
#endif
}