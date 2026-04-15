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
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
      "CoincidenceStudies/PairFiles/"
      "Pair_NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp.root";

  digiAnalysis::Analysis an(fname, 0000, 000, 0);
  std::cout << "getting the vector from an" << std::endl;

  std::vector<std::unique_ptr<digiAnalysis::Pair>> &vecOfPairs =
      an.GetPairsVec();
  int nentries = vecOfPairs.size();
  std::cout << "got the vector from an: " << nentries << std::endl;
  int nPairs = nentries;

  ULong64_t tNear, tFar;
  double Energy1, Energy2;
  Int_t delT;
  TH1 *hDelT = new TH1F("hDelT", "hDelT", 100000, -100, 100);
  TH2 *hDTE = new TH2F("hDTF", "hDTF", 2000, -100, 100, 1000, 0, 1000);
  for (int iter = 0; iter < nentries; iter++) {
    tNear = vecOfPairs[iter]->GetPairHitTime(0);
    tFar = vecOfPairs[iter]->GetPairHitTime(1);
    delT = tFar - tNear;
    vecOfPairs[iter]->GetPairHitEnergy(0) > 694
        ? Energy1 = vecOfPairs[iter]->GetPairHitEnergy(0) * 0.09465 - 5.7613
        : Energy1 = vecOfPairs[iter]->GetPairHitEnergy(0) * 0.08696 -
                    0.4222; // Calibration to get the energy
                            // 1900V
    Energy2 = 0.98215 * vecOfPairs[iter]->GetPairHitEnergy(1) - 15.47;
    hDelT->Fill(delT / 1E3);
    hDTE->Fill(delT / 1E3, Energy1 + Energy2);
  }
  TCanvas *c1 = new TCanvas("c1", "tFar-tNear", 800, 600);
  hDelT->Draw("HIST");
  TCanvas *c2 = new TCanvas("c2", "DelTvsE", 800, 600);
  hDTE->Draw("COLZ");
  fApp->Run();
}