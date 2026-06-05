#include "Analysis.h"
#include "TMath.h"
#include "WaveForm.h"
#include "globals.h"
#include "includes.hh"
#include "singleHits.h"
#include <TApplication.h>
#include <iostream>
#include <string>
#include <vector>
int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);

  std::string fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
                      "CoincidenceStudies/01JuneNoSrc/"
                      "NaI1342_June26_1750_1345_1350_1350_NoSrc_Thresh_15-30_"
                      "300_WAVES_Coinc_144ns_LeadPit_Sum_BLCorrected.root";

  // std::string fname =
  // "/media/kirtikesh/UbuntuFiles/SDataF_NaI1342_02June26_1750_1345_1350_1350_NoSrc_Thresh_30_300_WAVES_Coinc_144ns_LeadPit_1.root";

  // Read to singleHits
  digiAnalysis::Analysis an(fname, 0, 0, 0);
  an.CreatePairs();
  std::vector<std::unique_ptr<digiAnalysis::Pair>> &vecOfPairs =
      an.GetPairsVec();
  std::cout << "Number of pairs found: " << vecOfPairs.size() << std::endl;

  double energy0 = 0., energyOther = 0.;
  TH2 *hE1E2 = new TH2I("hE1E2", "hE1E2", 500, 0, 500, 2600, 0, 2600);

  for (int i = 0; i < vecOfPairs.size(); i++) {
    if (vecOfPairs[i]->GetPairHitCh(0) == 0) {
      energy0 = vecOfPairs[i]->GetPairHitEnergy(0) * 0.09032; //- 3.385;
      switch (vecOfPairs[i]->GetPairHitCh(1)) {
      case 2:
        energyOther = vecOfPairs[i]->GetPairHitEnergy(1) * 0.49444 - 14.93;
        break;
      case 3:
        energyOther = vecOfPairs[i]->GetPairHitEnergy(1) * 0.47714 + 6.75;
        break;
      case 4:
        energyOther = vecOfPairs[i]->GetPairHitEnergy(1) * 0.51929 - 9.29;
        break;

      default:
        break;
      }
      if (vecOfPairs[i]->GetHitPtr(0)->GetMeanTime() > 2.08)
        hE1E2->Fill(energy0, energyOther);
    }
  }
  TCanvas *c1 = new TCanvas("c1", "E1E2", 800, 600);
  c1->SetFrameFillColor(kBlack);
  hE1E2->Draw("COLZ");
  c1->Update();

  fApp->Run();
  return 0;
}