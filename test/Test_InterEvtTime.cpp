#include "Analysis.h"
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
#include <iostream>
#include <ratio>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;
  std::string fname =
      "/home/kirtikesh/Analysis/DATA/extCoincSep/1800V/Data/"
      "NaI134_NoSrc_WAVES_SplitSignal_Gain_6_Acquisition_4_ExtTrig_Threshold_"
      "6mV_160nsPromt_240nsDelay_800nsLong_1200nsCoinc_30SepSum."
      "root";
  digiAnalysis::Analysis an(0, fname, 0, 450000, 0);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  uint nentries = hitsVector.size();
  TH2 *hDelT_1 =
      new TH2F("hDelT_1", "Energy vs log(DelT)", 8000, 0, 80, 1300, 0, 13);
  ULong_t dt = 0;
  double e = 0;
  for (uint iter = 0; iter < nentries - 1; iter++) {
    dt =
        hitsVector[iter + 1]->GetTimestamp() - hitsVector[iter]->GetTimestamp();
    e = hitsVector[iter]->GetEvalEnergy() * 0.0062915 + 0.3047;
    hDelT_1->Fill(e, TMath::Log10(dt));
  }
  hDelT_1->Draw("COLZ");
  fApp->Run();
}