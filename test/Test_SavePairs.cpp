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
#include <iostream>
#include <ratio>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {
  std::cout << "hello DigiAnalysis..." << std::endl;
  std::string fname =
      "/home/kirtikesh/Analysis/DATA/"
      "extCoinc/"
      "NaI3124_15-17Jun-13Jul26_NoSrc_1350V_2000V_1350V_1350V_Gain2_NoSplit_"
      "ExtTrig_"
      "Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_1000nsCoinc_"
      "2Vpp_Thresh_100lsb_WAVES_Sum_BLCorrected.root";

  digiAnalysis::Analysis an(fname, 1000000, 2000000, 0);
  std::cout << "getting the vector from an" << std::endl;

  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int nentries = hitsVector.size();
  std::cout << "got the vector from an: " << nentries << std::endl;

  // create pairs to identitfy true photon events
  an.CreatePairs();
  std::vector<std::unique_ptr<digiAnalysis::Pair>> &vecOfPairs =
      an.GetPairsVec();
  int nPairs = vecOfPairs.size();
  std::cout << nPairs << " Pairs were formed in the data." << std::endl;

  std::string outfname =
      "/home/kirtikesh/Analysis/DATA/extCoinc/PairFiles/"
      "Pair_NaI3124_15-17Jun-13Jul26_NoSrc_1350V_2000V_1350V_1350V_Gain2_"
      "NoSplit_"
      "ExtTrig_"
      "Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_1000nsCoinc_"
      "2Vpp_Thresh_100lsb_WAVES_Sum_BLCorrected_1000000_2000000.root";

  TFile *fout = TFile::Open(outfname.c_str(), "RECREATE");
  TTree *t = new TTree("Data_Pair", "Data_Pair");

  digiAnalysis::Pair pairObj;
  t->Branch("pair", "digiAnalysis::Pair",
            &pairObj); // or similar depending on your class setup

  for (const std::unique_ptr<digiAnalysis::Pair> &p : vecOfPairs) {
    pairObj.ClearPair();
    pairObj.SetPair(*p->GetHitPtr(0), *p->GetHitPtr(1));
    t->Fill();
  }

  fout->Write();
  fout->Close();
}