#include "Analysis.h"
#include "Pair.h"
#include "WaveForm.h"
#include "singleHits.h"
#include <RtypesCore.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TString.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;
  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //     "Directionality/"
  //     "NaI_13_CsSrc_LinearConf_HV_1900V_1345V_50cm_Coinc_96ns_Run_CFD_WAVES/"
  //     "FILTERED/"
  //     "SDataF_NaI_13_CsSrc_LinearConf_HV_1900V_1345V_50cm_Coinc_96ns_Run_CFD_"
  //     "WAVES.root";

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
      "CoincidenceStudies/PairFiles/"
      "Pair_NaI13_12May26_1900_1345_Cs_Coinc144ns_35cm_NoCollimation_1.root";

  digiAnalysis::Analysis an(fname, 0000, 00, 1);
  std::cout << "getting the vector from an" << std::endl;

  std::string outfname =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
      "SPEFiles/"
      "SPE_Ch0_NaI13_12May26_1900_1345_Cs_Coinc144ns_35cm_NoCollimation_1.root";
  TFile *fout = TFile::Open(outfname.c_str(), "RECREATE");
  TTree *t = new TTree("SPE_WF", "SPE_WF");
  digiAnalysis::WaveForm WFSPE;
  t->Branch("SPE", "digiAnalysis::WaveForm", &WFSPE);

  // an.CreatePairs();

  std::vector<std::unique_ptr<digiAnalysis::Pair>> &vecOfPairs =
      an.GetPairsVec();
  int nentries = vecOfPairs.size();
  std::cout << "got the vector from an: " << nentries << std::endl;
  int nPairs = nentries;

  double Energy1 = 0;
  double Energy2 = 0;
  double PSD = 0, MT = 0;
  digiAnalysis::singleHits *hit;
  digiAnalysis::WaveForm *WF = nullptr;
  for (int iter = 0; iter < nPairs; iter++) {

    hit = vecOfPairs[iter]->GetHitPtr(0);

    if (iter % 10000 == 0) {
      std::cout << hit->GetEvNum() << " : " << hit->GetTimestamp() << std::endl;
    }
    WF = nullptr;
    WF = hit->GetWFPtr();
    // WF->SetSmooth(40);
    auto results = WF->DetectPeakValleys(10);
    // std::cout << iter << ": DetectedNumPeaks: " << results.first.size()
    //           << std::endl;
    int iterPeaks = 0;
    int isolationRange = 150;
    int saveRange = isolationRange;
    while (iterPeaks < results.first.size()) {
      //   std::cout << iter << " : PeakNum: " << iterPeaks << std::endl;
      int peakPos = results.first[iterPeaks];
      if ((peakPos > 2500 and peakPos < 4500) and
          (peakPos - results.first[iterPeaks - 1] > isolationRange)) {
        if ((iterPeaks + 1 < results.first.size() and
             (results.first[iterPeaks + 1] - peakPos) > 2 * isolationRange) ||
            (iterPeaks + 1 == results.first.size())) {
          double postBL = WF->EvalBaseLine(peakPos + 100, 50);
          double preBL = WF->EvalBaseLine(peakPos - 100, 50);
          if (fabs(preBL - postBL) < 2.0) {
            // WFSPE.Clear();
            WFSPE.SetWaveForm(*WF, peakPos - saveRange,
                              peakPos + 2.0 * saveRange, saveRange - 100, 50);
            // std::cout << iter << "WaveFormSet" << std::endl;
            t->Fill();
          }
        }
      }
      iterPeaks += 1;
    }
  }
  fout->Write();
  fout->Close();
}