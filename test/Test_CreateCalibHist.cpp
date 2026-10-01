// ################################################################# //
//  This code is meant to create a Histogram from the BLCorrected    //
//                   Waveforms for Calibration                       //
// ################################################################# //
#include "Analysis.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1.h>
#include <filesystem>
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::string fpath = "/home/kirtikesh/Analysis/DATA/extCoincSep/1800V/Calib/";
  std::string finit =
      "NaI1_AmSrc_1800_1337_1350_WAVES_FILTERED_NoSplitSignal_Gain_12_"
      "Acquisition_4_ExtTrig_Threshold1LSB_160nsPromt_240nsDelay_800nsLong_"
      "1496nsCoinc_FreeWrites_Ch248Singles_";
  for (int fiter = 14; fiter <= 15; fiter++) {
    std::string fname = fpath + finit + std::to_string(fiter) +
                        "/UNFILTERED/Data_" + finit + std::to_string(fiter) +
                        "_BLCorrected.root";
    if (std::filesystem::exists(fname)) {
      digiAnalysis::Analysis an(fname, 00000, 00000, 0);
      std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
          an.GetSingleHitsVec();
      std::string writefname = fpath + finit + std::to_string(fiter) +
                               "/UNFILTERED/CalibHist_HighGainCh0_" + finit +
                               std::to_string(fiter);
      TFile *fout = new TFile(writefname.c_str(), "RECREATE");
      TH1 *hECh0 =
          new TH1F("hECh0", "Ch0 BLCorrected Energies", 16384, 0, 16384);
      TH1 *hECh1 = new TH1F("hECh1", "Ch4 Energies", 16384, 0, 16384);
      TH1 *hECh2 = new TH1F("hECh2", "Ch8 Energies", 16384, 0, 16384);
      for (int hititer = 0; hititer < hitsVector.size(); hititer++) {
        switch (hitsVector[hititer]->GetChNum()) {
        case 0:
          hECh0->Fill(hitsVector[hititer]->GetEvalEnergy());
        case 4:
          hECh1->Fill(hitsVector[hititer]->GetEnergy());
        case 8:
          hECh2->Fill(hitsVector[hititer]->GetEnergy());
        }
      }
      hECh0->Write();
      hECh1->Write();
      hECh2->Write();
      fout->Close();
      delete fout;
      fout = nullptr;
    } else {
      std::cout << "ERROR: No File \"" << fname << "\"" << std::endl;
    }
  }
}