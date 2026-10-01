// ################################################################# //
//  This code is meant to reduce the file size by removing the WF    //
//                 from all channels but highgain                    //
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
  std::string fpath = "/home/kirtikesh/Analysis/DATA/extCoincSep/1800V/Data/";
  std::string finit =
      "NaI134_NoSrc_30Sep_1800_1337_1350_WAVES_NoSplitSignal_Gain_12_"
      "Acquisition_4_ExtTrig_Threshold1LSB_160nsPromt_240nsDelay_800nsLong_"
      "1496nsCoinc_FreeWrites_Ch248Singles_";
  for (int fiter = 26; fiter <= 40; fiter++) {
    std::string fname = fpath + finit + std::to_string(fiter) +
                        "/UNFILTERED/Data_" + finit + std::to_string(fiter) +
                        "_BLCorrected.root";
    if (std::filesystem::exists(fname)) {
      digiAnalysis::Analysis an(fname, 00000, 00000, 0);
      std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
          an.GetSingleHitsVec();

      // to write the subtracted BL to file
      UShort_t Channel;
      ULong64_t Timestamp;
      UShort_t Board;
      UShort_t Energy;
      UShort_t EnergyShort;
      TArrayS *Samples = nullptr;
      std::string writefname =
          fpath + finit + std::to_string(fiter) + "/UNFILTERED/Data_" + finit +
          std::to_string(fiter) + "_BLCorrected_Compressed.root";
      TFile *fout = new TFile(writefname.c_str(), "RECREATE");
      TTree *Data_F = new TTree("Data_F", "Filtered Data");
      Data_F->Branch("Channel", &Channel, "Channel/s");
      Data_F->Branch("Timestamp", &Timestamp, "Timestamp/l");
      Data_F->Branch("Board", &Board, "Board/s");
      Data_F->Branch("Energy", &Energy, "Energy/s");
      Data_F->Branch("EnergyShort", &EnergyShort, "EnergyShort/s");
      Data_F->Branch("Samples", &Samples);

      for (int hititer = 0; hititer < hitsVector.size(); hititer++) {
        Channel = hitsVector[hititer]->GetChNum();
        Timestamp = hitsVector[hititer]->GetTimestamp();
        Board = hitsVector[hititer]->GetBoard();
        Energy = hitsVector[hititer]->GetEnergy();
        EnergyShort = hitsVector[hititer]->GetEnergyShort();
        std::vector<double> traceother =
            hitsVector[hititer]->GetWFPtr()->GetTraces();
        if (hititer == 10) {
          hitsVector[hititer]->GetWFPtr()->Plot(traceother);
        }
        switch (hitsVector[hititer]->GetChNum()) {
        case 0:
          Samples = new TArrayS(traceother.size());
          for (size_t i = 0; i < traceother.size(); i++) {
            (*Samples)[i] = static_cast<Short_t>(std::round(traceother[i]));
          }
          Data_F->Fill();
          delete Samples;
          Samples = nullptr;
          break;
        default:
          Samples = new TArrayS(1);
          for (size_t i = 0; i < 1; i++) {
            (*Samples)[i] = static_cast<Short_t>(std::round(traceother[i]));
          }
          Data_F->Fill();
          delete Samples;
          Samples = nullptr;
          break;
        }
      }
      Data_F->Write();
      fout->Close();
      delete fout;
      fout = nullptr;
    } else {
      std::cout << "ERROR: No File \"" << fname << "\"" << std::endl;
    }
  }
}