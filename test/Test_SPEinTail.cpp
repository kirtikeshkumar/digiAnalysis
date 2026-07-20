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
#include <memory>
#include <ratio>
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

  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
  //     "CoincidenceStudies/PairFiles/"
  //     "Pair_NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp.root";

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/extCoinc/"
      "NaI3124_16Jun26_AmSrc_1350V_2000V_1350V_1350V_Gain2_NoSplit_ExtTrig_"
      "Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_1000nsCoinc_"
      "2Vpp_Thresh_100lsb_WAVES_21/FILTERED/"
      "DataF_NaI3124_16Jun26_AmSrc_1350V_2000V_1350V_1350V_Gain2_NoSplit_"
      "ExtTrig_Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_"
      "1000nsCoinc_2Vpp_Thresh_100lsb_WAVES_21_BLCorrected.root";

  digiAnalysis::Analysis an(fname, 0000, 200000, 0);
  std::cout << "getting the vector from an" << std::endl;

  // an.CreatePairs();

  // std::vector<std::unique_ptr<digiAnalysis::Pair>> &vecOfPairs =
  //     an.GetPairsVec();
  // int nentries = vecOfPairs.size();
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &vecOfHits =
      an.GetSingleHitsVec();
  int nentries = vecOfHits.size();
  std::cout << "got the vector from an: " << nentries << std::endl;
  int nPairs = nentries;

  double Energy1 = 0;
  double Energy2 = 0;
  double PSD = 0, MT = 0;
  digiAnalysis::singleHits *hit;
  digiAnalysis::WaveForm *WF = nullptr;
  std::vector<digiAnalysis::WaveForm> waveformVector;
#ifdef WAVES
  TH1 *hSPE = new TH1F("hSPE", "hSPE", 1000, 0, 5000);
  TH2 *hESPE = new TH2F("hESPE", "hESPE", 16384, 0, 16384, 1000, 0, 100);
  TH1 *hpreVal = new TH1F("hpreVal", "distance to pre valley", 2000, 0, 2000);
  TH1 *hpostVal =
      new TH1F("hpostVal", "distance to post valley", 2000, 0, 2000);
  bool keepGoing = true;
  std::string userInput;
  double intSPE, intWave;
  int wfSz;

  for (int iter = 0; iter < nPairs && keepGoing; iter++) {

    // hit = vecOfPairs[iter]->GetHitPtr(0);
    hit = vecOfHits[iter].get();

    if (iter % 10000 == 0) {
      std::cout << hit->GetEvNum() << " : " << hit->GetTimestamp() << std::endl;
    }

    // hit->GetEnergy() > 694
    //     ? Energy1 = vecOfPairs[iter]->GetPairHitEnergy(0) * 0.09465 - 5.7613
    //     : Energy1 = vecOfPairs[iter]->GetPairHitEnergy(0) * 0.08696 -
    //                 0.4222; // Calibration to get the energy
    //                         // 1900V
    Energy1 = hit->GetEnergy() * 0.01911 - 0.3;
    if (Energy1 > 55 and Energy1 < 64) {
      WF = nullptr;
      WF = hit->GetWFPtr();
      WF->SetSmooth(65);
      // WF->Plot();
      auto results = WF->DetectPeakValleys(6);
      // std::cout << "size of peaks: " << results.first.size() << std::endl;
      // std::cout << "size of valleys: " << results.second.size() << std::endl;
      // //
      // // Print the identified peaks
      // int iter1 = 0;
      // while (iter1 < results.first.size())
      // {
      //   std::cout
      //       << iter1 << ": Positions: " << results.second[2 * iter1] << ":"
      //       << results.first[iter1] << ":"
      //       << results.second[2 * iter1 + 1]
      //       << std::endl;
      //   std::cout
      //       << iter1 << ": Amplitudes: "
      //       << WF->GetTraces()[results.second[2 * iter1]]
      //       << ":" << WF->GetTraces()[results.first[iter1]] << ":"
      //       << WF->GetTraces()[results.second[2 * iter1 + 1]]
      //       << std::endl;

      //   iter1 += 1;
      //   // if (iter >= results.first.size())
      //   //     break;
      // }

      // Select isolated SPE peaks and integrate to get charge
      digiAnalysis::WaveForm WFSPE;
      int iterPeaks = 0;
      int isolationRange = 50;
      int saveRange = 250; // isolationRange - 50;
      while (iterPeaks < results.first.size() && keepGoing) {
        int peakPos = results.first[iterPeaks];
        if ((peakPos > 2000 and peakPos < 4500) and
            (peakPos - results.first[iterPeaks - 1] > isolationRange)) {
          if ((iterPeaks + 1 < results.first.size() and
               (results.first[iterPeaks + 1] - peakPos) > isolationRange) ||
              (iterPeaks + 1 == results.first.size())) {
            double postBL = WF->EvalBaseLine(peakPos + 100, 50);
            double preBL = WF->EvalBaseLine(peakPos - 100, 50);
            if (fabs(preBL - postBL) < 2.0) {
              WFSPE.SetWaveForm(*WF, peakPos - saveRange, peakPos + saveRange,
                                saveRange - 100, 50);
              // WFSPE.SetBaseLine(50, 50);
              wfSz = WFSPE.GetSize();
              WFSPE.SetSmooth(25);
              waveformVector.push_back(WFSPE);
              // std::cout << iterPeaks << ":" << peakPos << std::endl;
              intSPE = WFSPE.IntegrateWaveForm(saveRange - 50, saveRange + 100);
              intWave = WF->IntegrateWaveForm(digiAnalysis::GateStart,
                                              digiAnalysis::GateStart +
                                                  digiAnalysis::GateLenLong);
              hSPE->Fill(intSPE);
              hESPE->Fill(Energy1, intWave / intSPE / Energy1);

              // WFSPE.Plot();
              // std::cout << "Do you want to see the next waveform? (y/n): ";
              // std::getline(std::cin, userInput);
              // if (userInput != "y" && userInput != "Y") {
              //   keepGoing = false;
              // }
            }
          }
        }

        iterPeaks += 1;
      }

      // Integrate between valleys to get charge of peaks and plot against
      // energy to see if we can identify the SPE peak digiAnalysis::WaveForm
      // WFSPE;
      iterPeaks = 0;
      while (iterPeaks < results.first.size()) {
        int peakPos = results.first[iterPeaks];
        int valley1Pos = results.second[2 * iterPeaks];
        int valley2Pos = results.second[2 * iterPeaks + 1];
        if (peakPos > 2000 and peakPos < 4500) {
          // intSPE = WF->IntegrateWaveForm(valley1Pos, valley2Pos);
          // intWave = hit->GetEvalEnergy();
          // hSPE->Fill(intSPE);
          // hESPE->Fill(Energy1, intWave / intSPE / Energy1);
          hpreVal->Fill(peakPos - valley1Pos);
          hpostVal->Fill(valley2Pos - peakPos);
        }
        iterPeaks += 1;
      }
      iterPeaks = 0;

      // std::cout << std::endl
      //           << "VALLEYS:____________" << std::endl;
      // while (iterPeaks < results.second.size())
      // {
      //     std::cout << iterPeaks << ":" << results.second[iterPeaks]
      //<< ":"
      //     << traces[results.second[iterPeaks]]
      //               << std::endl;
      //     iterPeaks += 1;
      //     // if (iterPeaks >= results.second.size())
      //     //     break;
      // }

      //   WF->Plot();
      //   std::cout << "Energy is approx: " << Energy1 << std::endl;
      // std::cout << "Do you want to see the next waveform? (y/n): ";
      // std::getline(std::cin, userInput);
      // if (userInput != "y" && userInput != "Y") {
      //   keepGoing = false;
      // }
    }
  }

  std::cout << "now showing average SPE waveform (" << waveformVector.size()
            << ")" << std::endl;

  digiAnalysis::WaveForm WFAveraged(wfSz, waveformVector);
  // WFAveraged.SetSmooth(150);
  WFAveraged.SetTracesFFT("orig");
  WFAveraged.Plot();

  TCanvas *c1 = new TCanvas("c1", "SPECharge", 800, 600);
  hSPE->Draw("HIST");
  TCanvas *c2 = new TCanvas("c2", "EvsNSPE", 800, 600);
  hESPE->Draw("COLZ");
  TCanvas *c3 = new TCanvas("c3", "preVal dist", 800, 600);
  hpreVal->Draw("HIST");
  TCanvas *c4 = new TCanvas("c4", "postVal distance", 800, 600);
  hpostVal->Draw("HIST");
  fApp->Run();
#endif
}