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
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;
  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //     "NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp/FILTERED/"
  //     "SDataF_NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp.root";

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
      "CoincidenceStudies/PairFiles/"
      "Pair_NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp.root";

  digiAnalysis::Analysis an(fname, 0, 10, 0);
  std::cout << "getting the vector from an" << std::endl;

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
  std::vector<digiAnalysis::WaveForm> waveformVector;
#ifdef WAVES
  TH1 *hSPE = new TH1F("hSPE", "hSPE", 1500, 0, 1500);
  TH2 *hESPE = new TH2F("hESPE", "hESPE", 80, 20, 100, 100, 0, 100);
  bool keepGoing = true;
  std::string userInput;
  double intSPE, intWave;
  int wfSz;
  for (int iter = 0; iter < nPairs && keepGoing; iter++) {
    if (iter % 1 == 0) {
      hit = vecOfPairs[iter]->GetHit(0);
      std::cout << hit->GetEvNum() << " : " << hit->GetTimestamp() << std::endl;
    }

    hit->GetEnergy() > 694
        ? Energy1 = vecOfPairs[iter]->GetPairHitEnergy(0) * 0.09465 - 5.7613
        : Energy1 = vecOfPairs[iter]->GetPairHitEnergy(0) * 0.08696 -
                    0.4222; // Calibration to get the energy
                            // 1900V
    if (Energy1 > 90 and Energy1 < 100) {
      WF = nullptr;
      WF = hit->GetWFPtr();
      WF->SetSmooth(65);
      auto results = WF->DetectPeakValleys(3);
      // std::cout << "size of peaks: " << results.first.size() << std::endl;
      // std::cout << "size of valleys: " << results.second.size() << std::endl;
      //
      // Print the identified peaks
      //   int iter1 = 0;
      //   while (iter1 < results.first.size()) {
      //     std::cout
      //         << iter1 << ":"
      //         << results.first[iter1] //<< ":" << traces[results.first[iter]]
      //         << std::endl;
      //     iter1 += 1;
      //     // if (iter >= results.first.size())
      //     //     break;
      //   }

      // Select isolated SPE peaks and integrate to get charge
      int iterPeaks = 0;
      int isolationRange = 250;
      int saveRange = isolationRange - 50;
      while (iterPeaks < results.first.size()) {
        int peakPos = results.first[iterPeaks];
        if ((peakPos > 2500 and peakPos < 4500) and
            (peakPos - results.first[iterPeaks - 1] > isolationRange)) {
          if ((iterPeaks + 1 < results.first.size() and
               (results.first[iterPeaks + 1] - peakPos) > isolationRange) ||
              (iterPeaks + 1 == results.first.size())) {
            digiAnalysis::WaveForm WFSPE;
            WFSPE.SetWaveForm(*WF, peakPos - saveRange, peakPos + saveRange,
                              saveRange - 100, 50);
            // WFSPE.SetBaseLine(50, 50);
            wfSz = WFSPE.GetSize();
            waveformVector.push_back(WFSPE);
            // std::cout << iterPeaks << ":" << peakPos << std::endl;
            intSPE = WFSPE.IntegrateWaveForm(saveRange - 50, saveRange + 100);
            intWave = WF->IntegrateWaveForm(digiAnalysis::GateStart,
                                            digiAnalysis::GateStart +
                                                digiAnalysis::GateLenLong);
            // hSPE->Fill(intSPE / 150.0 * digiAnalysis::EvalNormFactor);
            hESPE->Fill(Energy1, intWave / intSPE / Energy1);
            // WFSPE.Plot();
            // std::cout << "Do you want to see the next waveform? (y/n): ";
            // std::getline(std::cin, userInput);
            // if (userInput != "y" && userInput != "Y") {
            //   keepGoing = false;
            // }
          }
        }

        iterPeaks += 1;
      }

      // iterPeaks = 0;

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
  // WFAveraged.SetTracesFFT("smooth");
  WFAveraged.Plot();

  // TCanvas *c1 = new TCanvas("c1", "SPECharge", 800, 600);
  // hSPE->Draw("HIST");
  TCanvas *c2 = new TCanvas("c2", "EvsNSPE", 800, 600);
  hESPE->Draw("COLZ");
  fApp->Run();
#endif
}