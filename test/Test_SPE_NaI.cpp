#include "Analysis.h"
#include "WaveForm.h"
#include "globals.h"
#include "singleHits.h"
#include <TApplication.h>
#include <TH1.h>
#include <iostream>
int main(int argc, char *argv[]) {
#ifdef WAVES
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;

  std::string fname = "/home/kirtikesh/analysisSSD/DATA/WCu_Test/"
                      "run_Cs_Direct_NaI0_1350V_18Nov/"
                      "FILTERED/DataF_run_Cs_Direct_NaI0_1350V_18Nov.root";

  digiAnalysis::Analysis an(fname, 0, 100000, 0);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int nentries = hitsVector.size();
  std::cout << "Got the vector from an: " << nentries << std::endl;
  int evi = 0;
  std::string userInput;
  bool keepGoing = true;
  digiAnalysis::WaveForm *WF = nullptr;
  std::vector<double> tracesSmooth;
  int sbSize = 30;

  int peakDist = 100; // peak seperation in ns
  peakDist /= 2;      // to convert to samples in 500MSPS digitizer

  double peakThreshold = 20.0; // threshold for integration peak
  double threshold = 10.0;     // threshold for peak detection above baseline

  UShort_t checkStartTime =
      digiAnalysis::GateStart +
      digiAnalysis::GateLenLong; // look for SPE after gate end
  int integrateRange = 20;

  //   while (keepGoing) {
  //     hitsVector[evi]->SetSmoothWF(sbSize);
  //     WF = nullptr;
  //     WF = hitsVector[evi]->GetWFPtr();
  //     WF->Plot();

  //     std::cout << "Do you want to see the next waveform? (y/n): ";
  //     std::getline(std::cin, userInput);
  //     if (userInput != "y" && userInput != "Y") {
  //       keepGoing = false;
  //     }
  //     evi += 1;
  //   }

  TH1 *hSPE = new TH1F("hSPE", "SPE Energy", 16384, 0, 16384);
  for (evi = 0; evi < nentries; evi++) {
    hitsVector[evi]->SetSmoothWF(sbSize);
    auto pvVec = hitsVector[evi]->DetectPeakValleys(threshold);
    WF = nullptr;
    tracesSmooth = {};
    WF = hitsVector[evi]->GetWFPtr();
    tracesSmooth = WF->GetTracesSmooth();
    auto peaks = pvVec.first;
    int iter = 0;
    int prevPeakPos = 0;
    int pos = peaks[iter];
    double SPEInt = 0;
    while (iter < peaks.size() - 1) {
      iter += 1;
      prevPeakPos = pos;
      pos = peaks[iter];
      if (pos > checkStartTime and tracesSmooth[pos] > peakThreshold) {
        if (pos - prevPeakPos > peakDist and peaks[iter + 1] - pos > peakDist) {
          SPEInt = WF->IntegrateSmoothWaveForm(pos - integrateRange,
                                               pos + integrateRange);
          hSPE->Fill(SPEInt);
        }
      }
    }
  }

  TCanvas *canvas1 = new TCanvas("canvas1", "Energy Hists", 1600, 1000);
  canvas1->cd();
  hSPE->Draw("HIST");
  canvas1->Update();

  fApp->Run();
#endif
  return 0;
}
