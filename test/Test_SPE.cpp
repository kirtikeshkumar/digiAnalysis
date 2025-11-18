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

  std::string fname = "/home/kirtikesh/analysisSSD/DATA/SPE/"
                      "run_Nov17_1GSPS_Direct_FreeWrites_CFD_3lsb_Delay_20ns_"
                      "Dark_Fiber_Bias1800V/FILTERED/"
                      "DataF_run_Nov17_1GSPS_Direct_FreeWrites_CFD_3lsb_Delay_"
                      "20ns_Dark_Fiber_Bias1800V.root";

  // std::string fname = "/home/kirtikesh/analysisSSD/DATA/SPE/"
  //                     "run_Nov11_1GSPS_Direct_FreeWrites_CFD_3lsb_Delay_20ns_"
  //                     "BlueLED_PicoOFF_Fiber_Bias1800V/FILTERED/"
  //                     "DataF_run_Nov11_1GSPS_Direct_FreeWrites_CFD_3lsb_Delay_"
  //                     "20ns_BlueLED_PicoOFF_Fiber_Bias1800V.root";

  // std::string fname = "/home/kirtikesh/analysisSSD/DATA/SPE/"
  //                     "run_Nov07_Direct_FreeWrites_CFD_2lsb_Delay_20ns_BlueLED_"
  //                     "Pico25_Fiber_Bias1800V/FILTERED/"
  //                     "DataF_run_Nov07_Direct_FreeWrites_CFD_2lsb_Delay_20ns_"
  //                     "BlueLED_Pico25_Fiber_Bias1800V.root";

  digiAnalysis::Analysis an(fname, 0, 0, 0);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int nentries = hitsVector.size();
  std::cout << "Got the vector from an: " << nentries << std::endl;
  int evi = 0;
  std::string userInput;
  bool keepGoing = true;
  digiAnalysis::WaveForm *WF = nullptr;
  std::vector<digiAnalysis::WaveForm> waveformVector;
  double peakThreshold = 3.5;
  double peakFWHM = 3.; // From PMT datasheet
  int sbSize = static_cast<int>(
      std::ceil(10. * peakFWHM / 2.35)); // assuming gaussian peak

  // Waveform integration of single PE peaks
  TH1 *hSPE = new TH1F("hSPE", "SPE Energy", 16384, 0, 16384);
  std::cout << "Integrating single photoelectrons from waveform" << std::endl;
  for (evi = 0; evi < nentries && keepGoing; ++evi) {
    if (evi % 100000 == 0) {
      std::cout << "Analyzing Event: " << evi << std::endl;
    }
    if (hitsVector[evi]->GetTimestamp() / 1E12 <= 1600 or
        hitsVector[evi]->GetTimestamp() / 1E12 >= 3200) {
      continue;
    }
    WF = nullptr;
    hitsVector[evi]->SetSmoothWF(sbSize);
    WF = hitsVector[evi]->GetWFPtr();
    auto pvVec = WF->DetectPeakValleys(peakThreshold);
    // Now check if peak to be integrated over
    // valley of previous and next peak to be > 20 from the peak
    // for first peak, location to be > 100
    // for last peak, location to be < wfsize-100
    int iter = 0;
    auto peaks = pvVec.first;
    auto valleys = pvVec.second;
    int startCheck = 100;
    int stopCheck = WF->GetSize() - 100;
    int distToValleys = 20;
    double integratedVal;
    int integrateRange = 20;
    double lowEnergy = 7500, hiEnergy = 8500;
    while (iter < peaks.size()) {
      int peakPos = peaks[iter];
      if (peakPos >= WF->GetSize() - startCheck) {
        break;
      }
      if (peakPos > startCheck and peakPos < stopCheck and
          ((peakPos - valleys[2 * iter - 1]) > distToValleys) and
          (valleys[2 * iter + 2] - peakPos > distToValleys)) {
        // Integrate this peak in +- 20
        integratedVal = WF->IntegrateSmoothWaveForm(peakPos - integrateRange,
                                                    peakPos + integrateRange);
        hSPE->Fill(integratedVal);

        // if (integratedVal >= lowEnergy and integratedVal <= hiEnergy) {
        if (hitsVector[evi]->GetEnergy() >= lowEnergy and
            hitsVector[evi]->GetEnergy() <= hiEnergy and
            hitsVector[evi]->GetMeanTime() < 2.34) {
          std::cout << evi << ": Peak Position = " << peakPos << std::endl;
          hitsVector[evi]->Print();
          WF->Plot();
          std::cout << "Do you want to see the next waveform? (y/n): ";
          std::getline(std::cin, userInput);
          if (userInput != "y" && userInput != "Y") {
            keepGoing = false;
          }
          std::cout << "\n\n";
        }
      }
      iter += 1;
    }
  }
  TCanvas *canvas1 = new TCanvas("canvas1", "Energy Hists", 1600, 1000);
  canvas1->cd();
  hSPE->Draw("HIST");
  canvas1->Update();

  // Waveform Plotting
  // for (evi = 0; evi < nentries && keepGoing; ++evi) {
  //   // if (evi % 1 == 0) {
  //   //   std::cout << evi << " : " << keepGoing << std::endl;
  //   // }
  //   if (fabs(hitsVector[evi]->GetEnergy() - 300) <
  //       100 // and
  //           //   hitsVector[evi]->GetTimestamp() / 1e12 > 30000 and
  //           //  fabs(hitsVector[evi]->GetMeanTime() - 2.7) < 0.1
  //   ) {
  //     WF = nullptr;
  //     WF = hitsVector[evi]->GetWFPtr();
  //     //   WF->SetSmooth(80, "Gauss");
  //     hitsVector[evi]->Print();
  //     WF->SetTracesFFT();
  //     WF->Plot();
  //     if (WF) {
  //       waveformVector.push_back(*WF);
  //     }

  //     std::cout << "Do you want to see the next waveform? (y/n): ";
  //     std::getline(std::cin, userInput);
  //     if (userInput != "y" && userInput != "Y") {
  //       keepGoing = false;
  //     }
  //   }
  // }
  // UShort_t wfSz = WF->GetSize();
  // std::cout << "Got size of waveforms" << wfSz << std::endl;
  // digiAnalysis::WaveForm *WFAveraged =
  //     new digiAnalysis::WaveForm(wfSz, waveformVector);
  // //   WFAveraged->FitExponential(280, 1100);
  // WFAveraged->SetSmooth(80, "Gauss");
  // WFAveraged->SetTracesFFT();
  // WFAveraged->Plot();

  // // Energy after moving baseline correction
  // TH1 *hE = new TH1F("hE", "Energy", 16384, 0, 16384);
  // TH1 *hEEval = new TH1F("hEEval", "Energy Eval", 16384, 0, 16384);
  // for (evi = 0; evi < nentries && keepGoing; ++evi) {
  //   hE->Fill(hitsVector[evi]->GetEnergy());
  //   WF = nullptr;
  //   WF = hitsVector[evi]->GetWFPtr();
  //   WF->SetTracesMovBLCorr();
  //   hitsVector[evi]->SetEvalEnergy();
  //   hEEval->Fill(hitsVector[evi]->GetEvalEnergy());
  // }
  // TCanvas *canvas1 = new TCanvas("canvas1", "Energy Hists", 1600, 1000);
  // canvas1->cd();
  // hE->Draw("HIST");
  // hEEval->SetLineColor(kBlack);
  // hEEval->Draw("HIST SAME");
  // canvas1->Update();
  // //   WFAveraged->Plot();
  std::cout << "Ending the run" << std::endl;
  fApp->Run();

#endif
  return 0;
}