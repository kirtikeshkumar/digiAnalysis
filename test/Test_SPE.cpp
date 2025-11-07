#include "Analysis.h"
#include "WaveForm.h"
#include "globals.h"
#include "singleHits.h"
#include <TApplication.h>
#include <iostream>
int main(int argc, char *argv[]) {
#ifdef WAVES
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;

  std::string fname =
      "/home/kirtikesh/analysisSSD/DATA/SPE/"
      "run_Nov06_Direct_SelfTrigger_10lsb_BlueLED_Pico15_"
      "Bias2100V/FILTERED/"
      "DataF_run_Nov06_Direct_SelfTrigger_10lsb_BlueLED_Pico15_Bias2100V.root";

  digiAnalysis::Analysis an(fname, 0, 00000, 0);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int nentries = hitsVector.size();
  std::cout << "Got the vector from an: " << nentries << std::endl;
  int evi = 0;
  std::string userInput;
  bool keepGoing = true;
  digiAnalysis::WaveForm *WF = nullptr;
  std::vector<digiAnalysis::WaveForm> waveformVector;

  // Waveform Plotting

  for (evi = 0; evi < nentries && keepGoing; ++evi) {
    // if (evi % 1 == 0) {
    //   std::cout << evi << " : " << keepGoing << std::endl;
    // }
    if (fabs(hitsVector[evi]->GetEnergy() - 8300) <
        100 // and
            //   hitsVector[evi]->GetTimestamp() / 1e12 > 30000 and
            //  fabs(hitsVector[evi]->GetMeanTime() - 2.7) < 0.1
    ) {
      WF = nullptr;
      WF = hitsVector[evi]->GetWFPtr();
      //   WF->SetSmooth(80, "Gauss");
      hitsVector[evi]->Print();
      WF->SetTracesFFT();
      WF->Plot();
      if (WF) {
        waveformVector.push_back(*WF);
      }

      std::cout << "Do you want to see the next waveform? (y/n): ";
      std::getline(std::cin, userInput);
      if (userInput != "y" && userInput != "Y") {
        keepGoing = false;
      }
    }
  }
  UShort_t wfSz = WF->GetSize();
  std::cout << "Got size of waveforms" << wfSz << std::endl;
  digiAnalysis::WaveForm *WFAveraged =
      new digiAnalysis::WaveForm(wfSz, waveformVector);
  //   WFAveraged->FitExponential(280, 1100);
  WFAveraged->SetSmooth(80, "Gauss");
  WFAveraged->SetTracesFFT();
  WFAveraged->Plot();

  // Energy after moving baseline correction
  TH1 *hE = new TH1F("hE", "Energy", 16384, 0, 16384);
  TH1 *hEEval = new TH1F("hEEval", "Energy Eval", 16384, 0, 16384);
  for (evi = 0; evi < nentries && keepGoing; ++evi) {
    hE->Fill(hitsVector[evi]->GetEnergy());
    WF = nullptr;
    WF = hitsVector[evi]->GetWFPtr();
    WF->SetTracesMovBLCorr();
    hitsVector[evi]->SetEvalEnergy();
    hEEval->Fill(hitsVector[evi]->GetEvalEnergy());
  }
  TCanvas *canvas1 = new TCanvas("canvas1", "Energy Hists", 1600, 1000);
  canvas1->cd();
  hE->Draw("HIST");
  hEEval->SetLineColor(kBlack);
  hEEval->Draw("HIST SAME");
  canvas1->Update();
  //   WFAveraged->Plot();
  std::cout << "Ending the run" << std::endl;
  fApp->Run();

#endif
  return 0;
}