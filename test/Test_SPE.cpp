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

  std::string fname = "/home/kirtikesh/analysisSSD/DATA/SPE/run_noSource_19Sep/"
                      "FILTERED/DataF_run_noSource_19Sep.root";

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

  for (evi = 0; evi < nentries && keepGoing; ++evi) {
    // if (evi % 1 == 0) {
    //   std::cout << evi << " : " << keepGoing << std::endl;
    // }
    if (fabs(hitsVector[evi]->GetEnergy() - 4000) <
        2000 // and
             // fabs(hitsVector[evi]->GetMeanTime() - 2.15) < 0.25
    ) {
      WF = nullptr;
      WF = hitsVector[evi]->GetWFPtr();
      //   WF->SetSmooth(80, "Gauss");
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
  std::cout << "Ending the run" << std::endl;
  fApp->Run();

#endif
  return 0;
}