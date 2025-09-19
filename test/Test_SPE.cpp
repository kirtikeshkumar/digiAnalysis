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
  int evi = 0;
  std::string userInput;
  bool keepGoing = true;
  digiAnalysis::WaveForm *WF = nullptr;
  for (evi = 0; evi < nentries && keepGoing; ++evi) {
    if (fabs(hitsVector[evi]->GetEnergy() - 4000) < 200) {
      WF = nullptr;
      WF = hitsVector[evi]->GetWFPtr();
      WF->SetSmooth(80, "Gauss");
      WF->Plot();

      std::cout << "Do you want to see the next waveform? (y/n): ";
      std::getline(std::cin, userInput);
      if (userInput != "y" && userInput != "Y") {
        keepGoing = false;
      }
    }
  }
#endif
  return 0;
}