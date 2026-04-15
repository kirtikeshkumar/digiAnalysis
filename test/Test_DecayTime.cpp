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
#include <ratio>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {
#ifdef WAVES
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
      "NaI_12_CoincidenceStudies_Na_HV_1900V_1365V_96nsCoinc_60min_2Vpp_WAVES/"
      "FILTERED/SDataF_NaI_12_CoincidenceStudies_Na_HV_1900V_1365V_96nsCoinc_"
      "60min_2Vpp_WAVES.root";

  digiAnalysis::Analysis an(0, fname, 0000, 00010, 0);
  std::cout << "getting the vector from an" << std::endl;

  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int nentries = hitsVector.size();
  std::cout << "got the vector from an: " << nentries << std::endl;

  bool keepGoing = true;
  std::string userInput;
  for (int iter = 0; iter < nentries && keepGoing; iter++) {
    digiAnalysis::WaveForm *wfPtr = hitsVector[iter]->GetWFPtr();
    wfPtr->SetSmooth(100);
    wfPtr->FitExponential(1, 1050, 3400);
    wfPtr->Plot();
    std::cout << "Do you want to see the next waveform? (y/n): ";
    std::getline(std::cin, userInput);
    if (userInput != "y" && userInput != "Y") {
      keepGoing = false;
    }
  }

  fApp->Run();
#else
  std::cout << "Build with WAVES set to ON" << std::endl;
#endif
}