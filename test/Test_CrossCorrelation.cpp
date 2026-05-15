// ################################################################# //
//  This code is meant to create a template using shifted averaging  //
//           with cross correlation used to identify shift           //
// ################################################################# //
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
  std::string fname =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
      "NaI1_15May26_1900_NoSrc_WAVES/UNFILTERED/"
      "Data_NaI1_15May26_1900_NoSrc_WAVES.root";
  // Read to singleHits
  digiAnalysis::Analysis an(0, fname, 0, 000, 0);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();

  digiAnalysis::WaveForm *WF = hitsVector[0]->GetWFPtr();
  WF->SetSmooth(500);
  std::vector<double> tracePrimary = WF->GetTracesSmooth();
  //   for (int iWF = 60; iWF < WF->GetSize(); iWF++) {
  //     if (abs(tracePrimary[iWF] - tracePrimary[iWF - 20]) > 1.5 or
  //         tracePrimary[iWF] > 5)
  //       tracePrimary[iWF] = tracePrimary[iWF - 50];
  //   }
  //   WF->Plot(WF->GetTracesSmooth(), tracePrimary); // this line is only for
  //   setting the cutoff for removing the single pulses
  WF->SetTracesFFT(tracePrimary);
  std::vector<double> trFFT_AmpPrim = WF->GetTracesFFT();
  std::vector<double> trFFT_PhsPrim = WF->GetTracesFFTPhase();

  digiAnalysis::WaveForm *WF1 = hitsVector[10]->GetWFPtr();
  WF1->SetSmooth(500);
  std::vector<double> traceSecondary = WF1->GetTracesSmooth();
  //   for (int iWF = 60; iWF < WF1->GetSize(); iWF++) {
  //     if (abs(traceSecondary[iWF] - traceSecondary[iWF - 20]) > 1.5 or
  //         traceSecondary[iWF] > 5)
  //       traceSecondary[iWF] = traceSecondary[iWF - 50];
  //   }
  WF1->SetTracesFFT(traceSecondary);
  std::vector<double> trFFT_AmpSecn = WF1->GetTracesFFT();
  std::vector<double> trFFT_PhsSecn = WF1->GetTracesFFTPhase();

  std::vector<double> trFFT_AmpCorr;
  std::vector<double> trFFT_PhsCorr;
  for (int iterFFT = 0; iterFFT < trFFT_AmpSecn.size(); iterFFT++) {
    trFFT_AmpCorr.push_back(trFFT_AmpPrim[iterFFT] * trFFT_AmpSecn[iterFFT]);
    trFFT_PhsCorr.push_back(trFFT_PhsPrim[iterFFT] - trFFT_PhsSecn[iterFFT]);
  }
  trFFT_AmpCorr[0] = 0;
  std::vector<double> trCorr = WF->EvalIFFT(trFFT_AmpCorr, trFFT_PhsCorr);
  int index = 0;
  auto maxIt = std::max_element(trCorr.begin(), trCorr.end());
  if (maxIt != trCorr.end())
    index = std::distance(trCorr.begin(), maxIt);

  int val = 0;
  int shift = WF->GetSize() - index;
  std::cout << "shift is: " << shift << std::endl;
  WF1->SetSmooth(16, "MovA");
  traceSecondary = WF1->GetTracesSmooth();
  WF->SetSmooth(16, "MovA");
  tracePrimary = WF->GetTracesSmooth();
  for (int iterWF = 0; iterWF < WF->GetSize(); iterWF++) {
    traceSecondary[iterWF] -= trFFT_AmpSecn[0] / WF->GetSize();
    tracePrimary[iterWF] -= trFFT_AmpPrim[0] / WF->GetSize();
  }
  for (int iterWF = 0; iterWF < WF->GetSize(); iterWF++) {
    val = iterWF - shift > 0 ? iterWF - shift : WF->GetSize() + iterWF - shift;
    traceSecondary[iterWF] = traceSecondary[iterWF] - tracePrimary[val];
  }
  WF1->SetTracesFFT(traceSecondary);
  //   WF->Plot(traceSecondary, tracePrimary);
  //   WF->Plot(traceSecondary, trCorr);
  //   WF->Plot(traceSecondary, WF1->GetTracesFFT());
  //   WF->Plot(tracePrimary, trFFT_AmpPrim);
  WF1->Plot(WF1->GetTracesSmooth(), traceSecondary);
  //   WF1->Plot(traceSecondary, "SAME_KGreen");
  fApp->Run();
}