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
      "SPEFiles/"
      "SPE_Ch0_NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp.root";
  TFile *fp = new TFile(fname.c_str(), "READ");
  TTree *tr = (TTree *)fp->Get("SPE_WF");

  if (tr) {
    digiAnalysis::WaveForm *SPEWF = new digiAnalysis::WaveForm();
    tr->SetBranchAddress("SPE", &SPEWF);
    Long64_t nentries = tr->GetEntries();
    std::cout << "Got " << nentries << " entries" << std::endl;
    std::vector<std::unique_ptr<digiAnalysis::WaveForm>> WFVec;
    for (int iter = 0; iter < 10; iter++) {
      tr->GetEntry(iter);
      WFVec.push_back(std::make_unique<digiAnalysis::WaveForm>(*SPEWF));
    }

    // Check data loading by plotting
    WFVec[0]->SetTracesFFT();
    // WFVec[0]->Plot();

    int wfSz = WFVec[0]->GetSize();

    // Defining the filter for smoothing the SPE Waveform
    int filterSz = wfSz / 2 + 1;
    int filterCutOff = 100; // this corresponds in frequency to filterCutOff *
                            // (500/NSampleSPE) MHz
    int filterFlatRange = 50;
    int filterGaussSigma = (filterCutOff - filterFlatRange) / 3;
    std::vector<Double_t> filter(filterSz);
    for (int iter = 0; iter < filterSz; iter++) {
      if (iter < filterFlatRange) {
        filter[iter] = 1.0;
      } else if (iter - filterFlatRange < 3 * filterGaussSigma) {
        filter[iter] = TMath::Gaus(iter, filterFlatRange, filterGaussSigma);
      } else {
        filter[iter] = 0;
      }
      //   std::cout << iter << " : " << filter[iter] << std::endl;
    }

    std::vector<double> SPEFFT_Amp = WFVec[0]->GetTracesFFT();
    // std::vector<double> SPEFFT_Amp1 = WFVec[0]->GetTracesFFT();
    std::vector<double> SPEFFT_Phase = WFVec[0]->GetTracesFFTPhase();
    std::cout << "FFT size: " << SPEFFT_Amp.size() << std::endl;
    std::cout << "Filter size: " << filter.size() << std::endl;
    for (int iter = 0; iter < SPEFFT_Amp.size(); iter++) {
      SPEFFT_Amp[iter] *= filter[iter];
    }
    WFVec[0]->ReSetTracesFFT(SPEFFT_Amp, SPEFFT_Phase);
    std::cout << "Reset FFT done" << std::endl;
    std::vector<double> trace = WFVec[0]->EvalIFFT(SPEFFT_Amp, SPEFFT_Phase);
    std::cout << "IFFT Evaluated" << std::endl;
    // WFVec[0]->SetWaveForm(trace);
    // WFVec[0]->Plot(SPEFFT_Amp, SPEFFT_Amp1);
    WFVec[0]->Plot(WFVec[0]->GetTraces(), trace);

    fApp->Run();
    return 0;
  } else {
    std::cout << "Could Not Read The Tree" << std::endl;
    return 1;
  }
}