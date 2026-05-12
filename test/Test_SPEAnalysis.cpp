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
      "SPE_Ch0_NaI13_12May26_1900_1345_Cs_Coinc144ns_35cm_NoCollimation_1.root";
  TFile *fp = new TFile(fname.c_str(), "READ");
  TTree *tr = (TTree *)fp->Get("SPE_WF");

  std::string fnameWF =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
      "CoincidenceStudies/PairFiles/"
      "Pair_NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp.root";

  digiAnalysis::Analysis an(fnameWF, 0000, 1000, 0);
  std::vector<std::unique_ptr<digiAnalysis::Pair>> &vecOfPairs =
      an.GetPairsVec();
  int pairVecSz = vecOfPairs.size();

  if (tr) {
    bool keepGoing = true;
    std::string userInput;

    digiAnalysis::WaveForm *SPEWF = new digiAnalysis::WaveForm();
    tr->SetBranchAddress("SPE", &SPEWF);
    Long64_t nentries = tr->GetEntries();
    std::cout << "Got " << nentries << " entries" << std::endl;

    std::vector<std::unique_ptr<digiAnalysis::WaveForm>> WFVec;
    for (int iter = 0; iter < nentries; iter++) {
      tr->GetEntry(iter);
      WFVec.push_back(std::make_unique<digiAnalysis::WaveForm>(*SPEWF));
    }

    std::cout << "Added " << WFVec.size() << " traces for analysis."
              << std::endl;

    // Check data loading by plotting
    // WFVec[0]->SetTracesFFT();
    // WFVec[0]->Plot();

    int wfSz = WFVec[0]->GetSize();

    // Defining the filter for smoothing the SPE Waveform
    int filterSz = wfSz / 2 + 1;
    std::cout << "Filter size: " << filterSz << std::endl;
    int filterCutOff = 90; // this corresponds in frequency to filterCutOff *
                           // (500/NSampleSPE) MHz
    int filterFlatRange = 60;
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

      // if ((iter > 500 and iter < 800) or (iter > 1100 and iter < 1400) or
      //     (iter > 1760 and iter < 1960)) {
      //   filter[iter] = 0.;
      // }
      //   std::cout << iter << " : " << filter[iter] << std::endl;
    }

    // std::vector<double> SPEFFT_Amp = WFVec[0]->GetTracesFFT();
    // // std::vector<double> SPEFFT_Amp1 = WFVec[0]->GetTracesFFT();
    // std::vector<double> SPEFFT_Phase = WFVec[0]->GetTracesFFTPhase();
    // std::cout << "FFT size: " << SPEFFT_Amp.size() << std::endl;
    // std::cout << "Filter size: " << filter.size() << std::endl;
    // for (int iter = 0; iter < SPEFFT_Amp.size(); iter++) {
    //   SPEFFT_Amp[iter] *= filter[iter];
    // }
    // WFVec[0]->ReSetTracesFFT(SPEFFT_Amp, SPEFFT_Phase);
    // std::cout << "Reset FFT done" << std::endl;
    // std::vector<double> trace = WFVec[0]->EvalIFFT(SPEFFT_Amp, SPEFFT_Phase);
    // std::cout << "IFFT Evaluated" << std::endl;
    // // WFVec[0]->SetWaveForm(trace);
    // // WFVec[0]->Plot(SPEFFT_Amp, SPEFFT_Amp1);
    // WFVec[0]->Plot(WFVec[0]->GetTraces(), trace);

    // Try fitting to extract central SPE Gaussian
    // Assume single Gaussian + polynomial(1)
    // SPE has integrated charge of ~172
    //

    std::vector<digiAnalysis::WaveForm> selWF;
    TH1 *gaussAmp = new TH1F("gaussAmp", "Amplitude", 10000, 0, 1000);
    TH1 *gaussMean = new TH1F("gaussMean", "Mean", 1000, 200, 350);
    TH1 *gaussSig = new TH1F("gaussSig", "Sigma", 1000, 0, 10);
    TH1 *gaussInt = new TH1F("gaussInt", "Integrated Gaussian", 10000, 0, 1000);
    TH1 *polSlope = new TH1F("polSlope", "Slope", 10000, -1, 1);
    TH2 *hsigInt = new TH2F("hsigInt", "Sigma vs Energy", 2000, -1000, 1000,
                            1000, 0.0, 10);
    TH2 *hsigAmp =
        new TH2F("hsigAmp", "Sigma vs Amplitude", 100, 0, 100, 1000, 0.0, 10);
    TH2 *hAmpInt = new TH2F("hAmpInt", "Amplitude vs Energy", 2000, -1000, 1000,
                            100, 0.0, 100);

    std::string function;
    function = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])) + "
               "([3]+[4]*x)";
    std::vector<double> parLimits = {0,  100, 120.0, 170.0, 0.01,
                                     10, -5,  5,     -1,    1};
    double intVal = 0;
    double amplitude, sigma;
    int count = 0;
    for (int iter = 0; iter < WFVec.size() && keepGoing && count < 300;
         iter++) {
      if (iter % 1000 == 0) {
        std::cout << iter << std::endl;
      }
      WFVec[iter]->FitFunction(function, parLimits, 5, wfSz - 5);
      intVal = WFVec[iter]->IntegrateWaveForm(
          100, 250); // WFVec[iter]->GetFitPar(0) * WFVec[iter]->GetFitPar(2) *
                     // TMath::Sqrt(2 * TMath::Pi());
      amplitude = WFVec[iter]->GetFitPar(0);
      sigma = WFVec[iter]->GetFitPar(2);

      gaussAmp->Fill(amplitude);
      gaussMean->Fill(WFVec[iter]->GetFitPar(1));
      gaussSig->Fill(sigma);
      if (sigma > 1) {
        gaussInt->Fill(intVal);
      }
      polSlope->Fill(WFVec[iter]->GetFitPar(4));
      hAmpInt->Fill(intVal, amplitude);
      hsigInt->Fill(intVal, sigma);
      hsigAmp->Fill(amplitude, sigma);

      // if ((amplitude > 3 and amplitude < 4) or (sigma > 8 and sigma < 9))
      // { std::cout << "Amplitude: " << amplitude << " | sigma: " << sigma
      //           << std::endl;
      // WFVec[iter]->Plot();
      // std::cout << "Do you want to see the next waveform? (y/n): ";
      // std::getline(std::cin, userInput);
      // if (userInput != "y" && userInput != "Y") {
      //   keepGoing = false;
      // }
      // }

      if ((sigma > 1 and sigma < 2) and
          (intVal > 0 and intVal < 1000)
          // ){
          and count < 5000) {

        if (count == 0)
          WFVec[iter]->Plot(WFVec[iter]->GetTraces());
        if (count != 0 and count < 100)
          WFVec[iter]->Plot(WFVec[iter]->GetTraces(), "SAME");

        // std::cout << "Count: " << count << " | Add this WaveForm? (y/n): ";
        // std::getline(std::cin, userInput);
        // if (userInput == "y" or userInput == "Y")
        // {
        // std::cout << "Count: " << count << std::endl;
        selWF.push_back(digiAnalysis::WaveForm(*WFVec[iter]));
        count++;
        // }
      }
    }

    TCanvas *c1 = new TCanvas("c1", "Amplitude", 800, 600);
    gaussAmp->Draw("Hist");
    TCanvas *c2 = new TCanvas("c2", "Mean", 800, 600);
    gaussMean->Draw("Hist");
    TCanvas *c3 = new TCanvas("c3", "Sigma", 800, 600);
    gaussSig->Draw("Hist");
    TCanvas *c4 = new TCanvas("c4", "IntVal", 800, 600);
    gaussInt->Draw("Hist");
    TCanvas *c5 = new TCanvas("c5", "Slope", 800, 600);
    polSlope->Draw("Hist");
    TCanvas *c6 = new TCanvas("c6", "Amplitude vs Energy", 800, 600);
    hAmpInt->Draw("COLZ");
    TCanvas *c7 = new TCanvas("c7", "Sigma vs Energy", 800, 600);
    hsigInt->Draw("COLZ");
    TCanvas *c8 = new TCanvas("c8", "Sigma vs Amplitude", 800, 600);
    hsigAmp->Draw("COLZ");

    std::cout << selWF.size() << " Waveforms to be averaged" << std::endl;
    digiAnalysis::WaveForm *WFAveraged =
        new digiAnalysis::WaveForm(wfSz, selWF);
    WFAveraged->SetTracesFFT();
    std::vector<double> SPEFFT_Amp = WFAveraged->GetTracesFFT();
    std::vector<double> SPEFFT_Phase = WFAveraged->GetTracesFFTPhase();
    // WFAveraged->Plot();
    for (int iter = 0; iter < SPEFFT_Amp.size(); iter++) {
      SPEFFT_Amp[iter] *= filter[iter];
    }
    WFAveraged->ReSetTracesFFT(SPEFFT_Amp, SPEFFT_Phase);
    std::vector<double> trace = WFAveraged->EvalIFFT(SPEFFT_Amp, SPEFFT_Phase);
    // WFAveraged->Plot(WFAveraged->GetTraces(), trace);
    WFAveraged->Plot(WFAveraged->GetTraces(), "SAME_kGreen_4");

    // WFAveraged->Plot();

    // Use Average WF to evaluate the NPE of the SPE

    // std::vector<double> trace = WFVec[101]->GetTraces();
    // WFVec[101]->SetTracesFFT();
    // std::vector<double> trFFT_Amp =
    //     WFAveraged->GetTracesFFT(); // WFVec[101]->GetTracesFFT();
    // std::vector<double> trFFT_Phase =
    //     WFAveraged->GetTracesFFTPhase(); // WFVec[101]->GetTracesFFTPhase();
    // for (int iter = 0; iter < wfSz / 2 + 1; iter++) {
    //   trFFT_Amp[iter] /= (0.5 * SPEFFT_Amp[iter]);
    //   trFFT_Amp[iter] *= filter[iter];
    //   // trFFT_Phase[iter] = trFFT_Phase[iter] - SPEFFT_Phase[iter];
    // }
    // std::vector<double> deconvolvedTrace =
    //     WFAveraged->EvalIFFT(trFFT_Amp, trFFT_Phase);
    // WFVec[101]->Plot(deconvolvedTrace, trFFT_Amp);

    // for (int iter = 0; iter < WFVec.size(); iter++)
    // {
    //   /* code */
    // }

    // digiAnalysis::WaveForm *highEWF;
    // for (int wfiter = 30; wfiter < pairVecSz; wfiter++) {
    //   if (vecOfPairs[wfiter]->GetPairHitEnergy(0) < 50) {
    //     vecOfPairs[wfiter]->Print();
    //     highEWF = vecOfPairs[wfiter]->GetHitPtr(0)->GetWFPtr();
    //     break;
    //   }
    // }

    // std::vector<double> wfAvTrace = WFAveraged->GetTraces();
    // WFAveraged->SetBaseLine(wfAvTrace, 400, 200);
    // std::vector<double> padtrace(highEWF->GetSize(),
    // WFAveraged->GetBaseLine()); std::copy(wfAvTrace.begin(), wfAvTrace.end(),
    // padtrace.begin()); WFAveraged->SetWaveForm(padtrace, 0, padtrace.size() -
    // 1, 400, 200); WFAveraged->SetTracesFFT(); std::vector<double> SPEFFT_Amp
    // = WFAveraged->GetTracesFFT(); std::vector<double> SPEFFT_Phase =
    // WFAveraged->GetTracesFFTPhase(); std::cout << "SPE FFT size: " <<
    // SPEFFT_Amp.size() << std::endl;

    // highEWF->SetTracesFFT();
    // int evt = 105;
    // std::vector<double> speAvTrace = WFVec[evt]->GetTraces();
    // WFVec[evt]->SetBaseLine(speAvTrace, 400, 200);
    // std::vector<double> spepadtrace(highEWF->GetSize(),
    //                                 WFVec[evt]->GetBaseLine());
    // std::copy(speAvTrace.begin(), speAvTrace.end(), spepadtrace.begin());
    // WFVec[evt]->SetWaveForm(spepadtrace, 0, spepadtrace.size() - 1, 400,
    // 200); WFVec[evt]->SetTracesFFT(); std::vector<double> trFFT_Amp =
    // highEWF->GetTracesFFT();
    // // WFAveraged->GetTracesFFT(); //WFVec[evt]->GetTracesFFT(); //
    // std::vector<double> trFFT_Phase = highEWF->GetTracesFFTPhase();
    // // WFAveraged->GetTracesFFTPhase(); //WFVec[evt]->GetTracesFFTPhase(); //
    // std::cout << "highE FFT size: " << SPEFFT_Amp.size() << std::endl;
    // // highEWF->Plot();
    // for (int iter = 0; iter < trFFT_Amp.size(); iter++) {
    //   trFFT_Amp[iter] /= (SPEFFT_Amp[iter]);
    //   trFFT_Amp[iter] *= filter[iter];
    //   // trFFT_Phase[iter] = trFFT_Phase[iter] - SPEFFT_Phase[iter];
    // }
    // std::vector<double> deconvolvedTrace =
    //     highEWF->EvalIFFT(trFFT_Amp, trFFT_Phase);

    // std::transform(deconvolvedTrace.begin(), deconvolvedTrace.end(),
    //                deconvolvedTrace.begin(), [](double x) { return 5.0 * x;
    //                });

    // // highEWF->Plot(deconvolvedTrace, trFFT_Amp);
    // highEWF->Plot(highEWF->GetTraces(), deconvolvedTrace);

    fApp->Run();
    return 0;
  } else {
    std::cout << "Could Not Read The Tree" << std::endl;
    return 1;
  }
}