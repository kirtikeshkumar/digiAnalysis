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

int main(int argc, char *argv[])
{
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;

  std::string fname =
      "/media/kirtikesh/UbuntuFiles/NaI/LowEnergy/SPE/"
      "SPE_Ch0_NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp.root";
  TFile *fp = new TFile(fname.c_str(), "READ");
  TTree *tr = (TTree *)fp->Get("SPE_WF");

  if (tr)
  {
    bool keepGoing = true;
    std::string userInput;

    digiAnalysis::WaveForm *SPEWF = new digiAnalysis::WaveForm();
    tr->SetBranchAddress("SPE", &SPEWF);
    Long64_t nentries = tr->GetEntries();
    std::cout << "Got " << nentries << " entries" << std::endl;

    std::vector<std::unique_ptr<digiAnalysis::WaveForm>> WFVec;
    for (int iter = 15000; iter < 15500; iter++)
    {
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
    for (int iter = 0; iter < filterSz; iter++)
    {
      if (iter < filterFlatRange)
      {
        filter[iter] = 1.0;
      }
      else if (iter - filterFlatRange < 3 * filterGaussSigma)
      {
        filter[iter] = TMath::Gaus(iter, filterFlatRange, filterGaussSigma);
      }
      else
      {
        filter[iter] = 0;
      }
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
    std::string function;
    function = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])) + "
               "([3]+[4]*x)";
    std::vector<double> parLimits = {0, 100, 220.0, 270.0, 1, 10, -5, 5, -1, 1};
    double intVal = 0;
    double amplitude, sigma;
    int count = 0;
    for (int iter = 0; iter < WFVec.size() && keepGoing && count < 100; iter++)
    {
      if (iter % 1000 == 0)
      {
        std::cout << iter << std::endl;
      }
      WFVec[iter]->FitFunction(function, parLimits, 5, wfSz - 5);
      intVal = WFVec[iter]->IntegrateWaveForm(200, 350); // WFVec[iter]->GetFitPar(0) * WFVec[iter]->GetFitPar(2) *
                                                         // TMath::Sqrt(2 * TMath::Pi());
      // amplitude = WFVec[iter]->GetFitPar(0);
      sigma = WFVec[iter]->GetFitPar(2);

      // gaussAmp->Fill(amplitude);
      // gaussMean->Fill(WFVec[iter]->GetFitPar(1));
      // gaussSig->Fill(sigma);
      gaussInt->Fill(intVal);
      // polSlope->Fill(WFVec[iter]->GetFitPar(4));

      // if ((amplitude > 3 and amplitude < 4) or (sigma > 8 and sigma < 9)) {
      // std::cout << "Amplitude: " << amplitude << " | sigma: " << sigma
      //           << std::endl;
      // WFVec[iter]->Plot();
      // std::cout << "Do you want to see the next waveform? (y/n): ";
      // std::getline(std::cin, userInput);
      // if (userInput != "y" && userInput != "Y") {
      //   keepGoing = false;
      // }
      // }

      if ((sigma > 3 and sigma < 5) and (intVal > 80 and intVal < 400)
          // ){
          and count < 100)
      {

        // if (count == 0)
        //   WFVec[iter]->Plot(WFVec[iter]->GetTraces());
        // if (count != 0)
        //   WFVec[iter]->Plot(WFVec[iter]->GetTraces(), "SAME");

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

    // TCanvas *c1 = new TCanvas("c1", "Amplitude", 800, 600);
    // gaussAmp->Draw("Hist");
    // TCanvas *c2 = new TCanvas("c2", "Mean", 800, 600);
    // gaussMean->Draw("Hist");
    // TCanvas *c3 = new TCanvas("c3", "Sigma", 800, 600);
    // gaussSig->Draw("Hist");
    // TCanvas *c4 = new TCanvas("c4", "IntVal", 800, 600);
    // gaussInt->Draw("Hist");
    // TCanvas *c5 = new TCanvas("c5", "Slope", 800, 600);
    // polSlope->Draw("Hist");

    std::cout << selWF.size() << " Waveforms to be averaged" << std::endl;
    digiAnalysis::WaveForm *WFAveraged =
        new digiAnalysis::WaveForm(wfSz, selWF);
    WFAveraged->SetTracesFFT();
    std::vector<double> SPEFFT_Amp = WFAveraged->GetTracesFFT();
    std::vector<double> SPEFFT_Phase = WFAveraged->GetTracesFFTPhase();
    // for (int iter = 0; iter < SPEFFT_Amp.size(); iter++)
    // {
    //   SPEFFT_Amp[iter] *= filter[iter];
    // }
    // WFAveraged->ReSetTracesFFT(SPEFFT_Amp, SPEFFT_Phase);
    // std::vector<double> trace = WFAveraged->EvalIFFT(SPEFFT_Amp, SPEFFT_Phase);
    // WFAveraged->Plot(WFAveraged->GetTraces(), trace);
    // WFAveraged->Plot(WFAveraged->GetTraces(), "SAME_kGreen_4");

    // WFAveraged->Plot();

    // Use Average WF to evaluate the NPE of the SPE

    std::vector<double> trace = WFVec[101]->GetTraces();
    WFVec[101]->SetTracesFFT();
    std::vector<double> trFFT_Amp = WFAveraged->GetTracesFFT();        // WFVec[101]->GetTracesFFT();
    std::vector<double> trFFT_Phase = WFAveraged->GetTracesFFTPhase(); // WFVec[101]->GetTracesFFTPhase();
    for (int iter = 0; iter < wfSz / 2 + 1; iter++)
    {
      trFFT_Amp[iter] /= (0.5 * SPEFFT_Amp[iter]);
      trFFT_Amp[iter] *= filter[iter];
      // trFFT_Phase[iter] = trFFT_Phase[iter] - SPEFFT_Phase[iter];
    }
    std::vector<double> deconvolvedTrace = WFAveraged->EvalIFFT(trFFT_Amp, trFFT_Phase);
    WFVec[101]->Plot(deconvolvedTrace, trFFT_Amp);

    // for (int iter = 0; iter < WFVec.size(); iter++)
    // {
    //   /* code */
    // }

    fApp->Run();
    return 0;
  }
  else
  {
    std::cout << "Could Not Read The Tree" << std::endl;
    return 1;
  }
}