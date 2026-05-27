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
#include <numeric>
#include <ratio>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
      "SPEFiles/"
      "SPE_Ch0_NaI1_21May26_1900_Cs_WAVES_2_BLCorrected.root";
  TFile *fp = new TFile(fname.c_str(), "READ");
  TTree *tr = (TTree *)fp->Get("SPE_WF");

  // this file is for estimating npe
  // std::string fnameWF =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
  //     "CoincidenceStudies/PairFiles/"
  //     "Pair_NaI13_12May26_1900_1345_Cs_Coinc144ns_35cm_NoCollimation_1.root";

  // digiAnalysis::Analysis an(fnameWF, 0000, 1000, 0);
  // std::vector<std::unique_ptr<digiAnalysis::Pair>> &vecOfPairs =
  //     an.GetPairsVec();
  // int pairVecSz = vecOfPairs.size();

  std::string fnameWF =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
      "NaI1_21May26_1900_Cs_WAVES_2/FILTERED/"
      "DataF_NaI1_21May26_1900_Cs_WAVES_2_BLCorrected.root";
  digiAnalysis::Analysis an(fnameWF, 0000, 1000, 0);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int hitsVecSz = hitsVector.size();

  if (tr) {
    bool keepGoing = true;
    std::string userInput;

    digiAnalysis::WaveForm *SPEWF = new digiAnalysis::WaveForm();
    tr->SetBranchAddress("SPE", &SPEWF);
    Long64_t nentries = tr->GetEntries();
    std::cout << "Got " << nentries << " entries" << std::endl;

    std::vector<std::unique_ptr<digiAnalysis::WaveForm>> WFVec;
    for (int iter = 0; iter < 5000; iter++) { // nentries
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
    int filterSz = 2501;
    std::cout << "Filter size: " << filterSz << std::endl;
    int filterCutOff = 1100; // this corresponds in frequency to filterCutOff *
                             // (500/NSampleSPE) MHz
    int filterFlatRange = 1000;
    int filterGaussSigma = (filterCutOff - filterFlatRange) / 3;
    std::vector<Double_t> filter(filterSz);
    for (int iter = 0; iter < filterSz; iter++) {
      if (iter < filterFlatRange) {
        filter[iter] = 1.0;
      } else if (iter - filterFlatRange < 5 * filterGaussSigma) {
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
    //

    std::vector<digiAnalysis::WaveForm> selWF1, selWF2, selWF3, selWF4, selWF;
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
    for (int iter = 0; iter < WFVec.size();
         iter++) { // && keepGoing && count < 300
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
      // if (sigma > 1) {
      gaussInt->Fill(intVal);
      // }
      polSlope->Fill(WFVec[iter]->GetFitPar(4));
      hAmpInt->Fill(intVal, amplitude);
      hsigInt->Fill(intVal, sigma);
      hsigAmp->Fill(amplitude, sigma);

      // if ((amplitude > 3 and amplitude < 4) or (sigma > 8 and sigma < 9))
      // { std::cout << "Amplitude: " << amplitude << " | sigma: " << sigma
      //           << std::endl;
      // WFVec[iter]->SetSmooth(16, "MovA");
      // WFVec[iter]->Plot();
      // std::cout << "Do you want to see the next waveform? (y/n): ";
      // std::getline(std::cin, userInput);
      // if (userInput != "y" && userInput != "Y") {
      //   keepGoing = false;
      // }
      // }

      if ( //(sigma > 1.2 and sigma < 1.6) and
          (intVal > 400 and intVal < 500)
          // ){
          and count < 50000) {

        // if (count == 0)
        //   WFVec[iter]->Plot(WFVec[iter]->GetTraces());
        // if (count != 0 and count < 100)
        //   WFVec[iter]->Plot(WFVec[iter]->GetTraces(), "SAME");

        // std::cout << "Count: " << count << " | Add this WaveForm? (y/n): ";
        // std::getline(std::cin, userInput);
        // if (userInput == "y" or userInput == "Y")
        // {
        // std::cout << "Count: " << count << std::endl;
        selWF1.push_back(digiAnalysis::WaveForm(*WFVec[iter]));
        count++;
        // }
      }

      if ((sigma > 0.4 and sigma < 0.8) and
          (intVal > 0 and intVal < 1000)
          // ){
          and count < 50000) {
        selWF2.push_back(digiAnalysis::WaveForm(*WFVec[iter]));
        // WFVec[iter]->Plot(WFVec[iter]->GetTraces(), "SAME_kBlue");
        count++;
      }

      if ((sigma > 0.22 and sigma < 0.33) and
          (intVal > 0 and intVal < 1000)
          // ){
          and count < 50000) {
        selWF3.push_back(digiAnalysis::WaveForm(*WFVec[iter]));
        count++;
      }
      if ((sigma > 0. and sigma < 0.1) and
          (intVal > 0 and intVal < 1000)
          // ){
          and count < 50000) {
        selWF.push_back(digiAnalysis::WaveForm(*WFVec[iter]));
        count++;
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

    std::cout << selWF1.size() << " Waveforms to be averaged" << std::endl;
    digiAnalysis::WaveForm *WFAveraged =
        new digiAnalysis::WaveForm(wfSz, selWF1);
    // WFAveraged->SetTracesFFT();
    // std::vector<double> SPEFFT_Amp = WFAveraged->GetTracesFFT();
    // std::vector<double> SPEFFT_Phase = WFAveraged->GetTracesFFTPhase();
    // WFAveraged1->Plot(WFAveraged1->GetTraces(), "kBlue_3");
    // for (int iter = 0; iter < SPEFFT_Amp.size(); iter++) {
    //   SPEFFT_Amp[iter] *= filter[iter];
    // }
    // WFAveraged->ReSetTracesFFT(SPEFFT_Amp, SPEFFT_Phase);
    // std::vector<double> trace = WFAveraged->EvalIFFT(SPEFFT_Amp,
    // SPEFFT_Phase); WFAveraged->Plot(WFAveraged->GetTraces(), trace);
    // WFAveraged->Plot(WFAveraged->GetTraces(), "SAME_kGreen_4");

    digiAnalysis::WaveForm *WFAveraged2 =
        new digiAnalysis::WaveForm(wfSz, selWF2);
    // WFAveraged2->Plot(WFAveraged2->GetTraces(), "SAME_kGreen_3");
    digiAnalysis::WaveForm *WFAveraged3 =
        new digiAnalysis::WaveForm(wfSz, selWF3);
    // WFAveraged1->Plot(WFAveraged3->GetTraces(), "SAME_kRed_3");
    digiAnalysis::WaveForm *WFAveraged4 =
        new digiAnalysis::WaveForm(wfSz, selWF4);
    // WFAveraged4->Plot(WFAveraged4->GetTraces(), "SAME_kMagenta_3");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Fitting the SPE waveform
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double start = 10;
    // std::string funcSPE = Form("[0] + ([1]*exp([2]*(x-%d)))", start);
    // std::vector<double> parLimSPE = {-1, 1, -20, -5, -0.061, -0.05};
    std::string funcSPE =
        Form("[0] - ([1]*exp([2]*(x-[3])))/(exp([4]*(x-[5]))+1)");
    std::vector<double> parLimSPE = {-1,  1,     0.5, 15,   -0.1, -0.0,
                                     140, 170.0, -3,  -0.0, 140,  160.0};
    funcSPE += " + [6]*exp(-0.5*((x-[7])/[8])*((x-[7])/[8]))";
    parLimSPE.insert(parLimSPE.end(), {0.9, 200, 148, 152, 0.9, 2});
    funcSPE += "- [9]*exp(-0.5*((x-[10])/[11])*((x-[10])/[11]))";
    parLimSPE.insert(parLimSPE.end(), {0.1, 100, 150, 156, 0.5, 5});
    int numParSet = 12;
    int numGauss = 8;
    double gaussMid = 150;
    for (int inum = 0; inum < numGauss; inum++) {
      funcSPE += Form(" + [%d]*exp(-0.5*((x-[%d])/[%d])*((x-[%d])/[%d]))",
                      (numParSet), (numParSet + 1), (numParSet + 2),
                      (numParSet + 1), (numParSet + 2));
      numParSet += 3;
      std::cout << inum << " : " << gaussMid << std::endl;
      parLimSPE.insert(parLimSPE.end(),
                       {0.1, 200, gaussMid - 8, gaussMid + 8, 1, 4});

      // below is the function in the paper 2025 JINST 20 P03019
      // funcSPE += Form(
      //     " + ((x>[%d]) ?
      //     [%d]*exp(-0.5*pow(TMath::Log((x-[%d])/[%d])/[%d],2)) "
      //     ": 0)",
      //     4 + inum * 4, 3 + inum * 4, 4 + inum * 4, 5 + inum * 4, 6 + inum *
      //     4);
      // parLimSPE.insert(parLimSPE.end(),
      //                  {0.1, 6, gaussMid - 5, gaussMid + 10, 1, 8, 0.1, 3 });
      gaussMid += 12;
    }
    // funcSPE += Form(" + [%d]*exp(-0.5*((x-[%d])/[%d])*((x-[%d])/[%d]))",
    //                 (numParSet), (numParSet + 1), (numParSet + 2),
    //                 (numParSet + 1), (numParSet + 2));
    // numParSet += 3;
    // parLimSPE.insert(parLimSPE.end(), {0.0, 80, 250, 400, 5, 100});
    std::cout << funcSPE << std::endl;
    WFAveraged->FitFunction(funcSPE, parLimSPE, start, wfSz - 5);
    // WFAveraged->Plot();
    WFAveraged->PrintFitParameters();

    std::vector<double> traceSPE;
    double baseline = WFAveraged->GetFitPar(0);
    for (int i = 0; i < wfSz; ++i) {
      traceSPE.push_back(WFAveraged->GetFitAt(i) - baseline);
    }
    digiAnalysis::WaveForm WFFit(traceSPE);
    WFFit.SetTracesFFT();
    std::vector<double> SPEFFT_Amp = WFFit.GetTracesFFT();
    std::vector<double> SPEFFT_Phase = WFFit.GetTracesFFTPhase();
    double intSPEVal = WFFit.IntegrateWaveForm(100, 200);
    std::cout << "Integrated Waveform Value: " << intSPEVal << std::endl;
    // WFFit.Plot(WFFit.GetTracesFFT(), WFAveraged->GetTracesFFT());
    // WFFit.Plot(WFFit.GetTracesFFT());
    // WFAveraged->Plot(WFAveraged->GetTracesFFT(), "SAME_kGreen");

    WFAveraged->Plot();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Plot the result of division of SPE waveform from average waveform

    // count = 0;
    // for (int iter = 0; iter < selWF.size() and count < 300; iter++) {
    //   std::vector<double> trace = selWF[iter].GetTraces();
    //   selWF[iter].SetTracesFFT();
    //   std::vector<double> trFFT_Amp =
    //       selWF[iter].GetTracesFFT(); // WFFit.GetTracesFFT(); //
    //   std::vector<double> trFFT_Phase =
    //       selWF[iter].GetTracesFFTPhase(); // WFFit.GetTracesFFTPhase();
    //                                        //
    //   // std::cout << trFFT_Amp.size() << std::endl;
    //   for (int iter1 = 0; iter1 < trFFT_Amp.size(); iter1++) {
    //     trFFT_Amp[iter1] /=
    //         (3 *
    //          SPEFFT_Amp[iter1]); // the constant value needs to be normalized
    //     trFFT_Amp[iter1] *= filter[iter1];
    //     // trFFT_Phase[iter] = trFFT_Phase[iter] - SPEFFT_Phase[iter];
    //   }

    //   std::vector<double> deconvolvedTrace =
    //       WFFit.EvalIFFT(trFFT_Amp, trFFT_Phase);
    //   // WFFit.Plot(deconvolvedTrace);
    //   if (count == 0)
    //     selWF[iter].Plot(deconvolvedTrace);
    //   if (count != 0 and count < 100)
    //     selWF[iter].Plot(deconvolvedTrace, "SAME");
    //   count++;
    // }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Use Average WF to evaluate the NPE of the SPE

    // digiAnalysis::WaveForm *highEWF;
    // double intWFVal, en;
    // TH2 *hENPE = new TH2F("hENPE", "E vs NPE", 4096, 0, 4096, 110, -10, 100);
    // // for (int wfiter = 0; wfiter < pairVecSz; wfiter++) {
    // //   intWFVal =
    // // vecOfPairs[wfiter]->GetHitPtr(0)->GetWFPtr()->IntegrateWaveForm();
    // //   en = vecOfPairs[wfiter]->GetPairHitEnergy(0) * 0.052966 - 4.547;
    // //   hENPE->Fill(en, intWFVal / intSPEVal / en);
    // // }
    // // TCanvas *c9 = new TCanvas("c9", "E vs NPE", 800, 600);
    // // hENPE->Draw("COLZ");

    // // for (int wfiter = 0; wfiter < pairVecSz; wfiter++) {
    // //   if (vecOfPairs[wfiter]->GetPairHitEnergy(0) * 0.052966 - 4.547 < -0
    // and
    // //       vecOfPairs[wfiter]->GetPairHitEnergy(0) * 0.052966 - 4.547 >
    // -10) {
    // //     vecOfPairs[wfiter]->Print();
    // //     highEWF = vecOfPairs[wfiter]->GetHitPtr(0)->GetWFPtr();
    // //     break;
    // //   }
    // // }

    // int wfiter = 0;
    // for (int wfiter = 0; wfiter < hitsVecSz; wfiter++) {
    //   intWFVal = hitsVector[wfiter]->GetWFPtr()->IntegrateWaveForm();
    //   en = hitsVector[wfiter]->GetEnergy() * 0.052966 - 4.547;
    //   hENPE->Fill(en, intWFVal / intSPEVal / en);
    // }
    // TCanvas *c9 = new TCanvas("c9", "E vs NPE", 800, 600);
    // hENPE->Draw("COLZ");

    // // for (int wfiter = 0; wfiter < hitsVecSz; wfiter++) {
    // //   if (hitsVector[wfiter]->GetEnergy() * 0.052966 - 4.547 < -0 and
    // //       hitsVector[wfiter]->GetEnergy() * 0.052966 - 4.547 > -10) {
    // //     hitsVector[wfiter]->Print();
    // //     highEWF = hitsVector[wfiter]->GetWFPtr();
    // //     break;
    // //   }
    // // }

    // highEWF = hitsVector[4]->GetWFPtr();
    // // pad the spe waveform to get same size as full waveform
    // std::vector<double> padtrace(highEWF->GetSize(), 0);
    // std::copy(traceSPE.begin(), traceSPE.end(), padtrace.begin());
    // WFAveraged->SetWaveForm(padtrace);
    // WFAveraged->SetTracesFFT();
    // SPEFFT_Amp = WFAveraged->GetTracesFFT();
    // SPEFFT_Phase = WFAveraged->GetTracesFFTPhase();
    // std::cout << "SPE FFT size: " << SPEFFT_Amp.size() << std::endl;
    // // SPEWF->Plot(SPEFFT_Amp);

    // // int evt = 105;
    // // std::vector<double> speAvTrace = WFVec[evt]->GetTraces();
    // // WFVec[evt]->SetBaseLine(speAvTrace, 400, 200);
    // // std::vector<double> spepadtrace(highEWF->GetSize(),
    // //                                 WFVec[evt]->GetBaseLine());
    // // std::copy(speAvTrace.begin(), speAvTrace.end(), spepadtrace.begin());
    // // WFVec[evt]->SetWaveForm(spepadtrace, 0, spepadtrace.size() - 1, 400,
    // // 200); WFVec[evt]->SetTracesFFT();

    // highEWF->SetTracesFFT();
    // std::vector<double> trFFT_Amp = highEWF->GetTracesFFT();
    // // WFAveraged->GetTracesFFT(); //WFVec[evt]->GetTracesFFT(); //
    // std::vector<double> trFFT_Phase = highEWF->GetTracesFFTPhase();
    // // WFAveraged->GetTracesFFTPhase(); //WFVec[evt]->GetTracesFFTPhase();
    // std::cout << "highE FFT size: " << trFFT_Amp.size() << std::endl;
    // // // highEWF->Plot();
    // for (int iter = 0; iter < trFFT_Amp.size(); iter++) {
    //   trFFT_Amp[iter] /= (0.6 * SPEFFT_Amp[iter]);
    //   trFFT_Amp[iter] *= filter[iter];
    //   trFFT_Phase[iter] = trFFT_Phase[iter] - SPEFFT_Phase[iter];
    // }
    // std::vector<double> deconvolvedTrace =
    //     highEWF->EvalIFFT(trFFT_Amp, trFFT_Phase);
    // std::cout << "trace 0: " << deconvolvedTrace[0] << std::endl;
    // double mean = std::accumulate(deconvolvedTrace.begin(),
    //                               deconvolvedTrace.begin() + 300, 0.0) /
    //               300;

    // std::cout << "mean: " << mean << std::endl;

    // for (auto &x : deconvolvedTrace)
    //   x -= mean;
    // std::cout << "trace 0: " << deconvolvedTrace[0] << std::endl;
    // // std::transform(deconvolvedTrace.begin(), deconvolvedTrace.end(),
    // //                deconvolvedTrace.begin(), [](double x) { return 5.0 *x;
    // //                });

    // // highEWF->Plot(deconvolvedTrace, trFFT_Amp);
    // // highEWF->SetSmooth(40);
    // highEWF->Plot(highEWF->GetTraces(), deconvolvedTrace);

    fApp->Run();
    return 0;
  } else {
    std::cout << "Could Not Read The Tree" << std::endl;
    return 1;
  }
}