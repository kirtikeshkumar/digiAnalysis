/*
**	Filename : Test_Dummy.cpp
**	2024-12-03
**	username : kirtikeshkumar
*/

#include "Analysis.h"
#include "TMath.h"
#include "WaveForm.h"
#include "globals.h"
#include "includes.hh"
#include "singleHits.h"
#include <TApplication.h>
#include <iostream>
#include <string>
#include <vector>
int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);

  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/NaI_1_Checks/Calibration/"
  //     "DataF_NaI_1_Na_Source_Gain_Calibration_HV1900_Waves_160FC.root";

  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //     "NaI13_20May26_1900_1345_Cs_Coinc144_WAVES_2/FILTERED/"
  //     "SDataF_NaI13_20May26_1900_1345_Cs_Coinc144_WAVES_2_BLCorrected.root";

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
      "NaI1_21May26_1900_Cs_WAVES_2/FILTERED/"
      "DataF_NaI1_21May26_1900_Cs_WAVES_2_BLCorrected.root";

  // Read to singleHits
  digiAnalysis::Analysis an(0, fname, 0, 000, 0);

  // Get the vector
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
#ifdef WAVES
  int spectralsize = 4096; // 8192
  double WFMT = 0;         // MeanTime parameter for flat waveform
  TH2 *hMTPlot = new TH2F("MTPlot", "Energy vs MeanTime", spectralsize, 0,
                          spectralsize, 500, -4, 4);
  TH2 *hLamPlot = new TH2F("LamPlot", "Energy vs Lambda", spectralsize, 0,
                           spectralsize, 2000, -10.0, 10.0);
  TH3 *hLamMTPlot = new TH3F("LamMTPlot", "Energy vs Lambda vs MT", 128, 0,
                             4096, 80, -10, -2, 50, 2.2, 3.2);
  TH2 *hMTLam =
      new TH2F("MTLam", "MeanTime vs Lambda", 500, -4, 4, 2000, -10.0, 10.0);
  TH2 *hPSDPlot = new TH2F("PSDPlot", "Energy vs PSD", spectralsize, 0,
                           spectralsize, 500, -8, 2);
  TH2 *hEPlot = new TH2F("EPlot", "Energy vs EvalEnergy", 4096, 0, spectralsize,
                         4096, 0, spectralsize);
  TH2 *hESPlot = new TH2F("ESPlot", "EnergyShort vs EvalEnergyShort", 410, 0,
                          spectralsize, 410, 0, spectralsize);
  TH2 *hPSDEvalPlot =
      new TH2F("PSDEvalPlot", "PSD vs PSDEval", 160, -8, 8, 160, -8, 8);
  TH1F *hEEvalRatio =
      new TH1F("hEEValRatio", "Ratio of EEval vs E", 1000, 0, 200);
  TH2 *hEdiffPlot = new TH2F("EdiffPlot", "Energy vs Energy-EnergyShort", 410,
                             0, spectralsize, 820, 0.0, 10);
  TH2 *hEdiffEvalPlot =
      new TH2F("EdiffEvalPlot", "EvalEnergy vs EvalEnergy-EvalEnergyShort", 410,
               0, spectralsize, 820, -0.0, 10.0);
  TH1F *hESpectra =
      new TH1F("hESpectra", "Energy Spectra", spectralsize, 0, spectralsize);
  TH1F *hEEvalSpectra = new TH1F("hEEvalSpectra", "Eval Energy Spectra",
                                 spectralsize, 0, spectralsize);
  int nentries = hitsVector.size();
  double psd = 0;
  double evalEnergy = 0;
  double evalEnergyShort = 0;
  double energy = 0;
  double energyShort = 0;
  double shortPSD = 0;
  double avQ1 = 0, avQ2 = 0, avT1 = 0, avT2 = 0, netQ = 0;
  double newLam = 0;
  bool keepGoing = true;
  std::string userInput;
  digiAnalysis::WaveForm *WF = nullptr;
  std::vector<double> trFFT;
  std::vector<double> trFFT_Phase;
  std::vector<digiAnalysis::WaveForm> waveformVector;
  for (int i = 0; i < nentries; i++) {
    if (i % 10000 == 0) {
      std::cout << i << std::endl;
    }
    if (hitsVector[i]->GetChNum() ==
        0 // and
          // hitsVector[i]->GetTimestamp() / 1E12 > 1600 and
          // hitsVector[i]->GetTimestamp() / 1E12 < 3200
          // and hitsVector[i]->GetMeanTime() < 3.8
    ) {   // (hitsVector[i]->GetPSD() > 0.0 and
      // hitsVector[i]->GetChNum() == 0) {
      energy = hitsVector[i]->GetEnergy();           // * 0.052966 - 4.547;
      energyShort = hitsVector[i]->GetEnergyShort(); // * 0.052966 - 4.547;
      WF = hitsVector[i]->GetWFPtr();
      // WF->SetTracesMovBLCorr();
      // WF->SetMeanTime();
      WFMT = TMath::Log10(WF->GetSize() / 2.0);
      avQ1 = 0;
      avQ2 = 0;
      avT1 = 0;
      avT2 = 0;
      netQ = 0;
      std::vector<double> trace = WF->GetTraces();
      int sz = digiAnalysis::GateLenLong; // trace.size();
      for (int iter = digiAnalysis::GateStart;
           iter < digiAnalysis::GateStart + sz / 2; iter++) {
        avT1 += trace[iter] * iter;
        avT2 += trace[iter + sz / 2] * (iter + sz / 2);
        avQ1 += trace[iter];
        avQ2 += trace[iter + sz / 2];
        netQ += (trace[iter] + trace[iter + sz / 2]);
      }
      newLam = TMath::Log(
          -1.0 * TMath::Log(avQ2 / avQ1) /
          (avT2 / avQ2 - avT1 / avQ1)); // /(avT2 / avQ2 - avT1 / avQ1)

      evalEnergy =
          hitsVector[i]->GetEvalEnergy(); // WF->IntegrateWaveForm(290, 1390);
      evalEnergyShort =
          hitsVector[i]->GetEvalEnergyShort(); // WF->IntegrateWaveForm(290,
                                               // 440);
      hLamPlot->Fill(evalEnergy, newLam);
      hLamMTPlot->Fill(evalEnergy, newLam, hitsVector[i]->GetMeanTime());
      hMTPlot->Fill(evalEnergy, hitsVector[i]->GetMeanTime());
      hMTLam->Fill(hitsVector[i]->GetMeanTime(), newLam);
      // shortPSD =
      //     WF->IntegrateWaveForm(440, 600) / WF->IntegrateWaveForm(290, 600);
      psd = 1.0 - evalEnergyShort * 1.0 / evalEnergy;
      hPSDPlot->Fill(energy, psd);
      hEPlot->Fill(energy, evalEnergy);
      hESPlot->Fill(energyShort, evalEnergyShort);
      hPSDEvalPlot->Fill(hitsVector[i]->GetPSD(), psd);
      hEEvalRatio->Fill(evalEnergyShort * 1.0 / energyShort);
      hEdiffPlot->Fill(energy, energyShort / energy);
      hEdiffEvalPlot->Fill(evalEnergy, evalEnergyShort / evalEnergy);
      hESpectra->Fill(energy);
      hEEvalSpectra->Fill(evalEnergy);

      // ########################################################### //
      //             This part plots selected waveforms              //
      // ########################################################### //

      if (keepGoing and
          // newLam > -6.5 and newLam < -6 and
          energy > 600 and
          energy < 710 // and
                       // hitsVector[i]->GetMeanTime() < 2.5 and
                       // hitsVector[i]->GetMeanTime() > 2.3
      ) {
        hitsVector[i]->Print();
        std::cout << "lam: " << newLam << " : " << avQ1 << " : " << avQ2
                  << " : " << avT1 << " : " << avT2 << " : " << netQ
                  << std::endl;
        std::cout << "Energy: " << energy << "\t | Lambda = " << newLam
                  << std::endl
                  << std::endl;
        WF->SetSmooth(16, "MovA");
        waveformVector.push_back(*WF);
        // WF->SetTracesFFT();
        // trFFT.clear();
        // trFFT_Phase.clear();
        // trFFT = WF->GetTracesFFT();
        // trFFT_Phase = WF->GetTracesFFTPhase();
        // int cutoff = 50;
        // std::fill(trFFT.begin() + cutoff, trFFT.end(), 0.0);
        // std::fill(trFFT_Phase.begin() + cutoff, trFFT_Phase.end(), 0.0);
        // WF->Plot(WF->GetTracesSmooth(), WF->EvalIFFT(trFFT,
        //                                              trFFT_Phase)); //
        WF->Plot(WF->GetTracesSmooth(), WF->GetTraces());
        std::cout << "Do you want to see the next waveform? (y/n): ";
        std::getline(std::cin, userInput);
        if (userInput != "y" && userInput != "Y") {
          keepGoing = false;
        }
      }
      // ########################################################### //
    }
  }
  TCanvas *c1 = new TCanvas("c1", "Energy vs MeanTime", 800, 600);
  hMTPlot->Draw("COLZ");
  TCanvas *c2 = new TCanvas("c2", "Energy vs PSD", 800, 600);
  hPSDPlot->Draw("COLZ");
  TCanvas *c3 = new TCanvas("c3", "Energy vs evalEnergy", 800, 600);
  hEPlot->Draw("COLZ");
  // hEEvalRatio->Draw("HIST");
  TCanvas *c4 = new TCanvas("c4", "EnergyShort vs evalEnergyShort", 800, 600);
  hESPlot->Draw("COLZ");
  TCanvas *c5 = new TCanvas("c5", "PSD vs  PSDEval", 800, 600);
  hPSDEvalPlot->Draw("COLZ");
  TCanvas *c6 = new TCanvas("c6", "Energy vs EnergyShort/Energy", 800, 600);
  hEdiffPlot->Draw("COLZ");
  TCanvas *c7 =
      new TCanvas("c7", "evalEnergy vs evalEnergyShort/evalEnergy", 800, 600);
  hEdiffEvalPlot->Draw("COLZ");
  TCanvas *c8 = new TCanvas("c8", "Energy Spectra", 800, 600);
  hESpectra->SetLineColor(kRed);
  hEEvalSpectra->SetLineColor(kGreen);
  hESpectra->Draw("HIST");
  hEEvalSpectra->Draw("HISTSAME");
  TCanvas *c9 = new TCanvas("c9", "Energy vs Lam", 800, 600);
  hLamPlot->Draw("COLZ");
  TCanvas *c10 = new TCanvas("c10", "MT vs Lam", 800, 600);
  hMTLam->Draw("COLZ");
  TCanvas *c11 = new TCanvas("c11", "E vs Lam vs MT", 800, 600);
  hLamMTPlot->Draw("");

  UShort_t wfSz = WF->GetSize();
  digiAnalysis::WaveForm WFAveraged(wfSz, waveformVector);
  WFAveraged.SetSmooth(40);
  WFAveraged.Plot();

  fApp->Run();
#endif
  return 0;
}
