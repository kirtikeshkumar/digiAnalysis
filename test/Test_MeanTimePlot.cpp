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
#include <RtypesCore.h>
#include <TApplication.h>
#include <iostream>
#include <string>
#include <vector>
int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/extCoinc/"
      "NaI3124_17Jun26_NoSrc_1350V_2000V_1350V_1350V_Gain2_NoSplit_ExtTrig_"
      "Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_1000nsCoinc_"
      "2Vpp_Thresh_100lsb_WAVES_Sum_BLCorrected.root";

  // Read to singleHits
  digiAnalysis::Analysis an(2, fname, 0, 20000, 0);

  // Get the vector
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  std::cout << "Got " << hitsVector.size() << " hits from file" << std::endl;
#ifdef WAVES
  int spectralsize = 80; // 8192
  double WFMT = 0;       // MeanTime parameter for flat waveform
  TH2 *hMTPlot = new TH2F("MTPlot", "Energy vs MeanTime", spectralsize, 0,
                          spectralsize, 500, -4, 4);
  TH2 *hLamPlot = new TH2F("LamPlot", "Energy vs Lambda", spectralsize, 0,
                           spectralsize, 2000, -10.0, 10.0);
  TH2 *hPSDLamPlot =
      new TH2F("PSDLamPlot", "PSD vs Lambda", 100, 0, 1, 2000, -10.0, 10.0);
  TH3 *hLamMTPlot =
      new TH3F("LamMTPlot", "Energy vs Lambda vs MT", spectralsize, 0,
               spectralsize, 80, -10, -2, 200, 0, 4);
  TH3 *hPSDLamMTPlot = new TH3F("PSDLamMTPlot", "PSD vs Lambda vs MT", 200, -1,
                                1, 60, -10, -4, 200, -4, 4);
  TH3 *hPSDLamEPlot =
      new TH3F("PSDLamEPlot", "PSD vs Lambda vs E", spectralsize, 0,
               spectralsize, 100, 0, 1, 80, -10, -2);
  TH2 *hMTLam =
      new TH2F("MTLam", "MeanTime vs Lambda", 500, -4, 4, 2000, -10.0, 10.0);
  TH2 *hPSDPlot = new TH2F("PSDPlot", "Energy vs PSD", spectralsize, 0,
                           spectralsize, 500, -8, 2);
  TH2 *hEPlot = new TH2F("EPlot", "Energy vs EvalEnergy", 4096, 0, spectralsize,
                         spectralsize, 0, spectralsize);
  TH2 *hESPlot = new TH2F("ESPlot", "EnergyShort vs EvalEnergyShort", 410, 0,
                          spectralsize, 410, 0, spectralsize);
  TH2 *hPSDEvalPlot =
      new TH2F("PSDEvalPlot", "PSD vs PSDEval", 160, -8, 8, 160, -8, 8);
  TH2 *hMTPSD = new TH2F("MTPSD", "MT vs PSD", 500, -4, 4, 40, -2, 2);
  TH1F *hEEvalRatio =
      new TH1F("hEEValRatio", "Ratio of EEval vs E", 1000, 0, 200);
  TH2 *hEdiffPlot = new TH2F("EdiffPlot", "Energy vs Energy-EnergyShort", 410,
                             0, spectralsize, 820, 0.0, 10);
  TH2 *hEdiffEvalPlot =
      new TH2F("EdiffEvalPlot", "EvalEnergy vs EvalEnergy-EvalEnergyShort", 410,
               0, spectralsize, 820, -0.0, 10.0);
  TH1F *hESpectra = new TH1F("hESpectra", "Energy Spectra", 5 * spectralsize, 0,
                             spectralsize);
  TH1F *hEEvalSpectra = new TH1F("hEEvalSpectra", "Eval Energy Spectra",
                                 5 * spectralsize, 0, spectralsize);
  TH2 *hELL = new TH2F("hELL", "hELL", 800, 0, 80, 1000, -5, 5);
  TH2 *hMTLL = new TH2F("hMTLL", "hMTLL", 500, -4, 4, 1000, -5, 5);
  TH2 *hX1X2 = new TH2F("hX1X2", "hX1X2", 60, -0.2, 1, 50, -1, 1);
  int nentries = hitsVector.size();
  double psd = 0;
  double evalEnergy = 0;
  double evalEnergyShort = 0;
  double energy = 0;
  double energyShort = 0;
  double shortPSD = 0;
  double Q1 = 0, Q2 = 0, avT1 = 0, avT2 = 0, netQ = 0;
  double X1 = 0, X2 = 0, XTot = 0;
  double newLam = 0;
  bool keepGoing = true;
  std::string userInput;
  digiAnalysis::WaveForm *WF = nullptr;
  std::vector<double> trFFT;
  std::vector<double> trFFT_Phase;
  std::vector<digiAnalysis::WaveForm> waveformVector;
  double meanTime;
  int badwf = 0;

  // for loglikelihood
  TFile f("output.root", "READ");
  std::vector<double> *traceAv = nullptr;
  f.GetObject("trace", traceAv);
  digiAnalysis::WaveForm WFNormAv(*traceAv);
  double intValAv = WFNormAv.IntegrateWaveForm();

  double X1LowCut = 0.06, X1HiCut = 0.16;
  double X2LowCut = 0.4, X2HiCut = 0.6;

  for (int i = 0; i < nentries; i++) {
    if (i % 10000 == 0) {
      std::cout << i << std::endl;
    }
    if (hitsVector[i]->GetChNum() ==
        2 // and
          // hitsVector[i]->GetTimestamp() / 1E12 > 1600 and
          // hitsVector[i]->GetTimestamp() / 1E12 < 3200
          // and hitsVector[i]->GetMeanTime() < 3.8
    ) {   // (hitsVector[i]->GetPSD() > 0.0 and
      // hitsVector[i]->GetChNum() == 0) {
      energy = hitsVector[i]->GetEvalEnergy() * 0.01911 -
               0.301; // * 0.052966 - 4.547;
      energyShort = hitsVector[i]->GetEvalEnergyShort() * 0.01911 -
                    0.301; // * 0.052966 - 4.547;
      WF = hitsVector[i]->GetWFPtr();
      WF->SetSmooth(100);
      std::vector<double> trace = WF->GetTracesSmooth();
      // WF->SetTracesMovBLCorr();
      // WF->SetMeanTime();

      double preInt = WF->IntegrateWaveForm(0, digiAnalysis::GateStart);

      WFMT = TMath::Log10(WF->GetSize() / 2.0);
      newLam = -1;
      if (preInt < 2.0) {
        Q1 = 0;
        Q2 = 0;
        X1 = 0;
        X2 = 0;
        XTot = 0;
        avT1 = 0;
        avT2 = 0;
        netQ = 0;
        int sz = 700;        // digiAnalysis::GateLenLong; // trace.size();
        int startVal = 1080; // digiAnalysis::GateStart
        int startVal1 = 980;
        for (int iter = startVal; iter < startVal + sz / 2; iter++) {
          avT1 += abs(trace[iter]) * iter;
          avT2 += abs(trace[iter + sz / 2]) * (iter + sz / 2);
          Q1 += abs(trace[iter]);
          Q2 += abs(trace[iter + sz / 2]);
          netQ += (trace[iter] + trace[iter + sz / 2]);
        }
        for (int iter = startVal1; iter < startVal1 + 301; iter++) {
          if (iter - startVal1 < 25) {
            X1 += trace[iter];
          }
          if (iter - startVal1 > 80 and iter - startVal1 < 300) {
            X2 += trace[iter];
          }
          if (iter - startVal1 < 300) {
            XTot += trace[iter];
          }
        }
        X1 = X1 / XTot;
        X2 = X2 / XTot;
      }
      newLam = TMath::Log(-1.0 * TMath::Log(Q2 / Q1) / (avT2 / Q2 - avT1 / Q1));
      evalEnergy = hitsVector[i]->GetEvalEnergy() * 0.01911 -
                   0.301; // WF->IntegrateWaveForm(290, 1390);
      evalEnergyShort = hitsVector[i]->GetEvalEnergyShort() * 0.01911 -
                        0.301; // WF->IntegrateWaveForm(290,
                               // 440);
      psd = 1.0 - energyShort * 1.0 / energy;
      meanTime = hitsVector[i]->GetMeanTime();
      // energy = energy * 0.09032 - 3.3849;
      // evalEnergy = evalEnergy * 0.09032 - 3.3849;
      if (psd > 0.2 and psd < 0.5 and preInt < 2.0) {
        hX1X2->Fill(X1, X2);
        hESpectra->Fill(energy);
        hLamPlot->Fill(evalEnergy, newLam);
      }

      if (psd > 0.2 and psd < 0.5 and preInt < 2.0 and X1 > X1LowCut and
          X1 < X1HiCut and X2 < X2HiCut and X2 > X2LowCut) {
        // if (meanTime > 1.8 and meanTime < 2.1)
        //   hLamPlot->Fill(evalEnergy, newLam);
        hMTPlot->Fill(energy, meanTime);
        hLamMTPlot->Fill(evalEnergy, newLam, meanTime);
        hMTLam->Fill(meanTime, newLam);
        // shortPSD =
        //     WF->IntegrateWaveForm(440, 600) / WF->IntegrateWaveForm(290,
        //     600);
        hPSDLamMTPlot->Fill(psd, newLam, meanTime);
        hPSDLamEPlot->Fill(evalEnergy, psd, newLam);
        hPSDLamPlot->Fill(psd, newLam);
        hPSDPlot->Fill(energy, psd);
        hEPlot->Fill(energy, evalEnergy);
        hESPlot->Fill(energyShort, evalEnergyShort);
        hPSDEvalPlot->Fill(hitsVector[i]->GetPSD(), psd);
        hEEvalRatio->Fill(evalEnergyShort * 1.0 / energyShort);
        hEdiffPlot->Fill(energy, energyShort / energy);
        hEdiffEvalPlot->Fill(evalEnergy, evalEnergyShort / evalEnergy);
        hMTPSD->Fill(meanTime, psd);
      }

      if (meanTime > 1.8 and meanTime < 2.1 and psd < 0.5 and psd > 0.01 and
          preInt < 2.0 and X1 > X1LowCut and X1 < X1HiCut and X2 < X2HiCut and
          X2 > X2LowCut) { // newLam > -5.4 and newLam < -4.5 and
        hEEvalSpectra->Fill(evalEnergy);
      }
      if (preInt >= 2.0) {
        badwf++;
      }

      double logLikelihood = 0;
      // ########################################################### //
      // This part evaluates the loglikelihood
      // ########################################################### //

      // if (newLam > -5.5 and newLam < -4.5 and psd < 0.5 and psd > 0.12 and
      //     meanTime < 2.1 and meanTime > 1.9 and preInt < 2.0) {
      //   digiAnalysis::WaveForm WFNorm(WF->ScaleWaveForm(1.0 / fabs(energy)));
      //   trace = WFNorm.GetTraces();

      //   for (int j = digiAnalysis::GateStart;
      //        j < digiAnalysis::GateStart + digiAnalysis::GateLenLong; j++) {
      //     logLikelihood += WFNorm.GetTraces()[j] *
      //                      TMath::Log(abs(trace[j] / traceAv->at(j)));
      //   }
      //   logLikelihood = logLikelihood + intValAv -
      //   WFNorm.IntegrateWaveForm(); hELL->Fill(energy, logLikelihood /
      //   intValAv); hMTLL->Fill(meanTime, logLikelihood / intValAv);
      // }
      // ########################################################### //

      // ########################################################### //
      //             This part plots selected waveforms              //
      // ########################################################### //

      // energy = energy * 0.09032 - 3.3849;
      // evalEnergy = evalEnergy * 0.09032 - 3.3849;
      // if (keepGoing and
      //     // newLam > -6.5 and newLam < -6 and
      //     energy > 1 and energy < 5 and newLam > -5.5 and newLam < -4.5 and
      //     psd < 0.5 and psd > 0.2 and meanTime < 2.1 and meanTime > 1.9 and
      //     preInt < 2.0 and X1 > X1LowCut and X1 < X1HiCut and X2 < X2HiCut
      //     and X2 > X2LowCut) { // and logLikelihood / 10000 < 0.3 and
      //     logLikelihood
      //                      // >0.
      //   hitsVector[i]->Print();
      //   std::cout << "lam: " << newLam << " : " << Q1 << " : " << Q2 << " : "
      //             << avT1 << " : " << avT2 << " : " << netQ << " : "
      //             << WF->IntegrateWaveForm(1500, 4900) / 3400.0 << std::endl;
      //   std::cout << "Energy: " << energy << " : " << evalEnergy
      //             << "\t | PSD = " << psd << " | InitEnergy: "
      //             << WF->IntegrateWaveForm(0, digiAnalysis::GateStart) * 1.0
      //             /
      //                    digiAnalysis::GateStart
      //             << std::endl
      //             << "X1: " << X1 << " : X2: " << X2 << std::endl
      //             << std::endl;
      //   std::cout << "log likelihood: " << logLikelihood << std::endl;
      //   WF->SetSmooth(16, "MovA");
      //   waveformVector.push_back(*WF);
      //   // WF->SetTracesFFT();
      //   // trFFT.clear();
      //   // trFFT_Phase.clear();
      //   // trFFT = WF->GetTracesFFT();
      //   // trFFT_Phase = WF->GetTracesFFTPhase();
      //   // int cutoff = 50;
      //   // std::fill(trFFT.begin() + cutoff, trFFT.end(), 0.0);
      //   // std::fill(trFFT_Phase.begin() + cutoff, trFFT_Phase.end(), 0.0);
      //   // WF->Plot(WF->GetTracesSmooth(), WF->EvalIFFT(trFFT,
      //   //                                              trFFT_Phase)); //
      //   WF->Plot();
      //   std::cout << "Energy: " << energy << " LL: " << logLikelihood
      //             << std::endl;
      //   std::cout << "Do you want to see the next waveform? (y/n): ";
      //   std::getline(std::cin, userInput);
      //   if (userInput != "y" && userInput != "Y") {
      //     keepGoing = false;
      //   }
      // }
      // ########################################################### //
    }
  }

  std::cout << "badWF: " << badwf << std::endl;
  TCanvas *c1 = new TCanvas("c1", "Energy vs MeanTime", 800, 600);
  hMTPlot->Draw("COLZ");
  TCanvas *c2 = new TCanvas("c2", "Energy vs PSD", 800, 600);
  hPSDPlot->Draw("COLZ");
  // TCanvas *c3 = new TCanvas("c3", "Energy vs evalEnergy", 800, 600);
  // hEPlot->Draw("COLZ");
  // hEEvalRatio->Draw("HIST");
  // TCanvas *c4 = new TCanvas("c4", "EnergyShort vs evalEnergyShort", 800,
  // 600); hESPlot->Draw("COLZ"); TCanvas *c5 = new TCanvas("c5", "PSD vs
  // PSDEval ", 800, 600); hPSDEvalPlot->Draw("COLZ"); TCanvas *c6 = new
  // TCanvas("c6", "Energy vs EnergyShort/Energy", 800, 600);
  // hEdiffPlot->Draw("COLZ");
  // TCanvas *c7 =
  //     new TCanvas("c7", "evalEnergy vs evalEnergyShort/evalEnergy", 800,
  //     600);
  // hEdiffEvalPlot->Draw("COLZ");
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
  TCanvas *c12 = new TCanvas("c12", "PSD vs Lam vs MT", 800, 600);
  hPSDLamMTPlot->Draw("");
  TCanvas *c13 = new TCanvas("c13", "PSD vs Lam", 800, 600);
  hPSDLamPlot->Draw("COLZ");
  // TCanvas *c14 = new TCanvas("c14", "Energy vs PSD vs Lam", 800, 600);
  // hPSDLamEPlot->Draw("");
  TCanvas *c15 = new TCanvas("c15", "MT vs PSD", 800, 600);
  hMTPSD->Draw("COLZ");
  TCanvas *c16 = new TCanvas("c16", "E vs LL", 800, 600);
  hELL->Draw("COLZ");
  TCanvas *c17 = new TCanvas("c17", "MT vs LL", 800, 600);
  hMTLL->Draw("COLZ");
  TCanvas *c18 = new TCanvas("c18", "X1 vs X2", 800, 600);
  hX1X2->Draw("LEGO");

  UShort_t wfSz = WF->GetSize();
  std::cout << " Averaging " << waveformVector.size()
            << " waveforms for plotting" << std::endl;
  digiAnalysis::WaveForm WFAveraged(wfSz, waveformVector);
  WFAveraged.SetSmooth(40);
  WFAveraged.SetTracesFFT();
  WFAveraged.Plot();

  fApp->Run();
#endif
  return 0;
}
