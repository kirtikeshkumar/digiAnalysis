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

  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //     "NaI1_21May26_1900_Cs_WAVES_2/FILTERED/"
  //     "DataF_NaI1_21May26_1900_Cs_WAVES_2_BLCorrected.root";

  //   std::string fname =
  //       "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //       "NaI31_26May26_1345_1750_NoSrc_Thresh50_WAVES_NoCoinc_LeadPit_5/FILTERED/"
  //       "DataF_NaI31_26May26_1345_1750_NoSrc_Thresh50_WAVES_NoCoinc_LeadPit_5_"
  //       "BLCorrected.root";

  //   std::string fname =
  //       "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //       "NaI31_01June26_1345_1750_Cs_Thresh_300_30_WAVES_Coinc_144ns_LeadPit/"
  //       "FILTERED/"
  //       "SDataF_NaI31_01June26_1345_1750_Cs_Thresh_300_30_WAVES_Coinc_144ns_"
  //       "LeadPit_BLCorrected.root";

  //   std::string fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
  //                       "CoincidenceStudies/01JuneNoSrc/"
  //                       "NaI1342_June26_1750_1345_1350_1350_NoSrc_Thresh_2_30_"
  //                       "300_WAVES_Coinc_144ns_LeadPit_Sum_BLCorrected.root";

  //   std::string fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
  //                       "CoincidenceStudies/01JuneNoSrc/"
  //                       "NaI1342_04June26_1750_1345_1350_1350_NoSrc_Thresh_120_"
  //                       "300_WAVES_Singles_LeadPit_45/FILTERED/"
  //                       "DataF_NaI1342_04June26_1750_1345_1350_1350_NoSrc_Thresh_"
  //                       "120_300_WAVES_Singles_LeadPit_45_BLCorrected.root";

  //   std::string fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
  //                       "CoincidenceStudies/01JuneNoSrc/"
  //                       "NaI1342_June26_1750_1345_1350_1350_NoSrc_Thresh_15-30_"
  //                       "300_WAVES_Coinc_144ns_LeadPit_Sum_BLCorrected.root";

  std::string fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
                      "CoincidenceStudies/01JuneNoSrc/CalibrationFiles/"
                      "DataF_NaI1_05June26_1750_CsSrc_Thresh_120_300_WAVES_"
                      "Singles_LeadPit_68_0-500k_BLCorrected.root";

  // Read to singleHits
  digiAnalysis::Analysis an(fname, 0, 200000, 0);

  // Get the vector
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
#ifdef WAVES
  int spectralsize = 500; // 8192
  double WFMT = 0;        // MeanTime parameter for flat waveform
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
  double Q1 = 0, Q2 = 0, avT1 = 0, avT2 = 0, netQ = 0;
  double newLam = 0;
  bool keepGoing = true;
  std::string userInput;
  digiAnalysis::WaveForm *WF = nullptr;
  std::vector<double> trFFT;
  std::vector<double> trFFT_Phase;
  std::vector<digiAnalysis::WaveForm> waveformVector;
  double meanTime;
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
      Q1 = 0;
      Q2 = 0;
      avT1 = 0;
      avT2 = 0;
      netQ = 0;
      std::vector<double> trace = WF->GetTraces();
      int sz = 700;       // digiAnalysis::GateLenLong; // trace.size();
      int startVal = 450; // digiAnalysis::GateStart
      for (int iter = startVal; iter < startVal + sz / 2; iter++) {
        avT1 += abs(trace[iter]) * iter;
        avT2 += abs(trace[iter + sz / 2]) * (iter + sz / 2);
        Q1 += abs(trace[iter]);
        Q2 += abs(trace[iter + sz / 2]);
        netQ += (trace[iter] + trace[iter + sz / 2]);
      }
      newLam =
          TMath::Log(-1.0 * TMath::Log(Q2 / Q1) /
                     (avT2 / Q2 - avT1 / Q1)); // /(avT2 / avQ2 - avT1 / avQ1)

      evalEnergy =
          hitsVector[i]->GetEvalEnergy(); // WF->IntegrateWaveForm(290, 1390);
      evalEnergyShort =
          hitsVector[i]->GetEvalEnergyShort(); // WF->IntegrateWaveForm(290,
                                               // 440);
      psd = 1.0 - energyShort * 1.0 / energy;
      meanTime = hitsVector[i]->GetMeanTime();
      energy = energy * 0.09032 - 3.3849;
      evalEnergy = evalEnergy * 0.09032 - 3.3849;
      hLamPlot->Fill(evalEnergy, newLam);
      hLamMTPlot->Fill(evalEnergy, newLam, meanTime);
      hMTPlot->Fill(energy, meanTime);
      hMTLam->Fill(meanTime, newLam);
      // shortPSD =
      //     WF->IntegrateWaveForm(440, 600) / WF->IntegrateWaveForm(290, 600);
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
      if (meanTime > 2.2) {
        hESpectra->Fill(energy);
        hEEvalSpectra->Fill(evalEnergy);
      }

      // ########################################################### //
      //             This part plots selected waveforms              //
      // ########################################################### //

      //   energy = energy * 0.09032 - 3.3849;
      //   evalEnergy = evalEnergy * 0.09032 - 3.3849;
      //   if (keepGoing and
      //       // newLam > -6.5 and newLam < -6 and
      //       energy > 0 and energy < 2 and newLam > -6 and psd > 0 and
      //       meanTime < 2. and meanTime > 0.) {
      //     hitsVector[i]->Print();
      //     std::cout << "lam: " << newLam << " : " << Q1 << " : " << Q2 << " :
      //     "
      //               << avT1 << " : " << avT2 << " : " << netQ << " : "
      //               << WF->IntegrateWaveForm(1500, 4900) / 3400.0 <<
      //               std::endl;
      //     std::cout << "Energy: " << energy << " : " << evalEnergy
      //               << "\t | PSD = " << psd << std::endl
      //               << std::endl;
      //     WF->SetSmooth(16, "MovA");
      //     waveformVector.push_back(*WF);
      //     // WF->SetTracesFFT();
      //     // trFFT.clear();
      //     // trFFT_Phase.clear();
      //     // trFFT = WF->GetTracesFFT();
      //     // trFFT_Phase = WF->GetTracesFFTPhase();
      //     // int cutoff = 50;
      //     // std::fill(trFFT.begin() + cutoff, trFFT.end(), 0.0);
      //     // std::fill(trFFT_Phase.begin() + cutoff, trFFT_Phase.end(), 0.0);
      //     // WF->Plot(WF->GetTracesSmooth(), WF->EvalIFFT(trFFT,
      //     //                                              trFFT_Phase)); //
      //     WF->Plot();
      //     std::cout << "Do you want to see the next waveform? (y/n): ";
      //     std::getline(std::cin, userInput);
      //     if (userInput != "y" && userInput != "Y") {
      //       keepGoing = false;
      //     }
      //   }
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
  TCanvas *c5 = new TCanvas("c5", "PSD vs PSDEval ", 800, 600);
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
  TCanvas *c12 = new TCanvas("c12", "PSD vs Lam vs MT", 800, 600);
  hPSDLamMTPlot->Draw("");
  TCanvas *c13 = new TCanvas("c13", "PSD vs Lam", 800, 600);
  hPSDLamPlot->Draw("COLZ");
  TCanvas *c14 = new TCanvas("c14", "Energy vs PSD vs Lam", 800, 600);
  hPSDLamEPlot->Draw("");

  UShort_t wfSz = WF->GetSize();
  std::cout << " Averaging " << waveformVector.size()
            << " waveforms for plotting" << std::endl;
  digiAnalysis::WaveForm WFAveraged(wfSz, waveformVector);
  WFAveraged.SetSmooth(40);
  WFAveraged.Plot();

  fApp->Run();
#endif
  return 0;
}
