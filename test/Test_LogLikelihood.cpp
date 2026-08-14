/*
Copied from Test_MultiChannelCoincidence for the single file part
Noise WF is taken from singles file by implementing X1X2 cut
Signal is taken from coincidence file
*/
#include "Analysis.h"
#include "TMath.h"
#include "WaveForm.h"
#include "globals.h"
#include "includes.hh"
#include "singleHits.h"
#include <Rtypes.h>
#include <TApplication.h>
#include <TH1.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>
int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);

  bool keepGoing = true;
  std::string userInput;

  std::vector<double> cluster0;
  std::vector<double> cluster1;
  double energy0 = 0.;

  double X1LowCut = 0.19, X1HiCut = 0.4;
  double mX = -1.096, mC_low = 0.764, mC_hi = 0.834;
  double X1 = 0, X2 = 0, XTot = 0;
  digiAnalysis::WaveForm *WF = nullptr, *WFNorm = nullptr;
  std::vector<digiAnalysis::WaveForm> waveformVector, wfVecNoise;

  TH2 *hEMT = new TH2I("hEMT", "hEMT", 800, 0, 80, 500, -4, 4);
  TH2 *hX1X2 = new TH2F("hX1X2", "hX1X2", 600, -0.2, 1, 500, -1, 1);

  // Defining the signal Template as cluster 1 from the X1X2 and meantime cut in
  // the pair data

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/extCoinc/PairFiles/"
      "Pair_NaI3124_15-17Jun26_NoSrc_1350V_2000V_1350V_1350V_Gain2_"
      "NoSplit_"
      "ExtTrig_"
      "Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_1000nsCoinc_"
      "2Vpp_Thresh_100lsb_WAVES_Sum_BLCorrected.root";
  digiAnalysis::Analysis an(fname, 0, 100000, 0);
  std::vector<std::unique_ptr<digiAnalysis::Pair>> &vecOfPairs =
      an.GetPairsVec();
  std::cout << "Number of pairs found: " << vecOfPairs.size() << std::endl;

  int badWFCount = 0;
  double MT, Lam;
  double preint;

  for (int i = 0; i < vecOfPairs.size(); i++) {
    X1 = 0;
    X2 = 0;
    XTot = 0;
    WF = vecOfPairs[i]->GetPairHitCh(0) == 2
             ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()
             : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr();
    preint = WF->IntegrateWaveForm(0, digiAnalysis::GateStart);
    preint /= digiAnalysis::GateStart;
    if (preint < 2.0) {
      X1 = WF->IntegrateWaveForm(980, 1025);
      X2 = WF->IntegrateWaveForm(1060, 1280);
      XTot = WF->IntegrateWaveForm(980, 1280);
      X1 = X1 / XTot;
      X2 = X2 / XTot;
      energy0 = vecOfPairs[i]->GetPairHitCh(0) == 2
                    ? vecOfPairs[i]->GetPairHitEvalEnergy(0) * 0.01911 - 0.301
                    : vecOfPairs[i]->GetPairHitEvalEnergy(1) * 0.01911 - 0.301;
      MT = vecOfPairs[i]->GetPairHitCh(0) == 2
               ? vecOfPairs[i]->GetHitPtr(0)->GetMeanTime()
               : vecOfPairs[i]->GetHitPtr(1)->GetMeanTime();
      hEMT->Fill(energy0, MT);
      //   if (MT > 1.9 and MT < 2.05)
      hX1X2->Fill(X1, X2);
      if (energy0 > 5 and energy0 < 15 and MT > 1.9 and MT < 2.05 and
          X1 > X1LowCut and X1 < X1HiCut and X2 < mX * X1 + mC_hi and
          X2 > mX * X1 + mC_low) {
        WFNorm = new digiAnalysis::WaveForm(WF->ScaleWaveForm(1.0 / energy0));
        waveformVector.push_back(*WFNorm);
      }
    }
  }

  UShort_t wfSz = WF->GetSize();
  std::cout << " Averaging " << waveformVector.size()
            << " Signal waveforms for plotting" << std::endl;
  digiAnalysis::WaveForm WFAveraged(wfSz, waveformVector);
  WFAveraged.SetSmooth(40);
  cluster1 = WFAveraged.GetTracesSmooth();

  // Define noise WF as cluster0 using X1X2 cuts
  std::string fname1 =
      "/home/kirtikesh/Analysis/DATA/extCoinc/"
      "NaI3124_17Jun26_NoSrc_1350V_2000V_1350V_1350V_Gain2_NoSplit_ExtTrig_"
      "Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_1000nsCoinc_"
      "2Vpp_Thresh_100lsb_WAVES_Sum_BLCorrected.root";
  digiAnalysis::Analysis an1(2, fname1, 0, 100000, 0);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVec1 =
      an1.GetSingleHitsVec();

  double X1LowCut_1 = 0.0, X1HiCut_1 = 0.8;
  double X2LowCut_1 = 0.08, X2HiCut_1 = 0.95;
  double mX_N = -1.0875, mC_N = 0.95;
  double X1_1, X2_1, XTot_1;
  for (int i = 0; i < std::min(ULong_t(20000), hitsVec1.size()); i++) {
    if (hitsVec1[i]->GetChNum() == 2) {
      WF = hitsVec1[i]->GetWFPtr();
      X1_1 = WF->IntegrateWaveForm(980, 1025);
      X2_1 = WF->IntegrateWaveForm(1060, 1280);
      XTot_1 = WF->IntegrateWaveForm(980, 1280);
      X1_1 = X1_1 / XTot_1;
      X2_1 = X2_1 / XTot_1;
      preint = WF->IntegrateWaveForm(0, digiAnalysis::GateStart) /
               digiAnalysis::GateStart;
      energy0 = hitsVec1[i]->GetEvalEnergy() * 0.01911 - 0.301;
      MT = hitsVec1[i]->GetMeanTime();
      if (preint < 2.0 and energy0 > 1. and energy0 < 10 and
          !(MT > 1.8 and MT < 2.1) and
          !(X1_1 > X1LowCut_1 and X1_1 < X1HiCut_1 and X2_1 > X2LowCut_1 and
            X2_1 < mX_N * X1_1 + mC_N)) {
        WFNorm = new digiAnalysis::WaveForm(WF->ScaleWaveForm(1.0 / energy0));
        wfVecNoise.push_back(*WFNorm);
        if (keepGoing) {
          WFNorm->Plot();
          std::cout << "Energy: " << energy0 << " X1: " << X1_1
                    << " X2: " << X2_1 << std::endl;
          std::cout << "Do you want to see the next waveform? (y/n): ";
          std::getline(std::cin, userInput);
          if (userInput != "y" && userInput != "Y") {
            keepGoing = false;
          }
        }
      }
    }
  }
  std::cout << " Averaging " << wfVecNoise.size()
            << " Noise waveforms for plotting" << std::endl;
  digiAnalysis::WaveForm WFAvNoise(wfSz, wfVecNoise);
  WFAvNoise.SetSmooth(40);
  cluster0 = WFAvNoise.GetTracesSmooth();
  WFAveraged.Plot(WFAveraged.NormWaveForm().second);
  WFAveraged.Plot(WFAvNoise.NormWaveForm().second, "SAME_kRed");

  // Now test log likelihood on the pair file
  double logLikelihood0, logLikelihood1, netLikelihood;
  double intValNoise, intValSig, intWF;
  intValSig = WFAveraged.IntegrateWaveForm();
  intValNoise = WFAvNoise.IntegrateWaveForm();

  std::vector<double> trace = WFNorm->GetTraces();
  TH2 *hELL = new TH2F("hELL", "hELL", 800, 0, 80, 1000, -10, 10);
  TH3 *hEMTLL = new TH3F("hEMTLL", "hEMTLL", 320, 0, 80, 50, 1, 3, 70, -4, 3);
  std::vector<digiAnalysis::WaveForm> *WFVec = nullptr;
  for (int i = 0; i < vecOfPairs.size(); i++) {
    if (i % 1000 == 0) {
      std::cout << i << std::endl;
    }
    WF = vecOfPairs[i]->GetPairHitCh(0) == 2
             ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()
             : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr();
    intWF = vecOfPairs[i]->GetPairHitCh(0) == 2
                ? vecOfPairs[i]->GetPairHitEvalEnergy(0)
                : vecOfPairs[i]->GetPairHitEvalEnergy(1);
    MT = vecOfPairs[i]->GetPairHitCh(0) == 2
             ? vecOfPairs[i]->GetHitPtr(0)->GetMeanTime()
             : vecOfPairs[i]->GetHitPtr(1)->GetMeanTime();
    energy0 = intWF * 0.01911 - 0.301;
    preint = WF->IntegrateWaveForm(0, digiAnalysis::GateStart) /
             digiAnalysis::GateStart;
    if (preint < 2.0) {
      WFNorm = new digiAnalysis::WaveForm(WF->ScaleWaveForm(1.0 / energy0));
      WFNorm->SetSmooth(100);
      trace = WFNorm->GetTracesSmooth();
      logLikelihood0 = 0;
      logLikelihood1 = 0;
      netLikelihood = 0;
      for (int j = digiAnalysis::GateStart + 30;
           j < digiAnalysis::GateStart + 250; j++) {
        logLikelihood0 += trace[j] * TMath::Log(abs(trace[j] / cluster0[j]));
        logLikelihood1 += trace[j] * TMath::Log(abs(trace[j] / cluster1[j]));
      }
      logLikelihood0 = logLikelihood0 + intValNoise - intWF;
      logLikelihood1 = logLikelihood1 + intValSig - intWF;
      netLikelihood =
          (logLikelihood0 - logLikelihood1) / (logLikelihood0 + logLikelihood1);
      hEMTLL->Fill(energy0, MT, netLikelihood);
      hELL->Fill(energy0, netLikelihood);
    }
  }

  TH2 *hELL_1 = new TH2F("hELL1", "hELL1", 800, 0, 80, 1000, -5, 5);
  TH2 *hEMT_1 = new TH2F("hEMT1", "hEMT1", 800, 0, 80, 1000, -5, 5);
  TH1F *hESpectra = new TH1F("hESpectra1", "Energy Spectra", 800, 0, 80);
  TH1F *hEEvalSpectra =
      new TH1F("hEEvalSpectra1", "Eval Energy Spectra", 800, 0, 80);
  TH1F *hELLCut = new TH1F("hELLCut1", "Energy Spectra", 800, 0, 80);
  TH2 *hMTLam_1 =
      new TH2F("MTLam_1", "MeanTime1 vs Lambda", 500, -4, 4, 2000, -10.0, 10.0);
  TH2 *hX1X2_Test =
      new TH2F("hX1X2_Test", "hX1X2_Test", 600, -0.2, 1, 500, -1, 1);
  TH2 *hMTLL1 =
      new TH2F("MTLL1", "MeanTime1 vs LL", 500, -4, 4, 1000, -10.0, 10.0);

  keepGoing = true;
  std::vector<digiAnalysis::WaveForm> waveformVector1;

  for (int i = 0; i < hitsVec1.size(); i++) {
    if (i % 1000 == 0) {
      std::cout << i << std::endl;
    }
    WF = hitsVec1[i]->GetWFPtr();
    // WF->SetSmooth(100);
    intWF = hitsVec1[i]->GetEvalEnergy();
    MT = hitsVec1[i]->GetMeanTime();
    energy0 = intWF * 0.01911 - 0.301;
    double preint = WF->IntegrateWaveForm(0, digiAnalysis::GateStart) /
                    digiAnalysis::GateStart;

    double lam = WF->EvalNoisePar2(1050, 1780);
    double psd = hitsVec1[i]->GetEvalPSD();

    X1_1 = WF->IntegrateWaveForm(980, 1025);
    X2_1 = WF->IntegrateWaveForm(1060, 1280);
    XTot_1 = WF->IntegrateWaveForm(980, 1280);
    X1_1 = X1_1 / XTot_1;
    X2_1 = X2_1 / XTot_1;
    hX1X2_Test->Fill(X1_1, X2_1);

    if (psd > 0.12 and psd < 1.0) {
      hESpectra->Fill(energy0);
    }

    digiAnalysis::WaveForm *WFNorm1 = new digiAnalysis::WaveForm(
        WF->ScaleWaveForm(1.0 / energy0)); // WF->NormWaveForm().second
    WFNorm1->SetSmooth(100);
    trace = WFNorm1->GetTracesSmooth();
    logLikelihood0 = 0;
    logLikelihood1 = 0;
    netLikelihood = 0;
    for (int j = digiAnalysis::GateStart + 30;
         j < digiAnalysis::GateStart + 250; j++) {
      logLikelihood0 += trace[j] * TMath::Log(abs(trace[j] / cluster0[j]));
      logLikelihood1 += trace[j] * TMath::Log(abs(trace[j] / cluster1[j]));
    }
    logLikelihood0 = logLikelihood0 + intValNoise - intWF;
    logLikelihood1 = logLikelihood1 + intValSig - intWF;
    netLikelihood =
        (logLikelihood0 - logLikelihood1) / (logLikelihood0 + logLikelihood1);
    if (psd > 0.12 and psd < 1.0 and preint < 2 and X1_1 > X1LowCut_1 and
        X1_1 < X1HiCut_1 and X2_1 < mX * X1_1 + mC_hi and
        X2_1 > mX * X1_1 + mC_low) {
      hMTLL1->Fill(MT, netLikelihood);
      hMTLam_1->Fill(MT, lam);
      hEMT_1->Fill(energy0, MT);
    }
    if (MT > 1.9 and MT < 2.1 and preint < 2.0 and psd > 0.12 and psd < 0.5 and
        lam > -5.4) {
      if (X1_1 > X1LowCut_1 and X1_1 < X1HiCut_1 and
          X2_1 < mX * X1_1 + mC_hi and X2_1 > mX * X1_1 + mC_low and
          lam > -5.5) {
        hEEvalSpectra->Fill(energy0);
        // } else
        // if (netLikelihood > 0.1 and keepGoing and energy0 < 5) {
        //   WF->Plot();
        //   std::cout << "Energy: " << energy0 << " LL: " << netLikelihood
        //             << " X1: " << X1_1 << " X2: " << X2_1 << std::endl;
        //   std::cout << "Do you want to see the next waveform? (y/n): ";
        //   std::getline(std::cin, userInput);
        //   if (userInput != "y" && userInput != "Y") {
        //     keepGoing = false;
        //   }
        // }
      }
      hELL_1->Fill(energy0, netLikelihood);

      // if (netLikelihood > -10.2 and netLikelihood < 0.0 and keepGoing and
      //     energy0 < 3) {
      //   WF->Plot();
      //   waveformVector1.push_back(*WF);
      //   std::cout << "Energy: " << energy0 << " LL: " << netLikelihood
      //             << std::endl;
      //   std::cout << "Do you want to see the next waveform? (y/n): ";
      //   std::getline(std::cin, userInput);
      //   if (userInput != "y" && userInput != "Y") {
      //     keepGoing = false;
      //   }
      // }
      if (netLikelihood > 0.0) {
        hELLCut->Fill(energy0);
      }
    }
  }

  // std::cout << " Averaging " << waveformVector1.size()
  //           << " waveforms for plotting" << std::endl;
  // digiAnalysis::WaveForm WFAveraged1(5000, waveformVector1);
  // WFAveraged1.SetSmooth(40);
  // WFAveraged1.SetTracesFFT();
  // WFAveraged1.Plot();
  // WFSig->Plot(cluster0, cluster1);

  std::cout << "Number of improper WF: " << badWFCount << std::endl;

  TCanvas *c3 = new TCanvas("c3", "MT", 800, 600);
  hEMT->Draw("COLZ");
  c3->Update();

  TCanvas *c5 = new TCanvas("c5", "E vs LL", 800, 600);
  hELL->Draw("COLZ");
  c5->Update();
  TCanvas *c7 = new TCanvas("c7", "E vs LL", 800, 600);
  hELL_1->Draw("COLZ");
  c7->Update();
  TCanvas *c8 = new TCanvas("c8", "MT vs LL", 800, 600);
  hMTLL1->Draw("COLZ");
  c8->Update();
  TCanvas *c9 = new TCanvas("c9", "E vs MT vs LL", 800, 600);
  hEMTLL->Draw();
  c9->Update();
  TCanvas *c10 = new TCanvas("c10", "Energy Spectra", 800, 600);
  hESpectra->SetLineColor(kRed);
  hEEvalSpectra->SetLineColor(kGreen);
  hELLCut->SetLineColor(kBlue);
  hESpectra->Draw("HIST");
  hEEvalSpectra->Draw("HISTSAME");
  hELLCut->Draw("HISTSAME");
  TCanvas *c11 = new TCanvas("c11", "MT vs Lam", 800, 600);
  hMTLam_1->Draw("COLZ");
  TCanvas *c12 = new TCanvas("c12", "E vs MT1", 800, 600);
  hEMT_1->Draw("COLZ");
  TCanvas *c14 = new TCanvas("c14", "X1 vs X2", 800, 600);
  hX1X2->Draw(); // "LEGO2"
  TCanvas *c15 = new TCanvas("c15", "X1 vs X2 test", 800, 600);
  hX1X2_Test->Draw();

  fApp->Run();
  return 0;
}