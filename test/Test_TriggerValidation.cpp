#include "Analysis.h"
#include "WaveForm.h"
#include "globals.h"
#include "includes.hh"
#include "singleHits.h"
#include <RtypesCore.h>
#include <TApplication.h>
#include <TH2.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

std::vector<double> FindThresholdCrossings(const std::vector<double> &WF,
                                           double thresh, int start = 1) {
  std::vector<double> inVec;
  if (WF.empty())
    return inVec;
  // Handle case where waveform starts above threshold
  if (WF[0] >= thresh)
    inVec.push_back(0);
  for (int i = start; i < WF.size(); ++i) {
    if (WF[i - 1] < thresh && WF[i] >= thresh)
      inVec.push_back(static_cast<double>(i));
  }
  return inVec;
}

std::vector<double> MakeGateVector(const std::vector<double> &inVec,
                                   UShort_t gateLen, UShort_t delay = 0) {
  constexpr size_t outLen = 5000;
  std::vector<double> gateVec(outLen, 0);
  for (UShort_t trigger : inVec) {
    size_t start = static_cast<size_t>(trigger) + delay;
    if (start >= outLen)
      continue;
    size_t end = std::min(start + gateLen, outLen);
    for (size_t i = start; i < end; ++i)
      gateVec[i] = 1;
  }
  return gateVec;
}

int main() {
#ifdef WAVES
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::string fname;
  //   fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
  //           "CoincidenceStudies/01JuneNoSrc/"
  //           "NaI1342_June26_1750_1345_1350_1350_NoSrc_Thresh_15-30_"
  //           "300_WAVES_Coinc_144ns_LeadPit_Sum_BLCorrected.root";

  fname =
      "/home/kirtikesh/Analysis/DATA/extCoinc/"
      "NaI3124_17Jun26_NoSrc_1350V_2000V_1350V_1350V_Gain2_NoSplit_ExtTrig_"
      "Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_1000nsCoinc_"
      "2Vpp_Thresh_100lsb_WAVES_Sum_BLCorrected.root";

  //   fname =
  //       "/home/kirtikesh/Analysis/DATA/extCoinc/"
  //       "NaI1_11Jun26_NoSrc_1900V_ExtTrig_Thresh2_DelayCoincLogic_1600nsCoinc_"
  //       "2Vpp_Thresh_12lsb_WAVES_4/FILTERED/"
  //       "DataF_NaI1_11Jun26_NoSrc_1900V_ExtTrig_Thresh2_DelayCoincLogic_"
  //       "1600nsCoinc_2Vpp_Thresh_12lsb_WAVES_4_BLCorrected.root";

  // fname =
  //     "/home/kirtikesh/Analysis/DATA/extCoinc/"
  //     "NaI3124_15-17Jun26_NoSrc_1350V_2000V_1350V_1350V_Gain2_NoSplit_"
  //     "ExtTrig_"
  //     "Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_1000nsCoinc_"
  //     "2Vpp_Thresh_100lsb_WAVES_Sum_BLCorrected.root";

  //   digiAnalysis::GateStart = 400;

  digiAnalysis::Analysis an(2, fname, 0, 200000, 1);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int nentries = hitsVector.size();
  std::cout << "Got " << nentries << " from file" << std::endl;

  int spectralsize = 4096;
  TH2 *hMTPlot = new TH2F("MTPlot", "Energy vs MeanTime", spectralsize, 0,
                          spectralsize, 500, -4, 4);
  TH2 *hLamPlot = new TH2F("LamPlot", "Energy vs Lambda", spectralsize, 0,
                           spectralsize, 2000, -10.0, 10.0);
  TH2 *hMTLam =
      new TH2F("MTLam", "MeanTime vs Lambda", 500, -4, 4, 2000, -10.0, 10.0);
  TH2 *hPSDPlot = new TH2F("PSDPlot", "Energy vs PSD", spectralsize, 0,
                           spectralsize, 500, -8, 2);
  TH1 *hEPlot = new TH1F("EPlot", "Energy", 4096, 0, spectralsize);

  TH2 *hMTPlotSel =
      new TH2F("MTPlotSel", "Energy vs MeanTime of Selected Events",
               spectralsize, 0, spectralsize, 500, -4, 4);
  TH2 *hLamPlotSel =
      new TH2F("LamPlotSel", "Energy vs Lambda of Selected Events",
               spectralsize, 0, spectralsize, 2000, -10.0, 10.0);
  TH2 *hMTLamSel = new TH2F("MTLamSel", "MeanTime vs Lambda of Selected Events",
                            500, -4, 4, 2000, -10.0, 10.0);
  TH2 *hPSDPlotSel = new TH2F("PSDPlotSel", "Energy vs PSD of Selected Events",
                              spectralsize, 0, spectralsize, 500, -8, 2);
  TH1 *hEPlotSel =
      new TH1F("EPlotSel", "Energy  of Selected Events", 4096, 0, spectralsize);

  TH2 *hMTPlotDel =
      new TH2F("MTPlotDel", "Energy vs MeanTime of Deleted Events",
               spectralsize, 0, spectralsize, 500, -4, 4);
  TH2 *hLamPlotDel =
      new TH2F("LamPlotDel", "Energy vs Lambda of Deleted Events", spectralsize,
               0, spectralsize, 2000, -10.0, 10.0);
  TH2 *hMTLamDel = new TH2F("MTLamDel", "MeanTime vs Lambda of Deleted Events",
                            500, -4, 4, 2000, -10.0, 10.0);
  TH2 *hPSDPlotDel = new TH2F("PSDPlotDel", "Energy vs PSD of Deleted Events",
                              spectralsize, 0, spectralsize, 500, -8, 2);
  TH1 *hEPlotDel =
      new TH1F("EPlotDel", "Energy  of Deleted Events", 4096, 0, spectralsize);

  TH1 *hCoincDelay =
      new TH1F("hCoincDelay", "Delay of trigger and coinc", 5000, 0, 5000);

  std::vector<double> trigVec, gateOrig, gateDelay;
  std::vector<double> trace;
  int gateLen = 80;
  int delay = 120;
  int gateDelayLen = 300;
  int threshold = 17;
  bool accFlag = false;
  double energy, MT, PSD, Lam, coincDelay;

  for (int iter = 0; iter < nentries; iter++) {
    if (iter % 10000 == 0) {
      std::cout << "Processing: " << iter << std::endl;
    }

    energy = hitsVector[iter]->GetEnergy();
    PSD = hitsVector[iter]->GetEvalPSD();
    MT = hitsVector[iter]->GetMeanTime();
    Lam = hitsVector[iter]->GetWFPtr()->EvalNoisePar2(
        digiAnalysis::GateStart, digiAnalysis::GateStart + 700);

    hEPlot->Fill(energy);
    hPSDPlot->Fill(energy, PSD);
    hMTPlot->Fill(energy, MT);
    hLamPlot->Fill(energy, Lam);
    hMTLam->Fill(MT, Lam);

    accFlag = false;
    trace = hitsVector[iter]->GetWFPtr()->GetTraces();
    trigVec = FindThresholdCrossings(trace, threshold, digiAnalysis::GateStart);
    gateOrig = MakeGateVector(trigVec, gateLen);
    gateDelay = MakeGateVector(trigVec, gateDelayLen, delay);
    int iterTr = 0;
    for (iterTr = digiAnalysis::GateStart + 10;
         iterTr < gateDelay.size() - delay; iterTr++) {
      // iterDelay = iterTr + delay;
      if (gateOrig[iterTr] == 1 and gateDelay[iterTr] == 1) {
        accFlag = true;
        break;
      }
    }
    coincDelay = iterTr - digiAnalysis::GateStart;
    if (accFlag and coincDelay < 650) {
      if (Lam > -6.5) {
        hEPlotSel->Fill(energy);
      }
      hPSDPlotSel->Fill(energy, PSD);
      hMTPlotSel->Fill(energy, MT);
      hLamPlotSel->Fill(energy, Lam);
      hMTLamSel->Fill(MT, Lam);
      hCoincDelay->Fill(coincDelay);

      //   if (abs(coincDelay - 210) < 10) {
      //     hitsVector[iter]->Print();
      //     hitsVector[iter]->GetWFPtr()->Plot();
      //     std::this_thread::sleep_for(std::chrono::seconds(1));
      //     hitsVector[iter]->GetWFPtr()->Plot(gateOrig);
      //     hitsVector[iter]->GetWFPtr()->Plot(gateDelay, "SAME_kRed");
      //     std::this_thread::sleep_for(std::chrono::seconds(1));
      //   }
    } else {
      hEPlotDel->Fill(energy);
      hPSDPlotDel->Fill(energy, PSD);
      hMTPlotDel->Fill(energy, MT);
      hLamPlotDel->Fill(energy, Lam);
      hMTLamDel->Fill(MT, Lam);

      //   if (abs(coincDelay - 210) < 10) {
      //     hitsVector[iter]->Print();
      //     hitsVector[iter]->GetWFPtr()->Plot();
      //     std::this_thread::sleep_for(std::chrono::seconds(1));
      //     hitsVector[iter]->GetWFPtr()->Plot(gateOrig);
      //     hitsVector[iter]->GetWFPtr()->Plot(gateDelay, "SAME_kRed");
      //     std::this_thread::sleep_for(std::chrono::seconds(1));
      //   }
    }
  }

  TCanvas *c2 = new TCanvas("c2", "Energy vs MeanTime", 800, 600);
  hMTPlot->Draw("COLZ");
  TCanvas *c3 = new TCanvas("c3", "Energy vs MeanTime Selected", 800, 600);
  hMTPlotSel->Draw("COLZ");
  TCanvas *c31 = new TCanvas("c31", "Energy vs MeanTime Deleted", 800, 600);
  hMTPlotDel->Draw("COLZ");

  //   TCanvas *c4 = new TCanvas("c4", "Energy vs PSD", 800, 600);
  //   hPSDPlot->Draw("COLZ");
  //   TCanvas *c5 = new TCanvas("c5", "Energy vs PSD Selected", 800, 600);
  //   hPSDPlotSel->Draw("COLZ");
  //   TCanvas *c51 = new TCanvas("c51", "Energy vs PSD Deleted", 800, 600);
  //   hPSDPlotDel->Draw("COLZ");

  //   TCanvas *c6 = new TCanvas("c6", "Energy vs Lam", 800, 600);
  //   hLamPlot->Draw("COLZ");
  //   TCanvas *c7 = new TCanvas("c7", "Energy vs Lam Selected", 800, 600);
  //   hLamPlotSel->Draw("COLZ");
  //   TCanvas *c71 = new TCanvas("c71", "Energy vs Lam Deleted", 800, 600);
  //   hLamPlotDel->Draw("COLZ");

  TCanvas *c8 = new TCanvas("c8", "MT vs Lam", 800, 600);
  hMTLam->Draw("COLZ");
  TCanvas *c9 = new TCanvas("c9", "MT vs Lam Selected", 800, 600);
  hMTLamSel->Draw("COLZ");
  //   TCanvas *c4 = new TCanvas("c4", "Energy vs PSD", 800, 600);
  //   hPSDPlot->Draw("COLZ");
  //   TCanvas *c5 = new TCanvas("c5", "Energy vs PSD Selected", 800, 600);
  //   hPSDPlotSel->Draw("COLZ");
  //   TCanvas *c51 = new TCanvas("c51", "Energy vs PSD Deleted", 800, 600);
  //   hPSDPlotDel->Draw("COLZ");

  //   TCanvas *c6 = new TCanvas("c6", "Energy vs Lam", 800, 600);
  //   hLamPlot->Draw("COLZ");
  //   TCanvas *c7 = new TCanvas("c7", "Energy vs Lam Selected", 800, 600);
  //   hLamPlotSel->Draw("COLZ");
  //   TCanvas *c71 = new TCanvas("c71", "Energy vs Lam Deleted", 800, 600);
  //   hLamPlotDel->Draw("COLZ");
  TCanvas *c91 = new TCanvas("c91", "MT vs Lam Deleted", 800, 600);
  hMTLamDel->Draw("COLZ");

  TCanvas *c10 =
      new TCanvas("c10", "Delay of Trigger and CoincTrigger", 800, 600);
  hCoincDelay->Draw("HIST");

  TCanvas *c1 = new TCanvas("c1", "Energy Spectra", 800, 600);
  hEPlot->SetLineColor(kBlack);
  hEPlotSel->SetLineColor(kGreen);
  hEPlotDel->SetLineColor(kRed);
  hEPlot->Draw("HIST");
  hEPlotSel->Draw("HISTSAME");
  hEPlotDel->Draw("HISTSAME");

  fApp->Run();
#endif
  return 0;
}