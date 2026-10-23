#include "Analysis.h"
#include "TMath.h"
#include "WaveForm.h"
#include "globals.h"
#include "includes.hh"
#include "singleHits.h"
#include <Rtypes.h>
#include <TApplication.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>
int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/extCoincSep/1800V/Data/PairFiles/"
      "pair_PreTrig1000.root";

  std::vector<double> cluster0;
  std::vector<double> cluster1;

  // Read to singleHits
  digiAnalysis::Analysis an(fname, 0, 100000, 0);
  // an.CreatePairs();
  std::vector<std::unique_ptr<digiAnalysis::Pair>> &vecOfPairs =
      an.GetPairsVec();
  std::cout << "Number of pairs found: " << vecOfPairs.size() << std::endl;

  double energy0 = 0., energy4 = 0., energy8;
  TH2 *hE1E2 = new TH2I("hE1E2", "hE1E2", 2600, 0, 5200, 2600, 0, 5200);
  TH1 *hDelT = new TH1F("hDelT", "Delta T (ns)", 1000, -1000, 1000);
  TH2 *hEMT = new TH2I("hEMT", "hEMT", 320, 0, 80, 500, -4, 4);
  TH2 *hEMTCut = new TH2I("hEMTCut", "hEMTCut", 320, 0, 80, 500, -4, 4);
  TH2 *hMTLam =
      new TH2F("MTLam", "MeanTime vs Lambda", 500, -4, 4, 1000, -10.0, 0);
  TH2 *hMTLL =
      new TH2F("MTLL", "MeanTime vs LL", 500, -4, 4, 1000, -10.0, 10.0);
  TH2 *hMTLL1 =
      new TH2F("MTLL1", "MeanTime1 vs LL", 500, -4, 4, 1000, -10.0, 10.0);
  TH2 *hLamLL =
      new TH2F("LamLL", "Lam vs LL", 1000, -10.0, 0, 1000, -10.0, 10.0);
  TH2 *hX1X2 = new TH2F("hX1X2", "hX1X2", 600, -0.2, 1, 500, -1, 1);
  TH1 *hECuts = new TH1F("hECuts", "ECuts", 320, 0, 80);
  TH1 *hELikeCuts = new TH1F("hELLCuts", "ELLCuts", 320, 0, 80);
  TH1 *hEOrig = new TH1F("hEOrig", "EOrig", 320, 0, 80);
  double X1LowCut = 0.2, X1HiCut = 0.4;
  double X2LowCut = 0.4, X2HiCut = 0.6;
  double X1 = 0, X2 = 0, XTot = 0;

  bool keepGoing = true;
  std::string userInput;
  int badWFCount = 0;
  double MT, Lam;
  digiAnalysis::WaveForm *WF = nullptr, *WFNorm = nullptr;
  std::vector<digiAnalysis::WaveForm> waveformVector, wfVecNoise;
  for (int i = 0; i < vecOfPairs.size(); i++) {
    X1 = 0;
    X2 = 0;
    XTot = 0;
    switch (vecOfPairs[i]->GetPairHitCh(1)) {
    // case 0: // for 2kV Data
    //   energyOther = vecOfPairs[i]->GetPairHitEnergy(1) * 0.47835 - 11.71;
    //   break;
    case 0: // 2 for 2kV Data
      // vecOfPairs[i]->GetHitPtr(1)->SetEvalEnergy();
      energy0 = vecOfPairs[i]->GetPairHitEvalEnergy(1) * 0.0062915 +
                0.3047; // * 0.01911 - 0.301 for 2kV Data;
      MT = vecOfPairs[i]->GetHitPtr(1)->GetMeanTime();
      vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()->SetSmooth(100);
      Lam = vecOfPairs[i]->GetHitPtr(1)->GetWFPtr()->EvalNoisePar2(1050, 1650);
      break;
    case 4: // 3 for 2kV Data
      energy4 = vecOfPairs[i]->GetPairHitEnergy(1) * 0.52958 -
                18.96; //* 0.54222 - 10.3 for 2kV Data;
      break;
    case 8: // 4 for 2kV Data
      energy8 = vecOfPairs[i]->GetPairHitEnergy(1) * 0.50226 -
                11.89; //* 0.48345 + 2.395 for 2kV Data;
      break;

    default:
      break;
    }
    switch (vecOfPairs[i]->GetPairHitCh(0)) {
      // case 0: // for 2kV Data
    //   energyOther = vecOfPairs[i]->GetPairHitEnergy(1) * 0.47835 - 11.71;
    //   break;
    case 0: // 2 for 2kV Data
      // vecOfPairs[i]->GetHitPtr(1)->SetEvalEnergy();
      energy0 = vecOfPairs[i]->GetPairHitEvalEnergy(0) * 0.0062915 +
                0.3047; // * 0.01911 - 0.301 for 2kV Data;
      // MT = vecOfPairs[i]->GetHitPtr(0)->GetMeanTime();
      // vecOfPairs[i]->GetHitPtr(1)->GetWFPtr()->SetSmooth(100);
      // Lam = vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()->EvalNoisePar2(1050,
      // 1650);
      break;
    case 4: // 3 for 2kV Data
      energy4 = vecOfPairs[i]->GetPairHitEnergy(0) * 0.52958 -
                18.96; //* 0.54222 - 10.3 for 2kV Data;
      break;
    case 8: // 4 for 2kV Data
      energy8 = vecOfPairs[i]->GetPairHitEnergy(0) * 0.50226 -
                11.89; //* 0.48345 + 2.395 for 2kV Data;
      break;

    default:
      break;
    }

    // if (vecOfPairs[i]->GetHitPtr(0)->GetMeanTime() > 2.08)
    if (vecOfPairs[i]->GetPairHitCh(0) != 0 and
        vecOfPairs[i]->GetPairHitCh(1) != 0) {
      // double preint =
      //     vecOfPairs[i]->GetPairHitCh(0) == 0
      //         ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()->IntegrateWaveForm(
      //               0, digiAnalysis::GateStart)
      //         : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr()->IntegrateWaveForm(
      //               0, digiAnalysis::GateStart);
      // if (preint / digiAnalysis::GateStart < 2.0) {
      hE1E2->Fill(energy4, energy8);
      // WF = vecOfPairs[i]->GetPairHitCh(0) == 0
      //          ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()
      //          : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr();
      // X1 = WF->IntegrateWaveForm(980, 1025);
      // X2 = WF->IntegrateWaveForm(1060, 1280);
      // XTot = WF->IntegrateWaveForm(980, 1280);
      // X1 = X1 / XTot;
      // X2 = X2 / XTot;
      hDelT->Fill(vecOfPairs[i]->GetPairDelTime() / 1E3);
      // hEMT->Fill(energy0, MT);
      // // if (energy0 > 15) {
      // hMTLam->Fill(MT, Lam);
      // }
      // hX1X2->Fill(X1, X2);
      // } else {
      //   badWFCount++;
      // }
      // std::cout << "E0: " << energy0 << " : EO: " << energyOther <<
      // std::endl;

      // if (keepGoing and energy0 > 5 and energy0 < 15 and MT > 1.9 and
      //     MT < 2.05 and preint < 2.1 and Lam > -5.4 and X1 > X1LowCut and
      //     X1 < X1HiCut and X2 < X2HiCut and X2 > X2LowCut) {

      //   // std::cout << "Energy: " << energy0 << " | MeanTime : " << MT
      //   //           << std::endl;
      //   WF = vecOfPairs[i]->GetPairHitCh(0) == 2
      //            ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()
      //            : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr();
      //   WFNorm = new digiAnalysis::WaveForm(WF->ScaleWaveForm(1.0 /
      //   energy0));
      //   // WF->NormWaveForm().second
      //   waveformVector.push_back(*WFNorm);
      //   // WF->Plot();
      //   // std::cout << "Do you want to see the next waveform? (y/n): ";
      //   // std::getline(std::cin, userInput);
      //   // if (userInput != "y" && userInput != "Y") {
      //   //   keepGoing = false;
      //   // }
      // }
      // if (energy0 > 0 and energy0 < 4 and !(MT > 1.75) and preint < 2.1 and
      //     Lam > -5.4 and !(X1 > X1LowCut and X1 < X1HiCut) and
      //     !(X2 < X2HiCut and X2 > X2LowCut)) { // or MT > 2.1

      //   // std::cout << "Energy: " << energy0 << " | MeanTime : " << MT
      //   //           << std::endl;
      //   WF = vecOfPairs[i]->GetPairHitCh(0) == 2
      //            ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()
      //            : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr();
      //   WFNorm = new digiAnalysis::WaveForm(
      //       WF->ScaleWaveForm(1.0 / energy0)); // WF->NormWaveForm().second

      //   wfVecNoise.push_back(*WFNorm);
      //   // if (keepGoing) {
      //   //   WF->Plot();
      //   //   std::cout << "Do you want to see the next waveform? (y/n): ";
      //   //   std::getline(std::cin, userInput);
      //   //   if (userInput != "y" && userInput != "Y") {
      //   //     keepGoing = false;
      //   //   }
      //   // }
      // }
    }
  }

  std::cout << "Number of improper WF: " << badWFCount << std::endl;
  TCanvas *c1 = new TCanvas("c1", "E1E2", 800, 600);
  c1->SetFrameFillColor(kBlack);
  hE1E2->Draw("COLZ");
  c1->Update();
  TCanvas *c2 = new TCanvas("c2", "DelT", 800, 600);
  hDelT->Draw("HIST");
  c2->Update();
  // TCanvas *c3 = new TCanvas("c3", "MT", 800, 600);
  // hEMT->Draw("COLZ");
  // c3->Update();
  // TCanvas *c4 = new TCanvas("c4", "MT vs Lam", 800, 600);
  // hMTLam->Draw("COLZ");
  // c4->Update();
  // TCanvas *c13 = new TCanvas("c13", "MT", 800, 600);
  // hEMTCut->Draw("COLZ");
  // c13->Update();

  // TCanvas *c14 = new TCanvas("c14", "X1 vs X2", 800, 600);
  // hX1X2->Draw(); // "LEGO2"

  fApp->Run();
  return 0;
}