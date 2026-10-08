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
int main(int argc, char *argv[])
{
  TApplication *fApp = new TApplication("TEST", NULL, NULL);

  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/extCoinc/PairFiles/"
  //     "Pair_NaI3124_15-17Jun26_NoSrc_1350V_2000V_1350V_1350V_Gain2_"
  //     "NoSplit_"
  //     "ExtTrig_"
  //     "Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_1000nsCoinc_"
  //     "2Vpp_Thresh_100lsb_WAVES_Sum_BLCorrected.root";

  std::string fname =
      "/media/kirtikesh/KKBlack/ToCopy/PairFiles/"
      "Pair_NaI134_NoSrc_01Oct_1800_1337_1350_WAVES_NoSplitSignal_Gain_12_Acquisition_4_ExtTrig_Threshold1LSB_160nsPromt_240nsDelay_800nsLong_1496nsCoinc_FreeWrites_Ch248Singles_68_BLCorrected.root";

  // std::string wfname =
  //     "/home/kirtikesh/Analysis/MLTests/AverageWaveforms.root"; // file
  //     for
  //                                                               //
  //                                                               different
  //                                                               // type
  //                                                               of
  //                                                               //
  //                                                               waveforms
  //                                                               //
  //                                                               identified
  //                                                               by
  //                                                               // PCA
  // TFile *wffile = TFile::Open(wfname.c_str());
  // TTree *wfTree = (TTree *)wffile->Get("AverageWaveforms");
  // int clusterID;
  // std::vector<double> *wfPCA = nullptr;
  // wfTree->SetBranchAddress("ClusterID", &clusterID);
  // wfTree->SetBranchAddress("Waveform", &wfPCA);

  std::vector<double> cluster0;
  std::vector<double> cluster1;
  // digiAnalysis::WaveForm *WFNoise, *WFSig;
  // Long64_t nentries = wfTree->GetEntries();

  // for (Long64_t i = 0; i < nentries; i++) {

  //   wfTree->GetEntry(i);

  //   if (clusterID == -1) {
  //     cluster0 = *wfPCA;
  //     WFNoise = new digiAnalysis::WaveForm(cluster0);
  //   }

  //   if (clusterID == 0) {
  //     cluster1 = *wfPCA;
  //     WFSig = new digiAnalysis::WaveForm(cluster1);
  //   }
  // }
  // std::string fname =
  // "/media/kirtikesh/UbuntuFiles/SDataF_NaI1342_02June26_1750_1345_1350_1350_NoSrc_Thresh_30_300_WAVES_Coinc_144ns_LeadPit_1.root";

  // Read to singleHits
  digiAnalysis::Analysis an(fname, 0, 100000, 0);
  // an.CreatePairs();
  std::vector<std::unique_ptr<digiAnalysis::Pair>> &vecOfPairs =
      an.GetPairsVec();
  std::cout << "Number of pairs found: " << vecOfPairs.size() << std::endl;

  double energy0 = 0., energyOther = 0.;
  TH2 *hE1E2 = new TH2I("hE1E2", "hE1E2", 320, 0, 80, 2600, 0, 2600);
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
  for (int i = 0; i < vecOfPairs.size(); i++)
  {
    X1 = 0;
    X2 = 0;
    XTot = 0;
    switch (vecOfPairs[i]->GetPairHitCh(1))
    {
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
      energyOther = vecOfPairs[i]->GetPairHitEnergy(1) * 0.52958 -
                    18.96; //* 0.54222 - 10.3 for 2kV Data;
      break;
    case 8: // 4 for 2kV Data
      energyOther = vecOfPairs[i]->GetPairHitEnergy(1) * 0.50226 -
                    11.89; //* 0.48345 + 2.395 for 2kV Data;
      break;

    default:
      break;
    }
    switch (vecOfPairs[i]->GetPairHitCh(0))
    {
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
      energyOther = vecOfPairs[i]->GetPairHitEnergy(0) * 0.52958 -
                    18.96; //* 0.54222 - 10.3 for 2kV Data;
      break;
    case 8: // 4 for 2kV Data
      energyOther = vecOfPairs[i]->GetPairHitEnergy(0) * 0.50226 -
                    11.89; //* 0.48345 + 2.395 for 2kV Data;
      break;

    default:
      break;
    }

    // if (vecOfPairs[i]->GetHitPtr(0)->GetMeanTime() > 2.08)
    if (vecOfPairs[i]->GetPairHitCh(0) == 0 or
        vecOfPairs[i]->GetPairHitCh(1) == 0)
    {
      double preint =
          vecOfPairs[i]->GetPairHitCh(0) == 0
              ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()->IntegrateWaveForm(
                    0, digiAnalysis::GateStart)
              : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr()->IntegrateWaveForm(
                    0, digiAnalysis::GateStart);
      if (preint / digiAnalysis::GateStart < 2.0)
      {
        hE1E2->Fill(energy0, energyOther);
        WF = vecOfPairs[i]->GetPairHitCh(0) == 0
                 ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()
                 : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr();
        X1 = WF->IntegrateWaveForm(980, 1025);
        X2 = WF->IntegrateWaveForm(1060, 1280);
        XTot = WF->IntegrateWaveForm(980, 1280);
        X1 = X1 / XTot;
        X2 = X2 / XTot;
        hDelT->Fill(vecOfPairs[i]->GetPairDelTime() / 1E3);
        hEMT->Fill(energy0, MT);
        // // if (energy0 > 15) {
        // hMTLam->Fill(MT, Lam);
        // }
        hX1X2->Fill(X1, X2);
      }
      else
      {
        badWFCount++;
      }
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

  // UShort_t wfSz = WF->GetSize();
  // std::cout << " Averaging " << waveformVector.size()
  //           << " Signal like waveforms" << std::endl;
  // digiAnalysis::WaveForm WFAveraged(wfSz, waveformVector);
  // WFAveraged.SetSmooth(40);
  // WFAveraged.SetTracesFFT();
  // WFAveraged.Plot(WFAveraged.NormWaveForm().second);
  // double intValAv = WFAveraged.IntegrateWaveForm();
  // cluster1 = WFAveraged.ScaleWaveForm(
  //     1.0 / (intValAv * 0.01911 - 0.301)); // ;
  //     WFAveraged.NormWaveForm().second
  // std::cout << "intValAv: " << intValAv * 0.01911 - 0.301 << std::endl;

  // // std::cout << " Averaging " << wfVecNoise.size() << " Noise like
  // waveforms"
  // //           << std::endl;
  // // digiAnalysis::WaveForm WFAvNoise(wfSz, wfVecNoise);
  // // WFAvNoise.SetSmooth(40);
  // // double intValAv2 = WFAvNoise.IntegrateWaveForm();
  // // cluster0 = WFAvNoise.ScaleWaveForm(
  // //     1.0 / (intValAv2 * 0.01911 - 0.301)); //
  // //     WFAvNoise.NormWaveForm().second;
  // //                                           // //
  // //                                           WFAvNoise.GetTracesSmooth();
  // // WFAveraged.Plot(WFAvNoise.NormWaveForm().second, "SAME_kRed");

  // TFile f("output.root", "RECREATE");
  // std::vector<double> waveformtrace;
  // waveformtrace = WFAveraged.GetTraces();
  // f.WriteObject(&waveformtrace, "trace");
  // f.Close();

  // // Testing the likelihood w.r.t averaged waveform
  // double logLikelihood0, logLikelihood1, netLikelihood;
  // double intValNoise, intValSig, intWF;
  // // intValNoise = WFNoise->IntegrateWaveForm();
  // // intValSig = WFSig->IntegrateWaveForm();
  // intValSig = intValAv; //
  // digiAnalysis::WaveForm(cluster1).IntegrateWaveForm();
  // // intValNoise =      intValAv2; //
  // // digiAnalysis::WaveForm(cluster0).IntegrateWaveForm();

  // keepGoing = true;
  // std::vector<double> trace = WFNorm->GetTraces();
  // TH2 *hELL = new TH2F("hELL", "hELL", 320, 0, 80, 1000, -10, 10);
  // TH3 *hEMTLL = new TH3F("hEMTLL", "hEMTLL", 320, 0, 80, 50, 1, 3, 70, -4,
  // 3);

  // // Log likelihood test on the pair file shifted to end

  // // Testing the log likelihood on non-pair files
  // keepGoing = true;
  // TH2 *hELL_1 = new TH2F("hELL1", "hELL1", 320, 0, 80, 1000, -5, 5);
  // TH2 *hEMT_1 = new TH2F("hEMT1", "hEMT1", 320, 0, 80, 1000, -5, 5);
  // TH1F *hESpectra = new TH1F("hESpectra1", "Energy Spectra", 320, 0, 80);
  // TH1F *hEEvalSpectra =
  //     new TH1F("hEEvalSpectra1", "Eval Energy Spectra", 320, 0, 80);
  // TH1F *hELLCut = new TH1F("hELLCut1", "Energy Spectra", 320, 0, 80);
  // TH2 *hMTLam_1 =
  //     new TH2F("MTLam_1", "MeanTime1 vs Lambda", 500, -4, 4, 2000,
  //     -10.0, 10.0);
  // TH2 *hX1X2_Test =
  //     new TH2F("hX1X2_Test", "hX1X2_Test", 600, -0.2, 1, 500, -1, 1);
  // std::vector<digiAnalysis::WaveForm> waveformVector1;

  // std::string fname1 =
  //     "/home/kirtikesh/Analysis/DATA/extCoinc/"
  //     "NaI3124_17Jun26_NoSrc_1350V_2000V_1350V_1350V_Gain2_NoSplit_ExtTrig_"
  //     "Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_1000nsCoinc_"
  //     "2Vpp_Thresh_100lsb_WAVES_Sum_BLCorrected.root";
  // digiAnalysis::Analysis an1(2, fname1, 0, 200000, 0);
  // std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVec1 =
  //     an1.GetSingleHitsVec();

  // wfVecNoise.clear();
  // std::cout << "The noise waveform vector is of size: " << wfVecNoise.size()
  //           << std::endl;
  // double X1LowCut_1 = 0.195, X1HiCut_1 = 0.35;
  // double mX = -1.096, mC_low = 0.764, mC_hi = 0.834;
  // double X2LowCut_1 = 0.4, X2HiCut_1 = 0.6;
  // double X1_1 = 0, X2_1 = 0, XTot_1 = 0;
  // bool likelihoodFlag = false, cutFlag = false;

  // std::vector<double> vecLam, vecpreInt, vecpsd, vecX1, vecX2;
  // for (int i = 0; i < hitsVec1.size(); i++) {
  //   likelihoodFlag = false, cutFlag = false;
  //   if (i % 1000 == 0) {
  //     std::cout << i << std::endl;
  //   }
  //   WF = hitsVec1[i]->GetWFPtr();
  //   // WF->SetSmooth(100);
  //   intWF = hitsVec1[i]->GetEvalEnergy();
  //   MT = hitsVec1[i]->GetMeanTime();
  //   energy0 = intWF * 0.01911 - 0.301;
  //   double preint = WF->IntegrateWaveForm(0, digiAnalysis::GateStart) /
  //                   digiAnalysis::GateStart;

  //   double lam = WF->EvalNoisePar2(1050, 1780);
  //   double psd = hitsVec1[i]->GetEvalPSD();

  //   vecLam.push_back(lam);
  //   vecpreInt.push_back(preint);
  //   vecpsd.push_back(psd);

  //   X1_1 = WF->IntegrateWaveForm(980, 1025);
  //   X2_1 = WF->IntegrateWaveForm(1060, 1280);
  //   XTot_1 = WF->IntegrateWaveForm(980, 1280);
  //   X1_1 = X1_1 / XTot_1;
  //   X2_1 = X2_1 / XTot_1;
  //   hX1X2_Test->Fill(X1_1, X2_1);

  //   vecX1.push_back(X1_1);
  //   vecX2.push_back(X2_1);

  //   if (psd > 0.12 and psd < 1.0 and preint < 2) {
  //     hESpectra->Fill(energy0);
  //     hEMT_1->Fill(energy0, MT);
  //   }

  //   if (psd > 0.12 and psd < 1.0 and preint < 2 and X1_1 > X1LowCut_1 and
  //       X1_1 < X1HiCut_1 and X2_1 < mX * X1_1 + mC_hi and
  //       X2_1 > mX * X1_1 + mC_low) {
  //     hMTLL1->Fill(MT, netLikelihood);
  //     hMTLam_1->Fill(MT, lam);
  //   }
  //   if (MT > 1.9 and MT < 2.1 and preint < 2.0 and psd > 0.12 and psd < 0.5
  //   and
  //       lam > -5.4) {
  //     if (X1_1 > X1LowCut_1 and X1_1 < X1HiCut_1 and
  //         X2_1 < mX * X1_1 + mC_hi and X2_1 > mX * X1_1 + mC_low and
  //         lam > -5.5) {
  //       hEEvalSpectra->Fill(energy0);
  //     }
  //   }

  //   if (energy0 > 0 and energy0 < 30 and !(MT > 1.5) and preint < 2.1 and
  //       lam > -5.4 and (X2_1 < 0.2 or X2_1 > (-1.0 * X1_1 + 0.95))) {
  //     // and !(X1_1 > X1LowCut_1 and X1_1 < X1HiCut_1)
  //     // and !(X2_1 < mX * X1_1 + mC_hi and X2_1 > mX * X1_1 + mC_low)
  //     WFNorm = new digiAnalysis::WaveForm(
  //         WF->ScaleWaveForm(1.0 / energy0)); // WF->NormWaveForm().second
  //     wfVecNoise.push_back(*WFNorm);
  //   }
  // }
  // std::cout << " Averaging " << wfVecNoise.size() << " Noise like waveforms"
  //           << std::endl;
  // digiAnalysis::WaveForm WFAvNoise1(wfSz, wfVecNoise);
  // WFAvNoise1.SetSmooth(40);
  // double intValAv2 = WFAvNoise1.IntegrateWaveForm();
  // cluster0 = WFAvNoise1.ScaleWaveForm(
  //     1.0 / (intValAv2 * 0.01911 - 0.301)); //
  //     WFAvNoise.NormWaveForm().second;
  //                                           // //
  //                                           WFAvNoise.GetTracesSmooth();
  // intValNoise = intValAv2;
  // WFAvNoise1.Plot(cluster0);
  // WFAvNoise1.Plot(cluster1, "SAME_kBlack");

  // TH2 *hBadLL = new TH2F("hBadLL", "BadLL", 600, -0.2, 1, 250, 0, 1);
  // TH2 *hGoodLL = new TH2F("hGoodLL", "GoodLL", 600, -0.2, 1, 250, 0, 1);

  // // for (int i = 0; i < hitsVec1.size(); i++) {
  // //   likelihoodFlag = false, cutFlag = false;
  // //   if (i % 1000 == 0) {
  // //     std::cout << i << std::endl;
  // //   }
  // //   WF = hitsVec1[i]->GetWFPtr();
  // //   intWF = hitsVec1[i]->GetEvalEnergy();
  // //   MT = hitsVec1[i]->GetMeanTime();
  // //   energy0 = intWF * 0.01911 - 0.301;
  // //   double preint = vecpreInt[i];

  // //   double lam = vecLam[i];
  // //   double psd = vecpsd[i];

  // //   X1_1 = vecX1[i];
  // //   X2_1 = vecX2[i];

  // //   digiAnalysis::WaveForm *WFNorm1 = new digiAnalysis::WaveForm(
  // //       WF->ScaleWaveForm(1.0 / energy0)); // WF->NormWaveForm().second
  // //   WFNorm1->SetSmooth(100);
  // //   trace = WFNorm1->GetTracesSmooth();
  // //   logLikelihood0 = 0;
  // //   logLikelihood1 = 0;
  // //   netLikelihood = 0;
  // //   for (int j = digiAnalysis::GateStart + 25;
  // //        j < digiAnalysis::GateStart + 470; j++) {
  // //     logLikelihood0 += trace[j] * TMath::Log(abs(trace[j] /
  // cluster0[j]));
  // //     logLikelihood1 += trace[j] * TMath::Log(abs(trace[j] /
  // cluster1[j]));
  // //   }
  // //   logLikelihood0 = logLikelihood0 + intValNoise - intWF;
  // //   logLikelihood1 = logLikelihood1 + intValSig - intWF;
  // //   netLikelihood =
  // //       (logLikelihood0 - logLikelihood1) / (logLikelihood0 +
  // //       logLikelihood1);
  // //   if (psd > 0.12 and psd < 1.0 and preint < 2 and X1_1 > X1LowCut_1 and
  // //       X1_1 < X1HiCut_1 and X2_1 < mX * X1_1 + mC_hi and
  // //       X2_1 > mX * X1_1 + mC_low) {
  // //     hMTLL1->Fill(MT, netLikelihood);
  // //   }
  // //   if (MT > 1.9 and MT < 2.1 and preint < 2.0 and psd > 0.12 and psd <
  // 0.5
  // //   and
  // //       lam > -5.4) {
  // //     if (X1_1 > X1LowCut_1 and X1_1 < X1HiCut_1 and
  // //         X2_1 < mX * X1_1 + mC_hi and X2_1 > mX * X1_1 + mC_low and
  // //         lam > -5.5) {
  // //       cutFlag = true;
  // //     }
  // //     hELL_1->Fill(energy0, netLikelihood);

  // //     if (netLikelihood > 0.14 and
  // //         X2_1 < (-1.031 * X1_1 +
  // //                 0.99)) { // and X1_1 > X1LowCut_1 and X1_1 < X1HiCut_1
  // and
  // //                          // X2_1<mX * X1_1 + mC_hi and X2_1> mX *X1_1 +
  // //                          mC_low
  // //       hELLCut->Fill(energy0);
  // //       likelihoodFlag = true;
  // //       hGoodLL->Fill(X1_1, X2_1);
  // //     }
  // //     if (!likelihoodFlag) {

  // //       // !(likelihoodFlag and cutFlag) and (likelihoodFlag or cutFlag)
  // and
  // //       // energy0 < 10 and keepGoing

  // //       // WF->Plot();
  // //       // std::cout << "Energy: " << energy0 << " LL: " <<
  // //       // netLikelihood
  // //       //           << " X1: " << X1_1 << " X2: " << X2_1 << std::endl;
  // //       // std::cout << "Do you want to see the next waveform?
  // //       // (yy/yn/n): "; std::getline(std::cin, userInput); if
  // //       // (userInput == "yy") {
  // //       hBadLL->Fill(X1_1, X2_1);
  // //       // }
  // //       // if (userInput != "yy" && userInput != "yn" && userInput != "y")
  // {
  // //       //   keepGoing = false;
  // //       // }
  // //     }
  // //   }
  // // }

  // for (int i = 0; i < vecOfPairs.size(); i++) {
  //   if (i % 1000 == 0) {
  //     std::cout << i << std::endl;
  //   }
  //   WF = vecOfPairs[i]->GetPairHitCh(0) == 2
  //            ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()
  //            : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr();
  //   intWF = vecOfPairs[i]->GetPairHitCh(0) == 2
  //               ? vecOfPairs[i]->GetPairHitEvalEnergy(0)
  //               : vecOfPairs[i]->GetPairHitEvalEnergy(1);
  //   MT = vecOfPairs[i]->GetPairHitCh(0) == 2
  //            ? vecOfPairs[i]->GetHitPtr(0)->GetMeanTime()
  //            : vecOfPairs[i]->GetHitPtr(1)->GetMeanTime();

  //   X1 = WF->IntegrateWaveForm(980, 1025);
  //   X2 = WF->IntegrateWaveForm(1060, 1280);
  //   XTot = WF->IntegrateWaveForm(980, 1280);
  //   X1 = X1 / XTot;
  //   X2 = X2 / XTot;

  //   energy0 = intWF * 0.01911 - 0.301;
  //   WFNorm = new digiAnalysis::WaveForm(
  //       WF->ScaleWaveForm(1.0 / energy0)); // WF->NormWaveForm().second

  //   WFNorm->SetSmooth(100);
  //   trace = WFNorm->GetTracesSmooth();
  //   logLikelihood0 = 0;
  //   logLikelihood1 = 0;
  //   netLikelihood = 0;
  //   for (int j = digiAnalysis::GateStart + 25;
  //        j < digiAnalysis::GateStart + 470; j++) {
  //     // logLikelihood0 += trace[j] * TMath::Log(abs(trace[j] /
  //     // waveformtrace[j]));
  //     logLikelihood0 += trace[j] * TMath::Log(abs(trace[j] / cluster0[j]));
  //     logLikelihood1 += trace[j] * TMath::Log(abs(trace[j] / cluster1[j]));
  //   }
  //   logLikelihood0 = logLikelihood0 + intValNoise - intWF;
  //   logLikelihood1 = logLikelihood1 + intValSig - intWF;
  //   netLikelihood =
  //       (logLikelihood0 - logLikelihood1) / (logLikelihood0 +
  //       logLikelihood1);

  //   double preint =
  //       vecOfPairs[i]->GetPairHitCh(0) == 2
  //           ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()->IntegrateWaveForm(
  //                 0, digiAnalysis::GateStart)
  //           : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr()->IntegrateWaveForm(
  //                 0, digiAnalysis::GateStart);
  //   preint /= digiAnalysis::GateStart;
  //   int chSel = vecOfPairs[i]->GetPairHitCh(1) == 2 ? 1 : 0;
  //   Lam =
  //       vecOfPairs[i]->GetHitPtr(chSel)->GetWFPtr()->EvalNoisePar2(1050,
  //       1650);
  //   if (preint < 2.0) {
  //     hEOrig->Fill(energy0);
  //     if (MT > 1.8 and MT < 2.1 and X1 > X1LowCut and X1 < X1HiCut and
  //         X2 < X2HiCut and X2 > X2LowCut)
  //       hELL->Fill(energy0, netLikelihood);
  //     hMTLL->Fill(MT, netLikelihood);
  //     hEMTLL->Fill(energy0, MT, netLikelihood);
  //   }
  //   if (preint < 2.0 and MT > 1.8 and MT < 2.1 and X1 > X1LowCut and
  //       X1 < X1HiCut and X2 < X2HiCut and
  //       X2 > X2LowCut) { // and netLikelihood < 0.0
  //     if (vecOfPairs[i]->GetPairDelTime() / 1E3 > -20)
  //       hE1E2->Fill(energy0, energyOther);
  //     hDelT->Fill(vecOfPairs[i]->GetPairDelTime() / 1E3);
  //     // if (energy0 > 15) {
  //     hMTLam->Fill(MT, Lam);
  //     hECuts->Fill(energy0);
  //     // }
  //   }
  //   if (preint < 2.0 and netLikelihood > 0.1) {
  //     hEMTCut->Fill(energy0, MT);
  //     hELikeCuts->Fill(energy0);
  //   }
  // }

  // std::cout << " Averaging " << waveformVector1.size()
  //           << " waveforms for plotting" << std::endl;
  // digiAnalysis::WaveForm WFAveraged1(5000, waveformVector1);
  // WFAveraged1.SetSmooth(40);
  // WFAveraged1.SetTracesFFT();
  // WFAveraged1.Plot();
  // WFSig->Plot(cluster0, cluster1);

  std::cout << "Number of improper WF: " << badWFCount << std::endl;
  TCanvas *c1 = new TCanvas("c1", "E1E2", 800, 600);
  c1->SetFrameFillColor(kBlack);
  hE1E2->Draw("COLZ");
  c1->Update();
  TCanvas *c2 = new TCanvas("c2", "DelT", 800, 600);
  hDelT->Draw("HIST");
  c2->Update();
  TCanvas *c3 = new TCanvas("c3", "MT", 800, 600);
  hEMT->Draw("COLZ");
  c3->Update();
  TCanvas *c4 = new TCanvas("c4", "MT vs Lam", 800, 600);
  hMTLam->Draw("COLZ");
  c4->Update();
  TCanvas *c13 = new TCanvas("c13", "MT", 800, 600);
  hEMTCut->Draw("COLZ");
  c13->Update();

  // TCanvas *c5 = new TCanvas("c5", "E vs LL", 800, 600);
  // hELL->Draw("COLZ");
  // c5->Update();
  // TCanvas *c6 = new TCanvas("c6", "MT vs LL", 800, 600);
  // hMTLL->Draw("COLZ");
  // c6->Update();
  // TCanvas *c7 = new TCanvas("c7", "E vs LL", 800, 600);
  // hELL_1->Draw("COLZ");
  // c7->Update();
  // TCanvas *c8 = new TCanvas("c8", "MT vs LL", 800, 600);
  // hMTLL1->Draw("COLZ");
  // c8->Update();
  // TCanvas *c9 = new TCanvas("c9", "E vs MT vs LL", 800, 600);
  // hEMTLL->Draw();
  // c9->Update();
  // TCanvas *c10 = new TCanvas("c10", "Energy Spectra", 800, 600);
  // hESpectra->SetLineColor(kRed);
  // hEEvalSpectra->SetLineColor(kBlue);
  // hELLCut->SetLineColor(kGreen);
  // hESpectra->SetLineWidth(4);
  // hEEvalSpectra->SetLineWidth(4);
  // hELLCut->SetLineWidth(2);
  // hESpectra->GetXaxis()->SetTitle("Energy (keV)");
  // hESpectra->GetYaxis()->SetTitle("Counts");
  // hESpectra->Draw("HIST");
  // hELLCut->Draw("HISTSAME");
  // hEEvalSpectra->Draw("HISTSAME");
  // TLegend *legend = new TLegend(0.65, 0.70, 0.88, 0.88);
  // legend->AddEntry(hESpectra, "Original spectra", "l");
  // legend->AddEntry(hEEvalSpectra, "Parameter cuts", "l");
  // legend->AddEntry(hELLCut, "Likelihood cut", "l");
  // legend->Draw();

  // TCanvas *c11 = new TCanvas("c11", "MT vs Lam", 800, 600);
  // hMTLam_1->Draw("COLZ");
  // TCanvas *c12 = new TCanvas("c12", "E vs MT1", 800, 600);
  // hEMT_1->Draw("COLZ");
  TCanvas *c14 = new TCanvas("c14", "X1 vs X2", 800, 600);
  hX1X2->Draw(); // "LEGO2"
  // TCanvas *c15 = new TCanvas("c15", "X1 vs X2 test", 800, 600);
  // hX1X2_Test->Draw();
  // TCanvas *c16 = new TCanvas("c16", "Energy Spectra Pair", 800, 600);
  // hEOrig->SetLineColor(kRed);
  // hECuts->SetLineColor(kBlue);
  // hELikeCuts->SetLineColor(kGreen);
  // hEOrig->GetXaxis()->SetTitle("Energy (keV)");
  // hEOrig->GetYaxis()->SetTitle("Counts");
  // hEOrig->Draw("HIST");
  // hELikeCuts->Draw("HISTSAME");
  // hECuts->Draw("HISTSAME");
  // TLegend *leg = new TLegend(0.65, 0.70, 0.88, 0.88);
  // leg->AddEntry(hEOrig, "Acquired Data", "l");
  // leg->AddEntry(hECuts, "Parameter Cuts", "l");
  // leg->AddEntry(hELikeCuts, "Likelihood Cut", "l");
  // leg->Draw();
  // TCanvas *c17 = new TCanvas("c17", "Bad X1-X2 passing LL", 800, 600);
  // hBadLL->Draw();
  // TCanvas *c18 = new TCanvas("c18", "Good X1-X2 passing LL", 800, 600);
  // hGoodLL->Draw();
  fApp->Run();
  return 0;
}