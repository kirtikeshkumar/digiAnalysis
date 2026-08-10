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
#include <iostream>
#include <ostream>
#include <string>
#include <vector>
int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);

  // std::string fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
  //                     "CoincidenceStudies/01JuneNoSrc/"
  //                     "NaI1342_June26_1750_1345_1350_1350_NoSrc_Thresh_15-30_"
  //                     "300_WAVES_Coinc_144ns_LeadPit_Sum_BLCorrected.root";

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/extCoinc/PairFiles/"
      "Pair_NaI3124_15-17Jun26_NoSrc_1350V_2000V_1350V_1350V_Gain2_"
      "NoSplit_"
      "ExtTrig_"
      "Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_1000nsCoinc_"
      "2Vpp_Thresh_100lsb_WAVES_Sum_BLCorrected.root";

  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/extCoinc/"
  //     "NaI1234_19Jun26_NoSrc_2000V_1350V_1350V_1350V_Gain2_NoSplit_Ch5_ExtTrig_"
  //     "Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_"
  //     "2000nsCoincCh0AndAny_2Vpp_Thresh_30lsb_WAVES_Sum_BLCorrected.root";

  std::string wfname =
      "/home/kirtikesh/Analysis/MLTests/AverageWaveforms.root"; // file for
                                                                // different
                                                                // type of
                                                                // waveforms
                                                                // identified by
                                                                // PCA
  TFile *wffile = TFile::Open(wfname.c_str());
  TTree *wfTree = (TTree *)wffile->Get("AverageWaveforms");
  int clusterID;
  std::vector<double> *wfPCA = nullptr;
  wfTree->SetBranchAddress("ClusterID", &clusterID);
  wfTree->SetBranchAddress("Waveform", &wfPCA);

  std::vector<double> cluster0;
  std::vector<double> cluster1;
  digiAnalysis::WaveForm *WFNoise, *WFSig;
  Long64_t nentries = wfTree->GetEntries();

  for (Long64_t i = 0; i < nentries; i++) {

    wfTree->GetEntry(i);

    if (clusterID == -1) {
      cluster0 = *wfPCA;
      WFNoise = new digiAnalysis::WaveForm(cluster0);
    }

    if (clusterID == 0) {
      cluster1 = *wfPCA;
      WFSig = new digiAnalysis::WaveForm(cluster1);
    }
  }
  // std::string fname =
  // "/media/kirtikesh/UbuntuFiles/SDataF_NaI1342_02June26_1750_1345_1350_1350_NoSrc_Thresh_30_300_WAVES_Coinc_144ns_LeadPit_1.root";

  // Read to singleHits
  digiAnalysis::Analysis an(fname, 0, 100000, 0);
  // an.CreatePairs();
  std::vector<std::unique_ptr<digiAnalysis::Pair>> &vecOfPairs =
      an.GetPairsVec();
  std::cout << "Number of pairs found: " << vecOfPairs.size() << std::endl;

  double energy0 = 0., energyOther = 0.;
  TH2 *hE1E2 = new TH2I("hE1E2", "hE1E2", 800, 0, 80, 2600, 0, 2600);
  TH1 *hDelT = new TH1F("hDelT", "Delta T (ns)", 1000, -1000, 1000);
  TH2 *hEMT = new TH2I("hEMT", "hEMT", 800, 0, 80, 500, -4, 4);
  TH2 *hEMTCut = new TH2I("hEMTCut", "hEMTCut", 800, 0, 80, 500, -4, 4);
  TH2 *hMTLam =
      new TH2F("MTLam", "MeanTime vs Lambda", 500, -4, 4, 1000, -10.0, 0);
  TH2 *hMTLL =
      new TH2F("MTLL", "MeanTime vs LL", 500, -4, 4, 1000, -10.0, 10.0);
  TH2 *hMTLL1 =
      new TH2F("MTLL1", "MeanTime1 vs LL", 500, -4, 4, 1000, -10.0, 10.0);
  TH2 *hLamLL =
      new TH2F("LamLL", "Lam vs LL", 1000, -10.0, 0, 1000, -10.0, 10.0);
  TH2 *hX1X2 = new TH2F("hX1X2", "hX1X2", 60, -0.2, 1, 50, -1, 1);
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
    case 0:
      energyOther = vecOfPairs[i]->GetPairHitEnergy(1) * 0.47835 - 11.71;
      break;
    case 2:
      // vecOfPairs[i]->GetHitPtr(1)->SetEvalEnergy();
      energy0 = vecOfPairs[i]->GetPairHitEvalEnergy(1) * 0.01911 - 0.301;
      MT = vecOfPairs[i]->GetHitPtr(1)->GetMeanTime();
      Lam = vecOfPairs[i]->GetHitPtr(1)->GetWFPtr()->EvalNoisePar2(950, 1650);
      break;
    case 3:
      energyOther = vecOfPairs[i]->GetPairHitEnergy(1) * 0.54222 - 10.3;
      break;
    case 4:
      energyOther = vecOfPairs[i]->GetPairHitEnergy(1) * 0.48345 + 2.395;
      break;

    default:
      break;
    }
    switch (vecOfPairs[i]->GetPairHitCh(0)) {
    case 0:
      energyOther = vecOfPairs[i]->GetPairHitEnergy(0) * 0.47835 - 11.71;
      break;
    case 2:
      // vecOfPairs[i]->GetHitPtr(0)->SetEvalEnergy();
      energy0 = vecOfPairs[i]->GetPairHitEvalEnergy(0) * 0.01911 - 0.301;
      MT = vecOfPairs[i]->GetHitPtr(0)->GetMeanTime();
      vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()->SetSmooth(100);
      Lam = vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()->EvalNoisePar2(950, 1650);
      break;
    case 3:
      energyOther = vecOfPairs[i]->GetPairHitEnergy(0) * 0.54222 - 10.3;
      break;
    case 4:
      energyOther = vecOfPairs[i]->GetPairHitEnergy(0) * 0.48345 + 2.395;
      break;

    default:
      break;
    }

    // if (vecOfPairs[i]->GetHitPtr(0)->GetMeanTime() > 2.08)
    if (vecOfPairs[i]->GetPairHitCh(0) == 2 or
        vecOfPairs[i]->GetPairHitCh(1) == 2) {
      double preint =
          vecOfPairs[i]->GetPairHitCh(0) == 2
              ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()->IntegrateWaveForm(
                    0, digiAnalysis::GateStart)
              : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr()->IntegrateWaveForm(
                    0, digiAnalysis::GateStart);
      if (preint / digiAnalysis::GateStart < 2.0) {
        WF = vecOfPairs[i]->GetPairHitCh(0) == 2
                 ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()
                 : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr();
        X1 = WF->IntegrateWaveForm(980, 1025);
        X2 = WF->IntegrateWaveForm(1060, 1280);
        XTot = WF->IntegrateWaveForm(980, 1280);
        X1 = X1 / XTot;
        X2 = X2 / XTot;
        // hE1E2->Fill(energy0, energyOther);
        // hDelT->Fill(vecOfPairs[i]->GetPairDelTime() / 1E3);
        hEMT->Fill(energy0, MT);
        // // if (energy0 > 15) {
        // hMTLam->Fill(MT, Lam);
        // }
        hX1X2->Fill(X1, X2);
      } else {
        badWFCount++;
      }
      // std::cout << "E0: " << energy0 << " : EO: " << energyOther <<
      // std::endl;

      if (keepGoing and energy0 > 1 and energy0 < 10 and MT > 1.9 and
          MT < 2.05 and preint < 2.1 and Lam > -5.4 and X1 > X1LowCut and
          X1 < X1HiCut and X2 < X2HiCut and X2 > X2LowCut) {

        // std::cout << "Energy: " << energy0 << " | MeanTime : " << MT
        //           << std::endl;
        WF = vecOfPairs[i]->GetPairHitCh(0) == 2
                 ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()
                 : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr();
        WFNorm = new digiAnalysis::WaveForm(
            WF->NormWaveForm().second); // WF->ScaleWaveForm(1.0 / energy0)
        waveformVector.push_back(*WFNorm);
        // WF->Plot();
        // std::cout << "Do you want to see the next waveform? (y/n): ";
        // std::getline(std::cin, userInput);
        // if (userInput != "y" && userInput != "Y") {
        //   keepGoing = false;
        // }
      }
      if (energy0 > 0 and energy0 < 4 and (MT < 1.0) and preint < 2.1 and
          Lam > -5.4 and !(X1 > X1LowCut and X1 < X1HiCut) and
          !(X2 < X2HiCut and X2 > X2LowCut)) { // or MT > 2.1

        // std::cout << "Energy: " << energy0 << " | MeanTime : " << MT
        //           << std::endl;
        WF = vecOfPairs[i]->GetPairHitCh(0) == 2
                 ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()
                 : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr();
        WFNorm = new digiAnalysis::WaveForm(
            WF->NormWaveForm().second); // WF->ScaleWaveForm(1.0 / energy0)

        wfVecNoise.push_back(*WFNorm);
        // if (keepGoing) {
        //   WF->Plot();
        //   std::cout << "Do you want to see the next waveform? (y/n): ";
        //   std::getline(std::cin, userInput);
        //   if (userInput != "y" && userInput != "Y") {
        //     keepGoing = false;
        //   }
        // }
      }
    }
  }

  UShort_t wfSz = WF->GetSize();
  std::cout << " Averaging " << waveformVector.size()
            << " waveforms for plotting" << std::endl;
  digiAnalysis::WaveForm WFAveraged(wfSz, waveformVector);
  WFAveraged.SetSmooth(40);
  WFAveraged.SetTracesFFT();
  WFAveraged.Plot(WFAveraged.NormWaveForm().second);
  cluster1 = WFAveraged.GetTracesSmooth();

  digiAnalysis::WaveForm WFAvNoise(wfSz, wfVecNoise);
  WFAvNoise.SetSmooth(40);
  cluster0 = WFAvNoise.GetTracesSmooth();
  double intValAv = WFAveraged.IntegrateWaveForm();
  std::cout << "intValAv: " << intValAv << std::endl;
  WFAveraged.Plot(WFAvNoise.NormWaveForm().second, "SAME_kRed");

  TFile f("output.root", "RECREATE");
  std::vector<double> waveformtrace;
  waveformtrace = WFAveraged.GetTraces();
  f.WriteObject(&waveformtrace, "trace");
  f.Close();

  // Testing the likelihood w.r.t averaged waveform
  double logLikelihood0, logLikelihood1, netLikelihood;
  double intValNoise, intValSig, intWF;
  // intValNoise = WFNoise->IntegrateWaveForm();
  // intValSig = WFSig->IntegrateWaveForm();
  intValSig = intValAv;
  intValNoise = WFAvNoise.IntegrateWaveForm();

  keepGoing = true;
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

    X1 = WF->IntegrateWaveForm(980, 1025);
    X2 = WF->IntegrateWaveForm(1060, 1280);
    XTot = WF->IntegrateWaveForm(980, 1280);
    X1 = X1 / XTot;
    X2 = X2 / XTot;

    energy0 = intWF * 0.01911 - 0.301;
    WFNorm = new digiAnalysis::WaveForm(
        WF->NormWaveForm().second); // WF->ScaleWaveForm(1.0 / energy0)

    WFNorm->SetSmooth(100);
    trace = WFNorm->GetTracesSmooth();
    logLikelihood0 = 0;
    logLikelihood1 = 0;
    netLikelihood = 0;
    for (int j = digiAnalysis::GateStart;
         j < digiAnalysis::GateStart + digiAnalysis::GateLenLong; j++) {
      // logLikelihood0 += trace[j] * TMath::Log(abs(trace[j] /
      // waveformtrace[j]));
      logLikelihood0 += trace[j] * TMath::Log(abs(trace[j] / cluster0[j]));
      logLikelihood1 += trace[j] * TMath::Log(abs(trace[j] / cluster1[j]));
    }
    // logLikelihood = logLikelihood + intValAv -
    //                 intWF * digiAnalysis::GateLenLong /
    //                     digiAnalysis::EvalNormFactor / energy0;
    logLikelihood0 = logLikelihood0 + intValNoise - intWF;
    logLikelihood1 = logLikelihood1 + intValSig - intWF;
    netLikelihood =
        (logLikelihood0 - logLikelihood1) / (logLikelihood0 + logLikelihood1);

    // std::cout << netLikelihood << std::endl;
    // if (i < 100) {
    //   std::cout << energy0 << " : " << logLikelihood / 10000.0 << " : "
    //             << WFNorm->IntegrateWaveForm() << " : "
    //             << intWF * digiAnalysis::GateLenLong /
    //                    digiAnalysis::EvalNormFactor / energy0
    //             << std::endl;
    // }

    double preint =
        vecOfPairs[i]->GetPairHitCh(0) == 2
            ? vecOfPairs[i]->GetHitPtr(0)->GetWFPtr()->IntegrateWaveForm(
                  0, digiAnalysis::GateStart)
            : vecOfPairs[i]->GetHitPtr(1)->GetWFPtr()->IntegrateWaveForm(
                  0, digiAnalysis::GateStart);
    preint /= digiAnalysis::GateStart;
    int chSel = vecOfPairs[i]->GetPairHitCh(0) == 2 ? 1 : 0;
    switch (vecOfPairs[i]->GetPairHitCh(chSel)) {
    case 0:
      energyOther = vecOfPairs[i]->GetPairHitEnergy(chSel) * 0.47835 - 11.71;
      break;
    case 2:
      break;
    case 3:
      energyOther = vecOfPairs[i]->GetPairHitEnergy(chSel) * 0.54222 - 10.3;
      break;
    case 4:
      energyOther = vecOfPairs[i]->GetPairHitEnergy(chSel) * 0.48345 + 2.395;
      break;

    default:
      break;
    }
    chSel = vecOfPairs[i]->GetPairHitCh(1) == 2 ? 1 : 0;
    Lam = vecOfPairs[i]->GetHitPtr(chSel)->GetWFPtr()->EvalNoisePar2(950, 1650);
    if (preint < 2.0) {
      if (MT > 1.8 and MT < 2.1 and X1 > X1LowCut and X1 < X1HiCut and
          X2 < X2HiCut and X2 > X2LowCut)
        hELL->Fill(energy0, netLikelihood);
      hMTLL->Fill(MT, netLikelihood);
      hEMTLL->Fill(energy0, MT, netLikelihood);
    }
    if (preint < 2.0 and MT > 1.8 and MT < 2.1 and X1 > X1LowCut and
        X1 < X1HiCut and X2 < X2HiCut and
        X2 > X2LowCut) { // and netLikelihood < 0.0
      if (vecOfPairs[i]->GetPairDelTime() / 1E3 > -20)
        hE1E2->Fill(energy0, energyOther);
      hDelT->Fill(vecOfPairs[i]->GetPairDelTime() / 1E3);
      // if (energy0 > 15) {
      hMTLam->Fill(MT, Lam);
      // }
    }
    if (preint < 2.0 and netLikelihood < 0.0) {
      hEMTCut->Fill(energy0, MT);
    }

    // if (MT > 1.8 and MT < 2.1 and netLikelihood > 0.0 and
    //     netLikelihood < 1.0 and keepGoing and energy0 < 3 and energy0 > 0 and
    //     X1 > X1LowCut and X1 < X1HiCut and X2 < X2HiCut and X2 > X2LowCut) {
    //   WF->Plot();
    //   // WFVec->push_back(*WF);
    //   std::cout << "Energy: " << energy0 << " LL: " << netLikelihood
    //             << std::endl;
    //   std::cout << "Do you want to see the next waveform? (y/n): ";
    //   std::getline(std::cin, userInput);
    //   if (userInput != "y" && userInput != "Y") {
    //     keepGoing = false;
    //   }
    // }
  }

  // WFAveraged

  // Testing the log likelihood on non-pair files
  keepGoing = true;
  TH2 *hELL_1 = new TH2F("hELL1", "hELL1", 800, 0, 80, 1000, -5, 5);
  TH2 *hEMT_1 = new TH2F("hEMT1", "hEMT1", 800, 0, 80, 1000, -5, 5);
  TH1F *hESpectra = new TH1F("hESpectra1", "Energy Spectra", 800, 0, 80);
  TH1F *hEEvalSpectra =
      new TH1F("hEEvalSpectra1", "Eval Energy Spectra", 800, 0, 80);
  TH1F *hELLCut = new TH1F("hELLCut1", "Energy Spectra", 800, 0, 80);
  TH2 *hMTLam_1 =
      new TH2F("MTLam_1", "MeanTime1 vs Lambda", 500, -4, 4, 2000, -10.0, 10.0);
  TH2 *hX1X2_Test =
      new TH2F("hX1X2_Test", "hX1X2_Test", 60, -0.2, 1, 50, -1, 1);
  std::vector<digiAnalysis::WaveForm> waveformVector1;

  std::string fname1 =
      "/home/kirtikesh/Analysis/DATA/extCoinc/"
      "NaI3124_17Jun26_NoSrc_1350V_2000V_1350V_1350V_Gain2_NoSplit_ExtTrig_"
      "Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_1000nsCoinc_"
      "2Vpp_Thresh_100lsb_WAVES_Sum_BLCorrected.root";
  digiAnalysis::Analysis an1(2, fname1, 0, 100000, 0);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVec1 =
      an1.GetSingleHitsVec();

  double X1LowCut_1 = 0.2, X1HiCut_1 = 0.4;
  double mX = -1.6, mC = 1.04;
  double X2LowCut_1 = 0.4, X2HiCut_1 = 0.6;
  double X1_1 = 0, X2_1 = 0, XTot_1 = 0;
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

    double lam = WF->EvalNoisePar2(1080, 1780);
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
        WF->NormWaveForm().second); // WF->ScaleWaveForm(1.0 / energy0)

    WFNorm1->SetSmooth(100);
    trace = WFNorm1->GetTracesSmooth();
    logLikelihood0 = 0;
    logLikelihood1 = 0;
    netLikelihood = 0;
    for (int j = digiAnalysis::GateStart;
         j < digiAnalysis::GateStart + digiAnalysis::GateLenLong; j++) {
      logLikelihood0 += trace[j] * TMath::Log(abs(trace[j] / cluster0[j]));
      logLikelihood1 += trace[j] * TMath::Log(abs(trace[j] / cluster1[j]));
    }
    logLikelihood0 = logLikelihood0 + intValNoise - intWF;
    logLikelihood1 = logLikelihood1 + intValSig - intWF;
    netLikelihood =
        (logLikelihood0 - logLikelihood1) / (logLikelihood0 + logLikelihood1);
    if (psd > 0.12 and psd < 1.0 and preint < 2 and X1_1 > X1LowCut_1 and
        X1_1 < X1HiCut_1 and X2_1 < mX * X1_1 + mC) {
      hMTLL1->Fill(MT, netLikelihood);
      hMTLam_1->Fill(MT, lam);
      hEMT_1->Fill(energy0, MT);
    }
    if (MT > 1.9 and MT < 2.1 and preint < 2.0 and psd > 0.12 and psd < 0.5 and
        lam > -5.4 and X1_1 > X1LowCut_1 and X1_1 < X1HiCut_1 and
        X2_1 < mX * X1_1 + mC and lam > -5.5) {
      hEEvalSpectra->Fill(energy0);
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
      if (netLikelihood < 0.0) {
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

  TCanvas *c5 = new TCanvas("c5", "E vs LL", 800, 600);
  hELL->Draw("COLZ");
  c5->Update();
  TCanvas *c6 = new TCanvas("c6", "MT vs LL", 800, 600);
  hMTLL->Draw("COLZ");
  c6->Update();
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
  hX1X2->Draw("LEGO2");
  TCanvas *c15 = new TCanvas("c15", "X1 vs X2 test", 800, 600);
  hX1X2_Test->Draw("LEGO2");

  fApp->Run();
  return 0;
}