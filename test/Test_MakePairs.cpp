#include "Analysis.h"
#include "Pair.h"
#include "TF1.h"
#include "TMath.h"
#include "WaveForm.h"
#include "includes.hh"
#include "singleHits.h"
#include <TApplication.h>
#include <cmath>
#include <iostream>
int main(int argc, char *argv[])
{
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;

  std::string fname = "/media/kirtikesh/KKBlack/NaI/"
                      "run_Cs_FAGain_2_10_CFDTHR_15_10_Mode_EXT_TRG_FREEWRITE_"
                      "SignalDelay_50ns_Aug26/FILTERED/"
                      "DataF_run_Cs_FAGain_2_10_CFDTHR_15_10_Mode_EXT_TRG_"
                      "FREEWRITE_SignalDelay_50ns_Aug26.root";

  digiAnalysis::Analysis an(fname, 0, 100000, 0);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  an.CreatePairs(5, 9);

  // std::vector<digiAnalysis::singleHits *> hitsVectorCh0 =
  //     an.GetSingleHitsVec(7);
  // std::vector<digiAnalysis::singleHits *> hitsVectorCh1 =
  //     an.GetSingleHitsVec(5);

  ULong64_t i = 0;
  ULong64_t j = 0;

  digiAnalysis::WaveForm *WF = nullptr;
  std::vector<digiAnalysis::WaveForm> waveformVector;
  std::vector<digiAnalysis::Pair *> vecOfPairs = an.GetPairsVec();

  // Long64_t Ch0Time = hitsVectorCh0[0]->GetTimestamp();
  // Long64_t Ch1Time = hitsVectorCh1[0]->GetTimestamp();
  // Long64_t TDiffCurr = Ch0Time - Ch1Time;
  // Long64_t TDiffMin = 0;
  // int EnergyNearestEvt = 0;
  // digiAnalysis::Pair currPair;
  // UShort_t ChEnergy = 0;

  TH1F *hICDelT =
      new TH1F("hICDelT", "inter channel nearest event time difference", 80000,
               -4000000, 4000000);
  TH1F *hECh0 =
      new TH1F("hECh0", "Energy Spectra of lower channel", 8192, 0, 8192);
  TH1F *hECh1 =
      new TH1F("hECh1", "Energy Spectra of higher channel", 8192, 0, 8192);
  TH1F *hE =
      new TH1F("hE", "Sum Energy Spectra of Coinc Events", 16384, 0, 16384);

  // while (i < hitsVectorCh0.size()) {
  //   Ch0Time = hitsVectorCh0[i]->GetTimestamp();
  //   TDiffMin = 100000000000000;
  //   while (j < hitsVectorCh1.size()) {
  //     Ch1Time = hitsVectorCh1[j]->GetTimestamp();
  //     TDiffCurr = Ch1Time - Ch0Time;
  //     if (fabs(TDiffCurr) < fabs(TDiffMin)) {
  //       EnergyNearestEvt = hitsVectorCh1[j]->GetEnergy();
  //       TDiffMin = TDiffCurr;
  //       j = j + 1;
  //     } else {
  //       hICDelT->Fill(TDiffMin / 1000.0);
  //       if (fabs(TDiffMin / 1000) < 30) {
  //         hECh0->Fill(hitsVectorCh0[i]->GetEnergy());
  //         hECh1->Fill(EnergyNearestEvt);

  //         currPair.ClearPair();
  //         currPair.SetPair(hitsVectorCh0[i], hitsVectorCh1[j - 1]);
  //         vecOfPairs.push_back(new digiAnalysis::Pair(currPair));
  //         // if (fabs(currPair.GetPairHitEnergy(0) - 400) < 200 or
  //         //     fabs(currPair.GetPairHitEnergy(1) - 400) < 200) {
  //         //   ChEnergy = fabs(currPair.GetPairHitEnergy(1) - 400) < 200
  //         //                  ? currPair.GetPairHitEnergy(0)
  //         //                  : currPair.GetPairHitEnergy(1);
  //         //   hE->Fill(ChEnergy);
  //         // }
  //         hE->Fill(currPair.GetPairEnergy());
  //       }
  //       i = i + 1;
  //       j = (j >= 2) ? j - 2 : j;
  //       break;
  //     }
  //   }
  //   if (i >= hitsVectorCh0.size() - 1 or j >= hitsVectorCh1.size() - 1) {
  //     break;
  //   }
  // }

  std::cout << "Loop Exited, " << vecOfPairs.size() << " pairs found"
            << std::endl;

  vecOfPairs[0]->Print();
  vecOfPairs[1]->Print();

  bool plotnext = true;
  i = 0;
  UShort_t CheckEnergy = 4000;
  UShort_t CheckWidth = 8000;
  std::string userInput;
  digiAnalysis::singleHits *hit;

  while (i < vecOfPairs.size())
  {
    if (fabs(vecOfPairs[i]->GetPairHitEnergy(1) - CheckEnergy) < CheckWidth)
    {
      hECh0->Fill(vecOfPairs[i]->GetPairHitEnergy(0));
      hECh1->Fill(vecOfPairs[i]->GetPairHitEnergy(1));
      if (plotnext)
      {
        hit = vecOfPairs[i]->GetHit(1);
        WF = hit->GetWFPtr();
        if (WF)
        {
          waveformVector.push_back(*WF);
        }
        vecOfPairs[i]->Print();
        WF->SetSmooth(16);
        WF->Plot();
        std::cout << "Do you want to see the next waveform? (y/n): ";
        std::getline(std::cin, userInput);
        if (userInput != "y" && userInput != "Y")
        {
          plotnext = false;
        }
      }
    }

    hE->Fill(vecOfPairs[i]->GetPairEnergy());
    hICDelT->Fill(vecOfPairs[i]->GetPairDelTime());
    i++;
  }
  UShort_t wfSz = WF->GetSize();
  std::cout << "Got size of waveforms" << wfSz << std::endl;
  if (waveformVector.size() > 0)
  {
    digiAnalysis::WaveForm *WFAveraged =
        new digiAnalysis::WaveForm(wfSz, waveformVector);
    WFAveraged->Plot();
  }

  TCanvas *c1 = new TCanvas("c1", "timeDiff", 800, 600);
  hICDelT->Draw("HIST");
  TCanvas *c2 = new TCanvas("c2", "LowerChannel Energy", 800, 600);
  hECh0->Draw("HIST");
  TCanvas *c3 = new TCanvas("c3", "HigherChannel Energy", 800, 600);
  hECh1->Draw("HIST");
  TCanvas *c4 = new TCanvas("c4", "Sum Energy", 800, 600);
  hE->Draw("HIST");

  fApp->Run();
  return 0;
}
