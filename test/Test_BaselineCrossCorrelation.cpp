// ################################################################# //
//  This code is meant to create a template using shifted averaging  //
//           with cross correlation used to identify shift           //
// ################################################################# //
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
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <ratio>
#include <string>
#include <thread>
#include <vector>

int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;
  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //     "NaI13_20May26_1900_1345_Cs_Coinc144_WAVES_2/FILTERED/"
  //     "SDataF_NaI13_20May26_1900_1345_Cs_Coinc144_WAVES_2.root";

  //   std::string fname =
  //   "/media/kirtikesh/UbuntuFiles/NaI/SPE/May2026/NaI1_15May26_1900_NoSrc_WAVES/UNFILTERED/"
  //                       "Data_NaI1_15May26_1900_NoSrc_WAVES.root";

  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //     "NaI31_01June26_1345_1750_Cs_Thresh_300_30_WAVES_Coinc_144ns_LeadPit/"
  //     "FILTERED/"
  //     "SDataF_NaI31_01June26_1345_1750_Cs_Thresh_300_30_WAVES_Coinc_144ns_"
  //     "LeadPit.root";

  // std::string fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
  //                     "CoincidenceStudies/01JuneNoSrc/"
  //                     "NaI1342_02June26_1750_1345_1350_1350_NoSrc_Thresh_30_"
  //                     "300_WAVES_Coinc_144ns_LeadPit_Sum.root";

  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //     "NaI31_26May26_1345_1750_NoSrc_Thresh50_WAVES_NoCoinc_LeadPit_5/FILTERED/"
  //     "DataF_NaI31_26May26_1345_1750_NoSrc_Thresh50_WAVES_NoCoinc_LeadPit_5."
  //     "root";

  // std::string fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
  //                     "CoincidenceStudies/01JuneNoSrc/"
  //                     "NaI1342_04June26_1750_1345_1350_1350_NoSrc_Thresh_120_"
  //                     "300_WAVES_Singles_LeadPit_45/FILTERED/"
  //                     "DataF_NaI1342_04June26_1750_1345_1350_1350_NoSrc_Thresh_"
  //                     "120_300_WAVES_Singles_LeadPit_45.root";

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/extCoinc/"
      "NaI3124_16Jun26_AmSrc_1350V_2000V_1350V_1350V_Gain2_NoSplit_ExtTrig_"
      "Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_1000nsCoinc_"
      "2Vpp_Thresh_100lsb_WAVES_21/FILTERED/"
      "DataF_NaI3124_16Jun26_AmSrc_1350V_2000V_1350V_1350V_Gain2_NoSplit_"
      "ExtTrig_Thresh75_DelayCoincLogic_PGate160ns_Delay240ns_DGate600ns_"
      "1000nsCoinc_2Vpp_Thresh_100lsb_WAVES_21.root";

  // Read to singleHits
  digiAnalysis::Analysis an(2, fname, 0, 200000, 1);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();

  std::vector<double> accumulateTrace(hitsVector[0]->GetWFPtr()->GetSize(),
                                      0.0);
  int numTraces = 0;

  int iterStart = 0;
  while (iterStart < hitsVector.size()) {
    if (hitsVector[iterStart]->GetEnergy() > 00 and
        hitsVector[iterStart]->GetEnergy() < 30 and
        hitsVector[iterStart]->GetPSD() > 0. and
        // hitsVector[iterStart]->GetPSD() > 0.5 and
        hitsVector[iterStart]->GetMeanTime() < 1.6 and
        hitsVector[iterStart]->GetMeanTime() > 0.0 and
        hitsVector[iterStart]->GetWFPtr()->IntegrateWaveForm(
            0, digiAnalysis::GateStart) *
                1.0 / digiAnalysis::GateStart <
            2)
      break;
    iterStart++;
  }

  digiAnalysis::WaveForm *WF = hitsVector[iterStart]->GetWFPtr();
  // WF->SetSmooth(500);
  WF->SetSmooth(16, "MovA");
  std::vector<double> tracePrimary = WF->GetTracesSmooth();
  // for (int iWF = 60; iWF < WF->GetSize(); iWF++) {
  //   if (abs(tracePrimary[iWF] - tracePrimary[iWF - 4]) > 3 or
  //       abs(tracePrimary[iWF] - tracePrimary[iWF + 4]) > 3 or
  //       tracePrimary[iWF] > 6)
  //     tracePrimary[iWF] = tracePrimary[iWF - 30];
  // }
  for (int iWF = WF->GetSize() - 50; iWF >= 0; iWF--) {
    if (abs(tracePrimary[iWF] - tracePrimary[iWF - 4]) > 3 or
        abs(tracePrimary[iWF] - tracePrimary[iWF + 4]) > 3 or
        abs(tracePrimary[iWF]) > 12) {
      tracePrimary[iWF] = tracePrimary[iWF + 30];
    }
  }
  WF->Plot(WF->GetTracesSmooth(), tracePrimary); // this line is only for
  //   setting the cutoff for removing the single pulses
  WF->SetTracesFFT(tracePrimary);
  std::vector<double> trFFT_AmpPrim = WF->GetTracesFFT();
  std::vector<double> trFFT_PhsPrim = WF->GetTracesFFTPhase();

  std::vector<double> traceSecondary;
  std::vector<double> traceSecondaryOrig, trCorr;
  std::vector<double> trFFT_AmpSecn, trFFT_PhsSecn;
  std::vector<double> trFFT_AmpCorr, trFFT_PhsCorr;
  digiAnalysis::WaveForm *WF1 = nullptr;
  for (int hititer = 0; hititer < hitsVector.size();
       hititer++) // hitsVector.size()
  {
    if (hititer % 1000 == 0)
      std::cout << "Processing hit " << hititer << std::endl;

    // choose events whose waveforms is expected to give good baseline variation
    if (hitsVector[hititer]->GetEnergy() > 00 and
        hitsVector[hititer]->GetEnergy() < 100 and
        hitsVector[hititer]->GetPSD() > 0. and
        // hitsVector[hititer]->GetPSD() > 0.5 and
        hitsVector[hititer]->GetMeanTime() < 1.6 and
        hitsVector[hititer]->GetMeanTime() > 0.0 and
        hitsVector[iterStart]->GetWFPtr()->IntegrateWaveForm(
            0, digiAnalysis::GateStart) *
                1.0 / digiAnalysis::GateStart <
            2) {
      WF1 = hitsVector[hititer]->GetWFPtr();
      // WF1->SetSmooth(500);
      WF1->SetSmooth(16, "MovA");
      traceSecondary = WF1->GetTracesSmooth();
      traceSecondaryOrig = WF1->GetTraces();
      // for (int iWF = 60; iWF < WF1->GetSize(); iWF++) {
      //   if (abs(traceSecondary[iWF] - traceSecondary[iWF - 4]) > 3 or
      //       abs(traceSecondary[iWF] - traceSecondary[iWF + 4]) > 3 or
      //       traceSecondary[iWF] > 6) {
      //     traceSecondary[iWF] = traceSecondary[iWF - 30];
      //     traceSecondaryOrig[iWF] = traceSecondaryOrig[iWF - 30];
      //   }
      // }

      // clean it to remove any fast fluctuations or large deviations from
      // baseline
      for (int iWF = WF1->GetSize() - 50; iWF >= 0; iWF--) {
        if (abs(traceSecondary[iWF] - traceSecondary[iWF - 4]) > 3 or
            abs(traceSecondary[iWF] - traceSecondary[iWF + 4]) > 3 or
            abs(traceSecondary[iWF]) > 12) {
          traceSecondary[iWF] = traceSecondary[iWF + 30];
          traceSecondaryOrig[iWF] = traceSecondaryOrig[iWF + 30];
        }
      }
      // WF1->Plot(traceSecondary);
      // std::this_thread::sleep_for(std::chrono::milliseconds(200));
      WF1->SetTracesFFT(traceSecondary);
      trFFT_AmpSecn = WF1->GetTracesFFT();
      trFFT_PhsSecn = WF1->GetTracesFFTPhase();

      for (int iterWF = 0; iterWF < WF->GetSize(); iterWF++) {
        // ensuring that the integrated waveform is 0
        traceSecondary[iterWF] -= trFFT_AmpSecn[0] / WF->GetSize();
        if (hititer == 1)
          tracePrimary[iterWF] -= trFFT_AmpPrim[0] / WF->GetSize();
      }

      // Evaluate the cross correlation and find the maxima to get the shift
      trFFT_AmpCorr.clear();
      trFFT_PhsCorr.clear();
      for (int iterFFT = 0; iterFFT < trFFT_AmpSecn.size(); iterFFT++) {
        trFFT_AmpCorr.push_back(trFFT_AmpPrim[iterFFT] *
                                trFFT_AmpSecn[iterFFT]);
        trFFT_PhsCorr.push_back(trFFT_PhsPrim[iterFFT] -
                                trFFT_PhsSecn[iterFFT]);
      }
      trFFT_AmpCorr[0] = 0;
      trCorr = WF1->EvalIFFT(trFFT_AmpCorr, trFFT_PhsCorr);
      int index = 0;
      auto maxIt = std::max_element(trCorr.begin(), trCorr.end());
      if (maxIt != trCorr.end())
        index = std::distance(trCorr.begin(), maxIt);
      int val = 0;
      int shift = WF->GetSize() - index;
      // std::cout << "shift is: " << shift << std::endl;

      // subtract the shifted waveform to check the removal of features.
      // // WF1->SetSmooth(16, "MovA");
      // traceSecondary = WF1->GetTraces();
      // // WF->SetSmooth(16, "MovA");
      // tracePrimary = WF->GetTraces();
      // for (int iterWF = 0; iterWF < WF->GetSize(); iterWF++)
      // {
      //   val = iterWF - shift > 0 ? iterWF - shift : WF->GetSize() + iterWF -
      //   shift; traceSecondary[iterWF] = traceSecondary[iterWF] -
      //   tracePrimary[val];
      // }
      // WF1->SetTracesFFT(traceSecondary);

      // // implementation of high pass filter to remove low frequency noise
      // in the subtracted waveform std::vector<double> trFFT_AmpRes =
      // WF1->GetTracesFFT(); std::vector<double> trFFT_AmpPhs =
      // WF1->GetTracesFFTPhase(); int noisestart = 500; int noisestop =
      // 1000;
      // int sigstop = trFFT_AmpRes.size() - 1;
      // int highpass = 5;
      // double sum = std::accumulate(trFFT_AmpRes.begin() + 500,
      // trFFT_AmpRes.begin() + 1000, 0.0); for (int iterFFT = 0; iterFFT <
      // sigstop; iterFFT++)
      // {
      //   if (iterFFT < highpass)
      //     trFFT_AmpRes[iterFFT] = 0;
      //   else
      //     trFFT_AmpRes[iterFFT] -= sum / (noisestop - noisestart);
      // }

      // std::fill(trFFT_AmpRes.begin() + sigstop, trFFT_AmpRes.end(), 0);
      // WF1->Plot(WF1->EvalIFFT(trFFT_AmpRes, trFFT_AmpPhs), trFFT_AmpRes);

      //   WF->Plot(traceSecondary, tracePrimary);
      //   WF->Plot(traceSecondary, trCorr);
      // WF->Plot(traceSecondary, WF1->GetTracesFFT());
      //   WF->Plot(tracePrimary, trFFT_AmpPrim);
      // WF1->Plot(WF1->GetTracesSmooth(), traceSecondary);
      //   WF1->Plot(traceSecondary, "SAME_KGreen");

      // For accumulation, use the original waveforms, not the smoothed
      // ones. traceSecondary = WF1->GetTraces();
      for (int iterWF = 0; iterWF < WF1->GetSize(); iterWF++) {
        val = iterWF - shift > 0 ? iterWF - shift
                                 : WF1->GetSize() + iterWF - shift;
        accumulateTrace[val] += traceSecondaryOrig[iterWF];
      }
      // if (numTraces % 50 == 0)
      // WF1->Plot(accumulateTrace, WF1->GetTraces());
      // std::this_thread::sleep_for(std::chrono::milliseconds(1000));
      // if (numTraces != 0 and numTraces < 100)
      //   WF1->Plot(accumulateTrace, "SAME");
      numTraces++;
      // WF1->Plot(accumulateTrace);
    }
  }
  std::transform(accumulateTrace.begin(), accumulateTrace.end(),
                 accumulateTrace.begin(),
                 [numTraces](double val) { return val / numTraces; });
  WF1->Plot(accumulateTrace);
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  digiAnalysis::WaveForm AccumulatedWF(accumulateTrace);
  AccumulatedWF.SetTracesFFT();
  // AccumulatedWF.Plot(AccumulatedWF.GetTraces(),
  // AccumulatedWF.GetTracesFFT());

  int filterSz = WF1->GetSize() / 2 + 1;
  int filterCutOff = 200; // this corresponds in frequency to filterCutOff *
                          // (500/NSampleSPE) MHz
  int filterFlatRange = 120;
  int filterGaussSigma = (filterCutOff - filterFlatRange) / 5;
  std::vector<Double_t> filter(filterSz);
  for (int iter = 0; iter < filterSz; iter++) {
    if (iter < filterFlatRange) {
      filter[iter] = 1.0;
    } else if (iter - filterFlatRange < 5 * filterGaussSigma) {
      filter[iter] = TMath::Gaus(iter, filterFlatRange, filterGaussSigma);
    } else {
      filter[iter] = 0;
    }
  }
  std::vector<double> trFFT_AmpAccum = AccumulatedWF.GetTracesFFT();
  std::transform(trFFT_AmpAccum.begin(), trFFT_AmpAccum.end(), filter.begin(),
                 trFFT_AmpAccum.begin(), std::multiplies<>());
  AccumulatedWF.Plot(
      AccumulatedWF.EvalIFFT(trFFT_AmpAccum, AccumulatedWF.GetTracesFFTPhase()),
      trFFT_AmpAccum);
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  std::vector<double> trFullRange =
      AccumulatedWF.EvalIFFT(trFFT_AmpAccum, AccumulatedWF.GetTracesFFTPhase());
  // std::vector<double> traceOnePeriod(trFullRange.begin() + 1428,
  //                                    trFullRange.begin() + 6428);
  std::vector<double> traceOnePeriod(trFullRange.begin(), trFullRange.end());
  digiAnalysis::WaveForm AccumulatedWF_5k(traceOnePeriod);
  AccumulatedWF_5k.SetTracesFFT();
  AccumulatedWF_5k.Plot();
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  // AccumulatedWF_5k.GetTracesFFT());

  std::string outfname = "/home/kirtikesh/Analysis/DATA/extCoinc/"
                         "Baseline_17June_singlePeriod_Ch2.root";
  TFile *fout = TFile::Open(outfname.c_str(), "RECREATE");
  TTree *t = new TTree("baseline", "b1aseline");
  digiAnalysis::WaveForm WFSPE;
  t->Branch("baseline", &traceOnePeriod);
  t->Fill();
  t->Write();
  fout->Close();

  fApp->Run();
}