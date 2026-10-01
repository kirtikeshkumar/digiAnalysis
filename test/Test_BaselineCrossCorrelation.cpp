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

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/extCoincSep/1800V/Calib/"
      "NaI1_AmSrc_1800_1337_1350_WAVES_FILTERED_NoSplitSignal_Gain_12_"
      "Acquisition_4_ExtTrig_Threshold1LSB_160nsPromt_240nsDelay_800nsLong_"
      "1496nsCoinc_FreeWrites_Ch248Singles_15/FILTERED/"
      "DataF_NaI1_AmSrc_1800_1337_1350_WAVES_FILTERED_NoSplitSignal_Gain_12_"
      "Acquisition_4_ExtTrig_Threshold1LSB_160nsPromt_240nsDelay_800nsLong_"
      "1496nsCoinc_FreeWrites_Ch248Singles_15.root";

  // Read to singleHits
  digiAnalysis::Analysis an(0, fname, 00000, 200000, 1);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();

  std::vector<double> accumulateTrace(hitsVector[0]->GetWFPtr()->GetSize(),
                                      0.0);
  int numTraces = 0;

  int iterStart = 0;
  while (iterStart < hitsVector.size()) {
    if (hitsVector[iterStart]->GetEnergy() > 00 and
        hitsVector[iterStart]->GetEnergy() < 200 and
        hitsVector[iterStart]->GetPSD() > 0.0 and
        // hitsVector[iterStart]->GetPSD() > 0.5 and
        hitsVector[iterStart]->GetMeanTime() < 1.6 and
        hitsVector[iterStart]->GetMeanTime() > 0.0 and
        fabs(hitsVector[iterStart]->GetWFPtr()->IntegrateWaveForm(0, 100) *
             1.0 / (100)) < 2.0)
      break;
    iterStart++;
  }

  digiAnalysis::WaveForm *WF = hitsVector[iterStart]->GetWFPtr();
  // WF->SetSmooth(500);
  WF->SetSmooth(16, "MovA");
  std::vector<double> tracePrimary = WF->GetTracesSmooth();
  for (int iWF = 60; iWF < WF->GetSize(); iWF++) {
    if (abs(tracePrimary[iWF] - tracePrimary[iWF - 4]) > 20 or
        abs(tracePrimary[iWF] - tracePrimary[iWF + 4]) > 20 or
        tracePrimary[iWF] > 30)
      tracePrimary[iWF] = tracePrimary[iWF - 30];
  }
  // for (int iWF = WF->GetSize() - 50; iWF >= 0; iWF--) {
  //   if (abs(tracePrimary[iWF] - tracePrimary[iWF - 4]) > 20 or
  //       abs(tracePrimary[iWF] - tracePrimary[iWF + 4]) > 20 or
  //       (tracePrimary[iWF]) > 30) {
  //     tracePrimary[iWF] = tracePrimary[iWF + 30];
  //   }
  // }
  // WF->Plot(WF->GetTracesSmooth(), tracePrimary); // this line is only for
  //   setting the cutoff for removing the single pulses
  // std::this_thread::sleep_for(std::chrono::milliseconds(10000));
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
        hitsVector[hititer]->GetEnergy() < 200 and
        hitsVector[hititer]->GetPSD() > 0.0 and
        // hitsVector[hititer]->GetPSD() > 0.5 and
        hitsVector[hititer]->GetMeanTime() < 1.6 and
        hitsVector[hititer]->GetMeanTime() > 0.0 and
        fabs(hitsVector[iterStart]->GetWFPtr()->IntegrateWaveForm(0, 100) *
             1.0 / (100)) < 2.0) {
      WF1 = hitsVector[hititer]->GetWFPtr();
      // WF1->SetSmooth(500);
      WF1->SetSmooth(16, "MovA");
      traceSecondary = WF1->GetTracesSmooth();
      traceSecondaryOrig = WF1->GetTraces();
      for (int iWF = 60; iWF < WF1->GetSize(); iWF++) {
        if (abs(traceSecondary[iWF] - traceSecondary[iWF - 4]) > 20 or
            abs(traceSecondary[iWF] - traceSecondary[iWF + 4]) > 20 or
            traceSecondary[iWF] > 30) {
          traceSecondary[iWF] = traceSecondary[iWF - 30];
          traceSecondaryOrig[iWF] = traceSecondaryOrig[iWF - 30];
        }
      }

      // clean it to remove any fast fluctuations or large deviations from
      // baseline
      // for (int iWF = WF1->GetSize() - 50; iWF >= 0; iWF--) {
      //   if (abs(traceSecondary[iWF] - traceSecondary[iWF - 4]) > 20 or
      //       abs(traceSecondary[iWF] - traceSecondary[iWF + 4]) > 20 or
      //       (traceSecondary[iWF]) > 30) {
      //     traceSecondary[iWF] = traceSecondary[iWF + 30];
      //     traceSecondaryOrig[iWF] = traceSecondaryOrig[iWF + 30];
      //   }
      // }

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
      for (int iterWF = 0; iterWF < WF->GetSize(); iterWF++) {
        val = iterWF - shift > 0 ? iterWF - shift
                                 : WF->GetSize() + iterWF - shift;
        traceSecondary[val] = traceSecondaryOrig[iterWF];
      }
      digiAnalysis::WaveForm *WF2 = new digiAnalysis::WaveForm(traceSecondary);
      WF2->SetSmooth();
      if (hititer and numTraces % 5 == 0) {
        WF1->Plot(WF1->GetTraces(), traceSecondary);
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
      }

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
      // if (numTraces % 5 == 0) {
      //   WF1->Plot(accumulateTrace, WF2->GetTracesSmooth());
      //   std::this_thread::sleep_for(std::chrono::milliseconds(100));
      // }
      // if (numTraces != 0 and numTraces < 100)
      //   WF1->Plot(accumulateTrace, "SAME");
      numTraces++;
      // WF1->Plot(accumulateTrace);
    }
  }
  std::transform(accumulateTrace.begin(), accumulateTrace.end(),
                 accumulateTrace.begin(),
                 [numTraces](double val) { return val / numTraces; });
  // WF1->Plot(accumulateTrace);
  // std::this_thread::sleep_for(std::chrono::milliseconds(100));
  digiAnalysis::WaveForm AccumulatedWF(accumulateTrace);
  AccumulatedWF.SetTracesFFT();
  AccumulatedWF.Plot(AccumulatedWF.GetTraces(), AccumulatedWF.GetTracesFFT());

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
  // AccumulatedWF.Plot(
  //     AccumulatedWF.EvalIFFT(trFFT_AmpAccum,
  //     AccumulatedWF.GetTracesFFTPhase()), trFFT_AmpAccum);
  // std::this_thread::sleep_for(std::chrono::milliseconds(100));
  std::vector<double> trFullRange =
      AccumulatedWF.EvalIFFT(trFFT_AmpAccum, AccumulatedWF.GetTracesFFTPhase());
  // std::vector<double> traceOnePeriod(trFullRange.begin() + 1428,
  //                                    trFullRange.begin() + 6428);
  std::vector<double> traceOnePeriod(trFullRange.begin(), trFullRange.end());
  digiAnalysis::WaveForm AccumulatedWF_5k(traceOnePeriod);
  AccumulatedWF_5k.SetTracesFFT();
  // AccumulatedWF_5k.Plot();
  // std::this_thread::sleep_for(std::chrono::milliseconds(100));
  // AccumulatedWF_5k.GetTracesFFT());

  std::string outfname = "/home/kirtikesh/Analysis/DATA/extCoincSep/1800V/"
                         "Baseline_29Sep_singlePeriod_Ch0.root";
  TFile *fout = TFile::Open(outfname.c_str(), "RECREATE");
  TTree *t = new TTree("baseline", "b1aseline");
  digiAnalysis::WaveForm WFSPE;
  t->Branch("baseline", &traceOnePeriod);
  t->Fill();
  t->Write();
  fout->Close();

  fApp->Run();
}