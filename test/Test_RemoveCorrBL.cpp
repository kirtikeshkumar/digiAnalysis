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
  std::string fnameBL =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
      "Baseline_04June_singlePeriod_Ch0.root";
  TFile *fp = new TFile(fnameBL.c_str(), "READ");
  TTree *tr = (TTree *)fp->Get("baseline");
  std::vector<double> *BL = nullptr;
  digiAnalysis::WaveForm baselineWF;
  if (tr) {
    tr->SetBranchAddress("baseline", &BL);
    Long64_t nentries = tr->GetEntries();
    tr->GetEntry(0);
    std::cout << nentries << " Baseline retrieved. Length is " << BL->size()
              << std::endl;
    baselineWF.SetWaveForm(*BL);
  } else {
    return 1;
  }
  std::vector<double> trFFT_AmpPrim, trFFT_PhsPrim;
  baselineWF.SetTracesFFT();
  trFFT_AmpPrim = baselineWF.GetTracesFFT();
  trFFT_PhsPrim = baselineWF.GetTracesFFTPhase();
  trFFT_AmpPrim[0] = 0;
  baselineWF.SetWaveForm(baselineWF.EvalIFFT(trFFT_AmpPrim, trFFT_PhsPrim));

  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //     "NaI13_20May26_1900_1345_Cs_Coinc144_WAVES_2/FILTERED/"
  //     "SDataF_NaI13_20May26_1900_1345_Cs_Coinc144_WAVES_2.root";

  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //     "NaI31_01June26_1345_1750_Cs_Thresh_300_30_WAVES_Coinc_144ns_LeadPit/"
  //     "FILTERED/"
  //     "SDataF_NaI31_01June26_1345_1750_Cs_Thresh_300_30_WAVES_Coinc_144ns_"
  //     "LeadPit.root";

  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //     "NaI31_26May26_1345_1750_NoSrc_Thresh50_WAVES_NoCoinc_LeadPit_5/FILTERED/"
  //     "DataF_NaI31_26May26_1345_1750_NoSrc_Thresh50_WAVES_NoCoinc_LeadPit_5."
  //     "root";

  // std::string fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
  //                     "CoincidenceStudies/01JuneNoSrc/"
  //                     "NaI1342_June26_1750_1345_1350_1350_NoSrc_Thresh_2_30_"
  //                     "300_WAVES_Coinc_144ns_LeadPit_Sum.root";

  // std::string fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
  //                     "CoincidenceStudies/01JuneNoSrc/"
  //                     "NaI1342_June26_1750_1345_1350_1350_NoSrc_Thresh_15-30_"
  //                     "300_WAVES_Coinc_144ns_LeadPit_Sum.root";

  std::string fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
                      "CoincidenceStudies/01JuneNoSrc/CalibrationFiles/"
                      "DataF_NaI1_05June26_1900_CoEuSrc_Thresh_100_WAVES_"
                      "Singles_LeadPit_75.root";

  // std::string fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
  //                     "CoincidenceStudies/01JuneNoSrc/"
  //                     "NaI1342_04June26_1750_1345_1350_1350_NoSrc_Thresh_120_"
  //                     "300_WAVES_Singles_LeadPit_45/FILTERED/"
  //                     "DataF_NaI1342_04June26_1750_1345_1350_1350_NoSrc_Thresh_"
  //                     "120_300_WAVES_Singles_LeadPit_45.root";

  // Read to singleHits
  digiAnalysis::Analysis an(fname, 00000, 00000, 0);
  // digiAnalysis::Analysis an2(0, fname, 0, 00000, 0);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  // std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector2 =
  //     an2.GetSingleHitsVec();

  // to write the subtracted BL to file
  UShort_t Channel;
  ULong64_t Timestamp;
  UShort_t Board;
  UShort_t Energy;
  UShort_t EnergyShort;
  TArrayS *Samples = nullptr;
  std::string writefname = fname;
  size_t posRoot = writefname.find(".root");
  if (posRoot != std::string::npos) {
    writefname.insert(posRoot, "_BLCorrected");
  }
  TFile *fout = new TFile(writefname.c_str(), "RECREATE");
  TTree *Data_F = new TTree("Data_F", "Filtered Data");
  Data_F->Branch("Channel", &Channel, "Channel/s");
  Data_F->Branch("Timestamp", &Timestamp, "Timestamp/l");
  Data_F->Branch("Board", &Board, "Board/s");
  Data_F->Branch("Energy", &Energy, "Energy/s");
  Data_F->Branch("EnergyShort", &EnergyShort, "EnergyShort/s");
  Data_F->Branch("Samples", &Samples);
  double baselineVal;

  int blLen = BL->size();
  int WFLen = hitsVector[0]->GetWFPtr()->GetSize();
  int WFCutStart = 2500;
  std::cout << WFLen << " : " << blLen << std::endl;

  std::vector<double> zeroTraceStart(WFCutStart, 0.0);
  std::vector<double> zeroTraceEnd(blLen - WFLen, 0.0);
  int numTraces = 0;

  std::vector<double> traceBL, traceSecondary;
  std::vector<double> traceSegment, tracePadded;
  std::vector<double> traceSecondaryOrig, trCorr;
  std::vector<double> trFFT_AmpSecn, trFFT_PhsSecn;
  std::vector<double> trFFT_AmpCorr, trFFT_PhsCorr;
  digiAnalysis::WaveForm *WF1 = nullptr;
  for (int hititer = 0; hititer < hitsVector.size();
       hititer++) // hitsVector.size()
  {
    if (hititer % 10000 == 0)
      std::cout << "Processing hit " << hititer << std::endl;
    if (hitsVector[hititer]->GetChNum() ==
        0 // and
          //  hitsVector[hititer]->GetEnergy() < 2000
    ) {

      //     if (hitsVector[hititer]->GetEnergy() < 100) {
      WF1 = hitsVector[hititer]->GetWFPtr();
      //       // WF1->SetSmooth(500);
      WF1->SetSmooth(16, "MovA");
      traceSecondary = WF1->GetTracesSmooth();
      traceSecondaryOrig = WF1->GetTraces();
      for (int iWF = 60; iWF < WF1->GetSize() - 5; iWF++) {
        if (abs(traceSecondary[iWF] - traceSecondary[iWF - 4]) > 3 or
            abs(traceSecondary[iWF] - traceSecondary[iWF + 4]) > 3 or
            traceSecondary[iWF] > 6) {
          traceSecondary[iWF] = traceSecondary[iWF - 30];
          traceSecondaryOrig[iWF] = traceSecondaryOrig[iWF - 30];
        }
      }
      // tracePadded.clear();
      // tracePadded.insert(tracePadded.end(), zeroTraceStart.begin(),
      //                    zeroTraceStart.end());
      // tracePadded.insert(tracePadded.end(), traceSecondary.begin() +
      // WFCutStart,
      //                    traceSecondary.end());
      // tracePadded.insert(tracePadded.end(), zeroTraceEnd.begin(),
      //                    zeroTraceEnd.end());
      // std::cout << "Len Padded: " << tracePadded.size() << std::endl;
      // WF1->Plot(traceSecondary, WF1->GetTracesSmooth());
      // std::this_thread::sleep_for(std::chrono::milliseconds(1000));

      // std::cout << "trace size: " << traceSecondary.size() << std::endl;
      WF1->SetTracesFFT(traceSecondary);
      trFFT_AmpSecn.clear();
      trFFT_PhsSecn.clear();
      trFFT_AmpSecn = WF1->GetTracesFFT();
      trFFT_PhsSecn = WF1->GetTracesFFTPhase();
      trFFT_AmpSecn[0] = 0;
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
      auto minIt = std::min_element(trCorr.begin(), trCorr.end());
      if (minIt != trCorr.end())
        index = std::distance(trCorr.begin(), minIt);
      int val = 0;
      int shift = baselineWF.GetSize() - index;
      // std::cout << "shift is: " << shift << std::endl;

      // traceSecondary = WF1->GetTraces();
      traceBL = baselineWF.GetTraces();

      // Evaluate baseline shift before subtraction.
      int starttime = 2500, stoptime = 3500;
      double blstart = starttime - shift > 0
                           ? starttime - shift
                           : baselineWF.GetSize() + starttime - shift;
      double blstop = stoptime - shift > 0
                          ? stoptime - shift
                          : baselineWF.GetSize() + stoptime - shift;
      if (blstop < blstart and blstop > 500) {
        starttime = starttime + (baselineWF.GetSize() - blstart);
        blstart = 0;
      }
      if (blstop < blstart and blstop < 500) {
        stoptime = stoptime - blstop;
        blstop = baselineWF.GetSize() - 1;
      }
      double blintegrate = baselineWF.IntegrateWaveForm(blstart, blstop);
      double wfintegrate =
          WF1->IntegrateWaveForm(traceSecondary, starttime, stoptime);
      double blcorr = (blintegrate - wfintegrate) / (stoptime - starttime);
      traceSecondary = WF1->GetTraces();
      // std::cout << shift << " : " << blstart << " : " << blstop << " : " <<
      // blcorr
      //           << std::endl;

      for (int iterWF = 0; iterWF < traceSecondary.size(); iterWF++) {
        val = iterWF - shift > 0 ? iterWF - shift
                                 : baselineWF.GetSize() + iterWF - shift;
        traceSecondary[iterWF] = traceSecondary[iterWF] - traceBL[val] + blcorr;
      }
      // WF1->Plot(traceSecondary, WF1->GetTracesSmooth());

      Channel = hitsVector[hititer]->GetChNum();
      Timestamp = hitsVector[hititer]->GetTimestamp();
      Board = hitsVector[hititer]->GetBoard();
      Energy = hitsVector[hititer]->GetEnergy();
      EnergyShort = hitsVector[hititer]->GetEnergyShort();
      Samples = new TArrayS(traceSecondary.size());
      // baselineVal = hitsVector[hititer]->GetWFPtr()->GetBaseLine();
      for (size_t i = 0; i < traceSecondary.size(); i++) {
        (*Samples)[i] = static_cast<Short_t>(std::round(0 - traceSecondary[i]));
      }
      Data_F->Fill();
      delete Samples;
      Samples = nullptr;
    }

    // fill the other channel as is
    else {
      Channel = hitsVector[hititer]->GetChNum();
      Timestamp = hitsVector[hititer]->GetTimestamp();
      Board = hitsVector[hititer]->GetBoard();
      Energy = hitsVector[hititer]->GetEnergy();
      EnergyShort = hitsVector[hititer]->GetEnergyShort();
      std::vector<double> traceother =
          hitsVector[hititer]->GetWFPtr()->GetTraces();
      Samples = new TArrayS(traceother.size());
      baselineVal = hitsVector[hititer]->GetWFPtr()->GetBaseLine();
      for (size_t i = 0; i < traceother.size(); i++) {
        (*Samples)[i] = static_cast<Short_t>(std::round(0 - traceother[i]));
      }
      Data_F->Fill();
      delete Samples;
      Samples = nullptr;
    }
    // digiAnalysis::WaveForm WF2(traceSecondary);
    // WF2.SetSmooth(16, "MovA");
    // WF1->SetSmooth(16, "MovA");
    // WF1->Plot(WF1->GetTraces(), WF2.GetTraces());
    // std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }
  Data_F->Write();
  fout->Close();

  // fApp->Run();
}