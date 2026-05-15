/*
** Code to create coincident event pairs between any number of channels
** Steps:
** 1. Segregate the events in increasing order of channel + time
** 2. Save the start and stop points of every channel
** 3. Start with the event that is first in time.
** 4. Find the nearest event among all channels including self
** 5. Save into pair if events fall in different channels
** 6. Change the start event to the nearest event and redo from Step 3
*/

#include "Analysis.h"
#include "WaveForm.h"
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
#include <iostream>
#include <ratio>
#include <vector>

int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;

  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //     "NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp/FILTERED/"
  //     "SDataF_NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp.root";

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
      "NaI13_15May26_1900_1345_NoSrc_Coinc144ns_WAVES/FILTERED/"
      "SDataF_NaI13_15May26_1900_1345_NoSrc_Coinc144ns_WAVES.root";

  digiAnalysis::Analysis an(fname, 0100, 000, 0);
  std::cout << "getting the vector from an" << std::endl;

  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int nentries = hitsVector.size();
  std::cout << "got the vector from an: " << nentries << std::endl;
  std::vector<digiAnalysis::Pair *> vecOfPairs;

  an.SortHits("Channel", "Time");

  // locate and segregate the channels
  UShort_t numChannels = 0;
  UShort_t currChannel = -1;
  std::vector<UShort_t> channels; // stores the value of channel number
  std::vector<int> chStart; // stores the index of first entry from the channel
  for (int i = 0; i < nentries; i++) {
    if (hitsVector[i]->GetChNum() != currChannel) {
      currChannel = hitsVector[i]->GetChNum();
      numChannels += 1;
      channels.push_back(currChannel);
      chStart.push_back(i);
    }
  }

  std::cout << "There are " << numChannels << " channels (";
  for (int iter = 0; iter < numChannels; iter++) {
    std::cout << channels[iter] << ", ";
  }
  std::cout << ") in the data." << std::endl;

  for (int iter = 0; iter < numChannels; iter++) {
    std::cout << "Channel " << channels[iter] << " starts at: " << chStart[iter]
              << std::endl;
  }

  // Estimate the number of events in a given time window around an event
  int timeWindow = 1000000; // in ps
  TH1 *hDelT = new TH1F("hDelT", "hDelT", 10000, 0, timeWindow / 1E6);
  TH2 *hE1E2 = new TH2I("hE1E2", "hE1E2", 4000, 0, 4000, 4000, 0, 4000);
  TH1 *hETot = new TH1F("hETot", "hETot", 8000, 0, 8000);
  TH2 *hE1ETot = new TH2I("hE1ETot", "hE1ETot", 4000, 0, 4000, 8000, 0, 8000);
  TH2 *hE1PSD = new TH2I("hE1PSD", "hE1PSD", 4000, 0, 4000, 1000, -1, 1);
  TH2 *hE1MT = new TH2I("hE1MT", "hE1MT", 4000, 0, 4000, 10000, -5, 5);
  TH2 *hPSDMT = new TH2I("hPSDMT", "hPSDMT", 1000, -10, 10, 1000, -5, 5);

  // within same detector
  //

  // between detectors
  int startIndex = -1;
  std::vector<int> currIndex(numChannels, 0);
  int startChannel = -1, stopChannel = -1;
  std::vector<ULong64_t> currTime;
  ULong64_t minTime = std::numeric_limits<ULong64_t>::max();

  // initialize with the smallest start time
  for (int iter = 0; iter < numChannels; iter++) {
    currIndex[iter] = currIndex[iter] + chStart[iter];
    currTime.push_back(hitsVector[currIndex[iter]]->GetTimestamp());
    if (currTime[iter] < minTime) {
      minTime = hitsVector[currIndex[iter]]->GetTimestamp();
      startIndex = currIndex[iter];
      startChannel = iter;
    }
  }
  ULong64_t startTime = hitsVector[startIndex]->GetTimestamp();
  currIndex[startChannel] += 1;
  currTime[startChannel] = hitsVector[currIndex[startChannel]]->GetTimestamp();

  // now iterativley find the closest event
  bool isEnd = false;
  ULong64_t chTime = 0;
  std::vector<bool> isChEnd(numChannels, false);
  std::vector<ULong64_t> chDelT(numChannels, 0);
  digiAnalysis::Pair coincPair;

  while (!isEnd) {
    chDelT.clear();
    // evaluate the time difference between current event and the next event in
    // each channel
    for (int iterCh = 0; iterCh < numChannels; iterCh++) {
      !isChEnd[iterCh]
          ? chDelT.push_back(currTime[iterCh] - startTime)
          : chDelT.push_back(std::numeric_limits<ULong64_t>::max());
    }
    // find the nearest event
    auto it = std::min_element(chDelT.begin(), chDelT.end());
    minTime = *it;
    stopChannel = std::distance(chDelT.begin(), it);

    // create pair if within the timeframe
    if ((minTime < timeWindow) and (startChannel != stopChannel)) {
      coincPair.ClearPair();
      coincPair.SetPair(*hitsVector[startIndex].get(),
                        *hitsVector[currIndex[stopChannel]].get());
      vecOfPairs.push_back(new digiAnalysis::Pair(coincPair));
    }
    // update the entries for next run
    startChannel = stopChannel;
    startIndex = currIndex[stopChannel];
    stopChannel = -1;
    currIndex[startChannel] += 1;
    currIndex[startChannel] < nentries
        ? currTime[startChannel] =
              hitsVector[currIndex[startChannel]]->GetTimestamp()
        : currTime[startChannel] = std::numeric_limits<ULong64_t>::max();
    startTime = hitsVector[startIndex]->GetTimestamp();

    // check if channel has finished
    bool endcheck = true;
    for (int iterCh = 0; iterCh < numChannels; iterCh++) {
      if ((iterCh + 1 < numChannels &&
           currIndex[iterCh] == chStart[iterCh + 1]) ||
          currIndex[iterCh] == nentries) {
        isChEnd[iterCh] = true;
      }
      endcheck = endcheck && isChEnd[iterCh];
    }
    isEnd = endcheck;
  }

  int nPairs = vecOfPairs.size();
  std::cout << nPairs << " Pairs were formed in the data." << std::endl;
  double Energy1 = 0;
  double Energy2 = 0;
  double PSD = 0, MT = 0;
  for (int iter = 0; iter < nPairs; iter++) {
    Energy1 = vecOfPairs[iter]->GetPairHitEnergy(0) * 1.0174 -
              38.84; // Gain match case
    PSD = vecOfPairs[iter]->GetHitPtr(0)->GetPSD();
#ifdef WAVES
    MT = vecOfPairs[iter]->GetHitPtr(0)->GetMeanTime();
#endif
    vecOfPairs[iter]->GetPairHitEnergy(0) > 694
        ? Energy1 = vecOfPairs[iter]->GetPairHitEnergy(0) * 0.09465 - 5.7613
        : Energy1 = vecOfPairs[iter]->GetPairHitEnergy(0) * 0.08696 - 0.4222;
    // 1900V
    Energy2 = vecOfPairs[iter]->GetPairHitEnergy(1) * 0.98145 -
              14.6; //* 1.0973 - 58.91;//
    double ETot = Energy1 + Energy2;
    // if (ETot > 600 and ETot < 800)
    hDelT->Fill(vecOfPairs[iter]->GetPairDelTime() / 1E6);
    hE1E2->Fill(Energy1, Energy2);
    hETot->Fill(ETot);
    hE1ETot->Fill(Energy1, ETot);
    hE1PSD->Fill(Energy1, PSD);
#ifdef WAVES
    hE1MT->Fill(Energy1, MT);
    hPSDMT->Fill(PSD, MT);
#endif
  }
  TCanvas *c1 = new TCanvas("c1", "timeDiff", 800, 600);
  hDelT->Draw("HIST");
  TCanvas *c2 = new TCanvas("c2", "E1E2", 800, 600);
  hE1E2->Draw("COLZ");
  TCanvas *c3 = new TCanvas("c3", "ETotal", 800, 600);
  hETot->Draw("HIST");
  TCanvas *c4 = new TCanvas("c4", "E1ETotal", 800, 600);
  hE1ETot->Draw("COLZ");
  TCanvas *c5 = new TCanvas("c5", "E1PSD", 800, 600);
  hE1PSD->Draw("COLZ");
#ifdef WAVES
  TCanvas *c6 = new TCanvas("c6", "E1MT", 800, 600);
  hE1MT->Draw("COLZ");
  TCanvas *c7 = new TCanvas("c7", "PSD MT", 800, 600);
  hPSDMT->Draw("COLZ");
#endif

#ifdef WAVES
  bool keepGoing = true;
  std::string userInput;
  UShort_t wfSz;
  std::vector<digiAnalysis::WaveForm> waveformVector;
  for (int iter = 0; iter < nPairs && keepGoing; iter++) {
    vecOfPairs[iter]->GetPairHitEnergy(0) > 694
        ? Energy1 = vecOfPairs[iter]->GetPairHitEnergy(0) * 0.09465 - 5.7613
        : Energy1 =
              vecOfPairs[iter]->GetPairHitEnergy(0) * 0.08696 - 0.4222; // 1900V
    Energy2 = vecOfPairs[iter]->GetPairHitEnergy(1) * 0.98145 -
              14.6; //* 1.0973 - 58.91;
    MT = vecOfPairs[iter]->GetHitPtr(0)->GetMeanTime();
    if (MT < 2.85 and MT > 2.0 and Energy1 < 390 and
        Energy1 > 370) { // and // MT > 3.05 and MT < 3.06 and
      // (Energy2 + Energy1 > 580) and (Energy2 + Energy1 < 740)) {
      vecOfPairs[iter]->GetHitPtr(0)->GetWFPtr()->SetTracesFFT();
      vecOfPairs[iter]->GetHitPtr(0)->GetWFPtr()->Plot();
      wfSz = vecOfPairs[iter]->GetHitPtr(0)->GetWFPtr()->GetSize();
      waveformVector.push_back(*vecOfPairs[iter]->GetHitPtr(0)->GetWFPtr());
      std::cout << iter << " Energy of D1 is approx: " << Energy1 << " keV"
                << std::endl;
      std::cout << "Energy Total is approx: " << Energy1 + Energy2 << " keV"
                << std::endl;
      std::cout << "Do you want to see the next waveform? (y/n): ";
      std::getline(std::cin, userInput);
      if (userInput != "y" && userInput != "Y") {
        keepGoing = false;
      }
    }
  }

  digiAnalysis::WaveForm WFAveraged(wfSz, waveformVector);
  WFAveraged.SetSmooth(30);
  WFAveraged.SetTracesFFT("smooth");
  WFAveraged.Plot();
#endif
  fApp->Run();
  return 0;
}