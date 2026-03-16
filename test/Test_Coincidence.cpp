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
  //     "NaI_12_CoincidenceStudies_Na_HV_GainMatch_10min_2Vpp/FILTERED/"
  //     "DataF_NaI_12_CoincidenceStudies_Na_HV_GainMatch_10min_2Vpp.root";
  std::string fname =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
      "NaI_12_CoincidenceStudies_Na_HV_1900V_1365V_96nsCoinc_60min_2Vpp_WAVES/"
      "FILTERED/"
      "SDataF_NaI_12_CoincidenceStudies_Na_HV_1900V_1365V_96nsCoinc_60min_2Vpp_"
      "WAVES."
      "root";

  digiAnalysis::Analysis an(fname, 0, 0, 0);
  std::cout << "getting the vector from an" << std::endl;

  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int nentries = hitsVector.size();
  std::cout << "got the vector from an" << nentries << std::endl;
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
  TH2 *hE1E2 = new TH2I("hE1E2", "hE1E2", 8192, 0, 8192, 8192, 0, 8192);
  TH1 *hETot = new TH1F("hETot", "hETot", 16384, 0, 2000);

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
    // std::cout << "startCH: " << startChannel << " Index: " << startIndex
    //           << " Time: " << startTime << std::endl;
    for (int iterCh = 0; iterCh < numChannels; iterCh++) {
      !isChEnd[iterCh]
          ? chDelT.push_back(currTime[iterCh] - startTime)
          : chDelT.push_back(std::numeric_limits<ULong64_t>::max());
      // std::cout << "iterch: " << iterCh << " Index: " << currIndex[iterCh]
      //           << " Time: " << currTime[iterCh]
      //           << " chDelT: " << chDelT[iterCh] << std::endl;
    }
    // find the nearest event
    auto it = std::min_element(chDelT.begin(), chDelT.end());
    minTime = *it;
    stopChannel = std::distance(chDelT.begin(), it);
    // std::cout << "minTime of: " << minTime << " in stopChannel: " <<
    // stopChannel
    //           << std::endl;

    // create pair if within the timeframe
    if ((minTime < timeWindow) and (startChannel != stopChannel)) {
      coincPair.ClearPair();
      coincPair.SetPair(hitsVector[startIndex].get(),
                        hitsVector[currIndex[stopChannel]].get());
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
      // std::cout << " Ch " << iterCh << " End: " << isChEnd[iterCh];
      endcheck = endcheck && isChEnd[iterCh];
    }
    // std::cout << std::endl;
    isEnd = endcheck;
    // std::cout << "End: " << isEnd << std::endl;
  }

  int nPairs = vecOfPairs.size();
  std::cout << nPairs << " Pairs were formed in the data." << std::endl;
  double Energy1 = 0;
  double Energy2 = 0;
  for (int iter = 0; iter < nPairs; iter++) {
    // Energy1 = vecOfPairs[iter]->GetPairHitEnergy(0) * 1.0174 -
    //           38.84; // Gain match case
    Energy1 =
        vecOfPairs[iter]->GetPairHitEnergy(0) * 0.0907585 - 0.4613; // 1900V
    Energy2 = vecOfPairs[iter]->GetPairHitEnergy(1) * 1.0973 - 58.91;
    double ETot = Energy1 + Energy2;
    // if (ETot > 1700 and ETot < 2000 and Energy2 > 1200 and Energy2 < 1400)
    hDelT->Fill(vecOfPairs[iter]->GetPairDelTime() / 1E6);
    hE1E2->Fill(Energy1, Energy2);
    hETot->Fill(ETot);
  }
  TCanvas *c1 = new TCanvas("c1", "timeDiff", 800, 600);
  hDelT->Draw("HIST");
  TCanvas *c2 = new TCanvas("c2", "E1E2", 800, 600);
  hE1E2->Draw("COLZ");
  TCanvas *c3 = new TCanvas("c3", "ETotal", 800, 600);
  hETot->Draw("HIST");
  fApp->Run();
  return 0;
}