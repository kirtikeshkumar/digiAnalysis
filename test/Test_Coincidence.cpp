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

int main(int argc, char *argv[])
{
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;
  std::string fname =
      "/media/kirtikesh/UbuntuFiles1/"
      "DataF_NaI_12_CoincidenceStudies_Na_HV_GainMatch_10min_2Vpp.root";

  digiAnalysis::Analysis an(fname, 0, 0000, 0);
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
  std::vector<int> chStart;       // stores the index of first entry from the channel
  for (int i = 0; i < nentries; i++)
  {
    if (hitsVector[i]->GetChNum() != currChannel)
    {
      currChannel = hitsVector[i]->GetChNum();
      numChannels += 1;
      channels.push_back(currChannel);
      chStart.push_back(i);
    }
  }

  std::cout << "There are " << numChannels << " channels (";
  for (int iter = 0; iter < numChannels; iter++)
  {
    std::cout << channels[iter] << ", ";
  }
  std::cout << ") in the data." << std::endl;

  // Estimate the number of events in a given time window around an event
  int timeWindow = 100000000; // in ps

  // within same detector
  //

  // between detectors
  int startIndex = -1;
  std::vector<int> currIndex(numChannels, 0);
  int startChannel = -1, stopChannel = -1;
  std::vector<Long64_t> currTime;
  Long64_t minTime = std::numeric_limits<Long64_t>::max();
  // initialize with the smallest start time
  for (int iter = 0; iter < numChannels; iter++)
  {
    currIndex[iter] = currIndex[iter] + chStart[iter];
    currTime[iter] = hitsVector[currIndex[iter]]->GetTimestamp();
    if (currTime[iter] < minTime)
    {
      minTime = hitsVector[currIndex[iter]]->GetTimestamp();
      startIndex = currIndex[iter];
      startChannel = iter;
    }
  }
  currIndex[startChannel] += 1;
  currTime[startChannel] = hitsVector[currIndex[startChannel]]->GetTimestamp();

  // now iterativley find the closest event
  bool isEnd = false;
  Long64_t startTime = hitsVector[startIndex]->GetTimestamp();
  Long64_t chTime = 0;
  std::vector<bool> isChEnd(numChannels, false);
  std::vector<Long64_t> chDelT(numChannels, 0);
  digiAnalysis::Pair coincPair;

  while (!isEnd)
  {
    chDelT.clear();
    // evaluate the time difference between current event and the next event in each channel
    for (int iterCh = 0; iterCh < numChannels; iterCh++)
    {
      !isChEnd[iterCh] ? chDelT.push_back(currTime[iterCh] - startTime) : chDelT.push_back(std::numeric_limits<Long64_t>::max());
    }
    // find the nearest event
    auto it = std::min_element(chDelT.begin(), chDelT.end());
    minTime = *it;
    stopChannel = std::distance(chDelT.begin(), it);

    // create pair if within the timeframe
    if (minTime < timeWindow)
    {
      coincPair.ClearPair();
      coincPair.SetPair(hitsVector[startIndex].get(), hitsVector[currIndex[stopChannel]].get());
      vecOfPairs.push_back(new digiAnalysis::Pair(coincPair));
    }

    // update the entries for next run
    startChannel = stopChannel;
    startIndex = currIndex[stopChannel];
    stopChannel = -1;
    currIndex[startChannel] += 1;
    currTime[startChannel] = hitsVector[currIndex[startChannel]]->GetTimestamp();
    startTime = hitsVector[startIndex]->GetTimestamp();

    // check if channel has finished
    bool endcheck = true;
    for (int iterCh = 0; iterCh < numChannels; iterCh++)
    {
      if ((iterCh + 1 < numChannels && currIndex[iterCh] == chStart[iterCh + 1]) || currIndex[iterCh] == nentries)
      {
        isChEnd[iterCh] = true;
      }
      endcheck = endcheck && isChEnd[iterCh];
    }
    isEnd = endcheck;
  }

  fApp->Run();
  return 0;
}