#ifndef Analysis_h
#define Analysis_h

#include "Events.h"
#include "PSBar.h"
#include "WaveForm.h"
#include "includes.hh"
#include "singleHits.h"
#include <TArray.h>
#include <TClass.h>
#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <TTreeIndex.h>
#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace digiAnalysis {

class singleHits;
class Event;
class PSBar;
class WaveForm;
// class diffWaveForm;

class Analysis {
private:
  std::vector<std::unique_ptr<Event>> vecOfEvents;
  std::vector<std::unique_ptr<singleHits>> vecOfHits;
  std::string fDatafileName;
  double EnergyThreshold;

public:
  Analysis();
  Analysis(std::string datafilename, ULong64_t numOfEvents = 0,
           double EThreshold = 0);
  Analysis(std::string datafilename, ULong64_t start, ULong64_t numOfEvents = 0,
           double EThreshold = 0);
  ~Analysis();

  void LoadData(ULong64_t numOfEvents, double EThreshold);
  void LoadData(ULong64_t start, ULong64_t numOfEvents, double EThreshold);
  void SetSingleHit(ULong64_t hitIndx, std::unique_ptr<singleHits> hit);
  void SortHits(const std::string &major, const std::string &minor = "");
  // void CreateEvents();

  void DeleteHit(ULong64_t hitNum);
  void ResizeHitsVector();

  // TTree *GetTreeFromFile(const std::string &filename);
  std::vector<std::unique_ptr<singleHits>> &GetSingleHitsVec();
  std::vector<singleHits> GetSingleHitsVec(ushort channel);
  std::unique_ptr<singleHits> GetSingleHit(ULong64_t hitIndx);
};

} // namespace digiAnalysis

#endif