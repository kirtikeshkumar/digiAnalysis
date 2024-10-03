#ifndef Analysis_h
#define Analysis_h

#include "Events.h"
#include "Hits.h"
#include "WaveForm.h"
#include "includes.hh"
#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace digiAnalysis {

class singleHits;
class WaveForm;
class diffWaveForm;

class Analysis {
private:
  std::vector<Event *> vecOfEvents;
  std::vector<singleHits *> vecOfHits;
  std::string fDatafileName;
  unsigned int EnergyThreshold;

public:
  Analysis();
  Analysis(std::string datafilename, unsigned int numOfEvents = 0,
           double EThreshold = 0);
  ~Analysis(unsigned int numOfEvents, unsigned int EnergyThreshold);

  void LoadData(); // read hits/events from file based on cmake variable
  void CreateEvents();
};

} // namespace digiAnalysis

#endif