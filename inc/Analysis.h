/*
the vecOfHits is the vector of unique_ptr to singleHits.
It is the only place where the singleHits objects are stored currently.
All others are only referencing this.
Never delete any object in this vector.
*/

#ifndef Analysis_h
#define Analysis_h

#include "PSBar.h"
#include "WaveForm.h"
#include "includes.hh"
#include "singleHits.h"
#include "globals.h"
#include "Pair.h"
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

namespace digiAnalysis
{

  class singleHits;
  class PSBar;
  class WaveForm;
  class Pair;
  // class diffWaveForm;

  class Analysis
  {
  private:
    std::vector<std::unique_ptr<singleHits>> vecOfHits;
    std::vector<Pair *> vecOfPairs;
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

    void CreatePairs(UShort_t Ch1, UShort_t Ch2);
    void DeleteHit(ULong64_t hitNum);
    void ResizeHitsVector();

    // TTree *GetTreeFromFile(const std::string &filename);
    std::vector<std::unique_ptr<singleHits>> &GetSingleHitsVec();
    std::vector<singleHits *> GetSingleHitsVec(ushort channel);
    std::unique_ptr<singleHits> GetSingleHit(ULong64_t hitIndx);

    std::vector<Pair *> GetPairsVec();
  };

} // namespace digiAnalysis

#endif