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
    /**
     * @class Analysis
     * @brief Class for analyzing event data from files, applying energy thresholds, and managing hit and pair data.
     */
  private:
    std::vector<std::unique_ptr<singleHits>> vecOfHits;
    std::vector<Pair *> vecOfPairs;
    std::string fDatafileName;
    double EnergyThreshold;

  public:
    /**
     * @brief Default constructor for the Analysis class.
     */
    Analysis();

    /**
     * @brief Constructs an Analysis object for processing data from a specified file.
     *
     * @param datafilename The path to the data file to be analyzed.
     * @param numOfEvents The number of events to process from the data file. If set to 0, all events are processed. (Default: 0)
     * @param EThreshold The energy threshold to apply during analysis. Events below this threshold are ignored. (Default: 0)
     */
    Analysis(std::string datafilename, ULong64_t numOfEvents = 0, double EThreshold = 0);

    /**
     * @brief Constructs an Analysis object to process event data from a file.
     *
     * @param datafilename   The path to the data file to be analyzed.
     * @param start          The index of the first event to process.
     * @param numOfEvents    The number of events to process (default is 0, which means all events).
     * @param EThreshold     The energy threshold for event selection (default is 0).
     */
    Analysis(std::string datafilename, ULong64_t start, ULong64_t numOfEvents = 0, double EThreshold = 0);

    /**
     * @brief Destructor for the Analysis class. Cleans up resources.
     */
    ~Analysis();

    /**
     * @brief Loads data for analysis, applying an energy threshold.
     *
     * @param numOfEvents The number of events to load.
     * @param EThreshold The minimum energy threshold to apply when loading events.
     */
    void LoadData(ULong64_t numOfEvents, double EThreshold);

    /**
     * @brief Loads data for analysis from a specified starting event, applying an energy threshold.
     *
     * @param start The index of the first event to load.
     * @param numOfEvents The number of events to load.
     * @param EThreshold The minimum energy threshold to apply when loading events.
     */
    void LoadData(ULong64_t start, ULong64_t numOfEvents, double EThreshold);

    /**
     * @brief Sets a single hit at the specified index in the hit vector.
     *
     * @param hitIndx The index at which to set the hit.
     * @param hit The unique pointer to the singleHits object to set.
     */
    void SetSingleHit(ULong64_t hitIndx, std::unique_ptr<singleHits> hit);

    /**
     * @brief Sorts the vector of hits based on the specified major and optional minor criteria.
     *
     * @param major The primary sorting criterion.
     * @param minor The secondary sorting criterion (optional).
     */
    void SortHits(const std::string &major, const std::string &minor = "");

    /**
     * @brief Creates pairs of hits based on the specified channels.
     *
     * @param Ch1 The first channel number.
     * @param Ch2 The second channel number.
     */
    void CreatePairs(UShort_t Ch1, UShort_t Ch2);

    /**
     * @brief Deletes a hit from the hit vector at the specified index.
     *
     * @param hitNum The index of the hit to delete.
     */
    void DeleteHit(ULong64_t hitNum);

    /**
     * @brief Resizes the vector of hits to fit the current number of hits.
     */
    void ResizeHitsVector();

    /**
     * @brief Returns a reference to the vector of unique pointers to singleHits.
     *
     * @return Reference to the vector of unique pointers to singleHits.
     */
    std::vector<std::unique_ptr<singleHits>> &GetSingleHitsVec();

    /**
     * @brief Returns a vector of pointers to singleHits objects for a specific channel.
     *
     * @param channel The channel number to filter hits by.
     * @return Vector of pointers to singleHits objects for the specified channel.
     */
    std::vector<singleHits *> GetSingleHitsVec(ushort channel);

    /**
     * @brief Returns a unique pointer to a singleHits object at the specified index.
     *
     * @param hitIndx The index of the hit to retrieve.
     * @return Unique pointer to the singleHits object at the specified index.
     */
    std::unique_ptr<singleHits> GetSingleHit(ULong64_t hitIndx);

    /**
     * @brief Returns a vector of pointers to Pair objects.
     *
     * @return Vector of pointers to Pair objects.
     */
    std::vector<Pair *> GetPairsVec();
  };

} // namespace digiAnalysis

#endif
