#include "Analysis.h"
#include "RtypesCore.h"
#include "includes.hh"
#include "singleHits.h"
#include <iostream>
#include <memory>

using namespace std;

namespace digiAnalysis
{
  Analysis::Analysis()
  {
    EnergyThreshold = 0;
    fDatafileName = "";
  }

  Analysis::Analysis(std::string datafilename, ULong64_t numOfEvents,
                     double EThreshold)
  {
    if (datafilename.empty())
    {
      throw std::logic_error("Data file name is empty!");
    }
    fDatafileName = datafilename;
    EnergyThreshold = EThreshold;
    LoadData(numOfEvents, EnergyThreshold);
  }

  Analysis::Analysis(std::string datafilename, ULong64_t start,
                     ULong64_t numOfEvents, double EThreshold)
  {
    if (datafilename.empty())
    {
      throw std::logic_error("Data file name is empty!");
    }
    fDatafileName = datafilename;
    EnergyThreshold = EThreshold;
    LoadData(start, numOfEvents, EnergyThreshold);
  }

  void Analysis::LoadData(ULong64_t numOfEvents, double EThreshold)
  {
    LoadData(0, numOfEvents, EThreshold);
  }

  void Analysis::LoadData(ULong64_t start, ULong64_t numOfEvents,
                          double EThreshold)
  {

    TFile *fp = new TFile(fDatafileName.c_str(), "READ");

    if (!fp || fp->IsZombie())
    {
      std::cerr << "Error: Unable to open file " << fDatafileName << std::endl;
      return; // Exit or handle the error appropriately
    }

    TTree *tr = (TTree *)fp->Get("Data_F");

    // Check if the TTree was retrieved successfully
    if (!tr)
    {
      std::cerr << "Error: Unable to retrieve TTree 'Data_F' from file."
                << std::endl;
      fp->Close();
      return; // Exit or handle the error appropriately
    }
    // Declaration of leaves types
    UShort_t Channel;
    ULong64_t Timestamp;
    UShort_t Board;
    UShort_t Energy;
    UShort_t EnergyShort;
    // #ifdef WAVES
    TArrayS *Samples = nullptr;
    // #endif

    std::cout << "Loading data from: " << fDatafileName << std::endl;
    // TTree *tr = GetTreeFromFile(fDatafileName);

    // Set branch addresses.
    tr->SetBranchAddress("Channel", &Channel);
    tr->SetBranchAddress("Timestamp", &Timestamp);
    tr->SetBranchAddress("Board", &Board);
    tr->SetBranchAddress("Energy", &Energy);
    tr->SetBranchAddress("EnergyShort", &EnergyShort);
#ifdef WAVES
    tr->SetBranchAddress("Samples", &Samples);
    std::cout << "Branch waves set" << std::endl;
#endif

    Long64_t nentries = tr->GetEntries();
    Long64_t nbytes = 0;

    if (numOfEvents + start > nentries)
    {
      std::cout << "Warning in LoadData: Cannot fetch till" << start + numOfEvents
                << " entries from file containing " << nentries << " entries"
                << std::endl;
      numOfEvents = nentries - start;
    }

    if (start > nentries)
    {
      std::cout << "Warning in LoadData: Cannot fetch " << start
                << " th entry from file containing " << nentries << " entries"
                << std::endl;
      start = 0;
      numOfEvents = nentries;
    }

    if (numOfEvents == 0)
    {
      std::cout << "Warning in LoadData: 0 numOfEvents passed to LoadData, "
                   "fetching all entries after "
                << start << std::endl;
      numOfEvents = nentries - start;
    }

    //////////////////////////////////////////////////////////////////////////////
    /// SORTING
    //////////////////////////////////////////////////////////////////////////////

    std::cout << "Now sorting: " << nentries << " entries in order of timestamp"
              << std::endl;

    tr->BuildIndex("Timestamp", "0");
    TTreeIndex *index = (TTreeIndex *)tr->GetTreeIndex();
    if (!index)
    {
      std::cerr << "Error creating tree index!" << std::endl;
      gROOT->ProcessLine(".q");
      return;
    }
    Long64_t *indices = index->GetIndex();
    if (!indices)
    {
      std::cerr << "Error retrieving index array!" << std::endl;
      gROOT->ProcessLine(".q");
      return;
    }

    std::cout << "sorting done" << std::endl;
    //////////////////////////////////////////////////////////////////////////////
    /// READING
    //////////////////////////////////////////////////////////////////////////////
    for (ULong64_t iev = start; iev < start + numOfEvents; iev++)
    {
      if ((iev - start) % 500000 == 0)
      {
        std::cout << "Reading: " << indices[iev] << std::endl;
      }
      nbytes += tr->GetEntry(indices[iev]);
      if (Energy > EThreshold)
      {
#ifndef WAVES
        std::unique_ptr<singleHits> hit = std::make_unique<singleHits>(
            iev, Channel, Board, Timestamp, Energy, EnergyShort);
#else
        std::unique_ptr<singleHits> hit = std::make_unique<singleHits>(
            iev, Channel, Board, Timestamp, Energy, EnergyShort, Samples);
#endif
        vecOfHits.push_back(std::move(hit));
      }
    }
  }

  void Analysis::SetSingleHit(ULong64_t hitIndx,
                              std::unique_ptr<singleHits> hit)
  {
    if (hitIndx < vecOfHits.size())
    {
      if (vecOfHits[hitIndx] == nullptr)
      {
        vecOfHits[hitIndx] = std::move(hit);
      }
      else
      {
        std::cout << "Warning in Analysis::SetSingleHit. hitIndx not Empty"
                  << std::endl;
        vecOfHits[hitIndx] = std::move(hit);
      }
    }
    else if (hitIndx == vecOfHits.size())
    {
      vecOfHits.push_back(std::move(hit));
    }
    else
    {
      std::cout << "Error in Analysis::SetSingleHit. hitIndx "
                << hitIndx - vecOfHits.size() << " out of bounds" << std::endl;
      std::cout << "hit added at end of vector";
      vecOfHits.push_back(std::move(hit));
    }
  }

  void Analysis::DeleteHit(ULong64_t hitNum)
  {
    if (hitNum < vecOfHits.size())
    {
      vecOfHits.push_back(nullptr);
      std::swap(vecOfHits[hitNum], vecOfHits.back());
      vecOfHits.pop_back();
    }
  }

  void Analysis::ResizeHitsVector()
  {
    vecOfHits.erase(std::remove_if(vecOfHits.begin(), vecOfHits.end(),
                                   [](const std::unique_ptr<singleHits> &ptr)
                                   {
                                     return ptr == nullptr;
                                   }),
                    vecOfHits.end());
  }

  std::vector<std::unique_ptr<singleHits>> &Analysis::GetSingleHitsVec()
  {
    /*
      This returns a reference to vecOfHits in this class.
      Thus vecOfHits and the variable created by using this func are
      "NOT" independent
    */
    return vecOfHits;
  }

  std::vector<singleHits *> Analysis::GetSingleHitsVec(ushort channel)
  {
    std::vector<singleHits *> vecOfHitsChannel;
    vecOfHitsChannel.reserve(vecOfHits.size());

    for (const auto &hit : vecOfHits)
    {
      if (hit && hit->GetChNum() == channel)
      {
        vecOfHitsChannel.push_back(hit.get()); // Get the raw pointer
      }
    }
    return vecOfHitsChannel; // Return by value
  }

  std::unique_ptr<singleHits> Analysis::GetSingleHit(ULong64_t hitIndx)
  {
    /*Here the ownership of singleHit transferred*/
    if (hitIndx < vecOfHits.size())
    {
      return std::move(vecOfHits[hitIndx]);
    }
    else
    {
      std::cout << "err: Analysis::GetSingleHit(ULong64_t hitIndx): hitIndx out "
                   "of range"
                << std::endl;
      return nullptr;
    }
  }

  void Analysis::SortHits(const std::string &major, const std::string &minor)
  {
    std::sort(vecOfHits.begin(), vecOfHits.end(),
              [major, minor](const std::unique_ptr<singleHits> &a,
                             const std::unique_ptr<singleHits> &b)
              {
                // Major sorting
                if (major == "Energy")
                {
                  if (a->GetEnergy() != b->GetEnergy())
                    return a->GetEnergy() < b->GetEnergy();
                }
                else if (major == "Time")
                {
                  if (a->GetTimestamp() != b->GetTimestamp())
                    return a->GetTimestamp() < b->GetTimestamp();
                }
                else if (major == "Channel")
                {
                  if (a->GetChNum() != b->GetChNum())
                    return a->GetChNum() < b->GetChNum();
                }
                else if (major == "Board")
                {
                  if (a->GetBoard() != b->GetBoard())
                    return a->GetBoard() < b->GetBoard();
                }
                else if (major == "PSD")
                {
                  if (a->GetPSD() != b->GetPSD())
                    return a->GetPSD() < b->GetPSD();
                }
#ifdef WAVES
                else if (major == "EvalEnergy")
                {
                  if (a->GetEvalEnergy() != b->GetEvalEnergy())
                    return a->GetEvalEnergy() < b->GetEvalEnergy();
                }
#endif

                // Minor sorting if specified
                if (!minor.empty())
                {
                  if (minor == "Energy")
                  {
                    return a->GetEnergy() < b->GetEnergy();
                  }
                  else if (minor == "Time")
                  {
                    return a->GetTimestamp() < b->GetTimestamp();
                  }
                  else if (minor == "Channel")
                  {
                    return a->GetChNum() < b->GetChNum();
                  }
                  else if (minor == "Board")
                  {
                    return a->GetBoard() < b->GetBoard();
                  }
                  else if (minor == "PSD")
                  {
                    if (a->GetPSD() != b->GetPSD())
                      return a->GetPSD() < b->GetPSD();
                  }

#ifdef WAVES
                  else if (major == "EvalEnergy")
                  {
                    if (a->GetEvalEnergy() != b->GetEvalEnergy())
                      return a->GetEvalEnergy() < b->GetEvalEnergy();
                  }
                }
#endif

                return false; // In case both major and minor criteria are equal
              });
  }

  void Analysis::CreatePairs(UShort_t Ch1, UShort_t Ch2)
  {

    /*
    Make Pairs of singleHits from Ch1 and Ch2 within PairCoincWindow
    iterate over hits in Ch1 and for each hit find the closest hit in Ch2
    */
    if (Ch1 == Ch2)
    {
      std::cerr << "Error in Analysis::CreatePairs: Ch1 and Ch2 cannot be same"
                << std::endl;
      return;
    }

    std::vector<singleHits *> hitsVectorCh0 = GetSingleHitsVec(Ch1);
    std::vector<singleHits *> hitsVectorCh1 = GetSingleHitsVec(Ch2);

    if (hitsVectorCh0.empty() || hitsVectorCh1.empty())
    {
      std::cerr << "Error in Analysis::CreatePairs: One or both channels have "
                   "no hits"
                << std::endl;
      return;
    }

    Long64_t i = 0, j = 0;
    Long64_t Ch0Time = hitsVectorCh0[0]->GetTimestamp();
    Long64_t Ch1Time = hitsVectorCh1[0]->GetTimestamp();
    Long64_t TDiffCurr = Ch0Time - Ch1Time;
    Long64_t TDiffMin = 0;
    Pair currPair;

    while (i < hitsVectorCh0.size())
    {
      Ch0Time = hitsVectorCh0[i]->GetTimestamp();
      TDiffMin = 100000000000000;
      while (j < hitsVectorCh1.size())
      {
        Ch1Time = hitsVectorCh1[j]->GetTimestamp();
        TDiffCurr = Ch1Time - Ch0Time;
        if (fabs(TDiffCurr) < fabs(TDiffMin))
        {
          TDiffMin = TDiffCurr;
          j = j + 1;
        }
        else
        {
          if (fabs(TDiffMin / 1000) < PairCoincWindow) // since the TimeStamp is in ps
          {
            currPair.ClearPair();
            currPair.SetPair(hitsVectorCh0[i], hitsVectorCh1[j - 1]);
            vecOfPairs.push_back(&currPair);
          }
          i = i + 1;
          j = (j >= 2) ? j - 2 : j;
          break;
        }
      }
      if (i >= hitsVectorCh0.size() - 1 or j >= hitsVectorCh1.size() - 1)
      {
        break;
      }
    }
  }

  std::vector<Pair *> Analysis::GetPairsVec()
  {
    return vecOfPairs;
  }

  /*
  TTree *Analysis::GetTreeFromFile(const std::string &filename) {
    TFile *fp = new TFile(filename.c_str(), "READ");
    if (!fp || fp->IsZombie()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return nullptr;
    }

    TTree *tree = nullptr;
    // Loop over all objects in the file
    TIter next(fp->GetListOfKeys());
    TKey *key;

    // Find the first object that is a TTree
    while ((key = (TKey *)next())) {
      TObject *obj = key->ReadObj();
      if (obj->IsA()->InheritsFrom(TTree::Class())) {
        tree = (TTree *)obj;
        std::cout << "Found tree: " << obj->GetName() << std::endl;
        break; // Exit after finding the first tree
      }
    }

    // If no tree was found
    if (!tree) {
      std::cerr << "No TTree found in the file: " << filename << std::endl;
    }

    return tree;
  }*/

  Analysis::~Analysis() {}
} // namespace digiAnalysis