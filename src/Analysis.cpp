#include "Analysis.h"
#include "RtypesCore.h"
#include "includes.hh"
#include "singleHits.h"
#include <memory>

using namespace std;

namespace digiAnalysis {
Analysis::Analysis() {
  EnergyThreshold = 0;
  fDatafileName = "";
}
Analysis::Analysis(std::string datafilename, unsigned int numOfEvents,
                   double EThreshold) {
  if (datafilename.empty()) {
    throw std::logic_error("Data file name is empty!");
  }
  fDatafileName = datafilename;
  EnergyThreshold = EThreshold;
  LoadData(numOfEvents, EnergyThreshold);
}

void Analysis::LoadData(unsigned int numOfEvents, double EThreshold) {

  // Declaration of leaves types
  UShort_t Channel;
  ULong64_t Timestamp;
  UShort_t Board;
  UShort_t Energy;
  UShort_t EnergyShort;
#ifdef WAVES
  TArrayS *Samples;
#endif

  TFile *fp = new TFile(fDatafileName.c_str(), "READ");
  TTree *tr = (TTree *)fp->Get("Data_F");

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
#endif

  Long64_t nentries = tr->GetEntries();
  Long64_t nbytes = 0;

  if (numOfEvents > nentries) {
    std::cout << "Warning in LoadData: Cannot fetch " << numOfEvents
              << " entries from file containing " << nentries << " entries"
              << std::endl;
    numOfEvents = nentries;
  }
  if (numOfEvents == 0) {
    std::cout << "Warning in LoadData: 0 numOfEvents passed to LoadData, "
                 "fetching all entries"
              << std::endl;
    numOfEvents = nentries;
  }

  //////////////////////////////////////////////////////////////////////////////
  /// SORTING
  //////////////////////////////////////////////////////////////////////////////

  std::cout << "Now sorting: " << nentries << " entries in order of timestamp"
            << std::endl;

  tr->BuildIndex("Timestamp", "0");
  std::cout << "BuiltIndex" << std::endl;
  TTreeIndex *index = (TTreeIndex *)tr->GetTreeIndex();
  if (!index) {
    std::cerr << "Error creating tree index!" << std::endl;
    gROOT->ProcessLine(".q");
    return;
  }
  Long64_t *indices = index->GetIndex();
  if (!indices) {
    std::cerr << "Error retrieving index array!" << std::endl;
    gROOT->ProcessLine(".q");
    return;
  }

  for (Long64_t iev = 0; iev < numOfEvents; iev++) {
    nbytes += tr->GetEntry(indices[iev]);
    if (Energy > EThreshold) {
      std::unique_ptr<singleHits> hit = std::make_unique<singleHits>(
          iev, Channel, Board, Timestamp, Energy, EnergyShort);
#ifdef WAVES
      std::unique_ptr<singleHits> hit = std::make_unique<singleHits>(
          iev, Channel, Board, Timestamp, Energy, EnergyShort, &Samples);
#endif
      vecOfHits.push_back(std::move(hit));
    }
  }
}

std::vector<std::unique_ptr<singleHits>> &Analysis::GetSingleHits() {
  /*
    This returns a pointer of the vecOfHits in this class.
    Any modification to vecOfHits will reflect in the variable populated by this
    function
  */
  return vecOfHits;
}

void Analysis::SortHits(const std::string &major, const std::string &minor) {
  std::sort(vecOfHits.begin(), vecOfHits.end(),
            [major, minor](const std::unique_ptr<singleHits> &a,
                           const std::unique_ptr<singleHits> &b) {
              // Major sorting
              if (major == "Energy") {
                if (a->GetEnergy() != b->GetEnergy())
                  return a->GetEnergy() < b->GetEnergy();
              } else if (major == "Time") {
                if (a->GetTimestamp() != b->GetTimestamp())
                  return a->GetTimestamp() < b->GetTimestamp();
              } else if (major == "Channel") {
                if (a->GetChNum() != b->GetChNum())
                  return a->GetChNum() < b->GetChNum();
              } else if (major == "Board") {
                if (a->GetBoard() != b->GetBoard())
                  return a->GetBoard() < b->GetBoard();
              }

              // Minor sorting if specified
              if (!minor.empty()) {
                if (minor == "Energy") {
                  return a->GetEnergy() < b->GetEnergy();
                } else if (minor == "Time") {
                  return a->GetTimestamp() < b->GetTimestamp();
                } else if (minor == "Channel") {
                  return a->GetChNum() < b->GetChNum();
                } else if (minor == "Board") {
                  return a->GetBoard() < b->GetBoard();
                }
              }

              return false; // In case both major and minor criteria are equal
            });
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