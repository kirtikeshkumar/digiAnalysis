#ifndef PSBar
#define PSBar

#include "singleHits.h"
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <TROOT.h>

class TH1F;

namespace digiAnalysis {
class singleHits;

class PSBar {
private:
  ushort BarIndex;    // Index of Bar starting from 0
  UInt_t Qlong;       // Folded Charge first 16 bits near, next 16 bits farftr
  ULong64_t TimeNear; // DAQ timestamp of smaller channel
  ULong64_t Timestamp;
  Int_t Delt; // Time diff between left and right PMT

  // Storing the event no. corresponding to TTree
  ULong64_t EvNo;

  UInt_t Qnear;
  UInt_t Qfar;

#ifdef WAVES
  std::unique_ptr<WaveForm> nearWF;
  std::unique_ptr<smoothWaveForm> nearSmoothWF;
  std::unique_ptr<WaveForm> farWF;
  std::unique_ptr<smoothWaveForm> farSmoothWF;
#endif

public:
  PSBar();
  PSBar(unsigned int bIndex);
  // Required copy constructor
  PSBar(const PSBar &sbar);
  PSBar(singleHits *h1, singleHits *h2);

  // Getters
  ushort GetBarIndex() const;
  ULong64_t GetEvNo() const { return fEvNo; }
  UInt_t GetQNear();
  UInt_t GetQFar();
  Double_t GetQMean();
  Long_t GetDelT() const;
  Long_t GetTStampNear();
  Long_t GetTStampFar();
  Long_t GetTStampAverage();
  double GetHitPos();

  // Printer
  void Print();

  ~PSBar();
};
} // namespace digiAnalysis
