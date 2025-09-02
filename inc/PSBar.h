#ifndef PSBar_h
#define PSBar_h

#include "WaveForm.h"
#include "singleHits.h"
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <TROOT.h>

class TH1F;

namespace digiAnalysis {
class singleHits;
class WaveForm;

class PSBar {
private:
  ushort BarIndex;    // Index of Bar starting from 0
  ULong64_t TimeNear; // DAQ timestamp of smaller channel
  ULong64_t Timestamp;
  Int_t Delt; // Time diff between left and right PMT
  UInt_t QNear;
  UInt_t QFar;
  UInt_t QNearShort;
  UInt_t QFarShort;
  double QMean;
  double QMeanShort;
  double PSD;

#ifdef WAVES
  std::unique_ptr<WaveForm> WF; // combined Waveform containing first half as
                                // near and second half as far
#endif

public:
  PSBar();
  PSBar(unsigned int bIndex);
  // Required copy constructor
  PSBar(const PSBar &sbar);
  PSBar(singleHits *h1, singleHits *h2);

  // Getters
  ushort GetBarIndex() const { return BarIndex; }
  UInt_t GetQNear();
  UInt_t GetQFar();
  double GetQMean();
  double GetQMeanShort();
  double GetPSD();
  Long_t GetDelT() const;
  Long_t GetTStampNear();
  Long_t GetTStampFar();
  Long_t GetTStampAverage();
  double GetHitPos();

  // Printer
  void Print();

#ifdef WAVES
  std::unique_ptr<WaveForm> GetWF(); // combined Waveform containing first half
                                     // as near and second half as far
  WaveForm *GetWFPtr();
  void SetWF(const WaveForm &wf);
#endif

  ~PSBar();
};
} // namespace digiAnalysis
#endif
