#ifndef Events_h
#define Events_h

#include "PSBar.h"
#include "WaveForm.h"
#include "includes.hh"
#include "singleHits.h"
#pragma once
#include <variant>

namespace digiAnalysis {
class singleHits;
class PSBar;
class WaveForm;

class Event {
private:
  std::variant<singleHits, PSBar> data;

public:
  Event();
  Event(const singleHits &hit);
  Event(const PSBar &bar);
  void Print();

  bool IsSingleHitEvt();
  bool IsPSBar();
  ULong64_t GetTStamp();
  UShort_t GetChNum();
  std::string GetEventType();
  UShort_t GetEnergy();
  UShort_t GetEnergyShort();
  float GetPSD();
#ifdef WAVES
  WaveForm GetWaveForm();
#endif
  ~Event();
};
} // namespace digiAnalysis
#endif