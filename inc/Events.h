#ifndef Events_h
#define Events_h

#include "PSBar.h"
#include "singleHits.h"
#pragma once
#include <variant>

namespace digiAnalysis {
class singleHits;
class PSBar;

class Event {
private:
  std::variant<singleHits, PSBars> data;

public:
  Event();
  Event(const singleHits &hit);
  Event(const PSBar &bar);
  void Print() const;

  ULong64_t GetTStamp();
  UShort_t GetChNum();
  std::string GetEventType();
  UShort_t GetEnergy();
  Ushort_t GetEnergyShort();
  UShort_t GetWaveForm();

  ~Event();
}
} // namespace digiAnalysis
#endif