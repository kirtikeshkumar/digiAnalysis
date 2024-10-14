#include "Events.h"
#include "PSBar.h"
#include "includes.hh"
#include "singleHits.h"

using namespace std;

namespace digiAnalysis {
/*Constructors*/
Event::Event() //  default constructor
{
  data = singleHits();
}

Event::Event(const singleHits &hit) // Copy from singleHit
{
  data = hit;
}

Event::Event(const PSBar &bar) // Copy from singleHit
{
  data = bar;
}

void Event::Print() {
  if (IsSingleHitEvt()) {
    std::get<singleHits>(data).Print();
  } else if (IsPSBar()) {
    std::get<PSBar>(data).Print();
  }
}

bool Event::IsSingleHitEvt() {
  return std::holds_alternative<singleHits>(data);
}
bool Event::IsPSBar() { return std::holds_alternative<PSBar>(data); }

ULong64_t Event::GetTStamp() {
  if (IsSingleHitEvt()) {
    return std::get<singleHits>(data).GetTimestamp();
  } else if (IsPSBar()) {
    return std::get<PSBar>(data).GetTStampAverage();
  }

  throw std::runtime_error("Event does not contain a valid data.");
}
UShort_t Event::GetChNum() {
  if (IsSingleHitEvt()) {
    return std::get<singleHits>(data).GetChNum();
  } else if (IsPSBar()) {
    return std::get<PSBar>(data).GetBarIndex();
  }

  throw std::runtime_error("Event does not contain a valid data.");
}
std::string Event::GetEventType() {
  if (IsSingleHitEvt()) {
    return "singleHit";
  } else if (IsPSBar()) {
    return "PSBar";
  }

  throw std::runtime_error("Event does not contain a valid data.");
}
UShort_t Event::GetEnergy() {
  if (IsSingleHitEvt()) {
    return std::get<singleHits>(data).GetEnergy();
  } else if (IsPSBar()) {
    return std::get<PSBar>(data).GetQMean();
  }

  throw std::runtime_error("Event does not contain a valid data.");
}
UShort_t Event::GetEnergyShort() {
  if (IsSingleHitEvt()) {
    return std::get<singleHits>(data).GetEnergyShort();
  } else if (IsPSBar()) {
    return std::get<PSBar>(data).GetQMeanShort();
  }

  throw std::runtime_error("Event does not contain a valid data.");
}
float Event::GetPSD() {
  if (IsSingleHitEvt()) {
    return std::get<singleHits>(data).GetPSD();
  } else if (IsPSBar()) {
    return std::get<PSBar>(data).GetPSD();
  }

  throw std::runtime_error("Event does not contain a valid data.");
}

#ifdef WAVES
WaveForm Event::GetWaveForm() {
  if (IsSingleHitEvt()) {
    return std::get<singleHits>(data).GetWF();
  } else if (IsPSBar()) {
    return std::get<PSBar>(data).GetWF();
  }

  throw std::runtime_error("Event does not contain a valid data.");
}
#endif
} // namespace digiAnalysis