#include "Pair.h"

using namespace std;

namespace digiAnalysis {

Pair::Pair() : hit1(nullptr), hit2(nullptr), SumEnergy(0), PairDelTime(0) {}

Pair::Pair(const singleHits *h1, const singleHits *h2)
    : hit1(const_cast<singleHits *>(h1)), hit2(const_cast<singleHits *>(h2)) {
  if (hit1 && hit2) {
    SumEnergy = hit1->GetEnergy() + hit2->GetEnergy();
    PairDelTime = hit1->GetTimestamp() > hit2->GetTimestamp()
                      ? hit1->GetTimestamp() - hit2->GetTimestamp()
                      : hit2->GetTimestamp() - hit1->GetTimestamp();
  } else {
    std::cout << "Pair not created" << std::endl;
    SumEnergy = 0;
    PairDelTime = 0;
  }
}

Pair::~Pair() {}

void Pair::SetPair(const singleHits *h1, const singleHits *h2) {
  hit1 = const_cast<singleHits *>(h1);
  hit2 = const_cast<singleHits *>(h2);

  if (hit1 && hit2) {
    SumEnergy = hit1->GetEnergy() + hit2->GetEnergy();
    PairDelTime = (hit1->GetTimestamp() > hit2->GetTimestamp())
                      ? hit1->GetTimestamp() - hit2->GetTimestamp()
                      : hit2->GetTimestamp() - hit1->GetTimestamp();
  } else {
    SumEnergy = 0;
    PairDelTime = 0;
  }
}

void Pair::ClearPair() {
  hit1 = nullptr;
  hit2 = nullptr;
  SumEnergy = 0;
  PairDelTime = 0;
}

singleHits *Pair::GetHit(Short_t Sel) {
  if (!hit1 || !hit2) {
    std::cout << "ERROR : one or both hits in pair are empty" << std::endl;
    return nullptr;
  }

  switch (Sel) {
  case -1:
    return (hit1->GetTimestamp() <= hit2->GetTimestamp()) ? hit1 : hit2;
  case -2:
    return (hit1->GetTimestamp() > hit2->GetTimestamp()) ? hit1 : hit2;
  case 0:
    return (hit1->GetChNum() <= hit2->GetChNum()) ? hit1 : hit2;
  case 1:
    return (hit1->GetChNum() > hit2->GetChNum()) ? hit1 : hit2;
  default:
    return nullptr;
  }
}

UShort_t Pair::GetPairHitCh(Short_t Sel) {
  auto h = GetHit(Sel);
  return h ? h->GetChNum() : 0;
}

ULong64_t Pair::GetPairHitTime(Short_t Sel) {
  auto h = GetHit(Sel);
  return h ? h->GetTimestamp() : 0;
}

ULong64_t Pair::GetPairDelTime() { return PairDelTime; }

UShort_t Pair::GetPairHitEnergy(Short_t Sel) {
  auto h = GetHit(Sel);
  return h ? h->GetEnergy() : 0;
}

UShort_t Pair::GetPairEnergy() { return SumEnergy; }

UShort_t Pair::GetPairHitEnergyShort(Short_t Sel) {
  auto h = GetHit(Sel);
  return h ? h->GetEnergyShort() : 0;
}

UShort_t Pair::GetPairHitPSD(Short_t Sel) {
  auto h = GetHit(Sel);
  return h ? h->GetPSD() : 0;
}

float Pair::GetPairHitEvalEnergy(Short_t Sel) {
  auto h = GetHit(Sel);
  return h ? h->GetEvalEnergy() : 0.0f;
}

float Pair::GetPairHitEvalEnergyShort(Short_t Sel) {
  auto h = GetHit(Sel);
  return h ? h->GetEvalEnergyShort() : 0.0f;
}

float Pair::GetPairHitEvalPSD(Short_t Sel) {
  auto h = GetHit(Sel);
  return h ? h->GetEvalPSD() : 0.0f;
}

void Pair::Print() {
  // Ensure both hits exist
  if (!hit1 || !hit2) {
    std::cout << "ERROR : Invalid Pair (one or both hits are null)"
              << std::endl;
    return;
  }

  // Formatting setup
  std::cout << "\n\n" << std::endl;

  std::cout << std::string(60, '*') << std::endl;
  std::cout << std::left << std::setw(20) << "Variable" << std::setw(20)
            << "Pair 0" << std::setw(20) << "Pair 1" << std::endl;

  std::cout << std::string(60, '-') << std::endl;

  // Print each variable

  std::cout << std::setw(20) << "ChNum" << std::setw(20) << hit1->GetChNum()
            << std::setw(20) << hit2->GetChNum() << std::endl;

  std::cout << std::setw(20) << "Board" << std::setw(20) << hit1->GetBoard()
            << std::setw(20) << hit2->GetBoard() << std::endl;

  std::cout << std::setw(20) << "hitNum" << std::setw(20) << hit1->GetEvNum()
            << std::setw(20) << hit2->GetEvNum() << std::endl;

  std::cout << std::setw(20) << "Timestamp" << std::setw(20)
            << hit1->GetTimestamp() << std::setw(20) << hit2->GetTimestamp()
            << std::endl;

  std::cout << std::setw(20) << "Energy" << std::setw(20) << hit1->GetEnergy()
            << std::setw(20) << hit2->GetEnergy() << std::endl;

  std::cout << std::setw(20) << "EnergyShort" << std::setw(20)
            << hit1->GetEnergyShort() << std::setw(20) << hit2->GetEnergyShort()
            << std::endl;

  std::cout << std::setw(20) << "PSD" << std::setw(20) << hit1->GetPSD()
            << std::setw(20) << hit2->GetPSD() << std::endl;

  std::cout << std::string(60, '-') << std::endl;

  std::cout << std::setw(20) << "SumEnergy" << std::setw(20) << SumEnergy
            << std::endl;

  std::cout << std::setw(20) << "DeltaTime" << std::setw(20)
            << PairDelTime / 1000.0 << std::setw(20) << "ns" << std::endl;

  std::cout << std::string(60, '*') << std::endl;
}

} // namespace digiAnalysis