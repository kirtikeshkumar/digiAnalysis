#include "singleHits.h"
#include "WaveForm.h"
#include "diffWaveFrom.h"
#include "includes.hh"

// ClassImp(digiAnalysis::singleHits);

using namespace std;

namespace digiAnalysis {
/*Constructors*/
singleHits::singleHits() {
  hitNum = 0;
  ChNum = 0;
  Board = 0;
  Timestamp = 0;
  Energy = 0;
  EnergyShort = 0;
  PSD = 0;
}

singleHits::singleHits(const singleHits &other)
    : hitNum(other.hitNum), ChNum(other.ChNum), Board(other.Board),
      Timestamp(other.Timestamp), Energy(other.Energy),
      EnergyShort(other.EnergyShort), PSD(other.PSD) {
#ifdef WAVES
  // Deep copy of unique pointers
  if (other.WF) {
    WF = std::make_unique<WaveForm>(*other.WF);
  }
  // if (other.dWF) {
  //   dWF = std::make_unique<diffWaveForm>(*other.dWF);
  // }
#endif
}

singleHits::singleHits(ULong64_t EvNum, UShort_t channel, UShort_t board,
                       ULong64_t timestamp, UShort_t energy,
                       UShort_t energyshort) {
  hitNum = EvNum;
  ChNum = channel;
  Board = board;
  Timestamp = timestamp;
  Energy = energy;
  EnergyShort = energyshort;
  PSD = 1.0 - (EnergyShort * 1.0) / Energy;
}

#ifdef WAVES
singleHits::singleHits(ULong64_t EvNum, UShort_t channel, UShort_t board,
                       ULong64_t timestamp, UShort_t energy,
                       UShort_t energyshort, TArrayS *arr) {
  hitNum = EvNum;
  ChNum = channel;
  Board = board;
  Timestamp = timestamp;
  Energy = energy;
  EnergyShort = energyshort;
  PSD = 1.0 - (EnergyShort * 1.0) / Energy;
  WF = std::make_unique<WaveForm>(arr);
}
#endif

/*Destructor*/
singleHits::~singleHits() {}

/*Getters*/
ULong64_t singleHits::GetEvNum() { return hitNum; }
UShort_t singleHits::GetChNum() { return ChNum; }
UShort_t singleHits::GetBoard() { return Board; }
ULong64_t singleHits::GetTimestamp() { return Timestamp; }
UShort_t singleHits::GetEnergy() { return Energy; }
UShort_t singleHits::GetEnergyShort() { return EnergyShort; }
float singleHits::GetPSD() { return PSD; }

#ifdef WAVES
std::unique_ptr<WaveForm> singleHits::GetWF() { return std::move(WF); }
WaveForm *singleHits::GetWFPtr() { return WF.get(); }
float singleHits::GetMeanTime() { return WF->GetMeanTime(); }
// std::unique_ptr<diffWaveForm> singleHits::GetDiffWF() { return
// std::move(dWF); }

void singleHits::SetWF(const WaveForm &wf) { WF->SetWaveForm(wf); }
void singleHits::SetSmoothWF() { WF->SetSmooth(); }
// void singleHits::SetDiffWF() {}
#endif

void singleHits::SetPSD() { PSD = 1.0 - (EnergyShort * 1.0) / Energy; }

void singleHits::Print() {
  std::cout << "###############################################################"
            << std::endl;
  std::cout << "                    Single Hit Event                           "
            << std::endl;
  std::cout << "###############################################################"
            << std::endl;
  std::cout << "Hit Number    : " << hitNum << std::endl;
  std::cout << "Channel Number: " << ChNum << std::endl;
  std::cout << "Board         : " << Board << std::endl;
  std::cout << "TimeStamp     : " << Timestamp << std::endl;
  std::cout << "Energy        : " << Energy << std::endl;
  std::cout << "EnergyShort   : " << EnergyShort << std::endl;
  std::cout << "PSD           : " << PSD << std::endl;
#ifdef WAVES
  std::cout << "Waves         : " << "ON" << std::endl;
  std::cout << "WaveForm      : " << "WF" << std::endl;
#else
  std::cout << "Waves         : " << "OFF" << std::endl;
#endif
  std::cout << "###############################################################"
            << std::endl;
}
} // namespace digiAnalysis