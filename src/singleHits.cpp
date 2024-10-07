#include "singleHits.h"
#include "WaveForm.h"
#include "diffWaveFrom.h"
#include "includes.hh"

// ClassImp(digiAnalysis::singleHits);

using namespace std;

namespace digiAnalysis
{
    /*Constructors*/
    singleHits::singleHits()
    {
        hitNum = 0;
        ChNum = 0;
        Board = 0;
        Timestamp = 0;
        Energy = 0;
        EnergyShort = 0;
    }

    singleHits::singleHits(ULong64_t EvNum, UShort_t channel, UShort_t board,
                           ULong64_t timestamp, UShort_t energy, UShort_t energyshort)
    {
        hitNum = EvNum;
        ChNum = channel;
        Board = board;
        Timestamp = timestamp;
        Energy = energy;
        EnergyShort = energyshort;
    }

#ifdef WAVES
    singleHits::singleHits(ULong64_t EvNum, UShort_t channel, UShort_t board,
                           ULong64_t timestamp, UShort_t energy, UShort_t energyshort,
                           TArrayS *arr)
    {
        hitNum = EvNum;
        ChNum = channel;
        Board = board;
        Timestamp = timestamp;
        Energy = energy;
        EnergyShort = energyshort;
        WF = std::make_unique<WaveForm>(arr);
    }
#endif

    /*Destructor*/
    singleHits::~singleHits() {}

    /*Getters*/
    ULong64_t singleHits::GetEvNum()
    {
        return hitNum;
    }
    UShort_t singleHits::GetChNum()
    {
        return ChNum;
    }
    UShort_t singleHits::GetBoard()
    {
        return Board;
    }
    ULong64_t singleHits::GetTimestamp()
    {
        return Timestamp;
    }
    UShort_t singleHits::GetEnergy()
    {
        return Energy;
    }
    UShort_t singleHits::GetEnergyShort()
    {
        return EnergyShort;
    }

#ifdef WAVES
    std::unique_ptr<WaveForm> singleHits::GetWF()
    {
        return WF;
    }
    std::unique_ptr<diffWaveForm> singleHits::GetDiffWF()
    {
        return dWF;
    }
#endif
}