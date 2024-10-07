/*
This file saves individual hits as registered
*/

#ifndef singleHits_h
#define singleHits_h

#include "WaveForm.h"
#include "diffWaveFrom.h"
#include "includes.hh"
#pragma once

namespace digiAnalysis
{

  class WaveForm;
  class diffWaveForm;

  class singleHits
  {
  private:
    // Compulsary Member Variables
    ULong64_t hitNum;
    UShort_t ChNum;
    UShort_t Board;
    ULong64_t Timestamp;
    UShort_t Energy;
    UShort_t EnergyShort;

#ifdef WAVES
    // Needed for Waves Acquisition
    std::unique_ptr<WaveForm> WF;

    // Initialize depending on situation
    std::unique_ptr<diffWaveForm> dWF;
#endif

  public:
    singleHits();
#ifdef WAVES
    singleHits(ULong64_t EvNum, UShort_t ChNum, UShort_t Board,
               ULong64_t Timestamp, UShort_t Energy, UShort_t EnergyShort,
               TArrayS *arr);
#endif
    singleHits(ULong64_t EvNum, UShort_t ChNum, UShort_t Board,
               ULong64_t Timestamp, UShort_t Energy, UShort_t EnergyShort);
    ~singleHits();

    // Getters
    ULong64_t GetEvNum();
    UShort_t GetChNum();
    UShort_t GetBoard();
    ULong64_t GetTimestamp();
    UShort_t GetEnergy();
    UShort_t GetEnergyShort();

#ifdef WAVES
    std::unique_ptr<WaveForm> GetWF();
    std::unique_ptr<smoothWaveForm> GetSmoothWF();
    std::unique_ptr<diffWaveForm> GetDiffWF();
    std::unique_ptr<cfdWaveForm> GetCFD();

    // Setters
    void SetSmoothWF();
    void SetDiffWF();
    void SetCFD();
#endif

    // ClassDef(singleHits, 1);
  };
} // namespace digiAnalysis
#endif