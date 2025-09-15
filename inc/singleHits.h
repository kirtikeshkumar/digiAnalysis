/*
This file saves individual hits as registered
*/

#ifndef singleHits_h
#define singleHits_h

#include "WaveForm.h"
#include "globals.h"
// #include "diffWaveFrom.h"
#include "includes.hh"
#pragma once

namespace digiAnalysis {

class WaveForm;
// class diffWaveForm;

class singleHits {
private:
  // Compulsary Member Variables
  ULong64_t hitNum;
  UShort_t ChNum;
  UShort_t Board;
  ULong64_t Timestamp;
  UShort_t Energy;
  UShort_t EnergyShort;
  double PSD;

#ifdef WAVES
  // Needed for Waves Acquisition
  std::unique_ptr<WaveForm> WF;
  double evalEnergy;
  double evalEnergyShort;
  double evalPSD;

  // Initialize depending on situation
  // std::unique_ptr<diffWaveForm> dWF;
#endif

public:
  singleHits();
  singleHits(const singleHits &hit);
#ifdef WAVES
  singleHits(ULong64_t EvNum, UShort_t ChNum, UShort_t Board,
             ULong64_t Timestamp, UShort_t Energy, UShort_t EnergyShort,
             TArrayS *arr);
  singleHits(ULong64_t EvNum, UShort_t ChNum, UShort_t Board,
             ULong64_t Timestamp, UShort_t Energy, UShort_t EnergyShort,
             WaveForm *WF);
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
  double GetPSD();

#ifdef WAVES
  std::unique_ptr<WaveForm> GetWF();
  WaveForm *GetWFPtr();
  double GetMeanTime();
  double GetEvalEnergy();
  double GetEvalEnergyShort();
  double GetEvalPSD();
  // std::unique_ptr<diffWaveForm> GetDiffWF();

  // Setters
  void SetWF(const WaveForm &wf);
  void SetSmoothWF();
  void SetSmoothWF(UShort_t sBoxSz);
  void SetEvalEnergy();
  void SetEvalEnergyShort();
  void SetEvalPSD();
  // void SetDiffWF();
  void SetCFD();
  void SetMovBLCorr();
#endif
  void SetPSD();

  void Print();

  // ClassDef(singleHits, 1);
};
} // namespace digiAnalysis
#endif