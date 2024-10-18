/*
This File defines how to hande the waveform
class WaveForm: parent class contains how to handle waveform including plotting
*/

#ifndef WaveForm_h
#define WaveForm_h

#include "globals.h"
#include "includes.hh"
#include <sys/types.h>
#pragma once

namespace digiAnalysis {
class WaveForm {
protected:
  std::vector<int> traces;
  std::vector<int> tracesSmooth;
  std::vector<int> CFDtraces;
  float meantime;
  float baseline;

public:
  WaveForm();
  WaveForm(TArrayS *arr);
  WaveForm(const std::vector<int> tr);
  WaveForm(const WaveForm &wf);
  WaveForm(const WaveForm &wf1,
           const WaveForm &wf2); // constructor to concatenate waveforms
  WaveForm(UShort_t numWaveForm, UShort_t sizeOfWaveForms,
           const std::vector<WaveForm>
               vecOfWaveForm); // Constructor to take average of waveforms

  virtual ~WaveForm();

  std::vector<int> GetTraces();
  std::vector<int> GetTracesSmooth();
  float GetMeanTime();
  float GetBaseLine();
  uint GetSize();

  void SetWaveForm(const std::vector<int> tr);
  void SetWaveForm(const WaveForm &wf);
  void SetSmooth(UShort_t sBoxSz);
  void SetSmooth();
  void SetCFD();
  void SetMeanTime();
  void SetMeanTime(const std::vector<int> tr);
  void SetMeanTime(const std::vector<int> tr, UShort_t start, UShort_t stop);
  void SetBaseLine();
  void SetBaseLine(TArrayS *arr);
  void SetBaseLine(const std::vector<int> tr);

  void Plot();
  // std::vector<float> EvalDecayTime(UShort_t FitStart, UShort_t FitEnd,
  // UShort_t numDecayConst);
  Long64_t IntegrateWaveForm();
  Long64_t IntegrateWaveForm(int startTime, int stopTime);
  void ConcatenateWaveForms(const WaveForm &wf1, const WaveForm &wf2);
  void AverageWaveForms(UShort_t numWaveForm, UShort_t sizeOfWaveForms,
                        const std::vector<WaveForm> vecOfWaveForm);
  std::vector<std::unique_ptr<WaveForm>> SplitWaveForm(UShort_t numSplits);

  // static UShort_t nSampleBL;   // Number of baseline samples
  // static UShort_t smoothBoxSz; // Size of smoothing box

  //  ClassDef(WaveForm, 1);
};

} // namespace digiAnalysis
#endif /*WaveForm_h*/