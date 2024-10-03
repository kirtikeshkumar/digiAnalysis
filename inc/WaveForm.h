/*
This File defines how to hande the waveform
class WaveForm: parent class contains how to handle waveform including plotting
*/

#ifndef WaveForm_h
#define WaveForm_h

#include "includes.hh"
#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace digiAnalysis {
class WaveForm {
protected:
  std::vector<UShort_t> traces;
  std::vector<UShort_t> tracesSmooth;
  std::vector<UShort_t> CFDtraces;
  float meantime;

public:
  WaveForm();
  WaveForm(
      TArrayS *arr); // in definition for this, include if smoothing
                     // and cfd is to be done while storing using another option
  WaveForm(std::vector<UShort_t> tr);
  WaveForm(const WaveForm &wf);
  virtual ~WaveForm();
  void Plot();
  std::vector<UShort_t> GetTraces();
  float GetMeanTime();

  void SetWaveForm(std::vector<UShort_t> tr);
  void SetSmooth();
  void SetCFD();

  void EvalDecayTime(UShort_t FitStart, UShort_t FitEnd,
                     UShort_t numDecayConst);

  ClassDef(WaveForm, 1);
};

} // namespace digiAnalysis
#endif /*WaveForm_h*/