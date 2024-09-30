/*
This File defines how to hande the waveform

class WaveForm: parent class contains how to handle waveform including plotting

class smoothWaveForm: derived from WaveForm. Evaluates and Stores the smoothed
waveform

class diffWaveForm: derived from WaveForm. Evaluates and Stores derivative of
Waveform. Useful for Slow-Fast Coincidence

class cfdWaveForm: derived from WaveForm. Evaluates and Stores CFD.
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
private:
  std::vector<UShort_t> traces;
  UShort_t meantime;

public:
  WaveForm();
  WaveForm(TArrayS *arr);
  WaveForm(std::vector<UShort_t> tr);
  WaveForm(const WaveForm &wf);
  virtual ~WaveForm();
  void Plot();
  std::vector<UShort_t> GetTraces();
  UShort_t GetMeanTime();

  void SetWaveForm(std::vector<UShort_t> tr);

  ClassDef(WaveForm, 1);
};

class smoothWaveForm : public WaveForm {
public:
  // to smooth existing waveform
  smoothWaveForm(const WaveForm &wf, int windowSize);
  // to smooth from existing vector
  smoothWaveForm(std::vector<UShort_t *> tr, int windowSize);
  // to smooth without creating waveform object
  smoothWaveForm(TArrayS *arr, int windowSize);

  void movingAverage();

  ClassDef(smoothWaveForm, 1);
};

class diffWaveForm : public WaveForm {
private:
  std::vector<float> difftraces;

public:
  diffWaveForm(const WaveForm &wf);
  diffWaveForm(const smoothWaveForm &swf);
  ~diffWaveForm() override;

  void computeDerivative();

  std::vector<float> GetDerivative();

  ClassDef(diffWaveForm, 1);
};

class cfdWaveForm : public WaveForm {
private:
  std::vector<float> cfdtraces;

public:
  cfdWaveForm(const WaveForm &wf);
  diffWaveForm(const smoothWaveForm &swf);
  ~cfdWaveForm() override;

  void computeCFD();

  std::vector<float> GetCFD();

  ClassDef(cfdWaveForm, 1);
};

} // namespace digiAnalysis
#endif /*WaveForm_h*/