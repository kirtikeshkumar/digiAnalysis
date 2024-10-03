/*
class diffWaveForm: derived from WaveForm. Evaluates and Stores derivative of
Waveform. Useful for Slow-Fast Coincidence
*/

#ifndef diffWaveForm_h
#define diffWaveForm_h

#include "WaveForm.h"
#include "includes.hh"
#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace digiAnalysis {
class WaveForm;

class diffWaveForm : public WaveForm {
private:
  std::vector<float> difftraces;

public:
  diffWaveForm();
  ~diffWaveForm();

  void computeDerivative();

  std::vector<float> GetDerivative();

  ClassDef(diffWaveForm, 1);
};

} // namespace digiAnalysis
#endif /*diffWaveForm_h*/