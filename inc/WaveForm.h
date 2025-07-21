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
  std::vector<float> traces;
  std::vector<float> tracesSmooth;
  std::vector<float> CFDtraces;
  float meantime;
  float baseline;
  int blStart = GateStart - nSampleBL - 50 > 0 ? GateStart - nSampleBL - 10 : 0;
  TF1 *fitFunc = nullptr;

public:
  WaveForm();
  WaveForm(TArrayS *arr);
  WaveForm(const std::vector<float> tr);
  WaveForm(const WaveForm &wf);
  WaveForm(const WaveForm &wf1,
           const WaveForm &wf2); // constructor to concatenate waveforms
  WaveForm(UShort_t sizeOfWaveForms,
           const std::vector<WaveForm>
               vecOfWaveForm); // Constructor to take average of waveforms

  virtual ~WaveForm();

  std::vector<float> GetTraces();
  std::vector<float> GetTracesSmooth();
  float GetMeanTime();
  float GetBaseLine();
  UShort_t GetSize();
  bool IsFit();
  double GetFitPar(int val);
  double GetFitParError(int val);

  void SetWaveForm(const std::vector<float> tr);
  void SetWaveForm(const WaveForm &wf);
  void SetSmooth(UShort_t sBoxSz);
  void SetSmooth();
  void SetCFD();
  void SetMeanTime();
  void SetMeanTime(const std::vector<float> tr);
  void SetMeanTime(const std::vector<float> tr, UShort_t start, UShort_t stop);
  void SetBaseLine();
  void SetBaseLine(TArrayS *arr);
  void SetBaseLine(const std::vector<float> tr);

  void Plot();
  // std::vector<float> EvalDecayTime(UShort_t FitStart, UShort_t FitEnd,
  // UShort_t numDecayConst);
  void ShiftWaveForm(int BL);
  float IntegrateWaveForm();
  float IntegrateWaveForm(int startTime, int stopTime);
  void ConcatenateWaveForms(const WaveForm &wf1, const WaveForm &wf2);
  void AverageWaveForms(UShort_t sizeOfWaveForms,
                        const std::vector<WaveForm> vecOfWaveForm);
  void AverageWaveForms(ULong_t start, UShort_t numWaveForm,
                        UShort_t sizeOfWaveForms,
                        const std::vector<WaveForm> vecOfWaveForm);
  void ScaleWaveForm(double Scale);
  void AddWaveForm(const WaveForm &wf1);
  void FitExponential(int start, int stop);
  std::vector<std::unique_ptr<WaveForm>> SplitWaveForm(UShort_t numSplits);

  // static UShort_t nSampleBL;   // Number of baseline samples
  // static UShort_t smoothBoxSz; // Size of smoothing box

  //  ClassDef(WaveForm, 1);
};

} // namespace digiAnalysis
#endif /*WaveForm_h*/