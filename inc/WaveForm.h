/*
This File defines how to hande the waveform
class WaveForm: parent class contains how to handle waveform including plotting
*/

#ifndef WaveForm_h
#define WaveForm_h

#include "globals.h"
#include "includes.hh"
#include <TROOT.h>
#include <sys/types.h>
#pragma once

namespace digiAnalysis {
class WaveForm {
protected:
  std::vector<double> traces;
  std::vector<double> tracesSmooth;
  std::vector<double> tracesFFT;
  double meantime;
  double baseline;
  int blStart = GateStart - nSampleBL - 50 > 0 ? GateStart - nSampleBL - 10 : 0;
  TF1 *fitFunc = nullptr;

public:
  WaveForm();
  WaveForm(TArrayS *arr);
  WaveForm(const std::vector<double> tr);
  WaveForm(const WaveForm &wf);
  WaveForm(const WaveForm &wf1,
           const WaveForm &wf2); // constructor to concatenate waveforms
  WaveForm(UShort_t sizeOfWaveForms,
           const std::vector<WaveForm>
               vecOfWaveForm); // Constructor to take average of waveforms

  virtual ~WaveForm();

  std::vector<double> GetTraces();
  std::vector<double> GetTracesSmooth();
  std::vector<double> GetTracesFFT();
  double GetMeanTime();
  double GetBaseLine();
  UShort_t GetSize();
  bool IsFit();
  double GetFitPar(int val);
  double GetFitParError(int val);

  void SetWaveForm(const std::vector<double> tr);
  void SetWaveForm(const WaveForm &wf);
  /**
   * @brief Smooth Waveforms.
   *
   * Smooths the waveform saved in traces for the particular hit.
   * kernel size is 10 sigma for gaussian kernel
   * kernel types are
   *    "MovA" for flat moving average (Default)
   *    "Gauss" for gaussian kernel
   *
   * @param sBoxSz smoothing kernel size.
   * @param kernel type of kernel.
   */
  void SetSmooth(UShort_t sBoxSz, std::string kernel = "Gauss");
  void SetSmooth();
  void SetCFD();
  void SetMeanTime();
  void SetMeanTime(const std::vector<double> tr);
  void SetMeanTime(const std::vector<double> tr, UShort_t start, UShort_t stop);
  void SetBaseLine();
  void SetBaseLine(TArrayS *arr);
  void SetBaseLine(const std::vector<double> tr);
  void SetTracesFFT();
  void SetTracesFFT(std::string whichTrace);
  void SetTracesFFT(std::vector<double> trFFT);
  std::vector<double> EvalTracesFFT(std::vector<double> trFFT);

  void Plot();
  void ShiftWaveForm(int BL);
  double IntegrateWaveForm();
  double IntegrateWaveForm(int startTime, int stopTime);
  double IntegrateSmoothWaveForm(int startTime, int stopTime);
  void FitExponential(int start, int stop);
  std::vector<double> GenerateWaveFromFFT();

  void AverageWaveForms(UShort_t sizeOfWaveForms,
                        const std::vector<WaveForm> vecOfWaveForm);
  void AverageWaveForms(ULong_t start, UShort_t numWaveForm,
                        UShort_t sizeOfWaveForms,
                        const std::vector<WaveForm> vecOfWaveForm);
  void ScaleWaveForm(double Scale);
  void AddWaveForm(const WaveForm &wf1);
  void ConcatenateWaveForms(const WaveForm &wf1, const WaveForm &wf2);
  std::vector<std::unique_ptr<WaveForm>> SplitWaveForm(UShort_t numSplits);

  // static UShort_t nSampleBL;   // Number of baseline samples
  // static UShort_t smoothBoxSz; // Size of smoothing box

  //  ClassDef(WaveForm, 1);
};

} // namespace digiAnalysis
#endif /*WaveForm_h*/