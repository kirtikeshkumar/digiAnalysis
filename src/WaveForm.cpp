#include "WaveForm.h"
#include "TMath.h"
#include <iostream>

using namespace std;
namespace digiAnalysis {
WaveForm::WaveForm() { meantime = 0.0; }
WaveForm::WaveForm(TArrayS *arr) {
#ifdef SMOOTH
  UShort_t smoothBoxSz = 4;
  unsigned int movingSum;
#endif
  meantime = 0;
  float sampleSum = 0;
  unsigned int size = arr->GetSize();
  UShort_t firstVal = arr->At(0);
  for (unsigned int j = 0; j < size; j++) {
    traces.push_back(firstVal - arr->At(j));
    meantime = meantime + traces[j] * j;
    sampleSum = sampleSum + traces[j];
#ifdef SMOOTH
    if (j < smoothBoxSz) {
      tracesSmooth.push_back(0);
      movingSum = movingSum + traces[j];
    }
    movingSum = movingSum + traces[j] - traces[j - smoothBoxSz];
    tracesSmooth = movingSum / smoothBoxSz;
#endif
#ifdef CFD
#endif
  }
  meantime = TMath::Log10(meantime / sampleSum);
}
WaveForm::WaveForm(std::vector<UShort_t> tr) {}
WaveForm::WaveForm(const WaveForm &wf) {
  meantime = wf.meantime;
  if (!wf.traces.empty()) {
  } // copy traces if not empty
  if (!wf.tracesSmooth.empty()) {
  } // copy traces smooth if not empty
  if (!wf.CFDtraces.empty()) {
  } // copy traces CFD if not empty
}
WaveForm::~WaveForm() {}
void WaveForm::Plot() {}
std::vector<UShort_t> WaveForm::GetTraces() { return traces; }
float WaveForm::GetMeanTime() { return meantime; }

void WaveForm::SetWaveForm(std::vector<UShort_t> tr) {}
void WaveForm::SetSmooth() {}
void WaveForm::SetCFD() {}
void WaveForm::EvalDecayTime(UShort_t FitStart, UShort_t FitEnd,
                             UShort_t numDecayConst) {}
} // namespace digiAnalysis
