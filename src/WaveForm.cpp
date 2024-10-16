/*define meantime using tracesSmooth if SMOOTH is turned on*/

#include "WaveForm.h"
#include "includes.hh"

// ClassImp(digiAnalysis::WaveForm);

using namespace std;

namespace digiAnalysis {
// UShort_t WaveForm::nSampleBL = 16;
const char *env_var_baselinesamples = std::getenv("NUM_SAMPLE_BASELINE");
UShort_t WaveForm::nSampleBL =
    static_cast<UShort_t>(std::stoul(env_var_baselinesamples));

#ifdef SMOOTH
const char *envVar = std::getenv("SMOOTH_BOX_SIZE");
UShort_t WaveForm::smoothBoxSz = static_cast<UShort_t>(std::stoul(envVar));
#endif

/*Constructors*/
WaveForm::WaveForm() // Default Constructor initializes empty vectors
{
  meantime = 0.0;
  baseline = 0.0;
}

WaveForm::WaveForm(TArrayS *arr) // CoMPASS saves waveforms as TArrayS
{
#ifdef SMOOTH
  unsigned int movingSum = 0;
#endif
  meantime = 0;
  baseline = 0;
  SetBaseLine(arr); // Alwasy Set Base Line First
  float sampleSum = 0;
  unsigned int size = arr->GetSize();
  for (unsigned int j = 0; j < size; j++) {
    traces.push_back(baseline - arr->At(j));
    meantime = meantime + traces[j] * j;
    sampleSum = sampleSum + traces[j];
#ifdef SMOOTH
    if (WaveForm::smoothBoxSz == 1) {
      tracesSmooth.push_back(traces[j]);
    } else {
      if (j < WaveForm::smoothBoxSz) {
        tracesSmooth.push_back(0);
        movingSum = movingSum + traces[j];
      }
      movingSum = movingSum + traces[j] - traces[j - WaveForm::smoothBoxSz];
      tracesSmooth.push_back(movingSum / WaveForm::smoothBoxSz);
    }
#endif
  }
  meantime = TMath::Log10(meantime / sampleSum);
}

WaveForm::WaveForm(const std::vector<UShort_t> tr) // copy from other vector
{
  if (!tr.empty()) {
    SetBaseLine(tr);
    if (std::fabs(baseline - tr[0]) >
        100) // check if the copied traces have baseline subtracted
    {
      SetWaveForm(tr);
    } else {
      traces = tr;
      SetMeanTime();
#ifdef SMOOTH
      SetSmooth(WaveForm::smoothBoxSz);
#endif
    }
  } else {
    std::cout << "err WaveForm: input vector is empty" << std::endl;
  }
}

WaveForm::WaveForm(const WaveForm &wf) // copy consructor
{
  meantime = wf.meantime;
  baseline = wf.baseline;
  if (!wf.traces.empty()) // copy traces if not empty
  {
    traces = wf.traces;
  }
  if (!wf.tracesSmooth.empty()) // copy traces smooth if not empty
  {
    tracesSmooth = wf.tracesSmooth;
  }
  if (!wf.CFDtraces.empty()) // copy traces CFD if not empty
  {
    CFDtraces = wf.CFDtraces;
  }
}

WaveForm::WaveForm(const WaveForm &wf1, const WaveForm &wf2) {
  ConcatenateWaveForms(wf1, wf2);
}

WaveForm::WaveForm(UShort_t numWaveForm, UShort_t sizeOfWaveForms,
                   const std::vector<WaveForm> vecOfWaveForm) {
  AverageWaveForms(numWaveForm, sizeOfWaveForms, vecOfWaveForm);
}

/*Destructor*/
WaveForm::~WaveForm() {
  /*clear() Removes elements
  shrink_to_fit() frees unused memory*/
  /*
  traces.clear();
  traces.shrink_to_fit();

  tracesSmooth.clear();
  tracesSmooth.shrink_to_fit();

  CFDtraces.clear();
  CFDtraces.shrink_to_fit();
  */
}

void WaveForm::Plot() {}

/*Getters*/
std::vector<UShort_t> WaveForm::GetTraces() { return traces; }
std::vector<UShort_t> WaveForm::GetTracesSmooth() { return tracesSmooth; }
float WaveForm::GetMeanTime() { return meantime; }
float WaveForm::GetBaseLine() { return baseline; }

/*Setters*/
void WaveForm::SetWaveForm(std::vector<UShort_t> tr) {
#ifdef SMOOTH
  unsigned int movingSum = 0;
#endif
  if (!tr.empty()) {
    SetBaseLine(tr);
    meantime = 0;
    float sampleSum = 0;
    unsigned int size = tr.size();
    for (unsigned int j = 0; j < size; j++) {
      traces.push_back(baseline - tr[j]);
      meantime = meantime + traces[j] * j;
      sampleSum = sampleSum + traces[j];
#ifdef SMOOTH
      if (WaveForm::smoothBoxSz == 1) {
        tracesSmooth.push_back(traces[j]);
      } else {
        if (j < WaveForm::smoothBoxSz) {
          tracesSmooth.push_back(0);
          movingSum = movingSum + traces[j];
        }
        movingSum = movingSum + traces[j] - traces[j - WaveForm::smoothBoxSz];
        tracesSmooth.push_back(movingSum / WaveForm::smoothBoxSz);
      }
#endif
    }
    meantime = TMath::Log10(meantime / sampleSum);
  } else {
    std::cout << "err SetWaveForm: input vector is empty" << std::endl;
  }
}

void WaveForm::SetWaveForm(const WaveForm &wf) {
  traces = wf.traces;
  tracesSmooth = wf.tracesSmooth;
  CFDtraces = wf.CFDtraces;

  // Copying the scalar values
  meantime = wf.meantime;
  baseline = wf.baseline;
}

void WaveForm::SetSmooth(UShort_t smoothBoxSz) {
  if (!traces.empty()) {
    unsigned int movingSum = 0;
    SetBaseLine();
    unsigned int size = traces.size();
    if (smoothBoxSz == 1) {
      tracesSmooth = traces;
    } else {
      for (unsigned int j = 0; j < size; j++) {
        if (j < smoothBoxSz) {
          tracesSmooth.push_back(0);
          movingSum = movingSum + traces[j];
        }
        movingSum = movingSum + traces[j] - traces[j - smoothBoxSz];
        tracesSmooth.push_back(movingSum / smoothBoxSz);
      }
    }
  } else {
    std::cout << "err SetSmooth: Fill traces first" << std::endl;
  }
}

void WaveForm::SetCFD() {}

void WaveForm::SetMeanTime() {
  meantime = 0;
  float sampleSum = 0;
  if (!traces.empty()) {
    unsigned int size = traces.size();
    for (unsigned int j = 0; j < size; j++) {
      meantime = meantime + traces[j] * j;
      sampleSum = sampleSum + traces[j];
    }
    meantime = TMath::Log10(meantime / sampleSum);
  } else {
    std::cout << "err SetMeanTime: traces is empty" << std::endl;
  }
}

void WaveForm::SetMeanTime(const std::vector<UShort_t> tr) {
  meantime = 0;
  float sampleSum = 0;
  if (!tr.empty()) {
    unsigned int size = tr.size();
    for (unsigned int j = 0; j < size; j++) {
      meantime = meantime + tr[j] * j;
      sampleSum = sampleSum + tr[j];
    }
    meantime = TMath::Log10(meantime / sampleSum);
  } else {
    std::cout << "err SetMeanTime: input vector is empty" << std::endl;
  }
}

void WaveForm::SetMeanTime(const std::vector<UShort_t> tr, UShort_t start,
                           UShort_t stop) {
  meantime = 0;
  float sampleSum = 0;
  if (!tr.empty()) {
    if (start < stop) {
      if (tr.size() < stop) {
        std::cout << "Warning: stop is less than input vector size. Evaluating "
                     "till end"
                  << std::endl;
        stop = tr.size();
      }
      for (unsigned int j = start; j < stop; j++) {
        meantime = meantime + tr[j] * j;
        sampleSum = sampleSum + tr[j];
      }
      meantime = TMath::Log10(meantime / sampleSum);
    } else {
      std::cout << "error: start > stop, exiting witout meantime setting"
                << std::endl;
    }
  } else {
    std::cout << "err SetMeanTime: input vector is empty" << std::endl;
  }
}

void WaveForm::SetBaseLine() {
  baseline = 0;
  float sum = 0;
  if (!traces.empty()) {
    for (unsigned int j = 0; j < WaveForm::nSampleBL; j++) {
      sum = sum + traces[j];
    }
    baseline = sum / WaveForm::nSampleBL;
  } else {
    std::cout << "err SetBaseLine: traces is empty" << std::endl;
  }
}

void WaveForm::SetBaseLine(std::vector<UShort_t> tr) {
  baseline = 0;
  float sum = 0;
  if (!tr.empty()) {
    for (unsigned int j = 0; j < WaveForm::nSampleBL; j++) {
      sum = sum + tr[j];
    }
    baseline = sum / WaveForm::nSampleBL;
  } else {
    std::cout << "err SetBaseLine: input vector is empty" << std::endl;
  }
}

void WaveForm::SetBaseLine(TArrayS *arr) {
  baseline = 0;
  float sum = 0;
  if (arr && (arr->GetSize() > WaveForm::nSampleBL)) {
    for (unsigned int j = 0; j < WaveForm::nSampleBL; j++) {
      sum = sum + arr->At(j);
    }
    baseline = sum / WaveForm::nSampleBL;
  } else {
    std::cout << "err SetBaseLine: input array is empty or smaller than "
                 "NUM_SAMPLE_BASELINE"
              << std::endl;
  }
}

// std::vector<float> WaveForm::EvalDecayTime(UShort_t FitStart, UShort_t
// FitEnd, UShort_t numDecayConst) {}
Long64_t WaveForm::IntegrateWaveForm() {
  Long64_t sum = 0;
  if (!traces.empty()) {
    sum = IntegrateWaveForm(0, traces.size());
  } else {
    std::cout << "err IntegrateWaveForm: fill traces first" << std::endl;
  }
  return sum;
}

Long64_t WaveForm::IntegrateWaveForm(int startTime, int stopTime) {
  Long64_t sum = 0;
  if (!traces.empty() && (traces.size() > stopTime)) {
    if (startTime < stopTime) {
      for (unsigned int j = startTime; j < stopTime; j++) {
        sum = sum + traces[j];
      }
    } else {
      std::cout << "err IntegrateWaveForm: order the times properly "
                << std::endl;
    }
  } else {
    std::cout << "err IntegrateWaveForm: fill traces first" << std::endl;
  }
  return sum;
}

void WaveForm::ConcatenateWaveForms(const WaveForm &wf1, const WaveForm &wf2) {
  if (!wf1.traces.empty()) // copy traces if not empty
  {
    traces = wf1.traces;
  }
  if (!wf1.tracesSmooth.empty()) // copy traces smooth if not empty
  {
    tracesSmooth = wf1.tracesSmooth;
  }
  if (!wf1.CFDtraces.empty()) // copy traces CFD if not empty
  {
    CFDtraces = wf1.CFDtraces;
  }

  if (!wf2.traces.empty()) // concatenate traces if not empty
  {
    traces.insert(traces.end(), wf2.traces.begin(), wf2.traces.end());
  }
  if (!wf2.tracesSmooth.empty()) // concatenate traces smooth if not empty
  {
    tracesSmooth.insert(tracesSmooth.end(), wf2.tracesSmooth.begin(),
                        wf2.tracesSmooth.end());
  }
  if (!wf2.CFDtraces.empty()) // concatenate traces CFD if not empty
  {
    CFDtraces.insert(CFDtraces.end(), wf2.CFDtraces.begin(),
                     wf2.CFDtraces.end());
  }

  meantime =
      traces.size() +
      1000; // can use this to check if the waveform is single or concatenated
  baseline = wf1.baseline;
}

void WaveForm::AverageWaveForms(UShort_t numWaveForm, UShort_t sizeOfWaveForms,
                                const std::vector<WaveForm> vecOfWaveForm) {
  // currently only traces is averaged
  UShort_t sum = 0;
  baseline = 0;
  for (unsigned int i = 0; i < WaveForm::nSampleBL; i++) {
    for (unsigned int j = 0; j < numWaveForm; j++) {
      sum = sum + vecOfWaveForm[j].traces[i] / numWaveForm;
    }
    baseline = sum / WaveForm::nSampleBL;
  }

  for (unsigned int i = 0; i < sizeOfWaveForms; i++) {
    sum = 0;
    for (unsigned int j = 0; j < numWaveForm; j++) {
      sum = sum + vecOfWaveForm[j].traces[i];
    }
    sum = baseline - sum / numWaveForm;
    traces.push_back(sum);
  }
  SetMeanTime();
}

std::vector<std::unique_ptr<WaveForm>>
WaveForm::SplitWaveForm(UShort_t numSplits) {

  std::vector<std::unique_ptr<WaveForm>> splitWaveForms;

  int traceSize = traces.size();
  int splitSize = traceSize / numSplits;

  // Handle cases where the number of splits is greater than the size
  if (splitSize == 0) {
    // If numSplits is too large, return empty vector
    return splitWaveForms;
  }

  // Loop to create each split waveform
  for (int i = 0; i < numSplits; ++i) {
    // Create a new unique_ptr for each new WaveForm object
    auto newWaveForm = std::make_unique<WaveForm>();

    // Define the range for the current split
    int start = i * splitSize;
    int end = (i == numSplits - 1) ? traceSize : (i + 1) * splitSize;

    // Split the traces and assign to the new waveform
    if (!traces.empty()) {
      newWaveForm->traces.assign(traces.begin() + start, traces.begin() + end);
    } else {
      std::cout << "err SplitWaveForm:  Empty traces passed for Splitting"
                << std::endl;
      return splitWaveForms;
    }
    if (!tracesSmooth.empty()) {
      newWaveForm->tracesSmooth.assign(tracesSmooth.begin() + start,
                                       tracesSmooth.begin() + end);
    }
    if (!CFDtraces.empty()) {
      newWaveForm->CFDtraces.assign(CFDtraces.begin() + start,
                                    CFDtraces.begin() + end);
    }

    newWaveForm->SetBaseLine();
    newWaveForm->SetMeanTime();

    // Add the new waveform to the vector
    splitWaveForms.push_back(std::move(newWaveForm));
  }

  return splitWaveForms;
}
} // namespace digiAnalysis
