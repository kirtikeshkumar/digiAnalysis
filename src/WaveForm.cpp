/*define meantime using tracesSmooth if SMOOTH is turned on*/

#include "WaveForm.h"
#include "Rtypes.h"
#include "TMath.h"
#include "globals.h"
#include "includes.hh"

// ClassImp(digiAnalysis::WaveForm);

using namespace std;

namespace digiAnalysis {
/*// UShort_t WaveForm::nSampleBL = 16;
const char *env_var_baselinesamples = std::getenv("NUM_SAMPLE_BASELINE");
UShort_t WaveForm::nSampleBL =
    static_cast<UShort_t>(std::stoul(env_var_baselinesamples));

#ifdef SMOOTH
const char *envVar = std::getenv("SMOOTH_BOX_SIZE");
UShort_t WaveForm::smoothBoxSz = static_cast<UShort_t>(std::stoul(envVar));
#endif
*/

/*Constructors*/
WaveForm::WaveForm() // Default Constructor initializes empty vectors
{
  meantime = 0.0;
  baseline = 0.0;
}

WaveForm::WaveForm(TArrayS *arr) // CoMPASS saves waveforms as TArrayS
{
#ifdef SMOOTH
  float movingSum = 0;
#endif
  meantime = 0;
  baseline = 0;
  SetBaseLine(arr); // Alwasy Set Base Line First
  float sampleSum = 0;
  unsigned int size = arr->GetSize();
  for (unsigned int j = 0; j < size; j++) {
    traces.push_back(baseline - arr->At(j));
    if (j >= GateStart and j <= GateStart + GateLenLong) {
      meantime = meantime + traces[j] * j;
      sampleSum = sampleSum + traces[j];
    }
#ifdef SMOOTH
    if (smoothBoxSz == 1) {
      tracesSmooth.push_back(traces[j]);
    } else {
      if (j < smoothBoxSz) {
        // tracesSmooth.push_back(0);
        movingSum = movingSum + traces[j];
        tracesSmooth.push_back(movingSum / j);
      }
      movingSum = movingSum + traces[j] - traces[j - smoothBoxSz];
      tracesSmooth.push_back(movingSum / smoothBoxSz);
    }
#endif
  }
  meantime = (GateLenLong * 0.5 + GateStart) - meantime / sampleSum;
  if (meantime > 0.0) {
    meantime = TMath::Log10(meantime);
  } else {
    meantime = -1.0 * TMath::Log10(fabs(meantime));
  }
}

WaveForm::WaveForm(const std::vector<float> tr) // copy from other vector
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
      SetSmooth(smoothBoxSz);
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
  if (wf.fitFunc) // copy traces CFD if not empty
  {
    fitFunc = wf.fitFunc;
  }
}

WaveForm::WaveForm(const WaveForm &wf1, const WaveForm &wf2) {
  ConcatenateWaveForms(wf1, wf2);
}

WaveForm::WaveForm(UShort_t sizeOfWaveForms,
                   const std::vector<WaveForm> vecOfWaveForm) {
  AverageWaveForms(sizeOfWaveForms, vecOfWaveForm);
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

void WaveForm::Plot() {
  if (traces.empty() && tracesSmooth.empty()) {
    std::cout
        << "Both traces and tracesSmooth are empty. No plot will be created."
        << std::endl;
    return;
  }

  // Create a canvas
  TObject *obj = gROOT->FindObject("canvas");
  TCanvas *canvas = dynamic_cast<TCanvas *>(obj);
  if (canvas) {
    canvas->cd();    // make it current
    canvas->Clear(); // optional: clear previous plot
    std::cout << "Reusing existing canvas: " << "canvas" << std::endl;
  } else {
    canvas = new TCanvas("canvas", "WaveForm Plot", 800, 600);
    std::cout << "Created new canvas: " << "canvas" << std::endl;
  }
  // TCanvas *canvas = new TCanvas("canvas", "WaveForm Plot", 800, 600);

  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Create a legend
  int nTraces = traces.size();
  int nTracesSmooth = tracesSmooth.size();

  TGraph *graphTraces = nullptr;
  TGraph *graphTracesSmooth = nullptr;

  // Plot traces if not empty
  if (!traces.empty()) {

    graphTraces = new TGraph(nTraces);
    for (int i = 0; i < nTraces; ++i) {
      graphTraces->SetPoint(i, i, traces[i]);
    }
    graphTraces->SetLineColor(kBlue);
    graphTraces->SetLineWidth(2);
    graphTraces->SetTitle("Traces");
    graphTraces->Draw("AL");
    legend->AddEntry(graphTraces, "Traces", "l");
  }

  // Plot tracesSmooth if not empty
  if (!tracesSmooth.empty()) {
    graphTracesSmooth = new TGraph(nTracesSmooth);
    for (int i = 0; i < nTracesSmooth; ++i) {
      graphTracesSmooth->SetPoint(i, i, tracesSmooth[i]);
    }
    graphTracesSmooth->SetLineColor(kRed);
    graphTracesSmooth->SetLineWidth(2);
    graphTracesSmooth->SetTitle("Smoothed Traces");

    if (graphTraces) {
      graphTracesSmooth->Draw("L SAME");
    } else {
      graphTracesSmooth->Draw("AL");
    }
    legend->AddEntry(graphTracesSmooth, "Smoothed Traces", "l");
  }

  TLine *line = new TLine(0.0, 0.0, 5000, 0.0);
  line->SetLineColor(kBlack);
  line->SetLineWidth(2);
  line->Draw("LSAME");
  // Draw Fit if present
  if (fitFunc) {
    fitFunc->SetLineColor(kGreen);
    fitFunc->SetLineWidth(2);
    fitFunc->Draw("LSAME");
  }
  // Draw the legend
  legend->Draw();

  // Update the canvas
  canvas->Update();
}

/*Getters*/
std::vector<float> WaveForm::GetTraces() { return traces; }
std::vector<float> WaveForm::GetTracesSmooth() { return tracesSmooth; }
float WaveForm::GetMeanTime() { return meantime; }
float WaveForm::GetBaseLine() { return baseline; }
UShort_t WaveForm::GetSize() { return traces.size(); }
bool WaveForm::IsFit() {
  if (fitFunc) {
    return true;
  } else {
    return false;
  }
}
double WaveForm::GetFitPar(int val) { return fitFunc->GetParameter(val); }
double WaveForm::GetFitParError(int val) { return fitFunc->GetParError(val); }

/*Setters*/
void WaveForm::SetWaveForm(std::vector<float> tr) {
#ifdef SMOOTH
  float movingSum = 0;
#endif
  if (!tr.empty()) {
    SetBaseLine(tr);
    meantime = 0;
    float sampleSum = 0;
    unsigned int size = tr.size();
    for (unsigned int j = 0; j < size; j++) {
      traces.push_back(baseline - tr[j]);
      if (j >= GateStart and j <= GateStart + GateLenLong) {
        meantime = meantime + traces[j] * j;
        sampleSum = sampleSum + traces[j];
      }
#ifdef SMOOTH
      if (smoothBoxSz == 1) {
        tracesSmooth.push_back(traces[j]);
      } else {
        if (j < smoothBoxSz) {
          // tracesSmooth.push_back(0);
          movingSum = movingSum + traces[j];
          tracesSmooth.push_back(movingSum / j);
        }
        movingSum = movingSum + traces[j] - traces[j - smoothBoxSz];
        tracesSmooth.push_back(movingSum / smoothBoxSz);
      }
#endif
    }
    meantime = (GateStart + GateLenLong * 0.5) - (meantime / sampleSum);
    if (meantime > 0.0) {
      meantime = TMath::Log10(meantime);
    } else {
      meantime = -1.0 * TMath::Log10(fabs(meantime));
    }
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

void WaveForm::SetSmooth(UShort_t sBoxSz) {
  // std::cout << "Smoothing with boxSz: " << sBoxSz << std::endl;
  if (!traces.empty()) {
    float movingSum = 0;
    SetBaseLine();
    // std::cout << "Baseline set to: " << baseline << std::endl;
    unsigned int size = traces.size();
    // std::cout << "traces size: " << size << std::endl;
    if (sBoxSz == 1) {
      tracesSmooth = traces;
    } else {
      for (unsigned int j = 0; j < size; j++) {
        if (j < sBoxSz) {
          tracesSmooth.push_back(0);
          movingSum = movingSum + traces[j];
          // tracesSmooth.push_back(movingSum / j);
        } else {
          movingSum = movingSum + traces[j] - traces[j - sBoxSz];
          tracesSmooth.push_back(movingSum / sBoxSz);
        }
      }
    }
  } else {
    std::cout << "err SetSmooth: Fill traces first" << std::endl;
  }
}

void WaveForm::SetSmooth() {
  if (!traces.empty()) {
    float movingSum = 0;
    SetBaseLine();
    unsigned int size = traces.size();
    if (smoothBoxSz == 1) {
      tracesSmooth = traces;
    } else {
      for (unsigned int j = 0; j < size; j++) {
        if (j < smoothBoxSz) {
          tracesSmooth.push_back(0);
          movingSum = movingSum + traces[j];
          // tracesSmooth.push_back(movingSum / j);
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
    for (unsigned int j = GateStart; j < GateLenLong + GateStart; j++) {
      meantime = meantime + traces[j] * j;
      sampleSum = sampleSum + traces[j];
    }
    meantime = (GateLenLong * 0.5 + GateStart) - meantime / sampleSum;
    if (meantime > 0.0) {
      meantime = TMath::Log10(meantime);
    } else {
      meantime = -1.0 * TMath::Log10(fabs(meantime));
    }
  } else {
    std::cout << "err SetMeanTime: traces is empty" << std::endl;
  }
}

void WaveForm::SetMeanTime(const std::vector<float> tr) {
  if (!tr.empty()) {
    SetMeanTime(tr, GateStart, GateLenLong + GateStart);
  } else {
    std::cout << "err SetMeanTime: input vector is empty" << std::endl;
  }
}

void WaveForm::SetMeanTime(const std::vector<float> tr, UShort_t start,
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
      meantime = (GateLenLong * 0.5 + GateStart) - meantime / sampleSum;
      if (meantime > 0.0) {
        meantime = TMath::Log10(meantime);
      } else {
        meantime = -1.0 * TMath::Log10(fabs(meantime));
      }
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
    for (unsigned int j = blStart; j < nSampleBL + blStart; j++) {
      sum = sum + traces[j];
    }
    baseline = sum / nSampleBL;
  } else {
    std::cout << "err SetBaseLine: traces is empty" << std::endl;
  }
}

void WaveForm::SetBaseLine(std::vector<float> tr) {
  baseline = 0;
  float sum = 0;

  if (!tr.empty()) {
    for (unsigned int j = blStart; j < nSampleBL + blStart; j++) {
      sum = sum + tr[j];
    }
    baseline = sum / nSampleBL;
  } else {
    std::cout << "err SetBaseLine: input vector is empty" << std::endl;
  }
}

void WaveForm::SetBaseLine(TArrayS *arr) {
  baseline = 0;
  float sum = 0;

  if (arr && (arr->GetSize() > nSampleBL)) {
    for (unsigned int j = blStart; j < nSampleBL + blStart; j++) {
      sum = sum + arr->At(j);
    }
    baseline = sum / nSampleBL;
  } else {
    std::cout << "err SetBaseLine: input array is empty or smaller than "
                 "NUM_SAMPLE_BASELINE"
              << std::endl;
  }
}

// std::vector<float> WaveForm::EvalDecayTime(UShort_t FitStart, UShort_t
// FitEnd, UShort_t numDecayConst) {}

void WaveForm::ShiftWaveForm(int BL) {
  if (!traces.empty()) {
    for (unsigned int j = 0; j < traces.size(); j++) {
      traces[j] = traces[j] + BL;
    }
  }
}

float WaveForm::IntegrateWaveForm() {
  float sum = 0;
  if (!traces.empty()) {
    sum = IntegrateWaveForm(0, traces.size());
  } else {
    std::cout << "err IntegrateWaveForm: fill traces first" << std::endl;
  }
  return sum;
}

float WaveForm::IntegrateWaveForm(int startTime, int stopTime) {
  float sum = 0;
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

float WaveForm::IntegrateSmoothWaveForm(int startTime, int stopTime) {
  float sum = 0;
  if (!tracesSmooth.empty() && (tracesSmooth.size() > stopTime)) {
    if (startTime < stopTime) {
      for (unsigned int j = startTime; j < stopTime; j++) {
        sum = sum + tracesSmooth[j];
      }
    } else {
      std::cout << "err IntegrateSmoothWaveForm: order the times properly "
                << std::endl;
    }
  } else {
    std::cout << "err IntegrateSmoothWaveForm: fill traces Smooth first"
              << std::endl;
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

void WaveForm::AverageWaveForms(UShort_t sizeOfWaveForms,
                                const std::vector<WaveForm> vecOfWaveForm) {
  // currently only traces is averaged
  UShort_t numWaveForm = vecOfWaveForm.size();
  float sum = 0;
  baseline = 0;

  for (unsigned int i = blStart; i < nSampleBL + blStart; i++) {
    for (unsigned int j = 0; j < numWaveForm; j++) {
      sum = sum + vecOfWaveForm[j].traces[i] / numWaveForm;
    }
  }
  baseline = sum / nSampleBL;

  for (unsigned int i = 0; i < sizeOfWaveForms; i++) {
    sum = 0;
    for (unsigned int j = 0; j < numWaveForm; j++) {
      sum = sum + vecOfWaveForm[j].traces[i];
    }
    if (baseline > 10) {
      sum = baseline - sum / numWaveForm;
    } else {
      sum = sum / numWaveForm;
    }
    traces.push_back(sum);
  }
  SetMeanTime();
}

void WaveForm::AverageWaveForms(ULong_t start, UShort_t numWaveForm,
                                UShort_t sizeOfWaveForms,
                                const std::vector<WaveForm> vecOfWaveForm) {
  // currently only traces is averaged
  float sum = 0;
  baseline = 0;

  for (unsigned int i = blStart; i < nSampleBL + blStart; i++) {
    for (unsigned int j = start; j < start + numWaveForm; j++) {
      sum = sum + vecOfWaveForm[j].traces[i] / numWaveForm;
    }
  }
  baseline = sum / nSampleBL;

  for (unsigned int i = 0; i < sizeOfWaveForms; i++) {
    sum = 0;
    for (unsigned int j = start; j < start + numWaveForm; j++) {
      sum = sum + vecOfWaveForm[j].traces[i];
    }
    if (baseline > 10) {
      sum = baseline - sum / numWaveForm;
    } else {
      sum = sum / numWaveForm;
    }
    traces.push_back(sum);
  }
  SetMeanTime();
}

void WaveForm::ScaleWaveForm(double Scale) {
  std::transform(traces.begin(), traces.end(), traces.begin(),
                 [Scale](float val) { return val * Scale; });
  if (!tracesSmooth.empty()) {
    std::transform(tracesSmooth.begin(), tracesSmooth.end(),
                   tracesSmooth.begin(),
                   [Scale](float val) { return val * Scale; });
  }
}

void WaveForm::AddWaveForm(const WaveForm &wf1) {
  if (traces.size() == wf1.traces.size()) {
    // Add each element from traces and wf1.traces using std::transform
    std::transform(traces.begin(), traces.end(), wf1.traces.begin(),
                   traces.begin(), std::plus<float>());
  } else {
    std::cout << "traces must be the same size for adding! Nothing done!!!!!!!!"
              << std::endl;
  }
  if (tracesSmooth.size() == wf1.tracesSmooth.size()) {
    std::transform(tracesSmooth.begin(), tracesSmooth.end(),
                   wf1.tracesSmooth.begin(), tracesSmooth.begin(),
                   std::plus<float>());
  } else {
    std::cout << "tracesSmooth not of same size hence not adding!" << std::endl;
  }
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

void WaveForm::FitExponential(int start, int stop) {
  // Check valid range
  if (start < 0 || stop >= GetSize() || start >= stop) {
    std::cerr << "Invalid fit range: [" << start << ", " << stop << "]\n";
  }

  int nPoints = stop - start + 1;
  auto graph = std::make_unique<TGraph>(nPoints);

  for (int i = 0; i < nPoints; ++i) {
    graph->SetPoint(i, start + i,
                    traces[start + i]); // x = index (or time), y = value
  }

  fitFunc = new TF1("expFit", "[0]*exp(-x/[1])", start, stop);
  fitFunc->SetParameters(traces[start], 100); // Initial guesses

  graph->Fit(fitFunc, "QR"); // Q = quiet, R = respect fit range

  // double A = fitFunc->GetParameter(0);
  // double tau = fitFunc->GetParameter(1);
  // double Aerr = fitFunc->GetParError(0);
  // double tauErr = fitFunc->GetParError(1);

  // std::cout << "Fit results ([0]*exp(-x/[1]) on range [" << start << ", "
  //           << stop << "]):\n";
  // std::cout << "  Constant A  = " << A << " ± " << Aerr << "\n";
  // std::cout << "  Tau         = " << tau << " ± " << tauErr << "\n";
}

} // namespace digiAnalysis
