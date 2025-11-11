/*define meantime using tracesSmooth if SMOOTH is turned on*/

#include "WaveForm.h"
#include "Rtypes.h"
#include "TMath.h"
#include "globals.h"
#include "includes.hh"
#include <iostream>
#include <vector>

// ClassImp(digiAnalysis::WaveForm);

using namespace std;

namespace digiAnalysis
{ // namespace digiAnalysis
  /*// UShort_t WaveForm::nSampleBL = 16;
  const char *env_var_baselinesamples = std::getenv("NUM_SAMPLE_BASELINE");
  UShort_t WaveForm::nSampleBL =
      static_cast<UShort_t>(std::stoul(env_var_baselinesamples));

  #ifdef SMOOTH
  const char *envVar = std::getenv("SMOOTH_BOX_SIZE");
  UShort_t WaveForm::smoothBoxSz = static_cast<UShort_t>(std::stoul(envVar));
  #endif
  */

  TVirtualFFT *WaveForm::fft = nullptr;
  TVirtualFFT *WaveForm::ifft = nullptr;
  void WaveForm::InitFFT(int N)
  {
    if (!fft)
    {
      fft = TVirtualFFT::FFT(1, &N, "R2C M K");
    }
    if (!ifft)
    {
      ifft = TVirtualFFT::FFT(1, &N, "C2R M K");
    }
  }

  /*Constructors*/
  WaveForm::WaveForm() // Default Constructor initializes empty vectors
  {
    meantime = 0.0;
    baseline = 0.0;
  }

  WaveForm::WaveForm(TArrayS *arr) // CoMPASS saves waveforms as TArrayS
  {
#ifdef SMOOTH
    double movingSum = 0;
#endif
    meantime = 0;
    baseline = 0;
    SetBaseLine(arr); // Alwasy Set Base Line First
    double sampleSum = 0;
    unsigned int size = arr->GetSize();
    for (unsigned int j = 0; j < size; j++)
    {
      traces.push_back(baseline - arr->At(j));
      if (j >= GateStart and j <= GateStart + GateLenLong)
      {
        meantime = meantime + traces[j] * j;
        sampleSum = sampleSum + traces[j];
      }
#ifdef SMOOTH
      if (smoothBoxSz == 1)
      {
        tracesSmooth.push_back(traces[j]);
      }
      else
      {
        if (j < smoothBoxSz)
        {
          // tracesSmooth.push_back(0);
          movingSum = movingSum + traces[j];
          tracesSmooth.push_back(movingSum / j);
        }
        movingSum = movingSum + traces[j] - traces[j - smoothBoxSz];
        tracesSmooth.push_back(movingSum / smoothBoxSz);
      }
#endif
    }
    // meantime = (GateLenLong * 0.5 + GateStart) - meantime / sampleSum;
    meantime = meantime / sampleSum;
    if (meantime > 0.0)
    {
      meantime = TMath::Log10(meantime);
    }
    else
    {
      meantime = -1.0 * TMath::Log10(fabs(meantime));
    }
    InitFFT(traces.size());
  }

  WaveForm::WaveForm(const std::vector<double> tr) // copy from other vector
  {
    if (!tr.empty())
    {
      SetBaseLine(tr);
      if (std::fabs(baseline - tr[0]) >
          100) // check if the copied traces have baseline subtracted
      {
        SetWaveForm(tr);
      }
      else
      {
        traces = tr;
        SetMeanTime();
#ifdef SMOOTH
        SetSmooth(smoothBoxSz);
#endif
      }
      InitFFT(traces.size());
    }
    else
    {
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
    if (!wf.tracesFFT.empty()) // copy traces CFD if not empty
    {
      tracesFFT = wf.tracesFFT;
    }
    if (wf.fitFunc) // copy traces CFD if not empty
    {
      fitFunc = wf.fitFunc;
    }
    InitFFT(traces.size());
  }

  WaveForm::WaveForm(const WaveForm &wf1, const WaveForm &wf2)
  {
    ConcatenateWaveForms(wf1, wf2);
  }

  WaveForm::WaveForm(UShort_t sizeOfWaveForms,
                     const std::vector<WaveForm> vecOfWaveForm)
  {
    AverageWaveForms(sizeOfWaveForms, vecOfWaveForm);
  }

  /*Destructor*/
  WaveForm::~WaveForm()
  {
    /*clear() Removes elements
    shrink_to_fit() frees unused memory*/

    traces.clear();
    traces.shrink_to_fit();

    tracesSmooth.clear();
    tracesSmooth.shrink_to_fit();

    tracesFFT.clear();
    tracesFFT.shrink_to_fit();

    tracesMovBLCorr.clear();
    tracesMovBLCorr.shrink_to_fit();

    if (fitFunc)
    {
      delete fitFunc;
      fitFunc = nullptr;
    }

    if (WaveForm::fft)
    {
      delete WaveForm::fft;
      WaveForm::fft = nullptr;
    }
    if (WaveForm::ifft)
    {
      delete WaveForm::ifft;
      WaveForm::ifft = nullptr;
    }
  }

  void WaveForm::Plot()
  {
    if (traces.empty() && tracesSmooth.empty())
    {
      std::cout
          << "Both traces and tracesSmooth are empty. No plot will be created."
          << std::endl;
      return;
    }

    // Create a canvas
    TObject *obj = gROOT->FindObject("canvas");
    TCanvas *canvas = dynamic_cast<TCanvas *>(obj);
    if (canvas)
    {
      canvas->cd();    // make it current
      canvas->Clear(); // optional: clear previous plot
      std::cout << "Reusing existing canvas: " << "canvas" << std::endl;
    }
    else
    {
      canvas = new TCanvas("canvas", "WaveForm Plot", 1600, 1000);
      std::cout << "Created new canvas: " << "canvas" << std::endl;
    }
    // TCanvas *canvas = new TCanvas("canvas", "WaveForm Plot", 800, 600);

    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Create a legend
    int nTraces = traces.size();
    int nTracesSmooth = tracesSmooth.size();
    int nTracesMovBLCorr = tracesMovBLCorr.size();

    TGraph *graphTraces = nullptr;
    TGraph *graphTracesSmooth = nullptr;
    TGraph *graphTracesMovBLCorr = nullptr;

    if (!tracesFFT.empty())
    {
      canvas->Divide(2, 1);
      canvas->cd(1);
    }

    double peakHt = 0;

    // Plot traces if not empty
    if (!traces.empty())
    {

      graphTraces = new TGraph(nTraces);
      for (int i = 0; i < nTraces; ++i)
      {
        graphTraces->SetPoint(i, i, traces[i]);
        if (traces[i] > peakHt)
          peakHt = traces[i];
      }
      graphTraces->SetLineColor(kBlue);
      graphTraces->SetLineWidth(2);
      graphTraces->SetTitle("Traces");
      graphTraces->Draw("AL");
      legend->AddEntry(graphTraces, "Traces", "l");
    }

    // Plot tracesSmooth if not empty
    if (!tracesSmooth.empty())
    {
      graphTracesSmooth = new TGraph(nTracesSmooth);
      for (int i = 0; i < nTracesSmooth; ++i)
      {
        graphTracesSmooth->SetPoint(i, i, tracesSmooth[i]);
      }
      graphTracesSmooth->SetLineColor(kRed);
      graphTracesSmooth->SetLineWidth(2);
      graphTracesSmooth->SetTitle("Smoothed Traces");

      if (graphTraces)
      {
        graphTracesSmooth->Draw("L SAME");
      }
      else
      {
        graphTracesSmooth->Draw("AL");
      }
      legend->AddEntry(graphTracesSmooth, "Smoothed Traces", "l");
    }

    // Plot tracesMovBLCorr if not empty
    if (!tracesMovBLCorr.empty())
    {
      graphTracesMovBLCorr = new TGraph(nTracesMovBLCorr);
      for (int i = 0; i < nTracesMovBLCorr; ++i)
      {
        graphTracesMovBLCorr->SetPoint(i, i, tracesMovBLCorr[i]);
      }
      graphTracesMovBLCorr->SetLineColor(kGreen);
      graphTracesMovBLCorr->SetLineWidth(2);
      graphTracesMovBLCorr->SetTitle("BL Corrected Traces");

      if (graphTraces)
      {
        graphTracesMovBLCorr->Draw("L SAME");
      }
      else
      {
        graphTracesMovBLCorr->Draw("AL");
      }
      legend->AddEntry(graphTracesMovBLCorr, "BL Corrected Traces", "l");
    }

    // Draw the Baseline
    TLine *line = new TLine(0.0, 0.0, nTraces, 0.0);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->Draw("LSAME");

    // Draw the Gates
    // LongGate is Red
    // ShortGate is Green
    TPolyLine *shortGate = DrawGate(GateStart, GateStart + GateLenShort, peakHt);
    TPolyLine *longGate = DrawGate(GateStart, GateStart + GateLenLong, peakHt);
    shortGate->SetLineColor(kGreen);
    longGate->SetLineColor(kRed);
    shortGate->Draw("LSAME");
    longGate->Draw("LSAME");

    // Draw Fit if present
    if (fitFunc)
    {
      fitFunc->SetLineColor(kGreen);
      fitFunc->SetLineWidth(2);
      fitFunc->Draw("LSAME");
    }
    // Draw the legend
    legend->Draw();

    if (!tracesFFT.empty())
    {
      canvas->cd(2);
      TGraph *graphTracesFFT = new TGraph(tracesFFT.size());
      for (size_t i = 0; i < tracesFFT.size(); ++i)
      {
        graphTracesFFT->SetPoint(i, i, tracesFFT[i]);
      }
      graphTracesFFT->SetLineColor(kBlack);
      graphTracesFFT->SetLineWidth(2);
      graphTracesFFT->SetTitle("FFT");
      graphTracesFFT->Draw("AL");
    }

    // Update the canvas
    canvas->Update();
  }

  TPolyLine *WaveForm::DrawGate(UShort_t start, UShort_t stop, Double_t height)
  {
    TPolyLine *rect = new TPolyLine(4);

    rect->SetPoint(0, start, 0);
    rect->SetPoint(1, start, height);
    rect->SetPoint(2, stop, height);
    rect->SetPoint(3, stop, 0);

    rect->SetLineWidth(2);

    return rect;
  }

  void WaveForm::Plot(std::vector<double> tr)
  {
    if (tr.empty())
    {
      std::cout << "Both traces is empty. No plot will be created." << std::endl;
      return;
    }

    // Create a canvas
    TObject *obj = gROOT->FindObject("canvas");
    TCanvas *canvas = dynamic_cast<TCanvas *>(obj);
    if (canvas)
    {
      canvas->cd();    // make it current
      canvas->Clear(); // optional: clear previous plot
      std::cout << "Reusing existing canvas: " << "canvas" << std::endl;
    }
    else
    {
      canvas = new TCanvas("canvas", "WaveForm Plot", 1600, 1000);
      std::cout << "Created new canvas: " << "canvas" << std::endl;
    }
    // TCanvas *canvas = new TCanvas("canvas", "WaveForm Plot", 800, 600);

    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Create a legend
    int nTraces = tr.size();

    TGraph *graphTraces = nullptr;

    // Plot traces if not empty
    if (!tr.empty())
    {

      graphTraces = new TGraph(nTraces);
      for (int i = 0; i < nTraces; ++i)
      {
        graphTraces->SetPoint(i, i, tr[i]);
      }
      graphTraces->SetLineColor(kBlue);
      graphTraces->SetLineWidth(2);
      graphTraces->SetTitle("Traces");
      graphTraces->Draw("AL");
      legend->AddEntry(graphTraces, "Traces", "l");
    }
    // Draw the legend
    legend->Draw();

    // Update the canvas
    canvas->Update();
  }

  /*Getters*/
  std::vector<double> WaveForm::GetTraces() { return traces; }
  std::vector<double> WaveForm::GetTracesSmooth() { return tracesSmooth; }
  std::vector<double> WaveForm::GetTracesFFT() { return tracesFFT; }
  double WaveForm::GetMeanTime() { return meantime; }
  double WaveForm::GetBaseLine() { return baseline; }
  UShort_t WaveForm::GetSize() { return traces.size(); }
  bool WaveForm::IsFit()
  {
    if (fitFunc)
    {
      return true;
    }
    else
    {
      return false;
    }
  }
  double WaveForm::GetFitPar(int val) { return fitFunc->GetParameter(val); }
  double WaveForm::GetFitParError(int val) { return fitFunc->GetParError(val); }

  /*Setters*/
  void WaveForm::SetWaveForm(std::vector<double> tr)
  {
#ifdef SMOOTH
    double movingSum = 0;
#endif
    if (!tr.empty())
    {
      SetBaseLine(tr);
      meantime = 0;
      double sampleSum = 0;
      unsigned int size = tr.size();
      for (unsigned int j = 0; j < size; j++)
      {
        traces.push_back(baseline - tr[j]);
        if (j >= GateStart and j <= GateStart + GateLenShort)
        {
          meantime = meantime + traces[j] * j;
          sampleSum = sampleSum + traces[j];
        }
#ifdef SMOOTH
        if (smoothBoxSz == 1)
        {
          tracesSmooth.push_back(traces[j]);
        }
        else
        {
          if (j < smoothBoxSz)
          {
            // tracesSmooth.push_back(0);
            movingSum = movingSum + traces[j];
            tracesSmooth.push_back(movingSum / j);
          }
          movingSum = movingSum + traces[j] - traces[j - smoothBoxSz];
          tracesSmooth.push_back(movingSum / smoothBoxSz);
        }
#endif
      }
      meantime = (meantime / sampleSum);
      if (meantime > 0.0)
      {
        meantime = TMath::Log10(meantime);
      }
      else
      {
        meantime = -1.0 * TMath::Log10(fabs(meantime));
      }
    }
    else
    {
      std::cout << "err SetWaveForm: input vector is empty" << std::endl;
    }
  }

  void WaveForm::SetWaveForm(const WaveForm &wf)
  {
    traces = wf.traces;
    tracesSmooth = wf.tracesSmooth;
    tracesFFT = wf.tracesFFT;

    // Copying the scalar values
    meantime = wf.meantime;
    baseline = wf.baseline;
  }

  void WaveForm::SetTracesMovBLCorr()
  {
    if (!traces.empty())
    {
      if (!IsTracesSmoothSet())
      {
        SetSmooth(smoothBoxSz);
      }
      int N = traces.size();
      int start = smoothBoxSz / 2;
      // Evaluate a baseline from the first few entries
      // It is known that there is no signal here.
      double sum = 0;
      double baselineVal = 0;
      std::vector<double> tracesBL(N, 0.0);
      tracesMovBLCorr.assign(N, 0.0);

      int end = std::min(N, nSampleMovBL + start);
      for (int i = start; i < end; i++)
      {
        sum = sum + traces[i];
      }
      baselineVal = sum / nSampleMovBL;

      // now evaluate the baseline at each point
      int prevcount = 0;
      double blprev = 0, bltemp = 0;
      for (int i = start + nSampleMovBL; i < N; i++)
      {
        tracesBL[i - nSampleMovBL / 2] = baselineVal;
        sum = sum - traces[i - nSampleMovBL];
        sum = sum + traces[i];

        blprev = bltemp;
        bltemp = sum / nSampleMovBL;

        if (prevcount > 0)
        {
          prevcount -= 1;
        }

        if (fabs(bltemp - blprev) > BLError)
        {
          prevcount = nSampleMovBL;
        }
        else if (prevcount == 0)
        {
          baselineVal = bltemp;
        }
      }

      // TVirtualFFT *fft = TVirtualFFT::FFT(1, &N, "R2C");
      fft->SetPoints(tracesBL.data());
      fft->Transform();
      std::vector<double> re(N / 2 + 1), im(N / 2 + 1);
      for (int i = 0; i <= N / 2; i++)
      {
        fft->GetPointComplex(i, re[i], im[i]);
      }

      // TVirtualFFT *ifft = TVirtualFFT::FFT(1, &N, "C2R");
      int cutoff = 10;
      for (int i = 0; i <= N / 2; i++)
      {
        if (i < cutoff)
        {
          TComplex c(re[i], im[i]);
          // std::cout << "re: " << re[i] << " im: " << im[i] << " c: " << c <<
          // std::endl;
          ifft->SetPointComplex(i, c);
        }
        else
        {
          TComplex c(0., 0.);
          ifft->SetPointComplex(i, c);
        }
      }
      ifft->Transform();
      // std::vector<double> tracesBL(N);
      ifft->GetPoints(&tracesBL[0]);
      // std::this_thread::sleep_for(std::chrono::seconds(1));

      for (int i = 0; i < N; i++)
      {
        tracesBL[i] /= N;
      }

      meantime = 0;
      double sampleSum = 0;
      for (int i = 0; i < N; i++)
      {
        // tracesMovBLCorr[i] = tracesSmooth[i] - tracesBL[i];
        tracesMovBLCorr[i] = traces[i] - tracesBL[i];
        if (i >= GateStart and i <= GateStart + GateLenShort)
        {
          meantime = meantime + tracesMovBLCorr[i] * i;
          sampleSum = sampleSum + tracesMovBLCorr[i];
        }
      }
      meantime = (meantime / sampleSum);
      if (meantime > 0.0)
      {
        meantime = TMath::Log10(meantime);
      }
      else
      {
        meantime = -1.0 * TMath::Log10(fabs(meantime));
      }
    }
    else
    {
      std::cout << "err SetTracesMovBLCorr: Fill traces first" << std::endl;
    }

    // delete fft;
    // delete ifft;
  }

  void WaveForm::SetSmooth(UShort_t sBoxSz, std::string kernel)
  {
    // std::cout << "Smoothing with boxSz: " << sBoxSz << std::endl;
    if (!traces.empty())
    {
      if (kernel == "MovA")
      {
        double movingSum = 0;
        SetBaseLine();
        // std::cout << "Baseline set to: " << baseline << std::endl;
        unsigned int size = traces.size();
        // std::cout << "traces size: " << size << std::endl;
        if (sBoxSz == 1)
        {
          tracesSmooth = traces;
        }
        else
        {
          for (unsigned int j = 0; j < size; j++)
          {
            if (j < sBoxSz)
            {
              tracesSmooth.push_back(0);
              movingSum = movingSum + traces[j];
              // tracesSmooth.push_back(movingSum / j);
            }
            else
            {
              movingSum = movingSum + traces[j] - traces[j - sBoxSz];
              tracesSmooth.push_back(movingSum / sBoxSz);
            }
          }
        }
      }
      else if (kernel == "Gauss")
      {
        // Create the Gaussian kernel
        std::vector<double> gauss(sBoxSz);
        double sig = sBoxSz / 10;
        double sumgauss = 0;
        for (int i = 0; i < sBoxSz; i++)
        {
          gauss[i] = TMath::Exp(-0.5 * TMath::Power((i - (sBoxSz / 2)) / sig, 2));
          sumgauss += gauss[i];
        }

        // Smoothing
        for (int i = 0; i < GetSize(); i++)
        {
          double gaussSmooth = 0;
          if (fabs(i - GetSize() / 2) <= (GetSize() - sBoxSz) / 2)
          {
            for (int j = 0; j < sBoxSz; j++)
            {
              gaussSmooth += gauss[j] / sumgauss * traces[i - sBoxSz / 2 + j];
            }
          }
          tracesSmooth.push_back(gaussSmooth);
        }
      }
    }
    else
    {
      std::cout << "err SetSmooth: Fill traces first" << std::endl;
    }
  }

  void WaveForm::SetSmooth()
  {
    if (!traces.empty())
    {
      double movingSum = 0;
      SetBaseLine();
      unsigned int size = traces.size();
      if (smoothBoxSz == 1)
      {
        tracesSmooth = traces;
      }
      else
      {
        for (unsigned int j = 0; j < size; j++)
        {
          if (j < smoothBoxSz)
          {
            tracesSmooth.push_back(0);
            movingSum = movingSum + traces[j];
            // tracesSmooth.push_back(movingSum / j);
          }
          movingSum = movingSum + traces[j] - traces[j - smoothBoxSz];
          tracesSmooth.push_back(movingSum / smoothBoxSz);
        }
      }
    }
    else
    {
      std::cout << "err SetSmooth: Fill traces first" << std::endl;
    }
  }

  void WaveForm::SetCFD() {}

  void WaveForm::SetMeanTime()
  {
    meantime = 0;
    double sampleSum = 0;
    if (!tracesMovBLCorr.empty())
    {
      unsigned int size = tracesMovBLCorr.size();
      for (unsigned int j = GateStart; j < GateLenLong + GateStart; j++)
      {
        meantime = meantime + tracesMovBLCorr[j] * j;
        sampleSum = sampleSum + tracesMovBLCorr[j];
      }
      meantime = meantime / sampleSum;
      if (meantime > 0.0)
      {
        meantime = TMath::Log10(meantime);
      }
      else
      {
        meantime = -1.0 * TMath::Log10(fabs(meantime));
      }
    }
    else if (!tracesSmooth.empty())
    {
      unsigned int size = traces.size();
      for (unsigned int j = GateStart; j < GateLenLong + GateStart; j++)
      {
        meantime = meantime + tracesSmooth[j] * j;
        sampleSum = sampleSum + tracesSmooth[j];
      }
      meantime = meantime / sampleSum;
      if (meantime > 0.0)
      {
        meantime = TMath::Log10(meantime);
      }
      else
      {
        meantime = -1.0 * TMath::Log10(fabs(meantime));
      }
    }
    else if (!traces.empty())
    {
      unsigned int size = traces.size();
      for (unsigned int j = GateStart; j < GateLenLong + GateStart; j++)
      {
        meantime = meantime + traces[j] * j;
        sampleSum = sampleSum + traces[j];
      }
      meantime = meantime / sampleSum;
      if (meantime > 0.0)
      {
        meantime = TMath::Log10(meantime);
      }
      else
      {
        meantime = -1.0 * TMath::Log10(fabs(meantime));
      }
    }
    else
    {
      std::cout << "err SetMeanTime: traces is empty" << std::endl;
    }
  }

  void WaveForm::SetMeanTime(const std::vector<double> tr)
  {
    if (!tr.empty())
    {
      SetMeanTime(tr, GateStart, GateLenLong + GateStart);
    }
    else
    {
      std::cout << "err SetMeanTime: input vector is empty" << std::endl;
    }
  }

  void WaveForm::SetMeanTime(const std::vector<double> tr, UShort_t start,
                             UShort_t stop)
  {
    meantime = 0;
    double sampleSum = 0;
    if (!tr.empty())
    {
      if (start < stop)
      {
        if (tr.size() < stop)
        {
          std::cout << "Warning: stop is less than input vector size. Evaluating "
                       "till end"
                    << std::endl;
          stop = tr.size();
        }
        for (unsigned int j = start; j < stop; j++)
        {
          meantime = meantime + tr[j] * j;
          sampleSum = sampleSum + tr[j];
        }
        meantime = meantime / sampleSum;
        if (meantime > 0.0)
        {
          meantime = TMath::Log10(meantime);
        }
        else
        {
          meantime = -1.0 * TMath::Log10(fabs(meantime));
        }
      }
      else
      {
        std::cout << "error: start > stop, exiting witout meantime setting"
                  << std::endl;
      }
    }
    else
    {
      std::cout << "err SetMeanTime: input vector is empty" << std::endl;
    }
  }

  void WaveForm::SetBaseLine()
  {
    baseline = 0;
    double sum = 0;

    if (!traces.empty())
    {
      for (unsigned int j = blStart; j < nSampleBL + blStart; j++)
      {
        sum = sum + traces[j];
      }
      baseline = sum / nSampleBL;
    }
    else
    {
      std::cout << "err SetBaseLine: traces is empty" << std::endl;
    }
  }

  void WaveForm::SetBaseLine(std::vector<double> tr)
  {
    baseline = 0;
    double sum = 0;

    if (!tr.empty())
    {
      for (unsigned int j = blStart; j < nSampleBL + blStart; j++)
      {
        sum = sum + tr[j];
      }
      baseline = sum / nSampleBL;
    }
    else
    {
      std::cout << "err SetBaseLine: input vector is empty" << std::endl;
    }
  }

  void WaveForm::SetBaseLine(TArrayS *arr)
  {
    baseline = 0;
    double sum = 0;

    if (arr && (arr->GetSize() > nSampleBL))
    {
      for (unsigned int j = blStart; j < nSampleBL + blStart; j++)
      {
        sum = sum + arr->At(j);
      }
      baseline = sum / nSampleBL;
    }
    else
    {
      std::cout << "err SetBaseLine: input array is empty or smaller than "
                   "NUM_SAMPLE_BASELINE"
                << std::endl;
    }
  }

  void WaveForm::SetTracesFFT()
  {
    if (!traces.empty())
    {
      Int_t n = static_cast<int>(traces.size());
      fft = TVirtualFFT::FFT(1, &n, "R2C");
      // TVirtualFFT *fft1 = TVirtualFFT::FFT(1, &n, "R2C");
      fft->SetPoints(traces.data());
      fft->Transform();
      int nFreq = traces.size() / 2 + 1; // only positive frequencies
      for (int i = 0; i < nFreq; i++)
      {
        double re, im;
        fft->GetPointComplex(i, re, im);
        tracesFFT.push_back(std::sqrt(re * re + im * im));
      }
    }
    else
    {
      std::cout << "err SetTracesFFT: traces is empty" << std::endl;
    }
  }

  void WaveForm::SetTracesFFT(std::string whichTrace)
  {
    if (whichTrace == "orig")
    {
      SetTracesFFT();
    }
    else if (whichTrace == "smooth")
    {
      if (!tracesSmooth.empty())
      {
        Int_t n = static_cast<int>(tracesSmooth.size());
        fft = TVirtualFFT::FFT(1, &n, "R2C");
        // TVirtualFFT *fft1 = TVirtualFFT::FFT(1, &n, "R2C");
        fft->SetPoints(tracesSmooth.data());
        fft->Transform();

        int nFreq = tracesSmooth.size() / 2 + 1; // only positive frequencies

        for (int i = 0; i < nFreq; i++)
        {
          double re, im;
          fft->GetPointComplex(i, re, im);
          tracesFFT.push_back(std::sqrt(re * re + im * im));
        }
      }
      else
      {
        std::cout << "err SetTracesFFT: traces is empty" << std::endl;
      }
    }
    else
    {
      std::cout << "err SetTracesFFT: Pass proper value for whichTrace"
                << std::endl;
      std::cout << "Valid options are: orig, smooth" << std::endl;
    }
  }

  void WaveForm::SetTracesFFT(std::vector<double> trFFT)
  {
    if (!trFFT.empty())
    {
      tracesFFT.clear();
      Int_t n = static_cast<int>(trFFT.size());
      // TVirtualFFT *fft1 = TVirtualFFT::FFT(1, &n, "R2C");
      fft->SetPoints(trFFT.data());
      fft->Transform();

      int nFreq = trFFT.size() / 2 + 1; // only positive frequencies

      for (int i = 0; i < nFreq; i++)
      {
        double re, im;
        fft->GetPointComplex(i, re, im);
        tracesFFT.push_back(std::sqrt(re * re + im * im));
      }
    }
    else
    {
      std::cout << "err SetTracesFFT: pass non empty vector" << std::endl;
    }
  }

  std::vector<double> WaveForm::EvalTracesFFT(std::vector<double> trFFT)
  {
    std::vector<double> res;
    if (!trFFT.empty())
    {
      Int_t n = static_cast<int>(trFFT.size());
      // TVirtualFFT *fft1 = TVirtualFFT::FFT(1, &n, "R2C");
      fft = TVirtualFFT::FFT(1, &n, "R2C");
      fft->SetPoints(trFFT.data());
      fft->Transform();

      int nFreq = trFFT.size() / 2 + 1; // only positive frequencies

      for (int i = 0; i < nFreq; i++)
      {
        double re, im;
        fft->GetPointComplex(i, re, im);
        res.push_back(std::sqrt(re * re + im * im));
      }
    }
    else
    {
      std::cout << "err SetTracesFFT: pass non empty vector" << std::endl;
    }
    return res;
  }

  // std::vector<double> WaveForm::EvalDecayTime(UShort_t FitStart, UShort_t
  // FitEnd, UShort_t numDecayConst) {}

  void WaveForm::ShiftWaveForm(int BL)
  {
    if (!traces.empty())
    {
      for (unsigned int j = 0; j < traces.size(); j++)
      {
        traces[j] = traces[j] + BL;
      }
    }
  }

  double WaveForm::IntegrateWaveForm()
  {
    double sum = 0;
    if (!traces.empty())
    {
      sum = IntegrateWaveForm(0, traces.size());
    }
    else
    {
      std::cout << "err IntegrateWaveForm: fill traces first" << std::endl;
    }
    return sum;
  }

  double WaveForm::IntegrateWaveForm(int startTime, int stopTime)
  {
    double sum = 0;
    if (!traces.empty() && (traces.size() > stopTime))
    {
      if (startTime < stopTime)
      {
        for (unsigned int j = startTime; j < stopTime; j++)
        {
          sum = sum + traces[j];
        }
      }
      else
      {
        std::cout << "err IntegrateWaveForm: order the times properly "
                  << std::endl;
      }
    }
    else
    {
      std::cout << "err IntegrateWaveForm: fill traces first: " << traces.size()
                << " : " << stopTime << std::endl;
    }
    return sum;
  }

  double WaveForm::IntegrateSmoothWaveForm(int startTime, int stopTime)
  {
    double sum = 0;
    if (!tracesSmooth.empty() && (tracesSmooth.size() > stopTime))
    {
      if (startTime < stopTime)
      {
        for (unsigned int j = startTime; j < stopTime; j++)
        {
          sum = sum + tracesSmooth[j];
        }
      }
      else
      {
        std::cout << "err IntegrateSmoothWaveForm: order the times properly "
                  << std::endl;
      }
    }
    else
    {
      std::cout << "err IntegrateSmoothWaveForm: fill traces Smooth first"
                << std::endl;
    }
    return sum;
  }

  double WaveForm::IntegrateBLCorrWaveForm(int startTime, int stopTime)
  {
    double sum = 0;
    if (!tracesMovBLCorr.empty() && (tracesMovBLCorr.size() > stopTime))
    {
      if (startTime < stopTime)
      {
        for (unsigned int j = startTime; j < stopTime; j++)
        {
          sum = sum + tracesMovBLCorr[j];
        }
      }
      else
      {
        std::cout << "err IntegrateBLCorrWaveForm: order the times properly "
                  << std::endl;
      }
    }
    else
    {
      std::cout << "err IntegrateBLCorrWaveForm: fill traces Smooth first"
                << std::endl;
    }
    return sum;
  }

  void WaveForm::ConcatenateWaveForms(const WaveForm &wf1, const WaveForm &wf2)
  {
    if (!wf1.traces.empty()) // copy traces if not empty
    {
      traces = wf1.traces;
    }
    if (!wf1.tracesSmooth.empty()) // copy traces smooth if not empty
    {
      tracesSmooth = wf1.tracesSmooth;
    }
    if (!wf1.tracesFFT.empty()) // copy traces CFD if not empty
    {
      tracesFFT = wf1.tracesFFT;
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
    if (!wf2.tracesFFT.empty()) // concatenate traces CFD if not empty
    {
      tracesFFT.insert(tracesFFT.end(), wf2.tracesFFT.begin(),
                       wf2.tracesFFT.end());
    }

    meantime =
        traces.size() +
        1000; // can use this to check if the waveform is single or concatenated
    baseline = wf1.baseline;
  }

  void WaveForm::AverageWaveForms(UShort_t sizeOfWaveForms,
                                  const std::vector<WaveForm> vecOfWaveForm)
  {
    // currently only traces is averaged
    UShort_t numWaveForm = vecOfWaveForm.size();
    double sum = 0;
    baseline = 0;

    for (unsigned int i = blStart; i < nSampleBL + blStart; i++)
    {
      for (unsigned int j = 0; j < numWaveForm; j++)
      {
        sum = sum + vecOfWaveForm[j].traces[i] / numWaveForm;
      }
    }
    baseline = sum / nSampleBL;

    for (unsigned int i = 0; i < sizeOfWaveForms; i++)
    {
      sum = 0;
      for (unsigned int j = 0; j < numWaveForm; j++)
      {
        sum = sum + vecOfWaveForm[j].traces[i];
      }
      if (baseline > 10)
      {
        sum = baseline - sum / numWaveForm;
      }
      else
      {
        sum = sum / numWaveForm;
      }
      traces.push_back(sum);
    }
    SetMeanTime();
  }

  void WaveForm::AverageWaveForms(ULong_t start, UShort_t numWaveForm,
                                  UShort_t sizeOfWaveForms,
                                  const std::vector<WaveForm> vecOfWaveForm)
  {
    // currently only traces is averaged
    double sum = 0;
    baseline = 0;

    for (unsigned int i = blStart; i < nSampleBL + blStart; i++)
    {
      for (unsigned int j = start; j < start + numWaveForm; j++)
      {
        sum = sum + vecOfWaveForm[j].traces[i] / numWaveForm;
      }
    }
    baseline = sum / nSampleBL;

    for (unsigned int i = 0; i < sizeOfWaveForms; i++)
    {
      sum = 0;
      for (unsigned int j = start; j < start + numWaveForm; j++)
      {
        sum = sum + vecOfWaveForm[j].traces[i];
      }
      if (baseline > 10)
      {
        sum = baseline - sum / numWaveForm;
      }
      else
      {
        sum = sum / numWaveForm;
      }
      traces.push_back(sum);
    }
    SetMeanTime();
  }

  void WaveForm::ScaleWaveForm(double Scale)
  {
    std::transform(traces.begin(), traces.end(), traces.begin(),
                   [Scale](double val)
                   { return val * Scale; });
    if (!tracesSmooth.empty())
    {
      std::transform(tracesSmooth.begin(), tracesSmooth.end(),
                     tracesSmooth.begin(),
                     [Scale](double val)
                     { return val * Scale; });
    }
  }

  void WaveForm::AddWaveForm(const WaveForm &wf1)
  {
    if (traces.size() == wf1.traces.size())
    {
      // Add each element from traces and wf1.traces using std::transform
      std::transform(traces.begin(), traces.end(), wf1.traces.begin(),
                     traces.begin(), std::plus<double>());
    }
    else
    {
      std::cout << "traces must be the same size for adding! Nothing done!!!!!!!!"
                << std::endl;
    }
    if (tracesSmooth.size() == wf1.tracesSmooth.size())
    {
      std::transform(tracesSmooth.begin(), tracesSmooth.end(),
                     wf1.tracesSmooth.begin(), tracesSmooth.begin(),
                     std::plus<double>());
    }
    else
    {
      std::cout << "tracesSmooth not of same size hence not adding!" << std::endl;
    }
  }

  std::vector<std::unique_ptr<WaveForm>>
  WaveForm::SplitWaveForm(UShort_t numSplits)
  {

    std::vector<std::unique_ptr<WaveForm>> splitWaveForms;

    int traceSize = traces.size();
    int splitSize = traceSize / numSplits;

    // Handle cases where the number of splits is greater than the size
    if (splitSize == 0)
    {
      // If numSplits is too large, return empty vector
      return splitWaveForms;
    }

    // Loop to create each split waveform
    for (int i = 0; i < numSplits; ++i)
    {
      // Create a new unique_ptr for each new WaveForm object
      auto newWaveForm = std::make_unique<WaveForm>();

      // Define the range for the current split
      int start = i * splitSize;
      int end = (i == numSplits - 1) ? traceSize : (i + 1) * splitSize;

      // Split the traces and assign to the new waveform
      if (!traces.empty())
      {
        newWaveForm->traces.assign(traces.begin() + start, traces.begin() + end);
      }
      else
      {
        std::cout << "err SplitWaveForm:  Empty traces passed for Splitting"
                  << std::endl;
        return splitWaveForms;
      }
      if (!tracesSmooth.empty())
      {
        newWaveForm->tracesSmooth.assign(tracesSmooth.begin() + start,
                                         tracesSmooth.begin() + end);
      }
      if (!tracesFFT.empty())
      {
        newWaveForm->tracesFFT.assign(tracesFFT.begin() + start,
                                      tracesFFT.begin() + end);
      }

      newWaveForm->SetBaseLine();
      newWaveForm->SetMeanTime();

      // Add the new waveform to the vector
      splitWaveForms.push_back(std::move(newWaveForm));
    }

    return splitWaveForms;
  }

  void WaveForm::FitExponential(int start, int stop)
  {
    // Check valid range
    if (start < 0 || stop >= GetSize() || start >= stop)
    {
      std::cerr << "Invalid fit range: [" << start << ", " << stop << "]\n";
    }

    int nPoints = stop - start + 1;
    auto graph = std::make_unique<TGraph>(nPoints);

    for (int i = 0; i < nPoints; ++i)
    {
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
