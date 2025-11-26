#include "Analysis.h"
#include "TVirtualFFT.h"
#include "WaveForm.h"
#include "globals.h"
#include "singleHits.h"
#include <TApplication.h>
#include <TComplex.h>
#include <TH1.h>
#include <TMath.h>
#include <cmath>
#include <iostream>
int main(int argc, char *argv[]) {
#ifdef WAVES
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;

  // std::string fname = "/home/kirtikesh/analysisSSD/DATA/WCu_Test/"
  //                     "run_CsBare_NaI0_1350V_0pt5Vpp_21Nov/"
  //                     "FILTERED/DataF_run_CsBare_NaI0_1350V_0pt5Vpp_21Nov.root";

  std::string fname = "/home/kirtikesh/analysisSSD/DATA/SPE/"
                      "run_CsBare_NaI0_1450V_INTCFD_26Nov/FILTERED/"
                      "DataF_run_CsBare_NaI0_1450V_INTCFD_26Nov.root";

  digiAnalysis::Analysis an(fname, 0, 10000, 0);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int nentries = hitsVector.size();
  std::cout << "Got the vector from an: " << nentries << std::endl;
  int evi = 0;
  std::string userInput;
  bool keepGoing = true;
  digiAnalysis::WaveForm *WF = nullptr;
  std::vector<double> tracesSmooth;
  std::vector<double> traces;
  std::vector<double> tracesFFTRe;
  std::vector<double> tracesFFTIm;
  int sbSize = 30;

  int peakDist = 200; // peak seperation in ns
  peakDist /= 2;      // to convert to samples in 500MSPS digitizer

  double peakThreshold = 5.0; // threshold for integration peak
  double threshold = 5.0;     // threshold for peak detection above baseline

  UShort_t checkStartTime =
      digiAnalysis::GateStart +
      digiAnalysis::GateLenLong; // look for SPE after gate end
  UShort_t checkStopTime = 4000;
  int integrateRange = 57;

  // while (keepGoing) {
  //   if (fabs(hitsVector[evi]->GetEnergy() - 75) < 10) {
  //     hitsVector[evi]->SetSmoothWF(sbSize);
  //     WF = nullptr;
  //     WF = hitsVector[evi]->GetWFPtr();
  //     WF->SetTracesFFT();
  //     hitsVector[evi]->Print();
  //     WF->Plot();

  //     std::cout << "Do you want to see the next waveform? (y/n): ";
  //     std::getline(std::cin, userInput);
  //     if (userInput != "y" && userInput != "Y") {
  //       keepGoing = false;
  //     }
  //   }
  //   evi += 1;
  // }

  WF = hitsVector[0]->GetWFPtr();
  int kernelSz = WF->GetSize();
  TVirtualFFT *fft = TVirtualFFT::FFT(1, &kernelSz, "R2C M K");
  TVirtualFFT *ifft = TVirtualFFT::FFT(1, &kernelSz, "C2R M K");
  std::vector<double> tracesifft(kernelSz);
  TComplex c;
  double re, im;
  int target = 1124;
  int halfwidth = 40;
  double phase;
  double new_re;
  double new_im;
  UShort_t filterStart = 1000;
  double filterSig = 300;
  double w;

  bool isSPE = false;
  TH1 *hSPE = new TH1F("hSPE", "SPE Energy", 16384, 0, 16384);
  TH2 *hPEH =
      new TH2F("hPEH", "PE Height vs Energy", 4000, 0, 4000, 250, 0, 250);
  for (evi = 0; evi < nentries and keepGoing; evi++) {
    isSPE = false;
    if (evi % 100000 == 0)
      std::cout << evi << std::endl;
    if (hitsVector[evi]->GetEnergy() > 90) {
      WF = nullptr;
      WF = hitsVector[evi]->GetWFPtr();
      // WF->SetTracesFFT();
      // WF->SetSmooth(sbSize);
      // tracesFFT = WF->GetTracesFFT();
      traces = {};
      traces = WF->GetTraces();
      fft->SetPoints(traces.data());
      fft->Transform();
      // fft at target is set to the average in the vicinity of target
      // double sum_mag = 0.0;
      // int count = 0;
      // for (int k = target - halfwidth; k <= target + halfwidth; k++) {
      //   if (k == target)
      //     continue; // skip the bin we want to replace
      //   if (k < 0 || k > kernelSz / 2)
      //     continue;
      //   fft->GetPointComplex(k, re, im);
      //   double mag = sqrt(re * re + im * im);
      //   sum_mag += mag;
      //   count++;
      // }
      // double avg_mag = sum_mag / count;
      // fft->GetPointComplex(target, re, im);
      // phase = atan2(im, re);
      // new_re = avg_mag * cos(phase);
      // new_im = avg_mag * sin(phase);
      // std::cout << "re: " << re << " : " << new_re << std::endl;
      // std::cout << "im: " << im << " : " << new_im << std::endl;

      // here a filter is applied to the fft
      for (int fftiter = 0; fftiter < kernelSz / 2; fftiter++) {
        fft->GetPointComplex(fftiter, re, im);
        if (fftiter < filterStart)
          c = TComplex(re, im);
        else {
          w = TMath::Exp(-0.5 *
                         TMath::Power((fftiter - filterStart) / filterSig, 2));
          new_re = re * w;
          new_im = im * w;
          c = TComplex(new_re, new_im);
        }
        // if (fftiter != target)
        //   c = TComplex(re, im);
        // else
        //   c = TComplex(new_re, new_im);
        ifft->SetPointComplex(fftiter, c);
      }
      ifft->Transform();
      ifft->GetPoints(&tracesifft[0]);
      for (int i = 0; i < kernelSz; i++) {
        tracesifft[i] /= (-1.0 * kernelSz); // normalization
      }
      hitsVector[evi]->SetWF(tracesifft);
      hitsVector[evi]->SetSmoothWF(sbSize);
      auto pvVec = hitsVector[evi]->DetectPeakValleys(threshold);
      WF = nullptr;
      tracesSmooth = {};
      WF = hitsVector[evi]->GetWFPtr();
      std::vector<double> sub(tracesifft.begin() + digiAnalysis::GateStart +
                                  digiAnalysis::GateLenLong + 1,
                              tracesifft.end());
      int subSz = sub.size();
      TVirtualFFT *fft1 = TVirtualFFT::FFT(1, &subSz, "R2C M K");
      TVirtualFFT *ifft1 = TVirtualFFT::FFT(1, &subSz, "C2R M K");
      fft1->SetPoints(sub.data());
      fft1->Transform();
      for (int fft1iter = 0; fft1iter < subSz / 2; fft1iter++) {
        fft1->GetPointComplex(fft1iter, re, im);
        if (fft1iter <= 5)
          c = TComplex(0, 0);
        else
          c = TComplex(re, im);
        ifft1->SetPointComplex(fft1iter, c);
      }
      ifft1->Transform();
      ifft1->GetPoints(&tracesifft[0]);
      // WF->SetTracesFFT(sub);
      // WF->SetTracesFFT();
      tracesSmooth = WF->GetTracesSmooth();
      auto peaks = pvVec.first;
      if (peaks.size() == 0) {
        continue;
      }
      int iter = 0;
      int prevPeakPos = 0;
      int pos = peaks[iter];
      double SPEInt = 0;
      while (iter < (int(peaks.size()) - 2)) {
        iter += 1;
        prevPeakPos = pos;
        pos = peaks[iter];
        if (pos > checkStartTime and tracesSmooth[pos] > peakThreshold and
            pos < checkStopTime) {
          if (pos - prevPeakPos > peakDist and
              peaks[iter + 1] - pos > peakDist) {
            // std::cout << evi << " : " << pos - integrateRange << " - "
            //           << std::min(pos + 6 * integrateRange,
            //                       int(tracesSmooth.size()) - 1)
            //           << std::endl;
            SPEInt = WF->IntegrateSmoothWaveForm(
                pos - integrateRange,
                std::min(pos + integrateRange, int(tracesSmooth.size()) - 1));
            hSPE->Fill(SPEInt);
            hPEH->Fill(SPEInt, tracesSmooth[pos]);

            // if (fabs(SPEInt - 300) < 100 and tracesSmooth[pos] < 10) {
            //   std::cout << evi << " Position: " << pos
            //             << " IntegratedVal: " << SPEInt << std::endl;
            //   isSPE = true;
            // }
          }
        }
      }
      if (hitsVector[evi]->GetEnergy() < 150) {
        std::cout << evi << std::endl;
        // WF->Plot();
        WF->Plot(tracesifft);
        std::cout << "Do you want to see the next waveform? (y/n): ";
        std::getline(std::cin, userInput);
        if (userInput != "y" && userInput != "Y") {
          keepGoing = false;
        }
      }
      // if (isSPE) {
      //   hitsVector[evi]->Print();
      //   WF->Plot();
      //   std::cout << "Do you want to see the next waveform? (y/n): ";
      //   std::getline(std::cin, userInput);
      //   if (userInput != "y" && userInput != "Y") {
      //     keepGoing = false;
      //   }
      // }
    }
  }

  // TH1 *hNPE = new TH1F("hNPE", "NPE", 1000, 0, 100);
  // double evalE = 0;
  // double gamE = 662;
  // double SPEEnergy = 140;
  // for (evi = 0; evi < nentries; evi++) {
  //   if (evi % 100000 == 0)
  //     std::cout << evi << std::endl;
  //   if (fabs(hitsVector[evi]->GetEnergy() - 1350.0) < 10) {
  //     hitsVector[evi]->SetSmoothWF(sbSize);
  //     hitsVector[evi]->SetEvalEnergy();
  //     evalE = hitsVector[evi]->GetEvalEnergy() * digiAnalysis::GateLenLong
  //     /
  //             digiAnalysis::EvalNormFactor;
  //     // std::cout << evalE / 70.0 / gamE << std::endl;
  //     hNPE->Fill(evalE / SPEEnergy / gamE);
  //   }
  // }

  TCanvas *canvas1 = new TCanvas("canvas1", "Energy Hists", 1600, 1000);
  canvas1->cd();
  hSPE->Draw("HIST");
  // hNPE->Draw("HIST");
  canvas1->Update();

  TCanvas *canvas2 = new TCanvas("canvas2", "Peak Height", 1600, 1000);
  canvas2->cd();
  hPEH->Draw("COLZ");
  canvas2->Update();

  fApp->Run();
#endif
  return 0;
}
