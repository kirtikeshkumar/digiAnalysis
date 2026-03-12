#include "Analysis.h"
#include "WaveForm.h"
#include "includes.hh"
#include "singleHits.h"
#include <RtypesCore.h>
#include <TApplication.h>
#include <TMath.h>
#include <iostream>
#include <vector>

double readUserInput(double defaultValue) {
  std::string input;
  while (true) {
    std::getline(std::cin, input);

    // 5. If no value entered, use default
    if (input.empty()) {
      return defaultValue;
    }

    try {
      size_t pos;
      double value = std::stod(input, &pos);

      // 2. Check if input is purely a double (no extra characters)
      if (pos == input.size()) {
        return value; // 4. Read value into variable
      } else {
        throw std::invalid_argument("Extra characters");
      }
    } catch (const std::exception &) {
      // 3. Ask again if not a valid double
      std::cout
          << "Invalid input. Please enter a valid floating-point number.\n";
    }
  }
}

int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;
  std::string fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
                      "Bkg_Ch0ANDAny_CorrectGeometry_WAVES_Feb25_NaI_1234_"
                      "Square_NaI1_Gain16/FILTERED/"
                      "SDataF_Bkg_Ch0ANDAny_CorrectGeometry_WAVES_Feb25_NaI_"
                      "1234_Square_NaI1_Gain16.root";
  UShort_t channel = 0;
  digiAnalysis::Analysis an(channel, fname, 0, 0, 0);
  std::cout << "getting the vector from an" << std::endl;

  // test Getting
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  std::cout << "got the vector from an" << hitsVector.size() << std::endl;
  int nentries = hitsVector.size();

  // Read the baseline waveform

  std::string fnameBL = "baseline.root";
  TFile *fBL = new TFile(fnameBL.c_str(), "READ");
  TTree *tree = (TTree *)fBL->Get("tree");
  std::vector<double> *Baseline = nullptr;
  tree->SetBranchAddress("Baseline", &Baseline);
  //   std::cout << "reading the baseline file: " << tree->GetEntries() <<
  //   std::endl;
  int nbytes = tree->GetEntry(tree->GetEntries() - 1);
  if (nbytes <= 0) {
    std::cerr << "Failed to read entry! nbytes=" << nbytes << std::endl;
    return 1;
  }
  if (!Baseline || Baseline->empty()) {
    std::cerr << "Vector is null or empty after read!" << std::endl;
    return 1;
  }

  Baseline->assign(Baseline->begin() + 75, Baseline->end() - 75);
  int BaselineSize = Baseline->size();
  double slope =
      (Baseline->at(BaselineSize - 1) - Baseline->at(0)) / (BaselineSize);
  double con = Baseline->at(0);
  for (int iter = 0; iter < BaselineSize; iter++) {
    Baseline->at(iter) = Baseline->at(iter) - (slope * (iter) + con);
  }
  // Evaluate FFT of Baseline
  TVirtualFFT *fft_base = TVirtualFFT::FFT(1, &BaselineSize, "R2C");
  TVirtualFFT *fft_corr = TVirtualFFT::FFT(1, &BaselineSize, "C2R M K");
  fft_base->SetPoints(Baseline->data());
  fft_base->Transform();
  std::vector<double> re_BL, im_BL;
  int BaselineFFTSize = BaselineSize / 2 + 1;
  re_BL.resize(BaselineFFTSize);
  im_BL.resize(BaselineFFTSize);
  double *re_BLptr = re_BL.data();
  double *im_BLptr = im_BL.data();
  fft_base->GetPointsComplex(re_BLptr, im_BLptr);
  digiAnalysis::WaveForm *WFFFT = new digiAnalysis::WaveForm();
  //   WFFFT->Plot(re_BL);

  //   std::cout << "Set the fft of baseline" << std::endl;

  std::vector<double> re_ccor, im_ccor;
  re_ccor.resize(BaselineSize);
  im_ccor.resize(BaselineSize);
  double *re_ccorptr = re_ccor.data();
  double *im_ccorptr = im_ccor.data();

#ifdef WAVES

  // To PLot Based on Cuts
  TH2 *hPSDPlot =
      new TH2F("PSDPlot", "Energy vs PSD", 16384, 0, 16384, 4096, 0, 1);
  for (const auto &hit : hitsVector) {
    if (hit->GetChNum() == channel) {
      hPSDPlot->Fill(hit->GetEnergy(), hit->GetPSD());
    }
  }
  TCanvas *c2 = new TCanvas("c2", "Energy vs PSD", 800, 600);
  hPSDPlot->Draw("COLZ");
  c2->Update();

  int ECutMin = 0, ECutMax = 16384;
  double PSDCutMin = 0, PSDCutMax = 1.0;
  std::string userInput;
  std::cout << "\n ECutMin: ";
  ECutMin = readUserInput(ECutMin);
  std::cout << "\n ECutMax: ";
  ECutMax = readUserInput(ECutMax);
  std::cout << "\n PSDCutMin: ";
  PSDCutMin = readUserInput(PSDCutMin);
  std::cout << "\n PSDCutMax: ";
  PSDCutMax = readUserInput(PSDCutMax);

  int count = 0;
  digiAnalysis::WaveForm *wfptr = nullptr;
  std::vector<digiAnalysis::WaveForm> waveformVector;
  bool keepGoing = true;
  double re_sig, im_sig;
  for (int i = 0; i < nentries && keepGoing; i++) {
    if (hitsVector[i]->GetChNum() == channel &&
        hitsVector[i]->GetPSD() >= PSDCutMin &&
        hitsVector[i]->GetPSD() <= PSDCutMax &&
        hitsVector[i]->GetEnergy() >= ECutMin &&
        hitsVector[i]->GetEnergy() <= ECutMax) {
      if (count < 50) {
        count++;
        wfptr = nullptr;
        wfptr = hitsVector[i]->GetWFPtr();
        // if (wfptr) {
        //   waveformVector.push_back(*wfptr);
        // }
        hitsVector[i]->Print();
        wfptr->SetSmooth(150);
        std::vector<double> Sig = wfptr->GetTracesSmooth();
        Sig.assign(Sig.begin() + 75, Sig.end() - 75);
        double Sslope = (Sig[BaselineSize - 1] - Sig[0]) / (BaselineSize);
        double Scon = Sig[0];
        for (int iter = 0; iter < BaselineSize; iter++) {
          Sig[iter] = Sig[iter] - (Sslope * (iter) + Scon);
        }
        wfptr->SetTracesFFT(Sig);
        std::vector<double> fft_Sig = wfptr->GetTracesFFT();
        std::vector<double> fft_SigPhase = wfptr->GetTracesFFTPhase();
        int SignalFFTSize = fft_Sig.size();
        if (SignalFFTSize == BaselineFFTSize) {
          // Convert FFT of signal to Real and Imaginary parts and do cross
          // correlation
          std::vector<double> re_corr, im_corr;
          for (int i = 0; i < SignalFFTSize; i++) {
            re_sig = fft_Sig[i] * TMath::Cos(fft_SigPhase[i]);
            im_sig = fft_Sig[i] * TMath::Sin(fft_SigPhase[i]);
            re_corr.push_back(re_sig * re_BL[i] + im_sig * im_BL[i]);
            im_corr.push_back(im_sig * re_BL[i] - re_sig * im_BL[i]);
          }
          // IFFT to recover the convoluted signal
          fft_corr->SetPointsComplex(re_corr.data(), im_corr.data());
          fft_corr->Transform();
          fft_corr->GetPointsComplex(re_ccorptr, im_ccorptr);
          int size = re_ccor.size();
          std::cout << "Size of cross correlation: " << size << std::endl;
          std::transform(re_ccor.begin(), re_ccor.end(), re_ccor.begin(),
                         [BaselineSize](double x) { return x / BaselineSize; });
          std::transform(im_ccor.begin(), im_ccor.end(), im_ccor.begin(),
                         [BaselineSize](double x) { return x / BaselineSize; });
          auto max_it = std::max_element(re_ccor.begin(), re_ccor.end());
          int max_lag = std::distance(re_ccor.begin(), max_it);
          if (max_lag > BaselineSize / 2) {
            max_lag = max_lag - BaselineSize;
          }
          std::cout << "Baseline vs Signal Lag: " << max_lag << std::endl;
          double yshift = (*Baseline)[BaselineSize - 75 - max_lag];
          for (int iter = 0; iter < Sig.size(); iter++) {
            int shift = iter + max_lag;
            if (shift > BaselineSize) {
              shift = shift - BaselineSize;
            }
            Sig[iter] = Sig[iter] - (*Baseline)[shift] - yshift;
          }

          wfptr->Plot(Sig, *Baseline);

        } else {
          std::cout << "SignalFFT (" << SignalFFTSize
                    << " events) and BaselineFFT (" << BaselineFFTSize
                    << " events) sizes are different. No "
                       "Correlation based subtraction."
                    << std::endl;
        }

        // wfptr->Plot();
        std::cout << "Do you want to see the next waveform? (y/n): ";
        std::getline(std::cin, userInput);
        if (userInput != "y" && userInput != "Y") {
          keepGoing = false;
        }
      }
    }
  }
  if (!wfptr) {
    std::cout << "No Event Found" << std::endl;
  }

  //   UShort_t wfSz = wfptr->GetSize();
  //   digiAnalysis::WaveForm WFAveraged(wfSz, waveformVector);
  //   WFAveraged.SetSmooth(150);
  //   WFAveraged.SetTracesFFT("smooth");

  //   WFAveraged.Plot();

#endif
  fBL->Close();

  fApp->Run();
}