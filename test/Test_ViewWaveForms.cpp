#include "Analysis.h"
#include "WaveForm.h"
#include "includes.hh"
#include "singleHits.h"
#include <RtypesCore.h>
#include <TApplication.h>
#include <iostream>

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
  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //     "NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp/FILTERED/"
  //     "SDataF_NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp.root";

  std::string fname =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
      "NaI31_26May26_1345_1750_NoSrc_Thresh50_WAVES_NoCoinc_LeadPit_5/FILTERED/"
      "DataF_NaI31_26May26_1345_1750_NoSrc_Thresh50_WAVES_NoCoinc_LeadPit_5_"
      "BLCorrected.root";

  // std::string fname =
  //     "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
  //     "NaI31_26May26_1345_1750_NoSrc_Thresh50_WAVES_Coinc_144ns_LeadPit/"
  //     "FILTERED/"
  //     "SDataF_NaI31_26May26_1345_1750_NoSrc_Thresh50_WAVES_Coinc_144ns_LeadPit_"
  //     "BLCorrected.root";

  UShort_t channel = 2;
  digiAnalysis::Analysis an(2, fname, 0000, 00000, 1);
  std::cout << "getting the vector from an" << std::endl;

  // test Getting
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  std::cout << "got the vector from an" << hitsVector.size() << std::endl;

  // an.SortHits("PSD", "Energy");
  int nentries = hitsVector.size();
  std::cout << "hitsVector size = " << nentries << std::endl;

#ifdef WAVES

  // To PLot Based on Cuts
  UShort_t spectralsz = 4096;
  TH2 *hPSDPlot = new TH2F("PSDPlot", "Energy vs PSD", spectralsz, 0,
                           spectralsz, 100, 0, 1);
  TH2 *hMTPlot =
      new TH2F("MTPlot", "Energy vs MT", spectralsz, 0, spectralsz, 4096, 2, 4);
  TH2 *hLamPlot = new TH2F("LamPlot", "Energy vs Lam", spectralsz, 0,
                           spectralsz, 2000, -10, 10);
  TH2 *hPSDLamPlot =
      new TH2F("PSDLamPlot", "PSD vs Lam", 100, 0, 1, 1000, -10, 0);
  for (const auto &hit : hitsVector) {
    if (hit->GetChNum() == channel) {
      hPSDPlot->Fill(hit->GetEvalEnergy(), hit->GetEvalPSD());
      hMTPlot->Fill(hit->GetEvalEnergy(), hit->GetMeanTime());
      hLamPlot->Fill(hit->GetEvalEnergy(), hit->GetWFPtr()->EvalNoisePar2());
      hPSDLamPlot->Fill(hit->GetEvalPSD(), hit->GetWFPtr()->EvalNoisePar2());
    }
  }
  // TCanvas *c2 = new TCanvas("c2", "Energy vs PSD", 1500, 1000);
  // hPSDPlot->Draw("COLZ");
  // TCanvas *c2 = new TCanvas("c2", "Energy vs MeanTime", 1500, 1000);
  // hMTPlot->Draw("COLZ");
  // TCanvas *c2 = new TCanvas("c2", "Energy vs LamPar", 1500, 1000);
  // hLamPlot->Draw("COLZ");
  // TCanvas *c2 = new TCanvas("c2", "PSD vs LamPar", 1500, 1000);
  // hPSDLamPlot->Draw("COLZ");
  // c2->Update();

  TH1 *ESpect = new TH1F("ESpect", "Energy Spectra", spectralsz, 0, spectralsz);

  double Par1CutMin = 0, Par1CutMax = 16384;
  double Par2CutMin = 0, Par2CutMax = 1.0;
  std::string userInput;
  // std::cout << "\n Par1CutMin: ";
  // Par1CutMin = readUserInput(Par1CutMin);
  // std::cout << "\n Par1CutMax: ";
  // Par1CutMax = readUserInput(Par1CutMax);
  // std::cout << "\n Par2CutMin: ";
  // Par2CutMin = readUserInput(Par2CutMin);
  // std::cout << "\n Par2CutMax: ";
  // Par2CutMax = readUserInput(Par2CutMax);

  std::cout << "Parameter 1 cut: " << Par1CutMin << " -> " << Par1CutMax
            << std::endl;
  std::cout << "Parameter 2 cut: " << Par2CutMin << " -> " << Par2CutMax
            << std::endl;

  int count = 0;
  digiAnalysis::WaveForm *wfptr = nullptr;
  std::vector<digiAnalysis::WaveForm> waveformVector;
  bool keepGoing = true;
  for (int i = 0; i < nentries; i++) {
    if (hitsVector[i]->GetChNum() == channel &&
        // hitsVector[i]->GetEvalPSD() >= Par1CutMin &&
        // hitsVector[i]->GetEvalPSD() <= Par1CutMax &&
        // hitsVector[i]->GetMeanTime() >= 2.3 &&
        // hitsVector[i]->GetMeanTime() <= 2.5 &&
        // hitsVector[i]->GetWFPtr()->EvalNoisePar2() >= Par2CutMin &&
        // hitsVector[i]->GetWFPtr()->EvalNoisePar2() <= Par2CutMax &&
        hitsVector[i]->GetEvalEnergy() >= 00 &&
        hitsVector[i]->GetEvalEnergy() <= 100) {
      ESpect->Fill(hitsVector[i]->GetEvalEnergy());
      if (count < 500 && keepGoing) {
        count++;
        wfptr = nullptr;
        wfptr = hitsVector[i]->GetWFPtr();
        if (wfptr) {
          waveformVector.push_back(*wfptr);
        }
        hitsVector[i]->Print();
        wfptr->SetSmooth(40);
        wfptr->SetTracesFFT();
        wfptr->Plot();
        std::cout << "Do you want to see the next waveform? (y/n): ";
        std::getline(std::cin, userInput);
        if (userInput != "y" && userInput != "Y") {
          keepGoing = false;
        }
      }
      // else {
      //   break;
      // }
    }
  }
  // std::cout << "Now showing Average" << std::endl;

  TCanvas *c3 =
      new TCanvas("c3", "Eenrgy Spectra of Selected Region", 1500, 1000);
  ESpect->Draw("HIST");
  c3->Update();

  if (!waveformVector.empty()) {
    UShort_t wfSz = wfptr->GetSize();
    digiAnalysis::WaveForm WFAveraged(wfSz, waveformVector);
    WFAveraged.SetSmooth(40);
    WFAveraged.SetTracesFFT("smooth");
    // std::vector<double> Baseline = WFAveraged.GetTracesSmooth();
    // TFile *f = new TFile("baseline.root", "RECREATE");
    // TTree *t = new TTree("tree", "Baseline");
    // std::vector<double> *Baseline_ptr = &Baseline;
    // t->Branch("Baseline", &Baseline_ptr);
    // t->Fill();
    // t->Write();
    // f->Close();

    WFAveraged.Plot();
  }

  // To plot particular Events

  fApp->Run();
#endif
}