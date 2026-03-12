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
  std::string fname =
      "/home/kirtikesh/Analysis/DATA/LeadPit/NaI_1_Checks/"
      "NaI_1_Cs_Source_HV1800_Waves_160FC_0pt5Vpp_GroundToLead/FILTERED/"
      "DataF_NaI_1_Cs_Source_HV1800_Waves_160FC_0pt5Vpp_GroundToLead."
      "root";
  UShort_t channel = 0;
  digiAnalysis::Analysis an(fname, 65730, 0000, 10);
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
        if (wfptr) {
          waveformVector.push_back(*wfptr);
        }
        hitsVector[i]->Print();
        wfptr->SetSmooth(150);
        wfptr->SetTracesFFT("smooth");
        wfptr->Plot();
        std::cout << "Do you want to see the next waveform? (y/n): ";
        std::getline(std::cin, userInput);
        if (userInput != "y" && userInput != "Y") {
          keepGoing = false;
        }
      }
    }
  }
  UShort_t wfSz = wfptr->GetSize();
  digiAnalysis::WaveForm WFAveraged(wfSz, waveformVector);
  WFAveraged.SetSmooth(150);
  WFAveraged.SetTracesFFT("smooth");

  std::vector<double> Baseline = WFAveraged.GetTracesSmooth();
  TFile *f = new TFile("baseline.root", "RECREATE");
  TTree *t = new TTree("tree", "Baseline");
  std::vector<double> *Baseline_ptr = &Baseline;
  t->Branch("Baseline", &Baseline_ptr);
  t->Fill();
  t->Write();
  f->Close();

  WFAveraged.Plot();

  // To plot particular Events

  fApp->Run();
#endif
}