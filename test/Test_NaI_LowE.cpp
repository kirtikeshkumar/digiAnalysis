/*
**	Filename : Test_NaI_LowE.cpp
**	2025-07-05
**	username : kirtikesh
*/

#include "Analysis.h"
#include "WaveForm.h"
#include "globals.h"
#include "singleHits.h"
#include <TApplication.h>
#include <iostream>
int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;

  // std::string fname =
  // "/media/kirtikesh/Ventoy/GGAG/DataF_CH1@V1730_167_GGAG_2inch_insideHutch_3Oct_4hr.root";

  // std::string fname =
  // "/media/kirtikesh/Ventoy/GGAG/DataF_CH1@V1730_167_GGAG_2inch_insideHutch_3Oct_calib_Na.root
  // ";
  // std::string fname =
  //     "/media/kirtikesh/Ventoy/NaI/"
  //     "run_NaI1_CsSrAl_AnodeDynodeCoinc_AmpAnode4Dynode8_2Vpp/FILTERED/"
  //     "DataF_run_NaI1_CsSrAl_AnodeDynodeCoinc_AmpAnode4Dynode8_2Vpp.root";

  std::string fname =
      "/media/kirtikesh/KirtikeshSSD/DATA/NaI/"
      "Calib_Waves_NaI12_Na_Coinc_PSDCut0pt4_AmpAnode10_2Vpp/FILTERED/"
      "DataF_Calib_Waves_NaI12_Na_Coinc_PSDCut0pt4_AmpAnode10_2Vpp.root";

  // std::string fname =
  //     "/media/kirtikesh/KirtikeshSSD/DATA/NaI/"
  //     "run_NaI1_Cs_AnodeDynodeCoinc_AmpAnode4Dynode8_2Vpp/FILTERED/"
  //     "DataF_run_NaI1_Cs_AnodeDynodeCoinc_AmpAnode4Dynode8_2Vpp.root";

  // test reading to singleHits
  digiAnalysis::Analysis an(fname, 0, 0, 0);

  std::cout << "getting the vector from an" << std::endl;

  // test Getting
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();

  std::cout << "got the vector from an: " << hitsVector.size() << std::endl;

  // sorting by time
  // an.SortHits("Time", "Channel");

  int nentries = hitsVector.size();
  std::cout << "hitsVector size = " << nentries << std::endl;

  // Plot Waveform to check
  int evi = 0;

  std::string userInput;
  bool keepGoing = true;

  digiAnalysis::WaveForm *WF = nullptr;
  std::vector<digiAnalysis::WaveForm> waveformVector;

  for (evi = 0; evi < nentries && keepGoing; ++evi) {
    if (evi % 1000 == 0) {
      std::cout << evi << std::endl;
    }
    if (hitsVector[evi]->GetChNum() == 2 and
        // hitsVector[evi]->GetEvalEnergy() > 100 and
        hitsVector[evi]->GetEnergy() < 250 and
        hitsVector[evi]->GetEnergy() > 230 and
        fabs(hitsVector[evi]->GetMeanTime() - 2.6) < 0.2
        // hitsVector[evi]->GetPSD() < 0.6
        )
    // hitsVector[evi]->GetEnergy() - hitsVector[evi]->GetEnergyShort() < 0
    // and hitsVector[evi]->GetEvalEnergy() -
    // hitsVector[evi]->GetEvalEnergyShort() > 0 and
    // hitsVector[evi]->GetPSD() < 0.8
    {
      std::cout << "evi: " << evi << " hitNum: " << hitsVector[evi]->GetEvNum()
                << std::endl;

      // digiAnalysis::WaveForm *WF = hitsVector[evi]->GetWFPtr();
      WF = nullptr;
      WF = hitsVector[evi]->GetWFPtr();
      WF->SetSmooth(32);
      if (WF) {
        waveformVector.push_back(*WF);
      }
      std::cout << "Got the waveform with size" << WF->GetSize() << std::endl;
      std::cout << "Got the waveform with baseline" << WF->GetBaseLine()
                << std::endl;
      std::cout << "Got the waveform with meantime" << WF->GetMeanTime()
                << std::endl;
      // WF->SetSmooth(32);
      WF->FitExponential(340, 1100);
      hitsVector[evi]->Print();
      std::cout << "Smooth Integral: "
                << (WF->IntegrateSmoothWaveForm(
                       digiAnalysis::GateStart,
                       digiAnalysis::GateStart + digiAnalysis::GateLenLong)) /
                       digiAnalysis::GateLenLong * 4.4
                << std::endl;
      WF->Plot();

      std::cout << "Do you want to see the next waveform? (y/n): ";
      std::getline(std::cin, userInput);
      if (userInput != "y" && userInput != "Y") {
        keepGoing = false;
      }
    }
  }

  UShort_t wfSz = WF->GetSize();
  std::cout << "Got size of waveforms" << wfSz << std::endl;
  digiAnalysis::WaveForm *WFAveraged =
      new digiAnalysis::WaveForm(wfSz, waveformVector);
  WFAveraged->FitExponential(340, 1100);
  WFAveraged->SetSmooth(32);
  digiAnalysis::singleHits *hitAveraged =
      new digiAnalysis::singleHits(0, 0, 0, 0, 0, 0, WFAveraged);
  hitAveraged->Print();
  WFAveraged->Plot();

  // evi = 15;
  // digiAnalysis::WaveForm *WF = hitsVector[evi]->GetWFPtr();
  // WF->SetSmooth(32);
  // std::cout << "Got the waveform with size" << WF->GetSize() << std::endl;
  // std::cout << "Got the waveform with baseline" << WF->GetBaseLine()
  //           << std::endl;
  // std::cout << "Got the waveform with meantime" << WF->GetMeanTime()
  //           << std::endl;
  // WF->SetSmooth(32);
  // WF->FitExponential(340, 1100);
  // hitsVector[evi]->Print();
  // WF->Plot();
  fApp->Run();
  return 0;
}