/*
**	Filename : Test_NaI_LowE.cpp
**	2025-07-05
**	username : kirtikesh
*/

#include "Analysis.h"
#include "WaveForm.h"
#include "includes.hh"
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
  std::string fname = "/media/kirtikesh/KirtikeshSSD/DATA/NaI/"
                      "run_NaI0_AnodeDynodeCoinc_AmpAnode2Dynode8/FILTERED/"
                      "DataF_run_NaI0_AnodeDynodeCoinc_AmpDynode2.root";

  // std::string fname =
  //     "/media/kirtikesh/KirtikeshSSD/DATA/NaI/"
  //     "run_NaI1_Cs_AnodeDynodeCoinc_AmpAnode4Dynode8_2Vpp/FILTERED/"
  //     "DataF_run_NaI1_Cs_AnodeDynodeCoinc_AmpAnode4Dynode8_2Vpp.root";

  // test reading to singleHits
  digiAnalysis::Analysis an(fname, 0, 10000, 0);

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

  for (evi = 0; evi < nentries && keepGoing; ++evi) {
    if (evi % 100 == 0) {
      std::cout << evi << std::endl;
    }
    if (hitsVector[evi]->GetChNum() == 0 and
        // hitsVector[evi]->GetEvalEnergy() > 100 and
        // hitsVector[evi]->GetEvalEnergy() < 200 and
        // hitsVector[evi]->GetMeanTime() < 2.8
        hitsVector[evi]->GetPSD() < 0.6)
    // hitsVector[evi]->GetEnergy() - hitsVector[evi]->GetEnergyShort() < 0
    // and hitsVector[evi]->GetEvalEnergy() -
    // hitsVector[evi]->GetEvalEnergyShort() > 0 and
    // hitsVector[evi]->GetPSD() < 0.8
    {
      std::cout << "evi: " << evi << " hitNum: " << hitsVector[evi]->GetEvNum()
                << std::endl;

      digiAnalysis::WaveForm *WF = hitsVector[evi]->GetWFPtr();
      WF->SetSmooth(32);
      std::cout << "Got the waveform with size" << WF->GetSize() << std::endl;
      std::cout << "Got the waveform with baseline" << WF->GetBaseLine()
                << std::endl;
      std::cout << "Got the waveform with meantime" << WF->GetMeanTime()
                << std::endl;
      WF->SetSmooth(32);
      WF->FitExponential(340, 1100);
      hitsVector[evi]->Print();
      WF->Plot();

      std::cout << "Do you want to see the next waveform? (y/n): ";
      std::getline(std::cin, userInput);
      if (userInput != "y" && userInput != "Y") {
        keepGoing = false;
      }
    }
  }

  evi = 15;
  digiAnalysis::WaveForm *WF = hitsVector[evi]->GetWFPtr();
  WF->SetSmooth(32);
  std::cout << "Got the waveform with size" << WF->GetSize() << std::endl;
  std::cout << "Got the waveform with baseline" << WF->GetBaseLine()
            << std::endl;
  std::cout << "Got the waveform with meantime" << WF->GetMeanTime()
            << std::endl;
  WF->SetSmooth(32);
  WF->FitExponential(340, 1100);
  hitsVector[evi]->Print();
  WF->Plot();
  fApp->Run();
  return 0;
}