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
  for (evi = 155; evi < nentries; ++evi) {
    if (evi % 100 == 0) {
      std::cout << evi << std::endl;
    }
    if (hitsVector[evi]->GetChNum() == 0 and
        hitsVector[evi]->GetEnergy() - hitsVector[evi]->GetEnergyShort() < 0 and
        hitsVector[evi]->GetEvalEnergy() -
                hitsVector[evi]->GetEvalEnergyShort() >
            0 and
        hitsVector[evi]->GetPSD() < 0.8) {
      std::cout << "evi: " << evi << " hitNum: " << hitsVector[evi]->GetEvNum()
                << std::endl;
      hitsVector[evi]->Print();
      break;
    }
  }
  // evi = 100;
  if (evi < nentries) {
    // hitsVector[evi]->Print();
    digiAnalysis::WaveForm *WF = hitsVector[evi]->GetWFPtr();
    // psd = 1.0 - WF->IntegrateWaveForm(290, 440) * 1.0 /
    //                 WF->IntegrateWaveForm(290, 1390);
    // std::cout<<
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
  } else {
    std::cout << "E.O.F: no event of required type found" << std::endl;
  }

  // create anode-dynode pairs
  // 0 is anode, 1 is dynode
  // It is observed that in the noAmp file, removing energy < 5 ADC in channel 1
  // results in same number of events as in channel 0

  // TH2 *hDelT =
  //     new TH2F("hDelT", "DelT between channels", 8192, 0, 8192, 1000, 0, 4);

  // for (evi=0 : evi<nentries:evi++) {

  // }

  fApp->Run();
  return 0;
}