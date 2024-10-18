/*
**	Filename : Test_Dummy.cpp
**	2024-09-30
**	username : rsehgal
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

  std::string fname =
      "../data/DataF_CH1@V1730_167_Cs_CeBr_0pt5Vpp_40lsb_Calib_Na_CsI_"
      "0pt5Vpp_20lsb_Calib.root";

  // test reading to singleHits
  digiAnalysis::Analysis an(fname, 0, 0);

  std::cout << "getting the vector from an" << std::endl;

  // test Getting
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHits();

  std::cout << "got the vector from an" << hitsVector.size() << std::endl;

  // test Printing
  hitsVector[0]->Print();
  // hitsVector[2]->Print();
  std::cout << "_________________________________________________" << std::endl;
  std::cout << "_________________________________________________" << std::endl;
  std::cout << "_________________________________________________" << std::endl;
  std::cout << "_________________________________________________" << std::endl;
  std::cout << "_________________________________________________" << std::endl;
  std::cout << "_________________________________________________" << std::endl;

  // test sorting
  an.SortHits("Energy", "Time");
  // hitsVector = an.GetSingleHits();
  hitsVector[0]->Print();
  // hitsVector[2]->Print();

  int nentries = hitsVector.size();
  std::cout << "hitsVector size = " << nentries << std::endl;
  std::cout << "Energy: " << hitsVector[nentries - 2]->GetEnergy() << std::endl;
  digiAnalysis::WaveForm *WF = hitsVector[nentries - 2]->GetWFPtr();
  std::cout << "Got the waveform with size" << WF->GetSize() << std::endl;
  std::cout << "Got the waveform with baseline" << WF->GetBaseLine()
            << std::endl;
  std::cout << "Got the waveform with meantime" << WF->GetMeanTime()
            << std::endl;
  WF->Plot();

  fApp->Run();
  return 0;
}
