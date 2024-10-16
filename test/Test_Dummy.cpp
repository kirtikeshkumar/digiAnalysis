/*
**	Filename : Test_Dummy.cpp
**	2024-09-30
**	username : rsehgal
*/

#include "Analysis.h"
#include "singleHits.h"
#include <iostream>
int main(int argc, char *argv[]) {
  std::cout << "hello DigiAnalysis..." << std::endl;

  std::string fname =
      "../data/DataF_CH1@V1730_167_Cs_CeBr_0pt5Vpp_40lsb_Calib_Na_CsI_"
      "0pt5Vpp_20lsb_Calib.root";

  digiAnalysis::Analysis an(fname, 0, 0);

  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHits();

  hitsVector[0]->Print();

  return 0;
}
