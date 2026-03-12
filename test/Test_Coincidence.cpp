#include "Analysis.h"
#include "WaveForm.h"
#include "includes.hh"
#include "singleHits.h"
#include <RtypesCore.h>
#include <TApplication.h>
#include <iostream>

int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;
  std::string fname =
      "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/CoincidenceStudies/"
      "NaI_12_CoincidenceStudies_Na_HV_GainMatch_10min_2Vpp/FILTERED/"
      "DataF_NaI_12_CoincidenceStudies_Na_HV_GainMatch_10min_2Vpp.root";

  digiAnalysis::Analysis an(fname, 0, 0000, 10);
  std::cout << "getting the vector from an" << std::endl;

  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int nentries = hitsVector.size();
  std::cout << "got the vector from an" << nentries << std::endl;

  an.SortHits("Channel", "Time");

  return 0;
}