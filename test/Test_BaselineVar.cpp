#include "Analysis.h"
#include "Pair.h"
#include "TF1.h"
#include "TMath.h"
#include "WaveForm.h"
#include "includes.hh"
#include "singleHits.h"
#include <TApplication.h>
#include <cmath>
#include <iostream>

int main(int argc, char *argv[]) {
#ifdef WAVES
  TApplication *fApp = new TApplication("TEST", NULL, NULL);

  // Extract folder path and name from argv[1]
  std::string folderPath = std::string(argv[1]);
  // Remove trailing slash if present
  if (!folderPath.empty() && folderPath.back() == '/')
    folderPath.pop_back();
  // Find the last '/' to get the folder name
  size_t lastSlash = folderPath.find_last_of('/');
  std::string folderName = (lastSlash != std::string::npos)
                               ? folderPath.substr(lastSlash + 1)
                               : folderPath;
  // Construct the data file name
  std::string dataFileName =
      folderPath + "/FILTERED/DataF_" + folderName + ".root";

  ULong64_t numEvt = 0;
  std::sscanf(argv[2], "%lld", &numEvt);
  ushort chNum = 0;
  std::sscanf(argv[3], "%hd", &chNum);

  digiAnalysis::Analysis an(dataFileName, 0, numEvt, 0.);

  std::vector<digiAnalysis::singleHits *> hitsVector =
      an.GetSingleHitsVec(chNum);

  int nentries = hitsVector.size();
  std::cout << "hitsVector size = " << nentries << std::endl;

  TH1F *hBaseLine =
      new TH1F("hBaseLine", "Baseline Distribution", 100000, 0, 100);

  for (int i = 0; i < nentries; i++) {
    double baseline = hitsVector[i]->GetWFPtr()->GetBaseLine();
    double timeStamp = hitsVector[i]->GetTimestamp() / 1E12; // in  seconds
    // std::cout << "timeStamp: " << timeStamp << " baseline: " << baseline <<
    // std::endl;
    hBaseLine->Fill(timeStamp, baseline);
  }
  hBaseLine->Draw("HIST");
  fApp->Run();
  return 0;
#endif
}