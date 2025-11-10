#include "Analysis.h"
#include "WaveForm.h"
#include "globals.h"
#include "singleHits.h"
#include <RtypesCore.h>
#include <TApplication.h>
#include <iostream>
#include <vector>

int PeakLocation(digiAnalysis::WaveForm *WF, int peakWidth, int startPoint = -1,
                 int threshold = 0);

int main(int argc, char *argv[]) {
#ifdef WAVES
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::string fname =
      "/home/kirtikesh/analysisSSD/DATA/SPE/"
      "run_Nov06_Direct_SelfTrigger_10lsb_BlueLED_Pico15_"
      "Bias2100V/FILTERED/"
      "DataF_run_Nov06_Direct_SelfTrigger_10lsb_BlueLED_Pico15_Bias2100V.root";

  digiAnalysis::Analysis an(fname, 0, 1, 0);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int nentries = hitsVector.size();
  std::cout << "Got the vector from an: " << nentries << std::endl;
  int evi = 0;
  std::string userInput;
  bool keepGoing = true;
  digiAnalysis::WaveForm *WF = nullptr;
  WF = hitsVector[0]->GetWFPtr();
  int peak = PeakLocation(WF, 10, 25, 5);
  std::cout << "Peak Position: " << peak << std::endl;
  WF->Plot();
  fApp->Run();
#endif
}

int PeakLocation(digiAnalysis::WaveForm *WF, int peakWidth, int startPoint,
                 int threshold) {
  if (startPoint == -1)
    startPoint = peakWidth;

  if (startPoint >= WF->GetSize()) {
    std::cout << "Error: Start Point is outside waveform size" << std::endl;
    return 0;
  }

  int peakPos = startPoint;
  std::vector<double> traces = WF->GetTraces();
  double peakVal = -0.1;

  while (peakVal < threshold) {
    int iPlus = 1, iMinus = 1;
    peakVal = traces[peakPos];
    while ((iPlus < peakWidth or iMinus < peakWidth) and
           peakPos < WF->GetSize() - peakWidth) {

      bool checkPlus =
          iPlus < peakWidth ? (peakVal - traces[peakPos + iPlus] > 0) : 1;
      bool checkMinus =
          iMinus < peakWidth ? (peakVal - traces[peakPos - iMinus] > 0) : 1;
      std::cout << "Pos: " << peakPos << " checkPlus: " << checkPlus
                << " checkMinus: " << checkMinus << std::endl;
      std::cout << "iPlus: " << iPlus << " iMinus: " << iMinus << std::endl;
      std::cout << "currPeakVal: " << peakVal
                << " valPlus: " << traces[peakPos + std::min(iPlus, peakWidth)]
                << " at " << peakPos + std::min(iPlus, peakWidth)
                << " valMinus: "
                << traces[peakPos - std::min(iMinus, peakWidth)] << " at "
                << peakPos - std::min(iMinus, peakWidth) << std::endl;
      if (checkPlus and checkMinus) {
        iMinus = std::min(iMinus + 1, peakWidth);
        iPlus = std::min(iPlus + 1, peakWidth);
        continue;
      } else if (!checkPlus) {
        peakPos = peakPos + iPlus;
        peakVal = traces[peakPos];
        iMinus = std::min(iMinus + iPlus + 1, peakWidth);
        iPlus = 1;
      } else if (!checkMinus) {
        peakPos = peakPos - iMinus;
        peakVal = traces[peakPos];
        iPlus = std::min(iMinus + iPlus + 1, peakWidth);
        iMinus = 1;
      }
    }
    if (peakVal < threshold)
      peakPos = std::min(int(WF->GetSize()),
                         peakPos + (peakPos - startPoint) /
                                       std::abs(peakPos - startPoint) *
                                       peakWidth * 2);
  }
  return peakPos;
}