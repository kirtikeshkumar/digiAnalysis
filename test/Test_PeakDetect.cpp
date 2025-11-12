#include "Analysis.h"
#include "WaveForm.h"
#include "globals.h"
#include "singleHits.h"
#include <RtypesCore.h>
#include <TApplication.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

int PeakLocation(digiAnalysis::WaveForm *WF, int peakWidth,
                 int startPoint = -1);

std::pair<std::vector<int>, std::vector<int>>
DetectPeakValleys(digiAnalysis::WaveForm *WF, double threshold = 0);

int main(int argc, char *argv[]) {
#ifdef WAVES
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::string fname = "/home/kirtikesh/analysisSSD/DATA/SPE/"
                      "run_Nov11_1GSPS_Direct_FreeWrites_CFD_3lsb_Delay_20ns_"
                      "BlueLED_PicoOFFHalfON_Fiber_Bias1800V/FILTERED/"
                      "DataF_run_Nov11_1GSPS_Direct_FreeWrites_CFD_3lsb_Delay_"
                      "20ns_BlueLED_PicoOFFHalfON_Fiber_Bias1800V.root";

  digiAnalysis::Analysis an(fname, 30, 1, 0);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int nentries = hitsVector.size();
  std::cout << "Got the vector from an: " << nentries << std::endl;
  int evi = 0;
  std::string userInput;
  bool keepGoing = true;
  digiAnalysis::WaveForm *WF = nullptr;
  WF = hitsVector[0]->GetWFPtr();
  int Width = 10;
  int startpoint = 0;
  int peak = PeakLocation(WF, Width, startpoint);
  std::set<int> identifiedPeaks;
  std::set<int> selectedPeaks;
  auto res = identifiedPeaks.insert(peak);
  std::vector<double> traces = WF->GetTracesSmooth();
  double threshold = 4;
  // while ((traces[peak] - traces[peak + Width]) / traces[peak] <= 0.4 or
  //        traces[peak] < threshold) { // this is the value expected from a
  //                                    // gaussian peak with width=sigma
  //   std::cout << "The width of identified peak at: " << peak
  //             << " is broader than required ("
  //             << (traces[peak] - traces[peak + Width]) / traces[peak]
  //             << ") CONTINUING" << std::endl;
  //   if (!res.second) { // Here I check if the peak identified in the previous
  //                      // run was different from previously identified peaks.
  //                      // This is required since if the previous peak is at
  //                      // the same location, i need to move farther for the
  //                      // next peak
  //     startpoint += 2 * Width;
  //   } else {
  //     startpoint = peak + 2 * Width;
  //   }
  //   if (startpoint >= WF->GetSize()) {
  //     std::cout << "No Peak of required or smaller width could be found"
  //               << std::endl;
  //     break;
  //   } else {
  //     peak = PeakLocation(WF, Width, startpoint);
  //     res = identifiedPeaks.insert(peak);
  //     std::cout << peak << " " << res.second << std::endl;
  //   }
  // }

  std::cout << "Peak Position: " << peak << " ("
            << (traces[peak] - traces[peak + Width]) / traces[peak] << ")"
            << std::endl;

  WF->SetSmooth(13);
  traces = WF->GetTracesSmooth();
  auto results = WF->DetectPeakValleys(3.5);
  std::cout << "size of peaks: " << results.first.size() << std::endl;
  std::cout << "size of valleys: " << results.second.size() << std::endl;
  int iter = 0;
  std::cout << std::endl << "PEAKS:____________" << std::endl;
  while (iter < results.first.size()) {
    int peakPos = results.first[iter];
    double ratiohi = 1.0 - traces[peakPos + 1] / traces[peakPos];
    double ratiolo = 1.0 - traces[peakPos + 2] / traces[peakPos];
    std::cout << iter << ":" << peakPos << ":" << ratiolo << ":" << ratiohi
              << std::endl;
    iter += 1;
  }

  iter = 0;

  std::cout << std::endl << "VALLEYS:____________" << std::endl;
  while (iter < results.second.size()) {
    std::cout << iter << ":" << results.second[iter] << std::endl;
    iter += 1;
  }

  WF->Plot();
  fApp->Run();
#endif
}

int PeakLocation(digiAnalysis::WaveForm *WF, int peakWidth, int startPoint) {
  if (startPoint == -1)
    startPoint = digiAnalysis::GateStart;

  if (startPoint >= WF->GetSize()) {
    std::cout << "Error: Start Point is outside waveform size" << std::endl;
    return 0;
  }

  int peakPos = startPoint;
  int iPlus = 1, iMinus = 1;
  WF->SetSmooth(2.5 * peakWidth);
  std::vector<double> traces = WF->GetTracesSmooth();
  double peakVal = traces[peakPos];
  double valPlus, valMinus;
  while ((iPlus < peakWidth or iMinus < peakWidth) and
         peakPos < WF->GetSize() - peakWidth) {

    valPlus = peakPos + iPlus < WF->GetSize() ? traces[peakPos + iPlus]
                                              : traces[WF->GetSize() - 1];
    valMinus = peakPos - iMinus > 0 ? traces[peakPos - iMinus] : traces[0];
    bool checkPlus = iPlus < peakWidth ? (peakVal - valPlus > 0) : 1;
    bool checkMinus = iMinus < peakWidth ? (peakVal - valMinus > 0) : 1;
    // std::cout << "Pos: " << peakPos << " checkPlus: " << checkPlus
    //           << " checkMinus: " << checkMinus << std::endl;
    // std::cout << "iPlus: " << iPlus << " iMinus: " << iMinus << std::endl;
    // std::cout << "currPeakVal: " << peakVal << " valPlus: " << valPlus << "
    // at "
    //           << peakPos + std::min(iPlus, peakWidth)
    //           << " valMinus: " << valMinus << " at "
    //           << peakPos - std::min(iMinus, peakWidth) << std::endl;
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
  return peakPos;
}

std::pair<std::vector<int>, std::vector<int>>
DetectPeakValleys(digiAnalysis::WaveForm *WF, double threshold) {
  std::vector<double> traces = WF->GetTracesSmooth();
  std::vector<int> peak;
  std::vector<int> valley;
  std::vector<int> valleyTemp;
  int iter = 0, peakPos = 0, valleyPos = 0;
  bool findPeak = true, findValley = true, peakFound = false;
  double peakVal = traces[peakPos], valleyVal = traces[valleyPos], currVal;
  std::cout << "traces size: " << traces.size() << std::endl;
  while (iter < traces.size()) {
    // std::cout << iter << " : " << traces[iter] << " : " << peakPos << " : "
    //           << peakVal << " : " << valleyPos << " : " << valleyVal
    //           << std::endl;
    currVal = traces[iter];
    if (findPeak and currVal < peakVal) {
      findPeak = false;
      findValley = true;
      if (peakVal > threshold) {
        peak.push_back(peakPos);
        peakFound = true;
        if (!valleyTemp.empty()) {
          valley.push_back(valleyTemp.back());
        }
      }
    } else {
      peakPos = iter;
      peakVal = currVal;
    }

    if (findValley and currVal > valleyVal) {
      findPeak = true;
      findValley = false;
      valleyTemp.push_back(valleyPos);
      if (peakFound == true) {
        peakFound = false;
        valley.push_back(valleyTemp.back());
      }
    } else {
      valleyPos = iter;
      valleyVal = currVal;
    }

    iter += 1;
  }
  std::cout << "peak vector in function " << peak.size() << std::endl;
  return std::make_pair(peak, valley);
}