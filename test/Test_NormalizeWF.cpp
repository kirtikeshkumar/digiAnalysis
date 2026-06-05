#include "Analysis.h"
#include "TMath.h"
#include "WaveForm.h"
#include "globals.h"
#include "includes.hh"
#include "singleHits.h"
#include <TApplication.h>
#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <vector>
int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);

  std::string fname = "/home/kirtikesh/Analysis/DATA/LeadPit/CopperLining/"
                      "CoincidenceStudies/01JuneNoSrc/"
                      "NaI1342_02June26_1750_1345_1350_1350_NoSrc_Thresh_30_"
                      "300_WAVES_Coinc_144ns_LeadPit_Sum_BLCorrected.root";

  // Read to singleHits
  digiAnalysis::Analysis an(fname, 0, 00, 0);
  std::cout << "getting the vector from an" << std::endl;

  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  int nentries = hitsVector.size();
  std::cout << "got the vector from an: " << nentries << std::endl;

  double energy;
  digiAnalysis::WaveForm *wfptr = nullptr;
  std::vector<digiAnalysis::WaveForm> waveformVector;
  int min = 50, range = 10;
  std::vector<double> traceSmooth;
  std::vector<double> normtraces;
  std::vector<double> cfdtraces;
  std::pair<double, std::vector<double>> normPair;
  digiAnalysis::WaveForm *normWF = new digiAnalysis::WaveForm();
  int expectedPeakPos = digiAnalysis::GateStart + 80;
  int start = 0;
  int delay = 25;
  double fraction = 0.4, cfdthreshold = 5;

  TH1 *hscale = new TH1F("hscale", "hist of norm", 2000, 0, 500);

  for (int iter = 0; iter < nentries; iter++) {
    if (hitsVector[iter]->GetChNum() == 0) {
      energy = hitsVector[iter]->GetEnergy() * 0.09032 - 3.385;
      if (energy > min and energy < min + range) {
        // std::cout << iter << " : ";
        wfptr = hitsVector[iter]->GetWFPtr();
        wfptr->SetSmooth(40);
        normPair = wfptr->NormWaveForm();
        // std::cout << iter
        //           << " retrieved normPair with scale: " << normPair.first
        //           << std::endl;
        // now shift the waveforms properly so that the peaks line up.
        // find the peak position
        traceSmooth = wfptr->GetTracesSmooth();
        int it = 0;
        // auto max_it =
        //     std::max_element(traceSmooth.begin() +
        //     digiAnalysis::GateStart,
        //                      traceSmooth.begin() +
        //                      digiAnalysis::GateStart +
        //                          digiAnalysis::GateLenShort);
        // int it = std::distance(traceSmooth.begin(), max_it);
        // std::cout << iter << " peak Position is at: " << it << std::endl;

        cfdtraces.clear();
        for (int itertrace = 0; itertrace < traceSmooth.size(); itertrace++) {
          if (itertrace + delay < traceSmooth.size()) {
            cfdtraces.push_back(traceSmooth[itertrace] -
                                fraction * traceSmooth[itertrace + delay]);
          } else {
            cfdtraces.push_back(
                traceSmooth[itertrace] -
                fraction * traceSmooth[itertrace + delay - traceSmooth.size()]);
          }
        }
        bool ready = false;
        for (int itercfd = digiAnalysis::GateStart; itercfd < cfdtraces.size();
             itercfd++) {
          if (cfdtraces[itercfd - 1] < -1 * cfdthreshold) {
            ready = true;
          }
          if (ready && cfdtraces[itercfd] > 0) {
            it = itercfd;
            break;
          }
        }

        //   std::cout << iter << " Zero crossing at: " << it << std::endl;
        //   wfptr->Plot(traceSmooth, cfdtraces);
        //   break;

        // now create the waveform
        normtraces = normPair.second;
        // std::cout << iter << " trace size is: " << normtraces.size()
        //           << std::endl;
        start = it - expectedPeakPos;
        // std::cout << iter << " start point: " << start << std::endl;
        // normWF = nullptr;
        // std::cout << iter << " cleared" << std::endl;
        normWF->SetWaveForm(normtraces, start, start + normtraces.size());
        // std::cout << iter << " trace situation is: " << normWF->IsTracesSet()
        //           << ": size: " << normWF->GetTracesSize() << std::endl;
        waveformVector.push_back(*normWF);
        hscale->Fill(1.0 / normPair.first / (energy + 0.5 * range));
        // if (normWF->GetSize() < 1) {
        //   std::cout << iter << " failed to set WF: " << start << " : "
        //             << start + normtraces.size() << std::endl;
        //   std::cout << iter << " Zero crossing at: " << it
        //             << " normtrace size: " << normWF->GetSize()
        //             << " supposed to be: " << normtraces.size() << std::endl;
        //   //   wfptr->Plot(normtraces, normWF->GetTraces());
        //   //   break;
        // }
      }
    }
  }

  TCanvas *c1 = new TCanvas("c1", "hist of scale", 1500, 1000);
  hscale->Draw("HIST");
  c1->Update();

  if (!waveformVector.empty()) {
    UShort_t wfSz = wfptr->GetSize();
    digiAnalysis::WaveForm WFAveraged(wfSz, waveformVector);
    WFAveraged.SetSmooth(40);
    // WFAveraged.SetTracesFFT("smooth");
    WFAveraged.Plot();
  }
  fApp->Run();
}