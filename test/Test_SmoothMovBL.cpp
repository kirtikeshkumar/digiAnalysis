#include "Analysis.h"
#include "TComplex.h"
#include "TVirtualFFT.h"
#include "TVirtualPad.h"
#include "WaveForm.h"
#include "globals.h"
#include "singleHits.h"
#include <TApplication.h>
#include <iostream>
#include <string>
int main(int argc, char *argv[]) {
#ifdef WAVES
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;

  // std::string fname = "/home/kirtikesh/analysisSSD/DATA/NaI/"
  //                     "run_Na_NaI12_FAGain_10_10_CFDTHR_5_255_Mode_EXT_TRG_"
  //                     "FREEWRITE_SignalDelay_80ns_Sep08/FILTERED/"
  //                     "DataF_run_Na_NaI12_FAGain_10_10_CFDTHR_5_255_Mode_EXT_"
  //                     "TRG_FREEWRITE_SignalDelay_80ns_Sep08.root";

  // std::string fname =
  //     "/home/kirtikesh/analysisSSD/DATA/NaI/"
  //     "run_Cs_FAGain_2_10_CFDTHR_15_10_Mode_EXT_TRG_FREEWRITE_"
  //     "SignalDelay_50ns_Aug26/FILTERED/"
  //     "DataF_run_Cs_FAGain_2_10_CFDTHR_15_10_Mode_EXT_TRG_FREEWRITE_"
  //     "SignalDelay_50ns_Aug26.root";

  std::string fname = "/home/kirtikesh/analysisSSD/DATA/SPE/run_noSource_19Sep/"
                      "FILTERED/DataF_run_noSource_19Sep.root";

  // test reading to singleHits
  digiAnalysis::Analysis an(fname, 0, 100000, 0);

  std::cout << "getting the vector from an" << std::endl;

  // test Getting
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();

  std::cout << "got the vector from an: " << hitsVector.size() << std::endl;

  // sorting by time
  // an.SortHits("Time", "Channel");

  int nentries = hitsVector.size();
  std::cout << "hitsVector size = " << nentries << std::endl;

  digiAnalysis::WaveForm *WF = nullptr;
  int evi = 0;
  std::string userInput;
  bool keepGoing = true;

  int spectralSz = 2048;
  TH1 *hE = new TH1F("hE", "Energy", spectralSz, 0, spectralSz);
  TH1 *hEEval = new TH1F("hEEval", "Energy Eval", spectralSz, 0, spectralSz);
  TH1 *hERough = new TH1F("hERough", "Energy Rough", spectralSz, 0, spectralSz);
  TH1 *hESmooth =
      new TH1F("hESmooth", "Energy Smooth", spectralSz, 0, spectralSz);
  /*for (evi = 28; evi < 45; ++evi) {
    if (evi % 100000 == 0) {
      std::cout << evi << std::endl;
    }
    if (hitsVector[evi]->GetChNum() == 5 and
        fabs(hitsVector[evi]->GetMeanTime() - 2.6) < 0.2) {
      hE->Fill(hitsVector[evi]->GetEnergy());
      WF = hitsVector[evi]->GetWFPtr();
      WF->SetTracesMovBLCorr();
      std::cout << "plotting: " << evi << std::endl;
      WF->Plot();
      std::this_thread::sleep_for(std::chrono::seconds(5));
      hEEval->Fill(hitsVector[evi]->GetEvalEnergy() / 2.0);
      hESmooth->Fill(WF->IntegrateBLCorrWaveForm(
                         digiAnalysis::GateStart,
                         digiAnalysis::GateStart + digiAnalysis::GateLenLong) /
                     digiAnalysis::GateLenLong * 32.0);
    }
  }*/

  TCanvas *canvas = new TCanvas("canvas", "WaveForm Plot", 1600, 1000);
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  for (evi = 0; evi < nentries && keepGoing; ++evi) // && keepGoing
  {
    if (evi % 10000 == 0) {
      std::cout << evi << std::endl;
    }
    if (hitsVector[evi]->GetChNum() == 9 and
        // fabs(hitsVector[evi]->GetMeanTime() - 2.7) < 0.3 and
        fabs(hitsVector[evi]->GetEnergy() - 500) < 20) {
      hitsVector[evi]->Print();
      hE->Fill(hitsVector[evi]->GetEnergy());
      hEEval->Fill(hitsVector[evi]->GetEvalEnergy());
      WF = hitsVector[evi]->GetWFPtr();
      // WF->SetSmooth(80);
      // std::vector<double> traces = WF->GetTracesSmooth();
      std::vector<double> traces = WF->GetTraces();
      //   std::cout << traces[1000] << std::endl;
      //   hitsVector[evi]->Print();
      //   WF->Plot();

      int N = traces.size();
      std::vector<double> tracesBL(N);
      std::vector<double> tracesNew(N);
      double sum = 0;
      int nBL = 16;
      double baselineVal = 0;
      int start = 45;
      double BLError = 0.5;
      double trCut = 1;

      // Evaluate a baseline from the first few entries
      // It is known that there is no signal here.
      for (int i = start; i < nBL + start; i++) {
        sum = sum + traces[i];
      }
      baselineVal = sum / nBL;

      // now evaluate the baseline at each point
      int prevcount = 0;
      std::vector<double> BLDevstack;
      double blprev = 0, bltemp = 0;
      double bldev = 0, bldevsum = 0;
      double edgediff = 0;

      for (int i = start + nBL; i < N; i++) {
        tracesBL[i - nBL / 2] = baselineVal;
        sum = sum - traces[i - nBL];
        sum = sum + traces[i];
        edgediff = (traces[i] - traces[i - 1]) *
                   (traces[i - nBL + 2] - traces[i - nBL + 1]);

        blprev = bltemp;
        bltemp = sum / nBL;

        BLDevstack.push_back(pow(bltemp - blprev, 2));
        bldevsum = bldevsum + BLDevstack.back();
        if (BLDevstack.size() >= nBL) {
          bldevsum = bldevsum - BLDevstack[BLDevstack.size() - nBL];
        }
        bldev = bldevsum / TMath::Min(nBL, (int)BLDevstack.size());

        if (prevcount > 0) {
          prevcount -= 1;
        }

        // if (fabs(traces[i] - baselineVal) > BLError) {
        if (fabs(bltemp - blprev) > BLError) {
          prevcount = nBL;
        } else if (prevcount == 0) // edgediff > 0.001 //(fabs(traces[i] -
                                   // baselineVal) < trCut)
        {
          baselineVal = bltemp;
        }
        // if (i > 700 and i < 852)
        // {
        //   std::cout << "i: " << i - nBL / 2 + 1 << "\t bltemp: " << bltemp <<
        //   "\t blprev: " << blprev << "\t edgediff: " << edgediff << "\t ERCK:
        //   " << fabs((bltemp - blprev) * edgediff) << "\t PC: " << prevcount
        //   << std::endl;
        // }
        // if (prevcount == 0)
        // {
        //   baselineVal = bltemp;
        // }

        // if (i < 1050 and i > 900) {
        //   std::cout << i << "\t : " << prevcount << "\t : " << traces[i]
        //             << "\t : " << baselineVal << std::endl;
        // }
      }

      for (int i = 0; i < N; i++) {
        tracesNew[i] = (traces[i] - tracesBL[i]);
      }
      digiAnalysis::WaveForm *wfNew = new digiAnalysis::WaveForm(tracesNew);

      // std::cout << "########### AFTER BL CORRECTION ##############"
      //           << std::endl;
      double evalenergy =
          wfNew->IntegrateWaveForm(digiAnalysis::GateStart,
                                   digiAnalysis::GateStart +
                                       digiAnalysis::GateLenLong) /
          64;
      double evalenergyshort =
          wfNew->IntegrateWaveForm(digiAnalysis::GateStart,
                                   digiAnalysis::GateStart +
                                       digiAnalysis::GateLenShort) /
          64;
      // std::cout << "EnergyEval:      " << evalenergy << std::endl;
      // std::cout << "EnergyEvalShort: " << evalenergyshort << std::endl;
      // std::cout << "EvalPSD:         " << 1.0 - evalenergyshort / evalenergy
      //           << "\n"
      //           << std::endl;
      hERough->Fill(evalenergy);

      TVirtualFFT *fft = TVirtualFFT::FFT(1, &N, "R2C");
      fft->SetPoints(tracesBL.data());
      fft->Transform();
      std::vector<double> re(N / 2 + 1), im(N / 2 + 1);
      for (int i = 0; i <= N / 2; i++) {
        fft->GetPointComplex(i, re[i], im[i]);
      }
      TVirtualFFT *ifft = TVirtualFFT::FFT(1, &N, "C2R");
      int cutoff = 10;
      for (int i = 0; i <= N / 2; i++) {
        if (i < cutoff) {
          TComplex c(re[i], im[i]);
          // std::cout << "re: " << re[i] << " im: " << im[i] << " c: " << c <<
          // std::endl;
          ifft->SetPointComplex(i, c);
        } else {
          TComplex c(0., 0.);
          ifft->SetPointComplex(i, c);
        }
      }
      ifft->Transform();
      std::vector<double> tracesBLFFT(N);
      ifft->GetPoints(&tracesBLFFT[0]);
      for (int i = 0; i < N; i++) {
        tracesBLFFT[i] /= N;
      }

      std::vector<double> tracesNewFFTCleaned(N);
      for (int i = 0; i < N; i++) {
        tracesNewFFTCleaned[i] = (traces[i] - tracesBLFFT[i]);
      }

      wfNew->SetWaveForm(tracesNewFFTCleaned);
      evalenergy = wfNew->IntegrateWaveForm(digiAnalysis::GateStart,
                                            digiAnalysis::GateStart +
                                                digiAnalysis::GateLenLong) /
                   digiAnalysis::GateLenLong * digiAnalysis::EvalNormFactor;
      std::cout << "Eval Energy: " << evalenergy << std::endl;
      hESmooth->Fill(evalenergy);

      digiAnalysis::WaveForm wf;
      // std::vector<double> tracesBLFFT = wf.EvalTracesFFT(tracesBL);
      std::vector<double> tracesFFT = wf.EvalTracesFFT(traces);
      std::vector<double> tracesNewFFT = wf.EvalTracesFFT(tracesNew);

      TGraph *graphTraces = nullptr;
      TGraph *graphTracesBL = nullptr;
      TGraph *graphTracesNew = nullptr;
      TGraph *graphTracesNewFFTCleaned = nullptr;
      TGraph *graphTracesBLFFT = nullptr;
      TGraph *graphTracesFFT = nullptr;
      TGraph *graphTracesNewFFT = nullptr;
      TGraph *graphTracesdiffFFT = nullptr;

      canvas->Clear();
      canvas->Divide(2, 2);
      canvas->cd(1);

      graphTraces = new TGraph(N);
      for (int i = 0; i < N; ++i) {
        graphTraces->SetPoint(i, i, traces[i]);
        // graphTraces->SetPoint(i, i, 0);
      }
      graphTraces->SetLineColor(kBlue);
      graphTraces->SetLineWidth(2);
      graphTraces->SetTitle("Traces");
      graphTraces->Draw("AL");
      legend->AddEntry(graphTraces, "Traces", "l");

      graphTracesNew = new TGraph(N);
      for (int i = 0; i < N; ++i) {
        graphTracesNew->SetPoint(i, i, tracesNew[i]);
      }
      graphTracesNew->SetLineColor(kGreen);
      graphTracesNew->SetLineWidth(2);
      graphTracesNew->SetTitle("Traces New");
      graphTracesNew->Draw("L SAME");
      legend->AddEntry(graphTracesNew, "Traces New", "l");

      graphTracesNewFFTCleaned = new TGraph(N);
      for (int i = 0; i < N; ++i) {
        graphTracesNewFFTCleaned->SetPoint(i, i, tracesNewFFTCleaned[i]);
      }
      graphTracesNewFFTCleaned->SetLineColor(kRed);
      graphTracesNewFFTCleaned->SetLineWidth(2);
      graphTracesNewFFTCleaned->SetTitle("Traces New FFT Cleaned");
      graphTracesNewFFTCleaned->Draw("L SAME");
      legend->AddEntry(graphTracesNewFFTCleaned, "Traces New  FFT Cleaned",
                       "l");

      double xmin = graphTracesNew->GetXaxis()->GetXmin();
      double xmax = graphTracesNew->GetXaxis()->GetXmax();

      TLine *line = new TLine(xmin, 0, xmax, 0);
      line->SetLineColor(kBlack);
      line->SetLineStyle(1);
      line->Draw("SAME");

      canvas->cd(2);

      graphTracesBL = new TGraph(N);
      for (int i = 0; i < N; ++i) {
        graphTracesBL->SetPoint(i, i, tracesBL[i]);
      }
      graphTracesBL->SetLineColor(kRed);
      graphTracesBL->SetLineWidth(2);
      graphTracesBL->SetTitle("Traces BL");
      graphTracesBL->Draw("AL");
      legend->AddEntry(graphTracesBL, "Traces BL", "l");

      TLine *line1 = new TLine(xmin, 0, xmax, 0);
      line1->SetLineColor(kBlack);
      line1->SetLineStyle(1);
      line1->Draw("SAME");

      canvas->cd(3);
      // gPad->DrawFrame(0.1, 10, 100, 10000, "TtracesDiff");
      // N = 300;
      graphTracesFFT = new TGraph(N / 2);
      for (int i = 0; i < N / 2; ++i) {
        if (!tracesFFT.empty()) {
          graphTracesFFT->SetPoint(i, i, tracesFFT[i]);
        } else {
          graphTracesFFT->SetPoint(i, i, 0);
        }
      }
      graphTracesFFT->SetLineColor(kBlue);
      graphTracesFFT->SetLineWidth(2);
      graphTracesFFT->SetTitle("Traces FFT");
      graphTracesFFT->Draw("AL");
      legend->AddEntry(graphTracesFFT, "Traces FFT", "l");

      graphTracesNewFFT = new TGraph(N / 2);
      for (int i = 0; i < N / 2; ++i) {
        if (!tracesNewFFT.empty()) {
          graphTracesNewFFT->SetPoint(i, i, tracesNewFFT[i]);
        } else {
          graphTracesNewFFT->SetPoint(i, i, 0);
        }
      }
      graphTracesNewFFT->SetLineColor(kGreen);
      graphTracesNewFFT->SetLineWidth(2);
      graphTracesNewFFT->SetTitle("Traces New FFT");
      graphTracesNewFFT->Draw("L SAME");
      legend->AddEntry(graphTracesNewFFT, "Traces New FFT", "l");

      // graphTracesdiffFFT = new TGraph(N / 2);
      // for (int i = 0; i < N / 2; ++i)
      // {
      //   if (!tracesNewFFT.empty())
      //   {
      //     graphTracesdiffFFT->SetPoint(i, i, tracesNewFFT[i] - tracesFFT[i]);
      //   }
      //   else
      //   {
      //     graphTracesdiffFFT->SetPoint(i, i, 0);
      //   }
      // }
      // graphTracesdiffFFT->SetLineColor(kBlack);
      // graphTracesdiffFFT->SetLineWidth(2);
      // graphTracesdiffFFT->SetTitle("Traces diff FFT");
      // graphTracesdiffFFT->Draw("L SAME");
      // legend->AddEntry(graphTracesdiffFFT, "Traces diff FFT", "l");

      gPad->SetLogy();

      canvas->cd(4);
      // N = tracesNew.size();
      graphTracesBLFFT = new TGraph(N);
      for (int i = 0; i < N; ++i) {
        if (!tracesBLFFT.empty()) {
          graphTracesBLFFT->SetPoint(i, i, tracesBLFFT[i]);
        } else {
          graphTracesBLFFT->SetPoint(i, i, 0);
        }
      }
      graphTracesBLFFT->SetLineColor(kRed);
      graphTracesBLFFT->SetLineWidth(2);
      graphTracesBLFFT->SetTitle("Traces BL FFT");
      graphTracesBLFFT->Draw("AL");
      legend->AddEntry(graphTracesBLFFT, "Traces BL FFT", "l");

      // gPad->SetLogy();

      canvas->Update();
      std::cout << "Do you want to see the next waveform? (y/n): ";
      std::getline(std::cin, userInput);
      if (userInput != "y" && userInput != "Y") {
        keepGoing = false;
      }
    }
  }
  hE->Draw("HIST");
  hEEval->SetLineColor(kBlack);
  hEEval->Draw("HIST SAME");
  hERough->SetLineColor(kRed);
  hERough->Draw("HIST SAME");
  hESmooth->SetLineColor(kGreen);
  hESmooth->Draw("HIST SAME");
  hE->SetTitle("Energy Comparison;Energy (a.u.);Counts");
  fApp->Run();
#endif
  return 0;
}