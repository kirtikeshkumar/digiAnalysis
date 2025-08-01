/*
**	Filename : Test_Dummy.cpp
**	2024-12-03
**	username : kirtikeshkumar
*/

#include "Analysis.h"
#include "TMath.h"
#include "WaveForm.h"
#include "includes.hh"
#include "singleHits.h"
#include <TApplication.h>
#include <iostream>
int main(int argc, char *argv[]) {
  TApplication *fApp = new TApplication("TEST", NULL, NULL);

  // std::string fname =
  //     "/media/kirtikesh/KirtikeshSSD/Waves/CeBr_CsI_0pt5Vpp/"
  //     "Am_CsI_950V_Waves_CFD_Delay150ns_Fraction50_10lsb_0pt5Vpp_8IS_300PG_"
  //     "6000Gate_5696Holdoff_TConnector_1/FILTERED/"
  //     "DataF_CH2@V1730_167_Am_CsI_950V_Waves_CFD_Delay150ns_Fraction50_10lsb_"
  //     "0pt5Vpp_8IS_300PG_6000Gate_5696Holdoff_TConnector_1.root";

  // std::string fname =
  //     "/media/kirtikesh/Ventoy/NaI/"
  //     "run_NaI1_CsSrAl_AnodeDynodeCoinc_AmpAnode4Dynode8_2Vpp/FILTERED/"
  //     "DataF_run_NaI1_CsSrAl_AnodeDynodeCoinc_AmpAnode4Dynode8_2Vpp.root";

  std::string fname =
      "/media/kirtikesh/KirtikeshSSD/DATA/NaI/"
      "Calib_Waves_NaI12_Na_Coinc_PSDCut0pt4_AmpAnode10_2Vpp/FILTERED/"
      "DataF_Calib_Waves_NaI12_Na_Coinc_PSDCut0pt4_AmpAnode10_2Vpp.root";

  // std::string fname =
  //     "/media/kirtikesh/KirtikeshSSD/DATA/NaI/"
  //     "run_NaI1_Cs_AnodeDynodeCoinc_AmpAnode4Dynode8_2Vpp/FILTERED/"
  //     "DataF_run_NaI1_Cs_AnodeDynodeCoinc_AmpAnode4Dynode8_2Vpp.root";

  // Read to singleHits
  digiAnalysis::Analysis an(fname, 0, 0, 1);

  // Get the vector
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
#ifdef WAVES
  int spectralsize = 8192;
  double WFMT = 0; // MeanTime parameter for flat waveform
  TH2 *hMTPlot = new TH2F("MTPlot", "Energy vs MeanTime", spectralsize, 0,
                          spectralsize, 500, -4, 4);
  TH2 *hPSDPlot = new TH2F("PSDPlot", "Energy vs PSD", spectralsize, 0,
                           spectralsize, 500, -8, 2);
  TH2 *hEPlot = new TH2F("EPlot", "Energy vs EvalEnergy", 410, 0, spectralsize,
                         410, 0, spectralsize);
  TH2 *hESPlot = new TH2F("ESPlot", "EnergyShort vs EvalEnergyShort", 410, 0,
                          spectralsize, 410, 0, spectralsize);
  TH2 *hPSDEvalPlot =
      new TH2F("PSDEvalPlot", "PSD vs PSDEval", 160, -8, 8, 160, -8, 8);
  TH1F *hEEvalRatio =
      new TH1F("hEEValRatio", "Ratio of EEval vs E", 1000, 0, 200);
  TH2 *hEdiffPlot = new TH2F("EdiffPlot", "Energy vs Energy-EnergyShort", 410,
                             0, spectralsize, 820, -spectralsize, spectralsize);
  TH2 *hEdiffEvalPlot =
      new TH2F("EdiffEvalPlot", "EvalEnergy vs EvalEnergy-EvalEnergyShort", 410,
               0, spectralsize, 820, -spectralsize, spectralsize);
  TH1F *hESpectra =
      new TH1F("hESpectra", "Energy Spectra", spectralsize, 0, spectralsize);
  TH1F *hEEvalSpectra = new TH1F("hEEvalSpectra", "Eval Energy Spectra",
                                 spectralsize, 0, spectralsize);
  int nentries = hitsVector.size();
  double psd = 0;
  double evalEnergy = 0;
  double evalEnergyShort = 0;
  double energy = 0;
  double energyShort = 0;
  digiAnalysis::WaveForm *WF = nullptr;
  for (int i = 0; i < nentries; i++) {
    if (hitsVector[i]->GetChNum() == 3
        // and hitsVector[i]->GetMeanTime() < 3.8
    ) { // (hitsVector[i]->GetPSD() > 0.0 and
        // hitsVector[i]->GetChNum() == 0) {
      energy = hitsVector[i]->GetEnergy();
      energyShort = hitsVector[i]->GetEnergyShort();
      WF = hitsVector[i]->GetWFPtr();
      WFMT = TMath::Log10(WF->GetSize() / 2.0);
      hMTPlot->Fill(energy, hitsVector[i]->GetMeanTime());
      evalEnergy =
          hitsVector[i]->GetEvalEnergy(); // WF->IntegrateWaveForm(290, 1390);
      evalEnergyShort =
          hitsVector[i]
              ->GetEvalEnergyShort(); // WF->IntegrateWaveForm(290, 440);
      psd = 1.0 - evalEnergyShort * 1.0 / evalEnergy;
      hPSDPlot->Fill(energy, psd);
      hEPlot->Fill(energy, evalEnergy);
      hESPlot->Fill(energyShort, evalEnergyShort);
      hPSDEvalPlot->Fill(hitsVector[i]->GetPSD(), psd);
      hEEvalRatio->Fill(evalEnergyShort * 1.0 / energyShort);
      hEdiffPlot->Fill(energy, energy - energyShort);
      hEdiffEvalPlot->Fill(evalEnergy, evalEnergy - evalEnergyShort);
      hESpectra->Fill(energy);
      hEEvalSpectra->Fill(evalEnergy);
      if (i < 1300 and i > 1285) {
        std::cout << i << std::endl;
        hitsVector[i]->Print();
        std::cout << "Eval Energy     : " << evalEnergy << std::endl;
        std::cout << "Eval EnergyShort: " << evalEnergyShort << std::endl;
        std::cout << "Eval PSD        : " << psd << std::endl;
      }
    }
  }
  TCanvas *c1 = new TCanvas("c1", "Energy vs MeanTime", 800, 600);
  hMTPlot->Draw("COLZ");
  TCanvas *c2 = new TCanvas("c2", "Energy vs PSD", 800, 600);
  hPSDPlot->Draw("COLZ");
  TCanvas *c3 = new TCanvas("c3", "Energy vs evalEnergy", 800, 600);
  hEPlot->Draw("COLZ");
  // hEEvalRatio->Draw("HIST");
  TCanvas *c4 = new TCanvas("c4", "EnergyShort vs evalEnergyShort", 800, 600);
  hESPlot->Draw("COLZ");
  TCanvas *c5 = new TCanvas("c5", "PSD vs PSDEval", 800, 600);
  hPSDEvalPlot->Draw("COLZ");
  TCanvas *c6 = new TCanvas("c6", "Energy vs Energy-EnergyShort", 800, 600);
  hEdiffPlot->Draw("COLZ");
  TCanvas *c7 =
      new TCanvas("c7", "evalEnergy vs evalEnergy-evalEnergyShort", 800, 600);
  hEdiffEvalPlot->Draw("COLZ");
  TCanvas *c8 = new TCanvas("c8", "Energy Spectra", 800, 600);
  hESpectra->SetLineColor(kRed);
  hEEvalSpectra->SetLineColor(kGreen);
  hESpectra->Draw("HIST");
  hEEvalSpectra->Draw("HISTSAME");
  fApp->Run();
#endif
  return 0;
}
