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
int main(int argc, char *argv[])
{
  TApplication *fApp = new TApplication("TEST", NULL, NULL);
  std::cout << "hello DigiAnalysis..." << std::endl;

  // std::string fname =
  // "/media/kirtikesh/Ventoy/GGAG/DataF_CH1@V1730_167_GGAG_2inch_insideHutch_3Oct_4hr.root";

  // std::string fname =
  // "/media/kirtikesh/Ventoy/GGAG/DataF_CH1@V1730_167_GGAG_2inch_insideHutch_3Oct_calib_Na.root
  // ";
  std::string fname =
      "/media/kirtikesh/47AEF45B17DB643D1/NaI/Am/FILTERED_LowE_Am/"
      "DataF_CH0@V1751_1231_NaI_60025_00911-4_1400V_CeBr_600V_LowE_Am.root";

  // test reading to singleHits
  digiAnalysis::Analysis an(fname, 0, 1000000, 0);

  std::cout << "getting the vector from an" << std::endl;

  // test Getting
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();

  std::cout << "got the vector from an: " << hitsVector.size() << std::endl;

  // test Printing
  // hitsVector[0]->Print();
  // // hitsVector[2]->Print();
  // std::cout << "_________________________________________________" << std::endl;
  // std::cout << "_________________________________________________" << std::endl;
  // std::cout << "_________________________________________________" << std::endl;
  // std::cout << "_________________________________________________" << std::endl;
  // std::cout << "_________________________________________________" << std::endl;
  // std::cout << "_________________________________________________" << std::endl;

  // test sorting
  // an.SortHits("Energy", "PSD");
  // hitsVector[0]->Print();
  // hitsVector[2]->Print();

  int nentries = hitsVector.size();
  std::cout << "hitsVector size = " << nentries << std::endl;

  // test smoothing and plotting
  int evi = 0;
  for (evi = 2010; evi < nentries; evi++)
  {
    if (hitsVector[evi]->GetEnergy() > 1750 and hitsVector[evi]->GetEnergy() < 1800 and hitsVector[evi]->GetMeanTime() < 2.7) // and hitsVector[evi]->GetMeanTime() < 2.35)
    {
      std::cout << evi << std::endl;
      break;
    }
  }
  hitsVector[evi]->Print();
  std::cout << "Energy: " << hitsVector[evi]->GetEnergy() << std::endl;
  digiAnalysis::WaveForm *WF = hitsVector[evi]->GetWFPtr();
  // WF->SetSmooth(16);
  std::cout << "Got the waveform with size" << WF->GetSize() << std::endl;
  std::cout << "Got the waveform with baseline" << WF->GetBaseLine()
            << std::endl;
  std::cout << "Got the waveform with meantime" << WF->GetMeanTime()
            << std::endl;
  WF->SetSmooth(32);
  WF->Plot();

  // test averaging
  digiAnalysis::WaveForm *wfptr = nullptr;
  std::vector<digiAnalysis::WaveForm> waveformVector;
  UShort_t wfSz = 0;
  Long64_t start = 1;
  UShort_t count = 20000;
  for (Long64_t i = start; i < start + count; i++)
  {
    if (hitsVector[i]->GetPSD() < 1.0)
    {
      wfptr = hitsVector[i]->GetWFPtr();
      if (wfptr and (hitsVector[i]->GetEnergy() > 170 and hitsVector[i]->GetEnergy() < 180) and hitsVector[i]->GetMeanTime() < 2.44)
      {
        waveformVector.push_back(*wfptr);
      }
    }
    else
    {
      wfptr = nullptr;
    }
  }
  wfSz = wfptr->GetSize();
  digiAnalysis::WaveForm WFAveraged(wfSz, waveformVector);
  std::cout << "Plot Averaged of: " << waveformVector.size() << std::endl;
  WFAveraged.Plot();

  // test plot data
  TH2 *hEnPSD =
      new TH2F("hEnPSD", "Energy(x) vs PSD(y)", 8192, 0, 8192, 819, 0, 1);
  TH2 *hEnMT = new TH2F("hEnMT", "Energy(x) vs MeanTime(y)", 8192, 0, 8192, 1000, 0, 4);
  TH2 *hEnIEvT = new TH2F("hEnIEvT", "Energy(x) vs InterEventTime(y)", 8192, 0, 8192, 1000, 5, 13);
  // std::cout << "hEnPSD created" << std::endl;
  // TH2 *hNumPSD = new TH2F("hNumPSD", "Index(x) vs PSD(y)", nentries / 100, 0,
  // nentries, 819, 0, 1);
  // TH2 *hNumEn = new TH2F("hNumEn", "Index(x) vs Energy(y)", nentries / 100, 0,
  //                        nentries, 8192, 0, 8192);
  float delT = 0;
  for (int i = 0; i < nentries; i++)
  {
    hEnPSD->Fill(hitsVector[i]->GetEnergy(), hitsVector[i]->GetPSD());
    hEnMT->Fill(hitsVector[i]->GetEnergy(), hitsVector[i]->GetMeanTime());
    if (i > 1) // nentries - 1)
    {
      delT = hitsVector[i]->GetTimestamp() - hitsVector[i - 1]->GetTimestamp();
      hEnIEvT->Fill(hitsVector[i]->GetEnergy(), TMath::Log10(delT));
    }
    // hNumPSD->Fill(i, hitsVector[i]->GetPSD());
    // hNumEn->Fill(i, hitsVector[i]->GetEnergy());
  }

  TCanvas *c1 = new TCanvas("c1", "Energy vs PSD", 800, 600);
  // TCanvas *c2 = new TCanvas("c2", "Index vs PSD", 800, 600);
  TCanvas *c3 = new TCanvas("c3", "Energy vs MeanTime", 800, 600);
  TCanvas *c4 = new TCanvas("c4", "Energy vs InterEventTime", 800, 600);
  c1->cd();
  hEnPSD->Draw("COLZ");

  // c2->cd();
  // hNumPSD->Draw("COLZ");

  c3->cd();
  // hNumEn->Draw("COLZ");
  hEnMT->Draw("COLZ");

  c4->cd();
  hEnIEvT->Draw("COLZ");

  c1->Update();
  // c2->Update();
  c3->Update();
  c4->Update();

  // std::cout << "Fin Final" << std::endl;

  fApp->Run();
  return 0;
}
