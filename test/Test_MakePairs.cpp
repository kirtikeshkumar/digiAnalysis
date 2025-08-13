#include "Analysis.h"
#include "TF1.h"
#include "TMath.h"
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
  //     "/media/kirtikesh/UbuntuFiles/NaI/Calib_Waves_NaI12_Na_Coinc_PSDCut0pt4_AmpAnode10_2Vpp/FILTERED/DataF_Calib_Waves_NaI12_Na_Coinc_PSDCut0pt4_AmpAnode10_2Vpp.root";

  std::string fname =
      "/media/kirtikesh/UbuntuFiles/NaI/run_NaI1_Cs_AnodeDynodeCoinc_AmpAnode4Dynode8_2Vpp/FILTERED/DataF_run_NaI1_Cs_AnodeDynodeCoinc_AmpAnode4Dynode8_2Vpp.root";

  digiAnalysis::Analysis an(fname, 0, 00000, 0);
  std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
      an.GetSingleHitsVec();
  std::vector<digiAnalysis::singleHits> hitsVectorCh0 = an.GetSingleHitsVec(5);
  std::vector<digiAnalysis::singleHits> hitsVectorCh1 = an.GetSingleHitsVec(4);

  std::cout << "got the vector of channel 0 from an: " << hitsVectorCh0.size()
            << std::endl;
  std::cout << "got the vector of channel 1 from an: " << hitsVectorCh1.size()
            << std::endl;

  ULong64_t i = 0;
  ULong64_t j = 0;

  Long64_t Ch0Time = hitsVectorCh0[0].GetTimestamp();
  Long64_t Ch1Time = hitsVectorCh1[0].GetTimestamp();
  Long64_t TDiffCurr = Ch0Time - Ch1Time;
  Long64_t TDiffMin = 0;

  TH1F *hICDelT =
      new TH1F("hICDelT", "inter channel nearest event time difference", 80000,
               -40000, 40000);

  // for (i = 0; i < hitsVector.size() - 1; i++)
  // {
  //   TDiffMin = hitsVector[i + 1]->GetTimestamp() - hitsVector[i]->GetTimestamp();
  //   hICDelT->Fill(TDiffMin / 1000000.0);
  // }

  while (i < hitsVectorCh0.size())
  {
    // if (fabs(hitsVectorCh0[i].GetEnergy() - 5050) < 150)
    // {
    Ch0Time = hitsVectorCh0[i].GetTimestamp();
    // }
    // else
    // {
    //   i = i + 1;
    // }
    // std::cout << i << " Ch: " << hitsVectorCh0[i].GetChNum() << " T0: " << Ch0Time
    //           << std::endl;
    TDiffMin = 100000000000000;
    while (j < hitsVectorCh1.size())
    {
      // // if (fabs(hitsVectorCh1[j].GetEnergy() - 5050) < 150)
      // // {
      // //   Ch1Time = hitsVectorCh1[j].GetTimestamp();
      // // }
      // else
      // {
      //   j = j + 1;
      //   break;
      // }
      Ch1Time = hitsVectorCh1[j].GetTimestamp();
      // std::cout << j << " Ch: " << hitsVectorCh1[j].GetChNum() << " T1: " << Ch1Time
      //           << std::endl;
      TDiffCurr = Ch1Time - Ch0Time;
      // std::cout << i << " : " << j << " : Current " << TDiffCurr / 1000000.0
      //           << std::endl;
      if (fabs(TDiffCurr) < fabs(TDiffMin))
      {
        TDiffMin = TDiffCurr;
        // std::cout << i << " : " << j << " : Minimum " << TDiffMin / 1000000.0
        //           << std::endl;
        j = j + 1;
      }
      else
      {
        // hitsVectorCh0[i].Print();
        // hitsVectorCh1[j - 1].Print();
        hICDelT->Fill(TDiffMin / 1000.0);
        i = i + 1;
        j = (j >= 2) ? j - 2 : j;
        break;
      }
    }
    if (i >= hitsVectorCh0.size() - 1 or j >= hitsVectorCh1.size() - 1)
    {
      break;
    }
  }

  std::cout << "Loop Exited" << std::endl;

  hICDelT->Draw("HIST");

  fApp->Run();
  return 0;
}
