#include "Analysis.h"
#include "globals.h"
#include "WaveForm.h"
#include "includes.hh"
#include "singleHits.h"
#include <RtypesCore.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TString.h>
#include <algorithm>
#include <iostream>
#include <ratio>
#include <vector>

int main(int argc, char *argv[])
{
    TApplication *fApp = new TApplication("TEST", NULL, NULL);
    std::cout << "hello DigiAnalysis..." << std::endl;
    std::string fname =
        "/media/kirtikesh/UbuntuFiles/NaI/LowEnergy/"
        "NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp/FILTERED/"
        "SDataF_NaI_13_CoincidenceStudies_Cs_HV_1900V_1365V_240min_2Vpp.root";

    digiAnalysis::Analysis an(fname, 0, 00, 0);
    std::cout << "getting the vector from an" << std::endl;

    std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
        an.GetSingleHitsVec();
    int nentries = hitsVector.size();
    std::cout << "got the vector from an: " << nentries << std::endl;

    an.CreatePairs();
    std::vector<digiAnalysis::Pair *> vecOfPairs = an.GetPairsVec();
    int nPairs = vecOfPairs.size();
    std::cout << nPairs << " Pairs were formed in the data." << std::endl;

    double Energy1 = 0;
    double Energy2 = 0;
    double PSD = 0, MT = 0;
    digiAnalysis::singleHits *hit;
    digiAnalysis::WaveForm *WF = nullptr;
#ifdef WAVES
    TH1 *hSPE = new TH1F("hSPE", "hSPE", 15000, 0, 15000);
    bool keepGoing = true;
    std::string userInput;
    for (int iter = 0; iter < nPairs && keepGoing; iter++)
    {
        hit = vecOfPairs[iter]->GetHit(0);
        hit->GetEnergy() > 694
            ? Energy1 = vecOfPairs[iter]->GetPairHitEnergy(0) * 0.09465 - 5.7613
            : Energy1 = vecOfPairs[iter]->GetPairHitEnergy(0) * 0.08696 - 0.4222;
        // 1900V
        if (Energy1 > 30 and Energy1 < 100)
        {
            WF = nullptr;
            WF = hit->GetWFPtr();
            WF->SetSmooth(65);
            auto results = WF->DetectPeakValleys(3);
            // std::cout << "size of peaks: " << results.first.size() << std::endl;
            // std::cout << "size of valleys: " << results.second.size() << std::endl;
            int iter = 0;
            // std::cout << std::endl
            //           << "PEAKS:____________" << std::endl;
            while (iter < results.first.size())
            {
                int peakPos = results.first[iter];
                if ((peakPos > 2500 and peakPos < 4000) and (peakPos - results.first[iter - 1] > 150))
                {
                    if ((iter + 1 < results.first.size() and (results.first[iter + 1] - peakPos) > 150) || (iter + 1 == results.first.size()))
                    {
                        // std::cout << iter << ":" << peakPos << std::endl;
                        hSPE->Fill(WF->IntegrateWaveForm(peakPos - 50, peakPos + 100) / 150.0 * digiAnalysis::EvalNormFactor);
                    }
                }

                iter += 1;
            }

            // iter = 0;

            // std::cout << std::endl
            //           << "VALLEYS:____________" << std::endl;
            // while (iter < results.second.size())
            // {
            //     std::cout << iter << ":" << results.second[iter] //<< ":" << traces[results.second[iter]]
            //               << std::endl;
            //     iter += 1;
            //     // if (iter >= results.second.size())
            //     //     break;
            // }

            // WF->Plot();
            // std::cout << "Do you want to see the next waveform? (y/n): ";
            // std::getline(std::cin, userInput);
            // if (userInput != "y" && userInput != "Y")
            // {
            //     keepGoing = false;
            // }
        }
    }
    TCanvas *c1 = new TCanvas("c1", "SPECharge", 800, 600);
    hSPE->Draw("HIST");
    fApp->Run();
#endif
}