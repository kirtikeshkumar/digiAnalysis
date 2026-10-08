#include "Analysis.h"
#include "Pair.h"
#include "WaveForm.h"
#include "globals.h"
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
#include <string>
#include <vector>

int main(int argc, char *argv[])
{
    std::cout << "hello DigiAnalysis..." << std::endl;

    int startIndex = 17, endIndex = 58;
    std::string fpath = "/home/kirtikesh/Analysis/DATA/extCoincSep/1800V/Calib/";
    std::string finit =
        "NaI134_CsSrc_1800_1337_1350_WAVES_FILTERED_NoSplitSignal_Gain_12_"
        "Acquisition_4_ExtTrig_Threshold1LSB_160nsPromt_240nsDelay_800nsLong_"
        "1496nsCoinc_FreeWrites_Ch248Singles_";
    for (int fileiter = startIndex; fileiter < endIndex; fileiter++)
    {
        std::string fname = fpath + finit + std::to_string(fileiter) + "/FILTERED/DataF_" + finit + std::to_string(fileiter) +
                            "_BLCorrected.root";

        digiAnalysis::Analysis an(fname, 1000000, 2000000, 0);
        std::cout << "getting the vector from an" << std::endl;

        std::vector<std::unique_ptr<digiAnalysis::singleHits>> &hitsVector =
            an.GetSingleHitsVec();
        int nentries = hitsVector.size();
        std::cout << "got the vector from an: " << nentries << std::endl;

        // create pairs to identitfy true photon events
        an.CreatePairs();
        std::vector<std::unique_ptr<digiAnalysis::Pair>> &vecOfPairs =
            an.GetPairsVec();
        int nPairs = vecOfPairs.size();
        std::cout << nPairs << " Pairs were formed in the data." << std::endl;

        std::string outfname =
            fpath + "PairFiles/Pair_" + finit + std::to_string(fileiter) + ".root";

        TFile *fout = TFile::Open(outfname.c_str(), "RECREATE");
        TTree *t = new TTree("Data_Pair", "Data_Pair");

        digiAnalysis::Pair pairObj;
        t->Branch("pair", "digiAnalysis::Pair",
                  &pairObj); // or similar depending on your class setup

        for (const std::unique_ptr<digiAnalysis::Pair> &p : vecOfPairs)
        {
            pairObj.ClearPair();
            pairObj.SetPair(*p->GetHitPtr(0), *p->GetHitPtr(1));
            t->Fill();
        }

        fout->Write();
        fout->Close();
    }
}