#include "Analysis.h"
#include "TVirtualPad.h"
#include "WaveForm.h"
#include "globals.h"
#include "singleHits.h"
#include <TApplication.h>
#include <iostream>
#include "TComplex.h"
#include "TVirtualFFT.h"
#include "TSpectrum.h"
int main(int argc, char *argv[])
{
#ifdef WAVES
    TApplication *fApp = new TApplication("TEST", NULL, NULL);
    std::cout << "hello DigiAnalysis..." << std::endl;

    std::string fname = "/media/kirtikesh/UbuntuFiles/NaI/"
                        "run_Cs_FAGain_2_10_CFDTHR_15_10_Mode_EXT_TRG_FREEWRITE_"
                        "SignalDelay_50ns_Aug26/FILTERED/"
                        "DataF_run_Cs_FAGain_2_10_CFDTHR_15_10_Mode_EXT_TRG_"
                        "FREEWRITE_SignalDelay_50ns_Aug26.root";

    // test reading to singleHits
    digiAnalysis::Analysis an(fname, 0000, 20000, 0);

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

    int N = 0;
    int kernelSz = 100;
    int stepSz = 1;
    double sumKernel = 0;
    double sig = kernelSz / 10;
    std::vector<double> gauss(kernelSz);
    sumKernel = 0;
    for (int i = 0; i < kernelSz; i++)
    {
        gauss[i] = TMath::Exp(-0.5 * TMath::Power((i - (kernelSz / 2)) / sig, 2));
        sumKernel += gauss[i];
    }

    TCanvas *canvas = new TCanvas("canvas", "WaveForm Plot", 1900, 1000);
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    TVirtualFFT *fft = TVirtualFFT::FFT(1, &kernelSz, "R2C");
    TH1D *projectionX = nullptr;
    TSpectrum *s = new TSpectrum();
    Double_t *xpeaks = nullptr;
    Double_t *ypeaks = nullptr;
    double peakThreshSTFT = 1000;
    TH1F *hFinal = new TH1F("hFinal", "Final Waveform", 4000, 0, 4000);
    for (evi = 0; evi < nentries && keepGoing; ++evi) //
    {
        if (evi % 1000 == 0)
        {
            std::cout << evi << std::endl;
        }
        if (hitsVector[evi]->GetChNum() == 9 and
            fabs(hitsVector[evi]->GetMeanTime() - 2.6) < 0.2 and
            fabs(hitsVector[evi]->GetEnergy() - 200) < 100)
        {
            hitsVector[evi]->Print();
            WF = hitsVector[evi]->GetWFPtr();
            WF->SetSmooth(80);
            std::vector<double> traces = WF->GetTracesSmooth();
            N = traces.size();
            TH2 *hTracesShort = new TH2F("hTracesShort", "Short Time Traces", N - kernelSz - 1, 0, N - kernelSz - 1, kernelSz, 0, kernelSz);
            TH2 *hSTFT = new TH2F("hSTFT", "Short Time Fourier Transform", N - kernelSz - 1, 0, N - kernelSz - 1, kernelSz / 2 + 1, 0, kernelSz / 2 + 1);

            double mintrsh = 0, minstft = 0;
            for (int i = 0; i < N - kernelSz - 1; i += stepSz)
            {
                hTracesShort->Clear();
                hSTFT->Clear();
                std::vector<double> tracesShort(kernelSz);
                auto minIt = std::min_element(traces.begin() + i, traces.begin() + i + kernelSz);
                double minVal = *minIt;
                for (int j = 0; j < kernelSz; j++)
                {
                    tracesShort[j] = gauss[j] * (traces[i + j] - minVal) * 100 / sumKernel;
                    if (tracesShort[j] < mintrsh)
                        mintrsh = tracesShort[j];
                    hTracesShort->Fill(i + 0.01, j + 0.01, tracesShort[j] / 100);
                }
                fft->SetPoints(tracesShort.data());
                fft->Transform();
                // std::vector<double> re(kernelSz / 2 + 1), im(kernelSz / 2 + 1);
                double re = 0, im = 0, mag = 0;
                for (int j = 0; j <= kernelSz / 2; j++)
                {
                    fft->GetPointComplex(j, re, im);
                    mag = TMath::Sqrt(re * re + im * im);
                    // mag = (re * re + im * im);
                    if (mag < minstft)
                        minstft = mag;
                    hSTFT->Fill(i + 0.01, j + 0.01, mag);
                    // std::cout << mag << " : ";
                }
                // std::cout << std::endl;
            }

            if (projectionX)
                projectionX->Clear();
            std::cout << "HIST CLEANED" << std::endl;
            projectionX = hSTFT->ProjectionX("projX", 0, 5);
            std::cout << "PROJECTION DONE" << std::endl;
            Int_t nfound = s->Search(projectionX, 2, "", 0.05); // sigma=2 bins, threshold=5% of max
            std::cout << nfound << " PEAKS FOUND" << std::endl;
            xpeaks = s->GetPositionX();
            ypeaks = s->GetPositionY();
            for (int numpeak = 0; numpeak < nfound; numpeak++)
            {
                std::cout << xpeaks[numpeak] << " : " << ypeaks[numpeak] << std::endl;
            }

            canvas->Clear();
            canvas->Divide(2, 1);
            TGraph *graphTraces = nullptr;
            TGraph *graphTracesTH2 = nullptr;
            canvas->cd(1);
            graphTraces = new TGraph(N);
            graphTracesTH2 = new TGraph(N);
            for (int i = 0; i < N; ++i)
            {
                graphTraces->SetPoint(i, i, traces[i]);
                if (i < N - kernelSz)
                    graphTracesTH2->SetPoint(i, i, (traces[i + kernelSz / 2] + kernelSz / 2));
                else
                    graphTracesTH2->SetPoint(i, i, (kernelSz / 2));
                // graphTraces->SetPoint(i, i, 0);
            }
            // graphTraces->SetLineColor(kBlue);
            // graphTraces->SetLineWidth(2);
            // graphTraces->SetTitle("Traces");
            // graphTraces->Draw("AL");
            // legend->AddEntry(graphTraces, "Traces", "l");
            double xmin = graphTraces->GetXaxis()->GetXmin();
            double xmax = graphTraces->GetXaxis()->GetXmax();
            TLine *line = new TLine(xmin, kernelSz / 2, xmax, kernelSz / 2);
            line->SetLineColor(kBlack);
            line->SetLineStyle(1);

            // canvas->cd(3);
            hTracesShort->SetStats(0);
            hTracesShort->Draw("COLZ");
            hTracesShort->SetMinimum(mintrsh);
            graphTracesTH2->SetLineColor(kBlack);
            graphTracesTH2->SetLineWidth(2);
            graphTracesTH2->Draw("SAME");
            line->Draw("SAME");

            canvas->cd(2);
            hSTFT->SetStats(0);
            hSTFT->SetMinimum(minstft);
            hSTFT->Draw("COLZ");

            canvas->Update();

            std::cout << "Do you want to see the next waveform? (y/n): ";
            std::getline(std::cin, userInput);
            if (userInput != "y" && userInput != "Y")
            {
                keepGoing = false;
            }
        }
    }
    fApp->Run();
#endif
    return 0;
}