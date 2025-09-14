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
    digiAnalysis::Analysis an(fname, 0000, 0000, 0);

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
    int kernelSz = 40;
    int stepSz = 1;
    double sumKernel = 0;
    double sig = kernelSz / 8;
    std::vector<double> gauss(kernelSz);
    sumKernel = 0;
    for (int i = 0; i < kernelSz; i++)
    {
        gauss[i] = TMath::Exp(-0.5 * TMath::Power((i - (kernelSz / 2)) / sig, 2));
        sumKernel += gauss[i];
    }

    TCanvas *canvas = new TCanvas("canvas", "WaveForm Plot", 1900, 1000);
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    TVirtualFFT *fft = TVirtualFFT::FFT(1, &kernelSz, "R2C M K");
    TVirtualFFT *ifft = TVirtualFFT::FFT(1, &kernelSz, "C2R M K");
    TH1D *projectionX = nullptr;
    TSpectrum *s = new TSpectrum();
    Double_t *xpeaks = nullptr;
    Double_t *ypeaks = nullptr;
    double peakThreshSTFT = 5;
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
            std::vector<double> tracesMod(N, 0.0);
            std::vector<double> tracesPeaks(N, 0.0);
            std::vector<double> tracesPeaksUni(N, 0.0);
            std::vector<double> tracesPeaksUni1(N, 0.0);
            int rangeSTFT = N - kernelSz - 1;
            TH2 *hTracesShort = new TH2F("hTracesShort", "Short Time Traces", rangeSTFT, 0, rangeSTFT, kernelSz, 0, kernelSz);
            TH2 *hSTFT = new TH2F("hSTFT", "Short Time Fourier Transform", rangeSTFT, 0, rangeSTFT, kernelSz / 2 + 1, 0, kernelSz / 2 + 1);
            TH2 *hSTFT_Re = new TH2F("hSTFT_Re", "Short Time Fourier Transform Real Part", rangeSTFT, 0, rangeSTFT, kernelSz / 2 + 1, 0, kernelSz / 2 + 1);
            TH2 *hSTFT_Im = new TH2F("hSTFT_Im", "Short Time Fourier Transform Imaginary Part", rangeSTFT, 0, rangeSTFT, kernelSz / 2 + 1, 0, kernelSz / 2 + 1);
            std::cout << "HISTS CREATED" << std::endl;
            double mintrsh = 0, minstft = 0;
            for (int i = 0; i < rangeSTFT; i += stepSz)
            {
                hTracesShort->Clear();
                hSTFT->Clear();
                hSTFT_Re->Clear();
                hSTFT_Im->Clear();
                std::vector<double> tracesShort(kernelSz);
                auto minIt = std::min_element(traces.begin() + i, traces.begin() + i + kernelSz);
                double minVal = *minIt;
                for (int j = 0; j < kernelSz; j++)
                {
                    tracesShort[j] = gauss[j] * (traces[i + j] - minVal) * 100 / sumKernel; // the traces-minval makes the pulse to be only >0 and shifts the baseline accordingly
                    tracesMod[i + j] += tracesShort[j];
                    if (tracesShort[j] < mintrsh)
                        mintrsh = tracesShort[j];
                    hTracesShort->Fill(i + 0.01, j + 0.01, tracesShort[j] / 100);
                }
                fft->SetPoints(tracesShort.data());
                // std::cout << "FFT SET" << std::endl;
                fft->Transform();
                // std::cout << "FFT PERFORMED" << std::endl;
                // std::vector<double> re(kernelSz / 2 + 1), im(kernelSz / 2 + 1);
                double re = 0, im = 0, mag = 0;
                for (int j = 0; j <= kernelSz / 2; j++)
                {
                    fft->GetPointComplex(j, re, im);
                    // mag = TMath::Sqrt(re * re + im * im);
                    mag = TMath::Log10(re * re + im * im + 0.001);
                    // std::cout << mag << std::endl;
                    if (mag < minstft)
                        minstft = mag;
                    hSTFT->Fill(i + 0.01, j + 0.01, mag);
                    hSTFT_Re->Fill(i + 0.01, j + 0.01, re);
                    hSTFT_Im->Fill(i + 0.01, j + 0.01, im);
                    // std::cout << mag << " : ";
                }

                // std::cout << minstft << std::endl;
                // std::cout << std::endl;
            }
            // std::cout << "FFT DONE" << std::endl;

            if (projectionX)
                projectionX->Clear();
            // std::cout << "HIST CLEANED" << std::endl;
            projectionX = hSTFT->ProjectionX("projX", 0, 1); // here the projection is taken to get peak positions
            // std::cout << "PROJECTION DONE" << std::endl;
            Int_t nfound = s->Search(projectionX, 2, "", 0.05); // sigma=2 bins, threshold=5% of max
            std::cout << nfound << " PEAKS FOUND" << std::endl;
            xpeaks = s->GetPositionX();
            ypeaks = s->GetPositionY();
            double re = 0, im = 0;
            TComplex c;
            for (int numpeak = 0; numpeak < nfound; numpeak++)
            {
                // std::cout << xpeaks[numpeak] << " : " << ypeaks[numpeak] << std::endl;
                std::vector<double> ifftSeg(kernelSz);
                if (ypeaks[numpeak] > peakThreshSTFT)
                {
                    int peakBin = std::floor(xpeaks[numpeak]);
                    std::cout << xpeaks[numpeak] << " : " << ypeaks[numpeak] << " : " << peakBin << std::endl;
                    int tracesLoc = 0;
                    for (int timeval = std::max(peakBin - (kernelSz + 10), 0); timeval < std::min(peakBin + (kernelSz + 10), rangeSTFT); timeval++)
                    {
                        // std::cout << "TAKING IFFT FOR: " << timeval << std::endl;
                        for (int fftiter = 0; fftiter < kernelSz / 2; fftiter++)
                        {
                            re = hSTFT_Re->GetBinContent(timeval, fftiter + 1);
                            im = hSTFT_Im->GetBinContent(timeval, fftiter + 1);
                            c = (re, im);
                            ifft->SetPointComplex(fftiter, c);
                        }
                        ifft->Transform();
                        ifft->GetPoints(&ifftSeg[0]);
                        for (int iter = 0; iter < kernelSz; iter++)
                        {
                            ifftSeg[iter] /= kernelSz;
                            ifftSeg[iter] *= (gauss[iter] / sumKernel);
                            tracesLoc = std::max(timeval - kernelSz / 2, 0) + iter;
                            tracesPeaks[tracesLoc] += ifftSeg[iter];
                        }
                    }
                }
            }
            // IFFT here gives bipolar pulses (differential of the original pulse) hence integrating it to get unipolar pulse
            double CDFtracesPeaks = 0;
            for (int iter = 0; iter < N; iter++)
            {
                CDFtracesPeaks += tracesPeaks[iter];
                tracesPeaksUni[iter] = -1.0 * CDFtracesPeaks;
            }
            // moving baseline back to zero and pulses become unipolar
            for (int iter = 0; iter < N; iter++)
            {
                if (iter < rangeSTFT)
                {
                    std::vector<double> tracesShort(kernelSz);
                    auto minIt = std::min_element(tracesPeaksUni.begin() + iter, tracesPeaksUni.begin() + iter + kernelSz);
                    double minVal = *minIt;
                    for (int j = 0; j < kernelSz; j++)
                    {
                        tracesShort[j] = gauss[j] * (tracesPeaksUni[iter + j] - minVal) / sumKernel;
                        // tracesShort[j] = (tracesPeaksUni[iter + j] - minVal) / 100;
                        tracesPeaksUni1[iter + j] += tracesShort[j];
                    }
                }
                else
                    tracesPeaksUni[iter] = 0.;
            }

            canvas->Clear();
            canvas->Divide(2, 1);
            TGraph *graphTraces = nullptr;
            TGraph *graphTracesTH2 = nullptr;
            TGraph *graphTracesPeaksTH2 = nullptr;
            TGraph *graphTracesPeaksTH2Fin = nullptr;
            canvas->cd(1);
            graphTraces = new TGraph(N);
            graphTracesTH2 = new TGraph(N);
            graphTracesPeaksTH2 = new TGraph(N);
            graphTracesPeaksTH2Fin = new TGraph(N);
            for (int i = 0; i < N; ++i)
            {
                graphTraces->SetPoint(i, i, traces[i]);
                if (i < N - kernelSz)
                    // graphTracesTH2->SetPoint(i, i, (tracesMod[i + kernelSz / 2] / 100 + kernelSz / 2));
                    graphTracesTH2->SetPoint(i, i, (tracesMod[i] / 100));
                else
                    graphTracesTH2->SetPoint(i, i, 0);
                graphTracesPeaksTH2->SetPoint(i, i, tracesPeaksUni[i] + kernelSz / 4);
                graphTracesPeaksTH2Fin->SetPoint(i, i, tracesPeaksUni1[i] + kernelSz / 4);
            }
            graphTraces->SetLineColor(kBlue);
            graphTraces->SetLineWidth(2);
            graphTraces->SetTitle("Traces");
            graphTraces->Draw("AL");
            // legend->AddEntry(graphTraces, "Traces", "l");
            double xmin = graphTraces->GetXaxis()->GetXmin();
            double xmax = graphTraces->GetXaxis()->GetXmax();
            TLine *line = new TLine(xmin, kernelSz / 2, xmax, kernelSz / 2);
            line->SetLineColor(kBlack);
            line->SetLineStyle(1);

            // canvas->cd(3);
            // hTracesShort->SetStats(0);
            // hTracesShort->Draw("COLZ");
            // hTracesShort->SetMinimum(mintrsh);
            graphTracesTH2->SetLineColor(kBlack);
            graphTracesTH2->SetLineWidth(2);
            graphTracesTH2->Draw("SAME");
            line->Draw("SAME");

            canvas->cd(2);
            // hSTFT->SetStats(0);
            // hSTFT->SetMinimum(minstft);
            // hSTFT->Draw("COLZ");
            graphTracesPeaksTH2->SetLineColor(kBlack);
            graphTracesPeaksTH2->SetLineWidth(2);
            graphTracesPeaksTH2->Draw("AL");
            graphTracesPeaksTH2Fin->SetLineColor(kGreen);
            graphTracesPeaksTH2Fin->SetLineWidth(2);
            graphTracesPeaksTH2Fin->Draw("SAME");

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