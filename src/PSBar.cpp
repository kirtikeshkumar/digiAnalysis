#include "PSBar.h"
#include "WaveForm.h"
#include "diffWaveFrom.h"
#include "includes.hh"
#include "singleHits.h"

using namespace std;

namespace digiAnalysis {
PSBar::PSBar() {}
PSBar::PSBar(unsigned int bIndex) {}
PSBar::PSBar(const PSBar &sbar) {}
PSBar::PSBar(singleHits *h1, singleHits *h2) {}

// Getters
UInt_t PSBar::GetQNear() { return 0; }
UInt_t PSBar::GetQFar() { return 0; }
double PSBar::GetQMean() { return 0; }
double PSBar::GetQMeanShort() { return 0; }
double PSBar::GetPSD() { return 0; }
Long_t PSBar::GetDelT() const { return 0; }
Long_t PSBar::GetTStampNear() { return 0; }
Long_t PSBar::GetTStampFar() { return 0; }
Long_t PSBar::GetTStampAverage() { return 0; }
double PSBar::GetHitPos() { return 0; }

// Printer
void PSBar::Print() {}

#ifdef WAVES
std::unique_ptr<WaveForm> PSBar::GetWF() { return std::move(WF); }
WaveForm *PSBar::GetWFPtr() { return WF.get(); }
void PSBar::SetWF(const WaveForm &wf) { WF->SetWaveForm(wf); }
#endif

PSBar::~PSBar() {}
} // namespace digiAnalysis