#include "globals.h"

namespace digiAnalysis {
UShort_t nSampleBL = 64;
UShort_t smoothBoxSz = 80;
UShort_t GateStart = 213;      // NaI 500MSPS has *2 ns
UShort_t GateLenLong = 2000;   // NaI 500MSPS has *2 ns
UShort_t GateLenShort = 275;   // NaI 500MSPS has *2 ns
UShort_t PairCoincWindow = 96; // in ns
UShort_t nSampleMovBL = 16;
double EvalNormFactor = 64; // Factor to scale Eval energy with Digitizer energy
double BLError = 0.15;
} // namespace digiAnalysis