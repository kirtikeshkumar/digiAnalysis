#include "globals.h"

namespace digiAnalysis {
UShort_t nSampleBL = 16;
UShort_t smoothBoxSz = 16;
// UShort_t GateStart = 265;    // NaI 500MSPS has *2 ns
// UShort_t GateLenLong = 2500; // NaI 500MSPS has *2 ns
// UShort_t GateLenShort = 350; // NaI 500MSPS has *2 ns
// UShort_t GateStart = 213;    // NaI 500MSPS has *2 ns
// UShort_t GateLenLong = 2000; // NaI 500MSPS has *2 ns
// UShort_t GateLenShort = 275; // NaI 500MSPS has *2 ns

// UShort_t GateStart = 280;      // NaI 500MSPS has *2 ns
// UShort_t GateLenLong = 1125;   // NaI 500MSPS has *2 ns
// UShort_t GateLenShort = 150;   // NaI 500MSPS has *2 ns

// SPE
UShort_t GateStart = 75;    // NaI 500MSPS has *2 ns
UShort_t GateLenLong = 250; // NaI 500MSPS has *2 ns
UShort_t GateLenShort = 15; // NaI 500MSPS has *2 ns

UShort_t PairCoincWindow = 96; // in ns
UShort_t nSampleMovBL = 16;
double EvalNormFactor =
    31; // 31; // Factor to scale Eval energy with Digitizer energy
double BLError = 0.15;
} // namespace digiAnalysis