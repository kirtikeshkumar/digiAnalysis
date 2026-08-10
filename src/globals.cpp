#include "globals.h"

namespace digiAnalysis {
UShort_t nSampleBL = 128;
UShort_t smoothBoxSz = 100;
// UShort_t GateStart = 265;    // NaI 500MSPS has *2 ns
// UShort_t GateLenLong = 2500; // NaI 500MSPS has *2 ns
// UShort_t GateLenShort = 350; // NaI 500MSPS has *2 ns
// UShort_t GateStart = 213;    // NaI 500MSPS has *2 ns
// UShort_t GateLenLong = 2000; // NaI 500MSPS has *2 ns
// UShort_t GateLenShort = 275; // NaI 500MSPS has *2 ns

UShort_t GateStart = 900;    // NaI 500MSPS has *2 ns
UShort_t GateLenLong = 2500; // NaI 500MSPS has *2 ns
UShort_t GateLenShort = 220; // NaI 500MSPS has *2 ns
UShort_t GateMeanTime = 350;
UShort_t TriggerTime = 985;
UShort_t noiseLevelPP = 15;

// SPE
// UShort_t GateStart = 75;    // NaI 500MSPS has *2 ns
// UShort_t GateLenLong = 250; // NaI 500MSPS has *2 ns
// UShort_t GateLenShort = 15; // NaI 500MSPS has *2 ns

// UShort_t GateStart = 188;   // 1GSPS
// UShort_t GateLenLong = 200; // 1GSPS
// UShort_t GateLenShort = 65; // 1GSPS

UShort_t PairCoincWindow = 2000; // in ns
UShort_t nSampleMovBL = 64;
double EvalNormFactor =
    20; // 31; // Factor to scale Eval energy with Digitizer energy
double BLError = 0.6;
} // namespace digiAnalysis