#include "globals.h"

namespace digiAnalysis
{
    UShort_t nSampleBL = 64;
    UShort_t smoothBoxSz = 16;
    UShort_t GateStart = 290;      // NaI 500MSPS has *2 ns
    UShort_t GateLenLong = 1110;   // NaI 500MSPS has *2 ns
    UShort_t GateLenShort = 275;   // NaI 500MSPS has *2 ns
    UShort_t PairCoincWindow = 96; // in ns
    float EvalNormFactor =
        31.4; // Factor to scale Eval energy with Digitizer energy
} // namespace digiAnalysis