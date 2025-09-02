#ifndef Global_h
#define Global_h

#include "includes.hh"

namespace digiAnalysis {
extern UShort_t nSampleBL;
extern UShort_t smoothBoxSz;
extern UShort_t GateStart;
extern UShort_t GateLenLong;
extern UShort_t GateLenShort;
extern UShort_t PairCoincWindow; // in ns
extern double EvalNormFactor;
} // namespace digiAnalysis

#endif
