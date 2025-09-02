#ifndef Pair_h
#define Pair_h

#include "RtypesCore.h"
#include "WaveForm.h"
#include "includes.hh"
#include "singleHits.h"
#pragma once
#include <iomanip>
#include <iostream>
#include <variant>

namespace digiAnalysis {
class singleHits;
class WaveForm;

class Pair {
private:
  singleHits *hit1;
  singleHits *hit2;
  UShort_t SumEnergy;
  ULong64_t PairDelTime;

public:
  Pair();
  Pair(const singleHits *hit1, const singleHits *hit2);
  Pair(const Pair &other);
  ~Pair();

  /**
   * @brief Set the pair of hits.
   *
   * Assigns the provided hits to the pair.
   *
   * @param hit1 Pointer to the first singleHits object.
   * @param hit2 Pointer to the second singleHits object.
   */
  void SetPair(const singleHits *hit1, const singleHits *hit2);

  /**
   * @brief Clear the pair of hits.
   *
   * Resets the pair, removing references to the hits.
   */
  void ClearPair();

  /**
   * @brief Retrieve a HIT based on a channel selector.
   *
   * Special values of `Ch` are interpreted as follows:
   *   - `-1` : The first hit in time
   *   - `-2` : The second hit in time
   *   - `0`  : The hit with the lower channel number
   *   - `1`  : The hit with the higher channel number
   *
   * @param Ch  Selector for which hit to return.
   * @return The corresponding singleHits object.
   */
  singleHits *GetHit(Short_t Ch);

  /**
   * @brief Retrieve hit CHANNEL NUMBER based on a selector.
   *
   * Special values of `Sel` are interpreted as follows:
   *   - `-1` : The first hit in time
   *   - `-2` : The second hit in time
   *   - `0`  : The hit with the lower channel number
   *   - `1`  : The hit with the higher channel number
   *
   * @param Sel  Selector for which hit to return.
   * @return The corresponding Hit Channel Number.
   */
  UShort_t GetPairHitCh(Short_t Sel);

  /**
   * @brief Retrieve hit TIME based on a selector.
   *
   * Special values of `Sel` are interpreted as follows:
   *   - `-1` : The first hit in time
   *   - `-2` : The second hit in time
   *   - `0`  : The hit with the lower channel number
   *   - `1`  : The hit with the higher channel number
   *
   * @param Sel  Selector for which hit to return.
   * @return The corresponding hit Time.
   */
  ULong64_t GetPairHitTime(Short_t Sel);

  /**
   * @brief Retrieve pair DEL TIME
   * @return The corresponding pair |Time Difference|.
   */
  ULong64_t GetPairDelTime();

  /**
   * @brief Retrieve hit ENERGY based on a selector.
   *
   * Special values of `Sel` are interpreted as follows:
   *   - `-1` : The first hit in time
   *   - `-2` : The second hit in time
   *   - `0`  : The hit with the lower channel number
   *   - `1`  : The hit with the higher channel number
   *
   * @param Sel  Selector for which hit to return.
   * @return The corresponding Hit Energy.
   */
  UShort_t GetPairHitEnergy(Short_t Ch);

  /**
   * @brief Retrieve pair ENERGY.
   * @return The corresponding Pair Energy.
   */
  UShort_t GetPairEnergy();

  /**
   * @brief Retrieve pair ENERGY.
   * @return The corresponding Pair Energy.
   *
   * @param fac1 Scaling factor for the lower channels energy.
   * @param fac2 Scaling factor for the higher channels energy.
   */
  double GetPairEnergy(double fac1, double fac2);

  /**
   * @brief Retrieve hit ENERGY SHORT based on a selector.
   *
   * Special values of `Sel` are interpreted as follows:
   *   - `-1` : The first hit in time
   *   - `-2` : The second hit in time
   *   - `0`  : The hit with the lower channel number
   *   - `1`  : The hit with the higher channel number
   *
   * @param Sel  Selector for which hit to return.
   * @return The corresponding Hit Energy Short.
   */
  UShort_t GetPairHitEnergyShort(Short_t Ch);

  /**
   * @brief Retrieve hit PSD based on a selector.
   *
   * Special values of `Sel` are interpreted as follows:
   *   - `-1` : The first hit in time
   *   - `-2` : The second hit in time
   *   - `0`  : The hit with the lower channel number
   *   - `1`  : The hit with the higher channel number
   *
   * @param Sel  Selector for which hit to return.
   * @return The corresponding Hit PSD.
   */
  UShort_t GetPairHitPSD(Short_t Ch);

#ifdef WAVES
  /**
   * @brief Retrieve hit EVALUATED ENERGY based on a selector.
   *
   * Special values of `Sel` are interpreted as follows:
   *   - `-1` : The first hit in time
   *   - `-2` : The second hit in time
   *   - `0`  : The hit with the lower channel number
   *   - `1`  : The hit with the higher channel number
   *
   * @param Sel  Selector for which hit to return.
   * @return The corresponding Hit Evaluated Energy.
   */
  double GetPairHitEvalEnergy(Short_t Ch);

  /**
   * @brief Retrieve hit EVALUATED ENERGY SHORT based on a selector.
   *
   * Special values of `Sel` are interpreted as follows:
   *   - `-1` : The first hit in time
   *   - `-2` : The second hit in time
   *   - `0`  : The hit with the lower channel number
   *   - `1`  : The hit with the higher channel number
   *
   * @param Sel  Selector for which hit to return.
   * @return The corresponding Hit Evaluated Energy Short.
   */
  double GetPairHitEvalEnergyShort(Short_t Ch);

  /**
   * @brief Retrieve hit EVALUATED PSD based on a selector.
   *
   * Special values of `Sel` are interpreted as follows:
   *   - `-1` : The first hit in time
   *   - `-2` : The second hit in time
   *   - `0`  : The hit with the lower channel number
   *   - `1`  : The hit with the higher channel number
   *
   * @param Sel  Selector for which hit to return.
   * @return The corresponding Hit Evaluated PSD.
   */
  double GetPairHitEvalPSD(Short_t Ch);
#endif
  void Print();
};
} // namespace digiAnalysis
#endif