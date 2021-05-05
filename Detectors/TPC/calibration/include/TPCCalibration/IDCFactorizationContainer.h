// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file IDCFactorizationContainer.h
/// \brief This file provides the structs for storing the factorized IDC values to be stored in the CCDB
///
/// \author  Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>
/// \date Apr 30, 2021

#ifndef ALICEO2_TPC_IDCFACTORIZATIONCONTAINER_H_
#define ALICEO2_TPC_IDCFACTORIZATIONCONTAINER_H_

#include <array>
#include <vector>
#include "DataFormatsTPC/Defs.h"

namespace o2
{
namespace tpc
{

template <typename DataT = float>
struct IDCDelta { ///< struct containing the IDC delta values
  void setValue(const float value, const o2::tpc::Side side, const unsigned int index) { mIDCDelta[side][index] = value; }
  const float getValue(const o2::tpc::Side side, const unsigned int index) const { return mIDCDelta[side][index]; }
  float& getValue(const o2::tpc::Side side, const unsigned int index) { return mIDCDelta[side][index]; }

  std::array<std::vector<DataT>, o2::tpc::SIDES> mIDCDelta{}; ///< \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
};

struct IDCZeroOne { ///< struct containing the IDC0 and IDC1 values
  void setValueIDCZero(const float value, const o2::tpc::Side side, const unsigned int index) { mIDCZero[side][index] = value; }
  const float getValueIDCZero(const o2::tpc::Side side, const unsigned int index) const { return mIDCZero[side][index]; }
  float& getValueIDCZero(const o2::tpc::Side side, const unsigned int index) { return mIDCZero[side][index]; }

  void setValueIDCOne(const float value, const o2::tpc::Side side, const unsigned int index) { mIDCOne[side][index] = value; }
  const float getValueIDCOne(const o2::tpc::Side side, const unsigned int index) const { return mIDCOne[side][index]; }
  float& getValueIDCOne(const o2::tpc::Side side, const unsigned int index) { return mIDCOne[side][index]; }

  std::array<std::vector<float>, o2::tpc::SIDES> mIDCZero{}; ///< I_0(r,\phi) = <I(r,\phi,t)>_t
  std::array<std::vector<float>, o2::tpc::SIDES> mIDCOne{};  ///< I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
};

} // namespace tpc
} // namespace o2

#endif
