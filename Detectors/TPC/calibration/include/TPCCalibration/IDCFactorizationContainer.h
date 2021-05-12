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
#include <limits>
#include <math.h>
#include <complex>
#include "DataFormatsTPC/Defs.h"
#include "Framework/Logger.h"
#include "TPCCalibration/ParameterIDC.h"

namespace o2
{
namespace tpc
{

/// struct containing the IDC delta values
template <typename DataT>
struct IDCDelta {

  /// set idcDelta for given index
  /// \param idcDelta Delta IDC value which will be set
  /// \param side side of the TPC
  /// \param index index in the storage
  void setValue(const float idcDelta, const o2::tpc::Side side, const unsigned int index) { mIDCDelta[side][index] = compressValue(idcDelta, side); }

  /// set idcDelta ath the end of storage
  /// \param idcDelta Delta IDC value which will be set
  /// \param side side of the TPC
  void emplace_back(const float idcDelta, const o2::tpc::Side side) { mIDCDelta[side].emplace_back(compressValue(idcDelta, side)); }

  /// \return returns converted IDC value from float to new data type
  /// \param idcDelta Delta IDC value which will be set
  /// \param side side of the TPC
  DataT compressValue(const float idcDelta, const o2::tpc::Side side) const
  {
    const static auto& paramIDCGroup = ParameterIDCCompression::Instance();
    return (std::abs(idcDelta) >= paramIDCGroup.MaxIDCDeltaValue) ? static_cast<DataT>(std::copysign(paramIDCGroup.MaxIDCDeltaValue * mFactors[side] + 0.5f, idcDelta)) : static_cast<DataT>(idcDelta * mFactors[side] + std::copysign(0.5f, idcDelta));
  }

  /// \return returns stored Delta IDC value
  /// \param side side of the TPC
  /// \param index index in the storage
  float getValue(const o2::tpc::Side side, const unsigned int index) const { return (static_cast<float>(mIDCDelta[side][index]) / mFactors[side]); }

  /// set compression factor
  /// \param side side of the TPC
  /// \param factor factor which will be used for the compression
  void setFactor(const o2::tpc::Side side, const float factor) { mFactors[side] = factor; }

  std::array<std::vector<DataT>, o2::tpc::SIDES> mIDCDelta{}; ///< \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  std::array<float, o2::tpc::SIDES> mFactors{1.f, 1.f};       ///< compression factors for each TPC side
};

template <>
struct IDCDelta<float> {
  /// set idcDelta for given index
  /// \param idcDelta Delta IDC value which will be set
  /// \param side side of the TPC
  /// \param index index in the storage
  void setValue(const float idcDelta, const o2::tpc::Side side, const unsigned int index) { mIDCDelta[side][index] = idcDelta; }

  /// \return returns stored Delta IDC value
  /// \param side side of the TPC
  /// \param index index in the storage
  float getValue(const o2::tpc::Side side, const unsigned int index) const { return mIDCDelta[side][index]; }

  std::array<std::vector<float>, o2::tpc::SIDES> mIDCDelta{}; ///< \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
};

template <typename DataT>
class IDCDeltaCompressionHelper
{
 public:
  IDCDeltaCompressionHelper() = default;

  /// static method to get the compressed Delta IDCs from uncompressed Delta IDCs
  /// \return returns compressed Delta IDC values
  /// \param idcDeltaUncompressed uncompressed Delta IDC values
  static IDCDelta<DataT> getCompressedIDCs(const IDCDelta<float>& idcDeltaUncompressed)
  {
    IDCDelta<DataT> idcCompressed{};
    compress(idcDeltaUncompressed, idcCompressed, o2::tpc::Side::A);
    compress(idcDeltaUncompressed, idcCompressed, o2::tpc::Side::C);
    return std::move(idcCompressed);
  }

 private:
  static void compress(const IDCDelta<float>& idcDeltaUncompressed, IDCDelta<DataT>& idcCompressed, const o2::tpc::Side side)
  {
    const float factor = getCompressionFactor(idcDeltaUncompressed, side);
    idcCompressed.mIDCDelta[side].reserve(idcDeltaUncompressed.mIDCDelta[side].size());
    idcCompressed.setFactor(side, factor);

    for (auto& idc : idcDeltaUncompressed.mIDCDelta[side]) {
      idcCompressed.emplace_back(idc, side);
    }
  }

  /// \return returns the factor which is used during the compression
  static float getCompressionFactor(const IDCDelta<float>& idcDeltaUncompressed, const o2::tpc::Side side)
  {
    const float maxAbsIDC = getMaxValue(idcDeltaUncompressed.mIDCDelta[side]);
    const auto& paramIDCGroup = ParameterIDCCompression::Instance();
    const float maxIDC = paramIDCGroup.MaxIDCDeltaValue;
    return (maxAbsIDC > maxIDC && maxIDC > 0) ? (std::numeric_limits<DataT>::max() / maxIDC) : (std::numeric_limits<DataT>::max() / maxAbsIDC);
  }

  /// \returns returns maximum abs value in vector
  static float getMaxValue(const std::vector<float>& idcs)
  {
    return std::abs(*std::max_element(idcs.begin(), idcs.end(), [](const int a, const int b) -> bool { return (std::abs(a) < std::abs(b)); }));
  };
};

struct IDCZeroOne { ///< struct containing the IDC0 and IDC1 values

  /// set IDC zero for given index
  /// \param idcZero Delta IDC value which will be set
  /// \param side side of the TPC
  /// \param index index in the storage
  void setValueIDCZero(const float idcZero, const o2::tpc::Side side, const unsigned int index) { mIDCZero[side][index] = idcZero; }

  /// increase IDC zero for given index
  /// \param idcZero Delta IDC value which will be set
  /// \param side side of the TPC
  /// \param index index in the storage
  void fillValueIDCZero(const float idcZero, const o2::tpc::Side side, const unsigned int index) { mIDCZero[side][index] += idcZero; }

  /// \return returns stored IDC zero value
  /// \param side side of the TPC
  /// \param index index in the storage
  float getValueIDCZero(const o2::tpc::Side side, const unsigned int index) const { return mIDCZero[side][index]; }

  /// set IDC one for given index
  /// \param idcOne Delta IDC value which will be set
  /// \param side side of the TPC
  /// \param index index in the storage
  void setValueIDCOne(const float idcOne, const o2::tpc::Side side, const unsigned int index) { mIDCOne[side][index] = idcOne; }

  /// \return returns stored IDC one value
  /// \param side side of the TPC
  /// \param index index in the storage
  float getValueIDCOne(const o2::tpc::Side side, const unsigned int index) const { return mIDCOne[side][index]; }

  std::array<std::vector<float>, o2::tpc::SIDES> mIDCZero{}; ///< I_0(r,\phi) = <I(r,\phi,t)>_t
  std::array<std::vector<float>, o2::tpc::SIDES> mIDCOne{};  ///< I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
};

struct FourierCoeff { ///< struct containing the fourier coefficients calculated from IDC0
  unsigned long getNIntervals(const o2::tpc::Side side) const { return mFourierCoefficients[side].size(); } /// \return returns numbers of stored intervals
  unsigned long getNCoefficients(const unsigned long interval) const { return mFourierCoefficients[interval].size(); } /// \return returns numbers of stored fourier coefficients
  auto getFourierCoefficient(const o2::tpc::Side side, const unsigned int interval, const unsigned int coefficient) const { return mFourierCoefficients[side][interval][coefficient]; }

  /// \return returns the stored value
  /// \param side side of the TPC
  /// \param interval interval of which the coefficients were calculated
  /// \param coefficient index of coefficient
  const std::complex<float>& operator()(const o2::tpc::Side side, unsigned long interval, unsigned int coefficient) const { return mFourierCoefficients[side][interval][coefficient]; }

  /// \return returns the stored value
  /// \param side side of the TPC
  /// \param interval interval of which the coefficients were calculated
  /// \param coefficient index of coefficient
  std::complex<float>& operator()(const o2::tpc::Side side, const unsigned long interval, const unsigned int coefficient) { return mFourierCoefficients[side][interval][coefficient]; }

  std::array<std::vector<std::vector<std::complex<float>>>, o2::tpc::SIDES> mFourierCoefficients; ///< fourier coefficients. side -> interval -> coefficient
};

} // namespace tpc
} // namespace o2

#endif
