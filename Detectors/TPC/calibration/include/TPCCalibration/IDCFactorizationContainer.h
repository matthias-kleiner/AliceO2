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
struct IDCDeltaContainer {
  std::array<std::vector<DataT>, o2::tpc::SIDES> mIDCDelta{}; ///< \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
};

/// storage for the factor used to compress IDCDelta.
/// This factor is separated from the IDCDelta struct to able to store those values independently in the CCDB
struct IDCDeltaCompressionFactors {
  std::array<float, o2::tpc::SIDES> mFactors{1.f, 1.f}; ///< compression factors for each TPC side
};

/// struct to access and set Delta IDCs
template <typename DataT>
struct IDCDelta {

  /// set idcDelta for given index
  /// \param idcDelta Delta IDC value which will be set
  /// \param side side of the TPC
  /// \param index index in the storage
  void setValue(const float idcDelta, const o2::tpc::Side side, const unsigned int index) { mIDCDelta.mIDCDelta[side][index] = compressValue(idcDelta, side); }

  /// set idcDelta ath the end of storage
  /// \param idcDelta Delta IDC value which will be set
  /// \param side side of the TPC
  void emplace_back(const float idcDelta, const o2::tpc::Side side) { mIDCDelta.mIDCDelta[side].emplace_back(compressValue(idcDelta, side)); }

  /// \return returns converted IDC value from float to new data type
  /// \param idcDelta Delta IDC value which will be set
  /// \param side side of the TPC
  DataT compressValue(const float idcDelta, const o2::tpc::Side side) const
  {
    const static auto& paramIDCGroup = ParameterIDCCompression::Instance();
    return (std::abs(idcDelta) >= paramIDCGroup.MaxIDCDeltaValue) ? static_cast<DataT>(std::copysign(paramIDCGroup.MaxIDCDeltaValue * mCompressionFactor.mFactors[side] + 0.5f, idcDelta)) : static_cast<DataT>(idcDelta * mCompressionFactor.mFactors[side] + std::copysign(0.5f, idcDelta));
  }

  /// \return returns stored Delta IDC value
  /// \param side side of the TPC
  /// \param index index in the storage
  float getValue(const o2::tpc::Side side, const unsigned int index) const { return (static_cast<float>(mIDCDelta.mIDCDelta[side][index]) / mCompressionFactor.mFactors[side]); }

  /// set compression factor
  /// \param side side of the TPC
  /// \param factor factor which will be used for the compression
  void setFactor(const o2::tpc::Side side, const float factor) { mCompressionFactor.mFactors[side] = factor; }

  /// \return returns vector of Delta IDCs for given side
  /// \param side side of the TPC
  const auto& getIDCDelta(const o2::tpc::Side side) const { return mIDCDelta.mIDCDelta[side]; }

  /// \return returns vector of Delta IDCs for given side
  /// \param side side of the TPC
  auto& getIDCDelta(const o2::tpc::Side side) { return mIDCDelta.mIDCDelta[side]; }

  /// \return returns IDCDeltaContainer container
  const auto& getIDCDelta() const { return mIDCDelta; }

  /// \return returns IDCDeltaContainer container
  auto& getIDCDelta() { return mIDCDelta; }

  /// \return returns compression factors to uncompress Delta IDC
  const auto& getCompressionFactors() const { return mCompressionFactor; }

  /// \return returns compression factors to uncompress Delta IDC
  /// \param side side of the TPC
  auto getCompressionFactor(const o2::tpc::Side side) const { return mCompressionFactor.mFactors[side]; }

  IDCDeltaContainer<DataT> mIDCDelta{};            ///< storage for Delta IDCs
  IDCDeltaCompressionFactors mCompressionFactor{}; ///< compression factor for Delta IDCs
};

template <>
struct IDCDelta<float> {
  /// set idcDelta for given index
  /// \param idcDelta Delta IDC value which will be set
  /// \param side side of the TPC
  /// \param index index in the storage
  void setValue(const float idcDelta, const o2::tpc::Side side, const unsigned int index) { mIDCDelta.mIDCDelta[side][index] = idcDelta; }

  /// \return returns stored Delta IDC value
  /// \param side side of the TPC
  /// \param index index in the storage
  float getValue(const o2::tpc::Side side, const unsigned int index) const { return mIDCDelta.mIDCDelta[side][index]; }

  /// \return returns vector of Delta IDCs for given side
  /// \param side side of the TPC
  const auto& getIDCDelta(const o2::tpc::Side side) const { return mIDCDelta.mIDCDelta[side]; }

  /// \return returns vector of Delta IDCs for given side
  /// \param side side of the TPC
  auto& getIDCDelta(const o2::tpc::Side side) { return mIDCDelta.mIDCDelta[side]; }

  IDCDeltaContainer<float> mIDCDelta{}; ///< storage for Delta IDCs
};

/// helper class to compress Delta IDC values
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
    idcCompressed.getIDCDelta(side).reserve(idcDeltaUncompressed.getIDCDelta(side).size());
    idcCompressed.setFactor(side, factor);

    for (auto& idc : idcDeltaUncompressed.getIDCDelta(side)) {
      idcCompressed.emplace_back(idc, side);
    }
  }

  /// \return returns the factor which is used during the compression
  static float getCompressionFactor(const IDCDelta<float>& idcDeltaUncompressed, const o2::tpc::Side side)
  {
    const float maxAbsIDC = getMaxValue(idcDeltaUncompressed.getIDCDelta(side));
    const auto& paramIDCGroup = ParameterIDCCompression::Instance();
    const float maxIDC = paramIDCGroup.MaxIDCDeltaValue;
    // TODO catch division by zero (maxAbsIDC=0)?
    return (maxAbsIDC > maxIDC && maxIDC > 0) ? (std::numeric_limits<DataT>::max() / maxIDC) : (std::numeric_limits<DataT>::max() / maxAbsIDC);
  }

  /// \returns returns maximum abs value in vector
  static float getMaxValue(const std::vector<float>& idcs)
  {
    return std::abs(*std::max_element(idcs.begin(), idcs.end(), [](const int a, const int b) -> bool { return (std::abs(a) < std::abs(b)); }));
  };
};

///<struct containing the IDC0 and IDC1 values
struct IDCZeroOne {

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

/// struct containing the fourier coefficients calculated from IDC0
struct FourierCoeff {
  /// IDC types
  enum class CoeffType { REAL = 0, ///< real coefficient
                         IMAG = 1  ///< imag coefficient
  };

  /// \return returns total number of stored coefficients for given side and real/complex type
  /// \param side side
  /// \param type coefficient type CoeffType::REAL or CoeffType::IMAG
  unsigned long getNValues(const o2::tpc::Side side, const CoeffType type = CoeffType::REAL) const { return mFourierCoefficients[side][static_cast<unsigned int>(type)].size(); }

  /// \return returns all stored coefficients for given side and real/complex type
  /// \param side side
  /// \param type coefficient type CoeffType::REAL or CoeffType::IMAG
  const auto& getFourierCoefficients(const o2::tpc::Side side, const FourierCoeff::CoeffType type) const { return mFourierCoefficients[side][static_cast<int>(type)]; }

  /// \return returns the stored value
  /// \param side side of the TPC
  /// \param index index of the data. For calculation see IDCFourierTransform::getIndex()
  /// \param CoeffType real or imag coefficient
  float operator()(const o2::tpc::Side side, unsigned int index, const CoeffType type) const { return mFourierCoefficients[side][static_cast<unsigned int>(type)][index]; }

  /// \return returns the stored value
  /// \param side side of the TPC
  /// \param index index of the data. TFor calculation see IDCFourierTransform::getIndex()
  /// \param CoeffType real or imag coefficient
  float& operator()(const o2::tpc::Side side, unsigned int index, const CoeffType type) { return mFourierCoefficients[side][static_cast<unsigned int>(type)][index]; }

  unsigned long getNTotalCoefficients(const o2::tpc::Side side, const CoeffType type) const { return mFourierCoefficients[side][static_cast<unsigned int>(type)].size(); }

  std::array<std::array<std::vector<float>, 2>, o2::tpc::SIDES> mFourierCoefficients{}; ///< fourier coefficients. side -> real/imag -> coefficient
};

/// struct for storing the number of fourier coefficients per interval to access the fourier ceofficients
struct FourierCoeffParameters {

  /// constructor
  /// nCoefficientsPerInterval number of fourier coefficientes per interval
  FourierCoeffParameters(const unsigned int nCoefficientsPerInterval = 50) : mNCoefficientsPerInterval{static_cast<unsigned char>(nCoefficientsPerInterval)} {};

  /// \return returns number of fourier coefficients per interval
  unsigned int getNCoefficientsPerInterval() const { return static_cast<unsigned int>(mNCoefficientsPerInterval); }

  const unsigned char mNCoefficientsPerInterval{}; ///< number of real/imag fourier coefficients per integration interval
};

} // namespace tpc
} // namespace o2

#endif
