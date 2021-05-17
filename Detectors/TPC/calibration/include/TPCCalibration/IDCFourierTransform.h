// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file IDCFourierTransform.h
/// \brief class for calculating the fourier coefficients from 1D-IDCs
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>
/// \date May 11, 2021

#ifndef ALICEO2_IDCFOURIERTRANSFORM_H_
#define ALICEO2_IDCFOURIERTRANSFORM_H_

#include <vector>
#include "Rtypes.h"
#include "DataFormatsTPC/Defs.h"
#include "Framework/Logger.h"
#include "TPCCalibration/IDCFactorizationContainer.h"

namespace o2::tpc
{

/// class for fourier transform of 1D-IDCs

class IDCFourierTransform
{
 public:
  /// contructor
  /// \param rangeIntegrationIntervals number of IDCs for each interval which will be used to calculate the fourier coefficients
  /// \param shift number of IDCs of which the range will be shifted
  /// \param nFourierCoefficients number of fourier coefficients which will be calculated (cant be larger than rangeIntegrationIntervals)
  IDCFourierTransform(const unsigned int rangeIntegrationIntervals = 10, const unsigned int shift = 1, const unsigned int nFourierCoefficients = 10) : mRangeIntegrationIntervals{rangeIntegrationIntervals}, mShift{shift}, mNFourierCoefficients{nFourierCoefficients} {};

  /// \return returns number of IDCs for each interval which will be used to calculate the fourier coefficients
  unsigned int getRangeIntegrationIntervals() const { return mRangeIntegrationIntervals; }

  /// \return returns shifting parameter for the range
  unsigned int getShift() const { return mShift; }

  /// \return returns numbers of stored fourier coeffiecients
  unsigned int getNCoefficients() const { return mNFourierCoefficients; }

  /// \return returns number of 1D-IDCs
  /// \param side TPC side
  unsigned long getNIDCs(const o2::tpc::Side side) const { return mIDCsOne[side].size(); }

  /// \return returns number of intervals for which the coefficients are obtained
  /// \param side TPC side
  unsigned long getNIntervals(const o2::tpc::Side side) const { return mFourierCoefficients.getNValues(side) / mNFourierCoefficients; }

  /// \return returns foruier coefficient
  /// \param side TPC side
  /// \param interval index of interval
  /// \param coefficient index of coefficient
  /// \param type whether to get real or imag coefficient
  auto getFourierCoefficient(const o2::tpc::Side side, const unsigned int interval, const unsigned int coefficient, const FourierCoeff::CoeffType type) const { return mFourierCoefficients(side, getIndex(interval, coefficient), type); }

  /// \return returns struct holding all fourier coefficients
  const auto& getFourierCoefficients() const { return mFourierCoefficients; }

  /// \return returns vector of all fourier coefficients for given side and type
  /// \param side TPC side
  /// \param type whether to get real or imag coefficient
  const auto& getFourierCoefficients(const o2::tpc::Side side, const FourierCoeff::CoeffType type) const { return mFourierCoefficients.getFourierCoefficients(side, type); }

  /// \return returns index to fourier coefficient
  /// \param interval index of interval
  /// \param coefficient index of coefficient
  unsigned int getIndex(const unsigned int interval, const unsigned int coefficient) const { return interval * mNFourierCoefficients + coefficient; }

  /// set input 1D-IDCs which are used to calculate fourier coefficients
  /// \param idcsOne 1D-IDCs
  /// \param side TPC side
  void setIDCs(std::vector<float>&& idcsOne, const o2::tpc::Side side);

  /// set input 1D-IDCs which are used to calculate fourier coefficients
  /// \param idcsOne 1D-IDCs
  /// \param side TPC side
  void setIDCs(const std::vector<float>& idcsOne, const o2::tpc::Side side);

  /// calculate fourier coefficients
  /// \param side TPC side
  void calcFourierCoefficients(const o2::tpc::Side side);

  /// get IDC0 values from the inverse fourier transform. Can be used for debugging. std::vector<std::vector<float>>: first vector interval second vector IDC0 values
  /// \param side TPC side
  std::vector<std::vector<float>> inverseFourierTransform(const o2::tpc::Side side) const;

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName = "Fourier.root", const char* outName = "FourierCoefficients") const;

  /// create debug tree
  void dumpToTree() const;

 private:
  const unsigned int mRangeIntegrationIntervals{};           ///< number of IDCs used for the calculation of fourier coefficients
  const unsigned int mShift{};                               ///< shifting parameter for the range. Should be >=1
  std::array<std::vector<float>, o2::tpc::SIDES> mIDCsOne{}; ///< all 1D-IDCs which are used to calculate the fourier coefficients
  const unsigned int mNFourierCoefficients{};                ///< number of fourier coefficients which will be calculated
  FourierCoeff mFourierCoefficients;                         ///< fourier coefficients. side -> interval -> coefficient

  /// initialize mFourierCoefficients member
  void initCoefficients(const o2::tpc::Side side);

  unsigned int getLastIndex(unsigned long interval, const o2::tpc::Side side) const { return (interval == getNIntervals(side) - 1) ? getNIDCs(side) - mShift * interval : mRangeIntegrationIntervals; }

  ClassDefNV(IDCFourierTransform, 1)
};

} // namespace o2::tpc

#endif
