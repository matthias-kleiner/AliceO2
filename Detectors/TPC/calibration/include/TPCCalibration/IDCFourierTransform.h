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

#ifndef ALICEO2_IDCFOURIERTRANSFORM_H_
#define ALICEO2_IDCFOURIERTRANSFORM_H_

#include <vector>

#include "Rtypes.h"
#include "DataFormatsTPC/Defs.h"
#include "Framework/Logger.h"
#include "TPCCalibration/IDCFactorizationContainer.h"

namespace o2
{
namespace tpc
{
class IDCFourierTransform
{
 public:
  IDCFourierTransform(const unsigned int rangeIntegrationIntervals, const unsigned int shift, const unsigned int nFourierCoefficients) : mRangeIntegrationIntervals{rangeIntegrationIntervals}, mShift{shift}, mNFourierCoefficients{nFourierCoefficients} {};

  /// adding default constructor for ROOT I/O
  IDCFourierTransform() = default;

  unsigned int getRangeIntegrationIntervals() const { return mRangeIntegrationIntervals; }                         ///< \return returns get number of integration intervals which are used for calculation of the fourier coefficients
  unsigned int getShift() const { return mShift; }                                                                 ///< \return returns shifting parameter for the range
  unsigned int getNCoefficients() const { return mNFourierCoefficients; }                                          /// \return returns numbers of stored fourier coeffiecients
  unsigned long getNIDCs(const o2::tpc::Side side) const { return mIDCsOne[side].size(); }                         /// \return returns number of 1D-IDCs
  unsigned long getNIntervals(const o2::tpc::Side side) const { return mFourierCoefficients.getNIntervals(side); } /// \return returns number of intervals for which the coefficients are obtained

  // auto getFourierCoefficients(const o2::tpc::Side side) const { return mFourierCoefficients[side]; }

  auto getFourierCoefficient(const o2::tpc::Side side, const unsigned int interval, const unsigned int coefficient) const { return mFourierCoefficients.getFourierCoefficient(side, interval, coefficient); }

  const auto& getFourierCoefficients() const { return mFourierCoefficients; }

  const auto& getFourierCoefficients(const o2::tpc::Side side) const { return mFourierCoefficients.mFourierCoefficients[side]; }

  void setIDCs(std::vector<float>&& idcs, const o2::tpc::Side side)
  {
    mIDCsOne[side] = std::move(idcs);
    initCoefficients(side);
  }

  void setIDCs(const std::vector<float>& idcs, const o2::tpc::Side side)
  {
    mIDCsOne[side] = idcs;
    initCoefficients(side);
  }

  /// calculate fourier coefficient
  void calcFourierCoefficients(const o2::tpc::Side side);

  /// get IDC0 values from the inverse fourier transform
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

} // namespace tpc
} // namespace o2

#endif
