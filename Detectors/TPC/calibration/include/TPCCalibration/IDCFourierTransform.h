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
// #include <fftw3.h>

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
  /// \param rangeIDC number of IDCs for each interval which will be used to calculate the fourier coefficients
  /// \param nFourierCoefficients number of fourier coefficients which will be calculated (cant be larger than rangeIDC)
  /// \param timeFrames number of time frames which will be stored
  IDCFourierTransform(const unsigned int rangeIDC = 10, const unsigned int nFourierCoefficients = 10, const unsigned int timeFrames = 1) : mRangeIDC{rangeIDC}, mNFourierCoefficients{nFourierCoefficients}, mTimeFrames{timeFrames} {};

  /// \return returns number of IDCs for each interval which will be used to calculate the fourier coefficients
  unsigned int getrangeIDC() const { return mRangeIDC; }

  /// \return returns numbers of stored fourier coeffiecients
  unsigned int getNCoefficients() const { return mNFourierCoefficients; }

  /// \return returns number of 1D-IDCs
  /// \param side TPC side
  unsigned long getNIDCs(const o2::tpc::Side side) const { return mIDCOne[!mBufferIndex].mIDCOne[side].size(); }

  /// \return returns number of intervals for which the coefficients are obtained
  unsigned int getNIntervals() const { return mTimeFrames; }

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

  /// \return returns vector of stored IDC1 I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  /// \param side TPC side
  const std::vector<float>& getIDCOne(const o2::tpc::Side side) const { return mIDCOne[!mBufferIndex].mIDCOne[side]; }

  /// \return returns struct of stored IDC1 I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  const auto& getIDCOne() const { return mIDCOne[!mBufferIndex]; }

  /// set input 1D-IDCs which are used to calculate fourier coefficients
  /// \param idcsOne 1D-IDCs
  /// \param integrationIntervalsPerTF vector containg for each TF the number of IDCs
  void setIDCs(IDCOne&& idcsOne, std::vector<unsigned int>&& integrationIntervalsPerTF)
  {
    mIDCOne[mBufferIndex] = std::move(idcsOne);
    mIntegrationIntervalsPerTF[mBufferIndex] = std::move(integrationIntervalsPerTF);
  }

  /// set input 1D-IDCs which are used to calculate fourier coefficients
  /// \param idcsOne 1D-IDCs
  /// \param integrationIntervalsPerTF vector containg for each TF the number of IDCs
  void setIDCs(const IDCOne& idcsOne, const std::vector<unsigned int>& integrationIntervalsPerTF)
  {
    mIDCOne[mBufferIndex] = idcsOne;
    mIntegrationIntervalsPerTF[mBufferIndex] = integrationIntervalsPerTF;
  }

  /// calculate fourier coefficients
  void calcFourierCoefficients();

  /// get IDC0 values from the inverse fourier transform. Can be used for debugging. std::vector<std::vector<float>>: first vector interval second vector IDC0 values
  /// \param side TPC side
  std::vector<std::vector<float>> inverseFourierTransform(const o2::tpc::Side side) const;

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName = "Fourier.root", const char* outName = "FourierCoefficients") const;

  /// create debug tree
  void dumpToTree(const char* outFileName = "FourierTree.root") const;

 private:
  const unsigned int mRangeIDC{};                                         ///< number of IDCs used for the calculation of fourier coefficients
  const unsigned int mNFourierCoefficients{};                             ///< number of fourier coefficients which will be calculated
  const unsigned int mTimeFrames{};                                       ///< number of timeframes which for which teh fourier coefficients are stored
  std::array<IDCOne, 2> mIDCOne{};                                        ///< all 1D-IDCs which are used to calculate the fourier coefficients. A buffer for the last aggregation interval is used to calculate the fourier coefficients for the first TFs
  std::array<std::vector<unsigned int>, 2> mIntegrationIntervalsPerTF{};  ///< number of integration intervals per TF used to set the correct range of IDCs. A buffer is needed for the last aggregation interval.
  bool mBufferIndex{false};                                               ///< index for the buffer
  FourierCoeff mFourierCoefficients{mTimeFrames * mNFourierCoefficients}; ///< fourier coefficients. side -> interval -> coefficient

  /// calculate fourier coefficients
  /// \param side TPC side
  void calcFourierCoefficients(const o2::tpc::Side side);

  /// returns number of IDCs which are used to calculate the fourier coefficients. This value can be smaller than the specified mRangeIDC if the buffer from last aggregation interval is empty
  /// \param currentIndexEndInterval index of integration
  /// \param buffer buffer to IDCs
  int getRangeIDC(const int currentIndexEndInterval, const o2::tpc::Side side, const bool buffer) const { return mIDCOne[!buffer].empty(side) ? ((currentIndexEndInterval + 1 - static_cast<int>(mRangeIDC) < 0) ? currentIndexEndInterval + 1 : mRangeIDC) : mRangeIDC; }

  /// returns 1D-IDC
  float getIDCOne(const int index, const o2::tpc::Side side, const bool buffer) const { return (index < 0) ? mIDCOne[!buffer].getValueIDCOne(side, mIDCOne[!buffer].getNIDCs(side) + index) : mIDCOne[buffer].getValueIDCOne(side, index); }

  ClassDefNV(IDCFourierTransform, 1)
};

} // namespace o2::tpc

#endif
