// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
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
#include "TPCCalibration/IDCContainer.h"
#include "TPCCalibration/IDCFourierTransformBase.h"

namespace o2::tpc
{

/// class for fourier transform of 1D-IDCs
/// For example usage see testO2TPCIDCFourierTransform.cxx

/// \tparam Type type which can either be IDCFTType::IDCFourierTransformBaseEPN for synchronous reconstruction or IDCFTType::IDCFourierTransformBaseAggregator for aggregator
template <class Type> // do not use enum class as type to avoid problems with ROOT dictionary generation!
class IDCFourierTransform : public IDCFourierTransformBase<Type>
{
 public:
  /// constructor for IDCFTType::AGGREGATOR type
  /// \param rangeIDC number of IDCs for each interval which will be used to calculate the fourier coefficients
  /// \param timeFrames number of time frames which will be stored
  /// \param nFourierCoefficientsStore number of courier coefficients (real+imag) which will be stored (the maximum can be 'rangeIDC + 2', should be an even number when using naive FT). If less than maximum is setn the inverse fourier transform will not work.
  template <bool IsEnabled = true, typename std::enable_if<(IsEnabled && (std::is_same<Type, IDCFTType::IDCFourierTransformBaseAggregator>::value)), int>::type = 0>
  IDCFourierTransform(const unsigned int rangeIDC = 200, const unsigned int timeFrames = 2000, const unsigned int nFourierCoefficientsStore = 200 + 2) : IDCFTType::IDCFourierTransformAggregator(rangeIDC, timeFrames), mFourierCoefficients{timeFrames, nFourierCoefficientsStore} {};

  /// constructor for IDCFTType::EPN type
  /// \param rangeIDC number of IDCs for each interval which will be used to calculate the fourier coefficients
  /// \param nFourierCoefficientsStore number of courier coefficients (real+imag) which will be stored (the maximum can be 'rangeIDC + 2', should be an even number when using naive FT). If less than maximum is setn the inverse fourier transform will not work.
  template <bool IsEnabled = true, typename std::enable_if<(IsEnabled && (std::is_same<Type, IDCFTType::IDCFourierTransformBaseEPN>::value)), int>::type = 0>
  IDCFourierTransform(const unsigned int rangeIDC = 200, const unsigned int nFourierCoefficientsStore = 200 + 2) : IDCFTType::IDCFourierTransformEPN(rangeIDC), mFourierCoefficients{1, nFourierCoefficientsStore} {};

  /// set fast fourier transform using FFTW3
  /// \param fft use FFTW3 or not (naive approach)
  static void setFFT(const bool fft) { sFftw = fft; }

  /// \param nThreads set the number of threads used for calculation of the fourier coefficients
  static void setNThreads(const int nThreads) { sNThreads = nThreads; }

  /// calculate fourier coefficients
  void calcFourierCoefficients() { sFftw ? calcFourierCoefficientsFFTW3() : calcFourierCoefficientsNaive(); }

  /// get IDC0 values from the inverse fourier transform. Can be used for debugging. std::vector<std::vector<float>>: first vector interval second vector IDC0 values
  /// \param side TPC side
  std::vector<std::vector<float>> inverseFourierTransform(const o2::tpc::Side side) const { return sFftw ? inverseFourierTransformFFTW3(side) : inverseFourierTransformNaive(side); }

  /// \return returns number of IDCs for each interval which will be used to calculate the fourier coefficients
  unsigned int getrangeIDC() const { return this->mRangeIDC; }

  /// \return returns struct holding all fourier coefficients
  const auto& getFourierCoefficients() const { return mFourierCoefficients; }

  /// get type of used fourier transform
  static bool getFFT() { return sFftw; }

  /// get the number of threads used for calculation of the fourier coefficients
  static int getNThreads() { return sNThreads; }

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName = "Fourier.root", const char* outName = "FourierCoefficients") const;

  /// create debug tree
  /// \param outFileName name of the output tree
  void dumpToTree(const char* outFileName = "FourierTree.root") const;

  /// printing information about the algorithms which are used by FFTW for debugging e.g. seeing if SIMD instructions will be used
  void printFFTWPlan() const;

 private:
  FourierCoeff mFourierCoefficients; ///< fourier coefficients. side -> interval -> coefficient
  inline static int sFftw{1};        ///< using fftw or naive approach for calculation of fourier coefficients
  inline static int sNThreads{1};    ///< number of threads which are used during the calculation of the fourier coefficients

  /// calculate fourier coefficients
  void calcFourierCoefficientsNaive();

  /// calculate fourier coefficients
  /// \param side TPC side
  void calcFourierCoefficientsNaive(const o2::tpc::Side side);

  /// calculate fourier coefficients
  void calcFourierCoefficientsFFTW3();

  /// calculate fourier coefficients using FFTW3 package
  /// get IDC0 values from the inverse fourier transform. Can be used for debugging. std::vector<std::vector<float>>: first vector interval second vector IDC0 values
  /// \param side TPC side
  void calcFourierCoefficientsFFTW3(const o2::tpc::Side side);

  /// get IDC0 values from the inverse fourier transform. Can be used for debugging. std::vector<std::vector<float>>: first vector interval second vector IDC0 values
  /// \param side TPC side
  std::vector<std::vector<float>> inverseFourierTransformNaive(const o2::tpc::Side side) const;

  /// get IDC0 values from the inverse fourier transform using FFTW3. Can be used for debugging. std::vector<std::vector<float>>: first vector interval second vector IDC0 values
  /// \param side TPC side
  std::vector<std::vector<float>> inverseFourierTransformFFTW3(const o2::tpc::Side side) const;

  /// divide coefficients by number of IDCs used
  void normalizeCoefficients(const o2::tpc::Side side)
  {
    std::transform(mFourierCoefficients.mFourierCoefficients[side].begin(), mFourierCoefficients.mFourierCoefficients[side].end(), mFourierCoefficients.mFourierCoefficients[side].begin(), [norm = this->mRangeIDC](auto& val) { return val / norm; });
  };

  /// \return returns maximum numbers of stored real/imag fourier coeffiecients
  unsigned int getNMaxCoefficients() const { return this->mRangeIDC / 2 + 1; }

  ClassDefNV(IDCFourierTransform, 1)
};

} // namespace o2::tpc

#endif
