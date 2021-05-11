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
#include <cmath>
#include <complex>
#include "Rtypes.h"
#include "TFile.h"
#include "DataFormatsTPC/Defs.h"
#include "Framework/Logger.h"
#include "CommonConstants/MathConstants.h"

#include "CommonUtils/TreeStreamRedirector.h"

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

  unsigned int getRangeIntegrationIntervals() const { return mRangeIntegrationIntervals; }
  unsigned int getShift() const { return mShift; }
  unsigned int getNCoefficients() const { return mFourierCoefficients[o2::tpc::Side::A].size(); }             /// \return returns numbers of stored fourier coeffiecients
  unsigned int getNIDCs(const o2::tpc::Side side) const { return mIDCsOne[side].size(); }                     /// \return returns number of 1D-IDCs
  unsigned int getNIntervals(const o2::tpc::Side side) const { return mFourierCoefficients[side][0].size(); } /// \return returns number of intervals for which the coefficients are obtained
  std::array<std::vector<float>, o2::tpc::SIDES> getIDCsOne() const { return mIDCsOne; }

  void setIDCs(std::vector<float>&& idcs, const o2::tpc::Side side)
  {
    mIDCsOne[side] = std::move(idcs);
    for (auto& coeff : mFourierCoefficients[side]) {
      coeff.resize(std::ceil((getNIDCs(side) - mRangeIntegrationIntervals) / static_cast<float>(mShift)) + 1);
    }
  }

  void setIDCs(const std::vector<float>& idcs, const o2::tpc::Side side)
  {
    mIDCsOne[side] = idcs;
    for (auto& coeff : mFourierCoefficients[side]) {
      coeff.resize(std::ceil((getNIDCs(side) - mRangeIntegrationIntervals) / static_cast<float>(mShift)) + 1);
    }
  }

  /// calculate fourier coefficient
  void calFourierCoefficients(const o2::tpc::Side side)
  {
    // see: https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Definitiona
    using namespace std::complex_literals;

    // loop over all the intervals. For each interval the coefficients are calculated
    for (unsigned long interval = 0; interval < getNIntervals(side); ++interval) {
      // loop over coefficients which will be calculated
      for (unsigned int coeff = 0; coeff < getNCoefficients(); ++coeff) {
        const unsigned int lastIndex = (interval == getNIntervals(side) - 1) ? getNIDCs(side) - mShift * interval : mRangeIntegrationIntervals;
        for (unsigned int index = 0; index < lastIndex; ++index) {
          const float term = o2::constants::math::TwoPI * coeff * index / lastIndex;
          const float idc = mIDCsOne[side][index + mShift * interval];
          mFourierCoefficients[side][coeff][interval] += idc * std::complex<float>(std::cos(term), std::sin(-term));
        }
        // normalize coefficient to number of used points
        mFourierCoefficients[side][coeff][interval] /= lastIndex;
      }
    }
  }

  auto getFourierCoefficients(const o2::tpc::Side side) const { return mFourierCoefficients[side]; }

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName = "Fourier.root", const char* outName = "FourierCoefficients") const
  {
    TFile fOut(outFileName, "RECREATE");
    fOut.WriteObject(this, outName);
    fOut.Close();
  }

  void dumpToTree() const
  {
    o2::utils::TreeStreamRedirector pcstream("FourierTree.root", "RECREATE");
    pcstream.GetFile()->cd();

    for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
      const o2::tpc::Side side = iSide == 0 ? Side::A : Side::C;
      for (unsigned int coeff = 0; coeff < mFourierCoefficients[iSide].size(); ++coeff) {
        for (unsigned int interval = 0; interval < mFourierCoefficients[iSide][coeff].size(); ++interval) {
          float real = mFourierCoefficients[iSide][coeff][interval].real();
          float imag = mFourierCoefficients[iSide][coeff][interval].imag();

          const unsigned int lastIndex = (interval == getNIntervals(side) - 1) ? getNIDCs(side) - mShift * interval : mRangeIntegrationIntervals;
          std::vector<float> idcOne;
          idcOne.reserve(lastIndex);
          for (unsigned int index = 0; index < lastIndex; ++index) {
            idcOne.emplace_back(mIDCsOne[side][index + mShift * interval]);
          }

          pcstream << "tree"
                   << "side=" << iSide
                   << "interval=" << interval
                   << "coefficient=" << coeff
                   << "real=" << real
                   << "imag=" << imag
                   << "IDCOne.=" << idcOne
                   << "\n";
        }
      }
    }
  }

 private:
  const unsigned int mRangeIntegrationIntervals{};           ///< number of IDCs used for the calculation of fourier coefficients
  const unsigned int mShift{};                               ///< shifting parameter for the range. Should be >=1
  std::array<std::vector<float>, o2::tpc::SIDES> mIDCsOne{}; ///< all 1D-IDCs which are used to calculate the fourier coefficients
  const unsigned int mNFourierCoefficients{};                ///< number of fourier coefficients which will be calculated
  std::array<std::vector<std::vector<std::complex<float>>>, o2::tpc::SIDES> mFourierCoefficients{
    std::vector<std::vector<std::complex<float>>>(mNFourierCoefficients),
    std::vector<std::vector<std::complex<float>>>(mNFourierCoefficients)}; ///< fourier coefficients. side -> coefficient -> interval

  // }; ///< calculated fourier coefficients

  ClassDefNV(IDCFourierTransform, 1)
};

} // namespace tpc
} // namespace o2

#endif
