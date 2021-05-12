// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TPCCalibration/IDCFourierTransform.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "TFile.h"
#include "CommonConstants/MathConstants.h"
#include <cmath>

void o2::tpc::IDCFourierTransform::calcFourierCoefficients(const o2::tpc::Side side)
{
  // see: https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Definitiona

  // loop over all the intervals. For each interval the coefficients are calculated
  for (unsigned long interval = 0; interval < getNIntervals(side); ++interval) {
    // loop over coefficients which will be calculated
    for (unsigned int coeff = 0; coeff < getNCoefficients(); ++coeff) {
      const unsigned int lastIndex = getLastIndex(interval, side);
      const float term0 = o2::constants::math::TwoPI * coeff / lastIndex;
      for (unsigned int index = 0; index < lastIndex; ++index) {
        const float term = term0 * index;
        const float idc = mIDCsOne[side][index + mShift * interval];
        mFourierCoefficients(side,interval,coeff) += idc * std::complex<float>(std::cos(term), -std::sin(term));
      }
      // normalize coefficient to number of used points
      mFourierCoefficients(side,interval,coeff) /= lastIndex;
    }
  }
}

std::vector<std::vector<float>> o2::tpc::IDCFourierTransform::inverseFourierTransform(const o2::tpc::Side side) const
{
  // vector containing for each intervall the inverse fourier IDCs
  std::vector<std::vector<float>> inverse(getNIntervals(side));

  // loop over all the intervals. For each interval the coefficients are calculated
  for (unsigned long interval = 0; interval < getNIntervals(side); ++interval) {
    // loop over coefficients which will be calculated
    const unsigned int lastIndex = getLastIndex(interval, side);
    inverse[interval].resize(lastIndex);
    for (unsigned int index = 0; index < lastIndex; ++index) {
      const float term0 = o2::constants::math::TwoPI * index / lastIndex;
      for (unsigned int coeff = 0; coeff < getNCoefficients(); ++coeff) {
        const float term = term0 * coeff;
        inverse[interval][index] += mFourierCoefficients(side,interval,coeff).real() * std::cos(term) - mFourierCoefficients(side,interval,coeff).imag() * std::sin(term);
      }
    }
  }
  return inverse;
}

void o2::tpc::IDCFourierTransform::dumpToFile(const char* outFileName, const char* outName) const
{
  TFile fOut(outFileName, "RECREATE");
  fOut.WriteObject(this, outName);
  fOut.Close();
}

void o2::tpc::IDCFourierTransform::dumpToTree() const
{
  o2::utils::TreeStreamRedirector pcstream("FourierTree.root", "RECREATE");
  pcstream.GetFile()->cd();

  for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
    const o2::tpc::Side side = iSide == 0 ? Side::A : Side::C;
    const auto inverseFourier = inverseFourierTransform(side);
    for (unsigned int coeff = 0; coeff < getNCoefficients(); ++coeff) {
      for (unsigned int interval = 0; interval < getNIntervals(o2::tpc::Side::A); ++interval) {
        float real = mFourierCoefficients(side,interval,coeff).real();
        float imag = mFourierCoefficients(side,interval,coeff).imag();

        const unsigned int lastIndex = getLastIndex(interval, side);
        std::vector<float> idcOne;
        idcOne.reserve(lastIndex);
        for (unsigned int index = 0; index < lastIndex; ++index) {
          idcOne.emplace_back(mIDCsOne[side][index + mShift * interval]);
        }

        std::vector<float> idcOneInverse = inverseFourier[interval];
        pcstream << "tree"
                 << "side=" << iSide
                 << "interval=" << interval
                 << "coefficient=" << coeff
                 << "real=" << real
                 << "imag=" << imag
                 << "IDCOne.=" << idcOne
                 << "IDCOneinverseFF.=" << idcOneInverse
                 << "\n";
      }
    }
  }
}

void o2::tpc::IDCFourierTransform::initCoefficients(const o2::tpc::Side side)
{
  mFourierCoefficients.mFourierCoefficients[side].resize(std::ceil((getNIDCs(side) - mRangeIntegrationIntervals) / static_cast<float>(mShift)) + 1); // n integration intervals
  for (auto& interval : mFourierCoefficients.mFourierCoefficients[side]) {
    interval.resize(mNFourierCoefficients);
  }
}
