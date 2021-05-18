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
#include "CommonConstants/MathConstants.h"
#include "TFile.h"
#include <cmath>

void o2::tpc::IDCFourierTransform::calcFourierCoefficients(const o2::tpc::Side side)
{
  // see: https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Definitiona
  // loop over all the intervals. For each interval the coefficients are calculated
  for (unsigned int interval = 0; interval < getNIntervals(side); ++interval) {
    // loop over coefficients which will be calculated
    for (unsigned int coeff = 0; coeff < getNCoefficients(); ++coeff) {
      const unsigned int indexData = getIndex(interval, coeff);

      const unsigned int lastIndex = getLastIndex(interval, side);
      const float term0 = o2::constants::math::TwoPI * coeff / lastIndex;
      for (unsigned int index = 0; index < lastIndex; ++index) {
        const float term = term0 * index;
        const float idc = mIDCsOne[mBufferIndex][side][index + mShift * interval];

        mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::REAL) += idc * std::cos(term);
        mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::IMAG) -= idc * std::sin(term);
      }
      // normalize coefficient to number of used points
      mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::REAL) /= lastIndex;
      mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::IMAG) /= lastIndex;
    }
  }
  mBufferIndex = !mBufferIndex; // switch buffer index
}

std::vector<std::vector<float>> o2::tpc::IDCFourierTransform::inverseFourierTransform(const o2::tpc::Side side) const
{
  // vector containing for each intervall the inverse fourier IDCs
  std::vector<std::vector<float>> inverse(getNIntervals(side));

  // loop over all the intervals. For each interval the coefficients are calculated
  for (unsigned int interval = 0; interval < getNIntervals(side); ++interval) {
    // loop over coefficients which will be calculated
    const unsigned int lastIndex = getLastIndex(interval, side);
    inverse[interval].resize(lastIndex);
    for (unsigned int index = 0; index < lastIndex; ++index) {
      const float term0 = o2::constants::math::TwoPI * index / lastIndex;
      for (unsigned int coeff = 0; coeff < getNCoefficients(); ++coeff) {
        const unsigned int indexData = getIndex(interval, coeff);
        const float term = term0 * coeff;
        inverse[interval][index] += mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::REAL) * std::cos(term) - mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::IMAG) * std::sin(term);
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
        const unsigned int indexData = getIndex(interval, coeff);
        float real = mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::REAL);
        float imag = mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::IMAG);

        const unsigned int lastIndex = getLastIndex(interval, side);
        std::vector<float> idcOne;
        idcOne.reserve(lastIndex);
        for (unsigned int index = 0; index < lastIndex; ++index) {
          idcOne.emplace_back(mIDCsOne[!mBufferIndex][side][index + mShift * interval]);
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
  const unsigned int factor = std::ceil((getNIDCs(side) - mRangeIntegrationIntervals) / static_cast<float>(mShift)) + 1;
  mFourierCoefficients.mFourierCoefficients[side][static_cast<int>(FourierCoeff::CoeffType::REAL)].resize(mNFourierCoefficients * factor); // n integration intervals
  mFourierCoefficients.mFourierCoefficients[side][static_cast<int>(FourierCoeff::CoeffType::IMAG)].resize(mNFourierCoefficients * factor); // n integration intervals
}

/// set input 1D-IDCs which are used to calculate fourier coefficients
/// \param idcsOne 1D-IDCs
/// \param side TPC side
void o2::tpc::IDCFourierTransform::setIDCs(std::vector<float>&& idcsOne, const o2::tpc::Side side)
{
  mIDCsOne[mBufferIndex][side] = std::move(idcsOne);
  initCoefficients(side);
}

/// set input 1D-IDCs which are used to calculate fourier coefficients
/// \param idcsOne 1D-IDCs
/// \param side TPC side
void o2::tpc::IDCFourierTransform::setIDCs(const std::vector<float>& idcsOne, const o2::tpc::Side side)
{
  mIDCsOne[mBufferIndex][side] = idcsOne;
  initCoefficients(side);
}
