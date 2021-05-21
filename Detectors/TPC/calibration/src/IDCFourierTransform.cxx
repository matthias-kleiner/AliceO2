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

void o2::tpc::IDCFourierTransform::calcFourierCoefficients()
{
  calcFourierCoefficients(o2::tpc::Side::A);
  calcFourierCoefficients(o2::tpc::Side::C);
  mBufferIndex = !mBufferIndex; // switch buffer index
}

void o2::tpc::IDCFourierTransform::calcFourierCoefficients(const o2::tpc::Side side)
{
  // see: https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Definitiona
  // loop over all the intervals. For each interval the coefficients are calculated. one interval corresponds to one TF
  int currentIndexEndInterval = -1;
  for (unsigned int interval = 0; interval < getNIntervals(); ++interval) {
    currentIndexEndInterval += mIntegrationIntervalsPerTF[mBufferIndex][interval];

    // set smaller range of IDCS if the buffer is empty
    const int rangeIDC = getRangeIDC(currentIndexEndInterval, side, mBufferIndex);

    // loop over coefficients which will be calculated
    for (unsigned int coeff = 0; coeff < getNCoefficients(); ++coeff) {
      const unsigned int indexData = getIndex(interval, coeff); // index for storing fourier coefficient
      const float term0 = o2::constants::math::TwoPI * coeff / rangeIDC;

      // start from last index and go backwards
      int indexNumber = rangeIDC - 1;
      for (int index = currentIndexEndInterval; index > (currentIndexEndInterval - rangeIDC); --index) {
        const float term = term0 * indexNumber;
        --indexNumber;

        const float idc = getIDCOne(index, side, mBufferIndex) - 1;
        // LOGP(info, "interval: {}  coeff: {}   index: {}     indexNumber: {}    rangeIDC: {}   idc: {}", interval, coeff, index, indexNumber + 1, rangeIDC, idc);
        mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::REAL) += idc * std::cos(term);
        mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::IMAG) -= idc * std::sin(term);
      }
      // normalize coefficient to number of used points
      mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::REAL) /= rangeIDC;
      mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::IMAG) /= rangeIDC;
    }
  }
}

std::vector<std::vector<float>> o2::tpc::IDCFourierTransform::inverseFourierTransform(const o2::tpc::Side side) const
{
  // vector containing for each intervall the inverse fourier IDCs
  std::vector<std::vector<float>> inverse(getNIntervals());

  // loop over all the intervals. For each interval the coefficients are calculated
  int currentIndexEndInterval = -1;
  for (unsigned int interval = 0; interval < getNIntervals(); ++interval) {
    currentIndexEndInterval += mIntegrationIntervalsPerTF[!mBufferIndex][interval];
    const int rangeIDC = getRangeIDC(currentIndexEndInterval, side, !mBufferIndex);

    // loop over coefficients which will be calculated
    inverse[interval].resize(rangeIDC);
    for (unsigned int index = 0; index < rangeIDC; ++index) {
      const float term0 = o2::constants::math::TwoPI * index / rangeIDC;
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

void o2::tpc::IDCFourierTransform::dumpToTree(const char* outFileName) const
{
  o2::utils::TreeStreamRedirector pcstream(outFileName, "RECREATE");
  pcstream.GetFile()->cd();

  for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
    const o2::tpc::Side side = iSide == 0 ? Side::A : Side::C;
    const auto inverseFourier = inverseFourierTransform(side);
    for (unsigned int coeff = 0; coeff < getNCoefficients(); ++coeff) {
      int currentIndexEndInterval = -1;
      for (unsigned int interval = 0; interval < getNIntervals(); ++interval) {
        const unsigned int indexData = getIndex(interval, coeff);
        float real = mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::REAL);
        float imag = mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::IMAG);
        std::vector<float> idcOneInverse = inverseFourier[interval];

        // const unsigned int lastIndex = getLastIndex(interval, side);
        currentIndexEndInterval += mIntegrationIntervalsPerTF[!mBufferIndex][interval];
        const int rangeIDC = getRangeIDC(currentIndexEndInterval, side, !mBufferIndex);

        // get 1D-IDC values used for calculation of the fourier coefficients
        std::vector<float> idcOne;
        idcOne.reserve(rangeIDC);
        for (int index = (currentIndexEndInterval - rangeIDC) + 1; index <= currentIndexEndInterval; ++index) {
          idcOne.emplace_back(getIDCOne(index, side, !mBufferIndex));
        }

        pcstream << "tree"
                 << "side=" << iSide
                 << "interval=" << interval
                 << "coefficient=" << coeff
                 << "real=" << real
                 << "imag=" << imag
                 << "IDCOne.=" << idcOne
                 << "IDCOneInverseFF.=" << idcOneInverse
                 << "\n";
      }
    }
  }
}
