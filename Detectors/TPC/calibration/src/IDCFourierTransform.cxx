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
#include <fftw3.h>


#if (defined(WITH_OPENMP) || defined(_OPENMP)) && !defined(__CLING__)
#include <omp.h>
#endif

void o2::tpc::IDCFourierTransform::setIDCs(IDCOne&& idcsOne, std::vector<unsigned int>&& integrationIntervalsPerTF)
{
  mIDCOne[mBufferIndex] = std::move(idcsOne);
  mIntegrationIntervalsPerTF[mBufferIndex] = std::move(integrationIntervalsPerTF);
  mBufferIndex = !mBufferIndex; // switch buffer index
}

void o2::tpc::IDCFourierTransform::setIDCs(const IDCOne& idcsOne, const std::vector<unsigned int>& integrationIntervalsPerTF)
{
  mIDCOne[mBufferIndex] = idcsOne;
  mIntegrationIntervalsPerTF[mBufferIndex] = integrationIntervalsPerTF;
  mBufferIndex = !mBufferIndex; // switch buffer index
}

void o2::tpc::IDCFourierTransform::calcFourierCoefficients()
{
  if (isEmpty()) {
    LOGP(info, "returning! One of the IDC aggregation intervals is empty");
    return;
  }
  if (sFftw) {
    calcFourierCoefficientsFFTW3();
  } else {
    calcFourierCoefficientsNaive();
  }
}

void o2::tpc::IDCFourierTransform::calcFourierCoefficientsNaive()
{
  const std::vector<int> endIndex = getLastIntervals();
  calcFourierCoefficientsNaive(o2::tpc::Side::A, endIndex);
  calcFourierCoefficientsNaive(o2::tpc::Side::C, endIndex);
}

void o2::tpc::IDCFourierTransform::calcFourierCoefficientsFFTW3()
{
  const std::vector<int> endIndex = getLastIntervals();
  calcFourierCoefficientsFFTW3(o2::tpc::Side::A, endIndex);
  calcFourierCoefficientsFFTW3(o2::tpc::Side::C, endIndex);
}

void o2::tpc::IDCFourierTransform::calcFourierCoefficientsNaive(const o2::tpc::Side side, const std::vector<int>& endIndex)
{
  // see: https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Definitiona
#pragma omp parallel for num_threads(sNThreads)
  for (unsigned int interval = 0; interval < getNIntervals(); ++interval) {
    for (unsigned int coeff = 0; coeff < mNFourierCoefficients; ++coeff) {
      const unsigned int indexData = getIndex(interval, coeff); // index for storing fourier coefficient
      const float term0 = o2::constants::math::TwoPI * coeff / mRangeIDC;

      unsigned int counter = 0;
      for (int index = getStartIndex(endIndex[interval]); index <= endIndex[interval]; ++index) {
        const float term = term0 * counter;
        const float idc = getIDCOne(index, side);
        mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::REAL) += idc * std::cos(term);
        mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::IMAG) -= idc * std::sin(term);
        ++counter;
      }
      // normalize coefficient to number of used points
      mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::REAL) /= mRangeIDC;
      mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::IMAG) /= mRangeIDC;
    }
  }
}

void o2::tpc::IDCFourierTransform::calcFourierCoefficientsFFTW3(const o2::tpc::Side side, const std::vector<int>& endIndex)
{
  const auto sizeInp = sizeof(float) * mRangeIDC;
  const auto sizeOut = sizeof(fftwf_complex) * mNFourierCoefficients;

  // for FFTW and OMP see: https://stackoverflow.com/questions/15012054/fftw-plan-creation-using-openmp
#pragma omp parallel num_threads(sNThreads)
  {
    float* val1DIDCs = (float*)fftwf_malloc(sizeInp);
    fftwf_complex* coefficients = (fftwf_complex*)fftwf_malloc(sizeOut);
    fftwf_plan fftwPlan;

#pragma omp critical(make_plan)
    fftwPlan = fftwf_plan_dft_r2c_1d(mRangeIDC, val1DIDCs, coefficients, FFTW_ESTIMATE); // fftw_plan_dft_r2c_1d: real input and complex output

#pragma omp for
    for (unsigned int interval = 0; interval < getNIntervals(); ++interval) {
      unsigned int counter = 0;
      for (int index = getStartIndex(endIndex[interval]); index <= endIndex[interval]; ++index) {
        val1DIDCs[counter++] = getIDCOne(index, side);
      }

      fftwf_execute_dft_r2c(fftwPlan, val1DIDCs, coefficients);
      for (unsigned int coeff = 0; coeff < mNFourierCoefficients; ++coeff) {
        const unsigned int indexData = getIndex(interval, coeff);
        // normalize coefficient to number of used points
        mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::REAL) = coefficients[coeff][0] / mRangeIDC;
        mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::IMAG) = coefficients[coeff][1] / mRangeIDC;
      }
    }
    fftwf_destroy_plan(fftwPlan);
    fftwf_free(coefficients);
    fftwf_free(val1DIDCs);
  }
}

std::vector<std::vector<float>> o2::tpc::IDCFourierTransform::inverseFourierTransform(const o2::tpc::Side side) const
{
  // vector containing for each intervall the inverse fourier IDCs
  std::vector<std::vector<float>> inverse(getNIntervals());

  // loop over all the intervals. For each interval the coefficients are calculated
  for (unsigned int interval = 0; interval < getNIntervals(); ++interval) {
    inverse[interval].resize(mRangeIDC);
    for (unsigned int index = 0; index < mRangeIDC; ++index) {
      const float term0 = o2::constants::math::TwoPI * index / mRangeIDC;
      unsigned int coeffTmp = 0;
      int fac = 1; // if input data is real (and it is) the coefficient is mirrored https://dsp.stackexchange.com/questions/4825/why-is-the-fft-mirrored
      for (unsigned int coeff = 0; coeff < mRangeIDC; ++coeff) {
        const unsigned int indexData = getIndex(interval, coeffTmp);
        const float term = term0 * coeff;
        inverse[interval][index] += mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::REAL) * std::cos(term) - fac * mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::IMAG) * std::sin(term);
        if (coeff < mNFourierCoefficients - 1) {
          ++coeffTmp;
        } else {
          --coeffTmp;
          fac = -1;
        };
      }
    }
  }
  return inverse;
}

std::vector<std::vector<float>> o2::tpc::IDCFourierTransform::inverseFourierTransformFFTW3(const o2::tpc::Side side) const
{
  // vector containing for each intervall the inverse fourier IDCs
  std::vector<std::vector<float>> inverse(getNIntervals());

  // loop over all the intervals. For each interval the coefficients are calculated
  for (unsigned int interval = 0; interval < getNIntervals(); ++interval) {
    inverse[interval].resize(mRangeIDC);
    std::vector<std::array<float, 2>> val1DIDCs;
    val1DIDCs.reserve(mRangeIDC);
    for (unsigned int index = 0; index < getNCoefficients(); ++index) {
      const unsigned int indexData = getIndex(interval, index);
      val1DIDCs.emplace_back(std::array<float, 2>{mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::REAL), mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::IMAG)});
    }
    const fftwf_plan fftwPlan = fftwf_plan_dft_c2r_1d(mRangeIDC, reinterpret_cast<fftwf_complex*>(val1DIDCs.data()), inverse[interval].data(), FFTW_ESTIMATE); // fftw_plan_dft_r2c_1d: real input and complex output
    fftwf_execute(fftwPlan);
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
  const std::vector<int> endIndex = getLastIntervals();
  for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
    const o2::tpc::Side side = iSide == 0 ? Side::A : Side::C;
    const auto inverseFourier = inverseFourierTransform(side);
    const auto inverseFourierFFTW3 = inverseFourierTransformFFTW3(side);

    for (unsigned int coeff = 0; coeff < getNCoefficients(); ++coeff) {
      for (unsigned int interval = 0; interval < getNIntervals(); ++interval) {
        const unsigned int indexData = getIndex(interval, coeff);
        float real = mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::REAL);
        float imag = mFourierCoefficients(side, indexData, FourierCoeff::CoeffType::IMAG);
        std::vector<float> idcOneInverse = inverseFourier[interval];
        std::vector<float> idcOneInverseFFTW3 = inverseFourierFFTW3[interval];

        // get 1D-IDC values used for calculation of the fourier coefficients
        std::vector<float> idcOne;
        idcOne.reserve(mRangeIDC);
        for (int index = getStartIndex(endIndex[interval]); index <= endIndex[interval]; ++index) {
          idcOne.emplace_back(getIDCOne(index, side));
        }

        pcstream << "tree"
                 << "side=" << iSide
                 << "interval=" << interval
                 << "coefficient=" << coeff
                 << "real=" << real
                 << "imag=" << imag
                 << "IDCOne.=" << idcOne
                 << "IDCOneiDFT.=" << idcOneInverse
                 << "IDCOneiDFTFFTW3.=" << idcOneInverseFFTW3
                 << "\n";
      }
    }
  }
}

std::vector<int> o2::tpc::IDCFourierTransform::getLastIntervals() const
{
  std::vector<int> endIndex;
  endIndex.reserve(mTimeFrames);
  endIndex.emplace_back(mIntegrationIntervalsPerTF[!mBufferIndex][0] - 1);
  for (unsigned int interval = 1; interval < mTimeFrames; ++interval) {
    endIndex.emplace_back(endIndex[interval - 1] + mIntegrationIntervalsPerTF[!mBufferIndex][interval]);
  }
  return endIndex;
}
