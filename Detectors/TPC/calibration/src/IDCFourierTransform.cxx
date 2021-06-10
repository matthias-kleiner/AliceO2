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

#include "TPCCalibration/IDCFourierTransform.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/Logger.h"
#include "TFile.h"
#include <fftw3.h>

#if (defined(WITH_OPENMP) || defined(_OPENMP)) && !defined(__CLING__)
#include <omp.h>
#endif

template <class Type>
void o2::tpc::IDCFourierTransform<Type>::calcFourierCoefficientsNaive()
{
  if (mFourierCoefficients.getNCoefficientsPerTF() % 2) {
    LOGP(warning, "number of specified fourier coefficients is {}, but should be an even number! you can use FFTW3 method instead!", mFourierCoefficients.getNCoefficientsPerTF());
  }
  LOGP(info, "calculating fourier coefficients for current TF using naive approach using {} threads", sNThreads);
  calcFourierCoefficientsNaive(o2::tpc::Side::A);
  calcFourierCoefficientsNaive(o2::tpc::Side::C);
}

template <class Type>
void o2::tpc::IDCFourierTransform<Type>::calcFourierCoefficientsFFTW3()
{
  LOGP(info, "calculating fourier coefficients for current TF using fftw3 using {} threads", sNThreads);
  calcFourierCoefficientsFFTW3(o2::tpc::Side::A);
  calcFourierCoefficientsFFTW3(o2::tpc::Side::C);
}

template <class Type>
void o2::tpc::IDCFourierTransform<Type>::calcFourierCoefficientsNaive(const o2::tpc::Side side)
{
  // check if IDCs are present for current side
  if (this->getNIDCs(side) == 0) {
    LOGP(warning, "no 1D-IDCs found!");
    mFourierCoefficients.reset(side);
    return;
  }

  const std::vector<unsigned int> offsetIndex = this->getLastIntervals(side);

  // see: https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Definitiona
#pragma omp parallel for num_threads(sNThreads)
  for (unsigned int interval = 0; interval < this->getNIntervals(); ++interval) {
    const std::vector<float>& idcOneExpanded{this->getExpandedIDCOne(side)}; // 1D-IDC values which will be used for the FFT
    for (unsigned int coeff = 0; coeff < mFourierCoefficients.getNCoefficientsPerTF() / 2; ++coeff) {
      const unsigned int indexDataReal = mFourierCoefficients.getIndex(interval, 2 * coeff); // index for storing real fourier coefficient
      const unsigned int indexDataImag = indexDataReal + 1;                                  // index for storing complex fourier coefficient
      const float term0 = o2::constants::math::TwoPI * coeff / this->mRangeIDC;
      for (unsigned int index = 0; index < this->mRangeIDC; ++index) {
        const float term = term0 * index;
        const float idc0 = idcOneExpanded[index + offsetIndex[interval]];
        mFourierCoefficients(side, indexDataReal) += idc0 * std::cos(term);
        mFourierCoefficients(side, indexDataImag) -= idc0 * std::sin(term);
      }
    }
  }
  // normalize coefficient to number of used points
  normalizeCoefficients(side);
}

template <class Type>
void o2::tpc::IDCFourierTransform<Type>::calcFourierCoefficientsFFTW3(const o2::tpc::Side side)
{
  // check if IDCs are present for current side
  if (this->getNIDCs(side) == 0) {
    LOGP(warning, "no 1D-IDCs found!");
    mFourierCoefficients.reset(side);
    return;
  }

  const std::vector<unsigned int> offsetIndex = this->getLastIntervals(side);

  // for FFTW and OMP see: https://stackoverflow.com/questions/15012054/fftw-plan-creation-using-openmp
#pragma omp parallel num_threads(sNThreads)
  {
    const std::vector<float>& idcOneExpanded{this->getExpandedIDCOne(side)}; // 1D-IDC values which will be used for the FFT
    fftwf_plan fftwPlan = nullptr;
    float* val1DIDCs = fftwf_alloc_real(this->mRangeIDC);
    fftwf_complex* coefficients = fftwf_alloc_complex(getNMaxCoefficients());

#pragma omp critical(make_plan)
    fftwPlan = fftwf_plan_dft_r2c_1d(this->mRangeIDC, val1DIDCs, coefficients, FFTW_ESTIMATE);
#pragma omp for
    for (unsigned int interval = 0; interval < this->getNIntervals(); ++interval) {
      std::memcpy(val1DIDCs, &idcOneExpanded[offsetIndex[interval]], this->mRangeIDC * sizeof(float)); // copy IDCs to avoid seg fault when using SIMD instructions
      fftwf_execute_dft_r2c(fftwPlan, val1DIDCs, coefficients);
      std::memcpy(&(*(mFourierCoefficients.mFourierCoefficients[side].begin() + mFourierCoefficients.getIndex(interval, 0))), coefficients, mFourierCoefficients.getNCoefficientsPerTF() * sizeof(float)); // store coefficients
    }
    // free memory
    fftwf_free(coefficients);
    fftwf_free(val1DIDCs);
    fftwf_destroy_plan(fftwPlan);
  }
  normalizeCoefficients(side);
}

template <class Type>
std::vector<std::vector<float>> o2::tpc::IDCFourierTransform<Type>::inverseFourierTransformNaive(const o2::tpc::Side side) const
{
  // vector containing for each intervall the inverse fourier IDCs
  std::vector<std::vector<float>> inverse(this->getNIntervals());
  const float factor = o2::constants::math::TwoPI / this->mRangeIDC;

  // loop over all the intervals. For each interval the coefficients are calculated
  for (unsigned int interval = 0; interval < this->getNIntervals(); ++interval) {
    inverse[interval].resize(this->mRangeIDC);
    for (unsigned int index = 0; index < this->mRangeIDC; ++index) {
      const float term0 = factor * index;
      unsigned int coeffTmp = 0;
      int fac = 1; // if input data is real (and it is) the coefficients are mirrored https://dsp.stackexchange.com/questions/4825/why-is-the-fft-mirrored
      for (unsigned int coeff = 0; coeff < this->mRangeIDC; ++coeff) {
        const unsigned int indexDataReal = mFourierCoefficients.getIndex(interval, 2 * coeffTmp); // index for storing real fourier coefficient
        const unsigned int indexDataImag = indexDataReal + 1;                                     // index for storing complex fourier coefficient
        const float term = term0 * coeff;
        inverse[interval][index] += mFourierCoefficients(side, indexDataReal) * std::cos(term) - fac * mFourierCoefficients(side, indexDataImag) * std::sin(term);
        if (coeff < getNMaxCoefficients() - 1) {
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

template <class Type>
std::vector<std::vector<float>> o2::tpc::IDCFourierTransform<Type>::inverseFourierTransformFFTW3(const o2::tpc::Side side) const
{
  // vector containing for each intervall the inverse fourier IDCs
  std::vector<std::vector<float>> inverse(this->getNIntervals());

  // loop over all the intervals. For each interval the coefficients are calculated
  // this loop and execution of FFTW is not optimized as it is used only for debugging
  for (unsigned int interval = 0; interval < this->getNIntervals(); ++interval) {
    inverse[interval].resize(this->mRangeIDC);
    std::vector<std::array<float, 2>> val1DIDCs;
    val1DIDCs.reserve(this->mRangeIDC);
    for (unsigned int index = 0; index < getNMaxCoefficients(); ++index) {
      const unsigned int indexDataReal = mFourierCoefficients.getIndex(interval, 2 * index); // index for storing real fourier coefficient
      const unsigned int indexDataImag = indexDataReal + 1;                                  // index for storing complex fourier coefficient
      val1DIDCs.emplace_back(std::array<float, 2>{mFourierCoefficients(side, indexDataReal), mFourierCoefficients(side, indexDataImag)});
    }
    const fftwf_plan fftwPlan = fftwf_plan_dft_c2r_1d(this->mRangeIDC, reinterpret_cast<fftwf_complex*>(val1DIDCs.data()), inverse[interval].data(), FFTW_ESTIMATE);
    fftwf_execute(fftwPlan);
    fftwf_destroy_plan(fftwPlan);
  }
  return inverse;
}

template <class Type>
void o2::tpc::IDCFourierTransform<Type>::dumpToFile(const char* outFileName, const char* outName) const
{
  TFile fOut(outFileName, "RECREATE");
  fOut.WriteObject(this, outName);
  fOut.Close();
}

template <class Type>
void o2::tpc::IDCFourierTransform<Type>::dumpToTree(const char* outFileName) const
{
  o2::utils::TreeStreamRedirector pcstream(outFileName, "RECREATE");
  pcstream.GetFile()->cd();
  for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
    const o2::tpc::Side side = iSide == 0 ? Side::A : Side::C;
    const std::vector<unsigned int> offsetIndex = this->getLastIntervals(side);
    const auto idcOneExpanded = this->getExpandedIDCOne(side);
    const auto inverseFourier = inverseFourierTransformNaive(side);
    const auto inverseFourierFFTW3 = inverseFourierTransformFFTW3(side);

    for (unsigned int interval = 0; interval < this->getNIntervals(); ++interval) {
      std::vector<float> oneDIDCInverse = inverseFourier[interval];
      std::vector<float> oneDIDCInverseFFTW3 = inverseFourierFFTW3[interval];

      // get 1D-IDC values used for calculation of the fourier coefficients
      std::vector<float> oneDIDC;
      oneDIDC.reserve(this->mRangeIDC);
      for (unsigned int index = 0; index < this->mRangeIDC; ++index) {
        oneDIDC.emplace_back(idcOneExpanded[index + offsetIndex[interval]]);
      }

      for (unsigned int coeff = 0; coeff < mFourierCoefficients.getNCoefficientsPerTF(); ++coeff) {
        float coefficient = mFourierCoefficients(side, mFourierCoefficients.getIndex(interval, coeff));

        pcstream << "tree"
                 << "side=" << iSide
                 << "interval=" << interval
                 << "icoefficient=" << coeff      // index of ith coefficient
                 << "coefficient=" << coefficient // value for ith coefficient
                 << "1DIDC.=" << oneDIDC
                 << "1DIDCiDFT.=" << oneDIDCInverse
                 << "1DIDCiDFTFFTW3.=" << oneDIDCInverseFFTW3
                 << "\n";
      }
    }
  }
}

template <class Type>
void o2::tpc::IDCFourierTransform<Type>::printFFTWPlan() const
{
  float* val1DIDCs = fftwf_alloc_real(this->mRangeIDC);
  fftwf_complex* coefficients = fftwf_alloc_complex(getNMaxCoefficients());
  fftwf_plan fftwPlan = fftwf_plan_dft_r2c_1d(this->mRangeIDC, val1DIDCs, coefficients, FFTW_ESTIMATE);
  char* splan = fftwf_sprint_plan(fftwPlan);

  LOGP(info, "========= printing FFTW plan ========= \n {}", splan);
  double add = 0;
  double mul = 0;
  double fusedMultAdd = 0;
  fftwf_flops(fftwPlan, &add, &mul, &fusedMultAdd);
  LOGP(info, "additions: {}    multiplications: {}    fused multiply-add: {}    sum: {}", add, mul, fusedMultAdd, add + mul + fusedMultAdd);

  // free memory
  free(splan);
  fftwf_free(coefficients);
  fftwf_free(val1DIDCs);
  fftwf_destroy_plan(fftwPlan);
}

template class o2::tpc::IDCFourierTransform<o2::tpc::IDCFTType::IDCFourierTransformBaseEPN>;
template class o2::tpc::IDCFourierTransform<o2::tpc::IDCFTType::IDCFourierTransformBaseAggregator>;
