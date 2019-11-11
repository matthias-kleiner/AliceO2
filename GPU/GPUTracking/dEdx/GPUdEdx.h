// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUdEdx.h
/// \author David Rohr

#ifndef GPUDEDX_H
#define GPUDEDX_H

#include "GPUDef.h"
#include "GPUTPCGeometry.h"
#include "GPUCommonMath.h"
#include "GPUParam.h"
#include "GPUdEdxInfo.h"

#include "Fairlogger.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
#ifndef HAVE_O2HEADERS

class GPUdEdx
{
 public:
  GPUd() void clear() {}
  GPUd() void fillCluster(float qtot, float qmax, int padRow, float trackSnp, float trackTgl, const GPUParam& param, float pad, float zz, float rms0) {}
  GPUd() void fillSubThreshold(int padRow, const GPUParam& param) {}
  GPUd() void computedEdx(GPUdEdxInfo& output, const GPUParam& param) {}
};

#else

class GPUdEdx
{
 public:
  // The driver must call clear(), fill clusters row by row outside-in, then run computedEdx() to get the result
  GPUd() void clear();
  GPUd() void fillCluster(float qtot, float qmax, int padRow, float trackSnp, float trackTgl, const GPUParam& param, float pad, float zz, float rms0);
  GPUd() void fillSubThreshold(int padRow, const GPUParam& param);
  GPUd() void computedEdx(GPUdEdxInfo& output, const GPUParam& param);

 private:
  GPUd() float GetSortTruncMean(float* array, int count, int trunclow, int trunchigh);
  GPUd() void checkSubThresh(int roc);

  //========================== qmax calib =============================
  /// ky: y angle - tan(y) - dy/dx (cm/cm)
  /// kz: Float_t fTAngleZ;     ///< z angle - tan(z) - dz/dx (cm/cm)
  GPUd() std::array<float, 7> qmaxCorrection(const GPUParam& param, int padRow, float cpad, float ctime, float ky, float kz, float rmsy0, float rmsz0, float effPad, float effDiff, int type);
  GPUd() std::array<float, 7> qmaxCorrectionOneDim(const GPUParam& param, int padRow, float cpad, float ctime, float ky, float kz, float rmsy0, float rmsz0, float effPad, float effDiff, int type);

  GPUd() static float GaussConvolution(float x0, float x1, float k0, float k1, float s0, float s1);
  GPUd() static float ErfFast(float x) { return 1 - ErfcFast(x); } // Error function erf(x)
  GPUd() static float ErfcFast(float x);
  GPUd() float GaussConvolutionGamma4(float x0, float x1, float k0, float k1, float s0, float s1, float tau);

  static constexpr int MAX_NCL = GPUCA_ROW_COUNT; // Must fit in mNClsROC (unsigned char)!

  float mChargeTot[MAX_NCL]; // No need for default, just some memory
  float mChargeMax[MAX_NCL]; // No need for default, just some memory
  float mSubThreshMinTot = 0.f;
  float mSubThreshMinMax = 0.f;
  unsigned char mNClsROC[4] = {0};
  unsigned char mNClsROCSubThresh[4] = {0};
  unsigned char mCount = 0;
  unsigned char mLastROC = 255;
  char mNSubThresh = 0;

  o2::tpc::ClusterNative mClNative[MAX_NCL]; //for debugging
};

GPUdi() void GPUdEdx::checkSubThresh(int roc)
{
  if (roc != mLastROC) {
    if (mNSubThresh && mCount + mNSubThresh <= MAX_NCL) {
      for (int i = 0; i < mNSubThresh; i++) {
        mChargeTot[mCount] = mSubThreshMinTot;
        mChargeMax[mCount++] = mSubThreshMinMax;
      }
      mNClsROC[mLastROC] += mNSubThresh;
      mNClsROCSubThresh[mLastROC] += mNSubThresh;
    }
    mNSubThresh = 0;
    mSubThreshMinTot = 1e10f;
    mSubThreshMinMax = 1e10f;
  }

  mLastROC = roc;
}

GPUdi() void GPUdEdx::fillCluster(float qtot, float qmax, int padRow, float trackSnp, float trackTgl, const GPUParam& GPUrestrict() param, float pad, float zz, float rms0)
{
  if (mCount >= MAX_NCL) {
    return;
  }
  int roc = param.tpcGeometry.GetROC(padRow);
  checkSubThresh(roc);
  float snp2 = trackSnp * trackSnp;
  if (snp2 > GPUCA_MAX_SIN_PHI_LOW) {
    snp2 = GPUCA_MAX_SIN_PHI_LOW;
  }
  float factor = CAMath::Sqrt((1 - snp2) / (1 + trackTgl * trackTgl));
  factor /= param.tpcGeometry.PadHeight(padRow);
  qtot *= factor;
  qmax *= factor; // / param.tpcGeometry.PadWidth(padRow);

  //===================== begin qmax as a function of pad correction =======================
  // make qmax as a function of pad position correction
  if (roc > 2) {
    roc = 2; // TODO: Add type 3
  }

  // const float rmsy0 = CAMath::Sqrt(param.GetClusterRMS(0, roc, 0, 0)); // roc is 0=0.205,1=0.285,2=0.295
  const float rmsy0 = 0.f;
  // const float rmsy0 = 0.;                                              // GEM PRF is close to 0
  const float rmsz0 = CAMath::Sqrt(param.GetClusterRMS(1, roc, 0, 0)); // roc is 0:0.244, 1:0.2475, 2:0.2565
  // const float rmsz0 = 0.175f;
  // const float rmsz0 = 0.175f;
  // const float rmsz0 = 0.f;

  // const float effPad = 0.5; // FIX ME get correct value
  const float effPad = 1.; // FIX ME get correct value
  // const float effDiff = 1;  // FIX ME get correct value

  // if is IROC use lower diffusion
  const float effDiff = 1; //padRow<63 ? 0.7 : 1;

  // Float_t zres0 = parcl->GetRMS0(1,ipad,0,0)/param->GetZWidth();
  // Float_t yres0 = parcl->GetRMS0(0,ipad,0,0)/padWidth;

  // ty: y angle - tan(y) - dy/dx (cm/cm)
  const float ty = CAMath::Abs(CAMath::Tan(o2::gpu::GPUCommonMath::ASin(trackSnp)));
  // ty=dy/dx
  // sin(Phi) = dy / sqrt(dx*dx+dy*dy)
  // cos(Phi) = dx / sqrt(dx*dx+dy*dy)
  // tan(Phi) = dy / dx

  // tz: z angle - tan(z) - dz/dx (cm/cm)
  const float tz = CAMath::Abs(trackTgl * CAMath::Sqrt(1 + ty * ty));
  // tz = dz / dx
  // tz = tan(lambda) * sqrt(1 + ty * ty)
  // tz = tan(lambda) * sqrt(1 + dy*dy / dx*dx)
  // tz = dz / sqrt(dx*dx+dy*dy) * sqrt(1 + dy*dy / dx*dx)
  // tz = dz * sqrt( (1+dy*dy/dx*dx) / (dx*dx+dy*dy)  )
  // tz = dz * sqrt( (dx*dx+dy*dy) / (dx*dx+dy*dy) * 1/dx*dx )
  // tz = dz * sqrt( 1/dx*dx )
  // tz = dz / dx
  // Float_t tz = TMath::Abs(point->GetAngleZ() * TMath::Sqrt(1 + ty * ty));
  const float timeBin = (250.f - GPUCommonMath::Abs(zz)) * 1.93798f; // (1/0.516f) = 1.93798f;
  // std::array<float, 7> qmaxcorrGaus = qmaxCorrection(param, padRow, pad, timeBin, ty, tz, rmsy0, rmsz0, effPad, effDiff, 0); // final correction ALL EFFECTS

  std::array<float, 7> qmaxcorrGaus = qmaxCorrectionOneDim(param, padRow, pad, timeBin, ty, tz, rmsy0, rmsz0, effPad, effDiff, 0); // ONLY ONE DIMENSIONAL


  // const float qmaxcorrGamma = qmaxCorrection(param, padRow, pad, timeBin, ty, tz, rmsy0, rmsz0, effPad, effDiff, 1);
  const float qmaxcorrGamma = 0;
  // const float qmaxcorr = 1;

  //==============alternate diffusion corr=====================
  const float padWidth = param.tpcGeometry.PadWidth(padRow); // width of the pad
  // static const auto& gasPar = o2::tpc::ParameterGas::Instance();
  const float diffT1 = 0.0209f; //gasPar.DiffT; //;
  const float diffL1 = 0.0221f; //gasPar.DiffL; //;

  // o2::tpc::GPUCATracking g;
  const float zwidth = 0.516f; // g.getPseudoVDrift(); //0.516f;

  // static const auto& parDet = o2::tpc::ParameterDetector::Instance();
  // const float driftlength = parDet.TPClength - GPUCommonMath::Abs(zz);
  const float driftlength = 250.f - GPUCommonMath::Abs(zz);
  const float fac1T = std::sqrt(driftlength) * diffT1;
  const float fac2L = std::sqrt(driftlength) * diffL1;
  const float diff1T = std::sqrt(fac1T * fac1T) / zwidth;
  const float diff2L = std::sqrt(fac2L * fac2L) / padWidth;
  // const float corrValDiff = 1 / (2 * M_PI * diff1T * diff2L);
  const float corrValDiff = 1 / ( sqrt(2 * M_PI) * diff2L);

  // TF1 f("","1/sqrt(250-x)*0.0209",0,250);

  mClNative[mCount].setPad(pad);
  mClNative[mCount].padrow = padRow;
  mClNative[mCount].qTot = qtot;
  mClNative[mCount].qMax = qmax;
  mClNative[mCount].z = zz;
  // mClNative[mCount] . z = rmsz0;
  mClNative[mCount].ty = ty;
  mClNative[mCount].tz = tz;
  mClNative[mCount].corrVal1 = qmaxcorrGaus.at(6);
  mClNative[mCount].corrVal2 = qmaxcorrGamma;
  mClNative[mCount].corrVal3 = corrValDiff;
  mClNative[mCount].timeVal = timeBin;

  mClNative[mCount].py = qmaxcorrGaus.at(0);
  mClNative[mCount].pz = qmaxcorrGaus.at(1);
  mClNative[mCount].pky = qmaxcorrGaus.at(2);
  mClNative[mCount].pkz = qmaxcorrGaus.at(3);
  mClNative[mCount].sy = qmaxcorrGaus.at(4);
  mClNative[mCount].sz = qmaxcorrGaus.at(5);

  mClNative[mCount].setSigmaPad(rms0);
  //===========================================================================================================

  mChargeTot[mCount] = qmax;
  mChargeMax[mCount++] = qmax / qmaxcorrGaus.at(6);
  mNClsROC[roc]++;
  if (qtot < mSubThreshMinTot) {
    mSubThreshMinTot = qtot;
  }
  if (qmax < mSubThreshMinMax) {
    mSubThreshMinMax = qmax;
  }
}

GPUdi() void GPUdEdx::fillSubThreshold(int padRow, const GPUParam& GPUrestrict() param)
{
  const int roc = param.tpcGeometry.GetROC(padRow);
  checkSubThresh(roc);
  mNSubThresh++;
}

#endif // !HAVE_O2HEADERS
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
