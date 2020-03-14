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
#include "DataFormatsTPC/ClusterNative.h"
#include <array>

// #include "Fairlogger.h"
// #include "/Users/matthias/alice/O2/Detectors/TPC/base/include/TPCBase/PadPos.h"
// #include "/Users/matthias/alice/O2/Detectors/TPC/base/include/TPCBase/Mapper.h"
// #include "/Users/matthias/alice/O2/Detectors/TPC/base/include/TPCBase/PadPos.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
#ifndef HAVE_O2HEADERS

class GPUdEdx
{
 public:
  GPUd() void clear() {}
  GPUd() void fillCluster(float qtot, float qmax, int padRow, float trackSnp, float trackTgl, const GPUParam& param, float pad, float zz, float rms0, float time, const tpc::ClusterNative& clNat) {}
  GPUd() void fillSubThreshold(int padRow, const GPUParam& param) {}
  GPUd() void computedEdx(GPUdEdxInfo& output, const GPUParam& param) {}
};

#else

class GPUdEdx
{
 public:
  // The driver must call clear(), fill clusters row by row outside-in, then run computedEdx() to get the result
  GPUd() void clear();
  GPUd() void fillCluster(float qtot, float qmax, int padRow, float trackSnp, float trackTgl, const GPUParam& param, float pad, float zz, float rms0, float time, const tpc::ClusterNative& clNat);
  GPUd() void fillSubThreshold(int padRow, const GPUParam& param);
  GPUd() void computedEdx(GPUdEdxInfo& output, const GPUParam& param);

 private:
  GPUd() float GetSortTruncMean(float* array, int count, int trunclow, int trunchigh);
  GPUd() void checkSubThresh(int roc);
  GPUd() std::array<float, 7> qmaxCorrection(const GPUParam& param, int padRow, float cpad, float driftDistance, float ky, float kz, float rmsy0, float rmsz0, float effPad, float effDiff, float time, int correctionType);
  GPUd() std::array<float, 7> qmaxCorrectionOneDim(const GPUParam& param, int padRow, float cpad, float driftDistance, float ky, float kz, float rmsy0, float rmsz0, float effPad, float effDiff, float altTime);

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

  o2::tpc::ClusterNative mClNative[MAX_NCL]{}; //for debugging
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

GPUdi() void GPUdEdx::fillCluster(float qtot, float qmax, int padRow, float trackSnp, float trackTgl, const GPUParam& GPUrestrict() param, float pad, float zz, float rms0, float time, const tpc::ClusterNative& clNat)
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
  float tgl2 = trackTgl * trackTgl;
  float factor = CAMath::Sqrt((1 - snp2) / (1 + tgl2));
  factor /= param.tpcGeometry.PadHeight(padRow);
  // qtot *= factor;
  // qmax *= factor; // / param.tpcGeometry.PadWidth(padRow);

  //===================== begin qmax as a function of pad correction =======================
  // make qmax as a function of pad position correction
  if (roc > 2) {
    roc = 2; // TODO: Add type 3
  }

  const float sec2 = 1.f / (1.f - snp2);
  const float ty =  CAMath::Sqrt( snp2 * sec2 ); // fast
  // ty=dy/dx
  // sin(Phi) = dy / sqrt(dx*dx+dy*dy)
  // cos(Phi) = dx / sqrt(dx*dx+dy*dy)
  // tan(Phi) = dy / dx

  // tz: z angle - tan(z) - dz/dx (cm/cm)
  const float tz = CAMath::Sqrt(tgl2 * sec2); // fast
  // tz = dz / dx
  // tz = tan(lambda) * sqrt(1 + ty * ty)
  // tz = tan(lambda) * sqrt(1 + dy*dy / dx*dx)
  // tz = dz / sqrt(dx*dx+dy*dy) * sqrt(1 + dy*dy / dx*dx)
  // tz = dz * sqrt( (1+dy*dy/dx*dx) / (dx*dx+dy*dy)  )
  // tz = dz * sqrt( (dx*dx+dy*dy) / (dx*dx+dy*dy) * 1/dx*dx )
  // tz = dz * sqrt( 1/dx*dx )
  // tz = dz / dx

  const float rmsy0 = 0.f;
  const float rmsz0 = 0.16f / 2.35f; //this doesnt effect the shape of the qmax vs rel pad distr., but a lower value shifts the qmxa distr down
  const float effPad  = 1.; // this should be needed , o big difference
  const float effDiff = 1;

  const float driftDistance = (250.f - GPUCommonMath::Abs(zz));
  std::array<float, 7> qmaxcorrGaus  = qmaxCorrection(param, padRow, clNat.getPad(), driftDistance, ty, tz, rmsy0, rmsz0, effPad, effDiff, clNat.getTime(), 1);
  // std::array<float, 7> qmaxcorrGaus2 = qmaxCorrection(param, padRow, clNat.getPad(), driftDistance, ty, tz, rmsy0, rmsz0, effPad, effDiff, clNat.getTime(), 1);
  // std::array<float, 7> qmaxcorrGaus3 = qmaxCorrection(param, padRow, clNat.getPad(), driftDistance, ty, tz, rmsy0, rmsz0, effPad, effDiff, clNat.getTime(), 2);

  mClNative[mCount].setPad(pad);
  mClNative[mCount].setTime(time);
  mClNative[mCount].clpad = clNat.getPad() - static_cast<int>(clNat.getPad() + 0.5f);
  mClNative[mCount].cltime = clNat.getTime() - static_cast<int>(clNat.getTime() + 0.5f);
  mClNative[mCount].padrow = padRow;
  mClNative[mCount].qTot = qtot;
  mClNative[mCount].qMax = qmax;

  mClNative[mCount].z = zz; // z from prop model
  mClNative[mCount].ty = ty; // angle ty
  mClNative[mCount].tz = tz; // angle tz

  mClNative[mCount].corrVal1 = qmaxcorrGaus[6];
  // mClNative[mCount].corrVal2 = qmaxcorrGaus2[6];
  // mClNative[mCount].corrVal3 = qmaxcorrGaus3[6];
  mClNative[mCount].corrVal2 = 0;
  mClNative[mCount].corrVal3 = factor;
  mClNative[mCount].timeVal = 0;

  mClNative[mCount].py = qmaxcorrGaus[0];
  mClNative[mCount].pz = qmaxcorrGaus[1];
  mClNative[mCount].pky = qmaxcorrGaus[2];
  mClNative[mCount].pkz = qmaxcorrGaus[3];
  mClNative[mCount].sy = qmaxcorrGaus[4];
  mClNative[mCount].sz = qmaxcorrGaus[5];

  mClNative[mCount].setSigmaTime( clNat.getSigmaTime() );
  mClNative[mCount].setSigmaPad( clNat.getSigmaPad() );
  //===========================================================================================================

  mChargeTot[mCount] = qmax;
  mChargeMax[mCount++] = qmax / qmaxcorrGaus[6];
  mNClsROC[roc]++;
  if (qtot < mSubThreshMinTot) {
    mSubThreshMinTot = qtot;
  }
  if (qmax < mSubThreshMinMax) {
    mSubThreshMinMax = qmax;
  }
}



GPUdi() void GPUdEdx::fillSubThreshold(int padRow, const GPUParam& param)
{
  const int roc = param.tpcGeometry.GetROC(padRow);
  checkSubThresh(roc);
  mNSubThresh++;
}

#endif // !HAVE_O2HEADERS
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
