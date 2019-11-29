// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUdEdx.cxx
/// \author David Rohr

#include "GPUdEdx.h"
#include "GPUTPCGeometry.h"
#include "GPUdEdxInfo.h"
#include "GPUCommonAlgorithm.h"
#include "GPUParam.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"

using namespace GPUCA_NAMESPACE::gpu;

GPUd() void GPUdEdx::clear() { new (this) GPUdEdx; }

GPUd() void GPUdEdx::computedEdx(GPUdEdxInfo& GPUrestrict() output, const GPUParam& GPUrestrict() param)
{
  checkSubThresh(255);
  const int truncLow = param.rec.dEdxTruncLow;
  const int truncHigh = param.rec.dEdxTruncHigh;
  const int countIROC = mNClsROC[0];
  const int countOROC1 = mNClsROC[1];
  const int countOROC2 = mNClsROC[2];
  const int countOROC3 = mNClsROC[3];
  output.dEdxTotIROC = GetSortTruncMean(mChargeTot + countOROC3 + countOROC2 + countOROC1, countIROC, truncLow, truncHigh);
  output.dEdxTotOROC1 = GetSortTruncMean(mChargeTot + countOROC3 + countOROC2, countOROC1, truncLow, truncHigh);
  output.dEdxTotOROC2 = GetSortTruncMean(mChargeTot + countOROC3, countOROC2, truncLow, truncHigh);
  output.dEdxTotOROC3 = GetSortTruncMean(mChargeTot, countOROC3, truncLow, truncHigh);
  output.dEdxTotTPC = GetSortTruncMean(mChargeTot, mCount, truncLow, truncHigh);
  output.dEdxMaxIROC = GetSortTruncMean(mChargeMax + countOROC3 + countOROC2 + countOROC1, countIROC, truncLow, truncHigh);
  output.dEdxMaxOROC1 = GetSortTruncMean(mChargeMax + countOROC3 + countOROC2, countOROC1, truncLow, truncHigh);
  output.dEdxMaxOROC2 = GetSortTruncMean(mChargeMax + countOROC3, countOROC2, truncLow, truncHigh);
  output.dEdxMaxOROC3 = GetSortTruncMean(mChargeMax, countOROC3, truncLow, truncHigh);
  output.dEdxMaxTPC = GetSortTruncMean(mChargeMax, mCount, truncLow, truncHigh);
  output.NHitsIROC = countIROC - mNClsROCSubThresh[0];
  output.NHitsSubThresholdIROC = countIROC;
  output.NHitsOROC1 = countOROC1 - mNClsROCSubThresh[1];
  output.NHitsSubThresholdOROC1 = countOROC1;
  output.NHitsOROC2 = countOROC2 - mNClsROCSubThresh[2];
  output.NHitsSubThresholdOROC2 = countOROC2;
  output.NHitsOROC2 = countOROC3 - mNClsROCSubThresh[3];
  output.NHitsSubThresholdOROC2 = countOROC3;

  for (int i = 0; i < mCount; ++i) {
    output.clNat[i] = mClNative[i];
  }
}

GPUd() float GPUdEdx::GetSortTruncMean(float* GPUrestrict() array, int count, int trunclow, int trunchigh)
{
  trunclow = count * trunclow / 128;
  trunchigh = count * trunchigh / 128;
  if (trunclow >= trunchigh) {
    return (0.);
  }
  CAAlgo::sort(array, array + count);
  float mean = 0;
  for (int i = trunclow; i < trunchigh; i++) {
    mean += array[i];
  }
  return (mean / (trunchigh - trunclow));
}

//===================================== qmax calib =================================
GPUd() std::array<float, 7> GPUdEdx::qmaxCorrectionOneDim(const GPUParam& param, int padRow, float cpad, float driftDistance, float ky, float kz, float rmsy0, float rmsz0, float effPad, float effDiff, int type, float altTime)
{
  /*
    THIS CORRECTION WORKS ONLY FOR TRACKS WITHOUT AN ANGLE THETA OR AN ANGLE Phi
    ONLY RELATIVE PAD POSITION IS CORRECTED
    PLOT; QMAX VS REL PAD SHOULD BE FLAT
  */

  /// This function approximates the charge distribution by assuming a convolution of two gaussian functions in pad and time direction.
  /// gaussian function in pad direction assumes that the width of the gaussian is given by the pad-response-function and the transversal diffusion.
  ///
  /// cpad      - pad (y) coordinate
  /// ctime     - time(z) coordinate
  /// ky        - dy/dx
  /// kz        - dz/dx
  /// rmsy0     - RF width. Width is normalized in this function to pad units
  /// rmsz0     - RF width in time bin units. Is normalized in this function to zwidth
  /// effLength - contribution of PRF and diffusion
  /// effDiff   - overwrite diffusion
  // Response function aproximated by convolution of gaussian with angular effect (integral=1)
  //
  // Gaus width sy and sz is determined by RF width and diffusion
  // Integral of Q is equal 1
  // Q max is calculated at position cpad, ctime
  // Example function:
  //  TF1 f1("f1", "AliTPCClusterParamNEW::QmaxCorrection(0,0.5,x,0,0,0.5,0.6)",0,1000)
  //
  // AliTPCParam *param   = AliTPCcalibDB::Instance()->GetParameters();
  // const PadRegionInfo& pregion = mapper.getPadRegionInfo(region);
  const float padLength = param.tpcGeometry.PadHeight(padRow); // length of the pad
  const float padWidth = param.tpcGeometry.PadWidth(padRow);   // width of the pad
  // const float padWidth = 1;

  //=========================== get zwidth =========================
  // static const auto& gasPar = o2::tpc::ParameterGas::Instance();
  // const float fDriftVel  = gasPar.DriftV;
  // const float fSamplRate = 5.; // 10 MHz
  // const float zwidth = fDriftVel * 1/fSamplRate;

  // o2::tpc::GPUCATracking g;
  // g.getPseudoVDrift() // == 0.516f
  const float zwidth = 0.516f; // bin width in z direction
  // auto& elParam = ParameterElectronics::Instance();
  // float vzbin = (elParam.ZbinWidth * gasParam.DriftV); // zBin

  // normalize to pad and zwidth
  rmsy0 /= padWidth;
  rmsz0 /= zwidth;

  // const float wwPitch = 0.25; //param->GetWWPitch(0); FIX ME: value is from aliroot in cm
  // const float wwPitch = 0.028f; // GEM hole pitch in LP cm (is private in ModelGEM.cxx) IS THIS NEEDED ANYMORE?
  // =============== get diffusion constants
  // const float diffT1 = gasPar.DiffT;
  // const float diffL1 = gasPar.DiffL;
  const float diffT1 = 0.0209f;
  const float diffL1 = 0.0221f;

  const float fDriftLength = CAMath::Sqrt(driftDistance); // ctime*zwidth = driftlength

  // effLength: effective length in x direction which takes the
  // padLength: length of a pad
  // wwPitch: the wire geometry. +0.5 pad pitch width (distance between two anode wires) REASON OF USAGE UNSURE: at the beginning and end of a pad (wwPitch),
  //
  // wire geometry for old MWPC readout system
  //   ——————- anode wire
  //   |‾‾‾‾‾|———- anode wire
  //   | pad |———- anode wire
  //   |_____|———- anode wire
  //   ——————- anode wire
  //
  // fDriftLength * diffT1: transverse diffusion
  // effPad: a one dim PRF is assumed with no dependence on the padrow. A shorter integration range is introduced with the effPad~0.5 factor
  // const float effLength = padLength + (wwPitch + fDriftLength * diffT1) * effPad; // TODO check if this is correct
  const float effLength = padLength + (fDriftLength * diffT1) * effPad;

  // diffT: transverse diffusion normalized to padWidth to take into account the different pad types
  // diffL: longitudinal diffusion normalized to padWidth to take into account the different pad types
  const float diffT = fDriftLength * diffT1 / padWidth * effDiff; // effDiff==1
  const float diffL = fDriftLength * diffL1 / zwidth * effDiff;

  // transform angular effect to pad units
  // pky: for tracks with an angle!=0 the charge distribution is further smeared by the tracklength=pky in y-direction
  // pkz: for tracks with an angle!=0 the charge distribution is further smeared by the tracklength=pkz in z-direction
  // ky=dy/dx
  // kz=dz/dx
  float pky = ky * effLength / padWidth; // dy/dx * effLength/padWidth    -> widening of track
  float pkz = kz * effLength / zwidth;

  // position in pad/time units. A shift of 0.5 is added due to the shift produced in the HWClusterer.cxx
  const float py = cpad - static_cast<int>(cpad + 0.5f);       // relative pad position. e.g.: pad=1.7 -> py=0.2
  const float pz = altTime - static_cast<int>(altTime + 0.5f); // relative time position

  // if(py==0) return -1; // return in case of single pad cluster

  //
  // sy: total transversal smearing! response function + diffusion! normalized tp padwidth
  // sz: same as sy but in z (time) direction
  const float sy = CAMath::Sqrt(rmsy0 * rmsy0 + diffT * diffT);
  // const float sy = diffT;
  const float sz = CAMath::Sqrt(rmsz0 * rmsz0 + diffL * diffL); // USE THIS TODO
  // const float sz = diffL; // WARNING TESTING

  float corrVal = 0;

  if (type) {
    const float tau = 160e-3f; // float PeakingTime = 160e-3f;
    //corrVal = GaussConvolutionGamma4(py, pz, pky, pkz, sy, sz, tau); // not used due to performance?
  }
  // const float length = padLength * CAMath::Sqrt(1 + ky * ky + kz * kz); // this correction is already applied in the GPUdEdx.h

  // py,pz: realative cluster and pad position
  // pky,pkz (Ly, Lz): smearing of the width over the length Ly,Lz for tracks with an angle!=0
  // sy,sz (sigmay, sigmaz): sigma value of the gaussian induced charge distribution. PRF+diffusion
  // pkz = 0.01;
  // pky = 0.01;
  if (!type) {
    // corrVal = GaussConvolution(py, pz, pky, pkz, sy, sz);// * length;
    static const float twoPi = CAMath::TwoPi();
    //corrVal = TMath::Gaus(py, 0, sy) / (sy * twoPi);

    // TF1 fGausOne("fGausOne", "TMath::Gaus(x, [1], [0]) / ([0] * sqrt(2*3.141592653589))",-2,2); //sqrt(2*PI) TODO
    // fGausOne.SetParameters(sy,py);
    // corrVal = fGausOne.Integral(-0.5,0.5);

    // TF2 fGaus("fGaus", "TMath::Gaus(x, [0], [1]) * TMath::Gaus(y, [2], [3]) / ([1]*[3]*2*3.141592653589)", -2, 2, -2, 2);
    // static tpc::SAMPAProcessing& sampaProcessing = tpc::SAMPAProcessing::instance();
    // TF1 fGamma4("fGamma4","sampaProcessing.getGamma4(x,0.,1.)",0,1);

    // static TF2 fGaus("fGaus", "TMath::Gaus(x, [0], [1]) * TMath::Gaus(y, [2], [3]) / ([1]*[3]*2*3.141592653589)", -2, 2, -2, 2);
    // fGaus.SetParameters(py, sy, pz, sz);
    // corrVal = fGaus.Integral(-0.5, 0.5, -0.5, 0.5);

    // with angulat effect rel pad postiion ONLY
    // x is in principle pad length: integrate over x from -0.5 to 0.5 integrate over one pad length
    // y is integrated over padwidth from -0.5 to 0.5
    // [0] is relative pad position
    // [1] is sigma
    // [2] is pky

    //full functions
    // std::string relPadWithAngular = "std::exp( - ([0]+x*[2] - y)*([0]+x*[2] - y)/(2*[1]*[1]) - ([3]+x*[5] - z)*([3]+x*[5] - z)/(2*[4]*[4]) ) / (2*3.141592653589*[1]*[4])"

    // two dim
    /*
    static std::string relPadWithAngular = "std::exp( - ([0]+x*[2] - y)*([0]+x*[2] - y)/(2*[1]*[1]) ) / (sqrt(2*3.141592653589)*[1])";
    static TF2 fGausRelPadWithAngular("fGausRelPadWithAngular", relPadWithAngular.data(), -2, 2, -2, 2);
    fGausRelPadWithAngular.SetParameters(py, sy, pky);
    corrVal = std::abs( fGausRelPadWithAngular.Integral(-0.5, 0.5, -0.5, 0.5));
    */

    std::string formularQMax = "std::exp( - ([0]+x*[1] - y)*([0]+x*[1] - y)/(2*[2]*[2]) - ([3]+x*[4] - z)*([3]+x*[4] - z)/(2*[5]*[5]) ) / (2*3.141592653589*[2]*[5])";
    TF3 fQMax("qMax", formularQMax.data(), -2, 2, -2, 2, -2, 2); // should be the correct version
    //    parameters   ty,   pky,   sy,  tz,  pkz,  sz
    fQMax.SetParameters(py,   pky,  sy,  pz,  pkz,  sz); //straight track in pad and time direction and at the pad and time centers
    corrVal = std::abs( fQMax.Integral(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5));

    //take angles into account
    // const float normFac = sy*sz*2*3.141592653589f;
    // TF2 fGaus("fGaus", "exp(-([0]-[1]*x)*([0]-[1]*x)/(2*[2]*[2]) - ([3]-[4]*y)*([3]-[4]*y)/(2*[5]*[5])  )", -2, 2, -2, 2);
    // fGaus.SetParameters(py,pky,sy,pz,pkz,sz);
    // corrVal = fGaus.Integral(-0.5,0.5,-0.5,0.5) / normFac;
  }
  std::array<float, 7> arr{py, pz, pky, pkz, sy, sz, corrVal};
  return arr;
  // return length;
}

GPUd() std::array<float, 7> GPUdEdx::qmaxCorrection(const GPUParam& param, int padRow, float cpad, float driftDistance, float ky, float kz, float rmsy0, float rmsz0, float effPad, float effDiff, float ctime)
{
  /// This function approximates the charge distribution by assuming a convolution of two gaussian functions in pad and time direction.
  /// gaussian function in pad direction assumes that the width of the gaussian is given by the pad-response-function and the transversal diffusion.
  ///
  /// cpad      - pad (y) coordinate
  /// ctime     - time(z) coordinate
  /// ky        - dy/dx
  /// kz        - dz/dx
  /// rmsy0     - RF width. Width is normalized in this function to pad units
  /// rmsz0     - RF width in time bin units. Is normalized in this function to zwidth
  /// effLength - contribution of PRF and diffusion
  /// effDiff   - overwrite diffusion
  // Response function aproximated by convolution of gaussian with angular effect (integral=1)
  //
  // Gaus width sy and sz is determined by RF width and diffusion
  // Integral of Q is equal 1
  // Q max is calculated at position cpad, ctime
  // Example function:
  //  TF1 f1("f1", "AliTPCClusterParamNEW::QmaxCorrection(0,0.5,x,0,0,0.5,0.6)",0,1000)
  //
  // AliTPCParam *param   = AliTPCcalibDB::Instance()->GetParameters();
  // const PadRegionInfo& pregion = mapper.getPadRegionInfo(region);
  const float padLength = param.tpcGeometry.PadHeight(padRow); // length of the pad
  const float padWidth = param.tpcGeometry.PadWidth(padRow);   // width of the pad
  // const float padWidth = 1;

  //=========================== get zwidth =========================
  // static const auto& gasPar = o2::tpc::ParameterGas::Instance();
  // const float fDriftVel  = gasPar.DriftV;
  // const float fSamplRate = 5.; // 10 MHz
  // const float zwidth = fDriftVel * 1/fSamplRate;

  // o2::tpc::GPUCATracking g;
  // g.getPseudoVDrift() // == 0.516f
  const float zwidth = 0.516f; // bin width in z direction
  // auto& elParam = ParameterElectronics::Instance();
  // float vzbin = (elParam.ZbinWidth * gasParam.DriftV); // zBin

  // normalize to pad and zwidth
  rmsy0 /= padWidth;
  rmsz0 /= zwidth;

  // const float wwPitch = 0.25; //param->GetWWPitch(0); FIX ME: value is from aliroot in cm
  // const float wwPitch = 0.028f; // GEM hole pitch in LP cm (is private in ModelGEM.cxx) IS THIS NEEDED ANYMORE?
  // =============== get diffusion constants
  // const float diffT1 = gasPar.DiffT;
  // const float diffL1 = gasPar.DiffL;
  const float diffT1 = 0.0209f;
  const float diffL1 = 0.0221f;

  const float fDriftLength = CAMath::Sqrt(driftDistance); // ctime*zwidth = driftlength

  // effLength: effective length in x direction which takes the
  // padLength: length of a pad
  // wwPitch: the wire geometry. +0.5 pad pitch width (distance between two anode wires) REASON OF USAGE UNSURE: at the beginning and end of a pad (wwPitch),
  //
  // wire geometry for old MWPC readout system
  //   ——————- anode wire
  //   |‾‾‾‾‾|———- anode wire
  //   | pad |———- anode wire
  //   |_____|———- anode wire
  //   ——————- anode wire
  //
  // fDriftLength * diffT1: transverse diffusion
  // effPad: a one dim PRF is assumed with no dependence on the padrow. A shorter integration range is introduced with the effPad~0.5 factor
  // const float effLength = padLength + (wwPitch + fDriftLength * diffT1) * effPad; // TODO check if this is correct
  const float effLength = padLength + (fDriftLength * diffT1) * effPad;

  // diffT: transverse diffusion normalized to padWidth to take into account the different pad types
  // diffL: longitudinal diffusion normalized to padWidth to take into account the different pad types
  const float diffT = fDriftLength * diffT1 / padWidth * effDiff; // effDiff==1
  const float diffL = fDriftLength * diffL1 / zwidth * effDiff;

  // transform angular effect to pad units
  // pky: for tracks with an angle!=0 the charge distribution is further smeared by the tracklength=pky in y-direction
  // pkz: for tracks with an angle!=0 the charge distribution is further smeared by the tracklength=pkz in z-direction
  // ky=dy/dx
  // kz=dz/dx
  float pky = ky * effLength / padWidth; // dy/dx * effLength/padWidth    -> widening of track
  float pkz = kz * effLength / zwidth;

  // position in pad/time units. A shift of 0.5 is added due to the shift produced in the HWClusterer.cxx
  const float py = cpad - static_cast<int>(cpad + 0.5f);   // relative pad position. e.g.: pad=1.7 -> py=0.2
  const float pz = ctime - static_cast<int>(ctime + 0.5f); // relative time position

  // if(py==0) return -1; // return in case of single pad cluster

  //
  // sy: total transversal smearing! response function + diffusion! normalized tp padwidth
  // sz: same as sy but in z (time) direction
  const float sy = CAMath::Sqrt(rmsy0 * rmsy0 + diffT * diffT);
  // const float sy = diffT;
  const float sz = CAMath::Sqrt(rmsz0 * rmsz0 + diffL * diffL); // USE THIS TODO
  // const float sz = diffL; // WARNING TESTING

  float corrVal = 0;

  // if (type) {
  //   const float tau = 160e-3f;                                       // float PeakingTime = 160e-3f;
  //   corrVal = GaussConvolutionGamma4(py, pz, pky, pkz, sy, sz, tau); // not used due to performance?
  // }
  // const float length = padLength * CAMath::Sqrt(1 + ky * ky + kz * kz); // this correction is already applied in the GPUdEdx.h

  // py,pz: realative cluster and pad position
  // pky,pkz (Ly, Lz): smearing of the width over the length Ly,Lz for tracks with an angle!=0
  // sy,sz (sigmay, sigmaz): sigma value of the gaussian induced charge distribution. PRF+diffusion
  // pkz = 0.01;
  // pky = 0.01;
  // if (!type) {
  corrVal = GaussConvolution(py, pz, pky, pkz, sy, sz); // * length;
  // }
  std::array<float, 7> arr{py, pz, pky, pkz, sy, sz, corrVal};
  return arr;
  // return length;
}

GPUd() float GPUdEdx::GaussConvolution(float x0, float x1, float k0, float k1, float s0, float s1)
{
  /// 2 D gaus convoluted with angular effect
  /// See in mathematica:
  /// Simplify[Integrate[Exp[-(x0-k0*xd)*(x0-k0*xd)/(2*s0*s0)-(x1-k1*xd)*(x1-k1*xd)/(2*s1*s1)]/(s0*s1),{xd,-1/2,1/2}]]
  /// WolframAlpha: 'integrate exp[-(x0-k0*x)^2/(2*s0^2)-(x1-k1*x)^2/(2*s1^2)]/(s0*s1) dx from x=-0.5 to x=0.5'
  /// TF1 f1("f1","GPUdEdx::GaussConvolution(x,0,1,0,0.1,0.1)",-2,2)
  /// TF2 f2("f2","AliTPCClusterParamNEW::GaussConvolution(x,y,1,1,0.1,0.1)",-2,2,-2,2)
  /// checked this function should be correct :)

  static const float kEpsilon = 0.0001f;
  static const float twoPi = CAMath::TwoPi();
  static const float hnorm = 0.5f / CAMath::Sqrt(twoPi);
  static const float sqtwo = CAMath::Sqrt(2.);

  if ((CAMath::Abs(k0) + CAMath::Abs(k1)) < kEpsilon * (s0 + s1)) {
    // small angular effect
    const float val = TMath::Gaus(x0, 0, s0) * TMath::Gaus(x1, 0, s1) / (s0 * s1 * twoPi);
    return val;
  }
  const float sigma2 = k1 * k1 * s0 * s0 + k0 * k0 * s1 * s1;
  const float sigma = CAMath::Sqrt(sigma2);
  const float exp0 = std::exp(-(k1 * x0 - k0 * x1) * (k1 * x0 - k0 * x1) / (2.f * sigma2));
  //
  const float sigmaErf = 1.f / (2.f * s0 * s1 * sqtwo * sigma);
  const float k0s1s1 = 2.f * k0 * s1 * s1;
  const float k1s0s0 = 2.f * k1 * s0 * s0;
  const float erf0 = ErfFast((sigma2 - k0s1s1 * x0 - k1s0s0 * x1) * sigmaErf);
  const float erf1 = ErfFast((sigma2 + k0s1s1 * x0 + k1s0s0 * x1) * sigmaErf);
  const float norm = hnorm / sigma;
  const float val = norm * exp0 * (erf0 + erf1);
  return val;
}

GPUd() float GPUdEdx::ErfcFast(float x)
{
  // Fast implementation of the complementary error function
  // The error of the approximation is |eps(x)| < 5E-4
  // See Abramowitz and Stegun, p.299, 7.1.27
  // http://people.math.sfu.ca/~cbm/aands/page_299.htm

  const float z = CAMath::Abs(x);
  float ans = 1 + z * (0.278393f + z * (0.230389f + z * (0.000972f + z * 0.078108f)));
  ans = 1.0f / ans;
  ans *= ans;
  ans *= ans;

  return (x >= 0.0 ? ans : 2.0f - ans);
}

float GPUdEdx::GaussConvolutionGamma4(float x0, float x1, float k0, float k1, float s0, float s1, float tau)
{
  /// 2 D gaus convoluted with angular effect and exponential tail in z-direction
  /// tail integrated numerically
  /// Integral normalized to one
  /// Mean at 0
  ///
  /// TF1 f1g4("f1g4","AliTPCClusterParamNEW::GaussConvolutionGamma4(0,x,0,0,0.5,0.2,1.6)",-5,5)
  /// TF2 f2g4("f2g4","AliTPCClusterParamNEW::GaussConvolutionGamma4(y,x,0,0,0.5,0.2,1.6)",-5,5,-5,5)

  /*
  Vc::float_v adcValue = 55.f * ADC * Vc::exp(-4.f * (time - startTime) / eleParam.PeakingTime) *
                        (time - startTime) / eleParam.PeakingTime * (time - startTime) /
                        eleParam.PeakingTime * (time - startTime) / eleParam.PeakingTime *
                        (time - startTime) / eleParam.PeakingTime;
  */

  float sum = 0;
  float mean = 0;

  // the COG of G4
  for (float iexp = 0; iexp < 5; iexp += 0.2f) {
    const float g4 = std::exp(-4. * iexp / tau) * std::pow(iexp / tau, 4.);
    mean += iexp * g4;
    sum += g4;
  }
  mean /= sum;
  //
  sum = 0;
  float val = 0;
  for (float iexp = 0; iexp < 5; iexp += 0.2f) {
    float g4 = std::exp(-4. * iexp / tau) * std::pow(iexp / tau, 4.);
    val += GaussConvolution(x0, x1 + mean - iexp, k0, k1, s0, s1) * g4;
    sum += g4;
  }
  return val / sum;
}
