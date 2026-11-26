// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file V0FullOutputParams.h
/// \brief Selection and downsampling of the V0s written with their clusters by o2-sv-study-workflow --full-output

#ifndef O2_SV_STUDY_V0_FULL_OUTPUT_PARAMS_H
#define O2_SV_STUDY_V0_FULL_OUTPUT_PARAMS_H

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"

namespace o2::svstudy
{
/// Prong cuts are applied before the refit of a candidate (they need only its tracks), V0 cuts after it. A V0 passing the cuts of a type
/// (K0s, Lambda, anti-Lambda, photon) is assigned to one sampling class (in this order of priority); the types it passed are flagged in
/// V0WithClusters::selection. Per TF and class the V0s are then sampled flat in V0 pT: in each of nPtBins log bins between ptMin and ptMax
/// (plus one bin below and one above) at most maxPerBin* V0s are kept, chosen at random (reproducible: hash of the TF and the V0 index);
/// each kept V0 gets the weight (V0s of its class and bin) / (kept ones). Far V0s (bit Far) are sampled apart from the others of their
/// class, with maxPerBinFar (default: all kept), and are not subject to maxV0sPerTF.
struct V0FullOutputParams : o2::conf::ConfigurableParamHelper<V0FullOutputParams> {
  // prong cuts, before the refit
  bool requireOppositeCharges = true; ///< prongs of opposite charge
  bool requireTPCBoth = true;         ///< both prongs with a TPC track
  int minTPCClusters = 0;             ///< min. number of TPC clusters of each prong (if it has a TPC track)
  // V0 cuts for K0s and Lambda, after the refit
  float maxDCA = 0.01f;             ///< max. V0::getDCA (chi2 at the point of closest approach of the prongs)
  float minCosPA = 0.9998f;         ///< min. cosine of the pointing angle to the primary vertex (Pb-Pb data, sidebands: 0.999 -> 0.9998
                                    ///< keeps 98 % of the K0s / 95 % of the Lambdas, S/B K0s 44 -> 61, Lambda 13 -> 20)
  float minQt = 0.01f;              ///< min. Armenteros qT [GeV]
  float maxQt = 0.25f;              ///< max. Armenteros qT [GeV] (K0s endpoint 0.206 GeV: a cut at 0.2 removes ~24 % of the K0s)
  float minR = 5.f;                 ///< min. transverse decay radius [cm]
  float maxR = 35.f;                ///< max. transverse decay radius [cm]
  bool selK0 = true;                ///< select K0s
  float maxDMassK0 = 0.010f;        ///< max. |m(pi+ pi-) - m(K0s)| [GeV] (Pb-Pb data: the peak, sigma ~4 MeV with tails, has 74 % of
                                    ///< its signal within 5 MeV, 94 % within 10 MeV; S/B 27 at 10 MeV)
  float vetoDMassLambdaK0 = 0.005f; ///< K0s within maxR: not a K0s if |m(p pi-) - m(Lambda)| or |m(pbar pi+) - m(Lambda)| is below this [GeV]
                                    ///< (Lambda reflection at qT < 0.1 GeV, |alpha| > 0.6: ~5 % of the K0s in Pb-Pb data, ~1.5 % K0s lost;
                                    ///< < 0: no veto). Far K0s: minQtK0Far removes the reflection
  bool selLambda = true;            ///< select Lambda and anti-Lambda
  float maxDMassLambda = 0.005f;    ///< max. |m(p pi) - m(Lambda)| [GeV]
  float minAlphaLambda = 0.28f;     ///< min. Armenteros alpha of a Lambda (anti-Lambda: max. -minAlphaLambda)
  float maxAlphaLambda = 0.95f;     ///< max. Armenteros alpha of a Lambda (anti-Lambda: min. -maxAlphaLambda); 0.95 = the whole Lambda
                                    ///< ellipse (fast Lambdas reach 0.87); with vetoDMassK0Lambda: +46 % Lambdas vs 0.7, S/B 10.9 -> 12.6 (Pb-Pb data)
  float vetoDMassK0Lambda = 0.005f; ///< Lambda / anti-Lambda within maxR: not a Lambda if |m(pi+ pi-) - m(K0s)| is below this [GeV] (the K0s
                                    ///< reflection where the two ellipses cross, |alpha| 0.7-0.85; < 0: no veto)
  // far K0s and Lambda: beyond maxR, both legs without ITS (V0 finder's svertexer.farV0Mode), with these cuts instead of the ones above.
  // Defaults tuned on Pb-Pb data at 19 kHz with mass sidebands (signal / background in the peak: K0s ~11, Lambda ~20; before ~0.3 / ~1)
  bool selFar = true;              ///< select far K0s / Lambda
  float maxRFar = 250.f;           ///< max. transverse decay radius [cm] of far V0s
  float maxDCAFar = 0.1f;          ///< max. V0::getDCA of far V0s (MC: 90 % of the true ones below 0.07, half of the fakes above 0.5)
  float minQtFar = 0.06f;          ///< min. Armenteros qT [GeV] of far V0s (the background piles up at low qT)
  float minQtK0Far = 0.105f;       ///< min. Armenteros qT [GeV] of far K0s: above the Lambda maximum (0.101 GeV), removes the Lambda
                                   ///< reflection in the K0s window; keeps the K0s decays with |cos(theta*)| < ~0.86
  float minCosPAFar = 0.9998f;     ///< min. cosine of the pointing angle of far V0s
  int minTPCClustersFar = 80;      ///< min. number of TPC clusters of each leg of a far V0 (short background legs)
  float maxCtK0Far = 20.f;         ///< max. proper decay length L m / p [cm] of far K0s (L: 3D distance from the primary vertex; c tau 2.68 cm)
  float maxCtLambdaFar = 40.f;     ///< max. proper decay length L m / p [cm] of far Lambdas (c tau 7.89 cm)
  float maxDMassK0Far = 0.03f;     ///< max. |m(pi+ pi-) - m(K0s)| [GeV] of far V0s (peak in data: mean -3 MeV, sigma 9 MeV)
  float maxDMassLambdaFar = 0.01f; ///< max. |m(p pi) - m(Lambda)| [GeV] of far V0s (peak in data: sigma 3 MeV)
  float maxAlphaLambdaFar = 0.7f;  ///< max. Armenteros alpha of a far Lambda (anti-Lambda: min. -maxAlphaLambdaFar); the far cuts were tuned
                                   ///< with 0.7, larger |alpha| not measured for far V0s
  // photon conversions, after the refit (only these cuts apply to them)
  bool selGamma = true;          ///< select photon conversions: V0 finder's photon-only flag or m(e+ e-) < maxMassGamma
  float maxMassGamma = 0.05f;    ///< max. m(e+ e-) [GeV]
  float maxQtGamma = 0.01f;      ///< Armenteros cut as the O2Physics photon-conversion-builder: (alpha / maxAlphaGamma)^2 + (qT / maxQtGamma)^2
                                 ///< < 1, maxQtGamma [GeV] (Pb-Pb data, dE/dx-tagged legs: purity 0.65 -> 0.85, 87 % of the conversions kept)
  float maxAlphaGamma = 0.95f;   ///< alpha of that ellipse
  float minLegDCAXYGamma = 0.1f; ///< min. |DCA_xy| [cm] of both legs to the primary vertex, as in O2Physics (< 0: no cut; a leg whose propagation
                                 ///< to the vertex fails passes)
  float minRGamma = 1.f;         ///< min. transverse conversion radius [cm]
  float maxRGamma = 90.f;        ///< max. transverse conversion radius [cm] (TPC inner field cage)
  float minCosPAGamma = 0.99f;   ///< min. cosine of the pointing angle (loose: the photon pointing is an observable)
  float maxDCAGamma = -1.f;      ///< max. V0::getDCA (< 0: no cut beyond the V0 finder's)
  // downsampling flat in V0 pT, per class (< 0: keep all)
  int nPtBins = 20;           ///< log V0-pT bins between ptMin and ptMax
  float ptMin = 0.1f;         ///< [GeV]
  float ptMax = 10.f;         ///< [GeV]
  int maxPerBinK0 = 5;        ///< K0s kept per TF and pT bin
  int maxPerBinLambda = 5;    ///< Lambdas kept per TF and pT bin, and the same for the anti-Lambdas
  int maxPerBinGamma = 5;     ///< photons kept per TF and pT bin
  int maxPerBinFar = -1;      ///< far K0s / Lambdas (bit Far) kept per TF and pT bin, sampled apart from the other V0s of their class (< 0: all)
  int maxV0sPerTF = -1;       ///< after the sampling at most this many V0s per TF that are not far, chosen at random (< 0: no limit)
  float pvTimeNSigma = 4.f;   ///< nCompatPV: primary vertices within this many sigma of their time error of the prongs' time windows
  float legScoreWeight = 0.f; ///< weight of the pointing angle in legScore (1 - weight: leg distance), as deduplicationScoreWeight in O2Physics
                              ///< (0.5 there); 0: the leg distance alone picks the true photon most often (MC)

  O2ParamDef(V0FullOutputParams, "svFullOutput");
};
} // namespace o2::svstudy

#endif
