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
/// \brief Selection of the V0s written with their clusters by o2-sv-study-workflow --full-output

#ifndef O2_SV_STUDY_V0_FULL_OUTPUT_PARAMS_H
#define O2_SV_STUDY_V0_FULL_OUTPUT_PARAMS_H

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"

namespace o2::svstudy
{
/// Prong cuts are applied before the refit of a candidate (they need only its tracks), V0 cuts after it. A V0 is written if it passes
/// the cuts of at least one type (K0s, Lambda, anti-Lambda); the types it passed are flagged in V0WithClusters::selection.
struct V0FullOutputParams : o2::conf::ConfigurableParamHelper<V0FullOutputParams> {
  // prong cuts, before the refit
  bool requireOppositeCharges = true; ///< prongs of opposite charge
  bool requireTPCBoth = true;         ///< both prongs with a TPC track
  int minTPCClusters = 0;             ///< min. number of TPC clusters of each prong (if it has a TPC track)
  // V0 cuts for all types, after the refit
  float maxDCA = 0.01f;    ///< max. V0::getDCA (chi2 at the point of closest approach of the prongs)
  float minCosPA = 0.999f; ///< min. cosine of the pointing angle to the primary vertex
  float minQt = 0.01f;     ///< min. Armenteros qT [GeV]
  float maxQt = 0.2f;      ///< max. Armenteros qT [GeV]
  float minR = 5.f;        ///< min. transverse decay radius [cm]
  float maxR = 35.f;       ///< max. transverse decay radius [cm]
  // per type
  bool selK0 = true;             ///< select K0s
  float maxDMassK0 = 0.005f;     ///< max. |m(pi+ pi-) - m(K0s)| [GeV]
  bool selLambda = true;         ///< select Lambda and anti-Lambda
  float maxDMassLambda = 0.005f; ///< max. |m(p pi) - m(Lambda)| [GeV]
  float minAlphaLambda = 0.28f;  ///< min. Armenteros alpha of a Lambda (anti-Lambda: max. -minAlphaLambda)
  float maxAlphaLambda = 0.7f;   ///< max. Armenteros alpha of a Lambda (anti-Lambda: min. -maxAlphaLambda)
  int maxV0sPerTF = -1;          ///< write at most this many V0s per TF, in the order of the V0 finder (< 0: no limit)

  O2ParamDef(V0FullOutputParams, "svFullOutput");
};
} // namespace o2::svstudy

#endif
