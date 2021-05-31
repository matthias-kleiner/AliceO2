// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  testO2TPCIDCFourierTransform.cxx
/// \brief this task tests  the calculation of fourier coefficients
///
/// \author  Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#define BOOST_TEST_MODULE Test TPC O2TPCIDCFourierTransform class
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "TPCCalibration/IDCFourierTransform.h"
#include "TRandom.h"
#include <numeric>

namespace o2::tpc
{

static constexpr float ABSTOLERANCE = 0.01f; // absolute tolerance is taken at small values near 0
static constexpr float TOLERANCE = 0.2f;

o2::tpc::IDCOne getIDCsOne(const std::vector<unsigned int>& integrationIntervals)
{
  const int nIDCs = std::accumulate(integrationIntervals.begin(), integrationIntervals.end(), 0);
  o2::tpc::IDCOne idcsOut;
  for (int iside = 0; iside < 2; ++iside) {
    std::vector<float> idcs(nIDCs);
    for (auto& val : idcs) {
      val = gRandom->Gaus(1, 0.2);
    }
    idcsOut.mIDCOne[iside] = std::move(idcs);
  }
  return idcsOut;
}

std::vector<unsigned int> getIntegrationIntervalsPerTF(const unsigned int integrationIntervals, const unsigned int tfs)
{
  std::vector<unsigned int> intervals;
  intervals.reserve(tfs);
  for (unsigned int i = 0; i < tfs; ++i) {
    const unsigned int additionalInterval = (i % 3) ? 1 : 0; // in each integration inerval are either 10 or 11 values when having 128 orbits per TF and 12 orbits integration length
    intervals.emplace_back(integrationIntervals + additionalInterval);
  }
  return intervals;
}

BOOST_AUTO_TEST_CASE(IDCFourierTransform_test)
{
  const unsigned int integrationIntervals = 10; // number of integration intervals for first TF
  const unsigned int tfs = 200;                 // number of aggregated TFs
  const unsigned int rangeIDC = 200;            // number of IDCs used to calculate the fourier coefficients
  gRandom->SetSeed(1);

  for (int iType = 0; iType < 2; ++iType) {
    const bool fft = iType == 0 ? false : true;
    o2::tpc::IDCFourierTransform::setFFT(fft);
    o2::tpc::IDCFourierTransform idcFourierTransform{rangeIDC, tfs};
    const auto intervalsPerTF = getIntegrationIntervalsPerTF(integrationIntervals, tfs);
    idcFourierTransform.setIDCs(getIDCsOne(intervalsPerTF), intervalsPerTF);
    idcFourierTransform.setIDCs(getIDCsOne(intervalsPerTF), intervalsPerTF);
    idcFourierTransform.calcFourierCoefficients();

    const std::vector<unsigned int> offsetIndex = idcFourierTransform.getLastIntervals();
    for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
      const o2::tpc::Side side = iSide == 0 ? Side::A : Side::C;
      const auto idcOneExpanded = idcFourierTransform.getExpandedIDCOne(side);
      const auto inverseFourierFFTW3 = idcFourierTransform.inverseFourierTransformFFTW3(side);
      for (unsigned int coeff = 0; coeff < idcFourierTransform.getNCoefficients(); ++coeff) {
        for (unsigned int interval = 0; interval < idcFourierTransform.getNIntervals(); ++interval) {
          unsigned int indexIDFT = 0;
          for (unsigned int index = 0; index < rangeIDC; ++index) {
            const float origIDCOne = idcOneExpanded[index + offsetIndex[interval]];
            const float idftIDCOne = inverseFourierFFTW3[interval][indexIDFT++];
            if (std::fabs(origIDCOne) < ABSTOLERANCE) {
              BOOST_CHECK_SMALL(idftIDCOne - origIDCOne, ABSTOLERANCE);
            } else {
              BOOST_CHECK_CLOSE(idftIDCOne, origIDCOne, TOLERANCE);
            }
          }
        }
      }
    }
  }
}

} // namespace o2::tpc
