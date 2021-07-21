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

#include "TPCCalibration/IDCFactorizationBase.h"

o2::tpc::IDCFactorizationSector::IDCFactorizationBase(const std::array<unsigned char, Mapper::NREGIONS>& groupPads, const std::array<unsigned char, Mapper::NREGIONS>& groupRows, const std::array<unsigned char, Mapper::NREGIONS>& groupLastRowsThreshold, const std::array<unsigned char, Mapper::NREGIONS>& groupLastPadsThreshold, const unsigned int timeFrames, const unsigned int timeframesDeltaIDC)
  : IDCGroupHelperSector{groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold}, mTimeFrames{timeFrames}, mTimeFramesDeltaIDC{timeframesDeltaIDC}, mIDCDelta{timeFrames / timeframesDeltaIDC + (timeFrames % timeframesDeltaIDC != 0)}
{
  for (auto& idc : mIDCs) {
    idc.resize(mTimeFrames);
  }
}

unsigned long o2::tpc::IDCFactorizationSector::getNIntegrationIntervals(const unsigned int chunk) const
{
  std::size_t sum = 0;
  const auto firstTF = chunk * getNTFsPerChunk(0);
  for (unsigned int i = firstTF; i < firstTF + getNTFsPerChunk(chunk); ++i) {
    sum += getSizeIDCs(0, i);
  }
  return sum / getIDCsPerCRU(0);
}

unsigned long o2::tpc::IDCFactorizationSector::getNIntegrationIntervals() const
{
  std::size_t sum = 0;
  for (const auto& idcsTF : mIDCs[0]) {
    sum += idcsTF.size();
  }
  return sum / getIDCsPerCRU(0);
}

unsigned int o2::tpc::IDCFactorizationSector::getNTFsPerChunk(const unsigned int chunk) const
{
  const unsigned int remain = mTimeFrames % mTimeFramesDeltaIDC;
  return ((chunk == getNChunks() - 1) && remain) ? remain : mTimeFramesDeltaIDC;
}

void o2::tpc::IDCFactorizationSector::reset()
{
  for (auto& tf : mIDCs) {
    for (auto& idcs : tf) {
      idcs.clear();
    }
  }
}

unsigned long o2::tpc::IDCFactorizationRegion::getNIntegrationIntervals() const
{
  std::size_t sum = 0;
  for (const auto& idcsTF : mIDCs) {
    sum += idcsTF.size();
  }
  return sum / mNIDCsPerCRU;
}

void o2::tpc::IDCFactorizationSector::calcIDCZero()
{
  const unsigned int nIDCsSide = mNIDCsPerSector * o2::tpc::SECTORSPERSIDE;
  mIDCZero.clear();
  mIDCZero.resize(nIDCsSide);

#pragma omp parallel for num_threads(sNThreads)
  for (unsigned int cru = 0; cru < mIDCs.size(); ++cru) {
    const o2::tpc::CRU cruTmp(cru);
    const unsigned int region = cruTmp.region();
    const auto factorIndexGlob = mRegionOffs[region] + mNIDCsPerSector * cruTmp.sector();

    for (unsigned int timeframe = 0; timeframe < mTimeFrames; ++timeframe) {
      for (unsigned int idcs = 0; idcs < getSizeIDCs(cru, timeframe); ++idcs) {
        const unsigned int indexGlob = (idcs % getIDCsPerCRU(region)) + factorIndexGlob;
        mIDCZero.fillValueIDCZero(mIDCs[cru][timeframe][idcs], cruTmp.side(), indexGlob % nIDCsSide);
      }
    }
  }
  std::transform(mIDCZero.mIDCZero[Side::A].begin(), mIDCZero.mIDCZero[Side::A].end(), mIDCZero.mIDCZero[Side::A].begin(), [norm = getNIntegrationIntervals()](auto& val) { return val / norm; });
  std::transform(mIDCZero.mIDCZero[Side::C].begin(), mIDCZero.mIDCZero[Side::C].end(), mIDCZero.mIDCZero[Side::C].begin(), [norm = getNIntegrationIntervals()](auto& val) { return val / norm; });
}

void o2::tpc::IDCFactorizationRegion::calcIDCZero()
{
  mIDCZero.clear();
  mIDCZero.resize(mNIDCsPerCRU);
#pragma omp parallel for num_threads(sNThreads)
  for (unsigned int timeframe = 0; timeframe < mTimeFrames; ++timeframe) {
    for (unsigned int idcs = 0; idcs < getSizeIDCs(0, timeframe); ++idcs) {
      const unsigned int indexGlob = idcs % mNIDCsPerCRU;
      mIDCZero[indexGlob] += mIDCs[timeframe][idcs];
    }
  }
  std::transform(mIDCZero.begin(), mIDCZero.end(), mIDCZero.begin(), [norm = getNIntegrationIntervals()](auto& val) { return val / norm; });
}

void o2::tpc::IDCFactorizationSector::calcIDCOne()
{
  const unsigned int nIDCsSide = mNIDCsPerSector * SECTORSPERSIDE;
  const unsigned int integrationIntervals = getNIntegrationIntervals();
  mIDCOne.clear();
  mIDCOne.resize(integrationIntervals);
  const unsigned int crusPerSide = Mapper::NREGIONS * SECTORSPERSIDE;

#pragma omp parallel for num_threads(sNThreads)
  for (unsigned int cru = 0; cru < mIDCs.size(); ++cru) {
    const o2::tpc::CRU cruTmp(cru);
    const unsigned int region = cruTmp.region();
    const auto factorIDCOne = crusPerSide * getIDCsPerCRU(region);
    const auto side = cruTmp.side();
    const auto factorIndexGlob = mRegionOffs[region] + mNIDCsPerSector * cruTmp.sector();
    unsigned int integrationIntervallast = 0;
    for (unsigned int timeframe = 0; timeframe < mTimeFrames; ++timeframe) {
      for (unsigned int idcs = 0; idcs < getSizeIDCs(cru, timeframe); ++idcs) {
        const unsigned int integrationInterval = idcs / getIDCsPerCRU(region) + integrationIntervallast;
        const unsigned int indexGlob = (idcs % getIDCsPerCRU(region)) + factorIndexGlob;
        mIDCOne.mIDCOne[side][integrationInterval] += mIDCs[cru][timeframe][idcs] / (factorIDCOne * mIDCZero.mIDCZero[side][indexGlob % nIDCsSide]);
      }
      integrationIntervallast += getSizeIDCs(cru, timeframe) / getIDCsPerCRU(region);
    }
  }
}

void o2::tpc::IDCFactorizationRegion::calcIDCOne()
{
  const unsigned int nIDCs = mNIDCsPerCRU;
  mIDCOne.clear();
  mIDCOne.resize(getNIntegrationIntervals());
  unsigned int integrationIntervallast = 0;
  for (unsigned int timeframe = 0; timeframe < mTimeFrames; ++timeframe) {
    for (unsigned int idcs = 0; idcs < getSizeIDCs(0, timeframe); ++idcs) {
      const unsigned int integrationInterval = idcs / nIDCs + integrationIntervallast;
      const unsigned int indexGlob = idcs % nIDCs;
      mIDCOne[integrationInterval] += mIDCs[timeframe][idcs] / (nIDCs * mIDCZero[indexGlob]);
    }
    integrationIntervallast += getSizeIDCs(0, timeframe) / nIDCs;
  }
}

void o2::tpc::IDCFactorizationRegion::calcIDCDelta()
{
  const unsigned int nIDCs = mNIDCsPerCRU;
  mIDCDelta.resize(nIDCs * getNIntegrationIntervals());
  unsigned int integrationIntervallast = 0;
  for (unsigned int timeframe = 0; timeframe < mTimeFrames; ++timeframe) {
    for (unsigned int idcs = 0; idcs < getSizeIDCs(0, timeframe); ++idcs) {
      const unsigned int indexGlob = idcs % mNIDCsPerCRU;
      const auto idcZero = mIDCZero[indexGlob];
      const unsigned int intervallocal = idcs / nIDCs;
      const unsigned int integrationIntervalGlobal = intervallocal + integrationIntervallast;
      const auto idcOne = mIDCOne[integrationIntervalGlobal];
      const auto mult = idcZero * idcOne;
      const auto val = (mult > 0) ? mIDCs[timeframe][idcs] / mult : 0;
      mIDCDelta[indexGlob] = val - 1.f;
    }
    integrationIntervallast += getSizeIDCs(0, timeframe) / nIDCs;
  }
}

void o2::tpc::IDCFactorizationSector::calcIDCDelta()
{
  const unsigned int nIDCsSide = mNIDCsPerSector * SECTORSPERSIDE;
  for (unsigned int i = 0; i < getNChunks(); ++i) {
    const auto idcsSide = nIDCsSide * IDCFactorizationBase<IDCGroupHelperSector>::getNIntegrationIntervals(i);
    mIDCDelta[i].getIDCDelta(Side::A).resize(idcsSide);
    mIDCDelta[i].getIDCDelta(Side::C).resize(idcsSide);
  }

#pragma omp parallel for num_threads(sNThreads)
  for (unsigned int cru = 0; cru < mIDCs.size(); ++cru) {
    const o2::tpc::CRU cruTmp(cru);
    const unsigned int region = cruTmp.region();
    const auto side = cruTmp.side();
    const auto factorIndexGlob = mRegionOffs[region] + mNIDCsPerSector * cruTmp.sector();
    unsigned int integrationIntervallast = 0;
    unsigned int integrationIntervallastLocal = 0;
    unsigned int lastChunk = 0;

    for (unsigned int timeframe = 0; timeframe < mTimeFrames; ++timeframe) {
      const unsigned int chunk = getChunk(timeframe);
      if (lastChunk != chunk) {
        integrationIntervallastLocal = 0;
      }

      for (unsigned int idcs = 0; idcs < getSizeIDCs(cru, timeframe); ++idcs) {
        const unsigned int intervallocal = idcs / getIDCsPerCRU(region);
        const unsigned int integrationIntervalGlobal = intervallocal + integrationIntervallast;
        const unsigned int integrationIntervalLocal = intervallocal + integrationIntervallastLocal;
        const unsigned int indexGlob = (idcs % getIDCsPerCRU(region)) + factorIndexGlob;
        const unsigned int indexGlobMod = indexGlob % nIDCsSide;
        const auto idcZero = mIDCZero.mIDCZero[side][indexGlobMod];
        const auto idcOne = mIDCOne.mIDCOne[side][integrationIntervalGlobal];
        const auto mult = idcZero * idcOne;
        const auto val = (mult > 0) ? mIDCs[cru][timeframe][idcs] / mult : 0;
        mIDCDelta[chunk].setIDCDelta(side, indexGlobMod + integrationIntervalLocal * nIDCsSide, val - 1.f);
      }

      const unsigned int intervals = getSizeIDCs(cru, timeframe) / getIDCsPerCRU(region);
      integrationIntervallast += intervals;
      integrationIntervallastLocal += intervals;
      lastChunk = chunk;
    }
  }
}
