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

#include "TPCCalibration/IDCFactorization.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "Framework/Logger.h"
#include "TPCCalibration/IDCDrawHelper.h"
#include "TFile.h"
#include <functional>

#if (defined(WITH_OPENMP) || defined(_OPENMP)) && !defined(__CLING__)
#include <omp.h>
#endif

template <typename Type>
void o2::tpc::IDCFactorization<Type>::dumpToFile(const char* outFileName, const char* outName) const
{
  TFile fOut(outFileName, "RECREATE");
  fOut.WriteObject(this, outName);
  fOut.Close();
}

template <typename Type>
void o2::tpc::IDCFactorization<Type>::dumpToTree(int integrationIntervals, const char* outFileName) const
{
  const Mapper& mapper = Mapper::instance();
  o2::utils::TreeStreamRedirector pcstream(outFileName, "RECREATE");
  pcstream.GetFile()->cd();

  if (this->mIDCZero.empty()) {
    LOGP(warning, "call factorizeIDCs() first. returning");
    return;
  }

  if (integrationIntervals <= 0) {
    integrationIntervals = static_cast<int>(this->getNIntegrationIntervals());
  }

  std::vector<float> idcOneA;
  std::vector<float> idcOneC;
  if constexpr (std::is_same_v<Type, IDCGroupHelperSector>) {
    idcOneA = this->mIDCOne.mIDCOne[0];
    idcOneC = this->mIDCOne.mIDCOne[1];
  } else {
    // just set the IDCs to the same value...
    idcOneA = this->mIDCOne;
    idcOneC = idcOneA;
  }

  for (unsigned int integrationInterval = 0; integrationInterval < integrationIntervals; ++integrationInterval) {
    const unsigned int startSector = this->getFirstSectorDraw();
    const unsigned int endSector = this->getLastSectorDraw();
    const unsigned int startRegion = this->getFirstRegionDraw();
    const unsigned int endRegion = this->getLastRegionDraw();

    int padsInRegions = 0;
    for (int reg = startRegion; reg < endRegion; ++reg) {
      padsInRegions += Mapper::PADSPERREGION[reg];
    }

    const unsigned int nIDCs = padsInRegions * (endSector - startSector);
    std::vector<int> vRow(nIDCs);
    std::vector<int> vPad(nIDCs);
    std::vector<float> vXPos(nIDCs);
    std::vector<float> vYPos(nIDCs);
    std::vector<float> vGlobalXPos(nIDCs);
    std::vector<float> vGlobalYPos(nIDCs);
    std::vector<float> idcs(nIDCs);
    std::vector<float> idcsZero(nIDCs);
    std::vector<float> idcsDelta(nIDCs);
    std::vector<float> idcsDeltaMedium(nIDCs);
    std::vector<float> idcsDeltaHigh(nIDCs);
    std::vector<unsigned int> sectorv(nIDCs);

    unsigned int chunk = 0;
    unsigned int localintegrationInterval = 0;
    getLocalIntegrationInterval(0, integrationInterval, chunk, localintegrationInterval);
    IDCDelta<short> idcDeltaMedium;
    IDCDelta<char> idcDeltaHigh;

    if constexpr (std::is_same_v<Type, IDCGroupHelperSector>) {
      idcDeltaMedium = this->getIDCDeltaMediumCompressed(chunk);
      idcDeltaHigh = this->getIDCDeltaHighCompressed(chunk);
    }

    unsigned int index = 0;
    for (unsigned int sector = startSector; sector < endSector; ++sector) {
      for (unsigned int region = startRegion; region < endRegion; ++region) {
        for (unsigned int irow = 0; irow < Mapper::ROWSPERREGION[region]; ++irow) {
          for (unsigned int ipad = 0; ipad < Mapper::PADSPERROW[region][irow]; ++ipad) {
            const auto padNum = Mapper::getGlobalPadNumber(irow, ipad, region);
            const auto padTmp = (sector < SECTORSPERSIDE) ? ipad : (Mapper::PADSPERROW[region][irow] - ipad); // C-Side is mirrored
            const auto& padPosLocal = mapper.padPos(padNum);
            vRow[index] = padPosLocal.getRow();
            vPad[index] = padPosLocal.getPad();
            vXPos[index] = mapper.getPadCentre(padPosLocal).X();
            vYPos[index] = mapper.getPadCentre(padPosLocal).Y();
            const GlobalPosition2D globalPos = mapper.LocalToGlobal(LocalPosition2D(vXPos[index], vYPos[index]), Sector(sector));
            vGlobalXPos[index] = globalPos.X();
            vGlobalYPos[index] = globalPos.Y();
            idcs[index] = getIDCValUngrouped(sector, region, irow, padTmp, integrationInterval);
            idcsZero[index] = this->getIDCZeroVal(sector, region, irow, padTmp);
            idcsDelta[index] = this->getIDCDeltaVal(sector, region, irow, padTmp, chunk, localintegrationInterval);
            if constexpr (std::is_same_v<Type, IDCGroupHelperSector>) {
              idcsDeltaMedium[index] = idcDeltaMedium.getValue(Sector(sector).side(), this->getIndexUngrouped(sector, region, irow, padTmp, localintegrationInterval));
              idcsDeltaHigh[index] = idcDeltaHigh.getValue(Sector(sector).side(), this->getIndexUngrouped(sector, region, irow, padTmp, localintegrationInterval));
            }
            sectorv[index++] = sector;
          }
        }
      }
    }
    float idcOneATmp = idcOneA[integrationInterval];
    float idcOneCTmp = idcOneC[integrationInterval];
    unsigned int timeFrame = 0;
    unsigned int interval = 0;
    getTF(0, integrationInterval, timeFrame, interval);

    pcstream << "tree"
             << "integrationInterval=" << integrationInterval
             << "localintervalinTF=" << interval
             << "indexinchunk=" << localintegrationInterval
             << "chunk=" << chunk
             << "timeframe=" << timeFrame
             << "IDC.=" << idcs
             << "IDC0.=" << idcsZero
             << "IDC1A=" << idcOneATmp
             << "IDC1C=" << idcOneCTmp
             << "IDCDeltaNoComp.=" << idcsDelta
             << "IDCDeltaMediumComp.=" << idcsDeltaMedium
             << "IDCDeltaHighComp.=" << idcsDeltaHigh
             << "pad.=" << vPad
             << "row.=" << vRow
             << "lx.=" << vXPos
             << "ly.=" << vYPos
             << "gx.=" << vGlobalXPos
             << "gy.=" << vGlobalYPos
             << "sector.=" << sectorv
             << "\n";
  }
  pcstream.Close();
}

template <typename Type>
float o2::tpc::IDCFactorization<Type>::getIDCValUngrouped(const unsigned int sector, const unsigned int region, unsigned int urow, unsigned int upad, unsigned int integrationInterval) const
{
  unsigned int timeFrame = 0;
  unsigned int interval = 0;
  getTF(region, integrationInterval, timeFrame, interval);
  if (this->idcsEmpty(sector * Mapper::NREGIONS + region, timeFrame)) {
    return 0.f;
  }
  return IDCFactorizationBase<Type>::getIDCValUngrouped(sector, region, urow, upad, timeFrame, interval);
}

template <typename Type>
void o2::tpc::IDCFactorization<Type>::getTF(const unsigned int region, unsigned int integrationInterval, unsigned int& timeFrame, unsigned int& interval) const
{
  unsigned int nintervals = 0;
  unsigned int intervalTmp = 0;
  for (unsigned int tf = 0; tf < this->mTimeFrames; ++tf) {
    nintervals += this->getSizeIDCs(region, tf) / this->getIDCsPerCRU(region);
    if (integrationInterval < nintervals) {
      timeFrame = tf;
      interval = integrationInterval - intervalTmp;
      return;
    }
    intervalTmp = nintervals;
  }
}

template <typename Type>
void o2::tpc::IDCFactorization<Type>::factorizeIDCs()
{
  LOGP(info, "Using {} threads for factorization of IDCs", getNThreads());
  this->calcIDCZero();
  this->calcIDCOne();
  this->calcIDCDelta();
}

template <typename Type>
void o2::tpc::IDCFactorization<Type>::getLocalIntegrationInterval(const unsigned int region, const unsigned int integrationInterval, unsigned int& chunk, unsigned int& localintegrationInterval) const
{
  unsigned int nintervals = 0;
  unsigned int nintervalsChunk = 0;
  unsigned int globalTF = 0;
  for (unsigned int ichunk = 0; ichunk < this->getNChunks(); ++ichunk) {
    const auto nTFsPerChunk = this->getNTFsPerChunk(ichunk);
    for (unsigned int tf = 0; tf < nTFsPerChunk; ++tf) {
      nintervals += this->getSizeIDCs(region, globalTF) / this->getIDCsPerCRU(region);
      if (integrationInterval < nintervals) {
        chunk = this->getChunk(globalTF);
        localintegrationInterval = integrationInterval - nintervalsChunk;
        return;
      }
      ++globalTF;
    }
    nintervalsChunk = nintervals;
  }
}

template <typename Type>
std::vector<unsigned int> o2::tpc::IDCFactorization<Type>::getIntegrationIntervalsPerTF(const unsigned int region) const
{
  std::vector<unsigned int> integrationIntervalsPerTF;
  integrationIntervalsPerTF.reserve(this->mTimeFrames);
  for (unsigned int tf = 0; tf < this->mTimeFrames; ++tf) {
    integrationIntervalsPerTF.emplace_back(this->getSizeIDCs(region, tf) / this->getIDCsPerCRU(region));
  }
  return integrationIntervalsPerTF;
}

template <typename Type>
void o2::tpc::IDCFactorization<Type>::drawIDCDeltaHelper(const bool type, const Sector sector, const unsigned int integrationInterval, const IDCDeltaCompression compression, const std::string filename) const
{
  if constexpr (std::is_same_v<Type, IDCGroupHelperSector>) {
    std::function<float(const unsigned int, const unsigned int, const unsigned int, const unsigned int)> idcFunc;

    unsigned int chunk = 0;
    unsigned int localintegrationInterval = 0;
    getLocalIntegrationInterval(0, integrationInterval, chunk, localintegrationInterval);
    const std::string zAxisTitle = IDCDrawHelper::getZAxisTitle(IDCType::IDCDelta, compression);

    IDCDrawHelper::IDCDraw drawFun;
    switch (compression) {
      case IDCDeltaCompression::NO:
      default: {
        idcFunc = [this, chunk, localintegrationInterval](const unsigned int sector, const unsigned int region, const unsigned int irow, const unsigned int pad) {
          return this->getIDCDeltaVal(sector, region, irow, pad, chunk, localintegrationInterval);
        };
        drawFun.mIDCFunc = idcFunc;
        type ? IDCDrawHelper::drawSide(drawFun, sector.side(), zAxisTitle, filename) : IDCDrawHelper::drawSector(drawFun, this->getFirstRegionDraw(), this->getLastRegionDraw(), sector, zAxisTitle, filename);
        break;
      }
      case IDCDeltaCompression::MEDIUM: {
        const auto idcDeltaMedium = this->getIDCDeltaMediumCompressed(chunk);
        idcFunc = [this, &idcDeltaMedium, chunk, localintegrationInterval = localintegrationInterval](const unsigned int sector, const unsigned int region, const unsigned int irow, const unsigned int pad) {
          return idcDeltaMedium.getValue(Sector(sector).side(), this->getIndexUngrouped(sector, region, irow, pad, localintegrationInterval));
        };
        drawFun.mIDCFunc = idcFunc;
        type ? IDCDrawHelper::drawSide(drawFun, sector.side(), zAxisTitle, filename) : IDCDrawHelper::drawSector(drawFun, this->getFirstRegionDraw(), this->getLastRegionDraw(), sector, zAxisTitle, filename);
        break;
      }
      case IDCDeltaCompression::HIGH: {
        const auto idcDeltaHigh = this->getIDCDeltaHighCompressed(chunk);
        idcFunc = [this, &idcDeltaHigh, chunk, localintegrationInterval](const unsigned int sector, const unsigned int region, const unsigned int irow, const unsigned int pad) {
          return idcDeltaHigh.getValue(Sector(sector).side(), this->getIndexUngrouped(sector, region, irow, pad, localintegrationInterval));
        };
        drawFun.mIDCFunc = idcFunc;
        type ? IDCDrawHelper::drawSide(drawFun, sector.side(), zAxisTitle, filename) : IDCDrawHelper::drawSector(drawFun, this->getFirstRegionDraw(), this->getLastRegionDraw(), sector, zAxisTitle, filename);
        break;
      }
    }
  }
}

template <typename Type>
void o2::tpc::IDCFactorization<Type>::drawIDCZeroHelper(const bool type, const Sector sector, const std::string filename) const
{
  std::function<float(const unsigned int, const unsigned int, const unsigned int, const unsigned int)> idcFunc = [this](const unsigned int sector, const unsigned int region, const unsigned int irow, const unsigned int pad) {
    return this->getIDCZeroVal(sector, region, irow, pad);
  };

  IDCDrawHelper::IDCDraw drawFun;
  drawFun.mIDCFunc = idcFunc;
  const std::string zAxisTitle = IDCDrawHelper::getZAxisTitle(IDCType::IDCZero);
  type ? IDCDrawHelper::drawSide(drawFun, sector.side(), zAxisTitle, filename) : IDCDrawHelper::drawSector(drawFun, this->getFirstRegionDraw(), this->getLastRegionDraw(), sector, zAxisTitle, filename);
}

template <typename Type>
void o2::tpc::IDCFactorization<Type>::drawIDCHelper(const bool type, const Sector sector, const unsigned int integrationInterval, const std::string filename) const
{
  std::function<float(const unsigned int, const unsigned int, const unsigned int, const unsigned int)> idcFunc = [this, integrationInterval](const unsigned int sector, const unsigned int region, const unsigned int irow, const unsigned int pad) {
    return this->getIDCValUngrouped(sector, region, irow, pad, integrationInterval);
  };

  IDCDrawHelper::IDCDraw drawFun;
  drawFun.mIDCFunc = idcFunc;

  const std::string zAxisTitleDraw = IDCDrawHelper::getZAxisTitle(IDCType::IDC);
  type ? IDCDrawHelper::drawSide(drawFun, sector.side(), zAxisTitleDraw, filename) : IDCDrawHelper::drawSector(drawFun, this->getFirstRegionDraw(), this->getLastRegionDraw(), sector, zAxisTitleDraw, filename);
}

template class o2::tpc::IDCFactorization<o2::tpc::IDCGroupHelperSector>;
template class o2::tpc::IDCFactorization<o2::tpc::IDCGroupHelperRegion>;
