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

#include "TPCCalibration/IDCCCDBHelper.h"
#include "TPCCalibration/IDCDrawHelper.h"
#include "TPCCalibration/IDCGroupHelperSector.h"
#include "TPCBase/Mapper.h"

template <typename DataT>
float o2::tpc::IDCCCDBHelper<DataT>::getIDCZeroVal(const unsigned int sector, const unsigned int region, unsigned int urow, unsigned int upad) const
{
  return mIDCZero->getValueIDCZero(Sector(sector).side(), mHelperSector->getIndexUngrouped(sector, region, urow, upad, 0));
}

template <typename DataT>
float o2::tpc::IDCCCDBHelper<DataT>::getIDCDeltaVal(const unsigned int sector, const unsigned int region, unsigned int urow, unsigned int upad, unsigned int localintegrationInterval) const
{
  return mIDCDelta->getValue(Sector(sector).side(), mHelperSector->getIndexUngrouped(sector, region, urow, upad, localintegrationInterval));
}

template <typename DataT>
void o2::tpc::IDCCCDBHelper<DataT>::loadGroupingParameter()
{
  mHelperSector = std::make_unique<IDCGroupHelperSector>(IDCGroupHelperSector{*mCCDBManager.get<o2::tpc::ParameterIDCGroupCCDB>("TPC/Calib/IDC/GROUPINGPAR")});
}

template <typename DataT>
void o2::tpc::IDCCCDBHelper<DataT>::drawIDCZeroHelper(const bool type, const o2::tpc::Sector sector, const std::string filename) const
{
  std::function<float(const unsigned int, const unsigned int, const unsigned int, const unsigned int)> idcFunc = [this](const unsigned int sector, const unsigned int region, const unsigned int irow, const unsigned int pad) {
    return this->getIDCZeroVal(sector, region, irow, pad);
  };

  IDCDrawHelper::IDCDraw drawFun;
  drawFun.mIDCFunc = idcFunc;
  const std::string zAxisTitle = IDCDrawHelper::getZAxisTitle(IDCType::IDCZero);
  type ? IDCDrawHelper::drawSide(drawFun, sector.side(), zAxisTitle, filename) : IDCDrawHelper::drawSector(drawFun, 0, Mapper::NREGIONS, sector, zAxisTitle, filename);
}

template <typename DataT>
void o2::tpc::IDCCCDBHelper<DataT>::drawIDCDeltaHelper(const bool type, const Sector sector, const unsigned int integrationInterval, const std::string filename) const
{
  std::function<float(const unsigned int, const unsigned int, const unsigned int, const unsigned int)> idcFunc = [this, integrationInterval](const unsigned int sector, const unsigned int region, const unsigned int irow, const unsigned int pad) {
    return this->getIDCDeltaVal(sector, region, irow, pad, integrationInterval);
  };

  IDCDrawHelper::IDCDraw drawFun;
  drawFun.mIDCFunc = idcFunc;
  const std::string zAxisTitle = IDCDrawHelper::getZAxisTitle(IDCType::IDCDelta);
  type ? IDCDrawHelper::drawSide(drawFun, sector.side(), zAxisTitle, filename) : IDCDrawHelper::drawSector(drawFun, 0, Mapper::NREGIONS, sector, zAxisTitle, filename);
}

/// load IDC-Delta, 0D-IDCs, grouping parameter
template <typename DataT>
void o2::tpc::IDCCCDBHelper<DataT>::loadAll()
{
  loadIDCDelta();
  loadIDCZero();
  loadGroupingParameter();
}

template class o2::tpc::IDCCCDBHelper<float>;
template class o2::tpc::IDCCCDBHelper<short>;
template class o2::tpc::IDCCCDBHelper<char>;
