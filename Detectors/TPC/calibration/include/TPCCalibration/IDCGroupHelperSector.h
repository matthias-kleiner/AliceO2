// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file IDCGroupHelperSector.h
/// \brief class for storing grouped IDCs
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_TPC_IDCGROUPHELPERSECTOR_H_
#define ALICEO2_TPC_IDCGROUPHELPERSECTOR_H_

#include <vector>
#include <numeric>
#include "Rtypes.h"
#include "TPCBase/Mapper.h"
#include "TPCCalibration/IDCGroupHelperRegion.h"
#include "TPCCalibration/IDCFactorizationContainer.h"

namespace o2::tpc
{

/// Class to hold grouped IDC values for one CRU for one TF

class IDCGroupHelperSector
{
 public:
  /// constructor
  /// \param groupPads number of pads in pad direction which will be grouped
  /// \param groupRows number of pads in row direction which will be grouped
  /// \param groupLastRowsThreshold minimum number of pads in row direction for the last group in row direction
  /// \param groupLastPadsThreshold minimum number of pads in pad direction for the last group in pad direction
  IDCGroupHelperSector(const std::array<unsigned char, Mapper::NREGIONS>& groupPads, const std::array<unsigned char, Mapper::NREGIONS>& groupRows, const std::array<unsigned char, Mapper::NREGIONS>& groupLastRowsThreshold, const std::array<unsigned char, Mapper::NREGIONS>& groupLastPadsThreshold)
    : mGroupPads{groupPads}, mGroupRows{groupRows}, mGroupLastRowsThreshold{groupLastRowsThreshold}, mGroupLastPadsThreshold{groupLastPadsThreshold}
  {
    for (unsigned int reg = 0; reg < Mapper::NREGIONS; ++reg) {
      const IDCGroupHelperRegion groupTmp(mGroupPads[reg], mGroupRows[reg], mGroupLastRowsThreshold[reg], mGroupLastPadsThreshold[reg], reg);
      mNIDCsPerCRU[reg] = groupTmp.getNIDCsPerIntegrationInterval();
      mRows[reg] = groupTmp.getNRows();
      mPadsPerRow[reg] = groupTmp.getPadsPerRow();
      mOffsRow[reg] = groupTmp.getRowOffset();
      if (reg > 0) {
        const unsigned int lastInd = reg - 1;
        mRegionOffs[reg] = mRegionOffs[lastInd] + mNIDCsPerCRU[lastInd];
      }
    }
    mNIDCsPerSector = static_cast<unsigned int>(std::accumulate(mNIDCsPerCRU.begin(), mNIDCsPerCRU.end(), 0));
  };

  /// default constructor for ROOT I/O
  IDCGroupHelperSector() = default;

  /// \return returns index to the data
  /// \param row row of the grouped IDCs
  /// \param pad pad of the grouped IDCs
  unsigned int getIndexGrouped(const unsigned int sector, const unsigned int region, const unsigned int row, const unsigned int pad, unsigned int integrationInterval) const { return mNIDCsPerSector * (integrationInterval * SECTORSPERSIDE + sector) + mRegionOffs[region] + mOffsRow[region][row] + pad; }

  /// \return returns offsets for rows to calculate data index
    /// \param grow grouped row in the region of the grouped IDCs
  const unsigned int getRowOffset(const unsigned int region, const unsigned int grow) const { return mOffsRow[region][grow]; }

  /// \return returns grouped pad for ungrouped row and pad
  /// \param region region
  /// \param urow local ungrouped row in a region
  /// \param upad ungrouped pad
  unsigned int getGroupedPad(const unsigned int region, unsigned int urow, unsigned int upad) const { return IDCGroupHelperRegion::getGroupedPad(upad, urow, region, mGroupPads[region], mGroupRows[region], mRows[region], mPadsPerRow[region]); }

  /// \return returns the row of the group from the local ungrouped row in a region
  /// \param region region
  /// \param urow local ungrouped row in a region
  unsigned int getGroupedRow(const unsigned int region, unsigned int urow) const { return IDCGroupHelperRegion::getGroupedRow(urow, mGroupRows[region], mRows[region]); }

  /// \return returns the index to the grouped data with ungrouped inputs
  /// \param sector sector
  /// \param region TPC region
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  unsigned int getIndexUngrouped(const unsigned int sector, const unsigned int region, unsigned int urow, unsigned int upad, unsigned int integrationInterval) const { return getIndexGrouped(sector % o2::tpc::SECTORSPERSIDE, region, getGroupedRow(region, urow), getGroupedPad(region, urow, upad), integrationInterval); }

  /// \returns grouping definition in pad direction (How many pads are grouped)
  /// \param region TPC region
  unsigned int getGroupPads(const unsigned int region) const { return mGroupPads[region]; }

  /// \returns grouping definition in pad direction (How many pads are grouped)
  const auto& getGroupPads() const { return mGroupPads; }

  /// \returns grouping definition in row direction (How many rows are grouped)
  /// \param region TPC region
  unsigned int getGroupRows(const unsigned int region) const { return mGroupRows[region]; }

  /// \returns grouping definition in row direction (How many rows are grouped)
  const auto& getGroupRows() const { return mGroupRows; }

  /// \returns grouping threshold for last group in pad direction
  /// \param region TPC region
  unsigned int getPadThreshold(const unsigned int region) const { return mGroupLastPadsThreshold[region]; }

  /// \returns grouping threshold for last group in pad direction
  const auto& getPadThreshold() const { return mGroupLastPadsThreshold; }

  /// \returns grouping threshold for last group in row direction
  /// \param region TPC region
  unsigned int getRowThreshold(const unsigned int region) const { return mGroupLastRowsThreshold[region]; }

  /// \returns grouping threshold for last group in row direction
  const auto& getRowThreshold() const { return mGroupLastRowsThreshold; }

  /// \return returns total number of IDCs per region per integration interval
  /// \param region TPC region
  unsigned int getIDCsPerCRU(const unsigned int region) const { return mNIDCsPerCRU[region]; }

  /// draw IDC zero I_0(r,\phi) = <I(r,\phi,t)>_t
  /// \param side side which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  // void drawIDCZeroSide(const o2::tpc::Side side, const std::string filename = "IDCZeroSide.pdf") const { drawSide(IDCType::IDCZero, side, 0, filename); }
  // void drawSide(IDCDelta<float> idcdelta);

protected:
  const std::array<unsigned char, Mapper::NREGIONS> mGroupPads{};              ///< grouping definition in pad direction (How many pads are grouped)
  const std::array<unsigned char, Mapper::NREGIONS> mGroupRows{};              ///< grouping definition in row direction (How many rows are grouped)
  const std::array<unsigned char, Mapper::NREGIONS> mGroupLastRowsThreshold{}; ///< if the last group (region edges) consists in row direction of less then mGroupLastRowsThreshold pads then it will be grouped into the previous group
  const std::array<unsigned char, Mapper::NREGIONS> mGroupLastPadsThreshold{}; ///< if the last group (sector edges) consists in pad direction of less then mGroupLastPadsThreshold pads then it will be grouped into the previous group
  std::array<unsigned int, Mapper::NREGIONS> mNIDCsPerCRU{1};                 ///< total number of IDCs per region per integration interval
  unsigned int mNIDCsPerSector{};                                             ///< number of grouped IDCs per sector
  std::array<unsigned int, Mapper::NREGIONS> mRows{};                         ///< number of grouped rows per region
  std::array<unsigned int, Mapper::NREGIONS> mRegionOffs{};                   ///< offset for the region per region
  std::array<std::vector<unsigned int>, Mapper::NREGIONS> mPadsPerRow{};      ///< number of pads per row per region
  std::array<std::vector<unsigned int>, Mapper::NREGIONS> mOffsRow{};         ///< offset to calculate the index in the data from row and pad per region

  ClassDefNV(IDCGroupHelperSector, 1)
};

} // namespace o2::tpc

#endif
