// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file IDCGroup.h
/// \brief class for storing grouped IDCs
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_TPC_IDCGROUP_H_
#define ALICEO2_TPC_IDCGROUP_H_

#include <vector>
#include <numeric>
#include "Rtypes.h"
#include "TPCBase/IDCHelper.h"

namespace o2
{
namespace tpc
{

/// Class to hold grouped and averaged IDC values for one CRU for one TF
/// Usage:
/// 1. set the number of rows setRows()
/// 2. set the number of pads for each row setPadsPerRow()
/// 3. initialize the storage and calculate offsets for data access initStorage()

class IDCGroup
{
 public:
  IDCGroup(const unsigned int groupPads = 4, const unsigned int groupRows = 4, const unsigned int groupLastRowsThreshold = 2, const unsigned int groupLastPadsThreshold = 2, const unsigned int region = 0) : mGroupPads{groupPads}, mGroupRows{groupRows}, mGroupLastRowsThreshold{groupLastRowsThreshold}, mGroupLastPadsThreshold{groupLastPadsThreshold}, mRegion{region}
  {
    initIDCGroup();
  }

  void setRows(const unsigned int nRows)
  {
    mRows = nRows;
    mPadsPerRow.resize(mRows);
    mOffsRow.resize(mRows);
  }

  /// \param row row of the grouped IDCs
  /// \param nPads number of pads of the grouped IDCs for given row
  void setPadsPerRow(const unsigned int row, const unsigned int nPads) { mPadsPerRow[row] = nPads; }

  /// \return returns the stored value
  /// \param row row of the grouped IDCs
  /// \param pad pad number of the grouped IDCs
  /// \param integrationInterval integration interval
  const float& operator()(unsigned int row, unsigned int pad, unsigned int integrationInterval) const { return mIDCsGrouped[getIndex(row, pad, integrationInterval)]; }

  /// \return returns the stored value
  /// \param row row of the grouped IDCs
  /// \param pad pad number of the grouped IDCs
  /// \param integrationInterval integration interval
  float& operator()(unsigned int row, unsigned int pad, unsigned int integrationInterval) { return mIDCsGrouped[getIndex(row, pad, integrationInterval)]; }

  /// \return returns the global pad number for given local pad row and pad
  /// \param lrow local row in a region
  /// \param lrow local row in a region
  unsigned int getGlobalPadNumber(const unsigned int lrow, const unsigned int pad) const { return IDCHelper::GLOBALPADOFFSET[mRegion] + IDCHelper::OFFSETCRULOCAL[mRegion][lrow] + pad; }

  /// \return returns index to the data
  /// \param row row of the grouped IDCs
  /// \param pad pad of the grouped IDCs
  unsigned int getIndex(const unsigned int row, const unsigned int pad, unsigned int integrationInterval) const { return mNIDCsPerCRU * integrationInterval + mOffsRow[row] + pad; }

  unsigned int getNRows() const { return mRows; }

  /// \param row row of the grouped IDCs
  unsigned int getPadsPerRow(const unsigned int row) const { return mPadsPerRow[row]; }

  /// initialize the member containing the grouped IDCs
  void initStorage();

  const auto& getData() const { return mIDCsGrouped; }
  auto& getData() { return mIDCsGrouped; }

  /// extend the size of the grouped and averaged IDC values corresponding to the number of integration intervals.
  /// without using this function the object can hold only one integration interval
  /// \param nIntegrationIntervals number of ontegration intervals for which teh IDCs are stored
  void resize(const unsigned int nIntegrationIntervals) { mIDCsGrouped.resize(mNIDCsPerCRU * nIntegrationIntervals); }

  /// dump the IDCs to a tree
  /// \param outname name of the output file
  void dumpToTree(const char* outname = "IDCGroup.root") const;

  int getNIntegrationIntervals() const { return mIDCsGrouped.size() / mNIDCsPerCRU; }

  /// \return returns the number pads in pad direction which are grouped
  unsigned int getGroupPads() const { return mGroupPads; }

  /// \return returns the number pads in row direction which are grouped
  unsigned int getGroupRows() const { return mGroupRows; }

  /// \return returns threshold for grouping the last group in row direction
  unsigned int getGroupLastRowsThreshold() const { return mGroupLastRowsThreshold; }

  /// \return returns threshold for grouping the last group in pad direction
  unsigned int getGroupLastPadsThreshold() const { return mGroupLastPadsThreshold; }

  /// \return returns the row of the group from the local row in a region
  /// \param lrow local row in a region
  unsigned int getGroupedRow(const unsigned int lrow) const;

  /// \return returns the region for which the IDCs are stored
  unsigned int getRegion() const { return mRegion; }

  /// \return returns the grouped pad index from ungrouped pad and row
  /// \param pad ungrouped pad
  /// \param lrow local ungrouped row in a region
  int getGroupedPad(const unsigned int pad, const unsigned int lrow) const;

  /// \return returns last ungrouped row index for mRegion
  int getLastRow() const;

  /// \return returns last ungrouped pad index for given flobal row
  int getLastPad(const int row) const;

  /// draw grouped IDCs
  /// \param integrationInterval integration interval for which the IDCs will be drawn
  void draw(const int integrationInterval = 0) const;

  void initIDCGroup();

 private:
  const unsigned int mGroupPads{4};              ///< group 4 pads
  const unsigned int mGroupRows{4};              ///< group 4 pads -> 4x4
  const unsigned int mGroupLastRowsThreshold{2}; ///< if the last group (region edges) consists in row direction less then mGroupLastRowsThreshold pads then it will be grouped into the previous group
  const unsigned int mGroupLastPadsThreshold{2}; ///< if the last group (sector edges) consists in pad direction less then mGroupLastPadsThreshold pads then it will be grouped into the previous group
  const unsigned int mRegion{};                  ///< region of input IDCs
  unsigned int mNIDCsPerCRU{1};                  ///< total number of IDCs per CRU per intgeration interval
  unsigned int mRows{};                          ///< number of rows
  std::vector<unsigned int> mPadsPerRow{};       ///< number of pads per row
  std::vector<unsigned int> mOffsRow{};          ///< offset to calculate the index in the data from row and pad
  std::vector<float> mIDCsGrouped{};             ///< grouped and averaged IDC values

  ClassDefNV(IDCGroup, 1)
};

} // namespace tpc
} // namespace o2

#endif
