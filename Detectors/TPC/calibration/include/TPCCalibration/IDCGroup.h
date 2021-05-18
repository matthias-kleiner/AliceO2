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
#include "TPCBase/Mapper.h"

namespace o2::tpc
{

/// Class to hold grouped IDC values for one CRU for one TF

class IDCGroup
{
 public:
  /// constructor
  /// \param groupPads number of pads in pad direction which will be grouped
  /// \param groupRows number of pads in row direction which will be grouped
  /// \param groupLastRowsThreshold minimum number of pads in row direction for the last group in row direction
  /// \param groupLastPadsThreshold minimum number of pads in pad direction for the last group in pad direction
  /// \param region region of the TPC
  IDCGroup(const unsigned int groupPads = 4, const unsigned int groupRows = 4, const unsigned int groupLastRowsThreshold = 2, const unsigned int groupLastPadsThreshold = 2, const unsigned int region = 0)
    : mGroupPads{groupPads}, mGroupRows{groupRows}, mGroupLastRowsThreshold{groupLastRowsThreshold}, mGroupLastPadsThreshold{groupLastPadsThreshold}, mRegion{region}
  {
    initIDCGroup();
  }

  /// extend the size of the grouped and averaged IDC values corresponding to the number of integration intervals.
  /// without using this function the object can hold only one integration interval
  /// \param nIntegrationIntervals number of ontegration intervals for which teh IDCs are stored
  void resize(const unsigned int nIntegrationIntervals) { mIDCsGrouped.resize(mNIDCsPerCRU * nIntegrationIntervals); }

  /// \return returns the stored value
  /// \param row row of the grouped IDCs
  /// \param pad pad number of the grouped IDCs
  /// \param integrationInterval integration interval
  float operator()(unsigned int row, unsigned int pad, unsigned int integrationInterval) const { return mIDCsGrouped[getIndex(row, pad, integrationInterval)]; }

  /// \return returns the stored value
  /// \param row row of the grouped IDCs
  /// \param pad pad number of the grouped IDCs
  /// \param integrationInterval integration interval
  float& operator()(unsigned int row, unsigned int pad, unsigned int integrationInterval) { return mIDCsGrouped[getIndex(row, pad, integrationInterval)]; }

  /// \return returns the stored value for local ungrouped pad row and ungrouped pad
  /// \param urow local row in region of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  float& setValUngrouped(unsigned int urow, unsigned int upad, unsigned int integrationInterval) { return mIDCsGrouped[getIndexUngrouped(urow, upad, integrationInterval)]; }

  /// \return returns the stored value for local ungrouped pad row and ungrouped pad
  /// \param urow local row in region of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  float getValUngrouped(unsigned int urow, unsigned int upad, unsigned int integrationInterval) const { return mIDCsGrouped[getIndexUngrouped(urow, upad, integrationInterval)]; }

  /// \return returns the stored value for local ungrouped pad row and ungrouped pad
  /// \param grow global row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  float getValUngroupedGlobal(unsigned int grow, unsigned int upad, unsigned int integrationInterval) const { return mIDCsGrouped[getIndexUngrouped(Mapper::getLocalRowFromGlobalRow(grow), upad, integrationInterval)]; }

  /// \return returns number of grouped rows
  unsigned int getNRows() const { return mRows; }

  /// \return returns number of grouped pads
  /// \param row grouped row
  unsigned int getPadsPerRow(const unsigned int row) const { return mPadsPerRow[row]; }

  /// \return returns number of grouped pads for all rows
  const std::vector<unsigned int>& getPadsPerRow() const { return mPadsPerRow; }

  /// \return returns offsets for rows to calculate data index
  const std::vector<unsigned int>& getRowOffset() const { return mOffsRow; }

  /// \return returns grouped and averaged IDC values
  const auto& getData() const { return mIDCsGrouped; }
  auto& getData() { return mIDCsGrouped; }

  /// \return returns number of stored integration intervals
  unsigned int getNIntegrationIntervals() const { return mIDCsGrouped.size() / mNIDCsPerCRU; }

  /// \return returns the number of pads in pad direction which are grouped
  unsigned int getGroupPads() const { return mGroupPads; }

  /// \return returns the number of pads in row direction which are grouped
  unsigned int getGroupRows() const { return mGroupRows; }

  /// \return returns threshold for grouping the last group in row direction
  unsigned int getGroupLastRowsThreshold() const { return mGroupLastRowsThreshold; }

  /// \return returns threshold for grouping the last group in pad direction
  unsigned int getGroupLastPadsThreshold() const { return mGroupLastPadsThreshold; }

  /// \return returns the region for which the IDCs are stored
  unsigned int getRegion() const { return mRegion; }

  /// \returns returns number of IDCS per integration interval
  unsigned int getNIDCsPerIntegrationInterval() const { return mNIDCsPerCRU; }

  /// \return returns the row of the group from the local ungrouped row in a region
  /// \param lrow local ungrouped row in a region
  /// \param groupRows grouping parameter for number of pads in row direction which are grouped
  /// \param groupedrows number of grouped rows
  static unsigned int getGroupedRow(const unsigned int lrow, const unsigned int groupRows, const unsigned int groupedrows);

  /// \return returns the row of the group from the local ungrouped row in a region
  /// \param lrow local ungrouped row in a region
  unsigned int getGroupedRow(const unsigned int lrow) const { return getGroupedRow(lrow, mGroupRows, mRows); }

  /// \return returns the grouped pad index from ungrouped pad and row
  /// \param pad ungrouped pad
  /// \param lrow local ungrouped row in a region
  /// \param region region
  /// \param groupPads grouping parameter for number of pads in pad direction which are grouped
  /// \param groupRows grouping parameter for number of pads in row direction which are grouped
  /// \param groupedrows number of grouped rows
  /// \param padsPerRow vector containing the number of pads per row
  static unsigned int getGroupedPad(const unsigned int pad, const unsigned int lrow, const unsigned int region, const unsigned int groupPads, const unsigned int groupRows, const unsigned int groupedrows, const std::vector<unsigned int>& padsPerRow);

  /// \return returns the grouped pad index from ungrouped pad and row
  /// \param pad ungrouped pad
  /// \param lrow local ungrouped row in a region
  unsigned int getGroupedPad(const unsigned int pad, const unsigned int lrow) const { return getGroupedPad(pad, lrow, mRegion, mGroupPads, mGroupRows, mRows, mPadsPerRow); };

  /// \return returns last ungrouped row
  unsigned int getLastRow() const;

  /// \return returns last ungrouped pad for given global row
  /// \param row ungrouped row
  unsigned int getLastPad(const unsigned int row) const;

  /// dump the IDCs to a tree
  /// \param outname name of the output file
  void dumpToTree(const char* outname = "IDCGroup.root") const;

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName = "IDCGroup.root", const char* outName = "IDCGroup") const;

  /// draw grouped IDCs
  /// \param integrationInterval integration interval for which the IDCs will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void draw(const unsigned int integrationInterval = 0, const std::string filename = "IDCsGrouped.pdf") const;

 private:
  const unsigned int mGroupPads{4};              ///< grouping parameter in pad direction (how many pads in pad direction are grouped)
  const unsigned int mGroupRows{4};              ///< grouping parameter in pad direction (how many pads in pad direction are grouped)
  const unsigned int mGroupLastRowsThreshold{2}; ///< if the last group (region edges) consists in row direction less then mGroupLastRowsThreshold pads then it will be grouped into the previous group
  const unsigned int mGroupLastPadsThreshold{2}; ///< if the last group (sector edges) consists in pad direction less then mGroupLastPadsThreshold pads then it will be grouped into the previous group
  const unsigned int mRegion{};                  ///< region of input IDCs
  unsigned int mNIDCsPerCRU{1};                  ///< total number of IDCs per CRU per integration interval
  unsigned int mRows{};                          ///< number of grouped rows
  std::vector<unsigned int> mPadsPerRow{};       ///< number of pads per row
  std::vector<unsigned int> mOffsRow{};          ///< offset to calculate the index in the data from row and pad
  std::vector<float> mIDCsGrouped{};             ///< grouped and averaged IDC values

  /// set number of grouped rows
  void setRows(const unsigned int nRows);

  /// initialize members
  void initIDCGroup();

  /// initialize the member containing the grouped IDCs
  void initStorage();

  /// \return returns index to the data
  /// \param row row of the grouped IDCs
  /// \param pad pad of the grouped IDCs
  /// \param integrationInterval integration interval
  unsigned int getIndex(const unsigned int row, const unsigned int pad, unsigned int integrationInterval) const { return mNIDCsPerCRU * integrationInterval + mOffsRow[row] + pad; }

  /// \return returns index to the data
  /// \param urow ungrouped row
  /// \param upad ungrouped pad
  /// \param integrationInterval integration interval
  unsigned int getIndexUngrouped(const unsigned int urow, const unsigned int upad, unsigned int integrationInterval) const { return getIndex(getGroupedRow(urow), getGroupedPad(upad, urow), integrationInterval); }

  /// \return returns the global pad number for given local pad row and pad
  /// \param lrow local ungrouped row in a region
  /// \param pad ungrouped pad in row
  unsigned int getGlobalPadNumber(const unsigned int lrow, const unsigned int pad) const { return Mapper::getGlobalPadNumber(lrow, pad, mRegion); }

  ClassDefNV(IDCGroup, 1)
};

} // namespace o2::tpc

#endif
