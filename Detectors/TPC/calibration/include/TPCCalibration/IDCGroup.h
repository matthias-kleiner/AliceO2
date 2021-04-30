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
#include "TPCCalibration/ParameterIDCGroup.h"

namespace o2
{
namespace tpc
{

/// Class to hold grouped and averaged IDC values for one CRU for one TF
/// Usage: TODO
/// 1. set the number of integration intervals after creation the object if the number of integration intervals which will be stored are larger than 1 by using resize(nIntegrationIntervals)

class IDCGroup
{
 public:
  IDCGroup(const unsigned int groupPads = 4, const unsigned int groupRows = 4, const unsigned int groupLastRowsThreshold = 2, const unsigned int groupLastPadsThreshold = 2, const unsigned int region = 0)
    : mGroupPads{groupPads}, mGroupRows{groupRows}, mGroupLastRowsThreshold{groupLastRowsThreshold}, mGroupLastPadsThreshold{groupLastPadsThreshold}, mRegion{region}
  {
    initIDCGroup();
  }

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

  /// \return returns the stored value for local ungrouped pad row and ungrouped pad
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  float& getVal(unsigned int urow, unsigned int upad, unsigned int integrationInterval) { return mIDCsGrouped[getIndex(getGroupedRow(urow), getGroupedPad(upad, urow), integrationInterval)]; }

  /// \return returns the stored value for local ungrouped pad row and ungrouped pad
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  const float& getVal(unsigned int urow, unsigned int upad, unsigned int integrationInterval) const { return mIDCsGrouped[getIndex(getGroupedRow(urow), getGroupedPad(upad, urow), integrationInterval)]; }

  /// \return returns the global pad number for given local pad row and pad
  /// \param lrow local row in a region
  /// \param pad pad in row
  unsigned int getGlobalPadNumber(const unsigned int lrow, const unsigned int pad) const { return Mapper::GLOBALPADOFFSET[mRegion] + Mapper::OFFSETCRULOCAL[mRegion][lrow] + pad; }

  /// \return returns the global pad number for given local pad row and pad
  /// \param lrow local row in a region
  /// \param pad pad in row
  unsigned int static getGlobalPadNumber(const unsigned int lrow, const unsigned int pad, const unsigned int region) { return Mapper::GLOBALPADOFFSET[region] + Mapper::OFFSETCRULOCAL[region][lrow] + pad; }

  unsigned int getNRows() const { return mRows; }

  /// \param row row of the grouped IDCs
  unsigned int getPadsPerRow(const unsigned int row) const { return mPadsPerRow[row]; }
  const std::vector<unsigned int>& getPadsPerRow() const { return mPadsPerRow; }

  const std::vector<unsigned int>& getRowOffset() const { return mOffsRow; }

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

  /// \return returns the region for which the IDCs are stored
  unsigned int getRegion() const { return mRegion; }

  /// \return returns the row of the group from the local row in a region
  /// \param lrow local ungrouped row in a region
  static unsigned int getGroupedRow(const unsigned int lrow, const unsigned int groupRows, const unsigned int rows);

  unsigned int getGroupedRow(const unsigned int lrow) const { return getGroupedRow(lrow, mGroupRows, mRows); }

  /// \return returns the grouped pad index from ungrouped pad and row
  /// \param pad ungrouped pad
  /// \param lrow local ungrouped row in a region
  static int getGroupedPad(const unsigned int pad, const unsigned int lrow, const unsigned int region, const unsigned int groupPads, const unsigned int groupRows, const unsigned int rows, const std::vector<unsigned int>& padsPerRow);

  /// \return returns the grouped pad index from ungrouped pad and row
  /// \param pad ungrouped pad
  /// \param lrow local ungrouped row in a region
  int getGroupedPad(const unsigned int pad, const unsigned int lrow) const { return getGroupedPad(pad, lrow, mRegion, mGroupPads, mGroupRows, mRows, mPadsPerRow); };

  /// \return returns last ungrouped row index for mRegion
  int getLastRow() const;

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName = "IDCGroup.root", const char* outName = "IDCGroup") const;

  /// \return returns last ungrouped pad index for given global row
  int getLastPad(const int row) const;

  /// draw grouped IDCs
  /// \param integrationInterval integration interval for which the IDCs will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void draw(const unsigned int integrationInterval = 0, const std::string filename = "IDCsGrouped.pdf") const;

  unsigned int getNIDCsPerCRU() const { return mNIDCsPerCRU; }

 private:
  const unsigned int mGroupPads{4};              ///< group 4 pads
  const unsigned int mGroupRows{4};              ///< group 4 pads -> 4x4
  const unsigned int mGroupLastRowsThreshold{2}; ///< if the last group (region edges) consists in row direction less then mGroupLastRowsThreshold pads then it will be grouped into the previous group
  const unsigned int mGroupLastPadsThreshold{2}; ///< if the last group (sector edges) consists in pad direction less then mGroupLastPadsThreshold pads then it will be grouped into the previous group
  const unsigned int mRegion{};                  ///< region of input IDCs
  unsigned int mNIDCsPerCRU{1};                  ///< total number of IDCs per CRU per integration interval
  unsigned int mRows{};                          ///< number of grouped rows
  std::vector<unsigned int> mPadsPerRow{};       ///< number of pads per row
  std::vector<unsigned int> mOffsRow{};          ///< offset to calculate the index in the data from row and pad
  std::vector<float> mIDCsGrouped{};             ///< grouped and averaged IDC values

  void setRows(const unsigned int nRows);

  void initIDCGroup();

  /// initialize the member containing the grouped IDCs
  void initStorage();

  /// \return returns index to the data
  /// \param row row of the grouped IDCs
  /// \param pad pad of the grouped IDCs
  unsigned int getIndex(const unsigned int row, const unsigned int pad, unsigned int integrationInterval) const { return mNIDCsPerCRU * integrationInterval + mOffsRow[row] + pad; }

  ClassDefNV(IDCGroup, 1)
};

} // namespace tpc
} // namespace o2

#endif
