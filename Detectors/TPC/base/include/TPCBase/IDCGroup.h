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
/// \brief class for averaging/merging and storing to IDCs
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_TPC_IDCGROUP_H_
#define ALICEO2_TPC_IDCGROUP_H_

#include <vector>
#include <numeric>
#include "Rtypes.h"

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
  IDCGroup() = default;

  void setRows(const uint32_t nRows)
  {
    mRows = nRows;
    mPadsPerRow.resize(mRows);
    mOffsRow.resize(mRows);
  }

  /// \param row row of the grouped IDCs
  /// \param nPads number of pads of the grouped IDCs for given row
  void setPadsPerRow(const uint32_t row, const uint32_t nPads) { mPadsPerRow[row] = nPads; }

  /// \return returns the stored value
  /// \param row row of the grouped IDCs
  /// \param pad pad number of the grouped IDCs
  const float& operator()(uint32_t row, uint32_t pad) const { return mIDCsGrouped[getIndex(row, pad)]; }

  /// \return returns the stored value
  /// \param row row of the grouped IDCs
  /// \param pad pad number of the grouped IDCs
  float& operator()(uint32_t row, uint32_t pad) { return mIDCsGrouped[getIndex(row, pad)]; }

  /// \return returns index to the data
  /// \param row row of the grouped IDCs
  /// \param pad pad of the grouped IDCs
  uint32_t getIndex(const uint32_t row, const uint32_t pad) const { return mOffsRow[row] + pad; }

  uint32_t getNRows() const { return mRows; }

  /// \param row row of the grouped IDCs
  uint32_t getPadsPerRow(const uint32_t row) { return mPadsPerRow[row]; }

  /// initialize the member containing the grouped IDCs
  void initStorage();

  const auto& getData() const { return mIDCsGrouped; }
  auto& getData() { return mIDCsGrouped; }

  /// dump the IDCs to a tree
  /// \param outname name of the output file
  void dumpToTree(const char* outname = "IDCGroup.root") const;

 private:
  uint32_t mRows{};                    ///< number of rows
  std::vector<uint32_t> mPadsPerRow{}; ///< number of pads per row
  std::vector<uint32_t> mOffsRow{};    ///< offset to calculate the index in the data from row and pad
  std::vector<float> mIDCsGrouped{};   ///< grouped and averaged IDC values

  ClassDefNV(IDCGroup, 1)
};

} // namespace tpc
} // namespace o2

#endif
