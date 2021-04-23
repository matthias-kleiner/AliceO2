// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file IDCAverageGroup.h
/// \brief class for averaging/grouping and storing to IDCs
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_IDCAVERAGEGROUP_H_
#define ALICEO2_IDCAVERAGEGROUP_H_

#include <vector>
#include "TPCBase/IDCGroup.h"
#include "TPCBase/IDCHelper.h"

#include "Framework/Logger.h"
namespace o2
{
namespace tpc
{
class IDCAverageGroup
{

 public:
  IDCAverageGroup(const unsigned int groupPads = 4, const unsigned int groupRows = 4, const unsigned int groupLastRowsThreshold = 2, const unsigned int groupLastPadsThreshold = 2, const unsigned int cru = 0)
    : mGroupPads{groupPads}, mGroupRows{groupRows}, mGroupLastRowsThreshold{groupLastRowsThreshold}, mGroupLastPadsThreshold{groupLastPadsThreshold}, mCRU{cru}
  {
    initIDCGroup();
  }

  /// \return returns the row of the group from the local row in a region
  /// \param lrow local row in a region
  unsigned int getGroupedRow(const unsigned int lrow) const;

  /// \return returns the grouped pad index from ungrouped pad and row
  /// \param pad ungrouped pad
  /// \param lrow local ungrouped row in a region
  int getGroupedPad(const unsigned int pad, const unsigned int lrow) const;

  /// \return returns the global pad number for given local pad row and pad
  /// \param lrow local row in a region
  /// \param lrow local row in a region
  unsigned int getGlobalPadNumber(const unsigned int lrow, const unsigned int pad) const { return IDCHelper::GLOBALPADOFFSET[mRegion] + IDCHelper::OFFSETCRULOCAL[mRegion][lrow] + pad; }

  /// \param IDCs vector containg the IDCs
  void setIDCs(const std::vector<float>& idcs)
  {
    mIDCs = idcs;
    mIDCsGrouped.resize(getNIntegrationIntervals());
  }

  /// \param IDCs vector containg the IDCs
  void setIDCs(std::vector<float>&& idcs)
  {
    mIDCs = std::move(idcs);
    mIDCsGrouped.resize(getNIntegrationIntervals());
  }

  /// \return returns number of integration intervalls stored in this object
  unsigned int getNIntegrationIntervals() const { return mIDCs.size() / IDCHelper::PADSPERREGION[mRegion]; }

  /// grouping and averaging of IDCs
  /// \param debug
  void processIDCs();

  /// \return returns grouped IDC object
  const auto& getIDCGroup() const { return mIDCsGrouped; }

  /// \return returns the number pads in pad direction which are grouped
  unsigned int getGroupPads() const { return mGroupPads; }

  /// \return returns the number pads in row direction which are grouped
  unsigned int getGroupRows() const { return mGroupRows; }

  /// \return returns threshold for grouping the last group in row direction
  unsigned int getGroupLastRowsThreshold() const { return mGroupLastRowsThreshold; }

  /// \return returns threshold for grouping the last group in pad direction
  unsigned int getGroupLastPadsThreshold() const { return mGroupLastPadsThreshold; }

  /// \return returns the region of the CRU for which the IDCs are grouped
  unsigned int getRegion() const { return mRegion; }

  /// \return returns the number of the CRU for which the IDCs are grouped
  uint32_t getCRU() const { return mCRU; }

  /// draw grouped IDCs
  /// \param integrationInterval integration interval for which the IDCs will be drawn
  void draw(const int integrationInterval = 0) const;

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName, const char* outName = "IDCGroup") const;

 private:
  const unsigned int mGroupPads{4};                                            ///< group 4 pads
  const unsigned int mGroupRows{4};                                            ///< group 4 pads -> 4x4
  const unsigned int mGroupLastRowsThreshold{2};                               ///< if the last group (region edges) consists in row direction less then mGroupLastRowsThreshold pads then it will be grouped into the previous group
  const unsigned int mGroupLastPadsThreshold{2};                               ///< if the last group (sector edges) consists in pad direction less then mGroupLastPadsThreshold pads then it will be grouped into the previous group
  const unsigned int mCRU{};                                                   ///< CRU of input IDCs
  const unsigned int mRegion{mCRU % IDCHelper::NREGIONS};                                 ///< region of input IDCs
  std::vector<float> mIDCs{};                                                  ///< integrated IDC values
  IDCGroup mIDCsGrouped{};                                                     ///< grouped and averaged IDC values

  int getLastRow() const;
  int getLastPad(const int row) const;
  void initIDCGroup();
};

} // namespace tpc
} // namespace o2

#endif
