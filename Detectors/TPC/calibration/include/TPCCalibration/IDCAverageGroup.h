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
/// \brief class for averaging/grouping and storing of IDCs
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_IDCAVERAGEGROUP_H_
#define ALICEO2_IDCAVERAGEGROUP_H_

#include <vector>
#include "TPCCalibration/IDCGroup.h"
#include "TPCBase/Mapper.h"

namespace o2
{
namespace tpc
{
class IDCAverageGroup
{

 public:
  IDCAverageGroup(const unsigned int groupPads = 4, const unsigned int groupRows = 4, const unsigned int groupLastRowsThreshold = 2, const unsigned int groupLastPadsThreshold = 2, const unsigned int region = 0)
    : mIDCsGrouped{groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, region} {}

  /// \param IDCs vector containing the IDCs
  void setIDCs(const std::vector<float>& idcs)
  {
    mIDCsUngrouped = idcs;
    mIDCsGrouped.resize(getNIntegrationIntervals());
  }

  /// \param IDCs vector containing the IDCs
  void setIDCs(std::vector<float>&& idcs)
  {
    mIDCsUngrouped = std::move(idcs);
    mIDCsGrouped.resize(getNIntegrationIntervals());
  }

  /// \return returns number of integration intervalls stored in this object
  unsigned int getNIntegrationIntervals() const { return mIDCsUngrouped.size() / Mapper::PADSPERREGION[mIDCsGrouped.getRegion()]; }

  /// grouping and averaging of IDCs
  void processIDCs();

  /// \return returns grouped IDC object
  const auto& getIDCGroup() const { return mIDCsGrouped; }

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName = "IDCAverageGroup.root", const char* outName = "IDCAverageGroup") const;

  /// draw ungrouped IDCs
  /// \param integrationInterval integration interval for which the IDCs will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawUngroupedIDCs(const unsigned int integrationInterval = 0, const std::string filename = "IDCsUngrouped.pdf") const;

  /// draw grouped IDCs
  /// \param integrationInterval integration interval for which the IDCs will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawGroupedIDCs(const unsigned int integrationInterval = 0, const std::string filename = "IDCsGrouped.pdf") const { mIDCsGrouped.draw(integrationInterval, filename); }

  /// \return returns the stored ungrouped IDC value for local ungrouped pad row and ungrouped pad
  /// \param urow ungrouped local row in region
  /// \param upad ungrouped pad in pad direction
  /// \param integrationInterval integration interval for which the IDCs will be returned
  const float& getUngroupedIDCVal(const unsigned int urow, const unsigned int upad, const unsigned int integrationInterval) const
  {
    const unsigned int indexIDC = integrationInterval * Mapper::PADSPERREGION[mIDCsGrouped.getRegion()] + Mapper::OFFSETCRULOCAL[mIDCsGrouped.getRegion()][urow] + upad;
    return mIDCsUngrouped[indexIDC];
  }

  /// \return returns the stored ungrouped IDC value for local ungrouped pad row and ungrouped pad
  /// \param urow ungrouped local row in region
  /// \param upad ungrouped pad in pad direction
  /// \param integrationInterval integration interval for which the IDCs will be returned
  float& getUngroupedIDCVal(const unsigned int urow, const unsigned int upad, const unsigned int integrationInterval)
  {
    const unsigned int indexIDC = integrationInterval * Mapper::PADSPERREGION[mIDCsGrouped.getRegion()] + Mapper::OFFSETCRULOCAL[mIDCsGrouped.getRegion()][urow] + upad;
    return mIDCsUngrouped[indexIDC];
  }

  /// \return returns the stored grouped IDC value for local ungrouped pad row and ungrouped pad
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  float& getGroupedIDCVal(unsigned int urow, unsigned int upad, unsigned int integrationInterval) { return mIDCsGrouped.getVal(urow, upad, integrationInterval); }

  /// \return returns the stored grouped IDC value for local ungrouped pad row and ungrouped pad
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  const float& getGroupedIDCVal(unsigned int urow, unsigned int upad, unsigned int integrationInterval) const { return mIDCsGrouped.getVal(urow, upad, integrationInterval); }


 private:
  std::vector<float> mIDCsUngrouped{}; ///< integrated ungrouped IDC values per pad
  IDCGroup mIDCsGrouped{};             ///< grouped and averaged IDC values
};

} // namespace tpc
} // namespace o2

#endif
