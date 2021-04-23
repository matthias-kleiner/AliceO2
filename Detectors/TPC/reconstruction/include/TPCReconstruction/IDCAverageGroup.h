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

namespace o2
{
namespace tpc
{
class IDCAverageGroup
{

 public:
  IDCAverageGroup(const unsigned int groupPads = 4, const unsigned int groupRows = 4, const unsigned int groupLastRowsThreshold = 2, const unsigned int groupLastPadsThreshold = 2, const unsigned int cru = 0)
    : mIDCsGrouped{groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, cru % IDCHelper::NREGIONS} {}

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
  unsigned int getNIntegrationIntervals() const { return mIDCs.size() / IDCHelper::PADSPERREGION[mIDCsGrouped.getRegion()]; }

  /// grouping and averaging of IDCs
  void processIDCs();

  /// \return returns grouped IDC object
  const auto& getIDCGroup() const { return mIDCsGrouped; }

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName, const char* outName = "IDCGroup") const;

 private:
  std::vector<float> mIDCs{}; ///< integrated IDC values
  IDCGroup mIDCsGrouped{};    ///< grouped and averaged IDC values
};

} // namespace tpc
} // namespace o2

#endif
