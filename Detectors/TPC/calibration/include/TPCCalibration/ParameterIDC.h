// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file ParameterIDCGroup.h
/// \brief Definition of the parameter for the grouping of the IDCs
/// \author Matthias Kleiner, mkleiner@ikf.uni-frankfurt.de

#ifndef ALICEO2_TPC_PARAMETERIDCGROUP_H_
#define ALICEO2_TPC_PARAMETERIDCGROUP_H_

#include <array>
#include <cmath>
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include "TPCBase/Mapper.h"

namespace o2
{
namespace tpc
{
/// struct for setting the parameters for the grouping of IDCs
struct ParameterIDCGroup : public o2::conf::ConfigurableParamHelper<ParameterIDCGroup> {
  unsigned int GroupPads[Mapper::NREGIONS]{7, 7, 7, 7, 6, 6, 6, 6, 5, 5};              ///< group 4 pads
  unsigned int GroupRows[Mapper::NREGIONS]{5, 5, 5, 5, 4, 4, 4, 4, 3, 3};              ///< group 4 pads -> 4x4
  unsigned int GroupLastRowsThreshold[Mapper::NREGIONS]{3, 3, 3, 3, 2, 2, 2, 2, 2, 2}; ///< if the last group (region edges) consists in row direction less then mGroupLastRowsThreshold pads then it will be grouped into the previous group
  unsigned int GroupLastPadsThreshold[Mapper::NREGIONS]{3, 3, 3, 3, 2, 2, 2, 2, 1, 1}; ///< if the last group (sector edges) consists in pad direction less then mGroupLastPadsThreshold pads then it will be grouped into the previous group
  O2ParamDef(ParameterIDCGroup, "TPCIDCGroupParam");
};

/// struct for storing the parameters for the grouping of IDCs to CCDB
struct ParameterIDCGroupCCDB {

  /// contructor
  /// \param groupPads number of pads in pad direction which are grouped
  /// \param groupRows number of pads in row direction which are grouped
  /// \param groupLastRowsThreshold minimum number of pads in row direction for the last group in row direction
  /// \param groupLastPadsThreshold minimum number of pads in pad direction for the last group in pad direction
  ParameterIDCGroupCCDB(const std::array<unsigned int, Mapper::NREGIONS>& groupPads, const std::array<unsigned int, Mapper::NREGIONS>& groupRows, const std::array<unsigned int, Mapper::NREGIONS>& groupLastRowsThreshold, const std::array<unsigned int, Mapper::NREGIONS>& groupLastPadsThreshold)
  {
    setParameters(groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold);
  }

  ParameterIDCGroupCCDB() = default;

  void setParameters(const std::array<unsigned int, Mapper::NREGIONS>& groupPads, const std::array<unsigned int, Mapper::NREGIONS>& groupRows, const std::array<unsigned int, Mapper::NREGIONS>& groupLastRowsThreshold, const std::array<unsigned int, Mapper::NREGIONS>& groupLastPadsThreshold)
  {
    for (unsigned int i = 0; i < Mapper::NREGIONS; ++i) {
      GroupPads[i] = static_cast<unsigned char>(groupPads[i]);
      GroupRows[i] = static_cast<unsigned char>(groupRows[i]);
      GroupLastRowsThreshold[i] = static_cast<unsigned char>(groupLastRowsThreshold[i]);
      GroupLastPadsThreshold[i] = static_cast<unsigned char>(groupLastPadsThreshold[i]);
    }
  }

  /// \return returns number of pads in pad direction which are grouped
  /// \parameter region TPC region
  unsigned int getGroupPads(const unsigned int region) const { return static_cast<unsigned int>(GroupPads[region]); }

  /// \return returns number of pads in row direction which are grouped
  /// \parameter region TPC region
  unsigned int getGroupRows(const unsigned int region) const { return static_cast<unsigned int>(GroupRows[region]); }

  /// \return returns minimum number of pads in row direction for the last group in row direction
  /// \parameter region TPC region
  unsigned int getGroupLastRowsThreshold(const unsigned int region) const { return static_cast<unsigned int>(GroupLastRowsThreshold[region]); }

  /// \return returns minimum number of pads in pad direction for the last group in pad direction
  /// \parameter region TPC region
  unsigned int getGroupLastPadsThreshold(const unsigned int region) const { return static_cast<unsigned int>(GroupLastPadsThreshold[region]); }

  /// \return returns number of pads in pad direction which are grouped for all regions
  std::array<unsigned int, Mapper::NREGIONS> getGroupPads() const
  {
    std::array<unsigned int, Mapper::NREGIONS> groupPadsTmp{};
    for (unsigned int i = 0; i < Mapper::NREGIONS; ++i) {
      groupPadsTmp[i] = static_cast<unsigned int>(GroupPads[i]);
    }
    return groupPadsTmp;
  }

  /// \return returns number of pads in row direction which are grouped for all regions
  std::array<unsigned int, Mapper::NREGIONS> getGroupRows() const
  {
    std::array<unsigned int, Mapper::NREGIONS> groupRowsTmp{};
    for (unsigned int i = 0; i < Mapper::NREGIONS; ++i) {
      groupRowsTmp[i] = static_cast<unsigned int>(GroupRows[i]);
    }
    return groupRowsTmp;
  }

  /// \return returns minimum number of pads in row direction for the last group in row direction for all regions
  std::array<unsigned int, Mapper::NREGIONS> getGroupLastRowsThreshold() const
  {
    std::array<unsigned int, Mapper::NREGIONS> groupLastRowsThresholdTmp{};
    for (unsigned int i = 0; i < Mapper::NREGIONS; ++i) {
      groupLastRowsThresholdTmp[i] = static_cast<unsigned int>(GroupLastRowsThreshold[i]);
    }
    return groupLastRowsThresholdTmp;
  }

  /// \return returns minimum number of pads in pad direction for the last group in pad direction for all regions
  std::array<unsigned int, Mapper::NREGIONS> getGroupLastPadsThreshold() const
  {
    std::array<unsigned int, Mapper::NREGIONS> groupLastPadsThresholdTmp{};
    for (unsigned int i = 0; i < Mapper::NREGIONS; ++i) {
      groupLastPadsThresholdTmp[i] = static_cast<unsigned int>(GroupLastPadsThreshold[i]);
    }
    return groupLastPadsThresholdTmp;
  }

  std::array<unsigned char, Mapper::NREGIONS> GroupPads{};              ///< grouping parameter in pad direction (how many pads in pad direction are grouped)
  std::array<unsigned char, Mapper::NREGIONS> GroupRows{};              ///< grouping parameter in pad direction (how many pads in pad direction are grouped)
  std::array<unsigned char, Mapper::NREGIONS> GroupLastRowsThreshold{}; ///< if the last group (region edges) consists in row direction less then mGroupLastRowsThreshold pads then it will be grouped into the previous group
  std::array<unsigned char, Mapper::NREGIONS> GroupLastPadsThreshold{}; ///< if the last group (sector edges) consists in pad direction less then mGroupLastPadsThreshold pads then it will be grouped into the previous group
};

struct ParameterIDCCompression : public o2::conf::ConfigurableParamHelper<ParameterIDCCompression> {
  float MaxIDCDeltaValue = 0.3f; ///< maximum Delta IDC
  O2ParamDef(ParameterIDCCompression, "TPCIDCCompressionParam");
};

} // namespace tpc
} // namespace o2

#endif // ALICEO2_TPC_ParameterGEM_H_
