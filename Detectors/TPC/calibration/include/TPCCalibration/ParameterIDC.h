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
struct ParameterIDCGroup : public o2::conf::ConfigurableParamHelper<ParameterIDCGroup> {
  unsigned int GroupPads[Mapper::NREGIONS]{7, 7, 7, 7, 6, 6, 6, 6, 5, 5};              ///< group 4 pads
  unsigned int GroupRows[Mapper::NREGIONS]{5, 5, 5, 5, 4, 4, 4, 4, 3, 3};              ///< group 4 pads -> 4x4
  unsigned int GroupLastRowsThreshold[Mapper::NREGIONS]{3, 3, 3, 3, 2, 2, 2, 2, 2, 2}; ///< if the last group (region edges) consists in row direction less then mGroupLastRowsThreshold pads then it will be grouped into the previous group
  unsigned int GroupLastPadsThreshold[Mapper::NREGIONS]{3, 3, 3, 3, 2, 2, 2, 2, 1, 1}; ///< if the last group (sector edges) consists in pad direction less then mGroupLastPadsThreshold pads then it will be grouped into the previous group
  O2ParamDef(ParameterIDCGroup, "TPCIDCGroupParam");
};

struct ParameterIDCCompression : public o2::conf::ConfigurableParamHelper<ParameterIDCCompression> {
  float MaxIDCDeltaValue = 0.3f; ///< maximum Delta IDC
  O2ParamDef(ParameterIDCCompression, "TPCIDCCompressionParam");
};

} // namespace tpc
} // namespace o2

#endif // ALICEO2_TPC_ParameterGEM_H_
