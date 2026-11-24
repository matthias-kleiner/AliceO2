// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   V0FullWriterSpec.cxx
/// @brief  Writes the V0s with the raw clusters of their prongs (o2-sv-study-workflow --full-output) to o2_v0s_full.root, one entry per TF

#include <vector>
#include <cstdint>
#include "GlobalTrackingStudy/SVStudy.h"
#include "DPLUtils/MakeRootTreeWriterSpec.h"
#include "DataFormatsGlobalTracking/V0WithClusters.h"
#include "CommonDataFormat/TFIDInfo.h"

using namespace o2::framework;

namespace o2::svstudy
{

template <typename T>
using BranchDefinition = MakeRootTreeWriterSpec::BranchDefinition<T>;

DataProcessorSpec getV0FullWriterSpec()
{
  auto logger = [](std::vector<o2::dataformats::V0WithClusters> const& v0s) {
    LOG(info) << "Writing " << v0s.size() << " V0s with clusters";
  };
  return MakeRootTreeWriterSpec("v0s-full-writer",
                                "o2_v0s_full.root",
                                "v0sFull",
                                -1,  // do not limit number of events to store
                                100, // periodically autosave
                                BranchDefinition<std::vector<o2::dataformats::V0WithClusters>>{InputSpec{"v0s", "GLO", "V0FULL", 0}, "v0s", 1, logger},
                                BranchDefinition<std::vector<o2::dataformats::V0Prong>>{InputSpec{"prongs", "GLO", "V0FULLPRONG", 0}, "prongs", 1},
                                BranchDefinition<std::vector<o2::dataformats::V0TPCCluster>>{InputSpec{"cltpc", "GLO", "V0FULLCLTPC", 0}, "clTPC", 1},
                                BranchDefinition<std::vector<o2::dataformats::V0ITSCluster>>{InputSpec{"clits", "GLO", "V0FULLCLITS", 0}, "clITS", 1},
                                BranchDefinition<std::vector<uint8_t>>{InputSpec{"itspatt", "GLO", "V0FULLITSPATT", 0}, "itsPatterns", 1},
                                BranchDefinition<std::vector<o2::dataformats::V0TOFCluster>>{InputSpec{"cltof", "GLO", "V0FULLCLTOF", 0}, "clTOF", 1},
                                BranchDefinition<std::vector<o2::dataformats::V0TRDTracklet>>{InputSpec{"trdtrkl", "GLO", "V0FULLTRDTRKL", 0}, "trdTracklets", 1},
                                BranchDefinition<o2::dataformats::V0sTFInfo>{InputSpec{"tfinfo", "GLO", "V0FULLTF", 0}, "tfInfo", 1},
                                BranchDefinition<std::vector<unsigned int>>{InputSpec{"occmap", "GLO", "V0FULLOCCMAP", 0}, "tpcOccMap", 1},
                                BranchDefinition<o2::dataformats::TFIDInfo>{InputSpec{"tfid", "GLO", "V0FULLTFID", 0}, "tfID", 1})();
}

} // namespace o2::svstudy
