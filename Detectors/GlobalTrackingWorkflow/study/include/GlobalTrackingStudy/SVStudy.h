// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_SV_STUDY_H
#define O2_SV_STUDY_H

#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "Framework/Task.h"
#include "Framework/DataProcessorSpec.h"
#include "ReconstructionDataFormats/Track.h"
#include "MathUtils/detail/Bracket.h"
#include "DataFormatsTPC/ClusterNative.h"

namespace o2::svstudy
{
/// create a processor spec; fullOutput: also send the selected V0s with the raw clusters of their prongs (V0WithClusters.h, see
/// V0FullOutputParams for the selection), written by getV0FullWriterSpec
o2::framework::DataProcessorSpec getSVStudySpec(o2::dataformats::GlobalTrackID::mask_t srcTracks, o2::dataformats::GlobalTrackID::mask_t srcCls, bool useMC, bool fullOutput = false);

/// writer of the --full-output V0s to o2_v0s_full.root
o2::framework::DataProcessorSpec getV0FullWriterSpec();

} // namespace o2::svstudy

#endif
