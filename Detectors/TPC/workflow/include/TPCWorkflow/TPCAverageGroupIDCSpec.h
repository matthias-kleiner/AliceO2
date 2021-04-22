// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TPCAverageGroupIDCSpec.h
/// \brief TPC merging and averaging of IDCs
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>
/// \date Apr 16, 2021

#ifndef O2_TPCAVERAGEGROUPIDCSPEC_H
#define O2_TPCAVERAGEGROUPIDCSPEC_H

#include <vector>
#include <fmt/format.h>

#include "Framework/Task.h"
#include "Framework/ControlService.h"
#include "Framework/Logger.h"
#include "Framework/DataProcessorSpec.h"
#include "Headers/DataHeader.h"
#include "TPCReconstruction/IDCAverageGroup.h"
#include "TPCWorkflow/TPCIntegrateIDCSpec.h"

using namespace o2::framework;
using o2::header::gDataOriginTPC;
using namespace o2::tpc;

namespace o2
{
namespace tpc
{

/*
class for grouping the IDCs and averaging them in those groups
*/

class TPCAverageGroupIDCDevice : public o2::framework::Task
{
 public:
  TPCAverageGroupIDCDevice(const uint32_t lane, const std::vector<uint32_t>& crus, const TPCIntegrateIDCDevice::IDCFormat inputFormat, const int groupPads, const int groupRows, const int groupLastRowsThreshold, const int groupLastPadsThreshold, const bool debug = false)
    : mLane{lane}, mCRUs{crus}, mInputFormat{inputFormat}, mGroupPads{groupPads}, mGroupRows{groupRows}, mGroupLastRowsThreshold{groupLastRowsThreshold}, mGroupLastPadsThreshold{groupLastPadsThreshold}, mDebug{debug} {}

  void init(o2::framework::InitContext& ic) final
  {
    for (const auto& cru : mCRUs) {
      mIDCs.emplace(cru, IDCAverageGroup(mGroupPads, mGroupRows, mGroupLastRowsThreshold, mGroupLastPadsThreshold, cru));
    }
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    const int nCRUs = mCRUs.size();
    for (int i = 0; i < nCRUs; ++i) {
      const DataRef ref = pc.inputs().getByPos(i);
      auto const* tpcCRUHeader = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
      const int cru = tpcCRUHeader->subSpecification >> 7;
      LOGP(info, "averaging and grouping digits for one TF for CRU: {} ", cru);
      mIDCs[cru].setIDCs(pc.inputs().get<std::vector<float>>(ref));
      mIDCs[cru].processIDCs(mDebug, mLane);

      // send the output for one CRU for one TF
      sendOutput(pc.outputs(), cru);
    }
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
    LOGP(info, "endOfStream");
    ec.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  }

 private:
  const uint32_t mLane{0};                                   ///< lane number of processor
  const std::vector<uint32_t> mCRUs{};                       ///< CRUs to process in this instance
  const TPCIntegrateIDCDevice::IDCFormat mInputFormat{};     ///< type of the input format. Sim=simulation, Real=realistic format
  std::unordered_map<unsigned int, IDCAverageGroup> mIDCs{}; ///< object for averaging and grouping the IDCs
  const int mGroupPads{4};                                   ///< group 4 pads
  const int mGroupRows{4};                                   ///< group 4 pads -> 4x4
  const int mGroupLastRowsThreshold{2};                      ///< if the last group (region edges) consists in row direction less then mGroupLastRowsThreshold pads then it will be grouped into the previous group
  const int mGroupLastPadsThreshold{2};                      ///< if the last group (sector edges) consists in pad direction less then mGroupLastPadsThreshold pads then it will be grouped into the previous group
  const bool mDebug{false};                                  ///< dump IDCs to tree for debugging

  void sendOutput(DataAllocator& output, const uint32_t cru)
  {
    const header::DataHeader::SubSpecificationType subSpec{cru << 7};
    output.snapshot(Output{gDataOriginTPC, "IDCGROUP", subSpec, Lifetime::Timeframe}, mIDCs[cru].getIDCGroup());
  }
};

DataProcessorSpec getTPCAverageGroupIDCSpec(const uint32_t ilane = 0, const std::vector<uint32_t>& crus = {}, const TPCIntegrateIDCDevice::IDCFormat inputFormat = {}, const int groupPads = {}, const int groupRows = {}, const int groupLastRowsThreshold = {}, const int groupLastPadsThreshold = {}, const bool debug = false)
{
  std::vector<OutputSpec> outputSpecs;
  outputSpecs.reserve(crus.size());

  std::vector<InputSpec> inputSpecs;
  inputSpecs.reserve(crus.size());

  for (const auto& cru : crus) {
    const header::DataHeader::SubSpecificationType subSpec{cru << 7};
    inputSpecs.emplace_back(InputSpec{"idcs", gDataOriginTPC, TPCIntegrateIDCDevice::getDataDescription(inputFormat), subSpec, Lifetime::Timeframe});
    outputSpecs.emplace_back(ConcreteDataMatcher{gDataOriginTPC, "IDCGROUP", subSpec});
  }

  const auto id = fmt::format("tpc-averagegroup-idc-{:02}", ilane);
  return DataProcessorSpec{
    id.data(),
    inputSpecs,
    outputSpecs,
    AlgorithmSpec{adaptFromTask<TPCAverageGroupIDCDevice>(ilane, crus, inputFormat, groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, debug)},
  }; // end DataProcessorSpec
}

} // namespace tpc
} // namespace o2

#endif
