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
#include "TPCCalibration/IDCAverageGroup.h"
#include "TPCWorkflow/TPCIntegrateIDCSpec.h"
#include "TPCBase/Mapper.h"
#include "TPCCalibration/ParameterIDC.h"

using namespace o2::framework;
using o2::header::gDataOriginTPC;
using namespace o2::tpc;

namespace o2::tpc
{

class TPCAverageGroupIDCDevice : public o2::framework::Task
{
 public:
  TPCAverageGroupIDCDevice(const int lane, const std::vector<uint32_t>& crus, const TPCIntegrateIDCDevice::IDCFormat inputFormat, const bool debug = false)
    : mLane{lane}, mCRUs{crus}, mInputFormat{inputFormat}, mDebug{debug}
  {
    auto& paramIDCGroup = ParameterIDCGroup::Instance();
    for (const auto& cru : mCRUs) {
      const CRU cruTmp(cru);
      const unsigned int reg = cruTmp.region();
      mIDCs.emplace(cru, IDCAverageGroup(paramIDCGroup.GroupPads[reg], paramIDCGroup.GroupRows[reg], paramIDCGroup.GroupLastRowsThreshold[reg], paramIDCGroup.GroupLastPadsThreshold[reg], reg, cruTmp.sector()));
    }
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    LOGP(info, "averaging and grouping IDCs for one TF for CRUs {} to {} using {} threads", mCRUs.front(), mCRUs.back(), mIDCs.begin()->second.getNThreads());

    for (int i = 0; i < mCRUs.size(); ++i) {
      const DataRef ref = pc.inputs().getByPos(i);
      auto const* tpcCRUHeader = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
      const int cru = tpcCRUHeader->subSpecification >> 7;
      mIDCs[cru].setIDCs(pc.inputs().get<std::vector<float>>(ref));
      mIDCs[cru].processIDCs();

      // send the output for one CRU for one TF
      sendOutput(pc.outputs(), cru);
    }

    if (mDebug) {
      TFile fOut(fmt::format("IDCGroup_{}_tf_{}.root", mLane, o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(pc.inputs().getByPos(0))->tfCounter).data(), "RECREATE");
      for (int i = 0; i < mCRUs.size(); ++i) {
        const DataRef ref = pc.inputs().getByPos(i);
        auto const* tpcCRUHeader = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
        const int cru = tpcCRUHeader->subSpecification >> 7;
        fOut.WriteObject(&mIDCs[cru], fmt::format("CRU_{}", cru).data());
      }
    }
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
    LOGP(info, "endOfStream");
    ec.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  }

 private:
  const int mLane{0};                                        ///< lane number of processor
  const std::vector<uint32_t> mCRUs{};                       ///< CRUs to process in this instance
  const TPCIntegrateIDCDevice::IDCFormat mInputFormat{};     ///< type of the input format. Sim=simulation, Real=realistic format
  const bool mDebug{};                                       ///< dump IDCs to tree for debugging
  std::unordered_map<unsigned int, IDCAverageGroup> mIDCs{}; ///< object for averaging and grouping the IDCs

  void sendOutput(DataAllocator& output, const uint32_t cru)
  {
    const header::DataHeader::SubSpecificationType subSpec{cru << 7};
    output.snapshot(Output{gDataOriginTPC, "IDCGROUP", subSpec, Lifetime::Timeframe}, mIDCs[cru].getIDCGroup().getData());
  }
};

DataProcessorSpec getTPCAverageGroupIDCSpec(const int ilane = 0, const std::vector<uint32_t>& crus = {}, const TPCIntegrateIDCDevice::IDCFormat inputFormat = {}, const bool debug = false)
{
  std::vector<OutputSpec> outputSpecs;
  std::vector<InputSpec> inputSpecs;
  outputSpecs.reserve(crus.size());
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
    AlgorithmSpec{adaptFromTask<TPCAverageGroupIDCDevice>(ilane, crus, inputFormat, debug)},
  }; // end DataProcessorSpec
}

} // namespace o2::tpc

#endif