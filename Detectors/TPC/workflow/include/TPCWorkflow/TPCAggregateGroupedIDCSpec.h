// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TPCAggregateGroupedIDCSpec.h
/// \brief TPC aggregation of grouped IDCs
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>
/// \date Apr 22, 2021

#ifndef O2_TPCAGGREGATEGROUPIDCSPEC_H
#define O2_TPCAGGREGATEGROUPIDCSPEC_H

#include <vector>
#include <fmt/format.h>

#include "Framework/Task.h"
#include "Framework/ControlService.h"
#include "Framework/Logger.h"
#include "Framework/DataProcessorSpec.h"
#include "CommonUtils/MemFileHelper.h"
#include "Headers/DataHeader.h"
#include "TPCBase/CDBInterface.h"
#include "TPCBase/IDCGroup.h"

using namespace o2::framework;
using o2::header::gDataOriginTPC;
using namespace o2::tpc;

namespace o2
{
namespace tpc
{

/*
class for aggregation of grouped IDCs
*/

class TPCAggregateGroupedIDCSpec : public o2::framework::Task
{
 public:
  TPCAggregateGroupedIDCSpec(const int lane, const std::vector<uint32_t>& crus, const bool debug = false) : mLane{lane}, mCRUs{crus}, mDebug{debug} {}

  void init(o2::framework::InitContext& ic) final
  {
    LOGP(info, "init TPCAggregateGroupedIDCSpec");
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    LOGP(info, "RUN aggreating IDCs for CRU ");

    const int nCRUs = mCRUs.size();
    for (int i = 0; i < nCRUs; ++i) {
      const DataRef ref = pc.inputs().getByPos(i);
      auto const* tpcCRUHeader = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
      const int cru = tpcCRUHeader->subSpecification >> 7;
      auto idcGroup = pc.inputs().get<IDCGroup*>(ref);

      const int nrows = idcGroup->getNRows();
      LOGP(info, "nrows: {} ", nrows);
    }
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
    LOGP(info, "endOfStream");
    ec.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  }

 private:
  const int mLane{0};                                 ///< lane number of processor
  const std::vector<uint32_t> mCRUs{};                ///< CRUs to process in this instance
  std::unordered_map<unsigned int, IDCGroup> mIDCs{}; ///< object for averaging and grouping the IDCs
  const bool mDebug{false};                           ///< dump IDCs to tree for debugging

  void sendOutput(DataAllocator& output)
  {
    for (const auto& [cru, idcs] : mIDCs) {
      const header::DataHeader::SubSpecificationType subSpec{cru << 7};
      // auto image = o2::utils::MemFileHelper::createFileImage(idcs.getIDCGroup(), typeid(*idcs.getIDCGroup()), cal->getName(), "data");
      // output.snapshot(Output{gDataOriginTPC, "IDCGroup", subSpec, Lifetime::Timeframe}, idcs.getIDCGroup());
    }
  }
};

DataProcessorSpec getTPCAggregateGroupedIDCSpec(const int ilane = 0, const std::vector<uint32_t>& crus = {}, const bool debug = false)
{
  std::vector<OutputSpec> outputSpecs;
  outputSpecs.reserve(crus.size());

  std::vector<InputSpec> inputSpecs;
  inputSpecs.reserve(crus.size());

  for (const auto& cru : crus) {
    const header::DataHeader::SubSpecificationType subSpec{cru << 7};
    inputSpecs.emplace_back(InputSpec{"idcsgroup", gDataOriginTPC, "IDCGROUP", subSpec, Lifetime::Timeframe});
    outputSpecs.emplace_back(ConcreteDataMatcher{gDataOriginTPC, "IDCAGGREGATED", subSpec});
  }

  const auto id = fmt::format("tpc-aggregate-IDC-{:02}", ilane);
  return DataProcessorSpec{
    id.data(),
    inputSpecs,
    outputSpecs,
    AlgorithmSpec{adaptFromTask<TPCAggregateGroupedIDCSpec>(ilane, crus, debug)},
  }; // end DataProcessorSpec
}

} // namespace tpc
} // namespace o2

#endif
