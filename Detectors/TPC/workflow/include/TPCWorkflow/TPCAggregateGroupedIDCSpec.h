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
#include "TPCCalibration/IDCFactorization.h"
#include "CCDB/CcdbApi.h"
#include "Framework/ConfigParamRegistry.h"

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
  TPCAggregateGroupedIDCSpec(const int lane, const std::vector<uint32_t>& crus, const unsigned int timeframes, std::array<unsigned int, Mapper::NREGIONS> groupPads, std::array<unsigned int, Mapper::NREGIONS> groupRows,
                             std::array<unsigned int, Mapper::NREGIONS> groupLastRowsThreshold, std::array<unsigned int, Mapper::NREGIONS> groupLastPadsThreshold, const bool debug = false) : mLane{lane}, mCRUs{crus}, mIDCs{groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, timeframes}, mDebug{debug} {};

  void init(o2::framework::InitContext& ic) final
  {
    mDBapi.init(ic.options().get<std::string>("ccdb-uri")); // or http://localhost:8080 for a local installation
    mWriteToDB = mDBapi.isHostReachable() ? true : false;
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    const int nCRUs = mCRUs.size();
    for (int i = 0; i < nCRUs; ++i) {
      const DataRef ref = pc.inputs().getByPos(i);
      auto const* tpcCRUHeader = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
      const int cru = tpcCRUHeader->subSpecification >> 7;
      mIDCs.setIDCs(pc.inputs().get<std::vector<float>>(ref), cru, mProcessedTFs);
    }
    ++mProcessedTFs;

    if (mProcessedTFs == mIDCs.getNTimeframes()) {
      mProcessedTFs = 0;
      mIDCs.factorizeIDCs();
      if (mWriteToDB) {
        mDBapi.storeAsTFileAny(&mIDCs.getIDCZeroOne(), "TPC/Calib/IDC", mMetadata);
        mDBapi.storeAsTFileAny(&mIDCs.getIDCDelta(), "TPC/Calib/IDC", mMetadata);
      }
      if (mDebug) {
        LOGP(info, "dumping aggregated IDCS to file");
        const DataRef ref = pc.inputs().getByPos(0);
        auto const* tpcCRUHeader = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
        const auto tf = tpcCRUHeader->tfCounter;
        mIDCs.dumpToFile(Form("IDCFactorized_%i.root", tf));
      }
    }
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
    LOGP(info, "endOfStream");
    ec.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  }

 private:
  const int mLane{0};                           ///< lane number of processor
  const std::vector<uint32_t> mCRUs{};          ///< CRUs to process in this instance
  int mProcessedTFs{0};                         ///< number of processed time frames to keep track of when the writing to CCDB will be done
  IDCFactorization mIDCs{};                     ///< object aggregating the IDCs and performing the factorization of the IDCs
  const bool mDebug{false};                     ///< dump IDCs to tree for debugging
  o2::ccdb::CcdbApi mDBapi;                     ///< object for storing the IDCs at CCDB
  std::map<std::string, std::string> mMetadata; ///< meta data of the stored object in CCDB
  bool mWriteToDB{};                            ///< flag if writing to CCDB will be done
  // void sendOutput(DataAllocator& output) {}
};

DataProcessorSpec getTPCAggregateGroupedIDCSpec(const int ilane, const std::vector<uint32_t>& crus, const int timeframes, const bool debug = false)
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

  const auto& paramIDCGroup = ParameterIDCGroup::Instance();
  std::array<unsigned int, Mapper::NREGIONS> groupPads{};
  std::array<unsigned int, Mapper::NREGIONS> groupRows{};
  std::array<unsigned int, Mapper::NREGIONS> groupLastRowsThreshold{};
  std::array<unsigned int, Mapper::NREGIONS> groupLastPadsThreshold{};
  std::copy(std::begin(paramIDCGroup.GroupPads), std::end(paramIDCGroup.GroupPads), std::begin(groupPads));
  std::copy(std::begin(paramIDCGroup.GroupRows), std::end(paramIDCGroup.GroupRows), std::begin(groupRows));
  std::copy(std::begin(paramIDCGroup.GroupLastRowsThreshold), std::end(paramIDCGroup.GroupLastRowsThreshold), std::begin(groupLastRowsThreshold));
  std::copy(std::begin(paramIDCGroup.GroupLastPadsThreshold), std::end(paramIDCGroup.GroupLastPadsThreshold), std::begin(groupLastPadsThreshold));

  const auto id = fmt::format("tpc-aggregate-IDC-{:02}", ilane);
  return DataProcessorSpec{
    id.data(),
    inputSpecs,
    outputSpecs,
    AlgorithmSpec{adaptFromTask<TPCAggregateGroupedIDCSpec>(ilane, crus, timeframes, groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, debug)},
    Options{{"ccdb-uri", VariantType::String, "http://ccdb-test.cern.ch:8080", {"URI for the CCDB access."}}}}; // end DataProcessorSpec
}

} // namespace tpc
} // namespace o2

#endif
