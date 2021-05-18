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
#include "TPCCalibration/ParameterIDC.h"

using namespace o2::framework;
using o2::header::gDataOriginTPC;
using namespace o2::tpc;

namespace o2::tpc
{

class TPCAggregateGroupedIDCSpec : public o2::framework::Task
{
 public:
  TPCAggregateGroupedIDCSpec(const std::vector<uint32_t>& crus, const unsigned int timeframes, std::array<unsigned int, Mapper::NREGIONS> groupPads,
                             std::array<unsigned int, Mapper::NREGIONS> groupRows, std::array<unsigned int, Mapper::NREGIONS> groupLastRowsThreshold,
                             std::array<unsigned int, Mapper::NREGIONS> groupLastPadsThreshold, const IDCDeltaCompression compression, const bool debug = false)
    : mCRUs{crus}, mIDCs{groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, timeframes}, mCompressionDeltaIDC{compression}, mDebug{debug} {};

  void init(o2::framework::InitContext& ic) final
  {
    // mDBapi.init(ic.options().get<std::string>("ccdb-uri")); // or http://localhost:8080 for a local installation
    mDBapi.init("http://localhost:8080"); // or http://localhost:8080 for a local installation
    mWriteToDB = mDBapi.isHostReachable() ? true : false;

    if (mWriteToDB) {
      // write struct containing grouping parameters to access grouped IDCs to CCDB
      const ParameterIDCGroupCCDB parGrouping(mIDCs.getGroupPads(), mIDCs.getGroupRows(), mIDCs.getPadThreshold(), mIDCs.getRowThreshold());
      mDBapi.storeAsTFileAny(&parGrouping, "TPC/Calib/IDC/GROUPINGPAR", mMetadata);
    }
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

      sendOutput(pc.outputs());

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
  const std::vector<uint32_t> mCRUs{};              ///< CRUs to process in this instance
  int mProcessedTFs{0};                             ///< number of processed time frames to keep track of when the writing to CCDB will be done
  IDCFactorization mIDCs{};                         ///< object aggregating the IDCs and performing the factorization of the IDCs
  const IDCDeltaCompression mCompressionDeltaIDC{}; ///< compression type of IDC delta
  const bool mDebug{false};                         ///< dump IDCs to tree for debugging
  o2::ccdb::CcdbApi mDBapi;                         ///< object for storing the IDCs at CCDB
  std::map<std::string, std::string> mMetadata;     ///< meta data of the stored object in CCDB
  bool mWriteToDB{};                                ///< flag if writing to CCDB will be done

  void sendOutput(DataAllocator& output)
  {
    if (mWriteToDB) {
      // store IDC Zero One in CCDB
      mDBapi.storeAsTFileAny(&mIDCs.getIDCZeroOne(), "TPC/Calib/IDC/IDCZEROONE", mMetadata);
    }

    for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
      const o2::tpc::Side side = iSide ? Side::C : Side::A;
      const header::DataHeader::SubSpecificationType subSpec{iSide};
      output.snapshot(Output{gDataOriginTPC, "IDCZERO", subSpec, Lifetime::Timeframe}, mIDCs.getIDCZero(side));
      output.snapshot(Output{gDataOriginTPC, "IDCONE", subSpec, Lifetime::Timeframe}, mIDCs.getIDCOne(side));
    }

    switch (mCompressionDeltaIDC) {
      case IDCDeltaCompression::MEDIUM: {
        auto idcDeltaMediumCompressed = mIDCs.getIDCDeltaMediumCompressed();
        if (mWriteToDB) {
          mDBapi.storeAsTFileAny(&idcDeltaMediumCompressed, "TPC/Calib/IDC/IDCDELTA", mMetadata);
        }
        for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
          const o2::tpc::Side side = iSide ? Side::C : Side::A;
          const header::DataHeader::SubSpecificationType subSpec{iSide};
          output.snapshot(Output{gDataOriginTPC, "IDCDELTA", subSpec, Lifetime::Timeframe}, idcDeltaMediumCompressed.mIDCDelta[side]);
        }
        break;
      }
      case IDCDeltaCompression::HIGH: {
        auto idcDeltaHighCompressed = mIDCs.getIDCDeltaHighCompressed();
        if (mWriteToDB) {
          mDBapi.storeAsTFileAny(&idcDeltaHighCompressed, "TPC/Calib/IDC/IDCDELTA", mMetadata);
        }
        for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
          const o2::tpc::Side side = iSide ? Side::C : Side::A;
          const header::DataHeader::SubSpecificationType subSpec{iSide};
          output.snapshot(Output{gDataOriginTPC, "IDCDELTA", subSpec, Lifetime::Timeframe}, idcDeltaHighCompressed.mIDCDelta[side]);
        }
        break;
      }
      case IDCDeltaCompression::NO:
      default:
        if (mWriteToDB) {
          mDBapi.storeAsTFileAny(&mIDCs.getIDCDeltaUncompressed(), "TPC/Calib/IDC/IDCDELTA", mMetadata);
        }
        for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
          const o2::tpc::Side side = iSide ? Side::C : Side::A;
          const header::DataHeader::SubSpecificationType subSpec{iSide};
          output.snapshot(Output{gDataOriginTPC, "IDCDELTA", subSpec, Lifetime::Timeframe}, mIDCs.getIDCDeltaUncompressed(side));
        }
        break;
    }
  }
};

DataProcessorSpec getTPCAggregateGroupedIDCSpec(const std::vector<uint32_t>& crus, const int timeframes, const IDCDeltaCompression compression, const bool debug = false)
{
  std::vector<OutputSpec> outputSpecs;
  for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
    const header::DataHeader::SubSpecificationType subSpec{iSide};
    outputSpecs.emplace_back(ConcreteDataMatcher{gDataOriginTPC, "IDCZERO", subSpec});
    outputSpecs.emplace_back(ConcreteDataMatcher{gDataOriginTPC, "IDCONE", subSpec});
    outputSpecs.emplace_back(ConcreteDataMatcher{gDataOriginTPC, "IDCDELTA", subSpec});
  }

  std::vector<InputSpec> inputSpecs;
  inputSpecs.reserve(crus.size());
  for (const auto& cru : crus) {
    const header::DataHeader::SubSpecificationType subSpec{cru << 7};
    inputSpecs.emplace_back(InputSpec{"idcsgroup", gDataOriginTPC, "IDCGROUP", subSpec, Lifetime::Timeframe});
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
  const auto id = fmt::format("tpc-aggregate-idc");
  return DataProcessorSpec{
    id.data(),
    inputSpecs,
    outputSpecs,
    AlgorithmSpec{adaptFromTask<TPCAggregateGroupedIDCSpec>(crus, timeframes, groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, compression, debug)},
    Options{{"ccdb-uri", VariantType::String, "http://ccdb-test.cern.ch:8080", {"URI for the CCDB access."}}}}; // end DataProcessorSpec
}

} // namespace o2::tpc

#endif
