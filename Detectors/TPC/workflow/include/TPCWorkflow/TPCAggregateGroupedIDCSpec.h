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
#include "TPCCalibration/IDCFourierTransform.h"

using namespace o2::framework;
using o2::header::gDataOriginTPC;
using namespace o2::tpc;

namespace o2::tpc
{

class TPCAggregateGroupedIDCSpec : public o2::framework::Task
{
 public:
  TPCAggregateGroupedIDCSpec(const std::vector<uint32_t>& crus, const unsigned int timeframes, const unsigned int timeframesDeltaIDC, std::array<unsigned int, Mapper::NREGIONS> groupPads,
                             std::array<unsigned int, Mapper::NREGIONS> groupRows, std::array<unsigned int, Mapper::NREGIONS> groupLastRowsThreshold,
                             std::array<unsigned int, Mapper::NREGIONS> groupLastPadsThreshold, const unsigned int rangeIDC, const unsigned int fourierCoefficients, const IDCDeltaCompression compression, const bool debug = false)
    : mCRUs{crus}, mIDCs{groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, timeframes, timeframesDeltaIDC}, mIDCFourierTransform{rangeIDC, fourierCoefficients, timeframes}, mCompressionDeltaIDC{compression}, mDebug{debug} {};

  void init(o2::framework::InitContext& ic) final
  {
    // mDBapi.init(ic.options().get<std::string>("ccdb-uri")); // or http://localhost:8080 for a local installation
    mDBapi.init("http://localhost:8080"); // or http://localhost:8080 for a local installation
    mWriteToDB = mDBapi.isHostReachable() ? true : false;

    if (mWriteToDB) {
      // write struct containing grouping parameters to access grouped IDCs to CCDB
      const ParameterIDCGroupCCDB parGrouping(mIDCs.getGroupPads(), mIDCs.getGroupRows(), mIDCs.getPadThreshold(), mIDCs.getRowThreshold());
      mDBapi.storeAsTFileAny<o2::tpc::ParameterIDCGroupCCDB>(&parGrouping, "TPC/Calib/IDC/GROUPINGPAR", mMetadata);
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

      if (mDebug) {
        // dump this first to also store 1D-IDCs! otherwise they will be moved out and are not available
        LOGP(info, "dumping aggregated IDCS to file");
        const DataRef ref = pc.inputs().getByPos(0);
        auto const* tpcCRUHeader = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
        const auto tf = tpcCRUHeader->tfCounter;
        mIDCs.dumpToFile(fmt::format("IDCFactorized_{:02}.root", tf).data());
      }

      mIDCFourierTransform.setIDCs(std::move(mIDCs).getIDCOne(), mIDCs.getIntegrationIntervalsPerTF()); // using move semantics here
      mIDCFourierTransform.calcFourierCoefficients();

      LOGP(info, "sendOutput calcFourierCoefficients");
      sendOutput(pc.outputs());

      if (mDebug) {
        LOGP(info, "dumping fourier transform  to file");
        const DataRef ref = pc.inputs().getByPos(0);
        auto const* tpcCRUHeader = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
        const auto tf = tpcCRUHeader->tfCounter;
        mIDCFourierTransform.dumpToFile(fmt::format("Fourier_{:02}.root", tf).data());
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
  IDCFourierTransform mIDCFourierTransform{};       ///< object for performing the fourier transform of 1D-IDCs
  const IDCDeltaCompression mCompressionDeltaIDC{}; ///< compression type of IDC delta
  const bool mDebug{false};                         ///< dump IDCs to tree for debugging
  o2::ccdb::CcdbApi mDBapi;                         ///< object for storing the IDCs at CCDB
  std::map<std::string, std::string> mMetadata;     ///< meta data of the stored object in CCDB
  bool mWriteToDB{};                                ///< flag if writing to CCDB will be done

  void sendOutput(DataAllocator& output)
  {
    if (mWriteToDB) {
      // store IDC Zero One in CCDB
      // const long timeStampStart = 33;
      // const long timeStampStart = 33;
      mDBapi.storeAsTFileAny<o2::tpc::IDCZero>(&mIDCs.getIDCZero(), "TPC/Calib/IDC/IDC0", mMetadata);
      mDBapi.storeAsTFileAny<o2::tpc::IDCOne>(&mIDCFourierTransform.getIDCOne(), "TPC/Calib/IDC/IDC1", mMetadata);
      mDBapi.storeAsTFileAny<o2::tpc::FourierCoeff>(&mIDCFourierTransform.getFourierCoefficients(), "TPC/Calib/IDC/FOURIER", mMetadata);
    }

    for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
      const o2::tpc::Side side = iSide ? Side::C : Side::A;
      const header::DataHeader::SubSpecificationType subSpec{iSide};
      output.snapshot(Output{gDataOriginTPC, "IDC0", subSpec, Lifetime::Timeframe}, mIDCs.getIDCZero(side));
      output.snapshot(Output{gDataOriginTPC, "IDC1", subSpec, Lifetime::Timeframe}, mIDCFourierTransform.getIDCOne(side));
      output.snapshot(Output{gDataOriginTPC, "FOURIERREAL", subSpec, Lifetime::Timeframe}, mIDCFourierTransform.getFourierCoefficients(side, FourierCoeff::CoeffType::REAL));
      output.snapshot(Output{gDataOriginTPC, "FOURIERIMAG", subSpec, Lifetime::Timeframe}, mIDCFourierTransform.getFourierCoefficients(side, FourierCoeff::CoeffType::IMAG));
    }

    switch (mCompressionDeltaIDC) {
      case IDCDeltaCompression::MEDIUM: {
        for (int iChunk = 0; iChunk < mIDCs.getNChunks(); ++iChunk) {
          auto idcDeltaMediumCompressed = mIDCs.getIDCDeltaMediumCompressed(iChunk);
          if (mWriteToDB) {
            if (iChunk == 0) {
              mDBapi.storeAsTFileAny<o2::tpc::IDCDeltaCompressionFactors>(&idcDeltaMediumCompressed.getCompressionFactors(), "TPC/Calib/IDC/IDCDELTA/COMPFACTOR", mMetadata);
            }
            mDBapi.storeAsTFileAny<o2::tpc::IDCDeltaContainer<short>>(&idcDeltaMediumCompressed.getIDCDelta(), "TPC/Calib/IDC/IDCDELTA/CONTAINER", mMetadata);
          }
          for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
            const o2::tpc::Side side = iSide ? Side::C : Side::A;
            const header::DataHeader::SubSpecificationType subSpec{iSide};
            output.snapshot(Output{gDataOriginTPC, "IDCDELTA", subSpec, Lifetime::Timeframe}, idcDeltaMediumCompressed.getIDCDelta(side));
            // output.snapshot(Output{gDataOriginTPC, "IDCDELTACOMP", subSpec, Lifetime::Timeframe}, idcDeltaMediumCompressed.getCompressionFactor(side));
          }
        }
        break;
      }
      case IDCDeltaCompression::HIGH: {
        for (int iChunk = 0; iChunk < mIDCs.getNChunks(); ++iChunk) {
          auto idcDeltaHighCompressed = mIDCs.getIDCDeltaHighCompressed(iChunk);
          if (mWriteToDB) {
            if (iChunk == 0) {
              mDBapi.storeAsTFileAny<o2::tpc::IDCDeltaCompressionFactors>(&idcDeltaHighCompressed.getCompressionFactors(), "TPC/Calib/IDC/IDCDELTA/COMPFACTOR", mMetadata);
            }
            mDBapi.storeAsTFileAny<o2::tpc::IDCDeltaContainer<char>>(&idcDeltaHighCompressed.getIDCDelta(), "TPC/Calib/IDC/IDCDELTA/CONTAINER", mMetadata);
          }
          for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
            const o2::tpc::Side side = iSide ? Side::C : Side::A;
            const header::DataHeader::SubSpecificationType subSpec{iSide};
            output.snapshot(Output{gDataOriginTPC, "IDCDELTA", subSpec, Lifetime::Timeframe}, idcDeltaHighCompressed.getIDCDelta(side));
            // output.snapshot(Output{gDataOriginTPC, "IDCDELTACOMP", subSpec, Lifetime::Timeframe}, idcDeltaHighCompressed.getCompressionFactor(side));
          }
        }
        break;
      }
      case IDCDeltaCompression::NO:
      default:
        for (int iChunk = 0; iChunk < mIDCs.getNChunks(); ++iChunk) {
          if (mWriteToDB) {
            mDBapi.storeAsTFileAny<o2::tpc::IDCDeltaContainer<float>>(&mIDCs.getIDCDeltaUncompressed(iChunk).mIDCDelta, "TPC/Calib/IDC/IDCDELTA/CONTAINER", mMetadata);
          }
          for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
            const o2::tpc::Side side = iSide ? Side::C : Side::A;
            const header::DataHeader::SubSpecificationType subSpec{iSide};
            output.snapshot(Output{gDataOriginTPC, "IDCDELTA", subSpec, Lifetime::Timeframe}, mIDCs.getIDCDeltaUncompressed(side, iChunk));
          }
        }
        break;
    }
  }
};

DataProcessorSpec getTPCAggregateGroupedIDCSpec(const std::vector<uint32_t>& crus, const unsigned int timeframes, const unsigned int timeframesDeltaIDC, const unsigned int rangeIDC, const unsigned int fourierCoefficients, const IDCDeltaCompression compression, const bool debug = false)
{
  std::vector<OutputSpec> outputSpecs;
  for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
    const header::DataHeader::SubSpecificationType subSpec{iSide};
    outputSpecs.emplace_back(ConcreteDataMatcher{gDataOriginTPC, "IDC0", subSpec});
    outputSpecs.emplace_back(ConcreteDataMatcher{gDataOriginTPC, "IDC1", subSpec});
    outputSpecs.emplace_back(ConcreteDataMatcher{gDataOriginTPC, "IDCDELTA", subSpec});
    outputSpecs.emplace_back(ConcreteDataMatcher{gDataOriginTPC, "FOURIERREAL", subSpec});
    outputSpecs.emplace_back(ConcreteDataMatcher{gDataOriginTPC, "FOURIERIMAG", subSpec});
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
  // const auto id = "tpc-aggregate-idc";
  return DataProcessorSpec{
    "tpc-aggregate-idc",
    inputSpecs,
    outputSpecs,
    AlgorithmSpec{adaptFromTask<TPCAggregateGroupedIDCSpec>(crus, timeframes, timeframesDeltaIDC, groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, rangeIDC, fourierCoefficients, compression, debug)},
    Options{{"ccdb-uri", VariantType::String, "http://ccdb-test.cern.ch:8080", {"URI for the CCDB access."}}}}; // end DataProcessorSpec
}

} // namespace o2::tpc

#endif
