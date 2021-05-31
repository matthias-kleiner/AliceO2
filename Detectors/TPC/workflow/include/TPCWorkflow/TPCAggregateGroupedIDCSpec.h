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
#include <limits>

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
  TPCAggregateGroupedIDCSpec(const std::vector<uint32_t>& crus, const unsigned int timeframes, const unsigned int timeframesDeltaIDC, std::array<unsigned char, Mapper::NREGIONS> groupPads,
                             std::array<unsigned char, Mapper::NREGIONS> groupRows, std::array<unsigned char, Mapper::NREGIONS> groupLastRowsThreshold,
                             std::array<unsigned char, Mapper::NREGIONS> groupLastPadsThreshold, const unsigned int rangeIDC, const IDCDeltaCompression compression, const bool debug = false)
    : mCRUs{crus}, mIDCs{groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, timeframes, timeframesDeltaIDC}, mIDCFourierTransform{rangeIDC, timeframes}, mCompressionDeltaIDC{compression}, mDebug{debug} {};

  void init(o2::framework::InitContext& ic) final
  {
    // mDBapi.init(ic.options().get<std::string>("ccdb-uri")); // or http://localhost:8080 for a local installation
    mDBapi.init("http://localhost:8080"); // or http://localhost:8080 for a local installation
    mWriteToDB = mDBapi.isHostReachable() ? true : false;
    mUpdateGroupingPar = ic.options().get<bool>("update-not-grouping-parameter");
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    // set the TF for first aggregated TF
    if (mProcessedTFs == 0) {
      mTFRange[0] = getCurrentTF(pc);

      if (mWriteToDB && !mUpdateGroupingPar) {
        // write struct containing grouping parameters to access grouped IDCs to CCDB
        const ParameterIDCGroupCCDB parGrouping(mIDCs.getGroupPads(), mIDCs.getGroupRows(), mIDCs.getRowThreshold(), mIDCs.getPadThreshold());
        // validity for grouping parameters is from first TF to some really large TF (until it is updated)
        mDBapi.storeAsTFileAny<o2::tpc::ParameterIDCGroupCCDB>(&parGrouping, "TPC/Calib/IDC/GROUPINGPAR", mMetadata, getFirstTF(), std::numeric_limits<uint32_t>::max());
        mUpdateGroupingPar = true;
      }
    }

    const int nCRUs = mCRUs.size();
    for (int i = 0; i < nCRUs; ++i) {
      const DataRef ref = pc.inputs().getByPos(i);
      auto const* tpcCRUHeader = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
      const int cru = tpcCRUHeader->subSpecification >> 7;
      mIDCs.setIDCs(pc.inputs().get<std::vector<float>>(ref), cru, mProcessedTFs);
    }
    ++mProcessedTFs;

    if (mProcessedTFs == mIDCs.getNTimeframes()) {
      // set the TF for last aggregated TF
      mTFRange[1] = getCurrentTF(pc);
      mProcessedTFs = 0;

      // calculate DeltaIDC, 0D-IDC, 1D-IDC
      mIDCs.factorizeIDCs();

      if (mDebug) {
        // dump this first to also store 1D-IDCs! otherwise they will be moved out and are not available
        LOGP(info, "dumping aggregated IDCS to file");
        mIDCs.dumpToFile(fmt::format("IDCFactorized_{:02}.root", getCurrentTF(pc)).data());
      }

      // perform fourier transform of 0D-IDCs
      mIDCFourierTransform.setIDCs(std::move(mIDCs).getIDCOne(), mIDCs.getIntegrationIntervalsPerTF());
      mIDCFourierTransform.calcFourierCoefficients();

      // storing to CCDB
      sendOutput(pc.outputs());

      if (mDebug) {
        LOGP(info, "dumping fourier transform  to file");
        mIDCFourierTransform.dumpToFile(fmt::format("Fourier_{:02}.root", getCurrentTF(pc)).data());
      }
    }
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
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
  std::array<uint32_t, 2> mTFRange{};               ///< storing of first and last TF used when setting the validity of the objects when w riting to CCDB
  bool mUpdateGroupingPar{true};                    ///< flag to set if grouping parameters should be updated or not

  uint32_t getCurrentTF(o2::framework::ProcessingContext& pc) const
  {
    const DataRef ref = pc.inputs().getByPos(0);
    auto const* tpcCRUHeader = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
    return tpcCRUHeader->tfCounter;
  }

  uint32_t getFirstTF() const { return mTFRange[0]; }
  uint32_t getLastTF() const { return mTFRange[1]; }
  unsigned int getFirstTFDeltaIDC(const unsigned int iChunk) const { return getFirstTF() + iChunk * mIDCs.getTimeFramesDeltaIDC(); }
  unsigned int getLastTFDeltaIDC(const unsigned int iChunk) const { return (iChunk == mIDCs.getNChunks() - 1) ? (mIDCs.getNTimeframes() - 1 + getFirstTF()) : (getFirstTFDeltaIDC(iChunk) + mIDCs.getTimeFramesDeltaIDC() - 1); }

  void sendOutput(DataAllocator& output)
  {
    if (mWriteToDB) {
      // store IDC Zero One in CCDB
      const long timeStampStart = getFirstTF();
      const long timeStampEnd = getLastTF();
      mDBapi.storeAsTFileAny<o2::tpc::IDCZero>(&mIDCs.getIDCZero(), "TPC/Calib/IDC/IDC0", mMetadata, timeStampStart, timeStampEnd);
      mDBapi.storeAsTFileAny<o2::tpc::IDCOne>(&mIDCFourierTransform.getIDCOne(), "TPC/Calib/IDC/IDC1", mMetadata, timeStampStart, timeStampEnd);
      mDBapi.storeAsTFileAny<o2::tpc::FourierCoeff>(&mIDCFourierTransform.getFourierCoefficients(), "TPC/Calib/IDC/FOURIER", mMetadata, timeStampStart, timeStampEnd);

      for (unsigned int iChunk = 0; iChunk < mIDCs.getNChunks(); ++iChunk) {
        switch (mCompressionDeltaIDC) {
          case IDCDeltaCompression::MEDIUM: {
            auto idcDeltaMediumCompressed = mIDCs.getIDCDeltaMediumCompressed(iChunk);
            mDBapi.storeAsTFileAny<o2::tpc::IDCDelta<short>>(&idcDeltaMediumCompressed, "TPC/Calib/IDC/IDCDELTA", mMetadata, getFirstTFDeltaIDC(iChunk), getLastTFDeltaIDC(iChunk));
            break;
          }
          case IDCDeltaCompression::HIGH: {
            auto idcDeltaHighCompressed = mIDCs.getIDCDeltaHighCompressed(iChunk);
            mDBapi.storeAsTFileAny<o2::tpc::IDCDelta<char>>(&idcDeltaHighCompressed, "TPC/Calib/IDC/IDCDELTA", mMetadata, getFirstTFDeltaIDC(iChunk), getLastTFDeltaIDC(iChunk));
            break;
          }
          case IDCDeltaCompression::NO:
          default:
            mDBapi.storeAsTFileAny<o2::tpc::IDCDelta<float>>(&mIDCs.getIDCDeltaUncompressed(iChunk), "TPC/Calib/IDC/IDCDELTA", mMetadata, getFirstTFDeltaIDC(iChunk), getLastTFDeltaIDC(iChunk));
            break;
        }
      }
    }
    // reseting aggregated IDCs. This is done for safety, but if all data is received in the next aggregation interval it isnt necessary... remove it?
    mIDCs.reset();
  }
};

DataProcessorSpec getTPCAggregateGroupedIDCSpec(const std::vector<uint32_t>& crus, const unsigned int timeframes, const unsigned int timeframesDeltaIDC, const unsigned int rangeIDC, const IDCDeltaCompression compression, const bool debug = false)
{
  std::vector<OutputSpec> outputSpecs;
  std::vector<InputSpec> inputSpecs;
  inputSpecs.reserve(crus.size());
  for (const auto& cru : crus) {
    const header::DataHeader::SubSpecificationType subSpec{cru << 7};
    inputSpecs.emplace_back(InputSpec{"idcsgroup", gDataOriginTPC, "IDCGROUP", subSpec, Lifetime::Timeframe});
  }

  const auto& paramIDCGroup = ParameterIDCGroup::Instance();
  std::array<unsigned char, Mapper::NREGIONS> groupPads{};
  std::array<unsigned char, Mapper::NREGIONS> groupRows{};
  std::array<unsigned char, Mapper::NREGIONS> groupLastRowsThreshold{};
  std::array<unsigned char, Mapper::NREGIONS> groupLastPadsThreshold{};
  std::copy(std::begin(paramIDCGroup.GroupPads), std::end(paramIDCGroup.GroupPads), std::begin(groupPads));
  std::copy(std::begin(paramIDCGroup.GroupRows), std::end(paramIDCGroup.GroupRows), std::begin(groupRows));
  std::copy(std::begin(paramIDCGroup.GroupLastRowsThreshold), std::end(paramIDCGroup.GroupLastRowsThreshold), std::begin(groupLastRowsThreshold));
  std::copy(std::begin(paramIDCGroup.GroupLastPadsThreshold), std::end(paramIDCGroup.GroupLastPadsThreshold), std::begin(groupLastPadsThreshold));

  return DataProcessorSpec{
    "tpc-aggregate-idc",
    inputSpecs,
    outputSpecs,
    AlgorithmSpec{adaptFromTask<TPCAggregateGroupedIDCSpec>(crus, timeframes, timeframesDeltaIDC, groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, rangeIDC, compression, debug)},
    Options{{"ccdb-uri", VariantType::String, "http://ccdb-test.cern.ch:8080", {"URI for the CCDB access."}},
            {"update-not-grouping-parameter", VariantType::Bool, false, {"Update/Writing grouping parameters to CCDB."}}}}; // end DataProcessorSpec
}

} // namespace o2::tpc

#endif
