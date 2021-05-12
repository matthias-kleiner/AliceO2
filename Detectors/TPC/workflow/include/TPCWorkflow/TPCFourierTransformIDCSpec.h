// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TPCFourierTransformIDCSpec.h
/// \brief TPC fourier transform of 1D-IDCs device
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>
/// \date May 10, 2021

#ifndef O2_TPCFOURIERTRANSFORMIDCSPEC_H
#define O2_TPCFOURIERTRANSFORMIDCSPEC_H

#include <vector>
#include <fmt/format.h>

#include "Framework/Task.h"
#include "Framework/ControlService.h"
#include "Framework/Logger.h"
#include "Framework/DataProcessorSpec.h"
#include "CommonUtils/MemFileHelper.h"
#include "Headers/DataHeader.h"
#include "DataFormatsTPC/Defs.h"
// #include "TPCBase/CDBInterface.h"
#include "CCDB/CcdbApi.h"
// #include "TPCCalibration/IDCFactorization.h"
// #include "CCDB/CcdbApi.h"
#include "Framework/ConfigParamRegistry.h"
#include "TPCCalibration/IDCFourierTransform.h"

using namespace o2::framework;
using o2::header::gDataOriginTPC;
using namespace o2::tpc;

namespace o2
{
namespace tpc
{

class TPCFourierTransformIDCSpec : public o2::framework::Task
{
 public:
  TPCFourierTransformIDCSpec(const unsigned int rangeIntegrationIntervals, const unsigned int shift, const unsigned int nFourierCoefficients, const bool debug)
    : mIDCFourierTransform{rangeIntegrationIntervals, shift, nFourierCoefficients}, mDebug{debug} {};

  void init(o2::framework::InitContext& ic) final
  {
    mDBapi.init(ic.options().get<std::string>("ccdb-uri")); // or http://localhost:8080 for a local installation
    mWriteToDB = false;                                     //mDBapi.isHostReachable() ? true : false;
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
      const o2::tpc::Side side = iSide ? Side::C : Side::A;
      const DataRef ref = pc.inputs().getByPos(iSide);
      mIDCFourierTransform.setIDCs(pc.inputs().get<std::vector<float>>(ref), side);
      mIDCFourierTransform.calcFourierCoefficients(side);
    }

    mIDCFourierTransform.getFourierCoefficients();

    if (mDebug) {
      const DataRef ref = pc.inputs().getByPos(0);
      auto const* header = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
      const auto tf = header->tfCounter;
      mIDCFourierTransform.dumpToFile(Form("Fourier_%i.root", tf));
    }
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
    LOGP(info, "endOfStream");
    ec.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  }

 private:
  IDCFourierTransform mIDCFourierTransform{};
  const bool mDebug{false};                     ///< dump IDCs to tree for debugging
  o2::ccdb::CcdbApi mDBapi;                     ///< object for storing the IDCs at CCDB
  std::map<std::string, std::string> mMetadata; ///< meta data of the stored object in CCDB
  bool mWriteToDB{};                            ///< flag if writing to CCDB will be done

  void sendOutput(DataAllocator& output)
  {
    if (mWriteToDB) {
      // store IDC Zero One in CCDB
      mDBapi.storeAsTFileAny(&mIDCFourierTransform.getFourierCoefficients(), "TPC/Calib/IDC", mMetadata);
    }

    for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
      const o2::tpc::Side side = iSide ? Side::C : Side::A;
      const header::DataHeader::SubSpecificationType subSpec{iSide};
      // output.snapshot(Output{gDataOriginTPC, "IDCFOURIERCOEFF", subSpec, Lifetime::Timeframe}, mIDCFourierTransform.getFourierCoefficients(side));
    }
  }
};

DataProcessorSpec getTPCFourierTransformIDCSpec(const unsigned int rangeIntegrationIntervals, const unsigned int shift, const unsigned int nFourierCoefficients, const bool debug = false)
{
  std::vector<OutputSpec> outputSpecs;
  for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
    const header::DataHeader::SubSpecificationType subSpec{iSide};
    outputSpecs.emplace_back(ConcreteDataMatcher{gDataOriginTPC, "IDCCOEFFREAL", subSpec});
    outputSpecs.emplace_back(ConcreteDataMatcher{gDataOriginTPC, "IDCCOEFFIMAG", subSpec});
  }

  std::vector<InputSpec> inputSpecs;
  inputSpecs.reserve(o2::tpc::SIDES);
  for (unsigned int iSide = 0; iSide < o2::tpc::SIDES; ++iSide) {
    const header::DataHeader::SubSpecificationType subSpec{iSide};
    inputSpecs.emplace_back(InputSpec{"idcone", gDataOriginTPC, "IDCONE", subSpec, Lifetime::Timeframe});
  }

  const auto id = fmt::format("tpc-fourier-IDC");
  return DataProcessorSpec{
    id.data(),
    inputSpecs,
    outputSpecs,
    AlgorithmSpec{adaptFromTask<TPCFourierTransformIDCSpec>(rangeIntegrationIntervals, shift, nFourierCoefficients, debug)},
    Options{{"ccdb-uri", VariantType::String, "http://ccdb-test.cern.ch:8080", {"URI for the CCDB access."}}}}; // end DataProcessorSpec
}

} // namespace tpc
} // namespace o2

#endif
