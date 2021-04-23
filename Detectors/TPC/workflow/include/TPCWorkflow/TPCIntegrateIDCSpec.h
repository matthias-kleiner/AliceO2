// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TPCIntegrateIDCSpec.h
/// \brief TPC integration of IDCs processor
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>
/// \date Apr 16, 2021

#ifndef O2_TPCINTEGRATEIDCSPEC_H
#define O2_TPCINTEGRATEIDCSPEC_H

#include <vector>
#include <fmt/format.h>

#include "Framework/Task.h"
#include "Framework/ControlService.h"
#include "Framework/Logger.h"
#include "Framework/DataProcessorSpec.h"
#include "CommonUtils/MemFileHelper.h"
#include "Headers/DataHeader.h"
#include "TPCSimulation/IDCSim.h"
#include "DataFormatsTPC/TPCSectorHeader.h"
#include "DataFormatsTPC/Digit.h"
#include "TPCBase/IDCHelper.h"

using namespace o2::framework;
using o2::header::gDataOriginTPC;
using namespace o2::tpc;

namespace o2
{
namespace tpc
{

class TPCIntegrateIDCDevice : public o2::framework::Task
{
 public:
  // enum for the output format of the integrated IDCs
  enum class IDCFormat : int {
    Sim = 0, // output format of simulation for faster processing
    Real = 1 // output format of real CRUs
  };

  TPCIntegrateIDCDevice(const int lane, const std::vector<unsigned int>& sectors, const int nOrbitsPerIDCIntervall, const IDCFormat outputFormat, const bool debug) : mLane{lane}, mSectors{sectors}, mIDCFormat{outputFormat}, mDebug{debug}
  {
    for (const auto& sector : mSectors) {
      mIDCs.emplace(sector, IDCSim(sector, nOrbitsPerIDCIntervall));
    }
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    // loop over sectors
    for (int iSec = 0; iSec < mSectors.size(); ++iSec) {
      const DataRef ref = pc.inputs().getByPos(iSec);
      auto const* tpcSectorHeader = o2::framework::DataRefUtils::getHeader<o2::tpc::TPCSectorHeader*>(ref);
      const int sector = tpcSectorHeader->sector();

      // integrate digits for given sector
      const gsl::span<const o2::tpc::Digit> digits = pc.inputs().get<gsl::span<o2::tpc::Digit>>(ref);
      mIDCs[sector].integrateDigitsForOneTF(digits);

      if (mDebug) {
        // mIDCs[sector].dumpIDCs(0);
        mIDCs[sector].createDebugTree(0);
      }

      // send the output for one sector for one TF
      sendOutput(pc.outputs(), sector);
    }
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
    LOGP(info, "endOfStream");
    ec.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  }

  /// return the kind of the output for given type.
  /// \param idcFormat type of the IDC format
  static header::DataDescription getDataDescription(const IDCFormat idcFormat)
  {
    return (idcFormat == IDCFormat::Sim) ? header::DataDescription{"IDCSIM"} : header::DataDescription{"IDC"};
  }

 private:
  const int mLane{};                                ///< lane number of processor
  const std::vector<unsigned int> mSectors{};       ///< sectors to process in this instance
  const IDCFormat mIDCFormat{IDCFormat::Sim};       ///< type of the output format. Sim=simulation, Real=realistic format
  const bool mDebug{false};                         ///< dump IDCs to tree for debugging
  std::unordered_map<unsigned int, IDCSim> mIDCs{}; ///< integrated IDCs for one TF for all specified sectors

  // send output for one sector
  void sendOutput(DataAllocator& output, const int sector)
  {
    uint32_t cru = sector * o2::tpc::IDCHelper::NREGIONS;
    for (const auto& idcs : mIDCs[sector].get()) {
      if (mIDCFormat == IDCFormat::Sim) {
        const header::DataHeader::SubSpecificationType subSpec{cru << 7};
        output.snapshot(Output{gDataOriginTPC, getDataDescription(mIDCFormat), subSpec, Lifetime::Timeframe}, idcs);

      } else {
        // TODO
        // convert to format from thorsten here
        // send.......
        // DUMMY FOR NOW
        // const TPCCRUHeader cruheader{cru, mIntegrationIntervalsPerTF};
        const header::DataHeader::SubSpecificationType subSpec{cru << 7};
        output.snapshot(Output{gDataOriginTPC, getDataDescription(mIDCFormat), subSpec, Lifetime::Timeframe}, idcs);
      }
      ++cru;
    }
  }
};

DataProcessorSpec getTPCIntegrateIDCSpec(const int ilane = 0, const std::vector<unsigned int>& sectors = {}, const int nOrbits = 12, const TPCIntegrateIDCDevice::IDCFormat outputFormat = TPCIntegrateIDCDevice::IDCFormat::Sim, const bool debug = false)
{
  // TPCIntegrateIDCDevice(const int lane, const std::vector<int>& sectors, const int nOrbitsPerIDCIntervall, const IDCFormat outputFormat, const bool debug) : mLane{lane}, mSectors{sectors}, mIDCFormat{outputFormat}, mDebug{debug}

  std::vector<InputSpec> inputSpecs;
  inputSpecs.reserve(sectors.size());

  std::vector<OutputSpec> outputSpecs;
  outputSpecs.reserve(sectors.size() * IDCHelper::NREGIONS);

  // define input and output specs
  for (const auto& sector : sectors) {
    inputSpecs.emplace_back(InputSpec{"digits", gDataOriginTPC, "DIGITS", sector, Lifetime::Timeframe});

    // output spec
    unsigned int cru = sector * IDCHelper::NREGIONS;
    for (int iRegion = 0; iRegion < IDCHelper::NREGIONS; ++iRegion) {
      const header::DataHeader::SubSpecificationType subSpec{cru << 7};
      outputSpecs.emplace_back(ConcreteDataMatcher{gDataOriginTPC, TPCIntegrateIDCDevice::getDataDescription(outputFormat), subSpec});
      ++cru;
    }
  }

  const auto id = fmt::format("tpc-integrate-idc-{:02}", ilane);
  return DataProcessorSpec{
    id.data(),
    inputSpecs,
    outputSpecs,
    AlgorithmSpec{adaptFromTask<TPCIntegrateIDCDevice>(ilane, sectors, nOrbits, outputFormat, debug)}}; // end DataProcessorSpec
}

} // namespace tpc
} // namespace o2

#endif
