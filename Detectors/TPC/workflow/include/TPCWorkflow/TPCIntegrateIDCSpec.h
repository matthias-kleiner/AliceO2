// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_CALIBRATION_TPCINTEGRATEIDCSPEC_H
#define O2_CALIBRATION_TPCINTEGRATEIDCSPEC_H

/// @file   TPCIntegrateIDCSpec.h
/// @brief  TPC integration of IDCs processor

#include <vector>
#include <string>
#include <chrono>
#include <fmt/format.h>

#include "Framework/Task.h"
#include "Framework/ControlService.h"
#include "Framework/Logger.h"
#include "Framework/DataProcessorSpec.h"

#include "CommonUtils/MemFileHelper.h"
#include "Headers/DataHeader.h"

#include "TPCBase/CalDet.h"
#include "TPCBase/CDBInterface.h"
#include "TPCBase/CRU.h"
#include "DataFormatsTPC/TPCSectorHeader.h"

#include "TPCCalibration/CalibTreeDump.h"

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
  TPCIntegrateIDCDevice(uint32_t lane, const std::vector<uint32_t>& sectors, uint32_t publishAfterTFs) : mLane{lane}, mSectors(sectors), mPublishAfter(publishAfterTFs) {}

  void init(o2::framework::InitContext& ic) final
  {
    LOGP(info, "init TPCIntegrateIDCDevice");
    mIDCs = CalPad("IDC", PadSubset::ROC);
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    LOGP(info, "run TPCIntegrateIDCDevice");

    for (int i = 0; i < mSectors.size(); ++i) {
      DataRef ref = pc.inputs().getByPos(i);
      auto const* tpcSectorHeader = o2::framework::DataRefUtils::getHeader<o2::tpc::TPCSectorHeader*>(ref);
      auto inDigits = pc.inputs().get<gsl::span<o2::tpc::Digit>>(ref);

      LOG(INFO) << "received " << inDigits.size() << " digits";

      for (const auto& digit : inDigits) {
        const auto timeBin = digit.getTimeStamp();

        const auto row = digit.getRow(); // global pad row
        const auto pad = digit.getPad(); // pad position
        const float charge = digit.getChargeFloat();
        o2::tpc::CRU cru = digit.getCRU();
        const auto sector = cru.sector();
        const float tmpCharge = mIDCs.getValue(sector, row, pad);
        mIDCs.setValue(sector, row, pad, tmpCharge + charge);
      }
    }
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
    LOGP(info, "endOfStream");
    dumpCalibData();
    sendOutput(ec.outputs());
    ec.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  }

 private:
  // const int nTimeBins = 200; ///< 200 zBins TODO FIX THIS
  CalPad mIDCs;
  bool mIDCsSet = false;
  // uint32_t mMaxEvents{0};           ///< maximum number of events to process
  uint32_t mPublishAfter{0};        ///< number of events after which to dump the calibration
  uint32_t mLane{0};                ///< lane number of processor
  std::vector<uint32_t> mSectors{}; ///< sectors to process in this instance
  // bool mReadyToQuit{false};         ///< if processor is ready to quit
  // bool mCalibDumped{false};         ///< if calibration object already dumped
  // bool mForceQuit{false};           ///< for quit after processing finished
  // bool mDirectFileDump{false};      ///< directly dump the calibration data to file

  //____________________________________________________________________________
  void sendOutput(DataAllocator& output)
  {
    LOGP(info, "sendOutput");
    auto image = o2::utils::MemFileHelper::createFileImage(&mIDCs, typeid(mIDCs), mIDCs.getName(), "data");
    CDBType dataType{CDBType::CalPedestal};
    int type = int(dataType);

    header::DataHeader::SubSpecificationType subSpec{(header::DataHeader::SubSpecificationType)((mLane << 4))};
    output.snapshot(Output{gDataOriginTPC, "CLBPART", subSpec}, *image.get());
    output.snapshot(Output{gDataOriginTPC, "CLBPARTINFO", subSpec}, type);
  }

  //____________________________________________________________________________
  void dumpCalibData()
  {
    LOGP(info, "Dumping output");

    CalibTreeDump dump;
    dump.add(&mIDCs);
    dump.dumpToFile(fmt::format("idcs{:02}.root", mLane));
    mIDCsSet = true;

    const std::string name = fmt::format("idcs_obj_{:02}.root", mLane);
    TFile f(name.data(), "recreate");
    f.WriteObject(&mIDCs, mIDCs.getName().data());
    f.Close();
  }
};

DataProcessorSpec getTPCIntegrateIDCSpec(uint32_t ilane = 0, std::vector<uint32_t> sectors = {}, uint32_t publishAfterTFs = 0)
{
  std::vector<o2::framework::OutputSpec> outputs{
    ConcreteDataTypeMatcher{gDataOriginTPC, "CLBPART"},
    ConcreteDataTypeMatcher{gDataOriginTPC, "CLBPARTINFO"}};

  std::vector<InputSpec> inputSpecs;
  inputSpecs.reserve(sectors.size());
  for (auto& sector : sectors) {
    inputSpecs.emplace_back(InputSpec{"digits", gDataOriginTPC, "DIGITS", sector, Lifetime::Timeframe});
  }

  const auto id = fmt::format("calib-tpc-integrateidc-{:02}", ilane);

  return DataProcessorSpec{
    id.data(),
    inputSpecs,
    outputs,
    AlgorithmSpec{adaptFromTask<TPCIntegrateIDCDevice>(ilane, sectors, publishAfterTFs)},
    Options{
      {"max-events", VariantType::Int, 0, {"maximum number of events to process"}},
      {"force-quit", VariantType::Bool, false, {"force quit after max-events have been reached"}},
      {"direct-file-dump", VariantType::Bool, false, {"directly dump calibration to file"}},
    } // end Options
  };  // end DataProcessorSpec
}

} // namespace tpc
} // namespace o2

#endif
