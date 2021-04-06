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
#include "CommonUtils/TreeStreamRedirector.h"

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
  TPCIntegrateIDCDevice(const uint32_t lane, const std::vector<uint32_t>& sectors, const uint32_t nTimeStamps) : mLane{lane}, mSectors(sectors), mTimeStamps(nTimeStamps) {}

  void init(o2::framework::InitContext& ic) final
  {
    LOGP(info, "init TPCIntegrateIDCDevice");
    // initialize map
    for (int iSec = 0; iSec < mSectors.size(); ++iSec) {
      const int cruOff = mSectors[iSec] * mCRUS;
      for (int iCRU = 0; iCRU < mCRUS; ++iCRU) {
        mIDCs.emplace(cruOff + iCRU, std::vector<float>(mPadsPerCRU[iCRU]));
      }
    }
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    LOGP(info, "run TPCIntegrateIDCDevice");
    ++mEvents;

    for (int i = 0; i < mSectors.size(); ++i) {
      DataRef ref = pc.inputs().getByPos(i);
      auto const* tpcSectorHeader = o2::framework::DataRefUtils::getHeader<o2::tpc::TPCSectorHeader*>(ref);
      auto inDigits = pc.inputs().get<gsl::span<o2::tpc::Digit>>(ref);
      const int nDigits = inDigits.size();

      LOG(INFO) << "received " << nDigits << " digits";
      LOG(INFO) << "sectors " << mSectors.size() << " digits";

      for (const auto& digit : inDigits) {
        auto timeBin = digit.getTimeStamp();
        // send every mNTimeBins the integrated IDCs
        if (timeBin % mTimeStamps == 0 && i == mSectors.size() - 1) {
          sendOutput(pc.outputs());
          // resetIDCs(); // set IDCs back to 0
        }
        const auto row = digit.getRow(); // global pad row
        const auto pad = digit.getPad(); // pad in row
        const float charge = digit.getChargeFloat();
        const int cru = digit.getCRU();
        int indexPad = getPadIndex(row, pad); // pad index for row and pad
        mIDCs[cru][indexPad] += charge;
      }
    }
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
    LOGP(info, "endOfStream");
    LOGP(info, "number of events processed {}", mEvents);

    dumpCalibData();
    sendOutput(ec.outputs());
    ec.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  }

 private:
  std::unordered_map<unsigned int, std::vector<float>> mIDCs;                                                           ///< key CRU value is vectors of IDCs per CRU
  static const int mCRUS = 10;                                                                                          ///< number of CRUs per sector
  static const int mPadRows = 152;                                                                                      ///< total number of pad rows
  const int mPadsPerCRU[mCRUS]{1200, 1200, 1440, 1440, 1440, 1440, 1600, 1600, 1600, 1600};                             ///< number of pads per CRU
  const int mGlobalPadOffs[mCRUS]{0, 1200, 2400, 3840, 5280, 6720, 8160, 9760, 11360, 12960};                           ///< offset of number of pads for region used for debugging only
  const int mOffs[mPadRows]{                                                                                            ///< row offset in cru for given global pad row
                            0, 66, 132, 198, 266, 334, 402, 472, 542, 612, 684, 756, 828, 902, 976, 1050, 1124,         // region 0
                            0, 76, 152, 228, 306, 384, 462, 542, 622, 702, 784, 866, 948, 1032, 1116,                   // region 1
                            0, 86, 172, 258, 346, 434, 522, 612, 702, 792, 882, 974, 1066, 1158, 1252, 1346,            // region 2
                            0, 92, 184, 276, 370, 464, 558, 654, 750, 846, 944, 1042, 1140, 1240, 1340,                 // region 3
                            0, 76, 152, 228, 304, 382, 460, 538, 618, 698, 778, 858, 940, 1022, 1104, 1188, 1272, 1356, // region 4
                            0, 86, 172, 258, 346, 434, 522, 612, 702, 792, 882, 974, 1066, 1158, 1252, 1346,            // region 5
                            0, 94, 190, 286, 382, 480, 578, 676, 776, 876, 978, 1080, 1182, 1286, 1390, 1494,           // region 6
                            0, 110, 220, 332, 444, 556, 670, 784, 898, 1014, 1130, 1246, 1364, 1482,                    // region 7
                            0, 118, 236, 356, 476, 598, 720, 844, 968, 1092, 1218, 1344, 1472,                          // region 8
                            0, 128, 258, 388, 520, 652, 784, 918, 1052, 1188, 1324, 1462};                              // region 9
  // bool mIDCsSet = false;
  // uint32_t mPublishAfter{0};        ///< number of events after which to dump the calibration
  const uint32_t mLane{0};                ///< lane number of processor
  const std::vector<uint32_t> mSectors{}; ///< sectors to process in this instance
  int mEvents{0};                         ///< number of processed events
  const uint32_t mTimeStamps{2000};       ///< number of time stamps for each integration interval
  // bool mReadyToQuit{false};         ///< if processor is ready to quit
  // bool mCalibDumped{false};         ///< if calibration object already dumped
  // bool mForceQuit{false};           ///< for quit after processing finished
  // bool mDirectFileDump{false};      ///< directly dump the calibration data to file

  /// \param row global pad row
  /// \param pad pad in row
  /// \return returns pad number for region
  unsigned int getPadIndex(const int row, const int pad) const
  {
    return mOffs[row] + pad;
  }

  /// set all IDC values to 0
  void resetIDCs()
  {
    for (auto& [cru, value] : mIDCs) {
      std::fill(value.begin(), value.end(), 0);
    }
  }

  //____________________________________________________________________________
  void sendOutput(DataAllocator& output)
  {
    LOGP(info, "sendOutput");
    for (const auto& [cru, value] : mIDCs) {
      header::DataHeader::SubSpecificationType subSpec{cru << 7};
      output.snapshot(Output{gDataOriginTPC, "IDC", subSpec}, value);
    }
  }

  //____________________________________________________________________________
  void dumpCalibData()
  {
    LOGP(info, "Dumping output");
    const static Mapper& mapper = Mapper::instance();

    const std::string name = fmt::format("idcs_obj_{:02}.root", mLane);
    TFile f(name.data(), "recreate");
    for (const auto& [cru, value] : mIDCs) {
      f.WriteObject(&value, Form("cru_%i", cru));
    }
    f.Close();

    o2::utils::TreeStreamRedirector pcstream(name.data(), "UPDATE");
    pcstream.GetFile()->cd();

    for (const auto& [cru, value] : mIDCs) {
      // loop over pads in CRU
      const o2::tpc::CRU cruTmp(cru);
      unsigned int region = cruTmp.region();
      const int padsPerCRU = mPadsPerCRU[region];
      std::vector<int> vRow(padsPerCRU);
      std::vector<int> vPad(padsPerCRU);
      std::vector<float> vXPos(padsPerCRU);
      std::vector<float> vYPos(padsPerCRU);

      for (int iPad = 0; iPad < padsPerCRU; ++iPad) {
        const GlobalPadNumber globalNum = mGlobalPadOffs[region] + iPad;
        const auto& padPos = mapper.padPos(globalNum);
        vRow[iPad] = padPos.getRow();
        vPad[iPad] = padPos.getPad();
        vXPos[iPad] = mapper.getPadCentre(padPos).X();
        vYPos[iPad] = mapper.getPadCentre(padPos).Y();
      }

      int cruiTmp = cru;
      int sector = cruTmp.sector();
      std::vector<float> valueTmp = value;
      pcstream << "tree"
               << "cru=" << cruiTmp
               << "sector=" << sector
               << "region=" << region
               << "value=" << valueTmp
               << "pad=" << vPad
               << "row=" << vRow
               << "x=" << vXPos
               << "y=" << vYPos
               << "\n";
    }
    pcstream.Close();
  }
};

DataProcessorSpec getTPCIntegrateIDCSpec(uint32_t ilane = 0, std::vector<uint32_t> sectors = {}, uint32_t publishAfterTFs = 0)
{
  std::vector<o2::framework::OutputSpec> outputs{ConcreteDataTypeMatcher{gDataOriginTPC, "IDC"}};

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
