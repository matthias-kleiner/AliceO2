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
#include <fmt/format.h>

#include "Framework/Task.h"
#include "Framework/ControlService.h"
#include "Framework/Logger.h"
#include "Framework/DataProcessorSpec.h"
#include "CommonUtils/MemFileHelper.h"
#include "Headers/DataHeader.h"
#include "TPCBase/CDBInterface.h"
#include "TPCBase/CRU.h"
#include "DataFormatsTPC/TPCSectorHeader.h"

#include "CommonUtils/TreeStreamRedirector.h" // for debugging

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
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    LOGP(info, "run TPCIntegrateIDCDevice");

    const int nSectors = mSectors.size();
    for (int i = 0; i < nSectors; ++i) {
      const DataRef ref = pc.inputs().getByPos(i);
      auto const* tpcSectorHeader = o2::framework::DataRefUtils::getHeader<o2::tpc::TPCSectorHeader*>(ref);
      const int sector = tpcSectorHeader->sector();

      // set CRU
      const int cruOff = sector * mCRUS;
      for (int iRegion = 0; iRegion < mCRUS; ++iRegion) {
        mIDCs[iRegion].first = cruOff + iRegion;
      }

      // reset vectors if necesseary
      if (i > 0) {
        resetIDCs();
      }

      // integrate digits for given sector
      const gsl::span<const o2::tpc::Digit> inDigits = pc.inputs().get<gsl::span<o2::tpc::Digit>>(ref);
      processSector(inDigits, pc);
    }
  }

  void processSector(const gsl::span<const o2::tpc::Digit>& inDigits, o2::framework::ProcessingContext& pc)
  {
    LOG(INFO) << "received " << inDigits.size() << " digits";

    unsigned int timeBinStop = mTimeStamps;
    for (const auto& digit : inDigits) {
      auto timeBin = digit.getTimeStamp();

      // sendoutput if maximum number of timebins is reached
      if (timeBin > timeBinStop) {
        sendOutput(pc.outputs());
        if (mDebug) {
          // dump integrated digits to a tree
          const o2::tpc::CRU cru(digit.getCRU());
          dumpIDCs(cru.sector(), (timeBinStop - mTimeStamps) / mTimeStamps);
        }
        resetIDCs(); // set IDCs
        timeBinStop += mTimeStamps;
      }

      const auto row = digit.getRow(); // global pad row
      const auto pad = digit.getPad(); // pad in row
      const float charge = digit.getChargeFloat();
      const int indexPad = getPadIndex(row, pad); // pad index for row and pad
      const o2::tpc::CRU cru(digit.getCRU());
      unsigned int region = cru.region();
      mIDCs[region].second[indexPad] += charge;
    }
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
    LOGP(info, "endOfStream");

    // dumpIDCs();
    sendOutput(ec.outputs());
    ec.services().get<ControlService>().readyToQuit(QuitRequest::Me);

    if (mDebug) {
      createDebugTree();
    }
  }

 private:
  static const int mCRUS{10};                                                                                           ///< number of CRUs per sector
  static const int mPadRows{152};                                                                                       ///< total number of pad rows
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
  std::array<std::pair<unsigned int, std::vector<float>>, mCRUS> mIDCs{
    std::make_pair(0, std::vector<float>(mPadsPerCRU[0])),
    std::make_pair(1, std::vector<float>(mPadsPerCRU[1])),
    std::make_pair(2, std::vector<float>(mPadsPerCRU[2])),
    std::make_pair(3, std::vector<float>(mPadsPerCRU[3])),
    std::make_pair(4, std::vector<float>(mPadsPerCRU[4])),
    std::make_pair(5, std::vector<float>(mPadsPerCRU[5])),
    std::make_pair(6, std::vector<float>(mPadsPerCRU[6])),
    std::make_pair(7, std::vector<float>(mPadsPerCRU[7])),
    std::make_pair(8, std::vector<float>(mPadsPerCRU[8])),
    std::make_pair(9, std::vector<float>(mPadsPerCRU[9]))}; ///< IDCs for one sector. The key value corresponds to the CRU. The index of the array to the partition.
  const uint32_t mLane{0};                                  ///< lane number of processor
  const std::vector<uint32_t> mSectors{};                   ///< sectors to process in this instance
  const uint32_t mTimeStamps{2000};                         ///< number of time stamps for each integration interval
  const bool mDebug{true};                                  ///< dump IDCs to tree for debugging

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
    for (const auto& [cru, value] : mIDCs) {
      const header::DataHeader::SubSpecificationType subSpec{cru << 7};
      output.snapshot(Output{gDataOriginTPC, "IDC", subSpec}, value);
    }
  }

  //____________________________________________________________________________
  void dumpIDCs(const int sector, const int timeBin)
  {
    const std::string name = fmt::format("idcs_obj_{:02}.root", sector);
    TFile f(name.data(), "UPDATE");
    for (const auto& [cru, value] : mIDCs) {
      f.WriteObject(&value, Form("cru_%i_%i", cru, timeBin));
    }
    f.Close();
  }

  void createDebugTree()
  {
    const static Mapper& mapper = Mapper::instance();

    const std::string nameTree = fmt::format("idcs_tree_{:02}.root", mLane);
    o2::utils::TreeStreamRedirector pcstream(nameTree.data(), "RECREATE");
    pcstream.GetFile()->cd();

    for (int iSec = 0; iSec < mSectors.size(); ++iSec) {
      const std::string nameObj = fmt::format("idcs_obj_{:02}.root", mSectors[iSec]);
      TFile fObj(nameObj.data(), "READ");

      TIter next(fObj.GetListOfKeys());
      TKey* key = nullptr;
      while ((key = (TKey*)next())) {
        const std::string nameObj = key->GetName();
        // std::cout << Form("key: %s points to an object of class: %s at %lld", nameObj.data(), key->GetClassName(), key->GetSeekKey()) << std::endl;
        const std::string delim = "_";
        unsigned first_delim_pos = nameObj.find(delim);
        unsigned end_pos_of_first_delim = first_delim_pos + delim.length();
        unsigned last_delim_pos = nameObj.find_first_of(delim, end_pos_of_first_delim);
        const int cru = std::stoi(nameObj.substr(end_pos_of_first_delim, last_delim_pos - end_pos_of_first_delim));
        int timeBin = std::stoi(nameObj.substr(last_delim_pos + delim.length(), nameObj.length()));

        std::vector<float>* idcs = (std::vector<float>*)fObj.Get(nameObj.data());

        const o2::tpc::CRU cruTmp(cru);
        unsigned int region = cruTmp.region();
        const int padsPerCRU = mPadsPerCRU[region];
        std::vector<int> vRow(padsPerCRU);
        std::vector<int> vPad(padsPerCRU);
        std::vector<float> vXPos(padsPerCRU);
        std::vector<float> vYPos(padsPerCRU);
        std::vector<float> vGlobalXPos(padsPerCRU);
        std::vector<float> vGlobalYPos(padsPerCRU);
        int sectorTmp = cruTmp.sector();

        for (int iPad = 0; iPad < padsPerCRU; ++iPad) {
          const GlobalPadNumber globalNum = mGlobalPadOffs[region] + iPad;
          const auto& padPosLocal = mapper.padPos(globalNum);
          vRow[iPad] = padPosLocal.getRow();
          vPad[iPad] = padPosLocal.getPad();
          vXPos[iPad] = mapper.getPadCentre(padPosLocal).X();
          vYPos[iPad] = mapper.getPadCentre(padPosLocal).Y();

          const GlobalPosition2D globalPos = mapper.LocalToGlobal(LocalPosition2D(vXPos[iPad], vYPos[iPad]), cruTmp.sector());
          vGlobalXPos[iPad] = globalPos.X();
          vGlobalYPos[iPad] = globalPos.Y();
        }

        int cruiTmp = cru;

        pcstream << "tree"
                 << "cru=" << cruiTmp
                 << "sector=" << sectorTmp
                 << "region=" << region
                 << "timeBin=" << timeBin
                 << "value.=" << *idcs
                 << "pad.=" << vPad
                 << "row.=" << vRow
                 << "lx.=" << vXPos
                 << "ly.=" << vYPos
                 << "gx.=" << vGlobalXPos
                 << "gy.=" << vGlobalYPos
                 << "\n";
        delete idcs;
      }
      fObj.Close();
    }
    pcstream.Close();
  }
};

DataProcessorSpec getTPCIntegrateIDCSpec(uint32_t ilane = 0, std::vector<uint32_t> sectors = {}, uint32_t nTimeBins = 2000)
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
    AlgorithmSpec{adaptFromTask<TPCIntegrateIDCDevice>(ilane, sectors, nTimeBins)},
    Options{
      // TODO SET THESE
      {"max-events", VariantType::Int, 0, {"maximum number of events to process"}},
      {"force-quit", VariantType::Bool, false, {"force quit after max-events have been reached"}},
      {"direct-file-dump", VariantType::Bool, false, {"directly dump calibration to file"}},
    } // end Options
  };  // end DataProcessorSpec
}

} // namespace tpc
} // namespace o2

#endif
