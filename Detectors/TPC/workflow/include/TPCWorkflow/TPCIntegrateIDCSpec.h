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
#include "DataFormatsTPC/Digit.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsTPC/Constants.h"

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
  // enum for the output format of the integrated IDCs
  enum class IDCFormat : int {
    Sim = 0, // output format of simulation for faster processing
    Real = 1 // output format of real CRUs
  };

  TPCIntegrateIDCDevice(const uint32_t lane, const std::vector<uint32_t>& sectors, const uint32_t nOrbitsPerIDCIntervall, const IDCFormat outputFormat, const bool debug) : mLane{lane}, mSectors{sectors}, mNOrbits{nOrbitsPerIDCIntervall}, mIDCFormat{outputFormat}, mDebug{debug} {}

  void run(o2::framework::ProcessingContext& pc) final
  {
    const int nSectors = mSectors.size();
    for (int i = 0; i < nSectors; ++i) {
      const DataRef ref = pc.inputs().getByPos(i);
      auto const* tpcSectorHeader = o2::framework::DataRefUtils::getHeader<o2::tpc::TPCSectorHeader*>(ref);
      const int sector = tpcSectorHeader->sector();

      // integrate digits for given sector
      const gsl::span<const o2::tpc::Digit> digits = pc.inputs().get<gsl::span<o2::tpc::Digit>>(ref);
      LOG(INFO) << "received " << digits.size() << " digits";
      resetIDCs();

      // loop over digits from one sector for one Time Frame
      for (const auto& digit : digits) {
        const auto timeBin = digit.getTimeStamp();
        const o2::tpc::CRU cru(digit.getCRU());
        const unsigned int region = cru.region();
        const unsigned int indexPadOffs = static_cast<int>(timeBin / mTimeStamps) * mPadsPerRegion[region]; // offset of the pad index for current integration interval
        const auto row = digit.getRow();                                                                    // global pad row
        const auto pad = digit.getPad();                                                                    // pad in row
        const unsigned int indexPad = indexPadOffs + getPadIndex(row, pad);                                 // pad index for row and pad
        const float charge = digit.getChargeFloat();
        mIDCs[region][indexPad] += charge;
      }

      // send the output for one sector for one TF
      sendOutput(pc.outputs(), sector);
      if (mDebug) {
        dumpIDCs(sector);
      }
    }
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
    LOGP(info, "endOfStream");
    ec.services().get<ControlService>().readyToQuit(QuitRequest::Me);
    if (mDebug) {
      createDebugTree();
    }
  }

  /// get the total number of regions in the TPC
  static int getNRegions() { return mRegions; } ///< returns total number of regions for one sector

  /// get pad offset to calculate global pad number from pad number in cru
  /// \param region region of the tpc
  static int getGlobalPadOffset(const int region) { return mGlobalPadOffs[region]; }

  /// get the nummber of pads for given region
  /// \param region region of the tpc
  static int getPadsPerRegion(const int region) { return mPadsPerRegion[region]; }

  /// get number of intergration intervals per TF
  /// \param idcs vector containg the idcs.
  /// \param region region of the idcs
  static int getIntegrationIntervalls(const std::vector<float>& idcs, const int region)
  {
    return idcs.size() / mPadsPerRegion[region];
  }

  /// return the kind of the output for given type.
  /// \param idcFormat type of the IDC format
  static header::DataDescription getDataDescription(const IDCFormat idcFormat)
  {
    return (idcFormat == IDCFormat::Sim) ? header::DataDescription{"IDCSIM"} : header::DataDescription{"IDC"};
  }

  static int getNIntegrationIntervals(const int nDigitsCRU, const int region)
  {
    return nDigitsCRU / mPadsPerRegion[region];
  }

 private:
  static constexpr int mRegions{10};                                                                                               ///< number of regions per sector
  static constexpr int mPadRows{152};                                                                                              ///< total number of pad rows
  static constexpr int mPadsPerRegion[mRegions]{1200, 1200, 1440, 1440, 1440, 1440, 1600, 1600, 1600, 1600};                       ///< number of pads per CRU
  static constexpr int mGlobalPadOffs[mRegions]{0, 1200, 2400, 3840, 5280, 6720, 8160, 9760, 11360, 12960};                        ///< offset of number of pads for region used for debugging only
  static constexpr int mOffs[mPadRows]{                                                                                            ///< row offset in cru for given global pad row
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
  const uint32_t mLane{};                                                                                                          ///< lane number of processor
  const std::vector<uint32_t> mSectors{};                                                                                          ///< sectors to process in this instance
  static constexpr const uint32_t mLengthOfTF{256};                                                                                ///< length of one TF in units of orbits
  const uint32_t mNOrbits{12};                                                                                                     ///< integration of IDCs in units of orbits
  const uint32_t mTimeStamps{o2::constants::lhc::LHCMaxBunches / o2::tpc::constants::LHCBCPERTIMEBIN * mNOrbits};                  ///< number of time stamps for each integration interval
  const uint32_t mIntegrationIntervalsPerTF{mLengthOfTF / mNOrbits + 1};                                                           ///< number of integration intervals per TF. Add 1: 256/12=21.333
  const IDCFormat mIDCFormat{IDCFormat::Sim};                                                                                      ///< type of the output format. Sim=simulation, Real=realistic format
  const bool mDebug{false};                                                                                                        ///< dump IDCs to tree for debugging
  std::array<std::vector<float>, mRegions> mIDCs{                                                                               ///< IDCs for one sector. The index of the array to the region.
                                                    std::vector<float>(mPadsPerRegion[0] * mIntegrationIntervalsPerTF),            // region 0
                                                    std::vector<float>(mPadsPerRegion[1] * mIntegrationIntervalsPerTF),            // region 1
                                                    std::vector<float>(mPadsPerRegion[2] * mIntegrationIntervalsPerTF),            // region 2
                                                    std::vector<float>(mPadsPerRegion[3] * mIntegrationIntervalsPerTF),            // region 3
                                                    std::vector<float>(mPadsPerRegion[4] * mIntegrationIntervalsPerTF),            // region 4
                                                    std::vector<float>(mPadsPerRegion[5] * mIntegrationIntervalsPerTF),            // region 5
                                                    std::vector<float>(mPadsPerRegion[6] * mIntegrationIntervalsPerTF),            // region 6
                                                    std::vector<float>(mPadsPerRegion[7] * mIntegrationIntervalsPerTF),            // region 7
                                                    std::vector<float>(mPadsPerRegion[8] * mIntegrationIntervalsPerTF),            // region 8
                                                    std::vector<float>(mPadsPerRegion[9] * mIntegrationIntervalsPerTF)};           // region 9

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
    for (auto& idcs : mIDCs) {
      std::fill(idcs.begin(), idcs.end(), 0);
    }
  }

  void sendOutput(DataAllocator& output, const int sector)
  {
    unsigned int cru = sector * mRegions;
    for (const auto& idcs : mIDCs) {
      if (mIDCFormat == IDCFormat::Sim) {
        // const TPCCRUHeader cruheader{cru, mIntegrationIntervalsPerTF};
        const header::DataHeader::SubSpecificationType subSpec{cru << 7};
        // output.snapshot(Output{gDataOriginTPC, getDataDescription(mIDCFormat), subSpec, Lifetime::Timeframe, cruheader}, idcs);
        output.snapshot(Output{gDataOriginTPC, getDataDescription(mIDCFormat), subSpec, Lifetime::Timeframe}, idcs);

      } else {
        // TODO
        // convert to format from thorsten here
        // send.......
        // DUMMY FOR NOW
        // const TPCCRUHeader cruheader{cru, mIntegrationIntervalsPerTF};
        const header::DataHeader::SubSpecificationType subSpec{cru << 7};
        // output.snapshot(Output{gDataOriginTPC, getDataDescription(mIDCFormat), subSpec, Lifetime::Timeframe, cruheader}, idcs);
        output.snapshot(Output{gDataOriginTPC, getDataDescription(mIDCFormat), subSpec, Lifetime::Timeframe}, idcs);
      }
      ++cru;
    }
  }

  void dumpIDCs(const int sector)
  {
    const std::string name = fmt::format("idcs_obj_{:02}.root", sector);
    TFile fOut(name.data(), "RECREATE");
    int cru = sector * mRegions;
    for (const auto& idcs : mIDCs) {
      fOut.WriteObject(&idcs, Form("cru_%i", cru));
      ++cru;
    }
    fOut.Close();
  }

  void createDebugTree()
  {
    const static Mapper& mapper = Mapper::instance();

    const std::string nameTree = fmt::format("idcs_tree_{:02}.root", mLane);
    o2::utils::TreeStreamRedirector pcstream(nameTree.data(), "RECREATE");
    pcstream.GetFile()->cd();

    for (int iSec = 0; iSec < mSectors.size(); ++iSec) {
      const std::string nameFile = fmt::format("idcs_obj_{:02}.root", mSectors[iSec]);
      TFile fObj(nameFile.data(), "READ");
      TIter next(fObj.GetListOfKeys());
      TKey* key = nullptr;
      while ((key = (TKey*)next())) {
        const std::string nameObj = key->GetName();
        std::vector<float>* idcs = (std::vector<float>*)fObj.Get(nameObj.data());           // IDCs for all time bins
        const int cru = std::stoi(nameObj.substr(nameObj.find("_") + 1, nameObj.length())); // extract cru from object name

        const o2::tpc::CRU cruTmp(cru);
        unsigned int region = cruTmp.region();
        const int padsPerCRU = mPadsPerRegion[region];
        std::vector<int> vRow(padsPerCRU);
        std::vector<int> vPad(padsPerCRU);
        std::vector<float> vXPos(padsPerCRU);
        std::vector<float> vYPos(padsPerCRU);
        std::vector<float> vGlobalXPos(padsPerCRU);
        std::vector<float> vGlobalYPos(padsPerCRU);
        std::vector<float> idcsPerTimeBin(padsPerCRU); // idcs for one time bin

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

        for (int iTimeBin = 0; iTimeBin < mIntegrationIntervalsPerTF; ++iTimeBin) {
          for (int iPad = 0; iPad < padsPerCRU; ++iPad) {
            idcsPerTimeBin[iPad] = (*idcs)[iPad + iTimeBin * mPadsPerRegion[region]];
          }

          pcstream << "tree"
                   << "cru=" << cruiTmp
                   << "sector=" << sectorTmp
                   << "region=" << region
                   << "timeBin=" << iTimeBin
                   << "IDCs.=" << idcsPerTimeBin
                   << "pad.=" << vPad
                   << "row.=" << vRow
                   << "lx.=" << vXPos
                   << "ly.=" << vYPos
                   << "gx.=" << vGlobalXPos
                   << "gy.=" << vGlobalYPos
                   << "\n";
        }
        delete idcs;
      }
      delete key;
      fObj.Close();
    }
    pcstream.Close();
  }
};

DataProcessorSpec getTPCIntegrateIDCSpec(const uint32_t ilane = 0, const std::vector<uint32_t>& sectors = {}, const uint32_t nOrbits = 22, const TPCIntegrateIDCDevice::IDCFormat outputFormat = TPCIntegrateIDCDevice::IDCFormat::Sim, const bool debug = false)
{
  std::vector<InputSpec> inputSpecs;
  inputSpecs.reserve(sectors.size());

  std::vector<OutputSpec> outputSpecs;
  const int nRegions = TPCIntegrateIDCDevice::getNRegions();
  outputSpecs.reserve(sectors.size() * nRegions);

  // define input and output specs
  for (const auto& sector : sectors) {
    inputSpecs.emplace_back(InputSpec{"digits", gDataOriginTPC, "DIGITS", sector, Lifetime::Timeframe});

    // output spec
    unsigned int cru = sector * nRegions;
    for (int iRegion = 0; iRegion < nRegions; ++iRegion) {
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
