// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_CALIBRATION_TPCAverageMergeIDCSpec_H
#define O2_CALIBRATION_TPCAverageMergeIDCSpec_H

/// @file   TPCAverageMergeIDCSpec.h
/// @brief  TPC merging and averaging of IDCs

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
#include "TPCWorkflow/TPCIntegrateIDCSpec.h"

#include "CommonUtils/TreeStreamRedirector.h" // for debugging
#include "TPCBase/Mapper.h"

using namespace o2::framework;
using o2::header::gDataOriginTPC;
using namespace o2::tpc;

namespace o2
{
namespace tpc
{

class TPCAverageMergeIDCDevice : public o2::framework::Task
{
 public:
  TPCAverageMergeIDCDevice(const uint32_t lane, const std::vector<uint32_t>& crus, const TPCIntegrateIDCDevice::IDCFormat inputFormat) : mLane{lane}, mCRUs{crus}, mInputFormat{inputFormat} {}

  void init(o2::framework::InitContext& ic) final
  {
    const static auto& mapper = Mapper::instance();

    if (mDebug) {
      const std::string name = fmt::format("SECOND_idcs_obj_{:02}.root", mLane);
      TFile fOut(name.data(), "RECREATE");
      fOut.Close();
    }

    for (int iCRU = 0; iCRU < mCRUs.size(); ++iCRU) {
      const int numberOfPads = 333; // total number of pads per CRU after averaging and merging
      mIDCs[iCRU] = std::vector<float>(numberOfPads);
    }
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    const int nCRUs = mCRUs.size();
    LOG(INFO) << "Lane: " << mLane;
    LOG(INFO) << "nCRUs: " << nCRUs;
    for (int i = 0; i < nCRUs; ++i) {
      const DataRef ref = pc.inputs().getByPos(i);
      std::vector<float> inIDCs = pc.inputs().get<std::vector<float>>(ref);
      unsigned int nDigits = inIDCs.size();
      auto const* tpcCRUHeader = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
      const int cru = tpcCRUHeader->subSpecification >> 7;
      const o2::tpc::CRU cruTmp(cru);
      unsigned int region = cruTmp.region();

      const auto nIntegrationIntervals = TPCIntegrateIDCDevice::getNIntegrationIntervals(nDigits, region);

      std::vector<std::vector<float>> mIDCs{};

      LOG(INFO) << "received " << nDigits << " IDCs";
      LOG(INFO) << "cru " << cru;
      LOG(INFO) << "nIntegrationIntervals " << nIntegrationIntervals;

      const int nTotalPads = TPCIntegrateIDCDevice::getPadsPerRegion(region); // total number of pads for current region
      const int facRows = mRowsPerRegion[region] / mMergeNRows;
      // mMergeNPads
      // mIDCs.push_back();

      dumpIDCs(cru, inIDCs);

      // for(auto& IDC : inIDCs){
      // LOG(INFO) << "received " << IDC;
      // }
    }
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
    LOGP(info, "endOfStream");
    ec.services().get<ControlService>().readyToQuit(QuitRequest::Me);
    createDebugTree();
    // using namespace std::chrono_literals;
    // std::this_thread::sleep_for(50s);
  }

 private:
  static constexpr int mRegions{10};
  const uint32_t mLane{0};                                                               ///< lane number of processor
  const std::vector<uint32_t> mCRUs{};                                                   ///< CRUs to process in this instance
  const bool mDebug{true};                                                               ///< dump IDCs to tree for debugging
  const TPCIntegrateIDCDevice::IDCFormat mInputFormat{};                                 ///< type of the input format. Sim=simulation, Real=realistic format
  const int mMergeNRows{2};                                                              ///< number of rows which will be merged
  const int mMergeNPads{2};                                                              ///< number of pads which will be merged
  std::vector<std::vector<float>> mIDCs{mCRUs.size()};                                   ///< outer vector crus, inner vector 'pads per CRU * integration intervals'
  static constexpr int mRowsPerRegion[mRegions]{17, 15, 16, 15, 18, 16, 16, 14, 13, 12}; ///< number of pad rows per region

  void sendOutput(DataAllocator& output)
  {
    // for (const auto& [cru, value] : mIDCs) {
    //   const header::DataHeader::SubSpecificationType subSpec{cru << 7};
    //   output.snapshot(Output{gDataOriginTPC, "IDC", subSpec}, value);
    // }
  }

  void dumpIDCs(const int cru, const std::vector<float>& value)
  {
    const std::string name = fmt::format("SECOND_idcs_obj_{:02}.root", mLane);
    TFile fOut(name.data(), "UPDATE");
    fOut.WriteObject(&value, Form("cru_%i", cru));
    fOut.Close();
  }

  void createDebugTree()
  {
    const static Mapper& mapper = Mapper::instance();

    const std::string nameTree = fmt::format("SECOND_idcs_tree_{:02}.root", mLane);
    o2::utils::TreeStreamRedirector pcstream(nameTree.data(), "RECREATE");
    pcstream.GetFile()->cd();

    const std::string nameFile = fmt::format("SECOND_idcs_obj_{:02}.root", mLane);
    TFile fObj(nameFile.data(), "READ");
    TIter next(fObj.GetListOfKeys());
    TKey* key = nullptr;
    while ((key = (TKey*)next())) {
      const std::string nameObj = key->GetName();
      std::vector<float>* idcs = (std::vector<float>*)fObj.Get(nameObj.data());           // IDCs for all time bins
      const int cru = std::stoi(nameObj.substr(nameObj.find("_") + 1, nameObj.length())); // extract cru from object name

      const o2::tpc::CRU cruTmp(cru);
      unsigned int region = cruTmp.region();
      const int padsPerCRU = TPCIntegrateIDCDevice::getPadsPerRegion(region); //mPadsPerCRU[region];
      std::vector<int> vRow(padsPerCRU);
      std::vector<int> vPad(padsPerCRU);
      std::vector<float> vXPos(padsPerCRU);
      std::vector<float> vYPos(padsPerCRU);
      std::vector<float> vGlobalXPos(padsPerCRU);
      std::vector<float> vGlobalYPos(padsPerCRU);
      std::vector<float> idcsPerTimeBin(padsPerCRU); // idcs for one time bin

      int sectorTmp = cruTmp.sector();

      for (int iPad = 0; iPad < padsPerCRU; ++iPad) {
        const GlobalPadNumber globalNum = TPCIntegrateIDCDevice::getGlobalPadOffset(region) + iPad; //mGlobalPadOffs[region] + iPad;
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

      for (int iTimeBin = 0; iTimeBin < TPCIntegrateIDCDevice::getIntegrationIntervalls(*idcs, region); ++iTimeBin) {
        for (int iPad = 0; iPad < padsPerCRU; ++iPad) {
          idcsPerTimeBin[iPad] = (*idcs)[iPad + iTimeBin * TPCIntegrateIDCDevice::getPadsPerRegion(region)];
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
    pcstream.Close();
  }
};

DataProcessorSpec getTPCAverageMergeIDCSpec(const uint32_t ilane = 0, const std::vector<uint32_t>& crus = {}, const TPCIntegrateIDCDevice::IDCFormat inputFormat = {})
{
  std::vector<OutputSpec> outputSpecs;
  outputSpecs.reserve(crus.size());

  std::vector<InputSpec> inputSpecs;
  inputSpecs.reserve(crus.size());

  for (const auto& cru : crus) {
    const header::DataHeader::SubSpecificationType subSpec{cru << 7};
    inputSpecs.emplace_back(InputSpec{"idcs", gDataOriginTPC, TPCIntegrateIDCDevice::getDataDescription(inputFormat), subSpec, Lifetime::Timeframe});
    outputSpecs.emplace_back(ConcreteDataMatcher{gDataOriginTPC, "IDCAVERAGE", subSpec});
  }

  const auto id = fmt::format("calib-tpc-averagemerge-{:02}", ilane);
  return DataProcessorSpec{
    id.data(),
    inputSpecs,
    outputSpecs,
    AlgorithmSpec{adaptFromTask<TPCAverageMergeIDCDevice>(ilane, crus, inputFormat)},
  }; // end DataProcessorSpec
}

} // namespace tpc
} // namespace o2

#endif
