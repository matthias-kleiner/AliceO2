// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TPCPressureTemperatureSpec.cxx
/// \brief device for providing pressure and temperature values
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>
/// \date Jun 4, 2025

#include "TPCWorkflow/TPCPressureTemperatureSpec.h"
#include "Framework/Task.h"
#include "Framework/Logger.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/CCDBParamSpec.h"
#include "TPCBase/CDBInterface.h"
#include "Framework/ConfigParamRegistry.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "CommonDataFormat/Pair.h"
#include "DataFormatsTPC/DCS.h"
#include "CommonUtils/TreeStreamRedirector.h"

using namespace o2::framework;

namespace o2
{
namespace tpc
{

class PressureTemperatureDevice : public o2::framework::Task
{
 public:
  PressureTemperatureDevice(std::shared_ptr<o2::base::GRPGeomRequest> req) : mCCDBRequest(req) {};
  void init(o2::framework::InitContext& ic) final
  {
    o2::base::GRPGeomHelper::instance().setRequest(mCCDBRequest);
    const int intInterval = ic.options().get<int>("fit-interval");
    mFitIntervalMS = intInterval * 1000;
    const bool enableDebugTree = ic.options().get<bool>("enable-root-output");
    if (enableDebugTree) {
      mStreamer = std::make_unique<o2::utils::TreeStreamRedirector>("pt.root", "recreate");
    }
  };

  void endOfStream(EndOfStreamContext& eos) final
  {
    if (mStreamer) {
      mStreamer->Close();
    }
  }

  /// \brief interpolate input values for given timestamp
  /// \param timestamps time stamps of the data
  /// \param values data points
  /// \param timestamp time where to interpolate the values
  float interpolate(const std::vector<uint64_t>& timestamps, const std::vector<float>& values, uint64_t timestamp)
  {
    if (auto idxClosest = o2::math_utils::findClosestIndices(timestamps, timestamp)) {
      auto [idxLeft, idxRight] = *idxClosest;
      if (idxRight > idxLeft) {
        const uint64_t x0 = timestamps[idxLeft];
        const uint64_t x1 = timestamps[idxRight];
        const float y0 = values[idxLeft];
        const float y1 = values[idxRight];
        const float y = (y0 * (x1 - timestamp) + y1 * (timestamp - x0)) / (x1 - x0);
        return y;
      } else {
        return values[idxLeft];
      }
    }
    return 0; // this should never happen
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    o2::base::GRPGeomHelper::instance().checkUpdates(pc);
    pc.inputs().get<dcs::Pressure*>("pressure");
    pc.inputs().get<dcs::Temperature*>("temperature");
    const auto orbitResetTimeMS = o2::base::GRPGeomHelper::instance().getOrbitResetTimeMS();
    const auto firstTFOrbit = pc.services().get<o2::framework::TimingInfo>().firstTForbit;
    const uint64_t timestamp = orbitResetTimeMS + firstTFOrbit * o2::constants::lhc::LHCOrbitMUS * 0.001;

    // find closest temperature and pressure
    const float pressure = interpolate(mPressure.second, mPressure.first, timestamp);
    const float tempA = interpolate(mTemperatureA.second, mTemperatureA.first, timestamp);
    const float tempC = interpolate(mTemperatureC.second, mTemperatureC.first, timestamp);

    if (mStreamer) {
      (*mStreamer) << "pt"
                   << "pressure=" << pressure
                   << "temperatureA=" << tempA
                   << "temperatureC=" << tempC
                   << "time=" << timestamp
                   << "\n";
    }

    LOGP(info, "Sending pressure {}, temperature A {} and temperature C {} for timestamp {}", pressure, tempA, tempC, timestamp);
    pc.outputs().snapshot(Output{o2::header::gDataOriginTPC, o2::tpc::getDataDescriptionTemperature()}, dataformats::Pair<float, float>{tempA, tempC});
    pc.outputs().snapshot(Output{o2::header::gDataOriginTPC, o2::tpc::getDataDescriptionPressure()}, pressure);
  }

  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) final
  {
    o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj);
    if (matcher == ConcreteDataMatcher(o2::header::gDataOriginTPC, "PRESSURECCDB", 0)) {
      LOGP(info, "Updating pressure");
      const auto& pressure = ((dcs::Pressure*)obj);
      mPressure.second = pressure->robustPressure.time;
      mPressure.first = pressure->robustPressure.robustPressure;

      if (mStreamer) {
        (*mStreamer) << "p"
                     << "pressureCCDB=" << pressure
                     << "pressure=" << mPressure
                     << "\n";
      }
    }

    if (matcher == ConcreteDataMatcher(o2::header::gDataOriginTPC, "TEMPERATURECCDB", 0)) {
      LOGP(info, "Updating temperature");
      auto temp = *(dcs::Temperature*)obj;
      temp.fitTemperature(o2::tpc::Side::A, mFitIntervalMS, false);
      temp.fitTemperature(o2::tpc::Side::C, mFitIntervalMS, false);

      mTemperatureA.first.clear();
      mTemperatureC.first.clear();
      mTemperatureA.second.clear();
      mTemperatureC.second.clear();

      for (const auto& dp : temp.statsA.data) {
        mTemperatureA.first.emplace_back(dp.value.mean);
        mTemperatureA.second.emplace_back(dp.time);
      }

      for (const auto& dp : temp.statsC.data) {
        mTemperatureC.first.emplace_back(dp.value.mean);
        mTemperatureC.second.emplace_back(dp.time);
      }

      if (mStreamer) {
        (*mStreamer) << "t"
                     << "temperatureCCDB=" << temp
                     << "temperatureA=" << mTemperatureA
                     << "temperatureC=" << mTemperatureC
                     << "\n";
      }
    }
  }

 private:
  std::shared_ptr<o2::base::GRPGeomRequest> mCCDBRequest;             ///< info for CCDB request
  std::pair<std::vector<float>, std::vector<uint64_t>> mPressure;     ///< pressure values for both measurements
  std::pair<std::vector<float>, std::vector<uint64_t>> mTemperatureA; ///< temperature values A-side
  std::pair<std::vector<float>, std::vector<uint64_t>> mTemperatureC; ///< temperature values C-side
  std::unique_ptr<o2::utils::TreeStreamRedirector> mStreamer;         ///< debug streamer
  int mFitIntervalMS{5 * 60 * 1000};                                  ///< fit interval for the temperature
};

o2::framework::DataProcessorSpec getTPCPressureTemperatureSpec()
{
  std::vector<InputSpec> inputs;
  std::vector<OutputSpec> outputs;
  o2::header::DataDescription dataDescription;

  inputs.emplace_back("pressure", o2::header::gDataOriginTPC, "PRESSURECCDB", 0, Lifetime::Condition, ccdbParamSpec(o2::tpc::CDBTypeMap.at(o2::tpc::CDBType::CalPressure), {}, 1));          // time-dependent
  inputs.emplace_back("temperature", o2::header::gDataOriginTPC, "TEMPERATURECCDB", 0, Lifetime::Condition, ccdbParamSpec(o2::tpc::CDBTypeMap.at(o2::tpc::CDBType::CalTemperature), {}, 1)); // time-dependent

  outputs.emplace_back(o2::header::gDataOriginTPC, o2::tpc::getDataDescriptionPressure(), 0, Lifetime::Timeframe);
  outputs.emplace_back(o2::header::gDataOriginTPC, o2::tpc::getDataDescriptionTemperature(), 0, Lifetime::Timeframe);

  auto ccdbRequest = std::make_shared<o2::base::GRPGeomRequest>(true,                           // orbitResetTime
                                                                false,                          // GRPECS=true for nHBF per TF
                                                                false,                          // GRPLHCIF
                                                                false,                          // GRPMagField
                                                                false,                          // askMatLUT
                                                                o2::base::GRPGeomRequest::None, // geometry
                                                                inputs);
  return DataProcessorSpec{
    "tpc-pressure-temperature",
    inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<PressureTemperatureDevice>(ccdbRequest)},
    Options{
      {"enable-root-output", VariantType::Bool, false, {"Enable root-files output writers"}},
      {"fit-interval", VariantType::Int, 300, {"interval in seconds for which to e.g. perform fits of the temperature sensors"}}} // end Options
  }; // end DataProcessorSpec
}

} // end namespace tpc
} // end namespace o2
