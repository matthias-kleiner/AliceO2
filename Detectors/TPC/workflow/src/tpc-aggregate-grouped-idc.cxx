// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <fmt/format.h>
#include "Algorithm/RangeTokenizer.h"
#include "Framework/WorkflowSpec.h"
#include "Framework/Logger.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/CompletionPolicy.h"
#include "Framework/CompletionPolicyHelpers.h"
#include "Headers/DataHeader.h"
#include "Headers/RAWDataHeader.h"
#include "CommonUtils/ConfigurableParam.h"
#include <vector>
#include <string>
#include "TPCWorkflow/TPCAggregateGroupedIDCSpec.h"
#include "TPCCalibration/IDCFactorization.h"
#include "TPCCalibration/IDCFourierTransform.h"

using namespace o2::framework;

// customize the completion policy
void customize(std::vector<o2::framework::CompletionPolicy>& policies)
{
  using o2::framework::CompletionPolicy;
  policies.push_back(CompletionPolicyHelpers::defineByName("tpc-idc-aggregate.*", CompletionPolicy::CompletionOp::Consume));
}

// we need to add workflow options before including Framework/runDataProcessing
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  const std::string cruDefault = "0-" + std::to_string(o2::tpc::CRU::MaxCRU - 1);

  std::vector<ConfigParamSpec> options{
    {"configFile", VariantType::String, "o2tpcaveragegroupidc_configuration.ini", {"configuration file for configurable parameters"}},
    {"timeframes", VariantType::Int, 10, {"Number of TFs which will be aggregated."}},
    {"timeframesDeltaIDC", VariantType::Int, 2, {"Number of TFs used for storing DeltaIDC in CCDB."}},
    {"rangeIDC", VariantType::Int, 10, {"Number 1D-IDCs which will be used for the calculation of the fourier coefficients."}},
    {"nthreads", VariantType::Int, 1, {"Number of threads which will be used during factorization of the IDCs."}},
    {"debug", VariantType::Bool, false, {"create debug files"}},
    {"fft", VariantType::Bool, false, {"using FFT or naive approach"}},
    {"crus", VariantType::String, cruDefault.c_str(), {"List of CRUs, comma separated ranges, e.g. 0-3,7,9-15"}},
    {"compression", VariantType::Int, 1, {"compression of DeltaIDC: 0 -> No, 1 -> Medium (data compression ratio 2), 2 -> High (data compression ratio ~6)"}},
    {"configKeyValues", VariantType::String, "", {"Semicolon separated key=value strings (e.g.: 'TPCIDCCompressionParam.MaxIDCDeltaValue=0.5;')"}}};

  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& config)
{
  using namespace o2::tpc;

  // set up configuration
  o2::conf::ConfigurableParam::updateFromFile(config.options().get<std::string>("configFile"));
  o2::conf::ConfigurableParam::updateFromString(config.options().get<std::string>("configKeyValues"));
  o2::conf::ConfigurableParam::writeINI("o2tpcaggregateidc_configuration.ini");

  const auto tpcCRUs = o2::RangeTokenizer::tokenize<int>(config.options().get<std::string>("crus"));
  const auto nCRUs = tpcCRUs.size();
  const auto timeframes = static_cast<unsigned int>(config.options().get<int>("timeframes"));
  const auto timeframesDeltaIDC = static_cast<unsigned int>(config.options().get<int>("timeframesDeltaIDC"));
  const auto debug = config.options().get<bool>("debug");
  const auto fft = config.options().get<bool>("fft");
  const auto rangeIDC = static_cast<unsigned int>(config.options().get<int>("rangeIDC"));
  const auto nthreads = static_cast<unsigned long>(config.options().get<int>("nthreads"));
  IDCFactorization::setNThreads(nthreads);
  IDCFourierTransform::setFFT(fft);

  const int compressionTmp = config.options().get<int>("compression");
  IDCDeltaCompression compression;
  switch (compressionTmp) {
    case static_cast<int>(IDCDeltaCompression::NO):
    case static_cast<int>(IDCDeltaCompression::MEDIUM):
    case static_cast<int>(IDCDeltaCompression::HIGH):
      compression = static_cast<IDCDeltaCompression>(compressionTmp);
      break;
    default:
      LOGP(error, "wrong compression type set. Setting compression to medium compression");
      compression = static_cast<IDCDeltaCompression>(IDCDeltaCompression::MEDIUM);
      break;
  }

  WorkflowSpec workflow;
  const auto first = tpcCRUs.begin();
  const auto last = std::min(tpcCRUs.end(), first + nCRUs);
  const std::vector<uint32_t> rangeCRUs(first, last);
  workflow.emplace_back(getTPCAggregateGroupedIDCSpec(rangeCRUs, timeframes, timeframesDeltaIDC, rangeIDC, compression, debug));

  return workflow;
}
