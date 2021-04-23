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
#include "TPCWorkflow/TPCAverageGroupIDCSpec.h"
#include "TPCWorkflow/PublisherSpec.h"
#include "TPCBase/CRU.h"

using namespace o2::framework;

// customize the completion policy
void customize(std::vector<o2::framework::CompletionPolicy>& policies)
{
  using o2::framework::CompletionPolicy;
  policies.push_back(CompletionPolicyHelpers::defineByName("tpc-idc-averagegroup.*", CompletionPolicy::CompletionOp::Consume));
}

// we need to add workflow options before including Framework/runDataProcessing
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  const std::string cruDefault = "0-" + std::to_string(o2::tpc::CRU::MaxCRU - 1);
  const int defaultlanes = std::max(1u, std::thread::hardware_concurrency() / 2);

  std::vector<ConfigParamSpec> options{
    {"inputFormat", VariantType::String, "Sim", {"setting the input format type: 'Sim'=IDC simulation format, 'Real'=real output format of CRUs"}},
    {"configFile", VariantType::String, "", {"configuration file for configurable parameters"}},
    {"debug", VariantType::Bool, false, {"create debug files"}},
    {"lanes", VariantType::Int, defaultlanes, {"Number of parallel processing lanes."}},
    {"crus", VariantType::String, cruDefault.c_str(), {"List of CRUs, comma separated ranges, e.g. 0-3,7,9-15"}},
    {"groupPads", VariantType::Int, 4, {"number of pads in a row which will be grouped"}},
    {"groupRows", VariantType::Int, 4, {"number of pads in row direction which will be grouped"}},
    {"groupLastRowsThreshold", VariantType::Int, 2, {"set threshold in row direction for merging the last group to the previous group"}},
    {"groupLastPadsThreshold", VariantType::Int, 2, {"set threshold in pad direction for merging the last group to the previous group"}},
  };

  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& config)
{
  using namespace o2::tpc;

  // set up configuration
  o2::conf::ConfigurableParam::updateFromFile(config.options().get<std::string>("configFile"));
  o2::conf::ConfigurableParam::writeINI("o2tpcaveragegroupidc_configuration.ini");

  const auto inputFormatStr = config.options().get<std::string>("inputFormat");
  const TPCIntegrateIDCDevice::IDCFormat inputFormat = inputFormatStr.compare("Sim") ? TPCIntegrateIDCDevice::IDCFormat::Real : TPCIntegrateIDCDevice::IDCFormat::Sim;

  const auto tpcCRUs = o2::RangeTokenizer::tokenize<int>(config.options().get<std::string>("crus"));
  const auto nCRUs = tpcCRUs.size();
  const auto nLanes = std::min(static_cast<unsigned long>(config.options().get<int>("lanes")), nCRUs);
  const auto crusPerLane = nCRUs / nLanes + ((nCRUs % nLanes) != 0);

  const auto groupPads = config.options().get<int>("groupPads");
  const auto groupRows = config.options().get<int>("groupRows");
  const auto groupLastRowsThreshold = config.options().get<int>("groupLastRowsThreshold");
  const auto groupLastPadsThreshold = config.options().get<int>("groupLastPadsThreshold");
  const auto debug = config.options().get<bool>("debug");

  WorkflowSpec workflow;
  if (nLanes <= 0) {
    return workflow;
  }

  for (int ilane = 0; ilane < nLanes; ++ilane) {
    const auto first = tpcCRUs.begin() + ilane * crusPerLane;
    if (first >= tpcCRUs.end()) {
      break;
    }
    const auto last = std::min(tpcCRUs.end(), first + crusPerLane);
    const std::vector<uint32_t> rangeCRUs(first, last);
    workflow.emplace_back(getTPCAverageGroupIDCSpec(ilane, rangeCRUs, inputFormat, groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, debug));
  }

  return workflow;
}
