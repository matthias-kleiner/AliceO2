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
#include "TPCWorkflow/TPCIntegrateIDCSpec.h"
#include "TPCWorkflow/PublisherSpec.h"

using namespace o2::framework;

// customize the completion policy
void customize(std::vector<o2::framework::CompletionPolicy>& policies)
{
  using o2::framework::CompletionPolicy;
  policies.push_back(CompletionPolicyHelpers::defineByName("calib-tpc-integrateidc.*", CompletionPolicy::CompletionOp::Consume));
}

// we need to add workflow options before including Framework/runDataProcessing
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  std::string sectorDefault = "0-" + std::to_string(o2::tpc::Sector::MAXSECTOR - 1);
  int defaultlanes = std::max(1u, std::thread::hardware_concurrency() / 2);

  std::vector<ConfigParamSpec> options{
    {"nOrbits", VariantType::Int, 12, {"number of orbits for which the IDCs are integrated"}},
    {"outputFormat", VariantType::String, "Sim", {"setting the output type: 'Sim'=IDC simulation format, 'Real'=real output format of CRUs"}},
    {"configFile", VariantType::String, "", {"configuration file for configurable parameters"}},
    {"lanes", VariantType::Int, defaultlanes, {"Number of parallel processing lanes."}},
    {"sectors", VariantType::String, sectorDefault.c_str(), {"List of TPC sectors, comma separated ranges, e.g. 0-3,7,9-15"}},
  };

  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& config)
{
  using namespace o2::tpc;

  // set up configuration
  o2::conf::ConfigurableParam::updateFromFile(config.options().get<std::string>("configFile"));
  o2::conf::ConfigurableParam::writeINI("o2tpcidc_configuration.ini");

  const auto nOrbits = (uint32_t)config.options().get<int>("nOrbits");
  const auto outputFormatStr = config.options().get<std::string>("outputFormat");
  const TPCIntegrateIDCDevice::OutputFormat outputFormat = outputFormatStr.compare("Sim") ? TPCIntegrateIDCDevice::OutputFormat::Real : TPCIntegrateIDCDevice::OutputFormat::Sim;

  const auto tpcsectors = o2::RangeTokenizer::tokenize<int>(config.options().get<std::string>("sectors"));
  const auto nSectors = (uint32_t)tpcsectors.size();
  const auto nLanes = std::min((uint32_t)config.options().get<int>("lanes"), nSectors);
  const auto sectorsPerLane = nSectors / nLanes + ((nSectors % nLanes) != 0);

  WorkflowSpec workflow;
  if (nLanes <= 0) {
    return workflow;
  }

  for (int ilane = 0; ilane < nLanes; ++ilane) {
    const auto first = tpcsectors.begin() + ilane * sectorsPerLane;
    if (first >= tpcsectors.end()) {
      break;
    }
    const auto last = std::min(tpcsectors.end(), first + sectorsPerLane);
    const std::vector<uint32_t> range(first, last);
    workflow.emplace_back(getTPCIntegrateIDCSpec(ilane, range, nOrbits, outputFormat));
  }

  return workflow;
}
