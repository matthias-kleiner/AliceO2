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
#include "TPCBase/Mapper.h"
#include "Framework/Variant.h"
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

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

  // std::string groupPads{"4,4,4,4,4,4,4,4,4,4"};
  // std::string groupRows{"4,4,4,4,4,4,4,4,4,4"};
  // std::string groupLastRowsThreshold{"2,2,2,2,2,2,2,2,2,2"};
  // std::string groupLastPadsThreshold{"2,2,2,2,2,2,2,2,2,2"};

  std::vector<ConfigParamSpec> options{
    {"inputFormat", VariantType::String, "Sim", {"setting the input format type: 'Sim'=IDC simulation format, 'Real'=real output format of CRUs"}},
    {"configFile", VariantType::String, "", {"configuration file for configurable parameters"}},
    {"debug", VariantType::Bool, false, {"create debug files"}},
    {"lanes", VariantType::Int, defaultlanes, {"Number of parallel processing lanes."}},
    {"crus", VariantType::String, cruDefault.c_str(), {"List of CRUs, comma separated ranges, e.g. 0-3,7,9-15"}},
    {"groupPads", VariantType::String, "4,4,4,4,4,4,4,4,4,4", {"number of pads in a row which will be grouped per region"}},
    {"groupRows", VariantType::String, "4,4,4,4,4,4,4,4,4,4", {"number of pads in row direction which will be grouped per region"}},
    {"groupLastRowsThreshold", VariantType::String, "2,2,2,2,2,2,2,2,2,2", {"set threshold in row direction for merging the last group to the previous group per region"}},
    {"groupLastPadsThreshold", VariantType::String, "2,2,2,2,2,2,2,2,2,2", {"set threshold in pad direction for merging the last group to the previous group per region"}}};

  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& config)
{
  using namespace o2::tpc;

  const auto inputFormatStr = config.options().get<std::string>("inputFormat");
  const TPCIntegrateIDCDevice::IDCFormat inputFormat = inputFormatStr.compare("Sim") ? TPCIntegrateIDCDevice::IDCFormat::Real : TPCIntegrateIDCDevice::IDCFormat::Sim;

  const auto tpcCRUs = o2::RangeTokenizer::tokenize<int>(config.options().get<std::string>("crus"));
  const auto nCRUs = tpcCRUs.size();
  const auto nLanes = std::min(static_cast<unsigned long>(config.options().get<int>("lanes")), nCRUs);
  const auto crusPerLane = nCRUs / nLanes + ((nCRUs % nLanes) != 0);

  const std::string sgroupPads = config.options().get<std::string>("groupPads");
  const std::string sgroupRows = config.options().get<std::string>("groupRows");
  const std::string sgroupLastRowsThreshold = config.options().get<std::string>("groupLastRowsThreshold");
  const std::string sgroupLastPadsThreshold = config.options().get<std::string>("groupLastPadsThreshold");

  const boost::char_separator<char> sep(","); /// char separator for the tokenizer
  const boost::tokenizer<boost::char_separator<char>> tgroupPads(sgroupPads, sep);
  const boost::tokenizer<boost::char_separator<char>> tgroupRows(sgroupRows, sep);
  const boost::tokenizer<boost::char_separator<char>> tgroupLastRowsThreshold(sgroupLastRowsThreshold, sep);
  const boost::tokenizer<boost::char_separator<char>> tgroupLastPadsThreshold(sgroupLastPadsThreshold, sep);

  std::vector<unsigned int> vgroupPads;
  std::vector<unsigned int> vgroupRows;
  std::vector<unsigned int> vgroupLastRowsThreshold;
  std::vector<unsigned int> vgroupLastPadsThreshold;
  vgroupPads.reserve(Mapper::NREGIONS);
  vgroupRows.reserve(Mapper::NREGIONS);
  vgroupLastRowsThreshold.reserve(Mapper::NREGIONS);
  vgroupLastPadsThreshold.reserve(Mapper::NREGIONS);

  std::transform(tgroupPads.begin(), tgroupPads.end(), std::back_inserter(vgroupPads), &boost::lexical_cast<int, std::string>);
  std::transform(tgroupRows.begin(), tgroupRows.end(), std::back_inserter(vgroupRows), &boost::lexical_cast<int, std::string>);
  std::transform(tgroupLastRowsThreshold.begin(), tgroupLastRowsThreshold.end(), std::back_inserter(vgroupLastRowsThreshold), &boost::lexical_cast<int, std::string>);
  std::transform(tgroupLastPadsThreshold.begin(), tgroupLastPadsThreshold.end(), std::back_inserter(vgroupLastPadsThreshold), &boost::lexical_cast<int, std::string>);

  if (vgroupPads.size() == 1) {
    vgroupPads = std::vector<unsigned int>(Mapper::NREGIONS, vgroupPads.front());
  } else if (vgroupPads.size() != Mapper::NREGIONS) {
    LOGP(error, "wrong number of parameters inserted for groupPads (n={}). Number should be 1 or {}", vgroupPads.size(), Mapper::NREGIONS);
  }

  if (vgroupRows.size() == 1) {
    vgroupRows = std::vector<unsigned int>(Mapper::NREGIONS, vgroupRows.front());
  } else if (vgroupRows.size() != Mapper::NREGIONS) {
    LOGP(error, "wrong number of parameters inserted for groupRows (n={}). Number should be 1 or {}", vgroupRows.size(), Mapper::NREGIONS);
  }

  if (vgroupLastRowsThreshold.size() == 1) {
    vgroupLastRowsThreshold = std::vector<unsigned int>(Mapper::NREGIONS, vgroupLastRowsThreshold.front());
  } else if (vgroupLastRowsThreshold.size() != Mapper::NREGIONS) {
    LOGP(error, "wrong number of parameters inserted for groupLastRowsThreshold (n={}). Number should be 1 or {}", vgroupLastRowsThreshold.size(), Mapper::NREGIONS);
  }

  if (vgroupLastPadsThreshold.size() == 1) {
    vgroupLastPadsThreshold = std::vector<unsigned int>(Mapper::NREGIONS, vgroupLastPadsThreshold.front());
  } else if (vgroupLastPadsThreshold.size() != Mapper::NREGIONS) {
    LOGP(error, "wrong number of parameters inserted for groupLastPadsThreshold (n={}). Number should be 1 or {}", vgroupLastPadsThreshold.size(), Mapper::NREGIONS);
  }

  for (int i = 0; i < vgroupPads.size(); ++i) {
    LOGP(info, "vgroupPads {} {}", i, vgroupPads[i]);
  }

  for (int i = 0; i < Mapper::NREGIONS; ++i) {
    o2::conf::ConfigurableParam::setValue<unsigned int>("TPCIDCGroupParam", fmt::format("GroupPads[{}]", i).data(), vgroupPads[i]);
    o2::conf::ConfigurableParam::setValue<unsigned int>("TPCIDCGroupParam", fmt::format("GroupRows[{}]", i).data(), vgroupRows[i]);
    o2::conf::ConfigurableParam::setValue<unsigned int>("TPCIDCGroupParam", fmt::format("GroupLastRowsThreshold[{}]", i).data(), vgroupLastRowsThreshold[i]);
    o2::conf::ConfigurableParam::setValue<unsigned int>("TPCIDCGroupParam", fmt::format("GroupLastPadsThreshold[{}]", i).data(), vgroupLastPadsThreshold[i]);
  }

  // set up configuration
  o2::conf::ConfigurableParam::updateFromFile(config.options().get<std::string>("configFile"));
  o2::conf::ConfigurableParam::writeINI("o2tpcaveragegroupidc_configuration.ini");

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
    workflow.emplace_back(getTPCAverageGroupIDCSpec(ilane, rangeCRUs, inputFormat, debug));
  }

  return workflow;
}
