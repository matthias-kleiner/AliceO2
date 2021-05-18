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
#include "TPCWorkflow/TPCFourierTransformIDCSpec.h"
#include "TPCWorkflow/PublisherSpec.h"
#include "TPCCalibration/IDCFactorization.h"

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
    {"range", VariantType::Int, 200, {"Number of integration intervals (IDC values) used for the calculation of one fourier coefficient."}},
    {"shiftrange", VariantType::Int, 10, {"Number of integration intervals of which the range will be shifted after the calculation of one fourier coefficient."}},
    {"fourierCoeff", VariantType::Int, 10, {"Number of fourier coefficients which will be calculated for each interval."}},
    {"debug", VariantType::Bool, false, {"create debug files"}}};

  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& config)
{
  using namespace o2::tpc;

  // set up configuration
  o2::conf::ConfigurableParam::updateFromFile(config.options().get<std::string>("configFile"));
  o2::conf::ConfigurableParam::writeINI("o2tpcfouriertransformidc_configuration.ini");

  const auto rangeIntegrationIntervals = static_cast<unsigned int>(config.options().get<int>("range"));
  const auto shift = static_cast<unsigned int>(config.options().get<int>("shiftrange"));
  const auto nFourierCoefficients = static_cast<unsigned int>(config.options().get<int>("fourierCoeff"));
  const auto debug = config.options().get<bool>("debug");

  WorkflowSpec workflow;
  workflow.emplace_back(getTPCFourierTransformIDCSpec(rangeIntegrationIntervals, shift, nFourierCoefficients, debug));

  return workflow;
}
