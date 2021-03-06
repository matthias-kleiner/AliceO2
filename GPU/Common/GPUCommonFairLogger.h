// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUCommonFairLogger.h
/// \author David Rohr

#ifndef GPUCOMMONFAIRLOGGER_H
#define GPUCOMMONFAIRLOGGER_H

#if defined(GPUCA_STANDALONE) || defined(GPUCA_ALIROOT_LIB) || defined(GPUCA_GPULIBRARY)

#include <iostream>
#define LOG(type) std::cout
namespace FairLogger
{
static constexpr const char* endl = "\n";
}

#else

#include <FairLogger.h>

#endif

#endif
