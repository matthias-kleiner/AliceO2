// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_TPCCRUHEADER_H
#define O2_TPCCRUHEADER_H

#include "Headers/DataHeader.h"
#include "DataFormatsTPC/Constants.h"

namespace o2
{
namespace tpc
{

/// @struct TPCCRUHeader
/// TPC specific header to be transported on the header stack
struct TPCCRUHeader : public o2::header::BaseHeader {
  // Required to do the lookup
  constexpr static const o2::header::HeaderType sHeaderType = "TPCCRUH";
  static const uint32_t sVersion = 1;
  static constexpr int NCRUs = o2::tpc::constants::MAXCRU;

  // TPCCRUHeader(int cru) : BaseHeader(sizeof(TPCCRUHeader), sHeaderType, o2::header::gSerializationMethodNone, sVersion), cruBits(((uint64_t)0x1) << cru) {}
  TPCCRUHeader(const uint32_t cru) : BaseHeader(sizeof(TPCCRUHeader), sHeaderType, o2::header::gSerializationMethodNone, sVersion), cruBits(cru << 7) {}

  int cru() const
  {
    // for (int cru = 0; cru < NCRUs; ++cru) {
    //   if ((cruBits >> cru) == 0x1) {
    //     return cru;
    //   }
    // }
    // if (cruBits != 0) {
    //   return NCRUs;
    // }
    // return -1;
    // return cruBits;
    return cruBits >> 7;
  }

  const uint32_t cruBits{};
  // union {
  //   uint64_t activeCRUsFlags = 0;
  //   struct {
  //     // uint64_t activeCRUs : NCRUs; // bit field
  //     uint64_t activeCRUs; // bit field
  //     uint64_t unused : 12;
  //     uint64_t flags : 16;
  //   };
  // };
};
} // namespace tpc
} // namespace o2

#endif // O2_TPCCRUHeader_H
