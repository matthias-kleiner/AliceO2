// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_TPC_IDCCRUSim_H_
#define ALICEO2_TPC_IDCCRUSim_H_

#include <vector>
#include <array>
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsTPC/Constants.h"

namespace o2
{
namespace tpc
{

class IDCCRUSim
{
 public:

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
  static constexpr const uint32_t mLengthOfTF{256};                                                                                ///< length of one TF in units of orbits
  const uint32_t mNOrbits{12};                                                                                                     ///< integration of IDCs in units of orbits
  const uint32_t mTimeStamps{o2::constants::lhc::LHCMaxBunches / o2::tpc::constants::LHCBCPERTIMEBIN * mNOrbits};                  ///< number of time stamps for each integration interval
  const uint32_t mIntegrationIntervalsPerTF{mLengthOfTF / mNOrbits + 1};                                                           ///< number of integration intervals per TF. Add 1: 256/12=21.333
  std::array<std::vector<float>, mRegions> mIDCs{                                                                                  ///< IDCs for one sector. The index of the array to the region.
                                                 std::vector<float>(mPadsPerRegion[0] * mIntegrationIntervalsPerTF),               // region 0
                                                 std::vector<float>(mPadsPerRegion[1] * mIntegrationIntervalsPerTF),               // region 1
                                                 std::vector<float>(mPadsPerRegion[2] * mIntegrationIntervalsPerTF),               // region 2
                                                 std::vector<float>(mPadsPerRegion[3] * mIntegrationIntervalsPerTF),               // region 3
                                                 std::vector<float>(mPadsPerRegion[4] * mIntegrationIntervalsPerTF),               // region 4
                                                 std::vector<float>(mPadsPerRegion[5] * mIntegrationIntervalsPerTF),               // region 5
                                                 std::vector<float>(mPadsPerRegion[6] * mIntegrationIntervalsPerTF),               // region 6
                                                 std::vector<float>(mPadsPerRegion[7] * mIntegrationIntervalsPerTF),               // region 7
                                                 std::vector<float>(mPadsPerRegion[8] * mIntegrationIntervalsPerTF),               // region 8
                                                 std::vector<float>(mPadsPerRegion[9] * mIntegrationIntervalsPerTF)};              // region 9
};

} // namespace tpc
} // namespace o2

#endif
