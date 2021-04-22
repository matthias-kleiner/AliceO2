// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file IDCAverageGroup.h
/// \brief class for averaging/grouping and storing to IDCs
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_IDCAVERAGEGROUP_H_
#define ALICEO2_IDCAVERAGEGROUP_H_

#include <vector>
#include "TPCBase/IDCGroup.h"

namespace o2
{
namespace tpc
{
class IDCAverageGroup
{

 public:
  IDCAverageGroup(const int groupPads = 4, const int groupRows = 4, const int groupLastRowsThreshold = 2, const int groupLastPadsThreshold = 2, const uint32_t cru = 0, const bool debug = false)
    : mGroupPads{groupPads}, mGroupRows{groupRows}, mGroupLastRowsThreshold{groupLastRowsThreshold}, mGroupLastPadsThreshold{groupLastPadsThreshold}, mCRU{cru}, mDebug{debug}
  {
    initIDCGroup();
  }

  static uint32_t getGlobalPadIndex(const uint32_t region, const uint32_t lrow, const uint32_t pad) { return mGlobalRowOff[region] + mOffs[region][lrow] + pad; }

  /// \param IDCs vector containg the IDCs
  void setIDCs(const std::vector<float>& idcs)
  {
    mIDCs = idcs;
    mIntegrationIntervals = mIDCs.size() / mPadsPerRegion[mRegion];
  }

  /// grouping and averaging of IDCs
  /// \param debug
  void processIDCs(const bool debug = false, const int lane = 0);

  const auto& getIDCGroup() const { return mIDCsGrouped; }

 private:
  static constexpr uint32_t mRegions{10};
  static constexpr uint32_t mPadsPerRegion[mRegions]{1200, 1200, 1440, 1440, 1440, 1440, 1600, 1600, 1600, 1600}; ///< number of pads per CRU
  static constexpr uint32_t mGlobalRowOff[mRegions]{0, 1200, 2400, 3840, 5280, 6720, 8160, 9760, 11360, 12960};   ///< global row offset
  static constexpr uint32_t mRowsPerRegion[mRegions]{17, 15, 16, 15, 18, 16, 16, 14, 13, 12};                     ///< number of pad rows for region
  const int mGroupPads{4};                                                                                        ///< group 4 pads
  const int mGroupRows{4};                                                                                        ///< group 4 pads -> 4x4
  const int mGroupLastRowsThreshold{2};                                                                           ///< if the last group (region edges) consists in row direction less then mGroupLastRowsThreshold pads then it will be grouped into the previous group
  const int mGroupLastPadsThreshold{2};                                                                           ///< if the last group (sector edges) consists in pad direction less then mGroupLastPadsThreshold pads then it will be grouped into the previous group
  const uint32_t mCRU{};                                                                                          ///< CRU of input IDCs
  const bool mDebug{};                                                                                            ///< dump IDCs to tree for debugging
  const uint32_t mRegion{mCRU % mRegions};                                                                        ///< region of input IDCs
  inline static std::vector<uint32_t> mAddPadsPerRow[mRegions]{
    {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5},    // region 0
    {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4},          // region 1
    {0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4},       // region 2
    {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4},          // region 3
    {0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4}, // region 4
    {0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4},       // region 5
    {0, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6},       // region 6
    {0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4},             // region 7
    {0, 0, 1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5},                // region 8
    {0, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5}                    // region 9
  };                                                        ///< additional pads per row on compared to first row
  inline static std::vector<uint32_t> mOffs[mRegions]{
    {0, 66, 132, 198, 266, 334, 402, 472, 542, 612, 684, 756, 828, 902, 976, 1050, 1124},         // region 0
    {0, 76, 152, 228, 306, 384, 462, 542, 622, 702, 784, 866, 948, 1032, 1116},                   // region 1
    {0, 86, 172, 258, 346, 434, 522, 612, 702, 792, 882, 974, 1066, 1158, 1252, 1346},            // region 2
    {0, 92, 184, 276, 370, 464, 558, 654, 750, 846, 944, 1042, 1140, 1240, 1340},                 // region 3
    {0, 76, 152, 228, 304, 382, 460, 538, 618, 698, 778, 858, 940, 1022, 1104, 1188, 1272, 1356}, // region 4
    {0, 86, 172, 258, 346, 434, 522, 612, 702, 792, 882, 974, 1066, 1158, 1252, 1346},            // region 5
    {0, 94, 190, 286, 382, 480, 578, 676, 776, 876, 978, 1080, 1182, 1286, 1390, 1494},           // region 6
    {0, 110, 220, 332, 444, 556, 670, 784, 898, 1014, 1130, 1246, 1364, 1482},                    // region 7
    {0, 118, 236, 356, 476, 598, 720, 844, 968, 1092, 1218, 1344, 1472},                          // region 8
    {0, 128, 258, 388, 520, 652, 784, 918, 1052, 1188, 1324, 1462}                                // region 9
  };                                                                                              ///< row offset in cru for given local pad row
  inline static std::vector<uint32_t> mPadsPerRow[mRegions]{
    {66, 66, 66, 68, 68, 68, 70, 70, 70, 72, 72, 72, 74, 74, 74, 74, 76},      // region 0
    {76, 76, 76, 78, 78, 78, 80, 80, 80, 82, 82, 82, 84, 84, 84},              // region 1
    {86, 86, 86, 88, 88, 88, 90, 90, 90, 90, 92, 92, 92, 94, 94, 94},          // region 2
    {92, 92, 92, 94, 94, 94, 96, 96, 96, 98, 98, 98, 100, 100, 100},           // region 3
    {76, 76, 76, 76, 78, 78, 78, 80, 80, 80, 80, 82, 82, 82, 84, 84, 84, 84},  // region 4
    {86, 86, 86, 88, 88, 88, 90, 90, 90, 90, 92, 92, 92, 94, 94, 94},          // region 5
    {94, 96, 96, 96, 98, 98, 98, 100, 100, 102, 102, 102, 104, 104, 104, 106}, // region 6
    {110, 110, 112, 112, 112, 114, 114, 114, 116, 116, 116, 118, 118, 118},    // region 7
    {118, 118, 120, 120, 122, 122, 124, 124, 124, 126, 126, 128, 128},         // region 8
    {128, 130, 130, 132, 132, 132, 134, 134, 136, 136, 138, 138}               // region 9
  };                                                                           ///< number of pads per row in region
  std::vector<float> mIDCs{};                                                  ///< integrated IDC values.
  IDCGroup mIDCsGrouped{};                                                     ///< group and averaged IDC values.
  int mIntegrationIntervals{};                                                 ///< number of integration intervals

  int getLastRow() const;
  int getLastPad(const int row) const;
  void initIDCGroup();
};

} // namespace tpc
} // namespace o2

#endif
