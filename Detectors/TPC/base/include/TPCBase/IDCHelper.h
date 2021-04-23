// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file IDCGroup.h
/// \brief class for storing grouped IDCs
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_TPC_IDCHELPER_H_
#define ALICEO2_TPC_IDCHELPER_H_

#include <vector>

namespace o2
{
namespace tpc
{
struct IDCHelper {
  static constexpr unsigned int NREGIONS{10};                                                                                                                      ///< total number of regions in one sector
  static constexpr unsigned int PADROWS{152};                                                                                                                      ///< total number of pad rows
  static constexpr unsigned int PADSPERREGION[NREGIONS]{1200, 1200, 1440, 1440, 1440, 1440, 1600, 1600, 1600, 1600};                                               ///< number of pads per CRU
  static constexpr unsigned int GLOBALPADOFFSET[NREGIONS]{0, 1200, 2400, 3840, 5280, 6720, 8160, 9760, 11360, 12960};                                              ///< offset of number of pads for region used for debugging only
  static constexpr unsigned int ROWSPERREGION[NREGIONS]{17, 15, 16, 15, 18, 16, 16, 14, 13, 12};                                                                   ///< number of pad rows for region
  static constexpr float PADAREA[NREGIONS]{1 / 0.312f, 1 / 0.315f, 1 / 0.315f, 1 / 0.327f, 1 / 0.6f, 1 / 0.6f, 1 / 0.7296f, 1 / 0.7056f, 1 / 0.906f, 1 / 0.9105f}; ///< inverse size of the pad area padwidth*padLength
  inline static std::vector<unsigned int> ADDITIONALPADSPERROW[NREGIONS]{
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
  };                                                        ///< additional pads per row compared to first row
  inline static std::vector<unsigned int> OFFSETCRULOCAL[NREGIONS]{
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
  inline static std::vector<unsigned int> PADSPERROW[NREGIONS]{
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
  static constexpr unsigned int OFFSETCRUGLOBAL[PADROWS]{
    0, 66, 132, 198, 266, 334, 402, 472, 542, 612, 684, 756, 828, 902, 976, 1050, 1124,         // region 0
    0, 76, 152, 228, 306, 384, 462, 542, 622, 702, 784, 866, 948, 1032, 1116,                   // region 1
    0, 86, 172, 258, 346, 434, 522, 612, 702, 792, 882, 974, 1066, 1158, 1252, 1346,            // region 2
    0, 92, 184, 276, 370, 464, 558, 654, 750, 846, 944, 1042, 1140, 1240, 1340,                 // region 3
    0, 76, 152, 228, 304, 382, 460, 538, 618, 698, 778, 858, 940, 1022, 1104, 1188, 1272, 1356, // region 4
    0, 86, 172, 258, 346, 434, 522, 612, 702, 792, 882, 974, 1066, 1158, 1252, 1346,            // region 5
    0, 94, 190, 286, 382, 480, 578, 676, 776, 876, 978, 1080, 1182, 1286, 1390, 1494,           // region 6
    0, 110, 220, 332, 444, 556, 670, 784, 898, 1014, 1130, 1246, 1364, 1482,                    // region 7
    0, 118, 236, 356, 476, 598, 720, 844, 968, 1092, 1218, 1344, 1472,                          // region 8
    0, 128, 258, 388, 520, 652, 784, 918, 1052, 1188, 1324, 1462                                // region 9
  };                                                                                            ///< row offset in cru for given global pad row
};

} // namespace tpc
} // namespace o2

#endif
