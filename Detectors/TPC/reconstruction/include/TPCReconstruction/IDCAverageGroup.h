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

#include "Framework/Logger.h"
namespace o2
{
namespace tpc
{
class IDCAverageGroup
{

 public:
  IDCAverageGroup(const unsigned int groupPads = 4, const unsigned int groupRows = 4, const unsigned int groupLastRowsThreshold = 2, const unsigned int groupLastPadsThreshold = 2, const unsigned int cru = 0)
    : mGroupPads{groupPads}, mGroupRows{groupRows}, mGroupLastRowsThreshold{groupLastRowsThreshold}, mGroupLastPadsThreshold{groupLastPadsThreshold}, mCRU{cru}
  {
    initIDCGroup();
  }

  /// \return returns the row of the group from the local row in a region
  /// \param lrow local row in a region
  unsigned int getGroupedRow(const unsigned int lrow) const;

  /// \return returns the grouped pad index from ungrouped pad and row
  /// \param pad ungrouped pad
  /// \param lrow local ungrouped row in a region
  int getGroupedPad(const unsigned int pad, const unsigned int lrow) const;

  /// \return returns the global pad number for given local pad row and pad
  /// \param lrow local row in a region
  /// \param lrow local row in a region
  unsigned int getGlobalPadNumber(const unsigned int lrow, const unsigned int pad) const { return mGlobalRowOff[mRegion] + mOffs[mRegion][lrow] + pad; }

  /// \param IDCs vector containg the IDCs
  void setIDCs(const std::vector<float>& idcs)
  {
    mIDCs = idcs;
    mIDCsGrouped.resize(getNIntegrationIntervals());
  }

  /// \param IDCs vector containg the IDCs
  void setIDCs(std::vector<float>&& idcs)
  {
    mIDCs = std::move(idcs);
    mIDCsGrouped.resize(getNIntegrationIntervals());
  }

  /// \return returns number of integration intervalls stored in this object
  unsigned int getNIntegrationIntervals() const { return mIDCs.size() / mPadsPerRegion[mRegion]; }

  /// grouping and averaging of IDCs
  /// \param debug
  void processIDCs();

  /// \return returns grouped IDC object
  const auto& getIDCGroup() const { return mIDCsGrouped; }

  /// \return returns the number pads in pad direction which are grouped
  unsigned int getGroupPads() const { return mGroupPads; }

  /// \return returns the number pads in row direction which are grouped
  unsigned int getGroupRows() const { return mGroupRows; }

  /// \return returns threshold for grouping the last group in row direction
  unsigned int getGroupLastRowsThreshold() const { return mGroupLastRowsThreshold; }

  /// \return returns threshold for grouping the last group in pad direction
  unsigned int getGroupLastPadsThreshold() const { return mGroupLastPadsThreshold; }

  /// \return returns the region of the CRU for which the IDCs are grouped
  unsigned int getRegion() const { return mRegion; }

  /// \return returns the number of the CRU for which the IDCs are grouped
  uint32_t getCRU() const { return mCRU; }

  /// draw grouped IDCs
  /// \param integrationInterval integration interval for which the IDCs will be drawn
  void draw(const int integrationInterval = 0) const;

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName, const char* outName = "IDCGroup") const;

 private:
  static constexpr unsigned int mRegions{10};
  static constexpr unsigned int mPadsPerRegion[mRegions]{1200, 1200, 1440, 1440, 1440, 1440, 1600, 1600, 1600, 1600};                                               ///< number of pads per CRU
  static constexpr unsigned int mGlobalRowOff[mRegions]{0, 1200, 2400, 3840, 5280, 6720, 8160, 9760, 11360, 12960};                                                 ///< global row offset
  static constexpr unsigned int mRowsPerRegion[mRegions]{17, 15, 16, 15, 18, 16, 16, 14, 13, 12};                                                                   ///< number of pad rows for region
  static constexpr float mPadArea[mRegions]{1 / 0.312f, 1 / 0.315f, 1 / 0.315f, 1 / 0.327f, 1 / 0.6f, 1 / 0.6f, 1 / 0.7296f, 1 / 0.7056f, 1 / 0.906f, 1 / 0.9105f}; ///< inverse size of the pad area padwidth*padLength
  inline static std::vector<unsigned int> mAddPadsPerRow[mRegions]{
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
  inline static std::vector<unsigned int> mOffs[mRegions]{
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
  inline static std::vector<unsigned int> mPadsPerRow[mRegions]{
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
  const unsigned int mGroupPads{4};                                            ///< group 4 pads
  const unsigned int mGroupRows{4};                                            ///< group 4 pads -> 4x4
  const unsigned int mGroupLastRowsThreshold{2};                               ///< if the last group (region edges) consists in row direction less then mGroupLastRowsThreshold pads then it will be grouped into the previous group
  const unsigned int mGroupLastPadsThreshold{2};                               ///< if the last group (sector edges) consists in pad direction less then mGroupLastPadsThreshold pads then it will be grouped into the previous group
  const unsigned int mCRU{};                                                   ///< CRU of input IDCs
  const unsigned int mRegion{mCRU % mRegions};                                 ///< region of input IDCs
  std::vector<float> mIDCs{};                                                  ///< integrated IDC values
  IDCGroup mIDCsGrouped{};                                                     ///< grouped and averaged IDC values

  int getLastRow() const;
  int getLastPad(const int row) const;
  void initIDCGroup();
};

} // namespace tpc
} // namespace o2

#endif
