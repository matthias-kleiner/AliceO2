// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file IDCAverageGroup.h
/// \brief class for averaging and grouping of IDCs
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_IDCAVERAGEGROUP_H_
#define ALICEO2_IDCAVERAGEGROUP_H_

#include <vector>
#include "TPCCalibration/IDCGroup.h"
#include "Rtypes.h"
#include "TPCBase/Sector.h"
#include "TPCCalibration/RobustAverage.h"
#include "TPCBase/CalDet.h"

namespace o2::utils
{
class TreeStreamRedirector;
}

class TH2Poly;

namespace o2::tpc
{

template <typename T>
struct Enable_enum_class_bitfield {
  static constexpr bool value = false;
};

// operator overload for allowing bitfiedls with enum
template <typename T>
typename std::enable_if<std::is_enum<T>::value && Enable_enum_class_bitfield<T>::value, T>::type
  operator&(T lhs, T rhs)
{
  typedef typename std::underlying_type<T>::type integer_type;
  return static_cast<T>(static_cast<integer_type>(lhs) & static_cast<integer_type>(rhs));
}

template <typename T>
typename std::enable_if<std::is_enum<T>::value && Enable_enum_class_bitfield<T>::value, T>::type
  operator|(T lhs, T rhs)
{
  typedef typename std::underlying_type<T>::type integer_type;
  return static_cast<T>(static_cast<integer_type>(lhs) | static_cast<integer_type>(rhs));
}

// TODO: https://dalzhim.github.io/2016/02/16/enum-class-bitfields/
enum class PadFlags : unsigned short {
  GOOD = 1 << 0,     ///< flag for a good pad binary 0001
  DEAD = 1 << 1,     ///< flag for a dead pad binary 0010
  UNKNOWN = 1 << 2,  ///< flag for unknown status binary 0100
  SATURATED = 1 << 3 ///< flag for unknown status binary 0100
};

template <>
struct Enable_enum_class_bitfield<PadFlags> {
  static constexpr bool value = true;
};

/// class for averaging and grouping IDCs
/// usage:
/// 1. Define grouping parameters
/// const int region = 3;
/// IDCAverageGroup idcaverage(6, 4, 3, 2, region);
/// 2. set the ungrouped IDCs for one CRU
/// const int nIntegrationIntervals = 3;
/// std::vector<float> idcsungrouped(nIntegrationIntervals*Mapper::PADSPERREGION[region], 11.11); // vector containing IDCs for one region
/// idcaverage.setIDCs(idcsungrouped)
/// 3. perform the averaging and grouping
/// idcaverage.processIDCs();
/// 4. draw IDCs
/// idcaverage.drawUngroupedIDCs(0)
/// idcaverage.drawGroupedIDCs(0)

class IDCAverageGroup
{
  template <class Type>
  class IDCAverageGroupType;

  class Group;
  class Draw;

  template <>
  class IDCAverageGroupType<Group>
  {
   public:
    IDCAverageGroupType(const unsigned char groupPads, const unsigned char groupRows, const unsigned char groupLastRowsThreshold, const unsigned char groupLastPadsThreshold, const unsigned int region, const int nThreads)
      : mIDCsGrouped{groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, region}, mRobustAverage(nThreads){};

    IDCGroup mIDCsGrouped{};                   ///< grouped and averaged IDC values
    std::vector<RobustAverage> mRobustAverage; ///<! object for averaging (each thread will get his one object)
  };

  template <>
  class IDCAverageGroupType<Draw>
  {
   public:
    IDCAverageGroupType(const unsigned char groupPads, const unsigned char groupRows, const unsigned char groupLastRowsThreshold, const unsigned char groupLastPadsThreshold, const unsigned int region, const unsigned int nPads, const PadRegionInfo& padInf, TH2Poly& poly)
      : mIDCsGrouped{groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, region}, mCountDraw(nPads), mPadInf{padInf}, mPoly{poly} {};

    IDCGroup mIDCsGrouped{};     ///< grouped and averaged IDC values
    std::vector<int> mCountDraw; ///< counter to keep track of the already drawn pads
    const PadRegionInfo& mPadInf;
    TH2Poly& mPoly;
    int mGroupCounter = 0;
    int mCol = 0;
    std::array<int, 4> mColors{1, 2, 3, 4};
  };

 public:
  /// constructor
  /// \param groupPads number of pads in pad direction which will be grouped
  /// \param groupRows number of pads in row direction which will be grouped
  /// \param groupLastRowsThreshold minimum number of pads in row direction for the last group in row direction
  /// \param groupLastPadsThreshold minimum number of pads in pad direction for the last group in pad direction
  /// \param region region of the TPC
  /// \param sigma maximum accepted standard deviation for filtering outliers: sigma*stdev
  /// \param overlapRows define parameter for additional overlapping pads in row direction
  /// \param overlapPads define parameter for additional overlapping pads in pad direction
  IDCAverageGroup(const unsigned char groupPads = 4, const unsigned char groupRows = 4, const unsigned char groupLastRowsThreshold = 2, const unsigned char groupLastPadsThreshold = 2, const unsigned int region = 0, const Sector sector = Sector{0}, const float sigma = 3, const unsigned char overlapRows = 0, const unsigned char overlapPads = 0);

  /// Update pad flag map from CCDB
  void updatePadStatusMap();

  /// Set pad flag map directly
  /// \param padStatus CalDet containing for each pad the status flag
  void setPadStatusMap(const CalDet<PadFlags>& padStatus) { mPadStatus = std::make_unique<CalDet<PadFlags>>(padStatus); }

  /// set the IDCs which will be averaged and grouped
  /// \param idcs vector containing the IDCs
  void setIDCs(const std::vector<float>& idcs);

  /// set the IDCs which will be averaged and grouped using move operator
  /// \param IDCs vector containing the IDCs
  void setIDCs(std::vector<float>&& idcs);

  /// \return returns number of integration intervalls stored in this object
  unsigned int getNIntegrationIntervals() const;

  /// grouping and averaging of IDCs
  void processIDCs();

  /// draw plot with information about the performed grouping
  void drawGrouping();

  /// \return returns grouped IDC object
  const auto& getIDCGroup() const { return mGroup.mIDCsGrouped; }

  /// \return returns grouped IDC object
  auto getIDCGroupData() && { return std::move(mGroup.mIDCsGrouped).getData(); }

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName = "IDCAverageGroup.root", const char* outName = "IDCAverageGroup") const;

  /// load ungrouped and grouped IDCs from File
  bool setFromFile(const char* fileName = "IDCAverageGroup.root", const char* name = "IDCAverageGroup");

  /// draw ungrouped IDCs
  /// \param integrationInterval integration interval for which the IDCs will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawUngroupedIDCs(const unsigned int integrationInterval = 0, const std::string filename = "IDCsUngrouped.pdf") const;

  /// draw grouped IDCs
  /// \param integrationInterval integration interval for which the IDCs will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawGroupedIDCs(const unsigned int integrationInterval = 0, const std::string filename = "IDCsGrouped.pdf") const { mGroup.mIDCsGrouped.draw(integrationInterval, filename); }

  /// \return returns the stored ungrouped IDC value for local ungrouped pad row and ungrouped pad
  /// \param ulrow ungrouped local row in region
  /// \param upad ungrouped pad in pad direction
  /// \param integrationInterval integration interval for which the IDCs will be returned
  float getUngroupedIDCValLocal(const unsigned int ulrow, const unsigned int upad, const unsigned int integrationInterval) const { return mIDCsUngrouped[getUngroupedIndex(ulrow, upad, integrationInterval)]; }

  /// \return returns the stored ungrouped IDC value for global ungrouped pad row and ungrouped pad
  /// \param ugrow ungrouped global row
  /// \param upad ungrouped pad in pad direction
  /// \param integrationInterval integration interval for which the IDCs will be returned
  float getUngroupedIDCValGlobal(const unsigned int ugrow, const unsigned int upad, const unsigned int integrationInterval) const { return mIDCsUngrouped[getUngroupedIndexGlobal(ugrow, upad, integrationInterval)]; }

  /// \return returns the stored ungrouped IDC value for local pad number
  /// \param localPadNumber local pad number for region
  /// \param integrationInterval integration interval for which the IDCs will be returned
  float getUngroupedIDCVal(const unsigned int localPadNumber, const unsigned int integrationInterval) const;

  /// \return returns the stored grouped IDC value for local ungrouped pad row and ungrouped pad
  /// \param ulrow local row in region of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  float getGroupedIDCValLocal(unsigned int ulrow, unsigned int upad, unsigned int integrationInterval) const { return mGroup.mIDCsGrouped.getValUngrouped(ulrow, upad, integrationInterval); }

  /// \return returns the stored grouped IDC value for local ungrouped pad row and ungrouped pad
  /// \param ugrow global ungrouped row
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  float getGroupedIDCValGlobal(unsigned int ugrow, unsigned int upad, unsigned int integrationInterval) const { return mGroup.mIDCsGrouped.getValUngroupedGlobal(ugrow, upad, integrationInterval); }

  /// get the number of threads used for some of the calculations
  static int getNThreads() { return sNThreads; }

  /// \return returns sector of which the IDCs are averaged and grouped
  Sector getSector() const { return mSector; }

  /// \return returns ungrouped IDCs
  const auto& getIDCsUngrouped() const { return mIDCsUngrouped; }

  /// \return returns region
  unsigned int getRegion() const { return mGroup.mIDCsGrouped.getRegion(); }

  /// \return returns sigma used for filtering
  float getSigma() const { return mSigma; }

  /// set the number of threads used for some of the calculations
  static void setNThreads(const int nThreads) { sNThreads = nThreads; }

  /// for debugging: creating debug tree
  /// \param nameFile name of the output file
  void createDebugTree(const char* nameFile) const;

  /// for debugging: creating debug tree for integrated IDCs for all objects which are in the same file
  /// \param nameFile name of the output file
  /// \param filename name of the input file containing all objects
  static void createDebugTreeForAllCRUs(const char* nameFile, const char* filename);

 private:
  inline static int sNThreads{1};                                                                                                 ///< number of threads which are used during the calculations
  std::vector<float> mIDCsUngrouped{};                                                                                            ///< integrated ungrouped IDC values per pad
  const Sector mSector{};                                                                                                         ///< sector of averaged and grouped IDCs (used for debugging)
  const float mSigma{};                                                                                                           ///< sigma cut for outlier filtering
  const unsigned char mOverlapRows{0};                                                                                            ///< additional/overlapping pads in row direction
  const unsigned char mOverlapPads{0};                                                                                            ///< additional/overlapping pads in pad direction
  IDCAverageGroupType<Group> mGroup;                                                                                              ///< object for averaging (grouping infornation and robust averaging)
  std::vector<float> mWeightsPad;                                                                                                 ///< storage of the weights in pad direction
  std::vector<float> mWeightsRow;                                                                                                 ///< storage of the weights in row direction
  std::unique_ptr<CalDet<PadFlags>> mPadStatus{std::make_unique<CalDet<PadFlags>>(CalDet<PadFlags>("flags", PadSubset::Region))}; ///< status flag for each pad (i.e. if the pad is dead)

  /// \return returns index to data from ungrouped pad and row
  /// \param ulrow ungrouped local row in region
  /// \param upad ungrouped pad in pad direction
  unsigned int getUngroupedIndex(const unsigned int ulrow, const unsigned int upad, const unsigned int integrationInterval) const;

  /// \return returns index to data from ungrouped pad and row
  /// \param ugrow ungrouped global row
  /// \param upad ungrouped pad in pad direction
  unsigned int getUngroupedIndexGlobal(const unsigned int ugrow, const unsigned int upad, const unsigned int integrationInterval) const;

  /// called from createDebugTreeForAllCRUs()
  static void createDebugTree(const IDCAverageGroup& idcavg, o2::utils::TreeStreamRedirector& pcstream);

  /// normal distribution used for weighting overlapping pads
  /// \param x distance to the center of the normal distribution
  /// \param sigma sigma of the normal distribution
  static float normal_dist(const float x, const float sigma);

  /// perform the loop over the groups bi etiher perform the grouping or the drawing
  /// \param threadNum thread (not used for drawing)
  /// \param integrationInterval integration interval which is averaged (not used for drawing)
  /// \type object containing necessary methods for either perform the grouping or the drawing
  template <class Type>
  void loopOverGroups(const unsigned int threadNum, const unsigned int integrationInterval, Type& type);

  ClassDefNV(IDCAverageGroup, 1)
};

} // namespace o2::tpc

#endif
