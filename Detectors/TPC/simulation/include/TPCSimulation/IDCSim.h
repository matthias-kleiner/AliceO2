// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_TPC_IDCSim_H_
#define ALICEO2_TPC_IDCSim_H_

#include <vector>
#include <array>
#include <fmt/format.h>
#include "Framework/Logger.h"
#include "DataFormatsTPC/Digit.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsTPC/Constants.h"
#include "CommonUtils/TreeStreamRedirector.h" // for debugging
#include "TFile.h"
#include <gsl/span>

namespace o2
{
namespace tpc
{

/// \class IDCSim
/// This class is for the integration of IDCs for one sector.
/// The input has to be provided for one TF.
/// The IDCs are stored per CRU for all integration intervals.

class IDCSim
{
 public:
  /// constructor
  /// \param nOrbits length of integration intervals
  IDCSim(const uint32_t sector = 0, const uint32_t nOrbits = 12) : mSector{sector}, mNOrbits{nOrbits} {}

  // integrate IDCs for one TF
  /// \param digits digits for one sector for one Time Frame
  void integrateDigitsForOneTF(const gsl::span<const o2::tpc::Digit>& digits);

  /// \return returns the total number of regions in the TPC for one sector
  static int getNRegions() { return mRegions; }

  /// \return returns pad offset to calculate global pad number from pad number in cru
  /// \param region region of the tpc
  static int getGlobalPadOffset(const int region) { return mGlobalPadOffs[region]; }

  /// \return returns the nummber of pads for given region
  /// \param region region of the tpc
  static int getPadsPerRegion(const int region) { return mPadsPerRegion[region]; }

  /// \return returns the total number of integration intervals for one TF
  /// \param nIDCs total number of IDCs for one TF for region
  /// \param region region of the tpc
  static int getNIntegrationIntervals(const int nIDCs, const int region) { return nIDCs / mPadsPerRegion[region]; }

  /// \param row global pad row
  /// \param pad pad in row
  /// \return returns local pad number in region
  static unsigned int getPadIndex(const int row, const int pad) { return mOffs[row] + pad; }

  /// set number of orbits per TF which is used to determine the size of the vectors etc.
  /// \param nOrbitsPerTF number of orbits per TF
  static void setNOrbitsPerTF(const uint32_t nOrbitsPerTF) { mOrbitsPerTF = nOrbitsPerTF; }

  /// \return returns the number of orbits for one TF
  static uint32_t getNOrbitsPerTF() { return mOrbitsPerTF; }

  /// for debugging: dumping IDCs to ROOT file
  void dumpIDCs(const int timeframe);

  /// for debugging: creating debug tree for integrated IDCs
  void createDebugTree(const int timeframe);

  /// return return the IDCs for all sector
  auto& get() { return mIDCsTmp[!mBufferIndex]; }
  const auto& get() const { return mIDCsTmp[!mBufferIndex]; }

  /// \return returns the sector for which the IDCs are integrated
  uint32_t getSector() const { return mSector; }

 private:
  inline static uint32_t mOrbitsPerTF{256};                                                                                               ///< length of one TF in units of orbits
  static constexpr int mRegions{10};                                                                                                      ///< number of regions per sector
  static constexpr int mPadRows{152};                                                                                                     ///< total number of pad rows
  static constexpr int mPadsPerRegion[mRegions]{1200, 1200, 1440, 1440, 1440, 1440, 1600, 1600, 1600, 1600};                              ///< number of pads per CRU
  static constexpr int mGlobalPadOffs[mRegions]{0, 1200, 2400, 3840, 5280, 6720, 8160, 9760, 11360, 12960};                               ///< offset of number of pads for region used for debugging only
  static constexpr int mOffs[mPadRows]{                                                                                                   ///< row offset in cru for given global pad row
                                       0, 66, 132, 198, 266, 334, 402, 472, 542, 612, 684, 756, 828, 902, 976, 1050, 1124,                // region 0
                                       0, 76, 152, 228, 306, 384, 462, 542, 622, 702, 784, 866, 948, 1032, 1116,                          // region 1
                                       0, 86, 172, 258, 346, 434, 522, 612, 702, 792, 882, 974, 1066, 1158, 1252, 1346,                   // region 2
                                       0, 92, 184, 276, 370, 464, 558, 654, 750, 846, 944, 1042, 1140, 1240, 1340,                        // region 3
                                       0, 76, 152, 228, 304, 382, 460, 538, 618, 698, 778, 858, 940, 1022, 1104, 1188, 1272, 1356,        // region 4
                                       0, 86, 172, 258, 346, 434, 522, 612, 702, 792, 882, 974, 1066, 1158, 1252, 1346,                   // region 5
                                       0, 94, 190, 286, 382, 480, 578, 676, 776, 876, 978, 1080, 1182, 1286, 1390, 1494,                  // region 6
                                       0, 110, 220, 332, 444, 556, 670, 784, 898, 1014, 1130, 1246, 1364, 1482,                           // region 7
                                       0, 118, 236, 356, 476, 598, 720, 844, 968, 1092, 1218, 1344, 1472,                                 // region 8
                                       0, 128, 258, 388, 520, 652, 784, 918, 1052, 1188, 1324, 1462};                                     // region 9
  const uint32_t mSector{};                                                                                                               ///< sector for which the IDCs are integrated
  const uint32_t mNOrbits{12};                                                                                                            ///< integration intervals of IDCs in units of orbits
  const uint32_t mTimeStampsPerIntegrationInterval{(o2::constants::lhc::LHCMaxBunches * mNOrbits) / o2::tpc::constants::LHCBCPERTIMEBIN}; ///< number of time stamps for each integration interval (5346)
  const bool mAddInterval{mOrbitsPerTF % mNOrbits > 0 ? true : false};                                                                    /// if the division has a reminder 256/12=21.333 then add an additional integration interval
  const uint32_t mIntegrationIntervalsPerTF{mOrbitsPerTF / mNOrbits + mAddInterval};                                                      ///< number of integration intervals per TF. Add 1: 256/12=21.333
  const uint32_t mTimeStampsReminder{mTimeStampsPerIntegrationInterval * (mOrbitsPerTF % mNOrbits) / mNOrbits};                           ///< number time stamps which remain in one TF and will be buffered to the next TF
  int mTimeBinsOff{};                                                                                                                     ///< offset from last time bin
  int mBufferIndex{};                                                                                                                     ///< index for the buffer
  std::array<std::vector<float>, mRegions> mIDCsTmp[2]{                                                                                   ///< IDCs for one sector. The array is needed to buffer the IDCs for the last integration interval
                                                       {std::vector<float>(mPadsPerRegion[0] * mIntegrationIntervalsPerTF),               // region 0
                                                        std::vector<float>(mPadsPerRegion[1] * mIntegrationIntervalsPerTF),               // region 1
                                                        std::vector<float>(mPadsPerRegion[2] * mIntegrationIntervalsPerTF),               // region 2
                                                        std::vector<float>(mPadsPerRegion[3] * mIntegrationIntervalsPerTF),               // region 3
                                                        std::vector<float>(mPadsPerRegion[4] * mIntegrationIntervalsPerTF),               // region 4
                                                        std::vector<float>(mPadsPerRegion[5] * mIntegrationIntervalsPerTF),               // region 5
                                                        std::vector<float>(mPadsPerRegion[6] * mIntegrationIntervalsPerTF),               // region 6
                                                        std::vector<float>(mPadsPerRegion[7] * mIntegrationIntervalsPerTF),               // region 7
                                                        std::vector<float>(mPadsPerRegion[8] * mIntegrationIntervalsPerTF),               // region 8
                                                        std::vector<float>(mPadsPerRegion[9] * mIntegrationIntervalsPerTF)},
                                                       {std::vector<float>(mPadsPerRegion[0] * mIntegrationIntervalsPerTF),   // region 0
                                                        std::vector<float>(mPadsPerRegion[1] * mIntegrationIntervalsPerTF),   // region 1
                                                        std::vector<float>(mPadsPerRegion[2] * mIntegrationIntervalsPerTF),   // region 2
                                                        std::vector<float>(mPadsPerRegion[3] * mIntegrationIntervalsPerTF),   // region 3
                                                        std::vector<float>(mPadsPerRegion[4] * mIntegrationIntervalsPerTF),   // region 4
                                                        std::vector<float>(mPadsPerRegion[5] * mIntegrationIntervalsPerTF),   // region 5
                                                        std::vector<float>(mPadsPerRegion[6] * mIntegrationIntervalsPerTF),   // region 6
                                                        std::vector<float>(mPadsPerRegion[7] * mIntegrationIntervalsPerTF),   // region 7
                                                        std::vector<float>(mPadsPerRegion[8] * mIntegrationIntervalsPerTF),   // region 8
                                                        std::vector<float>(mPadsPerRegion[9] * mIntegrationIntervalsPerTF)}}; // region 9

  unsigned int getLastTimeBinForSwitch() const
  {
    const int totaloffs = mTimeBinsOff + static_cast<int>(mTimeStampsReminder);
    return totaloffs >= mTimeStampsPerIntegrationInterval ? mIntegrationIntervalsPerTF * mTimeStampsPerIntegrationInterval : (mIntegrationIntervalsPerTF - mAddInterval) * mTimeStampsPerIntegrationInterval - mTimeBinsOff;
  }

  int getNewOffset() const
  {
    const int totaloffs = mTimeBinsOff + static_cast<int>(mTimeStampsReminder);
    return totaloffs >= mTimeStampsPerIntegrationInterval ? (totaloffs - static_cast<int>(mTimeStampsPerIntegrationInterval)) : totaloffs;
  }

  /// set all IDC values to 0
  void resetIDCs()
  {
    for (auto& idcs : mIDCsTmp[!mBufferIndex]) {
      std::fill(idcs.begin(), idcs.end(), 0);
    }
  }

  /// return orbit for given timeStamp
  unsigned int getOrbit(const unsigned int timeStamp) const
  {
    return static_cast<unsigned int>((timeStamp + mTimeBinsOff) / mTimeStampsPerIntegrationInterval);
  }

  /// \return returns index in the vector of the mIDCs member
  /// \param timeStamp timeStamp for which the index is calculated
  /// \param region region in the sector
  /// \param row global pad row
  /// \param pad pad in row
  unsigned int getIndex(const int timeStamp, const int region, const int row, const int pad) const
  {
    return getOrbit(timeStamp) * mPadsPerRegion[region] + getPadIndex(row, pad);
  }
};

} // namespace tpc
} // namespace o2

#endif