// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file IDCSim.h
/// \brief Integration of IDCs
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>
/// \date Apr 16, 2021

#ifndef ALICEO2_TPC_IDCSim_H_
#define ALICEO2_TPC_IDCSim_H_

#include <vector>
#include <array>
#include "DataFormatsTPC/Digit.h"
#include "DataFormatsTPC/Constants.h"
#include "CommonConstants/LHCConstants.h"
#include "TPCBase/Mapper.h"
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
  /// \param sector sector for which the data is processed
  /// \param nOrbits length of integration intervals
  IDCSim(const unsigned int sector = 0, const unsigned int nOrbits = 12) : mSector{sector}, mNOrbits{nOrbits} {}

  // integrate IDCs for one TF
  /// \param digits digits for one sector for one Time Frame
  void integrateDigitsForOneTF(const gsl::span<const o2::tpc::Digit>& digits);

  /// \return returns the total number of integration intervals for one TF
  /// \param nIDCs total number of IDCs for one TF for region
  /// \param region region of the tpc
  static int getNIntegrationIntervals(const int nIDCs, const int region) { return nIDCs / Mapper::PADSPERREGION[region]; }

  /// \param row global pad row
  /// \param pad pad in row
  /// \return returns local pad number in region
  static unsigned int getPadIndex(const int row, const int pad) { return Mapper::OFFSETCRUGLOBAL[row] + pad; }

  /// set number of orbits per TF which is used to determine the size of the vectors etc.
  /// \param nOrbitsPerTF number of orbits per TF
  static void setNOrbitsPerTF(const unsigned int nOrbitsPerTF) { mOrbitsPerTF = nOrbitsPerTF; }

  /// \return returns the number of orbits for one TF
  static unsigned int getNOrbitsPerTF() { return mOrbitsPerTF; }

  /// for debugging: dumping IDCs to ROOT file
  /// \param timeframe timeframe of the IDCs to avoid overwriting the output file for same TF
  void dumpIDCs(const int timeframe);

  /// for debugging: creating debug tree for integrated IDCs
  /// \param timeframe timeframe of the IDCs to avoid overwriting the output file for same TF
  void createDebugTree(const int timeframe);

  /// return return the IDCs for all sector
  auto& get() { return mIDCs[!mBufferIndex]; }
  const auto& get() const { return mIDCs[!mBufferIndex]; }

  /// \return returns the sector for which the IDCs are integrated
  unsigned int getSector() const { return mSector; }

 private:
  inline static unsigned int mOrbitsPerTF{256};                                                                                               ///< length of one TF in units of orbits
  const unsigned int mSector{};                                                                                                               ///< sector for which the IDCs are integrated
  const unsigned int mNOrbits{12};                                                                                                            ///< integration intervals of IDCs in units of orbits
  const unsigned int mTimeStampsPerIntegrationInterval{(o2::constants::lhc::LHCMaxBunches * mNOrbits) / o2::tpc::constants::LHCBCPERTIMEBIN}; ///< number of time stamps for each integration interval (5346)
  const bool mAddInterval{mOrbitsPerTF % mNOrbits > 0 ? true : false};                                                                        ///< if the division has a reminder 256/12=21.333 then add an additional integration interval
  const unsigned int mIntegrationIntervalsPerTF{mOrbitsPerTF / mNOrbits + mAddInterval};                                                      ///< number of integration intervals per TF. Add 1: 256/12=21.333
  const unsigned int mTimeStampsReminder{mTimeStampsPerIntegrationInterval * (mOrbitsPerTF % mNOrbits) / mNOrbits};                           ///< number time stamps which remain in one TF and will be buffered to the next TF
  int mTimeBinsOff{};                                                                                                                         ///< offset from last time bin
  int mBufferIndex{};                                                                                                                         ///< index for the buffer
  const std::array<unsigned int, Mapper::NREGIONS> mMaxIDCs{
    Mapper::PADSPERREGION[0] * mIntegrationIntervalsPerTF, // region 0
    Mapper::PADSPERREGION[1] * mIntegrationIntervalsPerTF, // region 1
    Mapper::PADSPERREGION[2] * mIntegrationIntervalsPerTF, // region 2
    Mapper::PADSPERREGION[3] * mIntegrationIntervalsPerTF, // region 3
    Mapper::PADSPERREGION[4] * mIntegrationIntervalsPerTF, // region 4
    Mapper::PADSPERREGION[5] * mIntegrationIntervalsPerTF, // region 5
    Mapper::PADSPERREGION[6] * mIntegrationIntervalsPerTF, // region 6
    Mapper::PADSPERREGION[7] * mIntegrationIntervalsPerTF, // region 7
    Mapper::PADSPERREGION[8] * mIntegrationIntervalsPerTF, // region 8
    Mapper::PADSPERREGION[9] * mIntegrationIntervalsPerTF  // region 9
  };                                                                   ///< maximum number of IDCs per region
  std::array<std::vector<float>, Mapper::NREGIONS> mIDCs[2]{
    {std::vector<float>(mMaxIDCs[0]),  // region 0
     std::vector<float>(mMaxIDCs[1]),  // region 1
     std::vector<float>(mMaxIDCs[2]),  // region 2
     std::vector<float>(mMaxIDCs[3]),  // region 3
     std::vector<float>(mMaxIDCs[4]),  // region 4
     std::vector<float>(mMaxIDCs[5]),  // region 5
     std::vector<float>(mMaxIDCs[6]),  // region 6
     std::vector<float>(mMaxIDCs[7]),  // region 7
     std::vector<float>(mMaxIDCs[8]),  // region 8
     std::vector<float>(mMaxIDCs[9])}, // region 9
    {std::vector<float>(mMaxIDCs[0]),  // region 0
     std::vector<float>(mMaxIDCs[1]),  // region 1
     std::vector<float>(mMaxIDCs[2]),  // region 2
     std::vector<float>(mMaxIDCs[3]),  // region 3
     std::vector<float>(mMaxIDCs[4]),  // region 4
     std::vector<float>(mMaxIDCs[5]),  // region 5
     std::vector<float>(mMaxIDCs[6]),  // region 6
     std::vector<float>(mMaxIDCs[7]),  // region 7
     std::vector<float>(mMaxIDCs[8]),  // region 8
     std::vector<float>(mMaxIDCs[9])}  // region 9
  };                                   ///< IDCs for one sector. The array is needed to buffer the IDCs for the last integration interval

  /// \return returns the last time bin after which the buffer is switched
  unsigned int getLastTimeBinForSwitch() const;
  int getNewOffset() const;

  /// set all IDC values to 0
  void resetIDCs();

  /// return orbit for given timeStamp
  unsigned int getOrbit(const unsigned int timeStamp) const { return static_cast<unsigned int>((timeStamp + mTimeBinsOff) / mTimeStampsPerIntegrationInterval); }

  /// \return returns index in the vector of the mIDCs member
  /// \param timeStamp timeStamp for which the index is calculated
  /// \param region region in the sector
  /// \param row global pad row
  /// \param pad pad in row
  unsigned int getIndex(const int timeStamp, const int region, const int row, const int pad) const { return getOrbit(timeStamp) * Mapper::PADSPERREGION[region] + getPadIndex(row, pad); }
};

} // namespace tpc
} // namespace o2

#endif
