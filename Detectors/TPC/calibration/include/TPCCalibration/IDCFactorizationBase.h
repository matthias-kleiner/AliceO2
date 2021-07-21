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

/// \file IDCFactorizationBase.h
/// \brief base class for holding factorization of IDCs for for a TPC and one cru
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>
/// \date Jun 10, 2021

#ifndef ALICEO2_IDCFACTORIZATIONBASE_H_
#define ALICEO2_IDCFACTORIZATIONBASE_H_

#include <vector>
#include "Rtypes.h"
#include "TPCCalibration/IDCGroupHelperSector.h"
#include "TPCCalibration/IDCGroupHelperRegion.h"
#include "TPCCalibration/IDCContainer.h"

namespace o2::tpc
{

template <class Type>
class IDCFactorizationBase;

using IDCFactorizationSector = o2::tpc::IDCFactorizationBase<IDCGroupHelperSector>;
using IDCFactorizationRegion = o2::tpc::IDCFactorizationBase<IDCGroupHelperRegion>;

/// template specialization to perform factorization for the full TPC (all sectors)
template <>
class IDCFactorizationBase<IDCGroupHelperSector> : public IDCGroupHelperSector
{
 public:
  /// constructor
  /// \param groupPads number of pads in pad direction which will be grouped for all regions
  /// \param groupRows number of pads in row direction which will be grouped for all regions
  /// \param groupLastRowsThreshold minimum number of pads in row direction for the last group in row direction for all regions
  /// \param groupLastPadsThreshold minimum number of pads in pad direction for the last group in pad direction for all regions
  /// \param timeFrames number of time frames which will be stored
  /// \param timeframesDeltaIDC number of time frames stored for each DeltaIDC object
  IDCFactorizationBase(const std::array<unsigned char, Mapper::NREGIONS>& groupPads, const std::array<unsigned char, Mapper::NREGIONS>& groupRows, const std::array<unsigned char, Mapper::NREGIONS>& groupLastRowsThreshold, const std::array<unsigned char, Mapper::NREGIONS>& groupLastPadsThreshold, const unsigned int timeFrames, const unsigned int timeframesDeltaIDC);

  /// default constructor for ROOT I/O
  IDCFactorizationBase() = default;

  /// set the IDC data
  /// \param idcs vector containing the IDCs
  /// \param cru CRU of the data
  /// \param timeframe time frame of the IDCs
  void setIDCs(std::vector<float>&& idcs, const unsigned int cru, const unsigned int timeframe) { mIDCs[cru][timeframe] = std::move(idcs); }

  /// calculate I_0(r,\phi) = <I(r,\phi,t)>_t
  void calcIDCZero();

  /// calculate I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  void calcIDCOne();

  /// calculate \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  void calcIDCDelta();

  /// \return returns number of IDC values for region
  /// \param region region of the TPC
  auto getIDCsPerCRU(const unsigned int region) const { return mNIDCsPerCRU[region]; }

  // get number of TFs in which the DeltaIDCs are split/stored
  unsigned int getTimeFramesDeltaIDC() const { return mTimeFramesDeltaIDC; }

  /// \return returns the stored grouped and integrated IDC
  /// \param sector sector
  /// \param region region
  /// \param urow row in the region of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param tf time frame
  /// \param interval integration interval
  float getIDCValUngrouped(const unsigned int sector, const unsigned int region, const unsigned int urow, const unsigned int upad, const unsigned int tf, const unsigned int interval) const
  {
    return mIDCs[sector * Mapper::NREGIONS + region][tf][interval * getIDCsPerCRU(region) + mOffsRow[region][getGroupedRow(region, urow)] + getGroupedPad(region, urow, upad)];
  }

  /// \return returns the stored IDC0 value for local ungrouped pad row and ungrouped pad
  /// \param sector sector
  /// \param region region
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  float getIDCZeroVal(const unsigned int sector, const unsigned int region, unsigned int urow, unsigned int upad) const
  {
    return mIDCZero.getValueIDCZero(Sector(sector).side(), getIndexUngrouped(sector, region, urow, upad, 0));
  }

  /// \return returns the stored DeltaIDC value for local ungrouped pad row and ungrouped pad
  /// \param sector sector
  /// \param region region
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param chunk chunk of the Delta IDC (can be obtained with getLocalIntegrationInterval())
  /// \param localintegrationInterval local integration interval for chunk (can be obtained with getLocalIntegrationInterval())
  float getIDCDeltaVal(const unsigned int sector, const unsigned int region, unsigned int urow, unsigned int upad, unsigned int chunk, unsigned int localintegrationInterval) const
  {
    return mIDCDelta[chunk].getValue(Sector(sector).side(), getIndexUngrouped(sector, region, urow, upad, localintegrationInterval));
  }

  /// \return returns total number of integration interval
  unsigned long getNIntegrationIntervals() const;

  /// \return returns stored IDCDelta \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  /// \param side TPC side
  /// \param chunk chunk of Delta IDC
  const std::vector<float>& getIDCDeltaUncompressed(const o2::tpc::Side side, const unsigned int chunk) const { return this->mIDCDelta[chunk].getIDCDelta(side); }

  /// \return returns stored IDCDelta \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  /// \param chunk chunk of Delta IDC
  const auto& getIDCDeltaUncompressed(const unsigned int chunk) const { return this->mIDCDelta[chunk]; }

  /// \return creates and returns medium compressed IDCDelta \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  /// \param chunk chunk of Delta IDC
  auto getIDCDeltaMediumCompressed(const unsigned int chunk) const { return IDCDeltaCompressionHelper<short>::getCompressedIDCs(this->mIDCDelta[chunk]); }

  /// \return creates and returns high compressed IDCDelta \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  /// \param chunk chunk of Delta IDC
  auto getIDCDeltaHighCompressed(const unsigned int chunk) const { return IDCDeltaCompressionHelper<char>::getCompressedIDCs(this->mIDCDelta[chunk]); }

  /// \param maxIDCDeltaValue maximum IDC delta value for compressed IDC delta
  static void setMaxCompressedIDCDelta(const float maxIDCDeltaValue) { o2::conf::ConfigurableParam::setValue<float>("TPCIDCCompressionParam", "MaxIDCDeltaValue", maxIDCDeltaValue); }

  /// resetting stored IDC values
  void reset();

  /// \return returns number of chunks for Delta IDCs
  unsigned int getNChunks() const { return mIDCDelta.size(); }

 protected:
  /// check if IDCs are empty
  bool idcsEmpty(const unsigned int cru, const unsigned int tf) const { return this->mIDCs[cru][tf].empty(); }

  /// \return returns the number of stored integration intervals for given Delta IDC chunk
  /// \param chunk chunk of Delta IDC
  unsigned long getNIntegrationIntervals(const unsigned int chunk) const;

  /// \return returns number of IDCs for region and time frame
  auto getSizeIDCs(const unsigned int region, const unsigned int tf) const { return mIDCs[region][tf].size(); }

  /// \return returns index of chunk for DeltaIDC
  /// \param tf time frame
  unsigned int getChunk(const unsigned int tf) const { return tf / mTimeFramesDeltaIDC; }

  /// \return returns number of time frames for chunk
  unsigned int getNTFsPerChunk(const unsigned int chunk) const;

  /// get time frame and index of integrationInterval in the TF
  void getTF(const unsigned int region, unsigned int integrationInterval, unsigned int& timeFrame, unsigned int& interval) const;

  /// \return returns first sector which will be drawn
  unsigned int getFirstSectorDraw() const { return 0; }

  /// \return returns last sector which will be drawn
  unsigned int getLastSectorDraw() const { return Mapper::NSECTORS; }

  /// \return returns first region which will be drawn
  unsigned int getFirstRegionDraw() const { return 0; }

  /// \return returns last region which will be drawn
  unsigned int getLastRegionDraw() const { return Mapper::NREGIONS; }

  const unsigned int mTimeFrames{};                                 ///< number of timeframes which are stored
  const unsigned int mTimeFramesDeltaIDC{};                         ///< number of timeframes of which Delta IDCs are stored
  std::array<std::vector<std::vector<float>>, CRU::MaxCRU> mIDCs{}; ///< grouped and integrated IDCs for the whole TPC. CRU -> time frame -> IDCs
  IDCZero mIDCZero{};                                               ///< I_0(r,\phi) = <I(r,\phi,t)>_t
  IDCOne mIDCOne{};                                                 ///< I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  std::vector<IDCDelta<float>> mIDCDelta{};                         ///< uncompressed: chunk -> Delta IDC: \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  inline static int sNThreads{1};                                   ///< number of threads which are used during the calculations

  ClassDefNV(IDCFactorizationBase, 1)
};

/// template specialization for perform the factorization for one CRU
template <>
class IDCFactorizationBase<IDCGroupHelperRegion> : public IDCGroupHelperRegion
{
 public:
  /// constructor
  /// \param groupPads number of pads in pad direction which will be grouped for all regions
  /// \param groupRows number of pads in row direction which will be grouped for all regions
  /// \param groupLastRowsThreshold minimum number of pads in row direction for the last group in row direction for all regions
  /// \param groupLastPadsThreshold minimum number of pads in pad direction for the last group in pad direction for all regions
  /// \param timeFrames number of time frames which will be stored
  IDCFactorizationBase(const unsigned char groupPads, const unsigned char groupRows, const unsigned char groupLastRowsThreshold, const unsigned char groupLastPadsThreshold, const unsigned int timeFrames, const unsigned int cru) : IDCGroupHelperRegion{groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, cru % Mapper::NREGIONS}, mTimeFrames{timeFrames}, mCRU{cru}, mIDCs(timeFrames){};

  /// default constructor for ROOT I/O
  IDCFactorizationBase() = default;

  /// set the IDC data
  /// \param idcs vector containing the IDCs
  /// \param timeframe time frame of the IDCs
  void setIDCs(std::vector<float>&& idcs, const unsigned int timeframe) { mIDCs[timeframe] = std::move(idcs); }

  /// calculate I_0(r,\phi) = <I(r,\phi,t)>_t
  void calcIDCZero();

  /// calculate I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  void calcIDCOne();

  /// calculate \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  void calcIDCDelta();

  /// \return returns number of IDC values
  auto getIDCsPerCRU(const unsigned int) const { return mNIDCsPerCRU; }

  /// \return returns the stored DeltaIDC value for local ungrouped pad row and ungrouped pad
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param localintegrationInterval local integration interval for chunk (can be obtained with getLocalIntegrationInterval())
  float getIDCDeltaVal(unsigned int urow, unsigned int upad, unsigned int localintegrationInterval) const { return mIDCDelta[getIndexUngrouped(urow, upad, localintegrationInterval)]; }

  /// \return returns the stored IDC0 value for local ungrouped pad row and ungrouped pad
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  float getIDCZeroVal(unsigned int urow, unsigned int upad) const { return mIDCZero[getIndexUngrouped(urow, upad, 0)]; }

  /// \return returns region
  unsigned int getRegion() const { return mCRU % Mapper::NREGIONS; }

  /// \return returns sector
  unsigned int getSector() const { return mCRU / Mapper::NSECTORS; }

  /// \return returns the stored grouped and integrated IDC
  /// \param ugrow row in the region of the ungrouped IDCs
  /// \param gpad pad number of the grouped IDCs
  /// \param interval integration interval
  float getIDCValUngrouped(const unsigned int ugrow, const unsigned int upad, const unsigned int tf, const unsigned int interval) const { return mIDCs[tf][getIndexUngrouped(ugrow, upad, interval)]; }

  /// \return returns total number of integration interval
  unsigned long getNIntegrationIntervals() const;

 protected:
  /// check if IDC values are empty
  bool idcsEmpty(const unsigned int, const unsigned int tf) const { return this->mIDCs[tf].empty(); }

  /// \return returns the stored DeltaIDC value for local ungrouped pad row and ungrouped pad
  float getIDCDeltaVal(const unsigned int, const unsigned int, unsigned int urow, unsigned int upad, unsigned int, unsigned int localintegrationInterval) const { return getIDCDeltaVal(urow, upad, localintegrationInterval); }

  /// \return returns the stored IDC0 value for local ungrouped pad row and ungrouped pad
  float getIDCZeroVal(const unsigned int, const unsigned int, unsigned int urow, unsigned int upad) const { return getIDCZeroVal(urow, upad); }

  /// \return returns the stored grouped and integrated IDC
  float getIDCValUngrouped(const unsigned int, const unsigned int, const unsigned int ugrow, const unsigned int upad, const unsigned int tf, const unsigned int interval) const { return getIDCValUngrouped(ugrow, upad, tf, interval); }

  /// \return returns index of chunk for DeltaIDC (as there are no chunks always return 0)
  unsigned int getChunk(const unsigned int) const { return 0; }

  /// \return returns number of time frames for chunk (as there are no chunks always return the total number of time frames)
  unsigned int getNTFsPerChunk(const unsigned int) const { return mTimeFrames; }

  /// \return returns number of IDCs for region and time frame
  auto getSizeIDCs(const unsigned int, const unsigned int tf) const { return mIDCs[tf].size(); }

  /// \return returns number of chunks for Delta IDCs
  unsigned int getNChunks() const { return 1; }

  /// get time frame and index of integrationInterval in the TF
  void getTF(const unsigned int region, unsigned int integrationInterval, unsigned int& timeFrame, unsigned int& interval) const;

  unsigned int getFirstSectorDraw() const { return getSector(); }
  unsigned int getLastSectorDraw() const { return getFirstSectorDraw() + 1; }
  unsigned int getFirstRegionDraw() const { return getRegion(); }
  unsigned int getLastRegionDraw() const { return getFirstRegionDraw() + 1; }

  const unsigned int mTimeFrames{};        ///< number of timeframes which are stored
  const unsigned int mCRU{};               ///< CRU of which the IDCs are stored
  std::vector<std::vector<float>> mIDCs{}; ///< grouped and integrated IDCs for the whole TPC. CRU -> time frame -> IDCs
  std::vector<float> mIDCZero{};           ///< I_0(r,\phi) = <I(r,\phi,t)>_t
  std::vector<float> mIDCOne{};            ///< I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  std::vector<float> mIDCDelta{};          ///< uncompressed: Delta IDC: \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  inline static int sNThreads{1};          ///< number of threads which are used during the calculations

  ClassDefNV(IDCFactorizationBase, 1)
};

} // namespace o2::tpc

#endif
