// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file IDCFactorization.h
/// \brief class for aggregating IDCs for the full TPC (all sectors) and factorization of aggregated IDCs
///
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>
/// \date Apr 30, 2021

#ifndef ALICEO2_IDCFACTORIZATION_H_
#define ALICEO2_IDCFACTORIZATION_H_

#include <vector>
#include "Rtypes.h"
#include "TPCBase/Mapper.h"
#include "TPCCalibration/IDCGroup.h"
#include "TPCCalibration/IDCFactorizationContainer.h"
#include "DataFormatsTPC/Defs.h"
#include "Framework/Logger.h"

#if (defined(WITH_OPENMP) || defined(_OPENMP)) && !defined(__CLING__)
#include <omp.h>
#endif

namespace o2::tpc
{

/// IDC Delta IDC Compression types
enum class IDCDeltaCompression { NO = 0,     ///< no compression using floats
                                 MEDIUM = 1, ///< medium compression using short (data compression ratio 2 when stored in CCDB)
                                 HIGH = 2    ///< high compression using char (data compression ratio ~5.5 when stored in CCDB)
};

class IDCFactorization
{
 public:
  /// constructor
  /// \param groupPads number of pads in pad direction which will be grouped for all regions
  /// \param groupRows number of pads in row direction which will be grouped for all regions
  /// \param groupLastRowsThreshold minimum number of pads in row direction for the last group in row direction for all regions
  /// \param groupLastPadsThreshold minimum number of pads in pad direction for the last group in pad direction for all regions
  /// \param timeFrames number of timeframes which will be stored
  IDCFactorization(const std::array<unsigned int, Mapper::NREGIONS>& groupPads, const std::array<unsigned int, Mapper::NREGIONS>& groupRows, const std::array<unsigned int, Mapper::NREGIONS>& groupLastRowsThreshold, const std::array<unsigned int, Mapper::NREGIONS>& groupLastPadsThreshold, const unsigned int timeFrames = 0);

  /// default constructor for ROOT I/O
  IDCFactorization() = default;

  /// IDC types
  enum class IDCType { IDC = 0,     ///< integrated and grouped IDCs
                       IDCZero = 1, ///< IDC0: I_0(r,\phi) = <I(r,\phi,t)>_t
                       IDCOne = 2,  ///< IDC1: I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
                       IDCDelta = 3 ///< IDCDelta: \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  };

  /// calculate I_0(r,\phi) = <I(r,\phi,t)>_t
  /// calculate I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  /// calculate \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  void factorizeIDCs();

  /// \return returns the stored grouped and integrated IDC
  /// \param sector sector
  /// \param region region
  /// \param grow row in the region of the grouped IDCs
  /// \param gpad pad number of the grouped IDCs
  /// \param integrationInterval integration interval
  float operator()(const unsigned int sector, const unsigned int region, unsigned int grow, unsigned int gpad, unsigned int integrationInterval) const { return mIDCs[sector * Mapper::NREGIONS + region][integrationInterval][mOffsRow[region][grow] + gpad]; }

  /// \return returns grouped pad for ungrouped row and pad
  /// \param region region
  /// \param urow local ungrouped row in a region
  /// \param upad ungrouped pad
  unsigned int getGroupedPad(const unsigned int region, unsigned int urow, unsigned int upad) const { return IDCGroup::getGroupedPad(upad, urow, region, mGroupPads[region], mGroupRows[region], mRows[region], mPadsPerRow[region]); }

  /// \return returns the row of the group from the local ungrouped row in a region
  /// \param region region
  /// \param urow local ungrouped row in a region
  unsigned int getGroupedRow(const unsigned int region, unsigned int urow) const { return IDCGroup::getGroupedRow(urow, mGroupRows[region], mRows[region]); }

  /// \return returns the stored value for local ungrouped pad row and ungrouped pad
  /// \param sector sector
  /// \param region region
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  float getIDCVal(const unsigned int sector, const unsigned int region, unsigned int urow, unsigned int upad, unsigned int integrationInterval) const;

  /// \return returns the stored IDC0 value for local ungrouped pad row and ungrouped pad
  /// \param sector sector
  /// \param region region
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  float getIDCZeroVal(const unsigned int sector, const unsigned int region, unsigned int urow, unsigned int upad) const { return mIDCZeroOne.getValueIDCZero(Sector(sector).side(), getIndexUngrouped(sector, region, urow, upad, 0)); }

  /// \return returns the stored DeltaIDC value for local ungrouped pad row and ungrouped pad
  /// \param sector sector
  /// \param region region
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  float getIDCDeltaVal(const unsigned int sector, const unsigned int region, unsigned int urow, unsigned int upad, unsigned int integrationInterval) const { return mIDCDelta.getValue(Sector(sector).side(), getIndexUngrouped(sector, region, urow, upad, integrationInterval)); }

  /// \return returns the index to the grouped data with ungrouped inputs
  /// \param sector sector
  /// \param region TPC region
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  unsigned int getIndexUngrouped(const unsigned int sector, const unsigned int region, unsigned int urow, unsigned int upad, unsigned int integrationInterval) const { return getIndexGrouped(sector % o2::tpc::SECTORSPERSIDE, region, getGroupedRow(region, urow), getGroupedPad(region, urow, upad), integrationInterval); }

  /// \returns grouping definition in pad direction (How many pads are grouped)
  /// \param region TPC region
  unsigned int getGroupPads(const unsigned int region) const { return mGroupPads[region]; }

  /// \returns grouping definition in pad direction (How many pads are grouped)
  const auto& getGroupPads() const { return mGroupPads; }

  /// \returns grouping definition in row direction (How many rows are grouped)
  /// \param region TPC region
  unsigned int getGroupRows(const unsigned int region) const { return mGroupRows[region]; }

  /// \returns grouping definition in row direction (How many rows are grouped)
  const auto& getGroupRows() const { return mGroupRows; }

  /// \returns grouping threshold for last group in pad direction
  /// \param region TPC region
  unsigned int getPadThreshold(const unsigned int region) const { return mGroupLastPadsThreshold[region]; }

  /// \returns grouping threshold for last group in pad direction
  const auto& getPadThreshold() const { return mGroupLastPadsThreshold; }

  /// \returns grouping threshold for last group in row direction
  /// \param region TPC region
  unsigned int getRowThreshold(const unsigned int region) const { return mGroupLastRowsThreshold[region]; }

  /// \returns grouping threshold for last group in row direction
  const auto& getRowThreshold() const { return mGroupLastRowsThreshold; }

  /// \return returns number of timeframes for which the IDCs are stored
  unsigned int getNTimeframes() const { return mTimeFrames; }

  /// \return returns total number of IDCs per region per integration interval
  /// \param region TPC region
  unsigned int getIDCsPerCRU(const unsigned int region) const { return mNIDCsPerCRU[region]; }

  /// \return returns number of grouped IDCs per sector
  unsigned int getIDCsPerSector() const { return mNIDCsPerSector; }

  /// \return returns the number of stored integration intervals
  unsigned long getNIntegrationIntervals() const;

  /// \return returns stored IDC0 I_0(r,\phi) = <I(r,\phi,t)>_t
  /// \param side TPC side
  const std::vector<float>& getIDCZero(const o2::tpc::Side side) const { return mIDCZeroOne.mIDCZero[side]; }

  /// \return returns stored IDC1 I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  /// \param side TPC side
  const std::vector<float>& getIDCOne(const o2::tpc::Side side) const { return mIDCZeroOne.mIDCOne[side]; }

  /// \return returns stored IDCDelta \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  /// \param side TPC side
  const std::vector<float>& getIDCDeltaUncompressed(const o2::tpc::Side side) const { return mIDCDelta.mIDCDelta[side]; }

  /// \return returns stored IDCDelta \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  const auto& getIDCDeltaUncompressed() const { return mIDCDelta; }

  /// \return creates and returns medium compressed IDCDelta \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  auto getIDCDeltaMediumCompressed() const { return IDCDeltaCompressionHelper<short>::getCompressedIDCs(mIDCDelta); }

  /// \return creates and returns high compressed IDCDelta \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  auto getIDCDeltaHighCompressed() const { return IDCDeltaCompressionHelper<char>::getCompressedIDCs(mIDCDelta); }

  /// \return returns struct containing IDC0 and IDC1
  const auto& getIDCZeroOne() const { return mIDCZeroOne; }

  /// \return returns grouped IDCs
  const auto& getIDCs() const { return mIDCs; }

  /// \return returns the number of threads used for some of the calculations
  static int getNThreads() { return sNThreads; }

  /// set the IDC data
  /// \param idcs vector containing the IDCs
  /// \param cru CRU
  /// \param timeframe time frame of the IDCs
  void setIDCs(std::vector<float>&& idcs, const unsigned int cru, const unsigned int timeframe) { mIDCs[cru][timeframe] = std::move(idcs); }

  /// set the number of threads used for some of the calculations
  /// \param nThreads number of threads
  static void setNThreads(const int nThreads) { sNThreads = nThreads; }

  /// \param maxIDCDeltaValue maximum IDC delta value for compressed IDC delta
  static void setMaxCompressedIDCDelta(const float maxIDCDeltaValue) { o2::conf::ConfigurableParam::setValue<float>("TPCIDCCompressionParam", "MaxIDCDeltaValue", maxIDCDeltaValue); }

  /// draw IDCs for one sector for one integration interval
  /// \param sector sector which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawIDCsSector(const unsigned int sector, const unsigned int integrationInterval, const std::string filename = "IDCsSector.pdf") const { drawSector(IDCType::IDC, sector, integrationInterval, filename); }

  /// draw IDC zero I_0(r,\phi) = <I(r,\phi,t)>_t
  /// \param sector sector which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawIDCZeroSector(const unsigned int sector, const std::string filename = "IDCZeroSector.pdf") const { drawSector(IDCType::IDCZero, sector, 0, filename); }

  /// draw IDCDelta for one sector for one integration interval
  /// \param sector sector which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param compression compression of Delta IDCs. (setMaxCompressedIDCDelta() should be called first in case of non standard compression parameter)
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawIDCDeltaSector(const unsigned int sector, const unsigned int integrationInterval, const IDCDeltaCompression compression, const std::string filename = "IDCDeltaSector.pdf") const { drawSector(IDCType::IDCDelta, sector, integrationInterval, filename, compression); }

  /// draw IDCs for one side for one integration interval
  /// \param side side which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawIDCsSide(const o2::tpc::Side side, const unsigned int integrationInterval, const std::string filename = "IDCsSide.pdf") const { drawSide(IDCType::IDC, side, integrationInterval, filename); }

  /// draw IDC zero I_0(r,\phi) = <I(r,\phi,t)>_t
  /// \param side side which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawIDCZeroSide(const o2::tpc::Side side, const std::string filename = "IDCZeroSide.pdf") const { drawSide(IDCType::IDCZero, side, 0, filename); }

  /// draw IDCDelta for one side for one integration interval
  /// \param side side which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param compression compression of Delta IDCs. (setMaxCompressedIDCDelta() should be called first in case of non standard compression parameter)
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawIDCDeltaSide(const o2::tpc::Side side, const unsigned int integrationInterval, const IDCDeltaCompression compression, const std::string filename = "IDCDeltaSide.pdf") const { drawSide(IDCType::IDCDelta, side, integrationInterval, filename, compression); }

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName = "IDCFactorized.root", const char* outName = "IDCFactorized") const;

  /// \param integrationIntervals number of integration intervals which will be dumped to the tree (-1: all integration intervalls)
  void dumpIDCsToTree(int integrationIntervals = -1) const;

 private:
  const std::array<unsigned int, Mapper::NREGIONS> mGroupPads{};              ///< grouping definition in pad direction (How many pads are grouped)
  const std::array<unsigned int, Mapper::NREGIONS> mGroupRows{};              ///< grouping definition in row direction (How many rows are grouped)
  const std::array<unsigned int, Mapper::NREGIONS> mGroupLastRowsThreshold{}; ///< if the last group (region edges) consists in row direction of less then mGroupLastRowsThreshold pads then it will be grouped into the previous group
  const std::array<unsigned int, Mapper::NREGIONS> mGroupLastPadsThreshold{}; ///< if the last group (sector edges) consists in pad direction of less then mGroupLastPadsThreshold pads then it will be grouped into the previous group
  const unsigned int mTimeFrames{};                                           ///< number of timeframes which are stored
  std::array<unsigned int, Mapper::NREGIONS> mNIDCsPerCRU{1};                 ///< total number of IDCs per region per integration interval
  std::array<std::vector<std::vector<float>>, CRU::MaxCRU> mIDCs{};           ///< grouped and integrated IDCs for the whole TPC. CRU -> time frame -> IDCs
  IDCZeroOne mIDCZeroOne{};                                                   ///< I_0(r,\phi) = <I(r,\phi,t)>_t and I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  IDCDelta<float> mIDCDelta{};                                                ///< uncompressed: \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  unsigned int mNIDCsPerSector{};                                             ///< number of grouped IDCs per sector
  std::array<unsigned int, Mapper::NREGIONS> mRows{};                         ///< number of grouped rows per region
  std::array<unsigned int, Mapper::NREGIONS> mRegionOffs{};                   ///< offset for the region per region
  std::array<std::vector<unsigned int>, Mapper::NREGIONS> mPadsPerRow{};      ///< number of pads per row per region
  std::array<std::vector<unsigned int>, Mapper::NREGIONS> mOffsRow{};         ///< offset to calculate the index in the data from row and pad per region
  inline static int sNThreads{1};                                             ///< number of threads which are used during the calculations

  /// calculate I_0(r,\phi) = <I(r,\phi,t)>_t
  void calcIDCZero();

  /// calculate I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  void calcIDCOne();

  /// calculate \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  void calcIDCDelta();

  /// \return returns index to the data
  /// \param row row of the grouped IDCs
  /// \param pad pad of the grouped IDCs
  unsigned int getIndexGrouped(const unsigned int sector, const unsigned int region, const unsigned int row, const unsigned int pad, unsigned int integrationInterval) const { return mNIDCsPerSector * (integrationInterval * SECTORSPERSIDE + sector) + mRegionOffs[region] + mOffsRow[region][row] + pad; }

  /// draw IDCs for one sector for one integration interval
  /// \param sector sector which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawSector(const IDCType type, const unsigned int sector, const unsigned int integrationInterval, const std::string filename, const IDCDeltaCompression compression = IDCDeltaCompression::NO) const;

  /// draw IDCs for one side for one integration interval
  /// \param Side side which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawSide(const IDCType type, const o2::tpc::Side side, const unsigned int integrationInterval, const std::string filename, const IDCDeltaCompression compression = IDCDeltaCompression::NO) const;

  /// get z axis title for given IDC type and compression type
  std::string getZAxisTitle(const IDCType type, const IDCDeltaCompression compression) const;

  void getTF(const unsigned int region, unsigned int integrationInterval, unsigned int& timeFrame, unsigned int& interval) const;

  ClassDefNV(IDCFactorization, 1)
};

} // namespace o2::tpc

#endif
