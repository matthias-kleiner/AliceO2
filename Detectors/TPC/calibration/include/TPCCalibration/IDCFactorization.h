// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file IDCDelta.h
/// \brief class for storing the aggregated IDCs for the full TPC
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

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

namespace o2
{
namespace tpc
{
class IDCFactorization
{
 public:
  IDCFactorization(const std::array<unsigned int, Mapper::NREGIONS>& groupPads, const std::array<unsigned int, Mapper::NREGIONS>& groupRows, const std::array<unsigned int, Mapper::NREGIONS>& groupLastRowsThreshold, const std::array<unsigned int, Mapper::NREGIONS>& groupLastPadsThreshold, const unsigned int timeFrames = 0);

  /// adding default constructor for ROOT I/O
  IDCFactorization() = default;

  /// IDC types
  enum class IDCType { IDC = 0,     ///< integrated and grouped IDCs
                       IDCZero = 1, ///< IDC0: I_0(r,\phi) = <I(r,\phi,t)>_t
                       IDCOne = 2,  ///< IDC1: I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
                       IDCDelta = 3 ///< IDCDelta: \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  };

  /// \return returns the stored IDC value
  /// \param sector sector
  /// \param region region
  /// \param row row of the grouped IDCs
  /// \param pad pad number of the grouped IDCs
  /// \param integrationInterval integration interval
  const float& operator()(const unsigned int sector, const unsigned int region, unsigned int row, unsigned int pad, unsigned int integrationInterval) const
  {
    return mIDCs[sector * Mapper::NREGIONS + region][integrationInterval][mOffsRow[region][row] + pad];
  }

  int getGroupedPad(const unsigned int region, unsigned int urow, unsigned int upad) const
  {
    return IDCGroup::getGroupedPad(upad, urow, region, mGroupPads[region], mGroupRows[region], mRows[region], mPadsPerRow[region]);
  }

  unsigned int getGroupedRow(const unsigned int region, unsigned int urow) const
  {
    return IDCGroup::getGroupedRow(urow, mGroupRows[region], mRows[region]);
  }

  /// \return returns the stored value for local ungrouped pad row and ungrouped pad
  /// \param sector sector
  /// \param region region
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  const float& getIDCVal(const unsigned int sector, const unsigned int region, unsigned int urow, unsigned int upad, unsigned int integrationInterval) const
  {
    unsigned int timeFrame = 0;
    unsigned int interval = 0;
    unsigned int nintervals = 0;

    for (unsigned int tf = 0; tf < mTimeFrames; ++tf) {
      timeFrame = tf;
      nintervals += mIDCs[region][tf].size() / mNIDCsPerCRU[region];
      if (integrationInterval < nintervals) {
        timeFrame = tf;
        interval = integrationInterval - interval;
        break;
      }
      interval += nintervals;
    }
    return mIDCs[sector * Mapper::NREGIONS + region][timeFrame][interval * mNIDCsPerCRU[region] + mOffsRow[region][getGroupedRow(region, urow)] + getGroupedPad(region, urow, upad)];
  }

  /// \return returns the stored value for local ungrouped pad row and ungrouped pad
  /// \param cru cru
  /// \param timeframe timeframe
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationIntervalinTF integration interval in timeframe
  // const float& getIDCVal(const unsigned int cru, const unsigned int timeframe, unsigned int urow, unsigned int upad, unsigned int integrationIntervalinTF) const
  // {
  //   const unsigned int region = cru % Mapper::NREGIONS;
  //   return mIDCs[cru][timeframe][integrationIntervalinTF * mNIDCsPerCRU[region] + mOffsRow[region][getGroupedRow(region, urow)] + getGroupedPad(region, urow, upad)];
  // }

  /// \return returns the stored value for local ungrouped pad row and ungrouped pad
  /// \param sector sector
  /// \param region region
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  const float& getIDCZeroVal(const unsigned int sector, const unsigned int region, unsigned int urow, unsigned int upad) const
  {
    return mIDCZeroOne.getValueIDCOne(getSide(sector), getIndex(sector % o2::tpc::SECTORSPERSIDE, region, getGroupedRow(region, urow), getGroupedPad(region, urow, upad), 0));
  }

  /// \return returns the stored value for local ungrouped pad row and ungrouped pad
  /// \param sector sector
  /// \param region region
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  const float& getIDCDeltaVal(const unsigned int sector, const unsigned int region, unsigned int urow, unsigned int upad, unsigned int integrationInterval) const
  {
    return mIDCDelta.getValue(getSide(sector), getIndex(sector % o2::tpc::SECTORSPERSIDE, region, getGroupedRow(region, urow), getGroupedPad(region, urow, upad), integrationInterval));
  }

  /// \param idcs vector containing the IDCs
  /// \param cru CRU
  /// \param timeframe time frame of the IDCs
  void setIDCs(std::vector<float>&& idcs, const unsigned int cru, const unsigned int timeframe)
  {
    mIDCs[cru][timeframe] = std::move(idcs);
  }

  /// calculate I_0(r,\phi) = <I(r,\phi,t)>_t
  /// calculate I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  /// calculate \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  void factorizeIDCs()
  {
    LOGP(info, "Using {} threads for factorization of IDCs", sNThreads);
    calcIDCZero();
    calcIDCOne();
    calcIDCDelta();
  }

  /// \returns grouping definition in pad direction (How many pads are grouped)
  unsigned int getGroupPads(const unsigned int region) const { return mGroupPads[region]; }

  /// \returns grouping definition in row direction (How many rows are grouped)
  unsigned int getGroupRows(const unsigned int region) const { return mGroupRows[region]; }

  /// \returns grouping threshold for last group in pad direction
  unsigned int getPadThreshold(const unsigned int region) const { return mGroupLastPadsThreshold[region]; }

  /// \returns grouping threshold for last group in row direction
  unsigned int getRowThreshold(const unsigned int region) const { return mGroupLastRowsThreshold[region]; }

  /// \return returns number of timeframes for which the IDCs are stored
  unsigned int getNTimeframes() const { return mTimeFrames; }

  /// \return returns total number of IDCs per region per integration interval
  unsigned int getIDCsPerCRU(const unsigned int region) const { return mNIDCsPerCRU[region]; }

  /// \return returns number of grouped IDCs per sector
  unsigned int getIDCsPerSector() const { return mNIDCsPerSector; }

  /// \return returns the number of integration intervals stored
  unsigned long getNIntegrationIntervals() const
  {
    std::size_t sum = 0;
    for (auto&& idcsTF : mIDCs[0]) {
      sum += idcsTF.size();
    }
    return sum / mNIDCsPerCRU[0];
  }

  /// \return returns stored IDC0 I_0(r,\phi) = <I(r,\phi,t)>_t
  const std::vector<float>& getIDCZero(const o2::tpc::Side side) const { return mIDCZeroOne.mIDCZero[side]; }

  /// \return returns stored IDC1 I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  const std::vector<float>& getIDCOne(const o2::tpc::Side side) const { return mIDCZeroOne.mIDCOne[side]; }

  /// \return returns stored IDCDelta \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  const auto& getIDCDelta(const o2::tpc::Side side) const { return mIDCDelta.mIDCDelta[side]; }

  const auto& getIDCZeroOne() const { return mIDCZeroOne; }
  const auto& getIDCDelta() const { return mIDCDelta; }

  /// \return returns grouped IDCs
  const auto& getIDCs() const { return mIDCs; }

  /// get the number of threads used for some of the calculations
  static int getNThreads() { return sNThreads; }

  /// set the number of threads used for some of the calculations
  static void setNThreads(const int nThreads)
  {
    sNThreads = nThreads;
  }

  /// draw IDCs for one sector for one integration interval
  /// \param sector sector which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawIDCsSector(const unsigned int sector, const unsigned int integrationInterval, const std::string filename = "IDCsSector.pdf") const
  {
    drawSector(IDCType::IDC, sector, integrationInterval, filename);
  }

  /// draw IDC zero I_0(r,\phi) = <I(r,\phi,t)>_t
  /// \param sector sector which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawIDCZeroSector(const unsigned int sector, const std::string filename = "IDCZeroSector.pdf") const
  {
    drawSector(IDCType::IDCZero, sector, 0, filename);
  }

  /// draw IDCDelta for one sector for one integration interval
  /// \param sector sector which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawIDCDeltaSector(const unsigned int sector, const unsigned int integrationInterval, const std::string filename = "IDCDeltaSector.pdf") const
  {
    drawSector(IDCType::IDCDelta, sector, integrationInterval, filename);
  }

  /// draw IDCs for one side for one integration interval
  /// \param Side side which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawIDCsSide(const o2::tpc::Side side, const unsigned int integrationInterval, const std::string filename = "IDCsSide.pdf") const
  {
    drawSide(IDCType::IDC, side, integrationInterval, filename);
  }

  /// draw IDC zero I_0(r,\phi) = <I(r,\phi,t)>_t
  /// \param Side side which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawIDCZeroSide(const o2::tpc::Side side, const std::string filename = "IDCZeroSide.pdf") const
  {
    drawSide(IDCType::IDCZero, side, 0, filename);
  }

  /// draw IDCDelta for one side for one integration interval
  /// \param Side side which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawIDCDeltaSide(const o2::tpc::Side side, const unsigned int integrationInterval, const std::string filename = "IDCZeroSide.pdf") const
  {
    drawSide(IDCType::IDCDelta, side, integrationInterval, filename);
  }

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName = "IDCFactorized.root", const char* outName = "IDCFactorized") const;

  void dumpIDCsToTree(int integrationIntervals) const;

 private:
  const std::array<unsigned int, Mapper::NREGIONS> mGroupPads{};                            ///< grouping definition in pad direction (How many pads are grouped)
  const std::array<unsigned int, Mapper::NREGIONS> mGroupRows{};                            ///< grouping definition in row direction (How many rows are grouped)
  const std::array<unsigned int, Mapper::NREGIONS> mGroupLastRowsThreshold{};               ///< if the last group (region edges) consists in row direction of less then mGroupLastRowsThreshold pads then it will be grouped into the previous group
  const std::array<unsigned int, Mapper::NREGIONS> mGroupLastPadsThreshold{};               ///< if the last group (sector edges) consists in pad direction of less then mGroupLastPadsThreshold pads then it will be grouped into the previous group
  const unsigned int mTimeFrames{};                                                         ///< number of timeframes which are stored
  std::array<unsigned int, Mapper::NREGIONS> mNIDCsPerCRU{1};                               ///< total number of IDCs per region per integration interval
  std::array<std::vector<std::vector<float>>, Mapper::NSECTORS * Mapper::NREGIONS> mIDCs{}; ///< grouped and IDCs for the whole TPC. sector -> region -> time frame -> IDCs
  // std::array<std::vector<float>, o2::tpc::SIDES> mIDCZero{};                                ///< I_0(r,\phi) = <I(r,\phi,t)>_t
  // std::array<std::vector<float>, o2::tpc::SIDES> mIDCOne{};                                 ///< I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  // std::array<std::vector<float>, o2::tpc::SIDES> mIDCDelta{};                               ///< \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  IDCZeroOne mIDCZeroOne{};                                              ///< I_0(r,\phi) = <I(r,\phi,t)>_t and I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  IDCDelta mIDCDelta{};                                                  ///< \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  unsigned int mNIDCsPerSector{};                                        ///< number of grouped IDCs per sector
  std::array<unsigned int, Mapper::NREGIONS> mRows{};                    ///< number of grouped rows per region
  std::array<unsigned int, Mapper::NREGIONS> mRegionOffs{};              ///< offset for the region per region
  std::array<std::vector<unsigned int>, Mapper::NREGIONS> mPadsPerRow{}; ///< number of pads per row per region
  std::array<std::vector<unsigned int>, Mapper::NREGIONS> mOffsRow{};    ///< offset to calculate the index in the data from row and pad per region
  inline static int sNThreads{1};                                        ///< number of threads which are used during the calculations

  /// calculate I_0(r,\phi) = <I(r,\phi,t)>_t
  void calcIDCZero();

  /// calculate I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  void calcIDCOne();

  /// calculate \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  void calcIDCDelta();

  /// \return returns index to the data
  /// \param row row of the grouped IDCs
  /// \param pad pad of the grouped IDCs
  unsigned int getIndex(const unsigned int sector, const unsigned int region, const unsigned int row, const unsigned int pad, unsigned int integrationInterval) const { return mNIDCsPerSector * (integrationInterval * Mapper::NSECTORS + sector) + mRegionOffs[region] + mOffsRow[region][row] + pad; }

  /// draw IDCs for one sector for one integration interval
  /// \param sector sector which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawSector(const IDCType type, const unsigned int sector, const unsigned int integrationInterval, const std::string filename) const;

  /// draw IDCs for one side for one integration interval
  /// \param Side side which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawSide(const IDCType type, const o2::tpc::Side side, const unsigned int integrationInterval, const std::string filename) const;

  Side getSide(const unsigned int sector) const { return sector < SECTORSPERSIDE ? Side::A : Side::C; }

  ClassDefNV(IDCFactorization, 1)
};

} // namespace tpc
} // namespace o2

#endif
