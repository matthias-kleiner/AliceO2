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

/// \file IDCFactorization.h
/// \brief class for factorization of IDCs for the full TPC (all sectors) or one CRU
///
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>
/// \date Apr 30, 2021

#ifndef ALICEO2_IDCFACTORIZATION_H_
#define ALICEO2_IDCFACTORIZATION_H_

#include <vector>
#include "Rtypes.h"
#include "TPCBase/Mapper.h"
#include "TPCCalibration/IDCFactorizationBase.h"
#include "DataFormatsTPC/Defs.h"

namespace o2::tpc
{

template <typename Type = IDCGroupHelperSector>
class IDCFactorization : public IDCFactorizationBase<Type>
{

 public:
  /// constructor
  /// \param groupPads number of pads in pad direction which will be grouped for all regions
  /// \param groupRows number of pads in row direction which will be grouped for all regions
  /// \param groupLastRowsThreshold minimum number of pads in row direction for the last group in row direction for all regions
  /// \param groupLastPadsThreshold minimum number of pads in pad direction for the last group in pad direction for all regions
  /// \param timeFrames number of time frames which will be stored
  /// \param timeframesDeltaIDC number of time frames stored for each DeltaIDC object
  template <bool IsEnabled = true, typename std::enable_if<(IsEnabled && (std::is_same<Type, IDCGroupHelperSector>::value)), int>::type = 0>
  IDCFactorization(const std::array<unsigned char, Mapper::NREGIONS>& groupPads, const std::array<unsigned char, Mapper::NREGIONS>& groupRows, const std::array<unsigned char, Mapper::NREGIONS>& groupLastRowsThreshold, const std::array<unsigned char, Mapper::NREGIONS>& groupLastPadsThreshold, const unsigned int timeFrames, const unsigned int timeframesDeltaIDC) : IDCFactorizationBase<IDCGroupHelperSector>(groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, timeFrames, timeframesDeltaIDC){};

  /// constructor
  /// \param groupPads number of pads in pad direction which will be grouped for all regions
  /// \param groupRows number of pads in row direction which will be grouped for all regions
  /// \param groupLastRowsThreshold minimum number of pads in row direction for the last group in row direction for all regions
  /// \param groupLastPadsThreshold minimum number of pads in pad direction for the last group in pad direction for all regions
  /// \param timeFrames number of time frames which will be stored
  template <bool IsEnabled = true, typename std::enable_if<(IsEnabled && (std::is_same<Type, IDCGroupHelperRegion>::value)), int>::type = 0>
  IDCFactorization(const unsigned char groupPads, const unsigned char groupRows, const unsigned char groupLastRowsThreshold, const unsigned char groupLastPadsThreshold, const unsigned int timeFrames, const unsigned int cru) : IDCFactorizationBase<IDCGroupHelperRegion>(groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, timeFrames, cru){};

  /// default constructor for ROOT I/O
  IDCFactorization() = default;

  /// calculate I_0(r,\phi) = <I(r,\phi,t)>_t
  /// calculate I_1(t) = <I(r,\phi,t) / I_0(r,\phi)>_{r,\phi}
  /// calculate \Delta I(r,\phi,t) = I(r,\phi,t) / ( I_0(r,\phi) * I_1(t) )
  void factorizeIDCs();

  /// \return returns the stored value for local ungrouped pad row and ungrouped pad
  /// \param sector sector
  /// \param region region
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  float getIDCValUngrouped(const unsigned int sector, const unsigned int region, unsigned int urow, unsigned int upad, unsigned int integrationInterval) const;

  /// \return returns index of integration interval in the chunk from global integration interval
  /// \param region TPC region
  /// \param integrationInterval integration interval
  /// \param chunk which will be set in the function
  /// \param localintegrationInterval local integration interval for chunk which will be set in the function
  void getLocalIntegrationInterval(const unsigned int region, const unsigned int integrationInterval, unsigned int& chunk, unsigned int& localintegrationInterval) const;

  /// \return returns number of timeframes for which the IDCs are stored
  unsigned int getNTimeframes() const { return this->mTimeFrames; }

  /// \return returns struct containing IDC0
  const auto& getIDCZero() const { return this->mIDCZero; }

  /// \return returns struct containing IDC1
  const auto& getIDCOne() const& { return this->mIDCOne; }

  /// \return returns struct containing IDC1 using move semantics
  auto getIDCOne() && { return std::move(this->mIDCOne); }

  /// \return returns grouped IDCs
  const auto& getIDCs() const { return this->mIDCs; }

  /// \return returns the number of threads used for some of the calculations
  static int getNThreads() { return o2::tpc::IDCFactorizationBase<Type>::sNThreads; }

  /// set the number of threads used for some of the calculations
  /// \param nThreads number of threads
  static void setNThreads(const int nThreads) { o2::tpc::IDCFactorizationBase<Type>::sNThreads = nThreads; }

  /// draw IDCs for one sector for one integration interval
  /// \param sector sector which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  template <bool IsEnabled = true, typename std::enable_if<(IsEnabled && (std::is_same<Type, IDCGroupHelperSector>::value)), int>::type = 0>
  void drawIDCsSector(const unsigned int sector, const unsigned int integrationInterval, const std::string filename = "IDCsSector.pdf") const
  {
    drawIDCHelper(false, Sector(sector), integrationInterval, filename);
  }

  /// draw IDCs for one region for one integration interval
  /// \param integrationInterval which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  template <bool IsEnabled = true, typename std::enable_if<(IsEnabled && (std::is_same<Type, IDCGroupHelperRegion>::value)), int>::type = 0>
  void drawIDCsRegion(const unsigned int integrationInterval, const std::string filename = "IDCsRegion.pdf") const
  {
    drawIDCHelper(false, this->getSector(), integrationInterval, filename);
  }

  /// draw IDC zero I_0(r,\phi) = <I(r,\phi,t)>_t
  /// \param sector sector which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  template <bool IsEnabled = true, typename std::enable_if<(IsEnabled && (std::is_same<Type, IDCGroupHelperSector>::value)), int>::type = 0>
  void drawIDCZeroSector(const unsigned int sector, const std::string filename = "IDCZeroSector.pdf") const
  {
    drawIDCZeroHelper(false, Sector(sector), filename);
  }

  /// draw IDC zero I_0(r,\phi) = <I(r,\phi,t)>_t
  /// \param filename name of the output file. If empty the canvas is drawn.
  template <bool IsEnabled = true, typename std::enable_if<(IsEnabled && (std::is_same<Type, IDCGroupHelperRegion>::value)), int>::type = 0>
  void drawIDCZeroRegion(const std::string filename = "IDCZeroRegion.pdf") const
  {
    drawIDCZeroHelper(false, this->getSector(), filename);
  }

  /// draw IDCDelta for one sector for one integration interval
  /// \param sector sector which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param compression compression of Delta IDCs. (setMaxCompressedIDCDelta() should be called first in case of non standard compression parameter)
  /// \param filename name of the output file. If empty the canvas is drawn.
  template <bool IsEnabled = true, typename std::enable_if<(IsEnabled && (std::is_same<Type, IDCGroupHelperSector>::value)), int>::type = 0>
  void drawIDCDeltaSector(const unsigned int sector, const unsigned int integrationInterval, const IDCDeltaCompression compression, const std::string filename = "IDCDeltaSector.pdf") const
  {
    drawIDCDeltaHelper(false, Sector(sector), integrationInterval, compression, filename);
  }

  /// draw IDCDelta for one sector for one integration interval
  /// \param integrationInterval which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  template <bool IsEnabled = true, typename std::enable_if<(IsEnabled && (std::is_same<Type, IDCGroupHelperRegion>::value)), int>::type = 0>
  void drawIDCDeltaRegion(const unsigned int integrationInterval, const std::string filename = "IDCDeltaRegion.pdf") const
  {
    drawIDCDeltaHelper(false, integrationInterval, IDCDeltaCompression::NO, filename);
  }

  /// draw IDCs for one side for one integration interval
  /// \param side side which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  template <bool IsEnabled = true, typename std::enable_if<(IsEnabled && (std::is_same<Type, IDCGroupHelperSector>::value)), int>::type = 0>
  void drawIDCsSide(const o2::tpc::Side side, const unsigned int integrationInterval, const std::string filename = "IDCsSide.pdf") const
  {
    drawIDCHelper(true, side == Side::A ? Sector(0) : Sector(Sector::MAXSECTOR - 1), integrationInterval, filename);
  }

  /// draw IDC zero I_0(r,\phi) = <I(r,\phi,t)>_t
  /// \param side side which will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  template <bool IsEnabled = true, typename std::enable_if<(IsEnabled && (std::is_same<Type, IDCGroupHelperSector>::value)), int>::type = 0>
  void drawIDCZeroSide(const o2::tpc::Side side, const std::string filename = "IDCZeroSide.pdf") const
  {
    drawIDCZeroHelper(true, side == Side::A ? Sector(0) : Sector(Sector::MAXSECTOR - 1), filename);
  }

  /// draw IDCDelta for one side for one integration interval
  /// \param side side which will be drawn
  /// \param integrationInterval which will be drawn
  /// \param compression compression of Delta IDCs. (setMaxCompressedIDCDelta() should be called first in case of non standard compression parameter)
  /// \param filename name of the output file. If empty the canvas is drawn.
  template <bool IsEnabled = true, typename std::enable_if<(IsEnabled && (std::is_same<Type, IDCGroupHelperSector>::value)), int>::type = 0>
  void drawIDCDeltaSide(const o2::tpc::Side side, const unsigned int integrationInterval, const IDCDeltaCompression compression, const std::string filename = "IDCDeltaSide.pdf") const
  {
    drawIDCDeltaHelper(true, side == Side::A ? Sector(0) : Sector(Sector::MAXSECTOR - 1), integrationInterval, compression, filename);
  }

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName = "IDCFactorized.root", const char* outName = "IDCFactorized") const;

  /// \param integrationIntervals number of integration intervals which will be dumped to the tree (-1: all integration intervalls)
  /// \param outFileName name of the output file
  void dumpToTree(int integrationIntervals = -1, const char* outFileName = "IDCTree.root") const;

  /// \returns vector containing the number of integration intervals for each stored TF
  std::vector<unsigned int> getIntegrationIntervalsPerTF(const unsigned int region = 0) const;

 private:
  /// helper function for drawing IDCDelta
  void drawIDCDeltaHelper(const bool type, const Sector sector, const unsigned int integrationInterval, const IDCDeltaCompression compression, const std::string filename) const;

  /// helper function for drawing IDCs
  void drawIDCHelper(const bool type, const Sector sector, const unsigned int integrationInterval, const std::string filename) const;

  /// helper function for drawing IDCZero
  void drawIDCZeroHelper(const bool type, const Sector sector, const std::string filename) const;

  /// get time frame and index of integrationInterval in the TF
  void getTF(const unsigned int region, unsigned int integrationInterval, unsigned int& timeFrame, unsigned int& interval) const;

  ClassDefNV(IDCFactorization, 1)
};

} // namespace o2::tpc

#endif
