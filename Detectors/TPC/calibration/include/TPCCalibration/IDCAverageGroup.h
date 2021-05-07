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
/// \brief class for averaging/grouping and storing of IDCs
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_IDCAVERAGEGROUP_H_
#define ALICEO2_IDCAVERAGEGROUP_H_

#include <vector>
#include "TPCCalibration/IDCGroup.h"
#include "TPCBase/Mapper.h"

#if (defined(WITH_OPENMP) || defined(_OPENMP)) && !defined(__CLING__)
#include <omp.h>
#endif

namespace o2
{
namespace tpc
{
class IDCAverageGroup
{
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

 public:
  /// constructor
  /// \param groupPads number of pads in pad direction which will be grouped
  /// \param groupRows number of pads in row direction which will be grouped
  /// \param groupLastRowsThreshold minimum number of pads in row direction for the last group in row direction
  /// \param groupLastPadsThreshold minimum number of pads in pad direction for the last group in pad direction
  /// \param region region of the TPC
  IDCAverageGroup(const unsigned int groupPads = 4, const unsigned int groupRows = 4, const unsigned int groupLastRowsThreshold = 2, const unsigned int groupLastPadsThreshold = 2, const unsigned int region = 0)
    : mIDCsGrouped{groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, region} {}

  /// set the IDCs which will be averaged and grouped
  /// \param idcs vector containing the IDCs
  void setIDCs(const std::vector<float>& idcs)
  {
    mIDCsUngrouped = idcs;
    mIDCsGrouped.resize(getNIntegrationIntervals());
  }

  /// set the IDCs which will be averaged and grouped using move operator
  /// \param IDCs vector containing the IDCs
  void setIDCs(std::vector<float>&& idcs)
  {
    mIDCsUngrouped = std::move(idcs);
    mIDCsGrouped.resize(getNIntegrationIntervals());
  }

  /// \return returns number of integration intervalls stored in this object
  unsigned int getNIntegrationIntervals() const { return mIDCsUngrouped.size() / Mapper::PADSPERREGION[mIDCsGrouped.getRegion()]; }

  /// grouping and averaging of IDCs
  void processIDCs();

  /// \return returns grouped IDC object
  const auto& getIDCGroup() const { return mIDCsGrouped; }

  /// dump object to disc
  /// \param outFileName name of the output file
  /// \param outName name of the object in the output file
  void dumpToFile(const char* outFileName = "IDCAverageGroup.root", const char* outName = "IDCAverageGroup") const;

  /// draw ungrouped IDCs
  /// \param integrationInterval integration interval for which the IDCs will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawUngroupedIDCs(const unsigned int integrationInterval = 0, const std::string filename = "IDCsUngrouped.pdf") const;

  /// draw grouped IDCs
  /// \param integrationInterval integration interval for which the IDCs will be drawn
  /// \param filename name of the output file. If empty the canvas is drawn.
  void drawGroupedIDCs(const unsigned int integrationInterval = 0, const std::string filename = "IDCsGrouped.pdf") const { mIDCsGrouped.draw(integrationInterval, filename); }

  /// \return returns the stored ungrouped IDC value for local ungrouped pad row and ungrouped pad
  /// \param urow ungrouped local row in region
  /// \param upad ungrouped pad in pad direction
  /// \param integrationInterval integration interval for which the IDCs will be returned
  const float& getUngroupedIDCVal(const unsigned int urow, const unsigned int upad, const unsigned int integrationInterval) const { return mIDCsUngrouped[getUngroupedIndex(urow, upad, integrationInterval)]; }

  /// \return returns the stored ungrouped IDC value for local ungrouped pad row and ungrouped pad
  /// \param urow ungrouped local row in region
  /// \param upad ungrouped pad in pad direction
  /// \param integrationInterval integration interval for which the IDCs will be returned
  float& getUngroupedIDCVal(const unsigned int urow, const unsigned int upad, const unsigned int integrationInterval) { return mIDCsUngrouped[getUngroupedIndex(urow, upad, integrationInterval)]; }

  /// \return returns the stored grouped IDC value for local ungrouped pad row and ungrouped pad
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  float& getGroupedIDCVal(unsigned int urow, unsigned int upad, unsigned int integrationInterval) { return mIDCsGrouped.getValUngrouped(urow, upad, integrationInterval); }

  /// \return returns the stored grouped IDC value for local ungrouped pad row and ungrouped pad
  /// \param urow row of the ungrouped IDCs
  /// \param upad pad number of the ungrouped IDCs
  /// \param integrationInterval integration interval
  const float& getGroupedIDCVal(unsigned int urow, unsigned int upad, unsigned int integrationInterval) const { return mIDCsGrouped.getValUngrouped(urow, upad, integrationInterval); }

  /// get the number of threads used for some of the calculations
  static int getNThreads() { return sNThreads; }

  /// set the number of threads used for some of the calculations
  static void setNThreads(const int nThreads) { sNThreads = nThreads; }

 private:
  inline static int sNThreads{1};      ///< number of threads which are used during the calculations
  std::vector<float> mIDCsUngrouped{}; ///< integrated ungrouped IDC values per pad
  IDCGroup mIDCsGrouped{};             ///< grouped and averaged IDC values

  unsigned int getUngroupedIndex(const unsigned int urow, const unsigned int upad, const unsigned int integrationInterval) const
  {
    return integrationInterval * Mapper::PADSPERREGION[mIDCsGrouped.getRegion()] + Mapper::OFFSETCRULOCAL[mIDCsGrouped.getRegion()][urow] + upad;
  }
};

} // namespace tpc
} // namespace o2

#endif
