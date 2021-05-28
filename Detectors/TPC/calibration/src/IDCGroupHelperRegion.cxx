// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TPCCalibration/IDCGroupHelperRegion.h"
#include "TPCBase/Mapper.h"
#include "TFile.h"
#include <numeric>

void o2::tpc::IDCGroupHelperRegion::initStorage()
{
  mNIDCsPerCRU = std::accumulate(mPadsPerRow.begin(), mPadsPerRow.end(), decltype(mPadsPerRow)::value_type(0));
  for (unsigned int i = 1; i < mRows; ++i) {
    const unsigned int lastInd = i - 1;
    mOffsRow[i] = mOffsRow[lastInd] + mPadsPerRow[lastInd];
  }
}

unsigned int o2::tpc::IDCGroupHelperRegion::getGroupedRow(const unsigned int lrow, const unsigned int groupRows, const unsigned int groupedrows)
{
  const unsigned int row = lrow / groupRows;
  return (row >= groupedrows) ? (groupedrows - 1) : row;
}

unsigned int o2::tpc::IDCGroupHelperRegion::getGroupedPad(const unsigned int pad, const unsigned int lrow, const unsigned int region, const unsigned int groupPads, const unsigned int groupRows, const unsigned int groupedrows, const std::vector<unsigned int>& padsPerRow)
{
  const int relPadHalf = static_cast<int>(std::floor((pad - 0.5f * Mapper::PADSPERROW[region][lrow]) / groupPads));
  const unsigned int nGroupedPads = padsPerRow[getGroupedRow(lrow, groupRows, groupedrows)];
  const unsigned int nGroupedPadsHalf = (nGroupedPads / 2);
  if (std::abs(relPadHalf) >= nGroupedPadsHalf) {
    return std::signbit(relPadHalf) ? 0 : nGroupedPads - 1;
  }
  return static_cast<unsigned int>(static_cast<int>(nGroupedPadsHalf) + relPadHalf);
}

void o2::tpc::IDCGroupHelperRegion::setRows(const unsigned int nRows)
{
  mRows = nRows;
  mPadsPerRow.resize(mRows);
  mOffsRow.resize(mRows);
}

unsigned int o2::tpc::IDCGroupHelperRegion::getLastRow() const
{
  const unsigned int nTotRows = Mapper::ROWSPERREGION[mRegion];
  const unsigned int rowsReminder = nTotRows % mGroupRows;
  unsigned int lastRow = nTotRows - rowsReminder;
  if (rowsReminder <= mGroupLastRowsThreshold) {
    lastRow -= mGroupRows;
  }
  return lastRow;
}

unsigned int o2::tpc::IDCGroupHelperRegion::getLastPad(const unsigned int row) const
{
  const unsigned int nPads = Mapper::PADSPERROW[mRegion][row] / 2;
  const unsigned int padsReminder = nPads % mGroupPads;
  int unsigned lastPad = (padsReminder == 0) ? nPads - mGroupPads : nPads - padsReminder;
  if (padsReminder && padsReminder <= mGroupLastPadsThreshold) {
    lastPad -= mGroupPads;
  }
  return lastPad;
}

void o2::tpc::IDCGroupHelperRegion::initIDCGroupHelperRegion()
{
  const unsigned int lastRow = getLastRow();
  const unsigned int nRows = lastRow / mGroupRows + 1;
  setRows(nRows);
  for (unsigned int irow = 0; irow < nRows; ++irow) {
    const unsigned int row = irow * mGroupRows;
    mPadsPerRow[irow] = 2 * (getLastPad(row) / mGroupPads + 1);
  }
  initStorage();
}

void o2::tpc::IDCGroupHelperRegion::dumpToFile(const char* outFileName, const char* outName) const
{
  TFile fOut(outFileName, "UPDATE");
  fOut.WriteObject(this, outName);
  fOut.Close();
}
