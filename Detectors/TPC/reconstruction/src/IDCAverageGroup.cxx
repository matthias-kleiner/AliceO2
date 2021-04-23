// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TPCReconstruction/IDCAverageGroup.h"
#include "TFile.h"

void o2::tpc::IDCAverageGroup::processIDCs()
{
  for (int integrationInterval = 0; integrationInterval < getNIntegrationIntervals(); ++integrationInterval) {
    const int lastRow = mIDCsGrouped.getLastRow();
    int rowGrouped = 0;
    for (int iRow = 0; iRow <= lastRow; iRow += mIDCsGrouped.getGroupRows()) {
      for (int iYLocalSide = 0; iYLocalSide < 2; ++iYLocalSide) {
        const int region = mIDCsGrouped.getRegion();
        const int nPads = 0.5 * IDCHelper::PADSPERROW[region][iRow];
        const int lastPad = mIDCsGrouped.getLastPad(iRow);
        const int endPads = lastPad + nPads;

        const int halfPadsInRow = 0.5 * mIDCsGrouped.getPadsPerRow(rowGrouped);
        int padGrouped = !iYLocalSide ? halfPadsInRow - 1 : halfPadsInRow;
        for (int ipad = nPads; ipad <= endPads; ipad += mIDCsGrouped.getGroupPads()) {
          int endRows = iRow == lastRow ? (IDCHelper::ROWSPERREGION[region] - iRow) : mIDCsGrouped.getGroupRows();

          float idcAverage = 0;
          float idcCounter = 0;
          for (int iRowMerge = 0; iRowMerge < endRows; ++iRowMerge) {
            const int iRowTmp = iRow + iRowMerge;
            const int offs = IDCHelper::ADDITIONALPADSPERROW[region][iRowTmp] - IDCHelper::ADDITIONALPADSPERROW[region][iRow];
            const int padStart = ipad == 0 ? 0 : offs;
            const int endPadsTmp = ipad == endPads ? (IDCHelper::PADSPERROW[region][iRowTmp] - ipad) : mIDCsGrouped.getGroupPads() + offs;
            for (int ipadMerge = padStart; ipadMerge < endPadsTmp; ++ipadMerge) {
              const int iPadTmp = ipad + ipadMerge;
              const int nPadsInRow = IDCHelper::PADSPERROW[region][iRowTmp];
              const int iPadSide = iYLocalSide == 0 ? nPadsInRow - iPadTmp - 1 : iPadTmp;
              const int indexIDC = integrationInterval * IDCHelper::PADSPERREGION[region] + IDCHelper::OFFSETCRULOCAL[region][iRowTmp] + iPadSide;

              idcAverage += mIDCs[indexIDC] * IDCHelper::PADAREA[region];
              ++idcCounter;
            }
          }
          mIDCsGrouped(rowGrouped, padGrouped, integrationInterval) = idcAverage / idcCounter;
          !iYLocalSide ? --padGrouped : ++padGrouped;
        }
      }
      ++rowGrouped;
    }
  }
}

void o2::tpc::IDCAverageGroup::dumpToFile(const char* outFileName, const char* outName) const
{
  TFile fOut(outFileName, "UPDATE");
  fOut.WriteObject(this, outName);
  fOut.Close();
}
