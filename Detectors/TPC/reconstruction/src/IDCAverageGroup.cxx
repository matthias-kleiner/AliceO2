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

#include "TH2Poly.h"
#include "TLatex.h"
#include "TPCBase/Painter.h"
#include "TCanvas.h"
#include "TFile.h"
#include "Framework/Logger.h"

void o2::tpc::IDCAverageGroup::processIDCs()
{
  for (int integrationInterval = 0; integrationInterval < getNIntegrationIntervals(); ++integrationInterval) {
    const int lastRow = getLastRow();
    int rowGrouped = 0;
    for (int iRow = 0; iRow <= lastRow; iRow += mGroupRows) {
      for (int iYLocalSide = 0; iYLocalSide < 2; ++iYLocalSide) {
        const int nPads = 0.5 * mPadsPerRow[mRegion][iRow];
        const int lastPad = getLastPad(iRow);
        const int endPads = lastPad + nPads;

        const int halfPadsInRow = 0.5 * mIDCsGrouped.getPadsPerRow(rowGrouped);
        int padGrouped = !iYLocalSide ? halfPadsInRow - 1 : halfPadsInRow;
        for (int ipad = nPads; ipad <= endPads; ipad += mGroupPads) {
          int endRows = iRow == lastRow ? (mRowsPerRegion[mRegion] - iRow) : mGroupRows;

          float idcAverage = 0;
          float idcCounter = 0;
          for (int iRowMerge = 0; iRowMerge < endRows; ++iRowMerge) {
            const int iRowTmp = iRow + iRowMerge;
            const int offs = mAddPadsPerRow[mRegion][iRowTmp] - mAddPadsPerRow[mRegion][iRow];
            const int padStart = ipad == 0 ? 0 : offs;
            const int endPadsTmp = ipad == endPads ? (mPadsPerRow[mRegion][iRowTmp] - ipad) : mGroupPads + offs;
            for (int ipadMerge = padStart; ipadMerge < endPadsTmp; ++ipadMerge) {
              const int iPadTmp = ipad + ipadMerge;
              const int nPadsInRow = mPadsPerRow[mRegion][iRowTmp];
              const int iPadSide = iYLocalSide == 0 ? nPadsInRow - iPadTmp - 1 : iPadTmp;
              const int indexIDC = integrationInterval * mPadsPerRegion[mRegion] + mOffs[mRegion][iRowTmp] + iPadSide;

              idcAverage += mIDCs[indexIDC] * mPadArea[mRegion];
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

unsigned int o2::tpc::IDCAverageGroup::getGroupedRow(const unsigned int lrow) const
{
  const unsigned int row = lrow / mGroupRows;
  if (row >= mIDCsGrouped.getNRows()) {
    return mIDCsGrouped.getNRows() - 1;
  }
  return row;
}

int o2::tpc::IDCAverageGroup::getGroupedPad(const unsigned int pad, const unsigned int lrow) const
{
  const int relPadHalf = static_cast<int>(std::floor((pad - 0.5f * mPadsPerRow[mRegion][lrow]) / mGroupPads));
  const int nGroupedPads = mIDCsGrouped.getPadsPerRow(getGroupedRow(lrow));
  const int nGroupedPadsHalf = static_cast<int>(0.5f * nGroupedPads);
  if (std::abs(relPadHalf) >= nGroupedPadsHalf) {
    if (std::signbit(relPadHalf)) {
      return 0; // if pad is negative
    } else {
      return nGroupedPads - 1;
    }
  }
  return nGroupedPadsHalf + relPadHalf;
}

void o2::tpc::IDCAverageGroup::draw(const int integrationInterval) const
{
  const auto coords = o2::tpc::painter::getPadCoordinatesSector();
  TH2Poly* poly = o2::tpc::painter::makeSectorHist("hSector", "Sector;local #it{x} (cm);local #it{y} (cm); #it{IDC}");
  poly->SetContour(255);
  poly->SetTitle(0);
  poly->GetYaxis()->SetTickSize(0.002f);
  poly->GetYaxis()->SetTitleOffset(0.7f);
  poly->GetZaxis()->SetTitleOffset(1.3f);
  poly->SetStats(0);

  TCanvas* can = new TCanvas("can", "can", 2000, 1400);
  can->SetRightMargin(0.14f);
  can->SetLeftMargin(0.06f);
  can->SetTopMargin(0.04f);

  TLatex lat;
  lat.SetTextFont(63);
  lat.SetTextSize(2);

  poly->Draw("colz");
  for (unsigned int irow = 0; irow < mRowsPerRegion[mRegion]; ++irow) {
    for (unsigned int ipad = 0; ipad < mPadsPerRow[mRegion][irow]; ++ipad) {
      const GlobalPadNumber padNum = getGlobalPadNumber(irow, ipad);
      const auto coordinate = coords[padNum];
      const float yPos = -0.5 * (coordinate.yVals[0] + coordinate.yVals[2]); // local coordinate system is mirrored
      const float xPos = 0.5 * (coordinate.xVals[0] + coordinate.xVals[2]);
      poly->Fill(xPos, yPos, mIDCsGrouped(getGroupedRow(irow), getGroupedPad(ipad, irow), integrationInterval));
      lat.SetTextAlign(12);
      lat.DrawLatex(xPos, yPos, Form("%i", ipad));
    }
  }
}

int o2::tpc::IDCAverageGroup::getLastRow() const
{
  const int nTotRows = mRowsPerRegion[mRegion];
  const int rowsReminder = nTotRows % mGroupRows;
  int lastRow = nTotRows - rowsReminder;
  if (rowsReminder <= mGroupLastRowsThreshold) {
    lastRow -= mGroupRows;
  }
  return lastRow;
}

int o2::tpc::IDCAverageGroup::getLastPad(const int row) const
{
  const int nPads = 0.5 * mPadsPerRow[mRegion][row];
  const int padsReminder = nPads % mGroupPads;
  int lastPad = padsReminder == 0 ? nPads - mGroupPads : nPads - padsReminder;
  if (padsReminder && padsReminder <= mGroupLastPadsThreshold) {
    lastPad -= mGroupPads;
  }
  return lastPad;
}

void o2::tpc::IDCAverageGroup::initIDCGroup()
{
  const int lastRow = getLastRow();
  const int nRows = lastRow / mGroupRows + 1;
  mIDCsGrouped.setRows(nRows);
  for (int irow = 0; irow < nRows; ++irow) {
    const int row = irow * mGroupRows;
    const int nPadsInRow = 2 * (getLastPad(row) / mGroupPads + 1);
    mIDCsGrouped.setPadsPerRow(irow, nPadsInRow);
  }
  mIDCsGrouped.initStorage();
}

void o2::tpc::IDCAverageGroup::dumpToFile(const char* outFileName, const char* outName) const
{
  TFile fOut(outFileName, "UPDATE");
  fOut.WriteObject(this, outName);
  fOut.Close();
}
