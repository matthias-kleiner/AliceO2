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

void o2::tpc::IDCAverageGroup::processIDCs(const bool debug, const int lane)
{
  for (int integrationInterval = 0; integrationInterval < mIntegrationIntervals; ++integrationInterval) {
    const auto coords = o2::tpc::painter::getPadCoordinatesSector();
    TH2Poly* poly = o2::tpc::painter::makeSectorHist("hSector", "Sector;local #it{x} (cm);local #it{y} (cm)");
    poly->SetContour(255);

    TCanvas can("can", "can", 2000, 1400);
    can.SetRightMargin(0.14);
    can.SetLeftMargin(0.06);
    can.SetTopMargin(0.04);
    poly->Draw("colz");
    poly->SetTitle(0);
    poly->GetYaxis()->SetTickSize(0.002);
    poly->GetYaxis()->SetTitleOffset(0.7);
    poly->SetStats(0);
    TLatex lat;
    lat.SetTextFont(63);
    lat.SetTextSize(2);

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

              const float idc = mIDCs[indexIDC];
              idcAverage += idc;
              ++idcCounter;
            }
          }

          mIDCsGrouped(rowGrouped, padGrouped) = idcAverage / idcCounter;
          !iYLocalSide ? --padGrouped : ++padGrouped;

          if (debug) {
            for (int iRowMerge = 0; iRowMerge < endRows; ++iRowMerge) {
              const int iRowTmp = iRow + iRowMerge;
              const int offs = mAddPadsPerRow[mRegion][iRowTmp] - mAddPadsPerRow[mRegion][iRow];
              const int padStart = ipad == 0 ? 0 : offs;
              const int endPadsTmp = ipad == endPads ? (mPadsPerRow[mRegion][iRowTmp] - ipad) : mGroupPads + offs;
              for (int ipadMerge = padStart; ipadMerge < endPadsTmp; ++ipadMerge) {
                const int iPadTmp = ipad + ipadMerge;
                const int nPadsInRow = mPadsPerRow[mRegion][iRowTmp];
                const int iPadSide = iYLocalSide == 0 ? nPadsInRow - iPadTmp - 1 : iPadTmp;

                const GlobalPadNumber padNum = getGlobalPadIndex(mRegion, iRowTmp, iPadSide);
                const auto coordinate = coords[padNum];
                const float yPos = -0.5 * (coordinate.yVals[0] + coordinate.yVals[2]); // local coordinate system is mirrored
                const float xPos = 0.5 * (coordinate.xVals[0] + coordinate.xVals[2]);

                poly->Fill(xPos, yPos, idcAverage / idcCounter);

                lat.SetTextAlign(12);
                lat.DrawLatex(xPos, yPos, Form("%i", iPadSide));
              }
            }
          }
        }
      }
      ++rowGrouped;
    }

    if (debug) {
      // TFile fOut( Form("grouped_IDCs_rows-%i_pads-%i_mRow-%i_mPad-%i.root", mGroupRows, mGroupPads, mGroupLastRowsThreshold, mGroupLastPadsThreshold), "UPDATE");
      TFile* fOut = TFile::Open(Form("grouped_IDCs_%i_rows-%i_pads-%i_mRow-%i_mPad-%i.root", lane, mGroupRows, mGroupPads, mGroupLastRowsThreshold, mGroupLastPadsThreshold), "UPDATE");
      can.Write(Form("CRU_can_%i_%i", mCRU, integrationInterval));
      poly->Write(Form("CRU_%i_%i", mCRU, integrationInterval));

      fOut->WriteObject(&mIDCsGrouped, Form("CRU_obj_%i_%i", mCRU, integrationInterval));
      fOut->Close();
      delete fOut;
    }
    delete poly;
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
