// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TPCCalibration/IDCAverageGroup.h"
#include "TPCCalibration/IDCGroup.h"
#include "TPCCalibration/RobustAverage.h"

#include "TFile.h"
#include "TPCBase/Painter.h"
#include "TH2Poly.h"
#include "TCanvas.h"
#include "TLatex.h"

void o2::tpc::IDCAverageGroup::processIDCs()
{
#pragma omp parallel for num_threads(sNThreads)
  for (int integrationInterval = 0; integrationInterval < getNIntegrationIntervals(); ++integrationInterval) {
    const int lastRow = mIDCsGrouped.getLastRow();
    int rowGrouped = 0;
    for (int iRow = 0; iRow <= lastRow; iRow += mIDCsGrouped.getGroupRows()) {
      // the sectors is divide in to two parts around ylocal=0 to get the same simmetric grouping around ylocal=0
      for (int iYLocalSide = 0; iYLocalSide < 2; ++iYLocalSide) {
        const unsigned int region = mIDCsGrouped.getRegion();
        const int nPads = 0.5f * Mapper::PADSPERROW[region][iRow];
        const int lastPad = mIDCsGrouped.getLastPad(iRow);
        const int endPads = lastPad + nPads;

        const int halfPadsInRow = 0.5 * mIDCsGrouped.getPadsPerRow(rowGrouped);
        int padGrouped = !iYLocalSide ? halfPadsInRow - 1 : halfPadsInRow;
        for (int ipad = nPads; ipad <= endPads; ipad += mIDCsGrouped.getGroupPads()) {
          int endRows = (iRow == lastRow) ? (Mapper::ROWSPERREGION[region] - iRow) : mIDCsGrouped.getGroupRows();

          // TODO Mapper::ADDITIONALPADSPERROW[region].back() factor is to large, but doesnt really matter
          const unsigned int maxGoup = (mIDCsGrouped.getGroupRows() + mIDCsGrouped.getGroupLastRowsThreshold()) * (mIDCsGrouped.getGroupPads() + mIDCsGrouped.getGroupLastPadsThreshold() + Mapper::ADDITIONALPADSPERROW[region].back());
          RobustAverage robustAverage(maxGoup);

          for (int iRowMerge = 0; iRowMerge < endRows; ++iRowMerge) {
            const int iRowTmp = iRow + iRowMerge;
            const auto offs = Mapper::ADDITIONALPADSPERROW[region][iRowTmp] - Mapper::ADDITIONALPADSPERROW[region][iRow];
            const auto padStart = (ipad == 0) ? 0 : offs;
            const int endPadsTmp = (ipad == endPads) ? (Mapper::PADSPERROW[region][iRowTmp] - ipad) : mIDCsGrouped.getGroupPads() + offs;
            for (int ipadMerge = padStart; ipadMerge < endPadsTmp; ++ipadMerge) {
              const int iPadTmp = ipad + ipadMerge;
              const int iPadSide = iYLocalSide == 0 ? Mapper::PADSPERROW[region][iRowTmp] - iPadTmp - 1 : iPadTmp;
              const int indexIDC = integrationInterval * Mapper::PADSPERREGION[region] + Mapper::OFFSETCRULOCAL[region][iRowTmp] + iPadSide;

              robustAverage.addValue(mIDCsUngrouped[indexIDC] * Mapper::PADAREA[region]);
            }
          }
          mIDCsGrouped(rowGrouped, padGrouped, integrationInterval) = robustAverage.getFilteredAverage();

          iYLocalSide ? ++padGrouped : --padGrouped;
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

void o2::tpc::IDCAverageGroup::drawUngroupedIDCs(const unsigned int integrationInterval, const std::string filename) const
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
  const unsigned int region = mIDCsGrouped.getRegion();
  for (unsigned int irow = 0; irow < Mapper::ROWSPERREGION[region]; ++irow) {
    for (unsigned int ipad = 0; ipad < Mapper::PADSPERROW[region][irow]; ++ipad) {
      const auto padNum = IDCGroup::getGlobalPadNumber(irow, ipad, region);
      const auto coordinate = coords[padNum];
      const float yPos = -0.5 * (coordinate.yVals[0] + coordinate.yVals[2]); // local coordinate system is mirrored
      const float xPos = 0.5 * (coordinate.xVals[0] + coordinate.xVals[2]);

      const unsigned int indexIDC = integrationInterval * Mapper::PADSPERREGION[region] + Mapper::OFFSETCRULOCAL[region][irow] + ipad;
      const float idc = mIDCsUngrouped[indexIDC] * Mapper::PADAREA[region];

      poly->Fill(xPos, yPos, idc);

      lat.SetTextAlign(12);
      lat.DrawLatex(xPos, yPos, Form("%i", ipad));
    }
  }

  if (!filename.empty()) {
    can->SaveAs(filename.data());
    delete poly;
    delete can;
  }
}
