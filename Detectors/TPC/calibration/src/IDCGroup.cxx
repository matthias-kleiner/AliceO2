// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TPCCalibration/IDCGroup.h"
#include "TPCBase/Mapper.h"
#include "CommonUtils/TreeStreamRedirector.h" // for debugging
#include "TPCBase/Painter.h"
#include "TH2Poly.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"

void o2::tpc::IDCGroup::initStorage()
{
  mNIDCsPerCRU = std::accumulate(mPadsPerRow.begin(), mPadsPerRow.end(), decltype(mPadsPerRow)::value_type(0));
  mIDCsGrouped.resize(mNIDCsPerCRU);

  for (unsigned int i = 1; i < mRows; ++i) {
    const unsigned int lastInd = i - 1;
    mOffsRow[i] = mOffsRow[lastInd] + mPadsPerRow[lastInd];
  }
}

void o2::tpc::IDCGroup::dumpToTree(const char* outname) const
{
  o2::utils::TreeStreamRedirector pcstream(outname, "RECREATE");
  pcstream.GetFile()->cd();
  for (int integrationInterval = 0; integrationInterval < getNIntegrationIntervals(); ++integrationInterval) {
    for (int irow = 0; irow < mRows; ++irow) {
      for (int ipad = 0; ipad < mPadsPerRow[irow]; ++ipad) {
        float idc = (*this)(irow, ipad, integrationInterval);
        pcstream << "idcs"
                 << "row=" << irow
                 << "pad=" << ipad
                 << "IDC=" << idc
                 << "\n";
      }
    }
  }
  pcstream.Close();
}

unsigned int o2::tpc::IDCGroup::getGroupedRow(const unsigned int lrow, const unsigned int groupRows, const unsigned int rows)
{
  const unsigned int row = lrow / groupRows;
  if (row >= rows) {
    return rows - 1;
  }
  return row;
}

int o2::tpc::IDCGroup::getGroupedPad(const unsigned int pad, const unsigned int lrow, const unsigned int region, const unsigned int groupPads, const unsigned int groupRows, const unsigned int rows, const std::vector<unsigned int>& padsPerRow)
{
  const int relPadHalf = static_cast<int>(std::floor((pad - 0.5f * Mapper::PADSPERROW[region][lrow]) / groupPads));
  const int nGroupedPads = padsPerRow[getGroupedRow(lrow, groupRows, rows)];
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

void o2::tpc::IDCGroup::setRows(const unsigned int nRows)
{
  mRows = nRows;
  mPadsPerRow.resize(mRows);
  mOffsRow.resize(mRows);
}

void o2::tpc::IDCGroup::draw(const unsigned int integrationInterval, const std::string filename) const
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
  for (unsigned int irow = 0; irow < Mapper::ROWSPERREGION[mRegion]; ++irow) {
    for (unsigned int ipad = 0; ipad < Mapper::PADSPERROW[mRegion][irow]; ++ipad) {
      const auto padNum = getGlobalPadNumber(irow, ipad);
      const auto coordinate = coords[padNum];
      const float yPos = -0.5 * (coordinate.yVals[0] + coordinate.yVals[2]); // local coordinate system is mirrored
      const float xPos = 0.5 * (coordinate.xVals[0] + coordinate.xVals[2]);
      poly->Fill(xPos, yPos, (*this)(getGroupedRow(irow, mGroupRows, mRows), getGroupedPad(ipad, irow, mRegion, mGroupPads, mGroupRows, mRows, mPadsPerRow), integrationInterval));
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

int o2::tpc::IDCGroup::getLastRow() const
{
  const int nTotRows = Mapper::ROWSPERREGION[mRegion];
  const int rowsReminder = nTotRows % mGroupRows;
  int lastRow = nTotRows - rowsReminder;
  if (rowsReminder <= mGroupLastRowsThreshold) {
    lastRow -= mGroupRows;
  }
  return lastRow;
}

int o2::tpc::IDCGroup::getLastPad(const int row) const
{
  const int nPads = 0.5 * Mapper::PADSPERROW[mRegion][row];
  const int padsReminder = nPads % mGroupPads;
  int lastPad = padsReminder == 0 ? nPads - mGroupPads : nPads - padsReminder;
  if (padsReminder && padsReminder <= mGroupLastPadsThreshold) {
    lastPad -= mGroupPads;
  }
  return lastPad;
}

void o2::tpc::IDCGroup::initIDCGroup()
{
  const int lastRow = getLastRow();
  const int nRows = lastRow / mGroupRows + 1;
  setRows(nRows);
  for (int irow = 0; irow < nRows; ++irow) {
    const int row = irow * mGroupRows;
    mPadsPerRow[irow] = 2 * (getLastPad(row) / mGroupPads + 1);
  }
  initStorage();
}

void o2::tpc::IDCGroup::dumpToFile(const char* outFileName, const char* outName) const
{
  TFile fOut(outFileName, "UPDATE");
  fOut.WriteObject(this, outName);
  fOut.Close();
}
