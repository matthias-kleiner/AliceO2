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
  for (unsigned int integrationInterval = 0; integrationInterval < getNIntegrationIntervals(); ++integrationInterval) {
    for (unsigned int irow = 0; irow < mRows; ++irow) {
      for (unsigned int ipad = 0; ipad < mPadsPerRow[irow]; ++ipad) {
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

unsigned int o2::tpc::IDCGroup::getGroupedRow(const unsigned int lrow, const unsigned int groupRows, const unsigned int groupedrows)
{
  const unsigned int row = lrow / groupRows;
  return (row >= groupedrows) ? (groupedrows - 1) : row;
}

unsigned int o2::tpc::IDCGroup::getGroupedPad(const unsigned int pad, const unsigned int lrow, const unsigned int region, const unsigned int groupPads, const unsigned int groupRows, const unsigned int groupedrows, const std::vector<unsigned int>& padsPerRow)
{
  const int relPadHalf = static_cast<int>(std::floor((pad - 0.5f * Mapper::PADSPERROW[region][lrow]) / groupPads));
  const unsigned int nGroupedPads = padsPerRow[getGroupedRow(lrow, groupRows, groupedrows)];
  const unsigned int nGroupedPadsHalf = (nGroupedPads / 2);
  if (std::abs(relPadHalf) >= nGroupedPadsHalf) {
    return std::signbit(relPadHalf) ? 0 : nGroupedPads - 1;
  }
  return static_cast<unsigned int>(static_cast<int>(nGroupedPadsHalf) + relPadHalf);
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
      const float yPos = -0.5f * (coordinate.yVals[0] + coordinate.yVals[2]); // local coordinate system is mirrored
      const float xPos = 0.5f * (coordinate.xVals[0] + coordinate.xVals[2]);
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

unsigned int o2::tpc::IDCGroup::getLastRow() const
{
  const unsigned int nTotRows = Mapper::ROWSPERREGION[mRegion];
  const unsigned int rowsReminder = nTotRows % mGroupRows;
  unsigned int lastRow = nTotRows - rowsReminder;
  if (rowsReminder <= mGroupLastRowsThreshold) {
    lastRow -= mGroupRows;
  }
  return lastRow;
}

unsigned int o2::tpc::IDCGroup::getLastPad(const unsigned int row) const
{
  const unsigned int nPads = Mapper::PADSPERROW[mRegion][row] / 2;
  const unsigned int padsReminder = nPads % mGroupPads;
  int unsigned lastPad = (padsReminder == 0) ? nPads - mGroupPads : nPads - padsReminder;
  if (padsReminder && padsReminder <= mGroupLastPadsThreshold) {
    lastPad -= mGroupPads;
  }
  return lastPad;
}

void o2::tpc::IDCGroup::initIDCGroup()
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

void o2::tpc::IDCGroup::dumpToFile(const char* outFileName, const char* outName) const
{
  TFile fOut(outFileName, "UPDATE");
  fOut.WriteObject(this, outName);
  fOut.Close();
}
