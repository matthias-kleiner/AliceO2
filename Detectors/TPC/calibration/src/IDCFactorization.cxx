// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TPCCalibration/IDCFactorization.h"
#include "TPCCalibration/IDCGroup.h"
#include "TPCBase/Mapper.h"
#include "CommonUtils/TreeStreamRedirector.h" // for debugging
#include "TPCBase/Painter.h"
#include "TH2Poly.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "Framework/Logger.h"
#include <functional>

o2::tpc::IDCFactorization::IDCFactorization(const std::array<unsigned int, Mapper::NREGIONS>& groupPads, const std::array<unsigned int, Mapper::NREGIONS>& groupRows, const std::array<unsigned int, Mapper::NREGIONS>& groupLastRowsThreshold, const std::array<unsigned int, Mapper::NREGIONS>& groupLastPadsThreshold, const unsigned int timeFrames)
  : mGroupPads{groupPads}, mGroupRows{groupRows}, mGroupLastRowsThreshold{groupLastRowsThreshold}, mGroupLastPadsThreshold{groupLastPadsThreshold}, mTimeFrames{timeFrames}
{
  for (unsigned int cru = 0; cru < Mapper::NSECTORS * Mapper::NREGIONS; ++cru) {
    mIDCs[cru].resize(mTimeFrames);
  }

  for (unsigned int reg = 0; reg < Mapper::NREGIONS; ++reg) {
    const IDCGroup groupTmp(mGroupPads[reg], mGroupRows[reg], mGroupLastRowsThreshold[reg], mGroupLastPadsThreshold[reg], reg);
    mNIDCsPerCRU[reg] = groupTmp.getNIDCsPerCRU();
    mRows[reg] = groupTmp.getNRows();
    mPadsPerRow[reg] = groupTmp.getPadsPerRow();
    mOffsRow[reg] = groupTmp.getRowOffset();
    if (reg > 0) {
      const unsigned int lastInd = reg - 1;
      mRegionOffs[reg] = mRegionOffs[lastInd] + mNIDCsPerCRU[lastInd];
    }
  }
  mNIDCsPerSector = std::accumulate(mNIDCsPerCRU.begin(), mNIDCsPerCRU.end(), decltype(mNIDCsPerCRU)::value_type(0));
}

void o2::tpc::IDCFactorization::drawSector(const IDCType type, const unsigned int sector, const unsigned int integrationInterval, const std::string filename) const
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
  for (unsigned int region = 0; region < Mapper::NREGIONS; ++region) {
    for (unsigned int irow = 0; irow < Mapper::ROWSPERREGION[region]; ++irow) {
      for (unsigned int ipad = 0; ipad < Mapper::PADSPERROW[region][irow]; ++ipad) {
        const auto padNum = IDCGroup::getGlobalPadNumber(irow, ipad, region);
        const auto coordinate = coords[padNum];
        const float yPos = -0.5f * (coordinate.yVals[0] + coordinate.yVals[2]); // local coordinate system is mirrored
        const float xPos = 0.5f * (coordinate.xVals[0] + coordinate.xVals[2]);

        switch (type) {
          case IDCType::IDC:
            poly->Fill(xPos, yPos, getIDCVal(sector, region, irow, ipad, integrationInterval));
            break;
          case IDCType::IDCZero:
            poly->Fill(xPos, yPos, getIDCZeroVal(sector, region, irow, ipad));
            break;
          case IDCType::IDCDelta:
            poly->Fill(xPos, yPos, getIDCDeltaVal(sector, region, irow, ipad, integrationInterval));
            break;
          case IDCType::IDCOne:
            break;
        }

        lat.SetTextAlign(12);
        lat.DrawLatex(xPos, yPos, Form("%i", ipad));
      }
    }
  }
  if (!filename.empty()) {
    can->SaveAs(filename.data());
    delete poly;
    delete can;
  }
}

void o2::tpc::IDCFactorization::drawSide(const IDCType type, const o2::tpc::Side side, const unsigned int integrationInterval, const std::string filename) const
{
  const auto coords = o2::tpc::painter::getPadCoordinatesSector();
  TH2Poly* poly = o2::tpc::painter::makeSideHist(side);
  poly->SetContour(255);
  poly->SetTitle(0);
  poly->GetXaxis()->SetTitleOffset(1.2f);
  poly->GetYaxis()->SetTitleOffset(1.3f);
  poly->GetZaxis()->SetTitleOffset(1.2f);
  poly->GetZaxis()->SetTitle("#it{IDC}");
  poly->GetZaxis()->SetMaxDigits(3); // force exponential axis
  poly->SetStats(0);

  TCanvas* can = new TCanvas("can", "can", 650, 600);
  can->SetTopMargin(0.04f);
  can->SetRightMargin(0.14f);
  can->SetLeftMargin(0.1f);

  poly->Draw("colz");
  unsigned int sectorStart = (side == Side::A) ? 0 : o2::tpc::SECTORSPERSIDE;
  unsigned int sectorEnd = (side == Side::A) ? o2::tpc::SECTORSPERSIDE : Mapper::NSECTORS;
  for (unsigned int sector = sectorStart; sector < sectorEnd; ++sector) {
    for (unsigned int region = 0; region < Mapper::NREGIONS; ++region) {
      for (unsigned int irow = 0; irow < Mapper::ROWSPERREGION[region]; ++irow) {
        for (unsigned int ipad = 0; ipad < Mapper::PADSPERROW[region][irow]; ++ipad) {
          const auto padNum = IDCGroup::getGlobalPadNumber(irow, ipad, region);
          const float angDeg = 10.f + sector * 20;
          auto coordinate = coords[padNum];
          coordinate.rotate(angDeg);
          const float yPos = 0.25f * (coordinate.yVals[0] + coordinate.yVals[1] + coordinate.yVals[2] + coordinate.yVals[3]);
          const float xPos = 0.25f * (coordinate.xVals[0] + coordinate.xVals[1] + coordinate.xVals[2] + coordinate.xVals[3]);
          const auto padTmp = (side == Side::A) ? ipad : (Mapper::PADSPERROW[region][irow] - ipad); // C-Side is mirrored
          switch (type) {
            case IDCType::IDC:
              poly->Fill(xPos, yPos, getIDCVal(sector, region, irow, padTmp, integrationInterval));
              break;
            case IDCType::IDCZero:
              poly->Fill(xPos, yPos, getIDCZeroVal(sector, region, irow, padTmp));
              break;
            case IDCType::IDCDelta:
              poly->Fill(xPos, yPos, getIDCDeltaVal(sector, region, irow, padTmp, integrationInterval));
              break;
            case IDCType::IDCOne:
              break;
          }
        }
      }
    }
  }
  if (!filename.empty()) {
    can->SaveAs(filename.data());
    delete poly;
    delete can;
  }
}

void o2::tpc::IDCFactorization::dumpToFile(const char* outFileName, const char* outName) const
{
  TFile fOut(outFileName, "RECREATE");
  fOut.WriteObject(this, outName);
  fOut.Close();
}

void o2::tpc::IDCFactorization::dumpIDCsToTree(int integrationIntervals) const
{
  const static Mapper& mapper = Mapper::instance();
  o2::utils::TreeStreamRedirector pcstream("IDCTree.root", "RECREATE");
  pcstream.GetFile()->cd();

  if (integrationIntervals <= 0) {
    integrationIntervals = getNIntegrationIntervals();
  }

  std::vector<float> idcOneA = mIDCZeroOne.mIDCOne[0];
  std::vector<float> idcOneC = mIDCZeroOne.mIDCOne[1];
  for (unsigned int integrationInterval = 0; integrationInterval < integrationIntervals; ++integrationInterval) {
    const unsigned int nIDCsSector = Mapper::getPadsInSector() * Mapper::NSECTORS;
    std::vector<int> vRow(nIDCsSector);
    std::vector<int> vPad(nIDCsSector);
    std::vector<float> vXPos(nIDCsSector);
    std::vector<float> vYPos(nIDCsSector);
    std::vector<float> vGlobalXPos(nIDCsSector);
    std::vector<float> vGlobalYPos(nIDCsSector);
    std::vector<float> idcs(nIDCsSector);
    std::vector<float> idcsZero(nIDCsSector);
    std::vector<float> idcsDelta(nIDCsSector);
    std::vector<unsigned int> sectorv(nIDCsSector);

    unsigned int index = 0;
    for (unsigned int sector = 0; sector < Mapper::NSECTORS; ++sector) {
      for (unsigned int region = 0; region < Mapper::NREGIONS; ++region) {
        for (unsigned int irow = 0; irow < Mapper::ROWSPERREGION[region]; ++irow) {
          for (unsigned int ipad = 0; ipad < Mapper::PADSPERROW[region][irow]; ++ipad) {
            const auto padNum = IDCGroup::getGlobalPadNumber(irow, ipad, region);
            const auto padTmp = (sector < SECTORSPERSIDE) ? ipad : (Mapper::PADSPERROW[region][irow] - ipad); // C-Side is mirrored
            const auto& padPosLocal = mapper.padPos(padNum);
            vRow[index] = padPosLocal.getRow();
            vPad[index] = padPosLocal.getPad();
            vXPos[index] = mapper.getPadCentre(padPosLocal).X();
            vYPos[index] = mapper.getPadCentre(padPosLocal).Y();
            const GlobalPosition2D globalPos = mapper.LocalToGlobal(LocalPosition2D(vXPos[index], vYPos[index]), Sector(sector));
            vGlobalXPos[index] = globalPos.X();
            vGlobalYPos[index] = globalPos.Y();
            idcs[index] = getIDCVal(sector, region, irow, padTmp, integrationInterval);
            idcsZero[index] = getIDCZeroVal(sector, region, irow, padTmp);
            idcsDelta[index] = getIDCDeltaVal(sector, region, irow, padTmp, integrationInterval);
            sectorv[index] = sector;
            ++index;
          }
        }
      }
    }
    pcstream << "tree"
             << "integrationInterval=" << integrationInterval
             << "IDC.=" << idcs
             << "IDC0.=" << idcsZero
             << "IDC1A.=" << idcOneA
             << "IDC1C.=" << idcOneC
             << "IDCDelta.=" << idcsDelta
             << "pad.=" << vPad
             << "row.=" << vRow
             << "lx.=" << vXPos
             << "ly.=" << vYPos
             << "gx.=" << vGlobalXPos
             << "gy.=" << vGlobalYPos
             << "sector.=" << sectorv
             << "\n";
  }
  pcstream.Close();
}

void o2::tpc::IDCFactorization::calcIDCZero()
{
  const unsigned int nIDCsSide = mNIDCsPerSector * o2::tpc::SECTORSPERSIDE;
  mIDCZeroOne.mIDCZero[Side::A].resize(nIDCsSide);
  mIDCZeroOne.mIDCZero[Side::C].resize(nIDCsSide);

#pragma omp parallel for num_threads(sNThreads)
  for (unsigned int cru = 0; cru < mIDCs.size(); ++cru) {
    const o2::tpc::CRU cruTmp(cru);
    const unsigned int region = cruTmp.region();
    for (unsigned int timeframe = 0; timeframe < mTimeFrames; ++timeframe) {
      for (unsigned int idcs = 0; idcs < mIDCs[cru][timeframe].size(); ++idcs) {
        const unsigned int indexGlob = (idcs % mNIDCsPerCRU[region]) + mRegionOffs[region] + mNIDCsPerSector * cruTmp.sector();
        mIDCZeroOne.getValueIDCZero(cruTmp.side(), indexGlob % nIDCsSide) += mIDCs[cru][timeframe][idcs];
      }
    }
  }
  std::transform(mIDCZeroOne.mIDCZero[Side::A].begin(), mIDCZeroOne.mIDCZero[Side::A].end(), mIDCZeroOne.mIDCZero[Side::A].begin(), std::bind(std::divides<float>(), std::placeholders::_1, getNIntegrationIntervals()));
  std::transform(mIDCZeroOne.mIDCZero[Side::C].begin(), mIDCZeroOne.mIDCZero[Side::C].end(), mIDCZeroOne.mIDCZero[Side::C].begin(), std::bind(std::divides<float>(), std::placeholders::_1, getNIntegrationIntervals()));
}

void o2::tpc::IDCFactorization::calcIDCOne()
{
  const unsigned int nIDCsSide = mNIDCsPerSector * SECTORSPERSIDE;
  const unsigned int integrationIntervals = getNIntegrationIntervals();
  mIDCZeroOne.mIDCOne[Side::A].resize(integrationIntervals);
  mIDCZeroOne.mIDCOne[Side::C].resize(integrationIntervals);

#pragma omp parallel for num_threads(sNThreads)
  for (unsigned int cru = 0; cru < mIDCs.size(); ++cru) {
    const o2::tpc::CRU cruTmp(cru);
    const unsigned int region = cruTmp.region();
    for (unsigned int timeframe = 0; timeframe < mTimeFrames; ++timeframe) {
      for (unsigned int idcs = 0; idcs < mIDCs[cru][timeframe].size(); ++idcs) {
        const unsigned int integrationInterval = idcs / mNIDCsPerCRU[region];
        const auto side = cruTmp.side();
        const unsigned int indexGlob = (idcs % mNIDCsPerCRU[region]) + mRegionOffs[region] + mNIDCsPerSector * cruTmp.sector();
        mIDCZeroOne.mIDCOne[side][integrationInterval] += mIDCs[cru][timeframe][idcs] / mIDCZeroOne.mIDCZero[side][indexGlob % nIDCsSide];
      }
    }
  }

  std::transform(mIDCZeroOne.mIDCOne[Side::A].begin(), mIDCZeroOne.mIDCOne[Side::A].end(), mIDCZeroOne.mIDCOne[Side::A].begin(), std::bind(std::divides<float>(), std::placeholders::_1, nIDCsSide));
  std::transform(mIDCZeroOne.mIDCOne[Side::C].begin(), mIDCZeroOne.mIDCOne[Side::C].end(), mIDCZeroOne.mIDCOne[Side::C].begin(), std::bind(std::divides<float>(), std::placeholders::_1, nIDCsSide));
}

void o2::tpc::IDCFactorization::calcIDCDelta()
{
  const unsigned int nIDCsSide = mNIDCsPerSector * SECTORSPERSIDE;
  const unsigned int idcsSide = nIDCsSide * getNIntegrationIntervals();
  mIDCDelta.mIDCDelta[Side::A].resize(idcsSide);
  mIDCDelta.mIDCDelta[Side::C].resize(idcsSide);

#pragma omp parallel for num_threads(sNThreads)
  for (unsigned int cru = 0; cru < mIDCs.size(); ++cru) {
    const o2::tpc::CRU cruTmp(cru);
    const unsigned int region = cruTmp.region();
    for (unsigned int timeframe = 0; timeframe < mTimeFrames; ++timeframe) {
      for (unsigned int idcs = 0; idcs < mIDCs[cru][timeframe].size(); ++idcs) {
        const unsigned int integrationInterval = idcs / mNIDCsPerCRU[region];
        const auto side = cruTmp.side();
        const unsigned int indexGlob = (idcs % mNIDCsPerCRU[region]) + mRegionOffs[region] + mNIDCsPerSector * cruTmp.sector();
        const unsigned int indexGlobMod = indexGlob % nIDCsSide;
        mIDCDelta.mIDCDelta[side][indexGlobMod + integrationInterval * nIDCsSide] += mIDCs[cru][timeframe][idcs] / (mIDCZeroOne.mIDCZero[side][indexGlobMod] * mIDCZeroOne.mIDCOne[side][integrationInterval]);
      }
    }
  }
}
