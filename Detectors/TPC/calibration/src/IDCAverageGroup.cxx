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

#include "TPCCalibration/IDCAverageGroup.h"
#include "TPCCalibration/IDCGroup.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "TPCCalibration/IDCGroupingParameter.h"
#include "TPCBase/Mapper.h"
#include "CommonConstants/MathConstants.h"
#include "TPCCalibration/IDCDrawHelper.h"

#include "TFile.h"
#include "TKey.h"
#include "TPCBase/Painter.h"
#include "TH2Poly.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TKey.h"
#include "Framework/Logger.h"
#include "TStyle.h"

#if (defined(WITH_OPENMP) || defined(_OPENMP)) && !defined(__CLING__)
#include <omp.h>
#else
static inline int omp_get_thread_num() { return 0; }
#endif

o2::tpc::IDCAverageGroup::IDCAverageGroup(const unsigned char groupPads, const unsigned char groupRows, const unsigned char groupLastRowsThreshold, const unsigned char groupLastPadsThreshold, const unsigned int region, const Sector sector, const float sigma, const unsigned char overlapRows, const unsigned char overlapPads)
  : mSector{sector}, mSigma{sigma}, mGroup(groupPads, groupRows, groupLastRowsThreshold, groupLastPadsThreshold, region, sNThreads), mOverlapRows{overlapRows}, mOverlapPads{overlapPads}
{
  unsigned int maxValues = 0;
  for (unsigned int i = 0; i < Mapper::NREGIONS; ++i) {
    const unsigned int maxGroup = (mGroup.mIDCsGrouped.getGroupRows() + mGroup.mIDCsGrouped.getGroupLastRowsThreshold()) * (mGroup.mIDCsGrouped.getGroupPads() + mGroup.mIDCsGrouped.getGroupLastPadsThreshold() + Mapper::ADDITIONALPADSPERROW[i].back());
    if (maxGroup > maxValues) {
      maxValues = maxGroup;
    }
  }

  for (auto& rob : mGroup.mRobustAverage) {
    rob.reserve(maxValues);
  }

  // init weights
  const float sigmaEdge = 1.f;
  mWeightsPad.reserve(mOverlapPads);
  for (int i = 0; i < mOverlapPads; ++i) {
    const float groupPadsHalf = groupPads / 2.f;
    const float sigmaPad = groupPadsHalf / sigmaEdge; // assume 3-sigma at the edge of the last pad
    mWeightsPad.emplace_back(normal_dist(groupPadsHalf + i, sigmaPad));
  }

  mWeightsRow.reserve(mOverlapRows);
  for (int i = 0; i < mOverlapRows; ++i) {
    const float groupRowsHalf = groupRows / 2.f;
    const float sigmaRow = groupRowsHalf / sigmaEdge; // assume 3-sigma at the edge of the last pad
    mWeightsRow.emplace_back(normal_dist(groupRowsHalf + i, sigmaRow));
  }
}

void o2::tpc::IDCAverageGroup::updatePadStatusMap()
{
  /// TODO USE CCDB instead of TFILE!!!
  TFile file("padflags.root");
  CalDet<PadFlags>* padStatus{nullptr};
  file.GetObject("flags", padStatus);
  if (!padStatus) {
    LOG(FATAL) << "No valid pedestal object was loaded";
  }

  const auto type = padStatus->getPadSubset();
  if (type != PadSubset::Region) {
    LOG(FATAL) << "Wrong pad subset type";
  }

  mPadStatus.reset(padStatus);
}

float o2::tpc::IDCAverageGroup::normal_dist(const float x, const float sigma)
{
  const float fac = x / sigma;
  return std::exp(-0.5f * fac * fac);
}

void o2::tpc::IDCAverageGroup::processIDCs()
{
#pragma omp parallel for num_threads(sNThreads)
  for (unsigned int integrationInterval = 0; integrationInterval < getNIntegrationIntervals(); ++integrationInterval) {
    const unsigned int threadNum = omp_get_thread_num();
    loopOverGroups(threadNum, integrationInterval, mGroup);
  }
}

void o2::tpc::IDCAverageGroup::drawGrouping()
{
  const auto& mapper = Mapper::instance();
  TH2Poly* poly = o2::tpc::painter::makeSectorHist("hSector", "Sector;#it{x} (cm);#it{y} (cm)");
  poly->SetContour(255);
  gStyle->SetNumberContours(255);

  TCanvas can("can", "can", 2000, 1400);
  can.SetRightMargin(0.01f);
  can.SetLeftMargin(0.06f);
  can.SetTopMargin(0.04f);
  poly->SetTitle(0);
  poly->GetYaxis()->SetTickSize(0.002f);
  poly->GetYaxis()->SetTitleOffset(0.7f);
  poly->SetStats(0);

  for (unsigned int i = 0; i < Mapper::NREGIONS; ++i) {
    IDCAverageGroupType<Draw> type(mGroup.mIDCsGrouped.getGroupPads(), mGroup.mIDCsGrouped.getGroupRows(), mGroup.mIDCsGrouped.getGroupLastRowsThreshold(), mGroup.mIDCsGrouped.getGroupLastPadsThreshold(), i, Mapper::PADSPERREGION[i], mapper.getPadRegionInfo(i), *poly);
    loopOverGroups(0, 0, type);
  }

  poly->Draw("col");
  painter::drawSectorLocalPadNumberPoly(kBlack);
  painter::drawSectorInformationPoly(kRed, kRed);
  can.SaveAs(Form("grouping_rows-%i_pads-%i_rowThr-%i_padThr-%i_ovRows-%i_ovPads-%i.pdf", mGroup.mIDCsGrouped.getGroupPads(), mGroup.mIDCsGrouped.getGroupRows(), mGroup.mIDCsGrouped.getGroupLastRowsThreshold(), mGroup.mIDCsGrouped.getGroupLastPadsThreshold(), mOverlapRows, mOverlapPads));
  delete poly;
}

template <class Type>
void o2::tpc::IDCAverageGroup::loopOverGroups(const unsigned int threadNum, const unsigned int integrationInterval, Type& type)
{
  const unsigned int region = type.mIDCsGrouped.getRegion();
  const int groupRows = static_cast<int>(type.mIDCsGrouped.getGroupRows());
  const int groupPads = static_cast<int>(type.mIDCsGrouped.getGroupPads());
  const int lastRow = static_cast<int>(type.mIDCsGrouped.getLastRow());
  unsigned int rowGrouped = 0;

  // loop over ungrouped row
  for (int iRow = 0; iRow <= lastRow; iRow += groupRows) {
    const bool bNotLastrow = iRow != lastRow;

    // the sectors is divide in to two parts around ylocal=0 to get the same simmetric grouping around ylocal=0
    for (int iYLocalSide = 0; iYLocalSide < 2; ++iYLocalSide) {
      if constexpr (std::is_same_v<Type, IDCAverageGroupType<Draw>>) {
        type.mCol = region + iRow / groupRows + iYLocalSide;
      }
      unsigned int padGrouped = iYLocalSide ? type.mIDCsGrouped.getPadsPerRow(rowGrouped) / 2 : type.mIDCsGrouped.getPadsPerRow(rowGrouped) / 2 - 1; // grouped pad in pad direction
      const int nPadsStart = Mapper::PADSPERROW[region][iRow] / 2;                                                                                   // first ungrouped pad in pad direction
      const int nPadsEnd = type.mIDCsGrouped.getLastPad(iRow) + nPadsStart;                                                                          // last grouped pad in pad direction

      // loop over ungrouped pads
      for (int iPad = nPadsStart; iPad <= nPadsEnd; iPad += groupPads) {
        if constexpr (std::is_same_v<Type, IDCAverageGroupType<Group>>) {
          type.mRobustAverage[threadNum].clear();
        }

        const int startRow = ((iRow - mOverlapRows) < 0) ? 0 : -mOverlapRows;                                                                                                          // first row in this group
        const int endRow = ((iRow + groupRows + mOverlapRows) >= Mapper::ROWSPERREGION[region] || !bNotLastrow) ? (Mapper::ROWSPERREGION[region] - iRow) : (mOverlapRows + groupRows); // last row in this group
        for (int iRowMerge = startRow; iRowMerge < endRow; ++iRowMerge) {
          const bool bOverlapRowRight = iRowMerge >= groupRows;
          const unsigned int ungroupedRow = iRow + iRowMerge;
          const int offsPad = static_cast<int>(Mapper::ADDITIONALPADSPERROW[region][ungroupedRow]) - static_cast<int>(Mapper::ADDITIONALPADSPERROW[region][iRow]); // offset due to additional pads in pad direction in the current row compared to the first row in the group

          const bool lastPad = iPad == nPadsEnd;
          const int padEnd = lastPad ? (static_cast<int>(Mapper::PADSPERROW[region][ungroupedRow]) - iPad) : (groupPads + offsPad + mOverlapPads); // last ungrouped pad in pad direction
          const int padStart = offsPad - mOverlapPads;                                                                                             // first ungrouped pad in pad direction

          for (int ipadMerge = padStart; ipadMerge < padEnd; ++ipadMerge) {
            const unsigned int ungroupedPad = iYLocalSide ? (iPad + ipadMerge) : Mapper::PADSPERROW[region][ungroupedRow] - (iPad + ipadMerge) - 1;
            const unsigned int padInRegion = Mapper::OFFSETCRULOCAL[region][ungroupedRow] + ungroupedPad;

            // averaging and grouping
            if constexpr (std::is_same_v<Type, IDCAverageGroupType<Group>>) {
              // check status flag
              auto flag = mPadStatus->getCalArray(region).getValue(padInRegion);

              // TODO add more cases...
              if (flag == PadFlags::DEAD) {
                // dead pad. just skip this
                continue;
              }

              float weight = 1;
              // set weight for outer pads which are not in the main group
              if (mOverlapRows && mOverlapPads) {
                if (iRowMerge < 0) {
                  // everything on the left border
                  const int relPosRow = std::abs(iRowMerge);
                  if (ipadMerge < offsPad) {
                    const int relPosPad = std::abs(ipadMerge - offsPad);
                    weight = mWeightsPad[std::max(relPosRow, relPosPad)];
                  } else if (!lastPad && ipadMerge >= (groupPads + offsPad)) {
                    const int relPosPad = std::abs(1 + ipadMerge - (groupPads + offsPad));
                    weight = mWeightsPad[std::max(relPosRow, relPosPad)];
                  } else {
                    weight = mWeightsPad[relPosRow];
                  }
                } else if (bNotLastrow && bOverlapRowRight) {
                  const int relPosRow = std::abs(1 + iRowMerge - (groupRows));
                  if (ipadMerge < offsPad) {
                    const int relPosPad = std::abs(ipadMerge - offsPad);
                    weight = mWeightsPad[std::max(relPosRow, relPosPad)];
                  } else if (!lastPad && ipadMerge >= (groupPads + offsPad)) {
                    const int relPosPad = std::abs(1 + ipadMerge - (groupPads + offsPad));
                    weight = mWeightsPad[std::max(relPosRow, relPosPad)];
                  } else {
                    weight = mWeightsPad[relPosRow];
                  }
                } else if (ipadMerge < offsPad) {
                  // bottom
                  const int relPadPos = std::abs(ipadMerge - offsPad);
                  weight = mWeightsPad[relPadPos];
                } else if (!lastPad && ipadMerge >= (groupPads + offsPad)) {
                  const int relPadPos = std::abs(1 + ipadMerge - (groupPads + offsPad));
                  weight = mWeightsPad[relPadPos];
                } else {
                }
              }
              const unsigned int indexIDC = integrationInterval * Mapper::PADSPERREGION[region] + padInRegion;
              type.mRobustAverage[threadNum].addValue(mIDCsUngrouped[indexIDC] * Mapper::PADAREA[region], weight);
            } else {
              // drawing
              const GlobalPadNumber padNum = o2::tpc::Mapper::getGlobalPadNumber(ungroupedRow, ungroupedPad, region);
              static auto coords = o2::tpc::painter::getPadCoordinatesSector();
              auto coordinate = coords[padNum];

              const float yPos = 0.5f * (coordinate.yVals[0] + coordinate.yVals[2]);
              const float xPos = 0.5f * (coordinate.xVals[0] + coordinate.xVals[2]);
              const int nCountDraw = type.mCountDraw[padInRegion]++;
              const float offsX = (nCountDraw % 2) * 0.6f * type.mPadInf.getPadHeight();
              const float offsY = (nCountDraw / 2) * 0.2f * type.mPadInf.getPadWidth();
              const float xPosDraw = xPos - 0.3f * type.mPadInf.getPadHeight() + offsX;
              const float yPosDraw = yPos - 0.4f * type.mPadInf.getPadWidth() + offsY;

              TLatex latex;
              latex.SetTextFont(63);
              latex.SetTextSize(1);
              const char* groupText = Form("#bf{#color[%lu]{%i}}", (type.mCol % type.mColors.size()) + 1, type.mGroupCounter);
              latex.DrawLatex(xPosDraw, yPosDraw, groupText);

              if (iRowMerge < 0 || (bNotLastrow && bOverlapRowRight) || (ipadMerge < offsPad) || (!lastPad && ipadMerge >= (groupPads + offsPad))) {
              } else {
                type.mPoly.Fill(xPos, yPos, type.mColors[type.mCol % type.mColors.size()]);
              }
            }
          }
        }

        if constexpr (std::is_same_v<Type, IDCAverageGroupType<Group>>) {
          const static auto& paramIDCGroup = ParameterIDCGroup::Instance();
          switch (paramIDCGroup.Method) {
            case o2::tpc::AveragingMethod::SLOW:
            default:
              type.mIDCsGrouped(rowGrouped, padGrouped, integrationInterval) = type.mRobustAverage[threadNum].getFilteredAverage(mSigma);
              break;
            case o2::tpc::AveragingMethod::FAST:
              type.mIDCsGrouped(rowGrouped, padGrouped, integrationInterval) = type.mRobustAverage[threadNum].getMean();
              break;
          }
        } else {
          ++type.mGroupCounter;
          ++type.mCol;
        }

        iYLocalSide ? ++padGrouped : --padGrouped;
      }
    }
    ++rowGrouped;
  }
}

void o2::tpc::IDCAverageGroup::dumpToFile(const char* outFileName, const char* outName) const
{
  TFile fOut(outFileName, "RECREATE");
  fOut.WriteObject(this, outName);
  fOut.Close();
}

bool o2::tpc::IDCAverageGroup::setFromFile(const char* fileName, const char* name)
{
  TFile inpf(fileName, "READ");
  IDCAverageGroup* idcAverageGroupTmp{nullptr};
  idcAverageGroupTmp = reinterpret_cast<IDCAverageGroup*>(inpf.GetObjectChecked(name, IDCAverageGroup::Class()));

  if (!idcAverageGroupTmp) {
    LOGP(ERROR, "Failed to load {} from {}", name, inpf.GetName());
    return false;
  }
  setIDCs(idcAverageGroupTmp->getIDCsUngrouped());

  delete idcAverageGroupTmp;
  return true;
}

void o2::tpc::IDCAverageGroup::drawUngroupedIDCs(const unsigned int integrationInterval, const std::string filename) const
{
  std::function<float(const unsigned int, const unsigned int, const unsigned int, const unsigned int)> idcFunc = [this, integrationInterval](const unsigned int, const unsigned int region, const unsigned int irow, const unsigned int pad) {
    const unsigned int indexIDC = integrationInterval * Mapper::PADSPERREGION[region] + Mapper::OFFSETCRULOCAL[region][irow] + pad;
    return mIDCsUngrouped[indexIDC] * Mapper::PADAREA[region];
  };

  IDCDrawHelper::IDCDraw drawFun;
  drawFun.mIDCFunc = idcFunc;
  const std::string zAxisTitle = IDCDrawHelper::getZAxisTitle(IDCType::IDC);
  IDCDrawHelper::drawSector(drawFun, mGroup.mIDCsGrouped.getRegion(), mGroup.mIDCsGrouped.getRegion() + 1, 0, zAxisTitle, filename);
}

/// for debugging: creating debug tree for integrated IDCs
/// \param nameFile name of the output file
void o2::tpc::IDCAverageGroup::createDebugTree(const char* nameFile) const
{
  o2::utils::TreeStreamRedirector pcstream(nameFile, "RECREATE");
  pcstream.GetFile()->cd();
  createDebugTree(*this, pcstream);
  pcstream.Close();
}

void o2::tpc::IDCAverageGroup::createDebugTreeForAllCRUs(const char* nameFile, const char* filename)
{
  o2::utils::TreeStreamRedirector pcstream(nameFile, "RECREATE");
  pcstream.GetFile()->cd();
  TFile fInp(filename, "READ");

  for (TObject* keyAsObj : *fInp.GetListOfKeys()) {
    const auto key = dynamic_cast<TKey*>(keyAsObj);
    LOGP(info, "Key name: {} Type: {}", key->GetName(), key->GetClassName());

    if (std::strcmp(o2::tpc::IDCAverageGroup::Class()->GetName(), key->GetClassName()) != 0) {
      LOGP(info, "skipping object. wrong class.");
      continue;
    }

    IDCAverageGroup* idcavg = (IDCAverageGroup*)fInp.Get(key->GetName());
    createDebugTree(*idcavg, pcstream);
    delete idcavg;
  }
  pcstream.Close();
}

void o2::tpc::IDCAverageGroup::createDebugTree(const IDCAverageGroup& idcavg, o2::utils::TreeStreamRedirector& pcstream)
{
  const Mapper& mapper = Mapper::instance();
  unsigned int sector = idcavg.getSector();
  unsigned int cru = sector * Mapper::NREGIONS + idcavg.getRegion();
  const o2::tpc::CRU cruTmp(cru);
  unsigned int region = cruTmp.region();

  for (unsigned int integrationInterval = 0; integrationInterval < idcavg.getNIntegrationIntervals(); ++integrationInterval) {
    const unsigned long padsPerCRU = Mapper::PADSPERREGION[region];
    std::vector<unsigned int> vRow(padsPerCRU);
    std::vector<unsigned int> vPad(padsPerCRU);
    std::vector<float> vXPos(padsPerCRU);
    std::vector<float> vYPos(padsPerCRU);
    std::vector<float> vGlobalXPos(padsPerCRU);
    std::vector<float> vGlobalYPos(padsPerCRU);
    std::vector<float> idcsPerIntegrationInterval(padsPerCRU);        // idcs for one time bin
    std::vector<float> groupedidcsPerIntegrationInterval(padsPerCRU); // idcs for one time bin
    std::vector<float> invPadArea(padsPerCRU);

    for (unsigned int iPad = 0; iPad < padsPerCRU; ++iPad) {
      const GlobalPadNumber globalNum = Mapper::GLOBALPADOFFSET[region] + iPad;
      const auto& padPosLocal = mapper.padPos(globalNum);
      vRow[iPad] = padPosLocal.getRow();
      vPad[iPad] = padPosLocal.getPad();
      vXPos[iPad] = mapper.getPadCentre(padPosLocal).X();
      vYPos[iPad] = mapper.getPadCentre(padPosLocal).Y();
      invPadArea[iPad] = Mapper::PADAREA[region];
      const GlobalPosition2D globalPos = mapper.LocalToGlobal(LocalPosition2D(vXPos[iPad], vYPos[iPad]), cruTmp.sector());
      vGlobalXPos[iPad] = globalPos.X();
      vGlobalYPos[iPad] = globalPos.Y();
      idcsPerIntegrationInterval[iPad] = idcavg.getUngroupedIDCVal(iPad, integrationInterval);
      groupedidcsPerIntegrationInterval[iPad] = idcavg.getGroupedIDCValGlobal(vRow[iPad], vPad[iPad], integrationInterval);
    }

    pcstream << "tree"
             << "cru=" << cru
             << "sector=" << sector
             << "region=" << region
             << "integrationInterval=" << integrationInterval
             << "IDCUngrouped.=" << idcsPerIntegrationInterval
             << "IDCGrouped.=" << groupedidcsPerIntegrationInterval
             << "invPadArea.=" << invPadArea
             << "pad.=" << vPad
             << "row.=" << vRow
             << "lx.=" << vXPos
             << "ly.=" << vYPos
             << "gx.=" << vGlobalXPos
             << "gy.=" << vGlobalYPos
             << "\n";
  }
}

void o2::tpc::IDCAverageGroup::setIDCs(const std::vector<float>& idcs)
{
  mIDCsUngrouped = idcs;
  mGroup.mIDCsGrouped.resize(getNIntegrationIntervals());
}

void o2::tpc::IDCAverageGroup::setIDCs(std::vector<float>&& idcs)
{
  mIDCsUngrouped = std::move(idcs);
  mGroup.mIDCsGrouped.resize(getNIntegrationIntervals());
}

unsigned int o2::tpc::IDCAverageGroup::getNIntegrationIntervals() const
{
  return mIDCsUngrouped.size() / Mapper::PADSPERREGION[mGroup.mIDCsGrouped.getRegion()];
}

float o2::tpc::IDCAverageGroup::getUngroupedIDCVal(const unsigned int localPadNumber, const unsigned int integrationInterval) const
{
  return mIDCsUngrouped[localPadNumber + integrationInterval * Mapper::PADSPERREGION[mGroup.mIDCsGrouped.getRegion()]];
}

unsigned int o2::tpc::IDCAverageGroup::getUngroupedIndex(const unsigned int ulrow, const unsigned int upad, const unsigned int integrationInterval) const
{
  return integrationInterval * Mapper::PADSPERREGION[mGroup.mIDCsGrouped.getRegion()] + Mapper::OFFSETCRULOCAL[mGroup.mIDCsGrouped.getRegion()][ulrow] + upad;
}

unsigned int o2::tpc::IDCAverageGroup::getUngroupedIndexGlobal(const unsigned int ugrow, const unsigned int upad, const unsigned int integrationInterval) const
{
  return integrationInterval * Mapper::PADSPERREGION[mGroup.mIDCsGrouped.getRegion()] + Mapper::OFFSETCRUGLOBAL[ugrow] + upad;
}

template void o2::tpc::IDCAverageGroup::loopOverGroups(const unsigned int, const unsigned int, IDCAverageGroupType<Group>&);
template void o2::tpc::IDCAverageGroup::loopOverGroups(const unsigned int, const unsigned int, IDCAverageGroupType<Draw>&);
