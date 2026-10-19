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

/// \file cosmicsdEdx.C
/// \brief Recalculate the TPC dE/dx of the cosmics in o2_cosmics_full.root (o2-cosmics-match-workflow --enable-cluster-output)
///        with o2::tpc::CalculatedEdx, per leg and for both legs combined in one truncated mean.
///
/// The TPC-only legs carry the z of the tracker's time0 guess. For cosmics with a TOF time (hasTOFTime()) every one-sided leg is
/// moved to that time before the dE/dx calculation: time0 = TOF time and z shifted by side * (time0 new - time0 old) * vDrift, so
/// the z-dependent corrections (topology, space charge) see the right drift length. Legs with clusters on both sides keep their
/// time0 (fixed by the central electrode). Cosmics without a TOF time are calculated at the reconstructed time.
///
/// Cluster sets:
///   A: the clusters attached to the leg's TPC track
///   B: A + road clusters (not attached to any track, same TPC side) filling row gaps of up to maxGapRows rows between two attached
///      clusters of the same sector: per gap row the road cluster closest to the straight-line interpolation of the two attached
///      clusters (raw pad-y and time, so the distortions cancel), within maxDy [cm] and maxDTime [time bins]
///
/// Output: friend tree "dedx" with one entry per TF (same order as cosmicsFull) and per cosmic one element in each vector:
///   totA, maxA                    combined truncated mean of both legs (qTot, qMax), -1 if no leg
///   totBotA, totTopA, maxBotA, maxTopA, nClBotA, nClTopA   per leg (-1 if the leg has no TPC part), nCl = accepted clusters
///   same with B, nGapFillBot / nGapFillTop = road clusters added in B
///   tofTime                       1 if the legs were moved to the TOF time
///
///   cosmicsFull->AddFriend("dedx", "cosmics_dedx.root");
///   cosmicsFull->Draw("dedx.totA : log(cosmics.cosmic.getP())>>g(200,-1,10,100,0,150)",
///                     "cosmics.scoreTOFPair>=0 && (cosmics.cosmic.mFlags&3)==0 && cosmics.duplicateOf<0", "colz");

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "Framework/Logger.h"
#include "DataFormatsGlobalTracking/TrackCosmicsExtended.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "DataFormatsTPC/ClusterNative.h"
#include "TPCBase/ParameterElectronics.h"
#include "GPUTPCGeometry.h"
#include "TPCCalibration/CalculatedEdx.h"
#endif

using namespace o2::tpc;

namespace
{
/// +1 A side, -1 C side, 0 clusters on both sides
int tpcSide(const TrackTPC& track)
{
  return track.hasASideClustersOnly() ? 1 : (track.hasCSideClustersOnly() ? -1 : 0);
}

/// truncated mean of the merged corrected charges [4] (all regions) of the legs
float combinedTruncMean(const CalculatedEdx& calc, const std::vector<const std::array<std::vector<float>, 5>*>& legCharges, float low, float high)
{
  std::vector<float> charges;
  for (const auto* legCharge : legCharges) {
    charges.insert(charges.end(), (*legCharge)[4].begin(), (*legCharge)[4].end());
  }
  if (charges.empty()) {
    return -1.f;
  }
  return calc.getTruncMean(charges, low, high);
}

struct LegResult {
  float tot = -1.f;
  float max = -1.f;
  int nCl = -1;
  int nGapFill = 0;
  std::array<std::vector<float>, 5> chargesTot;
  std::array<std::vector<float>, 5> chargesMax;
};

/// attached clusters (A) and, for B, the gap-filling road clusters in the track's row order
void buildClusters(const std::vector<o2::dataformats::CosmicTPCCluster>& legClusters, bool fillGaps, int maxGapRows, float maxDy, float maxDTime,
                   std::vector<ClusterNative>& clusters, ClInfoVec& clusterInfos, int& nGapFill)
{
  const o2::gpu::GPUTPCGeometry gpuGeom;
  std::vector<const o2::dataformats::CosmicTPCCluster*> attached;
  std::vector<const o2::dataformats::CosmicTPCCluster*> road;
  for (const auto& cl : legClusters) {
    if (cl.isAttached()) {
      attached.push_back(&cl);
    } else if (cl.isCorridor() && !cl.isUsed() && !cl.isAbsTime()) {
      road.push_back(&cl);
    }
  }

  // position along the attached list: attached cluster i at i, a gap cluster between i and i+1 proportionally to its row
  std::vector<std::pair<float, const o2::dataformats::CosmicTPCCluster*>> ordered;
  for (size_t i = 0; i < attached.size(); i++) {
    ordered.emplace_back(static_cast<float>(i), attached[i]);
  }
  nGapFill = 0;
  if (fillGaps) {
    for (size_t i = 0; i + 1 < attached.size(); i++) {
      const auto& first = *attached[i];
      const auto& second = *attached[i + 1];
      const int rowDiff = std::abs(static_cast<int>(second.row) - static_cast<int>(first.row));
      if (first.sector != second.sector || rowDiff < 2 || rowDiff > maxGapRows + 1) {
        continue;
      }
      const float xFirst = gpuGeom.Row2X(first.row);
      const float xSecond = gpuGeom.Row2X(second.row);
      const float yFirst = gpuGeom.LinearPad2Y(first.sector, first.row, first.cl.getPad());
      const float ySecond = gpuGeom.LinearPad2Y(second.sector, second.row, second.cl.getPad());
      const int rowStep = (second.row > first.row) ? 1 : -1;
      for (int row = first.row + rowStep; row != second.row; row += rowStep) {
        const float frac = (gpuGeom.Row2X(row) - xFirst) / (xSecond - xFirst);
        const float yExpected = yFirst + frac * (ySecond - yFirst);
        const float timeExpected = first.cl.getTime() + frac * (second.cl.getTime() - first.cl.getTime());
        const o2::dataformats::CosmicTPCCluster* best = nullptr;
        float bestDistance2 = 1.e9f;
        for (const auto* cand : road) {
          if (cand->sector != first.sector || cand->row != row) {
            continue;
          }
          const float dy = gpuGeom.LinearPad2Y(cand->sector, cand->row, cand->cl.getPad()) - yExpected;
          const float dTime = cand->cl.getTime() - timeExpected;
          if (std::abs(dy) > maxDy || std::abs(dTime) > maxDTime) {
            continue;
          }
          const float distance2 = (dy / maxDy) * (dy / maxDy) + (dTime / maxDTime) * (dTime / maxDTime);
          if (distance2 < bestDistance2) {
            bestDistance2 = distance2;
            best = cand;
          }
        }
        if (best) {
          ordered.emplace_back(static_cast<float>(i) + std::abs(row - static_cast<int>(first.row)) / static_cast<float>(rowDiff), best);
          nGapFill++;
        }
      }
    }
    std::stable_sort(ordered.begin(), ordered.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  }

  clusters.clear();
  clusterInfos.clear();
  for (const auto& [position, cl] : ordered) {
    clusters.push_back(cl->cl);
    ClInfo info;
    info.sectorIndex = cl->sector;
    info.rowIndex = cl->row;
    info.isShared = false; // sharing with other tracks is not stored; road clusters attached to another track are excluded above
    clusterInfos.push_back(info);
  }
}

LegResult calculateLeg(CalculatedEdx& calc, const TrackTPC& legTrack, const std::vector<o2::dataformats::CosmicTPCCluster>& legClusters, bool fillGaps,
                       float time0TOF, float vDrift, const dEdxSettings& settings, int maxGapRows, float maxDy, float maxDTime)
{
  LegResult result;
  if (legTrack.getNClusters() == 0) { // no TPC part
    return result;
  }
  std::vector<ClusterNative> clusters;
  ClInfoVec clusterInfos;
  buildClusters(legClusters, fillGaps, maxGapRows, maxDy, maxDTime, clusters, clusterInfos, result.nGapFill);
  if (clusters.empty()) {
    return result;
  }

  TrackTPC track = legTrack;
  const int side = tpcSide(track);
  if (time0TOF > -1.e8f && side != 0) {
    track.setZ(track.getZ() + side * (time0TOF - track.getTime0()) * vDrift);
    track.setTime0(time0TOF);
  }

  std::vector<dEdxInfo> outputs;
  AverageOccupancy averageOcc;
  calc.calculatedEdxMultipleSettings(track, clusters, clusterInfos, outputs, averageOcc, {settings});
  const auto& output = outputs.front();
  result.tot = output.dEdxTotTPC;
  result.max = output.dEdxMaxTPC;
  result.nCl = output.NHitsIROC + output.NHitsOROC1 + output.NHitsOROC2 + output.NHitsOROC3;
  result.chargesTot = calc.getClusterChargesTot().front();
  result.chargesMax = calc.getClusterChargesMax().front();
  return result;
}
} // namespace

/// \param inFile o2_cosmics_full.root (tree cosmicsFull)
/// \param outFile output file with the friend tree "dedx"
/// \param runNumberOrTimeStamp run number (calibrations at the middle of the run) or time stamp for the CCDB
/// \param firstEntry first TF (entry of cosmicsFull) to process
/// \param lastEntry process the TFs up to this entry, exclusive (-1: to the end); pieces of consecutive ranges merged with hadd in
///                  order form the friend tree of the whole file
/// \param maxGapRows B: longest row gap between two attached clusters filled with road clusters
/// \param maxDy B: max distance in pad direction of a road cluster to the interpolation [cm]
/// \param maxDTime B: max distance in time of a road cluster to the interpolation [time bins]
void cosmicsdEdx(const std::string inFile = "o2_cosmics_full.root",
                 const std::string outFile = "cosmics_dedx.root",
                 const long runNumberOrTimeStamp = 568041,
                 const long long firstEntry = 0,
                 const long long lastEntry = -1,
                 const int maxGapRows = 5,
                 const float maxDy = 0.3f,
                 const float maxDTime = 2.f)
{
  std::unique_ptr<TFile> fileIn(TFile::Open(inFile.data()));
  if (!fileIn || fileIn->IsZombie()) {
    LOGP(error, "Could not open {}", inFile);
    return;
  }
  auto* treeIn = static_cast<TTree*>(fileIn->Get("cosmicsFull"));
  if (!treeIn) {
    LOGP(error, "No tree cosmicsFull in {}", inFile);
    return;
  }
  std::vector<o2::dataformats::TrackCosmicsExtended>* cosmics = nullptr;
  o2::dataformats::CosmicsTFInfo* tfInfo = nullptr;
  treeIn->SetBranchAddress("cosmics", &cosmics);
  treeIn->SetBranchAddress("tfInfo", &tfInfo);

  dEdxSettings settings;
  settings.low = 0.015f;
  settings.high = 0.6f;
  settings.correctionMask = CorrectionFlags::TopologyPol | CorrectionFlags::dEdxResidual;
  settings.clusterMask = ClusterFlags::ExcludeEdgeCl;
  settings.subthresholdMethod = 0;
  settings.stackBoundaryMethod = 0;
  settings.sameRowClusterMethod = 2;

  CalculatedEdx calc;
  calc.loadCalibsFromCCDB(runNumberOrTimeStamp, false, false, false, false);
  calc.setPropagateTrack(true); // from the legs' own parameters (fitted with the distortion corrections)
  calc.setKeepClusterCharges(true);

  const float tbWidthMUS = ParameterElectronics::Instance().ZbinWidth;

  std::unique_ptr<TFile> fileOut(TFile::Open(outFile.data(), "RECREATE"));
  TTree treeOut("dedx", "recalculated dE/dx of the cosmics, friend of cosmicsFull");
  std::vector<float> totA;
  std::vector<float> maxA;
  std::vector<float> totBotA;
  std::vector<float> totTopA;
  std::vector<float> maxBotA;
  std::vector<float> maxTopA;
  std::vector<int> nClBotA;
  std::vector<int> nClTopA;
  std::vector<float> totB;
  std::vector<float> maxB;
  std::vector<float> totBotB;
  std::vector<float> totTopB;
  std::vector<float> maxBotB;
  std::vector<float> maxTopB;
  std::vector<int> nClBotB;
  std::vector<int> nClTopB;
  std::vector<int> nGapFillBot;
  std::vector<int> nGapFillTop;
  std::vector<int> tofTime;
  treeOut.Branch("totA", &totA);
  treeOut.Branch("maxA", &maxA);
  treeOut.Branch("totBotA", &totBotA);
  treeOut.Branch("totTopA", &totTopA);
  treeOut.Branch("maxBotA", &maxBotA);
  treeOut.Branch("maxTopA", &maxTopA);
  treeOut.Branch("nClBotA", &nClBotA);
  treeOut.Branch("nClTopA", &nClTopA);
  treeOut.Branch("totB", &totB);
  treeOut.Branch("maxB", &maxB);
  treeOut.Branch("totBotB", &totBotB);
  treeOut.Branch("totTopB", &totTopB);
  treeOut.Branch("maxBotB", &maxBotB);
  treeOut.Branch("maxTopB", &maxTopB);
  treeOut.Branch("nClBotB", &nClBotB);
  treeOut.Branch("nClTopB", &nClTopB);
  treeOut.Branch("nGapFillBot", &nGapFillBot);
  treeOut.Branch("nGapFillTop", &nGapFillTop);
  treeOut.Branch("tofTime", &tofTime);

  const long long endEntry = (lastEntry >= 0) ? std::min<long long>(lastEntry, treeIn->GetEntries()) : treeIn->GetEntries();
  long long nCosmics = 0;
  calc.resetDebugCounters();
  for (long long iEntry = firstEntry; iEntry < endEntry; iEntry++) {
    treeIn->GetEntry(iEntry);
    for (auto* vec : {&totA, &maxA, &totBotA, &totTopA, &maxBotA, &maxTopA, &totB, &maxB, &totBotB, &totTopB, &maxBotB, &maxTopB}) {
      vec->clear();
    }
    for (auto* vec : {&nClBotA, &nClTopA, &nClBotB, &nClTopB, &nGapFillBot, &nGapFillTop, &tofTime}) {
      vec->clear();
    }

    for (const auto& cosmic : *cosmics) {
      const float time0TOF = cosmic.hasTOFTime() ? cosmic.timeTOFMUS / tbWidthMUS : -1.e9f;
      tofTime.push_back(cosmic.hasTOFTime());
      for (const bool fillGaps : {false, true}) {
        const auto bottom = calculateLeg(calc, cosmic.tpcBottom, cosmic.clTPCBottom, fillGaps, time0TOF, tfInfo->vDrift, settings, maxGapRows, maxDy, maxDTime);
        const auto top = calculateLeg(calc, cosmic.tpcTop, cosmic.clTPCTop, fillGaps, time0TOF, tfInfo->vDrift, settings, maxGapRows, maxDy, maxDTime);
        std::vector<const std::array<std::vector<float>, 5>*> chargesTot;
        std::vector<const std::array<std::vector<float>, 5>*> chargesMax;
        for (const auto* leg : {&bottom, &top}) {
          if (leg->nCl > 0) {
            chargesTot.push_back(&leg->chargesTot);
            chargesMax.push_back(&leg->chargesMax);
          }
        }
        (fillGaps ? totB : totA).push_back(combinedTruncMean(calc, chargesTot, settings.low, settings.high));
        (fillGaps ? maxB : maxA).push_back(combinedTruncMean(calc, chargesMax, settings.low, settings.high));
        (fillGaps ? totBotB : totBotA).push_back(bottom.tot);
        (fillGaps ? totTopB : totTopA).push_back(top.tot);
        (fillGaps ? maxBotB : maxBotA).push_back(bottom.max);
        (fillGaps ? maxTopB : maxTopA).push_back(top.max);
        (fillGaps ? nClBotB : nClBotA).push_back(bottom.nCl);
        (fillGaps ? nClTopB : nClTopA).push_back(top.nCl);
        if (fillGaps) {
          nGapFillBot.push_back(bottom.nGapFill);
          nGapFillTop.push_back(top.nGapFill);
        }
      }
      nCosmics++;
    }
    treeOut.Fill();
    if ((iEntry - firstEntry) % 10000 == 0) {
      LOGP(info, "TF {} / [{}, {}): {} cosmics, rows {}, failed propagations {}", iEntry, firstEntry, endEntry, nCosmics, calc.getNRowsProcessed(), calc.getNPropagationFailed());
    }
  }
  LOGP(info, "Done: TFs [{}, {}), {} cosmics, rows {}, failed propagations {}", firstEntry, endEntry, nCosmics, calc.getNRowsProcessed(), calc.getNPropagationFailed());
  fileOut->cd();
  treeOut.Write();
}
