// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   CosmicsClusterCollectorSpec.cxx
/// @brief  Collects the raw clusters of matched cosmics (attached + road around the legs) for offline refits
///
/// The async reconstruction stores neither TPC tracks nor TPC clusters, so the leg references of TrackCosmics cannot be resolved
/// offline. For every cosmic this device stores the TPC tracks of both legs, the raw TPC clusters attached to them, all raw TPC clusters
/// in a road around each leg (gaps, split pieces, delta electrons; clusters attached to other tracks are flagged), the ITS clusters,
/// matched TOF clusters and TRD tracklets of the legs, and per TF the quantities of the TPC transformation.
///
/// Road: each leg's helix is propagated through all pad rows of all sectors on its TPC side, in the frame of its own time0 (the frame
/// in which its clusters are consistent). Points on the other side of the leg's closest approach to the beam line belong to the other
/// leg and are skipped. The predicted real (y, z) is mapped to nominal coordinates with the inverse correction, and clusters with
/// nominal coordinates within the corridor width of the predicted point (distance perpendicular to the track) are taken. The other TPC
/// side is searched only when the time of the cosmic is known (CE-crossing leg, or a small time error of the cosmic).
///
/// Roads in the other detectors (--road-detectors): TOF clusters and TRD tracklets close to the outward extrapolation of the legs, and
/// ITS clusters close to the trajectory near the beam line, within the time window of the cosmic (its time from the matcher, which is
/// precise for TPC-only legs on opposite sides; otherwise the time window of the legs, with a correspondingly loose z cut). They are
/// flagged as found on the road; hits of the legs' matched global tracks are flagged as matched.
///
/// Only raw detector data are stored (TPC ClusterNative, ITS compact clusters + patterns, raw TOF time, TRD tracklet words); calibrations,
/// the cluster dictionary and the geometry are applied offline.
///
/// --debug-tree writes cosmics_collector_debug.root for test runs: tree "cl" with every stored TPC cluster transformed (local x, y, z in the
/// frame used by the road, zCos in the common frame of the cosmic's time, global gx, gy), tree "road" with the predicted track points of
/// each leg (row by row, same frames), tree "its" with the local (with the ITS road also global) coordinates of the ITS clusters, tree
/// "tof" with the raw and calibrated TOF time and the cluster position, tree "trd" with the road tracklets and their residuals, and tree
/// "cosm" with one entry per cosmic.

#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include "TStopwatch.h"
#include "Framework/Task.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/DeviceSpec.h"
#include "GlobalTrackingWorkflow/CosmicsClusterCollectorSpec.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DataFormatsGlobalTracking/CosmicTrack.h"
#include "ReconstructionDataFormats/TrackCosmics.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "DataFormatsTPC/ClusterNative.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsITSMFT/ClusterPattern.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "DataFormatsTOF/Cluster.h"
#include "DataFormatsTRD/TrackTRD.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include "DataFormatsTRD/CalibratedTracklet.h"
#include "DataFormatsTRD/Constants.h"
#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "ITSBase/GeometryTGeo.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "CommonConstants/LHCConstants.h"
#include "CommonDataFormat/TFIDInfo.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/TFIDInfoHelper.h"
#include "TPCBase/ParameterElectronics.h"
#include "MathUtils/Utils.h"
#include "MathUtils/Primitive2D.h"
#include "TPCFastTransformPOD.h"
#include "CommonUtils/TreeStreamRedirector.h"

using namespace o2::framework;
using GTrackID = o2::dataformats::GlobalTrackID;
using TPCGeo = o2::gpu::TPCFastTransformGeoPOD;
using DetID = o2::detectors::DetID;

namespace o2::globaltracking
{

class CosmicsClusterCollectorSpec : public Task
{
 public:
  CosmicsClusterCollectorSpec(std::shared_ptr<DataRequest> dr, std::shared_ptr<o2::base::GRPGeomRequest> gr, bool useMC, DetID::mask_t roadDets) : mDataRequest(dr), mGGCCDBRequest(gr), mRoadDets(roadDets), mUseMC(useMC) {}
  ~CosmicsClusterCollectorSpec() override = default;
  void init(InitContext& ic) final;
  void run(ProcessingContext& pc) final;
  void endOfStream(EndOfStreamContext& ec) final;
  void finaliseCCDB(ConcreteDataMatcher& matcher, void* obj) final;

 private:
  /// which side of the closest approach to the beam line (transverse) a point is on, relative to the leg's own clusters
  struct LegBranch {
    bool isLine = false; ///< straight line (no field / very high pT) instead of a circle
    float centerX = 0.f; ///< circle centre
    float centerY = 0.f;
    float pcaX = 0.f; ///< point of closest approach to the beam line
    float pcaY = 0.f;
    float dirX = 0.f; ///< direction of the straight line
    float dirY = 0.f;
    int sign = 0; ///< side of the leg's clusters; 0: the leg spans both sides, accept everything
    int side(float x, float y) const
    {
      const float orientation = isLine ? (x - pcaX) * dirX + (y - pcaY) * dirY : (pcaX - centerX) * (y - centerY) - (pcaY - centerY) * (x - centerX);
      return orientation > 0.f ? 1 : -1;
    }
    bool accept(float x, float y) const { return sign == 0 || side(x, y) == sign; }
    void init(const o2::track::TrackPar& inner, const o2::track::TrackPar& outer, float bz);
  };

  /// time of the cosmic in TPC time bins
  struct CosmicTime {
    bool known = false; ///< known well enough to search the other TPC side of one-side legs
    float tb = 0.f;     ///< time [TB]
    float errTB = 0.f;  ///< its error [TB]
  };

  void updateTimeDependentParams(ProcessingContext& pc);
  void buildUsedMap(const RecoContainer& data, std::vector<std::pair<float, float>>& time0Windows);
  void addTPCAttached(const RecoContainer& data, const o2::tpc::TrackTPC& trk, std::vector<o2::dataformats::CosmicTPCCluster>& out, std::unordered_set<uint32_t>& taken) const;
  void addTPCCorridor(const RecoContainer& data, const o2::tpc::TrackTPC& trk, const CosmicTime& cosmicTime, std::vector<o2::dataformats::CosmicTPCCluster>& out, std::unordered_set<uint32_t>& taken) const;
  void searchRow(const o2::tpc::ClusterNativeAccess& clusters, int sector, int row, float y, float z, float snp, float tgl, float vertexTime, float zTolerance, uint8_t flag,
                 std::vector<o2::dataformats::CosmicTPCCluster>& out, std::unordered_set<uint32_t>& taken) const;
  struct ITSPattRequest {
    int cosmic;  ///< entry of the cosmic in the output
    int cluster; ///< entry of the cluster in its clITS
    int index;   ///< index of the cluster in the TF
  };
  void addITS(const RecoContainer& data, GTrackID gid, uint8_t leg, std::vector<o2::dataformats::CosmicITSCluster>& out, int icosm, std::vector<ITSPattRequest>& requests, std::unordered_set<int>& matched) const;
  void fillITSPatterns(const RecoContainer& data, std::vector<ITSPattRequest>& requests, std::vector<o2::dataformats::CosmicTrack>& cosmics) const;
  void addTOF(const RecoContainer& data, GTrackID gid, uint8_t leg, std::vector<o2::dataformats::CosmicTOFCluster>& out, std::unordered_set<int>& matched) const;
  void addTRD(const RecoContainer& data, GTrackID gid, uint8_t leg, std::vector<o2::dataformats::CosmicTRDTracklet>& out, std::unordered_set<int>& matched) const;
  std::pair<float, float> timeWindowMUS(const CosmicTime& cosmicTime) const;
  bool predictOutward(const o2::tpc::TrackTPC& leg, const CosmicTime& cosmicTime, int sector, float x, float& y, float& z) const;
  void roadTOF(const RecoContainer& data, const o2::tpc::TrackTPC* const* legs, const CosmicTime& cosmicTime, int icosm, std::vector<o2::dataformats::CosmicTOFCluster>& out, const std::unordered_set<int>& matched, float& timeTOFMUS) const;
  void roadTRD(const RecoContainer& data, const o2::tpc::TrackTPC* const* legs, const CosmicTime& cosmicTime, int icosm, std::vector<o2::dataformats::CosmicTRDTracklet>& out, const std::unordered_set<int>& matched) const;
  void roadITS(const RecoContainer& data, const o2::dataformats::TrackCosmics& cosm, const CosmicTime& cosmicTime, int legsSide, int icosm, std::vector<o2::dataformats::CosmicITSCluster>& out, std::vector<ITSPattRequest>& requests,
               const std::unordered_set<int>& matched) const;
  void cacheITSChipCentres();
  void writeDebug(const o2::dataformats::CosmicTrack& cosm, int icosm) const;
  void writeDebugTOF(const o2::tof::Cluster& c, int icosm, int leg, uint8_t flags) const;

  std::shared_ptr<DataRequest> mDataRequest;
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  const o2::gpu::TPCFastTransformPOD* mCorrMap = nullptr;
  const o2::itsmft::TopologyDictionary* mITSDict = nullptr;
  std::vector<bool> mUsed;                                     ///< TPC clusters attached to any TPC track in this TF
  o2::InteractionRecord mTFStart{};                            ///< first BC of the TF
  float mCorridor = 1.f;                                       ///< road half-width [cm]
  float mMaxAbsTimeErr = 0.5f;                                 ///< max. time error of the cosmic [mus] to search the other TPC side of a one-side leg
  size_t mMaxCosmicsPerTF = 100;                               ///< cosmics processed per TF at most (protection against fake-dominated settings)
  DetID::mask_t mRoadDets{};                                   ///< detectors searched along the road besides the TPC
  float mRoadTOF = 5.f;                                        ///< road half-width at the TOF [cm]
  float mTOFFlightTol = 2.f;                                   ///< max. deviation of the top/bottom TOF time difference from the flight time [ns]
  float mTOFTimeErr = 0.1f;                                    ///< error of the TOF time of a cosmic for the later roads [mus] (covers TPC vs TOF offsets)
  float mRoadTRD = 5.f;                                        ///< road half-width at the TRD [cm] (in z plus the pad length)
  float mRoadITS = 1.5f;                                       ///< road half-width in the ITS [cm]
  std::vector<o2::math_utils::Point3D<float>> mITSChipCentres; ///< global positions of the ITS chip centres (aligned geometry)
  float mTPCTBinMUS = 0.2f;                                    ///< TPC time bin [mus]
  float mBz = 0.f;
  bool mUseMC = false;
  mutable bool mStaggeredWarned = false;
  mutable bool mNoDictWarned = false;
  std::unique_ptr<o2::utils::TreeStreamRedirector> mDebugOut; ///< debug trees with transformed clusters and road points (--debug-tree)
  int mDbgTF = 0;                                             ///< context of the debug output: TF counter
  int mDbgCosmic = 0;                                         ///< entry of the cosmic in the TF
  int mDbgLeg = 0;                                            ///< leg (0 bottom, 1 top)
  float mDbgTauC = 0.f;                                       ///< time of the cosmic [TB]: common frame (zCos) of the debug output
  size_t mNCosmics = 0;                                       ///< cosmics processed
  size_t mNClAttached = 0;                                    ///< attached TPC clusters stored
  size_t mNClCorridor = 0;                                    ///< road TPC clusters stored
  TStopwatch mTimer;
};

void CosmicsClusterCollectorSpec::init(InitContext& ic)
{
  mTimer.Stop();
  mTimer.Reset();
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);
  mCorridor = ic.options().get<float>("corridor-width");
  mMaxAbsTimeErr = ic.options().get<float>("max-abs-time-err");
  mMaxCosmicsPerTF = ic.options().get<int>("max-cosmics-per-tf");
  mRoadTOF = ic.options().get<float>("tof-road-width");
  mTOFFlightTol = ic.options().get<float>("tof-flight-tolerance");
  mTOFTimeErr = ic.options().get<float>("tof-time-error");
  mRoadTRD = ic.options().get<float>("trd-road-width");
  mRoadITS = ic.options().get<float>("its-road-width");
  if (ic.options().get<bool>("debug-tree")) {
    const auto timesliceId = ic.services().get<const o2::framework::DeviceSpec>().inputTimesliceId;
    const std::string name = timesliceId == 0 ? "cosmics_collector_debug.root" : fmt::format("cosmics_collector_debug_{}.root", timesliceId);
    mDebugOut = std::make_unique<o2::utils::TreeStreamRedirector>(name.c_str(), "recreate");
  }
}

void CosmicsClusterCollectorSpec::run(ProcessingContext& pc)
{
  mTimer.Start(false);
  RecoContainer recoData;
  recoData.collectData(pc, *mDataRequest.get());
  updateTimeDependentParams(pc);

  o2::dataformats::TFIDInfo tfID;
  o2::base::TFIDInfoHelper::fillTFIDInfo(pc, tfID);
  mTFStart = {0, tfID.firstTForbit};
  o2::dataformats::CosmicsTFInfo tfInfo;
  tfInfo.vDrift = mCorrMap->getVDrift();
  tfInfo.t0 = mCorrMap->getT0();

  std::vector<o2::dataformats::CosmicTrack> cosmicsOut;
  const auto cosmics = recoData.getCosmicTracks();
  const size_t nCosmics = std::min(cosmics.size(), mMaxCosmicsPerTF);
  if (nCosmics < cosmics.size()) {
    LOGP(warning, "{} cosmics in TF {}, clusters collected for the first {} only", cosmics.size(), tfID.tfCounter, nCosmics);
  }
  if (nCosmics) {
    // only TPC tracks whose time0 is within two drift times of a leg's time0 can share clusters with the roads
    const float maxDistTB = 2.2f * TPCGeo::getTPCzLength() / mCorrMap->getVDrift();
    std::vector<std::pair<float, float>> time0Windows;
    for (size_t ic = 0; ic < nCosmics; ic++) {
      for (const auto leg : {cosmics[ic].getRefBottom(), cosmics[ic].getRefTop()}) {
        const auto refs = recoData.getSingleDetectorRefs(leg);
        if (refs[GTrackID::TPC].isIndexSet()) {
          const float time0 = recoData.getTPCTrack(refs[GTrackID::TPC]).getTime0();
          time0Windows.emplace_back(time0 - maxDistTB, time0 + maxDistTB);
        }
      }
    }
    buildUsedMap(recoData, time0Windows);
  }
  std::vector<ITSPattRequest> pattRequests; // ITS clusters whose patterns are not in the dictionary
  mDbgTF = tfID.tfCounter;
  for (size_t ic = 0; ic < nCosmics; ic++) {
    const auto& cosm = cosmics[ic];
    auto& out = cosmicsOut.emplace_back();
    out.cosmic = cosm;
    if (mUseMC) {
      out.label = recoData.getCosmicTrackMCLabel(ic);
    }
    const GTrackID legs[2] = {cosm.getRefBottom(), cosm.getRefTop()};
    const o2::tpc::TrackTPC* tpcLegs[2] = {nullptr, nullptr};
    std::vector<o2::dataformats::CosmicTPCCluster>* tpcCl[2] = {&out.clTPCBottom, &out.clTPCTop};
    std::unordered_set<int> matchedITS; // hits of the legs' matched tracks, not searched again on the road
    std::unordered_set<int> matchedTOF;
    std::unordered_set<int> matchedTRD;
    for (uint8_t leg = 0; leg < 2; leg++) {
      auto refs = recoData.getSingleDetectorRefs(legs[leg]);
      if (refs[GTrackID::TPC].isIndexSet()) {
        tpcLegs[leg] = &recoData.getTPCTrack(refs[GTrackID::TPC]);
        (leg == 0 ? out.tpcBottom : out.tpcTop) = *tpcLegs[leg];
      }
      if (refs[GTrackID::ITS].isIndexSet()) {
        addITS(recoData, refs[GTrackID::ITS], leg, out.clITS, ic, pattRequests, matchedITS);
      }
      if (refs[GTrackID::TOF].isIndexSet()) {
        addTOF(recoData, refs[GTrackID::TOF], leg, out.clTOF, matchedTOF);
        if (mDebugOut) {
          writeDebugTOF(recoData.getTOFClusters()[refs[GTrackID::TOF].getIndex()], ic, leg, o2::dataformats::HitMatched);
        }
      }
      if (refs[GTrackID::TRD].isIndexSet()) {
        addTRD(recoData, refs[GTrackID::TRD], leg, out.trdTracklets, matchedTRD);
      }
    }
    CosmicTime cosmicTime;
    cosmicTime.tb = cosm.getTimeMUS().getTimeStamp() / mTPCTBinMUS;
    cosmicTime.errTB = cosm.getTimeMUS().getTimeStampError() / mTPCTBinMUS;
    cosmicTime.known = cosm.getTimeMUS().getTimeStampError() < mMaxAbsTimeErr;
    mDbgCosmic = ic;
    // TOF road first: a top/bottom hit pair matching the muon's flight gives the cosmic's time to ~ns (also for one-side legs on the same
    // side, whose brackets leave it open by tens of mus); the TPC corridor of the other side and the TRD / ITS roads then use that time
    if (mRoadDets[DetID::TOF]) {
      roadTOF(recoData, tpcLegs, cosmicTime, ic, out.clTOF, matchedTOF, out.timeTOFMUS);
    }
    CosmicTime roadTime = cosmicTime;
    if (out.timeTOFMUS >= 0.f) {
      roadTime.tb = out.timeTOFMUS / mTPCTBinMUS;
      roadTime.errTB = mTOFTimeErr / mTPCTBinMUS;
      roadTime.known = true;
    }
    mDbgTauC = roadTime.tb;
    std::unordered_set<uint32_t> taken;
    for (int leg = 0; leg < 2; leg++) { // attached clusters of both legs first, so that a road never takes the other leg's clusters
      if (tpcLegs[leg]) {
        addTPCAttached(recoData, *tpcLegs[leg], *tpcCl[leg], taken);
        mNClAttached += tpcCl[leg]->size();
      }
    }
    for (int leg = 0; leg < 2; leg++) {
      if (tpcLegs[leg]) {
        const size_t nAttached = tpcCl[leg]->size();
        mDbgLeg = leg;
        addTPCCorridor(recoData, *tpcLegs[leg], roadTime, *tpcCl[leg], taken);
        mNClCorridor += tpcCl[leg]->size() - nAttached;
      }
    }
    if (mRoadDets[DetID::TRD]) {
      roadTRD(recoData, tpcLegs, roadTime, ic, out.trdTracklets, matchedTRD);
    }
    if (mRoadDets[DetID::ITS]) {
      const int side0 = tpcLegs[0] ? (tpcLegs[0]->hasASideClustersOnly() ? 1 : (tpcLegs[0]->hasCSideClustersOnly() ? -1 : 0)) : 0;
      const int side1 = tpcLegs[1] ? (tpcLegs[1]->hasASideClustersOnly() ? 1 : (tpcLegs[1]->hasCSideClustersOnly() ? -1 : 0)) : 0;
      roadITS(recoData, cosm, roadTime, side0 == side1 ? side0 : 0, ic, out.clITS, pattRequests, matchedITS);
    }
  }
  if (!pattRequests.empty()) {
    fillITSPatterns(recoData, pattRequests, cosmicsOut);
  }
  if (mDebugOut) {
    for (size_t ic = 0; ic < cosmicsOut.size(); ic++) {
      writeDebug(cosmicsOut[ic], ic);
    }
  }
  mNCosmics += cosmicsOut.size();
  LOGP(info, "Collected clusters for {} cosmics in TF {}", cosmicsOut.size(), tfID.tfCounter);
  pc.outputs().snapshot(Output{"GLO", "COSMFULL", 0}, cosmicsOut);
  pc.outputs().snapshot(Output{"GLO", "COSMFULLTF", 0}, tfInfo);
  pc.outputs().snapshot(Output{"GLO", "COSMFULLTFID", 0}, tfID);
  mTimer.Stop();
}

void CosmicsClusterCollectorSpec::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  mCorrMap = &o2::gpu::TPCFastTransformPOD::get(pc.inputs().get<const char*>("corrMap"));
  mTPCTBinMUS = o2::tpc::ParameterElectronics::Instance().ZbinWidth;
  mBz = o2::base::Propagator::Instance()->getNominalBz();
  if (mRoadDets[DetID::ITS] && mITSChipCentres.empty()) {
    cacheITSChipCentres();
  }
}

void CosmicsClusterCollectorSpec::cacheITSChipCentres()
{
  auto geom = o2::its::GeometryTGeo::Instance();
  geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  mITSChipCentres.resize(geom->getNumberOfChips());
  for (int chip = 0; chip < geom->getNumberOfChips(); chip++) {
    mITSChipCentres[chip] = geom->getMatrixL2G(chip) * o2::math_utils::Point3D<float>(0.f, 0.f, 0.f);
  }
}

void CosmicsClusterCollectorSpec::buildUsedMap(const RecoContainer& data, std::vector<std::pair<float, float>>& time0Windows)
{
  // merge the windows, then flag the clusters of the TPC tracks whose time0 falls into one of them
  std::sort(time0Windows.begin(), time0Windows.end());
  std::vector<std::pair<float, float>> merged;
  for (const auto& window : time0Windows) {
    if (!merged.empty() && window.first <= merged.back().second) {
      merged.back().second = std::max(merged.back().second, window.second);
    } else {
      merged.push_back(window);
    }
  }
  const auto& clusters = data.getTPCClusters();
  const auto tracks = data.getTPCTracks();
  const auto refs = data.getTPCTracksClusterRefs();
  mUsed.assign(clusters.nClustersTotal, false);
  for (const auto& trk : tracks) {
    const float time0 = trk.getTime0();
    const auto window = std::upper_bound(merged.begin(), merged.end(), time0, [](float t, const std::pair<float, float>& w) { return t < w.first; });
    if (window == merged.begin() || time0 > std::prev(window)->second) {
      continue;
    }
    for (int j = 0; j < trk.getNClusterReferences(); j++) {
      uint8_t sector = 0;
      uint8_t row = 0;
      uint32_t clIdx = 0;
      trk.getClusterReference(refs, j, sector, row, clIdx);
      mUsed[clusters.clusterOffset[sector][row] + clIdx] = true;
    }
  }
}

void CosmicsClusterCollectorSpec::addTPCAttached(const RecoContainer& data, const o2::tpc::TrackTPC& trk, std::vector<o2::dataformats::CosmicTPCCluster>& out, std::unordered_set<uint32_t>& taken) const
{
  const auto& clusters = data.getTPCClusters();
  const auto refs = data.getTPCTracksClusterRefs();
  for (int j = 0; j < trk.getNClusterReferences(); j++) {
    uint8_t sector = 0;
    uint8_t row = 0;
    uint32_t clIdx = 0;
    trk.getClusterReference(refs, j, sector, row, clIdx);
    if (!taken.insert(clusters.clusterOffset[sector][row] + clIdx).second) {
      continue;
    }
    auto& cl = out.emplace_back();
    cl.cl = clusters.clusters[sector][row][clIdx];
    cl.sector = sector;
    cl.row = row;
    cl.flags = o2::dataformats::CosmicTPCCluster::Attached | o2::dataformats::CosmicTPCCluster::Used;
  }
}

void CosmicsClusterCollectorSpec::LegBranch::init(const o2::track::TrackPar& inner, const o2::track::TrackPar& outer, float bz)
{
  if (std::abs(inner.getCurvature(bz)) < 1e-5f) { // radius > 1 km: straight line
    isLine = true;
    const auto point = inner.getXYZGlo();
    const float phi = inner.getPhi();
    dirX = std::cos(phi);
    dirY = std::sin(phi);
    const float proj = point.X() * dirX + point.Y() * dirY;
    pcaX = point.X() - proj * dirX;
    pcaY = point.Y() - proj * dirY;
  } else {
    o2::math_utils::CircleXYf_t circle;
    float sinAlpha = 0.f;
    float cosAlpha = 0.f;
    inner.getCircleParams(bz, circle, sinAlpha, cosAlpha);
    centerX = circle.xC;
    centerY = circle.yC;
    const float centerDist = std::sqrt(centerX * centerX + centerY * centerY);
    if (centerDist < 1e-3f) { // circle around the beam line: no closest approach
      sign = 0;
      return;
    }
    pcaX = centerX * (1.f - circle.rC / centerDist);
    pcaY = centerY * (1.f - circle.rC / centerDist);
  }
  // side of the leg from its end farther from the closest approach; a leg with ends on both sides (both > 10 cm away) spans the
  // closest approach and is accepted on both sides
  const auto pointIn = inner.getXYZGlo();
  const auto pointOut = outer.getXYZGlo();
  const float distIn2 = (pointIn.X() - pcaX) * (pointIn.X() - pcaX) + (pointIn.Y() - pcaY) * (pointIn.Y() - pcaY);
  const float distOut2 = (pointOut.X() - pcaX) * (pointOut.X() - pcaX) + (pointOut.Y() - pcaY) * (pointOut.Y() - pcaY);
  const int sideIn = side(pointIn.X(), pointIn.Y());
  const int sideOut = side(pointOut.X(), pointOut.Y());
  constexpr float MinDist2 = 10.f * 10.f;
  if (sideIn == sideOut) {
    sign = sideIn;
  } else if (std::min(distIn2, distOut2) < MinDist2) {
    sign = distIn2 > distOut2 ? sideIn : sideOut;
  } else {
    sign = 0;
  }
}

void CosmicsClusterCollectorSpec::addTPCCorridor(const RecoContainer& data, const o2::tpc::TrackTPC& trk, const CosmicTime& cosmicTime, std::vector<o2::dataformats::CosmicTPCCluster>& out, std::unordered_set<uint32_t>& taken) const
{
  constexpr int NSectorsA = TPCGeo::getNumberOfSectorsA();
  constexpr float TanSector = 0.17632698f; // tan(10 deg)
  const float zLength = TPCGeo::getTPCzLength();
  const float vDrift = mCorrMap->getVDrift();
  const auto& clusters = data.getTPCClusters();
  const float time0Leg = trk.getTime0();
  const int legSide = trk.hasASideClustersOnly() ? 1 : (trk.hasCSideClustersOnly() ? -1 : 0);

  struct Frame {
    float vertexTime; ///< vertex time used for the transformation [TB]
    float dz;         ///< shift of the leg's z into this frame
    float zTolerance; ///< extra z tolerance from the error of the vertex time
    int sectorMin;
    int sectorMax;
    uint8_t flag;
  };
  Frame frames[2];
  int nFrames = 0;
  if (legSide == 0) { // CE-crossing leg: its time0 is absolute
    frames[nFrames++] = {time0Leg, 0.f, 0.f, 0, 2 * NSectorsA, 0};
  } else {
    frames[nFrames++] = {time0Leg, 0.f, 0.f, legSide > 0 ? 0 : NSectorsA, legSide > 0 ? NSectorsA : 2 * NSectorsA, 0};
    if (cosmicTime.known) { // z of a one-side leg moves by side * vD * (t - time0) when its clusters are transformed with the time t
      frames[nFrames++] = {cosmicTime.tb, legSide * (cosmicTime.tb - time0Leg) * vDrift, cosmicTime.errTB * vDrift,
                           legSide > 0 ? NSectorsA : 0, legSide > 0 ? 2 * NSectorsA : NSectorsA, o2::dataformats::CosmicTPCCluster::AbsTime};
    }
  }

  const o2::track::TrackPar refPar[2] = {trk, trk.getParamOut()}; // rows below rMid: inner parameter, above: outer
  LegBranch branch;
  branch.init(refPar[0], refPar[1], mBz);
  const auto pointIn = refPar[0].getXYZGlo();
  const auto pointOut = refPar[1].getXYZGlo();
  const float rMid = 0.5f * (std::hypot(pointIn.X(), pointIn.Y()) + std::hypot(pointOut.X(), pointOut.Y()));
  auto outsideDriftVolume = [zLength](bool sideA, float z, float margin) {
    return sideA ? (z < -margin || z > zLength + margin) : (z > margin || z < -zLength - margin);
  };

  for (int iFrame = 0; iFrame < nFrames; iFrame++) {
    const auto& frame = frames[iFrame];
    for (int sector = frame.sectorMin; sector < frame.sectorMax; sector++) {
      const float alpha = o2::math_utils::sector2Angle(sector % NSectorsA);
      const float sinAlpha = std::sin(alpha);
      const float cosAlpha = std::cos(alpha);
      const bool sideA = sector < NSectorsA;
      for (int iPar = 0; iPar < 2; iPar++) {
        auto par = refPar[iPar];
        par.setZ(par.getZ() + frame.dz);
        if (!par.rotateParam(alpha)) { // the leg points away from this sector frame: same helix, opposite direction
          par.invertParam();
          if (!par.rotateParam(alpha)) {
            continue;
          }
        }
        for (int row = 0; row < o2::tpc::constants::MAXGLOBALPADROW; row++) {
          const float xRow = TPCGeo::getRowInfoX(row);
          if ((xRow < rMid) != (iPar == 0)) {
            continue;
          }
          auto parRow = par;
          if (!parRow.propagateParamTo(xRow, mBz)) {
            continue;
          }
          // coarse acceptance at the nominal x of the row before the more expensive move to its real x
          constexpr float CoarseMargin = 10.f; // [cm] covers the shift of y and z between the nominal and the real x
          if (std::abs(parRow.getY()) > xRow * TanSector + mCorridor + CoarseMargin || outsideDriftVolume(sideA, parRow.getZ(), mCorridor + frame.zTolerance + CoarseMargin)) {
            continue;
          }
          // the corrected clusters of this row lie at the real x of the row, x + dx(y, z), not at its nominal x
          float xReal = xRow;
          const float zInside = sideA ? std::clamp(parRow.getZ(), 0.f, zLength) : std::clamp(parRow.getZ(), -zLength, 0.f);
          mCorrMap->InverseTransformYZtoX(sector, row, parRow.getY(), zInside, xReal);
          if (!parRow.propagateParamTo(xReal, mBz)) {
            continue;
          }
          const float y = parRow.getY();
          const float z = parRow.getZ();
          if (std::abs(y) > xRow * TanSector + mCorridor) {
            continue;
          }
          if (outsideDriftVolume(sideA, z, mCorridor + frame.zTolerance)) {
            continue;
          }
          if (!branch.accept(xRow * cosAlpha - y * sinAlpha, xRow * sinAlpha + y * cosAlpha)) {
            continue;
          }
          if (mDebugOut) { // predicted point, z in the frame of the road and in the common frame of the cosmic's time
            const float zCos = z + (sideA ? 1.f : -1.f) * (mDbgTauC - frame.vertexTime) * vDrift;
            (*mDebugOut) << "road"
                         << "tf=" << mDbgTF << "cosm=" << mDbgCosmic << "leg=" << mDbgLeg << "frame=" << int(frame.flag != 0) << "sector=" << sector
                         << "row=" << row << "x=" << xReal << "y=" << y << "z=" << z << "zCos=" << zCos << "gx=" << xReal * cosAlpha - y * sinAlpha
                         << "gy=" << xReal * sinAlpha + y * cosAlpha << "snp=" << parRow.getSnp() << "tgl=" << parRow.getTgl() << "\n";
          }
          searchRow(clusters, sector, row, y, z, parRow.getSnp(), parRow.getTgl(), frame.vertexTime, frame.zTolerance, frame.flag, out, taken);
        }
      }
    }
  }
}

void CosmicsClusterCollectorSpec::searchRow(const o2::tpc::ClusterNativeAccess& clusters, int sector, int row, float y, float z, float snp, float tgl, float vertexTime, float zTolerance, uint8_t flag,
                                            std::vector<o2::dataformats::CosmicTPCCluster>& out, std::unordered_set<uint32_t>& taken) const
{
  const float zLength = TPCGeo::getTPCzLength();
  const float vDrift = mCorrMap->getVDrift();
  const float t0 = mCorrMap->getT0();
  // nominal (measured) coordinates of the predicted real point; the correction is evaluated inside the drift volume
  const float zInside = sector < TPCGeo::getNumberOfSectorsA() ? std::clamp(z, 0.f, zLength) : std::clamp(z, -zLength, 0.f);
  float yNominal = 0.f;
  float zNominal = 0.f;
  mCorrMap->InverseTransformYZtoNominalYZ(sector, row, y, zInside, yNominal, zNominal);
  zNominal += z - zInside;
  float padPred = 0.f;
  float driftLengthPred = 0.f;
  TPCGeo::convLocalToPadDriftLength(sector, row, yNominal, zNominal, padPred, driftLengthPred);
  const float timePred = driftLengthPred / vDrift + t0 + vertexTime;
  // the road is a cylinder around the track; its section with the pad-row plane is an ellipse with half-axes W / cos(phi) in y and
  // W * sqrt(cos^2(phi) + tgl^2) / cos(phi) in z
  const float cosPhi = std::sqrt((1.f - snp) * (1.f + snp));
  const float cosPhiWindow = std::max(cosPhi, 0.1f); // limits the window for tracks nearly parallel to the pad row
  const float norm = 1.f / std::sqrt(1.f + tgl * tgl);
  const float dirY = snp * norm; // track direction in the pad-row plane (y, z components)
  const float dirZ = tgl * norm;
  const float corridor2 = mCorridor * mCorridor;
  const float windowPad = mCorridor / (cosPhiWindow * TPCGeo::getRowInfoPadWidth(row));
  const float windowTime = (mCorridor * std::sqrt(cosPhiWindow * cosPhiWindow + tgl * tgl) / cosPhiWindow + zTolerance) / vDrift;
  const auto* rowClusters = clusters.clusters[sector][row];
  const uint32_t rowOffset = clusters.clusterOffset[sector][row];
  for (uint32_t k = 0; k < clusters.nClusters[sector][row]; k++) {
    const auto& c = rowClusters[k];
    const float time = c.getTime();
    const float pad = c.getPad();
    if (std::abs(time - timePred) > windowTime || std::abs(pad - padPred) > windowPad) {
      continue;
    }
    float yCluster = 0.f;
    float zCluster = 0.f;
    TPCGeo::convPadDriftLengthToLocal(sector, row, pad, (time - t0 - vertexTime) * vDrift, yCluster, zCluster);
    const float dy = yCluster - yNominal;
    float dz = zCluster - zNominal;
    if (zTolerance > 0.f) { // the vertex time of this frame is known within zTolerance / vD: allow that shift along z
      dz = std::abs(dz) > zTolerance ? dz - std::copysign(zTolerance, dz) : 0.f;
    }
    const float proj = dy * dirY + dz * dirZ;
    if (dy * dy + dz * dz - proj * proj > corridor2) { // distance perpendicular to the track
      continue;
    }
    if (!taken.insert(rowOffset + k).second) {
      continue;
    }
    auto& cl = out.emplace_back();
    cl.cl = c;
    cl.sector = sector;
    cl.row = row;
    cl.flags = o2::dataformats::CosmicTPCCluster::Corridor | flag | (mUsed[rowOffset + k] ? o2::dataformats::CosmicTPCCluster::Used : 0);
  }
}

void CosmicsClusterCollectorSpec::addITS(const RecoContainer& data, GTrackID gid, uint8_t leg, std::vector<o2::dataformats::CosmicITSCluster>& out, int icosm, std::vector<ITSPattRequest>& requests, std::unordered_set<int>& matched) const
{
  if (gid.getSource() != GTrackID::ITS) { // ITS-AB tracklets are not stored
    return;
  }
  if (data.getITSPerLayer()) {
    if (!mStaggeredWarned) {
      LOGP(warning, "Staggered ITS clusters are not supported, ITS clusters of cosmics are not stored");
      mStaggeredWarned = true;
    }
    return;
  }
  const auto& trk = data.getITSTrack(gid);
  const auto refs = data.getITSTracksClusterRefs();
  const auto clusters = data.getITSClusters();
  const auto rofs = data.getITSClustersROFRecords();
  for (int i = 0; i < trk.getNumberOfClusters(); i++) {
    const int idx = refs[trk.getFirstClusterEntry() + i];
    const auto& c = clusters[idx];
    auto& cl = out.emplace_back();
    cl.chipID = c.getSensorID();
    cl.row = c.getRow();
    cl.col = c.getCol();
    cl.pattID = c.getPatternID();
    cl.leg = leg;
    cl.flags = o2::dataformats::HitMatched;
    matched.insert(idx);
    if (c.getPatternID() == o2::itsmft::CompCluster::InvalidPatternID || (mITSDict && mITSDict->isGroup(c.getPatternID()))) {
      requests.push_back({icosm, int(out.size()) - 1, idx}); // the pattern is in the TF's pattern stream, copied by fillITSPatterns
    }
    auto rof = std::upper_bound(rofs.begin(), rofs.end(), idx, [](int v, const o2::itsmft::ROFRecord& r) { return v < r.getFirstEntry(); });
    if (rof != rofs.begin()) {
      cl.rofBC = std::prev(rof)->getBCData().differenceInBC(mTFStart);
    }
  }
}

void CosmicsClusterCollectorSpec::addTOF(const RecoContainer& data, GTrackID gid, uint8_t leg, std::vector<o2::dataformats::CosmicTOFCluster>& out, std::unordered_set<int>& matched) const
{
  const auto& c = data.getTOFClusters()[gid.getIndex()];
  auto& cl = out.emplace_back();
  cl.timeRaw = c.getTimeRaw();
  cl.tot = c.getTot();
  cl.channel = c.getMainContributingChannel();
  cl.leg = leg;
  cl.flags = o2::dataformats::HitMatched;
  matched.insert(gid.getIndex());
}

std::pair<float, float> CosmicsClusterCollectorSpec::timeWindowMUS(const CosmicTime& cosmicTime) const
{
  // a time fixed by z continuity is given with a 1 sigma error, otherwise the error is the half-width of the legs' time-bracket overlap
  const float timeMUS = cosmicTime.tb * mTPCTBinMUS;
  const float errMUS = cosmicTime.errTB * mTPCTBinMUS;
  const float halfWidth = cosmicTime.known ? 5.f * errMUS + 0.2f : errMUS;
  return {timeMUS - halfWidth, timeMUS + halfWidth};
}

bool CosmicsClusterCollectorSpec::predictOutward(const o2::tpc::TrackTPC& leg, const CosmicTime& cosmicTime, int sector, float x, float& y, float& z) const
{
  // outward continuation of a leg in the frame of a sector; z in the frame of the cosmic's time (a one-side TPC-only leg has z relative to
  // its time0)
  o2::track::TrackPar par = leg.getParamOut();
  if (!par.rotateParam(o2::math_utils::sector2Angle(sector % TPCGeo::getNumberOfSectorsA())) || !par.propagateParamTo(x, mBz)) {
    return false;
  }
  const int side = leg.hasASideClustersOnly() ? 1 : (leg.hasCSideClustersOnly() ? -1 : 0);
  y = par.getY();
  z = par.getZ() + side * (cosmicTime.tb - leg.getTime0()) * mCorrMap->getVDrift();
  return true;
}

void CosmicsClusterCollectorSpec::roadTOF(const RecoContainer& data, const o2::tpc::TrackTPC* const* legs, const CosmicTime& cosmicTime, int icosm, std::vector<o2::dataformats::CosmicTOFCluster>& out, const std::unordered_set<int>& matched, float& timeTOFMUS) const
{
  // TOF clusters along the legs' outward continuations (the road in PbPb also contains hits of collision tracks). Preferred: the top/bottom
  // pair whose time difference matches the muon's flight between them; its mean time fixes the cosmic's time, also when the legs' brackets
  // leave it open by tens of mus (one-side legs on the same side), and the z of one-side legs is shifted to it. Otherwise per leg the
  // cluster closest to its continuation.
  constexpr float MaxFlightMUS = 0.1f; // flight time of the muon between the TPC and the TOF, slow tails
  constexpr float CmPerNS = 29.9792458f;
  const auto window = timeWindowMUS(cosmicTime);
  const float zTimeTol = 0.5f * (window.second - window.first) / mTPCTBinMUS * mCorrMap->getVDrift(); // z uncertainty of one-side legs
  const float vDriftPerMUS = mCorrMap->getVDrift() / mTPCTBinMUS;
  const float cosmicTimeMUS = cosmicTime.tb * mTPCTBinMUS;
  struct Candidate {
    int index;
    double timeNS; // since the start of the TF
    float dy;
    float dzAtCosmicTime;
    int side; // TPC side of a one-side leg (z moves with the time), 0: z absolute
    float gx, gy, gz;
  };
  std::vector<Candidate> candidates[2];
  const auto clusters = data.getTOFClusters();
  int best[2] = {-1, -1};
  float bestScore[2] = {1.f, 1.f};
  for (int i = 0; i < (int)clusters.size(); i++) {
    const auto& c = clusters[i];
    const float timeMUS = c.getTime() * 1e-6f; // [ps] since the start of the TF
    if (timeMUS < window.first - MaxFlightMUS || timeMUS > window.second + MaxFlightMUS || matched.count(i)) {
      continue;
    }
    for (int leg = 0; leg < 2; leg++) {
      float y = 0.f;
      float z = 0.f;
      if (!legs[leg] || !predictOutward(*legs[leg], cosmicTime, c.getSector(), c.getX(), y, z)) {
        continue;
      }
      const float normY = (c.getY() - y) / mRoadTOF;
      const float normZ = (c.getZ() - z) / (mRoadTOF + zTimeTol);
      const float score = std::max(normY * normY, normZ * normZ);
      if (score >= 1.f) {
        continue;
      }
      if (score < bestScore[leg]) {
        bestScore[leg] = score;
        best[leg] = i;
      }
      const float alpha = o2::math_utils::sector2Angle(c.getSector());
      const int side = legs[leg]->hasASideClustersOnly() ? 1 : (legs[leg]->hasCSideClustersOnly() ? -1 : 0);
      candidates[leg].push_back(Candidate{i, c.getTime() * 1e-3, c.getY() - y, c.getZ() - z, side, c.getX() * std::cos(alpha) - c.getY() * std::sin(alpha),
                                          c.getX() * std::sin(alpha) + c.getY() * std::cos(alpha), c.getZ()});
    }
  }
  float bestPairScore = -1.f;
  for (const auto& c0 : candidates[0]) {
    for (const auto& c1 : candidates[1]) {
      if (c0.index == c1.index) {
        continue;
      }
      const auto& top = c0.gy > c1.gy ? c0 : c1;
      const auto& bottom = c0.gy > c1.gy ? c1 : c0;
      const float length = std::sqrt((top.gx - bottom.gx) * (top.gx - bottom.gx) + (top.gy - bottom.gy) * (top.gy - bottom.gy) + (top.gz - bottom.gz) * (top.gz - bottom.gz));
      const float flightDev = float(top.timeNS - bottom.timeNS) + length / CmPerNS; // the muon crosses the top TOF first
      if (std::abs(flightDev) > mTOFFlightTol) {
        continue;
      }
      const double pairTimeNS = 0.5 * (c0.timeNS + c1.timeNS);
      const float shiftMUS = float(pairTimeNS * 1e-3) - cosmicTimeMUS;
      const float dz0 = c0.dzAtCosmicTime - c0.side * shiftMUS * vDriftPerMUS;
      const float dz1 = c1.dzAtCosmicTime - c1.side * shiftMUS * vDriftPerMUS;
      if (std::abs(dz0) > mRoadTOF || std::abs(dz1) > mRoadTOF) {
        continue;
      }
      const float score = (c0.dy * c0.dy + c1.dy * c1.dy + dz0 * dz0 + dz1 * dz1) / (mRoadTOF * mRoadTOF) + flightDev * flightDev / (mTOFFlightTol * mTOFFlightTol);
      if (bestPairScore < 0.f || score < bestPairScore) {
        bestPairScore = score;
        best[0] = c0.index;
        best[1] = c1.index;
        timeTOFMUS = float(pairTimeNS * 1e-3);
      }
    }
  }
  const uint8_t flags = bestPairScore < 0.f ? o2::dataformats::HitRoad : (o2::dataformats::HitRoad | o2::dataformats::HitTOFFlight);
  for (int leg = 0; leg < 2; leg++) {
    if (best[leg] < 0) {
      continue;
    }
    const auto& c = clusters[best[leg]];
    auto& cl = out.emplace_back();
    cl.timeRaw = c.getTimeRaw();
    cl.tot = c.getTot();
    cl.channel = c.getMainContributingChannel();
    cl.leg = leg;
    cl.flags = flags;
    if (mDebugOut) {
      writeDebugTOF(c, icosm, leg, flags);
    }
  }
}

void CosmicsClusterCollectorSpec::roadTRD(const RecoContainer& data, const o2::tpc::TrackTPC* const* legs, const CosmicTime& cosmicTime, int icosm, std::vector<o2::dataformats::CosmicTRDTracklet>& out, const std::unordered_set<int>& matched) const
{
  // per leg and layer the TRD tracklet closest to the leg's outward continuation, from triggers whose readout window can contain the cosmic
  constexpr float ReadoutWindowMUS = 3.f; // a cosmic leaves tracklets only if it passes within the readout window after a trigger
  constexpr float PadLength = 10.f;       // [cm] longest TRD pads: the tracklet z is the pad-row centre
  constexpr int NLayers = o2::trd::constants::NLAYER;
  const auto window = timeWindowMUS(cosmicTime);
  const float zTimeTol = 0.5f * (window.second - window.first) / mTPCTBinMUS * mCorrMap->getVDrift();
  const auto tracklets = data.getTRDTracklets();
  const auto calibrated = data.getTRDCalibratedTracklets();
  if (calibrated.size() != tracklets.size()) {
    return;
  }
  struct Candidate {
    int tracklet = -1;
    int trigBC = 0;
    float score = 1.f;
    float dy = 0.f;
    float dz = 0.f;
  };
  Candidate best[2][NLayers];
  for (const auto& trig : data.getTRDTriggerRecords()) {
    const int trigBC = trig.getBCData().differenceInBC(mTFStart);
    const float trigMUS = trigBC * o2::constants::lhc::LHCBunchSpacingMUS;
    if (trigMUS < window.first - ReadoutWindowMUS || trigMUS > window.second) {
      continue;
    }
    for (int it = trig.getFirstTracklet(); it < trig.getFirstTracklet() + trig.getNumberOfTracklets(); it++) {
      if (matched.count(it)) {
        continue;
      }
      const int detector = tracklets[it].getDetector();
      const int sector = detector / o2::trd::constants::NCHAMBERPERSEC;
      const int layer = detector % NLayers;
      const auto& point = calibrated[it];
      for (int leg = 0; leg < 2; leg++) {
        float y = 0.f;
        float z = 0.f;
        if (!legs[leg] || !predictOutward(*legs[leg], cosmicTime, sector, point.getX(), y, z)) {
          continue;
        }
        const float normY = (point.getY() - y) / mRoadTRD;
        const float normZ = (point.getZ() - z) / (mRoadTRD + PadLength + zTimeTol);
        const float score = std::max(normY * normY, normZ * normZ);
        if (score < best[leg][layer].score) {
          best[leg][layer] = {it, trigBC, score, point.getY() - y, point.getZ() - z};
        }
      }
    }
  }
  for (int leg = 0; leg < 2; leg++) {
    for (int layer = 0; layer < NLayers; layer++) {
      const auto& cand = best[leg][layer];
      if (cand.tracklet < 0) {
        continue;
      }
      auto& tr = out.emplace_back();
      tr.word = tracklets[cand.tracklet].getTrackletWord();
      tr.trigBC = cand.trigBC;
      tr.layer = layer;
      tr.leg = leg;
      tr.flags = o2::dataformats::HitRoad;
      if (mDebugOut) {
        const auto& point = calibrated[cand.tracklet];
        (*mDebugOut) << "trd"
                     << "tf=" << mDbgTF << "cosm=" << icosm << "leg=" << leg << "layer=" << layer << "x=" << point.getX() << "y=" << point.getY()
                     << "z=" << point.getZ() << "dy=" << cand.dy << "dz=" << cand.dz << "trigMUS=" << cand.trigBC * float(o2::constants::lhc::LHCBunchSpacingMUS) << "\n";
      }
    }
  }
}

void CosmicsClusterCollectorSpec::roadITS(const RecoContainer& data, const o2::dataformats::TrackCosmics& cosm, const CosmicTime& cosmicTime, int legsSide, int icosm, std::vector<o2::dataformats::CosmicITSCluster>& out,
                                          std::vector<ITSPattRequest>& requests, const std::unordered_set<int>& matched) const
{
  if (data.getITSPerLayer() || mITSChipCentres.empty()) {
    return;
  }
  // trajectory near the beam line from the combined cosmic, sampled every cm along its direction (closest approach at local x = 0)
  constexpr float MaxRadius = 45.f;    // [cm] outer ITS barrel + margin
  constexpr float ChipHalfDiag = 1.7f; // [cm] half diagonal of an ITS chip
  o2::track::TrackPar par = cosm;
  if (!par.rotateParam(par.getPhi())) {
    return;
  }
  std::vector<o2::math_utils::Point3D<float>> points;
  for (float x = -MaxRadius - 5.f; x <= MaxRadius + 5.f; x += 1.f) {
    bool ok = false;
    const auto point = par.getXYZGloAt(x, mBz, ok);
    if (ok) {
      points.push_back(point);
    }
  }
  if (points.size() < 2) {
    return;
  }
  // closest segment in the transverse plane: distance, z of the trajectory there
  auto closest = [&points](float x, float y, float& dist2, float& zTraj) {
    dist2 = 1e10f;
    for (size_t i = 0; i + 1 < points.size(); i++) {
      const float segX = points[i + 1].X() - points[i].X();
      const float segY = points[i + 1].Y() - points[i].Y();
      const float len2 = segX * segX + segY * segY;
      const float frac = len2 > 0.f ? std::clamp(((x - points[i].X()) * segX + (y - points[i].Y()) * segY) / len2, 0.f, 1.f) : 0.f;
      const float dx = x - (points[i].X() + frac * segX);
      const float dy = y - (points[i].Y() + frac * segY);
      if (dx * dx + dy * dy < dist2) {
        dist2 = dx * dx + dy * dy;
        zTraj = points[i].Z() + frac * (points[i + 1].Z() - points[i].Z());
      }
    }
  };
  float minR2 = 1e10f;
  float pcaY = 0.f;
  for (const auto& point : points) {
    const float r2 = point.X() * point.X() + point.Y() * point.Y();
    if (r2 < minR2) {
      minR2 = r2;
      pcaY = point.Y();
    }
  }
  if (minR2 > MaxRadius * MaxRadius) { // the cosmic does not cross the ITS
    return;
  }
  std::vector<bool> candidateChip(mITSChipCentres.size(), false);
  bool anyChip = false;
  const float chipCut2 = (mRoadITS + ChipHalfDiag) * (mRoadITS + ChipHalfDiag);
  for (size_t chip = 0; chip < mITSChipCentres.size(); chip++) {
    float dist2 = 0.f;
    float zTraj = 0.f;
    closest(mITSChipCentres[chip].X(), mITSChipCentres[chip].Y(), dist2, zTraj);
    if (dist2 < chipCut2) {
      candidateChip[chip] = true;
      anyChip = true;
    }
  }
  if (!anyChip) {
    return;
  }
  const auto window = timeWindowMUS(cosmicTime);
  const float zTimeTol = 0.5f * (window.second - window.first) / mTPCTBinMUS * mCorrMap->getVDrift(); // z of the cosmic is in the frame of its time
  // the refitted cosmic has the z of its TPC time; with legs on one TPC side and a road time from the TOF its z moves by side * vD * dt
  const float zShift = legsSide * (cosmicTime.tb - cosm.getTimeMUS().getTimeStamp() / mTPCTBinMUS) * mCorrMap->getVDrift();
  const float rofLengthMUS = o2::itsmft::DPLAlpideParam<DetID::ITS>::Instance().roFrameLengthInBC * o2::constants::lhc::LHCBunchSpacingMUS;
  // per half of the cosmic and layer the ITS cluster closest to the trajectory (the road near the beam line also contains collision clusters)
  constexpr int NLayers = 7;
  struct Candidate {
    int index = -1;
    int rofBC = 0;
    float score = 1.f;
  };
  Candidate best[2][NLayers];
  const auto clusters = data.getITSClusters();
  auto geom = o2::its::GeometryTGeo::Instance();
  for (const auto& rof : data.getITSClustersROFRecords()) {
    const int rofBC = rof.getBCData().differenceInBC(mTFStart);
    const float rofMUS = rofBC * o2::constants::lhc::LHCBunchSpacingMUS;
    if (rofMUS + rofLengthMUS < window.first || rofMUS > window.second) {
      continue;
    }
    for (int idx = rof.getFirstEntry(); idx < rof.getFirstEntry() + rof.getNEntries(); idx++) {
      const auto& c = clusters[idx];
      if (!candidateChip[c.getSensorID()] || matched.count(idx)) {
        continue;
      }
      o2::math_utils::Point3D<float> local;
      if (mITSDict && c.getPatternID() != o2::itsmft::CompCluster::InvalidPatternID && !mITSDict->isGroup(c.getPatternID())) {
        local = mITSDict->getClusterCoordinates<float>(c);
      } else { // anchor pixel: good enough for a cm road
        o2::itsmft::SegmentationAlpide::detectorToLocalUnchecked(c.getRow(), c.getCol(), local);
      }
      const auto global = geom->getMatrixL2G(c.getSensorID()) * local;
      float dist2 = 0.f;
      float zTraj = 0.f;
      closest(global.X(), global.Y(), dist2, zTraj);
      const float normZ = (global.Z() - zTraj - zShift) / (mRoadITS + zTimeTol);
      const float score = std::max(dist2 / (mRoadITS * mRoadITS), normZ * normZ);
      const int half = global.Y() < pcaY ? 0 : 1; // bottom / top half of the cosmic
      const int layer = geom->getLayer(c.getSensorID());
      if (layer >= 0 && layer < NLayers && score < best[half][layer].score) {
        best[half][layer] = {idx, rofBC, score};
      }
    }
  }
  for (int half = 0; half < 2; half++) {
    for (int layer = 0; layer < NLayers; layer++) {
      const auto& cand = best[half][layer];
      if (cand.index < 0) {
        continue;
      }
      const auto& c = clusters[cand.index];
      auto& cl = out.emplace_back();
      cl.chipID = c.getSensorID();
      cl.row = c.getRow();
      cl.col = c.getCol();
      cl.pattID = c.getPatternID();
      cl.rofBC = cand.rofBC;
      cl.leg = half;
      cl.flags = o2::dataformats::HitRoad;
      if (c.getPatternID() == o2::itsmft::CompCluster::InvalidPatternID || (mITSDict && mITSDict->isGroup(c.getPatternID()))) {
        requests.push_back({icosm, int(out.size()) - 1, cand.index});
      }
    }
  }
}

void CosmicsClusterCollectorSpec::fillITSPatterns(const RecoContainer& data, std::vector<ITSPattRequest>& requests, std::vector<o2::dataformats::CosmicTrack>& cosmics) const
{
  // the TF's pattern stream holds, in cluster order, the patterns of the clusters with an invalid or a group pattern ID
  if (!mITSDict) {
    if (!mNoDictWarned) {
      LOGP(warning, "No ITS cluster dictionary: patterns of ITS clusters of cosmics are not stored");
      mNoDictWarned = true;
    }
    return;
  }
  std::sort(requests.begin(), requests.end(), [](const ITSPattRequest& a, const ITSPattRequest& b) { return a.index < b.index; });
  const auto clusters = data.getITSClusters();
  auto pattIt = data.getITSClustersPatterns().begin();
  size_t ir = 0;
  for (int k = 0; k < (int)clusters.size() && ir < requests.size(); k++) {
    const auto pattID = clusters[k].getPatternID();
    if (pattID != o2::itsmft::CompCluster::InvalidPatternID && !mITSDict->isGroup(pattID)) {
      continue;
    }
    const auto start = pattIt;
    o2::itsmft::ClusterPattern::skipPattern(pattIt);
    for (; ir < requests.size() && requests[ir].index == k; ir++) {
      auto& cosm = cosmics[requests[ir].cosmic];
      cosm.clITS[requests[ir].cluster].pattEntry = cosm.itsPatterns.size();
      cosm.itsPatterns.insert(cosm.itsPatterns.end(), start, pattIt);
    }
  }
}

void CosmicsClusterCollectorSpec::addTRD(const RecoContainer& data, GTrackID gid, uint8_t leg, std::vector<o2::dataformats::CosmicTRDTracklet>& out, std::unordered_set<int>& matched) const
{
  const auto& trk = data.getTrack<o2::trd::TrackTRD>(gid); // TPC-TRD or ITS-TPC-TRD track
  const auto tracklets = data.getTRDTracklets();
  const auto trigs = data.getTRDTriggerRecords();
  for (int layer = 0; layer < 6; layer++) {
    const int it = trk.getTrackletIndex(layer);
    if (it < 0) {
      continue;
    }
    auto& tr = out.emplace_back();
    tr.word = tracklets[it].getTrackletWord();
    tr.layer = layer;
    tr.leg = leg;
    tr.flags = o2::dataformats::HitMatched;
    matched.insert(it);
    auto trig = std::upper_bound(trigs.begin(), trigs.end(), it, [](int v, const o2::trd::TriggerRecord& t) { return v < t.getFirstTracklet(); });
    if (trig != trigs.begin()) {
      tr.trigBC = std::prev(trig)->getBCData().differenceInBC(mTFStart);
    }
  }
}

void CosmicsClusterCollectorSpec::writeDebug(const o2::dataformats::CosmicTrack& cosm, int icosm) const
{
  constexpr int NSectorsA = TPCGeo::getNumberOfSectorsA();
  const float timeCosmic = cosm.cosmic.getTimeMUS().getTimeStamp() / mTPCTBinMUS;
  const bool absTimeKnown = cosm.cosmic.getTimeMUS().getTimeStampError() < mMaxAbsTimeErr;
  int nAttached[2] = {0, 0};
  int nRoad[2] = {0, 0};
  int side[2] = {-2, -2};
  const o2::tpc::TrackTPC* legs[2] = {&cosm.tpcBottom, &cosm.tpcTop};
  const std::vector<o2::dataformats::CosmicTPCCluster>* legClusters[2] = {&cosm.clTPCBottom, &cosm.clTPCTop};
  for (int leg = 0; leg < 2; leg++) {
    if (legs[leg]->getNClusters() == 0) { // no TPC part
      continue;
    }
    side[leg] = legs[leg]->hasASideClustersOnly() ? 1 : (legs[leg]->hasCSideClustersOnly() ? -1 : 0);
    for (const auto& c : *legClusters[leg]) {
      // transform with the vertex time used by the road: the leg's time0, or the cosmic's time for clusters found on the other side
      const float vertexTime = c.isAbsTime() ? timeCosmic : legs[leg]->getTime0();
      float x = 0.f;
      float y = 0.f;
      float z = 0.f;
      mCorrMap->Transform(c.sector, c.row, c.cl.getPad(), c.cl.getTime(), x, y, z, vertexTime);
      float xCos = 0.f;
      float yCos = 0.f;
      float zCos = 0.f;
      mCorrMap->Transform(c.sector, c.row, c.cl.getPad(), c.cl.getTime(), xCos, yCos, zCos, timeCosmic);
      const float alpha = o2::math_utils::sector2Angle(c.sector % NSectorsA);
      const float sinAlpha = std::sin(alpha);
      const float cosAlpha = std::cos(alpha);
      (*mDebugOut) << "cl"
                   << "tf=" << mDbgTF << "cosm=" << icosm << "leg=" << leg << "flags=" << int(c.flags) << "sector=" << int(c.sector) << "row=" << int(c.row)
                   << "pad=" << c.cl.getPad() << "time=" << c.cl.getTime() << "qMax=" << float(c.cl.getQmax()) << "qTot=" << float(c.cl.getQtot())
                   << "x=" << x << "y=" << y << "z=" << z << "zCos=" << zCos << "gx=" << x * cosAlpha - y * sinAlpha << "gy=" << x * sinAlpha + y * cosAlpha << "\n";
      (c.isAttached() ? nAttached : nRoad)[leg]++;
    }
  }
  for (const auto& c : cosm.clITS) { // local coordinates on the chip from the dictionary or the stored pattern
    const o2::itsmft::CompClusterExt compCluster(c.row, c.col, c.pattID, c.chipID);
    o2::math_utils::Point3D<float> local;
    if (c.pattEntry >= 0) {
      o2::itsmft::ClusterPattern pattern;
      auto pattIt = cosm.itsPatterns.begin() + c.pattEntry;
      pattern.acquirePattern(pattIt);
      local = o2::itsmft::TopologyDictionary::getClusterCoordinates<float>(compCluster, pattern, c.pattID != o2::itsmft::CompCluster::InvalidPatternID);
    } else if (mITSDict && c.pattID != o2::itsmft::CompCluster::InvalidPatternID) {
      local = mITSDict->getClusterCoordinates<float>(compCluster);
    } else { // pattern not available: pixel centre
      o2::itsmft::SegmentationAlpide::detectorToLocalUnchecked(c.row, c.col, local);
    }
    o2::math_utils::Point3D<float> global(0.f, 0.f, 0.f);
    if (!mITSChipCentres.empty()) { // geometry loaded for the ITS road
      global = o2::its::GeometryTGeo::Instance()->getMatrixL2G(c.chipID) * local;
    }
    (*mDebugOut) << "its"
                 << "tf=" << mDbgTF << "cosm=" << icosm << "leg=" << int(c.leg) << "flags=" << int(c.flags) << "gx=" << global.X() << "gy=" << global.Y() << "gz=" << global.Z() << "chip=" << int(c.chipID) << "row=" << int(c.row) << "col=" << int(c.col)
                 << "pattID=" << int(c.pattID) << "hasPatt=" << int(c.pattEntry >= 0) << "xLoc=" << local.X() << "zLoc=" << local.Z() << "rofBC=" << c.rofBC << "\n";
  }
  const auto& time = cosm.cosmic.getTimeMUS();
  (*mDebugOut) << "cosm"
               << "tf=" << mDbgTF << "cosm=" << icosm << "t=" << time.getTimeStamp() << "tErr=" << time.getTimeStampError() << "absTime=" << int(absTimeKnown)
               << "chi2Match=" << cosm.cosmic.getChi2Match() << "chi2Refit=" << cosm.cosmic.getChi2Refit() << "q2pt=" << cosm.cosmic.getQ2Pt()
               << "tgl=" << cosm.cosmic.getTgl() << "side0=" << side[0] << "side1=" << side[1] << "nAtt0=" << nAttached[0] << "nAtt1=" << nAttached[1]
               << "nRoad0=" << nRoad[0] << "nRoad1=" << nRoad[1] << "nITS=" << int(cosm.clITS.size()) << "nTOF=" << int(cosm.clTOF.size())
               << "nTRD=" << int(cosm.trdTracklets.size()) << "tTOF=" << cosm.timeTOFMUS << "mcEvent=" << (cosm.label.isSet() ? cosm.label.getEventID() : -1)
               << "mcTrack=" << (cosm.label.isSet() ? cosm.label.getTrackID() : -1) << "\n";
}

void CosmicsClusterCollectorSpec::writeDebugTOF(const o2::tof::Cluster& c, int icosm, int leg, uint8_t flags) const
{
  // the TOF cluster position is in the frame of its sector
  const float alpha = o2::math_utils::sector2Angle(c.getSector());
  const float sinAlpha = std::sin(alpha);
  const float cosAlpha = std::cos(alpha);
  (*mDebugOut) << "tof"
               << "tf=" << mDbgTF << "cosm=" << icosm << "leg=" << leg << "channel=" << c.getMainContributingChannel() << "timeRaw=" << c.getTimeRaw()
               << "time=" << c.getTime() << "flags=" << int(flags) << "tot=" << c.getTot() << "x=" << c.getX() << "y=" << c.getY() << "z=" << c.getZ()
               << "gx=" << c.getX() * cosAlpha - c.getY() * sinAlpha << "gy=" << c.getX() * sinAlpha + c.getY() * cosAlpha << "\n";
}

void CosmicsClusterCollectorSpec::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "CLUSDICT", 0)) {
    mITSDict = (const o2::itsmft::TopologyDictionary*)obj;
    return;
  }
}

void CosmicsClusterCollectorSpec::endOfStream(EndOfStreamContext& ec)
{
  mDebugOut.reset();
  LOGP(info, "Cosmics cluster collector: {} cosmics, {} attached and {} road TPC clusters; Cpu: {:.3e} Real: {:.3e} s in {} slots",
       mNCosmics, mNClAttached, mNClCorridor, mTimer.CpuTime(), mTimer.RealTime(), mTimer.Counter() - 1);
}

DataProcessorSpec getCosmicsClusterCollectorSpec(GTrackID::mask_t src, bool useMC, bool itsStag, DetID::mask_t roadDets)
{
  auto dataRequest = std::make_shared<DataRequest>();
  dataRequest->setITSPerLayer(itsStag);
  dataRequest->requestTracks(src, false);
  dataRequest->requestClusters(src, false);
  if (roadDets[DetID::ITS]) {
    dataRequest->requestITSClusters(false);
  }
  if (roadDets[DetID::TOF]) {
    dataRequest->requestTOFClusters(false);
  }
  if (roadDets[DetID::TRD]) {
    dataRequest->requestTRDTracklets(false);
  }
  dataRequest->requestCoscmicTracks(useMC);
  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                                                                                     // orbitResetTime
                                                              false,                                                                                     // GRPECS
                                                              false,                                                                                     // GRPLHCIF
                                                              true,                                                                                      // GRPMagField
                                                              false,                                                                                     // askMatLUT
                                                              roadDets[DetID::ITS] ? o2::base::GRPGeomRequest::Aligned : o2::base::GRPGeomRequest::None, // ITS road: chip positions
                                                              dataRequest->inputs,
                                                              true);
  dataRequest->inputs.emplace_back("corrMap", o2::header::gDataOriginTPC, "TPCCORRMAP", 0, Lifetime::Timeframe);

  std::vector<OutputSpec> outputs;
  outputs.emplace_back("GLO", "COSMFULL", 0, Lifetime::Timeframe);
  outputs.emplace_back("GLO", "COSMFULLTF", 0, Lifetime::Timeframe);
  outputs.emplace_back("GLO", "COSMFULLTFID", 0, Lifetime::Timeframe);

  return DataProcessorSpec{
    "cosmics-cluster-collector",
    dataRequest->inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<CosmicsClusterCollectorSpec>(dataRequest, ggRequest, useMC, roadDets)},
    Options{
      {"corridor-width", VariantType::Float, 1.f, {"half-width of the road around each leg [cm]"}},
      {"max-abs-time-err", VariantType::Float, 0.5f, {"max. time error of a cosmic [mus] to search the other TPC side of its one-side legs"}},
      {"max-cosmics-per-tf", VariantType::Int, 100, {"collect the clusters of at most this many cosmics per TF"}},
      {"tof-road-width", VariantType::Float, 5.f, {"half-width of the road at the TOF [cm]"}},
      {"tof-flight-tolerance", VariantType::Float, 2.f, {"max. deviation of the top/bottom TOF time difference from the muon's flight time [ns]"}},
      {"tof-time-error", VariantType::Float, 0.1f, {"error of a cosmic's TOF time for the TPC other-side corridor and the TRD / ITS roads [mus]"}},
      {"trd-road-width", VariantType::Float, 5.f, {"half-width of the road at the TRD [cm] (in z plus the pad length)"}},
      {"its-road-width", VariantType::Float, 1.5f, {"half-width of the road in the ITS [cm]"}},
      {"debug-tree", VariantType::Bool, false, {"write cosmics_collector_debug.root with transformed clusters and road points (test runs)"}}}};
}

} // namespace o2::globaltracking
