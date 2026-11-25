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

#include <vector>
#include <TStopwatch.h>
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DataFormatsGlobalTracking/RecoContainerCreateTracksVariadic.h"
#include "DetectorsVertexing/SVertexerParams.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "SimulationDataFormat/MCEventLabel.h"
#include "SimulationDataFormat/MCUtils.h"
#include "SimulationDataFormat/MCTrack.h"
#include "CommonDataFormat/BunchFilling.h"
#include "CommonUtils/NameConf.h"
#include "DataFormatsFT0/RecPoints.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/CCDBParamSpec.h"
#include "FT0Reconstruction/InteractionTag.h"
#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "GlobalTrackingStudy/TrackingStudy.h"
#include "TPCBase/ParameterElectronics.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "DataFormatsFT0/RecPoints.h"
#include "DataFormatsITSMFT/TrkClusRef.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Steer/MCKinematicsReader.h"
#include "DCAFitter/DCAFitterN.h"
#include "MathUtils/fit.h"
#include "GlobalTrackingStudy/V0Ext.h"
#include "GlobalTrackingStudy/V0FullOutputParams.h"
#include "DataFormatsGlobalTracking/V0WithClusters.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsITSMFT/ClusterPattern.h"
#include "DataFormatsTOF/Cluster.h"
#include "DataFormatsTRD/TrackTRD.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include "CommonDataFormat/TFIDInfo.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DetectorsBase/TFIDInfoHelper.h"
#include "TPCFastTransformPOD.h"
#include "MathUtils/Utils.h"
#include "CommonDataFormat/InteractionRecord.h"
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <cmath>
#include "GPUO2InterfaceConfiguration.h"
// #include "GPUSettingsO2.h"
#include "GPUParam.h"
#include "GPUParam.inc"
#include "GPUTPCGeometry.h"
#include "GPUO2InterfaceRefit.h"
#include "GPUO2InterfaceUtils.h"

namespace o2::svstudy
{

using namespace o2::framework;
using DetID = o2::detectors::DetID;
using DataRequest = o2::globaltracking::DataRequest;

using PVertex = o2::dataformats::PrimaryVertex;
using V2TRef = o2::dataformats::VtxTrackRef;
using VTIndex = o2::dataformats::VtxTrackIndex;
using GTrackID = o2::dataformats::GlobalTrackID;
using TBracket = o2::math_utils::Bracketf_t;
using V0ID = o2::dataformats::V0Index;

using timeEst = o2::dataformats::TimeStampWithError<float, float>;

class SVStudySpec final : public Task
{
 public:
  SVStudySpec(std::shared_ptr<DataRequest> dr, std::shared_ptr<o2::base::GRPGeomRequest> gr, GTrackID::mask_t src, bool useTPCCl, bool useMC, bool fullOutput)
    : mDataRequest(dr), mGGCCDBRequest(gr), mTracksSrc(src), mUseTPCCl(useTPCCl), mUseMC(useMC), mFullOutput(fullOutput) {}
  ~SVStudySpec() final = default;
  void init(InitContext& ic) final;
  void run(ProcessingContext& pc) final;
  void endOfStream(EndOfStreamContext& ec) final;
  void finaliseCCDB(ConcreteDataMatcher& matcher, void* obj) final;
  void process(o2::globaltracking::RecoContainer& recoData);
  o2::dataformats::V0Ext processV0(int iv, o2::globaltracking::RecoContainer& recoData);

 private:
  /// ITS cluster of the full output whose pattern bytes are in the TF's pattern stream
  struct ITSPattRequest {
    int cluster; ///< entry of the cluster in the output ITS cluster table
    int index;   ///< index of the cluster in the TF
  };
  /// output of --full-output for one TF
  struct FullOutput {
    std::vector<o2::dataformats::V0WithClusters> v0s;
    std::vector<o2::dataformats::V0Prong> prongs;
    std::vector<o2::dataformats::V0TPCCluster> clTPC;
    std::vector<o2::dataformats::V0ITSCluster> clITS;
    std::vector<uint8_t> itsPatterns;
    std::vector<o2::dataformats::V0TOFCluster> clTOF;
    std::vector<o2::dataformats::V0TRDTracklet> trdTracklets;
    std::unordered_map<uint32_t, int> prongEntry; ///< prong table entry of a global track ID (GTrackID::getRawWOFlags)
    std::vector<ITSPattRequest> pattRequests;
  };

  void updateTimeDependentParams(ProcessingContext& pc);
  bool refitV0(const V0ID& id, o2::dataformats::V0& v0, o2::globaltracking::RecoContainer& recoData);
  const o2::dataformats::V0* getV0(int iv, o2::globaltracking::RecoContainer& recoData);
  void getTPCOccupancy(const PVertex& pv, float& occBef, float& occAft) const;
  void processFull(o2::globaltracking::RecoContainer& recoData, ProcessingContext& pc);
  bool passProngCuts(const V0ID& id, const o2::globaltracking::RecoContainer& recoData) const;
  int countCompatiblePVs(const V0ID& id, const o2::globaltracking::RecoContainer& recoData) const;
  int addProng(GTrackID gid, const o2::globaltracking::RecoContainer& recoData, FullOutput& out);
  void fillITSPatterns(const o2::globaltracking::RecoContainer& recoData, FullOutput& out) const;
  void writeFullDebug(const o2::dataformats::V0WithClusters& v0, const FullOutput& out, const PVertex& pv, int tfCounter) const;
  std::shared_ptr<DataRequest> mDataRequest;
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  bool mUseMC{false}; ///< MC flag
  std::unique_ptr<o2::utils::TreeStreamRedirector> mDBGOut;
  float mSelK0 = -1;
  bool mRefit = false;
  bool mUseTPCCl = false;
  bool mStudyTree = true;                                         ///< write svStudy.root
  bool mFullOutput = false;                                       ///< send the selected V0s with the raw clusters of their prongs
  std::vector<o2::dataformats::V0> mRefitV0s;                     ///< refitted V0s of the TF (index = V0 finder entry), shared by svStudy.root and the full output
  std::vector<int8_t> mRefitStatus;                               ///< per V0 finder entry: -1 not refitted yet, 0 refit failed, 1 refitted
  const o2::gpu::TPCFastTransformPOD* mCorrMap = nullptr;         ///< TPC transformation of the TF (full output)
  const o2::itsmft::TopologyDictionary* mITSDict = nullptr;       ///< ITS cluster dictionary (full output: which patterns to store)
  o2::InteractionRecord mTFStart{};                               ///< first BC of the TF
  std::unique_ptr<o2::utils::TreeStreamRedirector> mFullDebugOut; ///< --debug-tree: transformed clusters of the full-output V0s
  mutable bool mNoDictWarned = false;
  bool mMaxV0sWarned = false;
  size_t mNV0sFull = 0;    ///< V0s written with clusters
  size_t mNProngsFull = 0; ///< prongs written
  size_t mNClTPCFull = 0;  ///< TPC clusters written
  float mMaxEta = 0.8;
  float mBz = 0;
  int mNHBPerTF = 0;
  int mNTPCOccBinLength = 0; ///< TPC occ. histo bin length in TBs
  float mNTPCOccBinLengthInv = 0.f;
  float mTPCTBinMUSInv = 0.f;
  float mTPCTBinMUS = 0.f;
  GTrackID::mask_t mTracksSrc{};
  o2::vertexing::DCAFitterN<2> mFitterV0;
  std::vector<float> mTBinClOccAft, mTBinClOccBef;
  std::unique_ptr<o2::steer::MCKinematicsReader> mcReader; // reader of MC information
  std::shared_ptr<o2::gpu::GPUParam> mParam = nullptr;
};

void SVStudySpec::init(InitContext& ic)
{
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);
  mStudyTree = !ic.options().get<bool>("disable-study-tree");
  if (mStudyTree) {
    mDBGOut = std::make_unique<o2::utils::TreeStreamRedirector>("svStudy.root", "recreate");
  }
  if (mFullOutput && ic.options().get<bool>("debug-tree")) {
    mFullDebugOut = std::make_unique<o2::utils::TreeStreamRedirector>("v0s_full_debug.root", "recreate");
  }
  mRefit = ic.options().get<bool>("refit");
  mSelK0 = ic.options().get<float>("sel-k0");
  mMaxEta = ic.options().get<float>("max-eta");
  if (mUseMC) {
    mcReader = std::make_unique<o2::steer::MCKinematicsReader>("collisioncontext.root");
  }
}

void SVStudySpec::run(ProcessingContext& pc)
{
  o2::globaltracking::RecoContainer recoData;
  recoData.collectData(pc, *mDataRequest.get()); // select tracks of needed type, with minimal cuts, the real selected will be done in the vertexer
  updateTimeDependentParams(pc);                 // Make sure this is called after recoData.collectData, which may load some conditions

  size_t occupancyMapSizeBytes = o2::gpu::GPUO2InterfaceRefit::fillOccupancyMapGetSize(mNHBPerTF, mParam.get());
  gsl::span<const unsigned int> TPCRefitterOccMap = recoData.occupancyMapTPC;
  o2::gpu::GPUO2InterfaceUtils::paramUseExternalOccupancyMap(mParam.get(), mNHBPerTF, TPCRefitterOccMap.data(), occupancyMapSizeBytes);

  mTBinClOccBef.resize(1);
  mTBinClOccAft.resize(1);
  if (recoData.inputsTPCclusters && mUseTPCCl) {
    mNTPCOccBinLength = mParam->rec.tpc.occupancyMapTimeBins;
    mTBinClOccBef.clear();
    mTBinClOccAft.clear();
    // prepare TPC occupancy data
    if (mNTPCOccBinLength > 1 && recoData.occupancyMapTPC.size()) {
      mNTPCOccBinLengthInv = 1. / mNTPCOccBinLength;
      int nTPCBins = mNHBPerTF * o2::constants::lhc::LHCMaxBunches / 8, ninteg = 0;
      int nTPCOccBins = nTPCBins * mNTPCOccBinLengthInv, sumBins = std::max(1, int(o2::constants::lhc::LHCMaxBunches / 8 * mNTPCOccBinLengthInv));
      mTBinClOccAft.resize(nTPCOccBins);
      mTBinClOccBef.resize(nTPCOccBins);
      float sm = 0., tb = 0.5 * mNTPCOccBinLength;
      std::vector<float> mltHistTB(nTPCOccBins);
      for (int i = 0; i < nTPCOccBins; i++) {
        mltHistTB[i] = mParam->GetUnscaledMult(tb);
        tb += mNTPCOccBinLength;
      }
      for (int i = nTPCOccBins; i--;) {
        sm += mltHistTB[i];
        if (i + sumBins < nTPCOccBins) {
          sm -= mltHistTB[i + sumBins];
        }
        mTBinClOccAft[i] = sm;
      }
      sm = 0;
      for (int i = 0; i < nTPCOccBins; i++) {
        sm += mltHistTB[i];
        if (i - sumBins > 0) {
          sm -= mltHistTB[i - sumBins];
        }
        mTBinClOccBef[i] = sm;
      }
    }
  }

  // the V0 finder sends only the V0 indices unless svertexer.createFullV0s: then every V0 is refitted, once per TF for both outputs
  const auto nv0 = recoData.getV0sIdx().size();
  if (nv0 > recoData.getV0s().size()) {
    mRefit = true;
  }
  mRefitStatus.assign(nv0, -1);
  mRefitV0s.resize(nv0);

  if (mStudyTree) {
    process(recoData);
  }
  if (mFullOutput) {
    processFull(recoData, pc);
  }
}

void SVStudySpec::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  static bool initOnceDone = false;
  if (!initOnceDone) { // this params need to be queried only once
    initOnceDone = true;
    const auto& svparam = o2::vertexing::SVertexerParams::Instance();
    // Note: reading of the ITS AlpideParam needed for ITS timing is done by the RecoContainer
    mFitterV0.setBz(mBz);
    mFitterV0.setOldMode(svparam.oldDCAFitterMode);
    mFitterV0.setUseAbsDCA(svparam.useAbsDCA);
    mFitterV0.setPropagateToPCA(false);
    mFitterV0.setMaxR(svparam.maxRIni);
    mFitterV0.setMinParamChange(svparam.minParamChange);
    mFitterV0.setMinRelChi2Change(svparam.minRelChi2Change);
    mFitterV0.setMaxDZIni(svparam.maxDZIni);
    mFitterV0.setMaxDXYIni(svparam.maxDXYIni);
    mFitterV0.setMaxChi2(svparam.maxChi2);
    mFitterV0.setMatCorrType(o2::base::Propagator::MatCorrType(svparam.matCorr));
    mFitterV0.setUsePropagator(svparam.usePropagator);
    mFitterV0.setRefitWithMatCorr(svparam.refitWithMatCorr);
    mFitterV0.setMaxStep(svparam.maxStep);
    mFitterV0.setMaxSnp(svparam.maxSnp);
    mFitterV0.setMinXSeed(svparam.minXSeed);

    mNHBPerTF = o2::base::GRPGeomHelper::instance().getGRPECS()->getNHBFPerTF();
    if (!mParam) {
      // for occupancy estimator
      mParam = o2::gpu::GPUO2InterfaceUtils::getFullParamShared(0.f, mNHBPerTF);
    }
    auto& elParam = o2::tpc::ParameterElectronics::Instance();
    mTPCTBinMUS = elParam.ZbinWidth;
    mTPCTBinMUSInv = 1. / elParam.ZbinWidth; // 1./TPC bin in microseconds
  }
  mBz = o2::base::Propagator::Instance()->getNominalBz();
  mFitterV0.setBz(mBz);
  if (mFullOutput) {
    mCorrMap = &o2::gpu::TPCFastTransformPOD::get(pc.inputs().get<const char*>("corrMap"));
  }
}

o2::dataformats::V0Ext SVStudySpec::processV0(int iv, o2::globaltracking::RecoContainer& recoData)
{
  o2::dataformats::V0Ext v0ext;
  auto invalidate = [&v0ext]() { v0ext.v0ID.setVertexID(-1); };

  auto v0s = recoData.getV0s();
  auto v0IDs = recoData.getV0sIdx();
  static int tfID = 0;

  const auto& v0id = v0IDs[iv];
  if (mRefit) {
    const auto* v0refit = getV0(iv, recoData);
    if (!v0refit) {
      invalidate();
      return v0ext;
    }
    v0ext.v0 = *v0refit;
  }
  const auto& v0sel = mRefit ? v0ext.v0 : v0s[iv];
  if (mMaxEta < std::abs(v0sel.getEta())) {
    invalidate();
    return v0ext;
  }
  if (mSelK0 > 0 && std::abs(std::sqrt(v0sel.calcMass2AsK0()) - 0.497) > mSelK0) {
    invalidate();
    return v0ext;
  }
  if (!mRefit) {
    v0ext.v0 = v0sel;
  }
  v0ext.v0ID = v0id;
  const auto clRefs = recoData.getTPCTracksClusterRefs();
  o2::MCCompLabel lb[2];
  const o2::MCTrack* mcTrks[2];
  for (int ip = 0; ip < 2; ip++) {
    auto& prInfo = v0ext.prInfo[ip];
    auto gid = v0ext.v0ID.getProngID(ip);
    auto gidset = recoData.getSingleDetectorRefs(gid);
    lb[ip] = recoData.getTrackMCLabel(gid);
    if (lb[ip].isValid()) {
      prInfo.corrGlo = !lb[ip].isFake();
    }
    // get TPC tracks, if any
    if (gidset[GTrackID::TPC].isSourceSet()) {
      const auto& tpcTr = recoData.getTPCTrack(gidset[GTrackID::TPC]);
      prInfo.trackTPC = tpcTr;
      prInfo.nClTPC = tpcTr.getNClusters();
      lb[ip] = recoData.getTrackMCLabel(gidset[GTrackID::TPC]);
      if (lb[ip].isValid()) {
        prInfo.corrTPC = !lb[ip].isFake();
      }
      if (mParam && mUseTPCCl) {
        uint8_t clSect = 0, clRow = 0;
        uint32_t clIdx = 0;
        tpcTr.getClusterReference(clRefs, tpcTr.getNClusterReferences() - 1, clSect, clRow, clIdx);
        const auto& clus = recoData.getTPCClusters().clusters[clSect][clRow][clIdx];
        prInfo.lowestRow = clRow;
        int npads = o2::gpu::GPUTPCGeometry::NPads(clRow);
        prInfo.padFromEdge = uint8_t(clus.getPad());
        if (prInfo.padFromEdge > npads / 2) {
          prInfo.padFromEdge = npads - 1 - prInfo.padFromEdge;
        }
      }
    }
    // get ITS tracks, if any
    if (gid.includesDet(DetID::ITS)) {
      auto gidITS = recoData.getITSContributorGID(gid);
      if (gidset[GTrackID::ITS].isSourceSet()) {
        const auto& itsTr = recoData.getITSTrack(gidset[GTrackID::ITS]);
        prInfo.nClITS = itsTr.getNClusters();
        lb[ip] = recoData.getTrackMCLabel(gidset[GTrackID::ITS]);
        if (lb[ip].isValid()) {
          prInfo.corrITS = !lb[ip].isFake();
        }
        for (int il = 0; il < 7; il++) {
          if (itsTr.hasHitOnLayer(il)) {
            prInfo.pattITS |= 0x1 << il;
          }
        }
      } else {
        const auto& itsTrf = recoData.getITSABRefs()[gidset[GTrackID::ITSAB]];
        prInfo.nClITS = itsTrf.getNClusters();
        lb[ip] = recoData.getTrackMCLabel(gidset[GTrackID::ITSAB]);
        if (lb[ip].isValid()) {
          prInfo.corrITS = !lb[ip].isFake();
        }
        for (int il = 0; il < 7; il++) {
          if (itsTrf.hasHitOnLayer(il)) {
            prInfo.pattITS |= 0x1 << il;
          }
        }
        prInfo.pattITS |= 0x1 << 31; // flag AB
      }
      if (gidset[GTrackID::ITSTPC].isSourceSet()) {
        auto mtc = recoData.getTPCITSTrack(gidset[GTrackID::ITSTPC]);
        lb[ip] = recoData.getTrackMCLabel(gidset[GTrackID::ITSTPC]);
        prInfo.chi2ITSTPC = mtc.getChi2Match();
        if (lb[ip].isValid()) {
          prInfo.corrITSTPC = !lb[ip].isFake();
        }
      }
    }
    if (mUseMC && lb[ip].isValid()) { // temp store of mctrks
      mcTrks[ip] = mcReader->getTrack(lb[ip]);
    }
  }
  if (mUseMC && (mcTrks[0] != nullptr) && (mcTrks[1] != nullptr)) {
    // check majority vote on mother particle otherwise leave pdg -1
    if (lb[0].getSourceID() == lb[1].getSourceID() && lb[0].getEventID() == lb[1].getEventID() &&
        mcTrks[0]->getMotherTrackId() == mcTrks[1]->getMotherTrackId() && mcTrks[0]->getMotherTrackId() >= 0) {
      const auto mother = mcReader->getTrack(lb[0].getSourceID(), lb[0].getEventID(), mcTrks[0]->getMotherTrackId());
      v0ext.mcPID = mother->GetPdgCode();
    }
  }
  return v0ext;
}

void SVStudySpec::process(o2::globaltracking::RecoContainer& recoData)
{
  auto v0IDs = recoData.getV0sIdx();
  auto nv0 = v0IDs.size();
  std::map<int, std::vector<int>> pv2sv;
  static int tfID = 0;
  for (int iv = 0; iv < nv0; iv++) {
    const auto v0id = v0IDs[iv];
    pv2sv[v0id.getVertexID()].push_back(iv);
  }
  std::vector<o2::dataformats::V0Ext> v0extVec;
  for (auto it : pv2sv) {
    int pvID = it.first;
    auto& vv = it.second;
    if (pvID < 0 || vv.size() == 0) {
      continue;
    }
    v0extVec.clear();
    for (int iv0 : vv) {
      auto v0ext = processV0(iv0, recoData);
      if (v0ext.v0ID.getVertexID() < 0) {
        continue;
      }
      v0extVec.push_back(v0ext);
    }
    if (v0extVec.size()) {
      const auto& pv = recoData.getPrimaryVertex(pvID);
      float tpcOccBef = 0., tpcOccAft = 0.;
      getTPCOccupancy(pv, tpcOccBef, tpcOccAft);

      (*mDBGOut) << "v0"
                 << "orbit=" << recoData.startIR.orbit << "tfID=" << tfID << "tpcOccBef=" << tpcOccBef << "tpcOccAft=" << tpcOccAft
                 << "v0Ext=" << v0extVec
                 << "pv=" << pv
                 << "\n";
    }
  }
  tfID++;
}

bool SVStudySpec::refitV0(const V0ID& id, o2::dataformats::V0& v0, o2::globaltracking::RecoContainer& recoData)
{
  auto seedP = recoData.getTrackParam(id.getProngID(0));
  auto seedN = recoData.getTrackParam(id.getProngID(1));
  bool isTPConly = (id.getProngID(0).getSource() == GTrackID::TPC) || (id.getProngID(1).getSource() == GTrackID::TPC);
  const auto& svparam = o2::vertexing::SVertexerParams::Instance();
  if (svparam.mTPCTrackPhotonTune && isTPConly) {
    mFitterV0.setMaxDZIni(svparam.mTPCTrackMaxDZIni);
    mFitterV0.setMaxDXYIni(svparam.mTPCTrackMaxDXYIni);
    mFitterV0.setMaxChi2(svparam.mTPCTrackMaxChi2);
    mFitterV0.setCollinear(true);
  }
  int nCand = mFitterV0.process(seedP, seedN);
  if (svparam.mTPCTrackPhotonTune && isTPConly) { // restore
    // Reset immediately to the defaults
    mFitterV0.setMaxDZIni(svparam.maxDZIni);
    mFitterV0.setMaxDXYIni(svparam.maxDXYIni);
    mFitterV0.setMaxChi2(svparam.maxChi2);
    mFitterV0.setCollinear(false);
  }
  if (nCand == 0) { // discard this pair
    return false;
  }
  const int cand = 0;
  if (!mFitterV0.isPropagateTracksToVertexDone(cand) && !mFitterV0.propagateTracksToVertex(cand)) {
    return false;
  }
  const auto& trPProp = mFitterV0.getTrack(0, cand);
  const auto& trNProp = mFitterV0.getTrack(1, cand);
  std::array<float, 3> pP{}, pN{};
  trPProp.getPxPyPzGlo(pP);
  trNProp.getPxPyPzGlo(pN);
  std::array<float, 3> pV0 = {pP[0] + pN[0], pP[1] + pN[1], pP[2] + pN[2]};
  auto p2V0 = pV0[0] * pV0[0] + pV0[1] * pV0[1] + pV0[2] * pV0[2];
  const auto& pv = recoData.getPrimaryVertex(id.getVertexID());
  const auto v0XYZ = mFitterV0.getPCACandidatePos(cand);
  float dx = v0XYZ[0] - pv.getX(), dy = v0XYZ[1] - pv.getY(), dz = v0XYZ[2] - pv.getZ(), prodXYZv0 = dx * pV0[0] + dy * pV0[1] + dz * pV0[2];
  float cosPA = prodXYZv0 / std::sqrt((dx * dx + dy * dy + dz * dz) * p2V0);
  new (&v0) o2::dataformats::V0(v0XYZ, pV0, mFitterV0.calcPCACovMatrixFlat(cand), trPProp, trNProp);
  v0.setDCA(mFitterV0.getChi2AtPCACandidate(cand));
  v0.setCosPA(cosPA);
  return true;
}

const o2::dataformats::V0* SVStudySpec::getV0(int iv, o2::globaltracking::RecoContainer& recoData)
{
  // the V0 of the finder if it sent them (svertexer.createFullV0s), else the refit of this device, done once per TF
  if (!mRefit) {
    return &recoData.getV0s()[iv];
  }
  if (mRefitStatus[iv] < 0) {
    mRefitStatus[iv] = refitV0(recoData.getV0sIdx()[iv], mRefitV0s[iv], recoData) ? 1 : 0;
  }
  return mRefitStatus[iv] ? &mRefitV0s[iv] : nullptr;
}

void SVStudySpec::getTPCOccupancy(const PVertex& pv, float& occBef, float& occAft) const
{
  // TPC clusters before / after the time of the primary vertex (0 without the TPC occupancy map)
  occBef = occAft = 0.f;
  if (mTBinClOccBef.empty() || mTBinClOccAft.empty()) {
    return;
  }
  int tb = pv.getTimeStamp().getTimeStamp() * mTPCTBinMUSInv * mNTPCOccBinLengthInv;
  occBef = tb < 0 ? mTBinClOccBef[0] : (tb >= mTBinClOccBef.size() ? mTBinClOccBef.back() : mTBinClOccBef[tb]);
  occAft = tb < 0 ? mTBinClOccAft[0] : (tb >= mTBinClOccAft.size() ? mTBinClOccAft.back() : mTBinClOccAft[tb]);
}

bool SVStudySpec::passProngCuts(const V0ID& id, const o2::globaltracking::RecoContainer& recoData) const
{
  // cuts that need only the prong tracks: applied before the refit
  const auto& par = V0FullOutputParams::Instance();
  int charge[2] = {0, 0};
  for (int ip = 0; ip < 2; ip++) {
    const auto gid = id.getProngID(ip);
    charge[ip] = recoData.getTrackParam(gid).getCharge();
    const auto refs = recoData.getSingleDetectorRefs(gid);
    const bool hasTPC = refs[GTrackID::TPC].isIndexSet();
    if (par.requireTPCBoth && !hasTPC) {
      return false;
    }
    if (hasTPC && recoData.getTPCTrack(refs[GTrackID::TPC]).getNClusters() < par.minTPCClusters) {
      return false;
    }
  }
  return !par.requireOppositeCharges || charge[0] * charge[1] < 0;
}

int SVStudySpec::addProng(GTrackID gid, const o2::globaltracking::RecoContainer& recoData, FullOutput& out)
{
  // each prong is stored once per TF: the V0s sharing it point to the same entry
  const auto found = out.prongEntry.find(gid.getRawWOFlags());
  if (found != out.prongEntry.end()) {
    return found->second;
  }
  const int entry = out.prongs.size();
  out.prongEntry[gid.getRawWOFlags()] = entry;
  auto& prong = out.prongs.emplace_back();
  prong.gid = gid;
  prong.track = recoData.getTrackParam(gid);
  float t = 0.f, tErr = 0.f;
  recoData.getTrackTime(gid, t, tErr);
  if (gid.getSource() == GTrackID::TPC) { // TPC-only tracks give their time in time bins
    t *= mTPCTBinMUS;
    tErr *= mTPCTBinMUS;
  }
  prong.timeMUS = t;
  prong.timeErrMUS = tErr;
  const auto refs = recoData.getSingleDetectorRefs(gid);
  if (refs[GTrackID::ITSTPC].isIndexSet()) {
    prong.chi2ITSTPC = recoData.getTPCITSTrack(refs[GTrackID::ITSTPC]).getChi2Match();
  }
  // TPC: the attached clusters, in the order of the track's cluster references
  if (refs[GTrackID::TPC].isIndexSet()) {
    const auto& trk = recoData.getTPCTrack(refs[GTrackID::TPC]);
    prong.tpc = trk;
    prong.flags |= o2::dataformats::V0Prong::HasTPC;
    const auto& clusters = recoData.getTPCClusters();
    const auto clRefs = recoData.getTPCTracksClusterRefs();
    prong.tpcClusters.setFirstEntry(out.clTPC.size());
    for (int j = 0; j < trk.getNClusterReferences(); j++) {
      uint8_t sector = 0, row = 0;
      uint32_t clIdx = 0;
      trk.getClusterReference(clRefs, j, sector, row, clIdx);
      auto& cl = out.clTPC.emplace_back();
      cl.cl = clusters.clusters[sector][row][clIdx];
      cl.sector = sector;
      cl.row = row;
      if (!recoData.clusterShMapTPC.empty()) {
        cl.shMap = recoData.clusterShMapTPC[clusters.clusterOffset[sector][row] + clIdx];
      }
    }
    prong.tpcClusters.setEntries(out.clTPC.size() - prong.tpcClusters.getFirstEntry());
  }
  // ITS: clusters of the ITS track, or of the afterburner tracklet attached in the ITS-TPC matching
  std::vector<int> itsIndices;
  if (refs[GTrackID::ITS].isIndexSet()) {
    const auto& trk = recoData.getITSTrack(refs[GTrackID::ITS]);
    const auto clRefs = recoData.getITSTracksClusterRefs();
    for (int i = 0; i < trk.getNumberOfClusters(); i++) {
      itsIndices.push_back(clRefs[trk.getFirstClusterEntry() + i]);
    }
    prong.flags |= o2::dataformats::V0Prong::HasITS;
  } else if (refs[GTrackID::ITSAB].isIndexSet()) {
    const auto& ab = recoData.getITSABRef(refs[GTrackID::ITSAB]);
    const auto clRefs = recoData.getITSABClusterRefs();
    for (int i = 0; i < ab.getEntries(); i++) {
      itsIndices.push_back(clRefs[ab.getFirstEntry() + i]);
    }
    prong.flags |= o2::dataformats::V0Prong::HasITS | o2::dataformats::V0Prong::ITSAB;
  }
  const auto itsClusters = recoData.getITSClusters();
  const auto rofs = recoData.getITSClustersROFRecords();
  prong.itsClusters.setFirstEntry(out.clITS.size());
  for (const int idx : itsIndices) {
    if (idx < 0 || idx >= (int)itsClusters.size()) {
      continue;
    }
    const auto& c = itsClusters[idx];
    auto& cl = out.clITS.emplace_back();
    cl.chipID = c.getSensorID();
    cl.row = c.getRow();
    cl.col = c.getCol();
    cl.pattID = c.getPatternID();
    if (c.getPatternID() == o2::itsmft::CompCluster::InvalidPatternID || (mITSDict && mITSDict->isGroup(c.getPatternID()))) {
      out.pattRequests.push_back({int(out.clITS.size()) - 1, idx}); // pattern bytes in the TF's pattern stream, copied by fillITSPatterns
    }
    auto rof = std::upper_bound(rofs.begin(), rofs.end(), idx, [](int v, const o2::itsmft::ROFRecord& r) { return v < r.getFirstEntry(); });
    if (rof != rofs.begin()) {
      cl.rofBC = std::prev(rof)->getBCData().differenceInBC(mTFStart);
    }
  }
  prong.itsClusters.setEntries(out.clITS.size() - prong.itsClusters.getFirstEntry());
  // TOF: the matched cluster
  prong.tofClusters.setFirstEntry(out.clTOF.size());
  if (refs[GTrackID::TOF].isIndexSet()) {
    const auto tofClusters = recoData.getTOFClusters();
    const int idx = refs[GTrackID::TOF].getIndex();
    if (idx < (int)tofClusters.size()) {
      const auto& c = tofClusters[idx];
      auto& cl = out.clTOF.emplace_back();
      cl.timeRaw = c.getTimeRaw();
      cl.tot = c.getTot();
      cl.channel = c.getMainContributingChannel();
      prong.flags |= o2::dataformats::V0Prong::HasTOF;
    }
  }
  prong.tofClusters.setEntries(out.clTOF.size() - prong.tofClusters.getFirstEntry());
  // TRD: the tracklets of the TPC-TRD or ITS-TPC-TRD track
  prong.trdTracklets.setFirstEntry(out.trdTracklets.size());
  if (refs[GTrackID::TRD].isIndexSet()) {
    const auto& trk = recoData.getTrack<o2::trd::TrackTRD>(refs[GTrackID::TRD]);
    const auto tracklets = recoData.getTRDTracklets();
    const auto trigs = recoData.getTRDTriggerRecords();
    for (int layer = 0; layer < 6; layer++) {
      const int it = trk.getTrackletIndex(layer);
      if (it < 0 || it >= (int)tracklets.size()) {
        continue;
      }
      auto& tr = out.trdTracklets.emplace_back();
      tr.word = tracklets[it].getTrackletWord();
      tr.layer = layer;
      auto trig = std::upper_bound(trigs.begin(), trigs.end(), it, [](int v, const o2::trd::TriggerRecord& t) { return v < t.getFirstTracklet(); });
      if (trig != trigs.begin()) {
        tr.trigBC = std::prev(trig)->getBCData().differenceInBC(mTFStart);
      }
    }
    prong.flags |= o2::dataformats::V0Prong::HasTRD;
  }
  prong.trdTracklets.setEntries(out.trdTracklets.size() - prong.trdTracklets.getFirstEntry());
  return entry;
}

void SVStudySpec::fillITSPatterns(const o2::globaltracking::RecoContainer& recoData, FullOutput& out) const
{
  // the TF's pattern stream holds, in cluster order, the patterns of the clusters with an invalid or a group pattern ID
  if (!mITSDict) {
    if (!mNoDictWarned) {
      LOGP(warning, "No ITS cluster dictionary: patterns of ITS clusters of V0 prongs are not stored");
      mNoDictWarned = true;
    }
    return;
  }
  auto& requests = out.pattRequests;
  std::sort(requests.begin(), requests.end(), [](const ITSPattRequest& a, const ITSPattRequest& b) { return a.index < b.index; });
  const auto clusters = recoData.getITSClusters();
  auto pattIt = recoData.getITSClustersPatterns().begin();
  size_t ir = 0;
  for (int k = 0; k < (int)clusters.size() && ir < requests.size(); k++) {
    const auto pattID = clusters[k].getPatternID();
    if (pattID != o2::itsmft::CompCluster::InvalidPatternID && !mITSDict->isGroup(pattID)) {
      continue;
    }
    const auto start = pattIt;
    o2::itsmft::ClusterPattern::skipPattern(pattIt);
    for (; ir < requests.size() && requests[ir].index == k; ir++) {
      out.clITS[requests[ir].cluster].pattEntry = out.itsPatterns.size();
      out.itsPatterns.insert(out.itsPatterns.end(), start, pattIt);
    }
  }
}

namespace
{
/// reproducible pseudo-random number in [0, 1) of a V0 of a TF (splitmix64 of the TF counter and the V0 index)
double v0Random(uint32_t tfCounter, uint32_t iv)
{
  uint64_t z = (uint64_t(tfCounter) << 32 | iv) + 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  z = z ^ (z >> 31);
  return (z >> 11) * (1. / 9007199254740992.);
}
} // namespace

int SVStudySpec::countCompatiblePVs(const V0ID& id, const o2::globaltracking::RecoContainer& recoData) const
{
  // primary vertices whose time is compatible with both prongs: a TPC-only prong (on one side) by its TPC time window, else by its time
  const auto& par = V0FullOutputParams::Instance();
  float tMin = -1e9f, tMax = 1e9f;
  for (int ip = 0; ip < 2; ip++) {
    const auto gid = id.getProngID(ip);
    float lo = 0.f, hi = 0.f;
    if (gid.getSource() == GTrackID::TPC) {
      const auto& trk = recoData.getTPCTrack(gid);
      lo = (trk.getTime0() - trk.getDeltaTBwd()) * mTPCTBinMUS;
      hi = (trk.getTime0() + trk.getDeltaTFwd()) * mTPCTBinMUS;
    } else {
      float t = 0.f, tErr = 0.f;
      recoData.getTrackTime(gid, t, tErr);
      lo = t - par.pvTimeNSigma * tErr;
      hi = t + par.pvTimeNSigma * tErr;
    }
    tMin = std::max(tMin, lo);
    tMax = std::min(tMax, hi);
  }
  int n = 0;
  for (const auto& pv : recoData.getPrimaryVertices()) {
    const float t = pv.getTimeStamp().getTimeStamp(), err = par.pvTimeNSigma * pv.getTimeStamp().getTimeStampError();
    n += (t + err >= tMin && t - err <= tMax);
  }
  return n;
}

void SVStudySpec::processFull(o2::globaltracking::RecoContainer& recoData, ProcessingContext& pc)
{
  const auto& par = V0FullOutputParams::Instance();
  o2::dataformats::TFIDInfo tfID;
  o2::base::TFIDInfoHelper::fillTFIDInfo(pc, tfID);
  mTFStart = {0, tfID.firstTForbit};
  o2::dataformats::V0sTFInfo tfInfo;
  tfInfo.vDrift = mCorrMap->getVDrift();
  tfInfo.t0 = mCorrMap->getT0();
  tfInfo.lumi = mCorrMap->getLumi();
  tfInfo.idc = mCorrMap->getIDC();
  tfInfo.tpcOcc = recoData.occupancyMapTPC.empty() ? -1.f : recoData.occupancyMapTPC[0];
  tfInfo.nHBFPerTF = mNHBPerTF;

  using V0C = o2::dataformats::V0WithClusters;
  FullOutput out;
  const auto v0IDs = recoData.getV0sIdx();
  const float massK0 = o2::constants::physics::MassK0Short;
  const float massLambda = o2::constants::physics::MassLambda;
  // 1) selected candidates (no clusters yet), with their sampling class and pT bin
  struct Sel {
    int iv;
    int bin;
    double rnd;
    V0C cand;
  };
  std::vector<Sel> sel;
  const double logPtMin = std::log10(par.ptMin), logPtMax = std::log10(par.ptMax);
  // 0) the V0 finder builds a pair with TPC-only prongs once per compatible primary vertex (track shifted in z to the vertex time):
  //    per prong pair only the copy with the largest cosPA is kept; the number of copies and the second-best one are stored
  struct PairCopies {
    int best = -1;
    int second = -1;
    int n = 0;
  };
  std::unordered_map<uint64_t, PairCopies> pairs;
  std::vector<uint64_t> pairKey(v0IDs.size(), 0);
  size_t nCopies = 0;
  for (int iv = 0; iv < (int)v0IDs.size(); iv++) {
    const auto& id = v0IDs[iv];
    if (id.getVertexID() < 0 || !passProngCuts(id, recoData)) {
      continue;
    }
    const auto* v0 = getV0(iv, recoData);
    if (!v0) {
      continue;
    }
    pairKey[iv] = (uint64_t(id.getProngID(0).getRawWOFlags()) << 32) | id.getProngID(1).getRawWOFlags();
    auto& cp = pairs[pairKey[iv]]; // key != 0: the two prongs are different tracks
    cp.n++;
    nCopies++;
    const float cosPA = v0->getCosPA();
    if (cp.best < 0 || cosPA > getV0(cp.best, recoData)->getCosPA()) {
      cp.second = cp.best;
      cp.best = iv;
    } else if (cp.second < 0 || cosPA > getV0(cp.second, recoData)->getCosPA()) {
      cp.second = iv;
    }
  }
  for (int iv = 0; iv < (int)v0IDs.size(); iv++) {
    const auto& id = v0IDs[iv];
    if (!pairKey[iv]) { // no vertex, failed the prong cuts or the refit
      continue;
    }
    const auto& copies = pairs[pairKey[iv]];
    if (copies.best != iv) {
      continue;
    }
    const auto* v0 = getV0(iv, recoData);
    V0C cand;
    cand.nPairPV = copies.n;
    if (copies.second >= 0) {
      const auto& pv2 = recoData.getPrimaryVertex(v0IDs[copies.second].getVertexID());
      cand.cosPA2 = getV0(copies.second, recoData)->getCosPA();
      cand.pv2Z = pv2.getZ();
      cand.pv2TimeMUS = pv2.getTimeStamp().getTimeStamp();
    }
    cand.dca = v0->getDCA();
    cand.cosPA = v0->getCosPA();
    cand.qT = v0->calcAPQt();
    cand.alpha = v0->calcAPAlpha();
    cand.r = std::sqrt(v0->calcR2());
    cand.massK0 = std::sqrt(std::max(0.f, v0->calcMass2AsK0()));
    cand.massLambda = std::sqrt(std::max(0.f, v0->calcMass2AsLambda()));
    cand.massAntiLambda = std::sqrt(std::max(0.f, v0->calcMass2AsAntiLambda()));
    cand.massGamma = std::sqrt(std::max(0.f, v0->calcMass2AsPhoton()));
    // far V0s (beyond maxR, both legs without ITS) get their own cuts
    const bool far = par.selFar && cand.r > par.maxR && !id.getProngID(0).includesDet(DetID::ITS) && !id.getProngID(1).includesDet(DetID::ITS);
    const float maxDCA = far ? par.maxDCAFar : par.maxDCA, minCosPA = far ? par.minCosPAFar : par.minCosPA, maxR = far ? par.maxRFar : par.maxR;
    const float dmK0 = far ? par.maxDMassK0Far : par.maxDMassK0, dmL = far ? par.maxDMassLambdaFar : par.maxDMassLambda;
    const float minQt = far ? par.minQtFar : par.minQt;
    if (!(cand.dca > maxDCA || cand.cosPA < minCosPA || cand.qT < minQt || cand.qT > par.maxQt || cand.r < par.minR || cand.r > maxR)) {
      if (par.selK0 && std::abs(cand.massK0 - massK0) < dmK0) {
        cand.selection |= V0C::K0;
      }
      if (par.selLambda && cand.alpha > par.minAlphaLambda && cand.alpha < par.maxAlphaLambda && std::abs(cand.massLambda - massLambda) < dmL) {
        cand.selection |= V0C::Lambda;
      }
      if (par.selLambda && cand.alpha < -par.minAlphaLambda && cand.alpha > -par.maxAlphaLambda && std::abs(cand.massAntiLambda - massLambda) < dmL) {
        cand.selection |= V0C::AntiLambda;
      }
      if (far && cand.selection) {
        cand.selection |= V0C::Far;
      }
    }
    if (par.selGamma && (id.isPhotonOnly() || cand.massGamma < par.maxMassGamma) && cand.qT < par.maxQtGamma && cand.r > par.minRGamma && cand.r < par.maxRGamma &&
        cand.cosPA > par.minCosPAGamma && (par.maxDCAGamma < 0.f || cand.dca < par.maxDCAGamma)) {
      cand.selection |= V0C::Gamma;
    }
    if (!cand.selection) {
      continue;
    }
    cand.sampleClass = (cand.selection & V0C::K0) ? V0C::ClassK0 : ((cand.selection & V0C::Lambda) ? V0C::ClassLambda : ((cand.selection & V0C::AntiLambda) ? V0C::ClassAntiLambda : V0C::ClassGamma));
    if (id.isPhotonOnly()) {
      cand.selection |= V0C::PhotonOnly;
    }
    cand.v0 = *v0;
    cand.v0ID = id;
    const double logPt = std::log10(std::max(1e-6f, v0->getPt()));
    const int bin = logPt < logPtMin ? 0 : (logPt >= logPtMax ? par.nPtBins + 1 : 1 + int((logPt - logPtMin) / (logPtMax - logPtMin) * par.nPtBins));
    sel.push_back({iv, bin, v0Random(tfID.tfCounter, iv), cand});
  }
  // 1b) V0s sharing a track (before the downsampling, which may drop the competitors): count and greedy matching by legScore
  {
    std::unordered_map<uint32_t, std::vector<int>> byLeg;
    for (int k = 0; k < (int)sel.size(); k++) {
      auto& c = sel[k].cand;
      c.legScore = par.legScoreWeight * 60.f * std::acos(std::clamp(c.cosPA, -1.f, 1.f)) + (1.f - par.legScoreWeight) * std::sqrt(std::max(0.f, c.dca)) / 3.f;
      for (int ip = 0; ip < 2; ip++) {
        byLeg[c.v0ID.getProngID(ip).getRawWOFlags()].push_back(k);
      }
    }
    std::vector<int> order(sel.size());
    for (int k = 0; k < (int)sel.size(); k++) {
      order[k] = k;
      std::unordered_set<int> others;
      for (int ip = 0; ip < 2; ip++) {
        for (int o : byLeg[sel[k].cand.v0ID.getProngID(ip).getRawWOFlags()]) {
          if (o != k) {
            others.insert(o);
          }
        }
      }
      sel[k].cand.nSharedLeg = others.size();
    }
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return sel[a].cand.legScore < sel[b].cand.legScore; });
    std::unordered_set<uint32_t> usedLegs;
    for (int k : order) {
      const uint32_t l0 = sel[k].cand.v0ID.getProngID(0).getRawWOFlags();
      const uint32_t l1 = sel[k].cand.v0ID.getProngID(1).getRawWOFlags();
      if (!usedLegs.count(l0) && !usedLegs.count(l1)) {
        sel[k].cand.selection |= V0C::BestForLegs;
        usedLegs.insert(l0);
        usedLegs.insert(l1);
      }
    }
  }
  // 2) flat in V0 pT: per class and bin at most maxPerBin* V0s, those with the lowest random numbers; weight = candidates / kept
  const int maxPerBin[4] = {par.maxPerBinK0, par.maxPerBinLambda, par.maxPerBinLambda, par.maxPerBinGamma};
  std::map<std::pair<int, int>, std::vector<int>> groups;
  for (int k = 0; k < (int)sel.size(); k++) {
    groups[{sel[k].cand.sampleClass, sel[k].bin}].push_back(k);
  }
  std::vector<int> kept;
  for (auto& [key, idx] : groups) {
    const int cap = maxPerBin[key.first];
    if (cap >= 0 && (int)idx.size() > cap) {
      std::sort(idx.begin(), idx.end(), [&](int a, int b) { return sel[a].rnd < sel[b].rnd; });
    }
    const int nKeep = cap < 0 ? idx.size() : std::min<int>(cap, idx.size());
    for (int j = 0; j < nKeep; j++) {
      sel[idx[j]].cand.weight = float(idx.size()) / nKeep;
      kept.push_back(idx[j]);
    }
  }
  // 3) at most maxV0sPerTF, again at random
  if (par.maxV0sPerTF >= 0 && (int)kept.size() > par.maxV0sPerTF) {
    if (!mMaxV0sWarned) {
      LOGP(warning, "{} sampled V0s in TF {}: {} kept at random (svFullOutput.maxV0sPerTF; warned once)", kept.size(), tfID.tfCounter, par.maxV0sPerTF);
      mMaxV0sWarned = true;
    }
    std::sort(kept.begin(), kept.end(), [&](int a, int b) { return sel[a].rnd < sel[b].rnd; });
    const float scale = float(kept.size()) / par.maxV0sPerTF;
    kept.resize(par.maxV0sPerTF);
    for (int k : kept) {
      sel[k].cand.weight *= scale;
    }
  }
  std::sort(kept.begin(), kept.end()); // V0 finder order
  // 4) the kept V0s with their primary vertex, prongs and clusters
  for (int k : kept) {
    auto& cand = sel[k].cand;
    const auto& id = cand.v0ID;
    const auto& pv = recoData.getPrimaryVertex(id.getVertexID());
    getTPCOccupancy(pv, cand.tpcOccBef, cand.tpcOccAft);
    cand.pvX = pv.getX();
    cand.pvY = pv.getY();
    cand.pvZ = pv.getZ();
    cand.pvTimeMUS = pv.getTimeStamp().getTimeStamp();
    cand.pvTimeErrMUS = pv.getTimeStamp().getTimeStampError();
    cand.pvNContrib = pv.getNContributors();
    cand.nCompatPV = countCompatiblePVs(id, recoData);
    for (int ip = 0; ip < 2; ip++) {
      cand.prong[ip] = addProng(id.getProngID(ip), recoData, out);
    }
    out.v0s.push_back(cand);
    if (mFullDebugOut) {
      writeFullDebug(out.v0s.back(), out, pv, tfID.tfCounter);
    }
  }
  if (!out.pattRequests.empty()) {
    fillITSPatterns(recoData, out);
  }
  mNV0sFull += out.v0s.size();
  mNProngsFull += out.prongs.size();
  mNClTPCFull += out.clTPC.size();
  int nPerClass[4] = {0, 0, 0, 0};
  for (const auto& v : out.v0s) {
    nPerClass[v.sampleClass]++;
  }
  LOGP(info, "Full output of TF {}: {} V0s (K0s {}, Lambda {}, anti-Lambda {}, photons {}) of {} selected, {} candidates ({} after the prong cuts, {} prong pairs); {} prongs, {} TPC, {} ITS, {} TOF clusters, {} TRD tracklets",
       tfID.tfCounter, out.v0s.size(), nPerClass[0], nPerClass[1], nPerClass[2], nPerClass[3], sel.size(), v0IDs.size(), nCopies, pairs.size(), out.prongs.size(), out.clTPC.size(), out.clITS.size(),
       out.clTOF.size(), out.trdTracklets.size());
  pc.outputs().snapshot(Output{"GLO", "V0FULL", 0}, out.v0s);
  pc.outputs().snapshot(Output{"GLO", "V0FULLPRONG", 0}, out.prongs);
  pc.outputs().snapshot(Output{"GLO", "V0FULLCLTPC", 0}, out.clTPC);
  pc.outputs().snapshot(Output{"GLO", "V0FULLCLITS", 0}, out.clITS);
  pc.outputs().snapshot(Output{"GLO", "V0FULLITSPATT", 0}, out.itsPatterns);
  pc.outputs().snapshot(Output{"GLO", "V0FULLCLTOF", 0}, out.clTOF);
  pc.outputs().snapshot(Output{"GLO", "V0FULLTRDTRKL", 0}, out.trdTracklets);
  pc.outputs().snapshot(Output{"GLO", "V0FULLTF", 0}, tfInfo);
  pc.outputs().snapshot(Output{"GLO", "V0FULLOCCMAP", 0}, std::vector<unsigned int>(recoData.occupancyMapTPC.begin(), recoData.occupancyMapTPC.end()));
  pc.outputs().snapshot(Output{"GLO", "V0FULLTFID", 0}, tfID);
}

void SVStudySpec::writeFullDebug(const o2::dataformats::V0WithClusters& v0, const FullOutput& out, const PVertex& pv, int tfCounter) const
{
  // the stored TPC clusters of both prongs transformed with the map of the TF: z with the prong's time0 (the frame of its TPC track) and
  // zPV with the time of the primary vertex (absolute z, comparable to the V0 vertex)
  constexpr int NSectorsA = 18;
  const float pvTimeTB = pv.getTimeStamp().getTimeStamp() * mTPCTBinMUSInv;
  std::vector<int> clProng, clSector, clRow;
  std::vector<float> clPad, clTime, clX, clY, clZ, clZPV, clGx, clGy;
  float time0[2] = {-1.f, -1.f};
  int nClTPC[2] = {0, 0}, nClITS[2] = {0, 0}, prongFlags[2] = {0, 0};
  for (int ip = 0; ip < 2; ip++) {
    const auto& prong = out.prongs[v0.prong[ip]];
    prongFlags[ip] = prong.flags;
    nClTPC[ip] = prong.tpcClusters.getEntries();
    nClITS[ip] = prong.itsClusters.getEntries();
    if (!prong.hasTPC()) {
      continue;
    }
    time0[ip] = prong.tpc.getTime0();
    for (int i = prong.tpcClusters.getFirstEntry(); i < (int)prong.tpcClusters.getEntriesBound(); i++) {
      const auto& c = out.clTPC[i];
      float x = 0.f, y = 0.f, z = 0.f, xPV = 0.f, yPV = 0.f, zPV = 0.f;
      mCorrMap->Transform(c.sector, c.row, c.cl.getPad(), c.cl.getTime(), x, y, z, time0[ip]);
      mCorrMap->Transform(c.sector, c.row, c.cl.getPad(), c.cl.getTime(), xPV, yPV, zPV, pvTimeTB);
      const float alpha = o2::math_utils::sector2Angle(c.sector % NSectorsA);
      const float sinAlpha = std::sin(alpha), cosAlpha = std::cos(alpha);
      clProng.push_back(ip);
      clSector.push_back(c.sector);
      clRow.push_back(c.row);
      clPad.push_back(c.cl.getPad());
      clTime.push_back(c.cl.getTime());
      clX.push_back(x);
      clY.push_back(y);
      clZ.push_back(z);
      clZPV.push_back(zPV);
      clGx.push_back(x * cosAlpha - y * sinAlpha);
      clGy.push_back(x * sinAlpha + y * cosAlpha);
    }
  }
  int selection = v0.selection;
  float vx = v0.v0.getX(), vy = v0.v0.getY(), vz = v0.v0.getZ();
  float pvTimeMUS = pv.getTimeStamp().getTimeStamp();
  int nCompatPV = v0.nCompatPV, nPairPV = v0.nPairPV;
  float cosPA2 = v0.cosPA2, pv2Z = v0.pv2Z, pv2TimeMUS = v0.pv2TimeMUS;
  int nSharedLeg = v0.nSharedLeg;
  float legScore = v0.legScore;
  (*mFullDebugOut) << "v0"
                   << "tf=" << tfCounter << "selection=" << selection << "massK0=" << v0.massK0 << "massLambda=" << v0.massLambda
                   << "massAntiLambda=" << v0.massAntiLambda << "r=" << v0.r << "cosPA=" << v0.cosPA << "dca=" << v0.dca << "qT=" << v0.qT
                   << "alpha=" << v0.alpha << "vx=" << vx << "vy=" << vy << "vz=" << vz << "pvTimeMUS=" << pvTimeMUS
                   << "nCompatPV=" << nCompatPV << "nPairPV=" << nPairPV << "cosPA2=" << cosPA2 << "pv2Z=" << pv2Z << "pv2TimeMUS=" << pv2TimeMUS << "nSharedLeg=" << nSharedLeg << "legScore=" << legScore
                   << "time0P=" << time0[0] << "time0N=" << time0[1] << "nClTPCP=" << nClTPC[0] << "nClTPCN=" << nClTPC[1]
                   << "nClITSP=" << nClITS[0] << "nClITSN=" << nClITS[1] << "flagsP=" << prongFlags[0] << "flagsN=" << prongFlags[1]
                   << "clProng=" << clProng << "clSector=" << clSector << "clRow=" << clRow << "clPad=" << clPad << "clTime=" << clTime
                   << "clX=" << clX << "clY=" << clY << "clZ=" << clZ << "clZPV=" << clZPV << "clGx=" << clGx << "clGy=" << clGy
                   << "\n";
}

void SVStudySpec::endOfStream(EndOfStreamContext& ec)
{
  mDBGOut.reset();
  mFullDebugOut.reset();
  if (mFullOutput) {
    LOGP(info, "SV study full output: {} V0s, {} prongs, {} TPC clusters", mNV0sFull, mNProngsFull, mNClTPCFull);
  }
}

void SVStudySpec::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "CLUSDICT", 0)) {
    mITSDict = (const o2::itsmft::TopologyDictionary*)obj;
    return;
  }
}

DataProcessorSpec getSVStudySpec(GTrackID::mask_t srcTracks, GTrackID::mask_t srcCls, bool useMC, bool fullOutput)
{
  std::vector<OutputSpec> outputs;
  auto dataRequest = std::make_shared<DataRequest>();

  dataRequest->requestTracks(srcTracks, useMC);
  dataRequest->requestClusters(srcCls, false);
  dataRequest->requestPrimaryVertices(useMC);
  dataRequest->requestSecondaryVertices(useMC);
  dataRequest->inputs.emplace_back("meanvtx", "GLO", "MEANVERTEX", 0, Lifetime::Condition, ccdbParamSpec("GLO/Calib/MeanVertex", {}, 1));
  if (fullOutput) {
    dataRequest->inputs.emplace_back("corrMap", o2::header::gDataOriginTPC, "TPCCORRMAP", 0, Lifetime::Timeframe);
    outputs.emplace_back("GLO", "V0FULL", 0, Lifetime::Timeframe);
    outputs.emplace_back("GLO", "V0FULLPRONG", 0, Lifetime::Timeframe);
    outputs.emplace_back("GLO", "V0FULLCLTPC", 0, Lifetime::Timeframe);
    outputs.emplace_back("GLO", "V0FULLCLITS", 0, Lifetime::Timeframe);
    outputs.emplace_back("GLO", "V0FULLITSPATT", 0, Lifetime::Timeframe);
    outputs.emplace_back("GLO", "V0FULLCLTOF", 0, Lifetime::Timeframe);
    outputs.emplace_back("GLO", "V0FULLTRDTRKL", 0, Lifetime::Timeframe);
    outputs.emplace_back("GLO", "V0FULLTF", 0, Lifetime::Timeframe);
    outputs.emplace_back("GLO", "V0FULLOCCMAP", 0, Lifetime::Timeframe);
    outputs.emplace_back("GLO", "V0FULLTFID", 0, Lifetime::Timeframe);
  }
  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(true,                           // orbitResetTime
                                                              true,                           // GRPECS=true
                                                              false,                          // GRPLHCIF
                                                              true,                           // GRPMagField
                                                              true,                           // askMatLUT
                                                              o2::base::GRPGeomRequest::None, // geometry
                                                              dataRequest->inputs,
                                                              true);
  bool useTPCcl = srcCls[GTrackID::TPC];
  return DataProcessorSpec{
    "sv-study",
    dataRequest->inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<SVStudySpec>(dataRequest, ggRequest, srcTracks, useTPCcl, useMC, fullOutput)},
    Options{
      {"refit", VariantType::Bool, false, {"refit SVertices"}},
      {"sel-k0", VariantType::Float, -1.f, {"If positive, select K0s with this mass margin"}},
      {"max-eta", VariantType::Float, 1.2f, {"Cut on track eta"}},
      {"disable-study-tree", VariantType::Bool, false, {"do not write svStudy.root (e.g. with --full-output only)"}},
      {"debug-tree", VariantType::Bool, false, {"with --full-output: write v0s_full_debug.root with the transformed TPC clusters of the written V0s (test runs)"}},
    }};
}

} // namespace o2::svstudy
