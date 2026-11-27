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

/// \file V0WithClusters.h
/// \brief Selected V0s with the raw clusters of their prongs, for offline refits (re-transformation of the clusters with any TPC map)
///
/// One entry per TF: the selected V0s, each self-contained (as TrackCosmicsExtended): its two prongs with their own raw clusters. A track
/// used by several V0s of the TF is stored with each of them (V0WithClusters::nSharedLeg; prong gid identifies it within the TF). Per TF only
/// the TPC transformation info, the TPC occupancy map and the TF ID. Only raw detector data are stored (TPC ClusterNative, ITS compact
/// clusters + patterns, raw TOF time, TRD tracklet words); the TPC transformation, the ITS cluster dictionary and the geometry are applied offline. To reproduce the
/// reconstruction's TPC refit, the per-TF TPC occupancy map (cluster errors) and the cluster sharing-map bytes are stored as well.

#ifndef ALICEO2_V0_WITH_CLUSTERS_H
#define ALICEO2_V0_WITH_CLUSTERS_H

#include <cstdint>
#include <vector>
#include <Rtypes.h>
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "DataFormatsTPC/ClusterNative.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::dataformats
{

/// raw TPC cluster attached to a prong's TPC track, with its address
struct V0TPCCluster {
  o2::tpc::ClusterNative cl{}; ///< raw cluster: time, pad, widths, charges, flags
  uint8_t sector = 0;
  uint8_t row = 0;
  uint8_t shMap = 0; ///< byte of the reconstruction's TPC cluster sharing map (shared-cluster flags; the sharedmap input of the GPU refit,
                     ///< which changes the refitted track noticeably)
  ClassDefNV(V0TPCCluster, 1);
};

/// ITS cluster of a prong: the raw compact cluster (chip, anchor pixel, pattern ID) and its ROF; no coordinates
struct V0ITSCluster {
  uint16_t chipID = 0;
  uint16_t row = 0;
  uint16_t col = 0;
  uint16_t pattID = 0;
  int32_t pattEntry = -1; ///< start of the pattern bytes in the prong's V0Prong::itsPatterns (pattern not in the dictionary or group pattern), -1: none
  int32_t rofBC = 0;      ///< start of the cluster's ROF in BCs since the start of the TF
  ClassDefNV(V0ITSCluster, 1);
};

/// TOF cluster matched to a prong
struct V0TOFCluster {
  double timeRaw = 0.; ///< raw TOF time [ps] (the calibration is applied offline)
  float tot = 0.f;     ///< time over threshold
  int32_t channel = -1;
  ClassDefNV(V0TOFCluster, 1);
};

/// TRD tracklet attached to a prong's TRD track
struct V0TRDTracklet {
  uint64_t word = 0;  ///< raw Tracklet64 word
  int32_t trigBC = 0; ///< BC of its trigger since the start of the TF
  uint8_t layer = 0;
  ClassDefNV(V0TRDTracklet, 1);
};

/// prong track of a V0 with its raw clusters
struct V0Prong {
  enum Flags : uint8_t {
    HasTPC = 0x1, ///< has a TPC track (tpc is valid, TPC clusters stored)
    HasITS = 0x2, ///< has an ITS part (ITS clusters stored)
    ITSAB = 0x4,  ///< the ITS part is an afterburner tracklet (ITS clusters attached to the TPC track in the ITS-TPC matching), not an ITS track
    HasTOF = 0x8, ///< matched to a TOF cluster (stored)
    HasTRD = 0x10 ///< has a TRD track (its tracklets stored)
  };
  GlobalTrackID gid{};                     ///< global track ID of the prong (only valid within the TF)
  o2::track::TrackParCov track{};          ///< global track parameters as used by the V0 finder (inner point of the track)
  o2::tpc::TrackTPC tpc{};                 ///< TPC-only track (default if none); z refers to its time0; its cluster references are only valid
                                           ///< within the TF: the clusters are in clTPC
  float timeMUS = 0.f;                     ///< time of the global track [mus since the TF start]
  float timeErrMUS = 0.f;                  ///< its error [mus]
  float chi2ITSTPC = -1.f;                 ///< chi2 of the ITS-TPC match (< 0: no ITS-TPC match)
  std::vector<V0TPCCluster> clTPC;         ///< TPC clusters of the TPC track, in the order of its cluster references
  std::vector<V0ITSCluster> clITS;         ///< ITS clusters of the ITS track or of the afterburner tracklet
  std::vector<uint8_t> itsPatterns;        ///< pattern bytes (row span, column span, bitmap) of the ITS clusters that need them
  std::vector<V0TOFCluster> clTOF;         ///< the matched TOF cluster (0 or 1)
  std::vector<V0TRDTracklet> trdTracklets; ///< TRD tracklets of the TRD track
  uint8_t flags = 0;                       ///< Flags
  // TOF match of the reconstruction (HasTOF), for TOF PID: beta = tofLength / (c (tofSignal - V0WithClusters::evTimeTOF))
  double tofSignal = 0.;       ///< calibrated TOF time of the match [ps since the TF start] (MatchInfoTOF::getSignal)
  float tofLength = 0.f;       ///< integrated track length to the TOF cluster [cm] (from the track's start, as in the reconstruction)
  float tofExpTime[5] = {0.f}; ///< expected times of flight [ps] for e, mu, pi, K, p (o2::track::PID 0-4) over tofLength
  o2::MCCompLabel mcLabel{};   ///< MC label of the global track (MC input only; else not set)
  bool hasTPC() const { return flags & HasTPC; }
  bool hasITS() const { return flags & HasITS; }
  bool isITSAB() const { return flags & ITSAB; }
  bool hasTOF() const { return flags & HasTOF; }
  bool hasTRD() const { return flags & HasTRD; }
  ClassDefNV(V0Prong, 1);
};

/// selected V0 with its selection variables
struct V0WithClusters {
  enum Selection : uint8_t {
    K0 = 0x1,           ///< passed the K0s selection
    Lambda = 0x2,       ///< passed the Lambda selection (proton = positive prong)
    AntiLambda = 0x4,   ///< passed the anti-Lambda selection
    PhotonOnly = 0x8,   ///< the V0 finder flagged the candidate as photon only (V0Index::isPhotonOnly)
    Gamma = 0x10,       ///< passed the photon-conversion selection
    BestForLegs = 0x20, ///< best V0 for its tracks: greedy matching over the selected V0s of the TF, each track in at most one, lowest
                        ///< w 60 acos(cosPA) + (1 - w) sqrt(dca) / 3 cm first (w = svFullOutput.legScoreWeight)
    Far = 0x40          ///< selected with the far cuts (svFullOutput.*Far): beyond maxR, both legs without ITS
  };
  enum SampleClass : uint8_t { ClassK0 = 0,
                               ClassLambda = 1,
                               ClassAntiLambda = 2,
                               ClassGamma = 3 };
  V0 v0{};                    ///< the V0 finder's V0 if it sent them (svertexer.createFullV0s), else the study device's refit (DCAFitter with
                              ///< the V0 finder's settings); prongs at the decay vertex
  V0Prong prongPos{};         ///< positive prong with its clusters (gid: the V0 finder's prong ID)
  V0Prong prongNeg{};         ///< negative prong
  float dca = 0.f;            ///< V0::getDCA of the refit (chi2 at the point of closest approach of the prongs, as cut in the V0 finder)
  float cosPA = 0.f;          ///< cosine of the pointing angle to the primary vertex
  float qT = 0.f;             ///< Armenteros qT [GeV]
  float alpha = 0.f;          ///< Armenteros alpha
  float r = 0.f;              ///< transverse decay radius [cm]
  float massK0 = 0.f;         ///< invariant mass as pi+ pi- [GeV]
  float massLambda = 0.f;     ///< as p pi- [GeV]
  float massAntiLambda = 0.f; ///< as pbar pi+ [GeV]
  float massGamma = 0.f;      ///< as e+ e- [GeV]
  float tpcOccBef = 0.f;      ///< TPC clusters before the time of the primary vertex, as tpcOccBef in svStudy.root (0: no occupancy map)
  float tpcOccAft = 0.f;      ///< same after
  uint8_t selection = 0;      ///< Selection bits
  uint8_t sampleClass = 0;    ///< SampleClass of the downsampling (first passed type of K0s, Lambda, anti-Lambda, photon)
  float weight = 1.f;         ///< downsampling weight: V0s of its class (far ones apart: bit Far) and pT bin in the TF / kept ones (1: all kept)
  float pvX = 0.f;            ///< primary vertex the V0 finder assigned [cm]
  float pvY = 0.f;
  float pvZ = 0.f;
  float pvTimeMUS = 0.f;    ///< its time [mus since the TF start]
  float pvTimeErrMUS = 0.f; ///< and the time error
  int pvNContrib = 0;       ///< its number of contributors
  int nCompatPV = 0;        ///< primary vertices of the TF compatible in time with both prongs (TPC-only prong: its TPC time window;
                            ///< else its time +- svFullOutput.pvTimeNSigma errors): 1 = unique, e.g. the z of TPC-only photon legs
                            ///< from the vertex time
  int nPairPV = 1;          ///< primary vertices the V0 finder built this prong pair with (TPC-only prongs: one copy per vertex, shifted
                            ///< in z to its time); only the copy with the largest cosPA is written (this one)
  float cosPA2 = -2.f;      ///< cosPA of the copy with the second-largest cosPA (-2: no other copy)
  int nSharedLeg = 0;       ///< other selected V0s of the TF (any class, before the downsampling) using one of its tracks (stored with each)

  float pvCov[6] = {0.f};               ///< position covariance of the assigned primary vertex [cm^2] (VertexBase order: xx, xy, yy, xz, yz, zz)
  float legDCAxy[2] = {-999.f, -999.f}; ///< DCA of the positive / negative leg to that vertex in the bending plane [cm] (the V0's prong
                                        ///< propagated from the decay vertex, TPC-only legs in z shifted to the vertex time; -999: failed)
  float legDCAz[2] = {-999.f, -999.f};  ///< same along z [cm]
  double evTimeTOF = 0.;                ///< event time of that vertex from its TOF-matched tracks [ps since the TF start]: robust mean of
                                        ///< (TOF signal - expected pion time), the V0's own legs excluded (evTimeTOFN = 0: none)
  float evTimeTOFErr = 0.f;             ///< its error [ps] (spread / sqrt(N))
  int evTimeTOFN = 0;                   ///< TOF tracks used

  const V0Prong& prong(int ip) const { return ip == 0 ? prongPos : prongNeg; } ///< 0: positive, 1: negative (as V0::getProng)
  V0Prong& prong(int ip) { return ip == 0 ? prongPos : prongNeg; }

  // decay-tagged species (no data member: follows from selection and the prong order), e.g. TTree::Draw selections "v0s.isLambda()",
  // "v0s.legPID(0)==4" (protons of Lambda decays)
  static constexpr uint8_t TypeBits = K0 | Lambda | AntiLambda | Gamma;
  bool isK0() const { return (selection & TypeBits) == K0; }                 ///< K0s and no other type (far V0s included: isFar())
  bool isLambda() const { return (selection & TypeBits) == Lambda; }         ///< Lambda and no other type
  bool isAntiLambda() const { return (selection & TypeBits) == AntiLambda; } ///< anti-Lambda and no other type
  bool isGamma() const { return (selection & TypeBits) == Gamma; }           ///< photon conversion and no other type
  bool isFar() const { return selection & Far; }
  bool isBestForLegs() const { return selection & BestForLegs; }
  /// species of prong ip (0 positive, 1 negative) from the decay, as o2::track::PID::ID (0 e, 2 pi, 4 p): K0s pi+ pi-, Lambda p pi-,
  /// anti-Lambda pi+ pbar, photon e+ e-; -1 if not exactly one type bit is set (far V0s in the wide K0s and Lambda windows can have two)
  int legPID(int ip) const
  {
    switch (selection & TypeBits) {
      case K0:
        return o2::track::PID::Pion;
      case Lambda:
        return ip == 0 ? o2::track::PID::Proton : o2::track::PID::Pion;
      case AntiLambda:
        return ip == 0 ? o2::track::PID::Pion : o2::track::PID::Proton;
      case Gamma:
        return o2::track::PID::Electron;
      default:
        return -1;
    }
  }
  ClassDefNV(V0WithClusters, 1);
};

/// per-TF quantities: the TPC transformation used in the reconstruction and the TPC occupancy
struct V0sTFInfo {
  float vDrift = 0.f; ///< drift velocity of the transformation [cm/time bin]
  float t0 = 0.f;     ///< time offset of the transformation [time bins]
  float lumi = 0.f;   ///< luminosity estimator the map was scaled with (TPCFastTransformPOD::getLumi)
  float idc = 0.f;    ///< IDC estimator the map was scaled with (TPCFastTransformPOD::getIDC); TPC occupancy: branch tpcOccMap; the
                      ///< refit with it also needs the HBFs per TF (GRPECSObject::getNHBFPerTF of the run)
  ClassDefNV(V0sTFInfo, 1);
};

} // namespace o2::dataformats

#endif
