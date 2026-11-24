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
/// One entry per TF. The V0s of a TF share their prongs: every prong track is stored once in the prong table, its raw clusters once in the
/// per-TF cluster tables, which the prong addresses with index ranges. Only raw detector data are stored (TPC ClusterNative, ITS compact
/// clusters + patterns, raw TOF time, TRD tracklet words); the TPC transformation, the ITS cluster dictionary and the geometry are applied offline. To reproduce the
/// reconstruction's TPC refit, the per-TF TPC occupancy map (cluster errors) and the cluster sharing-map bytes are stored as well.

#ifndef ALICEO2_V0_WITH_CLUSTERS_H
#define ALICEO2_V0_WITH_CLUSTERS_H

#include <cstdint>
#include <Rtypes.h>
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"
#include "CommonDataFormat/RangeReference.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "DataFormatsTPC/ClusterNative.h"

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
  int32_t pattEntry = -1; ///< start of the pattern bytes in the TF's pattern table (pattern not in the dictionary or group pattern), -1: none
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

/// prong track of one or more V0s of the TF
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
                                           ///< within the TF: the clusters are in the TF's TPC cluster table at tpcClusters
  float timeMUS = 0.f;                     ///< time of the global track [mus since the TF start]
  float timeErrMUS = 0.f;                  ///< its error [mus]
  float chi2ITSTPC = -1.f;                 ///< chi2 of the ITS-TPC match (< 0: no ITS-TPC match)
  RangeReference<int, int> tpcClusters{};  ///< first entry and number of the prong's TPC clusters in the TF's TPC cluster table (in the
                                           ///< order of the track's cluster references)
  RangeReference<int, int> itsClusters{};  ///< same for the ITS clusters in the TF's ITS cluster table
  RangeReference<int, int> tofClusters{};  ///< same for the matched TOF cluster (0 or 1) in the TF's TOF cluster table
  RangeReference<int, int> trdTracklets{}; ///< same for the TRD tracklets in the TF's TRD tracklet table
  uint8_t flags = 0;                       ///< Flags
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
    K0 = 0x1,         ///< passed the K0s selection
    Lambda = 0x2,     ///< passed the Lambda selection (proton = positive prong)
    AntiLambda = 0x4, ///< passed the anti-Lambda selection
    PhotonOnly = 0x8  ///< the V0 finder flagged the candidate as photon only (V0Index::isPhotonOnly)
  };
  V0 v0{};                    ///< V0 refitted by the study device (DCAFitter with the V0 finder's settings), prongs at the decay vertex
  V0Index v0ID{};             ///< prong IDs and primary vertex ID of the V0 finder (only valid within the TF)
  int prong[2] = {-1, -1};    ///< entries of the positive / negative prong in the TF's prong table
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
  ClassDefNV(V0WithClusters, 1);
};

/// per-TF quantities: the TPC transformation used in the reconstruction and the TPC occupancy
struct V0sTFInfo {
  float vDrift = 0.f;  ///< drift velocity of the transformation [cm/time bin]
  float t0 = 0.f;      ///< time offset of the transformation [time bins]
  float lumi = 0.f;    ///< luminosity estimator the map was scaled with (TPCFastTransformPOD::getLumi)
  float idc = 0.f;     ///< IDC estimator the map was scaled with (TPCFastTransformPOD::getIDC)
  float tpcOcc = -1.f; ///< TPC clusters per orbit, average of the TF (TPC occupancy map) (< 0: no map)
  int nHBFPerTF = 0;   ///< HBFs per TF (needed with the occupancy map tpcOccMap in the TPC refit)
  ClassDefNV(V0sTFInfo, 1);
};

} // namespace o2::dataformats

#endif
