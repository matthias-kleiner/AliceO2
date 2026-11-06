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

/// \file ClusterBiasCorrection.h
/// \brief Per-cluster correction of the two TPC cluster-position biases of the reconstruction that leak into the space-charge
/// distortion map, for the map macros (staticMapCreatorCPM.C, fineMapCreatorGBL.C). Header only, no ROOT.
///
/// The biases are cluster-level effects, not distortions:
/// - the pad-phase S-curve: the centre of gravity of the pad charges is pulled towards the pad centre (pad discretisation of the
///   charge cloud), largest at short drift (narrow cloud) and in the OROC (wider pads);
/// - the sector-edge truncation: near the pad plane's edge the charge cloud is cut off (no pad beyond the last pad), so the cluster
///   is pushed inward, growing with the drift length.
/// Both are tabulated against what the map macros have per cluster: the measured cluster y SNAPPED to the 1/64-pad grid of
/// ClusterNative (snapClusterY), the pad row, the drift length 250 cm - |z|, the track slope tgSlp = dy/dx at the row and the
/// TPC cluster flags (UnbinnedResid::channel). The correction returned by dy() is ADDED to the snapped cluster y.
///
/// Coordinates: d = distance of the cluster from the NEAREST pad-plane edge in pads (pad-0 centre at d = 0.5), "inward" = towards
/// larger d. The pad phase is phase = p - round(p) with p = d - 0.5 the pad coordinate counted from that edge (pad centres at
/// integers), so both edges fold onto one table. A bias b > 0 (inward) moves the cluster towards the sector centre, so the
/// correction is dy = +sign(y) * b on the +y side and -b on the -y side: it flips between the two edges.
///
/// Table file (plain text, '#' starts a comment, one record per line, values in um, whitespace separated):
///   version 1 | 2                      2: with the occupancy term below
///   unit    stack | region             granularity of the records' first index (optional, default stack): 4 stacks or the 10 pad
///                                      regions (the IROC regions differ in pad width, 4.16-4.36 mm)
///   drift   <n> <e_0> ... <e_n>        drift-length bin edges [cm]; the table values belong to the bin centres
///   slope   <n> <e_0> ... <e_n>        |tgSlp| bin edges; values at the bin centres
///   phase   <n>                        S-curve: n phase nodes phase_j = -0.5 + j / n (j = 0 ... n-1); with n = 64 the nodes are the
///                                      1/64-pad grid of ClusterNative, so a snapped cluster takes the value of its own grid point
///   dist    <n> <dLo> <dHi>            edge part: n bins of d in [dLo, dHi) [pads]
///   S <unit> <iDrift> <iSlope> <class> <v_0> ... <v_{nPhase-1}>   S-curve: mean inward bias [um] at the phase nodes
///   E <unit> <iDrift> <iSlope> <class> <v_0> ... <v_{nDist-1}>    edge part: mean inward bias [um] per d bin, S-curve subtracted
///   mult    <n> <m_0> ... <m_n-1>      occupancy nodes (TrackDataCompact::multStack of the cluster's stack), increasing
///   M <unit> <iDrift> <class> <iMult> <v_0> ... <v_{nDist-1}>     occupancy term: what the edge part changes at node iMult [um] per d bin; added to E, no slope axis
/// unit: the stack (0 IROC, 1 OROC1, 2 OROC2, 3 OROC3) or the pad region (0-9). The S-curve is a cluster property and therefore odd
/// in the phase and even in the slope (only |tgSlp| is binned); the producer keeps only the odd part of the measured profile. class: 0 golden (no flag of 0x19), 1 split pad / shared (0x1 or 0x10 without 0x8),
/// 2 single pad (flag 0x8 AND exactly on the pad centre: a single-pad cluster sits there; 0x8 elsewhere is a single-time cluster,
/// several pads, treated as class 0); a class without records falls back to class 0. The S-curve part is taken everywhere, the
/// edge part only for d inside [dLo, dHi).
///
/// Interpolation: periodic and linear in the phase between the nodes; linear in the drift length and in |tgSlp| between the
/// bin centres, clamped to the first / last centre outside; linear in d between the bin centres, clamped below the first centre,
/// zero at and beyond dHi. The occupancy term: linear between the nodes, clamped outside, in drift and d as the edge part; only with
/// a mult >= 0 given to inwardBias.

#ifndef O2_TPC_SCD_CLUSTERBIASCORRECTION_H
#define O2_TPC_SCD_CLUSTERBIASCORRECTION_H

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace o2::tpc::cluster_bias
{

/// the TPC cluster flags of UnbinnedResid::channel used for the class (ClusterNative flagSplitPad 0x1 / flagSingle 0x8, merger
/// shared 0x10)
constexpr int FlagSplitPad = 0x1;
constexpr int FlagSingle = 0x8;
constexpr int FlagShared = 0x10;
constexpr int NStacks = 4;
constexpr int NRegions = 10;
constexpr int NClasses = 3;
constexpr int NRows = 152;

/// drift length [cm] of a cluster at z [cm] as the tables bin it: 250 cm - |z|
constexpr float TPCLength = 250.f;
inline float driftLength(float z)
{
  return TPCLength - std::abs(z);
}

/// stack of a pad row: 0 IROC (rows 0-62), 1 OROC1 (63-96), 2 OROC2 (97-126), 3 OROC3 (127-151)
inline int stackOfRow(int row)
{
  return row < 63 ? 0 : (row < 97 ? 1 : (row < 127 ? 2 : 3));
}

/// pad region of a pad row (Mapper::ROWSPERREGION 17, 15, 16, 15, 18, 16, 16, 14, 13, 12)
inline int regionOfRow(int row)
{
  constexpr int FirstRowNextRegion[NRegions] = {17, 32, 48, 63, 81, 97, 113, 127, 140, 152};
  int region = 0;
  while (region < NRegions - 1 && row >= FirstRowNextRegion[region]) {
    ++region;
  }
  return region;
}

/// class of a cluster from its flags (UnbinnedResid::channel) and whether it sits exactly on a pad centre; a negative flags value
/// (no flags stored) counts as golden. flagSingle is single pad OR single time: only on the pad centre it is a single-pad cluster
inline int classOfFlags(int flags, bool onPadCentre)
{
  if (flags < 0) {
    return 0;
  }
  if (flags & FlagSingle) {
    return onPadCentre ? 2 : 0;
  }
  if (flags & (FlagSplitPad | FlagShared)) {
    return 1;
  }
  return 0;
}

struct ClusterBiasCorrection {
  std::vector<double> driftEdges;        ///< [cm]
  std::vector<double> slopeEdges;        ///< |tgSlp|
  int nPhase = 0;                        ///< S-curve phase nodes -0.5 + j / nPhase
  int nDist = 0;                         ///< edge-part d bins in [distLo, distHi)
  double distLo = 0.;                    ///< [pads]
  double distHi = 0.;                    ///< [pads]
  std::vector<double> scurve;            ///< [unit][drift][slope][class][phase] in cm
  std::vector<double> edge;              ///< [unit][drift][slope][class][d] in cm
  std::vector<char> hasScurve;           ///< [unit][drift][slope][class]: record present
  std::vector<char> hasEdge;             ///< [unit][drift][slope][class]: record present
  std::vector<double> multNodes;         ///< occupancy term: nodes in multStack
  std::vector<double> occupancy;         ///< [unit][drift][class][mult][d] in cm
  std::array<double, NRows> padWidth{};  ///< [cm] per pad row (setRowGeometry)
  std::array<double, NRows> halfPlane{}; ///< nPads / 2 per pad row: the pad plane's edge in pads from y = 0
  std::array<int, NRows> regionByRow{};  ///< pad region per row (setRowGeometry)
  bool byRegion = false;                 ///< records per pad region instead of per stack
  bool loaded = false;

  /// pad geometry of a row (from the Mapper in the macro): pad width [cm], number of pads, pad region (< 0: from the row)
  void setRowGeometry(int row, double width, int nPads, int region = -1)
  {
    if (row >= 0 && row < NRows) {
      padWidth[row] = width;
      halfPlane[row] = 0.5 * nPads;
      regionByRow[row] = region >= 0 ? region : regionOfRow(row);
    }
  }

  int nUnits() const { return byRegion ? NRegions : NStacks; }
  /// the record unit (stack or pad region) of a pad row
  int unitOfRow(int row) const { return byRegion ? regionByRow[row] : stackOfRow(row); }

  bool isActive() const { return loaded; }
  bool hasOccupancy() const { return loaded && !multNodes.empty(); }
  int nDrift() const { return int(driftEdges.size()) - 1; }
  int nSlope() const { return int(slopeEdges.size()) - 1; }

  /// index of a (unit, drift bin, slope bin, class) record
  size_t record(int unit, int iDrift, int iSlope, int cls) const
  {
    return ((size_t(unit) * nDrift() + iDrift) * nSlope() + iSlope) * NClasses + cls;
  }

  /// read the table; false (and inactive) on any format error
  bool load(const std::string& fileName, std::string* error = nullptr)
  {
    loaded = false;
    byRegion = false;
    std::ifstream in(fileName);
    if (!in) {
      return fail(error, "cannot open " + fileName);
    }
    std::string line;
    int version = 0;
    std::vector<std::pair<char, std::vector<double>>> records;
    while (std::getline(in, line)) {
      const auto hash = line.find('#');
      if (hash != std::string::npos) {
        line.erase(hash);
      }
      std::istringstream ss(line);
      std::string key;
      if (!(ss >> key)) {
        continue;
      }
      if (key == "version") {
        ss >> version;
      } else if (key == "unit") {
        std::string unit;
        ss >> unit;
        if (unit != "stack" && unit != "region") {
          return fail(error, "unknown unit " + unit);
        }
        byRegion = unit == "region";
      } else if (key == "drift" || key == "slope") {
        int n = 0;
        ss >> n;
        std::vector<double> e(n + 1);
        for (auto& v : e) {
          ss >> v;
        }
        if (!ss || n < 1 || !std::is_sorted(e.begin(), e.end())) {
          return fail(error, "bad " + key + " line");
        }
        (key == "drift" ? driftEdges : slopeEdges) = e;
      } else if (key == "mult") {
        int n = 0;
        ss >> n;
        multNodes.resize(n);
        for (auto& v : multNodes) {
          ss >> v;
        }
        if (!ss || n < 1 || !std::is_sorted(multNodes.begin(), multNodes.end())) {
          return fail(error, "bad mult line");
        }
      } else if (key == "phase") {
        ss >> nPhase;
      } else if (key == "dist") {
        ss >> nDist >> distLo >> distHi;
      } else if (key == "S" || key == "E" || key == "M") {
        std::vector<double> v;
        double x = 0.;
        while (ss >> x) {
          v.push_back(x);
        }
        records.emplace_back(key[0], v);
      } else {
        return fail(error, "unknown record " + key);
      }
    }
    if ((version != 1 && version != 2) || driftEdges.size() < 2 || slopeEdges.size() < 2 || nPhase < 2 || nDist < 0 || (nDist > 0 && distHi <= distLo)) {
      return fail(error, "missing or bad header (version / drift / slope / phase / dist)");
    }
    const size_t nRec = size_t(nUnits()) * nDrift() * nSlope() * NClasses;
    scurve.assign(nRec * nPhase, 0.);
    edge.assign(nRec * std::max(nDist, 1), 0.);
    hasScurve.assign(nRec, 0);
    hasEdge.assign(nRec, 0);
    occupancy.assign(multNodes.empty() ? 0 : size_t(nUnits()) * nDrift() * NClasses * multNodes.size() * nDist, 0.);
    for (const auto& [type, v] : records) {
      if (type == 'M') {
        if (int(v.size()) != 4 + nDist || multNodes.empty()) {
          return fail(error, "record M: wrong number of values or no mult line");
        }
        const int unit = int(v[0]);
        const int iDrift = int(v[1]);
        const int cls = int(v[2]);
        const int iMult = int(v[3]);
        if (unit < 0 || unit >= nUnits() || iDrift < 0 || iDrift >= nDrift() || cls < 0 || cls >= NClasses || iMult < 0 || iMult >= int(multNodes.size())) {
          return fail(error, "record M: index out of range");
        }
        for (int i = 0; i < nDist; ++i) {
          occupancy[occupancyRecord(unit, iDrift, cls, iMult) + i] = 1.e-4 * v[4 + i]; // um -> cm
        }
        continue;
      }
      const int nVal = type == 'S' ? nPhase : nDist;
      if (int(v.size()) != 4 + nVal) {
        return fail(error, std::string("record ") + type + ": wrong number of values");
      }
      const int unit = int(v[0]);
      const int iDrift = int(v[1]);
      const int iSlope = int(v[2]);
      const int cls = int(v[3]);
      if (unit < 0 || unit >= nUnits() || iDrift < 0 || iDrift >= nDrift() || iSlope < 0 || iSlope >= nSlope() || cls < 0 || cls >= NClasses) {
        return fail(error, std::string("record ") + type + ": index out of range");
      }
      const size_t r = record(unit, iDrift, iSlope, cls);
      auto& dst = type == 'S' ? scurve : edge;
      for (int i = 0; i < nVal; ++i) {
        dst[r * nVal + i] = 1.e-4 * v[4 + i]; // um -> cm
      }
      (type == 'S' ? hasScurve : hasEdge)[r] = 1;
    }
    loaded = true;
    return true;
  }

  /// correction [cm] to ADD to the snapped cluster y (sector frame) of a TPC cluster on pad row `row`; drift = 250 cm - |z|,
  /// tSlope = dy/dx of the track at the row, flags = UnbinnedResid::channel. 0 if inactive or the row has no geometry.
  float dy(int row, float ySnapped, float drift, float tSlope, int flags) const
  {
    if (!loaded || row < 0 || row >= NRows || padWidth[row] <= 0.) {
      return 0.f;
    }
    const double side = ySnapped < 0.f ? -1. : 1.;
    const double d = halfPlane[row] - std::abs(double(ySnapped)) / padWidth[row]; // pads from the nearest edge
    const double pad = d - 0.5;
    const double phase = pad - std::round(pad);
    const int unit = unitOfRow(row);
    const int cls = classOfFlags(flags, std::abs(phase) < 0.25 / 64.); // snapped: on the pad centre to the float precision of y
    const double bias = interpolate(unit, cls, drift, std::abs(tSlope), [&](size_t r, bool& ok) { return scurveAt(r, phase, ok); }) +
                        interpolate(unit, cls, drift, std::abs(tSlope), [&](size_t r, bool& ok) { return edgeAt(r, d, ok); });
    return float(side * bias);
  }

  /// inward bias [cm] (S-curve + edge part, with mult >= 0 and a table with it the occupancy term) at a distance d [pads] from the nearest edge
  double inwardBias(int unit, int cls, double drift, double absSlope, double d, double mult = -1.) const
  {
    const double pad = d - 0.5;
    const double phase = pad - std::round(pad);
    return interpolate(unit, cls, drift, absSlope, [&](size_t r, bool& ok) { return scurveAt(r, phase, ok); }) +
           interpolate(unit, cls, drift, absSlope, [&](size_t r, bool& ok) { return edgeAt(r, d, ok); }) +
           (mult >= 0. ? occupancyAt(unit, cls, drift, d, mult) : 0.);
  }

 private:
  static bool fail(std::string* error, const std::string& message)
  {
    if (error) {
      *error = message;
    }
    return false;
  }

  /// start of the d values of an occupancy record
  size_t occupancyRecord(int unit, int iDrift, int cls, int iMult) const
  {
    return (((size_t(unit) * nDrift() + iDrift) * NClasses + cls) * multNodes.size() + iMult) * nDist;
  }

  /// occupancy term [cm]: linear in drift, mult (between the nodes, clamped outside) and d (as edgeAt)
  double occupancyAt(int unit, int cls, double drift, double d, double mult) const
  {
    if (multNodes.empty() || nDist < 1 || d >= distHi) {
      return 0.;
    }
    int iD = 0;
    double fD = 0.;
    bracket(driftEdges, drift, iD, fD);
    const int nM = int(multNodes.size());
    int iM = 0;
    double fM = 0.;
    if (nM > 1 && mult > multNodes[0]) {
      while (iM < nM - 2 && mult >= multNodes[iM + 1]) {
        ++iM;
      }
      fM = std::min(1., (mult - multNodes[iM]) / (multNodes[iM + 1] - multNodes[iM]));
    }
    const double w = (distHi - distLo) / nDist;
    const double x = std::max(0., (d - distLo) / w - 0.5);
    const int iX = std::min(int(std::floor(x)), nDist - 1);
    const double fX = iX == nDist - 1 ? 0. : x - iX;
    double sum = 0.;
    for (int a = 0; a < 2; ++a) {
      const double wD = a ? fD : 1. - fD;
      for (int b = 0; b < 2; ++b) {
        const double wM = b ? fM : 1. - fM;
        if (wD == 0. || wM == 0.) {
          continue;
        }
        const double* v = &occupancy[occupancyRecord(unit, iD + a, cls, iM + b)];
        sum += wD * wM * ((1. - fX) * v[iX] + (fX > 0. ? fX * v[iX + 1] : 0.));
      }
    }
    return sum;
  }

  /// position of x between the bin centres of edges: lower centre index and weight of the upper one, clamped at both ends
  static void bracket(const std::vector<double>& edges, double x, int& i0, double& f)
  {
    const int n = int(edges.size()) - 1;
    const double c0 = 0.5 * (edges[0] + edges[1]);
    const double cN = 0.5 * (edges[n - 1] + edges[n]);
    if (n == 1 || x <= c0) {
      i0 = 0;
      f = 0.;
      return;
    }
    if (x >= cN) {
      i0 = n - 2;
      f = 1.;
      return;
    }
    i0 = 0;
    while (i0 < n - 2 && x >= 0.5 * (edges[i0 + 1] + edges[i0 + 2])) {
      ++i0;
    }
    const double a = 0.5 * (edges[i0] + edges[i0 + 1]);
    const double b = 0.5 * (edges[i0 + 1] + edges[i0 + 2]);
    f = (x - a) / (b - a);
  }

  /// bilinear in (drift, |slope|) over the records of (unit, class), falling back to class 0 where a record is missing
  template <typename F>
  double interpolate(int unit, int cls, double drift, double absSlope, F&& value) const
  {
    int iD = 0;
    int iS = 0;
    double fD = 0.;
    double fS = 0.;
    bracket(driftEdges, drift, iD, fD);
    bracket(slopeEdges, absSlope, iS, fS);
    double sum = 0.;
    for (int a = 0; a < 2; ++a) {
      const double wD = a ? fD : 1. - fD;
      if (wD == 0. || (a && nDrift() == 1)) {
        continue;
      }
      for (int b = 0; b < 2; ++b) {
        const double wS = b ? fS : 1. - fS;
        if (wS == 0. || (b && nSlope() == 1)) {
          continue;
        }
        bool ok = false;
        double v = value(record(unit, iD + a, iS + b, cls), ok);
        if (!ok && cls != 0) {
          v = value(record(unit, iD + a, iS + b, 0), ok);
        }
        sum += wD * wS * (ok ? v : 0.);
      }
    }
    return sum;
  }

  /// S-curve of a record at a phase in [-0.5, 0.5], periodic linear interpolation between the nodes -0.5 + j / nPhase
  double scurveAt(size_t r, double phase, bool& ok) const
  {
    ok = hasScurve[r];
    if (!ok) {
      return 0.;
    }
    double x = (phase + 0.5) * nPhase; // in node units
    if (std::abs(x - std::round(x)) < 1.e-3) {
      x = std::round(x); // on a node (a snapped cluster with 64 nodes): its own value, no admixture from the float rounding of y
    }
    const double fl = std::floor(x);
    const double f = x - fl;
    const int i0 = ((int(fl) % nPhase) + nPhase) % nPhase;
    const int i1 = (i0 + 1) % nPhase;
    const double* v = &scurve[r * nPhase];
    return (1. - f) * v[i0] + f * v[i1];
  }

  /// edge part of a record at d [pads]: linear between the bin centres, clamped below the first centre, 0 from distHi on
  double edgeAt(size_t r, double d, bool& ok) const
  {
    ok = nDist > 0 && hasEdge[r];
    if (!ok || d >= distHi) {
      return 0.;
    }
    const double w = (distHi - distLo) / nDist;
    const double x = (d - distLo) / w - 0.5;
    const double* v = &edge[r * nDist];
    if (x <= 0.) {
      return v[0];
    }
    const int i0 = int(std::floor(x));
    if (i0 >= nDist - 1) {
      // between the last centre and distHi: towards 0 at distHi
      const double f = (d - (distHi - 0.5 * w)) / (0.5 * w);
      return (1. - f) * v[nDist - 1];
    }
    const double f = x - i0;
    return (1. - f) * v[i0] + f * v[i0 + 1];
  }
};

} // namespace o2::tpc::cluster_bias

#endif
