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
///   version 1
///   drift   <n> <e_0> ... <e_n>        drift-length bin edges [cm]; the table values belong to the bin centres
///   slope   <n> <e_0> ... <e_n>        |tgSlp| bin edges; values at the bin centres
///   phase   <n>                        S-curve: n phase nodes phase_j = -0.5 + j / n (j = 0 ... n-1); with n = 64 the nodes are the
///                                      1/64-pad grid of ClusterNative, so a snapped cluster takes the value of its own grid point
///   dist    <n> <dLo> <dHi>            edge part: n bins of d in [dLo, dHi) [pads]
///   S <stack> <iDrift> <iSlope> <class> <v_0> ... <v_{nPhase-1}>   S-curve: mean inward bias [um] at the phase nodes
///   E <stack> <iDrift> <iSlope> <class> <v_0> ... <v_{nDist-1}>    edge part: mean inward bias [um] per d bin, S-curve subtracted
/// stack: 0 IROC, 1 OROC1, 2 OROC2, 3 OROC3. class: 0 golden (no flag of 0x19), 1 split pad / shared (0x1 or 0x10 without 0x8),
/// 2 single pad or time (0x8); a class without records falls back to class 0. The S-curve part is taken everywhere, the edge
/// part only for d inside [dLo, dHi).
///
/// Interpolation: periodic and linear in the phase between the nodes; linear in the drift length and in |tgSlp| between the
/// bin centres, clamped to the first / last centre outside; linear in d between the bin centres, clamped below the first centre,
/// zero at and beyond dHi.

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

/// class of a cluster from its flags (UnbinnedResid::channel); a negative value (no flags stored) counts as golden
inline int classOfFlags(int flags)
{
  if (flags < 0) {
    return 0;
  }
  if (flags & FlagSingle) {
    return 2;
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
  std::vector<double> scurve;            ///< [stack][drift][slope][class][phase] in cm
  std::vector<double> edge;              ///< [stack][drift][slope][class][d] in cm
  std::vector<char> hasScurve;           ///< [stack][drift][slope][class]: record present
  std::vector<char> hasEdge;             ///< [stack][drift][slope][class]: record present
  std::array<double, NRows> padWidth{};  ///< [cm] per pad row (setRowGeometry)
  std::array<double, NRows> halfPlane{}; ///< nPads / 2 per pad row: the pad plane's edge in pads from y = 0
  bool loaded = false;

  /// pad geometry of a row (from the Mapper in the macro): pad width [cm] and number of pads
  void setRowGeometry(int row, double width, int nPads)
  {
    if (row >= 0 && row < NRows) {
      padWidth[row] = width;
      halfPlane[row] = 0.5 * nPads;
    }
  }

  bool isActive() const { return loaded; }
  int nDrift() const { return int(driftEdges.size()) - 1; }
  int nSlope() const { return int(slopeEdges.size()) - 1; }

  /// index of a (stack, drift bin, slope bin, class) record
  size_t record(int stack, int iDrift, int iSlope, int cls) const
  {
    return ((size_t(stack) * nDrift() + iDrift) * nSlope() + iSlope) * NClasses + cls;
  }

  /// read the table; false (and inactive) on any format error
  bool load(const std::string& fileName, std::string* error = nullptr)
  {
    loaded = false;
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
      } else if (key == "phase") {
        ss >> nPhase;
      } else if (key == "dist") {
        ss >> nDist >> distLo >> distHi;
      } else if (key == "S" || key == "E") {
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
    if (version != 1 || driftEdges.size() < 2 || slopeEdges.size() < 2 || nPhase < 2 || nDist < 0 || (nDist > 0 && distHi <= distLo)) {
      return fail(error, "missing or bad header (version / drift / slope / phase / dist)");
    }
    const size_t nRec = size_t(NStacks) * nDrift() * nSlope() * NClasses;
    scurve.assign(nRec * nPhase, 0.);
    edge.assign(nRec * std::max(nDist, 1), 0.);
    hasScurve.assign(nRec, 0);
    hasEdge.assign(nRec, 0);
    for (const auto& [type, v] : records) {
      const int nVal = type == 'S' ? nPhase : nDist;
      if (int(v.size()) != 4 + nVal) {
        return fail(error, std::string("record ") + type + ": wrong number of values");
      }
      const int stack = int(v[0]);
      const int iDrift = int(v[1]);
      const int iSlope = int(v[2]);
      const int cls = int(v[3]);
      if (stack < 0 || stack >= NStacks || iDrift < 0 || iDrift >= nDrift() || iSlope < 0 || iSlope >= nSlope() || cls < 0 || cls >= NClasses) {
        return fail(error, std::string("record ") + type + ": index out of range");
      }
      const size_t r = record(stack, iDrift, iSlope, cls);
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
    const int stack = stackOfRow(row);
    const int cls = classOfFlags(flags);
    const double bias = interpolate(stack, cls, drift, std::abs(tSlope), [&](size_t r, bool& ok) { return scurveAt(r, phase, ok); }) +
                        interpolate(stack, cls, drift, std::abs(tSlope), [&](size_t r, bool& ok) { return edgeAt(r, d, ok); });
    return float(side * bias);
  }

  /// inward bias [cm] (S-curve + edge part) at a distance d [pads] from the nearest edge: for tests and QA
  double inwardBias(int stack, int cls, double drift, double absSlope, double d) const
  {
    const double pad = d - 0.5;
    const double phase = pad - std::round(pad);
    return interpolate(stack, cls, drift, absSlope, [&](size_t r, bool& ok) { return scurveAt(r, phase, ok); }) +
           interpolate(stack, cls, drift, absSlope, [&](size_t r, bool& ok) { return edgeAt(r, d, ok); });
  }

 private:
  static bool fail(std::string* error, const std::string& message)
  {
    if (error) {
      *error = message;
    }
    return false;
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

  /// bilinear in (drift, |slope|) over the records of (stack, class), falling back to class 0 where a record is missing
  template <typename F>
  double interpolate(int stack, int cls, double drift, double absSlope, F&& value) const
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
        double v = value(record(stack, iD + a, iS + b, cls), ok);
        if (!ok && cls != 0) {
          v = value(record(stack, iD + a, iS + b, 0), ok);
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
