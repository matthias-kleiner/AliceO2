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

/// \file makeClusterBiasTable.C
/// \brief Producer of the ClusterBiasCorrection.h table (per-cluster correction of the pad-phase S-curve and the sector-edge
/// truncation) from CSV tables of mean inward cluster residuals (plain or gzip-compressed, read through gzip -dc).
///
/// - S-curve, from DATA (residuals to the track, clusters snapped to the 1/64-pad grid, far from the edge), one of two formats:
///   v2: region,stack,pad_width_cm,drift_lo_cm,drift_hi_cm,class,tOut_lo,tOut_hi,phase_node,count,mean_rIn_um,rms_rIn_um,
///       count_trim3mm,mean_trim3mm_um
///       per pad region (0-9) x drift x class x SIGNED slope tOut = tgSlp * sign(y) x phase node k / 64 (k = -32 ... 31)
///       => a table per pad region, the 64 nodes copied 1:1;
///   v1: stack,drift_lo_cm,drift_hi_cm,class,abs_tgSlp,phase_lo,phase_hi,count,mean_rIn_um,rms_rIn_um
///       per stack x drift x class x |tgSlp| (lt0.1 / 0.1-0.3 / ge0.3) x phase bin => a table per stack, interpolated onto the nodes.
///   The phase is p - round(p) of the pad coordinate p counted from the nearest edge; class golden / splitshared / single (flags as
///   in ClusterBiasCorrection.h). Only the ODD part [b(phase) - b(-phase)] / 2 of each profile is kept: it removes the smooth
///   distortion (a constant over the phase) and the even part, which is not a cluster bias but an artefact of the residual stream's
///   reference (odd in the signed slope, present in MC residuals, absent against the MC truth). The odd part is even in the slope,
///   so the signed slope bins of v2 are folded into |t| bins after the anti-symmetrisation, weighted by the inverse variance of the
///   odd part (node error = rms of the input profile / sqrt(count), the rms dominated by the distortion spread). Single-pad clusters
///   sit on the pad centre, where the odd part is 0: no record for them (class 0 is used, also 0 there). A record with fewer than
///   minCount clusters, an empty node, or an rms node error above maxNodeError (e.g. the innermost IROC rows at short drift: few
///   tracks, cm distortions) is treated as empty and takes the nearest accepted |t| bin, then (pad regions) the nearest region of
///   the same stack, then the nearest drift bin.
///   The plain means are used: the residuals contain the distortion (up to cm), so a truncation around 0 (the v2 trimmed columns,
///   useTrimmedMean) biases the profile wherever the distortion exceeds the cut.
/// - edge part, from a SIMULATION with the truth as reference (no distortion), one or several files (comma separated, e.g. two slope
///   arms; their sums add), in one of two formats:
///   full-MC truth (exportEdge.C; recognised by its mean_d_true column), the recommended source:
///       stack,drift_lo_cm,drift_hi_cm,class,abs_t,d_meas_lo,d_meas_hi,count,mean_dyIn_um,rms_dyIn_um,mean_d_true
///       class single (sigmaPad 0) / split (multi-pad, flag 0x1) / normal (multi-pad, not split) -> table class 2 / 1 / 0;
///       abs_t lt0.1 / ge0.1; drift in 50 cm bins: each table drift bin gets the profile linearly interpolated between the bin
///       centres (25, 75, ... cm; constant beyond the first / last one); d_meas 0.45-5.00 (use sMin = 2: the S-curve from d 2-5);
///   toy (toyEdgeTable.C):
///       stack,drift_lo_cm,drift_hi_cm,class,abs_t,d_meas_lo,d_meas_hi,count,mean_dyIn_um,rms_dyIn_um
///       class golden / split / singletime / single (singletime clusters spread over several pads and go to golden, as
///       ClusterBiasCorrection.h treats flag 0x8 off the pad centre); abs_t lt0.1 / 0.1-0.3 / ge0.3; a row fills every table drift
///       bin whose centre lies in its range.
///   d_meas = measured distance from the pad-plane edge in pads (pad-0 centre at 0.5, snapped). The simulation's own S-curve
///   (d_meas >= sMin, at the same phase) is subtracted, so that only the truncation is kept; beyond distHi the edge part is 0. d bins
///   with fewer than minBinCount clusters (or an S-curve phase bin below it) take the linear interpolation of their accepted
///   neighbours. Single-pad clusters (class 2): only the last pad (d < singleEdgeMax), nothing subtracted (in the bulk the single-pad
///   class of a simulation is dominated by noise clusters). Every pad region of a stack gets the stack's edge part; table bins still
///   empty take the nearest filled |t| bin, then the nearest filled drift bin, of the same unit and class.
///
/// smooth = true (default, v2 input): instead of the binned profiles, SMOOTH FITS (penalised B-spline least squares; the numerics in
/// ClusterBiasSmoothFit.h):
///   S per region x class (golden, splitshared) x |t| bin: penalised least squares (lambda 1e-4, second differences) of the odd part of
///     the v2 nodes (count >= 20, error rms / sqrt(n) (+) 5 um) on cubic B-splines in |u| (knots 0 .02 .05 .1 .15 .2 .3 .4 .5, the
///     value at 0+ free, S(0.5) = 0) x cubic B-splines in sqrt(drift) (knots at 2.5 10 25 50 100 160 237.5 cm), held at its 7.5 cm
///     value below 7.5 cm (the 0-5 cm bin is sparse); fewer than 40 points: the next lower |t| bin of the class, then golden;
///   E per stack x class x |t| arm of the simulation: the d bins with >= minBinCount clusters minus the simulation's S at the same
///     phase (d >= sMin), cubic B-splines in d (knots 0.5 ... 2.4, E = 0 from 2.4 on, constant below 0.5) x linear interpolation
///     between the simulation's drift-bin centres; single-pad clusters: the last-pad value (d <= 0.6) per drift centre, for d <
///     singleEdgeMax. Checked against an independent implementation of the same fit: equal to the printed 0.1 um.
/// The binned path (smooth = false, or a v1 input): nHarmonics > 0 replaces each S-curve by its sine series up to that order
/// (low-statistics inputs); 0 (default) keeps the profile. edgeScale scales the edge part (a simulated cloud known to differ from the
/// data).
///   root -l -b -q 'makeClusterBiasTable.C+("scurve_v2.csv.gz", "edge_truth_fullmc.csv", "clusterBias.txt", 4.0, 2.0)'   (smooth)

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "Framework/Logger.h"
#endif

#include "ClusterBiasSmoothFit.h"

namespace cluster_bias_table
{

constexpr double TwoPi = 6.283185307179586;
constexpr int NPhase = 64;
constexpr int NClasses = 3;
constexpr int NSelfPhase = 20;
constexpr double DistLo = 0.40;
constexpr double DistW = 0.05;

/// the lines of a plain or gzip-compressed (.gz) text file
std::vector<std::string> readLines(const std::string& fileName)
{
  std::vector<std::string> lines;
  const bool gz = fileName.size() > 3 && fileName.compare(fileName.size() - 3, 3, ".gz") == 0;
  if (gz) {
    FILE* pipe = popen(("gzip -dc '" + fileName + "'").c_str(), "r");
    if (!pipe) {
      LOGP(fatal, "makeClusterBiasTable: cannot read {}", fileName);
    }
    std::string line;
    std::array<char, 4096> buffer{};
    while (std::fgets(buffer.data(), buffer.size(), pipe)) {
      line += buffer.data();
      if (!line.empty() && line.back() == '\n') {
        line.pop_back();
        lines.push_back(line);
        line.clear();
      }
    }
    if (!line.empty()) {
      lines.push_back(line);
    }
    pclose(pipe);
  } else {
    std::ifstream in(fileName);
    if (!in) {
      LOGP(fatal, "makeClusterBiasTable: cannot open {}", fileName);
    }
    std::string line;
    while (std::getline(in, line)) {
      lines.push_back(line);
    }
  }
  if (lines.empty()) {
    LOGP(fatal, "makeClusterBiasTable: {} is empty", fileName);
  }
  return lines;
}

std::vector<std::string> splitCSV(const std::string& line)
{
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, ',')) {
    fields.push_back(field);
  }
  return fields;
}

/// column index by name in a CSV header (-1 if absent)
int column(const std::vector<std::string>& header, const std::string& name)
{
  const auto it = std::find(header.begin(), header.end(), name);
  return it == header.end() ? -1 : int(it - header.begin());
}

int stackIndex(const std::string& name)
{
  const char* names[] = {"IROC", "OROC1", "OROC2", "OROC3"};
  for (int i = 0; i < 4; ++i) {
    if (name == names[i]) {
      return i;
    }
  }
  return -1;
}

int stackOfRegion(int region)
{
  return region < 4 ? 0 : (region < 6 ? 1 : (region < 8 ? 2 : 3));
}

/// table class: 0 golden (with the single-time clusters of a simulation), 1 split pad / shared, 2 single pad
int classIndex(const std::string& name)
{
  if (name == "golden" || name == "singletime" || name == "normal") {
    return 0;
  }
  if (name == "splitshared" || name == "split") {
    return 1;
  }
  if (name == "single") {
    return 2;
  }
  return -1;
}

/// |t| range of the v1 / simulation slope labels
bool slopeLabelRange(const std::string& label, double maxSlope, double& lo, double& hi)
{
  if (label == "lt0.1") {
    lo = 0.;
    hi = 0.1;
  } else if (label == "0.1-0.3") {
    lo = 0.1;
    hi = 0.3;
  } else if (label == "ge0.3") {
    lo = 0.3;
    hi = maxSlope;
  } else if (label == "ge0.1") {
    lo = 0.1;
    hi = maxSlope;
  } else {
    return false;
  }
  return true;
}

/// count-weighted profile accumulator over fixed bins
struct Profile {
  std::vector<double> sum;
  std::vector<double> sum2;
  std::vector<double> cnt;
  explicit Profile(int n = 0) : sum(n, 0.), sum2(n, 0.), cnt(n, 0.) {}
  void fill(int bin, double count, double mean, double rms = 0.)
  {
    sum[bin] += count * mean;
    sum2[bin] += count * (rms * rms + mean * mean);
    cnt[bin] += count;
  }
  double value(int bin) const { return cnt[bin] > 0. ? sum[bin] / cnt[bin] : 0.; }
  double total() const
  {
    double c = 0.;
    for (const double v : cnt) {
      c += v;
    }
    return c;
  }
  /// variance of the entries pooled over all bins
  double variance() const
  {
    double n = 0.;
    double s = 0.;
    double s2 = 0.;
    for (size_t i = 0; i < cnt.size(); ++i) {
      n += cnt[i];
      s += sum[i];
      s2 += sum2[i];
    }
    return n > 0. ? std::max(s2 / n - (s / n) * (s / n), 0.) : 0.;
  }
};

/// error-weighted odd profile on the 64 nodes, summed over input profiles
struct OddAccumulator {
  std::vector<double> sumWV = std::vector<double>(NPhase, 0.);
  std::vector<double> sumW = std::vector<double>(NPhase, 0.);
  double count = 0.;
};

/// least-squares sine series up to order k of an odd periodic profile on nodes -0.5 + j / n
std::vector<double> sineSmooth(const std::vector<double>& v, int nHarmonics)
{
  const int n = int(v.size());
  std::vector<double> out(n, 0.);
  for (int k = 1; k <= nHarmonics && 2 * k < n; ++k) {
    double b = 0.;
    for (int j = 0; j < n; ++j) {
      b += 2. / n * v[j] * std::sin(TwoPi * k * (-0.5 + double(j) / n));
    }
    for (int j = 0; j < n; ++j) {
      out[j] += b * std::sin(TwoPi * k * (-0.5 + double(j) / n));
    }
  }
  return out;
}

int findBin(const std::vector<double>& edges, double x)
{
  const int i = int(std::upper_bound(edges.begin(), edges.end(), x) - edges.begin()) - 1;
  return (i < 0 || i >= int(edges.size()) - 1) ? -1 : i;
}

using RecordKey = std::tuple<int, int, int, int>; // unit, drift bin, |t| bin, class

/// records still empty take the nearest filled |t| bin of the same unit and drift bin, then (byRegion) the nearest filled pad region
/// of the same stack, then the nearest filled drift bin (same class throughout)
void fillEmpty(std::map<RecordKey, std::vector<double>>& records, int nUnits, int nDrift, int nSlope, bool byRegion)
{
  const int nPass = 3;
  for (int pass = 0; pass < nPass; ++pass) {
    if (pass == 1 && !byRegion) {
      continue;
    }
    // look-ups only in the records present at the start of the pass: the result does not depend on the loop order
    const auto filled = records;
    for (int unit = 0; unit < nUnits; ++unit) {
      for (int cls = 0; cls < NClasses; ++cls) {
        for (int iD = 0; iD < nDrift; ++iD) {
          for (int iS = 0; iS < nSlope; ++iS) {
            if (records.count({unit, iD, iS, cls})) {
              continue;
            }
            const int maxDist = pass == 0 ? nSlope : (pass == 1 ? nUnits : nDrift);
            for (int dist = 1; dist < maxDist; ++dist) {
              bool found = false;
              for (const int sgn : {-1, 1}) {
                const int jS = pass == 0 ? iS + sgn * dist : iS;
                const int jU = pass == 1 ? unit + sgn * dist : unit;
                const int jD = pass == 2 ? iD + sgn * dist : iD;
                if (pass == 1 && (jU < 0 || jU >= nUnits || stackOfRegion(jU) != stackOfRegion(unit))) {
                  continue;
                }
                const auto it = filled.find({jU, jD, jS, cls});
                if (it != filled.end()) {
                  records[{unit, iD, iS, cls}] = it->second;
                  found = true;
                  break;
                }
              }
              if (found) {
                break;
              }
            }
          }
        }
      }
    }
  }
}

// ===| smooth table (smooth = true) |====================================================================================================
// S per region x class (golden, splitshared) x |t| bin: penalised least squares of the odd part of the v2 nodes on cubic B-splines in
// |u| (dense near the pad centre, free value at 0+, S(0.5) = 0) x cubic B-splines in sqrt(drift), held constant below SMinDrift;
// E per stack x class x |t| arm of the simulation: B-splines in d on [0.5, 2.4] (E = 0 from 2.4 on, constant below 0.5) x linear
// interpolation between the drift-bin centres; single-pad clusters: the last-pad value per drift-bin centre, for d < 1.
const std::vector<double> SmoothUKnots = {0., 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5};
const std::vector<double> SmoothDriftKnots = {2.5, 10., 25., 50., 100., 160., 237.5}; // used as sqrt(drift)
const std::vector<double> SmoothDKnots = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.7, 2.0, 2.4};
constexpr double SmoothLambda = 1e-4;
constexpr double SMinDrift = 7.5;    // cm: S held at its value here below (the 0-5 cm bin is sparse: steep tracks at the readout)
constexpr double SMinNodeCount = 20; // v2 nodes with fewer clusters are not used
constexpr int SMinPoints = 40;       // fewer S points: fall back to the next lower |t| bin, then golden
constexpr double SErrFloor = 5.;     // um, added in quadrature to the S points
constexpr double EErrFloor = 10.;    // um, added in quadrature to the E points
constexpr double ESingleMaxD = 0.6;  // single-pad E points: the last pad only

/// |t| bin of a signed tOut bin [lo, hi) (the smaller |edge|, 0 for a bin across 0; the last bin above)
int absSlopeBin(double lo, double hi, const std::vector<double>& slopeEdges)
{
  const double a = lo * hi >= 0. ? std::min(std::abs(lo), std::abs(hi)) : 0.;
  for (size_t j = 0; j + 1 < slopeEdges.size(); ++j) {
    if (slopeEdges[j] - 1e-9 <= a && a < slopeEdges[j + 1] - 1e-9) {
      return int(j);
    }
  }
  return int(slopeEdges.size()) - 2;
}

void makeSmoothTable(const std::vector<std::string>& lines, const std::vector<std::string>& header, const std::vector<double>& driftEdges,
                     const std::vector<double>& slopeEdges, const std::string& scurveCSV, const std::string& edgeCSVs, const std::string& outFile,
                     double distHi, double sMin, double edgeScale, double maxSlope, double singleEdgeMax, double minBinCount,
                     const std::string& parameterLine)
{
  namespace sm = o2::tpc::cluster_bias::smooth;
  const int nDrift = int(driftEdges.size()) - 1;
  const int nSlope = int(slopeEdges.size()) - 1;
  const int nDist = int(std::round((distHi - DistLo) / DistW));
  std::vector<double> driftCentre(nDrift);
  for (int i = 0; i < nDrift; ++i) {
    driftCentre[i] = 0.5 * (driftEdges[i] + driftEdges[i + 1]);
  }
  std::vector<double> xKnots;
  for (const double d : SmoothDriftKnots) {
    xKnots.push_back(std::sqrt(d));
  }
  // ---- S points: per (region, class, drift bin, signed slope bin) the nodes k -> (mean, error)
  const int cUnit = column(header, "region");
  const int cDriftLo = column(header, "drift_lo_cm");
  const int cDriftHi = column(header, "drift_hi_cm");
  const int cClass = column(header, "class");
  const int cSlopeLo = column(header, "tOut_lo");
  const int cSlopeHi = column(header, "tOut_hi");
  const int cPhase = column(header, "phase_node");
  const int cCount = column(header, "count");
  const int cMean = column(header, "mean_rIn_um");
  const int cRms = column(header, "rms_rIn_um");
  std::map<std::tuple<int, int, int, double, double>, std::map<int, std::pair<double, double>>> raw;
  for (size_t i = 1; i < lines.size(); ++i) {
    const auto f = splitCSV(lines[i]);
    if (f.size() < header.size()) {
      continue;
    }
    const double n = std::stod(f[cCount]);
    const int cls = classIndex(f[cClass]);
    if (n < SMinNodeCount || cls < 0 || cls == 2) {
      continue;
    }
    const int iD = findBin(driftEdges, 0.5 * (std::stod(f[cDriftLo]) + std::stod(f[cDriftHi])));
    if (iD < 0) {
      continue;
    }
    const int k = int(std::lround(std::stod(f[cPhase]) * NPhase));
    raw[{std::stoi(f[cUnit]), cls, iD, std::stod(f[cSlopeLo]), std::stod(f[cSlopeHi])}][k] = {std::stod(f[cMean]), std::stod(f[cRms]) / std::sqrt(n)};
  }
  // odd part, signed slope folded into |t| (inverse variance)
  std::map<std::tuple<int, int, int, int>, std::map<int, std::pair<double, double>>> odd; // (region, class, |t| bin, drift bin) -> k -> (sum w y, sum w)
  for (const auto& [key, prof] : raw) {
    const auto [reg, cls, iD, tLo, tHi] = key;
    const int j = absSlopeBin(tLo, tHi, slopeEdges);
    for (int k = 1; k < NPhase / 2; ++k) {
      const auto ip = prof.find(k);
      const auto im = prof.find(-k);
      if (ip == prof.end() || im == prof.end()) {
        continue;
      }
      const double y = 0.5 * (ip->second.first - im->second.first);
      const double e = 0.5 * std::hypot(ip->second.second, im->second.second);
      const double w = 1. / (std::max(e, 1.) * std::max(e, 1.));
      auto& acc = odd[{reg, cls, j, iD}][k];
      acc.first += w * y;
      acc.second += w;
    }
  }
  std::map<std::tuple<int, int, int>, sm::SCurveFit> sfit;
  int nSFits = 0;
  for (int reg = 0; reg < 10; ++reg) {
    for (int cls = 0; cls < 2; ++cls) {
      for (int j = 0; j < nSlope; ++j) {
        std::vector<sm::SCurvePoint> pts;
        for (int iD = 0; iD < nDrift; ++iD) {
          const auto it = odd.find({reg, cls, j, iD});
          if (it == odd.end()) {
            continue;
          }
          for (const auto& [k, acc] : it->second) {
            pts.push_back({double(k) / NPhase, std::sqrt(driftCentre[iD]), acc.first / acc.second, std::sqrt(1. / acc.second + SErrFloor * SErrFloor)});
          }
        }
        if (int(pts.size()) < SMinPoints) {
          continue;
        }
        auto f = sm::fitSCurve(pts, SmoothUKnots, xKnots, SmoothLambda);
        if (f.valid()) {
          sfit[{reg, cls, j}] = f;
          ++nSFits;
        }
      }
    }
  }
  // ---- E points: per (stack, class, |t| arm, drift bin) the d bins (count >= minBinCount) -> (mean, error, count)
  struct EdgeBin {
    double dLo = 0.;
    double dHi = 0.;
    double n = 0.;
    double sum = 0.;
    double sumRms2 = 0.;
  };
  std::map<std::tuple<int, int, double, double, double, double>, std::map<int, EdgeBin>> erows;
  std::map<int, std::set<std::pair<double, double>>> stackDrifts;
  if (!edgeCSVs.empty()) {
    std::stringstream list(edgeCSVs);
    std::string name;
    while (std::getline(list, name, ',')) {
      const auto elines = readLines(name);
      const auto eh = splitCSV(elines[0]);
      const int cS = column(eh, "stack");
      const int cDL = column(eh, "drift_lo_cm");
      const int cDH = column(eh, "drift_hi_cm");
      const int cC = column(eh, "class");
      const int cT = column(eh, "abs_t");
      const int cXL = column(eh, "d_meas_lo");
      const int cXH = column(eh, "d_meas_hi");
      const int cN = column(eh, "count");
      const int cM = column(eh, "mean_dyIn_um");
      const int cR = column(eh, "rms_dyIn_um");
      if (std::min({cS, cDL, cDH, cC, cT, cXL, cXH, cN, cM, cR}) < 0) {
        LOGP(fatal, "makeClusterBiasTable: {} is not an edge table", name);
      }
      for (size_t i = 1; i < elines.size(); ++i) {
        const auto f = splitCSV(elines[i]);
        if (f.size() < eh.size()) {
          continue;
        }
        const int stack = stackIndex(f[cS]);
        const int cls = classIndex(f[cC]);
        double tLo = 0.;
        double tHi = 0.;
        const double n = std::stod(f[cN]);
        if (stack < 0 || cls < 0 || !slopeLabelRange(f[cT], maxSlope, tLo, tHi) || n < minBinCount) {
          continue;
        }
        const double dl = std::stod(f[cDL]);
        const double dh = std::stod(f[cDH]);
        stackDrifts[stack].insert({dl, dh});
        auto& b = erows[{stack, cls, tLo, tHi, dl, dh}][int(std::lround(std::stod(f[cXL]) * 100.))];
        const double rms = std::stod(f[cR]);
        b.dLo = std::stod(f[cXL]);
        b.dHi = std::stod(f[cXH]);
        b.n += n;
        b.sum += n * std::stod(f[cM]);
        b.sumRms2 += n * rms * rms;
      }
    }
  }
  // E = mean - the simulation's own S at the same phase (d >= sMin), fit per (stack, class, arm)
  std::map<std::tuple<int, int, double, double>, std::vector<sm::EdgePoint>> epts;
  for (const auto& [key, prof] : erows) {
    const auto [stack, cls, tLo, tHi, dl, dh] = key;
    for (const auto& [dk, b] : prof) {
      if (b.dLo >= SmoothDKnots.back()) {
        continue;
      }
      double sw = 0.;
      double sy = 0.;
      for (const auto& [dk2, b2] : prof) {
        if (b2.dLo >= sMin && (dk2 - dk) % 100 == 0) {
          sw += b2.n;
          sy += b2.sum;
        }
      }
      const double d = 0.5 * (b.dLo + b.dHi);
      if (sw <= 0. || d < SmoothDKnots.front() || (cls == 2 && d > ESingleMaxD)) {
        continue;
      }
      const double e = std::sqrt(b.sumRms2) / b.n;
      epts[{stack, cls, tLo, tHi}].push_back({d, 0.5 * (dl + dh), b.sum / b.n - sy / sw, std::sqrt(e * e + EErrFloor * EErrFloor)});
    }
  }
  std::map<std::tuple<int, int, double, double>, sm::EdgeFit> efit;
  for (const auto& [key, pts] : epts) {
    const auto [stack, cls, tLo, tHi] = key;
    if (pts.size() < 5) {
      continue;
    }
    std::vector<double> centres;
    for (const auto& [dl, dh] : stackDrifts[stack]) {
      centres.push_back(0.5 * (dl + dh));
    }
    auto f = sm::fitEdge(pts, centres, SmoothDKnots, SmoothLambda, cls == 2);
    if (f.valid()) {
      efit[key] = f;
    }
  }
  // ---- table
  FILE* out = std::fopen(outFile.c_str(), "w");
  if (!out) {
    LOGP(fatal, "makeClusterBiasTable: cannot write {}", outFile);
  }
  std::fprintf(out,
               "# ClusterBiasCorrection table from SMOOTH FITS (makeClusterBiasTable.C, smooth): S = odd part of %s "
               "(B-spline |u| x sqrt(drift)), E = %s minus its S (B-spline d x linear drift interpolation, x %.3f); inward bias [um]\n",
               scurveCSV.c_str(), edgeCSVs.empty() ? "-" : edgeCSVs.c_str(), edgeScale);
  std::fprintf(out, "%s", parameterLine.c_str());
  std::fprintf(out, "version 1\nunit region\ndrift %d", nDrift);
  for (const double e : driftEdges) {
    std::fprintf(out, " %g", e);
  }
  std::fprintf(out, "\nslope %d", nSlope);
  for (const double e : slopeEdges) {
    std::fprintf(out, " %g", e);
  }
  std::fprintf(out, "\nphase %d\ndist %d %g %g\n", NPhase, efit.empty() ? 0 : nDist, DistLo, distHi);
  int nS = 0;
  int nFallback = 0;
  int nE = 0;
  for (int reg = 0; reg < 10; ++reg) {
    const int stack = stackOfRegion(reg);
    for (int iD = 0; iD < nDrift; ++iD) {
      const double xS = std::sqrt(std::max(driftCentre[iD], SMinDrift));
      for (int j = 0; j < nSlope; ++j) {
        for (int cls = 0; cls < 2; ++cls) {
          auto it = sfit.find({reg, cls, j});
          for (int jj = j - 1; it == sfit.end() && jj >= 0; --jj) { // the nearest lower |t| bin, then golden
            it = sfit.find({reg, cls, jj});
          }
          if (it == sfit.end()) {
            it = sfit.find({reg, 0, j});
          }
          if (it == sfit.end()) {
            continue;
          }
          nFallback += (std::get<1>(it->first) != cls || std::get<2>(it->first) != j);
          std::fprintf(out, "S %d %d %d %d", reg, iD, j, cls);
          for (int n = 0; n < NPhase; ++n) {
            const double u = double(n - NPhase / 2) / NPhase;
            const double v = (n == NPhase / 2 || std::abs(u) >= 0.5) ? 0. : (u > 0. ? 1. : -1.) * sm::evalSCurve(it->second, SmoothUKnots, xKnots, std::abs(u), xS);
            std::fprintf(out, " %.1f", v);
          }
          std::fprintf(out, "\n");
          ++nS;
        }
        std::fprintf(out, "S %d %d %d 2", reg, iD, j); // single-pad clusters: on the pad centre, odd part 0
        for (int n = 0; n < NPhase; ++n) {
          std::fprintf(out, " 0");
        }
        std::fprintf(out, "\n");
        ++nS;
        const double tc = 0.5 * (slopeEdges[j] + slopeEdges[j + 1]);
        for (int cls = 0; cls < NClasses; ++cls) {
          for (const auto& [key, f] : efit) {
            const auto [st, c, tLo, tHi] = key;
            if (st != stack || c != cls || tc < tLo || tc >= tHi) {
              continue;
            }
            std::fprintf(out, "E %d %d %d %d", reg, iD, j, cls);
            for (int i = 0; i < nDist; ++i) {
              std::fprintf(out, " %.1f", edgeScale * sm::evalEdge(f, SmoothDKnots, DistLo + (i + 0.5) * DistW, driftCentre[iD], singleEdgeMax));
            }
            std::fprintf(out, "\n");
            ++nE;
            break;
          }
        }
      }
    }
  }
  std::fclose(out);
  LOGP(info, "makeClusterBiasTable (smooth): {} S fits, {} S records ({} from a fallback fit), {} E fits, {} E records written to {}", nSFits, nS,
       nFallback, efit.size(), nE, outFile);
}

} // namespace cluster_bias_table

void makeClusterBiasTable(std::string scurveCSV = "scurve_v2.csv.gz", // data S-curve table (v2 or v1, see above)
                          std::string edgeCSVs = "edge.csv",          // simulated edge tables (comma separated); "" = no edge part
                          std::string outFile = "clusterBias.txt",    // the ClusterBiasCorrection.h table
                          double distHi = 4.0,                        // edge part for d_meas < distHi [pads]
                          double sMin = 3.5,                          // the edge table's own S-curve from d_meas >= sMin
                          double edgeScale = 1.,                      // scale of the edge part
                          int nHarmonics = 0,                         // > 0: sine-smooth the S-curves to this order
                          double minCount = 1000.,                    // records with fewer clusters are treated as empty
                          double maxNodeError = 100.,                 // S-curve records with a larger rms node error [um] too
                          double maxSlope = 1.,                       // upper edge of the last |tgSlp| bin
                          double singleEdgeMax = 1.,                  // single-pad clusters: edge part for d_meas < singleEdgeMax only
                          bool useTrimmedMean = false,                // v2: the |rIn| < 3 mm truncated means instead of the plain ones
                          double minBinCount = 100.,                  // edge: d bins with fewer clusters are interpolated (smooth: left out)
                          bool smooth = true)                         // smooth B-spline fits (v2 input) instead of the binned profiles
{
  using namespace cluster_bias_table;
  const int nDist = int(std::round((distHi - DistLo) / DistW));

  // ===| S-curve input |=========================================================================================================
  const auto lines = readLines(scurveCSV);
  const auto header = splitCSV(lines[0]);
  const bool v2 = column(header, "phase_node") >= 0;
  const int cUnit = column(header, v2 ? "region" : "stack");
  const int cDriftLo = column(header, "drift_lo_cm");
  const int cDriftHi = column(header, "drift_hi_cm");
  const int cClass = column(header, "class");
  const int cSlopeLo = column(header, v2 ? "tOut_lo" : "abs_tgSlp");
  const int cSlopeHi = v2 ? column(header, "tOut_hi") : cSlopeLo;
  const int cPhase = column(header, v2 ? "phase_node" : "phase_lo");
  const int cCount = column(header, (v2 && useTrimmedMean) ? "count_trim3mm" : "count");
  const int cMean = column(header, (v2 && useTrimmedMean) ? "mean_trim3mm_um" : "mean_rIn_um");
  const int cRms = column(header, "rms_rIn_um"); // (the truncated means: the full rms, a conservative error)
  if (std::min({cUnit, cDriftLo, cDriftHi, cClass, cSlopeLo, cSlopeHi, cPhase, cCount, cMean, cRms}) < 0) {
    LOGP(fatal, "makeClusterBiasTable: {} is neither a v1 nor a v2 S-curve table", scurveCSV);
  }
  const int nUnits = v2 ? 10 : 4;
  // binning of the output: the input's drift bins; |t| bins from the non-negative signed edges (v2) or the v1 labels
  std::set<double> driftSet;
  std::set<double> slopeSet{0.};
  std::set<double> phaseSet;
  for (size_t i = 1; i < lines.size(); ++i) {
    const auto f = splitCSV(lines[i]);
    if (f.size() < header.size()) {
      continue;
    }
    driftSet.insert(std::stod(f[cDriftLo]));
    driftSet.insert(std::stod(f[cDriftHi]));
    phaseSet.insert(std::stod(f[cPhase]));
    if (v2) {
      for (const int c : {cSlopeLo, cSlopeHi}) {
        const double t = std::stod(f[c]);
        if (t >= 0. && t < maxSlope) {
          slopeSet.insert(t);
        }
      }
    }
  }
  const std::vector<double> driftEdges(driftSet.begin(), driftSet.end());
  std::vector<double> slopeEdges = v2 ? std::vector<double>(slopeSet.begin(), slopeSet.end()) : std::vector<double>{0., 0.1, 0.3};
  slopeEdges.push_back(maxSlope);
  const int nDrift = int(driftEdges.size()) - 1;
  const int nSlope = int(slopeEdges.size()) - 1;
  const int nIn = v2 ? NPhase : int(phaseSet.size());
  if (nDrift < 1 || nIn < 2) {
    LOGP(fatal, "makeClusterBiasTable: {} has no usable binning", scurveCSV);
  }
  char parameterLine[512];
  std::snprintf(parameterLine, sizeof(parameterLine),
                "# makeClusterBiasTable parameters: distHi %g sMin %g edgeScale %g nHarmonics %d minCount %g maxNodeError %g maxSlope %g "
                "singleEdgeMax %g useTrimmedMean %d minBinCount %g smooth %d\n",
                distHi, sMin, edgeScale, nHarmonics, minCount, maxNodeError, maxSlope, singleEdgeMax, int(useTrimmedMean), minBinCount, int(smooth));
  if (smooth && v2) {
    makeSmoothTable(lines, header, driftEdges, slopeEdges, scurveCSV, edgeCSVs, outFile, distHi, sMin, edgeScale, maxSlope, singleEdgeMax,
                    minBinCount, parameterLine);
    return;
  }
  if (smooth) {
    LOGP(warning, "makeClusterBiasTable: the smooth fit needs the v2 S-curve table; {} is v1, the binned profiles are used", scurveCSV);
  }
  // per (unit, drift, class, signed slope bin as its lower edge) the profile over the input phase bins
  std::map<std::tuple<int, int, int, double, double>, Profile> sprof;
  for (size_t i = 1; i < lines.size(); ++i) {
    const auto f = splitCSV(lines[i]);
    if (f.size() < header.size()) {
      continue;
    }
    const int unit = v2 ? std::stoi(f[cUnit]) : stackIndex(f[cUnit]);
    const int cls = classIndex(f[cClass]);
    const int iD = findBin(driftEdges, 0.5 * (std::stod(f[cDriftLo]) + std::stod(f[cDriftHi])));
    double tLo = 0.;
    double tHi = 0.;
    if (v2) {
      tLo = std::stod(f[cSlopeLo]);
      tHi = std::stod(f[cSlopeHi]);
    } else if (!slopeLabelRange(f[cSlopeLo], maxSlope, tLo, tHi)) {
      continue;
    }
    if (unit < 0 || unit >= nUnits || cls < 0 || iD < 0) {
      continue;
    }
    const double phase = std::stod(f[cPhase]);
    const int bin = v2 ? ((int(std::lround(phase * NPhase)) + NPhase / 2) % NPhase + NPhase) % NPhase
                       : std::clamp(int(std::floor((phase + 0.5) * nIn + 0.5)), 0, nIn - 1);
    sprof.try_emplace({unit, iD, cls, tLo, tHi}, Profile(nIn)).first->second.fill(bin, std::stod(f[cCount]), std::stod(f[cMean]), std::stod(f[cRms]));
  }
  // anti-symmetrise each profile, fold the signed slope into |t| (inverse-variance weighted), map onto the 64 nodes
  std::map<RecordKey, OddAccumulator> sodd;
  for (const auto& [key, p] : sprof) {
    const auto [unit, iD, cls, tLo, tHi] = key;
    if (cls == 2) {
      continue; // single-pad clusters: on the pad centre, odd part 0
    }
    const double absLo = tLo >= 0. ? tLo : -tHi;
    const int iS = findBin(slopeEdges, std::min(absLo, maxSlope - 1e-9) + 1e-9);
    if (iS < 0) {
      continue;
    }
    auto& acc = sodd[{unit, iD, iS, cls}];
    acc.count += p.total();
    // odd part and its variance per input bin; 0 where the bin or its mirror is empty
    const double var = std::max(p.variance(), 1.); // [um^2]
    std::vector<double> odd(nIn, 0.);
    std::vector<double> err2(nIn, 0.);
    for (int j = 0; j < nIn; ++j) {
      const int mirror = v2 ? (nIn - j) % nIn : nIn - 1 - j; // v2: node k <-> -k (k = j - 32); v1: bin <-> mirrored bin
      if (p.cnt[j] > 0. && p.cnt[mirror] > 0.) {
        odd[j] = 0.5 * (p.value(j) - p.value(mirror));
        err2[j] = 0.25 * var * (1. / p.cnt[j] + 1. / p.cnt[mirror]);
      }
    }
    for (int j = 0; j < NPhase; ++j) {
      double value = 0.;
      double e2 = 0.;
      if (v2) {
        value = odd[j];
        e2 = err2[j];
      } else {
        const double x = (double(j) / NPhase) * nIn - 0.5; // node phase in units of input bin centres
        const double fl = std::floor(x);
        const double fr = x - fl;
        const int i0 = ((int(fl) % nIn) + nIn) % nIn;
        const int i1 = (i0 + 1) % nIn;
        value = (1. - fr) * odd[i0] + fr * odd[i1];
        e2 = (err2[i0] > 0. && err2[i1] > 0.) ? (1. - fr) * (1. - fr) * err2[i0] + fr * fr * err2[i1] : 0.;
      }
      if (e2 > 0.) {
        acc.sumWV[j] += value / e2;
        acc.sumW[j] += 1. / e2;
      }
    }
  }

  FILE* out = std::fopen(outFile.c_str(), "w");
  if (!out) {
    LOGP(fatal, "makeClusterBiasTable: cannot write {}", outFile);
  }
  std::fprintf(out, "# ClusterBiasCorrection table: S-curve = odd part of %s (%s means), edge part from %s (x %.3f); inward bias [um]\n",
               scurveCSV.c_str(), useTrimmedMean ? "truncated" : "plain", edgeCSVs.empty() ? "-" : edgeCSVs.c_str(), edgeScale);
  std::fprintf(out, "%s", parameterLine);
  std::fprintf(out, "version 1\nunit %s\ndrift %d", v2 ? "region" : "stack", nDrift);
  for (const double e : driftEdges) {
    std::fprintf(out, " %g", e);
  }
  std::fprintf(out, "\nslope %d", nSlope);
  for (const double e : slopeEdges) {
    std::fprintf(out, " %g", e);
  }
  std::fprintf(out, "\nphase %d\ndist %d %g %g\n", NPhase, edgeCSVs.empty() ? 0 : nDist, DistLo, distHi);
  std::map<RecordKey, std::vector<double>> srec;
  int nLowCount = 0;
  int nNoisy = 0;
  for (const auto& [key, acc] : sodd) {
    if (acc.count < minCount) {
      ++nLowCount;
      continue;
    }
    std::vector<double> nodes(NPhase, 0.);
    double sumErr2 = 0.;
    bool hole = false;
    for (int j = 0; j < NPhase; ++j) {
      if (acc.sumW[j] > 0.) {
        nodes[j] = acc.sumWV[j] / acc.sumW[j];
        sumErr2 += 1. / acc.sumW[j];
      } else {
        hole = true;
      }
    }
    if (hole || std::sqrt(sumErr2 / NPhase) > maxNodeError) {
      ++nNoisy;
      continue;
    }
    if (nHarmonics > 0) {
      nodes = sineSmooth(nodes, nHarmonics);
    }
    srec[key] = nodes;
  }
  const size_t nSMeasured = srec.size();
  LOGP(info, "makeClusterBiasTable: S-curve records: {} accepted, {} below minCount {}, {} with an empty node or an rms node error above {} um",
       nSMeasured, nLowCount, minCount, nNoisy, maxNodeError);
  fillEmpty(srec, nUnits, nDrift, nSlope, v2); // e.g. the innermost IROC rows at short drift
  int nS = 0;
  for (const auto& [key, nodes] : srec) {
    const auto [unit, iD, iS, cls] = key;
    std::fprintf(out, "S %d %d %d %d", unit, iD, iS, cls);
    for (const double v : nodes) {
      std::fprintf(out, " %.1f", v);
    }
    std::fprintf(out, "\n");
    ++nS;
  }

  // ===| edge part |=============================================================================================================
  int nE = 0;
  if (!edgeCSVs.empty()) {
    // per (stack, drift lo, drift hi, class, |t| lo, |t| hi): the d profile and the table's own S-curve at d >= sMin
    using EdgeKey = std::tuple<int, double, double, int, double, double>;
    std::map<EdgeKey, Profile> eprof;
    std::map<EdgeKey, Profile> selfS;
    std::set<EdgeKey> interpolateDrift; // rows of the full-MC format
    std::stringstream list(edgeCSVs);
    std::string name;
    while (std::getline(list, name, ',')) {
      const auto elines = readLines(name);
      const auto eh = splitCSV(elines[0]);
      const int cS = column(eh, "stack");
      const int cDL = column(eh, "drift_lo_cm");
      const int cDH = column(eh, "drift_hi_cm");
      const int cC = column(eh, "class");
      const int cT = column(eh, "abs_t");
      const int cXL = column(eh, "d_meas_lo");
      const int cXH = column(eh, "d_meas_hi");
      const int cN = column(eh, "count");
      const int cM = column(eh, "mean_dyIn_um");
      if (std::min({cS, cDL, cDH, cC, cT, cXL, cXH, cN, cM}) < 0) {
        LOGP(fatal, "makeClusterBiasTable: {} is not an edge table", name);
      }
      const bool fullMC = column(eh, "mean_d_true") >= 0;
      LOGP(info, "makeClusterBiasTable: edge part from {} ({} format)", name, fullMC ? "full-MC truth" : "toy");
      for (size_t i = 1; i < elines.size(); ++i) {
        const auto f = splitCSV(elines[i]);
        if (f.size() < eh.size()) {
          continue;
        }
        const int stack = stackIndex(f[cS]);
        const int cls = classIndex(f[cC]);
        double tLo = 0.;
        double tHi = 0.;
        if (stack < 0 || cls < 0 || !slopeLabelRange(f[cT], maxSlope, tLo, tHi)) {
          continue;
        }
        const EdgeKey key{stack, std::stod(f[cDL]), std::stod(f[cDH]), cls, tLo, tHi};
        if (fullMC) {
          interpolateDrift.insert(key);
        }
        const double d = 0.5 * (std::stod(f[cXL]) + std::stod(f[cXH]));
        const double count = std::stod(f[cN]);
        const double mean = std::stod(f[cM]);
        if (d >= sMin) {
          const double p = (d - 0.5) - std::round(d - 0.5);
          selfS.try_emplace(key, Profile(NSelfPhase)).first->second.fill(std::clamp(int(std::floor((p + 0.5) * NSelfPhase)), 0, NSelfPhase - 1), count, mean);
        }
        if (d >= DistLo && d < distHi) {
          eprof.try_emplace(key, Profile(nDist)).first->second.fill(std::clamp(int(std::floor((d - DistLo) / DistW)), 0, nDist - 1), count, mean);
        }
      }
    }
    // [unit][drift][slope][class] -> edge values
    std::map<RecordKey, std::vector<double>> edge;
    // full-MC rows: per (stack, class, |t| range) the profiles at their drift-bin centres, interpolated below
    std::map<std::tuple<int, int, double, double>, std::map<double, std::vector<double>>> byDrift;
    int nThinBins = 0;
    for (const auto& [key, e] : eprof) {
      const auto [stack, dLo, dHi, cls, tLo, tHi] = key;
      if (e.total() < minCount) {
        continue;
      }
      std::vector<double> v(nDist, 0.);
      if (cls == 2) {
        // single-pad clusters sit exactly on d = n + 0.5, the boundary of two d bins: both bins get their common value; only the
        // last pad (d < singleEdgeMax), no S-curve to subtract
        for (int n = 0; n + 0.5 < std::min(distHi, singleEdgeMax); ++n) {
          const int b = int(std::lround((n + 0.5 - DistLo) / DistW));
          if (b >= 1 && b < nDist) {
            const double c = e.cnt[b - 1] + e.cnt[b];
            v[b - 1] = v[b] = c > 0. ? (e.sum[b - 1] + e.sum[b]) / c : 0.;
          }
        }
      } else {
        const auto self = selfS.find(key);
        std::vector<bool> good(nDist, false);
        for (int i = 0; i < nDist; ++i) {
          v[i] = e.value(i);
          good[i] = e.cnt[i] >= minBinCount;
          if (self != selfS.end()) {
            // subtract the simulation's own S-curve at the phase of the bin
            const double d = DistLo + (i + 0.5) * DistW;
            const double p = (d - 0.5) - std::round(d - 0.5);
            const int iP = std::clamp(int(std::floor((p + 0.5) * NSelfPhase)), 0, NSelfPhase - 1);
            v[i] -= self->second.value(iP);
            good[i] = good[i] && self->second.cnt[iP] >= minBinCount;
          }
        }
        // thin d bins: linear interpolation of the nearest accepted neighbours (the nearest one beyond the last)
        std::vector<int> ok;
        for (int i = 0; i < nDist; ++i) {
          if (good[i]) {
            ok.push_back(i);
          }
        }
        if (ok.empty()) {
          continue;
        }
        for (int i = 0; i < nDist; ++i) {
          if (good[i]) {
            continue;
          }
          ++nThinBins;
          const auto hi = std::upper_bound(ok.begin(), ok.end(), i);
          if (hi == ok.begin()) {
            v[i] = v[ok.front()];
          } else if (hi == ok.end()) {
            v[i] = v[ok.back()];
          } else {
            const int i0 = *(hi - 1);
            const int i1 = *hi;
            v[i] = v[i0] + (v[i1] - v[i0]) * double(i - i0) / double(i1 - i0);
          }
        }
      }
      for (auto& x : v) {
        x *= edgeScale;
      }
      if (interpolateDrift.count(key)) {
        byDrift[{stack, cls, tLo, tHi}][0.5 * (dLo + dHi)] = v;
        continue;
      }
      // toy rows: every output unit of the stack, drift and |t| bin whose centre lies in the row's range
      for (int unit = 0; unit < nUnits; ++unit) {
        if ((v2 ? stackOfRegion(unit) : unit) != stack) {
          continue;
        }
        for (int iD = 0; iD < nDrift; ++iD) {
          const double c = 0.5 * (driftEdges[iD] + driftEdges[iD + 1]);
          if (c < dLo || c >= dHi) {
            continue;
          }
          for (int iS = 0; iS < nSlope; ++iS) {
            const double t = 0.5 * (slopeEdges[iS] + slopeEdges[iS + 1]);
            if (t >= tLo && t < tHi) {
              edge[{unit, iD, iS, cls}] = v;
            }
          }
        }
      }
    }
    // full-MC rows: linear in drift between the bin centres, constant beyond the first / last one
    for (const auto& [gkey, profiles] : byDrift) {
      const auto [stack, cls, tLo, tHi] = gkey;
      for (int iD = 0; iD < nDrift; ++iD) {
        const double c = 0.5 * (driftEdges[iD] + driftEdges[iD + 1]);
        std::vector<double> v;
        const auto hi = profiles.lower_bound(c);
        if (hi == profiles.begin()) {
          v = hi->second;
        } else if (hi == profiles.end()) {
          v = std::prev(hi)->second;
        } else {
          const auto lo = std::prev(hi);
          const double fr = (c - lo->first) / (hi->first - lo->first);
          v.resize(nDist);
          for (int i = 0; i < nDist; ++i) {
            v[i] = (1. - fr) * lo->second[i] + fr * hi->second[i];
          }
        }
        for (int unit = 0; unit < nUnits; ++unit) {
          if ((v2 ? stackOfRegion(unit) : unit) != stack) {
            continue;
          }
          for (int iS = 0; iS < nSlope; ++iS) {
            const double t = 0.5 * (slopeEdges[iS] + slopeEdges[iS + 1]);
            if (t >= tLo && t < tHi) {
              edge[{unit, iD, iS, cls}] = v;
            }
          }
        }
      }
    }
    LOGP(info, "makeClusterBiasTable: edge part: {} d bins below minBinCount {} interpolated from their neighbours", nThinBins, minBinCount);
    fillEmpty(edge, nUnits, nDrift, nSlope, false); // (an edge row fills all regions of its stack already)
    for (const auto& [key, v] : edge) {
      const auto [unit, iD, iS, cls] = key;
      std::fprintf(out, "E %d %d %d %d", unit, iD, iS, cls);
      for (const double x : v) {
        std::fprintf(out, " %.1f", x);
      }
      std::fprintf(out, "\n");
      ++nE;
    }
  }
  std::fclose(out);
  LOGP(info,
       "makeClusterBiasTable: {} S-curve ({} measured, the rest from the nearest bin) and {} edge records ({} x {} drift x {} slope bins, "
       "{} phase nodes, {} d bins) written to {}",
       nS, nSMeasured, nE, v2 ? "10 pad regions" : "4 stacks", nDrift, nSlope, NPhase, nDist, outFile);
}
