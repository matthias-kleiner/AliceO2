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
/// truncation) from two CSV tables of mean inward cluster residuals:
///
/// - S-curve, from DATA (residuals to the track, clusters snapped to the 1/64-pad grid, far from the edge):
///     stack,drift_lo_cm,drift_hi_cm,class,abs_tgSlp,phase_lo,phase_hi,count,mean_rIn_um,rms_rIn_um
///   stack IROC / OROC1 / OROC2 / OROC3; class golden / splitshared / single (flags as ClusterBiasCorrection.h); abs_tgSlp lt0.1 /
///   0.1-0.3 / ge0.3; phase = p - round(p) of the pad coordinate p counted from the nearest edge, bins [phase_lo, phase_hi).
///   Per stack x drift x slope the golden mean over the phase is the smooth distortion: it is subtracted from all classes.
/// - edge part, from a SIMULATION with the truth as reference (no distortion):
///     stack,drift_lo_cm,drift_hi_cm,class,abs_t,d_meas_lo,d_meas_hi,count,mean_dyIn_um,rms_dyIn_um
///   class golden / split / singletime / single (single = single PAD: the table's class 2; singletime clusters spread over several
///   pads and are merged into golden, as ClusterBiasCorrection.h treats flag 0x8 off the pad centre); d_meas = measured distance
///   from the pad-plane edge in pads (pad-0 centre at 0.5, snapped). The S-curve of the same table (d_meas >= sMin) is subtracted
///   per phase, so that only the truncation is kept; beyond distHi the edge part is 0. Single-pad clusters: only the last pad
///   (d < singleEdgeMax) is taken and nothing is subtracted (in the bulk the single-pad class of a simulation is dominated by
///   noise clusters). The edge rows are matched to the S-curve's drift / slope bins by the bin centres.
///
/// The S-curve is written on 64 phase nodes (the 1/64-pad grid): copied for an input with 64 bins (one grid point per bin,
/// at the bin's lower edge), otherwise interpolated periodically between the bin centres. Single-pad clusters sit on the pad
/// centre and have no S-curve by symmetry: their record is 0 unless singleScurveFromData (then the data's value at phase 0, minus
/// the golden mean, everywhere: the data's single class mixes single-time clusters and has long tails). nHarmonics > 0 replaces each S-curve by its
/// Fourier series up to that order (even and odd terms) for low-statistics inputs; 0 (default) keeps the measured profile.
///   root -l -b -q 'makeClusterBiasTable.C+("scurve_phase.csv", "edge.csv", "clusterBias.txt")'
/// The edge part can be scaled (edgeScale) when the simulated cloud is known to differ from the data.

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "Framework/Logger.h"
#endif

namespace cluster_bias_table
{

constexpr double TwoPi = 6.283185307179586;

struct Row {
  int stack = -1;
  double driftLo = 0.;
  double driftHi = 0.;
  int cls = -1;
  int slope = -1;
  double xLo = 0.;
  double xHi = 0.;
  double count = 0.;
  double mean = 0.;
};

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

int slopeIndex(const std::string& name)
{
  if (name == "lt0.1") {
    return 0;
  }
  if (name == "0.1-0.3") {
    return 1;
  }
  if (name == "ge0.3") {
    return 2;
  }
  return -1;
}

/// table class: 0 golden (with the single-time clusters of a simulation), 1 split pad / shared, 2 single pad
int classIndex(const std::string& name)
{
  if (name == "golden" || name == "singletime") {
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

std::vector<Row> readCSV(const std::string& fileName)
{
  std::vector<Row> rows;
  std::ifstream in(fileName);
  if (!in) {
    LOGP(fatal, "makeClusterBiasTable: cannot open {}", fileName);
  }
  std::string line;
  std::getline(in, line); // header
  while (std::getline(in, line)) {
    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> f;
    while (std::getline(ss, field, ',')) {
      f.push_back(field);
    }
    if (f.size() < 10) {
      continue;
    }
    Row r;
    r.stack = stackIndex(f[0]);
    r.driftLo = std::stod(f[1]);
    r.driftHi = std::stod(f[2]);
    r.cls = classIndex(f[3]);
    r.slope = slopeIndex(f[4]);
    r.xLo = std::stod(f[5]);
    r.xHi = std::stod(f[6]);
    r.count = std::stod(f[7]);
    r.mean = std::stod(f[8]);
    if (r.stack < 0 || r.cls < 0 || r.slope < 0) {
      LOGP(warning, "makeClusterBiasTable: unknown stack / class / slope in {}: {}", fileName, line);
      continue;
    }
    rows.push_back(r);
  }
  LOGP(info, "makeClusterBiasTable: {} rows from {}", rows.size(), fileName);
  return rows;
}

/// count-weighted profile accumulator over fixed bins
struct Profile {
  std::vector<double> sum;
  std::vector<double> cnt;
  explicit Profile(int n = 0) : sum(n, 0.), cnt(n, 0.) {}
  void fill(int bin, double count, double mean)
  {
    sum[bin] += count * mean;
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
  double mean() const
  {
    double s = 0.;
    double c = 0.;
    for (size_t i = 0; i < sum.size(); ++i) {
      s += sum[i];
      c += cnt[i];
    }
    return c > 0. ? s / c : 0.;
  }
};

/// least-squares Fourier series (constant, cos and sin up to order k) of a periodic profile on nodes -0.5 + j / n
std::vector<double> fourierSmooth(const std::vector<double>& v, int nHarmonics)
{
  const int n = int(v.size());
  std::vector<double> out(n, 0.);
  double c0 = 0.;
  for (const double x : v) {
    c0 += x / n;
  }
  for (int j = 0; j < n; ++j) {
    out[j] = c0;
  }
  for (int k = 1; k <= nHarmonics && 2 * k < n; ++k) {
    double a = 0.;
    double b = 0.;
    for (int j = 0; j < n; ++j) {
      const double phase = -0.5 + double(j) / n;
      a += 2. / n * v[j] * std::cos(TwoPi * k * phase);
      b += 2. / n * v[j] * std::sin(TwoPi * k * phase);
    }
    for (int j = 0; j < n; ++j) {
      const double phase = -0.5 + double(j) / n;
      out[j] += a * std::cos(TwoPi * k * phase) + b * std::sin(TwoPi * k * phase);
    }
  }
  return out;
}

} // namespace cluster_bias_table

void makeClusterBiasTable(std::string scurveCSV = "scurve_phase.csv", // data S-curve table (format above)
                          std::string edgeCSV = "edge.csv",           // simulated edge table (format above); "" = no edge part
                          std::string outFile = "clusterBias.txt",    // the ClusterBiasCorrection.h table
                          double distHi = 4.0,                        // edge part for d_meas < distHi [pads]
                          double sMin = 3.5,                          // the edge table's own S-curve from d_meas >= sMin
                          double edgeScale = 1.,                      // scale of the edge part
                          int nHarmonics = 0,                         // > 0: Fourier-smooth the S-curves to this order
                          double minCount = 1000.,                    // records with fewer clusters are not written (class 0 is used)
                          double maxSlope = 1.,                       // upper edge of the last |tgSlp| bin
                          double singleEdgeMax = 1.,                  // single-pad clusters: edge part for d_meas < singleEdgeMax only
                          bool singleScurveFromData = false)          // single-pad clusters: the data's value at phase 0 instead of 0
{
  using namespace cluster_bias_table;
  constexpr int NPhase = 64;
  constexpr double DistLo = 0.40;
  constexpr double DistW = 0.05;
  const int nDist = int(std::round((distHi - DistLo) / DistW));

  // ===| S-curve |===============================================================================================================
  const auto srows = readCSV(scurveCSV);
  std::set<double> driftSet;
  std::set<double> phaseLoSet;
  for (const auto& r : srows) {
    driftSet.insert(r.driftLo);
    driftSet.insert(r.driftHi);
    phaseLoSet.insert(r.xLo);
  }
  const std::vector<double> driftEdges(driftSet.begin(), driftSet.end());
  const std::vector<double> slopeEdges{0., 0.1, 0.3, maxSlope};
  const int nDrift = int(driftEdges.size()) - 1;
  const int nIn = int(phaseLoSet.size());
  if (nDrift < 1 || nIn < 2) {
    LOGP(fatal, "makeClusterBiasTable: {} has no usable binning", scurveCSV);
  }
  auto driftBin = [&](double lo) { return int(std::lower_bound(driftEdges.begin(), driftEdges.end(), lo) - driftEdges.begin()); };
  auto phaseBin = [&](double lo) { return std::clamp(int(std::floor((lo + 0.5) * nIn + 0.5)), 0, nIn - 1); };
  // [stack][drift][slope][class]
  std::map<std::tuple<int, int, int, int>, Profile> sprof;
  for (const auto& r : srows) {
    auto& p = sprof.try_emplace({r.stack, driftBin(r.driftLo), r.slope, r.cls}, Profile(nIn)).first->second;
    p.fill(phaseBin(r.xLo), r.count, r.mean);
  }
  FILE* out = std::fopen(outFile.c_str(), "w");
  if (!out) {
    LOGP(fatal, "makeClusterBiasTable: cannot write {}", outFile);
  }
  std::fprintf(out, "# ClusterBiasCorrection table from %s (S-curve) and %s (edge part, x %.3f), inward bias [um]\n", scurveCSV.c_str(),
               edgeCSV.empty() ? "-" : edgeCSV.c_str(), edgeScale);
  std::fprintf(out, "version 1\ndrift %d", nDrift);
  for (const double e : driftEdges) {
    std::fprintf(out, " %g", e);
  }
  std::fprintf(out, "\nslope 3");
  for (const double e : slopeEdges) {
    std::fprintf(out, " %g", e);
  }
  std::fprintf(out, "\nphase %d\ndist %d %g %g\n", NPhase, edgeCSV.empty() ? 0 : nDist, DistLo, distHi);
  int nS = 0;
  for (const auto& [key, p] : sprof) {
    const auto [stack, iD, iS, cls] = key;
    if (p.total() < minCount) {
      continue;
    }
    // the smooth distortion: the golden mean over the phase of this stack / drift / slope
    const auto golden = sprof.find({stack, iD, iS, 0});
    const double offset = golden != sprof.end() ? golden->second.mean() : p.mean();
    std::vector<double> nodes(NPhase, 0.);
    if (cls == 2) {
      // single-pad clusters sit on the pad centre: no S-curve by symmetry (or the data's value there)
      std::fill(nodes.begin(), nodes.end(), singleScurveFromData ? p.value(phaseBin(0.)) - offset : 0.);
    } else if (nIn == NPhase) {
      for (int j = 0; j < NPhase; ++j) {
        nodes[j] = p.value(j) - offset;
      }
    } else {
      for (int j = 0; j < NPhase; ++j) {
        const double x = (double(j) / NPhase) * nIn - 0.5; // node phase in units of input bin centres
        const double fl = std::floor(x);
        const double f = x - fl;
        const int i0 = ((int(fl) % nIn) + nIn) % nIn;
        const int i1 = (i0 + 1) % nIn;
        nodes[j] = (1. - f) * p.value(i0) + f * p.value(i1) - offset;
      }
    }
    if (nHarmonics > 0 && cls != 2) {
      nodes = fourierSmooth(nodes, nHarmonics);
    }
    std::fprintf(out, "S %d %d %d %d", stack, iD, iS, cls);
    for (const double v : nodes) {
      std::fprintf(out, " %.1f", v);
    }
    std::fprintf(out, "\n");
    ++nS;
  }

  // ===| edge part |=============================================================================================================
  int nE = 0;
  if (!edgeCSV.empty()) {
    const auto erows = readCSV(edgeCSV);
    constexpr int NSelfPhase = 20;
    std::map<std::tuple<int, int, int, int>, Profile> eprof; // [stack][drift][slope][class] in d bins
    std::map<std::tuple<int, int, int, int>, Profile> selfS; // the same table's S-curve at d >= sMin, 20 phase bins
    for (const auto& r : erows) {
      const double driftCentre = 0.5 * (r.driftLo + r.driftHi);
      const int iD = int(std::upper_bound(driftEdges.begin(), driftEdges.end(), driftCentre) - driftEdges.begin()) - 1;
      if (iD < 0 || iD >= nDrift) {
        continue;
      }
      const double d = 0.5 * (r.xLo + r.xHi);
      const std::tuple<int, int, int, int> key{r.stack, iD, r.slope, r.cls};
      if (d >= sMin) {
        const double p = (d - 0.5) - std::round(d - 0.5);
        auto& s = selfS.try_emplace(key, Profile(NSelfPhase)).first->second;
        s.fill(std::clamp(int(std::floor((p + 0.5) * NSelfPhase)), 0, NSelfPhase - 1), r.count, r.mean);
      }
      if (d >= DistLo && d < distHi) {
        auto& e = eprof.try_emplace(key, Profile(nDist)).first->second;
        e.fill(std::clamp(int(std::floor((d - DistLo) / DistW)), 0, nDist - 1), r.count, r.mean);
      }
    }
    for (auto& [key, e] : eprof) {
      const auto [stack, iD, iS, cls] = key;
      if (e.total() < minCount) {
        continue;
      }
      const auto self = selfS.find(key);
      std::vector<double> v(nDist, 0.);
      for (int i = 0; i < nDist; ++i) {
        v[i] = e.value(i);
      }
      if (cls == 2) {
        // single-pad clusters sit exactly on d = n + 0.5, the boundary of two d bins: both bins get their common value; only the
        // last pad (d < singleEdgeMax), no S-curve to subtract
        std::fill(v.begin(), v.end(), 0.);
        for (int n = 0; n + 0.5 < std::min(distHi, singleEdgeMax); ++n) {
          const int b = int(std::lround((n + 0.5 - DistLo) / DistW));
          if (b >= 1 && b < nDist) {
            const double c = e.cnt[b - 1] + e.cnt[b];
            v[b - 1] = v[b] = c > 0. ? (e.sum[b - 1] + e.sum[b]) / c : 0.;
          }
        }
      }
      for (int i = 0; i < nDist; ++i) {
        if (cls != 2 && self != selfS.end()) {
          // subtract the table's own S-curve at the phase of the bin
          const double d = DistLo + (i + 0.5) * DistW;
          const double p = (d - 0.5) - std::round(d - 0.5);
          v[i] -= self->second.value(std::clamp(int(std::floor((p + 0.5) * NSelfPhase)), 0, NSelfPhase - 1));
        }
        v[i] *= edgeScale;
      }
      std::fprintf(out, "E %d %d %d %d", stack, iD, iS, cls);
      for (const double x : v) {
        std::fprintf(out, " %.1f", x);
      }
      std::fprintf(out, "\n");
      ++nE;
    }
  }
  std::fclose(out);
  LOGP(info, "makeClusterBiasTable: {} S-curve and {} edge records ({} drift x 3 slope bins, {} phase nodes, {} d bins) written to {}", nS, nE, nDrift,
       NPhase, nDist, outFile);
}
