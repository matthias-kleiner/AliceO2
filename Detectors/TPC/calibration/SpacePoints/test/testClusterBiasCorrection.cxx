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

/// \file testClusterBiasCorrection.cxx
/// \brief Standalone test of macro/ClusterBiasCorrection.h (no ROOT, no Boost): table reading, interpolation in phase / drift /
/// slope / distance, the sign on the two sides of the pad plane, clamping, the cluster classes and the inactive / error paths; and of
/// macro/ClusterBiasSmoothFit.h (the smooth-table numerics of makeClusterBiasTable.C): B-spline basis, solver, penalised fits.
///   g++ -std=c++17 -Wall -Wextra -I../macro testClusterBiasCorrection.cxx -o testClusterBiasCorrection && ./testClusterBiasCorrection

#include "ClusterBiasCorrection.h"
#include "ClusterBiasSmoothFit.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

namespace
{
int nFailed = 0;
int nChecks = 0;

void check(bool ok, const std::string& what)
{
  ++nChecks;
  if (!ok) {
    ++nFailed;
    std::printf("FAILED: %s\n", what.c_str());
  }
}

void checkClose(double a, double b, double tol, const std::string& what)
{
  const bool ok = std::abs(a - b) <= tol;
  if (!ok) {
    std::printf("  got %.9g expected %.9g (tol %.2g)\n", a, b, tol);
  }
  check(ok, what);
}

constexpr double TwoPi = 6.283185307179586;
const double AmpS[2][2] = {{-500., -300.}, {-200., -100.}}; // [drift bin][slope bin], um: inward bias -A sin(2 pi phase)
const double PadW = 0.6;
const int NPads = 80;
const int Row = 70; // OROC1

double scurveTrue(int iD, int iS, double phase) { return AmpS[iD][iS] * std::sin(TwoPi * phase); }
double edgeTrue(int iD, double dCentre) { return (iD + 1) * (1000. - 1000. * (dCentre - 0.4) / 0.5); }

std::string writeTable(const std::string& name)
{
  std::ofstream o(name);
  o << "# test table\nversion 1\ndrift 2 0 50 100\nslope 2 0 0.1 0.3\nphase 64\ndist 10 0.4 0.9\n";
  for (int iD = 0; iD < 2; ++iD) {
    for (int iS = 0; iS < 2; ++iS) {
      o << "S 1 " << iD << " " << iS << " 0";
      for (int j = 0; j < 64; ++j) {
        o << " " << scurveTrue(iD, iS, -0.5 + j / 64.);
      }
      o << "\nE 1 " << iD << " " << iS << " 0";
      for (int i = 0; i < 10; ++i) {
        o << " " << edgeTrue(iD, 0.4 + 0.05 * (i + 0.5));
      }
      o << "\n";
    }
  }
  o << "S 1 0 0 2";
  for (int j = 0; j < 64; ++j) {
    o << " 77";
  }
  o << "\n";
  return name;
}

/// snapped y of a cluster at distance d [pads] from the edge on the given side, as snapClusterY leaves it (float)
float yAt(double d, int side) { return float(side * (0.5 * NPads - d) * PadW); }
} // namespace

int main()
{
  using o2::tpc::cluster_bias::ClusterBiasCorrection;
  const std::string file = writeTable("testClusterBiasCorrection.table.txt");
  ClusterBiasCorrection corr;
  check(!corr.isActive(), "inactive before load");
  check(corr.dy(Row, 10.f, 25.f, 0.05f, 0) == 0.f, "inactive: no correction");
  std::string err;
  check(corr.load(file, &err), "load: " + err);
  check(corr.isActive(), "active after load");
  check(corr.dy(Row, 10.f, 25.f, 0.05f, 0) == 0.f, "row without geometry: no correction");
  corr.setRowGeometry(Row, PadW, NPads);
  const double um = 1.e-4;

  // 1. S-curve on a grid node (exact), both sides: dy = sign(y) * inward bias
  {
    const double d = 10.5 + 5. / 64.; // pad 10 + 5/64 from the edge: phase +5/64
    const double b = scurveTrue(0, 0, 5. / 64.) * um;
    checkClose(corr.dy(Row, yAt(d, +1), 25.f, 0.05f, 0), +b, 1e-7, "S-curve on a node, +y side");
    checkClose(corr.dy(Row, yAt(d, -1), 25.f, 0.05f, 0), -b, 1e-7, "S-curve on a node, -y side (sign flips)");
    // physics sign: phase > 0 (cluster inward of the pad centre) carries an outward pull (b < 0); the corrected cluster moves
    // further inward, i.e. |y| decreases on both sides
    const float yp = yAt(d, +1);
    const float ym = yAt(d, -1);
    check(std::abs(yp + corr.dy(Row, yp, 25.f, 0.05f, 0)) < std::abs(yp), "corrected cluster moves inward, +y side");
    check(std::abs(ym + corr.dy(Row, ym, 25.f, 0.05f, 0)) < std::abs(ym), "corrected cluster moves inward, -y side");
  }
  // 2. periodic phase interpolation between nodes, across phase = +-0.5
  {
    const double d = 10.5 + 0.5 / 64.; // halfway between the nodes 0 and 1/64
    const double b = 0.5 * (scurveTrue(0, 0, 0.) + scurveTrue(0, 0, 1. / 64.)) * um;
    checkClose(corr.dy(Row, yAt(d, +1), 25.f, 0.05f, 0), b, 1e-7, "phase interpolation between nodes");
    const double d2 = 10.5 + 31.5 / 64.; // between +31/64 and -32/64 (= the node -0.5)
    const double b2 = 0.5 * (scurveTrue(0, 0, 31. / 64.) + scurveTrue(0, 0, -0.5)) * um;
    checkClose(corr.dy(Row, yAt(d2, +1), 25.f, 0.05f, 0), b2, 1e-6, "periodic interpolation across phase 0.5"); // float y
  }
  // 3. drift and slope interpolation, clamping, |slope|
  {
    const double d = 20.5 + 16. / 64.; // phase 0.25: sin = 1
    checkClose(corr.dy(Row, yAt(d, +1), 50.f, 0.05f, 0), 0.5 * (AmpS[0][0] + AmpS[1][0]) * um, 1e-7, "drift interpolation");
    checkClose(corr.dy(Row, yAt(d, +1), 0.f, 0.05f, 0), AmpS[0][0] * um, 1e-7, "drift clamped below the first centre");
    checkClose(corr.dy(Row, yAt(d, +1), 400.f, 0.05f, 0), AmpS[1][0] * um, 1e-7, "drift clamped above the last centre");
    checkClose(corr.dy(Row, yAt(d, +1), 25.f, 0.125f, 0), 0.5 * (AmpS[0][0] + AmpS[0][1]) * um, 1e-7, "slope interpolation");
    checkClose(corr.dy(Row, yAt(d, +1), 25.f, -0.125f, 0), 0.5 * (AmpS[0][0] + AmpS[0][1]) * um, 1e-7, "slope sign ignored");
    checkClose(corr.dy(Row, yAt(d, +1), 25.f, 0.9f, 0), AmpS[0][1] * um, 1e-7, "slope clamped");
    const double bilin = 0.25 * (AmpS[0][0] + AmpS[0][1] + AmpS[1][0] + AmpS[1][1]) * um;
    checkClose(corr.dy(Row, yAt(d, +1), 50.f, 0.125f, 0), bilin, 1e-7, "bilinear drift x slope");
  }
  // 4. edge part: linear in d between the centres, clamped below, 0 beyond dHi; added to the S-curve
  {
    const double d = 0.5 + 6. / 64.; // between the centres 0.575 and 0.625
    const double f = (d - 0.575) / 0.05;
    const double e = (1. - f) * edgeTrue(0, 0.575) + f * edgeTrue(0, 0.625);
    checkClose(corr.inwardBias(1, 0, 25., 0.05, d), (e + scurveTrue(0, 0, 6. / 64.)) * um, 1e-7, "edge + S-curve");
    checkClose(corr.dy(Row, yAt(d, -1), 25.f, 0.05f, 0), -(e + scurveTrue(0, 0, 6. / 64.)) * um, 2e-6, "edge on the -y side");
    // (d on phase nodes, so that the S-curve part is exact and only the edge part is tested)
    checkClose(corr.inwardBias(1, 0, 25., 0.05, 0.5 - 6. / 64.) - scurveTrue(0, 0, -6. / 64.) * um, edgeTrue(0, 0.425) * um, 1e-7,
               "edge clamped below the first centre");
    checkClose(corr.inwardBias(1, 0, 25., 0.05, 0.5 + 26. / 64.) - scurveTrue(0, 0, 26. / 64.) * um, 0., 1e-7, "edge zero beyond dHi");
    const double dLast = 0.5 + 24. / 64.; // 0.875 = the last centre: full value; halfway to dHi: half of it
    checkClose(corr.inwardBias(1, 0, 25., 0.05, dLast) - scurveTrue(0, 0, 24. / 64.) * um, edgeTrue(0, 0.875) * um, 1e-7, "edge at the last centre");
    // halfway between the last centre (0.875) and dHi (0.9): half of the last value (S-curve part: linear between nodes 56 / 57)
    const double xs = 0.8875 * 64. - 56.;
    const double sLin = (1. - xs) * scurveTrue(0, 0, -0.5 + 56. / 64.) + xs * scurveTrue(0, 0, -0.5 + 57. / 64.);
    checkClose(corr.inwardBias(1, 0, 25., 0.05, 0.8875), (sLin + 0.5 * edgeTrue(0, 0.875)) * um, 1e-7, "edge falls to 0 towards dHi");
    checkClose(corr.inwardBias(1, 0, 25., 0.05, 3.25) - scurveTrue(0, 0, -0.25) * um, 0., 1e-9, "edge zero beyond dHi");
    checkClose(corr.inwardBias(1, 0, 75., 0.05, 0.875) - scurveTrue(1, 0, 0.375) * um, edgeTrue(1, 0.875) * um, 1e-7, "edge at the last centre");
  }
  // 5. classes from the flags: single pad (0x8 on the pad centre) has its own record where present, else class 0; 0x8 off the pad
  //    centre is a single-time cluster (class 0); split / shared fall back to class 0 without records
  {
    const double dCentre = 20.5;          // phase 0: single-pad clusters sit here
    const double dOff = 20.5 + 16. / 64.; // phase 0.25
    checkClose(corr.dy(Row, yAt(dCentre, +1), 25.f, 0.05f, 0x8), 77. * um, 1e-7, "single pad on the pad centre: class-2 record");
    checkClose(corr.dy(Row, yAt(dCentre, -1), 25.f, 0.05f, 0x9), -77. * um, 1e-7, "single wins over split pad, -y side");
    checkClose(corr.dy(Row, yAt(dCentre, +1), 25.f, 0.05f, 0), scurveTrue(0, 0, 0.) * um, 1e-7, "golden on the pad centre: class 0");
    checkClose(corr.dy(Row, yAt(dOff, +1), 25.f, 0.05f, 0x8), AmpS[0][0] * um, 1e-7, "0x8 off the pad centre (single time): class 0");
    checkClose(corr.dy(Row, yAt(dCentre, +1), 75.f, 0.05f, 0x8), scurveTrue(1, 0, 0.) * um, 1e-7, "single without record: class 0");
    checkClose(corr.dy(Row, yAt(dOff, +1), 25.f, 0.05f, 0x10), AmpS[0][0] * um, 1e-7, "shared: class 1 missing -> class 0");
    checkClose(corr.dy(Row, yAt(dOff, +1), 25.f, 0.05f, -1), AmpS[0][0] * um, 1e-7, "no flags stored (-1): golden");
  }
  // 6. a stack without records: no correction
  checkClose(corr.dy(10, 5.f, 25.f, 0.05f, 0), 0., 0., "IROC row without geometry");
  corr.setRowGeometry(10, 0.42, 66);
  checkClose(corr.dy(10, 5.f, 25.f, 0.05f, 0), 0., 0., "IROC row without records");
  // 7. per pad region: regionOfRow, and records of region 4 (rows 63-80) vs region 5 (rows 81-96)
  {
    using o2::tpc::cluster_bias::regionOfRow;
    check(regionOfRow(0) == 0 && regionOfRow(16) == 0 && regionOfRow(17) == 1 && regionOfRow(62) == 3 && regionOfRow(63) == 4 &&
            regionOfRow(81) == 5 && regionOfRow(127) == 8 && regionOfRow(140) == 9 && regionOfRow(151) == 9,
          "regionOfRow at the region boundaries");
    std::ofstream o("testClusterBiasCorrection.region.txt");
    o << "version 1\nunit region\ndrift 1 0 250\nslope 1 0 1\nphase 64\ndist 0 0 0\n";
    for (int region : {4, 5}) {
      o << "S " << region << " 0 0 0";
      for (int j = 0; j < 64; ++j) {
        o << " " << (region == 4 ? 100. : 300.) * std::sin(TwoPi * (-0.5 + j / 64.));
      }
      o << "\n";
    }
    o.close();
    ClusterBiasCorrection reg;
    check(reg.load("testClusterBiasCorrection.region.txt", &err) && reg.byRegion && reg.nUnits() == 10, "region table: " + err);
    reg.setRowGeometry(70, PadW, NPads);    // region from the row: 4
    reg.setRowGeometry(85, PadW, NPads, 5); // region given
    const double d = 20.5 + 16. / 64.;      // phase 0.25
    checkClose(reg.dy(70, yAt(d, +1), 100.f, 0.f, 0), 100. * um, 1e-7, "region 4 record for row 70");
    checkClose(reg.dy(85, yAt(d, -1), 100.f, 0.f, 0), -300. * um, 1e-7, "region 5 record for row 85, -y side");
    std::ofstream("testClusterBiasCorrection.badunit.txt") << "version 1\nunit sector\ndrift 1 0 50\nslope 1 0 1\nphase 4\n";
    check(!reg.load("testClusterBiasCorrection.badunit.txt", &err), "unknown unit rejected");
    check(corr.load(file, &err) && !corr.byRegion, "stack table again: unit reset to stack");
  }
  // 8. format errors
  {
    std::ofstream("testClusterBiasCorrection.bad.txt") << "version 1\ndrift 1 0 50\nslope 1 0 1\nphase 4\ndist 0 0 0\nS 0 0 0 0 1 2 3\n";
    ClusterBiasCorrection bad;
    check(!bad.load("testClusterBiasCorrection.bad.txt", &err) && !bad.isActive(), "wrong number of values rejected");
    check(!bad.load("does_not_exist.txt", &err), "missing file rejected");
    std::ofstream("testClusterBiasCorrection.bad2.txt") << "drift 1 0 50\nslope 1 0 1\nphase 4\n";
    check(!bad.load("testClusterBiasCorrection.bad2.txt", &err), "missing version rejected");
  }
  // 9. smooth-fit numerics (ClusterBiasSmoothFit.h)
  {
    namespace sm = o2::tpc::cluster_bias::smooth;
    const std::vector<double> knots = {0., 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5};
    bool unity = true;
    bool nonNeg = true;
    for (int i = 0; i <= 100; ++i) {
      const auto b = sm::bsplineBasis(0.5 * i / 100., knots);
      double sum = 0.;
      for (const double v : b) {
        sum += v;
        nonNeg = nonNeg && v >= -1e-15;
      }
      unity = unity && std::abs(sum - 1.) < 1e-12;
    }
    check(sm::bsplineBasis(0.1, knots).size() == knots.size() + 2, "cubic basis size = knots + 2");
    check(unity && nonNeg, "B-spline basis: partition of unity, non-negative (incl. the clipped right end)");
    checkClose(sm::bsplineBasis(0., knots)[0], 1., 1e-12, "clamped: first basis = 1 at the left knot");
    const auto x = sm::solve({{4., 1., 2.}, {1., 3., 0.}, {2., 0., 5.}}, {1., 2., 3.});
    checkClose(4. * x[0] + x[1] + 2. * x[2], 1., 1e-12, "solve: row 1");
    checkClose(x[0] + 3. * x[1], 2., 1e-12, "solve: row 2");
    checkClose(2. * x[0] + 5. * x[2], 3., 1e-12, "solve: row 3");
    check(sm::solve({{1., 2.}, {2., 4.}}, {1., 2.}).empty(), "solve: singular matrix -> empty");
    const auto h = sm::driftHats(60., {25., 75., 125.});
    checkClose(h[0], 0.3, 1e-12, "drift hats: weight of the lower centre");
    checkClose(h[1], 0.7, 1e-12, "drift hats: weight of the upper centre");
    checkClose(sm::driftHats(300., {25., 75., 125.})[2], 1., 1e-12, "drift hats: constant beyond the last centre");
    // S: a cubic in |u| vanishing at 0.5 times a linear function of x is in the spline space: recovered without penalty
    const std::vector<double> xKnots = {1., 2., 4., 8.};
    auto truth = [](double u, double xx) { return (u - 0.5) * (u + 0.3) * (u - 2.) * (100. + 10. * xx); };
    std::vector<sm::SCurvePoint> pts;
    for (int k = 1; k < 32; ++k) {
      for (const double xx : {1., 1.5, 2.5, 3., 5., 6., 7.9}) {
        pts.push_back({k / 64., xx, truth(k / 64., xx), 1.});
      }
    }
    const auto f = sm::fitSCurve(pts, knots, xKnots, 0.);
    check(f.valid() && f.nu == 10 && f.nx == 6, "S fit: valid, 10 x 6 coefficients");
    double maxDev = 0.;
    for (const double u : {0.013, 0.11, 0.27, 0.43}) {
      for (const double xx : {1.2, 3.7, 7.}) {
        maxDev = std::max(maxDev, std::abs(sm::evalSCurve(f, knots, xKnots, u, xx) - truth(u, xx)));
      }
    }
    checkClose(maxDev, 0., 0.02, "S fit: an exactly representable surface is recovered between the points (up to the 1e-9 ridge)");
    checkClose(sm::evalSCurve(f, knots, xKnots, 0.5, 3.), 0., 1e-9, "S fit: S(0.5) = 0 by construction");
    const auto fp = sm::fitSCurve(pts, knots, xKnots, 1e3);
    check(fp.valid() && fp.chi2 > f.chi2, "S fit: a strong penalty costs chi2");
    // E: linear in d (zero at the last knot) x linear in drift, recovered; single class = interpolation weights only
    const std::vector<double> dKnots = {0.5, 0.7, 1.0, 1.4, 2.4};
    const std::vector<double> centres = {25., 75., 125.};
    auto etruth = [](double d, double dr) { return (2.4 - d) * (500. + 4. * dr); };
    std::vector<sm::EdgePoint> ep;
    for (int i = 0; i < 38; ++i) {
      for (const double dr : centres) {
        ep.push_back({0.525 + 0.05 * i, dr, etruth(0.525 + 0.05 * i, dr), 10.});
      }
    }
    const auto e = sm::fitEdge(ep, centres, dKnots, 0., false);
    check(e.valid() && e.nd == 6, "E fit: valid, 6 d coefficients");
    checkClose(sm::evalEdge(e, dKnots, 1.23, 100.), etruth(1.23, 100.), 2e-3, "E fit: recovered between the drift centres");
    checkClose(sm::evalEdge(e, dKnots, 0.3, 75.), etruth(0.5, 75.), 2e-3, "E fit: constant below the first knot");
    checkClose(sm::evalEdge(e, dKnots, 2.5, 75.), 0., 1e-12, "E fit: 0 from the last knot on");
    std::vector<sm::EdgePoint> sp = {{0.525, 25., 3000., 10.}, {0.525, 75., 3400., 10.}, {0.525, 125., 3600., 10.}};
    const auto se = sm::fitEdge(sp, centres, dKnots, 0., true);
    checkClose(sm::evalEdge(se, dKnots, 0.5, 100., 1.), 3500., 2e-3, "E single: linear in drift between the centres");
    checkClose(sm::evalEdge(se, dKnots, 1.2, 100., 1.), 0., 1e-12, "E single: 0 from singleMax on");
  }
  std::printf("%d checks, %d failed\n", nChecks, nFailed);
  return nFailed == 0 ? 0 : 1;
}
