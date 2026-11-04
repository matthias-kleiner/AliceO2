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
/// slope / distance, the sign on the two sides of the pad plane, clamping, the cluster classes and the inactive / error paths.
///   g++ -std=c++17 -Wall -Wextra -I../macro testClusterBiasCorrection.cxx -o testClusterBiasCorrection && ./testClusterBiasCorrection

#include "ClusterBiasCorrection.h"

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
  std::printf("%d checks, %d failed\n", nChecks, nFailed);
  return nFailed == 0 ? 0 : 1;
}
