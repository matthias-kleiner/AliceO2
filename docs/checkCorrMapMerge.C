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

/// \file checkCorrMapMerge.C
/// \brief Checks the per-TF merge of TPC correction maps (TPCScalerSpec::buildMap, TPCFastSpaceChargeCorrectionHelper::addCorrections)
/// against the legacy runtime composition, and the inverse correction (initInverse) over the full pad-row area including the edges.
///
/// 1. Prints, per row, the spline grids of the two maps (knot positions in y and z, z where the fade starts). The merge is exact only
///    where the grids agree.
/// 2. merged = main + scale * extra with addCorrections (as buildMap does for the derivative, M-shape and sector-edge maps); compares
///    at random points with the legacy composition main(y, z) + scale * extra(y, z), for the direct and the inverse corrections,
///    split by region (inside / outside each map's grid).
/// 3. Inverse round trip on the main and on the merged map: real = measured + d(measured), inverse(real) vs d(measured), split by the
///    distance of the real point to the grid edge, and by what TPCFastSpaceChargeCorrection::testInverse checks (measured and real inside).
///
/// extra: "ccdb:<path>", a file with a TPCFastTransform "ccdb_object", or "test:<nKnotsY>,<nKnotsZ>" for a synthetic M-shape-like map
/// (createFromGlobalCorrection with uniform knots over the full row and z 0..250, as TPCScalerSpec builds the M-shape map).
///
/// root -b -q 'checkCorrMapMerge.C+("ccdb:TPC/Calib/CorrectionMapV2", "ccdb:TPC/Calib/CorrectionMapDerivativeV2", 1700000000000, 1.)'
/// root -b -q 'checkCorrMapMerge.C+("ccdb:TPC/Calib/CorrectionMapV2", "ccdb:TPC/Calib/CorrectionMapSecEdgeFluc", 1700000000000, 1.)'
/// root -b -q 'checkCorrMapMerge.C+("ccdb:TPC/Calib/CorrectionMapV2", "test:4,4", 1700000000000, 1.)'

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include <fmt/format.h>

#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/Logger.h"
#include "TPCCalibration/TPCFastSpaceChargeCorrectionHelper.h"
#include "TPCFastSpaceChargeCorrection.h"
#include "TPCFastTransform.h"
#endif

using o2::gpu::TPCFastSpaceChargeCorrection;
using o2::gpu::TPCFastTransform;

namespace check_merge
{

/// a TPCFastTransform from the CCDB ("ccdb:<path>") or from a file (object "ccdb_object"), as an own copy
std::unique_ptr<TPCFastTransform> loadMap(const std::string& source, long timestamp)
{
  std::unique_ptr<TPCFastTransform> copy(new TPCFastTransform);
  if (source.rfind("ccdb:", 0) == 0) {
    auto& ccdb = o2::ccdb::BasicCCDBManager::instance();
    auto* object = ccdb.getForTimeStamp<TPCFastTransform>(source.substr(5), timestamp);
    if (!object) {
      LOGP(fatal, "No {} at {}", source, timestamp);
    }
    object->rectifyAfterReadingFromFile();
    copy->cloneFromObject(*object, nullptr);
    return copy;
  }
  std::unique_ptr<TFile> file(TFile::Open(source.c_str()));
  auto* object = file ? file->Get<TPCFastTransform>("ccdb_object") : nullptr;
  if (!object) {
    LOGP(fatal, "No TPCFastTransform ccdb_object in {}", source);
  }
  object->rectifyAfterReadingFromFile();
  copy->cloneFromObject(*object, nullptr);
  return copy;
}

/// synthetic M-shape-like correction: smooth in z (two lobes over the drift), radial and z components, ~amplitude [cm]
std::unique_ptr<TPCFastSpaceChargeCorrection> makeTestMap(int nKnotsY, int nKnotsZ, double amplitude)
{
  auto* helper = o2::tpc::TPCFastSpaceChargeCorrectionHelper::instance();
  return helper->createFromGlobalCorrection(
    [amplitude](int sector, double gx, double gy, double gz, double& dgx, double& dgy, double& dgz) {
      const double az = std::min(1., std::abs(gz) / 250.);
      const double r = std::max(1., std::sqrt(gx * gx + gy * gy));
      const double f = amplitude * std::pow(std::sin(M_PI * az), 2) * std::cos(1.5 * M_PI * az) * (1. - (r - 85.) / 160.);
      dgx = f * gx / r;
      dgy = f * gy / r;
      dgz = 0.5 * f * (gz >= 0 ? 1. : -1.);
    },
    nKnotsY, nKnotsZ);
}

/// side-oriented z (as GridInfo::convLocalToGridUntruncated uses it)
float sideZ(int sector, float z) { return sector < 18 ? z : -z; }

/// inside the grid and before the z fade
bool insideGrid(const TPCFastSpaceChargeCorrection& c, int sector, int row, float y, float z, bool real)
{
  const auto& grid = real ? c.getRowInfo(row).gridReal : c.getRowInfo(row).gridMeasured;
  const bool inside = real ? c.isRealLocalInsideGrid(sector, row, y, z) : c.isLocalInsideGrid(sector, row, y, z);
  return inside && grid.getSpineScaleForZ(sideZ(sector, z)) >= 1.f;
}

/// per row: knots in y and z (local cm), and the z where the fade towards the readout starts
void printGrids(const char* name, const TPCFastSpaceChargeCorrection& c, const std::vector<int>& rows)
{
  for (int row : rows) {
    const auto& grid = c.getRowInfo(row).gridMeasured;
    const auto& spline = c.getSplineForRow(row);
    std::string knotsY, knotsZ;
    for (int i = 0; i < spline.getGridX1().getNumberOfKnots(); ++i) {
      knotsY += fmt::format("{:.2f} ", grid.getY0() + spline.getGridX1().getKnot(i).u / grid.getYscale());
    }
    for (int i = 0; i < spline.getGridX2().getNumberOfKnots(); ++i) {
      knotsZ += fmt::format("{:.1f} ", grid.getZ0() + spline.getGridX2().getKnot(i).u / grid.getZscale());
    }
    float fadeStart = -1.f;
    for (float z = 0.f; z <= 260.f; z += 0.1f) {
      if (grid.getSpineScaleForZ(z) < 1.f) {
        fadeStart = z;
        break;
      }
    }
    LOGP(info, "{} row {:3d}: y knots [cm] {}| z knots [cm] {}| fade from z = {:.1f} cm (-1: none)", name, row, knotsY, knotsZ, fadeStart);
  }
}

/// max and rms of |difference| per region and component
struct Stats {
  std::vector<std::string> labels;
  std::vector<std::array<double, 3>> max, sum2;
  std::vector<long> n;
  explicit Stats(std::vector<std::string> l) : labels(std::move(l)), max(labels.size()), sum2(labels.size()), n(labels.size(), 0) {}
  void fill(size_t region, const std::array<double, 3>& d)
  {
    ++n[region];
    for (int i = 0; i < 3; ++i) {
      max[region][i] = std::max(max[region][i], std::abs(d[i]));
      sum2[region][i] += d[i] * d[i];
    }
  }
  void print(const char* title) const
  {
    LOGP(info, "{} (um: max / rms for dx, dy, dz):", title);
    for (size_t r = 0; r < labels.size(); ++r) {
      if (!n[r]) {
        continue;
      }
      LOGP(info, "  {:58s} n {:8d}  dx {:8.1f} / {:7.1f}  dy {:8.1f} / {:7.1f}  dz {:8.1f} / {:7.1f}", labels[r], n[r],
           1e4 * max[r][0], 1e4 * std::sqrt(sum2[r][0] / n[r]), 1e4 * max[r][1], 1e4 * std::sqrt(sum2[r][1] / n[r]),
           1e4 * max[r][2], 1e4 * std::sqrt(sum2[r][2] / n[r]));
    }
  }
};

/// inverse round trip: for measured points over the full pad row, real = measured + d; inverse(real) - d
void roundTrip(const char* name, const TPCFastSpaceChargeCorrection& c, int nPoints, TRandom3& rnd)
{
  Stats stats({"measured and real inside the grid (what testInverse checks)", "  ... real point within 2 cm of the y grid edge",
               "  ... real point within 2 cm of the z grid edge", "measured or real outside the grid"});
  const auto& geo = c.getGeometry();
  for (int i = 0; i < nPoints; ++i) {
    const int sector = rnd.Integer(36);
    const int row = rnd.Integer(geo.getNumberOfRows());
    const auto [yMin, yMax] = geo.getRowInfo(row).getYrange();
    const auto [zMin, zMax] = geo.getZrange(sector);
    const float y = yMin + rnd.Rndm() * (yMax - yMin);
    const float z = zMin + rnd.Rndm() * (zMax - zMin);
    float dx, dy, dz;
    c.getCorrectionLocal(sector, row, y, z, dx, dy, dz);
    const float realY = y + dy;
    const float realZ = z + dz;
    float iy, iz;
    c.getCorrectionYZatRealYZ(sector, row, realY, realZ, iy, iz);
    const float ix = c.getCorrectionXatRealYZ(sector, row, realY, realZ);
    const std::array<double, 3> diff{ix - dx, iy - dy, iz - dz};
    if (c.isLocalInsideGrid(sector, row, y, z) && c.isRealLocalInsideGrid(sector, row, realY, realZ)) {
      stats.fill(0, diff);
      const auto& grid = c.getRowInfo(row).gridReal;
      const auto& spline = c.getSplineForRow(row);
      const float u = (realY - grid.getY0()) * grid.getYscale();
      const float v = (sideZ(sector, realZ) - grid.getZ0()) * grid.getZscale();
      if (std::min(u, spline.getGridX1().getUmax() - u) / grid.getYscale() < 2.f) {
        stats.fill(1, diff);
      }
      if (std::min(v, spline.getGridX2().getUmax() - v) / grid.getZscale() < 2.f) {
        stats.fill(2, diff);
      }
    } else {
      stats.fill(3, diff);
    }
  }
  stats.print(fmt::format("Inverse round trip of {}", name).c_str());
}

} // namespace check_merge

void checkCorrMapMerge(std::string mainSource = "ccdb:TPC/Calib/CorrectionMapV2",
                       std::string extraSource = "ccdb:TPC/Calib/CorrectionMapDerivativeV2", // or "test:<nKnotsY>,<nKnotsZ>"
                       long timestamp = -1, float scale = 1.f, int nPoints = 1000000, int nThreads = 4,
                       std::string ccdbUrl = "http://alice-ccdb.cern.ch")
{
  using namespace check_merge;
  if (timestamp < 0 && (mainSource.rfind("ccdb:", 0) == 0 || extraSource.rfind("ccdb:", 0) == 0)) {
    LOGP(fatal, "A timestamp [ms] is needed for the CCDB");
  }
  o2::ccdb::BasicCCDBManager::instance().setURL(ccdbUrl);
  auto* helper = o2::tpc::TPCFastSpaceChargeCorrectionHelper::instance();
  helper->setNthreads(nThreads);

  auto mainMap = loadMap(mainSource, timestamp);
  std::unique_ptr<TPCFastTransform> extraTransform;
  std::unique_ptr<TPCFastSpaceChargeCorrection> extraTest;
  if (extraSource.rfind("test:", 0) == 0) {
    int nKnotsY = 4;
    int nKnotsZ = 4;
    sscanf(extraSource.c_str() + 5, "%d,%d", &nKnotsY, &nKnotsZ);
    extraTest = makeTestMap(nKnotsY, nKnotsZ, 0.1);
  } else {
    extraTransform = loadMap(extraSource, timestamp);
  }
  const TPCFastSpaceChargeCorrection& main = mainMap->getCorrection();
  const TPCFastSpaceChargeCorrection& extra = extraTest ? *extraTest : extraTransform->getCorrection();

  // 1. the grids
  const std::vector<int> rows{0, 31, 62, 63, 96, 97, 126, 127, 151};
  printGrids("main ", main, rows);
  printGrids("extra", extra, rows);
  int nDifferentRows = 0;
  for (int row = 0; row < main.getGeometry().getNumberOfRows(); ++row) {
    const auto& gm = main.getRowInfo(row).gridMeasured;
    const auto& ge = extra.getRowInfo(row).gridMeasured;
    const auto& sm = main.getSplineForRow(row);
    const auto& se = extra.getSplineForRow(row);
    bool same = sm.getGridX1().getNumberOfKnots() == se.getGridX1().getNumberOfKnots() && sm.getGridX2().getNumberOfKnots() == se.getGridX2().getNumberOfKnots() &&
                std::abs(gm.getY0() - ge.getY0()) < 1e-4 && std::abs(gm.getYscale() - ge.getYscale()) < 1e-6 * std::abs(gm.getYscale()) &&
                std::abs(gm.getZ0() - ge.getZ0()) < 1e-4 && std::abs(gm.getZscale() - ge.getZscale()) < 1e-6 * std::abs(gm.getZscale());
    for (int i = 0; same && i < sm.getGridX1().getNumberOfKnots(); ++i) {
      same = sm.getGridX1().getKnot(i).u == se.getGridX1().getKnot(i).u;
    }
    for (int i = 0; same && i < sm.getGridX2().getNumberOfKnots(); ++i) {
      same = sm.getGridX2().getKnot(i).u == se.getGridX2().getKnot(i).u;
    }
    nDifferentRows += !same;
  }
  LOGP(info, "Rows with different spline grids (knots or y/z ranges): {} of {} -> the merge is {}", nDifferentRows, main.getGeometry().getNumberOfRows(),
       nDifferentRows ? "an approximation" : "exact up to float rounding (fade regions aside)");

  // 2. merged = main + scale * extra (addCorrections, as TPCScalerSpec::buildMap) vs the legacy runtime composition
  TPCFastTransform mergedTransform;
  mergedTransform.cloneFromObject(*mainMap, nullptr);
  helper->addCorrections(mergedTransform.getCorrection(), 1., {{&extra, double(scale)}});
  const TPCFastSpaceChargeCorrection& merged = mergedTransform.getCorrection();

  TRandom3 rnd(12345);
  Stats direct({"inside both grids", "inside the extra grid only (outside the main grid or in its fade)", "inside the main grid only",
                "outside both grids"});
  Stats inverse({"inside both real grids", "inside the extra real grid only", "inside the main real grid only", "outside both real grids"});
  const auto& geo = main.getGeometry();
  for (int i = 0; i < nPoints; ++i) {
    const int sector = rnd.Integer(36);
    const int row = rnd.Integer(geo.getNumberOfRows());
    const auto [yMin, yMax] = geo.getRowInfo(row).getYrange();
    const auto [zMin, zMax] = geo.getZrange(sector);
    const float y = yMin + rnd.Rndm() * (yMax - yMin);
    const float z = zMin + rnd.Rndm() * (zMax - zMin);
    float m[3], e[3], c[3];
    main.getCorrectionLocal(sector, row, y, z, m[0], m[1], m[2]);
    extra.getCorrectionLocal(sector, row, y, z, e[0], e[1], e[2]);
    merged.getCorrectionLocal(sector, row, y, z, c[0], c[1], c[2]);
    const bool inMain = insideGrid(main, sector, row, y, z, false);
    const bool inExtra = insideGrid(extra, sector, row, y, z, false);
    const size_t region = inMain && inExtra ? 0 : (inExtra ? 1 : (inMain ? 2 : 3));
    direct.fill(region, {c[0] - (m[0] + scale * e[0]), c[1] - (m[1] + scale * e[1]), c[2] - (m[2] + scale * e[2])});
    // inverse at the same (y, z) taken as a real position
    float my, mz, ey, ez, cy, cz;
    main.getCorrectionYZatRealYZ(sector, row, y, z, my, mz);
    extra.getCorrectionYZatRealYZ(sector, row, y, z, ey, ez);
    merged.getCorrectionYZatRealYZ(sector, row, y, z, cy, cz);
    const float mx = main.getCorrectionXatRealYZ(sector, row, y, z);
    const float ex = extra.getCorrectionXatRealYZ(sector, row, y, z);
    const float cx = merged.getCorrectionXatRealYZ(sector, row, y, z);
    const bool inMainReal = insideGrid(main, sector, row, y, z, true);
    const bool inExtraReal = insideGrid(extra, sector, row, y, z, true);
    const size_t regionReal = inMainReal && inExtraReal ? 0 : (inExtraReal ? 1 : (inMainReal ? 2 : 3));
    inverse.fill(regionReal, {cx - (mx + scale * ex), cy - (my + scale * ey), cz - (mz + scale * ez)});
  }
  direct.print("Direct correction, merged - (main + scale * extra)");
  inverse.print("Inverse correction, merged - (main + scale * extra)");

  // 3. inverse round trips over the full pad-row area
  roundTrip("the main map", main, nPoints, rnd);
  roundTrip("the merged map", merged, nPoints, rnd);
}
