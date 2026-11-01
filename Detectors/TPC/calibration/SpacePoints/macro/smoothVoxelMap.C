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

/// \file smoothVoxelMap.C
/// \brief Smoothing and extrapolation of a TPC voxel map as a penalised fit: one node (dX, dY, dZ) per voxel of the map-creation grid,
/// the voxels with data entering with their chi2 as a quadratic form at their position (bilinear between the nodes of their pad row),
/// and second-difference penalties along the pad rows, y/x and z/x with the prior rms tauRow, tauY2X, tauZ2X (none across the stack
/// boundaries) and optionally a ridge towards 0 (tau0). Nodes without data follow the penalties: a smooth continuation, linear in the
/// direction without data and coherent between the neighbouring columns, e.g. through the innermost IROC rows, which the large radial
/// distortion empties (trustFraction: those rows are no data). Replaces the Gaussian smoothing + per-column pol2 extrapolation of
/// SmoothingExtrapolate.C (a leave-rows-out test on PbPb 2025 maps: tauRow 100 um, tauY2X = tauZ2X = 1000 um predicted the hidden
/// rows best and is ~10-100 times smoother in the extrapolated rows).
///
/// Input (source):
/// - "voxResTree": D with the errors E (dX, dY correlated by EXYCorr) at the voxel's mean position (stat), for the fitted, unmasked voxels
///   with >= minEntries entries (staticMapCreator.C, staticMapCreatorCPM.C, fineMapCreatorGBL.C). errorMode 0: the input's E, EXYCorr;
///   1: E = s / sqrt(entries), EXYCorr = 0, with s per coordinate from the scatter of the rows' second differences in the IROC bulk (an
///   upper bound: it includes the curvature differences between columns); -1: 1 if the input has the placeholders of
///   staticMapCreatorCPM.C's CPM maps (E = 0.1, EXYCorr = 1), else 0;
/// - "fineSums": per voxel and charge the sums of the residual slope method (staticMapCreatorCPM.C writeSlopeSums, fineMapCreatorGBL.C
///   writeSums), solved here: D, E from the voxel solve; every voxel with entries enters with its chi2 scaled by max(1, chi2 / ndf) of its
///   solve (the median of the solved voxels below minEntries);
/// - "auto": fineSums if the file has them, else voxResTree.
/// trustFraction > 0: per column the innermost IROC rows with fewer entries than trustFraction x the column's median entries in rows
/// 10-30 (up to the first row at or above it) are no data. robust > 0: Tukey reweightings of the voxels by their excess chi2.
/// curvaturePrior: the taus are the prior rms of the curvature (um / cm^2) averaged over 1 cm^3 instead of the second difference per bin
/// (um): second derivatives from the bin centres' positions, each term weighted with its voxel's volume, so the prior does not depend
/// on the binning (coarse map, fine pass, refined bins) nor on the radius (the y/x and z/x bins grow with x).
///
/// Output: voxResTree with every voxel of the grid, DS = the nodes, SmoothDone, entries at least 1 (TPCFastSpaceChargeCorrectionHelper
/// takes voxels with entries as the map's points); D, E: the input's (voxResTree; with errorMode 1 the E used) or the voxel solve
/// (fineSums). User info: the input voxResTree's (else the binning and nTracks of fineSums) + "smoothing"; the input's other trees
/// (OrbitLumiInfo, MetaData) copied. joinIDC: meanIDC / medianIDC = mean / median of the TPC scalers (A + C) / 2 at the TFs of the
/// input's OrbitLumiInfo (timeMSsel; CCDB TPC/Calib/Scaler and ScalerWeights, as SmoothingExtrapolate.C), -1: only if the input has no
/// meanIDC > 0 (staticMapCreatorCPM.C leaves them 0), 1: always, 0: never.
///
/// root -b -q 'smoothVoxelMap.C+("voxRes.root", "voxResSmooth.root")'

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

#include <Eigen/Dense>
#include <fmt/format.h>

#include "TFile.h"
#include "TKey.h"
#include "TNamed.h"
#include "TTree.h"

#include "Algorithm/RangeTokenizer.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonUtils/ConfigurableParam.h"
#include "Framework/Logger.h"
#include "SpacePoints/TrackResiduals.h"
#include "TPCCalibration/TPCScaler.h"
#if __has_include("TPCBaseRecSim/CDBTypes.h")
#include "TPCBaseRecSim/CDBTypes.h"
#else
#include "TPCBase/CDBTypes.h"
#endif
#else
#error This macro must run in compiled mode
#endif

namespace voxel_smoother
{
using TrackResiduals = o2::tpc::TrackResiduals;

constexpr int NSectors = 36;
constexpr int NRows = 152;
constexpr int NIrocRows = 63;
constexpr int PlateauFirstRow = 10; ///< the rows whose median entries are a column's plateau (trustFraction)
constexpr int PlateauLastRow = 30;
constexpr double Micron = 1e-4;          ///< in cm
constexpr double Tukey = 4.685;          ///< Tukey's bisquare constant of the robust reweighting
constexpr double CGTolerance = 1e-8;     ///< relative residual norm of the conjugate gradients
constexpr int CGMaxIterations = 20000;   ///< per sector and cycle
constexpr float PlaceholderError = 0.1f; ///< staticMapCreatorCPM.C's E of every voxel

std::unique_ptr<TrackResiduals> makeBinning(const std::string& y2xBinning, const std::string& z2xBinning, float maxZ2X)
{
  o2::conf::ConfigurableParam::setValue("scdcalib.maxZ2X", std::to_string(maxZ2X));
  auto binning = std::make_unique<TrackResiduals>();
  binning->setY2XBinning(o2::RangeTokenizer::tokenize<float>(y2xBinning));
  binning->setZ2XBinning(o2::RangeTokenizer::tokenize<float>(z2xBinning));
  binning->init();
  return binning;
}

/// the mean TPC scalers (IDC) of the A and C sides at a TF time from the CCDB objects TPC/Calib/Scaler and ScalerWeights (as
/// SmoothingExtrapolate.C); false if missing or not positive
bool getScalerValues(o2::ccdb::BasicCCDBManager& manager, long timeMS, float& scalerA, float& scalerC)
{
  auto* scalerTree = manager.getForTimeStamp<TTree>(o2::tpc::CDBTypeMap.at(o2::tpc::CDBType::CalScaler), timeMS);
  auto* scalerWeights = manager.getForTimeStamp<o2::tpc::TPCScalerWeights>(o2::tpc::CDBTypeMap.at(o2::tpc::CDBType::CalScalerWeights), timeMS);
  if (!scalerTree || !scalerWeights) {
    LOGP(error, "IDC join: no TPC/Calib/Scaler or ScalerWeights for the time stamp {}", timeMS);
    return false;
  }
  o2::tpc::TPCScaler scaler;
  scaler.setFromTree(*scalerTree);
  scaler.setScalerWeights(*scalerWeights);
  scaler.useWeights(true);
  scaler.setIonDriftTimeMS(500);
  if (scaler.getRun() == 0) {
    return false; // the default object
  }
  scalerA = scaler.getMeanScaler(timeMS, o2::tpc::Side::A);
  scalerC = scaler.getMeanScaler(timeMS, o2::tpc::Side::C);
  return scalerA > 0.f && scalerC > 0.f;
}

/// the IDC join: mean and median of (A + C) / 2 over the TFs of the input's OrbitLumiInfo (timeMSsel, sorted: the CCDB cache hits while
/// consecutive TFs share a validity window); false without TFs or valid scalers
bool joinIDC(TFile& input, float& meanIDC, float& medianIDC)
{
  auto* orbitLumiTree = input.Get<TTree>("OrbitLumiInfo");
  std::vector<long>* times = nullptr;
  if (!orbitLumiTree || orbitLumiTree->GetEntries() == 0 || orbitLumiTree->SetBranchAddress("timeMSsel", &times) < 0) {
    LOGP(warning, "IDC join: no OrbitLumiInfo with timeMSsel in the input");
    return false;
  }
  orbitLumiTree->GetEntry(0);
  if (!times || times->empty()) {
    LOGP(warning, "IDC join: no TF times in OrbitLumiInfo");
    return false;
  }
  std::vector<long> sortedTimes(*times);
  std::sort(sortedTimes.begin(), sortedTimes.end());
  auto& manager = o2::ccdb::BasicCCDBManager::instance();
  manager.setCaching(true);
  manager.setFatalWhenNull(false);
  manager.setURL("http://alice-ccdb.cern.ch");
  std::vector<float> idcs;
  idcs.reserve(sortedTimes.size());
  for (const long timeMS : sortedTimes) {
    float scalerA = 0.f;
    float scalerC = 0.f;
    if (getScalerValues(manager, timeMS, scalerA, scalerC)) {
      idcs.push_back(0.5f * (scalerA + scalerC));
    }
  }
  if (idcs.empty()) {
    LOGP(warning, "IDC join: no valid scaler for any of the {} TFs", sortedTimes.size());
    return false;
  }
  meanIDC = std::accumulate(idcs.begin(), idcs.end(), 0.) / idcs.size();
  std::nth_element(idcs.begin(), idcs.begin() + idcs.size() / 2, idcs.end());
  medianIDC = idcs[idcs.size() / 2];
  LOGP(info, "IDC join: {} of {} TFs with a valid scaler, meanIDC {}, medianIDC {}", idcs.size(), sortedTimes.size(), meanIDC, medianIDC);
  return true;
}

/// the stack (IROC, OROC1 / 2 / 3) of a pad row
int stackOfRow(int row)
{
  return row < 63 ? 0 : (row < 97 ? 1 : (row < 127 ? 2 : 3));
}

/// the bilinear interpolation between the bin centres of the pad row in y/x and the side-oriented z/x (z/x on the A side, -z/x on the C
/// side; constant beyond the outermost centres): up to 4 voxels (sector-local index) and their weights; returns their number
int bilinearNodes(const TrackResiduals& binning, int row, double y2x, double z2xSide, int voxels[4], double weights[4])
{
  auto bracket = [](int nBins, double value, auto centre, int& low, double& fraction) {
    low = 0;
    fraction = 0.;
    if (nBins < 2 || value <= centre(0)) {
      return;
    }
    if (value >= centre(nBins - 1)) {
      low = nBins - 1;
      return;
    }
    while (low < nBins - 2 && value >= centre(low + 1)) {
      ++low;
    }
    fraction = (value - centre(low)) / (centre(low + 1) - centre(low));
  };
  int binY2X = 0;
  int binZ2X = 0;
  double fractionY2X = 0.;
  double fractionZ2X = 0.;
  bracket(binning.getNY2XBins(), y2x, [&](int bin) { return double(binning.getY2X(row, bin)); }, binY2X, fractionY2X);
  bracket(binning.getNZ2XBins(), z2xSide, [&](int bin) { return double(binning.getZ2X(bin)); }, binZ2X, fractionZ2X);
  int nVoxels = 0;
  for (int stepY2X = 0; stepY2X < 2; ++stepY2X) {
    for (int stepZ2X = 0; stepZ2X < 2; ++stepZ2X) {
      const double weight = (stepY2X ? fractionY2X : 1. - fractionY2X) * (stepZ2X ? fractionZ2X : 1. - fractionZ2X);
      if (weight > 0.) {
        voxels[nVoxels] = int(binning.getGlbVoxBin(row, binY2X + stepY2X, binZ2X + stepZ2X));
        weights[nVoxels++] = weight;
      }
    }
  }
  return nVoxels;
}

/// per voxel and charge of the residual slope method: entries, position sums, per channel (y, z) the sums w, w s, w s^2, w r, w s r,
/// w r^2 (staticMapCreatorCPM.C, fineMapCreatorGBL.C)
struct ChargeSums {
  enum Sum { W,
             WS,
             WSS,
             WR,
             WSR,
             WRR,
             NSums };
  double entries = 0.;
  double x = 0.;
  double y2x = 0.;
  double z2x = 0.;
  std::array<double, NSums> y{};
  std::array<double, NSums> z{};
};

struct VoxelSolution {
  std::array<double, 3> value{}; ///< dX, dY, dZ
  std::array<double, 3> error{}; ///< scaled by sqrt(chi2Scale)
  double correlationXY = 0.;
  double rmsY = 0.;
  double rmsZ = 0.;
  double chi2Scale = 1.; ///< max(1, chi2 / ndf)
};

/// the voxel solve of the residual slope method: dX, dY from the y channel (2 x 2), then dZ = <r_z> + <s_z> dX; false below minEntries
/// or singular
bool solveVoxel(const ChargeSums& sums, double minEntries, VoxelSolution& solution)
{
  using S = ChargeSums;
  if (sums.entries < minEntries || sums.y[S::W] <= 0. || sums.z[S::W] <= 0.) {
    return false;
  }
  const double determinant = sums.y[S::WSS] * sums.y[S::W] - sums.y[S::WS] * sums.y[S::WS];
  if (!(determinant > 1e-12 * sums.y[S::WSS] * sums.y[S::W])) {
    return false;
  }
  const double covarianceXX = sums.y[S::W] / determinant;
  const double covarianceXY = sums.y[S::WS] / determinant;
  const double covarianceYY = sums.y[S::WSS] / determinant;
  const double dX = -covarianceXX * sums.y[S::WSR] + covarianceXY * sums.y[S::WR];
  const double dY = -covarianceXY * sums.y[S::WSR] + covarianceYY * sums.y[S::WR];
  const double meanSlopeZ = sums.z[S::WS] / sums.z[S::W];
  const double dZ = sums.z[S::WR] / sums.z[S::W] + meanSlopeZ * dX;
  const double chi2Y = sums.y[S::WRR] + dY * dY * sums.y[S::W] + dX * dX * sums.y[S::WSS] - 2. * dY * sums.y[S::WR] + 2. * dX * sums.y[S::WSR] - 2. * dX * dY * sums.y[S::WS];
  const double chi2Z = sums.z[S::WRR] + dZ * dZ * sums.z[S::W] + dX * dX * sums.z[S::WSS] - 2. * dZ * sums.z[S::WR] + 2. * dX * sums.z[S::WSR] - 2. * dX * dZ * sums.z[S::WS];
  const int ndf = int(2. * sums.entries + 0.5) - 3;
  const double chi2 = std::max(0., chi2Y) + std::max(0., chi2Z);
  const double scale = ndf > 0 ? std::max(1., chi2 / ndf) : 1.;
  solution.value = {dX, dY, dZ};
  solution.error[0] = std::sqrt(std::max(0., covarianceXX) * scale);
  solution.error[1] = std::sqrt(std::max(0., covarianceYY) * scale);
  solution.error[2] = std::sqrt(std::max(0., 1. / sums.z[S::W] + meanSlopeZ * meanSlopeZ * covarianceXX) * scale);
  solution.correlationXY = covarianceXX * covarianceYY > 0. ? covarianceXY / std::sqrt(covarianceXX * covarianceYY) : 0.;
  solution.rmsY = std::sqrt(std::max(0., chi2Y) / sums.y[S::W]);
  solution.rmsZ = std::sqrt(std::max(0., chi2Z) / sums.z[S::W]);
  solution.chi2Scale = scale;
  return std::isfinite(dX) && std::isfinite(dY) && std::isfinite(dZ);
}

/// one voxel with data: its chi2 theta^T a theta - 2 b^T theta + const (theta = (dX, dY, dZ)) at its position, and the nodes there
struct DataVoxel {
  Eigen::Matrix3d a;
  Eigen::Vector3d b;
  double weight = 1.; ///< robust weight
  int nNodes = 0;
  std::array<int, 4> node{}; ///< sector-local
  std::array<double, 4> nodeWeight{};
};

/// the per-coordinate noise s of a voxel with N entries, s / sqrt(N): the rows' second difference d_r - (d_{r-1} + d_{r+1}) / 2 over its
/// expected spread sqrt(1/N_r + (1/N_{r-1} + 1/N_{r+1}) / 4), its MAD in the IROC bulk (rows 10-55)
std::array<double, 3> estimateNoise(const TrackResiduals& binning, const std::vector<TrackResiduals::VoxRes>& results, double minEntries)
{
  const int voxelsPerSector = binning.getNVoxelsPerSector();
  std::array<std::vector<double>, 3> scaled;
  for (int sector = 0; sector < NSectors; ++sector) {
    for (int binZ2X = 0; binZ2X < binning.getNZ2XBins(); ++binZ2X) {
      for (int binY2X = 0; binY2X < binning.getNY2XBins(); ++binY2X) {
        for (int row = 10; row <= 55; ++row) {
          const size_t voxel = size_t(sector) * voxelsPerSector + binning.getGlbVoxBin(row, binY2X, binZ2X);
          const auto& previous = results[voxel - 1];
          const auto& current = results[voxel];
          const auto& next = results[voxel + 1];
          const double entriesPrevious = previous.stat[TrackResiduals::VoxV];
          const double entriesCurrent = current.stat[TrackResiduals::VoxV];
          const double entriesNext = next.stat[TrackResiduals::VoxV];
          if (entriesPrevious < minEntries || entriesCurrent < minEntries || entriesNext < minEntries) {
            continue;
          }
          const double spread = std::sqrt(1. / entriesCurrent + 0.25 * (1. / entriesPrevious + 1. / entriesNext));
          for (int coordinate = 0; coordinate < 3; ++coordinate) {
            scaled[coordinate].push_back((current.D[coordinate] - 0.5 * (previous.D[coordinate] + next.D[coordinate])) / spread);
          }
        }
      }
    }
  }
  auto median = [](std::vector<double>& values) {
    if (values.empty()) {
      return 0.;
    }
    std::nth_element(values.begin(), values.begin() + values.size() / 2, values.end());
    return values[values.size() / 2];
  };
  std::array<double, 3> noise{};
  for (int coordinate = 0; coordinate < 3; ++coordinate) {
    const double centre = median(scaled[coordinate]);
    for (auto& value : scaled[coordinate]) {
      value = std::abs(value - centre);
    }
    noise[coordinate] = 1.4826 * median(scaled[coordinate]);
  }
  return noise;
}

/// per voxel of the grid: 1 if it is one of the innermost IROC rows of its column below trustFraction x the column's plateau entries
std::vector<char> untrustedRows(const TrackResiduals& binning, const std::vector<TrackResiduals::VoxRes>& results, double trustFraction, long& nUntrusted)
{
  const int voxelsPerSector = binning.getNVoxelsPerSector();
  std::vector<char> untrusted(results.size(), 0);
  nUntrusted = 0;
  if (trustFraction <= 0.) {
    return untrusted;
  }
  std::vector<double> plateau;
  for (int sector = 0; sector < NSectors; ++sector) {
    for (int binZ2X = 0; binZ2X < binning.getNZ2XBins(); ++binZ2X) {
      for (int binY2X = 0; binY2X < binning.getNY2XBins(); ++binY2X) {
        const size_t first = size_t(sector) * voxelsPerSector + binning.getGlbVoxBin(0, binY2X, binZ2X);
        plateau.clear();
        for (int row = PlateauFirstRow; row <= PlateauLastRow; ++row) {
          plateau.push_back(results[first + row].stat[TrackResiduals::VoxV]);
        }
        std::nth_element(plateau.begin(), plateau.begin() + plateau.size() / 2, plateau.end());
        const double threshold = trustFraction * plateau[plateau.size() / 2];
        for (int row = 0; row < NIrocRows && results[first + row].stat[TrackResiduals::VoxV] < threshold; ++row) {
          untrusted[first + row] = 1;
          ++nUntrusted;
        }
      }
    }
  }
  return untrusted;
}

/// one term of the penalty: lambda (sum_k coefficient_k node(first + k step))^2 along an axis
struct Stencil {
  int axis = 0;
  int first = 0;
  std::array<double, 3> coefficient{};
  double lambda = 0.;
};

/// the penalty's terms, the same in every sector (none across the stack boundaries). Per bin (curvature = false): the second difference
/// with the prior rms tau[axis] (cm) per term. Curvature: the second derivative from the bin centres' positions in cm (x of the pad rows,
/// y = y/x x, z = z/x x; unequal spacing allowed) weighted with the centre voxel's volume in cm^3, tau[axis] = the prior rms (cm / cm^2)
/// of the curvature over 1 cm^3: the penalty ~ the integral of curvature^2 over the volume, the same prior for any binning and radius
std::vector<Stencil> makeStencils(const TrackResiduals& binning, const std::array<double, 3>& tau, bool curvature)
{
  const int nBinsX = binning.getNXBins();
  const int nBinsY2X = binning.getNY2XBins();
  const int nBinsZ2X = binning.getNZ2XBins();
  const int bins[3] = {nBinsX, nBinsY2X, nBinsZ2X};
  std::vector<Stencil> stencils;
  for (int axis = 0; axis < 3; ++axis) {
    if (!(tau[axis] > 0.)) {
      continue;
    }
    const double lambda = 1. / (tau[axis] * tau[axis]);
    for (int binZ2X = 0; binZ2X < nBinsZ2X; ++binZ2X) {
      for (int binY2X = 0; binY2X < nBinsY2X; ++binY2X) {
        for (int binX = 0; binX < nBinsX; ++binX) {
          const int bin[3] = {binX, binY2X, binZ2X};
          if (bin[axis] + 3 > bins[axis] || (axis == 0 && stackOfRow(binX) != stackOfRow(binX + 2))) {
            continue;
          }
          Stencil stencil;
          stencil.axis = axis;
          stencil.first = int(binning.getGlbVoxBin(binX, binY2X, binZ2X));
          if (!curvature) {
            stencil.coefficient = {1., -2., 1.};
            stencil.lambda = lambda;
          } else {
            // the positions along the axis of the three bin centres (y/x and z/x in the pad row of the stencil)
            std::array<double, 3> position{};
            for (int k = 0; k < 3; ++k) {
              const int row = axis == 0 ? binX + k : binX;
              position[k] = axis == 0 ? binning.getX(row)
                                      : binning.getX(row) * (axis == 1 ? binning.getY2X(row, binY2X + k) : binning.getZ2X(binZ2X + k));
            }
            const double step1 = position[1] - position[0];
            const double step2 = position[2] - position[1];
            const double scale = 2. / (step1 + step2);
            stencil.coefficient = {scale / step1, -scale * (1. / step1 + 1. / step2), scale / step2};
            // the centre voxel's volume
            const int centre[3] = {binX + (axis == 0), binY2X + (axis == 1), binZ2X + (axis == 2)};
            const double x = binning.getX(centre[0]);
            const double volume = (1. / binning.getDXI(centre[0])) * (x / binning.getDY2XI(centre[0], centre[1])) * (x / binning.getDZ2XI(centre[2]));
            stencil.lambda = lambda * volume;
          }
          stencils.push_back(stencil);
        }
      }
    }
  }
  return stencils;
}

/// the penalised fit per sector (see the file description); nodes: 3 per voxel of the grid (dX, dY, dZ); false if not finite
bool smoothNodes(const TrackResiduals& binning, std::vector<std::vector<DataVoxel>>& data, const std::array<double, 3>& tau, bool curvature,
                 double tau0, int robust, int nThreads, std::vector<double>& nodes)
{
  const int nBinsX = binning.getNXBins();
  const int nBinsY2X = binning.getNY2XBins();
  const int voxelsPerSector = binning.getNVoxelsPerSector();
  const int strides[3] = {1, nBinsX, nBinsX * nBinsY2X}; // getGlbVoxBin = row + nX * (y/x bin + nY2X * z/x bin)
  const double ridge = tau0 > 0. ? 1. / (tau0 * tau0) : 0.;

  // the penalty's terms and their diagonal per node
  const auto stencils = makeStencils(binning, tau, curvature);
  std::vector<double> penaltyDiagonal(voxelsPerSector, 0.);
  for (const auto& stencil : stencils) {
    for (int k = 0; k < 3; ++k) {
      penaltyDiagonal[stencil.first + k * strides[stencil.axis]] += stencil.lambda * stencil.coefficient[k] * stencil.coefficient[k];
    }
  }

  nodes.assign(size_t(NSectors) * voxelsPerSector * 3, 0.);
  std::array<int, NSectors> iterations{};
  std::array<double, NSectors> residuals{};
  std::array<bool, NSectors> finite{};
  std::array<long, NSectors> nLowWeight{};
  std::array<long, NSectors> nZeroWeight{};
  auto smoothSector = [&](int sector) {
    finite[sector] = true;
    auto& voxels = data[sector];
    if (voxels.empty()) {
      return;
    }
    const int n = 3 * voxelsPerSector;
    double* sectorNodes = nodes.data() + 3 * size_t(sector) * voxelsPerSector;
    auto nodeValue = [&](const DataVoxel& voxel, const double* x) {
      Eigen::Vector3d theta = Eigen::Vector3d::Zero();
      for (int j = 0; j < voxel.nNodes; ++j) {
        theta += voxel.nodeWeight[j] * Eigen::Map<const Eigen::Vector3d>(x + 3 * voxel.node[j]);
      }
      return theta;
    };
    for (int cycle = 0; cycle <= robust; ++cycle) {
      if (cycle > 0) {
        // Tukey's bisquare of sqrt(excess chi2 / rank): the voxel's chi2 at the solution minus its minimum -b^T a^+ b
        nLowWeight[sector] = 0;
        nZeroWeight[sector] = 0;
        for (auto& voxel : voxels) {
          const Eigen::Vector3d theta = nodeValue(voxel, sectorNodes);
          Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(voxel.a);
          const double largest = eigen.eigenvalues().maxCoeff();
          double minimum = 0.;
          int rank = 0;
          for (int k = 0; k < 3; ++k) {
            const double value = eigen.eigenvalues()(k);
            if (value > 1e-12 * largest) {
              const double projection = eigen.eigenvectors().col(k).dot(voxel.b);
              minimum -= projection * projection / value;
              ++rank;
            }
          }
          const double excess = theta.dot(voxel.a * theta) - 2. * voxel.b.dot(theta) - minimum;
          const double normalised = rank > 0 ? std::sqrt(std::max(0., excess) / rank) / Tukey : 0.;
          voxel.weight = normalised < 1. ? (1. - normalised * normalised) * (1. - normalised * normalised) : 0.;
          nLowWeight[sector] += voxel.weight < 0.5;
          nZeroWeight[sector] += voxel.weight == 0.;
        }
      }
      // the normal equations: ridge + penalty + sum over the voxels of w_j w_k a at their nodes; preconditioner = the inverse 3 x 3 diagonal blocks
      Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
      std::vector<Eigen::Matrix3d> block(voxelsPerSector);
      for (int l = 0; l < voxelsPerSector; ++l) {
        block[l] = (penaltyDiagonal[l] + ridge) * Eigen::Matrix3d::Identity();
      }
      std::vector<Eigen::Matrix3d> weighted(voxels.size());
      for (size_t i = 0; i < voxels.size(); ++i) {
        const auto& voxel = voxels[i];
        weighted[i] = voxel.weight * voxel.a;
        const Eigen::Vector3d r = voxel.weight * voxel.b;
        for (int j = 0; j < voxel.nNodes; ++j) {
          rhs.segment<3>(3 * voxel.node[j]) += voxel.nodeWeight[j] * r;
          block[voxel.node[j]] += voxel.nodeWeight[j] * voxel.nodeWeight[j] * weighted[i];
        }
      }
      for (auto& matrix : block) {
        // a node without data and penalty (no stencil, no ridge) stays 0
        matrix = matrix.determinant() != 0. ? matrix.inverse().eval() : Eigen::Matrix3d::Zero();
      }
      auto apply = [&](const Eigen::VectorXd& x, Eigen::VectorXd& y) {
        y = ridge * x;
        for (const auto& stencil : stencils) {
          const int step = strides[stencil.axis];
          for (int c = 0; c < 3; ++c) {
            double difference = 0.;
            for (int k = 0; k < 3; ++k) {
              difference += stencil.coefficient[k] * x(3 * (stencil.first + k * step) + c);
            }
            for (int k = 0; k < 3; ++k) {
              y(3 * (stencil.first + k * step) + c) += stencil.lambda * stencil.coefficient[k] * difference;
            }
          }
        }
        for (size_t i = 0; i < voxels.size(); ++i) {
          const auto& voxel = voxels[i];
          const Eigen::Vector3d u = weighted[i] * nodeValue(voxel, x.data());
          for (int j = 0; j < voxel.nNodes; ++j) {
            y.segment<3>(3 * voxel.node[j]) += voxel.nodeWeight[j] * u;
          }
        }
      };
      auto precondition = [&](const Eigen::VectorXd& r, Eigen::VectorXd& z) {
        for (int l = 0; l < voxelsPerSector; ++l) {
          z.segment<3>(3 * l) = block[l] * r.segment<3>(3 * l);
        }
      };
      // preconditioned conjugate gradients from the previous cycle's solution
      Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(sectorNodes, n);
      Eigen::VectorXd r(n);
      Eigen::VectorXd z(n);
      Eigen::VectorXd q(n);
      const double rhsNorm = rhs.norm();
      iterations[sector] = 0;
      residuals[sector] = 0.;
      if (rhsNorm == 0.) {
        std::fill(sectorNodes, sectorNodes + n, 0.);
        continue;
      }
      apply(x, q);
      r = rhs - q;
      precondition(r, z);
      Eigen::VectorXd p = z;
      double rz = r.dot(z);
      residuals[sector] = r.norm() / rhsNorm;
      while (residuals[sector] > CGTolerance && iterations[sector] < CGMaxIterations) {
        apply(p, q);
        const double alpha = rz / p.dot(q);
        x += alpha * p;
        r -= alpha * q;
        residuals[sector] = r.norm() / rhsNorm;
        ++iterations[sector];
        precondition(r, z);
        const double rzNew = r.dot(z);
        p = z + (rzNew / rz) * p;
        rz = rzNew;
      }
      finite[sector] = finite[sector] && x.allFinite();
      Eigen::Map<Eigen::VectorXd>(sectorNodes, n) = x;
    }
  };
  std::atomic<int> nextSector{0};
  std::vector<std::thread> threads;
  for (int thread = 0; thread < std::clamp(nThreads, 1, NSectors); ++thread) {
    threads.emplace_back([&]() {
      for (int sector = nextSector++; sector < NSectors; sector = nextSector++) {
        smoothSector(sector);
      }
    });
  }
  for (auto& thread : threads) {
    thread.join();
  }

  const int nNotConverged = std::count_if(residuals.begin(), residuals.end(), [](double residual) { return residual > CGTolerance; });
  const int nNotFinite = std::count(finite.begin(), finite.end(), false);
  LOGP(info, "Smoothing: {} stencils per sector; CG iterations per sector {} - {}, worst relative residual {:.2e}", stencils.size(),
       *std::min_element(iterations.begin(), iterations.end()), *std::max_element(iterations.begin(), iterations.end()),
       *std::max_element(residuals.begin(), residuals.end()));
  if (robust > 0) {
    LOGP(info, "Smoothing: {} Tukey reweightings (c = {}), at the last {} voxels below weight 0.5, {} at 0", robust, Tukey,
         std::accumulate(nLowWeight.begin(), nLowWeight.end(), 0L), std::accumulate(nZeroWeight.begin(), nZeroWeight.end(), 0L));
  }
  if (nNotConverged) {
    LOGP(warning, "Smoothing: the conjugate gradients did not reach the tolerance {} in {} sectors", CGTolerance, nNotConverged);
  }
  if (nNotFinite) {
    LOGP(error, "Smoothing: the solution is not finite in {} sectors", nNotFinite);
  }
  return nNotFinite == 0;
}
} // namespace voxel_smoother

void smoothVoxelMap(std::string inputFile = "voxRes.root",
                    std::string outputFile = "voxResSmooth.root",
                    float tauRow = 100.,       // prior rms (um) of the second difference along the pad rows (0: no penalty); curvaturePrior: um / cm^2
                    float tauY2X = 1000.,      // ... along y/x
                    float tauZ2X = 1000.,      // ... along z/x
                    float tau0 = 0.,           // prior rms (um) of a node: ridge towards 0 (0: none)
                    float trustFraction = 0.7, // the innermost IROC rows below this fraction of the column's plateau entries are no data (<= 0: off)
                    float minEntries = 20.,    // voxResTree: voxels with fewer entries are no data; fineSums: the voxel solve's minimum
                    int errorMode = -1,        // voxResTree: -1 auto, 0 the input's E and EXYCorr, 1 s / sqrt(entries) (see the file description)
                    int joinIDCs = -1,         // meanIDC / medianIDC from the CCDB scalers at the input's TFs: -1 if the input has none, 1 always, 0 never
                    int robust = 0,            // Tukey reweighting cycles
                    int nThreads = 8,
                    bool curvaturePrior = false, // the taus as the prior rms of the curvature over 1 cm^3 (um / cm^2): the same for any binning (makeStencils)
                    std::string source = "auto") // "voxResTree", "fineSums" or "auto" (fineSums if present)
{
  using namespace voxel_smoother;
  std::unique_ptr<TFile> input(TFile::Open(inputFile.data()));
  if (!input || input->IsZombie()) {
    LOGP(error, "Cannot open {}", inputFile);
    return;
  }
  auto* voxResTree = input->Get<TTree>("voxResTree");
  auto* sumsTree = input->Get<TTree>("fineSums");
  if (source == "auto") {
    source = sumsTree ? "fineSums" : "voxResTree";
  }
  const bool useSums = source == "fineSums";
  if ((useSums && !sumsTree) || (!useSums && (source != "voxResTree" || !voxResTree))) {
    LOGP(error, "{}: no input {} (source: voxResTree, fineSums or auto)", inputFile, source);
    return;
  }
  TTree* infoTree = useSums ? sumsTree : voxResTree;
  auto info = [&](const char* key) -> std::string {
    auto* object = infoTree->GetUserInfo()->FindObject(key);
    return object ? object->GetTitle() : "";
  };
  if (info("y2xBinning").empty() || info("z2xBinning").empty() || info("maxZ2X").empty()) {
    LOGP(error, "{}: {} without the binning in its user info (y2xBinning, z2xBinning, maxZ2X)", inputFile, source);
    return;
  }
  const auto binning = makeBinning(info("y2xBinning"), info("z2xBinning"), std::stof(info("maxZ2X")));
  const int voxelsPerSector = binning->getNVoxelsPerSector();
  const size_t nVoxels = size_t(NSectors) * voxelsPerSector;

  // the output's voxels: the input's (voxResTree) or the voxel solve (fineSums), positions at the centres where empty
  std::vector<TrackResiduals::VoxRes> results(nVoxels);
  for (int sector = 0; sector < NSectors; ++sector) {
    for (int binZ2X = 0; binZ2X < binning->getNZ2XBins(); ++binZ2X) {
      for (int binY2X = 0; binY2X < binning->getNY2XBins(); ++binY2X) {
        for (int binX = 0; binX < binning->getNXBins(); ++binX) {
          auto& result = results[size_t(sector) * voxelsPerSector + binning->getGlbVoxBin(binX, binY2X, binZ2X)];
          float centreX = 0.f;
          float centreY2X = 0.f;
          float centreZ2X = 0.f;
          binning->getVoxelCoordinates(sector, binX, binY2X, binZ2X, centreX, centreY2X, centreZ2X);
          result.stat[TrackResiduals::VoxX] = centreX;
          result.stat[TrackResiduals::VoxF] = centreY2X;
          result.stat[TrackResiduals::VoxZ] = centreZ2X;
          result.bvox[TrackResiduals::VoxZ] = binZ2X;
          result.bvox[TrackResiduals::VoxF] = binY2X;
          result.bvox[TrackResiduals::VoxX] = binX;
          result.bsec = sector;
        }
      }
    }
  }
  std::vector<std::array<ChargeSums, 2>> sums;
  if (useSums) {
    sums.resize(nVoxels);
    int voxelIndex = 0;
    float entries[2] = {};
    float position[2][3] = {};
    float channelY[2][ChargeSums::NSums] = {};
    float channelZ[2][ChargeSums::NSums] = {};
    sumsTree->SetBranchAddress("voxel", &voxelIndex);
    sumsTree->SetBranchAddress("entries", entries);
    sumsTree->SetBranchAddress("position", position);
    sumsTree->SetBranchAddress("y", channelY);
    sumsTree->SetBranchAddress("z", channelZ);
    for (long entry = 0; entry < sumsTree->GetEntries(); ++entry) {
      sumsTree->GetEntry(entry);
      if (voxelIndex < 0 || size_t(voxelIndex) >= nVoxels) {
        LOGP(error, "{}: fineSums voxel {} outside the grid", inputFile, voxelIndex);
        return;
      }
      for (int charge = 0; charge < 2; ++charge) {
        auto& chargeSums = sums[voxelIndex][charge];
        chargeSums.entries = entries[charge];
        chargeSums.x = position[charge][0];
        chargeSums.y2x = position[charge][1];
        chargeSums.z2x = position[charge][2];
        for (int sum = 0; sum < ChargeSums::NSums; ++sum) {
          chargeSums.y[sum] = channelY[charge][sum];
          chargeSums.z[sum] = channelZ[charge][sum];
        }
      }
    }
    sumsTree->ResetBranchAddresses();
    for (size_t voxel = 0; voxel < nVoxels; ++voxel) {
      ChargeSums total = sums[voxel][0];
      for (int sum = 0; sum < ChargeSums::NSums; ++sum) {
        total.y[sum] += sums[voxel][1].y[sum];
        total.z[sum] += sums[voxel][1].z[sum];
      }
      total.entries += sums[voxel][1].entries;
      total.x += sums[voxel][1].x;
      total.y2x += sums[voxel][1].y2x;
      total.z2x += sums[voxel][1].z2x;
      sums[voxel][0] = total; // the charges are not needed separately below
      auto& result = results[voxel];
      result.stat[TrackResiduals::VoxV] = total.entries;
      if (total.entries <= 0.) {
        continue;
      }
      result.stat[TrackResiduals::VoxX] = total.x / total.entries;
      result.stat[TrackResiduals::VoxF] = total.y2x / total.entries;
      result.stat[TrackResiduals::VoxZ] = total.z2x / total.entries;
      VoxelSolution solution;
      if (solveVoxel(total, minEntries, solution)) {
        for (int coordinate = 0; coordinate < 3; ++coordinate) {
          result.D[coordinate] = solution.value[coordinate];
          result.E[coordinate] = solution.error[coordinate];
        }
        result.EXYCorr = solution.correlationXY;
        result.dYSigMAD = solution.rmsY;
        result.dZSigLTM = solution.rmsZ;
        result.D[TrackResiduals::ResD] = solution.rmsY;
        result.flags = TrackResiduals::DistDone | TrackResiduals::DispDone;
      }
    }
  } else {
    TrackResiduals::VoxRes* result = nullptr;
    voxResTree->SetBranchAddress("voxRes", &result);
    for (long entry = 0; entry < voxResTree->GetEntries(); ++entry) {
      voxResTree->GetEntry(entry);
      results[size_t(result->bsec) * voxelsPerSector + binning->getGlbVoxBin(result->bvox)] = *result;
    }
    voxResTree->ResetBranchAddresses();
  }

  // the errors of a voxResTree: the input's, or s / sqrt(entries) (CPM placeholders)
  if (!useSums) {
    if (errorMode < 0) {
      bool placeholders = true;
      for (const auto& result : results) {
        if (result.stat[TrackResiduals::VoxV] > 0.f && (result.E[0] != PlaceholderError || result.E[1] != PlaceholderError || result.E[2] != PlaceholderError || result.EXYCorr != 1.f)) {
          placeholders = false;
          break;
        }
      }
      errorMode = placeholders ? 1 : 0;
    }
    if (errorMode == 1) {
      const auto noise = estimateNoise(*binning, results, minEntries);
      LOGP(info, "Errors: s / sqrt(entries), s = {:.4f} / {:.4f} / {:.4f} cm (dX / dY / dZ, from the rows' scatter in the IROC bulk)", noise[0], noise[1], noise[2]);
      for (auto& result : results) {
        const double entries = result.stat[TrackResiduals::VoxV];
        for (int coordinate = 0; coordinate < 3; ++coordinate) {
          result.E[coordinate] = entries > 0. ? noise[coordinate] / std::sqrt(entries) : 0.;
        }
        result.EXYCorr = 0.f;
      }
    } else {
      LOGP(info, "Errors: the input's E and EXYCorr");
    }
  }

  // the data per sector
  long nUntrusted = 0;
  const auto untrusted = untrustedRows(*binning, results, trustFraction, nUntrusted);
  double medianScale = 1.;
  if (useSums) {
    std::vector<double> scales;
    for (size_t voxel = 0; voxel < nVoxels; ++voxel) {
      VoxelSolution solution;
      if (solveVoxel(sums[voxel][0], minEntries, solution)) {
        scales.push_back(solution.chi2Scale);
      }
    }
    if (!scales.empty()) {
      std::nth_element(scales.begin(), scales.begin() + scales.size() / 2, scales.end());
      medianScale = scales[scales.size() / 2];
    }
  }
  std::vector<std::vector<DataVoxel>> data(NSectors);
  long nData = 0;
  for (size_t voxel = 0; voxel < nVoxels; ++voxel) {
    if (untrusted[voxel]) {
      continue;
    }
    const auto& result = results[voxel];
    const int sector = voxel / voxelsPerSector;
    const int row = result.bvox[TrackResiduals::VoxX];
    DataVoxel dataVoxel;
    double z2xSide = 0.;
    if (useSums) {
      using S = ChargeSums;
      const auto& total = sums[voxel][0];
      if (total.entries <= 0. || total.y[S::W] <= 0. || total.z[S::W] <= 0.) {
        continue;
      }
      VoxelSolution solution;
      const double scale = solveVoxel(total, minEntries, solution) ? solution.chi2Scale : medianScale;
      // the y channel (J = (-s_y, 1, 0)) and the z channel's mean W_z (<r_z> - dZ + <s_z> dX)^2
      const double meanSlopeZ = total.z[S::WS] / total.z[S::W];
      dataVoxel.a << total.y[S::WSS] + total.z[S::WS] * meanSlopeZ, -total.y[S::WS], -total.z[S::WS], -total.y[S::WS], total.y[S::W], 0., -total.z[S::WS], 0.,
        total.z[S::W];
      dataVoxel.b << -total.y[S::WSR] - meanSlopeZ * total.z[S::WR], total.y[S::WR], total.z[S::WR];
      dataVoxel.a /= scale;
      dataVoxel.b /= scale;
      z2xSide = (sector < NSectors / 2 ? 1. : -1.) * result.stat[TrackResiduals::VoxZ];
    } else {
      const bool fitted = (result.flags & TrackResiduals::DistDone) && !(result.flags & TrackResiduals::Masked) && result.stat[TrackResiduals::VoxV] >= std::max(1.f, minEntries);
      if (!fitted || !(result.E[0] > 0.f) || !(result.E[1] > 0.f) || !(result.E[2] > 0.f) ||
          !std::isfinite(result.D[0]) || !std::isfinite(result.D[1]) || !std::isfinite(result.D[2])) {
        continue;
      }
      Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
      covariance(0, 0) = result.E[0] * result.E[0];
      covariance(1, 1) = result.E[1] * result.E[1];
      covariance(0, 1) = std::clamp(double(result.EXYCorr), -0.99, 0.99) * result.E[0] * result.E[1];
      covariance(1, 0) = covariance(0, 1);
      covariance(2, 2) = result.E[2] * result.E[2];
      dataVoxel.a = covariance.inverse();
      dataVoxel.b = dataVoxel.a * Eigen::Vector3d(result.D[0], result.D[1], result.D[2]);
      z2xSide = std::abs(result.stat[TrackResiduals::VoxZ]);
    }
    int nodeIndices[4];
    double nodeWeights[4];
    dataVoxel.nNodes = bilinearNodes(*binning, row, result.stat[TrackResiduals::VoxF], z2xSide, nodeIndices, nodeWeights);
    for (int j = 0; j < dataVoxel.nNodes; ++j) {
      dataVoxel.node[j] = nodeIndices[j];
      dataVoxel.nodeWeight[j] = nodeWeights[j];
    }
    data[sector].push_back(dataVoxel);
    ++nData;
  }
  const std::string description = fmt::format(
    "penalised node fit (smoothVoxelMap.C): input {}, tau rows / y/x / z/x {} / {} / {} {}, tau0 {} um, trustFraction {} "
    "({} untrusted inner IROC voxels), minEntries {}, errors {}, robust {}",
    source, tauRow, tauY2X, tauZ2X, curvaturePrior ? "um / cm^2 (curvature over 1 cm^3)" : "um (second difference per bin)", tau0, trustFraction, nUntrusted, minEntries,
    useSums ? std::string("voxel solve") : (errorMode == 1 ? std::string("s / sqrt(entries)") : std::string("input")), robust);
  LOGP(info, "{}: {} of {} voxels are data; {}", inputFile, nData, nVoxels, description);

  // the IDC join (staticMapCreatorCPM.C leaves meanIDC / medianIDC 0)
  auto* treeWithInfo = voxResTree ? voxResTree : sumsTree;
  auto* storedMeanIDC = treeWithInfo->GetUserInfo()->FindObject("meanIDC");
  const bool hasIDC = storedMeanIDC && std::atof(storedMeanIDC->GetTitle()) > 0.;
  float meanIDC = 0.f;
  float medianIDC = 0.f;
  const bool joinedIDC = (joinIDCs > 0 || (joinIDCs < 0 && !hasIDC)) && joinIDC(*input, meanIDC, medianIDC);
  if (!joinedIDC && !hasIDC) {
    LOGP(warning, "No meanIDC in the input and none joined: the output's meanIDC / medianIDC stay as in the input");
  }

  std::vector<double> nodes;
  const bool finite = smoothNodes(*binning, data, {tauRow * Micron, tauY2X * Micron, tauZ2X * Micron}, curvaturePrior, tau0 * Micron, robust, nThreads, nodes);
  if (!finite) {
    return;
  }

  // the output: the voxels with DS = the nodes, the input's user info and other trees
  TFile output(outputFile.data(), "recreate");
  TTree tree("voxResTree", "voxRes map results and statistics, smoothed by smoothVoxelMap.C");
  TrackResiduals::VoxRes result;
  TrackResiduals::VoxRes* resultPointer = &result;
  tree.Branch("voxRes", &resultPointer);
  for (auto* object : *treeWithInfo->GetUserInfo()) {
    const std::string name = object->GetName();
    if (!joinedIDC || (name != "meanIDC" && name != "medianIDC")) {
      tree.GetUserInfo()->Add(object->Clone());
    }
  }
  if (joinedIDC) {
    tree.GetUserInfo()->Add(new TNamed("meanIDC", std::to_string(meanIDC).data()));
    tree.GetUserInfo()->Add(new TNamed("medianIDC", std::to_string(medianIDC).data()));
  }
  if (!voxResTree) {
    tree.GetUserInfo()->Add(new TNamed("nTracksProcessed", info("nTracks").data()));
  }
  tree.GetUserInfo()->Add(new TNamed("smoothing", description.data()));
  for (const char* alias : {"z2xBin:bvox[0]", "y2xBin:bvox[1]", "xBin:bvox[2]", "z2xAV:stat[0]", "y2xAV:stat[1]", "xAV:stat[2]", "entries:stat[3]",
                            "dX:D[0]", "dY:D[1]", "dZ:D[2]", "dXS:DS[0]", "dYS:DS[1]", "dZS:DS[2]", "dXE:E[0]", "dYE:E[1]", "dZE:E[2]"}) {
    const std::string pair(alias);
    tree.SetAlias(pair.substr(0, pair.find(':')).data(), pair.substr(pair.find(':') + 1).data());
  }
  for (int sector = 0; sector < NSectors; ++sector) {
    for (int binZ2X = 0; binZ2X < binning->getNZ2XBins(); ++binZ2X) {
      for (int binY2X = 0; binY2X < binning->getNY2XBins(); ++binY2X) {
        for (int binX = 0; binX < binning->getNXBins(); ++binX) {
          const size_t voxel = size_t(sector) * voxelsPerSector + binning->getGlbVoxBin(binX, binY2X, binZ2X);
          result = results[voxel];
          for (int coordinate = 0; coordinate < 3; ++coordinate) {
            result.DS[coordinate] = nodes[3 * voxel + coordinate];
          }
          result.flags |= TrackResiduals::SmoothDone;
          result.stat[TrackResiduals::VoxV] = std::max(result.stat[TrackResiduals::VoxV], 1.f);
          tree.Fill();
        }
      }
    }
  }
  tree.Write();
  for (auto* key : *input->GetListOfKeys()) {
    auto* object = static_cast<TKey*>(key)->ReadObj();
    const std::string name = object->GetName();
    if (object->InheritsFrom(TTree::Class()) && name != "voxResTree" && name != "fineSums") {
      output.cd();
      static_cast<TTree*>(object)->CloneTree(-1, "fast")->Write();
    }
  }
  LOGP(info, "{} written: {} voxels, DS = the smoothed nodes", outputFile, nVoxels);
}
