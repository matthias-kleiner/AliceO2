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
/// SmoothingExtrapolate.C. Defaults tuned on a PbPb 2025 map (568210, 22 kHz) by a leave-rows-out test (hide the first measured inner
/// rows, predict them) and a hold-out cross-validation of the bulk (hide 10 % of the voxels): the curvature prior 1000 / 2000 / 1000
/// um / cm^2 (rows / y/x / z/x) predicts within 2-3 % of the best scanned setting in both (hold-out rms 1.7-3.4 times smaller than
/// SmoothingExtrapolate.C's, 1.2-1.5 in the OROCs) and is ~10-100 times smoother in the extrapolated rows.
///
/// Input (source):
/// - "voxResTree": D with the errors E (dX, dY correlated by EXYCorr) at the voxel's mean position (stat), for the fitted, unmasked voxels
///   with >= minEntries entries (staticMapCreator.C, staticMapCreatorCPM.C, fineMapCreatorGBL.C). errorMode 0: the input's E, EXYCorr;
///   1: E = s / sqrt(entries), EXYCorr = 0, with s per coordinate from the scatter of the rows' second differences in the IROC bulk (an
///   upper bound: it includes the curvature differences between columns); 2: E, EXYCorr from the residual slope method's voxel solve of
///   the same file's fineSums (staticMapCreatorCPM.C voxelMethod 0 + writeSlopeSums: the CPM values weighted with the slope method's
///   errors; voxels it does not solve are no data); -1: 1 if the input has the placeholders of staticMapCreatorCPM.C's CPM maps (E = 0.1,
///   EXYCorr = 1), else 0;
/// - "fineSums": per voxel and charge the sums of the residual slope method (staticMapCreatorCPM.C writeSlopeSums, fineMapCreatorGBL.C
///   writeSums), solved here: D, E from the voxel solve; every voxel with entries enters with its chi2 scaled by max(1, chi2 / ndf) of its
///   solve (the median of the solved voxels below minEntries);
/// - "auto": voxResTree if the file has one, else fineSums (a file with both is the voxResTree's map: a CPM map with the slope sums,
///   or a two-pass fineMapCreatorGBL.C map, whose fineSums hold pass 1 only).
/// trustFraction > 0: per column the innermost IROC rows with fewer entries than trustFraction x the column's median entries in rows
/// 10-30 (up to the first row at or above it) are no data; -1 (default): 0.7 with the s / sqrt(entries) errors (errorMode 1: the large
/// radial distortion empties these rows and biases them, the stand-in errors do not weight them down: without the cut the map follows them
/// by up to ~3 mm), else 0 (with real errors they hardly pull: <= 0.3 mm in the rows with clusters, the bulk unchanged).
/// robust > 0: Tukey reweightings of the voxels by their excess chi2.
/// curvaturePrior (default): the taus are the prior rms of the curvature (um / cm^2) averaged over 1 cm^3: second derivatives from the
/// bin centres' positions, each term weighted with its voxel's volume, so the prior does not depend on the binning (coarse map, fine
/// pass, refined bins) nor on the radius (the y/x and z/x bins grow with x); false: the taus are the second difference per bin (um).
///
/// Output: voxResTree with every voxel of the grid, DS = the nodes interpolated bilinearly to the voxel's position stat for the voxels that
/// are data, the voxel's node and stat = the bin centre for the others, SmoothDone, entries at least 1 (TPCFastSpaceChargeCorrectionHelper
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
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

#ifndef SMOOTHVOXELMAP_SOLVER_ONLY // the solver alone (solveSector and its algebra) compiles without ROOT and O2
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
#endif
#else
#error This macro must run in compiled mode
#endif

namespace voxel_smoother
{
constexpr int NSectors = 36;
constexpr int NRows = 152;
constexpr int NIrocRows = 63;
constexpr int PlateauFirstRow = 10; ///< the rows whose median entries are a column's plateau (trustFraction)
constexpr int PlateauLastRow = 30;
constexpr double Micron = 1e-4;          ///< in cm
constexpr double Tukey = 4.685;          ///< Tukey's bisquare constant of the robust reweighting
constexpr double CGTolerance = 1e-8;     ///< relative residual norm of the conjugate gradients
constexpr int CGMaxIterations = 20000;   ///< per sector and cycle
constexpr int JacobiMaxSweeps = 50;      ///< of the 3 x 3 eigen-decomposition (robust reweighting of a voxel with a near-singular chi2)
constexpr float PlaceholderError = 0.1f; ///< staticMapCreatorCPM.C's E of every voxel

// the algebra of the fit with std::array / std::vector (no Eigen: it is not on the ROOT_INCLUDE_PATH of ACLiC on the GRID), in the operation
// order of Eigen 3.4 with SSE2 / NEON (packets of two doubles) where it is cheap
using Vector3 = std::array<double, 3>;
using Matrix3 = std::array<Vector3, 3>; ///< [row][column]

/// a^T b (from left to right as Eigen's Vector3d::dot)
double dot(const Vector3& a, const Vector3& b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/// matrix x vector: rows 0 and 1 from left to right, row 2 as m20 v0 + (m21 v1 + m22 v2) (Eigen's Matrix3d * Vector3d: a packet of rows 0
/// and 1, row 2 by its coefficient)
Vector3 multiply(const Matrix3& matrix, const Vector3& vector)
{
  Vector3 product;
  product[0] = matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2];
  product[1] = matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2];
  product[2] = matrix[2][0] * vector[0] + (matrix[2][1] * vector[1] + matrix[2][2] * vector[2]);
  return product;
}

/// the determinant (as Eigen's Matrix3d::determinant)
double determinant(const Matrix3& m)
{
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

/// the inverse of a symmetric matrix via the adjugate (as Eigen's Matrix3d::inverse: the cofactors of the first column give the
/// determinant; for a symmetric matrix Eigen's cofactors (i, j) and (j, i) are the same products)
Matrix3 inverseSymmetric(const Matrix3& m)
{
  const double cofactor00 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
  const double cofactor01 = m[2][1] * m[0][2] - m[2][2] * m[0][1];
  const double cofactor02 = m[0][1] * m[1][2] - m[0][2] * m[1][1];
  const double cofactor11 = m[2][2] * m[0][0] - m[2][0] * m[0][2];
  const double cofactor12 = m[0][2] * m[1][0] - m[0][0] * m[1][2];
  const double cofactor22 = m[0][0] * m[1][1] - m[0][1] * m[1][0];
  const double inverseDeterminant = 1. / (cofactor00 * m[0][0] + cofactor01 * m[1][0] + cofactor02 * m[2][0]);
  Matrix3 inverse;
  inverse[0] = {cofactor00 * inverseDeterminant, cofactor01 * inverseDeterminant, cofactor02 * inverseDeterminant};
  inverse[1] = {inverse[0][1], cofactor11 * inverseDeterminant, cofactor12 * inverseDeterminant};
  inverse[2] = {inverse[0][2], inverse[1][2], cofactor22 * inverseDeterminant};
  return inverse;
}

/// a^T b of two CG vectors in four interleaved partial sums (as Eigen's VectorXd::dot with two packets of two doubles)
double dot(const std::vector<double>& a, const std::vector<double>& b)
{
  const size_t size = a.size();
  const size_t end = size / 4 * 4;
  std::array<double, 4> partial{};
  for (size_t i = 0; i < end; i += 4) {
    for (size_t k = 0; k < 4; ++k) {
      partial[k] += a[i + k] * b[i + k];
    }
  }
  partial[0] += partial[2];
  partial[1] += partial[3];
  size_t i = end;
  if (size - end >= 2) {
    partial[0] += a[i] * b[i];
    partial[1] += a[i + 1] * b[i + 1];
    i += 2;
  }
  double sum = partial[0] + partial[1];
  for (; i < size; ++i) {
    sum += a[i] * b[i];
  }
  return sum;
}

/// the minimum -b^T a^+ b of theta^T a theta - 2 b^T theta for a symmetric positive semi-definite a, and the rank of a (its eigenvalues
/// above 1e-12 x the largest): with a^-1 if det(a) > 1e-12 trace(a)^3 (then all eigenvalues are above), else from the eigen-decomposition
/// by cyclic Jacobi rotations
double quadraticMinimum(const Matrix3& a, const Vector3& b, int& rank)
{
  const double trace = a[0][0] + a[1][1] + a[2][2];
  if (determinant(a) > 1e-12 * trace * trace * trace) {
    rank = 3;
    return -dot(b, multiply(inverseSymmetric(a), b));
  }
  Matrix3 rotated = a;    // V^T a V, diagonal at the end
  Matrix3 eigenvectors{}; // V: in its columns
  double squaredNorm = 0.;
  for (int row = 0; row < 3; ++row) {
    eigenvectors[row][row] = 1.;
    for (int column = 0; column < 3; ++column) {
      squaredNorm += a[row][column] * a[row][column];
    }
  }
  const double epsilon = std::numeric_limits<double>::epsilon();
  for (int sweep = 0; sweep < JacobiMaxSweeps; ++sweep) {
    const double offDiagonal = rotated[0][1] * rotated[0][1] + rotated[0][2] * rotated[0][2] + rotated[1][2] * rotated[1][2];
    if (!(offDiagonal > epsilon * epsilon * squaredNorm)) {
      break;
    }
    for (int first = 0; first < 2; ++first) {
      for (int second = first + 1; second < 3; ++second) {
        const double element = rotated[first][second];
        if (element == 0.) {
          continue;
        }
        // the rotation by atan(tangent) in the plane (first, second) that zeroes the element
        const int third = 3 - first - second;
        const double cotangentTwoAngle = (rotated[second][second] - rotated[first][first]) / (2. * element);
        const double tangent = (cotangentTwoAngle < 0. ? -1. : 1.) / (std::abs(cotangentTwoAngle) + std::hypot(cotangentTwoAngle, 1.));
        const double cosine = 1. / std::sqrt(tangent * tangent + 1.);
        const double sine = tangent * cosine;
        rotated[first][first] -= tangent * element;
        rotated[second][second] += tangent * element;
        rotated[first][second] = 0.;
        rotated[second][first] = 0.;
        const double thirdFirst = rotated[third][first];
        const double thirdSecond = rotated[third][second];
        rotated[third][first] = cosine * thirdFirst - sine * thirdSecond;
        rotated[first][third] = rotated[third][first];
        rotated[third][second] = sine * thirdFirst + cosine * thirdSecond;
        rotated[second][third] = rotated[third][second];
        for (int row = 0; row < 3; ++row) {
          const double vectorFirst = eigenvectors[row][first];
          const double vectorSecond = eigenvectors[row][second];
          eigenvectors[row][first] = cosine * vectorFirst - sine * vectorSecond;
          eigenvectors[row][second] = sine * vectorFirst + cosine * vectorSecond;
        }
      }
    }
  }
  const double largest = std::max({rotated[0][0], rotated[1][1], rotated[2][2]});
  double minimum = 0.;
  rank = 0;
  for (int k = 0; k < 3; ++k) {
    const double eigenvalue = rotated[k][k];
    if (eigenvalue > 1e-12 * largest) {
      const double projection = eigenvectors[0][k] * b[0] + eigenvectors[1][k] * b[1] + eigenvectors[2][k] * b[2];
      minimum -= projection * projection / eigenvalue;
      ++rank;
    }
  }
  return minimum;
}

/// one voxel with data: its chi2 theta^T a theta - 2 b^T theta + const (theta = (dX, dY, dZ)) at its position, and the nodes there
struct DataVoxel {
  Matrix3 a{};
  Vector3 b{};
  double weight = 1.; ///< robust weight
  int nNodes = 0;
  std::array<int, 4> node{}; ///< sector-local
  std::array<double, 4> nodeWeight{};
};

/// one term of the penalty: lambda (sum_k coefficient_k node(first + k step))^2 along an axis
struct Stencil {
  int axis = 0;
  int first = 0;
  std::array<double, 3> coefficient{};
  double lambda = 0.;
};

/// per sector: the solver statistics
struct SectorStatistics {
  int iterations = 0;
  double residual = 0.; ///< relative residual norm of the last CG
  bool finite = true;
  long nLowWeight = 0;  ///< robust: voxels below weight 0.5 at the last reweighting
  long nZeroWeight = 0; ///< ... at 0
};

// The preconditioner of solveSector: a multigrid V-cycle. The normal matrix couples the nodes within +-2 pad rows (row penalty), +-2 y/x
// bins (y/x penalty), +-2 z/x bins (z/x penalty) and the bilinear neighbours in y/x and z/x within a pad row (data); the row penalty is by
// far the stiffest term (the rows ~1 cm apart, the y/x and z/x bins several cm), then the y/x penalty. The smoother is a block
// Gauss-Seidel sweep over the columns (the nodes of one y/x and z/x bin) with the column's block of the matrix (the row penalty, the data
// at the nodes and the diagonal of the other terms) solved exactly by its block L D L^T: the rows are smoothed exactly. The coarse levels
// keep the rows and halve y/x and / or z/x (linear interpolation P; only the stiff direction(s), see coarseLevel) with the Galerkin
// matrix P^T A P, down to <= MultigridCoarsestColumns columns, solved by a banded Cholesky. The forward sweep before and the backward
// sweep after the coarse correction make it symmetric positive definite (a preconditioner for CG). Every level in the layout of the
// sector's nodes: node = row + nRows (y/x bin + nY2X z/x bin), 3 values (dX, dY, dZ) per node.
constexpr int MultigridCoarsestColumns = 4; ///< the coarsest level of the multigrid: at most this many columns
constexpr double MultigridStiffness = 0.25; ///< semi-coarsening: a direction whose penalty is >= this x the stiffer direction's is coarsened
constexpr int MaxNodeOffset = 2;            ///< the matrix couples the nodes within +-2 pad rows, y/x bins and z/x bins (on all levels)

/// the difference of two nodes' pad row, y/x bin and z/x bin
struct NodeOffset {
  int row = 0;
  int y2x = 0;
  int z2x = 0;
};

/// one level of the multigrid: the matrix on a grid of nRows x nY2X x nZ2X nodes (node = row + nRows (y/x bin + nY2X z/x bin), as in the
/// sector) per coupling offset (both directions and the diagonal) and node, the penalty part (the same for dX, dY, dZ) and the data part
/// (3 x 3); the block LDL^T of its columns (or, on the coarsest level, the banded Cholesky factor of the whole matrix); work vectors
struct MultigridLevel {
  int nRows = 0;
  int nY2X = 0;
  int nZ2X = 0;
  std::vector<NodeOffset> penaltyOffsets;
  std::vector<double> penalty; ///< [offset][node]: the coupling of the node to the node at the offset
  std::vector<NodeOffset> dataOffsets;
  std::vector<Matrix3> data;         ///< [offset][node]
  std::vector<Matrix3> pivotInverse; ///< per node: D^-1 of the column's L D L^T
  std::vector<Matrix3> lower1;       ///< per node: L(node, node 1 row before)
  int twoRowsBeforeSlot = -1;        ///< the penalty offset to the node 2 rows before
  int bandWidth = 0;                 ///< coarsest level: the band of the Cholesky factor (unknowns ordered by row, column, component)
  std::vector<double> band;          ///< [unknown][k]: L(unknown, unknown - k)
  std::vector<double> rhs;           ///< the right-hand side of the coarse correction (levels > 0)
  std::vector<double> solution;      ///< the coarse correction (levels > 0)
  std::vector<double> work;          ///< the residual (V-cycle) or the unknowns in the order of the band (coarsest)
  std::vector<double> columnRhs;     ///< the right-hand side of a column (Gauss-Seidel)

  int nNodes() const { return nRows * nY2X * nZ2X; }
  int nColumns() const { return nY2X * nZ2X; }
  int shift(const NodeOffset& offset) const { return offset.row + nRows * (offset.y2x + nY2X * offset.z2x); }
};

/// calls function(node, neighbour) for the nodes of the level whose neighbour at the offset is on the grid, along the rows
template <typename Function>
void forEachCoupling(const MultigridLevel& level, const NodeOffset& offset, Function function)
{
  const int shift = level.shift(offset);
  const int firstRow = std::max(0, -offset.row);
  const int endRow = std::min(level.nRows, level.nRows - offset.row);
  for (int z2x = std::max(0, -offset.z2x); z2x < std::min(level.nZ2X, level.nZ2X - offset.z2x); ++z2x) {
    for (int y2x = std::max(0, -offset.y2x); y2x < std::min(level.nY2X, level.nY2X - offset.y2x); ++y2x) {
      const int column = level.nRows * (y2x + level.nY2X * z2x);
      for (int node = column + firstRow; node < column + endRow; ++node) {
        function(node, node + shift);
      }
    }
  }
}

/// fills the level's offsets and couplings from addCouplings(addPenalty, addData), which calls addPenalty(node, offset, value) and
/// addData(node, offset, weight, matrix) for the terms of the matrix (value x 1, weight x matrix between the node and the node at the
/// offset); offsets beyond +-MaxNodeOffset are not represented
template <typename AddCouplings>
void assembleLevel(MultigridLevel& level, AddCouplings addCouplings)
{
  constexpr int Width = 2 * MaxNodeOffset + 1;
  auto key = [](const NodeOffset& offset) {
    if (std::abs(offset.row) > MaxNodeOffset || std::abs(offset.y2x) > MaxNodeOffset || std::abs(offset.z2x) > MaxNodeOffset) {
      return -1;
    }
    return offset.row + MaxNodeOffset + Width * (offset.y2x + MaxNodeOffset + Width * (offset.z2x + MaxNodeOffset));
  };
  std::array<int, Width * Width * Width> penaltySlot;
  std::array<int, Width * Width * Width> dataSlot;
  penaltySlot.fill(-1);
  dataSlot.fill(-1);
  auto addSlot = [&](const NodeOffset& offset, std::array<int, Width * Width * Width>& slots, std::vector<NodeOffset>& offsets) {
    const int index = key(offset);
    if (index >= 0 && slots[index] < 0) {
      slots[index] = int(offsets.size());
      offsets.push_back(offset);
    }
  };
  addCouplings([&](int, const NodeOffset& offset, double) { addSlot(offset, penaltySlot, level.penaltyOffsets); },
               [&](int, const NodeOffset& offset, double, const Matrix3&) { addSlot(offset, dataSlot, level.dataOffsets); });
  const size_t nNodes = level.nNodes();
  level.penalty.assign(level.penaltyOffsets.size() * nNodes, 0.);
  level.data.assign(level.dataOffsets.size() * nNodes, Matrix3{});
  addCouplings(
    [&](int node, const NodeOffset& offset, double value) {
      const int index = key(offset);
      if (index >= 0) {
        level.penalty[penaltySlot[index] * nNodes + node] += value;
      }
    },
    [&](int node, const NodeOffset& offset, double weight, const Matrix3& matrix) {
      const int index = key(offset);
      if (index >= 0) {
        Matrix3& coupling = level.data[dataSlot[index] * nNodes + node];
        for (int row = 0; row < 3; ++row) {
          for (int column = 0; column < 3; ++column) {
            coupling[row][column] += weight * matrix[row][column];
          }
        }
      }
    });
}

/// the coarse bins of a bin (1 or 2) and their weights in the linear interpolation: the fine bin 2 k is the coarse bin k, 2 k + 1 the mean
/// of k and k + 1; the same bin if not coarsened
int coarseBins(int bin, bool coarsened, std::array<int, 2>& coarse, std::array<double, 2>& weight)
{
  if (!coarsened) {
    coarse[0] = bin;
    weight[0] = 1.;
    return 1;
  }
  coarse = {bin / 2, bin / 2 + 1};
  weight = {bin % 2 ? 0.5 : 1., 0.5};
  return bin % 2 ? 2 : 1;
}

/// calls function(fine column, coarse column, weight) for the interpolation P of the coarse level's columns to the fine level's
template <typename Function>
void forEachInterpolation(const MultigridLevel& fine, const MultigridLevel& coarse, Function function)
{
  std::array<int, 2> binY2X;
  std::array<int, 2> binZ2X;
  std::array<double, 2> weightY2X;
  std::array<double, 2> weightZ2X;
  for (int z2x = 0; z2x < fine.nZ2X; ++z2x) {
    const int nZ2X = coarseBins(z2x, coarse.nZ2X != fine.nZ2X, binZ2X, weightZ2X);
    for (int y2x = 0; y2x < fine.nY2X; ++y2x) {
      const int nY2X = coarseBins(y2x, coarse.nY2X != fine.nY2X, binY2X, weightY2X);
      for (int j = 0; j < nY2X; ++j) {
        for (int k = 0; k < nZ2X; ++k) {
          function(y2x + fine.nY2X * z2x, binY2X[j] + coarse.nY2X * binZ2X[k], weightY2X[j] * weightZ2X[k]);
        }
      }
    }
  }
}

/// the coarse level: the rows kept, y/x and / or z/x halved; the matrix P^T A P. Semi-coarsening: only the direction(s) whose penalty is
/// >= MultigridStiffness x the stiffer one's (the penalty energy of the vector alternating along it and constant along the others, = -2 x
/// the sum of the penalty couplings at odd offsets along it): the error after the sweeps is smooth only along the stiff direction(s), e.g.
/// y/x alone first at the default binning (~20 x stiffer than z/x), whose stiffness relative to z/x drops by ~16 per halving
MultigridLevel coarseLevel(const MultigridLevel& fine)
{
  const size_t nNodes = fine.nNodes();
  double stiffnessY2X = 0.;
  double stiffnessZ2X = 0.;
  for (size_t slot = 0; slot < fine.penaltyOffsets.size(); ++slot) {
    const NodeOffset& offset = fine.penaltyOffsets[slot];
    double sum = 0.;
    forEachCoupling(fine, offset, [&](int node, int) { sum += fine.penalty[slot * nNodes + node]; });
    stiffnessY2X -= offset.y2x % 2 ? 2. * sum : 0.;
    stiffnessZ2X -= offset.z2x % 2 ? 2. * sum : 0.;
  }
  const double stiffest = std::max(stiffnessY2X, stiffnessZ2X);
  bool coarsenY2X = fine.nY2X > 2 && stiffnessY2X >= MultigridStiffness * stiffest;
  bool coarsenZ2X = fine.nZ2X > 2 && stiffnessZ2X >= MultigridStiffness * stiffest;
  if (!coarsenY2X && !coarsenZ2X) {
    coarsenY2X = fine.nY2X > 2;
    coarsenZ2X = fine.nZ2X > 2;
  }
  MultigridLevel coarse;
  coarse.nRows = fine.nRows;
  coarse.nY2X = coarsenY2X ? fine.nY2X / 2 + 1 : fine.nY2X;
  coarse.nZ2X = coarsenZ2X ? fine.nZ2X / 2 + 1 : fine.nZ2X;
  // a coupling of the fine nodes (row, column a) and (row + d, column b) adds w_a(A) w_b(B) x it to the coarse nodes (row, A) and (row + d,
  // B), for the coarse columns A, B of a, b and their interpolation weights: function(fine column, coarse column, coarse offset, weight)
  auto forEachCoarseCoupling = [&](const NodeOffset& offset, auto function) {
    std::array<int, 2> binY2X;
    std::array<int, 2> binZ2X;
    std::array<int, 2> neighbourY2X;
    std::array<int, 2> neighbourZ2X;
    std::array<double, 2> weightY2X;
    std::array<double, 2> weightZ2X;
    std::array<double, 2> neighbourWeightY2X;
    std::array<double, 2> neighbourWeightZ2X;
    for (int z2x = std::max(0, -offset.z2x); z2x < std::min(fine.nZ2X, fine.nZ2X - offset.z2x); ++z2x) {
      const int nZ2X = coarseBins(z2x, coarsenZ2X, binZ2X, weightZ2X);
      const int nNeighbourZ2X = coarseBins(z2x + offset.z2x, coarsenZ2X, neighbourZ2X, neighbourWeightZ2X);
      for (int y2x = std::max(0, -offset.y2x); y2x < std::min(fine.nY2X, fine.nY2X - offset.y2x); ++y2x) {
        const int nY2X = coarseBins(y2x, coarsenY2X, binY2X, weightY2X);
        const int nNeighbourY2X = coarseBins(y2x + offset.y2x, coarsenY2X, neighbourY2X, neighbourWeightY2X);
        for (int j = 0; j < nY2X; ++j) {
          for (int k = 0; k < nZ2X; ++k) {
            for (int l = 0; l < nNeighbourY2X; ++l) {
              for (int m = 0; m < nNeighbourZ2X; ++m) {
                const NodeOffset coarseOffset{offset.row, neighbourY2X[l] - binY2X[j], neighbourZ2X[m] - binZ2X[k]};
                const double weight = weightY2X[j] * weightZ2X[k] * neighbourWeightY2X[l] * neighbourWeightZ2X[m];
                function(y2x + fine.nY2X * z2x, binY2X[j] + coarse.nY2X * binZ2X[k], coarseOffset, weight);
              }
            }
          }
        }
      }
    }
  };
  const int nRows = fine.nRows;
  assembleLevel(coarse, [&](auto addPenalty, auto addData) {
    for (size_t slot = 0; slot < fine.penaltyOffsets.size(); ++slot) {
      const NodeOffset& offset = fine.penaltyOffsets[slot];
      forEachCoarseCoupling(offset, [&](int fineColumn, int coarseColumn, const NodeOffset& coarseOffset, double weight) {
        const double* coupling = fine.penalty.data() + slot * nNodes + nRows * fineColumn;
        for (int row = std::max(0, -offset.row); row < std::min(nRows, nRows - offset.row); ++row) {
          addPenalty(row + nRows * coarseColumn, coarseOffset, weight * coupling[row]);
        }
      });
    }
    for (size_t slot = 0; slot < fine.dataOffsets.size(); ++slot) {
      const NodeOffset& offset = fine.dataOffsets[slot];
      forEachCoarseCoupling(offset, [&](int fineColumn, int coarseColumn, const NodeOffset& coarseOffset, double weight) {
        const Matrix3* coupling = fine.data.data() + slot * nNodes + nRows * fineColumn;
        for (int row = std::max(0, -offset.row); row < std::min(nRows, nRows - offset.row); ++row) {
          addData(row + nRows * coarseColumn, coarseOffset, weight, coupling[row]);
        }
      });
    }
  });
  return coarse;
}

/// a b (3 x 3)
Matrix3 matrixProduct(const Matrix3& a, const Matrix3& b)
{
  Matrix3 product{};
  for (int row = 0; row < 3; ++row) {
    for (int column = 0; column < 3; ++column) {
      product[row][column] = a[row][0] * b[0][column] + a[row][1] * b[1][column] + a[row][2] * b[2][column];
    }
  }
  return product;
}

/// the block L D L^T of each column's block of the level's matrix: the diagonal blocks (penalty x 1 + data) and the couplings to the
/// nodes 1 and 2 rows before and after (penalty x 1; L(node, node - 2) = A(node, node - 2) D^-1(node - 2), not stored). A pivot whose
/// eigenvalues' geometric mean is below 1e-12 x the mean diagonal element of its block is left out: D^-1 = 0 (a node without data and
/// penalty, no stencil and no ridge, stays 0)
void factorColumns(MultigridLevel& level)
{
  const size_t nNodes = level.nNodes();
  int diagonalSlot = -1;
  int oneRowBeforeSlot = -1;
  for (size_t slot = 0; slot < level.penaltyOffsets.size(); ++slot) {
    const NodeOffset& offset = level.penaltyOffsets[slot];
    if (offset.y2x == 0 && offset.z2x == 0 && offset.row == 0) {
      diagonalSlot = int(slot);
    } else if (offset.y2x == 0 && offset.z2x == 0 && offset.row == -1) {
      oneRowBeforeSlot = int(slot);
    } else if (offset.y2x == 0 && offset.z2x == 0 && offset.row == -2) {
      level.twoRowsBeforeSlot = int(slot);
    }
  }
  int dataDiagonalSlot = -1;
  for (size_t slot = 0; slot < level.dataOffsets.size(); ++slot) {
    const NodeOffset& offset = level.dataOffsets[slot];
    if (offset.row == 0 && offset.y2x == 0 && offset.z2x == 0) {
      dataDiagonalSlot = int(slot);
    }
  }
  auto penalty = [&](int slot, int node) { return slot >= 0 ? level.penalty[slot * nNodes + node] : 0.; };
  level.pivotInverse.assign(nNodes, Matrix3{});
  level.lower1.assign(nNodes, Matrix3{});
  for (int column = 0; column < level.nColumns(); ++column) {
    for (int row = 0; row < level.nRows; ++row) {
      const int node = row + level.nRows * column;
      // D(node) = A(node, node) - sum over k = node - 2, node - 1 of L(node, k) D(k) L(node, k)^T, where L(node, node - 2) D(node - 2) =
      // A(node, node - 2) = coupling2 x 1 and L(node, node - 1) D(node - 1) = A(node, node - 1) - coupling2 L(node - 1, node - 2)^T =: reduced
      Matrix3 pivot = dataDiagonalSlot >= 0 ? level.data[dataDiagonalSlot * nNodes + node] : Matrix3{};
      for (int c = 0; c < 3; ++c) {
        pivot[c][c] += penalty(diagonalSlot, node);
      }
      const double scale = 1e-12 * (pivot[0][0] + pivot[1][1] + pivot[2][2]) / 3.;
      const double coupling2 = row >= 2 ? penalty(level.twoRowsBeforeSlot, node) : 0.;
      if (row >= 2) {
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            pivot[i][j] -= coupling2 * coupling2 * level.pivotInverse[node - 2][i][j];
          }
        }
      }
      if (row >= 1) {
        Matrix3 reduced{};
        for (int i = 0; i < 3; ++i) {
          reduced[i][i] = penalty(oneRowBeforeSlot, node);
          for (int j = 0; j < 3; ++j) {
            reduced[i][j] -= coupling2 * level.lower1[node - 1][j][i];
          }
        }
        level.lower1[node] = matrixProduct(reduced, level.pivotInverse[node - 1]);
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            pivot[i][j] -= level.lower1[node][i][0] * reduced[j][0] + level.lower1[node][i][1] * reduced[j][1] + level.lower1[node][i][2] * reduced[j][2];
          }
        }
      }
      level.pivotInverse[node] = determinant(pivot) > scale * scale * scale ? inverseSymmetric(pivot) : Matrix3{};
    }
  }
}

/// solution = the column's block of the matrix ^-1 rhs (the column's 3 nRows values)
void solveColumn(const MultigridLevel& level, int column, const double* rhs, double* solution)
{
  const int first = level.nRows * column;
  const double* coupling2 = level.twoRowsBeforeSlot >= 0 ? level.penalty.data() + level.twoRowsBeforeSlot * size_t(level.nNodes()) + first : nullptr;
  // L y = rhs, solution = D^-1 y; L(row, row - 2) y(row - 2) = A(row, row - 2) solution(row - 2)
  Vector3 previous{}; // y of the row before
  for (int row = 0; row < level.nRows; ++row) {
    const Vector3 term1 = multiply(level.lower1[first + row], previous);
    Vector3 y;
    for (int c = 0; c < 3; ++c) {
      y[c] = rhs[3 * row + c] - term1[c] - (row >= 2 && coupling2 ? coupling2[row] * solution[3 * (row - 2) + c] : 0.);
    }
    const Vector3 scaled = multiply(level.pivotInverse[first + row], y);
    std::copy(scaled.begin(), scaled.end(), solution + 3 * row);
    previous = y;
  }
  // L^T solution = D^-1 y; L(row + 2, row)^T = A(row + 2, row) D^-1(row)
  for (int row = level.nRows - 2; row >= 0; --row) {
    const Matrix3& lower = level.lower1[first + row + 1];
    Vector3 term2{};
    if (row + 2 < level.nRows && coupling2) {
      term2 = multiply(level.pivotInverse[first + row], {solution[3 * (row + 2)], solution[3 * (row + 2) + 1], solution[3 * (row + 2) + 2]});
    }
    for (int c = 0; c < 3; ++c) {
      const double term1 = lower[0][c] * solution[3 * (row + 1)] + lower[1][c] * solution[3 * (row + 1) + 1] + lower[2][c] * solution[3 * (row + 1) + 2];
      solution[3 * row + c] -= term1 + (row + 2 < level.nRows && coupling2 ? coupling2[row + 2] * term2[c] : 0.);
    }
  }
}

/// the couplings of a column's nodes to the other columns (laterOnly: to the columns after it in the order of the nodes) times x: y = the
/// column's 3 nRows values
void multiplyColumn(const MultigridLevel& level, int column, bool laterOnly, const double* x, double* y)
{
  const size_t nNodes = level.nNodes();
  const int y2x = column % level.nY2X;
  const int z2x = column / level.nY2X;
  const int first = level.nRows * column;
  std::fill(y, y + 3 * level.nRows, 0.);
  // the rows of the column whose neighbour at the offset is on the grid; false if none or not to be included
  auto rowRange = [&](const NodeOffset& offset, int& firstRow, int& endRow) {
    const int columnShift = offset.y2x + level.nY2X * offset.z2x;
    const bool included = laterOnly ? columnShift > 0 : columnShift != 0;
    const bool onGrid = y2x + offset.y2x >= 0 && y2x + offset.y2x < level.nY2X && z2x + offset.z2x >= 0 && z2x + offset.z2x < level.nZ2X;
    firstRow = std::max(0, -offset.row);
    endRow = std::min(level.nRows, level.nRows - offset.row);
    return included && onGrid && firstRow < endRow;
  };
  int firstRow = 0;
  int endRow = 0;
  for (size_t slot = 0; slot < level.penaltyOffsets.size(); ++slot) {
    if (!rowRange(level.penaltyOffsets[slot], firstRow, endRow)) {
      continue;
    }
    const double* coupling = level.penalty.data() + slot * nNodes + first;
    const double* neighbour = x + 3 * (first + level.shift(level.penaltyOffsets[slot]));
    for (int row = firstRow; row < endRow; ++row) {
      for (int c = 0; c < 3; ++c) {
        y[3 * row + c] += coupling[row] * neighbour[3 * row + c];
      }
    }
  }
  for (size_t slot = 0; slot < level.dataOffsets.size(); ++slot) {
    if (!rowRange(level.dataOffsets[slot], firstRow, endRow)) {
      continue;
    }
    const Matrix3* coupling = level.data.data() + slot * nNodes + first;
    const double* neighbour = x + 3 * (first + level.shift(level.dataOffsets[slot]));
    for (int row = firstRow; row < endRow; ++row) {
      const Vector3 product = multiply(coupling[row], {neighbour[3 * row], neighbour[3 * row + 1], neighbour[3 * row + 2]});
      for (int c = 0; c < 3; ++c) {
        y[3 * row + c] += product[c];
      }
    }
  }
}

/// a block Gauss-Seidel sweep over the level's columns, forward (in the order of the nodes) or backward: per column solution = its block
/// ^-1 (rhs - the couplings to the other columns at their current solution)
void gaussSeidelSweep(MultigridLevel& level, const double* rhs, double* solution, bool forward)
{
  const int nColumnValues = 3 * level.nRows;
  double* columnRhs = level.columnRhs.data();
  for (int step = 0; step < level.nColumns(); ++step) {
    const int column = forward ? step : level.nColumns() - 1 - step;
    multiplyColumn(level, column, false, solution, columnRhs);
    for (int k = 0; k < nColumnValues; ++k) {
      columnRhs[k] = rhs[nColumnValues * column + k] - columnRhs[k];
    }
    solveColumn(level, column, columnRhs, solution + nColumnValues * column);
  }
}

/// the coarsest level: the banded Cholesky factor of its matrix (unknowns ordered by row, column, component: the couplings within +-2
/// rows); a pivot below 1e-12 x its diagonal element is left out
void factorCoarsest(MultigridLevel& level)
{
  const int nColumns = level.nColumns();
  const int nUnknowns = 3 * level.nNodes();
  level.bandWidth = 3 * (MaxNodeOffset * nColumns + nColumns) - 1;
  const int width = level.bandWidth + 1;
  level.band.assign(size_t(nUnknowns) * width, 0.);
  level.work.assign(nUnknowns, 0.);
  auto unknown = [&](int node, int component) { return 3 * (node / level.nRows + nColumns * (node % level.nRows)) + component; };
  auto addLower = [&](int i, int j, double value) {
    if (j <= i) {
      level.band[size_t(i) * width + (i - j)] += value;
    }
  };
  const size_t nNodes = level.nNodes();
  for (size_t slot = 0; slot < level.penaltyOffsets.size(); ++slot) {
    const double* coupling = level.penalty.data() + slot * nNodes;
    forEachCoupling(level, level.penaltyOffsets[slot], [&](int node, int neighbour) {
      for (int c = 0; c < 3; ++c) {
        addLower(unknown(node, c), unknown(neighbour, c), coupling[node]);
      }
    });
  }
  for (size_t slot = 0; slot < level.dataOffsets.size(); ++slot) {
    const Matrix3* coupling = level.data.data() + slot * nNodes;
    forEachCoupling(level, level.dataOffsets[slot], [&](int node, int neighbour) {
      for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
          addLower(unknown(node, row), unknown(neighbour, column), coupling[node][row][column]);
        }
      }
    });
  }
  double* factor = level.band.data();
  for (int i = 0; i < nUnknowns; ++i) {
    const int first = std::max(0, i - level.bandWidth);
    for (int j = first; j <= i; ++j) {
      double sum = factor[size_t(i) * width + (i - j)];
      for (int k = std::max(first, j - level.bandWidth); k < j; ++k) {
        sum -= factor[size_t(i) * width + (i - k)] * factor[size_t(j) * width + (j - k)];
      }
      if (j < i) {
        const double pivot = factor[size_t(j) * width];
        factor[size_t(i) * width + (i - j)] = pivot > 0. ? sum / pivot : 0.;
      } else {
        const double diagonal = factor[size_t(i) * width];
        factor[size_t(i) * width] = sum > 1e-12 * diagonal ? std::sqrt(sum) : 0.;
      }
    }
  }
}

/// solution = the coarsest level's matrix ^-1 rhs
void solveCoarsest(MultigridLevel& level, const double* rhs, double* solution)
{
  const int nColumns = level.nColumns();
  const int nUnknowns = 3 * level.nNodes();
  const int width = level.bandWidth + 1;
  const double* factor = level.band.data();
  double* permuted = level.work.data();
  auto unknown = [&](int node) { return 3 * (node / level.nRows + nColumns * (node % level.nRows)); };
  for (int node = 0; node < level.nNodes(); ++node) {
    std::copy(rhs + 3 * node, rhs + 3 * node + 3, permuted + unknown(node));
  }
  for (int i = 0; i < nUnknowns; ++i) {
    double sum = permuted[i];
    for (int k = std::max(0, i - level.bandWidth); k < i; ++k) {
      sum -= factor[size_t(i) * width + (i - k)] * permuted[k];
    }
    permuted[i] = factor[size_t(i) * width] > 0. ? sum / factor[size_t(i) * width] : 0.;
  }
  for (int i = nUnknowns - 1; i >= 0; --i) {
    double sum = permuted[i];
    for (int k = i + 1; k <= std::min(nUnknowns - 1, i + level.bandWidth); ++k) {
      sum -= factor[size_t(k) * width + (k - i)] * permuted[k];
    }
    permuted[i] = factor[size_t(i) * width] > 0. ? sum / factor[size_t(i) * width] : 0.;
  }
  for (int node = 0; node < level.nNodes(); ++node) {
    std::copy(permuted + unknown(node), permuted + unknown(node) + 3, solution + 3 * node);
  }
}

/// the V-cycle from the level index on: solution = M^-1 rhs (M symmetric positive definite)
void vCycle(std::vector<MultigridLevel>& levels, size_t index, const double* rhs, double* solution)
{
  MultigridLevel& level = levels[index];
  if (index + 1 == levels.size()) {
    solveCoarsest(level, rhs, solution);
    return;
  }
  MultigridLevel& coarse = levels[index + 1];
  const int nValues = 3 * level.nNodes();
  const int nColumnValues = 3 * level.nRows;
  std::fill(solution, solution + nValues, 0.);
  gaussSeidelSweep(level, rhs, solution, true);
  // the coarse correction of the residual (P^T, solve, P): the forward sweep from 0 solved each column with the later columns at 0, its
  // residual rhs - A solution = -(the couplings to the later columns) x solution
  double* residual = level.work.data();
  for (int column = 0; column < level.nColumns(); ++column) {
    multiplyColumn(level, column, true, solution, residual + nColumnValues * column);
  }
  std::fill(coarse.rhs.begin(), coarse.rhs.end(), 0.);
  forEachInterpolation(level, coarse, [&](int fineColumn, int coarseColumn, double weight) {
    for (int k = 0; k < nColumnValues; ++k) {
      coarse.rhs[nColumnValues * coarseColumn + k] -= weight * residual[nColumnValues * fineColumn + k];
    }
  });
  vCycle(levels, index + 1, coarse.rhs.data(), coarse.solution.data());
  forEachInterpolation(level, coarse, [&](int fineColumn, int coarseColumn, double weight) {
    for (int k = 0; k < nColumnValues; ++k) {
      solution[nColumnValues * fineColumn + k] += weight * coarse.solution[nColumnValues * coarseColumn + k];
    }
  });
  gaussSeidelSweep(level, rhs, solution, false);
}

/// the multigrid hierarchy of a sector's normal matrix (see solveSector), down to <= MultigridCoarsestColumns columns (or 2 x 2 bins)
std::vector<MultigridLevel> makeMultigrid(const std::vector<Stencil>& stencils, const std::array<int, 3>& strides,
                                          const std::vector<double>& penaltyDiagonal, double ridge, const std::vector<DataVoxel>& voxels)
{
  MultigridLevel fine;
  fine.nRows = strides[1];
  fine.nY2X = strides[2] / strides[1];
  fine.nZ2X = int(penaltyDiagonal.size()) / strides[2];
  auto offsetBetween = [&](int node, int other) {
    NodeOffset offset;
    offset.row = other % fine.nRows - node % fine.nRows;
    offset.y2x = other / fine.nRows % fine.nY2X - node / fine.nRows % fine.nY2X;
    offset.z2x = other / strides[2] - node / strides[2];
    return offset;
  };
  assembleLevel(fine, [&](auto addPenalty, auto addData) {
    for (int node = 0; node < fine.nNodes(); ++node) {
      addPenalty(node, NodeOffset{}, penaltyDiagonal[node] + ridge);
    }
    for (const auto& stencil : stencils) {
      const int step = strides[stencil.axis];
      for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
          if (l != k) {
            NodeOffset offset;
            (stencil.axis == 0 ? offset.row : (stencil.axis == 1 ? offset.y2x : offset.z2x)) = l - k;
            addPenalty(stencil.first + k * step, offset, stencil.lambda * stencil.coefficient[k] * stencil.coefficient[l]);
          }
        }
      }
    }
    for (const auto& voxel : voxels) {
      for (int j = 0; j < voxel.nNodes; ++j) {
        for (int k = 0; k < voxel.nNodes; ++k) {
          addData(voxel.node[j], offsetBetween(voxel.node[j], voxel.node[k]), voxel.weight * voxel.nodeWeight[j] * voxel.nodeWeight[k], voxel.a);
        }
      }
    }
  });
  std::vector<MultigridLevel> levels;
  levels.push_back(std::move(fine));
  while (levels.back().nColumns() > MultigridCoarsestColumns && (levels.back().nY2X > 2 || levels.back().nZ2X > 2)) {
    factorColumns(levels.back());
    levels.back().work.assign(3 * size_t(levels.back().nNodes()), 0.);
    levels.back().columnRhs.assign(3 * size_t(levels.back().nRows), 0.);
    MultigridLevel coarse = coarseLevel(levels.back());
    coarse.rhs.assign(3 * size_t(coarse.nNodes()), 0.);
    coarse.solution.assign(3 * size_t(coarse.nNodes()), 0.);
    levels.push_back(std::move(coarse));
  }
  factorCoarsest(levels.back());
  return levels;
}

/// the penalised fit of one sector: minimises the voxels' chi2 (at their nodes, robust weights) + the penalty (stencils, strides of the
/// sector-local node index along rows / y/x / z/x, penaltyDiagonal = its diagonal per node) + ridge |x|^2 by conjugate gradients
/// preconditioned by a multigrid V-cycle (see above) from the given nodes (3 per node: dX, dY, dZ; the solution on return); robust > 0:
/// Tukey reweighting cycles
void solveSector(const std::vector<Stencil>& stencils, const std::array<int, 3>& strides, const std::vector<double>& penaltyDiagonal,
                 double ridge, int robust, std::vector<DataVoxel>& voxels, double* sectorNodes, SectorStatistics& statistics)
{
  const int voxelsPerSector = int(penaltyDiagonal.size());
  statistics = SectorStatistics{};
  if (voxels.empty()) {
    return;
  }
  const int n = 3 * voxelsPerSector;
  auto nodeValue = [&](const DataVoxel& voxel, const double* x) {
    Vector3 theta{};
    for (int j = 0; j < voxel.nNodes; ++j) {
      for (int c = 0; c < 3; ++c) {
        theta[c] += voxel.nodeWeight[j] * x[3 * voxel.node[j] + c];
      }
    }
    return theta;
  };
  for (int cycle = 0; cycle <= robust; ++cycle) {
    if (cycle > 0) {
      // Tukey's bisquare of sqrt(excess chi2 / rank): the voxel's chi2 at the solution minus its minimum -b^T a^+ b
      statistics.nLowWeight = 0;
      statistics.nZeroWeight = 0;
      for (auto& voxel : voxels) {
        const Vector3 theta = nodeValue(voxel, sectorNodes);
        int rank = 0;
        const double minimum = quadraticMinimum(voxel.a, voxel.b, rank);
        const double excess = dot(theta, multiply(voxel.a, theta)) - 2. * dot(voxel.b, theta) - minimum;
        const double normalised = rank > 0 ? std::sqrt(std::max(0., excess) / rank) / Tukey : 0.;
        voxel.weight = normalised < 1. ? (1. - normalised * normalised) * (1. - normalised * normalised) : 0.;
        statistics.nLowWeight += voxel.weight < 0.5;
        statistics.nZeroWeight += voxel.weight == 0.;
      }
    }
    // the normal equations: ridge + penalty + sum over the voxels of w_j w_k a at their nodes; preconditioner = a multigrid V-cycle
    std::vector<double> rhs(n, 0.);
    std::vector<Matrix3> weighted(voxels.size());
    for (size_t i = 0; i < voxels.size(); ++i) {
      const auto& voxel = voxels[i];
      Vector3 weightedB;
      for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
          weighted[i][row][column] = voxel.weight * voxel.a[row][column];
        }
        weightedB[row] = voxel.weight * voxel.b[row];
      }
      for (int j = 0; j < voxel.nNodes; ++j) {
        for (int row = 0; row < 3; ++row) {
          rhs[3 * voxel.node[j] + row] += voxel.nodeWeight[j] * weightedB[row];
        }
      }
    }
    std::vector<MultigridLevel> levels = makeMultigrid(stencils, strides, penaltyDiagonal, ridge, voxels);
    auto apply = [&](const std::vector<double>& x, std::vector<double>& y) {
      for (int i = 0; i < n; ++i) {
        y[i] = ridge * x[i];
      }
      for (const auto& stencil : stencils) {
        const int step = strides[stencil.axis];
        for (int c = 0; c < 3; ++c) {
          double difference = 0.;
          for (int k = 0; k < 3; ++k) {
            difference += stencil.coefficient[k] * x[3 * (stencil.first + k * step) + c];
          }
          for (int k = 0; k < 3; ++k) {
            y[3 * (stencil.first + k * step) + c] += stencil.lambda * stencil.coefficient[k] * difference;
          }
        }
      }
      for (size_t i = 0; i < voxels.size(); ++i) {
        const auto& voxel = voxels[i];
        const Vector3 u = multiply(weighted[i], nodeValue(voxel, x.data()));
        for (int j = 0; j < voxel.nNodes; ++j) {
          for (int c = 0; c < 3; ++c) {
            y[3 * voxel.node[j] + c] += voxel.nodeWeight[j] * u[c];
          }
        }
      }
    };
    auto precondition = [&](const std::vector<double>& r, std::vector<double>& z) { vCycle(levels, 0, r.data(), z.data()); };
    // preconditioned conjugate gradients from the previous cycle's solution
    std::vector<double> x(sectorNodes, sectorNodes + n);
    std::vector<double> r(n);
    std::vector<double> z(n);
    std::vector<double> q(n);
    const double rhsNorm = std::sqrt(dot(rhs, rhs));
    statistics.iterations = 0;
    statistics.residual = 0.;
    if (rhsNorm == 0.) {
      std::fill(sectorNodes, sectorNodes + n, 0.);
      continue;
    }
    apply(x, q);
    for (int i = 0; i < n; ++i) {
      r[i] = rhs[i] - q[i];
    }
    precondition(r, z);
    std::vector<double> p = z;
    double rz = dot(r, z);
    statistics.residual = std::sqrt(dot(r, r)) / rhsNorm;
    while (statistics.residual > CGTolerance && statistics.iterations < CGMaxIterations) {
      apply(p, q);
      const double alpha = rz / dot(p, q);
      for (int i = 0; i < n; ++i) {
        x[i] += alpha * p[i];
        r[i] -= alpha * q[i];
      }
      statistics.residual = std::sqrt(dot(r, r)) / rhsNorm;
      ++statistics.iterations;
      precondition(r, z);
      const double rzNew = dot(r, z);
      const double beta = rzNew / rz;
      for (int i = 0; i < n; ++i) {
        p[i] = z[i] + beta * p[i];
      }
      rz = rzNew;
    }
    statistics.finite = statistics.finite && std::all_of(x.begin(), x.end(), [](double value) { return std::isfinite(value); });
    std::copy(x.begin(), x.end(), sectorNodes);
  }
}
} // namespace voxel_smoother

#ifndef SMOOTHVOXELMAP_SOLVER_ONLY
namespace voxel_smoother
{
using TrackResiduals = o2::tpc::TrackResiduals;

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
  std::vector<long> sortedTimes(times ? *times : std::vector<long>());
  orbitLumiTree->ResetBranchAddresses();
  delete times;
  if (sortedTimes.empty()) {
    LOGP(warning, "IDC join: no TF times in OrbitLumiInfo");
    return false;
  }
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
  std::sort(idcs.begin(), idcs.end()); // the median as TMath::Median (SmoothingExtrapolate.C): the mean of the central two for an even number
  const size_t middle = idcs.size() / 2;
  medianIDC = idcs.size() % 2 ? idcs[middle] : 0.5f * (idcs[middle - 1] + idcs[middle]);
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
  const std::array<int, 3> strides = {1, nBinsX, nBinsX * nBinsY2X}; // getGlbVoxBin = row + nX * (y/x bin + nY2X * z/x bin)
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
  std::array<SectorStatistics, NSectors> statistics{};
  auto smoothSector = [&](int sector) {
    solveSector(stencils, strides, penaltyDiagonal, ridge, robust, data[sector], nodes.data() + 3 * size_t(sector) * voxelsPerSector, statistics[sector]);
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

  const int nNotConverged = std::count_if(statistics.begin(), statistics.end(), [](const auto& s) { return s.residual > CGTolerance; });
  const int nNotFinite = std::count_if(statistics.begin(), statistics.end(), [](const auto& s) { return !s.finite; });
  const auto [fewest, most] = std::minmax_element(statistics.begin(), statistics.end(), [](const auto& a, const auto& b) { return a.iterations < b.iterations; });
  const auto worst = std::max_element(statistics.begin(), statistics.end(), [](const auto& a, const auto& b) { return a.residual < b.residual; });
  LOGP(info, "Smoothing: {} stencils per sector; CG iterations per sector {} - {}, worst relative residual {:.2e}", stencils.size(),
       fewest->iterations, most->iterations, worst->residual);
  if (robust > 0) {
    LOGP(info, "Smoothing: {} Tukey reweightings (c = {}), at the last {} voxels below weight 0.5, {} at 0", robust, Tukey,
         std::accumulate(statistics.begin(), statistics.end(), 0L, [](long sum, const auto& s) { return sum + s.nLowWeight; }),
         std::accumulate(statistics.begin(), statistics.end(), 0L, [](long sum, const auto& s) { return sum + s.nZeroWeight; }));
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
                    float tauRow = 1000.,      // prior rms of the curvature (um / cm^2) along the pad rows (0: no penalty); !curvaturePrior: second difference (um)
                    float tauY2X = 2000.,      // ... along y/x
                    float tauZ2X = 1000.,      // ... along z/x
                    float tau0 = 0.,           // prior rms (um) of a node: ridge towards 0 (0: none)
                    float trustFraction = -1., // the innermost IROC rows below this fraction of the column's plateau entries are no data (0: off; -1: auto)
                    float minEntries = 20.,    // voxResTree: voxels with fewer entries are no data; fineSums: the voxel solve's minimum
                    int errorMode = -1,        // voxResTree: -1 auto, 0 the input's E and EXYCorr, 1 s / sqrt(entries), 2 from fineSums (see the file description)
                    int joinIDCs = -1,         // meanIDC / medianIDC from the CCDB scalers at the input's TFs: -1 if the input has none, 1 always, 0 never
                    int robust = 0,            // Tukey reweighting cycles
                    int nThreads = 8,
                    bool curvaturePrior = true,  // the taus as the prior rms of the curvature over 1 cm^3 (um / cm^2): the same for any binning (makeStencils)
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
    source = voxResTree ? "voxResTree" : "fineSums";
  }
  const bool useSums = source == "fineSums";
  if ((useSums && !sumsTree) || (!useSums && (source != "voxResTree" || !voxResTree))) {
    LOGP(error, "{}: no input {} (source: voxResTree, fineSums or auto)", inputFile, source);
    return;
  }
  if (errorMode == 2 && (useSums || !sumsTree)) {
    LOGP(error, "{}: errorMode 2 needs a voxResTree input and the fineSums in the same file", inputFile);
    return;
  }
  if (useSums && errorMode >= 0) {
    LOGP(warning, "{}: errorMode {} ignored: the fineSums input has the errors of its voxel solve", inputFile, errorMode);
  }
  TTree* infoTree = useSums ? sumsTree : voxResTree;
  auto info = [&](const char* key) -> std::string {
    auto* object = infoTree->GetUserInfo()->FindObject(key);
    return object ? object->GetTitle() : "";
  };
  // an empty binning string is the default binning (staticMapCreatorCPM.C writes the arguments as given)
  if (!infoTree->GetUserInfo()->FindObject("y2xBinning") || !infoTree->GetUserInfo()->FindObject("z2xBinning") || info("maxZ2X").empty()) {
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
  // the fineSums per voxel, both charges added
  std::vector<ChargeSums> sums;
  auto readSums = [&]() {
    sums.assign(nVoxels, ChargeSums());
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
        return false;
      }
      auto& total = sums[voxelIndex];
      for (int charge = 0; charge < 2; ++charge) {
        total.entries += entries[charge];
        total.x += position[charge][0];
        total.y2x += position[charge][1];
        total.z2x += position[charge][2];
        for (int sum = 0; sum < ChargeSums::NSums; ++sum) {
          total.y[sum] += channelY[charge][sum];
          total.z[sum] += channelZ[charge][sum];
        }
      }
    }
    sumsTree->ResetBranchAddresses();
    return true;
  };
  if (useSums) {
    if (!readSums()) {
      return;
    }
    for (size_t voxel = 0; voxel < nVoxels; ++voxel) {
      const auto& total = sums[voxel];
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
    std::vector<bool> read(nVoxels, false);
    long nDuplicates = 0;
    for (long entry = 0; entry < voxResTree->GetEntries(); ++entry) {
      voxResTree->GetEntry(entry);
      if (result->bsec >= NSectors || result->bvox[TrackResiduals::VoxX] >= binning->getNXBins() || result->bvox[TrackResiduals::VoxF] >= binning->getNY2XBins() ||
          result->bvox[TrackResiduals::VoxZ] >= binning->getNZ2XBins()) {
        LOGP(error, "{}: voxel (sector {}, bins {} / {} / {}) outside the binning of its user info", inputFile, int(result->bsec), int(result->bvox[TrackResiduals::VoxX]),
             int(result->bvox[TrackResiduals::VoxF]), int(result->bvox[TrackResiduals::VoxZ]));
        return;
      }
      const size_t voxel = size_t(result->bsec) * voxelsPerSector + binning->getGlbVoxBin(result->bvox);
      nDuplicates += read[voxel];
      read[voxel] = true;
      results[voxel] = *result;
    }
    voxResTree->ResetBranchAddresses();
    delete result;
    if (nDuplicates) {
      LOGP(warning, "{}: {} voxels appear more than once in the voxResTree (merged maps?): the last entry is used", inputFile, nDuplicates);
    }
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
    if (errorMode == 2) {
      if (!readSums()) {
        return;
      }
      long nSolved = 0;
      for (size_t voxel = 0; voxel < nVoxels; ++voxel) {
        auto& result = results[voxel];
        VoxelSolution solution;
        const bool solved = solveVoxel(sums[voxel], minEntries, solution);
        nSolved += solved;
        for (int coordinate = 0; coordinate < 3; ++coordinate) {
          result.E[coordinate] = solved ? solution.error[coordinate] : 0.;
        }
        result.EXYCorr = solved ? solution.correlationXY : 0.;
      }
      LOGP(info, "Errors: the residual slope method's voxel solve of the fineSums ({} voxels solved, the others no data)", nSolved);
    } else if (errorMode == 1) {
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

  // the data per sector (the trust cut only with the stand-in errors, see the file description)
  if (trustFraction < 0.f) {
    trustFraction = !useSums && errorMode == 1 ? 0.7f : 0.f;
  }
  long nUntrusted = 0;
  const auto untrusted = untrustedRows(*binning, results, trustFraction, nUntrusted);
  double medianScale = 1.;
  if (useSums) {
    std::vector<double> scales;
    for (size_t voxel = 0; voxel < nVoxels; ++voxel) {
      VoxelSolution solution;
      if (solveVoxel(sums[voxel], minEntries, solution)) {
        scales.push_back(solution.chi2Scale);
      }
    }
    if (!scales.empty()) {
      std::nth_element(scales.begin(), scales.begin() + scales.size() / 2, scales.end());
      medianScale = scales[scales.size() / 2];
    }
  }
  std::vector<std::vector<DataVoxel>> data(NSectors);
  std::vector<bool> isData(nVoxels, false);
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
      const auto& total = sums[voxel];
      if (total.entries <= 0. || total.y[S::W] <= 0. || total.z[S::W] <= 0.) {
        continue;
      }
      VoxelSolution solution;
      const double scale = solveVoxel(total, minEntries, solution) ? solution.chi2Scale : medianScale;
      // the y channel (J = (-s_y, 1, 0)) and the z channel's mean W_z (<r_z> - dZ + <s_z> dX)^2
      const double meanSlopeZ = total.z[S::WS] / total.z[S::W];
      dataVoxel.a[0] = {total.y[S::WSS] + total.z[S::WS] * meanSlopeZ, -total.y[S::WS], -total.z[S::WS]};
      dataVoxel.a[1] = {-total.y[S::WS], total.y[S::W], 0.};
      dataVoxel.a[2] = {-total.z[S::WS], 0., total.z[S::W]};
      dataVoxel.b = {-total.y[S::WSR] - meanSlopeZ * total.z[S::WR], total.y[S::WR], total.z[S::WR]};
      for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
          dataVoxel.a[row][column] /= scale;
        }
        dataVoxel.b[row] /= scale;
      }
      z2xSide = (sector < NSectors / 2 ? 1. : -1.) * result.stat[TrackResiduals::VoxZ];
    } else {
      const bool fitted = (result.flags & TrackResiduals::DistDone) && !(result.flags & TrackResiduals::Masked) && result.stat[TrackResiduals::VoxV] >= std::max(1.f, minEntries);
      if (!fitted || !(result.E[0] > 0.f) || !(result.E[1] > 0.f) || !(result.E[2] > 0.f) ||
          !std::isfinite(result.D[0]) || !std::isfinite(result.D[1]) || !std::isfinite(result.D[2])) {
        continue;
      }
      Matrix3 covariance{};
      covariance[0][0] = result.E[0] * result.E[0];
      covariance[1][1] = result.E[1] * result.E[1];
      covariance[0][1] = std::clamp(double(result.EXYCorr), -0.99, 0.99) * result.E[0] * result.E[1];
      covariance[1][0] = covariance[0][1];
      covariance[2][2] = result.E[2] * result.E[2];
      dataVoxel.a = inverseSymmetric(covariance);
      dataVoxel.b = multiply(dataVoxel.a, {result.D[0], result.D[1], result.D[2]});
      z2xSide = (sector < NSectors / 2 ? 1. : -1.) * result.stat[TrackResiduals::VoxZ];
    }
    int nodeIndices[4];
    double nodeWeights[4];
    dataVoxel.nNodes = bilinearNodes(*binning, row, result.stat[TrackResiduals::VoxF], z2xSide, nodeIndices, nodeWeights);
    for (int j = 0; j < dataVoxel.nNodes; ++j) {
      dataVoxel.node[j] = nodeIndices[j];
      dataVoxel.nodeWeight[j] = nodeWeights[j];
    }
    data[sector].push_back(dataVoxel);
    isData[voxel] = true;
    ++nData;
  }
  const std::string description = fmt::format(
    "penalised node fit (smoothVoxelMap.C): input {}, tau rows / y/x / z/x {} / {} / {} {}, tau0 {} um, trustFraction {} "
    "({} untrusted inner IROC voxels), minEntries {}, errors {}, robust {}",
    source, tauRow, tauY2X, tauZ2X, curvaturePrior ? "um / cm^2 (curvature over 1 cm^3)" : "um (second difference per bin)", tau0, trustFraction, nUntrusted, minEntries,
    useSums ? std::string("voxel solve") : (errorMode == 1 ? std::string("s / sqrt(entries)") : (errorMode == 2 ? std::string("from the fineSums voxel solve") : std::string("input"))), robust);
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

  // the output: the voxels with DS = the nodes at the voxel positions, the input's user info and other trees
  TFile output(outputFile.data(), "recreate");
  TTree tree("voxResTree", "voxRes map results and statistics, smoothed by smoothVoxelMap.C");
  TrackResiduals::VoxRes result;
  TrackResiduals::VoxRes* resultPointer = &result;
  tree.Branch("voxRes", &resultPointer);
  for (auto* object : *treeWithInfo->GetUserInfo()) {
    const std::string name = object->GetName();
    if (name != "smoothing" && (!joinedIDC || (name != "meanIDC" && name != "medianIDC"))) {
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
          // DS where the transform puts it, at stat: the nodes are the values at the bin centres. The voxels that are no data (empty, below
          // minEntries, untrusted) at their bin centre, DS = their node: a few clusters' mean position in a data-free region is arbitrary
          if (!isData[voxel]) {
            binning->getVoxelCoordinates(sector, binX, binY2X, binZ2X, result.stat[TrackResiduals::VoxX], result.stat[TrackResiduals::VoxF], result.stat[TrackResiduals::VoxZ]);
          }
          int nodeIndices[4];
          double nodeWeights[4];
          const double z2xSide = (sector < NSectors / 2 ? 1. : -1.) * result.stat[TrackResiduals::VoxZ];
          const int nNodes = bilinearNodes(*binning, binX, result.stat[TrackResiduals::VoxF], z2xSide, nodeIndices, nodeWeights);
          for (int coordinate = 0; coordinate < 3; ++coordinate) {
            double value = 0.;
            for (int j = 0; j < nNodes; ++j) {
              value += nodeWeights[j] * nodes[3 * (size_t(sector) * voxelsPerSector + nodeIndices[j]) + coordinate];
            }
            result.DS[coordinate] = value;
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
  LOGP(info, "{} written: {} voxels, DS = the smoothed map at the voxels' positions", outputFile, nVoxels);
}
#endif // SMOOTHVOXELMAP_SOLVER_ONLY
