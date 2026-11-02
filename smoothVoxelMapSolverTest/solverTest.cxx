// solverTest.cxx — stand-alone test of the per-sector solver of Detectors/TPC/calibration/SpacePoints/macro/smoothVoxelMap.C
// (TEST ONLY, branch tpc-fine-map-gbl-solver). Includes the macro itself with SMOOTHVOXELMAP_SOLVER_ONLY (no ROOT / O2), reads sector
// systems dumped by the macro (SMOOTHVOXELMAP_DUMP), solves each from zero with the macro's solveSector and compares with the macro's own
// solution in the file:
//   iterations, final relative residual, wall time; per component (dX, dY, dZ) rms and max |x - x_reference| over the nodes in um;
//   the objective F(x) = sum_voxels w (theta^T a theta - 2 b^T theta) + sum_stencils lambda (sum_k c_k x_k)^2 + ridge |x|^2 (theta = the
//   bilinear node value at the voxel), at both solutions: F(new) - F(reference) <= ~0 means at least as converged.
//
//   c++ -std=c++20 -O2 -o solverTest solverTest.cxx && ./solverTest data/t1_sector05.bin data/t30_sector05.bin ...

#define SMOOTHVOXELMAP_SOLVER_ONLY
#include "../Detectors/TPC/calibration/SpacePoints/macro/smoothVoxelMap.C"

#include <chrono>
#include <cstdio>

using namespace voxel_smoother;

double objective(const SectorSystem& system, const std::vector<double>& x)
{
  double sum = 0.;
  for (const auto& voxel : system.voxels) {
    Vector3 theta{};
    for (int j = 0; j < voxel.nNodes; ++j) {
      for (int c = 0; c < 3; ++c) {
        theta[c] += voxel.nodeWeight[j] * x[3 * voxel.node[j] + c];
      }
    }
    sum += voxel.weight * (dot(theta, multiply(voxel.a, theta)) - 2. * dot(voxel.b, theta));
  }
  for (const auto& stencil : system.stencils) {
    const int step = system.strides[stencil.axis];
    for (int c = 0; c < 3; ++c) {
      double difference = 0.;
      for (int k = 0; k < 3; ++k) {
        difference += stencil.coefficient[k] * x[3 * (stencil.first + k * step) + c];
      }
      sum += stencil.lambda * difference * difference;
    }
  }
  for (const double value : x) {
    sum += system.ridge * value * value;
  }
  return sum;
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    printf("usage: %s <sector system .bin> ...\n", argv[0]);
    return 1;
  }
  int nFailed = 0;
  for (int i = 1; i < argc; ++i) {
    SectorSystem system;
    if (!readSectorSystem(argv[i], system)) {
      printf("%s: cannot read\n", argv[i]);
      ++nFailed;
      continue;
    }
    const auto penaltyDiagonal = penaltyDiagonalOf(system.stencils, system.strides, system.voxelsPerSector);
    std::vector<DataVoxel> voxels = system.voxels;
    std::vector<double> x(3 * size_t(system.voxelsPerSector), 0.);
    SectorStatistics statistics;
    const auto start = std::chrono::steady_clock::now();
    solveSector(system.stencils, system.strides, penaltyDiagonal, system.ridge, system.robust, voxels, x.data(), statistics);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::array<double, 3> sumSquares{};
    std::array<double, 3> maximum{};
    for (size_t k = 0; k < x.size(); ++k) {
      const double difference = std::abs(x[k] - system.nodes[k]) / Micron;
      sumSquares[k % 3] += difference * difference;
      maximum[k % 3] = std::max(maximum[k % 3], difference);
    }
    const double nNodes = system.voxelsPerSector;
    const double objectiveNew = objective(system, x);
    const double objectiveReference = objective(system, system.nodes);
    printf("%s: %zu data voxels, %zu stencils | reference: %d iterations, residual %.2e | now: %d iterations, residual %.2e, %.2f s, finite %d\n",
           argv[i], system.voxels.size(), system.stencils.size(), system.iterations, system.residual, statistics.iterations, statistics.residual,
           seconds, int(statistics.finite));
    printf("   |x - reference| [um] rms / max: dX %.3f / %.3f, dY %.3f / %.3f, dZ %.3f / %.3f; F(now) - F(reference) = %.3e (relative %.2e)\n",
           std::sqrt(sumSquares[0] / nNodes), maximum[0], std::sqrt(sumSquares[1] / nNodes), maximum[1], std::sqrt(sumSquares[2] / nNodes), maximum[2],
           objectiveNew - objectiveReference, (objectiveNew - objectiveReference) / std::abs(objectiveReference));
    nFailed += !statistics.finite;
  }
  return nFailed ? 1 : 0;
}
