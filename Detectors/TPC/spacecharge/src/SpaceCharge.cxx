// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  SpaceCharge.cxx
/// \brief Definition of SpaceCharge class
///
/// \author  Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#include "TPCSpacecharge/SpaceCharge.h"
#include "fmt/core.h"
#include "Framework/Logger.h"
#include <chrono>

templateClassImp(o2::tpc::SpaceCharge);

using namespace o2::tpc;

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::calculateDistortionsCorrections(const o2::tpc::Side side, const GlobalDistType globalDistType, const GlobalDistCorrMethod globalDistCorrMethod)
{
  using timer = std::chrono::high_resolution_clock;
  using SC = o2::tpc::SpaceCharge<DataT, Nz, Nr, Nphi>;
  if (!mIsChargeSet[side]) {
    LOGP(ERROR, "the charge is not set!");
  }

  const std::array sglobalType{"local distortion/correction interpolator", "Electric fields"};
  const std::array sglobalDistType{"Standard method", "interpolation of global corrections"};

  LOGP(info, "====== starting calculation of distortions and corrections ======");
  LOGP(info, "Using {} threads", getNThreads());

  if (globalDistCorrMethod == SC::GlobalDistCorrMethod::LocalDistCorr) {
    LOGP(info, "calculation of global distortions and corrections are performed by using: {}", sglobalType[0]);
  } else {
    LOGP(info, "calculation of global distortions and corrections are performed by using: {}", sglobalType[1]);
  }

  if (globalDistType == SC::GlobalDistType::Standard) {
    LOGP(info, "calculation of global distortions performed by following method: {}", sglobalDistType[0]);
  } else {
    LOGP(info, "calculation of global distortions performed by following method: {}", sglobalDistType[1]);
  }

  auto startTotal = timer::now();

  auto start = timer::now();
  poissonSolver(side);
  auto stop = timer::now();
  std::chrono::duration<float> time = stop - start;
  LOGP(info, "Poisson Solver time: {}", time.count());

  start = timer::now();
  calcEField(side);
  stop = timer::now();
  time = stop - start;
  LOGP(info, "electric field calculation time: {}", time.count());

  const auto numEFields = getElectricFieldsInterpolator(side);
  if (globalDistType == SC::GlobalDistType::Standard) {
    start = timer::now();
    const auto dist = o2::tpc::SpaceCharge<DataT, Nz, Nr, Nphi>::Type::Distortions;
    calcLocalDistortionsCorrections(dist, numEFields); // local distortion calculation
    stop = timer::now();
    time = stop - start;
    LOGP(info, "local distortions time: {}", time.count());
  } else {
    LOGP(info, "skipping local distortions (not needed)");
  }

  start = timer::now();
  const auto corr = o2::tpc::SpaceCharge<DataT, Nz, Nr, Nphi>::Type::Corrections;
  calcLocalDistortionsCorrections(corr, numEFields); // local correction calculation
  stop = timer::now();
  time = stop - start;
  LOGP(info, "local corrections time: {}", time.count());

  start = timer::now();
  const auto lCorrInterpolator = getLocalCorrInterpolator(side);
  (globalDistCorrMethod == SC::GlobalDistCorrMethod::LocalDistCorr) ? calcGlobalCorrections(lCorrInterpolator) : calcGlobalCorrections(numEFields);
  stop = timer::now();
  time = stop - start;
  LOGP(info, "global corrections time: {}", time.count());

  start = timer::now();
  if (globalDistType == SC::GlobalDistType::Standard) {
    const auto lDistInterpolator = getLocalDistInterpolator(side);
    (globalDistCorrMethod == SC::GlobalDistCorrMethod::LocalDistCorr) ? calcGlobalDistortions(lDistInterpolator) : calcGlobalDistortions(numEFields);
  } else {
    const auto globalCorrInterpolator = getGlobalCorrInterpolator(side);
    calcGlobalDistWithGlobalCorrIterative(globalCorrInterpolator);
  }
  stop = timer::now();
  time = stop - start;
  LOGP(info, "global distortions time: {}", time.count());

  stop = timer::now();
  time = stop - startTotal;
  LOGP(info, "everything is done. Total Time: {}", time.count());
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
DataT SpaceCharge<DataT, Nz, Nr, Nphi>::regulateR(const DataT posR) const
{
  const DataT minR = getRMin() - 4 * getGridSpacingR();
  if (posR < minR) {
    return minR;
  }
  const DataT maxR = getRMax() + 2 * getGridSpacingR();
  if (posR > maxR) {
    return maxR;
  }
  return posR;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::setFromFile(TFile& file, const Side side)
{
  setDensityFromFile(file, side);
  setPotentialFromFile(file, side);
  setElectricFieldsFromFile(file, side);
  setLocalDistortionsFromFile(file, side);
  setLocalCorrectionsFromFile(file, side);
  setGlobalDistortionsFromFile(file, side);
  setGlobalCorrectionsFromFile(file, side);
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::setChargeDensityFromFormula(const AnalyticalFields<DataT>& formulaStruct)
{
  const Side side = formulaStruct.getSide();
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi = getPhiVertex(iPhi);
    for (size_t iR = 0; iR < Nr; ++iR) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz; ++iZ) {
        const DataT z = getZVertex(iZ);
        mDensity[side](iZ, iR, iPhi) = formulaStruct.evalDensity(z, radius, phi);
      }
    }
  }
  mIsChargeSet[side] = true;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::setPotentialFromFormula(const AnalyticalFields<DataT>& formulaStruct)
{
  const Side side = formulaStruct.getSide();
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi = getPhiVertex(iPhi);
    for (size_t iR = 0; iR < Nr; ++iR) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz; ++iZ) {
        const DataT z = getZVertex(iZ);
        mPotential[side](iZ, iR, iPhi) = formulaStruct.evalPotential(z, radius, phi);
      }
    }
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::setPotentialBoundaryFromFormula(const AnalyticalFields<DataT>& formulaStruct)
{
  const Side side = formulaStruct.getSide();
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi = getPhiVertex(iPhi);
    for (size_t iZ = 0; iZ < Nz; ++iZ) {
      const DataT z = getZVertex(iZ);
      const size_t iR = 0;
      const DataT radius = getRVertex(iR);
      mPotential[side](iZ, iR, iPhi) = formulaStruct.evalPotential(z, radius, phi);
    }
  }

  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi = getPhiVertex(iPhi);
    for (size_t iZ = 0; iZ < Nz; ++iZ) {
      const DataT z = getZVertex(iZ);
      const size_t iR = Nr - 1;
      const DataT radius = getRVertex(iR);
      mPotential[side](iZ, iR, iPhi) = formulaStruct.evalPotential(z, radius, phi);
    }
  }

  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi = getPhiVertex(iPhi);
    for (size_t iR = 0; iR < Nr; ++iR) {
      const DataT radius = getRVertex(iR);
      const size_t iZ = 0;
      const DataT z = getZVertex(iZ);
      mPotential[side](iZ, iR, iPhi) = formulaStruct.evalPotential(z, radius, phi);
    }
  }

  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi = getPhiVertex(iPhi);
    for (size_t iR = 0; iR < Nr; ++iR) {
      const DataT radius = getRVertex(iR);
      const size_t iZ = Nz - 1;
      const DataT z = getZVertex(iZ);
      mPotential[side](iZ, iR, iPhi) = formulaStruct.evalPotential(z, radius, phi);
    }
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::poissonSolver(const Side side, const int maxIteration, const DataT stoppingConvergence, const int symmetry)
{
  ASolv::setConvergenceError(stoppingConvergence);
  ASolv poissonSolver(mGrid3D);
  poissonSolver.poissonSolver3D(mPotential[side], mDensity[side], symmetry);
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::setEFieldFromFormula(const AnalyticalFields<DataT>& formulaStruct)
{
  const Side side = formulaStruct.getSide();
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    for (size_t iR = 0; iR < Nr; ++iR) {
      for (size_t iZ = 0; iZ < Nz; ++iZ) {
        const DataT radius = getRVertex(iR);
        const DataT z = getZVertex(iZ);
        const DataT phi = getPhiVertex(iPhi);
        mElectricFieldEr[side](iZ, iR, iPhi) = formulaStruct.evalEr(z, radius, phi);
        mElectricFieldEz[side](iZ, iR, iPhi) = formulaStruct.evalEz(z, radius, phi);
        mElectricFieldEphi[side](iZ, iR, iPhi) = formulaStruct.evalEphi(z, radius, phi);
      }
    }
  }
  mIsEfieldSet[side] = true;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::calcEField(const Side side)
{
#pragma omp parallel for num_threads(sNThreads)
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const int symmetry = 0;
    size_t tmpPlus = iPhi + 1;
    int signPlus = 1;
    int tmpMinus = static_cast<int>(iPhi - 1);
    int signMinus = 1;
    if (symmetry == 1 || symmetry == -1) { // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
      if (tmpPlus > Nphi - 1) {
        if (symmetry == -1) {
          signPlus = -1;
        }
        tmpPlus = Nphi - 2;
      }
      if (tmpMinus < 0) {
        tmpMinus = 1; // SHOULD IT BE =0?
        if (symmetry == -1) {
          signMinus = -1;
        }
      }
    } else { // No Symmetries in phi, no boundaries, the calculations is continuous across all phi
      if (tmpPlus > Nphi - 1) {
        tmpPlus = iPhi + 1 - Nphi;
      }
      if (tmpMinus < 0) {
        tmpMinus = static_cast<int>(iPhi - 1 + Nphi);
      }
    }

    // for non-boundary V
    for (size_t iR = 1; iR < Nr - 1; iR++) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 1; iZ < Nz - 1; iZ++) {
        mElectricFieldEr[side](iZ, iR, iPhi) = -1 * (mPotential[side](iZ, iR + 1, iPhi) - mPotential[side](iZ, iR - 1, iPhi)) * static_cast<DataT>(0.5) * getInvSpacingR();                                    // r direction
        mElectricFieldEz[side](iZ, iR, iPhi) = -1 * (mPotential[side](iZ + 1, iR, iPhi) - mPotential[side](iZ - 1, iR, iPhi)) * static_cast<DataT>(0.5) * getInvSpacingZ();                                    // z direction
        mElectricFieldEphi[side](iZ, iR, iPhi) = -1 * (signPlus * mPotential[side](iZ, iR, tmpPlus) - signMinus * mPotential[side](iZ, iR, tmpMinus)) * static_cast<DataT>(0.5) * getInvSpacingPhi() / radius; // phi direction
      }
    }

    // for boundary-r
    for (size_t iZ = 0; iZ < Nz; iZ++) {
      mElectricFieldEr[side](iZ, 0, iPhi) = -1 * (-static_cast<DataT>(0.5) * mPotential[side](iZ, 2, iPhi) + 2 * mPotential[side](iZ, 1, iPhi) - static_cast<DataT>(1.5) * mPotential[side](iZ, 0, iPhi)) * getInvSpacingR();                    // forward difference
      mElectricFieldEr[side](iZ, Nr - 1, iPhi) = -1 * (static_cast<DataT>(1.5) * mPotential[side](iZ, Nr - 1, iPhi) - 2 * mPotential[side](iZ, Nr - 2, iPhi) + static_cast<DataT>(0.5) * mPotential[side](iZ, Nr - 3, iPhi)) * getInvSpacingR(); // backward difference
    }

    for (size_t iR = 0; iR < Nr; iR += Nr - 1) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 1; iZ < Nz - 1; iZ++) {
        mElectricFieldEz[side](iZ, iR, iPhi) = -1 * (mPotential[side](iZ + 1, iR, iPhi) - mPotential[side](iZ - 1, iR, iPhi)) * static_cast<DataT>(0.5) * getInvSpacingZ();                                    // z direction
        mElectricFieldEphi[side](iZ, iR, iPhi) = -1 * (signPlus * mPotential[side](iZ, iR, tmpPlus) - signMinus * mPotential[side](iZ, iR, tmpMinus)) * static_cast<DataT>(0.5) * getInvSpacingPhi() / radius; // phi direction
      }
    }

    // for boundary-z
    for (size_t iR = 0; iR < Nr; ++iR) {
      mElectricFieldEz[side](0, iR, iPhi) = -1 * (-static_cast<DataT>(0.5) * mPotential[side](2, iR, iPhi) + 2 * mPotential[side](1, iR, iPhi) - static_cast<DataT>(1.5) * mPotential[side](0, iR, iPhi)) * getInvSpacingZ();
      mElectricFieldEz[side](Nz - 1, iR, iPhi) = -1 * (static_cast<DataT>(1.5) * mPotential[side](Nz - 1, iR, iPhi) - 2 * mPotential[side](Nz - 2, iR, iPhi) + static_cast<DataT>(0.5) * mPotential[side](Nz - 3, iR, iPhi)) * getInvSpacingZ();
    }

    for (size_t iR = 1; iR < Nr - 1; ++iR) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz; iZ += Nz - 1) {
        mElectricFieldEr[side](iZ, iR, iPhi) = -1 * (mPotential[side](iZ, iR + 1, iPhi) - mPotential[side](iZ, iR - 1, iPhi)) * static_cast<DataT>(0.5) * getInvSpacingR();                                    // r direction
        mElectricFieldEphi[side](iZ, iR, iPhi) = -1 * (signPlus * mPotential[side](iZ, iR, tmpPlus) - signMinus * mPotential[side](iZ, iR, tmpMinus)) * static_cast<DataT>(0.5) * getInvSpacingPhi() / radius; // phi direction
      }
    }

    // corner points for EPhi
    for (size_t iR = 0; iR < Nr; iR += Nr - 1) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz; iZ += Nz - 1) {
        mElectricFieldEphi[side](iZ, iR, iPhi) = -1 * (signPlus * mPotential[side](iZ, iR, tmpPlus) - signMinus * mPotential[side](iZ, iR, tmpMinus)) * static_cast<DataT>(0.5) * getInvSpacingPhi() / radius; // phi direction
      }
    }
  }
  mIsEfieldSet[side] = true;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::calcGlobalDistWithGlobalCorrIterative(const DistCorrInterpolator<DataT, Nz, Nr, Nphi>& globCorr, const int maxIter, const DataT approachZ, const DataT approachR, const DataT approachPhi, const DataT diffCorr)
{
  // for nearest neighbour search see: https://doc.cgal.org/latest/Spatial_searching/Spatial_searching_2searching_with_point_with_info_8cpp-example.html
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Point_3 = Kernel::Point_3;
  using Point_Index = boost::tuple<Point_3, std::array<unsigned int, 3>>;
  using Traits_base = CGAL::Search_traits_3<Kernel>;
  using Traits = CGAL::Search_traits_adapter<Point_Index, CGAL::Nth_of_tuple_property_map<0, Point_Index>, Traits_base>;
  using K_neighbor_search = CGAL::Orthogonal_k_neighbor_search<Traits>;
  using Tree = K_neighbor_search::Tree;

  const Side side = globCorr.getSide();
  const int nPoints = Nz * Nr * Nphi; // maximum number of points

  std::vector<Point_3> points{};
  points.reserve(nPoints);
  std::vector<std::array<unsigned int, 3>> indices{};
  indices.reserve(nPoints);

  // loop over global corrections and calculate the positions of the global correction
  for (unsigned int iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi = getPhiVertex(iPhi);
    for (unsigned int iR = 0; iR < Nr; ++iR) {
      const DataT radius = getRVertex(iR);
      for (unsigned int iZ = 0; iZ < Nz; ++iZ) {
        const DataT z = getZVertex(iZ);

        const DataT globalCorrR = mGlobalCorrdR[side](iZ, iR, iPhi);
        const DataT globalCorrRPhi = mGlobalCorrdRPhi[side](iZ, iR, iPhi);
        const DataT globalCorrZ = mGlobalCorrdZ[side](iZ, iR, iPhi);

        const DataT posRCorr = radius + globalCorrR; // position of global correction
        const DataT posPhiCorr = regulatePhi(phi + globalCorrRPhi / radius);
        const DataT posZCorr = z + globalCorrZ;

        // check if the postion lies in the TPC volume
        if (posRCorr >= getRMin() && posRCorr <= getRMax() && posZCorr >= getZMin() && posZCorr <= getZMax()) {
          // normalize coordinates to gridspacing for searching the nearest neighbour. Otherwise one would have to convert the coordinates to x,y,z
          points.emplace_back(Point_3(posZCorr * getInvSpacingZ(), posRCorr * getInvSpacingR(), posPhiCorr * getInvSpacingPhi()));
          indices.emplace_back(std::array<unsigned int, 3>{iZ, iR, iPhi});
        }
      }
    }
  }

  // Insert data points in the tree
  Tree tree(boost::make_zip_iterator(boost::make_tuple(points.begin(), indices.begin())), boost::make_zip_iterator(boost::make_tuple(points.end(), indices.end())));

#pragma omp parallel for num_threads(sNThreads)
  for (unsigned int iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi = getPhiVertex(iPhi);
    for (unsigned int iR = 0; iR < Nr; ++iR) {
      const DataT radius = getRVertex(iR);
      for (unsigned int iZ = 0; iZ < Nz; ++iZ) {
        const DataT z = getZVertex(iZ);

        // find nearest neighbour
        const Point_3 query(z * getInvSpacingZ(), radius * getInvSpacingR(), phi * getInvSpacingPhi());
        const int neighbors = 1; // we need only one nearest neighbour
        const K_neighbor_search search(tree, query, neighbors);

        const K_neighbor_search::iterator it = search.begin();
        const DataT nearestZ = static_cast<DataT>(boost::get<0>(it->first)[0]) * getGridSpacingZ();
        const DataT nearestR = static_cast<DataT>(boost::get<0>(it->first)[1]) * getGridSpacingR();
        const DataT nearestPhi = static_cast<DataT>(boost::get<0>(it->first)[2]) * getGridSpacingPhi();

        const unsigned int nearestiZ = boost::get<1>(it->first)[0];
        const unsigned int nearestiR = boost::get<1>(it->first)[1];
        const unsigned int nearestiPhi = boost::get<1>(it->first)[2];

        //
        //==========================================================================================
        //==== start algorithm: use tricubic upsampling to numerically approach the query point ====
        //==========================================================================================
        //
        // 1. calculate difference from nearest point to query point with stepwidth factor x
        // and approach the new point
        //
        DataT stepR = (radius - nearestR) * approachR;
        DataT stepZ = (z - nearestZ) * approachZ;
        DataT stepPhi = (phi - nearestPhi) * approachPhi;

        // needed to check for convergence
        DataT lastCorrdR = std::numeric_limits<DataT>::max();
        DataT lastCorrdZ = std::numeric_limits<DataT>::max();
        DataT lastCorrdRPhi = std::numeric_limits<DataT>::max();

        // interpolated global correction
        DataT corrdR = 0;
        DataT corrdRPhi = 0;
        DataT corrdZ = 0;

        for (int iter = 0; iter < maxIter; ++iter) {
          // 2. get new point coordinates
          const DataT rCurrPos = getRVertex(nearestiR) + stepR;
          const DataT zCurrPos = getZVertex(nearestiZ) + stepZ;
          const DataT phiCurrPos = getPhiVertex(nearestiPhi) + stepPhi;

          // interpolate global correction at new point and calculate position of global correction
          // corrdR = globCorr.evalSparsedR(zCurrPos, rCurrPos, phiCurrPos);
          corrdR = globCorr.evaldR(zCurrPos, rCurrPos, phiCurrPos);
          const DataT rNewPos = rCurrPos + corrdR;

          // const DataT corrPhi = globCorr.evalSparsedRPhi(zCurrPos, rCurrPos, phiCurrPos) / rCurrPos;
          const DataT corrPhi = globCorr.evaldRPhi(zCurrPos, rCurrPos, phiCurrPos) / rCurrPos;
          corrdRPhi = corrPhi * rNewPos; // normalize to new r coordinate
          const DataT phiNewPos = phiCurrPos + corrPhi;

          // corrdZ = globCorr.evalSparsedZ(zCurrPos, rCurrPos, phiCurrPos);
          corrdZ = globCorr.evaldZ(zCurrPos, rCurrPos, phiCurrPos);
          const DataT zNewPos = zCurrPos + corrdZ;

          // approach desired coordinate
          stepR += (radius - rNewPos) * approachR;
          stepZ += (z - zNewPos) * approachZ;
          stepPhi += (phi - phiNewPos) * approachPhi;

          // check for convergence
          const DataT diffCorrdR = std::abs(corrdR - lastCorrdR);
          const DataT diffCorrdRZ = std::abs(corrdZ - lastCorrdZ);
          const DataT diffCorrdRPhi = std::abs(corrdRPhi - lastCorrdRPhi);

          // stop algorithm if converged
          if (diffCorrdR < diffCorr && diffCorrdRZ < diffCorr && diffCorrdRPhi < diffCorr) {
            break;
          }

          lastCorrdR = corrdR;
          lastCorrdZ = corrdZ;
          lastCorrdRPhi = corrdRPhi;
        }
        // set global distortions if algorithm converged or iterations exceed max numbers of iterations
        mGlobalDistdR[side](iZ, iR, iPhi) = -corrdR;
        mGlobalDistdRPhi[side](iZ, iR, iPhi) = -corrdRPhi;
        mGlobalDistdZ[side](iZ, iR, iPhi) = -corrdZ;
      }
    }
  }
  // set flag that global distortions are set to true
  mIsGlobalDistSet[side] = true;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
NumericalFields<DataT, Nz, Nr, Nphi> SpaceCharge<DataT, Nz, Nr, Nphi>::getElectricFieldsInterpolator(const Side side) const
{
  if (!mIsEfieldSet[side]) {
    LOGP(warning, "============== E-Fields are not set! ==============\n");
  }
  NumericalFields<DataT, Nz, Nr, Nphi> numFields(mElectricFieldEr[side], mElectricFieldEz[side], mElectricFieldEphi[side], mGrid3D, side);
  return numFields;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
DistCorrInterpolator<DataT, Nz, Nr, Nphi> SpaceCharge<DataT, Nz, Nr, Nphi>::getLocalDistInterpolator(const Side side) const
{
  if (!mIsLocalDistSet[side]) {
    LOGP(warning, "============== local distortions not set! ==============\n");
  }
  DistCorrInterpolator<DataT, Nz, Nr, Nphi> numFields(mLocalDistdR[side], mLocalDistdZ[side], mLocalDistdRPhi[side], mGrid3D, side);
  return numFields;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
DistCorrInterpolator<DataT, Nz, Nr, Nphi> SpaceCharge<DataT, Nz, Nr, Nphi>::getLocalCorrInterpolator(const Side side) const
{
  if (!mIsLocalCorrSet[side]) {
    LOGP(warning, "============== local corrections not set!  ==============\n");
  }
  DistCorrInterpolator<DataT, Nz, Nr, Nphi> numFields(mLocalCorrdR[side], mLocalCorrdZ[side], mLocalCorrdRPhi[side], mGrid3D, side);
  return numFields;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
DistCorrInterpolator<DataT, Nz, Nr, Nphi> SpaceCharge<DataT, Nz, Nr, Nphi>::getGlobalDistInterpolator(const Side side) const
{
  if (!mIsGlobalDistSet[side]) {
    LOGP(warning, "============== global distortions not set ==============\n");
  }
  DistCorrInterpolator<DataT, Nz, Nr, Nphi> numFields(mGlobalDistdR[side], mGlobalDistdZ[side], mGlobalDistdRPhi[side], mGrid3D, side);
  return numFields;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
DistCorrInterpolator<DataT, Nz, Nr, Nphi> SpaceCharge<DataT, Nz, Nr, Nphi>::getGlobalCorrInterpolator(const Side side) const
{
  if (!mIsGlobalCorrSet[side]) {
    LOGP(warning, "============== global corrections not set ==============\n");
  }
  DistCorrInterpolator<DataT, Nz, Nr, Nphi> numFields(mGlobalCorrdR[side], mGlobalCorrdZ[side], mGlobalCorrdRPhi[side], mGrid3D, side);
  return numFields;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
int SpaceCharge<DataT, Nz, Nr, Nphi>::dumpElectricFields(TFile& outf, const Side side) const
{
  if (!mIsEfieldSet[side]) {
    LOGP(warning, "============== E-Fields are not set! returning ==============\n");
    return 0;
  }
  const std::string sideName = getSideName(side);
  const int er = mElectricFieldEr[side].writeToFile(outf, fmt::format("fieldEr_side{}", sideName).data());
  const int ez = mElectricFieldEz[side].writeToFile(outf, fmt::format("fieldEz_side{}", sideName).data());
  const int ephi = mElectricFieldEphi[side].writeToFile(outf, fmt::format("fieldEphi_side{}", sideName).data());
  return er + ez + ephi;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::setElectricFieldsFromFile(TFile& inpf, const Side side)
{
  const std::string sideName = getSideName(side);
  mElectricFieldEr[side].initFromFile(inpf, fmt::format("fieldEr_side{}", sideName).data());
  mElectricFieldEz[side].initFromFile(inpf, fmt::format("fieldEz_side{}", sideName).data());
  mElectricFieldEphi[side].initFromFile(inpf, fmt::format("fieldEphi_side{}", sideName).data());
  mIsEfieldSet[side] = true;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
int SpaceCharge<DataT, Nz, Nr, Nphi>::dumpGlobalDistortions(TFile& outf, const Side side) const
{
  if (!mIsGlobalDistSet[side]) {
    LOGP(warning, "============== global distortions are not set! returning ==============\n");
    return 0;
  }
  const std::string sideName = getSideName(side);
  const int er = mGlobalDistdR[side].writeToFile(outf, fmt::format("distR_side{}", sideName).data());
  const int ez = mGlobalDistdZ[side].writeToFile(outf, fmt::format("distZ_side{}", sideName).data());
  const int ephi = mGlobalDistdRPhi[side].writeToFile(outf, fmt::format("distRphi_side{}", sideName).data());
  return er + ez + ephi;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::setGlobalDistortionsFromFile(TFile& inpf, const Side side)
{
  mIsGlobalDistSet[side] = true;
  const std::string sideName = getSideName(side);
  mGlobalDistdR[side].initFromFile(inpf, fmt::format("distR_side{}", sideName).data());
  mGlobalDistdZ[side].initFromFile(inpf, fmt::format("distZ_side{}", sideName).data());
  mGlobalDistdRPhi[side].initFromFile(inpf, fmt::format("distRphi_side{}", sideName).data());
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
int SpaceCharge<DataT, Nz, Nr, Nphi>::dumpGlobalCorrections(TFile& outf, const Side side) const
{
  if (!mIsGlobalCorrSet[side]) {
    LOGP(warning, "============== global corrections are not set! returning ==============\n");
    return 0;
  }
  const std::string sideName = getSideName(side);
  const int er = mGlobalCorrdR[side].writeToFile(outf, fmt::format("corrR_side{}", sideName).data());
  const int ez = mGlobalCorrdZ[side].writeToFile(outf, fmt::format("corrZ_side{}", sideName).data());
  const int ephi = mGlobalCorrdRPhi[side].writeToFile(outf, fmt::format("corrRPhi_side{}", sideName).data());
  return er + ez + ephi;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::setGlobalCorrectionsFromFile(TFile& inpf, const Side side)
{
  mIsGlobalCorrSet[side] = true;
  const std::string sideName = getSideName(side);
  mGlobalCorrdR[side].initFromFile(inpf, fmt::format("corrR_side{}", sideName).data());
  mGlobalCorrdZ[side].initFromFile(inpf, fmt::format("corrZ_side{}", sideName).data());
  mGlobalCorrdRPhi[side].initFromFile(inpf, fmt::format("corrRPhi_side{}", sideName).data());
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
int SpaceCharge<DataT, Nz, Nr, Nphi>::dumpLocalCorrections(TFile& outf, const Side side) const
{
  if (!mIsLocalCorrSet[side]) {
    LOGP(warning, "============== local corrections are not set! returning ==============\n");
    return 0;
  }
  const std::string sideName = getSideName(side);
  const int lCorrdR = mLocalCorrdR[side].writeToFile(outf, fmt::format("lcorrR_side{}", sideName).data());
  const int lCorrdZ = mLocalCorrdZ[side].writeToFile(outf, fmt::format("lcorrZ_side{}", sideName).data());
  const int lCorrdRPhi = mLocalCorrdRPhi[side].writeToFile(outf, fmt::format("lcorrRPhi_side{}", sideName).data());
  return lCorrdR + lCorrdZ + lCorrdRPhi;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::setLocalCorrectionsFromFile(TFile& inpf, const Side side)
{
  const std::string sideName = getSideName(side);
  const bool lCorrdR = mLocalCorrdR[side].initFromFile(inpf, fmt::format("lcorrR_side{}", sideName).data());
  const bool lCorrdZ = mLocalCorrdZ[side].initFromFile(inpf, fmt::format("lcorrZ_side{}", sideName).data());
  const bool lCorrdRPhi = mLocalCorrdRPhi[side].initFromFile(inpf, fmt::format("lcorrRPhi_side{}", sideName).data());
  if (lCorrdR && lCorrdZ && lCorrdRPhi) {
    mIsLocalCorrSet[side] = true;
  } else {
    mIsLocalCorrSet[side] = false;
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
int SpaceCharge<DataT, Nz, Nr, Nphi>::dumpLocalDistortions(TFile& outf, const Side side) const
{
  if (!mIsLocalDistSet[side]) {
    LOGP(warning, "============== local distortions are not set! returning ==============\n");
    return 0;
  }
  const std::string sideName = getSideName(side);
  const int lDistdR = mLocalDistdR[side].writeToFile(outf, fmt::format("ldistR_side{}", sideName).data());
  const int lDistdZ = mLocalDistdZ[side].writeToFile(outf, fmt::format("ldistZ_side{}", sideName).data());
  const int lDistdRPhi = mLocalDistdRPhi[side].writeToFile(outf, fmt::format("ldistRPhi_side{}", sideName).data());
  return lDistdR + lDistdZ + lDistdRPhi;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::setLocalDistortionsFromFile(TFile& inpf, const Side side)
{
  const std::string sideName = getSideName(side);
  const bool lDistdR = mLocalDistdR[side].initFromFile(inpf, fmt::format("ldistR_side{}", sideName).data());
  const bool lDistdZ = mLocalDistdZ[side].initFromFile(inpf, fmt::format("ldistZ_side{}", sideName).data());
  const bool lDistdRPhi = mLocalDistdRPhi[side].initFromFile(inpf, fmt::format("ldistRPhi_side{}", sideName).data());

  if (lDistdR && lDistdZ && lDistdRPhi) {
    mIsLocalDistSet[side] = true;
  } else {
    mIsLocalDistSet[side] = false;
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::fillChargeDensityFromHisto(TFile& fInp, const char* name)
{
  TH3* hisSCDensity3D = (TH3*)fInp.Get(name);
  TH3D hRebin = rebinDensityHisto(*hisSCDensity3D);
  for (int side = Side::A; side < SIDES; ++side) {
    for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
      for (size_t iR = 0; iR < Nr; ++iR) {
        for (size_t iZ = 0; iZ < Nz; ++iZ) {
          const size_t zBin = side == Side::A ? Nz + iZ + 1 : Nz - iZ; // TODO CHECK THIS!
          mDensity[side](iZ, iR, iPhi) = hRebin.GetBinContent(iPhi, iR, zBin);
        }
      }
    }
    mIsChargeSet[side] = true;
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
TH3D SpaceCharge<DataT, Nz, Nr, Nphi>::rebinDensityHisto(const TH3& hOrig) const
{
  TH3D hRebin{};

  const int nBinsPhiNew = Nphi;
  const int nBinsRNew = Nr;
  const int nBinsZNew = 2 * Nz;

  const auto phiLow = hOrig.GetXaxis()->GetBinLowEdge(1);
  const auto phiUp = hOrig.GetXaxis()->GetBinUpEdge(hOrig.GetNbinsX());
  const auto rLow = hOrig.GetYaxis()->GetBinLowEdge(1);
  const auto rUp = hOrig.GetYaxis()->GetBinUpEdge(hOrig.GetNbinsY());
  const auto zLow = hOrig.GetZaxis()->GetBinLowEdge(1);
  const auto zUp = hOrig.GetZaxis()->GetBinUpEdge(hOrig.GetNbinsZ());
  hRebin.SetBins(nBinsPhiNew, phiLow, phiUp, nBinsRNew, rLow, rUp, nBinsZNew, zLow, zUp);

  for (int iBinPhi = 1; iBinPhi <= nBinsPhiNew; ++iBinPhi) {
    const auto phiLowEdge = hRebin.GetXaxis()->GetBinLowEdge(iBinPhi);
    const auto phiUpEdge = hRebin.GetXaxis()->GetBinUpEdge(iBinPhi);

    const int phiLowBinOrig = hOrig.GetXaxis()->FindBin(phiLowEdge);
    const int phiUpBinOrig = hOrig.GetXaxis()->FindBin(phiUpEdge);

    // calculate the weights (area of original bin lies in the new bin / binwidthOrig) of the first and last bins
    const auto binWidthPhiOrig = hOrig.GetXaxis()->GetBinWidth(phiLowBinOrig);
    const auto lowerBinWeightPhi = std::abs(phiLowEdge - hOrig.GetXaxis()->GetBinUpEdge(phiLowBinOrig)) / binWidthPhiOrig;
    const auto upperBinWeightPhi = std::abs(phiUpEdge - hOrig.GetXaxis()->GetBinLowEdge(phiUpBinOrig)) / binWidthPhiOrig;

    for (int iBinR = 1; iBinR <= nBinsRNew; ++iBinR) {
      const auto rLowEdge = hRebin.GetYaxis()->GetBinLowEdge(iBinR);
      const auto rUpEdge = hRebin.GetYaxis()->GetBinUpEdge(iBinR);

      const int rLowBinOrig = hOrig.GetYaxis()->FindBin(rLowEdge);
      const int rUpBinOrig = hOrig.GetYaxis()->FindBin(rUpEdge);

      // calculate the weights (area of original bin lies in the new bin / binwidthOrig) of the first and last bins
      const auto binWidthROrig = hOrig.GetYaxis()->GetBinWidth(rLowBinOrig);
      const auto lowerBinWeightR = std::abs(rLowEdge - hOrig.GetYaxis()->GetBinUpEdge(rLowBinOrig)) / binWidthROrig;
      const auto upperBinWeightR = std::abs(rUpEdge - hOrig.GetYaxis()->GetBinLowEdge(rUpBinOrig)) / binWidthROrig;

      for (int iBinZ = 1; iBinZ <= nBinsZNew; ++iBinZ) {
        const auto zLowEdge = hRebin.GetZaxis()->GetBinLowEdge(iBinZ);
        const auto zUpEdge = hRebin.GetZaxis()->GetBinUpEdge(iBinZ);
        const auto zCenter = hRebin.GetZaxis()->GetBinCenter(iBinZ);

        int zLowBinOrig = hOrig.GetZaxis()->FindBin(zLowEdge);
        int zUpBinOrig = hOrig.GetZaxis()->FindBin(zUpEdge);
        const int currside = getSide(zCenter); // set the side of the current z-bin
        // get the side of the lowest and uppest bin from the orig histo
        const int sideLowOrig = getSide(hOrig.GetZaxis()->GetBinCenter(zLowBinOrig));
        const int sideUpOrig = getSide(hOrig.GetZaxis()->GetBinCenter(zUpBinOrig));

        // make bounds/side check of the zLowBinOrig and zUpBinOrig bins. They must be on the same side as the currside!!!
        if (currside != sideLowOrig && zLowBinOrig != zUpBinOrig) {
          // if the lower bins from the orig histo are not on the same side as the rebinned increase the binnumber until they are on the same side
          bool notequal = true;
          do {
            zLowBinOrig += 1;
            if (zLowBinOrig > zUpBinOrig) {
              LOGP(WARNING, "SOMETHING WENT WRONG: SETTING BINS TO: {}", zUpBinOrig);
              zLowBinOrig = zUpBinOrig;
              notequal = false;
            }
            const int sideTmp = getSide(hOrig.GetZaxis()->GetBinCenter(zLowBinOrig));
            if (sideTmp == currside) {
              notequal = false;
            }
          } while (notequal);
        }

        if (currside != sideUpOrig && zLowBinOrig != zUpBinOrig) {
          // if the upper bins from the orig histo are not on the same side as the rebinned increase the binnumber until they are on the same side
          bool notequal = true;
          do {
            zUpBinOrig -= 1;
            if (zUpBinOrig < zLowBinOrig) {
              LOGP(WARNING, "SOMETHING WENT WRONG: SETTING BINS TO: {}", zLowBinOrig);
              zUpBinOrig = zLowBinOrig;
              notequal = false;
            }
            const int sideTmp = getSide(hOrig.GetZaxis()->GetBinCenter(zUpBinOrig));
            if (sideTmp == currside) {
              notequal = false;
            }
          } while (notequal);
        }

        const auto binWidthZOrig = hOrig.GetZaxis()->GetBinWidth(zLowBinOrig);
        const auto lowerBinWeightZ = std::abs(zLowEdge - hOrig.GetZaxis()->GetBinUpEdge(zLowBinOrig)) / binWidthZOrig;
        const auto upperBinWeightZ = std::abs(zUpEdge - hOrig.GetZaxis()->GetBinLowEdge(zUpBinOrig)) / binWidthZOrig;

        // get the mean value of the original histogram of the found bin range
        DataT sum = 0;
        DataT sumW = 0;
        for (int iPhi = phiLowBinOrig; iPhi <= phiUpBinOrig; ++iPhi) {
          DataT weightPhi = 1;
          if (iPhi == phiLowBinOrig) {
            weightPhi = lowerBinWeightPhi;
          } else if (iPhi == phiUpBinOrig) {
            weightPhi = upperBinWeightPhi;
          }

          for (int iR = rLowBinOrig; iR <= rUpBinOrig; ++iR) {
            DataT weightR = 1;
            if (iR == rLowBinOrig) {
              weightR = lowerBinWeightR;
            } else if (iR == rUpBinOrig) {
              weightR = upperBinWeightR;
            }

            for (int iZ = zLowBinOrig; iZ <= zUpBinOrig; ++iZ) {
              DataT weightZ = 1;
              if (iZ == zLowBinOrig) {
                weightZ = lowerBinWeightZ;
              } else if (iZ == zUpBinOrig) {
                weightZ = upperBinWeightZ;
              }
              const auto val = hOrig.GetBinContent(iPhi, iR, iZ);
              // if(val==0){
              // what to do now???
              // }
              const auto totalWeight = weightPhi * weightR * weightZ;
              sum += val * totalWeight;
              sumW += totalWeight;
            }
          }
        }
        sum /= sumW;
        hRebin.SetBinContent(iBinPhi, iBinR, iBinZ, sum);
      }
    }
  }
  return hRebin;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename ElectricFields>
void SpaceCharge<DataT, Nz, Nr, Nphi>::calcLocalDistortionsCorrections(const SpaceCharge<DataT, Nz, Nr, Nphi>::Type type, const ElectricFields& formulaStruct)
{
  const Side side = formulaStruct.getSide();
  // calculate local distortions/corrections for each vertex in the tpc
#pragma omp parallel for num_threads(sNThreads)
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi = getPhiVertex(iPhi);
    for (size_t iR = 0; iR < Nr; ++iR) {
      const DataT radius = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz - 1; ++iZ) {
        // set z coordinate depending on distortions or correction calculation
        const DataT z0 = type == Type::Corrections ? getZVertex(iZ + 1) : getZVertex(iZ);
        const DataT z1 = type == Type::Corrections ? getZVertex(iZ) : getZVertex(iZ + 1);

        DataT drTmp = 0;   // local distortion dR
        DataT dPhiTmp = 0; // local distortion dPhi (multiplication with R has to be done at the end)
        DataT dzTmp = 0;   // local distortion dZ

        const DataT stepSize = (z1 - z0) / mSteps; // the distortions are calculated by leting the elctron drift this distance in z direction
        for (int iter = 0; iter < mSteps; ++iter) {
          const DataT z0Tmp = (z0 + iter * stepSize + dzTmp); // starting z position
          const DataT z1Tmp = (z0Tmp + stepSize);             // electron drifts from z0Tmp to z1Tmp

          DataT ddR = 0;   // distortion dR for drift from z0Tmp to z1Tmp
          DataT ddPhi = 0; // distortion dPhi for drift from z0Tmp to z1Tmp
          DataT ddZ = 0;   // distortion dZ for drift from z0Tmp to z1Tmp

          const DataT radiusTmp = regulateR(radius + drTmp); // current radial position
          const DataT phiTmp = regulatePhi(phi + dPhiTmp);   // current phi position

          // calculate distortions/corrections
          calcDistCorr(radiusTmp, phiTmp, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct, true);

          // add temp distortions to local distortions
          drTmp += ddR;
          dPhiTmp += ddPhi;
          dzTmp += ddZ;
        }

        // store local distortions/corrections
        switch (type) {
          case Type::Corrections:
            mLocalCorrdR[side](iZ + 1, iR, iPhi) = drTmp;
            mLocalCorrdRPhi[side](iZ + 1, iR, iPhi) = dPhiTmp * radius;
            mLocalCorrdZ[side](iZ + 1, iR, iPhi) = dzTmp;
            break;

          case Type::Distortions:
            mLocalDistdR[side](iZ, iR, iPhi) = drTmp;
            mLocalDistdRPhi[side](iZ, iR, iPhi) = dPhiTmp * radius;
            mLocalDistdZ[side](iZ, iR, iPhi) = dzTmp;
            break;
        }
      }
      //extrapolate local distortion/correction to last/first bin using legendre polynoms with x0=0, x1=1, x2=2 and x=-1. This has to be done to ensure correct interpolation in the last,second last/first,second bin!
      switch (type) {
        case Type::Corrections:
          mLocalCorrdR[side](0, iR, iPhi) = 3 * (mLocalCorrdR[side](1, iR, iPhi) - mLocalCorrdR[side](2, iR, iPhi)) + mLocalCorrdR[side](3, iR, iPhi);
          mLocalCorrdRPhi[side](0, iR, iPhi) = 3 * (mLocalCorrdRPhi[side](1, iR, iPhi) - mLocalCorrdRPhi[side](2, iR, iPhi)) + mLocalCorrdRPhi[side](3, iR, iPhi);
          mLocalCorrdZ[side](0, iR, iPhi) = 3 * (mLocalCorrdZ[side](1, iR, iPhi) - mLocalCorrdZ[side](2, iR, iPhi)) + mLocalCorrdZ[side](3, iR, iPhi);
          break;

        case Type::Distortions:
          mLocalDistdR[side](Nz - 1, iR, iPhi) = 3 * (mLocalDistdR[side](Nz - 2, iR, iPhi) - mLocalDistdR[side](Nz - 3, iR, iPhi)) + mLocalDistdR[side](Nz - 4, iR, iPhi);
          mLocalDistdRPhi[side](Nz - 1, iR, iPhi) = 3 * (mLocalDistdRPhi[side](Nz - 2, iR, iPhi) - mLocalDistdRPhi[side](Nz - 3, iR, iPhi)) + mLocalDistdRPhi[side](Nz - 4, iR, iPhi);
          mLocalDistdZ[side](Nz - 1, iR, iPhi) = 3 * (mLocalDistdZ[side](Nz - 2, iR, iPhi) - mLocalDistdZ[side](Nz - 3, iR, iPhi)) + mLocalDistdZ[side](Nz - 4, iR, iPhi);
          break;
      }
    }
  }
  switch (type) {
    case Type::Corrections:
      mIsLocalCorrSet[side] = true;
      break;
    case Type::Distortions:
      mIsLocalDistSet[side] = true;
      break;
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename Fields>
void SpaceCharge<DataT, Nz, Nr, Nphi>::calcGlobalDistortions(const Fields& formulaStruct)
{
  const Side side = formulaStruct.getSide();
  const DataT stepSize = formulaStruct.getID() == 2 ? getGridSpacingZ() : getGridSpacingZ() / mSteps; // if one used local distortions then no smaller stepsize is needed. if electric fields are used then smaller stepsize can be used
  // loop over tpc volume and let the electron drift from each vertex to the readout of the tpc
#pragma omp parallel for num_threads(sNThreads)
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi0 = getPhiVertex(iPhi);
    for (size_t iR = 0; iR < Nr; ++iR) {
      const DataT r0 = getRVertex(iR);
      for (size_t iZ = 0; iZ < Nz - 1; ++iZ) {
        const DataT z0 = getZVertex(iZ); // the electron starts at z0, r0, phi0
        DataT drDist = 0.0;              // global distortion dR
        DataT dPhiDist = 0.0;            // global distortion dPhi (multiplication with R has to be done at the end)
        DataT dzDist = 0.0;              // global distortion dZ
        int iter = 0;

        for (;;) {
          const DataT z0Tmp = z0 + dzDist + iter * stepSize; // starting z position
          const DataT z1Tmp = regulateZ(z0Tmp + stepSize);   // electron drifts from z0Tmp to z1Tmp
          const DataT radius = regulateR(r0 + drDist);       // current radial position of the electron
          const DataT phi = regulatePhi(phi0 + dPhiDist);    // current phi position of the electron

          DataT ddR = 0;   // distortion dR for drift from z0Tmp to z1Tmp
          DataT ddPhi = 0; // distortion dPhi for drift from z0Tmp to z1Tmp
          DataT ddZ = 0;   // distortion dZ for drift from z0Tmp to z1Tmp

          // get the distortion from interpolation of local distortions or calculate distortions with the electric field
          processGlobalDistCorr(radius, phi, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct);

          // if one uses local distortions the interpolated value for the last bin has to be scaled.
          // This has to be done because of the interpolated value is defined for a drift length of one z bin, but in the last bin the distance to the readout can be smaller than one z bin.
          if (formulaStruct.getID() == 2 && z1Tmp >= getZMax()) {
            const DataT fac = (getZMax() - z0Tmp) * getInvSpacingZ();
            ddR *= fac;
            ddZ *= fac;
            ddPhi *= fac;
          }

          // add local distortions to global distortions
          drDist += ddR;
          dPhiDist += ddPhi;
          dzDist += ddZ;

          // set loop to exit if the readout is reached and approximate distortion of 'missing' (one never ends exactly on the readout: z1Tmp + ddZ != ZMAX) drift distance.
          // approximation is done by the current calculated values of the distortions and scaled linear to the 'missing' distance.
          if (z1Tmp >= getZMax()) {
            const DataT endPoint = z1Tmp + ddZ;
            const DataT deltaZ = getZMax() - endPoint; // distance from last point to read out
            const DataT diff = endPoint - z0Tmp;
            const DataT fac = diff != 0 ? deltaZ / diff : 0; // approximate the distortions for the 'missing' distance deltaZ
            drDist += ddR * fac;
            dPhiDist += ddPhi * fac;
            dzDist += ddZ * fac;
            break;
          }
          ++iter;
        }
        // store global distortions
        mGlobalDistdR[side](iZ, iR, iPhi) = drDist;
        mGlobalDistdRPhi[side](iZ, iR, iPhi) = dPhiDist * r0;
        mGlobalDistdZ[side](iZ, iR, iPhi) = dzDist;
      }
    }
  }
  // set flag that global distortions are set to true
  mIsGlobalDistSet[side] = true;
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename Formulas>
void SpaceCharge<DataT, Nz, Nr, Nphi>::calcGlobalCorrections(const Formulas& formulaStruct)
{
  const Side side = formulaStruct.getSide();
  const int iSteps = formulaStruct.getID() == 2 ? 1 : mSteps; // if one used local corrections no step width is needed. since it is already used for calculation of the local corrections
  const DataT stepSize = -getGridSpacingZ() / iSteps;
  // loop over tpc volume and let the electron drift from each vertex to the readout of the tpc
#pragma omp parallel for num_threads(sNThreads)
  for (size_t iPhi = 0; iPhi < Nphi; ++iPhi) {
    const DataT phi0 = getPhiVertex(iPhi);
    for (size_t iR = 0; iR < Nr; ++iR) {
      const DataT r0 = getRVertex(iR);
      DataT drCorr = 0;
      DataT dPhiCorr = 0;
      DataT dzCorr = 0;

      // start at the readout and follow electron towards central electrode
      for (size_t iZ = Nz - 1; iZ >= 1; --iZ) {
        const DataT z0 = getZVertex(iZ); // the electron starts at z0, r0, phi0
        // flag which is set when the central electrode is reached. if the central electrode is reached the calculation of the global corrections is aborted and the value set is the last calculated value.
        bool centralElectrodeReached = false;
        for (int iter = 0; iter < iSteps; ++iter) {
          if (centralElectrodeReached) {
            break;
          }
          const DataT radius = regulateR(r0 + drCorr);       // current radial position of the electron
          const DataT phi = regulatePhi(phi0 + dPhiCorr);    // current phi position of the electron
          const DataT z0Tmp = z0 + dzCorr + iter * stepSize; // starting z position
          const DataT z1Tmp = regulateZ(z0Tmp + stepSize);   // follow electron from z0Tmp to z1Tmp

          DataT ddR = 0;   // distortion dR for z0Tmp to z1Tmp
          DataT ddPhi = 0; // distortion dPhi for z0Tmp to z1Tmp
          DataT ddZ = 0;   // distortion dZ for z0Tmp to z1Tmp

          // get the distortion from interpolation of local distortions or calculate distortions with the electric field
          processGlobalDistCorr(radius, phi, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct);

          // if one uses local corrections the interpolated value for the first bin has to be scaled.
          // This has to be done because of the interpolated value is defined for a drift length of one z bin, but in the first bin the distance to the readout can be smaller than one z bin.
          if (formulaStruct.getID() == 2 && z1Tmp <= getZMin()) {
            const DataT fac = (z0Tmp - getZMin()) * getInvSpacingZ();
            ddR *= fac;
            ddZ *= fac;
            ddPhi *= fac;
          }

          // add local corrections to global corrections
          drCorr += ddR;
          dPhiCorr += ddPhi;
          dzCorr += ddZ;

          // set loop to exit if the central electrode is reached and approximate correction of 'missing' (one never ends exactly on the central electrode: z1Tmp + ddZ != ZMIN) distance.
          // approximation is done by the current calculated values of the corrections and scaled linear to the 'missing' distance deltaZ. (NOT TESTED)
          if (z1Tmp <= getZMin()) {
            const DataT endPoint = z1Tmp + ddZ;
            const DataT deltaZ = endPoint - getZMin();
            const DataT diff = z0Tmp - endPoint;
            const DataT fac = diff != 0 ? deltaZ / diff : 0; // approximate the distortions for the 'missing' distance deltaZ
            drCorr += ddR * fac;
            dPhiCorr += ddPhi * fac;
            dzCorr += ddZ * fac;
            centralElectrodeReached = true;
            break;
          }
        }
        // store global corrections
        mGlobalCorrdR[side](iZ - 1, iR, iPhi) = drCorr;
        mGlobalCorrdRPhi[side](iZ - 1, iR, iPhi) = dPhiCorr * r0;
        mGlobalCorrdZ[side](iZ - 1, iR, iPhi) = dzCorr;
      }
    }
  }
  // set flag that global corrections are set to true
  mIsGlobalCorrSet[side] = true;
}

using DataTD = double;
template class o2::tpc::SpaceCharge<DataTD, 17, 17, 90>;
template class o2::tpc::SpaceCharge<DataTD, 33, 33, 180>;
template class o2::tpc::SpaceCharge<DataTD, 65, 65, 180>;
template class o2::tpc::SpaceCharge<DataTD, 129, 129, 180>;

// 129*129*180
using NumFields129D = NumericalFields<DataTD, 129, 129, 180>;
using AnaFields129D = AnalyticalFields<DataTD>;
using DistCorrInterp129D = DistCorrInterpolator<DataTD, 129, 129, 180>;
using O2TPCSpaceCharge3DCalc129D = SpaceCharge<DataTD, 129, 129, 180>;

template void O2TPCSpaceCharge3DCalc129D::calcLocalDistortionsCorrections(const O2TPCSpaceCharge3DCalc129D::Type, const NumFields129D&);
template void O2TPCSpaceCharge3DCalc129D::calcLocalDistortionsCorrections(const O2TPCSpaceCharge3DCalc129D::Type, const AnaFields129D&);
template void O2TPCSpaceCharge3DCalc129D::calcGlobalCorrections(const NumFields129D&);
template void O2TPCSpaceCharge3DCalc129D::calcGlobalCorrections(const AnaFields129D&);
template void O2TPCSpaceCharge3DCalc129D::calcGlobalCorrections(const DistCorrInterp129D&);

template void O2TPCSpaceCharge3DCalc129D::calcGlobalDistortions(const NumFields129D&);
template void O2TPCSpaceCharge3DCalc129D::calcGlobalDistortions(const AnaFields129D&);
template void O2TPCSpaceCharge3DCalc129D::calcGlobalDistortions(const DistCorrInterp129D&);

// 33*33*180
using NumFields33D = NumericalFields<DataTD, 33, 33, 180>;
using AnaFields33D = AnalyticalFields<DataTD>;
using DistCorrInterp33D = DistCorrInterpolator<DataTD, 33, 33, 180>;
using O2TPCSpaceCharge3DCalc33D = SpaceCharge<DataTD, 33, 33, 180>;

template void O2TPCSpaceCharge3DCalc33D::calcLocalDistortionsCorrections(const O2TPCSpaceCharge3DCalc33D::Type, const NumFields33D&);
template void O2TPCSpaceCharge3DCalc33D::calcLocalDistortionsCorrections(const O2TPCSpaceCharge3DCalc33D::Type, const AnaFields33D&);
template void O2TPCSpaceCharge3DCalc33D::calcGlobalCorrections(const NumFields33D&);
template void O2TPCSpaceCharge3DCalc33D::calcGlobalCorrections(const AnaFields33D&);
template void O2TPCSpaceCharge3DCalc33D::calcGlobalCorrections(const DistCorrInterp33D&);

template void O2TPCSpaceCharge3DCalc33D::calcGlobalDistortions(const NumFields33D&);
template void O2TPCSpaceCharge3DCalc33D::calcGlobalDistortions(const AnaFields33D&);
template void O2TPCSpaceCharge3DCalc33D::calcGlobalDistortions(const DistCorrInterp33D&);

// 65*65*180
using NumFields65D = NumericalFields<DataTD, 65, 65, 180>;
using AnaFields65D = AnalyticalFields<DataTD>;
using DistCorrInterp65D = DistCorrInterpolator<DataTD, 65, 65, 180>;
using O2TPCSpaceCharge3DCalc65D = SpaceCharge<DataTD, 65, 65, 180>;

template void O2TPCSpaceCharge3DCalc65D::calcLocalDistortionsCorrections(const O2TPCSpaceCharge3DCalc65D::Type, const NumFields65D&);
template void O2TPCSpaceCharge3DCalc65D::calcLocalDistortionsCorrections(const O2TPCSpaceCharge3DCalc65D::Type, const AnaFields65D&);
template void O2TPCSpaceCharge3DCalc65D::calcGlobalCorrections(const NumFields65D&);
template void O2TPCSpaceCharge3DCalc65D::calcGlobalCorrections(const AnaFields65D&);
template void O2TPCSpaceCharge3DCalc65D::calcGlobalCorrections(const DistCorrInterp65D&);

template void O2TPCSpaceCharge3DCalc65D::calcGlobalDistortions(const NumFields65D&);
template void O2TPCSpaceCharge3DCalc65D::calcGlobalDistortions(const AnaFields65D&);
template void O2TPCSpaceCharge3DCalc65D::calcGlobalDistortions(const DistCorrInterp65D&);

// 17*17*90
using NumFields17D = NumericalFields<DataTD, 17, 17, 90>;
using AnaFields17D = AnalyticalFields<DataTD>;
using DistCorrInterp17D = DistCorrInterpolator<DataTD, 17, 17, 90>;
using O2TPCSpaceCharge3DCalc17D = SpaceCharge<DataTD, 17, 17, 90>;

template void O2TPCSpaceCharge3DCalc17D::calcLocalDistortionsCorrections(const O2TPCSpaceCharge3DCalc17D::Type, const NumFields17D&);
template void O2TPCSpaceCharge3DCalc17D::calcLocalDistortionsCorrections(const O2TPCSpaceCharge3DCalc17D::Type, const AnaFields17D&);
template void O2TPCSpaceCharge3DCalc17D::calcGlobalCorrections(const NumFields17D&);
template void O2TPCSpaceCharge3DCalc17D::calcGlobalCorrections(const AnaFields17D&);
template void O2TPCSpaceCharge3DCalc17D::calcGlobalCorrections(const DistCorrInterp17D&);

template void O2TPCSpaceCharge3DCalc17D::calcGlobalDistortions(const NumFields17D&);
template void O2TPCSpaceCharge3DCalc17D::calcGlobalDistortions(const AnaFields17D&);
template void O2TPCSpaceCharge3DCalc17D::calcGlobalDistortions(const DistCorrInterp17D&);
