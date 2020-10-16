// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  SpaceCharge.h
/// \brief This class contains the algorithms for calculation the distortions and corrections
///
/// \author  Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_TPC_SpaceCharge_H_
#define ALICEO2_TPC_SpaceCharge_H_

#include "TPCSpacecharge/TriCubic.h"
#include "TPCSpacecharge/PoissonSolver.h"
#include "TPCSpacecharge/SpaceChargeStructs.h"
#include "TPCSpacecharge/RegularGrid3D.h"
#include "TPCSpacecharge/DataContainer3D.h"

#include "DataFormatsTPC/Defs.h"
#include "Framework/Logger.h"

// Root includes
#include "TF1.h" /// for numerical intergration only
#include "TH3.h"

#ifdef WITH_OPENMP
#include <omp.h>
#endif

// for nearest neighbour search. needed for the calculation of the global distortions by iterative method using global corrections
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <boost/iterator/zip_iterator.hpp>

namespace o2
{
namespace tpc
{

/// \class SpaceCharge
/// this class provides the algorithms for calculating the global distortions and corrections from the space charge density.
/// The calculations can be done by a realistic space charge histogram as an input or by an analytical formula.
/// An example of of the usage can be found in 'macro/calculateDistortionsCorrections.C'

/// \tparam DataT the data type which is used during the calculations
/// \tparam Nz number of vertices in z direction
/// \tparam Nr number of vertices in r direction
/// \tparam Nphi number of vertices in phi direction
template <typename DataT = double, size_t Nz = 129, size_t Nr = 129, size_t Nphi = 180>
class SpaceCharge
{
  using RegularGrid = RegularGrid3D<DataT, Nz, Nr, Nphi>;
  using DataContainer = DataContainer3D<DataT, Nz, Nr, Nphi>;
  using GridProp = GridProperties<DataT, Nr, Nz, Nphi>;

 public:
  /// default constructor
  SpaceCharge() = default;

  /// numerical integration strategys
  enum class IntegrationStrategy { Trapezoidal = 0,     ///< trapezoidal integration (https://en.wikipedia.org/wiki/Trapezoidal_rule). straight electron drift line assumed: z0->z1, r0->r0, phi0->phi0
                                   Simpson = 1,         ///< simpon integration. see: https://en.wikipedia.org/wiki/Simpson%27s_rule. straight electron drift line assumed: z0->z1, r0->r0, phi0->phi0
                                   Root = 2,            ///< Root integration. straight electron drift line assumed: z0->z1, r0->r0, phi0->phi0
                                   SimpsonIterative = 3 ///< simpon integration, but using an iterative method to approximate the drift path. No straight electron drift line assumed: z0->z1, r0->r1, phi0->phi1
  };

  enum class Type {
    Distortions = 0, ///< distortions
    Corrections = 1  ///< corrections
  };

  /// step 0: set the charge density from TH3 histogram containing the space charge density
  /// \param fInp input file containing a histogram for the space charge density
  /// \param name the name of the space charge density histogram in the file
  void fillChargeDensityFromHisto(TFile& fInp, const char* name);

  /// step 0: this function fills the internal storage for the charge density using an analytical formula
  /// \param formulaStruct struct containing a method to evaluate the density
  void setChargeDensityFromFormula(const AnalyticalFields<DataT>& formulaStruct);

  /// step 0: this function fills the boundary of the potential using an analytical formula. The boundary is used in the PoissonSolver.
  /// \param formulaStruct struct containing a method to evaluate the potential
  void setPotentialBoundaryFromFormula(const AnalyticalFields<DataT>& formulaStruct);

  /// step 0: this function fills the potential using an analytical formula
  /// \param formulaStruct struct containing a method to evaluate the potential
  void setPotentialFromFormula(const AnalyticalFields<DataT>& formulaStruct);

  /// step 1: use the O2TPCPoissonSolver class to numerically calculate the potential with set space charge density and boundary conditions from potential
  /// \param side side of the TPC
  /// \param maxIteration maximum number of iterations used in the poisson solver
  /// \param stoppingConvergence stopping criterion used in the poisson solver
  /// \param symmetry use symmetry or not in the poisson solver
  void poissonSolver(const Side side, const int maxIteration = 300, const DataT stoppingConvergence = 1e-6, const int symmetry = 0);

  /// step 2: calculate numerically the electric field from the potential
  /// \param side side of the TPC
  void calcEField(const Side side);

  /// step 2a: set the electric field from an analytical formula
  /// \param formulaStruct struct containing a method to evaluate the electric fields
  void setEFieldFromFormula(const AnalyticalFields<DataT>& formulaStruct);

  /// step 3: calculate the local distortions and corrections with an electric field
  /// \param type calculate local corrections or local distortions: type = o2::tpc::SpaceCharge<>::Type::Distortions or o2::tpc::SpaceCharge<>::Type::Corrections
  /// \param formulaStruct struct containing a method to evaluate the electric field Er, Ez, Ephi (analytical formula or by TriCubic interpolator)
  template <typename ElectricFields = AnalyticalFields<DataT>>
  void calcLocalDistortionsCorrections(const Type type, const ElectricFields& formulaStruct);

  /// step 4: calculate global corrections by using the electric field or the local corrections
  /// \param formulaStruct struct containing a method to evaluate the electric field Er, Ez, Ephi or the local corrections
  template <typename Fields = AnalyticalFields<DataT>>
  void calcGlobalCorrections(const Fields& formulaStruct);

  /// step 5: calculate global distortions by using the electric field or the local distortions (SLOW)
  /// \param formulaStruct struct containing a method to evaluate the electric field Er, Ez, Ephi or the local distortions
  template <typename Fields = AnalyticalFields<DataT>>
  void calcGlobalDistortions(const Fields& formulaStruct);

  /// step 5: calculate global distortions using the global corrections (FAST)
  /// \param globCorr interpolator for global corrections
  /// \param maxIter maximum iterations per global distortion
  /// \param approachZ when the difference between the desired z coordinate and the position of the global correction is deltaZ, approach the desired z coordinate by deltaZ * \p approachZ.
  /// \param approachR when the difference between the desired r coordinate and the position of the global correction is deltaR, approach the desired r coordinate by deltaR * \p approachR.
  /// \param approachPhi when the difference between the desired phi coordinate and the position of the global correction is deltaPhi, approach the desired phi coordinate by deltaPhi * \p approachPhi.
  /// \param diffCorr if the absolute differences from the interpolated values for the global corrections from the last iteration compared to the current iteration is smaller than this value, set converged to true for current global distortion
  void calcGlobalDistWithGlobalCorrIterative(const DistCorrInterpolator<DataT, Nz, Nr, Nphi>& globCorr, const int maxIter = 100, const DataT approachZ = 0.5, const DataT approachR = 0.5, const DataT approachPhi = 0.5, const DataT diffCorr = 1e-6);

  /// set the density, potential, electric fields, local distortions/corrections, global distortions/corrections from a file. Missing objects in the file are ignored.
  /// \file file containing the stored values for the density, potential, electric fields, local distortions/corrections, global distortions/corrections
  /// \param side side of the TPC
  void setFromFile(TFile& file, const Side side);

  /// Get grid spacing in r direction
  static constexpr DataT getGridSpacingR() { return GridProp::GRIDSPACINGR; }

  /// Get grid spacing in z direction
  static constexpr DataT getGridSpacingZ() { return GridProp::GRIDSPACINGZ; }

  /// Get grid spacing in phi direction
  static constexpr DataT getGridSpacingPhi() { return GridProp::GRIDSPACINGPHI; }

  /// Get constant electric field
  static constexpr DataT getEzField() { return (TPCParameters<DataT>::cathodev - TPCParameters<DataT>::vg1t) / TPCParameters<DataT>::TPCZ0; }

  /// Get inner radius of tpc
  static constexpr DataT getRMin() { return GridProp::RMIN; }

  /// Get min z position which is used during the calaculations
  static constexpr DataT getZMin() { return GridProp::ZMIN; }

  /// Get min phi
  static constexpr DataT getPhiMin() { return GridProp::PHIMIN; }

  /// Get max r
  static constexpr DataT getRMax() { return GridProp::RMAX; };

  /// Get max z
  static constexpr DataT getZMax() { return GridProp::ZMAX; }

  /// Get max phi
  static constexpr DataT getPhiMax() { return GridProp::PHIMAX; }

  // get side of TPC for z coordinate TODO rewrite this
  int getSide(const DataT z) const { return z <= 0 ? 0 : 1; }

  /// Get the grid object
  const RegularGrid& getGrid3D() const { return mGrid3D; }

  /// Get struct containing interpolators for the electrical fields
  /// \param side side of the TPC
  NumericalFields<DataT, Nz, Nr, Nphi> getElectricFieldsInterpolator(const Side side) const;

  /// Get struct containing interpolators for local distortions dR, dZ, dPhi
  /// \param side side of the TPC
  DistCorrInterpolator<DataT, Nz, Nr, Nphi> getLocalDistInterpolator(const Side side) const;

  /// Get struct containing interpolators for local corrections dR, dZ, dPhi
  /// \param side side of the TPC
  DistCorrInterpolator<DataT, Nz, Nr, Nphi> getLocalCorrInterpolator(const Side side) const;

  /// Get struct containing interpolators for global distortions dR, dZ, dPhi
  /// \param side side of the TPC
  DistCorrInterpolator<DataT, Nz, Nr, Nphi> getGlobalDistInterpolator(const Side side) const;

  /// Get struct containing interpolators for global corrections dR, dZ, dPhi
  /// \param side side of the TPC
  DistCorrInterpolator<DataT, Nz, Nr, Nphi> getGlobalCorrInterpolator(const Side side) const;

  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  /// \return returns local distortion dR for given vertex
  DataT getLocalDistR(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mLocalDistdR[side](iz, ir, iphi); }

  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  /// \return returns local distortion dZ for given vertex
  DataT getLocalDistZ(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mLocalDistdZ[side](iz, ir, iphi); }

  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  /// \return returns local distortion dRPhi for given vertex
  DataT getLocalDistRPhi(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mLocalDistdRPhi[side](iz, ir, iphi); }

  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  /// \return returns local correction dR for given vertex
  DataT getLocalCorrR(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mLocalCorrdR[side](iz, ir, iphi); }

  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  /// \return returns local correction dZ for given vertex
  DataT getLocalCorrZ(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mLocalCorrdZ[side](iz, ir, iphi); }

  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  /// \return returns local correction dRPhi for given vertex
  DataT getLocalCorrRPhi(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mLocalCorrdRPhi[side](iz, ir, iphi); }

  /// Get global distortion dR for vertex
  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  DataT getGlobalDistR(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mGlobalDistdR[side](iz, ir, iphi); }

  /// Get global distortion dZ for vertex
  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  DataT getGlobalDistZ(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mGlobalDistdZ[side](iz, ir, iphi); }

  /// Get global distortion dRPhi for vertex
  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  DataT getGlobalDistRPhi(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mGlobalDistdRPhi[side](iz, ir, iphi); }

  /// Get global correction dR for vertex
  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  DataT getGlobalCorrR(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mGlobalCorrdR[side](iz, ir, iphi); }

  /// Get global correction dZ for vertex
  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  DataT getGlobalCorrZ(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mGlobalCorrdZ[side](iz, ir, iphi); }

  /// Get global correction dRPhi for vertex
  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  DataT getGlobalCorrRPhi(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mGlobalCorrdRPhi[side](iz, ir, iphi); }

  /// Get global electric Field Er for vertex
  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  DataT getEr(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mElectricFieldEr[side](iz, ir, iphi); }

  /// Get global electric Field Ez for vertex
  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  DataT getEz(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mElectricFieldEz[side](iz, ir, iphi); }

  /// Get global electric Field Ephi for vertex
  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  DataT getEphi(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mElectricFieldEphi[side](iz, ir, iphi); }

  /// Get density for vertex
  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  DataT getDensity(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mDensity[side](iz, ir, iphi); }

  /// Get potential for vertex
  /// \param vertex in iz dimension
  /// \param vertex in ir dimension
  /// \param vertex in iphi dimension
  /// \param side side of the TPC
  DataT getPotential(const size_t iz, const size_t ir, const size_t iphi, const Side side) const { return mPotential[side](iz, ir, iphi); }

  /// Get the step width which is used for the calculation of the correction/distortions in units of the z-bin
  int getStepWidth() const { return 1 / mSteps; }

  /// Get phi vertex position for index in phi direction
  /// \param indexPhi index in phi direction
  DataT getPhiVertex(const size_t indexPhi) const { return mGrid3D.getZVertex(indexPhi); }

  /// Get r vertex position for index in r direction
  /// \param indexR index in r direction
  DataT getRVertex(const size_t indexR) const { return mGrid3D.getYVertex(indexR); }

  /// Get z vertex position for index in z direction
  /// \param indexZ index in z direction
  DataT getZVertex(const size_t indexZ) const { return mGrid3D.getXVertex(indexZ); }

  /// \param omegaTau \omega \tau value
  /// \param t1 value for t1 see: ???
  /// \param t2 value for t2 see: ???
  void setOmegaTauT1T2(const DataT omegaTau, const DataT t1, const DataT t2)
  {
    const DataT wt0 = t2 * omegaTau;
    mC0 = 1 / (1 + wt0 * wt0);
    const DataT wt1 = t1 * omegaTau;
    mC1 = wt1 / (1 + wt1 * wt1);
  };

  /// \param c0 coefficient C0 (compare Jim Thomas's notes for definitions)
  /// \param c1 coefficient C1 (compare Jim Thomas's notes for definitions)
  void setC0C1(const DataT c0, const DataT c1)
  {
    mC0 = c0;
    mC1 = c1;
  }

  /// set number of steps used for calculation of distortions/corrections per z bin
  /// \param nSteps number of steps per z bin
  void setNStep(const int nSteps) { mSteps = nSteps; }

  /// set which kind of numerical integration is used for calcution of the integrals int Er/Ez dz, int Ephi/Ez dz, int Ez dz
  /// \param strategy numerical integration strategy. see enum IntegrationStrategy for the different types
  void setNumericalIntegrationStrategy(const IntegrationStrategy strategy) { mNumericalIntegrationStrategy = strategy; }

  void setSimpsonNIteratives(const int nIter) { mSimpsonNIteratives = nIter; }

  /// write electric fields to root file
  /// \param outf output file where the electrical fields will be written to
  /// \side side of the TPC
  int dumpElectricFields(TFile& outf, const Side side) const;

  /// set electric field from root file
  /// \param inpf input file where the electrical fields are stored
  /// \side side of the TPC
  void setElectricFieldsFromFile(TFile& inpf, const Side side);

  /// write potential to root file
  /// \param outf output file where the potential will be written to
  /// \side side of the TPC
  int dumpPotential(TFile& outf, const Side side) const { return mPotential[side].writeToFile(outf, Form("potential_side%s", getSideName(side).data())); }

  /// set potential from root file
  /// \param inpf input file where the potential is stored
  /// \side side of the TPC
  void setPotentialFromFile(TFile& inpf, const Side side) { mPotential[side].initFromFile(inpf, Form("potential_side%s", getSideName(side).data())); }

  /// write potential to root file
  /// \param outf output file where the charge density will be written to
  /// \side side of the TPC
  int dumpDensity(TFile& outf, const Side side) const { return mDensity[side].writeToFile(outf, Form("density_side%s", getSideName(side).data())); }

  /// set potential from root file
  /// \param inpf input file where the charge density is stored
  /// \side side of the TPC
  void setDensityFromFile(TFile& inpf, const Side side) { mDensity[side].initFromFile(inpf, Form("density_side%s", getSideName(side).data())); }

  /// write global distortions to root file
  /// \param outf output file where the global distortions will be written to
  /// \side side of the TPC
  int dumpGlobalDistortions(TFile& outf, const Side side) const;

  /// set global distortions from root file
  /// \param inpf input file where the global distortions are stored
  /// \side side of the TPC
  void setGlobalDistortionsFromFile(TFile& inpf, const Side side);

  /// write global corrections to root file
  /// \param outf output file where the global corrections will be written to
  /// \side side of the TPC
  int dumpGlobalCorrections(TFile& outf, const Side side) const;

  /// set global corrections from root file
  /// \param inpf input file where the global corrections are stored
  /// \side side of the TPC
  void setGlobalCorrectionsFromFile(TFile& inpf, const Side side);

  /// write local corrections to root file
  /// \param outf output file where the local corrections will be written to
  /// \side side of the TPC
  int dumpLocalCorrections(TFile& outf, const Side side) const;

  /// set local corrections from root file
  /// \param inpf input file where the local corrections are stored
  /// \side side of the TPC
  void setLocalCorrectionsFromFile(TFile& inpf, const Side side);

  /// write local distortions to root file
  /// \param outf output file where the local distortions will be written to
  /// \side side of the TPC
  int dumpLocalDistortions(TFile& outf, const Side side) const;

  /// set local distortions from root file
  /// \param inpf input file where the local distortions are stored
  /// \side side of the TPC
  void setLocalDistortionsFromFile(TFile& inpf, const Side side);

  /// set z coordinate between min z max z
  /// \param posZ z position which will be regulated if needed
  DataT regulateZ(const DataT posZ) const { return mGrid3D.clampToGrid(posZ, 0); }

  /// set r coordinate between 'RMIN - 4 * GRIDSPACINGR' and 'RMAX + 2 * GRIDSPACINGR'. the r coordinate is not clamped to RMIN and RMAX to ensure correct interpolation at the borders of the grid.
  DataT regulateR(const DataT posR) const;

  /// set phi coordinate between min phi max phi
  DataT regulatePhi(const DataT posPhi) const { return mGrid3D.clampToGridCircular(posPhi, 2); }

 private:
  using ASolv = o2::tpc::PoissonSolver<DataT, Nz, Nr, Nphi>;

  IntegrationStrategy mNumericalIntegrationStrategy = IntegrationStrategy::SimpsonIterative; ///< numerical integration strategy of integration of the E-Field: 0: trapezoidal, 1: Simpson, 2: Root (only for analytical formula case)
  int mSimpsonNIteratives = 3;                                                               ///< number of iterations which are performed in the iterative simpson calculation of distortions/corrections
  int mSteps = 1;                                                                            ///< during the calculation of the corrections/distortions it is assumed that the electron drifts on a line from deltaZ = z0 -> z1. The value sets the deltaZ width: 1: deltaZ=zBin/1, 5: deltaZ=zBin/5

  DataT mC0 = 0; ///< coefficient C0 (compare Jim Thomas's notes for definitions)
  DataT mC1 = 0; ///< coefficient C1 (compare Jim Thomas's notes for definitions)

  static constexpr int FNSIDES = SIDES; ///< number of sides of the TPC
  bool mIsEfieldSet[FNSIDES]{};         ///< flag if E-fields are set
  bool mIsLocalCorrSet[FNSIDES]{};      ///< flag if local corrections are set
  bool mIsLocalDistSet[FNSIDES]{};      ///< flag if local distortions are set
  bool mIsGlobalCorrSet[FNSIDES]{};     ///< flag if global corrections are set
  bool mIsGlobalDistSet[FNSIDES]{};     ///< flag if global distortions are set

  RegularGrid mGrid3D{getZMin(), getRMin(), getPhiMin(), getGridSpacingZ(), getGridSpacingR(), getGridSpacingPhi()}; ///< grid properties

  DataContainer mLocalDistdR[FNSIDES]{};    ///< data storage for local distortions dR
  DataContainer mLocalDistdZ[FNSIDES]{};    ///< data storage for local distortions dZ
  DataContainer mLocalDistdRPhi[FNSIDES]{}; ///< data storage for local distortions dRPhi

  DataContainer mLocalCorrdR[FNSIDES]{};    ///< data storage for local corrections dR
  DataContainer mLocalCorrdZ[FNSIDES]{};    ///< data storage for local corrections dZ
  DataContainer mLocalCorrdRPhi[FNSIDES]{}; ///< data storage for local corrections dRPhi

  DataContainer mGlobalDistdR[FNSIDES]{};    ///< data storage for global distortions dR
  DataContainer mGlobalDistdZ[FNSIDES]{};    ///< data storage for global distortions dZ
  DataContainer mGlobalDistdRPhi[FNSIDES]{}; ///< data storage for global distortions dRPhi

  DataContainer mGlobalCorrdR[FNSIDES]{};    ///< data storage for global corrections dR
  DataContainer mGlobalCorrdZ[FNSIDES]{};    ///< data storage for global corrections dZ
  DataContainer mGlobalCorrdRPhi[FNSIDES]{}; ///< data storage for global corrections dRPhi

  DataContainer mDensity[FNSIDES]{};   ///< data storage for space charge density
  DataContainer mPotential[FNSIDES]{}; ///< data storage for the potential

  DataContainer mElectricFieldEr[FNSIDES]{};   ///< data storage for the electric field Er
  DataContainer mElectricFieldEz[FNSIDES]{};   ///< data storage for the electric field Ez
  DataContainer mElectricFieldEphi[FNSIDES]{}; ///< data storage for the electric field Ephi

  /// rebin the input space charge density histogram to desired binning
  /// \param hOrig original histogram
  TH3D rebinDensityHisto(const TH3& hOrig) const;

  /// get inverse spacing in z direction
  DataT getInvSpacingZ() const { return mGrid3D.getInvSpacingX(); }

  /// get inverse spacing in r direction
  DataT getInvSpacingR() const { return mGrid3D.getInvSpacingY(); }

  /// get inverse spacing in phi direction
  DataT getInvSpacingPhi() const { return mGrid3D.getInvSpacingZ(); }

  std::string getSideName(const Side side) const { return side == Side::A ? "A" : "C"; }

  /// calculate distortions or corrections analytical with electric fields
  template <typename Fields = AnalyticalFields<DataT>>
  void calcDistCorr(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& ddR, DataT& ddPhi, DataT& ddZ, const Fields& formulaStruct, const bool localDistCorr) const;

  /// calculate distortions/corrections using the formulas proposed in https://edms.cern.ch/ui/file/1108138/1/ALICE-INT-2010-016.pdf page 7
  void langevinCylindrical(DataT& ddR, DataT& ddPhi, DataT& ddZ, const DataT radius, const DataT localIntErOverEz, const DataT localIntEPhiOverEz, const DataT localIntDeltaEz) const;

  /// integrate electrical fields using root integration method
  template <typename Fields = AnalyticalFields<DataT>>
  void integrateEFieldsRoot(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, const Fields& formulaStruct) const;

  /// integrate electrical fields using trapezoidal integration method
  template <typename Fields = AnalyticalFields<DataT>>
  void integrateEFieldsTrapezoidal(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, const Fields& formulaStruct) const;

  /// integrate electrical fields using simpson integration method
  template <typename Fields = AnalyticalFields<DataT>>
  void integrateEFieldsSimpson(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, const Fields& formulaStruct) const;

  /// integrate electrical fields using simpson integration method with non straight drift of electrons
  template <typename Fields = AnalyticalFields<DataT>>
  void integrateEFieldsSimpsonIterative(const DataT p1r, const DataT p2r, const DataT p1phi, const DataT p2phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, const Fields& formulaStruct) const;

  /// calculate distortions/corrections using analytical electric fields
  void processGlobalDistCorr(const DataT radius, const DataT phi, const DataT z0Tmp, const DataT z1Tmp, DataT& ddR, DataT& ddPhi, DataT& ddZ, const AnalyticalFields<DataT>& formulaStruct) const
  {
    calcDistCorr(radius, phi, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct, false);
  }

  /// calculate distortions/corrections using electric fields from tricubic interpolator
  void processGlobalDistCorr(const DataT radius, const DataT phi, const DataT z0Tmp, const DataT z1Tmp, DataT& ddR, DataT& ddPhi, DataT& ddZ, const NumericalFields<DataT, Nz, Nr, Nphi>& formulaStruct) const
  {
    calcDistCorr(radius, phi, z0Tmp, z1Tmp, ddR, ddPhi, ddZ, formulaStruct, false);
  }

  /// calculate distortions/corrections by interpolation of local distortions/corrections
  void processGlobalDistCorr(const DataT radius, const DataT phi, const DataT z0Tmp, [[maybe_unused]] const DataT z1Tmp, DataT& ddR, DataT& ddPhi, DataT& ddZ, const DistCorrInterpolator<DataT, Nz, Nr, Nphi>& localDistCorr) const
  {
    ddR = localDistCorr.evaldR(z0Tmp, radius, phi);
    ddZ = localDistCorr.evaldZ(z0Tmp, radius, phi);
    ddPhi = localDistCorr.evaldRPhi(z0Tmp, radius, phi) / radius;
  }
};

///
/// ========================================================================================================
///                                Inline implementations of some methods
/// ========================================================================================================
///

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename Fields>
void SpaceCharge<DataT, Nz, Nr, Nphi>::integrateEFieldsRoot(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, const Fields& formulaStruct) const
{
  const DataT ezField = getEzField();
  TF1 fErOverEz("fErOverEz", [&](double* x, double* p) { (void)p; return static_cast<double>(formulaStruct.evalEr(static_cast<DataT>(x[0]), p1r, p1phi) / (formulaStruct.evalEz(static_cast<DataT>(x[0]), p1r, p1phi) + ezField)); }, p1z, p2z, 1);
  localIntErOverEz = static_cast<DataT>(fErOverEz.Integral(p1z, p2z));

  TF1 fEphiOverEz("fEPhiOverEz", [&](double* x, double* p) { (void)p; return static_cast<double>(formulaStruct.evalEphi(static_cast<DataT>(x[0]), p1r, p1phi) / (formulaStruct.evalEz(static_cast<DataT>(x[0]), p1r, p1phi) + ezField)); }, p1z, p2z, 1);
  localIntEPhiOverEz = static_cast<DataT>(fEphiOverEz.Integral(p1z, p2z));

  TF1 fEz("fEZOverEz", [&](double* x, double* p) { (void)p; return static_cast<double>(formulaStruct.evalEz(static_cast<DataT>(x[0]), p1r, p1phi) - ezField); }, p1z, p2z, 1);
  localIntDeltaEz = static_cast<DataT>(fEz.Integral(p1z, p2z));
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename Fields>
void SpaceCharge<DataT, Nz, Nr, Nphi>::integrateEFieldsTrapezoidal(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, const Fields& formulaStruct) const
{
  //========trapezoidal rule see: https://en.wikipedia.org/wiki/Trapezoidal_rule ==============
  const DataT fielder0 = formulaStruct.evalEr(p1z, p1r, p1phi);
  const DataT fieldez0 = formulaStruct.evalEz(p1z, p1r, p1phi);
  const DataT fieldephi0 = formulaStruct.evalEphi(p1z, p1r, p1phi);

  const DataT fielder1 = formulaStruct.evalEr(p2z, p1r, p1phi);
  const DataT fieldez1 = formulaStruct.evalEz(p2z, p1r, p1phi);
  const DataT fieldephi1 = formulaStruct.evalEphi(p2z, p1r, p1phi);

  const DataT ezField = getEzField();
  const DataT eZ0 = 1. / (ezField + fieldez0);
  const DataT eZ1 = 1. / (ezField + fieldez1);

  const DataT deltaX = 0.5 * (p2z - p1z);
  localIntErOverEz = deltaX * (fielder0 * eZ0 + fielder1 * eZ1);
  localIntEPhiOverEz = deltaX * (fieldephi0 * eZ0 + fieldephi1 * eZ1);
  localIntDeltaEz = deltaX * (fieldez0 + fieldez1);
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename Fields>
void SpaceCharge<DataT, Nz, Nr, Nphi>::integrateEFieldsSimpson(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, const Fields& formulaStruct) const
{
  //==========simpsons rule see: https://en.wikipedia.org/wiki/Simpson%27s_rule =============================
  const DataT fielder0 = formulaStruct.evalEr(p1z, p1r, p1phi);
  const DataT fieldez0 = formulaStruct.evalEz(p1z, p1r, p1phi);
  const DataT fieldephi0 = formulaStruct.evalEphi(p1z, p1r, p1phi);

  const DataT fielder1 = formulaStruct.evalEr(p2z, p1r, p1phi);
  const DataT fieldez1 = formulaStruct.evalEz(p2z, p1r, p1phi);
  const DataT fieldephi1 = formulaStruct.evalEphi(p2z, p1r, p1phi);

  const DataT deltaX = p2z - p1z;
  const DataT ezField = getEzField();
  const DataT xk2N = (p2z - static_cast<DataT>(0.5) * deltaX);
  const DataT ezField2 = formulaStruct.evalEz(xk2N, p1r, p1phi);
  const DataT ezField2Denominator = 1. / (ezField + ezField2);
  const DataT fieldSum2ErOverEz = formulaStruct.evalEr(xk2N, p1r, p1phi) * ezField2Denominator;
  const DataT fieldSum2EphiOverEz = formulaStruct.evalEphi(xk2N, p1r, p1phi) * ezField2Denominator;

  const DataT eZ0 = 1. / (ezField + fieldez0);
  const DataT eZ1 = 1. / (ezField + fieldez1);

  const DataT deltaXSimpsonSixth = deltaX / 6.;
  localIntErOverEz = deltaXSimpsonSixth * (4. * fieldSum2ErOverEz + fielder0 * eZ0 + fielder1 * eZ1);
  localIntEPhiOverEz = deltaXSimpsonSixth * (4. * fieldSum2EphiOverEz + fieldephi0 * eZ0 + fieldephi1 * eZ1);
  localIntDeltaEz = deltaXSimpsonSixth * (4. * ezField2 + fieldez0 + fieldez1);
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename Fields>
void SpaceCharge<DataT, Nz, Nr, Nphi>::integrateEFieldsSimpsonIterative(const DataT p1r, const DataT p2r, const DataT p1phi, const DataT p2phi, const DataT p1z, const DataT p2z, DataT& localIntErOverEz, DataT& localIntEPhiOverEz, DataT& localIntDeltaEz, const Fields& formulaStruct) const
{
  //==========simpsons rule see: https://en.wikipedia.org/wiki/Simpson%27s_rule =============================
  const DataT ezField = getEzField();
  const DataT p2phiSave = regulatePhi(p2phi);

  const DataT fielder0 = formulaStruct.evalEr(p1z, p1r, p1phi);
  const DataT fieldez0 = formulaStruct.evalEz(p1z, p1r, p1phi);
  const DataT fieldephi0 = formulaStruct.evalEphi(p1z, p1r, p1phi);

  const DataT fielder1 = formulaStruct.evalEr(p2z, p2r, p2phiSave);
  const DataT fieldez1 = formulaStruct.evalEz(p2z, p2r, p2phiSave);
  const DataT fieldephi1 = formulaStruct.evalEphi(p2z, p2r, p2phiSave);

  const DataT eZ0Inv = 1. / (ezField + fieldez0);
  const DataT eZ1Inv = 1. / (ezField + fieldez1);

  const DataT pHalfZ = 0.5 * (p1z + p2z);                        // dont needs to be regulated since p1z and p2z are already regulated
  const DataT pHalfPhiSave = regulatePhi(0.5 * (p1phi + p2phi)); // needs to be regulated since p2phi is not regulated
  const DataT pHalfR = 0.5 * (p1r + p2r);

  const DataT ezField2 = formulaStruct.evalEz(pHalfZ, pHalfR, pHalfPhiSave);
  const DataT eZHalfInv = 1. / (ezField + ezField2);
  const DataT fieldSum2ErOverEz = formulaStruct.evalEr(pHalfZ, pHalfR, pHalfPhiSave);
  const DataT fieldSum2EphiOverEz = formulaStruct.evalEphi(pHalfZ, pHalfR, pHalfPhiSave);

  const DataT deltaXSimpsonSixth = (p2z - p1z) / 6;
  localIntErOverEz = deltaXSimpsonSixth * (4 * fieldSum2ErOverEz * eZHalfInv + fielder0 * eZ0Inv + fielder1 * eZ1Inv);
  localIntEPhiOverEz = deltaXSimpsonSixth * (4 * fieldSum2EphiOverEz * eZHalfInv + fieldephi0 * eZ0Inv + fieldephi1 * eZ1Inv);
  localIntDeltaEz = deltaXSimpsonSixth * (4 * ezField2 + fieldez0 + fieldez1);
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
template <typename Fields>
void SpaceCharge<DataT, Nz, Nr, Nphi>::calcDistCorr(const DataT p1r, const DataT p1phi, const DataT p1z, const DataT p2z, DataT& ddR, DataT& ddPhi, DataT& ddZ, const Fields& formulaStruct, const bool localDistCorr) const
{
  // see: https://edms.cern.ch/ui/file/1108138/1/ALICE-INT-2010-016.pdf
  // needed for calculation of distortions/corrections
  DataT localIntErOverEz = 0;   // integral_p1z^p2z Er/Ez dz
  DataT localIntEPhiOverEz = 0; // integral_p1z^p2z Ephi/Ez dz
  DataT localIntDeltaEz = 0;    // integral_p1z^p2z Ez dz

  // there are differentnumerical integration strategys implements. for details see each function.
  switch (mNumericalIntegrationStrategy) {
    case IntegrationStrategy::SimpsonIterative:                                // iterative simpson integration (should be more precise at least for the analytical E-Field case but takes alot more time than normal simpson integration)
      for (int i = 0; i < mSimpsonNIteratives; ++i) {                          // TODO define a convergence criterion to abort the algorithm earlier for speed up.
        const DataT tmpZ = localDistCorr ? (p2z + ddZ) : regulateZ(p2z + ddZ); // dont regulate for local distortions/corrections! (to get same result as using electric field at last/first bin)
        integrateEFieldsSimpsonIterative(p1r, p1r + ddR, p1phi, p1phi + ddPhi, p1z, tmpZ, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, formulaStruct);
        langevinCylindrical(ddR, ddPhi, ddZ, (p1r + 0.5 * ddR), localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz); // using the mean radius '(p1r + 0.5 * ddR)' for calculation of distortions/corections
      }
      break;
    case IntegrationStrategy::Simpson: // simpson integration
      integrateEFieldsSimpson(p1r, p1phi, p1z, p2z, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, formulaStruct);
      langevinCylindrical(ddR, ddPhi, ddZ, p1r, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz);
      break;
    case IntegrationStrategy::Trapezoidal: // trapezoidal integration (fastest)
      integrateEFieldsTrapezoidal(p1r, p1phi, p1z, p2z, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, formulaStruct);
      langevinCylindrical(ddR, ddPhi, ddZ, p1r, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz);
      break;
    case IntegrationStrategy::Root: // using integration implemented in ROOT (slow)
      integrateEFieldsRoot(p1r, p1phi, p1z, p2z, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, formulaStruct);
      langevinCylindrical(ddR, ddPhi, ddZ, p1r, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz);
      break;
    default:
      LOGP(INFO, "no matching case: Using Simpson");
      integrateEFieldsSimpson(p1r, p1phi, p1z, p2z, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, formulaStruct);
      langevinCylindrical(ddR, ddPhi, ddZ, p1r, localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz);
  }
}

template <typename DataT, size_t Nz, size_t Nr, size_t Nphi>
void SpaceCharge<DataT, Nz, Nr, Nphi>::langevinCylindrical(DataT& ddR, DataT& ddPhi, DataT& ddZ, const DataT radius, const DataT localIntErOverEz, const DataT localIntEPhiOverEz, const DataT localIntDeltaEz) const
{
  // calculated distortions/correction with the formula described in https://edms.cern.ch/ui/file/1108138/1/ALICE-INT-2010-016.pdf page 7.
  ddR = mC0 * localIntErOverEz + mC1 * localIntEPhiOverEz;
  ddPhi = (mC0 * localIntEPhiOverEz - mC1 * localIntErOverEz) / radius;
  ddZ = -localIntDeltaEz * TPCParameters<DataT>::DVDE;
}

} // namespace tpc
} // namespace o2

#endif
