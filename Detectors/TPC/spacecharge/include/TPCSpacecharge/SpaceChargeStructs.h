// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file SpaceChargeStructs.h
/// \brief This file provides all necesseray structs which are used during the calcution of the distortions and corrections
///
/// \author  Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>
/// \date Aug 21, 2020

#ifndef ALICEO2_TPC_SpaceChargeStructs_H_
#define ALICEO2_TPC_SpaceChargeStructs_H_

#include <functional>
#include <cmath>
#include "TPCSpacecharge/TriCubic.h"
#include "DataFormatsTPC/Defs.h"
#include "TPCSpacecharge/SpaceChargeStructs.h"

namespace o2
{
namespace tpc
{

///
/// this struct contains an analytical description of the space charge, potential and the electric fields.
/// The analytical functions can be used to test the poisson solver and the caluclation of distortions/corrections.
///
template <typename DataT = double>
struct AnalyticalFields {

  DataT parA{1e-5};                     ///< parameter [0] of functions
  DataT parB{0.5};                      ///< parameter [1] of functions
  DataT parC{1e-4};                     ///< parameter [2] of functions
  o2::tpc::Side side{o2::tpc::Side::A}; ///< side of the TPC. Since the absolute value is taken during the calculations the choice of the side is arbitrary.

  AnalyticalFields(const o2::tpc::Side sideTmp = o2::tpc::Side::A) : side{sideTmp} {};

  o2::tpc::Side getSide() const { return side; }

  void setSide(const o2::tpc::Side sideTmp) { side = sideTmp; }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Er for given coordinate
  DataT evalEr(DataT z, DataT r, DataT phi) const
  {
    return erFunc(z, r, phi);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Ez for given coordinate
  DataT evalEz(DataT z, DataT r, DataT phi) const
  {
    return ezFunc(z, r, phi);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Ephi for given coordinate
  DataT evalEphi(DataT z, DataT r, DataT phi) const
  {
    return ephiFunc(z, r, phi);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the potential for given coordinate
  DataT evalPotential(DataT z, DataT r, DataT phi) const
  {
    return potentialFunc(z, r, phi);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the space charge density for given coordinate
  DataT evalDensity(DataT z, DataT r, DataT phi) const
  {
    return densityFunc(z, r, phi);
  }

  /// analytical potential
  std::function<DataT(DataT, DataT, DataT)> potentialFunc = [& parA = parA, &parB = parB, &parC = parC](DataT z, DataT r, DataT phi) {
    return -parA * (std::pow((-r + 254.5 + 83.5), 4) - 338.0 * std::pow((-r + 254.5 + 83.5), 3) + 21250.75 * std::pow((-r + 254.5 + 83.5), 2)) * std::cos(parB * phi) * std::cos(parB * phi) * std::exp(-1 * parC * (z - 125) * (z - 125));
  };

  /// analytical space charge - NOTE: if the space charge density is calculated analytical there would be a - sign in the formula (-parA)  - however since its an e- the sign is flipped (IS THIS CORRECT??? see for minus sign: AliTPCSpaceCharge3DCalc::SetPotentialBoundaryAndChargeFormula)-
  std::function<DataT(DataT, DataT, DataT)> densityFunc = [& parA = parA, &parB = parB, &parC = parC](DataT z, DataT r, DataT phi) {
    return parA * ((1 / r * 16 * (-3311250 + 90995.5 * r - 570.375 * r * r + r * r * r)) * std::cos(parB * phi) * std::cos(parB * phi) * std::exp(-1 * parC * (z - 125) * (z - 125)) +
                   (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * std::pow(-r + 254.5 + 83.5, 2)) / (r * r) * std::exp(-1 * parC * (z - 125) * (z - 125)) * -2 * parB * parB * std::cos(2 * parB * phi) +
                   (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * std::pow(-r + 254.5 + 83.5, 2)) * std::cos(parB * phi) * std::cos(parB * phi) * 2 * parC * std::exp(-1 * parC * (z - 125) * (z - 125)) * (2 * parC * (z - 125) * (z - 125) - 1));
  };

  /// analytical electric field Er
  std::function<DataT(DataT, DataT, DataT)> erFunc = [& parA = parA, &parB = parB, &parC = parC](DataT z, DataT r, DataT phi) {
    return parA * 4 * (r * r * r - 760.5 * r * r + 181991 * r - 1.3245 * std::pow(10, 7)) * std::cos(parB * phi) * std::cos(parB * phi) * std::exp(-1 * parC * (z - 125) * (z - 125));
  };

  /// analytical electric field Ephi
  std::function<DataT(DataT, DataT, DataT)> ephiFunc = [& parA = parA, &parB = parB, &parC = parC](DataT z, DataT r, DataT phi) {
    return parA * (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * (-r + 254.5 + 83.5) * (-r + 254.5 + 83.5)) / r * std::exp(-1 * parC * (z - 125) * (z - 125)) * -parB * std::sin(2 * parB * phi);
  };

  /// analytical electric field Ez
  std::function<DataT(DataT, DataT, DataT)> ezFunc = [& parA = parA, &parB = parB, &parC = parC](DataT z, DataT r, DataT phi) {
    return parA * (std::pow(-r + 254.5 + 83.5, 4) - 338.0 * std::pow(-r + 254.5 + 83.5, 3) + 21250.75 * (-r + 254.5 + 83.5) * (-r + 254.5 + 83.5)) * std::cos(parB * phi) * std::cos(parB * phi) * -2 * parC * (z - 125) * std::exp(-1 * parC * (z - 125) * (z - 125));
  };

  static constexpr unsigned int getID() { return ID; }

  // private:
  static constexpr unsigned int ID = 0; ///< needed to distinguish between the differrent structs
};

///
/// This struct gives tricubic interpolation of the electric fields and can be used to calculate the distortions/corrections.
/// The electric fields have to be calculated by the poisson solver or given by the analytical formula.
///
template <typename DataT = double, size_t Nr = 129, size_t Nz = 129, size_t Nphi = 180>
struct NumericalFields {
  using RegularGrid = o2::tpc::RegularGrid3D<DataT, Nz, Nr, Nphi>;
  using DataContainer = o2::tpc::DataContainer3D<DataT, Nz, Nr, Nphi>;
  using TriCubic = o2::tpc::TriCubicInterpolator<DataT, Nz, Nr, Nphi>;

  /// constructor
  /// \param dataEr container for the data of the electrical field Er
  /// \param dataEz container for the data of the electrical field Ez
  /// \param dataEphi container for the data of the electrical field Ephi
  /// \param gridProperties properties of the grid
  /// \param side side of the tpc
  NumericalFields(const DataContainer& dataErA, const DataContainer& dataEzA, const DataContainer& dataEphiA, const RegularGrid& gridPropertiesA, const o2::tpc::Side sideA) : dataEr{dataErA}, dataEz{dataEzA}, dataEphi{dataEphiA}, gridProperties{gridPropertiesA}, side{sideA} {};

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Er for given coordinate
  DataT evalEr(DataT z, DataT r, DataT phi) const
  {
    return interpolatorEr(z, r, phi, interpolType);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Ez for given coordinate
  DataT evalEz(DataT z, DataT r, DataT phi) const
  {
    return interpolatorEz(z, r, phi, interpolType);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for electric field Ephi for given coordinate
  DataT evalEphi(DataT z, DataT r, DataT phi) const
  {
    return interpolatorEphi(z, r, phi, interpolType);
  }

  o2::tpc::Side getSide() const
  {
    return side;
  }

  static constexpr unsigned int getID() { return ID; }

  /// set which kind of TriCubic interpolation algorithm is used
  void setInterpolationType(const int type)
  {
    interpolType = type;
  }

  // \return returns which kind of TriCubic interpolation algorithm is used
  int getInterpolationType() const
  {
    return interpolType;
  }

 private:
  const DataContainer& dataEr{};       ///< adress to the data container of the grid
  const DataContainer& dataEz{};       ///< adress to the data container of the grid
  const DataContainer& dataEphi{};     ///< adress to the data container of the grid
  const RegularGrid& gridProperties{}; ///< properties of the regular grid
  const o2::tpc::Side side{};          ///< side of the TPC

  TriCubic interpolatorEr{dataEr, gridProperties};                                         ///< TriCubic interpolator of the electric field Er
  TriCubic interpolatorEz{dataEz, gridProperties};                                         ///< TriCubic interpolator of the electric field Ez
  TriCubic interpolatorEphi{dataEphi, gridProperties};                                     ///< TriCubic interpolator of the electric field Ephi
  typename TriCubic::InterpolationType interpolType = TriCubic::InterpolationType::Sparse; ///< type of TriCubic interpolation
  static constexpr unsigned int ID = 1;                                                    ///< needed to distinguish between the different structs
};

///
/// This struct gives tricubic interpolation of the local distortions or corrections.
/// The the local distortions or corrections can be used to calculate the global distortions/corrections.
///
template <typename DataT = double, size_t Nr = 129, size_t Nz = 129, size_t Nphi = 180>
struct DistCorrInterpolator {
  using RegularGrid = o2::tpc::RegularGrid3D<DataT, Nz, Nr, Nphi>;
  using DataContainer = o2::tpc::DataContainer3D<DataT, Nz, Nr, Nphi>;
  using TriCubic = o2::tpc::TriCubicInterpolator<DataT, Nz, Nr, Nphi>;

  /// constructor
  /// \param dataDistCorrdR container for the data of the distortions dR
  /// \param dataDistCorrdZ container for the data of the distortions dZ
  /// \param dataDistCorrdRPhi container for the data of the distortions dPhi
  /// \param gridProperties properties of the grid
  /// \param side side of the tpc
  DistCorrInterpolator(const DataContainer& dataDistCorrdR, const DataContainer& dataDistCorrdZ, const DataContainer& dataDistCorrdRPhi, const RegularGrid& gridProperties, const o2::tpc::Side side) : dataDistCorrdR{dataDistCorrdR}, dataDistCorrdZ{dataDistCorrdZ}, dataDistCorrdRPhi{dataDistCorrdRPhi}, gridProperties{gridProperties}, side{side} {};

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the local distortion or correction dR for given coordinate
  DataT evaldR(const DataT z, const DataT r, const DataT phi) const
  {
    return interpolatorDistCorrdR(z, r, phi, interpolType);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the local distortion or correction dZ for given coordinate
  DataT evaldZ(const DataT z, const DataT r, const DataT phi) const
  {
    return interpolatorDistCorrdZ(z, r, phi, interpolType);
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the local distortion or correction dRPhi for given coordinate
  DataT evaldRPhi(const DataT z, const DataT r, const DataT phi) const
  {
    return interpolatorDistCorrdRPhi(z, r, phi, interpolType);
  }

  o2::tpc::Side getSide() const
  {
    return side;
  }

  static constexpr unsigned int getID() { return ID; }

  /// set which kind of TriCubic interpolation algorithm is used
  void setInterpolationType(const int type)
  {
    interpolType = type;
  }

  /// \return returns which kind of TriCubic interpolation algorithm is used
  int getInterpolationType() const
  {
    return interpolType;
  }

 private:
  const DataContainer& dataDistCorrdR{};    ///< adress to the data container of the grid
  const DataContainer& dataDistCorrdZ{};    ///< adress to the data container of the grid
  const DataContainer& dataDistCorrdRPhi{}; ///< adress to the data container of the grid
  const RegularGrid& gridProperties{};      ///< properties of the regular grid
  const o2::tpc::Side side{};               ///< side of the TPC.

  TriCubic interpolatorDistCorrdR{dataDistCorrdR, gridProperties};                         ///< TriCubic interpolator of distortion or correction dR
  TriCubic interpolatorDistCorrdZ{dataDistCorrdZ, gridProperties};                         ///< TriCubic interpolator of distortion or correction dZ
  TriCubic interpolatorDistCorrdRPhi{dataDistCorrdRPhi, gridProperties};                   ///< TriCubic interpolator of distortion or correction dRPhi
  typename TriCubic::InterpolationType interpolType = TriCubic::InterpolationType::Sparse; ///< type of TriCubic interpolation
  static constexpr unsigned int ID = 2;                                                    ///< needed to distinguish between the different structs
};

} // namespace tpc
} // namespace o2

#endif
