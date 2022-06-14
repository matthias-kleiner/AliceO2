// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  TPCFastSpaceChargeCorrection.h
/// \brief Definition of TPCFastSpaceChargeCorrection class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_TPCSLOWSPACECHARGECORRECTION_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_TPCSLOWSPACECHARGECORRECTION_H

#include "FlatObject.h"
#include "GPUCommonDef.h"
#include "TriCubic.h"
#include <iostream>

namespace GPUCA_NAMESPACE
{
namespace gpu
{
using namespace o2::tpc;

class TPCSlowSpaceChargeCorrection
{
 public:
  using DataT = float;

  // TPCSlowSpaceChargeCorrection() = default;

  /// convert x and y coordinates from cartesian to the radius in polar coordinates
  DataT getRadiusFromCartesian(const DataT x, const DataT y) const { return std::sqrt(x * x + y * y); }

  /// convert x and y coordinates from cartesian to phi in polar coordinates
  DataT getPhiFromCartesian(const DataT x, const DataT y) const { return std::atan2(y, x); }

  /// convert radius and phi coordinates from polar coordinates to x cartesian coordinates
  static DataT getXFromPolar(const DataT r, const DataT phi) { return r * std::cos(phi); }

  /// convert radius and phi coordinates from polar coordinates to y cartesian coordinate
  static DataT getYFromPolar(const DataT r, const DataT phi) { return r * std::sin(phi); }

  void getCorrections(const DataT x, const DataT y, const DataT z, const bool side, DataT& corrX, DataT& corrY, DataT& corrZ) const
  {
    // convert cartesian to polar
    const DataT radius = getRadiusFromCartesian(x, y);
    const DataT phi = getPhiFromCartesian(x, y);

    DataT corrR = evaldR(side, z, radius, phi);
    DataT corrRPhi = evaldRPhi(side, z, radius, phi);
    corrZ = evaldZ(side, z, radius, phi);

    // Calculate corrected position
    const DataT radiusCorr = radius + corrR;
    const DataT phiCorr = phi + corrRPhi / radius;

    corrX = getXFromPolar(radiusCorr, phiCorr) - x; // difference between corrected and original x coordinate
    corrY = getYFromPolar(radiusCorr, phiCorr) - y; // difference between corrected and original y coordinate
  }

  template <typename DataTIn = DataT>
  void initFromFile(TFile& inpf)
  {
    std::cout << "Init from file: " << std::endl;
    mGlobalCorrdR[0].template initFromFile<DataTIn>(inpf, "corrR_sideA");
    mGlobalCorrdZ[0].template initFromFile<DataTIn>(inpf, "corrZ_sideA");
    mGlobalCorrdRPhi[0].template initFromFile<DataTIn>(inpf, "corrRPhi_sideA");

    mGlobalCorrdR[1].template initFromFile<DataTIn>(inpf, "corrR_sideC");
    mGlobalCorrdZ[1].template initFromFile<DataTIn>(inpf, "corrZ_sideC");
    mGlobalCorrdRPhi[1].template initFromFile<DataTIn>(inpf, "corrRPhi_sideC");
  }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the local distortion or correction dR for given coordinate
  DataT evaldR(bool SideA, const DataT z, const DataT r, const DataT phi) const { return (SideA==0) ? interpolatorDistCorrdRA(z, r, phi, mInterpolType) : interpolatorDistCorrdRC(z, r, phi, mInterpolType); }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the local distortion or correction dZ for given coordinate
  DataT evaldZ(bool SideA, const DataT z, const DataT r, const DataT phi) const { return (SideA==0) ? interpolatorDistCorrdZA(z, r, phi, mInterpolType) : interpolatorDistCorrdZC(z, r, phi, mInterpolType); }

  /// \param r r coordinate
  /// \param phi phi coordinate
  /// \param z z coordinate
  /// \return returns the function value for the local distortion or correction dRPhi for given coordinate
  DataT evaldRPhi(bool SideA, const DataT z, const DataT r, const DataT phi) const { return (SideA==0) ? interpolatorDistCorrdRPhiA(z, r, phi, mInterpolType) : interpolatorDistCorrdRPhiC(z, r, phi, mInterpolType); }

 private:
  static constexpr int NSIDES = 2;
  using RegularGrid = o2::tpc::RegularGrid3DGPU<DataT>;
  using DataContainer = o2::tpc::DataContainer3DGPU<DataT>;
  using TriCubic = o2::tpc::TriCubicInterpolatorGPU<DataT>;
  using GridProp = GridPropertiesGPU<DataT>;

  ParameterSpaceChargeGPU mParams;

  DataContainer mGlobalCorrdR[NSIDES]{DataContainer(mParams.NZVertices, mParams.NRVertices, mParams.NPhiVertices), DataContainer(mParams.NZVertices, mParams.NRVertices, mParams.NPhiVertices)};    ///< data storage for local corrections dR
  DataContainer mGlobalCorrdZ[NSIDES]{DataContainer(mParams.NZVertices, mParams.NRVertices, mParams.NPhiVertices), DataContainer(mParams.NZVertices, mParams.NRVertices, mParams.NPhiVertices)};    ///< data storage for local corrections dZ
  DataContainer mGlobalCorrdRPhi[NSIDES]{DataContainer(mParams.NZVertices, mParams.NRVertices, mParams.NPhiVertices), DataContainer(mParams.NZVertices, mParams.NRVertices, mParams.NPhiVertices)}; ///< data storage for local corrections dRPhi

  const RegularGrid mParamGrid[NSIDES]{
    {GridProp::ZMIN, GridProp::RMIN, GridProp::PHIMIN, GridProp::getGridSpacingZ(mParams.NZVertices), GridProp::getGridSpacingR(mParams.NRVertices), GridProp::getGridSpacingPhi(mParams.NPhiVertices)},
    {GridProp::ZMIN, GridProp::RMIN, GridProp::PHIMIN, -GridProp::getGridSpacingZ(mParams.NZVertices), GridProp::getGridSpacingR(mParams.NRVertices), GridProp::getGridSpacingPhi(mParams.NPhiVertices)}}; ///< grid properties

  TriCubic interpolatorDistCorrdRA{mGlobalCorrdR[0], mParamGrid[0]};       ///< TriCubic interpolator of distortion or correction dR
  TriCubic interpolatorDistCorrdZA{mGlobalCorrdZ[0], mParamGrid[0]};       ///< TriCubic interpolator of distortion or correction dZ
  TriCubic interpolatorDistCorrdRPhiA{mGlobalCorrdRPhi[0], mParamGrid[0]}; ///< TriCubic interpolator of distortion or correction dRPhi

  TriCubic interpolatorDistCorrdRC{mGlobalCorrdR[1], mParamGrid[1]};                        ///< TriCubic interpolator of distortion or correction dR
  TriCubic interpolatorDistCorrdZC{mGlobalCorrdZ[1], mParamGrid[1]};                        ///< TriCubic interpolator of distortion or correction dZ
  TriCubic interpolatorDistCorrdRPhiC{mGlobalCorrdRPhi[1], mParamGrid[1]};                  ///< TriCubic interpolator of distortion or correction dRPhi
  typename TriCubic::InterpolationType mInterpolType = TriCubic::InterpolationType::Sparse; ///< type of TriCubic interpolation
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
