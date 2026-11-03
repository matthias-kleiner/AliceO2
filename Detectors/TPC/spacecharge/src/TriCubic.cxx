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

/// \file  TriCubic.cxx
/// \brief Definition of TriCubic class
///
/// \author  Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#include "TPCSpaceCharge/TriCubic.h"
#include "TPCSpaceCharge/DataContainer3D.h"
#include <algorithm>

using namespace o2::tpc;

/// default move constructor
template <typename DataT>
TriCubicInterpolator<DataT>::TriCubicInterpolator(TriCubicInterpolator<DataT>&& other)
{
  if (this != &other) {
    mExtrapolationType = other.mExtrapolationType;
    mExtraPolateValues = other.mExtraPolateValues;
    // do not change address of old pointers!
    other.mGridData = nullptr;
    other.mGridProperties = nullptr;
  }
}

template <typename DataT>
TriCubicInterpolator<DataT>& TriCubicInterpolator<DataT>::operator=(TriCubicInterpolator<DataT>&& other)
{
  if (this != &other) {
    mExtrapolationType = other.mExtrapolationType;
    mExtraPolateValues = other.mExtraPolateValues;
    // do not change address of old pointers!
    other.mGridData = nullptr;
    other.mGridProperties = nullptr;
  }
  return *this;
}

template <typename DataT>
DataT TriCubicInterpolator<DataT>::extrapolation(const DataT valk, const DataT valk1, const DataT valk2) const
{
  switch (mExtrapolationType) {
    case ExtrapolationType::Linear:
    default:
      return linearExtrapolation(valk, valk1);
      break;
    case ExtrapolationType::Parabola:
      return parabolExtrapolation(valk, valk1, valk2);
      break;
  }
}

template <typename DataT>
DataT TriCubicInterpolator<DataT>::linearExtrapolation(const DataT valk, const DataT valk1) const
{
  const DataT val = 2 * valk - valk1;
  return val;
}

template <typename DataT>
DataT TriCubicInterpolator<DataT>::parabolExtrapolation(const DataT valk, const DataT valk1, const DataT valk2) const
{
  const DataT val = 3 * (valk - valk1) + valk2; // legendre polynom with x0=0, x1=1, x2=2 and z=-1
  return val;
}

namespace
{
/// calculate the 4 weights of the 1D cubic (Catmull-Rom) interpolation for the relative position t in the cell
template <typename DataT>
inline void getWeights(const DataT t, DataT (&weights)[4])
{
  const DataT t2 = t * t;
  const DataT t3 = t2 * t;
  weights[0] = DataT(0.5) * (-t + 2 * t2 - t3);
  weights[1] = DataT(0.5) * (2 - 5 * t2 + 3 * t3);
  weights[2] = DataT(0.5) * (t + 4 * t2 - 3 * t3);
  weights[3] = DataT(0.5) * (t3 - t2);
}

/// sum_{slice, row, iz} weightsPhi[slice] * weightsR[row] * weightsZ[iz] * vals[sliceOffset[slice] + row * rowStride + iz]
/// \param vals values of the grid
/// \param sliceOffset offset of the 4 phi slices
/// \param rowStride offset between the rows in r direction (the 4 values in z direction are consecutive)
template <typename DataT>
inline DataT interpolateCell(const DataT* vals, const int (&sliceOffset)[4], const int rowStride, const DataT (&weightsZ)[4], const DataT (&weightsR)[4], const DataT (&weightsPhi)[4])
{
  DataT result = 0;
  for (int slice = 0; slice < 4; ++slice) {
    DataT sumR = 0;
    for (int row = 0; row < 4; ++row) {
      const DataT* val = vals + sliceOffset[slice] + row * rowStride;
      sumR += weightsR[row] * (weightsZ[0] * val[0] + weightsZ[1] * val[1] + weightsZ[2] * val[2] + weightsZ[3] * val[3]);
    }
    result += weightsPhi[slice] * sumR;
  }
  return result;
}
} // namespace

template <typename DataT>
TriCubicCell<DataT> TriCubicInterpolator<DataT>::getCell(const DataT z, const DataT r, const DataT phi) const
{
  TriCubicCell<DataT> cell;

  // relative position in the grid (in units of the grid spacing)
  DataT posRel[FDim]{(z - mGridProperties->getGridMinZ()) * mGridProperties->getInvSpacingZ(), (r - mGridProperties->getGridMinR()) * mGridProperties->getInvSpacingR(), (phi - mGridProperties->getGridMinPhi()) * mGridProperties->getInvSpacingPhi()};
  posRel[FPHI] = mGridProperties->clampToGridCircularRel(posRel[FPHI], FPHI);
  DataT posRelN[FDim]{posRel[FZ], posRel[FR], posRel[FPHI]};
  posRel[FZ] = mGridProperties->clampToGridRel(posRel[FZ], FZ);
  posRel[FR] = mGridProperties->clampToGridRel(posRel[FR], FR);
  if (!mExtraPolateValues) {
    // restrict the query position to the grid. The index of the cell is at most the second last vertex (see clampToGridRel()), but the position has to be clamped to the last vertex
    posRelN[FZ] = std::min(std::max(posRelN[FZ], DataT(0)), mGridProperties->getMaxIndexZ());
    posRelN[FR] = std::min(std::max(posRelN[FR], DataT(0)), mGridProperties->getMaxIndexR());
  }

  // index of the cell: the relative positions are >= 0 after the clamping, so the conversion to int is equal to floor()
  cell.indexZ = static_cast<int>(posRel[FZ]);
  cell.indexR = static_cast<int>(posRel[FR]);
  cell.indexPhi = static_cast<int>(posRel[FPHI]);

  // weights of the interpolation in each dimension
  getWeights<DataT>(posRelN[FZ] - cell.indexZ, cell.weightsZ);
  getWeights<DataT>(posRelN[FR] - cell.indexR, cell.weightsR);
  getWeights<DataT>(posRelN[FPHI] - cell.indexPhi, cell.weightsPhi);

  const GridPos location = findPos(cell.indexZ, cell.indexR, cell.indexPhi);
  cell.location = static_cast<int>(location);
  cell.inGrid = (location == GridPos::InnerVolume) || (location == GridPos::SideZRight) || (location == GridPos::SideZLeft);
  if (cell.inGrid) {
    // all 64 values are inside the grid (the phi direction is periodic): the values are directly taken from the grid.
    // The values in z direction are consecutive in memory (DataContainer3D::getDataIndex()): start at (iz - 1, ir - 1, iphi)
    int deltaPhi[3]{};
    getDataIndexCircularArray(cell.indexPhi, FPHI, deltaPhi);
    cell.sliceOffset[0] = deltaPhi[0];
    cell.sliceOffset[1] = 0;
    cell.sliceOffset[2] = deltaPhi[1];
    cell.sliceOffset[3] = deltaPhi[2];
    cell.rowStride = mGridProperties->getDeltaDataIndex(1, FR);
    cell.dataIndex = mGridData->getDataIndex(cell.indexZ, cell.indexR, cell.indexPhi) - cell.rowStride - 1;
  }
  return cell;
}

template <typename DataT>
DataT TriCubicInterpolator<DataT>::interpolateData(const DataContainer& data, const TriCubicCell<DataT>& cell) const
{
  // check if data is empty
  if (!data.getNDataPoints()) {
    return 0;
  }

  if (cell.inGrid) {
    return interpolateCell<DataT>(data.getData().data() + cell.dataIndex, cell.sliceOffset, cell.rowStride, cell.weightsZ, cell.weightsR, cell.weightsPhi);
  }

  // missing values outside of the grid are extrapolated
  DataT cVals[16][4];
  setValues(data, cell.indexZ, cell.indexR, cell.indexPhi, static_cast<GridPos>(cell.location), cVals);
  static constexpr int sliceOffset[4]{0, 16, 32, 48};
  return interpolateCell<DataT>(&cVals[0][0], sliceOffset, 4, cell.weightsZ, cell.weightsR, cell.weightsPhi);
}

template <typename DataT>
DataT TriCubicInterpolator<DataT>::interpolateSparse(const DataT z, const DataT r, const DataT phi) const
{
  // check if data is empty
  if (!mGridData->getNDataPoints()) {
    return 0;
  }
  return interpolateData(*mGridData, getCell(z, r, phi));
}

template <typename DataT>
void TriCubicInterpolator<DataT>::interpolate(const DataT z, const DataT r, const DataT phi, const TriCubicInterpolator<DataT>& other1, const TriCubicInterpolator<DataT>& other2, DataT& val0, DataT& val1, DataT& val2) const
{
  const size_t nPoints = mGridData->getNDataPoints();
  if (!nPoints || (other1.mGridData->getNDataPoints() != nPoints) || (other2.mGridData->getNDataPoints() != nPoints)) {
    // empty data or data with different grid: separate interpolation
    val0 = interpolateSparse(z, r, phi);
    val1 = other1.interpolateSparse(z, r, phi);
    val2 = other2.interpolateSparse(z, r, phi);
    return;
  }
  const TriCubicCell<DataT> cell = getCell(z, r, phi);
  val0 = interpolateData(*mGridData, cell);
  val1 = interpolateData(*other1.mGridData, cell);
  val2 = interpolateData(*other2.mGridData, cell);
}

// for perdiodic boundary condition
template <typename DataT>
void TriCubicInterpolator<DataT>::getDataIndexCircularArray(const int index0, const int dim, int arr[]) const
{
  const int delta_min1 = getRegulatedDelta(index0, -1, dim, mGridProperties->getN(dim) - 1);
  const int delta_plus1 = getRegulatedDelta(index0, +1, dim, 1 - mGridProperties->getN(dim));
  const int delta_plus2 = getRegulatedDelta(index0, +2, dim, 2 - mGridProperties->getN(dim));

  arr[0] = mGridProperties->getDeltaDataIndex(delta_min1, dim);
  arr[1] = mGridProperties->getDeltaDataIndex(delta_plus1, dim);
  arr[2] = mGridProperties->getDeltaDataIndex(delta_plus2, dim);
}

template <typename DataT>
typename TriCubicInterpolator<DataT>::GridPos TriCubicInterpolator<DataT>::findPos(const int iz, const int ir, const int iphi) const
{
  GridPos pos = GridPos::None;
  if (isInInnerVolume(iz, ir, iphi, pos)) {
    return pos;
  }

  if (findEdge(iz, ir, iphi, pos)) {
    return pos;
  }

  if (findLine(iz, ir, iphi, pos)) {
    return pos;
  }

  if (findSide(iz, ir, iphi, pos)) {
    return pos;
  }
  return GridPos::None;
}

template <typename DataT>
bool TriCubicInterpolator<DataT>::findEdge(const int iz, const int ir, const int iphi, GridPos& posType) const
{
  const int iR = 2;
  if (iz == 0 && ir == 0) {
    if (iphi == 0) {
      posType = GridPos::Edge0;
      return true;
    } else if (iphi == mGridData->getNPhi() - iR) {
      posType = GridPos::Edge4;
      return true;
    }
  } else if (iz == mGridData->getNZ() - iR && ir == 0) {
    if (iphi == 0) {
      posType = GridPos::Edge1;
      return true;
    } else if (iphi == mGridData->getNPhi() - iR) {
      posType = GridPos::Edge5;
      return true;
    }
  } else if (iz == 0 && ir == mGridData->getNR() - iR) {
    if (iphi == 0) {
      posType = GridPos::Edge2;
      return true;
    } else if (iphi == mGridData->getNPhi() - iR) {
      posType = GridPos::Edge6;
      return true;
    }
  } else if (iz == mGridData->getNZ() - iR && ir == mGridData->getNR() - iR) {
    if (iphi == 0) {
      posType = GridPos::Edge3;
      return true;
    } else if (iphi == mGridData->getNPhi() - iR) {
      posType = GridPos::Edge7;
      return true;
    }
  }
  return false;
}

template <typename DataT>
bool TriCubicInterpolator<DataT>::findLine(const int iz, const int ir, const int iphi, GridPos& posType) const
{
  const int iR = 2;
  // check line
  if (ir == 0) {
    if (iphi == 0) {
      posType = GridPos::LineA;
      return true;
    } else if (iphi == mGridData->getNPhi() - iR) {
      posType = GridPos::LineE;
      return true;
    }
    if (iz == 0) {
      posType = GridPos::LineI;
      return true;
    } else if (iz == mGridData->getNZ() - iR) {
      posType = GridPos::LineJ;
      return true;
    }
  } else if (ir == mGridData->getNR() - iR) {
    if (iphi == 0) {
      posType = GridPos::LineB;
      return true;
    } else if (iphi == mGridData->getNPhi() - iR) {
      posType = GridPos::LineF;
      return true;
    }
    if (iz == 0) {
      posType = GridPos::LineK;
      return true;
    } else if (iz == mGridData->getNZ() - iR) {
      posType = GridPos::LineL;
      return true;
    }
  } else if (iz == 0) {
    if (iphi == 0) {
      posType = GridPos::LineC;
      return true;
    } else if (iphi == mGridData->getNPhi() - iR) {
      posType = GridPos::LineG;
      return true;
    }
  } else if (iz == mGridData->getNZ() - iR) {
    if (iphi == 0) {
      posType = GridPos::LineD;
      return true;
    } else if (iphi == mGridData->getNPhi() - iR) {
      posType = GridPos::LineH;
      return true;
    }
  }
  return false;
}

template <typename DataT>
bool TriCubicInterpolator<DataT>::findSide(const int iz, const int ir, const int iphi, GridPos& posType) const
{
  if (isSideRight(iz, FZ)) {
    posType = GridPos::SideXRight;
    return true;
  } else if (isSideLeft(iz)) {
    posType = GridPos::SideXLeft;
    return true;
  }
  if (isSideRight(ir, FR)) {
    posType = GridPos::SideYRight;
    return true;
  } else if (isSideLeft(ir)) {
    posType = GridPos::SideYLeft;
    return true;
  }
  if (isSideRight(iphi, FPHI)) {
    posType = GridPos::SideZRight;
    return true;
  } else if (isSideLeft(iphi)) {
    posType = GridPos::SideZLeft;
    return true;
  }
  return false;
}

template <typename DataT>
bool TriCubicInterpolator<DataT>::isInInnerVolume(const int iz, const int ir, const int iphi, GridPos& posType) const
{
  if (iz >= 1 && iz < static_cast<int>(mGridData->getNZ() - 2) && ir >= 1 && ir < static_cast<int>(mGridData->getNR() - 2) && iphi >= 1 && iphi < static_cast<int>(mGridData->getNPhi() - 2)) {
    posType = GridPos::InnerVolume;
    return true;
  }
  return false;
}

template <typename DataT>
bool TriCubicInterpolator<DataT>::isSideRight(const int ind, const int dim) const
{
  if (ind == static_cast<int>(mGridProperties->getN(dim) - 2)) {
    return true;
  }
  return false;
}

template <typename DataT>
bool TriCubicInterpolator<DataT>::isSideLeft(const int ind) const
{
  if (ind == 0) {
    return true;
  }
  return false;
}

template <typename DataT>
void TriCubicInterpolator<DataT>::setValues(const DataContainer& data, const int iz, const int ir, const int iphi, const GridPos location, DataT (&cVals)[16][4]) const
{
  const int ii_x_y_z = data.getDataIndex(iz, ir, iphi);

  int deltaZ[3]{mGridProperties->getDeltaDataIndex(-1, 0), mGridProperties->getDeltaDataIndex(1, 0), mGridProperties->getDeltaDataIndex(2, 0)};
  int deltaR[3]{mGridProperties->getDeltaDataIndex(-1, 1), mGridProperties->getDeltaDataIndex(1, 1), mGridProperties->getDeltaDataIndex(2, 1)};
  int deltaPhi[3]{};
  getDataIndexCircularArray(iphi, FPHI, deltaPhi);

  const int i0 = 0;
  const int i1 = 1;
  const int i2 = 2;

  cVals[5][1] = data[ii_x_y_z];

  switch (location) {
    case GridPos::InnerVolume:
    case GridPos::SideZRight:
    case GridPos::SideZLeft:
    default: {
      const int ind[4][4][4]{
        {{ii_x_y_z + deltaPhi[i0] + deltaR[i0] + deltaZ[i0], ind[0][0][0] - deltaZ[i0], ind[0][0][1] - deltaZ[i0], ind[0][0][2] - deltaZ[i0]},
         {ind[0][0][0] - deltaR[i0], ind[0][1][0] - deltaZ[i0], ind[0][1][1] - deltaZ[i0], ind[0][1][2] - deltaZ[i0]},
         {ind[0][1][0] - deltaR[i0], ind[0][2][0] - deltaZ[i0], ind[0][2][1] - deltaZ[i0], ind[0][2][2] - deltaZ[i0]},
         {ind[0][2][0] - deltaR[i0], ind[0][3][0] - deltaZ[i0], ind[0][3][1] - deltaZ[i0], ind[0][3][2] - deltaZ[i0]}},
        {{ii_x_y_z + deltaR[i0] + deltaZ[i0], ind[1][0][0] - deltaZ[i0], ind[1][0][1] - deltaZ[i0], ind[1][0][2] - deltaZ[i0]},
         {ind[1][0][0] - deltaR[i0], ind[1][1][0] - deltaZ[i0], ind[1][1][1] - deltaZ[i0], ind[1][1][2] - deltaZ[i0]},
         {ind[1][1][0] - deltaR[i0], ind[1][2][0] - deltaZ[i0], ind[1][2][1] - deltaZ[i0], ind[1][2][2] - deltaZ[i0]},
         {ind[1][2][0] - deltaR[i0], ind[1][3][0] - deltaZ[i0], ind[1][3][1] - deltaZ[i0], ind[1][3][2] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i1] + deltaR[i0] + deltaZ[i0], ind[2][0][0] - deltaZ[i0], ind[2][0][1] - deltaZ[i0], ind[2][0][2] - deltaZ[i0]},
         {ind[2][0][0] - deltaR[i0], ind[2][1][0] - deltaZ[i0], ind[2][1][1] - deltaZ[i0], ind[2][1][2] - deltaZ[i0]},
         {ind[2][1][0] - deltaR[i0], ind[2][2][0] - deltaZ[i0], ind[2][2][1] - deltaZ[i0], ind[2][2][2] - deltaZ[i0]},
         {ind[2][2][0] - deltaR[i0], ind[2][3][0] - deltaZ[i0], ind[2][3][1] - deltaZ[i0], ind[2][3][2] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i2] + deltaR[i0] + deltaZ[i0], ind[3][0][0] - deltaZ[i0], ind[3][0][1] - deltaZ[i0], ind[3][0][2] - deltaZ[i0]},
         {ind[3][0][0] - deltaR[i0], ind[3][1][0] - deltaZ[i0], ind[3][1][1] - deltaZ[i0], ind[3][1][2] - deltaZ[i0]},
         {ind[3][1][0] - deltaR[i0], ind[3][2][0] - deltaZ[i0], ind[3][2][1] - deltaZ[i0], ind[3][2][2] - deltaZ[i0]},
         {ind[3][2][0] - deltaR[i0], ind[3][3][0] - deltaZ[i0], ind[3][3][1] - deltaZ[i0], ind[3][3][2] - deltaZ[i0]}}};

      cVals[0][0] = data[ind[0][0][0]];
      cVals[0][1] = data[ind[0][0][1]];
      cVals[0][2] = data[ind[0][0][2]];
      cVals[0][3] = data[ind[0][0][3]];
      cVals[1][0] = data[ind[0][1][0]];
      cVals[1][1] = data[ind[0][1][1]];
      cVals[1][2] = data[ind[0][1][2]];
      cVals[1][3] = data[ind[0][1][3]];
      cVals[2][0] = data[ind[0][2][0]];
      cVals[2][1] = data[ind[0][2][1]];
      cVals[2][2] = data[ind[0][2][2]];
      cVals[2][3] = data[ind[0][2][3]];
      cVals[3][0] = data[ind[0][3][0]];
      cVals[3][1] = data[ind[0][3][1]];
      cVals[3][2] = data[ind[0][3][2]];
      cVals[3][3] = data[ind[0][3][3]];
      cVals[4][0] = data[ind[1][0][0]];
      cVals[4][1] = data[ind[1][0][1]];
      cVals[4][2] = data[ind[1][0][2]];
      cVals[4][3] = data[ind[1][0][3]];
      cVals[5][2] = data[ind[1][1][2]];
      cVals[5][0] = data[ind[1][1][0]];
      cVals[5][3] = data[ind[1][1][3]];
      cVals[6][0] = data[ind[1][2][0]];
      cVals[6][1] = data[ind[1][2][1]];
      cVals[6][2] = data[ind[1][2][2]];
      cVals[6][3] = data[ind[1][2][3]];
      cVals[7][0] = data[ind[1][3][0]];
      cVals[7][1] = data[ind[1][3][1]];
      cVals[7][2] = data[ind[1][3][2]];
      cVals[7][3] = data[ind[1][3][3]];
      cVals[8][0] = data[ind[2][0][0]];
      cVals[8][1] = data[ind[2][0][1]];
      cVals[8][2] = data[ind[2][0][2]];
      cVals[8][3] = data[ind[2][0][3]];
      cVals[9][0] = data[ind[2][1][0]];
      cVals[9][1] = data[ind[2][1][1]];
      cVals[9][2] = data[ind[2][1][2]];
      cVals[9][3] = data[ind[2][1][3]];
      cVals[10][0] = data[ind[2][2][0]];
      cVals[10][1] = data[ind[2][2][1]];
      cVals[10][2] = data[ind[2][2][2]];
      cVals[10][3] = data[ind[2][2][3]];
      cVals[11][0] = data[ind[2][3][0]];
      cVals[11][1] = data[ind[2][3][1]];
      cVals[11][2] = data[ind[2][3][2]];
      cVals[11][3] = data[ind[2][3][3]];
      cVals[12][0] = data[ind[3][0][0]];
      cVals[12][1] = data[ind[3][0][1]];
      cVals[12][2] = data[ind[3][0][2]];
      cVals[12][3] = data[ind[3][0][3]];
      cVals[13][0] = data[ind[3][1][0]];
      cVals[13][1] = data[ind[3][1][1]];
      cVals[13][2] = data[ind[3][1][2]];
      cVals[13][3] = data[ind[3][1][3]];
      cVals[14][0] = data[ind[3][2][0]];
      cVals[14][1] = data[ind[3][2][1]];
      cVals[14][2] = data[ind[3][2][2]];
      cVals[14][3] = data[ind[3][2][3]];
      cVals[15][0] = data[ind[3][3][0]];
      cVals[15][1] = data[ind[3][3][1]];
      cVals[15][2] = data[ind[3][3][2]];
      cVals[15][3] = data[ind[3][3][3]];
    } break;

    case GridPos::SideXRight:
    case GridPos::LineD:
    case GridPos::LineH: {
      const int ind[4][4][3]{
        {{ii_x_y_z + deltaPhi[i0] + deltaR[i0] + deltaZ[i0], ind[0][0][0] - deltaZ[i0], ind[0][0][1] - deltaZ[i0]},
         {ind[0][0][0] - deltaR[i0], ind[0][1][0] - deltaZ[i0], ind[0][1][1] - deltaZ[i0]},
         {ind[0][1][0] - deltaR[i0], ind[0][2][0] - deltaZ[i0], ind[0][2][1] - deltaZ[i0]},
         {ind[0][2][0] - deltaR[i0], ind[0][3][0] - deltaZ[i0], ind[0][3][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaR[i0] + deltaZ[i0], ind[1][0][0] - deltaZ[i0], ind[1][0][1] - deltaZ[i0]},
         {ind[1][0][0] - deltaR[i0], ind[1][1][0] - deltaZ[i0], ind[1][1][1] - deltaZ[i0]},
         {ind[1][1][0] - deltaR[i0], ind[1][2][0] - deltaZ[i0], ind[1][2][1] - deltaZ[i0]},
         {ind[1][2][0] - deltaR[i0], ind[1][3][0] - deltaZ[i0], ind[1][3][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i1] + deltaR[i0] + deltaZ[i0], ind[2][0][0] - deltaZ[i0], ind[2][0][1] - deltaZ[i0]},
         {ind[2][0][0] - deltaR[i0], ind[2][1][0] - deltaZ[i0], ind[2][1][1] - deltaZ[i0]},
         {ind[2][1][0] - deltaR[i0], ind[2][2][0] - deltaZ[i0], ind[2][2][1] - deltaZ[i0]},
         {ind[2][2][0] - deltaR[i0], ind[2][3][0] - deltaZ[i0], ind[2][3][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i2] + deltaR[i0] + deltaZ[i0], ind[3][0][0] - deltaZ[i0], ind[3][0][1] - deltaZ[i0]},
         {ind[3][0][0] - deltaR[i0], ind[3][1][0] - deltaZ[i0], ind[3][1][1] - deltaZ[i0]},
         {ind[3][1][0] - deltaR[i0], ind[3][2][0] - deltaZ[i0], ind[3][2][1] - deltaZ[i0]},
         {ind[3][2][0] - deltaR[i0], ind[3][3][0] - deltaZ[i0], ind[3][3][1] - deltaZ[i0]}}};

      cVals[0][0] = data[ind[0][0][0]];
      cVals[0][1] = data[ind[0][0][1]];
      cVals[0][2] = data[ind[0][0][2]];
      cVals[0][3] = extrapolation(data[ind[0][0][2]], data[ind[0][0][1]], data[ind[0][0][0]]);
      cVals[1][0] = data[ind[0][1][0]];
      cVals[1][1] = data[ind[0][1][1]];
      cVals[1][2] = data[ind[0][1][2]];
      cVals[1][3] = extrapolation(data[ind[0][1][2]], data[ind[0][1][1]], data[ind[0][1][0]]);
      cVals[2][0] = data[ind[0][2][0]];
      cVals[2][1] = data[ind[0][2][1]];
      cVals[2][2] = data[ind[0][2][2]];
      cVals[2][3] = extrapolation(data[ind[0][2][2]], data[ind[0][2][1]], data[ind[0][2][0]]);
      cVals[3][0] = data[ind[0][3][0]];
      cVals[3][1] = data[ind[0][3][1]];
      cVals[3][2] = data[ind[0][3][2]];
      cVals[3][3] = extrapolation(data[ind[0][3][2]], data[ind[0][3][1]], data[ind[0][3][0]]);
      cVals[4][0] = data[ind[1][0][0]];
      cVals[4][1] = data[ind[1][0][1]];
      cVals[4][2] = data[ind[1][0][2]];
      cVals[4][3] = extrapolation(data[ind[1][0][2]], data[ind[1][0][1]], data[ind[1][0][0]]);
      cVals[5][0] = data[ind[1][1][0]];
      cVals[5][2] = data[ind[1][1][2]];
      cVals[5][3] = extrapolation(data[ind[1][1][2]], data[ii_x_y_z], data[ind[1][1][0]]);
      cVals[6][0] = data[ind[1][2][0]];
      cVals[6][1] = data[ind[1][2][1]];
      cVals[6][2] = data[ind[1][2][2]];
      cVals[6][3] = extrapolation(data[ind[1][2][2]], data[ind[1][2][1]], data[ind[1][2][0]]);
      cVals[7][0] = data[ind[1][3][0]];
      cVals[7][1] = data[ind[1][3][1]];
      cVals[7][2] = data[ind[1][3][2]];
      cVals[7][3] = extrapolation(data[ind[1][3][2]], data[ind[1][3][1]], data[ind[1][3][0]]);
      cVals[8][0] = data[ind[2][0][0]];
      cVals[8][1] = data[ind[2][0][1]];
      cVals[8][2] = data[ind[2][0][2]];
      cVals[8][3] = extrapolation(data[ind[2][0][2]], data[ind[2][0][1]], data[ind[2][0][0]]);
      cVals[9][0] = data[ind[2][1][0]];
      cVals[9][1] = data[ind[2][1][1]];
      cVals[9][2] = data[ind[2][1][2]];
      cVals[9][3] = extrapolation(data[ind[2][1][2]], data[ind[2][1][1]], data[ind[2][1][0]]);
      cVals[10][0] = data[ind[2][2][0]];
      cVals[10][1] = data[ind[2][2][1]];
      cVals[10][2] = data[ind[2][2][2]];
      cVals[10][3] = extrapolation(data[ind[2][2][2]], data[ind[2][2][1]], data[ind[2][2][0]]);
      cVals[11][0] = data[ind[2][3][0]];
      cVals[11][1] = data[ind[2][3][1]];
      cVals[11][2] = data[ind[2][3][2]];
      cVals[11][3] = extrapolation(data[ind[2][3][2]], data[ind[2][3][1]], data[ind[2][3][0]]);
      cVals[12][0] = data[ind[3][0][0]];
      cVals[12][1] = data[ind[3][0][1]];
      cVals[12][2] = data[ind[3][0][2]];
      cVals[13][0] = data[ind[3][1][0]];
      cVals[12][3] = extrapolation(data[ind[3][0][2]], data[ind[3][0][1]], data[ind[3][0][0]]);
      cVals[13][1] = data[ind[3][1][1]];
      cVals[13][2] = data[ind[3][1][2]];
      cVals[13][3] = extrapolation(data[ind[3][1][2]], data[ind[3][1][1]], data[ind[3][1][0]]);
      cVals[14][0] = data[ind[3][2][0]];
      cVals[14][1] = data[ind[3][2][1]];
      cVals[14][2] = data[ind[3][2][2]];
      cVals[14][3] = extrapolation(data[ind[3][2][2]], data[ind[3][2][1]], data[ind[3][2][0]]);
      cVals[15][0] = data[ind[3][3][0]];
      cVals[15][1] = data[ind[3][3][1]];
      cVals[15][2] = data[ind[3][3][2]];
      cVals[15][3] = extrapolation(data[ind[3][3][2]], data[ind[3][3][1]], data[ind[3][3][0]]);
    } break;

    case GridPos::SideYRight:
    case GridPos::LineB:
    case GridPos::LineF: {
      const int ind[4][3][4]{
        {{ii_x_y_z + deltaPhi[i0] + deltaR[i0] + deltaZ[i0], ind[0][0][0] - deltaZ[i0], ind[0][0][1] - deltaZ[i0], ind[0][0][2] - deltaZ[i0]},
         {ind[0][0][0] - deltaR[i0], ind[0][1][0] - deltaZ[i0], ind[0][1][1] - deltaZ[i0], ind[0][1][2] - deltaZ[i0]},
         {ind[0][1][0] - deltaR[i0], ind[0][2][0] - deltaZ[i0], ind[0][2][1] - deltaZ[i0], ind[0][2][2] - deltaZ[i0]}},
        {{ii_x_y_z + deltaR[i0] + deltaZ[i0], ind[1][0][0] - deltaZ[i0], ind[1][0][1] - deltaZ[i0], ind[1][0][2] - deltaZ[i0]},
         {ind[1][0][0] - deltaR[i0], ind[1][1][0] - deltaZ[i0], ind[1][1][1] - deltaZ[i0], ind[1][1][2] - deltaZ[i0]},
         {ind[1][1][0] - deltaR[i0], ind[1][2][0] - deltaZ[i0], ind[1][2][1] - deltaZ[i0], ind[1][2][2] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i1] + deltaR[i0] + deltaZ[i0], ind[2][0][0] - deltaZ[i0], ind[2][0][1] - deltaZ[i0], ind[2][0][2] - deltaZ[i0]},
         {ind[2][0][0] - deltaR[i0], ind[2][1][0] - deltaZ[i0], ind[2][1][1] - deltaZ[i0], ind[2][1][2] - deltaZ[i0]},
         {ind[2][1][0] - deltaR[i0], ind[2][2][0] - deltaZ[i0], ind[2][2][1] - deltaZ[i0], ind[2][2][2] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i2] + deltaR[i0] + deltaZ[i0], ind[3][0][0] - deltaZ[i0], ind[3][0][1] - deltaZ[i0], ind[3][0][2] - deltaZ[i0]},
         {ind[3][0][0] - deltaR[i0], ind[3][1][0] - deltaZ[i0], ind[3][1][1] - deltaZ[i0], ind[3][1][2] - deltaZ[i0]},
         {ind[3][1][0] - deltaR[i0], ind[3][2][0] - deltaZ[i0], ind[3][2][1] - deltaZ[i0], ind[3][2][2] - deltaZ[i0]}}};

      cVals[0][0] = data[ind[0][0][0]];
      cVals[0][1] = data[ind[0][0][1]];
      cVals[0][2] = data[ind[0][0][2]];
      cVals[0][3] = data[ind[0][0][3]];
      cVals[1][0] = data[ind[0][1][0]];
      cVals[1][1] = data[ind[0][1][1]];
      cVals[1][2] = data[ind[0][1][2]];
      cVals[1][3] = data[ind[0][1][3]];
      cVals[2][0] = data[ind[0][2][0]];
      cVals[2][1] = data[ind[0][2][1]];
      cVals[2][2] = data[ind[0][2][2]];
      cVals[2][3] = data[ind[0][2][3]];
      cVals[3][0] = extrapolation(data[ind[0][2][0]], data[ind[0][1][0]], data[ind[0][0][0]]);
      cVals[3][1] = extrapolation(data[ind[0][2][1]], data[ind[0][1][1]], data[ind[0][0][1]]);
      cVals[3][2] = extrapolation(data[ind[0][2][2]], data[ind[0][1][2]], data[ind[0][0][2]]);
      cVals[3][3] = extrapolation(data[ind[0][2][3]], data[ind[0][1][3]], data[ind[0][0][3]]);
      cVals[4][0] = data[ind[1][0][0]];
      cVals[4][1] = data[ind[1][0][1]];
      cVals[4][2] = data[ind[1][0][2]];
      cVals[4][3] = data[ind[1][0][3]];
      cVals[5][0] = data[ind[1][1][0]];
      cVals[5][2] = data[ind[1][1][2]];
      cVals[5][3] = data[ind[1][1][3]];
      cVals[6][0] = data[ind[1][2][0]];
      cVals[6][1] = data[ind[1][2][1]];
      cVals[6][2] = data[ind[1][2][2]];
      cVals[6][3] = data[ind[1][2][3]];
      cVals[7][0] = extrapolation(data[ind[1][2][0]], data[ind[1][1][0]], data[ind[1][0][0]]);
      cVals[7][1] = extrapolation(data[ind[1][2][1]], data[ii_x_y_z], data[ind[1][0][1]]);
      cVals[7][2] = extrapolation(data[ind[1][2][2]], data[ind[1][1][2]], data[ind[1][0][2]]);
      cVals[7][3] = extrapolation(data[ind[1][2][3]], data[ind[1][1][3]], data[ind[1][0][3]]);
      cVals[8][0] = data[ind[2][0][0]];
      cVals[8][1] = data[ind[2][0][1]];
      cVals[8][2] = data[ind[2][0][2]];
      cVals[8][3] = data[ind[2][0][3]];
      cVals[9][0] = data[ind[2][1][0]];
      cVals[9][1] = data[ind[2][1][1]];
      cVals[9][2] = data[ind[2][1][2]];
      cVals[9][3] = data[ind[2][1][3]];
      cVals[10][0] = data[ind[2][2][0]];
      cVals[10][1] = data[ind[2][2][1]];
      cVals[10][2] = data[ind[2][2][2]];
      cVals[10][3] = data[ind[2][2][3]];
      cVals[11][0] = extrapolation(data[ind[2][2][0]], data[ind[2][1][0]], data[ind[2][0][0]]);
      cVals[11][1] = extrapolation(data[ind[2][2][1]], data[ind[2][1][1]], data[ind[2][0][1]]);
      cVals[11][2] = extrapolation(data[ind[2][2][2]], data[ind[2][1][2]], data[ind[2][0][2]]);
      cVals[11][3] = extrapolation(data[ind[2][2][3]], data[ind[2][1][3]], data[ind[2][0][3]]);
      cVals[12][0] = data[ind[3][0][0]];
      cVals[12][1] = data[ind[3][0][1]];
      cVals[12][2] = data[ind[3][0][2]];
      cVals[12][3] = data[ind[3][0][3]];
      cVals[13][0] = data[ind[3][1][0]];
      cVals[13][1] = data[ind[3][1][1]];
      cVals[13][2] = data[ind[3][1][2]];
      cVals[13][3] = data[ind[3][1][3]];
      cVals[14][0] = data[ind[3][2][0]];
      cVals[14][1] = data[ind[3][2][1]];
      cVals[14][2] = data[ind[3][2][2]];
      cVals[14][3] = data[ind[3][2][3]];
      cVals[15][0] = extrapolation(data[ind[3][2][0]], data[ind[3][1][0]], data[ind[3][0][0]]);
      cVals[15][1] = extrapolation(data[ind[3][2][1]], data[ind[3][1][1]], data[ind[3][0][1]]);
      cVals[15][2] = extrapolation(data[ind[3][2][2]], data[ind[3][1][2]], data[ind[3][0][2]]);
      cVals[15][3] = extrapolation(data[ind[3][2][3]], data[ind[3][1][3]], data[ind[3][0][3]]);
    } break;

    case GridPos::SideYLeft:
    case GridPos::LineA:
    case GridPos::LineE: {
      const int ind[4][3][4]{
        {{ii_x_y_z + deltaPhi[i0] + deltaZ[i0], ind[0][0][0] - deltaZ[i0], ind[0][0][1] - deltaZ[i0], ind[0][0][2] - deltaZ[i0]},
         {ind[0][0][0] - deltaR[i0], ind[0][1][0] - deltaZ[i0], ind[0][1][1] - deltaZ[i0], ind[0][1][2] - deltaZ[i0]},
         {ind[0][1][0] - deltaR[i0], ind[0][2][0] - deltaZ[i0], ind[0][2][1] - deltaZ[i0], ind[0][2][2] - deltaZ[i0]}},
        {{ii_x_y_z + deltaZ[i0], ind[1][0][0] - deltaZ[i0], ind[1][0][1] - deltaZ[i0], ind[1][0][2] - deltaZ[i0]},
         {ind[1][0][0] - deltaR[i0], ind[1][1][0] - deltaZ[i0], ind[1][1][1] - deltaZ[i0], ind[1][1][2] - deltaZ[i0]},
         {ind[1][1][0] - deltaR[i0], ind[1][2][0] - deltaZ[i0], ind[1][2][1] - deltaZ[i0], ind[1][2][2] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i1] + deltaZ[i0], ind[2][0][0] - deltaZ[i0], ind[2][0][1] - deltaZ[i0], ind[2][0][2] - deltaZ[i0]},
         {ind[2][0][0] - deltaR[i0], ind[2][1][0] - deltaZ[i0], ind[2][1][1] - deltaZ[i0], ind[2][1][2] - deltaZ[i0]},
         {ind[2][1][0] - deltaR[i0], ind[2][2][0] - deltaZ[i0], ind[2][2][1] - deltaZ[i0], ind[2][2][2] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i2] + deltaZ[i0], ind[3][0][0] - deltaZ[i0], ind[3][0][1] - deltaZ[i0], ind[3][0][2] - deltaZ[i0]},
         {ind[3][0][0] - deltaR[i0], ind[3][1][0] - deltaZ[i0], ind[3][1][1] - deltaZ[i0], ind[3][1][2] - deltaZ[i0]},
         {ind[3][1][0] - deltaR[i0], ind[3][2][0] - deltaZ[i0], ind[3][2][1] - deltaZ[i0], ind[3][2][2] - deltaZ[i0]}}};

      cVals[0][0] = extrapolation(data[ind[0][0][0]], data[ind[0][1][0]], data[ind[0][2][0]]);
      cVals[0][1] = extrapolation(data[ind[0][0][1]], data[ind[0][1][1]], data[ind[0][2][1]]);
      cVals[0][2] = extrapolation(data[ind[0][0][2]], data[ind[0][1][2]], data[ind[0][2][2]]);
      cVals[0][3] = extrapolation(data[ind[0][0][3]], data[ind[0][1][3]], data[ind[0][2][3]]);
      cVals[1][0] = data[ind[0][0][0]];
      cVals[1][1] = data[ind[0][0][1]];
      cVals[1][2] = data[ind[0][0][2]];
      cVals[1][3] = data[ind[0][0][3]];
      cVals[2][0] = data[ind[0][1][0]];
      cVals[2][1] = data[ind[0][1][1]];
      cVals[2][2] = data[ind[0][1][2]];
      cVals[2][3] = data[ind[0][1][3]];
      cVals[3][0] = data[ind[0][2][0]];
      cVals[3][1] = data[ind[0][2][1]];
      cVals[3][2] = data[ind[0][2][2]];
      cVals[3][3] = data[ind[0][2][3]];
      cVals[4][0] = extrapolation(data[ind[1][0][0]], data[ind[1][1][0]], data[ind[1][2][0]]);
      cVals[4][1] = extrapolation(data[ii_x_y_z], data[ind[1][1][1]], data[ind[1][2][1]]);
      cVals[4][2] = extrapolation(data[ind[1][0][2]], data[ind[1][1][2]], data[ind[1][2][2]]);
      cVals[4][3] = extrapolation(data[ind[1][0][3]], data[ind[1][1][3]], data[ind[1][2][3]]);
      cVals[5][0] = data[ind[1][0][0]];
      cVals[5][2] = data[ind[1][0][2]];
      cVals[5][3] = data[ind[1][0][3]];
      cVals[6][0] = data[ind[1][1][0]];
      cVals[6][1] = data[ind[1][1][1]];
      cVals[6][2] = data[ind[1][1][2]];
      cVals[6][3] = data[ind[1][1][3]];
      cVals[7][0] = data[ind[1][2][0]];
      cVals[7][1] = data[ind[1][2][1]];
      cVals[7][2] = data[ind[1][2][2]];
      cVals[7][3] = data[ind[1][2][3]];
      cVals[8][0] = extrapolation(data[ind[2][0][0]], data[ind[2][1][0]], data[ind[2][2][0]]);
      cVals[8][1] = extrapolation(data[ind[2][0][1]], data[ind[2][1][1]], data[ind[2][2][1]]);
      cVals[8][2] = extrapolation(data[ind[2][0][2]], data[ind[2][1][2]], data[ind[2][2][2]]);
      cVals[8][3] = extrapolation(data[ind[2][0][3]], data[ind[2][1][3]], data[ind[2][2][3]]);
      cVals[9][0] = data[ind[2][0][0]];
      cVals[9][1] = data[ind[2][0][1]];
      cVals[9][2] = data[ind[2][0][2]];
      cVals[9][3] = data[ind[2][0][3]];
      cVals[10][0] = data[ind[2][1][0]];
      cVals[10][1] = data[ind[2][1][1]];
      cVals[10][2] = data[ind[2][1][2]];
      cVals[10][3] = data[ind[2][1][3]];
      cVals[11][0] = data[ind[2][2][0]];
      cVals[11][1] = data[ind[2][2][1]];
      cVals[11][2] = data[ind[2][2][2]];
      cVals[11][3] = data[ind[2][2][3]];
      cVals[12][0] = extrapolation(data[ind[3][0][0]], data[ind[3][1][0]], data[ind[3][2][0]]);
      cVals[12][1] = extrapolation(data[ind[3][0][1]], data[ind[3][1][1]], data[ind[3][2][1]]);
      cVals[12][2] = extrapolation(data[ind[3][0][2]], data[ind[3][1][2]], data[ind[3][2][2]]);
      cVals[12][3] = extrapolation(data[ind[3][0][3]], data[ind[3][1][3]], data[ind[3][2][3]]);
      cVals[13][0] = data[ind[3][0][0]];
      cVals[13][1] = data[ind[3][0][1]];
      cVals[13][2] = data[ind[3][0][2]];
      cVals[13][3] = data[ind[3][0][3]];
      cVals[14][0] = data[ind[3][1][0]];
      cVals[14][1] = data[ind[3][1][1]];
      cVals[14][2] = data[ind[3][1][2]];
      cVals[14][3] = data[ind[3][1][3]];
      cVals[15][0] = data[ind[3][2][0]];
      cVals[15][1] = data[ind[3][2][1]];
      cVals[15][2] = data[ind[3][2][2]];
      cVals[15][3] = data[ind[3][2][3]];
    } break;

    case GridPos::SideXLeft:
    case GridPos::LineC:
    case GridPos::LineG: {
      const int ind[4][4][3]{
        {{ii_x_y_z + deltaPhi[i0] + deltaR[i0], ind[0][0][0] - deltaZ[i0], ind[0][0][1] - deltaZ[i0]},
         {ind[0][0][0] - deltaR[i0], ind[0][1][0] - deltaZ[i0], ind[0][1][1] - deltaZ[i0]},
         {ind[0][1][0] - deltaR[i0], ind[0][2][0] - deltaZ[i0], ind[0][2][1] - deltaZ[i0]},
         {ind[0][2][0] - deltaR[i0], ind[0][3][0] - deltaZ[i0], ind[0][3][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaR[i0], ind[1][0][0] - deltaZ[i0], ind[1][0][1] - deltaZ[i0]},
         {ind[1][0][0] - deltaR[i0], ind[1][1][0] - deltaZ[i0], ind[1][1][1] - deltaZ[i0]},
         {ind[1][1][0] - deltaR[i0], ind[1][2][0] - deltaZ[i0], ind[1][2][1] - deltaZ[i0]},
         {ind[1][2][0] - deltaR[i0], ind[1][3][0] - deltaZ[i0], ind[1][3][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i1] + deltaR[i0], ind[2][0][0] - deltaZ[i0], ind[2][0][1] - deltaZ[i0]},
         {ind[2][0][0] - deltaR[i0], ind[2][1][0] - deltaZ[i0], ind[2][1][1] - deltaZ[i0]},
         {ind[2][1][0] - deltaR[i0], ind[2][2][0] - deltaZ[i0], ind[2][2][1] - deltaZ[i0]},
         {ind[2][2][0] - deltaR[i0], ind[2][3][0] - deltaZ[i0], ind[2][3][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i2] + deltaR[i0], ind[3][0][0] - deltaZ[i0], ind[3][0][1] - deltaZ[i0]},
         {ind[3][0][0] - deltaR[i0], ind[3][1][0] - deltaZ[i0], ind[3][1][1] - deltaZ[i0]},
         {ind[3][1][0] - deltaR[i0], ind[3][2][0] - deltaZ[i0], ind[3][2][1] - deltaZ[i0]},
         {ind[3][2][0] - deltaR[i0], ind[3][3][0] - deltaZ[i0], ind[3][3][1] - deltaZ[i0]}}};

      cVals[0][0] = extrapolation(data[ind[0][0][0]], data[ind[0][0][1]], data[ind[0][0][2]]);
      cVals[0][1] = data[ind[0][0][0]];
      cVals[0][2] = data[ind[0][0][1]];
      cVals[0][3] = data[ind[0][0][2]];
      cVals[1][0] = extrapolation(data[ind[0][1][0]], data[ind[0][1][1]], data[ind[0][1][2]]);
      cVals[1][1] = data[ind[0][1][0]];
      cVals[1][2] = data[ind[0][1][1]];
      cVals[1][3] = data[ind[0][1][2]];
      cVals[2][0] = extrapolation(data[ind[0][2][0]], data[ind[0][2][1]], data[ind[0][2][2]]);
      cVals[2][1] = data[ind[0][2][0]];
      cVals[2][2] = data[ind[0][2][1]];
      cVals[2][3] = data[ind[0][2][2]];
      cVals[3][0] = extrapolation(data[ind[0][3][0]], data[ind[0][3][1]], data[ind[0][3][2]]);
      cVals[3][1] = data[ind[0][3][0]];
      cVals[3][2] = data[ind[0][3][1]];
      cVals[3][3] = data[ind[0][3][2]];
      cVals[4][0] = extrapolation(data[ind[1][0][0]], data[ind[1][0][1]], data[ind[1][0][2]]);
      cVals[4][1] = data[ind[1][0][0]];
      cVals[4][2] = data[ind[1][0][1]];
      cVals[4][3] = data[ind[1][0][2]];
      cVals[5][0] = extrapolation(data[ii_x_y_z], data[ind[1][1][1]], data[ind[1][1][2]]);
      cVals[5][2] = data[ind[1][1][1]];
      cVals[5][3] = data[ind[1][1][2]];
      cVals[6][0] = extrapolation(data[ind[1][2][0]], data[ind[1][2][1]], data[ind[1][2][2]]);
      cVals[6][1] = data[ind[1][2][0]];
      cVals[6][2] = data[ind[1][2][1]];
      cVals[6][3] = data[ind[1][2][2]];
      cVals[7][0] = extrapolation(data[ind[1][3][0]], data[ind[1][3][1]], data[ind[1][3][2]]);
      cVals[7][1] = data[ind[1][3][0]];
      cVals[7][2] = data[ind[1][3][1]];
      cVals[7][3] = data[ind[1][3][2]];
      cVals[8][0] = extrapolation(data[ind[2][0][0]], data[ind[2][0][1]], data[ind[2][0][2]]);
      cVals[8][1] = data[ind[2][0][0]];
      cVals[8][2] = data[ind[2][0][1]];
      cVals[8][3] = data[ind[2][0][2]];
      cVals[9][0] = extrapolation(data[ind[2][1][0]], data[ind[2][1][1]], data[ind[2][1][2]]);
      cVals[9][1] = data[ind[2][1][0]];
      cVals[9][2] = data[ind[2][1][1]];
      cVals[9][3] = data[ind[2][1][2]];
      cVals[10][0] = extrapolation(data[ind[2][2][0]], data[ind[2][2][1]], data[ind[2][2][2]]);
      cVals[10][1] = data[ind[2][2][0]];
      cVals[10][2] = data[ind[2][2][1]];
      cVals[10][3] = data[ind[2][2][2]];
      cVals[11][0] = extrapolation(data[ind[2][3][0]], data[ind[2][3][1]], data[ind[2][3][2]]);
      cVals[11][1] = data[ind[2][3][0]];
      cVals[11][2] = data[ind[2][3][1]];
      cVals[11][3] = data[ind[2][3][2]];
      cVals[12][0] = extrapolation(data[ind[3][0][0]], data[ind[3][0][1]], data[ind[3][0][2]]);
      cVals[12][1] = data[ind[3][0][0]];
      cVals[12][2] = data[ind[3][0][1]];
      cVals[12][3] = data[ind[3][0][2]];
      cVals[13][0] = extrapolation(data[ind[3][1][0]], data[ind[3][1][1]], data[ind[3][1][2]]);
      cVals[13][1] = data[ind[3][1][0]];
      cVals[13][2] = data[ind[3][1][1]];
      cVals[13][3] = data[ind[3][1][2]];
      cVals[14][0] = extrapolation(data[ind[3][2][0]], data[ind[3][2][1]], data[ind[3][2][2]]);
      cVals[14][1] = data[ind[3][2][0]];
      cVals[14][2] = data[ind[3][2][1]];
      cVals[14][3] = data[ind[3][2][2]];
      cVals[15][0] = extrapolation(data[ind[3][3][0]], data[ind[3][3][1]], data[ind[3][3][2]]);
      cVals[15][1] = data[ind[3][3][0]];
      cVals[15][2] = data[ind[3][3][1]];
      cVals[15][3] = data[ind[3][3][2]];
    } break;

    case GridPos::Edge0:
    case GridPos::Edge4:
    case GridPos::LineI: {
      const int ind[4][3][3]{
        {{ii_x_y_z + deltaPhi[i0], ind[0][0][0] - deltaZ[i0], ind[0][0][1] - deltaZ[i0]},
         {ind[0][0][0] - deltaR[i0], ind[0][1][0] - deltaZ[i0], ind[0][1][1] - deltaZ[i0]},
         {ind[0][1][0] - deltaR[i0], ind[0][2][0] - deltaZ[i0], ind[0][2][1] - deltaZ[i0]}},
        {{ii_x_y_z, ind[1][0][0] - deltaZ[i0], ind[1][0][1] - deltaZ[i0]},
         {ind[1][0][0] - deltaR[i0], ind[1][1][0] - deltaZ[i0], ind[1][1][1] - deltaZ[i0]},
         {ind[1][1][0] - deltaR[i0], ind[1][2][0] - deltaZ[i0], ind[1][2][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i1], ind[2][0][0] - deltaZ[i0], ind[2][0][1] - deltaZ[i0]},
         {ind[2][0][0] - deltaR[i0], ind[2][1][0] - deltaZ[i0], ind[2][1][1] - deltaZ[i0]},
         {ind[2][1][0] - deltaR[i0], ind[2][2][0] - deltaZ[i0], ind[2][2][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i2], ind[3][0][0] - deltaZ[i0], ind[3][0][1] - deltaZ[i0]},
         {ind[3][0][0] - deltaR[i0], ind[3][1][0] - deltaZ[i0], ind[3][1][1] - deltaZ[i0]},
         {ind[3][1][0] - deltaR[i0], ind[3][2][0] - deltaZ[i0], ind[3][2][1] - deltaZ[i0]}}};

      cVals[0][0] = extrapolation(data[ind[0][0][0]], data[ind[0][1][1]], data[ind[0][2][2]]);
      cVals[0][1] = extrapolation(data[ind[0][0][0]], data[ind[0][1][0]], data[ind[0][2][0]]);
      cVals[0][2] = extrapolation(data[ind[0][0][1]], data[ind[0][1][1]], data[ind[0][2][1]]);
      cVals[0][3] = extrapolation(data[ind[0][0][2]], data[ind[0][1][2]], data[ind[0][2][2]]);
      cVals[1][0] = extrapolation(data[ind[0][0][0]], data[ind[0][0][1]], data[ind[0][0][2]]);
      cVals[1][1] = data[ind[0][0][0]];
      cVals[1][2] = data[ind[0][0][1]];
      cVals[1][3] = data[ind[0][0][2]];
      cVals[2][0] = extrapolation(data[ind[0][1][0]], data[ind[0][1][1]], data[ind[0][1][2]]);
      cVals[2][1] = data[ind[0][1][0]];
      cVals[2][2] = data[ind[0][1][1]];
      cVals[2][3] = data[ind[0][1][2]];
      cVals[3][0] = extrapolation(data[ind[0][2][0]], data[ind[0][2][1]], data[ind[0][2][2]]);
      cVals[3][1] = data[ind[0][2][0]];
      cVals[3][2] = data[ind[0][2][1]];
      cVals[3][3] = data[ind[0][2][2]];
      cVals[4][0] = extrapolation(data[ii_x_y_z], data[ind[1][1][1]], data[ind[1][2][2]]);
      cVals[4][1] = extrapolation(data[ii_x_y_z], data[ind[1][1][0]], data[ind[1][2][0]]);
      cVals[4][2] = extrapolation(data[ind[1][0][1]], data[ind[1][1][1]], data[ind[1][2][1]]);
      cVals[4][3] = extrapolation(data[ind[1][0][2]], data[ind[1][1][2]], data[ind[1][2][2]]);
      cVals[5][0] = extrapolation(data[ii_x_y_z], data[ind[1][0][1]], data[ind[1][0][2]]);
      cVals[5][2] = data[ind[1][0][1]];
      cVals[5][3] = data[ind[1][0][2]];
      cVals[6][0] = extrapolation(data[ind[1][1][0]], data[ind[1][1][1]], data[ind[1][1][2]]);
      cVals[6][1] = data[ind[1][1][0]];
      cVals[6][2] = data[ind[1][1][1]];
      cVals[6][3] = data[ind[1][1][2]];
      cVals[7][0] = extrapolation(data[ind[1][2][0]], data[ind[1][2][1]], data[ind[1][2][2]]);
      cVals[7][1] = data[ind[1][2][0]];
      cVals[7][2] = data[ind[1][2][1]];
      cVals[7][3] = data[ind[1][2][2]];
      cVals[8][0] = extrapolation(data[ind[2][0][0]], data[ind[2][1][1]], data[ind[2][2][2]]);
      cVals[8][1] = extrapolation(data[ind[2][0][0]], data[ind[2][1][0]], data[ind[2][2][0]]);
      cVals[8][2] = extrapolation(data[ind[2][0][1]], data[ind[2][1][1]], data[ind[2][2][1]]);
      cVals[8][3] = extrapolation(data[ind[2][0][2]], data[ind[2][1][2]], data[ind[2][2][2]]);
      cVals[9][0] = extrapolation(data[ind[2][0][0]], data[ind[2][0][1]], data[ind[2][0][2]]);
      cVals[9][1] = data[ind[2][0][0]];
      cVals[9][2] = data[ind[2][0][1]];
      cVals[9][3] = data[ind[2][0][2]];
      cVals[10][0] = extrapolation(data[ind[2][1][0]], data[ind[2][1][1]], data[ind[2][1][2]]);
      cVals[10][1] = data[ind[2][1][0]];
      cVals[10][2] = data[ind[2][1][1]];
      cVals[10][3] = data[ind[2][1][2]];
      cVals[11][0] = extrapolation(data[ind[2][2][0]], data[ind[2][2][1]], data[ind[2][2][2]]);
      cVals[11][1] = data[ind[2][2][0]];
      cVals[11][2] = data[ind[2][2][1]];
      cVals[11][3] = data[ind[2][2][2]];
      cVals[12][0] = extrapolation(data[ind[3][0][0]], data[ind[3][1][1]], data[ind[3][2][2]]);
      cVals[12][1] = extrapolation(data[ind[3][0][0]], data[ind[3][1][0]], data[ind[3][2][0]]);
      cVals[12][2] = extrapolation(data[ind[3][0][1]], data[ind[3][1][1]], data[ind[3][2][1]]);
      cVals[12][3] = extrapolation(data[ind[3][0][2]], data[ind[3][1][2]], data[ind[3][2][2]]);
      cVals[13][0] = extrapolation(data[ind[3][0][0]], data[ind[3][0][1]], data[ind[3][0][2]]);
      cVals[13][1] = data[ind[3][0][0]];
      cVals[13][2] = data[ind[3][0][1]];
      cVals[13][3] = data[ind[3][0][2]];
      cVals[14][0] = extrapolation(data[ind[3][1][0]], data[ind[3][1][1]], data[ind[3][1][2]]);
      cVals[14][1] = data[ind[3][1][0]];
      cVals[14][2] = data[ind[3][1][1]];
      cVals[14][3] = data[ind[3][1][2]];
      cVals[15][0] = extrapolation(data[ind[3][2][0]], data[ind[3][2][1]], data[ind[3][2][2]]);
      cVals[15][1] = data[ind[3][2][0]];
      cVals[15][2] = data[ind[3][2][1]];
      cVals[15][3] = data[ind[3][2][2]];
    } break;

    case GridPos::Edge1:
    case GridPos::Edge5:
    case GridPos::LineJ: {
      const int ind[4][3][3]{
        {{ii_x_y_z + deltaPhi[i0] + deltaZ[i0], ind[0][0][0] - deltaZ[i0], ind[0][0][1] - deltaZ[i0]},
         {ind[0][0][0] - deltaR[i0], ind[0][1][0] - deltaZ[i0], ind[0][1][1] - deltaZ[i0]},
         {ind[0][1][0] - deltaR[i0], ind[0][2][0] - deltaZ[i0], ind[0][2][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaZ[i0], ind[1][0][0] - deltaZ[i0], ind[1][0][1] - deltaZ[i0]},
         {ind[1][0][0] - deltaR[i0], ind[1][1][0] - deltaZ[i0], ind[1][1][1] - deltaZ[i0]},
         {ind[1][1][0] - deltaR[i0], ind[1][2][0] - deltaZ[i0], ind[1][2][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i1] + deltaZ[i0], ind[2][0][0] - deltaZ[i0], ind[2][0][1] - deltaZ[i0]},
         {ind[2][0][0] - deltaR[i0], ind[2][1][0] - deltaZ[i0], ind[2][1][1] - deltaZ[i0]},
         {ind[2][1][0] - deltaR[i0], ind[2][2][0] - deltaZ[i0], ind[2][2][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i2] + deltaZ[i0], ind[3][0][0] - deltaZ[i0], ind[3][0][1] - deltaZ[i0]},
         {ind[3][0][0] - deltaR[i0], ind[3][1][0] - deltaZ[i0], ind[3][1][1] - deltaZ[i0]},
         {ind[3][1][0] - deltaR[i0], ind[3][2][0] - deltaZ[i0], ind[3][2][1] - deltaZ[i0]}}};

      cVals[0][0] = extrapolation(data[ind[0][0][0]], data[ind[0][1][0]], data[ind[0][2][0]]);
      cVals[0][1] = extrapolation(data[ind[0][0][1]], data[ind[0][1][1]], data[ind[0][2][1]]);
      cVals[0][2] = extrapolation(data[ind[0][0][1]], data[ind[0][1][0]], data[ind[0][2][0] + deltaZ[i0]]);
      cVals[0][3] = extrapolation(data[ind[0][0][2]], data[ind[0][1][1]], data[ind[0][2][0]]);
      cVals[1][0] = data[ind[0][0][0]];
      cVals[1][1] = data[ind[0][0][1]];
      cVals[1][2] = data[ind[0][0][2]];
      cVals[1][3] = extrapolation(data[ind[0][0][2]], data[ind[0][0][1]], data[ind[0][0][0]]);
      cVals[2][0] = data[ind[0][1][0]];
      cVals[2][1] = data[ind[0][1][1]];
      cVals[2][2] = data[ind[0][1][2]];
      cVals[2][3] = extrapolation(data[ind[0][1][2]], data[ind[0][1][1]], data[ind[0][1][0]]);
      cVals[3][0] = data[ind[0][2][0]];
      cVals[3][1] = data[ind[0][2][1]];
      cVals[3][2] = data[ind[0][2][2]];
      cVals[3][3] = extrapolation(data[ind[0][2][2]], data[ind[0][2][1]], data[ind[0][2][0]]);
      cVals[4][0] = extrapolation(data[ind[1][0][0]], data[ind[1][1][0]], data[ind[1][2][0]]);
      cVals[4][1] = extrapolation(data[ii_x_y_z], data[ind[1][1][1]], data[ind[1][2][1]]);
      cVals[4][2] = extrapolation(data[ii_x_y_z], data[ind[1][1][0]], data[ind[1][2][0] + deltaZ[i0]]);
      cVals[4][3] = extrapolation(data[ind[1][0][2]], data[ind[1][1][1]], data[ind[1][2][0]]);
      cVals[5][0] = data[ind[1][0][0]];
      cVals[5][2] = data[ind[1][0][2]];
      cVals[5][3] = extrapolation(data[ind[1][0][2]], data[ii_x_y_z], data[ind[1][0][0]]);
      cVals[6][0] = data[ind[1][1][0]];
      cVals[6][1] = data[ind[1][1][1]];
      cVals[6][2] = data[ind[1][1][2]];
      cVals[6][3] = extrapolation(data[ind[1][1][2]], data[ind[1][1][1]], data[ind[1][1][0]]);
      cVals[7][0] = data[ind[1][2][0]];
      cVals[7][1] = data[ind[1][2][1]];
      cVals[7][2] = data[ind[1][2][2]];
      cVals[7][3] = extrapolation(data[ind[1][2][2]], data[ind[1][2][1]], data[ind[1][2][0]]);
      cVals[8][0] = extrapolation(data[ind[2][0][0]], data[ind[2][1][0]], data[ind[2][2][0]]);
      cVals[8][1] = extrapolation(data[ind[2][0][1]], data[ind[2][1][1]], data[ind[2][2][1]]);
      cVals[8][2] = extrapolation(data[ind[2][0][1]], data[ind[2][1][0]], data[ind[2][2][0] + deltaZ[i0]]);
      cVals[8][3] = extrapolation(data[ind[2][0][2]], data[ind[2][1][1]], data[ind[2][2][0]]);
      cVals[9][0] = data[ind[2][0][0]];
      cVals[9][1] = data[ind[2][0][1]];
      cVals[9][2] = data[ind[2][0][2]];
      cVals[9][3] = extrapolation(data[ind[2][0][2]], data[ind[2][0][1]], data[ind[2][0][0]]);
      cVals[10][0] = data[ind[2][1][0]];
      cVals[10][1] = data[ind[2][1][1]];
      cVals[10][2] = data[ind[2][1][2]];
      cVals[10][3] = extrapolation(data[ind[2][1][2]], data[ind[2][1][1]], data[ind[2][1][0]]);
      cVals[11][0] = data[ind[2][2][0]];
      cVals[11][1] = data[ind[2][2][1]];
      cVals[11][2] = data[ind[2][2][2]];
      cVals[11][3] = extrapolation(data[ind[2][2][2]], data[ind[2][2][1]], data[ind[2][2][0]]);
      cVals[12][0] = extrapolation(data[ind[3][0][0]], data[ind[3][1][0]], data[ind[3][2][0]]);
      cVals[12][1] = extrapolation(data[ind[3][0][1]], data[ind[3][1][1]], data[ind[3][2][1]]);
      cVals[12][2] = extrapolation(data[ind[3][0][1]], data[ind[3][1][0]], data[ind[3][2][0] + deltaZ[i0]]);
      cVals[12][3] = extrapolation(data[ind[3][0][2]], data[ind[3][1][1]], data[ind[3][2][0]]);
      cVals[13][0] = data[ind[3][0][0]];
      cVals[13][1] = data[ind[3][0][1]];
      cVals[13][2] = data[ind[3][0][2]];
      cVals[13][3] = extrapolation(data[ind[3][0][2]], data[ind[3][0][1]], data[ind[3][0][0]]);
      cVals[14][0] = data[ind[3][1][0]];
      cVals[14][1] = data[ind[3][1][1]];
      cVals[14][2] = data[ind[3][1][2]];
      cVals[14][3] = extrapolation(data[ind[3][1][2]], data[ind[3][1][1]], data[ind[3][1][0]]);
      cVals[15][0] = data[ind[3][2][0]];
      cVals[15][1] = data[ind[3][2][1]];
      cVals[15][2] = data[ind[3][2][2]];
      cVals[15][3] = extrapolation(data[ind[3][2][2]], data[ind[3][2][1]], data[ind[3][2][0]]);
    } break;

    case GridPos::Edge2:
    case GridPos::Edge6:
    case GridPos::LineK: {
      const int ind[4][3][3]{
        {{ii_x_y_z + deltaPhi[i0] + deltaR[i0], ind[0][0][0] - deltaZ[i0], ind[0][0][1] - deltaZ[i0]},
         {ind[0][0][0] - deltaR[i0], ind[0][1][0] - deltaZ[i0], ind[0][1][1] - deltaZ[i0]},
         {ind[0][1][0] - deltaR[i0], ind[0][2][0] - deltaZ[i0], ind[0][2][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaR[i0], ind[1][0][0] - deltaZ[i0], ind[1][0][1] - deltaZ[i0]},
         {ind[1][0][0] - deltaR[i0], ind[1][1][0] - deltaZ[i0], ind[1][1][1] - deltaZ[i0]},
         {ind[1][1][0] - deltaR[i0], ind[1][2][0] - deltaZ[i0], ind[1][2][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i1] + deltaR[i0], ind[2][0][0] - deltaZ[i0], ind[2][0][1] - deltaZ[i0]},
         {ind[2][0][0] - deltaR[i0], ind[2][1][0] - deltaZ[i0], ind[2][1][1] - deltaZ[i0]},
         {ind[2][1][0] - deltaR[i0], ind[2][2][0] - deltaZ[i0], ind[2][2][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i2] + deltaR[i0], ind[3][0][0] - deltaZ[i0], ind[3][0][1] - deltaZ[i0]},
         {ind[3][0][0] - deltaR[i0], ind[3][1][0] - deltaZ[i0], ind[3][1][1] - deltaZ[i0]},
         {ind[3][1][0] - deltaR[i0], ind[3][2][0] - deltaZ[i0], ind[3][2][1] - deltaZ[i0]}}};

      cVals[0][0] = extrapolation(data[ind[0][0][0]], data[ind[0][0][1]], data[ind[0][0][2]]);
      cVals[0][1] = data[ind[0][0][0]];
      cVals[0][2] = data[ind[0][0][1]];
      cVals[0][3] = data[ind[0][0][2]];
      cVals[1][0] = extrapolation(data[ind[0][1][0]], data[ind[0][1][1]], data[ind[0][1][2]]);
      cVals[1][1] = data[ind[0][1][0]];
      cVals[1][2] = data[ind[0][1][1]];
      cVals[1][3] = data[ind[0][1][2]];
      cVals[2][0] = extrapolation(data[ind[0][1][0]], data[ind[0][0][1]], data[ind[0][0][2] + deltaR[i0]]);
      cVals[2][1] = data[ind[0][2][0]];
      cVals[2][2] = data[ind[0][2][1]];
      cVals[2][3] = data[ind[0][2][2]];
      cVals[3][0] = extrapolation(data[ind[0][2][0]], data[ind[0][1][1]], data[ind[0][0][2]]);
      cVals[3][1] = extrapolation(data[ind[0][2][0]], data[ind[0][1][0]], data[ind[0][0][0]]);
      cVals[3][2] = extrapolation(data[ind[0][2][1]], data[ind[0][1][1]], data[ind[0][0][1]]);
      cVals[3][3] = extrapolation(data[ind[0][2][2]], data[ind[0][1][2]], data[ind[0][0][2]]);
      cVals[4][0] = extrapolation(data[ind[1][0][0]], data[ind[1][0][1]], data[ind[1][0][2]]);
      cVals[4][1] = data[ind[1][0][0]];
      cVals[4][2] = data[ind[1][0][1]];
      cVals[4][3] = data[ind[1][0][2]];
      cVals[5][0] = extrapolation(data[ii_x_y_z], data[ind[1][1][1]], data[ind[1][1][2]]);
      cVals[5][2] = data[ind[1][1][1]];
      cVals[5][3] = data[ind[1][1][2]];
      cVals[6][0] = extrapolation(data[ii_x_y_z], data[ind[1][0][1]], data[ind[1][0][2] + deltaR[i0]]);
      cVals[6][1] = data[ind[1][2][0]];
      cVals[6][2] = data[ind[1][2][1]];
      cVals[6][3] = data[ind[1][2][2]];
      cVals[7][0] = extrapolation(data[ind[1][2][0]], data[ind[1][1][1]], data[ind[1][0][2]]);
      cVals[7][1] = extrapolation(data[ind[1][2][0]], data[ii_x_y_z], data[ind[1][0][0]]);
      cVals[7][2] = extrapolation(data[ind[1][2][1]], data[ind[1][1][1]], data[ind[1][0][1]]);
      cVals[8][0] = extrapolation(data[ind[2][0][0]], data[ind[2][0][1]], data[ind[2][0][2]]);
      cVals[7][3] = extrapolation(data[ind[1][2][2]], data[ind[1][1][2]], data[ind[1][0][2]]);
      cVals[8][1] = data[ind[2][0][0]];
      cVals[8][2] = data[ind[2][0][1]];
      cVals[8][3] = data[ind[2][0][2]];
      cVals[9][0] = extrapolation(data[ind[2][1][0]], data[ind[2][1][1]], data[ind[2][1][2]]);
      cVals[9][1] = data[ind[2][1][0]];
      cVals[9][2] = data[ind[2][1][1]];
      cVals[9][3] = data[ind[2][1][2]];
      cVals[10][0] = extrapolation(data[ind[2][1][0]], data[ind[2][0][1]], data[ind[2][0][2] + deltaR[i0]]);
      cVals[10][1] = data[ind[2][2][0]];
      cVals[10][2] = data[ind[2][2][1]];
      cVals[10][3] = data[ind[2][2][2]];
      cVals[11][0] = extrapolation(data[ind[2][2][0]], data[ind[2][1][1]], data[ind[2][0][2]]);
      cVals[11][1] = extrapolation(data[ind[2][2][0]], data[ind[2][1][0]], data[ind[2][0][0]]);
      cVals[11][2] = extrapolation(data[ind[2][2][1]], data[ind[2][1][1]], data[ind[2][0][1]]);
      cVals[11][3] = extrapolation(data[ind[2][2][2]], data[ind[2][1][2]], data[ind[2][0][2]]);
      cVals[12][0] = extrapolation(data[ind[3][0][0]], data[ind[3][0][1]], data[ind[3][0][2]]);
      cVals[12][1] = data[ind[3][0][0]];
      cVals[12][2] = data[ind[3][0][1]];
      cVals[12][3] = data[ind[3][0][2]];
      cVals[13][0] = extrapolation(data[ind[3][1][0]], data[ind[3][1][1]], data[ind[3][1][2]]);
      cVals[13][1] = data[ind[3][1][0]];
      cVals[13][2] = data[ind[3][1][1]];
      cVals[13][3] = data[ind[3][1][2]];
      cVals[14][0] = extrapolation(data[ind[3][1][0]], data[ind[3][0][1]], data[ind[3][0][2] + deltaR[i0]]);
      cVals[14][1] = data[ind[3][2][0]];
      cVals[14][2] = data[ind[3][2][1]];
      cVals[14][3] = data[ind[3][2][2]];
      cVals[15][0] = extrapolation(data[ind[3][2][0]], data[ind[3][1][1]], data[ind[3][0][2]]);
      cVals[15][1] = extrapolation(data[ind[3][2][0]], data[ind[3][1][0]], data[ind[3][0][0]]);
      cVals[15][2] = extrapolation(data[ind[3][2][1]], data[ind[3][1][1]], data[ind[3][0][1]]);
      cVals[15][3] = extrapolation(data[ind[3][2][2]], data[ind[3][1][2]], data[ind[3][0][2]]);
    } break;

    case GridPos::Edge3:
    case GridPos::Edge7:
    case GridPos::LineL: {
      const int ind[4][3][3]{
        {{ii_x_y_z + deltaPhi[i0] + deltaR[i0] + deltaZ[i0], ind[0][0][0] - deltaZ[i0], ind[0][0][1] - deltaZ[i0]},
         {ind[0][0][0] - deltaR[i0], ind[0][1][0] - deltaZ[i0], ind[0][1][1] - deltaZ[i0]},
         {ind[0][1][0] - deltaR[i0], ind[0][2][0] - deltaZ[i0], ind[0][2][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaR[i0] + deltaZ[i0], ind[1][0][0] - deltaZ[i0], ind[1][0][1] - deltaZ[i0]},
         {ind[1][0][0] - deltaR[i0], ind[1][1][0] - deltaZ[i0], ind[1][1][1] - deltaZ[i0]},
         {ind[1][1][0] - deltaR[i0], ind[1][2][0] - deltaZ[i0], ind[1][2][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i1] + deltaR[i0] + deltaZ[i0], ind[2][0][0] - deltaZ[i0], ind[2][0][1] - deltaZ[i0]},
         {ind[2][0][0] - deltaR[i0], ind[2][1][0] - deltaZ[i0], ind[2][1][1] - deltaZ[i0]},
         {ind[2][1][0] - deltaR[i0], ind[2][2][0] - deltaZ[i0], ind[2][2][1] - deltaZ[i0]}},
        {{ii_x_y_z + deltaPhi[i2] + deltaR[i0] + deltaZ[i0], ind[3][0][0] - deltaZ[i0], ind[3][0][1] - deltaZ[i0]},
         {ind[3][0][0] - deltaR[i0], ind[3][1][0] - deltaZ[i0], ind[3][1][1] - deltaZ[i0]},
         {ind[3][1][0] - deltaR[i0], ind[3][2][0] - deltaZ[i0], ind[3][2][1] - deltaZ[i0]}}};

      cVals[0][0] = data[ind[0][0][0]];
      cVals[0][1] = data[ind[0][0][1]];
      cVals[0][2] = data[ind[0][0][2]];
      cVals[0][3] = extrapolation(data[ind[0][0][2]], data[ind[0][0][1]], data[ind[0][0][0]]);
      cVals[1][0] = data[ind[0][1][0]];
      cVals[1][1] = data[ind[0][1][1]];
      cVals[1][2] = data[ind[0][1][2]];
      cVals[1][3] = extrapolation(data[ind[0][1][2]], data[ind[0][1][1]], data[ind[0][1][0]]);
      cVals[2][0] = data[ind[0][2][0]];
      cVals[2][1] = data[ind[0][2][1]];
      cVals[2][2] = data[ind[0][2][2]];
      cVals[2][3] = extrapolation(data[ind[0][2][2]], data[ind[0][2][1]], data[ind[0][2][0]]);
      cVals[3][0] = extrapolation(data[ind[0][2][0]], data[ind[0][1][0]], data[ind[0][0][0]]);
      cVals[3][1] = extrapolation(data[ind[0][2][1]], data[ind[0][1][1]], data[ind[0][0][1]]);
      cVals[3][2] = extrapolation(data[ind[0][2][2]], data[ind[0][1][2]], data[ind[0][0][2]]);
      cVals[3][3] = extrapolation(data[ind[0][2][2]], data[ind[0][1][1]], data[ind[0][0][0]]);
      cVals[4][0] = data[ind[1][0][0]];
      cVals[4][1] = data[ind[1][0][1]];
      cVals[4][2] = data[ind[1][0][2]];
      cVals[4][3] = extrapolation(data[ind[1][0][2]], data[ind[1][0][1]], data[ind[1][0][0]]);
      cVals[5][0] = data[ind[1][1][0]];
      cVals[5][2] = data[ind[1][1][2]];
      cVals[5][3] = extrapolation(data[ind[1][1][2]], data[ii_x_y_z], data[ind[1][1][0]]);
      cVals[6][0] = data[ind[1][2][0]];
      cVals[6][1] = data[ind[1][2][1]];
      cVals[6][2] = data[ind[1][2][2]];
      cVals[6][3] = extrapolation(data[ind[1][2][2]], data[ind[1][2][1]], data[ind[1][2][0]]);
      cVals[7][0] = extrapolation(data[ind[1][2][0]], data[ind[1][1][0]], data[ind[1][0][0]]);
      cVals[7][1] = extrapolation(data[ind[1][2][1]], data[ii_x_y_z], data[ind[1][0][1]]);
      cVals[7][2] = extrapolation(data[ind[1][2][2]], data[ind[1][1][2]], data[ind[1][0][2]]);
      cVals[7][3] = extrapolation(data[ind[1][2][2]], data[ind[1][1][1]], data[ind[1][0][0]]);
      cVals[8][0] = data[ind[2][0][0]];
      cVals[8][1] = data[ind[2][0][1]];
      cVals[8][2] = data[ind[2][0][2]];
      cVals[8][3] = extrapolation(data[ind[2][0][2]], data[ind[2][0][1]], data[ind[2][0][0]]);
      cVals[9][0] = data[ind[2][1][0]];
      cVals[9][1] = data[ind[2][1][1]];
      cVals[9][2] = data[ind[2][1][2]];
      cVals[9][3] = extrapolation(data[ind[2][1][2]], data[ind[2][1][1]], data[ind[2][1][0]]);
      cVals[10][0] = data[ind[2][2][0]];
      cVals[10][1] = data[ind[2][2][1]];
      cVals[10][2] = data[ind[2][2][2]];
      cVals[10][3] = extrapolation(data[ind[2][2][2]], data[ind[2][2][1]], data[ind[2][2][0]]);
      cVals[11][0] = extrapolation(data[ind[2][2][0]], data[ind[2][1][0]], data[ind[2][0][0]]);
      cVals[11][1] = extrapolation(data[ind[2][2][1]], data[ind[2][1][1]], data[ind[2][0][1]]);
      cVals[11][2] = extrapolation(data[ind[2][2][2]], data[ind[2][1][2]], data[ind[2][0][2]]);
      cVals[11][3] = extrapolation(data[ind[2][2][2]], data[ind[2][1][1]], data[ind[2][0][0]]);
      cVals[12][0] = data[ind[3][0][0]];
      cVals[12][1] = data[ind[3][0][1]];
      cVals[12][2] = data[ind[3][0][2]];
      cVals[12][3] = extrapolation(data[ind[3][0][2]], data[ind[3][0][1]], data[ind[3][0][0]]);
      cVals[13][0] = data[ind[3][1][0]];
      cVals[13][1] = data[ind[3][1][1]];
      cVals[13][2] = data[ind[3][1][2]];
      cVals[13][3] = extrapolation(data[ind[3][1][2]], data[ind[3][1][1]], data[ind[3][1][0]]);
      cVals[14][0] = data[ind[3][2][0]];
      cVals[14][1] = data[ind[3][2][1]];
      cVals[14][2] = data[ind[3][2][2]];
      cVals[14][3] = extrapolation(data[ind[3][2][2]], data[ind[3][2][1]], data[ind[3][2][0]]);
      cVals[15][0] = extrapolation(data[ind[3][2][0]], data[ind[3][1][0]], data[ind[3][0][0]]);
      cVals[15][1] = extrapolation(data[ind[3][2][1]], data[ind[3][1][1]], data[ind[3][0][1]]);
      cVals[15][2] = extrapolation(data[ind[3][2][2]], data[ind[3][1][2]], data[ind[3][0][2]]);
      cVals[15][3] = extrapolation(data[ind[3][2][2]], data[ind[3][1][1]], data[ind[3][0][0]]);
    } break;
  }
}

template class o2::tpc::TriCubicInterpolator<double>;
template class o2::tpc::TriCubicInterpolator<float>;
