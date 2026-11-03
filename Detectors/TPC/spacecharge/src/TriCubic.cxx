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
#include "TPCSpaceCharge/Vector.h"
#include <algorithm>
#include <type_traits>

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
/// sum of four values with the same association of the terms as the horizontal sum of a Vc::Memory<Vc::Vector<DataT>, 4> (Vector<DataT, 4>),
/// which was used in the previous implementation (to get bit-identical results):
/// Vc scalar: ((x0 + x1) + x2) + x3, SSE float (4 lanes): (x0 + x3) + (x1 + x2), SSE double (2 x 2 lanes), AVX double (4 lanes), AVX float (8 lanes with zero padding): (x0 + x2) + (x1 + x3)
template <typename DataT>
inline DataT sum4(const DataT x0, const DataT x1, const DataT x2, const DataT x3)
{
  if constexpr (Vc::Vector<DataT>::Size == 1) {
    return ((x0 + x1) + x2) + x3;
  } else if constexpr (std::is_same_v<DataT, float> && (Vc::Vector<DataT>::Size == 4)) {
    return (x0 + x3) + (x1 + x2);
  } else {
    return (x0 + x2) + (x1 + x3);
  }
}

/// calculate the 4 weights of the 1D cubic (Catmull-Rom) interpolation for the relative positions t in the cell for all dimensions (the last one is not used)
/// \param t relative position in the cell for each dimension
/// \param weights weights[i][dim]
template <typename DataT>
inline void getWeights(const DataT (&t)[4], DataT (&weights)[4][4])
{
  static constexpr DataT matrixA[4][4]{{0, -0.5, 1, -0.5},
                                       {1, 0, -2.5, 1.5},
                                       {0, 0.5, 2., -1.5},
                                       {0, 0, -0.5, 0.5}};
  DataT t2[4];
  DataT t3[4];
  for (int dim = 0; dim < 4; ++dim) {
    t2[dim] = t[dim] * t[dim];
    t3[dim] = t[dim] * t2[dim];
  }
  for (int i = 0; i < 4; ++i) {
    for (int dim = 0; dim < 4; ++dim) {
      weights[i][dim] = sum4<DataT>(DataT(0) + matrixA[i][0] * DataT(1), DataT(0) + matrixA[i][1] * t[dim], DataT(0) + matrixA[i][2] * t2[dim], DataT(0) + matrixA[i][3] * t3[dim]);
    }
  }
}

/// sum_{slice, row, iz} weightsPhi[slice] * weightsR[row] * weightsZ[iz] * vals[sliceOffset[slice] + row * rowStride + iz]
/// \param vals values of the grid
/// \param sliceOffset offset of the 4 phi slices
/// \param rowStride offset between the rows in r direction (the 4 values in z direction are consecutive)
template <typename DataT>
inline DataT interpolateCell(const DataT* vals, const int (&sliceOffset)[4], const int rowStride, const DataT (&weightsZ)[4], const DataT (&weightsR)[4], const DataT (&weightsPhi)[4])
{
  using VDataT = Vc::Vector<DataT>;
  DataT result{};
  if constexpr (VDataT::Size <= 4) {
    // same operations as in the previous Vc implementation
    constexpr size_t NVec = 4 / VDataT::Size;
    VDataT vecWeightsZ[NVec];
    for (size_t j = 0; j < NVec; ++j) {
      vecWeightsZ[j] = VDataT(weightsZ + j * VDataT::Size, Vc::Unaligned);
    }
    for (int slice = 0; slice < 4; ++slice) {
      for (int row = 0; row < 4; ++row) {
        const DataT weightPhiR = weightsPhi[slice] * weightsR[row];
        const DataT* val = vals + sliceOffset[slice] + row * rowStride;
        VDataT sum = (weightPhiR * vecWeightsZ[0]) * VDataT(val, Vc::Unaligned);
        for (size_t j = 1; j < NVec; ++j) {
          sum += (weightPhiR * vecWeightsZ[j]) * VDataT(val + j * VDataT::Size, Vc::Unaligned);
        }
        result += sum.sum();
      }
    }
  } else {
    for (int slice = 0; slice < 4; ++slice) {
      for (int row = 0; row < 4; ++row) {
        const DataT weightPhiR = weightsPhi[slice] * weightsR[row];
        const DataT* val = vals + sliceOffset[slice] + row * rowStride;
        result += sum4<DataT>((weightPhiR * weightsZ[0]) * val[0], (weightPhiR * weightsZ[1]) * val[1], (weightPhiR * weightsZ[2]) * val[2], (weightPhiR * weightsZ[3]) * val[3]);
      }
    }
  }
  return result;
}
} // namespace

template <typename DataT>
DataT TriCubicInterpolator<DataT>::interpolateSparse(const DataT z, const DataT r, const DataT phi) const
{
  // check if data is empty
  if (!mGridData->getNDataPoints()) {
    return 0;
  }

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
  const int indexZ = static_cast<int>(posRel[FZ]);
  const int indexR = static_cast<int>(posRel[FR]);
  const int indexPhi = static_cast<int>(posRel[FPHI]);

  // weights of the interpolation in each dimension
  const DataT t[4]{posRelN[FZ] - indexZ, posRelN[FR] - indexR, posRelN[FPHI] - indexPhi, 0};
  DataT weights[4][4];
  getWeights<DataT>(t, weights);
  const DataT weightsZ[4]{weights[0][FZ], weights[1][FZ], weights[2][FZ], weights[3][FZ]};
  const DataT weightsR[4]{weights[0][FR], weights[1][FR], weights[2][FR], weights[3][FR]};
  const DataT weightsPhi[4]{weights[0][FPHI], weights[1][FPHI], weights[2][FPHI], weights[3][FPHI]};

  const GridPos location = findPos(indexZ, indexR, indexPhi);
  if ((location == GridPos::InnerVolume) || (location == GridPos::SideZRight) || (location == GridPos::SideZLeft)) {
    // all 64 values are inside the grid (the phi direction is periodic): directly use the values from the grid.
    // The values in z direction are consecutive in memory (DataContainer3D::getDataIndex()): start at (iz - 1, ir - 1, iphi)
    int deltaPhi[3]{};
    getDataIndexCircularArray(indexPhi, FPHI, deltaPhi);
    const int deltaPhiSlice[4]{deltaPhi[0], 0, deltaPhi[1], deltaPhi[2]};
    const int deltaR = mGridProperties->getDeltaDataIndex(1, FR);
    const DataT* data = mGridData->getData().data() + mGridData->getDataIndex(indexZ, indexR, indexPhi) - deltaR - 1;
    return interpolateCell<DataT>(data, deltaPhiSlice, deltaR, weightsZ, weightsR, weightsPhi);
  }

  // missing values outside of the grid are extrapolated
  DataT cVals[16][4];
  setValues(indexZ, indexR, indexPhi, location, cVals);
  static constexpr int sliceOffset[4]{0, 16, 32, 48};
  return interpolateCell<DataT>(&cVals[0][0], sliceOffset, 4, weightsZ, weightsR, weightsPhi);
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
void TriCubicInterpolator<DataT>::setValues(const int iz, const int ir, const int iphi, const GridPos location, DataT (&cVals)[16][4]) const
{
  const int ii_x_y_z = mGridData->getDataIndex(iz, ir, iphi);

  int deltaZ[3]{mGridProperties->getDeltaDataIndex(-1, 0), mGridProperties->getDeltaDataIndex(1, 0), mGridProperties->getDeltaDataIndex(2, 0)};
  int deltaR[3]{mGridProperties->getDeltaDataIndex(-1, 1), mGridProperties->getDeltaDataIndex(1, 1), mGridProperties->getDeltaDataIndex(2, 1)};
  int deltaPhi[3]{};
  getDataIndexCircularArray(iphi, FPHI, deltaPhi);

  const int i0 = 0;
  const int i1 = 1;
  const int i2 = 2;

  cVals[5][1] = (*mGridData)[ii_x_y_z];

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

      cVals[0][0] = (*mGridData)[ind[0][0][0]];
      cVals[0][1] = (*mGridData)[ind[0][0][1]];
      cVals[0][2] = (*mGridData)[ind[0][0][2]];
      cVals[0][3] = (*mGridData)[ind[0][0][3]];
      cVals[1][0] = (*mGridData)[ind[0][1][0]];
      cVals[1][1] = (*mGridData)[ind[0][1][1]];
      cVals[1][2] = (*mGridData)[ind[0][1][2]];
      cVals[1][3] = (*mGridData)[ind[0][1][3]];
      cVals[2][0] = (*mGridData)[ind[0][2][0]];
      cVals[2][1] = (*mGridData)[ind[0][2][1]];
      cVals[2][2] = (*mGridData)[ind[0][2][2]];
      cVals[2][3] = (*mGridData)[ind[0][2][3]];
      cVals[3][0] = (*mGridData)[ind[0][3][0]];
      cVals[3][1] = (*mGridData)[ind[0][3][1]];
      cVals[3][2] = (*mGridData)[ind[0][3][2]];
      cVals[3][3] = (*mGridData)[ind[0][3][3]];
      cVals[4][0] = (*mGridData)[ind[1][0][0]];
      cVals[4][1] = (*mGridData)[ind[1][0][1]];
      cVals[4][2] = (*mGridData)[ind[1][0][2]];
      cVals[4][3] = (*mGridData)[ind[1][0][3]];
      cVals[5][2] = (*mGridData)[ind[1][1][2]];
      cVals[5][0] = (*mGridData)[ind[1][1][0]];
      cVals[5][3] = (*mGridData)[ind[1][1][3]];
      cVals[6][0] = (*mGridData)[ind[1][2][0]];
      cVals[6][1] = (*mGridData)[ind[1][2][1]];
      cVals[6][2] = (*mGridData)[ind[1][2][2]];
      cVals[6][3] = (*mGridData)[ind[1][2][3]];
      cVals[7][0] = (*mGridData)[ind[1][3][0]];
      cVals[7][1] = (*mGridData)[ind[1][3][1]];
      cVals[7][2] = (*mGridData)[ind[1][3][2]];
      cVals[7][3] = (*mGridData)[ind[1][3][3]];
      cVals[8][0] = (*mGridData)[ind[2][0][0]];
      cVals[8][1] = (*mGridData)[ind[2][0][1]];
      cVals[8][2] = (*mGridData)[ind[2][0][2]];
      cVals[8][3] = (*mGridData)[ind[2][0][3]];
      cVals[9][0] = (*mGridData)[ind[2][1][0]];
      cVals[9][1] = (*mGridData)[ind[2][1][1]];
      cVals[9][2] = (*mGridData)[ind[2][1][2]];
      cVals[9][3] = (*mGridData)[ind[2][1][3]];
      cVals[10][0] = (*mGridData)[ind[2][2][0]];
      cVals[10][1] = (*mGridData)[ind[2][2][1]];
      cVals[10][2] = (*mGridData)[ind[2][2][2]];
      cVals[10][3] = (*mGridData)[ind[2][2][3]];
      cVals[11][0] = (*mGridData)[ind[2][3][0]];
      cVals[11][1] = (*mGridData)[ind[2][3][1]];
      cVals[11][2] = (*mGridData)[ind[2][3][2]];
      cVals[11][3] = (*mGridData)[ind[2][3][3]];
      cVals[12][0] = (*mGridData)[ind[3][0][0]];
      cVals[12][1] = (*mGridData)[ind[3][0][1]];
      cVals[12][2] = (*mGridData)[ind[3][0][2]];
      cVals[12][3] = (*mGridData)[ind[3][0][3]];
      cVals[13][0] = (*mGridData)[ind[3][1][0]];
      cVals[13][1] = (*mGridData)[ind[3][1][1]];
      cVals[13][2] = (*mGridData)[ind[3][1][2]];
      cVals[13][3] = (*mGridData)[ind[3][1][3]];
      cVals[14][0] = (*mGridData)[ind[3][2][0]];
      cVals[14][1] = (*mGridData)[ind[3][2][1]];
      cVals[14][2] = (*mGridData)[ind[3][2][2]];
      cVals[14][3] = (*mGridData)[ind[3][2][3]];
      cVals[15][0] = (*mGridData)[ind[3][3][0]];
      cVals[15][1] = (*mGridData)[ind[3][3][1]];
      cVals[15][2] = (*mGridData)[ind[3][3][2]];
      cVals[15][3] = (*mGridData)[ind[3][3][3]];
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

      cVals[0][0] = (*mGridData)[ind[0][0][0]];
      cVals[0][1] = (*mGridData)[ind[0][0][1]];
      cVals[0][2] = (*mGridData)[ind[0][0][2]];
      cVals[0][3] = extrapolation((*mGridData)[ind[0][0][2]], (*mGridData)[ind[0][0][1]], (*mGridData)[ind[0][0][0]]);
      cVals[1][0] = (*mGridData)[ind[0][1][0]];
      cVals[1][1] = (*mGridData)[ind[0][1][1]];
      cVals[1][2] = (*mGridData)[ind[0][1][2]];
      cVals[1][3] = extrapolation((*mGridData)[ind[0][1][2]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][1][0]]);
      cVals[2][0] = (*mGridData)[ind[0][2][0]];
      cVals[2][1] = (*mGridData)[ind[0][2][1]];
      cVals[2][2] = (*mGridData)[ind[0][2][2]];
      cVals[2][3] = extrapolation((*mGridData)[ind[0][2][2]], (*mGridData)[ind[0][2][1]], (*mGridData)[ind[0][2][0]]);
      cVals[3][0] = (*mGridData)[ind[0][3][0]];
      cVals[3][1] = (*mGridData)[ind[0][3][1]];
      cVals[3][2] = (*mGridData)[ind[0][3][2]];
      cVals[3][3] = extrapolation((*mGridData)[ind[0][3][2]], (*mGridData)[ind[0][3][1]], (*mGridData)[ind[0][3][0]]);
      cVals[4][0] = (*mGridData)[ind[1][0][0]];
      cVals[4][1] = (*mGridData)[ind[1][0][1]];
      cVals[4][2] = (*mGridData)[ind[1][0][2]];
      cVals[4][3] = extrapolation((*mGridData)[ind[1][0][2]], (*mGridData)[ind[1][0][1]], (*mGridData)[ind[1][0][0]]);
      cVals[5][0] = (*mGridData)[ind[1][1][0]];
      cVals[5][2] = (*mGridData)[ind[1][1][2]];
      cVals[5][3] = extrapolation((*mGridData)[ind[1][1][2]], (*mGridData)[ii_x_y_z], (*mGridData)[ind[1][1][0]]);
      cVals[6][0] = (*mGridData)[ind[1][2][0]];
      cVals[6][1] = (*mGridData)[ind[1][2][1]];
      cVals[6][2] = (*mGridData)[ind[1][2][2]];
      cVals[6][3] = extrapolation((*mGridData)[ind[1][2][2]], (*mGridData)[ind[1][2][1]], (*mGridData)[ind[1][2][0]]);
      cVals[7][0] = (*mGridData)[ind[1][3][0]];
      cVals[7][1] = (*mGridData)[ind[1][3][1]];
      cVals[7][2] = (*mGridData)[ind[1][3][2]];
      cVals[7][3] = extrapolation((*mGridData)[ind[1][3][2]], (*mGridData)[ind[1][3][1]], (*mGridData)[ind[1][3][0]]);
      cVals[8][0] = (*mGridData)[ind[2][0][0]];
      cVals[8][1] = (*mGridData)[ind[2][0][1]];
      cVals[8][2] = (*mGridData)[ind[2][0][2]];
      cVals[8][3] = extrapolation((*mGridData)[ind[2][0][2]], (*mGridData)[ind[2][0][1]], (*mGridData)[ind[2][0][0]]);
      cVals[9][0] = (*mGridData)[ind[2][1][0]];
      cVals[9][1] = (*mGridData)[ind[2][1][1]];
      cVals[9][2] = (*mGridData)[ind[2][1][2]];
      cVals[9][3] = extrapolation((*mGridData)[ind[2][1][2]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][1][0]]);
      cVals[10][0] = (*mGridData)[ind[2][2][0]];
      cVals[10][1] = (*mGridData)[ind[2][2][1]];
      cVals[10][2] = (*mGridData)[ind[2][2][2]];
      cVals[10][3] = extrapolation((*mGridData)[ind[2][2][2]], (*mGridData)[ind[2][2][1]], (*mGridData)[ind[2][2][0]]);
      cVals[11][0] = (*mGridData)[ind[2][3][0]];
      cVals[11][1] = (*mGridData)[ind[2][3][1]];
      cVals[11][2] = (*mGridData)[ind[2][3][2]];
      cVals[11][3] = extrapolation((*mGridData)[ind[2][3][2]], (*mGridData)[ind[2][3][1]], (*mGridData)[ind[2][3][0]]);
      cVals[12][0] = (*mGridData)[ind[3][0][0]];
      cVals[12][1] = (*mGridData)[ind[3][0][1]];
      cVals[12][2] = (*mGridData)[ind[3][0][2]];
      cVals[13][0] = (*mGridData)[ind[3][1][0]];
      cVals[12][3] = extrapolation((*mGridData)[ind[3][0][2]], (*mGridData)[ind[3][0][1]], (*mGridData)[ind[3][0][0]]);
      cVals[13][1] = (*mGridData)[ind[3][1][1]];
      cVals[13][2] = (*mGridData)[ind[3][1][2]];
      cVals[13][3] = extrapolation((*mGridData)[ind[3][1][2]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][1][0]]);
      cVals[14][0] = (*mGridData)[ind[3][2][0]];
      cVals[14][1] = (*mGridData)[ind[3][2][1]];
      cVals[14][2] = (*mGridData)[ind[3][2][2]];
      cVals[14][3] = extrapolation((*mGridData)[ind[3][2][2]], (*mGridData)[ind[3][2][1]], (*mGridData)[ind[3][2][0]]);
      cVals[15][0] = (*mGridData)[ind[3][3][0]];
      cVals[15][1] = (*mGridData)[ind[3][3][1]];
      cVals[15][2] = (*mGridData)[ind[3][3][2]];
      cVals[15][3] = extrapolation((*mGridData)[ind[3][3][2]], (*mGridData)[ind[3][3][1]], (*mGridData)[ind[3][3][0]]);
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

      cVals[0][0] = (*mGridData)[ind[0][0][0]];
      cVals[0][1] = (*mGridData)[ind[0][0][1]];
      cVals[0][2] = (*mGridData)[ind[0][0][2]];
      cVals[0][3] = (*mGridData)[ind[0][0][3]];
      cVals[1][0] = (*mGridData)[ind[0][1][0]];
      cVals[1][1] = (*mGridData)[ind[0][1][1]];
      cVals[1][2] = (*mGridData)[ind[0][1][2]];
      cVals[1][3] = (*mGridData)[ind[0][1][3]];
      cVals[2][0] = (*mGridData)[ind[0][2][0]];
      cVals[2][1] = (*mGridData)[ind[0][2][1]];
      cVals[2][2] = (*mGridData)[ind[0][2][2]];
      cVals[2][3] = (*mGridData)[ind[0][2][3]];
      cVals[3][0] = extrapolation((*mGridData)[ind[0][2][0]], (*mGridData)[ind[0][1][0]], (*mGridData)[ind[0][0][0]]);
      cVals[3][1] = extrapolation((*mGridData)[ind[0][2][1]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][0][1]]);
      cVals[3][2] = extrapolation((*mGridData)[ind[0][2][2]], (*mGridData)[ind[0][1][2]], (*mGridData)[ind[0][0][2]]);
      cVals[3][3] = extrapolation((*mGridData)[ind[0][2][3]], (*mGridData)[ind[0][1][3]], (*mGridData)[ind[0][0][3]]);
      cVals[4][0] = (*mGridData)[ind[1][0][0]];
      cVals[4][1] = (*mGridData)[ind[1][0][1]];
      cVals[4][2] = (*mGridData)[ind[1][0][2]];
      cVals[4][3] = (*mGridData)[ind[1][0][3]];
      cVals[5][0] = (*mGridData)[ind[1][1][0]];
      cVals[5][2] = (*mGridData)[ind[1][1][2]];
      cVals[5][3] = (*mGridData)[ind[1][1][3]];
      cVals[6][0] = (*mGridData)[ind[1][2][0]];
      cVals[6][1] = (*mGridData)[ind[1][2][1]];
      cVals[6][2] = (*mGridData)[ind[1][2][2]];
      cVals[6][3] = (*mGridData)[ind[1][2][3]];
      cVals[7][0] = extrapolation((*mGridData)[ind[1][2][0]], (*mGridData)[ind[1][1][0]], (*mGridData)[ind[1][0][0]]);
      cVals[7][1] = extrapolation((*mGridData)[ind[1][2][1]], (*mGridData)[ii_x_y_z], (*mGridData)[ind[1][0][1]]);
      cVals[7][2] = extrapolation((*mGridData)[ind[1][2][2]], (*mGridData)[ind[1][1][2]], (*mGridData)[ind[1][0][2]]);
      cVals[7][3] = extrapolation((*mGridData)[ind[1][2][3]], (*mGridData)[ind[1][1][3]], (*mGridData)[ind[1][0][3]]);
      cVals[8][0] = (*mGridData)[ind[2][0][0]];
      cVals[8][1] = (*mGridData)[ind[2][0][1]];
      cVals[8][2] = (*mGridData)[ind[2][0][2]];
      cVals[8][3] = (*mGridData)[ind[2][0][3]];
      cVals[9][0] = (*mGridData)[ind[2][1][0]];
      cVals[9][1] = (*mGridData)[ind[2][1][1]];
      cVals[9][2] = (*mGridData)[ind[2][1][2]];
      cVals[9][3] = (*mGridData)[ind[2][1][3]];
      cVals[10][0] = (*mGridData)[ind[2][2][0]];
      cVals[10][1] = (*mGridData)[ind[2][2][1]];
      cVals[10][2] = (*mGridData)[ind[2][2][2]];
      cVals[10][3] = (*mGridData)[ind[2][2][3]];
      cVals[11][0] = extrapolation((*mGridData)[ind[2][2][0]], (*mGridData)[ind[2][1][0]], (*mGridData)[ind[2][0][0]]);
      cVals[11][1] = extrapolation((*mGridData)[ind[2][2][1]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][0][1]]);
      cVals[11][2] = extrapolation((*mGridData)[ind[2][2][2]], (*mGridData)[ind[2][1][2]], (*mGridData)[ind[2][0][2]]);
      cVals[11][3] = extrapolation((*mGridData)[ind[2][2][3]], (*mGridData)[ind[2][1][3]], (*mGridData)[ind[2][0][3]]);
      cVals[12][0] = (*mGridData)[ind[3][0][0]];
      cVals[12][1] = (*mGridData)[ind[3][0][1]];
      cVals[12][2] = (*mGridData)[ind[3][0][2]];
      cVals[12][3] = (*mGridData)[ind[3][0][3]];
      cVals[13][0] = (*mGridData)[ind[3][1][0]];
      cVals[13][1] = (*mGridData)[ind[3][1][1]];
      cVals[13][2] = (*mGridData)[ind[3][1][2]];
      cVals[13][3] = (*mGridData)[ind[3][1][3]];
      cVals[14][0] = (*mGridData)[ind[3][2][0]];
      cVals[14][1] = (*mGridData)[ind[3][2][1]];
      cVals[14][2] = (*mGridData)[ind[3][2][2]];
      cVals[14][3] = (*mGridData)[ind[3][2][3]];
      cVals[15][0] = extrapolation((*mGridData)[ind[3][2][0]], (*mGridData)[ind[3][1][0]], (*mGridData)[ind[3][0][0]]);
      cVals[15][1] = extrapolation((*mGridData)[ind[3][2][1]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][0][1]]);
      cVals[15][2] = extrapolation((*mGridData)[ind[3][2][2]], (*mGridData)[ind[3][1][2]], (*mGridData)[ind[3][0][2]]);
      cVals[15][3] = extrapolation((*mGridData)[ind[3][2][3]], (*mGridData)[ind[3][1][3]], (*mGridData)[ind[3][0][3]]);
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

      cVals[0][0] = extrapolation((*mGridData)[ind[0][0][0]], (*mGridData)[ind[0][1][0]], (*mGridData)[ind[0][2][0]]);
      cVals[0][1] = extrapolation((*mGridData)[ind[0][0][1]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][2][1]]);
      cVals[0][2] = extrapolation((*mGridData)[ind[0][0][2]], (*mGridData)[ind[0][1][2]], (*mGridData)[ind[0][2][2]]);
      cVals[0][3] = extrapolation((*mGridData)[ind[0][0][3]], (*mGridData)[ind[0][1][3]], (*mGridData)[ind[0][2][3]]);
      cVals[1][0] = (*mGridData)[ind[0][0][0]];
      cVals[1][1] = (*mGridData)[ind[0][0][1]];
      cVals[1][2] = (*mGridData)[ind[0][0][2]];
      cVals[1][3] = (*mGridData)[ind[0][0][3]];
      cVals[2][0] = (*mGridData)[ind[0][1][0]];
      cVals[2][1] = (*mGridData)[ind[0][1][1]];
      cVals[2][2] = (*mGridData)[ind[0][1][2]];
      cVals[2][3] = (*mGridData)[ind[0][1][3]];
      cVals[3][0] = (*mGridData)[ind[0][2][0]];
      cVals[3][1] = (*mGridData)[ind[0][2][1]];
      cVals[3][2] = (*mGridData)[ind[0][2][2]];
      cVals[3][3] = (*mGridData)[ind[0][2][3]];
      cVals[4][0] = extrapolation((*mGridData)[ind[1][0][0]], (*mGridData)[ind[1][1][0]], (*mGridData)[ind[1][2][0]]);
      cVals[4][1] = extrapolation((*mGridData)[ii_x_y_z], (*mGridData)[ind[1][1][1]], (*mGridData)[ind[1][2][1]]);
      cVals[4][2] = extrapolation((*mGridData)[ind[1][0][2]], (*mGridData)[ind[1][1][2]], (*mGridData)[ind[1][2][2]]);
      cVals[4][3] = extrapolation((*mGridData)[ind[1][0][3]], (*mGridData)[ind[1][1][3]], (*mGridData)[ind[1][2][3]]);
      cVals[5][0] = (*mGridData)[ind[1][0][0]];
      cVals[5][2] = (*mGridData)[ind[1][0][2]];
      cVals[5][3] = (*mGridData)[ind[1][0][3]];
      cVals[6][0] = (*mGridData)[ind[1][1][0]];
      cVals[6][1] = (*mGridData)[ind[1][1][1]];
      cVals[6][2] = (*mGridData)[ind[1][1][2]];
      cVals[6][3] = (*mGridData)[ind[1][1][3]];
      cVals[7][0] = (*mGridData)[ind[1][2][0]];
      cVals[7][1] = (*mGridData)[ind[1][2][1]];
      cVals[7][2] = (*mGridData)[ind[1][2][2]];
      cVals[7][3] = (*mGridData)[ind[1][2][3]];
      cVals[8][0] = extrapolation((*mGridData)[ind[2][0][0]], (*mGridData)[ind[2][1][0]], (*mGridData)[ind[2][2][0]]);
      cVals[8][1] = extrapolation((*mGridData)[ind[2][0][1]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][2][1]]);
      cVals[8][2] = extrapolation((*mGridData)[ind[2][0][2]], (*mGridData)[ind[2][1][2]], (*mGridData)[ind[2][2][2]]);
      cVals[8][3] = extrapolation((*mGridData)[ind[2][0][3]], (*mGridData)[ind[2][1][3]], (*mGridData)[ind[2][2][3]]);
      cVals[9][0] = (*mGridData)[ind[2][0][0]];
      cVals[9][1] = (*mGridData)[ind[2][0][1]];
      cVals[9][2] = (*mGridData)[ind[2][0][2]];
      cVals[9][3] = (*mGridData)[ind[2][0][3]];
      cVals[10][0] = (*mGridData)[ind[2][1][0]];
      cVals[10][1] = (*mGridData)[ind[2][1][1]];
      cVals[10][2] = (*mGridData)[ind[2][1][2]];
      cVals[10][3] = (*mGridData)[ind[2][1][3]];
      cVals[11][0] = (*mGridData)[ind[2][2][0]];
      cVals[11][1] = (*mGridData)[ind[2][2][1]];
      cVals[11][2] = (*mGridData)[ind[2][2][2]];
      cVals[11][3] = (*mGridData)[ind[2][2][3]];
      cVals[12][0] = extrapolation((*mGridData)[ind[3][0][0]], (*mGridData)[ind[3][1][0]], (*mGridData)[ind[3][2][0]]);
      cVals[12][1] = extrapolation((*mGridData)[ind[3][0][1]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][2][1]]);
      cVals[12][2] = extrapolation((*mGridData)[ind[3][0][2]], (*mGridData)[ind[3][1][2]], (*mGridData)[ind[3][2][2]]);
      cVals[12][3] = extrapolation((*mGridData)[ind[3][0][3]], (*mGridData)[ind[3][1][3]], (*mGridData)[ind[3][2][3]]);
      cVals[13][0] = (*mGridData)[ind[3][0][0]];
      cVals[13][1] = (*mGridData)[ind[3][0][1]];
      cVals[13][2] = (*mGridData)[ind[3][0][2]];
      cVals[13][3] = (*mGridData)[ind[3][0][3]];
      cVals[14][0] = (*mGridData)[ind[3][1][0]];
      cVals[14][1] = (*mGridData)[ind[3][1][1]];
      cVals[14][2] = (*mGridData)[ind[3][1][2]];
      cVals[14][3] = (*mGridData)[ind[3][1][3]];
      cVals[15][0] = (*mGridData)[ind[3][2][0]];
      cVals[15][1] = (*mGridData)[ind[3][2][1]];
      cVals[15][2] = (*mGridData)[ind[3][2][2]];
      cVals[15][3] = (*mGridData)[ind[3][2][3]];
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

      cVals[0][0] = extrapolation((*mGridData)[ind[0][0][0]], (*mGridData)[ind[0][0][1]], (*mGridData)[ind[0][0][2]]);
      cVals[0][1] = (*mGridData)[ind[0][0][0]];
      cVals[0][2] = (*mGridData)[ind[0][0][1]];
      cVals[0][3] = (*mGridData)[ind[0][0][2]];
      cVals[1][0] = extrapolation((*mGridData)[ind[0][1][0]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][1][2]]);
      cVals[1][1] = (*mGridData)[ind[0][1][0]];
      cVals[1][2] = (*mGridData)[ind[0][1][1]];
      cVals[1][3] = (*mGridData)[ind[0][1][2]];
      cVals[2][0] = extrapolation((*mGridData)[ind[0][2][0]], (*mGridData)[ind[0][2][1]], (*mGridData)[ind[0][2][2]]);
      cVals[2][1] = (*mGridData)[ind[0][2][0]];
      cVals[2][2] = (*mGridData)[ind[0][2][1]];
      cVals[2][3] = (*mGridData)[ind[0][2][2]];
      cVals[3][0] = extrapolation((*mGridData)[ind[0][3][0]], (*mGridData)[ind[0][3][1]], (*mGridData)[ind[0][3][2]]);
      cVals[3][1] = (*mGridData)[ind[0][3][0]];
      cVals[3][2] = (*mGridData)[ind[0][3][1]];
      cVals[3][3] = (*mGridData)[ind[0][3][2]];
      cVals[4][0] = extrapolation((*mGridData)[ind[1][0][0]], (*mGridData)[ind[1][0][1]], (*mGridData)[ind[1][0][2]]);
      cVals[4][1] = (*mGridData)[ind[1][0][0]];
      cVals[4][2] = (*mGridData)[ind[1][0][1]];
      cVals[4][3] = (*mGridData)[ind[1][0][2]];
      cVals[5][0] = extrapolation((*mGridData)[ii_x_y_z], (*mGridData)[ind[1][1][1]], (*mGridData)[ind[1][1][2]]);
      cVals[5][2] = (*mGridData)[ind[1][1][1]];
      cVals[5][3] = (*mGridData)[ind[1][1][2]];
      cVals[6][0] = extrapolation((*mGridData)[ind[1][2][0]], (*mGridData)[ind[1][2][1]], (*mGridData)[ind[1][2][2]]);
      cVals[6][1] = (*mGridData)[ind[1][2][0]];
      cVals[6][2] = (*mGridData)[ind[1][2][1]];
      cVals[6][3] = (*mGridData)[ind[1][2][2]];
      cVals[7][0] = extrapolation((*mGridData)[ind[1][3][0]], (*mGridData)[ind[1][3][1]], (*mGridData)[ind[1][3][2]]);
      cVals[7][1] = (*mGridData)[ind[1][3][0]];
      cVals[7][2] = (*mGridData)[ind[1][3][1]];
      cVals[7][3] = (*mGridData)[ind[1][3][2]];
      cVals[8][0] = extrapolation((*mGridData)[ind[2][0][0]], (*mGridData)[ind[2][0][1]], (*mGridData)[ind[2][0][2]]);
      cVals[8][1] = (*mGridData)[ind[2][0][0]];
      cVals[8][2] = (*mGridData)[ind[2][0][1]];
      cVals[8][3] = (*mGridData)[ind[2][0][2]];
      cVals[9][0] = extrapolation((*mGridData)[ind[2][1][0]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][1][2]]);
      cVals[9][1] = (*mGridData)[ind[2][1][0]];
      cVals[9][2] = (*mGridData)[ind[2][1][1]];
      cVals[9][3] = (*mGridData)[ind[2][1][2]];
      cVals[10][0] = extrapolation((*mGridData)[ind[2][2][0]], (*mGridData)[ind[2][2][1]], (*mGridData)[ind[2][2][2]]);
      cVals[10][1] = (*mGridData)[ind[2][2][0]];
      cVals[10][2] = (*mGridData)[ind[2][2][1]];
      cVals[10][3] = (*mGridData)[ind[2][2][2]];
      cVals[11][0] = extrapolation((*mGridData)[ind[2][3][0]], (*mGridData)[ind[2][3][1]], (*mGridData)[ind[2][3][2]]);
      cVals[11][1] = (*mGridData)[ind[2][3][0]];
      cVals[11][2] = (*mGridData)[ind[2][3][1]];
      cVals[11][3] = (*mGridData)[ind[2][3][2]];
      cVals[12][0] = extrapolation((*mGridData)[ind[3][0][0]], (*mGridData)[ind[3][0][1]], (*mGridData)[ind[3][0][2]]);
      cVals[12][1] = (*mGridData)[ind[3][0][0]];
      cVals[12][2] = (*mGridData)[ind[3][0][1]];
      cVals[12][3] = (*mGridData)[ind[3][0][2]];
      cVals[13][0] = extrapolation((*mGridData)[ind[3][1][0]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][1][2]]);
      cVals[13][1] = (*mGridData)[ind[3][1][0]];
      cVals[13][2] = (*mGridData)[ind[3][1][1]];
      cVals[13][3] = (*mGridData)[ind[3][1][2]];
      cVals[14][0] = extrapolation((*mGridData)[ind[3][2][0]], (*mGridData)[ind[3][2][1]], (*mGridData)[ind[3][2][2]]);
      cVals[14][1] = (*mGridData)[ind[3][2][0]];
      cVals[14][2] = (*mGridData)[ind[3][2][1]];
      cVals[14][3] = (*mGridData)[ind[3][2][2]];
      cVals[15][0] = extrapolation((*mGridData)[ind[3][3][0]], (*mGridData)[ind[3][3][1]], (*mGridData)[ind[3][3][2]]);
      cVals[15][1] = (*mGridData)[ind[3][3][0]];
      cVals[15][2] = (*mGridData)[ind[3][3][1]];
      cVals[15][3] = (*mGridData)[ind[3][3][2]];
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

      cVals[0][0] = extrapolation((*mGridData)[ind[0][0][0]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][2][2]]);
      cVals[0][1] = extrapolation((*mGridData)[ind[0][0][0]], (*mGridData)[ind[0][1][0]], (*mGridData)[ind[0][2][0]]);
      cVals[0][2] = extrapolation((*mGridData)[ind[0][0][1]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][2][1]]);
      cVals[0][3] = extrapolation((*mGridData)[ind[0][0][2]], (*mGridData)[ind[0][1][2]], (*mGridData)[ind[0][2][2]]);
      cVals[1][0] = extrapolation((*mGridData)[ind[0][0][0]], (*mGridData)[ind[0][0][1]], (*mGridData)[ind[0][0][2]]);
      cVals[1][1] = (*mGridData)[ind[0][0][0]];
      cVals[1][2] = (*mGridData)[ind[0][0][1]];
      cVals[1][3] = (*mGridData)[ind[0][0][2]];
      cVals[2][0] = extrapolation((*mGridData)[ind[0][1][0]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][1][2]]);
      cVals[2][1] = (*mGridData)[ind[0][1][0]];
      cVals[2][2] = (*mGridData)[ind[0][1][1]];
      cVals[2][3] = (*mGridData)[ind[0][1][2]];
      cVals[3][0] = extrapolation((*mGridData)[ind[0][2][0]], (*mGridData)[ind[0][2][1]], (*mGridData)[ind[0][2][2]]);
      cVals[3][1] = (*mGridData)[ind[0][2][0]];
      cVals[3][2] = (*mGridData)[ind[0][2][1]];
      cVals[3][3] = (*mGridData)[ind[0][2][2]];
      cVals[4][0] = extrapolation((*mGridData)[ii_x_y_z], (*mGridData)[ind[1][1][1]], (*mGridData)[ind[1][2][2]]);
      cVals[4][1] = extrapolation((*mGridData)[ii_x_y_z], (*mGridData)[ind[1][1][0]], (*mGridData)[ind[1][2][0]]);
      cVals[4][2] = extrapolation((*mGridData)[ind[1][0][1]], (*mGridData)[ind[1][1][1]], (*mGridData)[ind[1][2][1]]);
      cVals[4][3] = extrapolation((*mGridData)[ind[1][0][2]], (*mGridData)[ind[1][1][2]], (*mGridData)[ind[1][2][2]]);
      cVals[5][0] = extrapolation((*mGridData)[ii_x_y_z], (*mGridData)[ind[1][0][1]], (*mGridData)[ind[1][0][2]]);
      cVals[5][2] = (*mGridData)[ind[1][0][1]];
      cVals[5][3] = (*mGridData)[ind[1][0][2]];
      cVals[6][0] = extrapolation((*mGridData)[ind[1][1][0]], (*mGridData)[ind[1][1][1]], (*mGridData)[ind[1][1][2]]);
      cVals[6][1] = (*mGridData)[ind[1][1][0]];
      cVals[6][2] = (*mGridData)[ind[1][1][1]];
      cVals[6][3] = (*mGridData)[ind[1][1][2]];
      cVals[7][0] = extrapolation((*mGridData)[ind[1][2][0]], (*mGridData)[ind[1][2][1]], (*mGridData)[ind[1][2][2]]);
      cVals[7][1] = (*mGridData)[ind[1][2][0]];
      cVals[7][2] = (*mGridData)[ind[1][2][1]];
      cVals[7][3] = (*mGridData)[ind[1][2][2]];
      cVals[8][0] = extrapolation((*mGridData)[ind[2][0][0]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][2][2]]);
      cVals[8][1] = extrapolation((*mGridData)[ind[2][0][0]], (*mGridData)[ind[2][1][0]], (*mGridData)[ind[2][2][0]]);
      cVals[8][2] = extrapolation((*mGridData)[ind[2][0][1]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][2][1]]);
      cVals[8][3] = extrapolation((*mGridData)[ind[2][0][2]], (*mGridData)[ind[2][1][2]], (*mGridData)[ind[2][2][2]]);
      cVals[9][0] = extrapolation((*mGridData)[ind[2][0][0]], (*mGridData)[ind[2][0][1]], (*mGridData)[ind[2][0][2]]);
      cVals[9][1] = (*mGridData)[ind[2][0][0]];
      cVals[9][2] = (*mGridData)[ind[2][0][1]];
      cVals[9][3] = (*mGridData)[ind[2][0][2]];
      cVals[10][0] = extrapolation((*mGridData)[ind[2][1][0]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][1][2]]);
      cVals[10][1] = (*mGridData)[ind[2][1][0]];
      cVals[10][2] = (*mGridData)[ind[2][1][1]];
      cVals[10][3] = (*mGridData)[ind[2][1][2]];
      cVals[11][0] = extrapolation((*mGridData)[ind[2][2][0]], (*mGridData)[ind[2][2][1]], (*mGridData)[ind[2][2][2]]);
      cVals[11][1] = (*mGridData)[ind[2][2][0]];
      cVals[11][2] = (*mGridData)[ind[2][2][1]];
      cVals[11][3] = (*mGridData)[ind[2][2][2]];
      cVals[12][0] = extrapolation((*mGridData)[ind[3][0][0]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][2][2]]);
      cVals[12][1] = extrapolation((*mGridData)[ind[3][0][0]], (*mGridData)[ind[3][1][0]], (*mGridData)[ind[3][2][0]]);
      cVals[12][2] = extrapolation((*mGridData)[ind[3][0][1]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][2][1]]);
      cVals[12][3] = extrapolation((*mGridData)[ind[3][0][2]], (*mGridData)[ind[3][1][2]], (*mGridData)[ind[3][2][2]]);
      cVals[13][0] = extrapolation((*mGridData)[ind[3][0][0]], (*mGridData)[ind[3][0][1]], (*mGridData)[ind[3][0][2]]);
      cVals[13][1] = (*mGridData)[ind[3][0][0]];
      cVals[13][2] = (*mGridData)[ind[3][0][1]];
      cVals[13][3] = (*mGridData)[ind[3][0][2]];
      cVals[14][0] = extrapolation((*mGridData)[ind[3][1][0]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][1][2]]);
      cVals[14][1] = (*mGridData)[ind[3][1][0]];
      cVals[14][2] = (*mGridData)[ind[3][1][1]];
      cVals[14][3] = (*mGridData)[ind[3][1][2]];
      cVals[15][0] = extrapolation((*mGridData)[ind[3][2][0]], (*mGridData)[ind[3][2][1]], (*mGridData)[ind[3][2][2]]);
      cVals[15][1] = (*mGridData)[ind[3][2][0]];
      cVals[15][2] = (*mGridData)[ind[3][2][1]];
      cVals[15][3] = (*mGridData)[ind[3][2][2]];
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

      cVals[0][0] = extrapolation((*mGridData)[ind[0][0][0]], (*mGridData)[ind[0][1][0]], (*mGridData)[ind[0][2][0]]);
      cVals[0][1] = extrapolation((*mGridData)[ind[0][0][1]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][2][1]]);
      cVals[0][2] = extrapolation((*mGridData)[ind[0][0][1]], (*mGridData)[ind[0][1][0]], (*mGridData)[ind[0][2][0] + deltaZ[i0]]);
      cVals[0][3] = extrapolation((*mGridData)[ind[0][0][2]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][2][0]]);
      cVals[1][0] = (*mGridData)[ind[0][0][0]];
      cVals[1][1] = (*mGridData)[ind[0][0][1]];
      cVals[1][2] = (*mGridData)[ind[0][0][2]];
      cVals[1][3] = extrapolation((*mGridData)[ind[0][0][2]], (*mGridData)[ind[0][0][1]], (*mGridData)[ind[0][0][0]]);
      cVals[2][0] = (*mGridData)[ind[0][1][0]];
      cVals[2][1] = (*mGridData)[ind[0][1][1]];
      cVals[2][2] = (*mGridData)[ind[0][1][2]];
      cVals[2][3] = extrapolation((*mGridData)[ind[0][1][2]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][1][0]]);
      cVals[3][0] = (*mGridData)[ind[0][2][0]];
      cVals[3][1] = (*mGridData)[ind[0][2][1]];
      cVals[3][2] = (*mGridData)[ind[0][2][2]];
      cVals[3][3] = extrapolation((*mGridData)[ind[0][2][2]], (*mGridData)[ind[0][2][1]], (*mGridData)[ind[0][2][0]]);
      cVals[4][0] = extrapolation((*mGridData)[ind[1][0][0]], (*mGridData)[ind[1][1][0]], (*mGridData)[ind[1][2][0]]);
      cVals[4][1] = extrapolation((*mGridData)[ii_x_y_z], (*mGridData)[ind[1][1][1]], (*mGridData)[ind[1][2][1]]);
      cVals[4][2] = extrapolation((*mGridData)[ii_x_y_z], (*mGridData)[ind[1][1][0]], (*mGridData)[ind[1][2][0] + deltaZ[i0]]);
      cVals[4][3] = extrapolation((*mGridData)[ind[1][0][2]], (*mGridData)[ind[1][1][1]], (*mGridData)[ind[1][2][0]]);
      cVals[5][0] = (*mGridData)[ind[1][0][0]];
      cVals[5][2] = (*mGridData)[ind[1][0][2]];
      cVals[5][3] = extrapolation((*mGridData)[ind[1][0][2]], (*mGridData)[ii_x_y_z], (*mGridData)[ind[1][0][0]]);
      cVals[6][0] = (*mGridData)[ind[1][1][0]];
      cVals[6][1] = (*mGridData)[ind[1][1][1]];
      cVals[6][2] = (*mGridData)[ind[1][1][2]];
      cVals[6][3] = extrapolation((*mGridData)[ind[1][1][2]], (*mGridData)[ind[1][1][1]], (*mGridData)[ind[1][1][0]]);
      cVals[7][0] = (*mGridData)[ind[1][2][0]];
      cVals[7][1] = (*mGridData)[ind[1][2][1]];
      cVals[7][2] = (*mGridData)[ind[1][2][2]];
      cVals[7][3] = extrapolation((*mGridData)[ind[1][2][2]], (*mGridData)[ind[1][2][1]], (*mGridData)[ind[1][2][0]]);
      cVals[8][0] = extrapolation((*mGridData)[ind[2][0][0]], (*mGridData)[ind[2][1][0]], (*mGridData)[ind[2][2][0]]);
      cVals[8][1] = extrapolation((*mGridData)[ind[2][0][1]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][2][1]]);
      cVals[8][2] = extrapolation((*mGridData)[ind[2][0][1]], (*mGridData)[ind[2][1][0]], (*mGridData)[ind[2][2][0] + deltaZ[i0]]);
      cVals[8][3] = extrapolation((*mGridData)[ind[2][0][2]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][2][0]]);
      cVals[9][0] = (*mGridData)[ind[2][0][0]];
      cVals[9][1] = (*mGridData)[ind[2][0][1]];
      cVals[9][2] = (*mGridData)[ind[2][0][2]];
      cVals[9][3] = extrapolation((*mGridData)[ind[2][0][2]], (*mGridData)[ind[2][0][1]], (*mGridData)[ind[2][0][0]]);
      cVals[10][0] = (*mGridData)[ind[2][1][0]];
      cVals[10][1] = (*mGridData)[ind[2][1][1]];
      cVals[10][2] = (*mGridData)[ind[2][1][2]];
      cVals[10][3] = extrapolation((*mGridData)[ind[2][1][2]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][1][0]]);
      cVals[11][0] = (*mGridData)[ind[2][2][0]];
      cVals[11][1] = (*mGridData)[ind[2][2][1]];
      cVals[11][2] = (*mGridData)[ind[2][2][2]];
      cVals[11][3] = extrapolation((*mGridData)[ind[2][2][2]], (*mGridData)[ind[2][2][1]], (*mGridData)[ind[2][2][0]]);
      cVals[12][0] = extrapolation((*mGridData)[ind[3][0][0]], (*mGridData)[ind[3][1][0]], (*mGridData)[ind[3][2][0]]);
      cVals[12][1] = extrapolation((*mGridData)[ind[3][0][1]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][2][1]]);
      cVals[12][2] = extrapolation((*mGridData)[ind[3][0][1]], (*mGridData)[ind[3][1][0]], (*mGridData)[ind[3][2][0] + deltaZ[i0]]);
      cVals[12][3] = extrapolation((*mGridData)[ind[3][0][2]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][2][0]]);
      cVals[13][0] = (*mGridData)[ind[3][0][0]];
      cVals[13][1] = (*mGridData)[ind[3][0][1]];
      cVals[13][2] = (*mGridData)[ind[3][0][2]];
      cVals[13][3] = extrapolation((*mGridData)[ind[3][0][2]], (*mGridData)[ind[3][0][1]], (*mGridData)[ind[3][0][0]]);
      cVals[14][0] = (*mGridData)[ind[3][1][0]];
      cVals[14][1] = (*mGridData)[ind[3][1][1]];
      cVals[14][2] = (*mGridData)[ind[3][1][2]];
      cVals[14][3] = extrapolation((*mGridData)[ind[3][1][2]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][1][0]]);
      cVals[15][0] = (*mGridData)[ind[3][2][0]];
      cVals[15][1] = (*mGridData)[ind[3][2][1]];
      cVals[15][2] = (*mGridData)[ind[3][2][2]];
      cVals[15][3] = extrapolation((*mGridData)[ind[3][2][2]], (*mGridData)[ind[3][2][1]], (*mGridData)[ind[3][2][0]]);
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

      cVals[0][0] = extrapolation((*mGridData)[ind[0][0][0]], (*mGridData)[ind[0][0][1]], (*mGridData)[ind[0][0][2]]);
      cVals[0][1] = (*mGridData)[ind[0][0][0]];
      cVals[0][2] = (*mGridData)[ind[0][0][1]];
      cVals[0][3] = (*mGridData)[ind[0][0][2]];
      cVals[1][0] = extrapolation((*mGridData)[ind[0][1][0]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][1][2]]);
      cVals[1][1] = (*mGridData)[ind[0][1][0]];
      cVals[1][2] = (*mGridData)[ind[0][1][1]];
      cVals[1][3] = (*mGridData)[ind[0][1][2]];
      cVals[2][0] = extrapolation((*mGridData)[ind[0][1][0]], (*mGridData)[ind[0][0][1]], (*mGridData)[ind[0][0][2] + deltaR[i0]]);
      cVals[2][1] = (*mGridData)[ind[0][2][0]];
      cVals[2][2] = (*mGridData)[ind[0][2][1]];
      cVals[2][3] = (*mGridData)[ind[0][2][2]];
      cVals[3][0] = extrapolation((*mGridData)[ind[0][2][0]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][0][2]]);
      cVals[3][1] = extrapolation((*mGridData)[ind[0][2][0]], (*mGridData)[ind[0][1][0]], (*mGridData)[ind[0][0][0]]);
      cVals[3][2] = extrapolation((*mGridData)[ind[0][2][1]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][0][1]]);
      cVals[3][3] = extrapolation((*mGridData)[ind[0][2][2]], (*mGridData)[ind[0][1][2]], (*mGridData)[ind[0][0][2]]);
      cVals[4][0] = extrapolation((*mGridData)[ind[1][0][0]], (*mGridData)[ind[1][0][1]], (*mGridData)[ind[1][0][2]]);
      cVals[4][1] = (*mGridData)[ind[1][0][0]];
      cVals[4][2] = (*mGridData)[ind[1][0][1]];
      cVals[4][3] = (*mGridData)[ind[1][0][2]];
      cVals[5][0] = extrapolation((*mGridData)[ii_x_y_z], (*mGridData)[ind[1][1][1]], (*mGridData)[ind[1][1][2]]);
      cVals[5][2] = (*mGridData)[ind[1][1][1]];
      cVals[5][3] = (*mGridData)[ind[1][1][2]];
      cVals[6][0] = extrapolation((*mGridData)[ii_x_y_z], (*mGridData)[ind[1][0][1]], (*mGridData)[ind[1][0][2] + deltaR[i0]]);
      cVals[6][1] = (*mGridData)[ind[1][2][0]];
      cVals[6][2] = (*mGridData)[ind[1][2][1]];
      cVals[6][3] = (*mGridData)[ind[1][2][2]];
      cVals[7][0] = extrapolation((*mGridData)[ind[1][2][0]], (*mGridData)[ind[1][1][1]], (*mGridData)[ind[1][0][2]]);
      cVals[7][1] = extrapolation((*mGridData)[ind[1][2][0]], (*mGridData)[ii_x_y_z], (*mGridData)[ind[1][0][0]]);
      cVals[7][2] = extrapolation((*mGridData)[ind[1][2][1]], (*mGridData)[ind[1][1][1]], (*mGridData)[ind[1][0][1]]);
      cVals[8][0] = extrapolation((*mGridData)[ind[2][0][0]], (*mGridData)[ind[2][0][1]], (*mGridData)[ind[2][0][2]]);
      cVals[7][3] = extrapolation((*mGridData)[ind[1][2][2]], (*mGridData)[ind[1][1][2]], (*mGridData)[ind[1][0][2]]);
      cVals[8][1] = (*mGridData)[ind[2][0][0]];
      cVals[8][2] = (*mGridData)[ind[2][0][1]];
      cVals[8][3] = (*mGridData)[ind[2][0][2]];
      cVals[9][0] = extrapolation((*mGridData)[ind[2][1][0]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][1][2]]);
      cVals[9][1] = (*mGridData)[ind[2][1][0]];
      cVals[9][2] = (*mGridData)[ind[2][1][1]];
      cVals[9][3] = (*mGridData)[ind[2][1][2]];
      cVals[10][0] = extrapolation((*mGridData)[ind[2][1][0]], (*mGridData)[ind[2][0][1]], (*mGridData)[ind[2][0][2] + deltaR[i0]]);
      cVals[10][1] = (*mGridData)[ind[2][2][0]];
      cVals[10][2] = (*mGridData)[ind[2][2][1]];
      cVals[10][3] = (*mGridData)[ind[2][2][2]];
      cVals[11][0] = extrapolation((*mGridData)[ind[2][2][0]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][0][2]]);
      cVals[11][1] = extrapolation((*mGridData)[ind[2][2][0]], (*mGridData)[ind[2][1][0]], (*mGridData)[ind[2][0][0]]);
      cVals[11][2] = extrapolation((*mGridData)[ind[2][2][1]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][0][1]]);
      cVals[11][3] = extrapolation((*mGridData)[ind[2][2][2]], (*mGridData)[ind[2][1][2]], (*mGridData)[ind[2][0][2]]);
      cVals[12][0] = extrapolation((*mGridData)[ind[3][0][0]], (*mGridData)[ind[3][0][1]], (*mGridData)[ind[3][0][2]]);
      cVals[12][1] = (*mGridData)[ind[3][0][0]];
      cVals[12][2] = (*mGridData)[ind[3][0][1]];
      cVals[12][3] = (*mGridData)[ind[3][0][2]];
      cVals[13][0] = extrapolation((*mGridData)[ind[3][1][0]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][1][2]]);
      cVals[13][1] = (*mGridData)[ind[3][1][0]];
      cVals[13][2] = (*mGridData)[ind[3][1][1]];
      cVals[13][3] = (*mGridData)[ind[3][1][2]];
      cVals[14][0] = extrapolation((*mGridData)[ind[3][1][0]], (*mGridData)[ind[3][0][1]], (*mGridData)[ind[3][0][2] + deltaR[i0]]);
      cVals[14][1] = (*mGridData)[ind[3][2][0]];
      cVals[14][2] = (*mGridData)[ind[3][2][1]];
      cVals[14][3] = (*mGridData)[ind[3][2][2]];
      cVals[15][0] = extrapolation((*mGridData)[ind[3][2][0]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][0][2]]);
      cVals[15][1] = extrapolation((*mGridData)[ind[3][2][0]], (*mGridData)[ind[3][1][0]], (*mGridData)[ind[3][0][0]]);
      cVals[15][2] = extrapolation((*mGridData)[ind[3][2][1]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][0][1]]);
      cVals[15][3] = extrapolation((*mGridData)[ind[3][2][2]], (*mGridData)[ind[3][1][2]], (*mGridData)[ind[3][0][2]]);
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

      cVals[0][0] = (*mGridData)[ind[0][0][0]];
      cVals[0][1] = (*mGridData)[ind[0][0][1]];
      cVals[0][2] = (*mGridData)[ind[0][0][2]];
      cVals[0][3] = extrapolation((*mGridData)[ind[0][0][2]], (*mGridData)[ind[0][0][1]], (*mGridData)[ind[0][0][0]]);
      cVals[1][0] = (*mGridData)[ind[0][1][0]];
      cVals[1][1] = (*mGridData)[ind[0][1][1]];
      cVals[1][2] = (*mGridData)[ind[0][1][2]];
      cVals[1][3] = extrapolation((*mGridData)[ind[0][1][2]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][1][0]]);
      cVals[2][0] = (*mGridData)[ind[0][2][0]];
      cVals[2][1] = (*mGridData)[ind[0][2][1]];
      cVals[2][2] = (*mGridData)[ind[0][2][2]];
      cVals[2][3] = extrapolation((*mGridData)[ind[0][2][2]], (*mGridData)[ind[0][2][1]], (*mGridData)[ind[0][2][0]]);
      cVals[3][0] = extrapolation((*mGridData)[ind[0][2][0]], (*mGridData)[ind[0][1][0]], (*mGridData)[ind[0][0][0]]);
      cVals[3][1] = extrapolation((*mGridData)[ind[0][2][1]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][0][1]]);
      cVals[3][2] = extrapolation((*mGridData)[ind[0][2][2]], (*mGridData)[ind[0][1][2]], (*mGridData)[ind[0][0][2]]);
      cVals[3][3] = extrapolation((*mGridData)[ind[0][2][2]], (*mGridData)[ind[0][1][1]], (*mGridData)[ind[0][0][0]]);
      cVals[4][0] = (*mGridData)[ind[1][0][0]];
      cVals[4][1] = (*mGridData)[ind[1][0][1]];
      cVals[4][2] = (*mGridData)[ind[1][0][2]];
      cVals[4][3] = extrapolation((*mGridData)[ind[1][0][2]], (*mGridData)[ind[1][0][1]], (*mGridData)[ind[1][0][0]]);
      cVals[5][0] = (*mGridData)[ind[1][1][0]];
      cVals[5][2] = (*mGridData)[ind[1][1][2]];
      cVals[5][3] = extrapolation((*mGridData)[ind[1][1][2]], (*mGridData)[ii_x_y_z], (*mGridData)[ind[1][1][0]]);
      cVals[6][0] = (*mGridData)[ind[1][2][0]];
      cVals[6][1] = (*mGridData)[ind[1][2][1]];
      cVals[6][2] = (*mGridData)[ind[1][2][2]];
      cVals[6][3] = extrapolation((*mGridData)[ind[1][2][2]], (*mGridData)[ind[1][2][1]], (*mGridData)[ind[1][2][0]]);
      cVals[7][0] = extrapolation((*mGridData)[ind[1][2][0]], (*mGridData)[ind[1][1][0]], (*mGridData)[ind[1][0][0]]);
      cVals[7][1] = extrapolation((*mGridData)[ind[1][2][1]], (*mGridData)[ii_x_y_z], (*mGridData)[ind[1][0][1]]);
      cVals[7][2] = extrapolation((*mGridData)[ind[1][2][2]], (*mGridData)[ind[1][1][2]], (*mGridData)[ind[1][0][2]]);
      cVals[7][3] = extrapolation((*mGridData)[ind[1][2][2]], (*mGridData)[ind[1][1][1]], (*mGridData)[ind[1][0][0]]);
      cVals[8][0] = (*mGridData)[ind[2][0][0]];
      cVals[8][1] = (*mGridData)[ind[2][0][1]];
      cVals[8][2] = (*mGridData)[ind[2][0][2]];
      cVals[8][3] = extrapolation((*mGridData)[ind[2][0][2]], (*mGridData)[ind[2][0][1]], (*mGridData)[ind[2][0][0]]);
      cVals[9][0] = (*mGridData)[ind[2][1][0]];
      cVals[9][1] = (*mGridData)[ind[2][1][1]];
      cVals[9][2] = (*mGridData)[ind[2][1][2]];
      cVals[9][3] = extrapolation((*mGridData)[ind[2][1][2]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][1][0]]);
      cVals[10][0] = (*mGridData)[ind[2][2][0]];
      cVals[10][1] = (*mGridData)[ind[2][2][1]];
      cVals[10][2] = (*mGridData)[ind[2][2][2]];
      cVals[10][3] = extrapolation((*mGridData)[ind[2][2][2]], (*mGridData)[ind[2][2][1]], (*mGridData)[ind[2][2][0]]);
      cVals[11][0] = extrapolation((*mGridData)[ind[2][2][0]], (*mGridData)[ind[2][1][0]], (*mGridData)[ind[2][0][0]]);
      cVals[11][1] = extrapolation((*mGridData)[ind[2][2][1]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][0][1]]);
      cVals[11][2] = extrapolation((*mGridData)[ind[2][2][2]], (*mGridData)[ind[2][1][2]], (*mGridData)[ind[2][0][2]]);
      cVals[11][3] = extrapolation((*mGridData)[ind[2][2][2]], (*mGridData)[ind[2][1][1]], (*mGridData)[ind[2][0][0]]);
      cVals[12][0] = (*mGridData)[ind[3][0][0]];
      cVals[12][1] = (*mGridData)[ind[3][0][1]];
      cVals[12][2] = (*mGridData)[ind[3][0][2]];
      cVals[12][3] = extrapolation((*mGridData)[ind[3][0][2]], (*mGridData)[ind[3][0][1]], (*mGridData)[ind[3][0][0]]);
      cVals[13][0] = (*mGridData)[ind[3][1][0]];
      cVals[13][1] = (*mGridData)[ind[3][1][1]];
      cVals[13][2] = (*mGridData)[ind[3][1][2]];
      cVals[13][3] = extrapolation((*mGridData)[ind[3][1][2]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][1][0]]);
      cVals[14][0] = (*mGridData)[ind[3][2][0]];
      cVals[14][1] = (*mGridData)[ind[3][2][1]];
      cVals[14][2] = (*mGridData)[ind[3][2][2]];
      cVals[14][3] = extrapolation((*mGridData)[ind[3][2][2]], (*mGridData)[ind[3][2][1]], (*mGridData)[ind[3][2][0]]);
      cVals[15][0] = extrapolation((*mGridData)[ind[3][2][0]], (*mGridData)[ind[3][1][0]], (*mGridData)[ind[3][0][0]]);
      cVals[15][1] = extrapolation((*mGridData)[ind[3][2][1]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][0][1]]);
      cVals[15][2] = extrapolation((*mGridData)[ind[3][2][2]], (*mGridData)[ind[3][1][2]], (*mGridData)[ind[3][0][2]]);
      cVals[15][3] = extrapolation((*mGridData)[ind[3][2][2]], (*mGridData)[ind[3][1][1]], (*mGridData)[ind[3][0][0]]);
    } break;
  }
}

template class o2::tpc::TriCubicInterpolator<double>;
template class o2::tpc::TriCubicInterpolator<float>;
