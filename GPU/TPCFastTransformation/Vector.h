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

/// \file  matrix.h
/// \brief Definition of Vector and MatrixGPU class
///
/// \author  Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_TPC_VECTORGPU_H_
#define ALICEO2_TPC_VECTORGPU_H_

#include <Vc/vector>
#include <Vc/Memory>

namespace o2
{
namespace tpc
{
template <typename DataT = float, size_t N = 64>
class MatrixGPU
{
  using VDataT = Vc::Vector<DataT>;

 public:
  /// constructor
  /// \param dataMatrix pointer to the data
  MatrixGPU(const Vc::Memory<VDataT, N>* dataMatrix) : mDataMatrix(dataMatrix) {}

  const Vc::Memory<VDataT, N>& operator[](size_t i) const { return mDataMatrix[i]; }

 private:
  const Vc::Memory<VDataT, N>* mDataMatrix{};
};

template <typename DataT = float, size_t N = 64>
class VectorGPU
{
  using VDataT = Vc::Vector<DataT>;

 public:
  /// default constructor
  VectorGPU() = default;

  /// constructor
  /// \param dataVector data which is assigned to the vector
  VectorGPU(const Vc::Memory<VDataT, N>& dataVector) : mDataVector(dataVector) {}

  /// operator access
  const DataT operator[](size_t i) const { return mDataVector.scalar(i); }
  DataT& operator[](size_t i) { return mDataVector.scalar(i); }

  /// sets the vector with index j
  void setVector(const size_t j, const VDataT& vector) { mDataVector.vector(j) = vector; }

  /// \return returns the vector with index j
  const VDataT getVector(const size_t j) const { return mDataVector.vector(j); }

  /// \return returns the number of Vc::VectorGPU<DataT> stored in the Vector
  size_t getvectorsCount() const { return mDataVector.vectorsCount(); }

  /// \return returns the number of entries stored in the Vector
  size_t getentriesCount() const { return mDataVector.entriesCount(); }

 private:
  // storage for the data
  Vc::Memory<VDataT, N> mDataVector{};
};

template <typename DataT, size_t N>
inline VectorGPU<DataT, N> operator*(const MatrixGPU<DataT, N>& a, const VectorGPU<DataT, N>& b)
{
  using V = Vc::Vector<DataT>;
  // resulting vector c
  VectorGPU<DataT, N> c;
  for (size_t i = 0; i < N; ++i) {
    V c_ij{};
    for (size_t j = 0; j < a[i].vectorsCount(); ++j) {
      c_ij += a[i].vector(j) * b.getVector(j);
    }
    c[i] = c_ij.sum();
  }
  return c;
}

template <typename DataT, size_t N>
inline VectorGPU<DataT, N> floorVEC(const VectorGPU<DataT, N>& a)
{
  VectorGPU<DataT, N> c;
  for (size_t j = 0; j < a.getvectorsCount(); ++j) {
    c.setVector(j, Vc::floor(a.getVector(j)));
  }
  return c;
}

template <typename DataT, size_t N>
inline VectorGPU<DataT, N> operator-(const VectorGPU<DataT, N>& a, const VectorGPU<DataT, N>& b)
{
  // resulting matrix c
  VectorGPU<DataT, N> c;
  for (size_t j = 0; j < a.getvectorsCount(); ++j) {
    c.setVector(j, a.getVector(j) - b.getVector(j));
  }
  return c;
}

template <typename DataT, size_t N>
inline VectorGPU<DataT, N> operator+(const VectorGPU<DataT, N>& a, const VectorGPU<DataT, N>& b)
{
  // resulting matrix c
  VectorGPU<DataT, N> c;
  for (size_t j = 0; j < a.getvectorsCount(); ++j) {
    c.setVector(j, a.getVector(j) + b.getVector(j));
  }
  return c;
}

template <typename DataT, size_t N>
inline VectorGPU<DataT, N> operator*(const DataT a, const VectorGPU<DataT, N>& b)
{
  // resulting matrix c
  VectorGPU<DataT, N> c;
  for (size_t j = 0; j < b.getvectorsCount(); ++j) {
    c.setVector(j, a * b.getVector(j));
  }
  return c;
}

// compute the sum of one Vector
template <typename DataT, size_t N>
inline DataT sum(const VectorGPU<DataT, N>& a)
{
  // resulting matrix c
  Vc::Vector<DataT> b = a.getVector(0);
  for (size_t j = 1; j < a.getvectorsCount(); ++j) {
    b += a.getVector(j);
  }
  return b.sum();
}

// multiply each row from a vector with the row from a second vector
template <typename DataT, size_t N>
inline VectorGPU<DataT, N> operator*(const VectorGPU<DataT, N>& a, const VectorGPU<DataT, N>& b)
{
  // resulting matrix c
  VectorGPU<DataT, N> c{};
  for (size_t j = 0; j < a.getvectorsCount(); ++j) {
    c.setVector(j, a.getVector(j) * b.getVector(j));
  }
  return c;
}

// check if all elements are equal
template <typename DataT, size_t N>
inline bool operator==(const VectorGPU<DataT, N>& a, const VectorGPU<DataT, N>& b)
{
  for (size_t j = 0; j < a.getvectorsCount(); ++j) {
    if (any_of(a.getVector(j) != b.getVector(j))) {
      return false;
    }
  }
  return true;
}

} // namespace tpc
} // namespace o2

#endif
