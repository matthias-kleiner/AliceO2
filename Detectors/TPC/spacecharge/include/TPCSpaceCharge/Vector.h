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

/// \file  Vector.h
/// \brief Definition of a simple fixed size vector which is used e.g. in RegularGrid3D
///
/// \author  Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_TPC_VECTOR_H_
#define ALICEO2_TPC_VECTOR_H_

#include <array>
#include <cstddef>

namespace o2
{
namespace tpc
{

/// simple fixed size vector
/// \tparam DataT type of the data
/// \tparam N number of elements
template <typename DataT, size_t N>
class Vector
{
 public:
  /// default constructor
  Vector() = default;

  /// constructor
  /// \param data values of the vector
  Vector(const std::array<DataT, N>& data) : mData{data} {}

  /// operator access
  DataT operator[](size_t i) const { return mData[i]; }
  DataT& operator[](size_t i) { return mData[i]; }

  /// \return returns the number of elements
  static constexpr size_t size() { return N; }

 private:
  std::array<DataT, N> mData{}; ///< storage for the data
};

} // namespace tpc
} // namespace o2

#endif
