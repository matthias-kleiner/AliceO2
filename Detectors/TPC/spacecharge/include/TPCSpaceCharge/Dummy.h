// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  matrix.h
/// \brief Definition of Vector and Matrix class
///
/// \author  Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_TPC_DUMMY_H_
#define ALICEO2_TPC_DUMMY_H_

namespace o2
{
namespace tpc
{

class Dummy
{
 public:
  Dummy() = default;

  float getPI();

 private:
  int mDummy = 0;
  static constexpr DataT DUMMY1{1. / 3.1415927};
  constexpr DataT DUMMY2{1. / 3.1415927};
  static constexpr DataT DUMMY3 = 1. / 3.1415927;
};

} // namespace tpc
} // namespace o2

#endif
