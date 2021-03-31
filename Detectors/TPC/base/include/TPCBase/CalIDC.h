// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_TPC_CalIDC_H_
#define ALICEO2_TPC_CalIDC_H_

#include <memory>
#include <vector>
#include <string>
#include <cassert>

#include "DataFormatsTPC/Defs.h"
#include "TPCBase/Mapper.h"
#include "TPCBase/ROC.h"
#include "TPCBase/Sector.h"

#ifndef GPUCA_ALIGPUCODE
#include <Framework/Logger.h>
#include <boost/format.hpp>
#include <boost/range/combine.hpp>
#endif

namespace o2
{
namespace tpc
{
/// Class to hold calibration data on a pad level
///
template <class DataT = float>
class CalIDC
{
 public:
  CalIDC(const CRU cru, const int mergeNRows, const int mergeNPads) : mCRU{cru}, mMergeNrows{mergeNRows}, mMergeNpads{mergeNPads} {
    const static Mapper& mapper = Mapper::instance();
    const auto regionInfo = mapper.getPadRegionInfo(cru.region());
    int rows = mapper.getNumberOfRowsRegion(cru.region());

    int nPadsinRow = mapper.getNumberOfPadsInRowRegion(cru.region(), rows);
  };

  CalIDC(CalIDC const&) = default;
  CalIDC& operator=(CalIDC const&) = default;
  ~CalIDC() = default;

  const DataT getValue(const int sec, const int globalPadInSector) const;
  void setValue(const int sec, const int globalPadInSector, const DataT& value);
  void setValue(const int sec, const int rowInSector, const int padInRow, const DataT& value);

  /// \todo return value of T& not possible if a default value should be returned, e.g. T{}:
  ///       warning: returning reference to temporary
  const DataT getValue(const ROC roc, const size_t row, const size_t pad) const;
  const DataT getValue(const CRU cru, const size_t row, const size_t pad) const;
  const DataT getValue(const Sector sec, const int rowInSector, const int padInRow) const;

  // void setName(const std::string_view name) { mName = name.data(); }
  // const std::string& getName() const { return mName; }

  const CalIDC& multiply(const DataT& val) { return *this *= val; }
  const CalIDC& operator+=(const CalIDC& other);
  const CalIDC& operator-=(const CalIDC& other);
  const CalIDC& operator*=(const CalIDC& other);
  const CalIDC& operator/=(const CalIDC& other);

  const CalIDC& operator+=(const DataT& val);
  const CalIDC& operator-=(const DataT& val);
  const CalIDC& operator*=(const DataT& val);
  const CalIDC& operator/=(const DataT& val);

  template <class U>
  friend CalIDC<U> operator+(const CalIDC<U>&, const CalIDC<U>&);

  template <class U>
  friend CalIDC<U> operator-(const CalIDC<U>&, const CalIDC<U>&);

 private:
  static constexpr int sSectors = 36;
  static constexpr int sRegions = 10;
  static constexpr int sRowsPerRegion[sRegions] = {17, 15, 16, 15, 18, 16, 16, 14, 13, 12};

  const int mCRU = 1;
  const int mMergeNrows = 1;
  const int mMergeNpads = 1;
  std::vector<DataT> mData;
};

} // namespace tpc
} // namespace o2

#endif
