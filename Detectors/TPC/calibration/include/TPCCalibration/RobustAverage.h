// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file RobustAverage.h
/// \brief class for performing robust averaging and outlier filtering
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_ROBUSTAVERAGE_H_
#define ALICEO2_ROBUSTAVERAGE_H_

#include <vector>
#include <numeric>
#include "Framework/Logger.h"

namespace o2
{
namespace tpc
{
/// class to perform filtering of outliers and robust averaging of a set of values.
/// This class is more or less a dummy for now... TODO add more sophisticated methods

class RobustAverage
{
 public:
  /// \param maxValues maximum number of values which will be averaged
  RobustAverage(const unsigned int maxValues) { mValues.reserve(maxValues); }

  /// clear the stored values
  void clear() { mValues.clear(); }

  /// \param value value which will be added to the list of stored values for averaging
  void addValue(const float value)
  {
    mValues.emplace_back(value);
  }

  /// returns the average value
  float getFilteredAverage(const float sigma = 3)
  {
    if (mValues.empty()) {
      return 0;
    }

    const float mean = getMean();
    const float stdev = getStdDev(mean);
    filterOutliers(mean, stdev, sigma);
    return getMean();
  }

  float getMean() const
  {
    return std::accumulate(mValues.begin(), mValues.end(), decltype(mValues)::value_type(0)) / mValues.size();
  }

  /// performing outlier filtering of the stored values
  float getStdDev(const float mean) const
  {
    std::vector<float> diff(mValues.size());
    std::transform(mValues.begin(), mValues.end(), diff.begin(), [mean](const float val) { return val - mean; });

    const float sqsum = std::inner_product(diff.begin(), diff.end(), diff.begin(), decltype(mValues)::value_type(0));
    const float stdev = std::sqrt(sqsum / diff.size());
    return stdev;
  }

  /// performing outlier filtering of the stored values by defining range of included values in terms of standard deviation
  /// \param mean mean of the stored values
  /// \param stdev standard deviation of the values
  /// \param sigma maximum accepted standard deviation: sigma*stdev
  void filterOutliers(const float mean, const float stdev, const float sigma = 3)
  {
    std::sort(mValues.begin(), mValues.end());
    const float minVal = mean - sigma * stdev;
    const float maxVal = mean + sigma * stdev;
    const auto upper = std::upper_bound(mValues.begin(), mValues.end(), maxVal);
    const auto lower = std::lower_bound(mValues.begin(), mValues.end(), minVal);
    mValues.erase(upper, mValues.end());
    mValues.erase(mValues.begin(), lower);
  }

  void print() const
  {
    LOGP(info, "PRINTING STORED VALUES");
    for (auto& val : mValues) {
      LOGP(info, "{}", val);
    }
  }

 private:
  std::vector<float> mValues{}; ///< values which will be averaged and filtered
};

} // namespace tpc
} // namespace o2

#endif
