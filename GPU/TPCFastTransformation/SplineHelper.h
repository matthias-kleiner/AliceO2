// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  SplineHelper.h
/// \brief Definition of SplineHelper class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINEHELPER_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINEHELPER_H

#include <cmath>
#include <vector>

#include "GPUCommonDef.h"
#include "Rtypes.h"
#include "TString.h"
#include "Spline1D.h"
#include "Spline.h"
#include "SplineHelper1D.h"
#include <functional>

namespace GPUCA_NAMESPACE
{
namespace gpu
{

///
/// The SplineHelper class is to initialize Spline* objects
///
template <typename DataT>
class SplineHelper
{
 public:
  /// _____________  Constructors / destructors __________________________

  /// Default constructor
  SplineHelper();

  /// Copy constructor: disabled
  SplineHelper(const SplineHelper&) CON_DELETE;

  /// Assignment operator: disabled
  SplineHelper& operator=(const SplineHelper&) CON_DELETE;

  /// Destructor
  ~SplineHelper() CON_DEFAULT;

  /// _______________  Main functionality  ________________________

  /// Create best-fit spline parameters for a given input function F
  template <bool isConsistentT>
  void approximateFunction(SplineBase<DataT, isConsistentT>& spline,
                           const DataT xMin[/* Xdim */], const DataT xMax[/* Xdim */],
                           std::function<void(std::vector<DataT> x[/* Xdim */], DataT f[/* Fdim */])> F, const int batchsize,
                           const int nAxiliaryDataPoints[/* Xdim */] = nullptr);

  /// _______________   Interface for a step-wise construction of the best-fit spline   ________________________

  /// precompute everything needed for the construction
  template <bool isConsistentT>
  int setSpline(const SplineBase<DataT, isConsistentT>& spline, const int nAxiliaryPoints[/* Xdim */]);

  /// approximate std::function, output in Fparameters
  void approximateFunction(
    DataT* Fparameters, const DataT xMin[/* mXdimensions */], const DataT xMax[/* mXdimensions */],
    std::function<void(std::vector<DataT> x[/* mXdimensions */], DataT f[/* mFdimensions */])> F, const int batchsize) const;

  /// approximate std::function, output in Fparameters. F calculates values for a batch of points.
  void approximateFunctionBatch(
    DataT* Fparameters, const DataT xMin[/* mXdimensions */], const DataT xMax[/* mXdimensions */],
    std::function<void(const std::vector<DataT> x[/* mXdimensions */], std::vector<DataT> f[/*mFdimensions*/])> F,
    unsigned int batchsize) const;

  /// approximate a function given as an array of values at data points
  void approximateFunction(
    DataT* Fparameters, const DataT DataPointF[/*getNumberOfDataPoints() x nFdim*/]) const;

  int getNumberOfDataPoints(int dimX) const { return mHelpers[dimX].getNumberOfDataPoints(); }

  int getNumberOfDataPoints() const { return mNumberOfDataPoints; }

  const SplineHelper1D<DataT>& getHelper(int dimX) const { return mHelpers[dimX]; }

  /// _______________  Utilities   ________________________

  ///  Gives error string
  const char* getLastError() const { return mError.Data(); }


  static int arraytopoints ( int point, int result[], const int numbers[],  int dim);

  static int pointstoarray ( const int indices[], const int numbers[], int dim);

 private:
  /// Stores an error message
  int storeError(Int_t code, const char* msg);

  TString mError = "";     ///< error string
  int mXdimensions;        ///< number of X dimensions
  int mFdimensions;        ///< number of F dimensions
  int mNumberOfParameters; ///< number of parameters
  int mNumberOfDataPoints; ///< number of data points
  std::vector<SplineHelper1D<DataT>> mHelpers;
};



template <typename DataT>
template <bool isConsistentT>
void SplineHelper<DataT>::approximateFunction(
  SplineBase<DataT, isConsistentT>& spline,
  const DataT xMin[/* Xdim */], const DataT xMax[/* Xdim */],
  std::function<void(std::vector<DataT> x[/* Xdim */], DataT f[/* Fdim */])> F, const int batchsize,
  const int nAxiliaryDataPoints[/* Xdim */])
{
  /// Create best-fit spline parameters for a given input function F
  if (spline.isConsistent()) {
    setSpline(spline, nAxiliaryDataPoints);
    approximateFunction(spline.getFparameters(), xMin, xMax, F, batchsize);
  }
  spline.setXrange(xMin, xMax);
}

template <typename DataT>
template <bool isConsistentT>
int SplineHelper<DataT>::setSpline(
  const SplineBase<DataT, isConsistentT>& spline, const int nAxiliaryPoints[/* Xdim */])
{
  // Prepare creation of an irregular spline
  // The should be at least one (better, two) axiliary measurements on each segnment between two knots and at least 2*nKnots measurements in total
  // Returns 0 when the spline can not be constructed with the given nAxiliaryPoints

  int ret = 0;
  mXdimensions = spline.getXdimensions();
  mFdimensions = spline.getFdimensions();
  mNumberOfParameters = spline.getNumberOfParameters();
  mNumberOfDataPoints = 1;
  mHelpers.clear();
  mHelpers.resize(mXdimensions);
  for (int i = 0; i < mXdimensions; i++) {
    int np = (nAxiliaryPoints != nullptr) ? nAxiliaryPoints[i] : 4;
    if (mHelpers[i].setSpline(spline.getGrid(i), mFdimensions, np) != 0) {
      ret = storeError(-2, "SplineHelper::setSpline: error by setting an axis");
    }
    mNumberOfDataPoints *= mHelpers[i].getNumberOfDataPoints();
  }

  return ret;
}

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
