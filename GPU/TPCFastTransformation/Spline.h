// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  Spline.h
/// \brief Definition of Spline class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINE_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINE_H

#include "Spline1D.h"
#include "FlatObject.h"
#include "GPUCommonDef.h"

#if !defined(__CINT__) && !defined(__ROOTCINT__) && !defined(GPUCA_GPUCODE) && !defined(GPUCA_NO_VC) && defined(__cplusplus) && __cplusplus >= 201703L
#include <Vc/Vc>
#include <Vc/SimdArray>
#endif

class TFile;

namespace GPUCA_NAMESPACE
{
namespace gpu
{
///
/// The Spline class performs a cubic spline interpolation on a multi-dimensional nonunifom grid.
/// The class is an extension of the Spline1D class.
/// See Spline1D.h for more details.
///
/// The spline S(x) approximates a function F(x):R^n->R^m,
/// where x is an n-dimentional vector which belongs to some area.
/// The F value may be also multi-dimensional.
///
/// --- Example of creating a spline ---
///
///  constexpr int nXdim=2, nFdim=1;
///  auto F = [&](float[nXdim] x, float f[nFdim] ) {
///   f[0] = 1.f + x[0] + x[1]*x[1]; // F(x)
///  };
///  const int nKnots[nXdim]={2,3};
///  const int knotsU1[nKnots[0]] = {0, 1};
///  const int knotsU2[nKnots[1]] = {0, 2, 5};
///  const int *knotsU[nXdim] = { knotsU1, knotsU2 };
///  Spline<float, nXdim, nFdim> spline(nKnots, knotsU ); // create a grid and allocate memory for 1-dimensional F
///  float xMin[nXdim] = {0.,0.};
///  float xMax[nXdim] = {1.,1.};
///  spline.approximateFunction(xMin, xMax, F); //initialize spline to approximate F on area [0., 1.]x[0., 1.]
///  float S = spline.interpolate( {.1, .3} ); // interpolated value at x=(.1,.3)
///
///  --- See also Spline::test();
//

/// Base class to store data members and non-inline methods
template <typename DataT, bool isConsistentT>
class SplineBase : public FlatObject
{
 public:
  /// _____________  Version control __________________________

  /// Version control
  GPUhd() static constexpr int getVersion() { return 1; }

    /// _____________  Constructors / destructors __________________________

#if !defined(GPUCA_GPUCODE)
  /// constructor && default constructor for the ROOT streamer
  SplineBase(int nXdim = 1, int nFdim = 1);
#else
  /// Disable constructors
  SplineBase() CON_DELETE;
#endif
  /// Disable constructors
  SplineBase(const SplineBase&) CON_DELETE;
  SplineBase& operator=(const SplineBase&) CON_DELETE;

  /// Destructor
  ~SplineBase() CON_DEFAULT;

  /// _______________  Construction interface  ________________________

#if !defined(GPUCA_GPUCODE)
  /// Constructor for a regular spline.
  void recreate(const int numberOfKnots[/* mXdim */]);

  /// Constructor for an irregular spline
  void recreate(const int numberOfKnots[/* mXdim */], const int* const knots[/* mXdim */]);

  /// approximate a function F with this spline.
  void approximateFunction(const DataT xMin[/* mXdim */], const DataT xMax[/* mXdim */],
                           std::function<void(std::vector<DataT> x[/* mXdim */], DataT f[/* mFdim */])> F, const int batchsize,
                           const int nAxiliaryDataPoints[/* mXdim */] = nullptr);
#endif

  /// _______________  IO   ________________________

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
  /// write a class object to the file
  int writeToFile(TFile& outf, const char* name);

  /// read a class object from the file
  static SplineBase* readFromFile(TFile& inpf, const char* name);
#endif

  /// _______________  Getters   ________________________

  /// Get number of F dimensions
  GPUhd() int getXdimensions() const { return mXdim; }

  /// Get number of F dimensions
  GPUhd() int getFdimensions() const { return mFdim; }

  ///
  GPUhd() static constexpr bool isConsistent() { return isConsistentT; }

  /// Get minimal required alignment for the spline parameters
  GPUhd() static constexpr size_t getParameterAlignmentBytes() { return 16; }

  /// Number of parameters
  GPUhd() int getNumberOfParametersPerKnot() const
  {
    return (1 << mXdim) * mFdim; // 2^mXdim parameters per F dimension
  }

  /// Number of parameters
  GPUhd() int getNumberOfParameters() const { return mNparameters; }

  /// Size of the parameter array in bytes
  GPUhd() size_t getSizeOfParameters() const { return sizeof(DataT) * getNumberOfParameters(); }

  /// Get 1-D grid for dim dimension
  GPUhd() const Spline1D<DataT>& getGrid(int dim) const { return mGrid[dim]; }

  /// Get u[] coordinate of i-th knot
  GPUhd() void getKnotU(int iKnot, DataT u[/* mXdim */]) const;

  /// Get index of a knot (iKnot1,iKnot2,..,iKnotN)
  GPUhd() int getKnotIndex(const int iKnot[/* mXdim */]) const;

  /// Get number of F parameters
  GPUhd() DataT* getFparameters() { return mFparameters; }

  /// Get number of F parameters
  GPUhd() const DataT* getFparameters() const { return mFparameters; }

  /// _______________  Technical stuff  ________________________

  /// Set X range
  void setXrange(const DataT xMin[/* mXdim */], const DataT xMax[/* mXdim */]);

  /// Print method
  void print() const;

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
  /// Test the class functionality
  static int test(const bool draw = 0, const bool drawDataPoints = 1);
#endif

  /// _____________  FlatObject functionality, see FlatObject class for description  ____________

  using FlatObject::getBufferAlignmentBytes;
  using FlatObject::getClassAlignmentBytes;

#if !defined(GPUCA_GPUCODE)
  void cloneFromObject(const SplineBase& obj, char* newFlatBufferPtr);
  void moveBufferTo(char* newBufferPtr);
#endif

  using FlatObject::releaseInternalBuffer;

  void destroy();
  void setActualBufferAddress(char* actualFlatBufferPtr);
  void setFutureBufferAddress(char* futureFlatBufferPtr);

 private:
  int mXdim;        ///< dimentionality of X
  int mFdim;        ///< dimentionality of F
  int mNparameters; ///< number of spline parameters

 protected:               /// _____________  Data members  ____________
  Spline1D<DataT>* mGrid; //! (transient!!) mXdim grids
  DataT* mFparameters;    //! (transient!!) F-dependent parameters of the spline

  ClassDefNV(SplineBase, 1);
};

///
/// The main Spline class. Contains constructors and interpolation.
///
/// F dimensions can be set as a template parameter (faster; the only option for the GPU)
/// or as a constructor parameter (slower; the only option for the ROOT interpretator).
/// In a compiled CPU code one can use both options.
///
template <typename DataT, int nXdimT = 1, int nFdimT = 1, bool isConsistentT = true>
class Spline : public SplineBase<DataT, isConsistentT>
{
 public:
  typedef SplineBase<DataT, isConsistentT> BaseT;

  /// _____________  Constructors / destructors __________________________

#if !defined(GPUCA_GPUCODE)

  /// Constructor for a regular spline && default constructor
  Spline(const int numberOfKnots[] = nullptr);

  /// Constructor for an irregular spline
  Spline(const int numberOfKnots[], const int* const knots[]);

  /// Copy constructor
  Spline(const Spline&);

  /// Assignment operator
  Spline& operator=(const Spline&);
#else
  /// Disable constructors for the GPU implementation
  Spline() CON_DELETE;
  Spline(const Spline&) CON_DELETE;
  Spline& operator=(const Spline&) CON_DELETE;
#endif

  /// Destructor
  ~Spline() CON_DEFAULT;

  /// _______________  Main functionality   ________________________

  /// Get interpolated value for F(x1,x2)
  GPUhd() void interpolate(GPUgeneric() const DataT x[], GPUgeneric() DataT S[]) const;

  /// Get interpolated value for the first dimension of F. (Simplified interface for 1D)
  GPUhd() DataT interpolate(GPUgeneric() const DataT x[]) const;

  /// ================ Expert tools   ================================

  /// Get interpolated value S for an nFdim-dimensional F(u) using spline parameters Fparameters.
  /// Fparameters can be created via SplineHelper.
  GPUhd() void interpolateU(GPUgeneric() const DataT Fparameters[],
                            GPUgeneric() const DataT u[], GPUgeneric() DataT S[]) const;

  /// _______________  IO   ________________________

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
  /// write a class object to the file
  using BaseT::writeToFile;

  /// read a class object from the file
  static Spline* readFromFile(TFile& inpf, const char* name)
  {
    return reinterpret_cast<Spline*>(BaseT::readFromFile(inpf, name));
  }
#endif

  /// _______________  Getters   ________________________

  /// Get number of X dimensions.
  GPUhd() static constexpr int getXdimensions() { return nXdimT; }

  /// Get number of F dimensions.
  GPUhd() static constexpr int getFdimensions() { return nFdimT; }

  using BaseT::getKnotIndex;
  using BaseT::mFparameters;
  using BaseT::mGrid;
};

///
/// ========================================================================================================
///       Inline implementations of some methods
/// ========================================================================================================
///

template <typename DataT, bool isConsistentT>
GPUhdi() void SplineBase<DataT, isConsistentT>::getKnotU(int iKnot, DataT u[/* mXdim */]) const
{
  /// Get u[] coordinate of i-th knot
  for (int dim = 0; dim < mXdim; dim++) {
    int n = mGrid[dim].getNumberOfKnots();
    u[dim] = mGrid[dim].getKnot(iKnot % n).u;
    iKnot /= n;
  }
}

template <typename DataT, bool isConsistentT>
GPUhdi() int SplineBase<DataT, isConsistentT>::getKnotIndex(const int iKnot[/* mXdim */]) const
{
  /// Get index of a knot (iKnot1,iKnot2,..,iKnotN)
  int ind = iKnot[0];
  int n = 1;
  for (int dim = 1; dim < mXdim; dim++) {
    n *= mGrid[dim - 1].getNumberOfKnots();
    ind += n * iKnot[dim];
  }
  return ind;
}

template <typename DataT, bool isConsistentT>
GPUhdi() void SplineBase<DataT, isConsistentT>::
  setXrange(const DataT xMin[/* mXdim */], const DataT xMax[/* mXdim */])
{
  /// Set X range
  for (int i = 0; i < mXdim; i++) {
    mGrid[i].setXrange(xMin[i], xMax[i]);
  }
}

#if !defined(GPUCA_GPUCODE)

template <typename DataT, int nXdimT, int nFdimT, bool isConsistentT>
GPUhdi() Spline<DataT, nXdimT, nFdimT, isConsistentT>::
  Spline(const int numberOfKnots[])
  : SplineBase<DataT, isConsistentT>(nXdimT, nFdimT)
{
  this->recreate(numberOfKnots);
}

template <typename DataT, int nXdimT, int nFdimT, bool isConsistentT>
GPUhdi() Spline<DataT, nXdimT, nFdimT, isConsistentT>::
  Spline(const int numberOfKnots[], const int* const knots[])
  : SplineBase<DataT, isConsistentT>(nXdimT, nFdimT)
{
  this->recreate(numberOfKnots, knots);
}

template <typename DataT, int nXdimT, int nFdimT, bool isConsistentT>
GPUhdi() Spline<DataT, nXdimT, nFdimT, isConsistentT>::
  Spline(const Spline& spline)
  : SplineBase<DataT, isConsistentT>(nXdimT, nFdimT)
{
  this->cloneFromObject(spline, nullptr);
}

template <typename DataT, int nXdimT, int nFdimT, bool isConsistentT>
GPUhdi() Spline<DataT, nXdimT, nFdimT, isConsistentT>& Spline<DataT, nXdimT, nFdimT, isConsistentT>::
  operator=(const Spline& spline)
{
  this->cloneFromObject(spline, nullptr);
  return *this;
}
#endif

template <typename DataT, int nXdimT, int nFdimT, bool isConsistentT>
GPUhdi() void Spline<DataT, nXdimT, nFdimT, isConsistentT>::
  interpolate(GPUgeneric() const DataT x[], GPUgeneric() DataT S[]) const
{
  /// Get interpolated value for F(x)
  assert(isConsistentT);
  if (isConsistentT) {
    DataT u[nXdimT];
    for (int i = 0; i < nXdimT; i++) {
      u[i] = mGrid[i].convXtoU(x[i]);
    }
    interpolateU(mFparameters, u, S);
  } else {
    for (int i = 0; i < nFdimT; i++) {
      S[i] = 0.;
    }
  }
}

template <typename DataT, int nXdimT, int nFdimT, bool isConsistentT>
GPUhdi() DataT Spline<DataT, nXdimT, nFdimT, isConsistentT>::
  interpolate(GPUgeneric() const DataT x[]) const
{
/// Simplified interface for 1D: get interpolated value for the first dimension of F(x)

#if defined(GPUCA_GPUCODE)
  DataT S[nFdimT]; // constexpr array size for the GPU compiler
#else
  DataT S[getFdimensions()];
#endif
  interpolate(x, S);
  return S[0];
}

#ifdef XXX
template <typename DataT, int nXdimT, int nFdimT, bool isConsistentT>
GPUhdi() void Spline<DataT, nXdimT, nFdimT, isConsistentT>::interpolateU(
  GPUgeneric() const DataT Fparameters[],
  GPUgeneric() const DataT u[], GPUgeneric() DataT S[]) const

/// Get interpolated value S for an nFdim-dimensional F(u) using spline parameters Fparameters.
/// Fparameters can be created via SplineHelper.
/*{//DUMMYVERSION:
  for (int i = 0; i < nFdimT; i++) {
    S[i] = 0.;

  }*/

{ //MY VERSION:
  //std::cout<<"INTERPOLATEU"<<std::endl;

  int nParameters = (int)pow(4.0, nXdimT); //total Nr of Parameters necessary for one interpolation
  //std::cout<<"nParameters "<<nParameters<<std::endl;
  int nKnotParameters = (int)pow(2.0, nXdimT); // Nr of Parameters per Knot
  DataT iParameters[nParameters * nFdimT];     // Array for all parameters

  int nrofInterpolations = ((int)nParameters / 4) * nFdimT;

  /// TO BE REMOVED TESTOUTPUT
  //std::cout<< "numberofInterpolations: " <<nrofInterpolations << std::endl;
  /// END TESTOUTPUT

  //TO BE REMOVED (testing fparameters)
  //std::cout << "Fparameters: " << std::endl;
  //for (int j = 0; j < 128; j++){
  //    std::cout << Fparameters[j]  <<",";
  //  }
  //  std::cout << std::endl;
  //END TESTOUTPUT

  //get the indices of the "most left" Knot:

  int indices[nXdimT]; //indices of the 'most left' knot
  for (int i = 0; i < nXdimT; i++) {
    indices[i] = mGrid[i].getKnotIndexU(u[i]);
  }
  // get all the needed parameters into one array iParameters[nParameters]:
  int indicestmp[nXdimT];
  for (int i = 0; i < (int)(pow(2.0, nXdimT)); i++) { // for every necessary Knot
    for (int k = 0; k < nXdimT; k++) {
      indicestmp[k] = indices[k] + (int)(i / pow(2.0, k)) % 2; //get the knot-indices in every dimension (mirrored order binary counting)
    }
    int index = BaseT::getKnotIndex(indicestmp); //get index of the current Knot

    for (int j = 0; j < nKnotParameters * nFdimT; j++) { //and fill the iparameter array with according parameters
      iParameters[i * nKnotParameters * nFdimT + j] = Fparameters[index * nFdimT * nKnotParameters + j];
    }
  }
  //now start with the interpolation loop:
  for (int d = 0; d < nXdimT; d++) { //for every dimension
    int nrofInterpolations = (int)pow(4.0, nXdimT - d - 1) * nFdimT;
    //std::cout<<"nrofInterpolations = "<<nrofInterpolations<<std::endl;
    int nrofKnots = (int)pow(2.0, nXdimT - d);
    DataT S0[nrofInterpolations];
    DataT D0[nrofInterpolations];
    DataT S1[nrofInterpolations];
    DataT D1[nrofInterpolations];

    DataT* pointer[4] = {S0, D0, S1, D1}; // pointers for interpolation arrays S0, D0, S1, D1 point to Arraystart

    for (int i = 0; i < (int)(pow(2.0, nXdimT - d)); i++) {   //for every knot
      for (int j = 0; j < (int)(pow(2.0, nXdimT - d)); j++) { // for every parametertype
        int pointernr = 2 * (i % 2) + (j % 2);                //to which array should it be delivered

        for (int k = 0; k < nFdimT; k++) {
          pointer[pointernr][0] = iParameters[i * (int)(pow(2.0, nXdimT - d)) * nFdimT + j * nFdimT + k];
          pointer[pointernr]++;
        }
      } // end for j (every parametertype)
    }   // end for i (every knot)

    const typename Spline1D<DataT>::Knot& knotL = mGrid[d].getKnot(indices[d]);
    int Fdim = (int)pow(4.0, nXdimT - d - 1) * nFdimT;
    DataT coordinate = u[d];

    mGrid[d].interpolateU(Fdim, knotL, S0, D0, S1, D1, coordinate, iParameters);

  } //end d (every dimension)

  //std::cout<<std::endl<< " ERGEBNIS INTERPOLATU =" << nFdimT << std::endl;
  for (int i = 0; i < nFdimT; i++) {
    S[i] = iParameters[i]; // write into result-array
    //std::cout<<iParameters[i] <<", ";
  }
  //std::cout<<std::endl;
} // end interpolateU

  // TODO: Implement
  /*
  int nu = mGridU1.getNumberOfKnots();
  int iu = mGridU1.getKnotIndexU(u);
  int iv = mGridU2.getKnotIndexU(v);

  const typename Spline1D<DataT>::Knot& knotU = mGridU1.getKnot(iu);
  const typename Spline1D<DataT>::Knot& knotV = mGridU2.getKnot(iv);

#if defined(GPUCA_GPUCODE) // constexpr array size for the GPU compiler
  const int nFdim = nFdimT;
#else
  const int nFdim = getFdimensions();
#endif

  const int nFdim2 = nFdimT * 2;
  const int nFdim4 = nFdimT * 4;

  // X:=Sx, Y:=Sy, Z:=Sz

  const DataT* par00 = Fparameters + (nu * iv + iu) * nFdim4; // values { {X,Y,Z}, {X,Y,Z}'v, {X,Y,Z}'u, {X,Y,Z}''vu } at {u0, v0}
  const DataT* par10 = par00 + nFdim4;                        // values { ... } at {u1, v0}
  const DataT* par01 = par00 + nFdim4 * nu;                   // values { ... } at {u0, v1}
  const DataT* par11 = par01 + nFdim4;                        // values { ... } at {u1, v1}

  DataT Su0[nFdim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u0
  DataT Du0[nFdim4]; // derivatives {}'_u  at u0
  DataT Su1[nFdim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u1
  DataT Du1[nFdim4]; // derivatives {}'_u  at u1

  for (int i = 0; i < nFdim2; i++) {
    Su0[i] = par00[i];
    Su0[nFdim2 + i] = par01[i];

    Du0[i] = par00[nFdim2 + i];
    Du0[nFdim2 + i] = par01[nFdim2 + i];

    Su1[i] = par10[i];
    Su1[nFdim2 + i] = par11[i];

    Du1[i] = par10[nFdim2 + i];
    Du1[nFdim2 + i] = par11[nFdim2 + i];
  }

  DataT parU[nFdim4]; // interpolated values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) } at u
  mGridU1.interpolateU(nFdim4, knotU, Su0, Du0, Su1, Du1, u, parU);

  const DataT* Sv0 = parU + 0;
  const DataT* Dv0 = parU + nFdim;
  const DataT* Sv1 = parU + nFdim2;
  const DataT* Dv1 = parU + nFdim2 + nFdim;

  mGridU2.interpolateU(nFdim, knotV, Sv0, Dv0, Sv1, Dv1, v, S);
  */
#endif

template <typename DataT, int nXdimT, int nFdimT, bool isConsistentT>
GPUhdi() void Spline<DataT, nXdimT, nFdimT, isConsistentT>::interpolateU(
  GPUgeneric() const DataT Fparameters[],
  GPUgeneric() const DataT u[], GPUgeneric() DataT S[]) const
{ //MY VERSION:

  /// Get interpolated value S for an nFdim-dimensional F(u) using spline parameters Fparameters.
  /// Fparameters can be created via SplineHelper.

  //std::cout<<"INTERPOLATEU"<<std::endl;

  constexpr int nParameters = 1 << (2 * nXdimT); //(int)pow(4.0, nXdimT); //total Nr of Parameters necessary for one interpolation
  //std::cout<<"nParameters "<<nParameters<<std::endl;
  constexpr int nKnotParameters = 1 << nXdimT; //(int)pow(2.0, nXdimT); // Nr of Parameters per Knot
  DataT iParameters[nParameters * nFdimT];     // Array for all parameters

  //constexpr int nrofInterpolations = ((int)nParameters / 4) * nFdimT;

  //get the indices of the "most left" Knot:

  int indices[nXdimT]; //indices of the 'most left' knot
  for (int i = 0; i < nXdimT; i++) {
    indices[i] = mGrid[i].getKnotIndexU(u[i]);
  }
  // get all the needed parameters into one array iParameters[nParameters]:
  int indicestmp[nXdimT];
  for (int i = 0; i < nKnotParameters; i++) { // for every necessary Knot
    for (int k = 0; k < nXdimT; k++) {
      indicestmp[k] = indices[k] + (i / (1 << k)) % 2; //get the knot-indices in every dimension (mirrored order binary counting)
    }
    int index = BaseT::getKnotIndex(indicestmp); //get index of the current Knot

    for (int j = 0; j < nKnotParameters * nFdimT; j++) { //and fill the iparameter array with according parameters
      iParameters[i * nKnotParameters * nFdimT + j] = Fparameters[index * nFdimT * nKnotParameters + j];
    }
  }
  //now start with the interpolation loop:
  int nrofInterpolations = (1 << (2 * nXdimT - 2)) * nFdimT;
  int nrofKnots = 1 << (nXdimT);

  DataT S0[nrofInterpolations];
  DataT D0[nrofInterpolations];
  DataT S1[nrofInterpolations];
  DataT D1[nrofInterpolations];

  for (int d = 0; d < nXdimT; d++) { //for every dimension
    //std::cout<<"nrofInterpolations = "<<nrofInterpolations<<std::endl;

    DataT* pointer[4] = {S0, D0, S1, D1}; // pointers for interpolation arrays S0, D0, S1, D1 point to Arraystart

    for (int i = 0; i < nrofKnots; i++) {      //for every knot
      for (int j = 0; j < nrofKnots; j++) {    // for every parametertype
        int pointernr = 2 * (i % 2) + (j % 2); //to which array should it be delivered
        for (int k = 0; k < nFdimT; k++) {
          pointer[pointernr][0] = iParameters[(i * nrofKnots + j) * nFdimT + k];
          pointer[pointernr]++;
        }
      } // end for j (every parametertype)
    }   // end for i (every knot)

    const typename Spline1D<DataT>::Knot& knotL = mGrid[d].getKnot(indices[d]);
    int Fdim = nrofInterpolations;
    DataT coordinate = u[d];

    mGrid[d].interpolateU(Fdim, knotL, S0, D0, S1, D1, coordinate, iParameters);

    nrofInterpolations = nrofInterpolations / 4;
    nrofKnots = nrofKnots / 2;
  } //end d (every dimension)

  //std::cout<<std::endl<< " ERGEBNIS INTERPOLATU =" << nFdimT << std::endl;
  for (int i = 0; i < nFdimT; i++) {
    S[i] = iParameters[i]; // write into result-array
    //std::cout<<iParameters[i] <<", ";
  }
  //std::cout<<std::endl;
} // end interpolateU

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
