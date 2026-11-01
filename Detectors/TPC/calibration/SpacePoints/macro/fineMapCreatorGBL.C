// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file fineMapCreatorGBL.C
/// \brief Fine level of the TPC space-charge distortion map from the unbinned residuals: per ITS-TPC track a TPC-only GBL fit, the
/// unbiased residuals of its clusters summed per voxel of the map-creation grid, and the voxel solve of the slope model.
///
/// Per track:
/// - the TPC clusters as the residual stream stores them (uncorrected, at the track's time), corrected with the base correction
///   d = d_map + lumiScale * d_derivative (both evaluated at the cluster; lumiScale from the TF's IDC as in the reconstruction,
///   LumiScaleMode::DerivativeMap) and optionally with a previous fine map (the iteration); along the motion (a cluster on a row not
///   beyond the previous one is dropped);
/// - clusters near the stack boundaries, on the first / last rows and at the sector edges (on the pad plane) stay in the fit with their
///   errors x ExcludedSigmaScale: they hardly pull the track, and their unbiased residuals are summed like the others;
/// - reference: the ITS track propagated cluster to cluster with the PropagatorD (field map, LUT energy loss, pion); one GBL point per
///   cluster with the step Jacobian (helix in the mean local Bz, central differences) and a thin scatterer every ScatterEvery points
///   with the x/X0 of its steps (Highland, the log term of the whole track); the seed takes GBL's correction at the first cluster every
///   pass; from FirstOutlierPass on, measurements beyond OutlierSigma of the previous pass's fit are dropped; a track stops once no new
///   outlier appears and GBL's correction is negligible;
/// - per cluster the unbiased residual r m^2 / e^2 (m: its error in the fit, e: GBL's residual error) with the variance sigma^2 + u,
///   u = v m^2 / (m^2 - v), v = m^2 - e^2 (sigma: the cluster's own error); the voxel from the best position y - r sigma^2 / (sigma^2 + u)
///   (cluster and prediction combined: binning by the measured position is Berkson); unbiased pull > OutlierSigma and the outermost
///   SkipEdgePads pads of the row (the best position moved onto the pad plane by the cluster's own correction) not summed.
/// Voxel model ((dX, dY, dZ) = measured - true, as voxResTree): r_y = dY - t dX, r_z = dZ - c dX with t = dy/dx, c = dz/dx of the track;
/// dX, dY from the y channel, then dZ = <r_z> + <c> dX. Errors scaled by sqrt(chi2 / ndf) if > 1.
///
/// Smoothing (smoothTau > 0): a discrete smoothing spline with one node per voxel and component at the bin centre, bilinear in y/x and z/x
/// within the pad row in between (as a fine map is applied), so a voxel's sums enter at their mean position through the bilinear weights.
/// Per voxel with entries the chi2 of its residuals as a quadratic form in (dX, dY, dZ) (the y channel, and the z channel's mean
/// W_z (<r_z> - dZ + <c> dX)^2 as in the voxel solve), divided by max(1, chi2 / ndf) of its solve (the median of the solved voxels for the
/// others); penalties: the second differences along the rows (none across a stack boundary), y/x and z/x with weight 1 / smoothTau^2, and
/// a ridge 1 / smoothTau0^2 towards 0 (no data: no fine correction). Per sector one sparse system, solved with conjugate gradients (3 x 3
/// block-Jacobi preconditioner); smoothRobust > 0: that many reweightings of the voxels with Tukey's bisquare of sqrt(excess chi2 / rank)
/// at the solution (the voxel's chi2 there minus its own minimum).
///
/// Two passes (nPasses = 2, needs the smoothing): the TFs split alternately in time; pass 1 on the even ones, smoothed: S1; pass 2 on the
/// odd ones with S1 subtracted from the clusters (after a fineMapInput), smoothed: S2. Every track is fitted once. Output DS = S1 + S2,
/// D = pass 2's voxel solve + S1 at its mean position (E and the entries pass 2's). Cross-validation (log and user info): the chi2 change
/// of pass 2's residuals with pass 1's raw voxels and with S1, relative to no fine correction.
///
/// ITS clusters and primary vertex in the fit (itsInFit, for tests; ITS-TPC tracks with at least MinITSClusters ITS clusters): the points
/// of the residual stream (reference + residual; ITS in the sensor's tracking frame of the ideal geometry for MC, the aligned one for data)
/// in front of the TPC clusters, a scatterer at every point before the TPC and scatterer points every GapStep in the gap to the TPC, the
/// reference from the vertex (or the innermost ITS cluster) on; only the TPC clusters' unbiased residuals are summed. The track's kinks are
/// then fitted: the reference does not follow them, so the slopes of the voxel model should be the fitted track's (fitSlopes).
///
/// TF time and IDC: the first orbit (records tree) and the CTPLumi of every TF are read first, with the time window cut of
/// staticMapCreatorCPM.C (firstTFTime / lastTFTime: the TFs of the coarse map of that window); the TFs are sorted in time and the TPC
/// scaler objects fetched from the CCDB in that order (each once); the IDC (mean of A and C side) and the lumi scale per TF go into a table.
/// Only this reader thread touches the input files and the CCDB (TGrid is not thread-safe): it reads one TF after the other into memory,
/// the fit threads take chunks of its tracks. With a derivative map, TFs without a valid scaler object are skipped; without one, the map
/// is applied as is (no IDC needed, e.g. the coarse map of the same time window).
///
/// Output: voxResTree (D = (dX, dY, dZ), E, flags DistDone | DispDone where solved; DS = the smoothed map in every voxel with SmoothDone,
/// without smoothing DS = D; binning, CTP / IDC means and medians in the user info) and the OrbitLumiInfo tree as staticMapCreatorCPM.C, with the per-TF IDC values; with writeSums also the tree fineSums
/// (per voxel and charge the sums of the solve, float) for adding outputs.
///
/// The GBL, Mille (GBL >= V04 includes its headers) and Eigen headers and libGBL are not on ROOT's paths in an O2 environment; e.g. in alienv:
///   GBL=$BASEDIR/$(tr ':' '\n' <<< "$LOADEDMODULES" | grep '^GBL/')
///   MILLE=$BASEDIR/$(tr ':' '\n' <<< "$LOADEDMODULES" | grep '^Mille/')
///   EIGEN=$(ls -d $BASEDIR/Eigen3/* | tail -1)
///   export ROOT_INCLUDE_PATH=$GBL/include:$MILLE/include:$EIGEN/include/eigen3:$ROOT_INCLUDE_PATH
///   export LD_LIBRARY_PATH=$GBL/lib:$MILLE/lib:$LD_LIBRARY_PATH
///   root -l -b -q 'fineMapCreatorGBL.C+("files.txt", 568041, "fineVoxRes.root", "corrMap.root", "corrMapDerivative.root")'

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <fstream>
#include <memory>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

#include <sys/resource.h>
#include <unistd.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "GblPoint.h"
#include "GblTrajectory.h"

#include "TFile.h"
#include "TGeoManager.h"
#include "TGrid.h"
#include "TNamed.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "Algorithm/RangeTokenizer.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/LHCConstants.h"
#include "CommonUtils/ConfigurableParam.h"
#include "CorrectionMapsHelper.h"
#include "DataFormatsCTP/LumiInfo.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/MatLayerCylSet.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "ITSBase/GeometryTGeo.h"
#include "MathUtils/Utils.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/TrackParametrization.h"
#include "SpacePoints/TrackInterpolation.h"
#include "SpacePoints/TrackResiduals.h"
#include "TPCBase/Mapper.h"
#include "TPCBaseRecSim/CDBTypes.h"
#include "TPCCalibration/TPCScaler.h"
#include "TPCFastTransform.h"
R__LOAD_LIBRARY(libGBL)
#else
#error This macro must run in compiled mode
#endif

namespace fine_map
{

using GID = o2::dataformats::GlobalTrackID;
using MatCorrType = o2::base::PropagatorD::MatCorrType;
using TrackParD = o2::track::TrackParametrization<double>;
using TrackResiduals = o2::tpc::TrackResiduals;

constexpr int NSectors = 36;
constexpr int NRows = 152;

// track fit, validated at the production grid (20 x 28) on data and MC
constexpr int MaxFitPasses = 7;            ///< at most; a track stops once converged
constexpr int FirstOutlierPass = 2;        ///< from this pass on, the previous pass's outliers are dropped
constexpr double OutlierSigma = 5.;        ///< outliers in the fit and of the unbiased residuals
constexpr int MinMeasurements = 30;        ///< clusters reached by the reference
constexpr double MaxSnp = 0.98;            ///< clusters beyond are dropped from the fit
constexpr int ScatterEvery = 5;            ///< a scatterer every n points with the x/X0 of its n steps
constexpr int MinFitClusters = 60;         ///< clusters not de-weighted
constexpr double ExcludedSigmaScale = 10.; ///< errors of the de-weighted clusters in the fit
// clusters de-weighted in the fit, on the pad plane
constexpr int StackRows = 2;     ///< rows on each side of a stack boundary
constexpr int EdgeRows = 4;      ///< the first / last rows
constexpr double EdgePads = 2.;  ///< pads from the pad plane's edge
constexpr double EdgeY2X = 0.14; ///< |y/x| above
// residuals not summed
constexpr double SkipEdgePads = 1.; ///< the outermost pads of the row, at the best position
// threading
constexpr size_t TFQueueDepth = 2;       ///< TFs read ahead of the fit
constexpr size_t TracksPerChunk = 16;    ///< tracks a fit thread takes at once
constexpr size_t SumsBatchSize = 100000; ///< contributions of a fit thread added to the shared sums at once
// ITS clusters and primary vertex in the track fit (itsInFit)
constexpr int RowVertex = 190;                         ///< the primary vertex in the residual stream (UnbinnedResid::isITS() is true for it too)
constexpr int MinITSClusters = 3;                      ///< for a track to be fitted
constexpr double ITSSigmaY = 0.0012;                   ///< [cm] ITS cluster errors
constexpr double ITSSigmaZ = 0.0036;                   ///< [cm]
constexpr double VertexSigmaY = 0.002;                 ///< [cm] primary vertex errors
constexpr double VertexSigmaZ = 0.005;                 ///< [cm]
constexpr double GapStep = 5.;                         ///< [cm] scatterers in the gap from the ITS to the TPC: one per step of at most this length
constexpr long IdealGeometryTimeStamp = 1704063600000; ///< GLO/Config/Geometry of MC (the ideal geometry)
// smoothing
constexpr double SmoothTukey = 4.685;        ///< Tukey's bisquare constant of the robust reweighting
constexpr double SmoothCGTolerance = 1e-8;   ///< relative residual norm of the conjugate gradients
constexpr int SmoothCGMaxIterations = 20000; ///< per sector and cycle

/// resident memory of the process [GB]: the current one where /proc/self/statm exists (Linux), else the peak (getrusage)
double residentMemoryGB()
{
  std::ifstream statm("/proc/self/statm");
  long size = 0;
  long resident = 0;
  if (statm >> size >> resident) {
    return resident * double(sysconf(_SC_PAGESIZE)) / 1e9;
  }
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
#ifdef __APPLE__
  return usage.ru_maxrss / 1e9;
#else
  return usage.ru_maxrss / 1e6;
#endif
}

template <typename T>
double calculateMean(const std::vector<T>& values)
{
  return values.empty() ? 0. : std::accumulate(values.begin(), values.end(), 0.) / values.size();
}

template <typename T>
double calculateMedian(std::vector<T> values)
{
  if (values.empty()) {
    return 0.;
  }
  const size_t middle = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + middle, values.end());
  if (values.size() % 2) {
    return values[middle];
  }
  const T upper = values[middle];
  std::nth_element(values.begin(), values.begin() + middle - 1, values.end());
  return (values[middle - 1] + upper) / 2.;
}

/// one residual file, or a .txt file with one per line
std::vector<std::string> getInputFileList(const std::string& fileInput)
{
  std::vector<std::string> files;
  if (fileInput.size() > 4 && fileInput.substr(fileInput.size() - 4) == ".txt") {
    std::ifstream list(fileInput);
    for (std::string line; std::getline(list, line);) {
      if (!line.empty() && line[0] != '#') {
        files.push_back(line);
      }
    }
  } else {
    files.push_back(fileInput);
  }
  return files;
}

/// the map-creation binning: one number = uniform, else the bin boundaries (scaled by maxZ2X in z/x)
std::unique_ptr<TrackResiduals> makeBinning(const std::string& y2xBinning, const std::string& z2xBinning, float maxZ2X)
{
  o2::conf::ConfigurableParam::setValue("scdcalib.maxZ2X", std::to_string(maxZ2X));
  auto binning = std::make_unique<TrackResiduals>();
  binning->setY2XBinning(o2::RangeTokenizer::tokenize<float>(y2xBinning));
  binning->setZ2XBinning(o2::RangeTokenizer::tokenize<float>(z2xBinning));
  binning->init();
  return binning;
}

/// pad-row x, pad width and the pad plane's edge (|y| of the outer edge of the outermost pad) per row
struct RowGeometry {
  std::array<double, NRows> x{};
  std::array<double, NRows> padWidth{};
  std::array<double, NRows> padPlaneEdge{};
  RowGeometry()
  {
    const auto& mapper = o2::tpc::Mapper::instance();
    for (int row = 0; row < NRows; ++row) {
      const int nPads = mapper.getNumberOfPadsInRowSector(row);
      x[row] = mapper.getPadCentre(o2::tpc::PadPos(row, 0)).X();
      padWidth[row] = mapper.getPadRegionInfo(o2::tpc::Mapper::REGION[row]).getPadWidth();
      padPlaneEdge[row] = std::max(std::abs(mapper.getPadCentre(o2::tpc::PadPos(row, 0)).Y()),
                                   std::abs(mapper.getPadCentre(o2::tpc::PadPos(row, nPads - 1)).Y())) +
                          0.5 * padWidth[row];
    }
  }
};

/// the clusters de-weighted in the track fit: stack boundaries (first rows of OROC1 / 2 / 3) and the first / last rows, within EdgePads of
/// the pad plane's edge, |y/x| above EdgeY2X
bool isExcludedFromFit(int row, double padY, const RowGeometry& geometry)
{
  for (int boundary : {63, 97, 127}) {
    if (row >= boundary - StackRows && row < boundary + StackRows) {
      return true;
    }
  }
  if (row < EdgeRows || row >= NRows - EdgeRows) {
    return true;
  }
  const double absY = std::abs(padY);
  return absY > geometry.padPlaneEdge[row] - EdgePads * geometry.padWidth[row] || absY > EdgeY2X * geometry.x[row];
}

/// TPC cluster errors (GPUParam's parametrisation)
void clusterErrors(int row, double z, double tgSlp, double tgl, double& sigmaY, double& sigmaZ)
{
  constexpr double ErrorsY[3] = {0.06, 0.24, 0.12};
  constexpr double ErrorsZ[3] = {0.06, 0.24, 0.15};
  constexpr double PadHeight[4] = {0.75, 1.0, 1.2, 1.5};
  constexpr double DriftLength = 250.;
  const double padHeight = PadHeight[row < 63 ? 0 : (row < 97 ? 1 : (row < 127 ? 2 : 3))];
  const double drift = DriftLength - std::abs(z);
  const double angleY2 = tgSlp * tgSlp;
  const double angleZ2 = tgl * tgl * (1. + tgSlp * tgSlp);
  sigmaY = std::sqrt(ErrorsY[0] * ErrorsY[0] + ErrorsY[1] * ErrorsY[1] * padHeight * angleY2 + ErrorsY[2] * ErrorsY[2] / DriftLength / padHeight * drift);
  sigmaZ = std::sqrt(ErrorsZ[0] * ErrorsZ[0] + ErrorsZ[1] * ErrorsZ[1] * padHeight * angleZ2 + ErrorsZ[2] * ErrorsZ[2] / DriftLength / padHeight * drift);
}

/// the stack (IROC, OROC1 / 2 / 3) of a pad row
int stackOfRow(int row)
{
  return row < 63 ? 0 : (row < 97 ? 1 : (row < 127 ? 2 : 3));
}

/// the bilinear interpolation between the bin centres of the pad row in y/x and the side-oriented z/x (z/x on the A side, -z/x on the C
/// side; constant beyond the outermost centres): up to 4 voxels (sector-local index) and their weights; returns their number. A fine
/// map is applied with it, and it defines the smoothed map between its nodes.
int bilinearNodes(const TrackResiduals& binning, int row, double y2x, double z2xSide, int voxels[4], double weights[4])
{
  auto bracket = [](int nBins, double value, auto centre, int& low, double& fraction) {
    low = 0;
    fraction = 0.;
    if (nBins < 2 || value <= centre(0)) {
      return;
    }
    if (value >= centre(nBins - 1)) {
      low = nBins - 1;
      return;
    }
    while (low < nBins - 2 && value >= centre(low + 1)) {
      ++low;
    }
    fraction = (value - centre(low)) / (centre(low + 1) - centre(low));
  };
  int binY2X = 0;
  int binZ2X = 0;
  double fractionY2X = 0.;
  double fractionZ2X = 0.;
  bracket(binning.getNY2XBins(), y2x, [&](int bin) { return double(binning.getY2X(row, bin)); }, binY2X, fractionY2X);
  bracket(binning.getNZ2XBins(), z2xSide, [&](int bin) { return double(binning.getZ2X(bin)); }, binZ2X, fractionZ2X);
  int nVoxels = 0;
  for (int stepY2X = 0; stepY2X < 2; ++stepY2X) {
    for (int stepZ2X = 0; stepZ2X < 2; ++stepZ2X) {
      const double weight = (stepY2X ? fractionY2X : 1. - fractionY2X) * (stepZ2X ? fractionZ2X : 1. - fractionZ2X);
      if (weight > 0.) {
        voxels[nVoxels] = int(binning.getGlbVoxBin(row, binY2X + stepY2X, binZ2X + stepZ2X));
        weights[nVoxels++] = weight;
      }
    }
  }
  return nVoxels;
}

/// a previous pass's voxResTree (DS where smoothed, else D where fitted, else 0), bilinear between the bin centres (bilinearNodes)
struct FineMapInput {
  std::unique_ptr<TrackResiduals> binning;
  int voxelsPerSector = 0;
  std::vector<float> nodes;

  bool load(const std::string& fileName)
  {
    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str()));
    auto* tree = file ? file->Get<TTree>("voxResTree") : nullptr;
    auto info = [&](const char* key) -> std::string {
      auto* object = tree ? tree->GetUserInfo()->FindObject(key) : nullptr;
      return object ? object->GetTitle() : "";
    };
    if (!tree || info("y2xBinning").empty() || info("z2xBinning").empty() || info("maxZ2X").empty()) {
      LOGP(error, "Fine map {}: no voxResTree with the binning in its user info", fileName);
      return false;
    }
    binning = makeBinning(info("y2xBinning"), info("z2xBinning"), std::stof(info("maxZ2X")));
    voxelsPerSector = binning->getNVoxelsPerSector();
    TrackResiduals::VoxRes* result = nullptr;
    tree->SetBranchAddress("voxRes", &result);
    nodes.assign(size_t(NSectors) * voxelsPerSector * 3, 0.f);
    for (long entry = 0; entry < tree->GetEntries(); ++entry) {
      tree->GetEntry(entry);
      const size_t voxel = size_t(result->bsec) * voxelsPerSector + binning->getGlbVoxBin(result->bvox);
      const bool smoothed = result->flags & TrackResiduals::SmoothDone;
      const bool fitted = result->flags & TrackResiduals::DistDone;
      for (int coordinate = 0; coordinate < 3; ++coordinate) {
        nodes[3 * voxel + coordinate] = smoothed ? result->DS[coordinate] : (fitted ? result->D[coordinate] : 0.f);
      }
    }
    tree->ResetBranchAddresses();
    LOGP(info, "Fine map {} is subtracted from the clusters", fileName);
    return true;
  }

  void apply(int sector, int row, double& x, double& y, double& z) const
  {
    int voxels[4];
    double weights[4];
    const int nVoxels = bilinearNodes(*binning, row, y / x, (sector < 18 ? z : -z) / x, voxels, weights);
    double correction[3] = {0., 0., 0.};
    for (int node = 0; node < nVoxels; ++node) {
      const size_t voxel = size_t(sector) * voxelsPerSector + voxels[node];
      for (int coordinate = 0; coordinate < 3; ++coordinate) {
        correction[coordinate] += weights[node] * nodes[3 * voxel + coordinate];
      }
    }
    x -= correction[0];
    y -= correction[1];
    z -= correction[2];
  }
};

/// the base correction at the cluster's nominal position: d_map + lumiScale * d_derivative, both evaluated there (the derivative-map
/// merging of the reconstruction, without resampling the derivative onto the map's knots)
struct BaseCorrection {
  std::unique_ptr<o2::gpu::TPCFastTransform> map;
  std::unique_ptr<o2::gpu::TPCFastTransform> derivative;

  void apply(int sector, int row, float lumiScale, float& x, float& y, float& z) const
  {
    float dx = 0.f;
    float dy = 0.f;
    float dz = 0.f;
    float derivativeX = 0.f;
    float derivativeY = 0.f;
    float derivativeZ = 0.f;
    map->getCorrection().getCorrectionLocal(sector, row, y, z, dx, dy, dz);
    if (derivative) {
      derivative->getCorrection().getCorrectionLocal(sector, row, y, z, derivativeX, derivativeY, derivativeZ);
    }
    x += dx + lumiScale * derivativeX;
    y += dy + lumiScale * derivativeY;
    z += dz + lumiScale * derivativeZ;
  }
};

/// a TPCFastTransform from a file (object "ccdb_object") or from the CCDB at timestamp: "ccdb" = defaultPath, "ccdb:<path>"
std::unique_ptr<o2::gpu::TPCFastTransform> loadTransform(const std::string& source, const std::string& defaultPath, long timestamp)
{
  if (source == "ccdb" || source.rfind("ccdb:", 0) == 0) {
    const std::string path = source == "ccdb" ? defaultPath : source.substr(5);
    auto* object = o2::ccdb::BasicCCDBManager::instance().getForTimeStamp<o2::gpu::TPCFastTransform>(path, timestamp);
    if (!object) {
      LOGP(fatal, "No {} in the CCDB at {}", path, timestamp);
    }
    object->rectifyAfterReadingFromFile();
    auto transform = std::make_unique<o2::gpu::TPCFastTransform>();
    transform->cloneFromObject(*object, nullptr);
    LOGP(info, "Correction {} at {} from the CCDB, IDC {}", path, timestamp, transform->getIDC());
    return transform;
  }
  std::unique_ptr<TFile> file(TFile::Open(source.c_str()));
  std::unique_ptr<o2::gpu::TPCFastTransform> transform((file && !file->IsZombie()) ? file->Get<o2::gpu::TPCFastTransform>("ccdb_object") : nullptr);
  if (!transform) {
    LOGP(fatal, "No TPCFastTransform ccdb_object in {}", source);
  }
  transform->rectifyAfterReadingFromFile();
  LOGP(info, "Correction {}, IDC {}", source, transform->getIDC());
  return transform;
}

/// one TF of the input
struct TFInfo {
  int file = 0;
  long entry = 0;
  uint32_t orbit = 0; ///< CTPLumi orbit
  float ctpLumi = 0.f;
  double timeMS = 0.;
  bool hasIDC = false;
  float idcA = 0.f;
  float idcC = 0.f;
  float lumiScale = 0.f;
};

/// the TFs of the input files with their time as staticMapCreatorCPM.C (orbit reset + firstTForbit of the records tree; in data the
/// CTPLumi orbit is firstTForbit + 1, in MC equal to it) and their CTPLumi orbit and luminosity; isMC: a file has the MC truth (trkMC)
std::vector<TFInfo> scanTFs(const std::vector<std::string>& files, long firstEntry, long entriesPerFile, double orbitResetMS, long firstTFTime, long lastTFTime,
                            long& nOutsideWindow, bool& isMC)
{
  std::vector<TFInfo> tfs;
  for (size_t fileIndex = 0; fileIndex < files.size(); ++fileIndex) {
    std::unique_ptr<TFile> file(TFile::Open(files[fileIndex].c_str()));
    auto* residualTree = (file && !file->IsZombie()) ? file->Get<TTree>("unbinnedResid") : nullptr;
    auto* trackTree = (file && !file->IsZombie()) ? file->Get<TTree>("trackData") : nullptr;
    auto* recordsTree = (file && !file->IsZombie()) ? file->Get<TTree>("records") : nullptr;
    if (!residualTree || !trackTree || !recordsTree || residualTree->GetEntries() != trackTree->GetEntries()) {
      LOGP(warning, "Skipping {}: no unbinnedResid / trackData / records or different entries", files[fileIndex]);
      continue;
    }
    isMC |= trackTree->GetBranch("trkMC") != nullptr;
    // the first orbits of the TFs, in TF order (MC files can have several records entries, one per merged chunk)
    std::vector<uint32_t> firstOrbits;
    TTreeReader recordsReader(recordsTree);
    TTreeReaderValue<std::vector<uint32_t>> recordOrbits(recordsReader, "firstTForbit");
    while (recordsReader.Next()) {
      firstOrbits.insert(firstOrbits.end(), recordOrbits->begin(), recordOrbits->end());
    }
    if (long(firstOrbits.size()) < residualTree->GetEntries()) {
      LOGP(warning, "Skipping {}: firstTForbit has fewer entries than the residual tree ({} vs {})", files[fileIndex], firstOrbits.size(), residualTree->GetEntries());
      continue;
    }
    const long endEntry = entriesPerFile > 0 ? std::min<long>(firstEntry + entriesPerFile, residualTree->GetEntries()) : residualTree->GetEntries();
    TTreeReader reader(residualTree);
    TTreeReaderValue<o2::ctp::LumiInfo> lumi(reader, "CTPLumi");
    for (long entry = firstEntry; entry < endEntry; ++entry) {
      // the time window with the integer arithmetic of staticMapCreatorCPM.C: the same TFs as the coarse map of the window
      const int64_t orbitResetTimeMS = static_cast<int64_t>(orbitResetMS);
      const int64_t tfTimeInMS = static_cast<int64_t>(orbitResetTimeMS + firstOrbits[entry] * o2::constants::lhc::LHCOrbitMUS * 1.e-3);
      if (firstTFTime > 0 && (tfTimeInMS < firstTFTime || tfTimeInMS > lastTFTime)) {
        ++nOutsideWindow;
        continue;
      }
      if (reader.SetEntry(entry) != TTreeReader::kEntryValid) {
        LOGP(warning, "Skipping entry {} of {}: cannot be read", entry, files[fileIndex]);
        continue;
      }
      const auto& lumiInfo = *lumi;
      if (lumi.GetReadStatus() != ROOT::Internal::TTreeReaderValueBase::kReadSuccess) {
        LOGP(warning, "Skipping entry {} of {}: no CTPLumi", entry, files[fileIndex]);
        continue;
      }
      TFInfo tf;
      tf.file = int(fileIndex);
      tf.entry = entry;
      tf.orbit = lumiInfo.orbit;
      tf.ctpLumi = lumiInfo.getLumi();
      tf.timeMS = orbitResetMS + double(firstOrbits[entry]) * o2::constants::lhc::LHCOrbitMUS * 1e-3;
      tfs.push_back(tf);
    }
  }
  return tfs;
}

/// the IDC per TF: the TFs in time order, a TPC/Calib/Scaler object fetched only when the current one does not cover the TF (each object
/// once); weights and integration time as the reconstruction's TPCScalerSpec (weights from TPC/Calib/ScalerWeights, at most 500 ms)
void fillIDCs(std::vector<TFInfo>& tfs, const o2::tpc::TPCScalerWeights* weights)
{
  std::vector<size_t> order(tfs.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](size_t first, size_t second) { return tfs[first].timeMS < tfs[second].timeMS; });
  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();
  std::unique_ptr<o2::tpc::TPCScaler> scaler;
  double retryAfterMS = 0.;
  long nObjects = 0;
  for (const size_t index : order) {
    auto& tf = tfs[index];
    const bool inside = scaler && tf.timeMS >= scaler->getStartTimeStampMS() && tf.timeMS <= scaler->getEndTimeStampMS(o2::tpc::Side::A);
    if (!inside && tf.timeMS >= retryAfterMS) {
      auto* tree = ccdb.getForTimeStamp<TTree>(o2::tpc::CDBTypeMap.at(o2::tpc::CDBType::CalScaler), long(tf.timeMS));
      auto next = std::make_unique<o2::tpc::TPCScaler>();
      if (tree) {
        next->setFromTree(*tree);
      }
      // a default object (no scaler at this time, e.g. after the end of the run) starts at 0: keep the current one, try again 10 s later
      if (!tree || next->getStartTimeStampMS() <= 0) {
        retryAfterMS = tf.timeMS + 1e4;
      } else {
        if (weights && weights->isValid()) {
          next->setScalerWeights(*weights);
          next->useWeights(true);
          const float integrationTimeMS = std::min(std::abs(weights->mFirstTimeStampMS), 500.f);
          if (integrationTimeMS > 0) {
            next->setIonDriftTimeMS(integrationTimeMS);
          }
        }
        scaler = std::move(next);
        ++nObjects;
        LOGP(info, "TPC scaler {:.0f} - {:.0f} ms (TF time {:.0f} ms), ion drift time {} ms", scaler->getStartTimeStampMS(),
             scaler->getEndTimeStampMS(o2::tpc::Side::A), tf.timeMS, scaler->getIonDriftTimeMS());
      }
    }
    if (scaler && tf.timeMS >= scaler->getStartTimeStampMS() - 1e3 && tf.timeMS <= scaler->getEndTimeStampMS(o2::tpc::Side::A) + 1e3) {
      tf.hasIDC = true;
      tf.idcA = scaler->getMeanScaler(tf.timeMS, o2::tpc::Side::A);
      tf.idcC = scaler->getMeanScaler(tf.timeMS, o2::tpc::Side::C);
    }
  }
  const long nWithIDC = std::count_if(tfs.begin(), tfs.end(), [](const TFInfo& tf) { return tf.hasIDC; });
  LOGP(info, "IDC for {} of {} TFs from {} scaler objects", nWithIDC, tfs.size(), nObjects);
}

/// one TF's tracks and residuals (copied out of the readers) and the next of its tracks to deal (guarded by the queue's lock)
struct TFTask {
  const TFInfo* info = nullptr;
  std::vector<o2::tpc::TrackDataCompact> trackRefs;
  std::vector<o2::tpc::TrackData> tracks;
  std::vector<o2::tpc::UnbinnedResid> residuals;
  size_t nextTrack = 0;
};

/// the TFs read and not yet dealt completely: a fit thread takes the next chunk of tracks of the oldest TF; a TF is released with its
/// last chunk. The reader waits while TFQueueDepth TFs are queued.
class TFQueue
{
 public:
  void push(std::shared_ptr<TFTask> task)
  {
    std::unique_lock<std::mutex> lock(mMutex);
    mNotFull.wait(lock, [this] { return mTasks.size() < TFQueueDepth; });
    mTasks.push_back(std::move(task));
    lock.unlock();
    mNotEmpty.notify_all();
  }

  void setDone()
  {
    {
      std::lock_guard<std::mutex> lock(mMutex);
      mDone = true;
    }
    mNotEmpty.notify_all();
  }

  /// tracks [firstTrack, lastTrack) of task; false once the reader is done and every TF is dealt
  bool next(std::shared_ptr<TFTask>& task, size_t& firstTrack, size_t& lastTrack)
  {
    std::unique_lock<std::mutex> lock(mMutex);
    mNotEmpty.wait(lock, [this] { return !mTasks.empty() || mDone; });
    if (mTasks.empty()) {
      return false;
    }
    task = mTasks.front();
    firstTrack = task->nextTrack;
    lastTrack = std::min(firstTrack + TracksPerChunk, task->trackRefs.size());
    task->nextTrack = lastTrack;
    if (lastTrack >= task->trackRefs.size()) {
      mTasks.pop_front();
      lock.unlock();
      mNotFull.notify_one();
    }
    return true;
  }

 private:
  std::mutex mMutex;
  std::condition_variable mNotEmpty;
  std::condition_variable mNotFull;
  std::deque<std::shared_ptr<TFTask>> mTasks;
  bool mDone = false;
};

/// one summed residual
struct Contribution {
  size_t voxel = 0;
  int charge = 0;
  double x = 0.;
  double y2x = 0.;
  double z2x = 0.;
  double slopeY = 0.;
  double residualY = 0.;
  double weightY = 0.;
  double slopeZ = 0.;
  double residualZ = 0.;
  double weightZ = 0.;
};

/// per voxel and charge: entries, position sums, and per channel (y, z) sum w, w s, w s^2, w r, w s r, w r^2 (s the track slope)
struct ChargeSums {
  enum Sum { W,
             WS,
             WSS,
             WR,
             WSR,
             WRR,
             NSums };
  double entries = 0.;
  double x = 0.;
  double y2x = 0.;
  double z2x = 0.;
  std::array<double, NSums> y{};
  std::array<double, NSums> z{};

  static void fill(std::array<double, NSums>& channel, double slope, double residual, double weight)
  {
    channel[W] += weight;
    channel[WS] += weight * slope;
    channel[WSS] += weight * slope * slope;
    channel[WR] += weight * residual;
    channel[WSR] += weight * slope * residual;
    channel[WRR] += weight * residual * residual;
  }

  void add(const Contribution& contribution)
  {
    entries += 1.;
    x += contribution.x;
    y2x += contribution.y2x;
    z2x += contribution.z2x;
    fill(y, contribution.slopeY, contribution.residualY, contribution.weightY);
    fill(z, contribution.slopeZ, contribution.residualZ, contribution.weightZ);
  }

  void add(const ChargeSums& other)
  {
    entries += other.entries;
    x += other.x;
    y2x += other.y2x;
    z2x += other.z2x;
    for (int sum = 0; sum < NSums; ++sum) {
      y[sum] += other.y[sum];
      z[sum] += other.z[sum];
    }
  }
};

struct VoxelSolution {
  std::array<double, 3> value{}; ///< dX, dY, dZ
  std::array<double, 3> error{}; ///< scaled by sqrt(chi2Scale)
  double correlationXY = 0.;
  double rmsY = 0.;
  double rmsZ = 0.;
  double chi2Scale = 1.; ///< max(1, chi2 / ndf)
};

/// the voxel solve: dX, dY from the y channel (2 x 2), then dZ = <r_z> + <c> dX; false below minEntries or singular
bool solveVoxel(const ChargeSums& sums, double minEntries, VoxelSolution& solution)
{
  using S = ChargeSums;
  if (sums.entries < minEntries || sums.y[S::W] <= 0. || sums.z[S::W] <= 0.) {
    return false;
  }
  Eigen::Matrix2d normalY;
  normalY << sums.y[S::WSS], -sums.y[S::WS], -sums.y[S::WS], sums.y[S::W];
  if (!(normalY.determinant() > 1e-12 * sums.y[S::WSS] * sums.y[S::W])) {
    return false;
  }
  const Eigen::Matrix2d covarianceY = normalY.inverse();
  const Eigen::Vector2d solutionY = covarianceY * Eigen::Vector2d(-sums.y[S::WSR], sums.y[S::WR]);
  const double meanSlopeZ = sums.z[S::WS] / sums.z[S::W];
  const double dX = solutionY(0);
  const double dY = solutionY(1);
  const double dZ = sums.z[S::WR] / sums.z[S::W] + meanSlopeZ * dX;
  const double chi2Y = sums.y[S::WRR] + dY * dY * sums.y[S::W] + dX * dX * sums.y[S::WSS] - 2. * dY * sums.y[S::WR] + 2. * dX * sums.y[S::WSR] - 2. * dX * dY * sums.y[S::WS];
  const double chi2Z = sums.z[S::WRR] + dZ * dZ * sums.z[S::W] + dX * dX * sums.z[S::WSS] - 2. * dZ * sums.z[S::WR] + 2. * dX * sums.z[S::WSR] - 2. * dX * dZ * sums.z[S::WS];
  const int ndf = int(2. * sums.entries + 0.5) - 3;
  const double chi2 = std::max(0., chi2Y) + std::max(0., chi2Z);
  const double scale = ndf > 0 ? std::max(1., chi2 / ndf) : 1.;
  solution.value = {dX, dY, dZ};
  solution.error[0] = std::sqrt(std::max(0., covarianceY(0, 0)) * scale);
  solution.error[1] = std::sqrt(std::max(0., covarianceY(1, 1)) * scale);
  solution.error[2] = std::sqrt(std::max(0., 1. / sums.z[S::W] + meanSlopeZ * meanSlopeZ * covarianceY(0, 0)) * scale);
  solution.correlationXY = covarianceY(0, 0) * covarianceY(1, 1) > 0. ? covarianceY(0, 1) / std::sqrt(covarianceY(0, 0) * covarianceY(1, 1)) : 0.;
  solution.rmsY = std::sqrt(std::max(0., chi2Y) / sums.y[S::W]);
  solution.rmsZ = std::sqrt(std::max(0., chi2Z) / sums.z[S::W]);
  solution.chi2Scale = scale;
  return std::isfinite(dX) && std::isfinite(dY) && std::isfinite(dZ);
}

/// one voxel with entries in the smoothing: the chi2 of its residuals theta^T A theta - 2 b^T theta + const at its mean position
/// (theta = (dX, dY, dZ)), divided by its chi2 scale, and the nodes there
struct SmoothVoxel {
  Eigen::Matrix3d a;
  Eigen::Vector3d b;
  double weight = 1.; ///< robust weight
  int nNodes = 0;
  std::array<int, 4> node{}; ///< sector-local
  std::array<double, 4> nodeWeight{};
};

/// the smoothing (see the file description): per sector the nodes (dX, dY, dZ) of every voxel, filled into nodes (3 per voxel); false if
/// a sector's solution is not finite
bool smoothMap(const TrackResiduals& binning, const std::vector<std::array<ChargeSums, 2>>& sums, double minEntries, double tau, double tau0,
               int robust, int nThreads, std::vector<double>& nodes)
{
  using S = ChargeSums;
  const int nBinsX = binning.getNXBins();
  const int nBinsY2X = binning.getNY2XBins();
  const int nBinsZ2X = binning.getNZ2XBins();
  const int voxelsPerSector = binning.getNVoxelsPerSector();
  const int bins[3] = {nBinsX, nBinsY2X, nBinsZ2X};
  const int strides[3] = {1, nBinsX, nBinsX * nBinsY2X}; // getGlbVoxBin = row + nX * (y/x bin + nY2X * z/x bin)
  const double lambda = 1. / (tau * tau);
  const double ridge = 1. / (tau0 * tau0);
  constexpr double Coefficients[3] = {1., -2., 1.};

  // the chi2 scale of the voxels without a solve: the median of the solved ones
  std::vector<double> scales;
  for (const auto& voxelSums : sums) {
    ChargeSums total = voxelSums[0];
    total.add(voxelSums[1]);
    VoxelSolution solution;
    if (solveVoxel(total, minEntries, solution)) {
      scales.push_back(solution.chi2Scale);
    }
  }
  double medianScale = 1.;
  if (!scales.empty()) {
    std::nth_element(scales.begin(), scales.begin() + scales.size() / 2, scales.end());
    medianScale = scales[scales.size() / 2];
  }

  // the second differences, the same in every sector: per stencil the axis and its first node; the penalty's diagonal per node
  std::vector<std::pair<int, int>> stencils;
  std::vector<double> penaltyDiagonal(voxelsPerSector, 0.);
  for (int axis = 0; axis < 3; ++axis) {
    for (int binZ2X = 0; binZ2X < nBinsZ2X; ++binZ2X) {
      for (int binY2X = 0; binY2X < nBinsY2X; ++binY2X) {
        for (int binX = 0; binX < nBinsX; ++binX) {
          const int bin[3] = {binX, binY2X, binZ2X};
          if (bin[axis] + 3 > bins[axis] || (axis == 0 && stackOfRow(binX) != stackOfRow(binX + 2))) {
            continue;
          }
          const int first = int(binning.getGlbVoxBin(binX, binY2X, binZ2X));
          stencils.emplace_back(axis, first);
          for (int k = 0; k < 3; ++k) {
            penaltyDiagonal[first + k * strides[axis]] += lambda * Coefficients[k] * Coefficients[k];
          }
        }
      }
    }
  }

  nodes.assign(sums.size() * 3, 0.);
  std::array<long, NSectors> nVoxels{};
  std::array<int, NSectors> iterations{};
  std::array<double, NSectors> residuals{};
  std::array<bool, NSectors> finite{};
  std::array<long, NSectors> nLowWeight{};
  std::array<long, NSectors> nZeroWeight{};
  auto smoothSector = [&](int sector) {
    finite[sector] = true;
    std::vector<SmoothVoxel> voxels;
    for (int local = 0; local < voxelsPerSector; ++local) {
      const auto& voxelSums = sums[size_t(sector) * voxelsPerSector + local];
      ChargeSums total = voxelSums[0];
      total.add(voxelSums[1]);
      if (total.entries <= 0. || total.y[S::W] <= 0. || total.z[S::W] <= 0.) {
        continue;
      }
      VoxelSolution solution;
      const double scale = solveVoxel(total, minEntries, solution) ? solution.chi2Scale : medianScale;
      // the y channel (J = (-t, 1, 0)) and the z channel's mean W_z (<r_z> - dZ + <c> dX)^2
      const double meanSlopeZ = total.z[S::WS] / total.z[S::W];
      SmoothVoxel voxel;
      voxel.a << total.y[S::WSS] + total.z[S::WS] * meanSlopeZ, -total.y[S::WS], -total.z[S::WS], -total.y[S::WS], total.y[S::W], 0., -total.z[S::WS], 0.,
        total.z[S::W];
      voxel.b << -total.y[S::WSR] - meanSlopeZ * total.z[S::WR], total.y[S::WR], total.z[S::WR];
      voxel.a /= scale;
      voxel.b /= scale;
      int nodeIndices[4];
      double nodeWeights[4];
      voxel.nNodes = bilinearNodes(binning, local % nBinsX, total.y2x / total.entries, (sector < 18 ? 1. : -1.) * total.z2x / total.entries, nodeIndices,
                                   nodeWeights);
      for (int j = 0; j < voxel.nNodes; ++j) {
        voxel.node[j] = nodeIndices[j];
        voxel.nodeWeight[j] = nodeWeights[j];
      }
      voxels.push_back(voxel);
    }
    nVoxels[sector] = voxels.size();
    if (voxels.empty()) {
      return;
    }
    const int n = 3 * voxelsPerSector;
    double* sectorNodes = nodes.data() + 3 * size_t(sector) * voxelsPerSector;
    auto nodeValue = [&](const SmoothVoxel& voxel, const double* x) {
      Eigen::Vector3d theta = Eigen::Vector3d::Zero();
      for (int j = 0; j < voxel.nNodes; ++j) {
        theta += voxel.nodeWeight[j] * Eigen::Map<const Eigen::Vector3d>(x + 3 * voxel.node[j]);
      }
      return theta;
    };
    for (int cycle = 0; cycle <= robust; ++cycle) {
      if (cycle > 0) {
        // Tukey's bisquare of sqrt(excess chi2 / rank): the voxel's chi2 at the solution minus its minimum -b^T A^+ b
        nLowWeight[sector] = 0;
        nZeroWeight[sector] = 0;
        for (auto& voxel : voxels) {
          const Eigen::Vector3d theta = nodeValue(voxel, sectorNodes);
          Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(voxel.a);
          const double largest = eigen.eigenvalues().maxCoeff();
          double minimum = 0.;
          int rank = 0;
          for (int k = 0; k < 3; ++k) {
            const double value = eigen.eigenvalues()(k);
            if (value > 1e-12 * largest) {
              const double projection = eigen.eigenvectors().col(k).dot(voxel.b);
              minimum -= projection * projection / value;
              ++rank;
            }
          }
          const double excess = theta.dot(voxel.a * theta) - 2. * voxel.b.dot(theta) - minimum;
          const double normalised = rank > 0 ? std::sqrt(std::max(0., excess) / rank) / SmoothTukey : 0.;
          voxel.weight = normalised < 1. ? (1. - normalised * normalised) * (1. - normalised * normalised) : 0.;
          nLowWeight[sector] += voxel.weight < 0.5;
          nZeroWeight[sector] += voxel.weight == 0.;
        }
      }
      // the normal equations: ridge + penalty + sum over the voxels of w_j w_k A at their nodes; preconditioner = the inverse 3 x 3 diagonal blocks
      Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
      std::vector<Eigen::Matrix3d> block(voxelsPerSector);
      for (int l = 0; l < voxelsPerSector; ++l) {
        block[l] = (penaltyDiagonal[l] + ridge) * Eigen::Matrix3d::Identity();
      }
      std::vector<Eigen::Matrix3d> weighted(voxels.size());
      for (size_t i = 0; i < voxels.size(); ++i) {
        const auto& voxel = voxels[i];
        weighted[i] = voxel.weight * voxel.a;
        const Eigen::Vector3d r = voxel.weight * voxel.b;
        for (int j = 0; j < voxel.nNodes; ++j) {
          rhs.segment<3>(3 * voxel.node[j]) += voxel.nodeWeight[j] * r;
          block[voxel.node[j]] += voxel.nodeWeight[j] * voxel.nodeWeight[j] * weighted[i];
        }
      }
      for (auto& matrix : block) {
        matrix = matrix.inverse().eval();
      }
      auto apply = [&](const Eigen::VectorXd& x, Eigen::VectorXd& y) {
        y = ridge * x;
        for (const auto& [axis, first] : stencils) {
          const int step = strides[axis];
          for (int c = 0; c < 3; ++c) {
            double difference = 0.;
            for (int k = 0; k < 3; ++k) {
              difference += Coefficients[k] * x(3 * (first + k * step) + c);
            }
            for (int k = 0; k < 3; ++k) {
              y(3 * (first + k * step) + c) += lambda * Coefficients[k] * difference;
            }
          }
        }
        for (size_t i = 0; i < voxels.size(); ++i) {
          const auto& voxel = voxels[i];
          const Eigen::Vector3d u = weighted[i] * nodeValue(voxel, x.data());
          for (int j = 0; j < voxel.nNodes; ++j) {
            y.segment<3>(3 * voxel.node[j]) += voxel.nodeWeight[j] * u;
          }
        }
      };
      auto precondition = [&](const Eigen::VectorXd& r, Eigen::VectorXd& z) {
        for (int l = 0; l < voxelsPerSector; ++l) {
          z.segment<3>(3 * l) = block[l] * r.segment<3>(3 * l);
        }
      };
      // preconditioned conjugate gradients from the previous cycle's solution
      Eigen::VectorXd x = Eigen::Map<Eigen::VectorXd>(sectorNodes, n);
      Eigen::VectorXd r(n);
      Eigen::VectorXd z(n);
      Eigen::VectorXd q(n);
      const double rhsNorm = rhs.norm();
      iterations[sector] = 0;
      residuals[sector] = 0.;
      if (rhsNorm == 0.) {
        std::fill(sectorNodes, sectorNodes + n, 0.);
        continue;
      }
      apply(x, q);
      r = rhs - q;
      precondition(r, z);
      Eigen::VectorXd p = z;
      double rz = r.dot(z);
      residuals[sector] = r.norm() / rhsNorm;
      while (residuals[sector] > SmoothCGTolerance && iterations[sector] < SmoothCGMaxIterations) {
        apply(p, q);
        const double alpha = rz / p.dot(q);
        x += alpha * p;
        r -= alpha * q;
        residuals[sector] = r.norm() / rhsNorm;
        ++iterations[sector];
        precondition(r, z);
        const double rzNew = r.dot(z);
        p = z + (rzNew / rz) * p;
        rz = rzNew;
      }
      finite[sector] = finite[sector] && x.allFinite();
      Eigen::Map<Eigen::VectorXd>(sectorNodes, n) = x;
    }
  };
  std::atomic<int> nextSector{0};
  std::vector<std::thread> threads;
  for (int thread = 0; thread < std::clamp(nThreads, 1, NSectors); ++thread) {
    threads.emplace_back([&]() {
      for (int sector = nextSector++; sector < NSectors; sector = nextSector++) {
        smoothSector(sector);
      }
    });
  }
  for (auto& thread : threads) {
    thread.join();
  }

  const long nWithEntries = std::accumulate(nVoxels.begin(), nVoxels.end(), 0L);
  const int nNotConverged = std::count_if(residuals.begin(), residuals.end(), [](double residual) { return residual > SmoothCGTolerance; });
  const int nNotFinite = std::count(finite.begin(), finite.end(), false);
  LOGP(info, "Smoothing: {} voxels with entries ({} solved, median chi2 / ndf scale {:.3f}), {} stencils per sector, tau {} um, tau0 {} um; CG iterations per sector {} - {}, worst relative residual {:.2e}",
       nWithEntries, scales.size(), medianScale, stencils.size(), tau * 1e4, tau0 * 1e4, *std::min_element(iterations.begin(), iterations.end()),
       *std::max_element(iterations.begin(), iterations.end()), *std::max_element(residuals.begin(), residuals.end()));
  if (robust > 0) {
    LOGP(info, "Smoothing: {} Tukey reweightings (c = {}), at the last {} voxels below weight 0.5, {} at 0", robust, SmoothTukey,
         std::accumulate(nLowWeight.begin(), nLowWeight.end(), 0L), std::accumulate(nZeroWeight.begin(), nZeroWeight.end(), 0L));
  }
  if (nNotConverged) {
    LOGP(warning, "Smoothing: the conjugate gradients did not reach the tolerance {} in {} sectors", SmoothCGTolerance, nNotConverged);
  }
  if (nNotFinite) {
    LOGP(error, "Smoothing: the solution is not finite in {} sectors", nNotFinite);
  }
  return nNotFinite == 0;
}

/// the cross-validation of pass 1 on the residuals of pass 2 (taken with pass 1's smoothed map S1 subtracted), per channel (y, z): the chi2
/// change of pass 2's residuals with pass 1's raw voxels (where solved, else none) and with S1, relative to no fine correction; the
/// corrections at the voxels' mean positions
struct CrossValidation {
  double entries = 0.;
  std::array<double, 2> raw{};
  std::array<double, 2> smoothed{};
};

CrossValidation crossValidate(const TrackResiduals& binning, const std::vector<std::array<ChargeSums, 2>>& sums, const std::vector<double>& passOne,
                              const std::vector<float>& rawPassOne)
{
  using S = ChargeSums;
  // chi2 of one channel's residuals at (dX, d) minus at (0, 0): sum w (r - d + s dX)^2 - sum w r^2
  auto channelChi2 = [](const std::array<double, S::NSums>& channel, double dX, double d) {
    return d * d * channel[S::W] + dX * dX * channel[S::WSS] - 2. * dX * d * channel[S::WS] - 2. * d * channel[S::WR] + 2. * dX * channel[S::WSR];
  };
  const int voxelsPerSector = binning.getNVoxelsPerSector();
  CrossValidation result;
  for (size_t voxel = 0; voxel < sums.size(); ++voxel) {
    ChargeSums total = sums[voxel][0];
    total.add(sums[voxel][1]);
    if (total.entries <= 0. || total.y[S::W] <= 0. || total.z[S::W] <= 0.) {
      continue;
    }
    const int sector = int(voxel / voxelsPerSector);
    const int local = int(voxel % voxelsPerSector);
    int nodes[4];
    double weights[4];
    const int nNodes = bilinearNodes(binning, local % binning.getNXBins(), total.y2x / total.entries, (sector < 18 ? 1. : -1.) * total.z2x / total.entries,
                                     nodes, weights);
    std::array<double, 3> smoothed{};
    for (int j = 0; j < nNodes; ++j) {
      for (int c = 0; c < 3; ++c) {
        smoothed[c] += weights[j] * passOne[3 * (size_t(sector) * voxelsPerSector + nodes[j]) + c];
      }
    }
    // relative to S1: no fine correction = -S1, the raw voxel = raw - S1, S1 = 0
    std::array<double, 3> raw{};
    for (int c = 0; c < 3; ++c) {
      raw[c] = rawPassOne[3 * voxel + c] - smoothed[c];
    }
    const double noneY = channelChi2(total.y, -smoothed[0], -smoothed[1]);
    const double noneZ = channelChi2(total.z, -smoothed[0], -smoothed[2]);
    result.entries += total.entries;
    result.raw[0] += channelChi2(total.y, raw[0], raw[1]) - noneY;
    result.raw[1] += channelChi2(total.z, raw[0], raw[2]) - noneZ;
    result.smoothed[0] -= noneY;
    result.smoothed[1] -= noneZ;
  }
  return result;
}

/// GBL's parameters (q/p, snp, tgl, y, z) of a state: a scattering kink changes only snp and tgl
std::array<double, 5> toGbl(const TrackParD& track)
{
  const double tgl = track.getTgl();
  return {track.getQ2Pt() / std::sqrt(1. + tgl * tgl), track.getSnp(), tgl, track.getY(), track.getZ()};
}

TrackParD fromGbl(const std::array<double, 5>& parameters, double alpha, double x)
{
  return TrackParD(x, alpha, {parameters[3], parameters[4], parameters[1], parameters[2], parameters[0] * std::sqrt(1. + parameters[2] * parameters[2])}, 1, o2::track::PID::Pion);
}

bool rotateTo(TrackParD& track, double alpha)
{
  return std::abs(std::remainder(alpha - track.getAlpha(), 2. * M_PI)) < 1e-12 || track.rotateParam(alpha);
}

/// d(GBL parameters at the plane x of frame alpha) / d(those of from): the helix in the uniform bz, central differences (TrackParD
/// computes in double precision)
bool stepJacobian(const TrackParD& from, double alpha, double x, double bz, Eigen::Matrix<double, 5, 5>& jacobian)
{
  const auto base = toGbl(from);
  const double steps[5] = {1e-6 * std::max(std::abs(base[0]), 1e-2), 1e-6, 1e-6, 1e-5, 1e-5};
  for (int column = 0; column < 5; ++column) {
    std::array<double, 5> moved[2];
    for (int sign = 0; sign < 2; ++sign) {
      auto parameters = base;
      parameters[column] += (sign ? -1. : 1.) * steps[column];
      TrackParD track = fromGbl(parameters, from.getAlpha(), from.getX());
      if (!rotateTo(track, alpha) || !track.propagateParamTo(x, bz)) {
        return false;
      }
      moved[sign] = toGbl(track);
    }
    for (int row = 0; row < 5; ++row) {
      jacobian(row, column) = (moved[0][row] - moved[1][row]) / (2. * steps[column]);
    }
  }
  return true;
}

double localBz(const TrackParD& track)
{
  double field[3];
  const auto position = track.getXYZGlo();
  o2::base::PropagatorD::Instance()->getFieldXYZ(o2::math_utils::Point3D<double>(position.X(), position.Y(), position.Z()), field);
  return field[2];
}

/// a point of the track fit: a TPC cluster, or with itsInFit an ITS cluster, the primary vertex or a scatterer without measurement
struct Cluster {
  enum Type { TPC,
              ITS,
              Vertex,
              Material };
  int type = TPC;
  int sector = 0;
  int row = 0;       ///< TPC pad row, ITS layer, -1 vertex
  double alpha = 0.; ///< the frame
  double x = 0.;     ///< TPC: corrected
  double y = 0.;
  double z = 0.;
  double padY = 0.; ///< measured y on the pad plane, before any correction
  double tgSlp = 0.;
  bool tgSlpValid = true;
  bool excluded = false; ///< de-weighted in the fit
  double sigmaY = 0.;    ///< the cluster's own error
  double sigmaZ = 0.;
};

/// the fit of a track: active = the clusters in the final pass, reference the states at their planes, outlier the last pass's flags
struct FitResult {
  std::vector<int> active;
  std::vector<TrackParD> reference;
  std::vector<char> outlier;
  std::unique_ptr<gbl::GblTrajectory> trajectory;
  double qOverPt = 0.;
  int passes = 0;
};

bool fitTrack(std::vector<Cluster>& clusters, TrackParD seed, FitResult& result)
{
  auto* propagator = o2::base::PropagatorD::Instance();
  result.active.resize(clusters.size());
  std::iota(result.active.begin(), result.active.end(), 0);
  std::vector<char> flagged(clusters.size(), 0);
  // with points before the TPC (itsInFit) the reference follows the kinks fitted there: per point the kink (snp, tgl) of the previous pass,
  // the reference continues with it (Gauss-Newton on the kinked track; a smooth reference misses ~kink^2 x lever arm, mrad kinks over the
  // gap: tens of um, proportional to the track slope)
  const bool kinkedReference = clusters.front().type != Cluster::TPC;
  std::vector<std::array<double, 2>> kinks(clusters.size(), {0., 0.});
  const double momentum = std::sqrt(1. + seed.getTgl() * seed.getTgl()) / std::abs(seed.getQ2Pt());
  const double mass = o2::track::PID::getMass(o2::track::PID::Pion);
  const double beta = momentum / std::sqrt(momentum * momentum + mass * mass);
  Eigen::VectorXd residuals(2);
  Eigen::VectorXd measurementErrors(2);
  for (int pass = 0; pass < MaxFitPasses; ++pass) {
    // reference from the seed through the clusters still active (not reached: dropped; flagged by the previous pass: dropped)
    std::vector<int> kept;
    std::vector<TrackParD> reference; ///< at each point, in front of its kink
    std::vector<TrackParD> leaving;   ///< at each point, after its kink
    TrackParD state = seed;
    for (int index : result.active) {
      if (pass >= FirstOutlierPass && flagged[index]) {
        continue;
      }
      TrackParD next = state;
      if (!rotateTo(next, clusters[index].alpha) || !propagator->PropagateToXBxByBz(next, clusters[index].x, 0.999, 2., MatCorrType::USEMatCorrLUT) ||
          std::abs(next.getSnp()) > MaxSnp) {
        continue;
      }
      kept.push_back(index);
      reference.push_back(next);
      if (kinkedReference && kept.size() > 1 && clusters[index].type != Cluster::TPC && (kinks[index][0] != 0. || kinks[index][1] != 0.)) {
        auto parameters = toGbl(next);
        parameters[1] += kinks[index][0];
        parameters[2] += kinks[index][1];
        next = fromGbl(parameters, next.getAlpha(), next.getX());
      }
      leaving.push_back(next);
      state = next;
    }
    result.active.swap(kept);
    const size_t nPoints = result.active.size();
    if (std::count_if(result.active.begin(), result.active.end(), [&](int index) { return clusters[index].type != Cluster::Material; }) < MinMeasurements) {
      return false;
    }
    // the points before the first TPC cluster (vertex, ITS, gap)
    const size_t nInner = std::find_if(result.active.begin(), result.active.end(), [&](int index) { return clusters[index].type == Cluster::TPC; }) -
                          result.active.begin();
    // the TPC cluster errors from the first reference (the stored tgSlp where valid)
    if (pass == 0) {
      for (size_t point = nInner; point < nPoints; ++point) {
        auto& cluster = clusters[result.active[point]];
        const double snp = reference[point].getSnp();
        clusterErrors(cluster.row, cluster.z, cluster.tgSlpValid ? cluster.tgSlp : snp / std::sqrt((1. - snp) * (1. + snp)), reference[point].getTgl(),
                      cluster.sigmaY, cluster.sigmaZ);
      }
    }
    // x/X0 of the step before each point, the Highland log term of the whole track
    std::vector<double> stepXOverX0(nPoints, 0.);
    double totalXOverX0 = 0.;
    for (size_t point = 1; point < nPoints; ++point) {
      stepXOverX0[point] = propagator->getMatBudget(MatCorrType::USEMatCorrLUT, reference[point - 1].getXYZGlo(), reference[point].getXYZGlo()).meanX2X0;
      totalXOverX0 += stepXOverX0[point];
    }
    const double theta0 = 0.0136 / (beta * momentum) * (1. + 0.038 * std::log(std::max(totalXOverX0, 1e-12)));
    // before the TPC a scatterer at every point with half of the steps on either side (the step after the first point to the second, the
    // step into the TPC to the last point before it); in the TPC a scatterer every ScatterEvery points, carrying the x/X0 of its
    // ScatterEvery steps (centred)
    std::vector<double> scattererXOverX0(nPoints, 0.);
    for (size_t point = 1; point <= nInner && point < nPoints; ++point) {
      if (point == 1) {
        scattererXOverX0[1] += stepXOverX0[1];
      } else if (point == nInner) {
        scattererXOverX0[point - 1] += stepXOverX0[point];
      } else {
        scattererXOverX0[point - 1] += 0.5 * stepXOverX0[point];
        scattererXOverX0[point] += 0.5 * stepXOverX0[point];
      }
    }
    constexpr size_t Every = ScatterEvery;
    const size_t nTPC = nPoints - nInner;
    for (size_t point = nInner + 1; point < nPoints; ++point) {
      const size_t step = point - nInner;
      const size_t scatterer = nInner + (Every == 1 ? step : std::clamp<size_t>(((step - 1) / Every) * Every + (Every + 1) / 2, 1, nTPC > 2 ? nTPC - 2 : 1));
      scattererXOverX0[scatterer] += stepXOverX0[point];
    }
    std::vector<gbl::GblPoint> points;
    points.reserve(nPoints);
    for (size_t point = 0; point < nPoints; ++point) {
      const auto& cluster = clusters[result.active[point]];
      Eigen::Matrix<double, 5, 5> jacobian = Eigen::Matrix<double, 5, 5>::Identity();
      if (point > 0 && !stepJacobian(leaving[point - 1], reference[point].getAlpha(), reference[point].getX(),
                                     0.5 * (localBz(leaving[point - 1]) + localBz(reference[point])), jacobian)) {
        return false;
      }
      gbl::GblPoint gblPoint(jacobian);
      if (cluster.type != Cluster::Material) {
        const double scale = cluster.excluded ? ExcludedSigmaScale : 1.;
        gblPoint.addMeasurement(Eigen::Vector2d(cluster.y - reference[point].getY(), cluster.z - reference[point].getZ()),
                                Eigen::Vector2d(1. / std::pow(cluster.sigmaY * scale, 2), 1. / std::pow(cluster.sigmaZ * scale, 2)));
      }
      // a kinked point before the TPC keeps its scatterer (residual = minus the reference's kink) even without material in its steps
      const bool kinked = kinkedReference && point < nInner;
      const double theta2 = theta0 * theta0 * (kinked ? std::max(scattererXOverX0[point], 1e-9) : scattererXOverX0[point]);
      if (point > 0 && point + 1 < nPoints && theta2 > 0.) {
        const double snp = reference[point].getSnp();
        const double secLambda2 = 1. + reference[point].getTgl() * reference[point].getTgl();
        const auto& kink = kinks[result.active[point]];
        gblPoint.addScatterer(kinked ? Eigen::Vector2d(-kink[0], -kink[1]) : Eigen::Vector2d::Zero(),
                              Eigen::Vector2d(1. / (theta2 * (1. - snp * snp) * secLambda2), 1. / (theta2 * secLambda2 * secLambda2)));
      }
      points.push_back(gblPoint);
    }
    result.trajectory = std::make_unique<gbl::GblTrajectory>(points, true);
    double chi2 = 0.;
    int ndf = 0;
    double lostWeight = 0.;
    if (!result.trajectory->isValid() || result.trajectory->fit(chi2, ndf, lostWeight) != 0) {
      return false;
    }
    // outliers: residuals to this pass's fit beyond OutlierSigma of the error in the fit (dropped at the start of the next pass)
    result.outlier.assign(nPoints, 0);
    bool newOutliers = false;
    if (pass + 1 >= FirstOutlierPass) {
      for (size_t point = 0; point < nPoints; ++point) {
        unsigned int numData = 0;
        result.trajectory->getMeasResults(point + 1, numData, residuals, measurementErrors);
        result.outlier[point] = numData >= 2 && (std::abs(residuals(0)) > OutlierSigma * measurementErrors(0) ||
                                                 std::abs(residuals(1)) > OutlierSigma * measurementErrors(1));
        flagged[result.active[point]] = result.outlier[point];
        newOutliers |= result.outlier[point] != 0;
      }
    }
    // the seed: the first point's reference plus GBL's correction there
    Eigen::VectorXd corrections;
    Eigen::MatrixXd covariance;
    if (result.trajectory->getResults(1, corrections, covariance) != 0) {
      return false;
    }
    auto parameters = toGbl(reference.front());
    for (int parameter = 0; parameter < 5; ++parameter) {
      parameters[parameter] += corrections(parameter);
    }
    seed = fromGbl(parameters, reference.front().getAlpha(), reference.front().getX());
    // the kinks before the TPC for the next reference: + the fitted change of the direction across each point
    double kinkUpdate = 0.;
    if (kinkedReference) {
      Eigen::VectorXd before;
      Eigen::VectorXd after;
      for (size_t point = 1; point < nInner && point + 1 < nPoints; ++point) {
        if (result.trajectory->getResults(-int(point + 1), before, covariance) != 0 || result.trajectory->getResults(int(point + 1), after, covariance) != 0) {
          return false;
        }
        auto& kink = kinks[result.active[point]];
        for (int direction = 0; direction < 2; ++direction) {
          kink[direction] += after(direction + 1) - before(direction + 1);
          kinkUpdate = std::max(kinkUpdate, std::abs(after(direction + 1) - before(direction + 1)));
        }
      }
    }
    result.reference.swap(reference);
    result.passes = pass + 1;
    // converged: from the first outlier pass on, no new outlier and GBL's correction at the first point negligible (1 um in y / z, 1e-6 in
    // snp / tgl, 1e-4 of q/p): the reference is the fit, further passes change nothing
    if (pass >= FirstOutlierPass && !newOutliers && std::abs(corrections(3)) < 1e-4 && std::abs(corrections(4)) < 1e-4 &&
        std::abs(corrections(1)) < 1e-6 && std::abs(corrections(2)) < 1e-6 && std::abs(corrections(0)) < 1e-4 * std::abs(parameters[0]) && kinkUpdate < 1e-6) {
      break;
    }
  }
  result.qOverPt = seed.getQ2Pt();
  return true;
}

enum Counter { Selected,
               Fitted,
               FitFailed,
               FitPasses,
               AtMaxPasses,
               OutOfRowOrder,
               Residuals,
               ResidualsExcluded,
               Outliers,
               EdgePadResiduals,
               OutsideGrid,
               NoResidual,
               NoInnerPoints,
               NCounters };

/// what the fit threads share read-only
struct FitContext {
  GID::mask_t sources;
  const BaseCorrection* baseCorrection = nullptr; ///< nullptr: no base correction
  const FineMapInput* fineMap = nullptr;          ///< nullptr: no previous fine map
  const FineMapInput* passOne = nullptr;          ///< pass 2: the smoothed map of pass 1, subtracted after fineMap
  const TrackResiduals* binning = nullptr;
  const RowGeometry* geometry = nullptr;
  float maxZ2X = 0.f;
  int voxelsPerSector = 0;
  int itsInFit = 0;       ///< 1: the ITS clusters in the fit, 2: also the primary vertex
  bool fitSlopes = false; ///< the slopes of the voxel model from the fitted track instead of the reference
};

/// itsInFit: the track's ITS clusters (2: and its primary vertex) in front of its TPC clusters, from the inside out, and scatterer points in
/// the gap between the last of them and the first TPC cluster (in that cluster's frame, centres of equal steps of at most GapStep); false
/// with fewer than MinITSClusters ITS clusters or without the vertex
bool addInnerPoints(const TFTask& task, const o2::tpc::TrackDataCompact& trackRef, int itsInFit, std::vector<Cluster>& clusters)
{
  auto decode = [](short value, float maximum) { return value * maximum / 0x7fff; };
  std::vector<Cluster> inner;
  const unsigned int firstExternal = trackRef.idxFirstResidual + trackRef.nResiduals;
  for (unsigned int residualIndex = firstExternal; residualIndex < firstExternal + trackRef.nExtDetResid; ++residualIndex) {
    const auto& residual = task.residuals[residualIndex];
    Cluster point;
    // reference + residual, in the frame of the sensor (ITS) or of the track at its closest approach (vertex)
    point.y = decode(residual.y, o2::tpc::param::MaxY) + decode(residual.dy, o2::tpc::param::MaxResid);
    point.z = decode(residual.z, o2::tpc::param::MaxZ) + decode(residual.dz, o2::tpc::param::MaxResid);
    if (residual.row == RowVertex) {
      if (itsInFit < 2) {
        continue;
      }
      point.type = Cluster::Vertex;
      point.row = -1;
      point.alpha = decode(residual.tgSlp, o2::tpc::param::MaxTgSlp) * M_PI;
      point.x = decode(residual.channel, o2::tpc::param::MaxVtxX);
      point.sigmaY = VertexSigmaY;
      point.sigmaZ = VertexSigmaZ;
    } else if (residual.isITS()) {
      point.type = Cluster::ITS;
      point.row = residual.getITSLayer();
      point.alpha = residual.getAlpha();
      point.x = residual.getX();
      point.sigmaY = ITSSigmaY;
      point.sigmaZ = ITSSigmaZ;
    } else {
      continue;
    }
    inner.push_back(point);
  }
  std::sort(inner.begin(), inner.end(), [](const Cluster& first, const Cluster& second) { return first.row < second.row || (first.row == second.row && first.x < second.x); });
  const long nITS = std::count_if(inner.begin(), inner.end(), [](const Cluster& point) { return point.type == Cluster::ITS; });
  if (nITS < MinITSClusters || (itsInFit >= 2 && inner.front().type != Cluster::Vertex)) {
    return false;
  }
  const Cluster& last = inner.back();
  const Cluster& first = clusters.front();
  const double start = last.x * std::cos(last.alpha - first.alpha) - last.y * std::sin(last.alpha - first.alpha);
  const int nSteps = std::max(1, int(std::ceil((first.x - start) / GapStep)));
  for (int step = 0; step < nSteps; ++step) {
    Cluster point;
    point.type = Cluster::Material;
    point.sector = first.sector;
    point.row = -1;
    point.alpha = first.alpha;
    point.x = start + (step + 0.5) * (first.x - start) / nSteps;
    inner.push_back(point);
  }
  clusters.insert(clusters.begin(), inner.begin(), inner.end());
  return true;
}

/// fit one track of a TF and append the contributions of its clusters
void processTrack(const TFTask& task, size_t trackIndex, const FitContext& context, std::array<long, NCounters>& counters,
                  std::vector<Cluster>& clusters, std::vector<Contribution>& batch)
{
  const auto& trackRef = task.trackRefs[trackIndex];
  const auto& itsTrack = task.tracks[trackIndex].par;
  if (!GID::includesSource(trackRef.sourceId, context.sources) || trackRef.filterFlag > 0 || std::abs(itsTrack.getQ2Pt()) < 1e-9 ||
      trackRef.idxFirstResidual + trackRef.nResiduals + (context.itsInFit ? trackRef.nExtDetResid : 0) > task.residuals.size()) {
    return;
  }
  auto decodeResidual = [](short value) { return value * o2::tpc::param::MaxResid / 0x7fff; };
  clusters.clear();
  int nFitClusters = 0;
  for (unsigned int residualIndex = trackRef.idxFirstResidual; residualIndex < trackRef.idxFirstResidual + trackRef.nResiduals; ++residualIndex) {
    const auto& residual = task.residuals[residualIndex];
    if (residual.rejected) {
      continue;
    }
    Cluster cluster;
    cluster.sector = residual.sec;
    cluster.row = residual.row;
    cluster.alpha = o2::math_utils::sector2Angle(cluster.sector % 18);
    if (!clusters.empty() && cluster.row <= clusters.back().row) {
      ++counters[OutOfRowOrder];
      continue;
    }
    cluster.x = residual.getX();
    cluster.y = residual.y * o2::tpc::param::MaxY / 0x7fff + decodeResidual(residual.dy);
    cluster.z = residual.z * o2::tpc::param::MaxZ / 0x7fff + decodeResidual(residual.dz);
    cluster.padY = cluster.y;
    // position-only cluster or saturated tgSlp: the reference's angle instead
    cluster.tgSlpValid = !residual.isPositionOnly() && !residual.isTgSlpClamped();
    cluster.tgSlp = residual.tgSlp * o2::tpc::param::MaxTgSlp / 0x7fff;
    if (context.baseCorrection) {
      float x = cluster.x;
      float y = cluster.y;
      float z = cluster.z;
      context.baseCorrection->apply(cluster.sector, cluster.row, task.info->lumiScale, x, y, z);
      cluster.x = x;
      cluster.y = y;
      cluster.z = z;
    }
    if (context.fineMap) {
      context.fineMap->apply(cluster.sector, cluster.row, cluster.x, cluster.y, cluster.z);
    }
    if (context.passOne) {
      context.passOne->apply(cluster.sector, cluster.row, cluster.x, cluster.y, cluster.z);
    }
    cluster.excluded = isExcludedFromFit(cluster.row, cluster.padY, *context.geometry);
    nFitClusters += !cluster.excluded;
    clusters.push_back(cluster);
  }
  if (nFitClusters < MinFitClusters || clusters.empty()) {
    return;
  }
  ++counters[Selected];
  if (context.itsInFit && !addInnerPoints(task, trackRef, context.itsInFit, clusters)) {
    ++counters[NoInnerPoints];
    return;
  }
  // the seed: the ITS track's outer parameters at the first point
  TrackParD seed(itsTrack.getX(), itsTrack.getAlpha(), {itsTrack.getY(), itsTrack.getZ(), itsTrack.getSnp(), itsTrack.getTgl(), itsTrack.getQ2Pt()}, 1, o2::track::PID::Pion);
  FitResult fit;
  if (!rotateTo(seed, clusters.front().alpha) ||
      !o2::base::PropagatorD::Instance()->PropagateToXBxByBz(seed, clusters.front().x, 0.999, 2., MatCorrType::USEMatCorrLUT) ||
      !fitTrack(clusters, seed, fit)) {
    ++counters[FitFailed];
    return;
  }
  ++counters[Fitted];
  counters[FitPasses] += fit.passes;
  counters[AtMaxPasses] += fit.passes >= MaxFitPasses;
  const int charge = fit.qOverPt > 0. ? 1 : 0;
  Eigen::VectorXd residuals(2);
  Eigen::VectorXd measurementErrors(2);
  Eigen::VectorXd residualErrors(2);
  Eigen::VectorXd downWeights(2);
  Eigen::VectorXd corrections;
  Eigen::MatrixXd covariance;
  for (size_t point = 0; point < fit.active.size(); ++point) {
    const auto& cluster = clusters[fit.active[point]];
    if (fit.outlier[point] || cluster.type != Cluster::TPC) {
      continue;
    }
    unsigned int numData = 0;
    fit.trajectory->getMeasResults(point + 1, numData, residuals, measurementErrors, residualErrors, downWeights);
    double unbiased[2] = {0., 0.};
    double weight[2] = {0., 0.};
    double bestPosition[2] = {cluster.y, cluster.z};
    bool valid = numData >= 2;
    for (int coordinate = 0; coordinate < 2 && valid; ++coordinate) {
      const double fitError2 = measurementErrors(coordinate) * measurementErrors(coordinate);
      const double residualError2 = residualErrors(coordinate) * residualErrors(coordinate);
      if (!(residualError2 > 1e-4 * fitError2)) {
        valid = false;
        break;
      }
      const double sigma = coordinate == 0 ? cluster.sigmaY : cluster.sigmaZ;
      const double withoutPoint = (fitError2 - residualError2) * fitError2 / residualError2;
      unbiased[coordinate] = residuals(coordinate) * fitError2 / residualError2;
      weight[coordinate] = 1. / (sigma * sigma + withoutPoint);
      bestPosition[coordinate] -= unbiased[coordinate] * sigma * sigma * weight[coordinate];
    }
    if (!valid) {
      ++counters[NoResidual];
      continue;
    }
    if (unbiased[0] * unbiased[0] * weight[0] > OutlierSigma * OutlierSigma || unbiased[1] * unbiased[1] * weight[1] > OutlierSigma * OutlierSigma) {
      ++counters[Outliers];
      continue;
    }
    const double bestPadY = cluster.padY + bestPosition[0] - cluster.y;
    if (std::abs(bestPadY) > context.geometry->padPlaneEdge[cluster.row] - SkipEdgePads * context.geometry->padWidth[cluster.row]) {
      ++counters[EdgePadResiduals];
      continue;
    }
    const double y2x = bestPosition[0] / cluster.x;
    const double z2x = bestPosition[1] / cluster.x;
    const int binY2X = std::abs(z2x) > context.maxZ2X ? -1 : context.binning->getY2XBinExact(y2x, cluster.row);
    const int binZ2X = context.binning->getZ2XBinExact(cluster.sector < 18 ? z2x : -z2x);
    if (binY2X < 0 || binY2X >= context.binning->getNY2XBins() || binZ2X < 0 || std::abs(z2x) > context.maxZ2X) {
      ++counters[OutsideGrid];
      continue;
    }
    // the slopes: the reference's, or the fitted track's arriving at the cluster: with the kinks fitted (the ITS in the fit), only these
    // are uncorrelated with the cluster's unbiased residual
    double snp = fit.reference[point].getSnp();
    double tgl = fit.reference[point].getTgl();
    if (context.fitSlopes && fit.trajectory->getResults(-int(point + 1), corrections, covariance) == 0) {
      snp += corrections(1);
      tgl += corrections(2);
    }
    const double cosPhi = std::sqrt(std::max(1e-6, (1. - snp) * (1. + snp)));
    Contribution contribution;
    contribution.voxel = size_t(cluster.sector) * context.voxelsPerSector + context.binning->getGlbVoxBin(cluster.row, binY2X, binZ2X);
    contribution.charge = charge;
    contribution.x = cluster.x;
    contribution.y2x = y2x;
    contribution.z2x = z2x;
    contribution.slopeY = snp / cosPhi;
    contribution.residualY = unbiased[0];
    contribution.weightY = weight[0];
    contribution.slopeZ = tgl / cosPhi;
    contribution.residualZ = unbiased[1];
    contribution.weightZ = weight[1];
    batch.push_back(contribution);
    ++counters[Residuals];
    counters[ResidualsExcluded] += cluster.excluded;
  }
}

/// per selected TF: the time, CTPLumi orbit and luminosity, the IDC and the lumi scale applied
struct SelectedTFs {
  std::vector<long> timeMS;
  std::vector<uint32_t> orbits;
  std::vector<float> ctpLumi;
  std::vector<float> idcA;
  std::vector<float> idcC;
  std::vector<float> lumiScale;
  std::vector<float> idc; ///< (A + C) / 2 of the TFs with IDC
};

/// voxResTree (with the binning, the CTP / IDC means and medians in the user info), OrbitLumiInfo and optionally fineSums. smoothed: the
/// nodes of the map (empty: DS = D); two passes: sums = pass 2's, passOne = S1 (added to pass 2's voxel solve at the mean position: D),
/// sumsPassOne = pass 1's (if writeSums: fineSums, pass 2's in fineSumsPass2); extraInfo: more user info
void writeOutput(const std::string& fileOutput, const TrackResiduals& binning, const std::string& y2xBinning, const std::string& z2xBinning, float maxZ2X,
                 const std::vector<std::array<ChargeSums, 2>>& sums, const std::vector<double>& smoothed, const std::vector<double>& passOne,
                 const std::vector<std::array<ChargeSums, 2>>& sumsPassOne, const std::vector<std::pair<std::string, std::string>>& extraInfo, long nTracks,
                 double minEntries, bool writeSums, int runNumber, double orbitResetMS, bool isPbPb, long firstTFTime, long lastTFTime, const SelectedTFs& selected)
{
  // CTP luminosity and IDC of the selected TFs, as staticMapCreatorCPM.C (PbPb: mean CTP luminosity x 2.414)
  double meanCTP = calculateMean(selected.ctpLumi) * (isPbPb ? 2.414 : 1.);
  float medianCTP = static_cast<float>(calculateMedian(selected.ctpLumi));
  float meanIDC = static_cast<float>(calculateMean(selected.idc));
  float medianIDC = static_cast<float>(calculateMedian(selected.idc));
  long meanTimeMS = static_cast<long>(calculateMean(selected.timeMS));
  long medianTimeMS = static_cast<long>(calculateMedian(selected.timeMS));

  // the trees are destroyed before the file (end of this scope)
  TFile output(fileOutput.c_str(), "recreate");
  TTree voxTree("voxResTree", "voxRes map results and statistics");
  voxTree.SetAlias("z2xBin", "bvox[0]");
  voxTree.SetAlias("y2xBin", "bvox[1]");
  voxTree.SetAlias("xBin", "bvox[2]");
  voxTree.SetAlias("z2xAV", "stat[0]");
  voxTree.SetAlias("y2xAV", "stat[1]");
  voxTree.SetAlias("xAV", "stat[2]");
  voxTree.SetAlias("fsector", "bsec+0.5+9.*(y2xAV)/pi");
  voxTree.SetAlias("phi", "(bsec%18+0.5+9.*(stat[1])/pi)/9*pi");
  voxTree.SetAlias("r", "stat[2]");
  voxTree.SetAlias("z", "z2xAV*xAV");
  voxTree.SetAlias("dX", "D[0]");
  voxTree.SetAlias("dY", "D[1]");
  voxTree.SetAlias("dZ", "D[2]");
  voxTree.SetAlias("dXS", "DS[0]");
  voxTree.SetAlias("dYS", "DS[1]");
  voxTree.SetAlias("dZS", "DS[2]");
  voxTree.SetAlias("dXE", "E[0]");
  voxTree.SetAlias("dYE", "E[1]");
  voxTree.SetAlias("dZE", "E[2]");
  voxTree.SetAlias("voxelIndex", Form("xBin + %i * (y2xBin + %i * z2xBin) + %i * bsec", binning.getNXBins(), binning.getNY2XBins(), binning.getNVoxelsPerSector()));
  voxTree.SetAlias("entries", "stat[3]");
  voxTree.SetAlias("fitOK", Form("(flags & %u) == %u", TrackResiduals::DistDone, TrackResiduals::DistDone));
  voxTree.SetAlias("dispOK", Form("(flags & %u) == %u", TrackResiduals::DispDone, TrackResiduals::DispDone));
  voxTree.SetAlias("smtOK", Form("(flags & %u) == %u", TrackResiduals::SmoothDone, TrackResiduals::SmoothDone));
  voxTree.SetAlias("masked", Form("(flags & %u) == %u", TrackResiduals::Masked, TrackResiduals::Masked));
  TrackResiduals::VoxRes result;
  TrackResiduals::VoxRes* resultPointer = &result;
  voxTree.Branch("voxRes", &resultPointer);
  auto* userInfo = voxTree.GetUserInfo();
  userInfo->Add(new TNamed("meanIDC", std::to_string(meanIDC).data()));
  userInfo->Add(new TNamed("meanCTP", std::to_string(meanCTP).data()));
  userInfo->Add(new TNamed("meanTimeMS", std::to_string(meanTimeMS).data()));
  userInfo->Add(new TNamed("medianIDC", std::to_string(medianIDC).data()));
  userInfo->Add(new TNamed("medianCTP", std::to_string(medianCTP).data()));
  userInfo->Add(new TNamed("medianTimeMS", std::to_string(medianTimeMS).data()));
  userInfo->Add(new TNamed("y2xBinning", y2xBinning.data()));
  userInfo->Add(new TNamed("z2xBinning", z2xBinning.data()));
  userInfo->Add(new TNamed("maxZ2X", std::to_string(maxZ2X).data()));
  userInfo->Add(new TNamed("nTracks", std::to_string(nTracks).data()));
  userInfo->Add(new TNamed("nTracksProcessed", std::to_string(nTracks).data()));
  userInfo->Add(new TNamed("runNumber", std::to_string(runNumber).data()));
  userInfo->Add(new TNamed("firstTFTime", std::to_string(firstTFTime).data()));
  userInfo->Add(new TNamed("lastTFTime", std::to_string(lastTFTime).data()));
  for (const auto& [key, value] : extraInfo) {
    userInfo->Add(new TNamed(key.data(), value.data()));
  }
  if (!selected.timeMS.empty()) {
    const auto [minTimeMS, maxTimeMS] = std::minmax_element(selected.timeMS.begin(), selected.timeMS.end());
    userInfo->Add(new TNamed("minTimeMS", std::to_string(*minTimeMS).data()));
    userInfo->Add(new TNamed("maxTimeMS", std::to_string(*maxTimeMS).data()));
  }
  const int voxelsPerSector = binning.getNVoxelsPerSector();
  long nSolved = 0;
  for (int sector = 0; sector < NSectors; ++sector) {
    for (int binZ2X = 0; binZ2X < binning.getNZ2XBins(); ++binZ2X) {
      for (int binY2X = 0; binY2X < binning.getNY2XBins(); ++binY2X) {
        for (int binX = 0; binX < binning.getNXBins(); ++binX) {
          const size_t voxel = size_t(sector) * voxelsPerSector + binning.getGlbVoxBin(binX, binY2X, binZ2X);
          ChargeSums total = sums[voxel][0];
          total.add(sums[voxel][1]);
          result = TrackResiduals::VoxRes();
          float centreX = 0.f;
          float centreY2X = 0.f;
          float centreZ2X = 0.f;
          binning.getVoxelCoordinates(sector, binX, binY2X, binZ2X, centreX, centreY2X, centreZ2X);
          result.stat[TrackResiduals::VoxX] = total.entries > 0. ? total.x / total.entries : centreX;
          result.stat[TrackResiduals::VoxF] = total.entries > 0. ? total.y2x / total.entries : centreY2X;
          result.stat[TrackResiduals::VoxZ] = total.entries > 0. ? total.z2x / total.entries : centreZ2X;
          result.stat[TrackResiduals::VoxV] = total.entries;
          result.bvox[TrackResiduals::VoxZ] = binZ2X;
          result.bvox[TrackResiduals::VoxF] = binY2X;
          result.bvox[TrackResiduals::VoxX] = binX;
          result.bsec = sector;
          VoxelSolution solution;
          if (solveVoxel(total, minEntries, solution)) {
            ++nSolved;
            if (!passOne.empty()) {
              int nodes[4];
              double weights[4];
              const int nNodes = bilinearNodes(binning, binX, result.stat[TrackResiduals::VoxF], (sector < 18 ? 1. : -1.) * result.stat[TrackResiduals::VoxZ],
                                               nodes, weights);
              for (int j = 0; j < nNodes; ++j) {
                for (int coordinate = 0; coordinate < 3; ++coordinate) {
                  solution.value[coordinate] += weights[j] * passOne[3 * (size_t(sector) * voxelsPerSector + nodes[j]) + coordinate];
                }
              }
            }
            for (int coordinate = 0; coordinate < 3; ++coordinate) {
              result.D[coordinate] = solution.value[coordinate];
              result.DS[coordinate] = solution.value[coordinate];
              result.E[coordinate] = solution.error[coordinate];
            }
            result.EXYCorr = solution.correlationXY;
            result.dYSigMAD = solution.rmsY;
            result.dZSigLTM = solution.rmsZ;
            result.D[TrackResiduals::ResD] = solution.rmsY;
            result.DS[TrackResiduals::ResD] = solution.rmsY;
            result.flags = TrackResiduals::DistDone | TrackResiduals::DispDone;
          }
          if (!smoothed.empty()) {
            for (int coordinate = 0; coordinate < 3; ++coordinate) {
              result.DS[coordinate] = smoothed[3 * voxel + coordinate];
            }
            result.flags |= TrackResiduals::SmoothDone;
          }
          voxTree.Fill();
        }
      }
    }
  }
  voxTree.Write();
  LOGP(info, "{}: {} voxels, {} solved (min entries {}), {}", fileOutput, sums.size(), nSolved, minEntries, smoothed.empty() ? "not smoothed" : "smoothed");

  TTree orbitLumiTree("OrbitLumiInfo", "Orbit and Lumi Info");
  int64_t orbitResetTimeMS = orbitResetMS;
  std::vector<long> timeMSsel = selected.timeMS;
  std::vector<uint32_t> orbitsSel = selected.orbits;
  std::vector<float> idcScalerASel = selected.idcA;
  std::vector<float> idcScalerCSel = selected.idcC;
  std::vector<float> ctpLumiSel = selected.ctpLumi;
  std::vector<float> lumiScaleSel = selected.lumiScale;
  orbitLumiTree.Branch("orbitResetTimeMS", &orbitResetTimeMS);
  orbitLumiTree.Branch("timeMSsel", &timeMSsel);
  orbitLumiTree.Branch("orbitsSel", &orbitsSel);
  orbitLumiTree.Branch("idcScalerASel", &idcScalerASel);
  orbitLumiTree.Branch("idcScalerCSel", &idcScalerCSel);
  orbitLumiTree.Branch("ctpLumiSel", &ctpLumiSel);
  orbitLumiTree.Branch("lumiScaleSel", &lumiScaleSel);
  orbitLumiTree.Branch("meanIDC", &meanIDC);
  orbitLumiTree.Branch("meanCTP", &meanCTP);
  orbitLumiTree.Branch("meanTimeMS", &meanTimeMS);
  orbitLumiTree.Branch("medianIDC", &medianIDC);
  orbitLumiTree.Branch("medianCTP", &medianCTP);
  orbitLumiTree.Branch("medianTimeMS", &medianTimeMS);
  orbitLumiTree.Fill();
  orbitLumiTree.Write();

  auto writeSumsTree = [&](const char* name, const char* title, const std::vector<std::array<ChargeSums, 2>>& treeSums) {
    TTree sumsTree(name, title);
    for (const char* key : {"y2xBinning", "z2xBinning", "maxZ2X", "nTracks"}) {
      sumsTree.GetUserInfo()->Add(userInfo->FindObject(key)->Clone());
    }
    int voxelIndex = 0;
    float entries[2] = {};
    float position[2][3] = {};
    float channelY[2][ChargeSums::NSums] = {};
    float channelZ[2][ChargeSums::NSums] = {};
    sumsTree.Branch("voxel", &voxelIndex);
    sumsTree.Branch("entries", entries, "entries[2]/F");
    sumsTree.Branch("position", position, "position[2][3]/F");
    sumsTree.Branch("y", channelY, "y[2][6]/F");
    sumsTree.Branch("z", channelZ, "z[2][6]/F");
    for (size_t voxel = 0; voxel < treeSums.size(); ++voxel) {
      if (treeSums[voxel][0].entries + treeSums[voxel][1].entries <= 0.) {
        continue;
      }
      voxelIndex = int(voxel);
      for (int charge = 0; charge < 2; ++charge) {
        const auto& chargeSums = treeSums[voxel][charge];
        entries[charge] = chargeSums.entries;
        position[charge][0] = chargeSums.x;
        position[charge][1] = chargeSums.y2x;
        position[charge][2] = chargeSums.z2x;
        for (int sum = 0; sum < ChargeSums::NSums; ++sum) {
          channelY[charge][sum] = chargeSums.y[sum];
          channelZ[charge][sum] = chargeSums.z[sum];
        }
      }
      sumsTree.Fill();
    }
    sumsTree.Write();
  };
  if (writeSums && sumsPassOne.empty()) {
    writeSumsTree("fineSums", "per voxel and charge the sums of the unbiased TPC-only residuals", sums);
  } else if (writeSums) {
    writeSumsTree("fineSums", "pass 1: per voxel and charge the sums of the unbiased TPC-only residuals", sumsPassOne);
    writeSumsTree("fineSumsPass2", "pass 2: per voxel and charge the sums of the unbiased TPC-only residuals with pass 1's smoothed map subtracted", sums);
  }

  const std::filesystem::path outputPath(fileOutput);
  const std::string infoFile = (outputPath.parent_path().empty() ? std::filesystem::path(".") : outputPath.parent_path()) / (outputPath.stem().string() + ".txt");
  std::ofstream info(infoFile);
  info << "meanIDC: " << meanIDC << "\n";
  info << "meanCTP: " << meanCTP << "\n";
  info << "meanTimeMS: " << meanTimeMS << "\n";
  info << "medianIDC: " << medianIDC << "\n";
  info << "medianCTP: " << medianCTP << "\n";
  info << "medianTimeMS: " << medianTimeMS << "\n";
  LOGP(info, "meanIDC {}, medianIDC {}, meanCTP {}, medianCTP {}, meanTimeMS {}, medianTimeMS {}", meanIDC, medianIDC, meanCTP, medianCTP, meanTimeMS, medianTimeMS);
}

/// the per-TF readers of the open file (destroyed before the file)
struct InputFile {
  std::unique_ptr<TFile> file;
  std::unique_ptr<TTreeReader> residualReader;
  std::unique_ptr<TTreeReader> trackReader;
  std::unique_ptr<TTreeReaderValue<std::vector<o2::tpc::UnbinnedResid>>> residuals;
  std::unique_ptr<TTreeReaderValue<std::vector<o2::tpc::TrackDataCompact>>> trackRefs;
  std::unique_ptr<TTreeReaderValue<std::vector<o2::tpc::TrackData>>> tracks;

  void close()
  {
    tracks.reset();
    trackRefs.reset();
    residuals.reset();
    trackReader.reset();
    residualReader.reset();
    file.reset();
  }

  /// externalDetectors: also what the ITS clusters and the vertex need
  bool open(const std::string& fileName, bool externalDetectors)
  {
    close();
    file.reset(TFile::Open(fileName.c_str()));
    if (!file || file->IsZombie()) {
      LOGP(warning, "Cannot open {}", fileName);
      return false;
    }
    const bool isAlien = fileName.rfind("alien://", 0) == 0;
    if (isAlien) {
      file->SetBufferSize(4000000);
    }
    residualReader = std::make_unique<TTreeReader>("unbinnedResid", file.get());
    trackReader = std::make_unique<TTreeReader>("trackData", file.get());
    if (!residualReader->GetTree() || !trackReader->GetTree()) {
      LOGP(warning, "No unbinnedResid / trackData in {}", fileName);
      return false;
    }
    // read only what the fit uses; a large TTreeCache against slow remote reads
    for (auto* tree : {residualReader->GetTree(), trackReader->GetTree()}) {
      tree->SetCacheSize(256LL * 1024 * 1024);
      tree->AddBranchToCache("*", true);
    }
    std::vector<const char*> unused{"trackInfo.multStack*", "detInfo*"};
    if (!externalDetectors) {
      unused.insert(unused.end(), {"res.channel", "trackInfo.nExtDetResid"});
    }
    for (const char* branch : unused) {
      residualReader->GetTree()->SetBranchStatus(branch, 0);
    }
    for (const char* branch : {"trk.gid*", "trk.chi2TRD", "trk.deltaTOF", "trk.nTrkltsTRD", "trk.clAvailTOF", "trk.TRDTrkltSlope*", "trk.nExtDetResid",
                               "trk.clIdx.*", "trk.multStack*"}) {
      trackReader->GetTree()->SetBranchStatus(branch, 0);
    }
    residuals = std::make_unique<TTreeReaderValue<std::vector<o2::tpc::UnbinnedResid>>>(*residualReader, "res");
    trackRefs = std::make_unique<TTreeReaderValue<std::vector<o2::tpc::TrackDataCompact>>>(*residualReader, "trackInfo");
    tracks = std::make_unique<TTreeReaderValue<std::vector<o2::tpc::TrackData>>>(*trackReader, "trk");
    return true;
  }

  /// the TF's data copied into task; false on a read error (a storage hiccup leaves the values unusable without an exception)
  bool read(long entry, TFTask& task)
  {
    if (residualReader->SetEntry(entry) != TTreeReader::kEntryValid || trackReader->SetEntry(entry) != TTreeReader::kEntryValid) {
      return false;
    }
    const auto readOK = [](ROOT::Internal::TTreeReaderValueBase& value) { return value.GetReadStatus() == ROOT::Internal::TTreeReaderValueBase::kReadSuccess; };
    task.trackRefs = **trackRefs;
    if (!readOK(*trackRefs)) {
      return false;
    }
    task.residuals = **residuals;
    if (!readOK(*residuals)) {
      return false;
    }
    task.tracks = **tracks;
    if (!readOK(*tracks)) {
      return false;
    }
    return task.trackRefs.size() == task.tracks.size();
  }
};

/// one pass over the TFs (in the order given; a file is reopened when the next TF is in another one): this thread reads them, the fit
/// threads take chunks of their tracks, the residuals summed per voxel and charge in batches; the TFs read are appended to selected;
/// returns the tracks fitted
long runPass(const std::vector<const TFInfo*>& passTFs, const std::vector<std::string>& files, const FitContext& context, int nThreads,
             std::vector<std::array<ChargeSums, 2>>& sums, SelectedTFs& selected, const std::string& label)
{
  std::mutex sumsMutex;
  std::array<std::atomic<long>, NCounters> counters{};
  TFQueue queue;
  auto worker = [&]() {
    std::vector<Contribution> batch;
    batch.reserve(SumsBatchSize);
    std::vector<Cluster> clusters;
    std::array<long, NCounters> localCounters{};
    auto flush = [&]() {
      std::lock_guard<std::mutex> lock(sumsMutex);
      for (const auto& contribution : batch) {
        sums[contribution.voxel][contribution.charge].add(contribution);
      }
      batch.clear();
    };
    std::shared_ptr<TFTask> task;
    size_t firstTrack = 0;
    size_t lastTrack = 0;
    while (queue.next(task, firstTrack, lastTrack)) {
      for (size_t track = firstTrack; track < lastTrack; ++track) {
        processTrack(*task, track, context, localCounters, clusters, batch);
        if (batch.size() >= SumsBatchSize) {
          flush();
        }
      }
      task.reset();
      for (int counter = 0; counter < NCounters; ++counter) {
        counters[counter] += localCounters[counter];
        localCounters[counter] = 0;
      }
    }
    flush();
  };
  std::vector<std::thread> threads;
  for (int thread = 0; thread < nThreads; ++thread) {
    threads.emplace_back(worker);
  }

  long nTFs = 0;
  long nTFsReadError = 0;
  InputFile input;
  int openFile = -1;
  bool fileOpen = false;
  for (const TFInfo* tf : passTFs) {
    if (tf->file != openFile) {
      openFile = tf->file;
      LOGP(info, "Reading {}", files[tf->file]);
      fileOpen = input.open(files[tf->file], context.itsInFit > 0);
    }
    auto task = std::make_shared<TFTask>();
    task->info = tf;
    if (!fileOpen || !input.read(tf->entry, *task)) {
      if (fileOpen) {
        LOGP(warning, "Skipping entry {} of {}: read error", tf->entry, files[tf->file]);
      }
      ++nTFsReadError;
      continue;
    }
    queue.push(std::move(task));
    ++nTFs;
    selected.timeMS.push_back(static_cast<long>(tf->timeMS));
    selected.orbits.push_back(tf->orbit);
    selected.ctpLumi.push_back(tf->ctpLumi);
    selected.idcA.push_back(tf->idcA);
    selected.idcC.push_back(tf->idcC);
    selected.lumiScale.push_back(tf->lumiScale);
    if (tf->hasIDC) {
      selected.idc.push_back(0.5f * (tf->idcA + tf->idcC));
    }
    if (nTFs % 20 == 0) {
      LOGP(info, "TF {}{}: {} tracks selected, {} fitted, {} residuals, resident memory {:.2f} GB", nTFs, label, counters[Selected].load(),
           counters[Fitted].load(), counters[Residuals].load(), residentMemoryGB());
    }
  }
  queue.setDone();
  for (auto& thread : threads) {
    thread.join();
  }
  input.close();

  const long nFitted = counters[Fitted];
  LOGP(info,
       "fineMapCreatorGBL{}: TFs {} (read errors {}); tracks selected {}, fitted {} (failed {}); clusters out of row order {}; "
       "residuals {} (of de-weighted clusters {}), unbiased pull > {} {}, edge pads {}, outside the grid {}, without a residual {}; without the ITS "
       "clusters / vertex for the fit {}",
       label, nTFs, nTFsReadError, counters[Selected].load(), nFitted, counters[FitFailed].load(), counters[OutOfRowOrder].load(), counters[Residuals].load(),
       counters[ResidualsExcluded].load(), OutlierSigma, counters[Outliers].load(), counters[EdgePadResiduals].load(), counters[OutsideGrid].load(),
       counters[NoResidual].load(), counters[NoInnerPoints].load());
  LOGP(info, "fineMapCreatorGBL{}: fit passes per track {:.2f} ({:.1f} % at the maximum {}), resident memory {:.2f} GB", label,
       nFitted ? double(counters[FitPasses]) / nFitted : 0., nFitted ? 100. * counters[AtMaxPasses] / nFitted : 0., MaxFitPasses, residentMemoryGB());
  return nFitted;
}

} // namespace fine_map

void fineMapCreatorGBL(std::string fileInput = "files.txt",                                          // one residual file, or a .txt file with one file per line (local or alien://)
                       int runNumber = -1,                                                           // field, material LUT, orbit reset, beam type and TPC scalers from the CCDB
                       std::string fileOutput = "fineVoxRes.root",                                   //
                       std::string corrMap = "",                                                     // base correction: "" none, a TPCFastTransform file (local or alien://), "ccdb" (TPC/Calib/CorrectionMapV2 at mid-run) or "ccdb:<path>"
                       std::string corrMapDerivative = "",                                           // its derivative map, as corrMap ("ccdb": TPC/Calib/CorrectionMapDerivativeV2), scaled per TF with the IDC; "": the map as is
                       std::string fineMapInput = "",                                                // voxResTree of a previous pass, subtracted after the base correction
                       std::string trackSources = "ITS-TPC,ITS-TPC-TRD,ITS-TPC-TOF,ITS-TPC-TRD-TOF", //
                       std::string y2xBinning = "20",                                                // single number: uniform binning; otherwise bin boundaries
                       std::string z2xBinning = "28",                                                // single number: uniform binning; otherwise bin boundaries (scaled by maxZ2X)
                       float maxZ2X = 1.4f,                                                          //
                       long firstTFTime = -1,                                                        // TF time window [ms] as staticMapCreatorCPM.C (firstTFTime <= 0: all TFs)
                       long lastTFTime = -1,                                                         //
                       long entriesPerFile = -1,                                                     // the TFs of each file from firstEntry on (-1: all)
                       long firstEntry = 0,                                                          //
                       int minEntries = 15,                                                          // residuals for a voxel to be solved
                       bool writeSums = false,                                                       // also write the per-voxel sums (tree fineSums)
                       int nThreads = 8,                                                             // fit threads; one more thread reads the input
                       float smoothTau = 1000.f,                                                     // [um] prior rms of the second differences of the smoothed map (DS); 0: no smoothing, DS = D
                       float smoothTau0 = 1000.f,                                                    // [um] prior rms of a node of the smoothed map (the ridge towards 0)
                       int smoothRobust = 0,                                                         // Tukey reweightings of the voxels in the smoothing
                       int nPasses = 1,                                                              // 2: pass 1 on the even TFs in time, pass 2 on the odd ones with pass 1's smoothed map subtracted
                       int itsInFit = 0,                                                             // tests: 1 the ITS clusters in the track fit, 2 also the primary vertex; 0 TPC clusters only
                       bool fitSlopes = false)                                                       // the slopes of the voxel model from the fitted track at the cluster instead of the reference
{
  using namespace fine_map;
  const auto startTime = std::chrono::steady_clock::now();
  ROOT::EnableThreadSafety();
  if (runNumber < 0) {
    LOGP(fatal, "A run number is required");
  }
  nThreads = std::max(1, nThreads);
  if (smoothTau > 0.f && !(smoothTau0 > 0.f)) {
    LOGP(fatal, "smoothTau0 must be > 0 (the ridge keeps the voxels without entries defined)");
  }
  if (nPasses != 1 && nPasses != 2) {
    LOGP(fatal, "nPasses {}: 1 or 2", nPasses);
  }
  if (itsInFit < 0 || itsInFit > 2) {
    LOGP(fatal, "itsInFit {}: 0, 1 or 2", itsInFit);
  }
  if (nPasses == 2 && !(smoothTau > 0.f)) {
    LOGP(fatal, "Two passes need the smoothing (smoothTau > 0): pass 2 subtracts pass 1's smoothed map");
  }
  // a local CCDB snapshot keeps one object per path: the time-dependent TPC scalers need the server
  if (gSystem->Getenv("ALICEO2_CCDB_LOCALCACHE")) {
    LOGP(warning, "Unsetting ALICEO2_CCDB_LOCALCACHE: its snapshot holds one TPC/Calib/Scaler object for all times");
    gSystem->Unsetenv("ALICEO2_CCDB_LOCALCACHE");
  }

  // CCDB objects at mid-run
  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();
  ccdb.setCaching(true);
  ccdb.setFatalWhenNull(false);
  ccdb.setURL("http://alice-ccdb.cern.ch");
  const auto runDuration = ccdb.getRunDuration(runNumber);
  const long midRunTime = runDuration.first + (runDuration.second - runDuration.first) / 2;
  ccdb.setTimestamp(midRunTime);
  const double orbitResetMS = (*ccdb.get<std::vector<int64_t>>("CTP/Calib/OrbitReset"))[0] * 1e-3;
  const auto* magField = ccdb.get<o2::parameters::GRPMagField>("GLO/Config/GRPMagField");
  const auto* matLut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb.get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
  o2::base::Propagator::initFieldFromGRP(magField);
  o2::base::Propagator::Instance()->setMatLUT(matLut);
  o2::base::PropagatorD::Instance()->setMatLUT(matLut);
  const auto* grpLHC = ccdb.get<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF");
  const bool isPbPb = grpLHC->getBeamZ(o2::constants::lhc::BeamA) == 82 && grpLHC->getBeamZ(o2::constants::lhc::BeamC) == 82;
  const auto* scalerWeights = ccdb.get<o2::tpc::TPCScalerWeights>(o2::tpc::CDBTypeMap.at(o2::tpc::CDBType::CalScalerWeights));
  LOGP(info, "Run {}: orbit reset {:.0f} ms, L3 {} A, dipole {} A, PbPb {}, scaler weights {}", runNumber, orbitResetMS, magField->getL3Current(),
       magField->getDipoleCurrent(), isPbPb, scalerWeights && scalerWeights->isValid());

  // input files and maps may be on AliEn
  const auto files = getInputFileList(fileInput);
  auto isAlien = [](const std::string& name) { return name.rfind("alien://", 0) == 0; };
  if ((std::any_of(files.begin(), files.end(), isAlien) || isAlien(corrMap) || isAlien(corrMapDerivative) || isAlien(fineMapInput)) && !gGrid &&
      !TGrid::Connect("alien://")) {
    LOGP(fatal, "Failed to open the AliEn connection");
  }

  // base correction, previous fine map, binning
  BaseCorrection baseCorrection;
  if (!corrMap.empty()) {
    baseCorrection.map = loadTransform(corrMap, o2::tpc::CDBTypeMap.at(o2::tpc::CDBType::CalCorrMap), midRunTime);
    if (!corrMapDerivative.empty()) {
      baseCorrection.derivative = loadTransform(corrMapDerivative, o2::tpc::CDBTypeMap.at(o2::tpc::CDBType::CalCorrDerivMap), midRunTime);
    }
  } else if (!corrMapDerivative.empty()) {
    LOGP(fatal, "A derivative map needs a correction map");
  }
  // lumi scale as the reconstruction (LumiScaleMode::DerivativeMap): (IDC - IDC of the map) / IDC of the derivative
  o2::gpu::CorrectionMapsHelper lumiScaling;
  if (baseCorrection.derivative) {
    if (!(baseCorrection.map->getIDC() > 0.f) || !(baseCorrection.derivative->getIDC() > 0.f)) {
      LOGP(fatal, "The correction map and its derivative need an IDC (map {}, derivative {})", baseCorrection.map->getIDC(), baseCorrection.derivative->getIDC());
    }
    lumiScaling.setLumiScaleMode(o2::tpc::LumiScaleMode::DerivativeMap);
    lumiScaling.setMeanLumi(baseCorrection.map->getIDC(), false);
    lumiScaling.setMeanLumiRef(baseCorrection.derivative->getIDC(), false);
  }
  FineMapInput fineMap;
  if (!fineMapInput.empty() && !fineMap.load(fineMapInput)) {
    LOGP(fatal, "Cannot load the fine map {}", fineMapInput);
  }
  auto binning = makeBinning(y2xBinning, z2xBinning, maxZ2X);
  const int voxelsPerSector = binning->getNVoxelsPerSector();
  const RowGeometry geometry;

  // the TFs, their time and IDC (all CCDB access before the fit threads start)
  long nTFsOutsideWindow = 0;
  bool isMC = false;
  auto tfs = scanTFs(files, std::max(0L, firstEntry), entriesPerFile, orbitResetMS, firstTFTime, lastTFTime, nTFsOutsideWindow, isMC);
  if (firstTFTime > 0) {
    LOGP(info, "TF time window [{}, {}] ms: {} TFs inside, {} outside", firstTFTime, lastTFTime, tfs.size(), nTFsOutsideWindow);
  }
  fillIDCs(tfs, scalerWeights);
  // the ITS sensors' tracking frames (UnbinnedResid::getX / getAlpha), built before the fit threads start: MC is reconstructed with the ideal
  // geometry, data with the aligned one
  if (itsInFit > 0) {
    if (isMC) {
      ccdb.getForTimeStamp<TGeoManager>("GLO/Config/Geometry", IdealGeometryTimeStamp);
    }
    o2::tpc::UnbinnedResid::init(midRunTime);
    if (!gGeoManager) {
      LOGP(fatal, "No geometry for the ITS sensors");
    }
    o2::its::GeometryTGeo::Instance();
    LOGP(info, "ITS clusters{} in the track fit, {} geometry; slopes from the {}", itsInFit > 1 ? " and the primary vertex" : "",
         isMC ? "ideal (MC)" : "aligned", fitSlopes ? "fitted track" : "reference");
  }
  long nTFsNoIDC = 0;
  for (auto& tf : tfs) {
    if (baseCorrection.derivative && tf.hasIDC) {
      lumiScaling.setInstLumi(0.5f * (tf.idcA + tf.idcC), false);
      lumiScaling.updateLumiScale(false);
      tf.lumiScale = lumiScaling.getLumiScale();
    }
    nTFsNoIDC += !tf.hasIDC;
  }
  const bool needIDC = baseCorrection.derivative != nullptr;
  LOGP(info, "fineMapCreatorGBL: {} files, {} TFs ({} without IDC{}); correction map {}, derivative {}, fine map {}; grid y/x {}, z/x {} up to {} ({} voxels per sector); {} fit threads, {} passes",
       files.size(), tfs.size(), nTFsNoIDC, needIDC ? ": skipped" : "", corrMap.empty() ? "none" : corrMap, corrMapDerivative.empty() ? "none" : corrMapDerivative,
       fineMapInput.empty() ? "none" : fineMapInput, y2xBinning, z2xBinning, maxZ2X, voxelsPerSector, nThreads, nPasses);

  // the TFs to fit (with a derivative map those with IDC) in input order; two passes: pass 1 the even TFs in time, pass 2 the odd ones
  std::vector<const TFInfo*> fitTFs;
  for (const auto& tf : tfs) {
    if (!needIDC || tf.hasIDC) {
      fitTFs.push_back(&tf);
    }
  }
  std::array<std::vector<const TFInfo*>, 2> passTFs;
  if (nPasses == 2) {
    std::vector<const TFInfo*> byTime = fitTFs;
    std::stable_sort(byTime.begin(), byTime.end(), [](const TFInfo* first, const TFInfo* second) { return first->timeMS < second->timeMS; });
    std::vector<char> isOdd(tfs.size(), 0);
    for (size_t rank = 1; rank < byTime.size(); rank += 2) {
      isOdd[byTime[rank] - tfs.data()] = 1;
    }
    for (const TFInfo* tf : fitTFs) {
      passTFs[isOdd[tf - tfs.data()]].push_back(tf);
    }
    LOGP(info, "Two passes: {} TFs (even in time) in pass 1, {} (odd) in pass 2", passTFs[0].size(), passTFs[1].size());
  } else {
    passTFs[0] = fitTFs;
  }

  FitContext context;
  context.sources = GID::getSourcesMask(trackSources);
  context.baseCorrection = baseCorrection.map ? &baseCorrection : nullptr;
  context.fineMap = fineMap.nodes.empty() ? nullptr : &fineMap;
  context.binning = binning.get();
  context.geometry = &geometry;
  context.maxZ2X = maxZ2X;
  context.voxelsPerSector = voxelsPerSector;
  context.itsInFit = itsInFit;
  context.fitSlopes = fitSlopes;
  std::vector<std::array<ChargeSums, 2>> sums(size_t(NSectors) * voxelsPerSector);
  SelectedTFs selected;
  long nFitted = runPass(passTFs[0], files, context, nThreads, sums, selected, nPasses == 2 ? " pass 1" : "");

  // the smoothing; two passes: S1 of pass 1, pass 2 with S1 subtracted, S2 of pass 2, DS = S1 + S2
  auto smooth = [&](const std::vector<std::array<ChargeSums, 2>>& passSums, std::vector<double>& nodes, const std::string& label) {
    const auto smoothStart = std::chrono::steady_clock::now();
    const bool ok = smoothMap(*binning, passSums, minEntries, smoothTau * 1e-4, smoothTau0 * 1e-4, smoothRobust, nThreads, nodes);
    LOGP(info, "fineMapCreatorGBL{}: smoothing {:.1f} s", label, std::chrono::duration<double>(std::chrono::steady_clock::now() - smoothStart).count());
    return ok;
  };
  std::vector<std::pair<std::string, std::string>> extraInfo;
  std::vector<double> smoothed;
  std::vector<double> passOne;
  std::vector<std::array<ChargeSums, 2>> sumsPassOne;
  if (nPasses == 2) {
    if (!smooth(sums, passOne, " pass 1")) {
      LOGP(fatal, "Smoothing of pass 1 failed");
    }
    // pass 1's raw voxels for the cross-validation; its smoothed map subtracted in pass 2
    std::vector<float> rawPassOne(3 * sums.size(), 0.f);
    for (size_t voxel = 0; voxel < sums.size(); ++voxel) {
      ChargeSums total = sums[voxel][0];
      total.add(sums[voxel][1]);
      VoxelSolution solution;
      if (solveVoxel(total, minEntries, solution)) {
        for (int coordinate = 0; coordinate < 3; ++coordinate) {
          rawPassOne[3 * voxel + coordinate] = solution.value[coordinate];
        }
      }
    }
    FineMapInput passOneMap;
    passOneMap.binning = makeBinning(y2xBinning, z2xBinning, maxZ2X);
    passOneMap.voxelsPerSector = voxelsPerSector;
    passOneMap.nodes.assign(passOne.begin(), passOne.end());
    context.passOne = &passOneMap;
    if (writeSums) {
      sumsPassOne = sums;
    }
    std::fill(sums.begin(), sums.end(), std::array<ChargeSums, 2>{});
    nFitted += runPass(passTFs[1], files, context, nThreads, sums, selected, " pass 2");
    std::vector<double> passTwo;
    if (!smooth(sums, passTwo, " pass 2")) {
      LOGP(fatal, "Smoothing of pass 2 failed");
    }
    smoothed.resize(passOne.size());
    for (size_t node = 0; node < smoothed.size(); ++node) {
      smoothed[node] = passOne[node] + passTwo[node];
    }
    const auto validation = crossValidate(*binning, sums, passOne, rawPassOne);
    const double perMillion = 1e6 / std::max(1., validation.entries);
    const std::string validationText = fmt::format("pass 1 on pass 2's {:.0f} residuals, chi2 change per 1e6 residuals y / z: raw voxels {:.1f} / {:.1f}, smoothed {:.1f} / {:.1f}",
                                                   validation.entries, validation.raw[0] * perMillion, validation.raw[1] * perMillion,
                                                   validation.smoothed[0] * perMillion, validation.smoothed[1] * perMillion);
    LOGP(info, "Cross-validation: {}", validationText);
    extraInfo.emplace_back("nPasses", "2");
    extraInfo.emplace_back("passTFs", fmt::format("{} {}", passTFs[0].size(), passTFs[1].size()));
    extraInfo.emplace_back("crossValidation", validationText);
  } else if (smoothTau > 0.f && !smooth(sums, smoothed, "")) {
    LOGP(error, "Smoothing failed: DS = D");
    smoothed.clear();
  }
  extraInfo.emplace_back("smoothing", smoothed.empty() ? std::string("none") : fmt::format("tau {} um, tau0 {} um, robust {}", smoothTau, smoothTau0, smoothRobust));
  writeOutput(fileOutput, *binning, y2xBinning, z2xBinning, maxZ2X, sums, smoothed, passOne, sumsPassOne, extraInfo, nFitted, minEntries, writeSums, runNumber,
              orbitResetMS, isPbPb, firstTFTime, lastTFTime, selected);
  LOGP(info, "fineMapCreatorGBL: {:.1f} s", std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count());
}
