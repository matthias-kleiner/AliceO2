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
/// - the TPC clusters as the residual stream stores them (uncorrected, at the track's time), y snapped back to the 1/64-pad cluster grid
///   and, with clusterBiasFile, corrected for the cluster-position biases of the reconstruction (ClusterBiasCorrection.h: pad-phase
///   S-curve, sector-edge truncation), then corrected with the base correction
///   d = d_map + lumiScale * d_derivative (both evaluated at the cluster; lumiScale from the TF's IDC as in the reconstruction,
///   LumiScaleMode::DerivativeMap) and optionally with a previous fine map (the iteration); along the motion (a cluster on a row not
///   beyond the previous one is dropped);
/// - clusters near the stack boundaries, on the first / last rows and at the sector edges (on the pad plane) stay in the fit with their
///   errors x ExcludedSigmaScale: they hardly pull the track, and their unbiased residuals are summed like the others;
/// - reference: the ITS track propagated cluster to cluster with the PropagatorD (field map, LUT energy loss, the mass hypothesis of the
///   tracking = the stored track's PID); one GBL point per cluster with the step Jacobian (helix in the mean local Bz, central
///   differences) and a thin scatterer every ScatterEvery points with the x/X0 of its steps (Highland with the PID's mass, the charge and
///   the momentum of the pass's seed, the log term of the whole track); the seed takes GBL's correction at the first cluster every
///   pass; from FirstOutlierPass on, measurements beyond OutlierSigma of the previous pass's fit are dropped; a track stops once no new
///   outlier appears and GBL's correction is negligible;
/// - per cluster the unbiased residual r m^2 / e^2 (m: its error in the fit, e: GBL's residual error) with the variance sigma^2 + u,
///   u = v m^2 / (m^2 - v), v = m^2 - e^2 (sigma: the cluster's own error); the voxel from the best position y - r sigma^2 / (sigma^2 + u)
///   (cluster and prediction combined: binning by the measured position is Berkson); unbiased pull > OutlierSigma and the outermost
///   skipEdgePads pads of the row (argument, default 1; the best position moved onto the pad plane by the cluster's own correction)
///   not summed;
/// - bestEdgeBiasFile: the edge part of the cluster bias as a function of the BEST position (ClusterBiasCorrection.h text table,
///   E records only, dist = d_best = pads of the best position from the pad-plane edge; built from the fine pass's own edge tree of
///   no-distortion MC with makeEdgeBestTable.C): per summed residual r_y += sign(y) * B(d_best, unit, drift, |t|, class) with t the
///   voxel model's slope, applied after the outlier and edge cuts (the table is the mean bias of exactly these clusters). The cluster
///   itself, its weight and its voxel are not changed.
///   edgeBiasInFit: instead, the CLUSTERS are corrected inside the track fit (the fitted track is then no longer pulled by the
///   biased clusters): after every pass each TPC cluster gets y = y_0 + sign(y) * B(d_pred, ...), d_pred = pads of the unbiased track
///   position (this pass's fit without the cluster) from the pad-plane edge, t from this pass's reference; the next pass fits the
///   corrected clusters, and the fit is converged only once the corrections no longer change (1 um). The table must then be indexed
///   by d_pred (makeEdgeBestTable.C byPred, from the edge tree's dPred).
/// Voxel model ((dX, dY, dZ) = measured - true, as voxResTree): r_y = dY - t dX, r_z = dZ - c dX with t = dy/dx, c = dz/dx of the track;
/// dX, dY from the y channel, then dZ = <r_z> + <c> dX. Errors scaled by sqrt(chi2 / ndf) if > 1.
///
/// No smoothing here: DS = D (the voxel solve). The map is smoothed from the per-voxel sums (writeSums: tree fineSums) with smoothVoxelMap.C
/// (source "fineSums"), which also adds the sums of several runs (e.g. the chunks of a split GRID stage 3). Iterations: a second run with
/// fineMapInput = the smoothed map of the first, subtracted from the clusters after the base correction; with addBackFineMap the
/// subtracted correction is added back to each residual (r_y + S_Y - t S_X, r_z + S_Z - c S_X, the voxel model), so the output (D, sums)
/// is the full fine map S + what the first one left, else only the latter (e.g. injection tests: the recovered part of an injected map).
///
/// ITS clusters and primary vertex in the fit (itsInFit, for tests; ITS-TPC tracks with at least MinITSClusters ITS clusters): the points
/// of the residual stream (reference + residual; ITS in the sensor's tracking frame of the ideal geometry for MC, the aligned one for data)
/// in front of the TPC clusters, a scatterer at every point before the TPC and scatterer points every GapStep in the gap to the TPC, the
/// reference from the vertex (or the innermost ITS cluster) on; only the TPC clusters' unbiased residuals are summed. The track's kinks are
/// then fitted: the reference does not follow them, so the slopes of the voxel model should be the fitted track's (fitSlopes).
/// The scattering before the TPC has non-Gaussian tails (single large-angle and nuclear scatters; in MC 3-5 % of the gap kinks beyond 3
/// sigma of Highland): from the first outlier pass on, a kink there with |residual / error| > KinkTailPull (the larger of its two
/// directions) gets its variance x (pull / KinkTailPull)^2 for this track and the fit continues until no kink is beyond. Otherwise a hard
/// scatter in the gap bends the TPC part of the track towards the ITS direction, a bias of the slope regressor (in no-dist MC a radial dX
/// tilt OROC1 - OROC3 of 16 um, 7 um with the inflation; a scaled Gaussian core in the gap on top changes it by ~1 um).
/// kinkQAEvery > 0 (with itsInFit): per sampled track the kink residuals of the scatterers before the TPC with their prior (Highland) and
/// post-fit errors: unit Gaussian pulls residual / error if the modelled scattering is right (tree kinkQA in <output stem>.kinkQA.root).
/// edgeTreePads > 0: the edge tree (tree edgeTree in <output stem>.edgeTree.root, for studies of the cluster biases at the pad plane's
/// edge; one entry per fitted track of every edgeTreeEvery-th track of a TF with a TPC cluster within edgeTreePads pads of the edge at the
/// measured, best or true position): the track's TF (file, entry, CTP luminosity, IDC, lumi scale), occupancy (TrackDataCompact::multStack
/// per stack: the reconstruction's TPC occupancy map over the track's time span) and fit, and per such cluster its row, sector, nominal x,
/// stack and occupancy of its stack (mult), measured pad-plane y before any correction (padY0, snapped), after the cluster-bias correction (padY), the best position on the pad
/// plane (bestPadY), z, slope (tgSlp: the residual stream's = ITS-extrapolated reference's; tgSlpFit: the fitted track's, as the voxel
/// model; NaN without a residual), flags, de-weighting, the unbiased residuals with the cluster's and the track's errors, the y/x bin and what
/// the fine pass did with it (status: 0 summed, 1 outlier of the fit's last pass, 2 flagged, 3 no residual, 4 unbiased pull > OutlierSigma,
/// 5 within SkipEdgePads at the best position (edgeByPred: at the track without the cluster), 6 outside the grid, 7 not in the fit's last pass, 8 beside a dead pad) and the distances from the edge in pads
/// (d = (pad plane edge - side * y) / pad width, side = the sign of padY0; pad-0 centre at d = 0.5) of the measured, best and true position.
/// MC (trkMC in the input, e.g. --enable-mc residuals): the true y (yTrue) at the cluster's row, from TrackDataMC::parTPCIn swept row by row
/// (sector frame, nominal row x, material LUT, true PID) for a usable truth (parTPCIn, not fake, distTPCRef <= EdgeTreeMaxDistTPCRef).
/// Tests of the fit near the edge on such tracks (MC, needs the edge tree; a track without a complete truth is dropped): truthClusters: every
/// TPC cluster y = the truth + Gaussian noise of the cluster's error (no cluster bias at all); edgeRemovePads > 0: the clusters within this
/// many pads of the edge at the true position practically out of the fit (errors x RemovedSigmaScale), their unbiased residuals still there;
/// edgeWeightByTruth: the de-weighting within EdgePads decided at the true position instead of the measured one.
/// edgeFitPads > 0 replaces the de-weighting within EdgePads, which is decided by the cluster's own measured position: near that line clusters
/// with an inward fluctuation keep their weight and those with an outward one lose it, a selection on the cluster's noise that pulls the fit
/// inward, worst for tracks riding along the edge. Instead the clusters within edgeFitPads of the edge at the track WITHOUT the cluster are
/// practically out of the fit (errors x RemovedSigmaScale), re-decided after every pass; the first pass by the measured position with a margin
/// of one pad, which with MinFitClusters also drops the tracks riding along the edge over most of their length.
/// deadPadFile ("sector row pad", pad = pad number in the row; the dead channel map): a dead pad cuts the charge cloud like the pad plane's edge, so
/// the clusters beside it are shifted by up to a pad. Clusters within deadPadMargin pads of a dead pad at the track WITHOUT the cluster are out of the fit
/// (with edgeFitPads) and not summed (status 8), and the voxels there are left to the smoothing.
/// edgeByPred: the skipEdgePads cut, the y/x bin and the lookup in bestEdgeBiasFile by the track WITHOUT the cluster (d_pred) instead of the
/// best position (d_best), so that nothing there depends on the cluster's own noise; z/x stays at the best position. The table must then be
/// indexed by d_pred, from a fine pass with the same setting. Binning by the prediction alone is Berkson through the prediction's error.
///
/// TF time and IDC: the first orbit (records tree) and the CTPLumi of every TF are read first, with the time window cut of
/// staticMapCreatorCPM.C (firstTFTime / lastTFTime: the TFs of the coarse map of that window); the TFs are sorted in time and the TPC
/// scaler objects fetched from the CCDB in that order (each once); the IDC (mean of A and C side) and the lumi scale per TF go into a table.
/// Only this reader thread touches the input files and the CCDB (TGrid is not thread-safe): it reads one TF after the other into memory,
/// the fit threads take chunks of its tracks. With a derivative map, TFs without a valid scaler object are skipped; without one, the map
/// is applied as is (no IDC needed, e.g. the coarse map of the same time window).
///
/// Output: voxResTree (D = (dX, dY, dZ), E, flags DistDone | DispDone where solved; DS = D; binning, CTP / IDC means and medians in the user info) and the OrbitLumiInfo tree as staticMapCreatorCPM.C, with the per-TF IDC values; with writeSums also the tree fineSums
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
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
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
#include "CommonUtils/TreeStreamRedirector.h"
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

#include "ClusterBiasCorrection.h"
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
constexpr double RemovedSigmaScale = 1000.; ///< errors of the clusters of the test edgeRemovePads in the fit
// clusters de-weighted in the fit, on the pad plane
constexpr int StackRows = 2;    ///< rows on each side of a stack boundary
constexpr int EdgeRows = 4;     ///< the first / last rows
constexpr double EdgePads = 2.; ///< pads from the pad plane's edge, in every row
// residuals not summed
constexpr double SkipEdgePads = 1.; ///< default of skipEdgePads: the outermost pads of the row, at the best position
// edge tree (edgeTreePads > 0)
constexpr float EdgeTreeMaxDistTPCRef = 1.f; ///< [cm] MC truth used up to this distance of parTPCIn from the track's TPC reference
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
constexpr double KinkTailPull = 3.;                    ///< a kink before the TPC beyond this pull gets its scattering variance x (pull / KinkTailPull)^2
constexpr double MaxKinkInflation = 1e4;               ///< at most
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

/// the y/x grid edge of the TrackResiduals this runs with: "padplane" (the pad plane's edge in each row; O2 with the y/x edge fix) or
/// "deadzone" (x tan 10 deg - 1.5 cm; before). Written to the user info "y2xEdge" of the maps (missing in older files: "deadzone")
static std::string y2xEdgeDefinition(const o2::tpc::TrackResiduals& binning)
{
  const auto& mapper = o2::tpc::Mapper::instance();
  const float x = binning.getX(0);
  int row = 0;
  for (int other = 1; other < o2::tpc::param::NPadRows; ++other) {
    if (std::abs(o2::tpc::param::RowX[other] - x) < std::abs(o2::tpc::param::RowX[row] - x)) {
      row = other;
    }
  }
  const double halfWidth = 0.5 * mapper.getNumberOfPadsInRowSector(row) * mapper.getPadRegionInfo(o2::tpc::Mapper::REGION[row]).getPadWidth();
  return std::abs(binning.getMaxY2X(0) * x - halfWidth) < 1e-3 ? "padplane" : "deadzone";
}

/// warns if a map was binned with another y/x grid edge than the TrackResiduals this runs with (its bin centres then differ by up to
/// about a pad at the sector edges: read it with an O2 version of its own definition)
static void checkY2XEdge(const TList* userInfo, const o2::tpc::TrackResiduals& binning, const std::string& file)
{
  const auto* stored = userInfo ? userInfo->FindObject("y2xEdge") : nullptr;
  const std::string storedEdge = stored ? stored->GetTitle() : "deadzone";
  const std::string currentEdge = y2xEdgeDefinition(binning);
  if (storedEdge != currentEdge) {
    LOGP(warning, "{}: binned with the y/x grid edge \"{}\"{}, this O2 has \"{}\": the bin centres differ by up to about a pad at the sector "
                  "edges -- use an O2 version with the map's definition",
         file, storedEdge, stored ? "" : " (no y2xEdge: an older file)", currentEdge);
  }
}

/// pad-row x, pad width and the pad plane's edge (|y| of the outer edge of the outermost pad) per row
struct RowGeometry {
  std::array<double, NRows> x{};
  std::array<double, NRows> padWidth{};
  std::array<double, NRows> padPlaneEdge{};
  std::array<double, NRows> padAtYZero{}; ///< pad coordinate of y = 0 ((nPads - 1) / 2), pad 0 = the pad at the +y edge
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
      padAtYZero[row] = 0.5 * (nPads - 1);
    }
  }
};

/// the rows de-weighted in the track fit: stack boundaries (first rows of OROC1 / 2 / 3) and the first / last rows
bool isExcludedRow(int row)
{
  for (int boundary : {63, 97, 127}) {
    if (row >= boundary - StackRows && row < boundary + StackRows) {
      return true;
    }
  }
  return row < EdgeRows || row >= NRows - EdgeRows;
}

/// the clusters de-weighted in the track fit: the rows above and within EdgePads of the pad plane's edge (counted in pads per row: the pad
/// plane's edge does not follow a constant y/x)
bool isExcludedFromFit(int row, double padY, const RowGeometry& geometry)
{
  if (isExcludedRow(row)) {
    return true;
  }
  const double absY = std::abs(padY);
  return absY > geometry.padPlaneEdge[row] - EdgePads * geometry.padWidth[row];
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
/// map is applied with it.
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
    // an empty binning string is the default binning (staticMapCreatorCPM.C writes the arguments as given)
    auto has = [&](const char* key) { return tree->GetUserInfo()->FindObject(key) != nullptr; };
    if (!tree || !has("y2xBinning") || !has("z2xBinning") || info("maxZ2X").empty()) {
      LOGP(error, "Fine map {}: no voxResTree with the binning in its user info", fileName);
      return false;
    }
    binning = makeBinning(info("y2xBinning"), info("z2xBinning"), std::stof(info("maxZ2X")));
    checkY2XEdge(tree->GetUserInfo(), *binning, fileName);
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

  /// subtracts the map from (x, y, z); the correction subtracted in applied (optional)
  void apply(int sector, int row, double& x, double& y, double& z, double* applied = nullptr) const
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
    if (applied) {
      std::copy(correction, correction + 3, applied);
    }
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
      tf.idcA = scaler->getMeanScaler(tf.timeMS, o2::tpc::Side::A);
      tf.idcC = scaler->getMeanScaler(tf.timeMS, o2::tpc::Side::C);
      tf.hasIDC = tf.idcA > 0.f && tf.idcC > 0.f; // getMeanScaler: -1 without values, 0 with zero weights
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
  std::vector<o2::tpc::TrackDataMC> tracksMC; ///< the MC truth per track (edge tree with trkMC in the input; empty otherwise)
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

/// GBL's parameters (q/p, snp, tgl, y, z) of a state: a scattering kink changes only snp and tgl
std::array<double, 5> toGbl(const TrackParD& track)
{
  const double tgl = track.getTgl();
  return {track.getQ2Pt() / std::sqrt(1. + tgl * tgl), track.getSnp(), tgl, track.getY(), track.getZ()};
}

/// a state from GBL's parameters with the charge and mass hypothesis of track
TrackParD fromGbl(const std::array<double, 5>& parameters, double alpha, double x, const TrackParD& track)
{
  return TrackParD(x, alpha, {parameters[3], parameters[4], parameters[1], parameters[2], parameters[0] * std::sqrt(1. + parameters[2] * parameters[2])},
                   track.getAbsCharge(), track.getPID());
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
      TrackParD track = fromGbl(parameters, from.getAlpha(), from.getX(), from);
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
  double padY = 0.; ///< measured y on the pad plane (snapped, cluster-bias corrected), before the map corrections
  double tgSlp = 0.;
  bool tgSlpValid = true;
  bool excluded = false; ///< de-weighted in the fit
  bool removed = false;  ///< test edgeRemovePads: practically out of the fit, the unbiased residual still there
  bool flagged = false;  ///< a TPC cluster flag of clusterFlagMask: de-weighted in the fit (excluded) and not summed
  double sigmaY = 0.;    ///< the cluster's own error
  double sigmaZ = 0.;
  double fineCorrection[3] = {0., 0., 0.}; ///< TPC: the fineMapInput correction subtracted from (x, y, z)
  // edge tree only
  int flags = 0;                                           ///< TPC cluster flags (UnbinnedResid::channel)
  double xNominal = 0.;                                    ///< the row's x
  double padY0 = 0.;                                       ///< measured y on the pad plane (snapped), before any correction
  double z0 = 0.;                                          ///< measured z, before any correction
  double yUncorrected = 0.;                                ///< y as fitted without the in-fit edge correction (after the map corrections)
  double edgeShift = 0.;                                   ///< edgeBiasInFit: the edge correction currently added to y [cm]
  double yTrue = std::numeric_limits<double>::quiet_NaN(); ///< MC: the true y at the row (sector frame, nominal x)
};

/// the fit of a track: active = the clusters in the final pass, reference the states at their planes, outlier the last pass's flags
struct FitResult {
  std::vector<int> active;
  std::vector<TrackParD> reference;
  std::vector<char> outlier;
  std::unique_ptr<gbl::GblTrajectory> trajectory;
  double qOverPt = 0.;
  int passes = 0;
  std::vector<float> kinkInflation; ///< per active point: the tail inflation of its scattering variance (1: none)
};

/// the dead channel map per sector and row, in pads from the -y edge of the pad plane
struct DeadPads {
  std::vector<char> dead;  ///< [sector][row][pad]
  int maxPads = 0;
  double margin = 0.;
  const RowGeometry* geometry = nullptr;
  bool active() const { return !dead.empty(); }
  /// the list "sector row pad", pad number in the row: pad 0 at the -y edge on the A side, at the +y edge on the C side (checked against the clusters per pad)
  void load(const std::string& fileName, double marginPads, const RowGeometry& rowGeometry)
  {
    geometry = &rowGeometry;
    margin = marginPads;
    maxPads = 0;
    for (int row = 0; row < NRows; ++row) {
      maxPads = std::max(maxPads, int(std::round(2. * rowGeometry.padPlaneEdge[row] / rowGeometry.padWidth[row])));
    }
    dead.assign(size_t(36) * NRows * maxPads, 0);
    std::ifstream in(fileName);
    int sector = 0;
    int row = 0;
    int pad = 0;
    long n = 0;
    while (in >> sector >> row >> pad) {
      const int nPads = int(std::round(2. * rowGeometry.padPlaneEdge[row] / rowGeometry.padWidth[row]));
      dead[(size_t(sector) * NRows + row) * maxPads + (sector < 18 ? pad : nPads - 1 - pad)] = 1;
      ++n;
    }
    if (n == 0) {
      LOGP(fatal, "no dead pads in {}", fileName);
    }
    LOGP(info, "{} dead pads from {}: clusters within {} pads of them at the track without the cluster out of the fit and not summed", n, fileName, margin);
  }
  /// within margin + extra pads of a dead pad at pad-plane y
  bool near(int sector, int row, double padY, double extra = 0.) const
  {
    if (dead.empty()) {
      return false;
    }
    const double width = geometry->padWidth[row];
    const double position = padY / width + geometry->padPlaneEdge[row] / width; // pads from the -y edge
    const int nPads = int(std::round(2. * geometry->padPlaneEdge[row] / width));
    const int first = std::max(0, int(std::floor(position - margin - extra)) - 1);
    const int last = std::min(nPads - 1, int(std::floor(position + margin + extra)) + 1);
    const char* flags = &dead[(size_t(sector) * NRows + row) * maxPads];
    for (int pad = first; pad <= last; ++pad) {
      if (flags[pad] && position >= pad - margin - extra && position < pad + 1 + margin + extra) {
        return true;
      }
    }
    return false;
  }
};

/// edgeBiasInFit: the edge part of the cluster bias vs the unbiased track position, applied to the clusters between the fit passes
struct EdgeFitCorrection {
  const o2::tpc::cluster_bias::ClusterBiasCorrection* table = nullptr;
  const RowGeometry* geometry = nullptr;
  double pads = 0.; ///< edgeFitPads: clusters within this many pads of the edge at the track without the cluster out of the fit
  const DeadPads* deadPads = nullptr; ///< with edgeFitPads: also the clusters beside a dead pad at the track without the cluster out of the fit
};

bool fitTrack(std::vector<Cluster>& clusters, TrackParD seed, FitResult& result, const EdgeFitCorrection* edgeFit = nullptr)
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
  std::vector<double> kinkInflation(clusters.size(), 1.);
  Eigen::VectorXd residuals(2);
  Eigen::VectorXd measurementErrors(2);
  // GBL fills these without resizing them: two directions of a thin scatterer
  Eigen::VectorXd kinkResiduals(2);
  Eigen::VectorXd kinkPriors(2);
  Eigen::VectorXd kinkErrors(2);
  Eigen::VectorXd kinkDownWeights(2);
  for (int pass = 0; pass < MaxFitPasses; ++pass) {
    // multiple scattering with the momentum of the current seed (pass 0: the ITS track, then GBL's fit), the mass hypothesis of the
    // tracking (the PID of the stored track) and the charge
    const double momentum = seed.getP();
    const double mass = seed.getPID().getMass();
    const double beta = momentum / std::sqrt(momentum * momentum + mass * mass);
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
        next = fromGbl(parameters, next.getAlpha(), next.getX(), next);
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
    const double theta0 = 0.0136 * seed.getAbsCharge() / (beta * momentum) * (1. + 0.038 * std::log(std::max(totalXOverX0, 1e-12)));
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
        const double scale = cluster.removed ? RemovedSigmaScale : (cluster.excluded ? ExcludedSigmaScale : 1.);
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
        const double inflatedTheta2 = theta2 * (point < nInner ? kinkInflation[result.active[point]] : 1.);
        gblPoint.addScatterer(kinked ? Eigen::Vector2d(-kink[0], -kink[1]) : Eigen::Vector2d::Zero(),
                              Eigen::Vector2d(1. / (inflatedTheta2 * (1. - snp * snp) * secLambda2), 1. / (inflatedTheta2 * secLambda2 * secLambda2)));
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
    // tail kinks before the TPC: variance x (pull / KinkTailPull)^2 for the next pass
    bool newInflation = false;
    if (pass + 1 >= FirstOutlierPass) {
      for (size_t point = 1; point < nInner && point + 1 < nPoints; ++point) {
        unsigned int numData = 0;
        if (result.trajectory->getScatResults(point + 1, numData, kinkResiduals, kinkPriors, kinkErrors, kinkDownWeights) != 0 || numData < 2) {
          continue;
        }
        double pull = 0.;
        for (int direction = 0; direction < 2; ++direction) {
          if (kinkErrors(direction) > 0.) {
            pull = std::max(pull, std::abs(kinkResiduals(direction) / kinkErrors(direction)));
          }
        }
        auto& inflation = kinkInflation[result.active[point]];
        if (pull > KinkTailPull && inflation < MaxKinkInflation) {
          inflation = std::min(MaxKinkInflation, inflation * (pull / KinkTailPull) * (pull / KinkTailPull));
          newInflation = true;
        }
      }
    }
    // edgeBiasInFit: the clusters' edge correction at this pass's unbiased track position, for the next pass
    bool newEdgeShift = false;
    if (edgeFit) {
      Eigen::VectorXd residualErrors(2);
      Eigen::VectorXd downWeights(2);
      const auto& geometry = *edgeFit->geometry;
      for (size_t point = nInner; point < nPoints; ++point) {
        auto& cluster = clusters[result.active[point]];
        if (cluster.type != Cluster::TPC) {
          continue;
        }
        unsigned int numData = 0;
        result.trajectory->getMeasResults(point + 1, numData, residuals, measurementErrors, residualErrors, downWeights);
        const double fitError2 = measurementErrors(0) * measurementErrors(0);
        const double residualError2 = residualErrors(0) * residualErrors(0);
        if (numData < 2 || !(residualError2 > 1e-4 * fitError2)) {
          continue;
        }
        const double yPredicted = cluster.y - residuals(0) * fitError2 / residualError2; // the track without the cluster
        const double predictedPadY = cluster.padY + (yPredicted - cluster.yUncorrected);
        const double width = geometry.padWidth[cluster.row];
        const double side = cluster.padY0 < 0. ? -1. : 1.;
        const double dPredicted = (geometry.padPlaneEdge[cluster.row] - side * predictedPadY) / width;
        if (edgeFit->pads > 0.) {
          const bool removed = dPredicted < edgeFit->pads || (edgeFit->deadPads && edgeFit->deadPads->near(cluster.sector, cluster.row, predictedPadY));
          newEdgeShift |= removed != cluster.removed;
          cluster.removed = removed;
        }
        if (!edgeFit->table) {
          continue;
        }
        const double padMeas = (geometry.padPlaneEdge[cluster.row] - std::abs(cluster.padY0)) / width - 0.5;
        const int cls = o2::tpc::cluster_bias::classOfFlags(cluster.flags, std::abs(padMeas - std::round(padMeas)) < 0.25 / 64.);
        const double snp = reference[point].getSnp();
        const double slope = snp / std::sqrt(std::max(1e-6, (1. - snp) * (1. + snp)));
        const double shift = side * edgeFit->table->inwardBias(edgeFit->table->unitOfRow(cluster.row), cls, o2::tpc::cluster_bias::driftLength(cluster.z0),
                                                               std::abs(slope), dPredicted);
        newEdgeShift |= std::abs(shift - cluster.edgeShift) > 1e-4;
        cluster.edgeShift = shift;
        cluster.y = cluster.yUncorrected + shift;
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
    seed = fromGbl(parameters, reference.front().getAlpha(), reference.front().getX(), reference.front());
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
    // converged: from the first outlier pass on, no new outlier or inflated kink and GBL's correction at the first point negligible (1 um
    // in y / z, 1e-6 in snp / tgl, 1e-4 of q/p): the reference is the fit, further passes change nothing
    if (pass >= FirstOutlierPass && !newOutliers && !newInflation && !newEdgeShift && std::abs(corrections(3)) < 1e-4 && std::abs(corrections(4)) < 1e-4 &&
        std::abs(corrections(1)) < 1e-6 && std::abs(corrections(2)) < 1e-6 && std::abs(corrections(0)) < 1e-4 * std::abs(parameters[0]) && kinkUpdate < 1e-6) {
      break;
    }
  }
  result.qOverPt = seed.getQ2Pt();
  result.kinkInflation.clear();
  for (int index : result.active) {
    result.kinkInflation.push_back(kinkInflation[index]);
  }
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
               DeadPadResiduals,
               FlaggedResiduals,
               OutsideGrid,
               NoResidual,
               NoInnerPoints,
               KinkQAWritten,
               KinkQASkipped,
               EdgeTreeTracks,
               EdgeTreeClusters,
               EdgeTreeTruth,
               EdgeTreeSweepFailed,
               EdgeBiasCorrected,
               NCounters };

/// what the fit threads share read-only
struct FitContext {
  GID::mask_t sources;
  const BaseCorrection* baseCorrection = nullptr; ///< nullptr: no base correction
  const FineMapInput* fineMap = nullptr;          ///< nullptr: no previous fine map
  bool addBackFineMap = false;                    ///< the fineMap correction added back to the summed residuals
  const TrackResiduals* binning = nullptr;
  const RowGeometry* geometry = nullptr;
  float maxZ2X = 0.f;
  int clusterFlagMask = 0; ///< TPC cluster flags (UnbinnedResid::channel) whose clusters are de-weighted in the fit and not summed
  const o2::tpc::cluster_bias::ClusterBiasCorrection* clusterBias = nullptr; ///< nullptr: no per-cluster bias correction
  int voxelsPerSector = 0;
  int itsInFit = 0;                                  ///< 1: the ITS clusters in the fit, 2: also the primary vertex
  bool fitSlopes = false;                            ///< the slopes of the voxel model from the fitted track instead of the reference
  o2::utils::TreeStreamRedirector* kinkQA = nullptr; ///< kink QA (itsInFit): one streamer for all fit threads
  std::mutex* kinkQAMutex = nullptr;
  int kinkQAEvery = 0;
  o2::utils::TreeStreamRedirector* edgeTree = nullptr; ///< edge tree: one streamer for all fit threads
  std::mutex* edgeTreeMutex = nullptr;
  double edgeTreePads = 0.; ///< clusters within this many pads of the pad plane's edge
  int edgeTreeEvery = 1;    ///< every n-th track of a TF
  const o2::tpc::cluster_bias::ClusterBiasCorrection* bestEdgeBias = nullptr; ///< edge part vs the best position, on the summed residual
  double skipEdgePads = SkipEdgePads;                                         ///< residuals within this many pads of the edge not summed
  EdgeFitCorrection edgeFit;                                                  ///< edgeBiasInFit: the table on the clusters inside the fit
  bool truthClusters = false;                                                 ///< test: the TPC clusters at the MC truth + noise
  double edgeRemovePads = 0.;                                                 ///< test: clusters within this many pads of the edge (truth) out of the fit
  bool edgeWeightByTruth = false;                                             ///< test: the de-weighting within EdgePads decided at the truth
  bool edgeByPred = false;                                                    ///< edge cut, y/x bin and table lookup by the track without the cluster
  DeadPads deadPads;                                                          ///< deadPadFile: clusters beside dead pads out
  bool lastRowsInFit = false;                                                 ///< the last EdgeRows rows keep their weight in the fit (with edgeFitPads)
};

/// kink QA: the scatterers before the TPC (ITS layers, gap) of a fitted track (x, z of each: the material crossed depends on z / x, not only
/// on tgl, with the vertex z spread), per scatterer and direction (snp, tgl) GBL's kink residual
/// (0 - fitted kink), its prior error (Highland) and its error after the fit: residual / error is unit Gaussian if the scattering model is
/// right; the summed (residual / prior)^2 of the gap and its expectation sum (error / prior)^2; the prior includes the tail inflation
/// (inflation per scatterer, 1: none). An entry is skipped while another thread writes.
void writeKinkQA(const FitResult& fit, const std::vector<Cluster>& clusters, const FitContext& context, std::array<long, NCounters>& counters)
{
  std::vector<float> x;
  std::vector<float> z;
  std::vector<int> type;
  std::vector<float> inflation;
  std::array<std::vector<float>, 2> residual;
  std::array<std::vector<float>, 2> prior;
  std::array<std::vector<float>, 2> error;
  double chi2Gap = 0.;
  double expectedGap = 0.;
  double chi2Inner = 0.;
  double expectedInner = 0.;
  // GBL fills these without resizing them: two directions of a thin scatterer
  Eigen::VectorXd residuals(2);
  Eigen::VectorXd priorErrors(2);
  Eigen::VectorXd errors(2);
  Eigen::VectorXd downWeights(2);
  for (size_t point = 1; point + 1 < fit.active.size(); ++point) {
    const auto& cluster = clusters[fit.active[point]];
    if (cluster.type == Cluster::TPC) {
      break;
    }
    unsigned int numData = 0;
    if (fit.trajectory->getScatResults(point + 1, numData, residuals, priorErrors, errors, downWeights) != 0 || numData < 2) {
      continue;
    }
    x.push_back(cluster.x);
    z.push_back(fit.reference[point].getZ());
    type.push_back(cluster.type);
    inflation.push_back(fit.kinkInflation[point]);
    for (int direction = 0; direction < 2; ++direction) {
      residual[direction].push_back(residuals(direction));
      prior[direction].push_back(priorErrors(direction));
      error[direction].push_back(errors(direction));
      const double chi2 = std::pow(residuals(direction) / priorErrors(direction), 2);
      const double expected = std::pow(errors(direction) / priorErrors(direction), 2);
      chi2Inner += chi2;
      expectedInner += expected;
      if (cluster.type == Cluster::Material) {
        chi2Gap += chi2;
        expectedGap += expected;
      }
    }
  }
  if (x.empty()) {
    return;
  }
  std::unique_lock<std::mutex> lock(*context.kinkQAMutex, std::try_to_lock);
  if (!lock.owns_lock()) {
    ++counters[KinkQASkipped];
    return;
  }
  const auto& first = fit.reference.front();
  float p = first.getP();
  float pt = first.getPt();
  float tgl = first.getTgl();
  int charge = fit.qOverPt > 0. ? 1 : -1;
  int pid = first.getPID().getID();
  (*context.kinkQA) << "kinkQA"
                    << "p=" << p << "pt=" << pt << "tgl=" << tgl << "charge=" << charge << "pid=" << pid << "x=" << x << "z=" << z << "inflation=" << inflation << "type=" << type
                    << "residualSnp=" << residual[0] << "residualTgl=" << residual[1] << "priorSnp=" << prior[0] << "priorTgl=" << prior[1]
                    << "errorSnp=" << error[0] << "errorTgl=" << error[1] << "chi2Gap=" << chi2Gap << "expectedGap=" << expectedGap
                    << "chi2Inner=" << chi2Inner << "expectedInner=" << expectedInner << "\n";
  ++counters[KinkQAWritten];
}

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

/// edge tree: per cluster what the fine pass did with it (the status codes of the edge tree), its best position and unbiased residuals
struct EdgeRecord {
  int status = 7; ///< not in the fit's last pass
  double bestPadY = std::numeric_limits<double>::quiet_NaN();
  double residualY = std::numeric_limits<double>::quiet_NaN();
  double residualZ = std::numeric_limits<double>::quiet_NaN();
  double trackSigmaY = std::numeric_limits<double>::quiet_NaN(); ///< the track's error at the cluster without the cluster
  double tgSlpFit = std::numeric_limits<double>::quiet_NaN();    ///< dy/dx of the fitted track at the cluster (the voxel model's slope)
  double edgeCorr = 0.;                                          ///< bestEdgeBias correction added to the summed residual [cm]
  double predPadY = std::numeric_limits<double>::quiet_NaN();    ///< the unbiased track position (fit without the cluster) on the pad plane
  double edgeShift = 0.;                                         ///< edgeBiasInFit: the correction of the cluster in the last fit pass [cm]
  int binY2X = -1;
};

/// edge tree: one entry per track with its TPC clusters within edgeTreePads pads of the pad plane's edge (measured, best or true position)
void writeEdgeTree(const TFTask& task, size_t trackIndex, const FitResult& fit, const std::vector<Cluster>& clusters, const std::vector<EdgeRecord>& records,
                   int truthState, bool sweepOK, const FitContext& context, std::array<long, NCounters>& counters)
{
  std::vector<int> row;
  std::vector<int> stack;
  std::vector<int> sector;
  std::vector<int> flags;
  std::vector<int> excluded;
  std::vector<int> status;
  std::vector<int> binY2X;
  std::vector<float> x;
  std::vector<float> padY0;
  std::vector<float> padY;
  std::vector<float> bestPadY;
  std::vector<float> yTrue;
  std::vector<float> z;
  std::vector<float> tgSlp;
  std::vector<float> tgSlpFit;
  std::vector<float> edgeCorr;
  std::vector<float> edgeShift;
  std::vector<float> predPadY;
  std::vector<float> dPred;
  std::vector<float> residualY;
  std::vector<float> residualZ;
  std::vector<float> sigmaY;
  std::vector<float> trackSigmaY;
  std::vector<float> dMeas;
  std::vector<float> dBest;
  std::vector<float> dTrue;
  std::vector<float> mult;
  const auto& trackRef = task.trackRefs[trackIndex];
  std::vector<float> multStack(4, 0.f);
  for (int iStack = 0; iStack < 4; ++iStack) {
    multStack[iStack] = trackRef.getMultStackPacked(iStack) > 0 ? trackRef.getMultStack(iStack) : 0.f;
  }
  const auto& geometry = *context.geometry;
  for (size_t index = 0; index < clusters.size(); ++index) {
    const auto& cluster = clusters[index];
    if (cluster.type != Cluster::TPC) {
      continue;
    }
    const auto& record = records[index];
    const double width = geometry.padWidth[cluster.row];
    const double edge = geometry.padPlaneEdge[cluster.row];
    const double side = cluster.padY0 < 0. ? -1. : 1.;
    const double distMeas = (edge - std::abs(cluster.padY0)) / width;
    const double distBest = (edge - side * record.bestPadY) / width; // NaN without a best position
    const double distTrue = (edge - side * cluster.yTrue) / width;   // NaN without the truth
    if (!(distMeas < context.edgeTreePads) && !(distBest < context.edgeTreePads) && !(distTrue < context.edgeTreePads)) {
      continue;
    }
    row.push_back(cluster.row);
    const int clusterStack = cluster.row < 63 ? 0 : (cluster.row < 97 ? 1 : (cluster.row < 127 ? 2 : 3));
    stack.push_back(clusterStack);
    mult.push_back(multStack[clusterStack]);
    sector.push_back(cluster.sector);
    flags.push_back(cluster.flags);
    excluded.push_back(cluster.excluded);
    status.push_back(record.status);
    binY2X.push_back(record.binY2X);
    x.push_back(cluster.xNominal);
    padY0.push_back(cluster.padY0);
    padY.push_back(cluster.padY);
    bestPadY.push_back(record.bestPadY);
    yTrue.push_back(cluster.yTrue);
    z.push_back(cluster.z0);
    tgSlp.push_back(cluster.tgSlpValid ? cluster.tgSlp : std::numeric_limits<double>::quiet_NaN());
    tgSlpFit.push_back(record.tgSlpFit);
    edgeCorr.push_back(record.edgeCorr);
    edgeShift.push_back(record.edgeShift);
    predPadY.push_back(record.predPadY);
    dPred.push_back((edge - side * record.predPadY) / width); // NaN without a residual
    residualY.push_back(record.residualY);
    residualZ.push_back(record.residualZ);
    sigmaY.push_back(cluster.sigmaY);
    trackSigmaY.push_back(record.trackSigmaY);
    dMeas.push_back(distMeas);
    dBest.push_back(distBest);
    dTrue.push_back(distTrue);
  }
  if (row.empty()) {
    return;
  }
  const auto& info = *task.info;
  int file = info.file;
  int entry = static_cast<int>(info.entry);
  int track = static_cast<int>(trackIndex);
  int source = trackRef.sourceId;
  float ctpLumi = info.ctpLumi;
  float idc = info.hasIDC ? 0.5f * (info.idcA + info.idcC) : -1.f;
  float lumiScale = info.lumiScale;
  const auto& first = fit.reference.front();
  float qOverPt = fit.qOverPt;
  float p = first.getP();
  float tgl = first.getTgl();
  int passes = fit.passes;
  int mc = truthState;
  int truthComplete = sweepOK;
  float pTrue = -1.f;
  float distTPCRef = -1.f;
  int pidTrue = -1;
  if (truthState > 0) {
    const auto& truthIn = task.tracksMC[trackIndex];
    pTrue = truthIn.parTPCIn.getP();
    distTPCRef = truthIn.distTPCRef;
    pidTrue = truthIn.parTPCIn.getPID().getID();
  }
  {
    std::lock_guard<std::mutex> lock(*context.edgeTreeMutex);
    (*context.edgeTree) << "edgeTree"
                        << "file=" << file << "entry=" << entry << "track=" << track << "source=" << source << "ctpLumi=" << ctpLumi << "idc=" << idc
                        << "lumiScale=" << lumiScale << "multStack=" << multStack << "qOverPt=" << qOverPt << "p=" << p << "tgl=" << tgl
                        << "passes=" << passes << "mc=" << mc << "truthComplete=" << truthComplete << "pTrue=" << pTrue << "distTPCRef=" << distTPCRef
                        << "pidTrue=" << pidTrue << "row=" << row << "stack=" << stack << "mult=" << mult << "sector=" << sector << "flags=" << flags << "excluded=" << excluded
                        << "status=" << status << "binY2X=" << binY2X << "x=" << x << "padY0=" << padY0 << "padY=" << padY << "bestPadY=" << bestPadY
                        << "yTrue=" << yTrue << "z=" << z << "tgSlp=" << tgSlp << "tgSlpFit=" << tgSlpFit << "edgeCorr=" << edgeCorr << "edgeShift=" << edgeShift << "predPadY=" << predPadY << "dPred=" << dPred << "residualY=" << residualY << "residualZ=" << residualZ
                        << "sigmaY=" << sigmaY << "trackSigmaY=" << trackSigmaY << "dMeas=" << dMeas << "dBest=" << dBest << "dTrue=" << dTrue << "\n";
  }
  ++counters[EdgeTreeTracks];
  counters[EdgeTreeClusters] += row.size();
  counters[EdgeTreeTruth] += truthState == 1;
  counters[EdgeTreeSweepFailed] += truthState == 1 && !sweepOK;
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
  // edge tree: a sampled track with a TPC cluster near the pad plane's edge (measured; one pad of margin for the best / true position)
  bool edgeTrack = false;
  if (context.edgeTree && trackIndex % context.edgeTreeEvery == 0) {
    for (unsigned int residualIndex = trackRef.idxFirstResidual; residualIndex < trackRef.idxFirstResidual + trackRef.nResiduals && !edgeTrack; ++residualIndex) {
      const auto& residual = task.residuals[residualIndex];
      if (residual.row >= NRows) {
        continue;
      }
      const double y = residual.y * o2::tpc::param::MaxY / 0x7fff + decodeResidual(residual.dy);
      edgeTrack = context.geometry->padPlaneEdge[residual.row] - std::abs(y) < (context.edgeTreePads + 1.) * context.geometry->padWidth[residual.row];
    }
  }
  // MC: the truth swept along the rows (edge tree only)
  int truthState = 0; // 0 no truth, 1 swept, 2 present but not usable
  TrackParD truth;
  if (edgeTrack && trackIndex < task.tracksMC.size()) {
    const auto& mc = task.tracksMC[trackIndex];
    if (mc.hasTPCIn() && !mc.isFakeITSTPC() && mc.label.isValid() && !mc.label.isFake() && mc.distTPCRef >= 0.f && mc.distTPCRef <= EdgeTreeMaxDistTPCRef) {
      const auto& in = mc.parTPCIn;
      truth = TrackParD(in.getX(), in.getAlpha(), {in.getY(), in.getZ(), in.getSnp(), in.getTgl(), in.getQ2Pt()}, in.getAbsCharge(), in.getPID());
      truthState = 1;
    } else {
      truthState = 2;
    }
  }
  bool sweepOK = truthState == 1;
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
    // the stored cluster y (track y + dy, ideal pad-plane coordinates) back on the ClusterNative grid of 1/64 pad: UnbinnedResid packed
    // y and dy into shorts by truncation towards zero (rounding since the O2 fix), which moves y + dy by up to 21 um depending on the
    // sign of dy; binned on sub-pad scales that splits the clusters of a grid point by the sign of dy. Exact while the packing error is
    // below half a grid step (>= 32 um). Before any correction (a corrected position is not on the grid)
    {
      const double padWidth = context.geometry->padWidth[cluster.row];
      const double padAtYZero = context.geometry->padAtYZero[cluster.row];
      const double pad = padAtYZero - cluster.y / padWidth;
      cluster.y = (padAtYZero - std::round(pad * 64.) / 64.) * padWidth;
    }
    cluster.z = residual.z * o2::tpc::param::MaxZ / 0x7fff + decodeResidual(residual.dz);
    cluster.flags = residual.channel;
    cluster.xNominal = cluster.x;
    cluster.padY0 = cluster.y;
    cluster.z0 = cluster.z;
    if (edgeTrack) {
      if (sweepOK) {
        sweepOK = rotateTo(truth, cluster.alpha) && o2::base::PropagatorD::Instance()->PropagateToXBxByBz(truth, cluster.x, 0.99, 2., MatCorrType::USEMatCorrLUT);
        if (sweepOK) {
          cluster.yTrue = truth.getY();
          if (context.truthClusters) {
            thread_local std::mt19937 generator{std::random_device{}()};
            thread_local std::normal_distribution<double> gauss;
            const double snp = truth.getSnp();
            double sigmaY = 0.;
            double sigmaZ = 0.;
            clusterErrors(cluster.row, cluster.z, snp / std::sqrt(std::max(1e-6, (1. - snp) * (1. + snp))), truth.getTgl(), sigmaY, sigmaZ);
            cluster.y = cluster.yTrue + sigmaY * gauss(generator);
            cluster.padY0 = cluster.y;
          }
          cluster.removed = context.geometry->padPlaneEdge[cluster.row] - std::abs(cluster.yTrue) < context.edgeRemovePads * context.geometry->padWidth[cluster.row];
        }
      }
    }
    // position-only cluster or saturated tgSlp: the reference's angle instead
    cluster.tgSlpValid = !residual.isPositionOnly() && !residual.isTgSlpClamped();
    cluster.tgSlp = residual.tgSlp * o2::tpc::param::MaxTgSlp / 0x7fff;
    // per-cluster correction of the cluster-position biases of the reconstruction (pad-phase S-curve, sector-edge truncation;
    // ClusterBiasCorrection.h) on the snapped position, before the edge cut and the binning
    if (context.clusterBias) {
      cluster.y += context.clusterBias->dy(cluster.row, cluster.y, o2::tpc::cluster_bias::driftLength(cluster.z),
                                           residual.isPositionOnly() ? 0. : cluster.tgSlp, residual.channel);
    }
    cluster.padY = cluster.y;
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
      context.fineMap->apply(cluster.sector, cluster.row, cluster.x, cluster.y, cluster.z, cluster.fineCorrection);
    }
    cluster.yUncorrected = cluster.y;
    // TPC cluster flags (ClusterNative split pad 0x1 / split time 0x2 / edge 0x4 / single 0x8, merger shared 0x10), stored by
    // TrackInterpolation in UnbinnedResid::channel
    cluster.flagged = context.clusterFlagMask != 0 && (residual.channel & context.clusterFlagMask);
    if (context.edgeFit.pads > 0.) {
      cluster.excluded = cluster.flagged || (isExcludedRow(cluster.row) && !(context.lastRowsInFit && cluster.row >= NRows - EdgeRows));
      cluster.removed = context.geometry->padPlaneEdge[cluster.row] - std::abs(cluster.padY) < (context.edgeFit.pads + 1.) * context.geometry->padWidth[cluster.row] ||
                        context.deadPads.near(cluster.sector, cluster.row, cluster.padY, 1.);
    } else {
      cluster.excluded = cluster.flagged || isExcludedFromFit(cluster.row, context.edgeWeightByTruth && std::isfinite(cluster.yTrue) ? cluster.yTrue : cluster.padY, *context.geometry);
    }
    nFitClusters += !cluster.excluded && !cluster.removed;
    clusters.push_back(cluster);
  }
  if (nFitClusters < MinFitClusters || clusters.empty()) {
    return;
  }
  if (edgeTrack && !sweepOK && (context.truthClusters || context.edgeRemovePads > 0. || context.edgeWeightByTruth)) {
    return; // the tests need the truth at every cluster
  }
  ++counters[Selected];
  if (context.itsInFit && !addInnerPoints(task, trackRef, context.itsInFit, clusters)) {
    ++counters[NoInnerPoints];
    return;
  }
  // the seed: the ITS track's outer parameters at the first point
  TrackParD seed(itsTrack.getX(), itsTrack.getAlpha(), {itsTrack.getY(), itsTrack.getZ(), itsTrack.getSnp(), itsTrack.getTgl(), itsTrack.getQ2Pt()},
                 itsTrack.getAbsCharge(), itsTrack.getPID());
  FitResult fit;
  if (!rotateTo(seed, clusters.front().alpha) ||
      !o2::base::PropagatorD::Instance()->PropagateToXBxByBz(seed, clusters.front().x, 0.999, 2., MatCorrType::USEMatCorrLUT) ||
      !fitTrack(clusters, seed, fit, context.edgeFit.table || context.edgeFit.pads > 0. ? &context.edgeFit : nullptr)) {
    ++counters[FitFailed];
    return;
  }
  ++counters[Fitted];
  if (context.kinkQA) {
    thread_local long nKinkQA = 0;
    if (++nKinkQA % context.kinkQAEvery == 0) {
      writeKinkQA(fit, clusters, context, counters);
    }
  }
  counters[FitPasses] += fit.passes;
  counters[AtMaxPasses] += fit.passes >= MaxFitPasses;
  const int charge = fit.qOverPt > 0. ? 1 : 0;
  Eigen::VectorXd residuals(2);
  Eigen::VectorXd measurementErrors(2);
  Eigen::VectorXd residualErrors(2);
  Eigen::VectorXd downWeights(2);
  Eigen::VectorXd corrections;
  Eigen::MatrixXd covariance;
  // edge tree: per cluster what the fine pass did with it, its best position and unbiased residuals
  std::vector<EdgeRecord> edgeRecords(edgeTrack ? clusters.size() : 0);
  auto mark = [&](size_t point, int status) {
    if (edgeTrack) {
      edgeRecords[fit.active[point]].status = status;
    }
  };
  for (size_t point = 0; point < fit.active.size(); ++point) {
    const auto& cluster = clusters[fit.active[point]];
    if (cluster.type != Cluster::TPC) {
      continue;
    }
    if (fit.outlier[point]) {
      mark(point, 1);
      continue;
    }
    if (cluster.flagged) {
      ++counters[FlaggedResiduals];
      mark(point, 2);
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
      mark(point, 3);
      continue;
    }
    const double bestPadY = cluster.padY + bestPosition[0] - cluster.yUncorrected; // = - cluster.y without edgeBiasInFit
    const double predY = cluster.y - unbiased[0];                                  // the track without the cluster
    const double predPadY = cluster.padY + predY - cluster.yUncorrected;
    const double edgePadY = context.edgeByPred ? predPadY : bestPadY; // the position of the edge cut and the table
    if (edgeTrack) {
      auto& record = edgeRecords[fit.active[point]];
      record.bestPadY = bestPadY;
      record.residualY = unbiased[0];
      record.residualZ = unbiased[1];
      record.trackSigmaY = std::sqrt(std::max(0., 1. / weight[0] - cluster.sigmaY * cluster.sigmaY));
      record.predPadY = predPadY;
      record.edgeShift = cluster.edgeShift;
      // the slope of the voxel model (as below): the reference's, with fitSlopes the fitted track's
      double snpFit = fit.reference[point].getSnp();
      if (context.fitSlopes && fit.trajectory->getResults(-int(point + 1), corrections, covariance) == 0) {
        snpFit += corrections(1);
      }
      record.tgSlpFit = snpFit / std::sqrt(std::max(1e-6, (1. - snpFit) * (1. + snpFit)));
    }
    if (unbiased[0] * unbiased[0] * weight[0] > OutlierSigma * OutlierSigma || unbiased[1] * unbiased[1] * weight[1] > OutlierSigma * OutlierSigma) {
      ++counters[Outliers];
      mark(point, 4);
      continue;
    }
    if (std::abs(edgePadY) > context.geometry->padPlaneEdge[cluster.row] - context.skipEdgePads * context.geometry->padWidth[cluster.row]) {
      ++counters[EdgePadResiduals];
      mark(point, 5);
      continue;
    }
    if (context.deadPads.near(cluster.sector, cluster.row, predPadY)) {
      ++counters[DeadPadResiduals];
      mark(point, 8);
      continue;
    }
    const double y2x = (context.edgeByPred ? predY : bestPosition[0]) / cluster.x;
    const double z2x = bestPosition[1] / cluster.x;
    const int binY2X = std::abs(z2x) > context.maxZ2X ? -1 : context.binning->getY2XBinExact(y2x, cluster.row);
    const int binZ2X = context.binning->getZ2XBinExact(cluster.sector < 18 ? z2x : -z2x);
    if (edgeTrack) {
      edgeRecords[fit.active[point]].binY2X = binY2X;
    }
    if (binY2X < 0 || binY2X >= context.binning->getNY2XBins() || binZ2X < 0 || std::abs(z2x) > context.maxZ2X) {
      ++counters[OutsideGrid];
      mark(point, 6);
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
    if (context.addBackFineMap) {
      // the residual of the uncorrected cluster: (dX, dY, dZ) = fine correction + what is left (r_y = dY - t dX, r_z = dZ - c dX)
      const double* correction = cluster.fineCorrection;
      contribution.residualY += correction[1] - contribution.slopeY * correction[0];
      contribution.residualZ += correction[2] - contribution.slopeZ * correction[0];
    }
    if (context.bestEdgeBias && !context.edgeFit.table) {
      // the edge part of the cluster bias at the table's conditioning variable (best position; edgeByPred: the track without the cluster), the voxel model's slope
      const double width = context.geometry->padWidth[cluster.row];
      const double side = cluster.padY0 < 0. ? -1. : 1.;
      const double dEdge = (context.geometry->padPlaneEdge[cluster.row] - side * edgePadY) / width;
      const double padMeas = (context.geometry->padPlaneEdge[cluster.row] - std::abs(cluster.padY0)) / width - 0.5;
      const int cls = o2::tpc::cluster_bias::classOfFlags(cluster.flags, std::abs(padMeas - std::round(padMeas)) < 0.25 / 64.);
      const auto& trackRefOcc = task.trackRefs[trackIndex];
      const int occStack = o2::tpc::cluster_bias::stackOfRow(cluster.row);
      const double occupancy = context.bestEdgeBias->hasOccupancy() ? (trackRefOcc.getMultStackPacked(occStack) > 0 ? trackRefOcc.getMultStack(occStack) : 0.) : -1.;
      const double bias = context.bestEdgeBias->inwardBias(context.bestEdgeBias->unitOfRow(cluster.row), cls, o2::tpc::cluster_bias::driftLength(cluster.z0),
                                                           std::abs(contribution.slopeY), dEdge, occupancy);
      contribution.residualY += side * bias;
      counters[EdgeBiasCorrected] += bias != 0.;
      if (edgeTrack) {
        edgeRecords[fit.active[point]].edgeCorr = side * bias;
      }
    }
    contribution.weightZ = weight[1];
    batch.push_back(contribution);
    ++counters[Residuals];
    counters[ResidualsExcluded] += cluster.excluded;
    mark(point, 0);
  }
  if (edgeTrack) {
    writeEdgeTree(task, trackIndex, fit, clusters, edgeRecords, truthState, sweepOK, context, counters);
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

/// voxResTree (DS = D; with the binning, the CTP / IDC means and medians in the user info), OrbitLumiInfo and optionally fineSums;
/// extraInfo: more user info
void writeOutput(const std::string& fileOutput, const TrackResiduals& binning, const std::string& y2xBinning, const std::string& z2xBinning, float maxZ2X,
                 const std::vector<std::array<ChargeSums, 2>>& sums, const std::vector<std::pair<std::string, std::string>>& extraInfo, long nTracks,
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
  userInfo->Add(new TNamed("y2xEdge", y2xEdgeDefinition(binning).data()));
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
          voxTree.Fill();
        }
      }
    }
  }
  voxTree.Write();
  LOGP(info, "{}: {} voxels, {} solved (min entries {}), not smoothed (DS = D)", fileOutput, sums.size(), nSolved, minEntries);

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
    for (const char* key : {"y2xBinning", "z2xBinning", "maxZ2X", "nTracks", "y2xEdge"}) {
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
  if (writeSums) {
    writeSumsTree("fineSums", "per voxel and charge the sums of the unbiased TPC-only residuals", sums);
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
  std::unique_ptr<TTreeReaderValue<std::vector<o2::tpc::TrackDataMC>>> tracksMC; ///< edge tree with trkMC in the input

  void close()
  {
    tracksMC.reset();
    tracks.reset();
    trackRefs.reset();
    residuals.reset();
    trackReader.reset();
    residualReader.reset();
    file.reset();
  }

  /// externalDetectors: also what the ITS clusters and the vertex need; edgeTree: also the occupancy and the MC truth (if present)
  bool open(const std::string& fileName, bool externalDetectors, bool clusterFlags = false, bool edgeTree = false)
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
    std::vector<const char*> unused{"detInfo*"};
    if (!edgeTree) {
      unused.push_back("trackInfo.multStack*");
    }
    if (!externalDetectors) {
      unused.push_back("trackInfo.nExtDetResid");
      if (!clusterFlags) {
        unused.push_back("res.channel"); // also the TPC cluster flags (clusterFlagMask)
      }
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
    if (edgeTree && trackReader->GetTree()->GetBranch("trkMC")) {
      tracksMC = std::make_unique<TTreeReaderValue<std::vector<o2::tpc::TrackDataMC>>>(*trackReader, "trkMC");
    }
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
    if (tracksMC) {
      task.tracksMC = **tracksMC;
      if (!readOK(*tracksMC) || task.tracksMC.size() != task.tracks.size()) {
        task.tracksMC.clear(); // the edge tree then without the truth for this TF
      }
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
      fileOpen = input.open(files[tf->file], context.itsInFit > 0,
                            context.clusterFlagMask != 0 || context.clusterBias != nullptr || context.edgeTree != nullptr || context.bestEdgeBias != nullptr,
                            context.edgeTree != nullptr || (context.bestEdgeBias != nullptr && context.bestEdgeBias->hasOccupancy()));
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
       "residuals {} (of de-weighted clusters {}), unbiased pull > {} {}, edge pads {}, dead pads {}, cluster flags 0x{:x} {}, outside the grid {}, without a residual {}; "
       "without the ITS clusters / vertex for the fit {}",
       label, nTFs, nTFsReadError, counters[Selected].load(), nFitted, counters[FitFailed].load(), counters[OutOfRowOrder].load(), counters[Residuals].load(),
       counters[ResidualsExcluded].load(), OutlierSigma, counters[Outliers].load(), counters[EdgePadResiduals].load(), counters[DeadPadResiduals].load(), context.clusterFlagMask,
       counters[FlaggedResiduals].load(), counters[OutsideGrid].load(),
       counters[NoResidual].load(), counters[NoInnerPoints].load());
  if (context.kinkQA) {
    LOGP(info, "fineMapCreatorGBL{}: kink QA entries {} (skipped while busy {})", label, counters[KinkQAWritten].load(), counters[KinkQASkipped].load());
  }
  if (context.edgeTree) {
    LOGP(info, "fineMapCreatorGBL{}: edge tree: {} tracks, {} clusters; MC truth swept for {} tracks (sweep failed {})", label, counters[EdgeTreeTracks].load(),
         counters[EdgeTreeClusters].load(), counters[EdgeTreeTruth].load(), counters[EdgeTreeSweepFailed].load());
  }
  if (context.bestEdgeBias) {
    LOGP(info, "fineMapCreatorGBL{}: best-position edge bias added to {} summed residuals", label, counters[EdgeBiasCorrected].load());
  }
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
                       int itsInFit = 0,                                                             // tests: 1 the ITS clusters in the track fit, 2 also the primary vertex; 0 TPC clusters only
                       bool fitSlopes = false,                                                       // the slopes of the voxel model from the fitted track at the cluster instead of the reference
                       int kinkQAEvery = 0,                                                          // with itsInFit: the kinks before the TPC of every n-th fitted track (per thread) to <output stem>.kinkQA.root; 0 off
                       bool addBackFineMap = false,                                                  // with fineMapInput: its correction added back to the summed residuals (an iteration: the output is the full fine map)
                       int clusterFlagMask = 0,                                                      // TPC clusters with any of these flags (UnbinnedResid::channel: 0x1 split pad, 0x2 split time, 0x4 edge,
                                                                                                     // 0x8 single pad or time, 0x10 shared) are de-weighted in the fit and not summed; 0: none. "Golden": 0x19
                       std::string clusterBiasFile = "",                                             // table of the per-cluster correction of the cluster-position biases (pad-phase S-curve,
                                                                                                     // sector-edge truncation; format in ClusterBiasCorrection.h), applied to the snapped cluster
                                                                                                     // y before the edge cut and the binning; empty: off
                       float edgeTreePads = 0.f,                                                     // > 0: the edge tree (<output stem>.edgeTree.root) with the clusters within this many pads of the
                                                                                                     // pad plane's edge (measured, best or MC true position); 0: off
                       int edgeTreeEvery = 1,                                                        // edge tree: every n-th track of a TF
                       std::string bestEdgeBiasFile = "",                                            // edge part of the cluster bias vs the BEST position (E-only ClusterBiasCorrection
                                                                                                     // table, dist = d_best), added to each summed residual; empty: off
                       float skipEdgePads = 1.f,                                                     // residuals with the best position within this many pads of the pad-plane edge not summed (= fine_map::SkipEdgePads)
                       bool edgeBiasInFit = false,                                                   // with bestEdgeBiasFile: correct the clusters inside the fit by the table at the unbiased
                                                                                                     // track position (table by d_pred) instead of the summed residuals
                       bool truthClusters = false,                                                   // test (MC, edge tree): the TPC clusters at the swept truth + Gaussian noise of their error
                       float edgeRemovePads = 0.f,                                                   // test (MC, edge tree): the clusters within this many pads of the edge (truth) out of the fit
                       bool edgeWeightByTruth = false,                                               // test (MC, edge tree): the de-weighting within EdgePads decided at the truth
                       float edgeFitPads = 0.f,                                                      // > 0: instead of that de-weighting the clusters within this many pads of the edge at the track
                                                                                                     // without the cluster out of the fit
                       bool edgeByPred = false,                                                      // the skipEdgePads cut, the y/x bin and the bestEdgeBiasFile lookup by the track without the cluster
                       std::string deadPadFile = "",                                                 // dead channel list ("sector row pad"): clusters beside a dead pad out of the fit and not summed; empty: off
                       float deadPadMargin = 2.f,                                                    // ... within this many pads of a dead pad at the track without the cluster
                       bool lastRowsInFit = false)                                                   // test (no-distortion MC, with edgeFitPads): the last rows are not de-weighted in the fit
{
  using namespace fine_map;
  const auto startTime = std::chrono::steady_clock::now();
  ROOT::EnableThreadSafety();
  if (runNumber < 0) {
    LOGP(fatal, "A run number is required");
  }
  nThreads = std::max(1, nThreads);
  if (itsInFit < 0 || itsInFit > 2) {
    LOGP(fatal, "itsInFit {}: 0, 1 or 2", itsInFit);
  }
  if (addBackFineMap && fineMapInput.empty()) {
    LOGP(fatal, "addBackFineMap needs a fineMapInput");
  }
  if (edgeTreePads > 0.f && edgeTreeEvery < 1) {
    LOGP(fatal, "edgeTreeEvery {}: at least 1", edgeTreeEvery);
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
  // per-cluster correction of the cluster-position biases (empty file name: off)
  o2::tpc::cluster_bias::ClusterBiasCorrection clusterBias;
  if (!clusterBiasFile.empty()) {
    std::string error;
    if (!clusterBias.load(clusterBiasFile, &error)) {
      LOGP(fatal, "cluster-bias correction: cannot use {}: {}", clusterBiasFile, error);
    }
    const auto& mapper = o2::tpc::Mapper::instance();
    for (int row = 0; row < NRows; ++row) {
      clusterBias.setRowGeometry(row, geometry.padWidth[row], mapper.getNumberOfPadsInRowSector(row), o2::tpc::Mapper::REGION[row]);
    }
    LOGP(info, "cluster-bias correction from {}: per {}, {} drift x {} slope bins, {} phase nodes, {} distance bins", clusterBiasFile,
         clusterBias.byRegion ? "pad region" : "stack", clusterBias.nDrift(), clusterBias.nSlope(), clusterBias.nPhase, clusterBias.nDist);
  }

  // the edge part of the cluster bias vs the best position, on the summed residuals (empty file name: off)
  o2::tpc::cluster_bias::ClusterBiasCorrection bestEdgeBias;
  if (!bestEdgeBiasFile.empty()) {
    std::string error;
    if (!bestEdgeBias.load(bestEdgeBiasFile, &error)) {
      LOGP(fatal, "best-position edge bias: cannot use {}: {}", bestEdgeBiasFile, error);
    }
    const auto& mapper = o2::tpc::Mapper::instance();
    for (int row = 0; row < NRows; ++row) {
      bestEdgeBias.setRowGeometry(row, geometry.padWidth[row], mapper.getNumberOfPadsInRowSector(row), o2::tpc::Mapper::REGION[row]);
    }
    LOGP(info, "best-position edge bias from {}: per {}, {} drift x {} slope bins, {} d_best bins in [{}, {}) pads, added to the summed residuals",
         bestEdgeBiasFile, bestEdgeBias.byRegion ? "pad region" : "stack", bestEdgeBias.nDrift(), bestEdgeBias.nSlope(), bestEdgeBias.nDist,
         bestEdgeBias.distLo, bestEdgeBias.distHi);
  }
  if (skipEdgePads != SkipEdgePads) {
    LOGP(info, "Residuals with the best position within {} pads of the pad-plane edge not summed (default {})", skipEdgePads, SkipEdgePads);
  }

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
    LOGP(info, "ITS clusters{} in the track fit, {} geometry; slopes from the {}; kinks before the TPC beyond {} sigma inflated",
         itsInFit > 1 ? " and the primary vertex" : "", isMC ? "ideal (MC)" : "aligned", fitSlopes ? "fitted track" : "reference", KinkTailPull);
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
  LOGP(info, "fineMapCreatorGBL: {} files, {} TFs ({} without IDC{}); correction map {}, derivative {}, fine map {}{}; grid y/x {}, z/x {} up to {} ({} voxels per sector); {} fit threads",
       files.size(), tfs.size(), nTFsNoIDC, needIDC ? ": skipped" : "", corrMap.empty() ? "none" : corrMap, corrMapDerivative.empty() ? "none" : corrMapDerivative,
       fineMapInput.empty() ? "none" : fineMapInput, addBackFineMap ? " (added back to the residuals)" : "", y2xBinning, z2xBinning, maxZ2X, voxelsPerSector, nThreads);

  // the TFs to fit (with a derivative map those with IDC) in input order
  std::vector<const TFInfo*> fitTFs;
  for (const auto& tf : tfs) {
    if (!needIDC || tf.hasIDC) {
      fitTFs.push_back(&tf);
    }
  }
  FitContext context;
  context.sources = GID::getSourcesMask(trackSources);
  context.baseCorrection = baseCorrection.map ? &baseCorrection : nullptr;
  context.fineMap = fineMap.nodes.empty() ? nullptr : &fineMap;
  context.addBackFineMap = addBackFineMap;
  context.binning = binning.get();
  context.geometry = &geometry;
  context.maxZ2X = maxZ2X;
  context.clusterFlagMask = clusterFlagMask;
  context.clusterBias = clusterBias.isActive() ? &clusterBias : nullptr;
  context.bestEdgeBias = bestEdgeBias.isActive() ? &bestEdgeBias : nullptr;
  if (edgeBiasInFit && context.bestEdgeBias) {
    context.edgeFit.table = context.bestEdgeBias;
    context.edgeFit.geometry = &geometry;
    LOGP(info, "best-position edge bias applied to the clusters inside the track fit (table by d_pred), not to the summed residuals");
  } else if (edgeBiasInFit) {
    LOGP(fatal, "edgeBiasInFit needs a bestEdgeBiasFile");
  }
  context.skipEdgePads = skipEdgePads;
  if ((truthClusters || edgeRemovePads > 0.f || edgeWeightByTruth) && edgeTreePads <= 0.f) {
    LOGP(fatal, "truthClusters / edgeRemovePads / edgeWeightByTruth need the edge tree (the MC truth is swept only for its tracks)");
  }
  context.truthClusters = truthClusters;
  context.edgeRemovePads = edgeRemovePads;
  context.edgeWeightByTruth = edgeWeightByTruth;
  context.edgeFit.geometry = &geometry;
  context.edgeFit.pads = edgeFitPads;
  context.lastRowsInFit = lastRowsInFit;
  if (lastRowsInFit) {
    LOGP(info, "The last {} rows keep their weight in the fit", EdgeRows);
  }
  if (!deadPadFile.empty()) {
    context.deadPads.load(deadPadFile, deadPadMargin, geometry);
    context.edgeFit.deadPads = &context.deadPads;
  }
  context.edgeByPred = edgeByPred;
  if (edgeByPred) {
    LOGP(info, "Edge cut, y/x bin and edge table by the track without the cluster instead of the best position");
  }
  if (edgeWeightByTruth) {
    LOGP(info, "TEST: de-weighting within {} pads of the edge at the MC truth", EdgePads);
  }
  if (edgeFitPads > 0.f) {
    LOGP(info, "Clusters within {} pads of the edge at the track without the cluster out of the fit (instead of the de-weighting by the measured position)", edgeFitPads);
  }
  if (truthClusters) {
    LOGP(info, "TEST: TPC clusters of the edge-tree tracks at the MC truth + noise");
  }
  if (edgeRemovePads > 0.f) {
    LOGP(info, "TEST: clusters within {} pads of the edge (MC truth) out of the fit", edgeRemovePads);
  }
  context.voxelsPerSector = voxelsPerSector;
  context.itsInFit = itsInFit;
  context.fitSlopes = fitSlopes;
  // the kink QA: one streamer for the fit threads, <output stem>.kinkQA.root
  std::unique_ptr<o2::utils::TreeStreamRedirector> kinkQA;
  std::mutex kinkQAMutex;
  if (kinkQAEvery > 0 && itsInFit > 0) {
    const std::filesystem::path outputPath(fileOutput);
    const std::string kinkQAFile = ((outputPath.parent_path().empty() ? std::filesystem::path(".") : outputPath.parent_path()) / (outputPath.stem().string() + ".kinkQA.root")).string();
    kinkQA = std::make_unique<o2::utils::TreeStreamRedirector>(kinkQAFile.c_str(), "recreate");
    context.kinkQA = kinkQA.get();
    context.kinkQAMutex = &kinkQAMutex;
    context.kinkQAEvery = kinkQAEvery;
    LOGP(info, "Kink QA of every {}th fitted track per fit thread: {}", kinkQAEvery, kinkQAFile);
  }
  // the edge tree: one streamer for the fit threads, <output stem>.edgeTree.root
  std::unique_ptr<o2::utils::TreeStreamRedirector> edgeTree;
  std::mutex edgeTreeMutex;
  if (edgeTreePads > 0.f) {
    const std::filesystem::path outputPath(fileOutput);
    const std::string edgeTreeFile = ((outputPath.parent_path().empty() ? std::filesystem::path(".") : outputPath.parent_path()) / (outputPath.stem().string() + ".edgeTree.root")).string();
    edgeTree = std::make_unique<o2::utils::TreeStreamRedirector>(edgeTreeFile.c_str(), "recreate");
    context.edgeTree = edgeTree.get();
    context.edgeTreeMutex = &edgeTreeMutex;
    context.edgeTreePads = edgeTreePads;
    context.edgeTreeEvery = edgeTreeEvery;
    LOGP(info, "Edge tree: clusters within {} pads of the pad plane's edge of every {}th track of a TF{}: {}", edgeTreePads, edgeTreeEvery,
         isMC ? ", with the MC truth" : "", edgeTreeFile);
  }
  std::vector<std::array<ChargeSums, 2>> sums(size_t(NSectors) * voxelsPerSector);
  SelectedTFs selected;
  const long nFitted = runPass(fitTFs, files, context, nThreads, sums, selected, "");
  std::vector<std::pair<std::string, std::string>> extraInfo;
  extraInfo.emplace_back("smoothing", "none");
  extraInfo.emplace_back("fineMapInput", fineMapInput.empty() ? std::string("none") : fineMapInput + (addBackFineMap ? " (added back)" : " (subtracted)"));
  if (kinkQA) {
    kinkQA->Close();
  }
  if (edgeTree) {
    edgeTree->Close();
  }
  writeOutput(fileOutput, *binning, y2xBinning, z2xBinning, maxZ2X, sums, extraInfo, nFitted, minEntries, writeSums, runNumber,
              orbitResetMS, isPbPb, firstTFTime, lastTFTime, selected);
  LOGP(info, "fineMapCreatorGBL: {:.1f} s", std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count());
}
