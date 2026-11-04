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

/// \file ClusterBiasSmoothFit.h
/// \brief Numerics of the SMOOTH ClusterBiasCorrection table (makeClusterBiasTable.C, smooth = true): clamped cubic B-splines,
/// second-difference penalties and a penalised weighted least-squares solve (small dense system), for
/// - the S-curve: a tensor product of B-splines in the phase |u| on (0, 0.5] (the last basis dropped, so S(0.5) = 0) and in a drift
///   coordinate (sqrt(drift)), fitted to the odd part of the data profiles;
/// - the edge part: B-splines in d (the last basis dropped, E = 0 at the last knot) times linear interpolation weights between the
///   drift-bin centres of the simulation, or (single-pad clusters) the interpolation weights alone.
/// Header only, no ROOT; the knots and the penalty strength are chosen by the caller.

#ifndef O2_TPC_CLUSTERBIASSMOOTHFIT_H
#define O2_TPC_CLUSTERBIASSMOOTHFIT_H

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace o2::tpc::cluster_bias::smooth
{

using Matrix = std::vector<std::vector<double>>;
using SparseRow = std::vector<std::pair<int, double>>; // (column, value) of the non-zero entries of a design-matrix row

/// clamped B-spline basis of the given degree on the knots (boundary knots repeated degree times), evaluated at x clipped to
/// [knots.front(), knots.back()); size knots.size() + degree - 1 (Cox-de Boor recursion)
inline std::vector<double> bsplineBasis(double x, const std::vector<double>& knots, int degree = 3)
{
  std::vector<double> t;
  t.insert(t.end(), degree, knots.front());
  t.insert(t.end(), knots.begin(), knots.end());
  t.insert(t.end(), degree, knots.back());
  const int nBasis = int(t.size()) - degree - 1;
  x = std::clamp(x, knots.front(), knots.back() - 1e-12);
  std::vector<double> b(t.size() - 1, 0.);
  for (size_t i = 0; i + 1 < t.size(); ++i) {
    b[i] = (t[i] <= x && x < t[i + 1]) ? 1. : 0.;
  }
  for (int k = 1; k <= degree; ++k) {
    std::vector<double> bn(t.size() - 1 - k, 0.);
    for (size_t i = 0; i < bn.size(); ++i) {
      const double a = t[i + k] > t[i] ? (x - t[i]) / (t[i + k] - t[i]) : 0.;
      const double c = t[i + k + 1] > t[i + 1] ? (t[i + k + 1] - x) / (t[i + k + 1] - t[i + 1]) : 0.;
      bn[i] = a * b[i] + c * b[i + 1];
    }
    b = bn;
  }
  b.resize(nBasis);
  return b;
}

/// D^T D of the second-difference matrix D ((n - 2) x n)
inline Matrix secondDifferencePenalty(int n)
{
  Matrix p(n, std::vector<double>(n, 0.));
  const double d[3] = {1., -2., 1.};
  for (int r = 0; r + 2 < n; ++r) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        p[r + i][r + j] += d[i] * d[j];
      }
    }
  }
  return p;
}

/// kron(a, I_m) (aIsFirst) or kron(I_m, a) of a square matrix a
inline Matrix kronIdentity(const Matrix& a, int m, bool aIsFirst)
{
  const int n = int(a.size());
  Matrix k(n * m, std::vector<double>(n * m, 0.));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int l = 0; l < m; ++l) {
        if (aIsFirst) {
          k[i * m + l][j * m + l] = a[i][j];
        } else {
          k[l * n + i][l * n + j] = a[i][j];
        }
      }
    }
  }
  return k;
}

/// solution of a x = b by Gaussian elimination with partial pivoting; empty if singular
inline std::vector<double> solve(Matrix a, std::vector<double> b)
{
  const int n = int(b.size());
  for (int c = 0; c < n; ++c) {
    int piv = c;
    for (int r = c + 1; r < n; ++r) {
      if (std::abs(a[r][c]) > std::abs(a[piv][c])) {
        piv = r;
      }
    }
    if (a[piv][c] == 0.) {
      return {};
    }
    std::swap(a[c], a[piv]);
    std::swap(b[c], b[piv]);
    for (int r = c + 1; r < n; ++r) {
      const double f = a[r][c] / a[c][c];
      if (f == 0.) {
        continue;
      }
      for (int k = c; k < n; ++k) {
        a[r][k] -= f * a[c][k];
      }
      b[r] -= f * b[c];
    }
  }
  std::vector<double> x(n, 0.);
  for (int r = n - 1; r >= 0; --r) {
    double s = b[r];
    for (int k = r + 1; k < n; ++k) {
      s -= a[r][k] * x[k];
    }
    x[r] = s / a[r][r];
  }
  return x;
}

/// penalised weighted least squares: minimise sum ((y - A c) / sig)^2 + lam c^T P c (+ 1e-9 |c|^2 against empty columns);
/// rows of A sparse; P empty = no penalty
inline std::vector<double> penalisedLeastSquares(const std::vector<SparseRow>& rows, const std::vector<double>& y, const std::vector<double>& sig,
                                                 const Matrix& p, double lam, int nPar)
{
  Matrix m(nPar, std::vector<double>(nPar, 0.));
  std::vector<double> rhs(nPar, 0.);
  for (size_t n = 0; n < rows.size(); ++n) {
    const double w = 1. / (sig[n] * sig[n]);
    for (const auto& [i, ai] : rows[n]) {
      rhs[i] += w * ai * y[n];
      for (const auto& [j, aj] : rows[n]) {
        m[i][j] += w * ai * aj;
      }
    }
  }
  for (int i = 0; i < nPar; ++i) {
    for (int j = 0; j < nPar && !p.empty(); ++j) {
      m[i][j] += lam * p[i][j];
    }
    m[i][i] += 1e-9;
  }
  return solve(m, rhs);
}

/// linear interpolation weights of x between the (sorted) centres, constant outside
inline std::vector<double> driftHats(double x, const std::vector<double>& centres)
{
  const int n = int(centres.size());
  std::vector<double> w(n, 0.);
  if (n == 1) {
    w[0] = 1.;
    return w;
  }
  x = std::clamp(x, centres.front(), centres.back());
  int i = int(std::upper_bound(centres.begin(), centres.end(), x) - centres.begin()) - 1;
  i = std::min(i, n - 2);
  const double f = (x - centres[i]) / (centres[i + 1] - centres[i]);
  w[i] = 1. - f;
  w[i + 1] = f;
  return w;
}

/// non-zero entries of a basis vector
inline SparseRow nonZero(const std::vector<double>& v, int offset = 0)
{
  SparseRow r;
  for (size_t i = 0; i < v.size(); ++i) {
    if (v[i] != 0.) {
      r.emplace_back(offset + int(i), v[i]);
    }
  }
  return r;
}

/// S-curve surface: coef[i * nx + j] of B_i(|u|) (last u basis dropped) x B_j(x)
struct SCurveFit {
  std::vector<double> coef;
  int nu = 0;
  int nx = 0;
  double chi2 = 0.;
  int nPoints = 0;
  bool valid() const { return !coef.empty(); }
};

struct SCurvePoint {
  double u;   // |phase| in (0, 0.5)
  double x;   // drift coordinate (sqrt(drift))
  double y;   // odd part [um]
  double err; // [um]
};

inline SCurveFit fitSCurve(const std::vector<SCurvePoint>& pts, const std::vector<double>& uKnots, const std::vector<double>& xKnots, double lam)
{
  SCurveFit f;
  f.nu = int(uKnots.size()) + 1; // knots + degree - 1, minus the dropped last one
  f.nx = int(xKnots.size()) + 2;
  const int nPar = f.nu * f.nx;
  std::vector<SparseRow> rows;
  std::vector<double> y;
  std::vector<double> e;
  for (const auto& p : pts) {
    auto bu = bsplineBasis(p.u, uKnots);
    bu.pop_back(); // S(0.5) = 0
    const auto bx = bsplineBasis(p.x, xKnots);
    SparseRow r;
    for (int i = 0; i < f.nu; ++i) {
      for (int j = 0; j < f.nx && bu[i] != 0.; ++j) {
        if (bx[j] != 0.) {
          r.emplace_back(i * f.nx + j, bu[i] * bx[j]);
        }
      }
    }
    rows.push_back(r);
    y.push_back(p.y);
    e.push_back(p.err);
  }
  auto pen = kronIdentity(secondDifferencePenalty(f.nu), f.nx, true);
  const auto px = kronIdentity(secondDifferencePenalty(f.nx), f.nu, false);
  for (int i = 0; i < nPar; ++i) {
    for (int j = 0; j < nPar; ++j) {
      pen[i][j] += px[i][j];
    }
  }
  f.coef = penalisedLeastSquares(rows, y, e, pen, lam, nPar);
  if (f.valid()) {
    for (size_t n = 0; n < rows.size(); ++n) {
      double v = 0.;
      for (const auto& [i, a] : rows[n]) {
        v += a * f.coef[i];
      }
      f.chi2 += (y[n] - v) * (y[n] - v) / (e[n] * e[n]);
    }
    f.nPoints = int(rows.size());
  }
  return f;
}

/// S at |u| (> 0) and the drift coordinate x
inline double evalSCurve(const SCurveFit& f, const std::vector<double>& uKnots, const std::vector<double>& xKnots, double absU, double x)
{
  auto bu = bsplineBasis(absU, uKnots);
  bu.pop_back();
  const auto bx = bsplineBasis(x, xKnots);
  double v = 0.;
  for (int i = 0; i < f.nu; ++i) {
    for (int j = 0; j < f.nx; ++j) {
      v += bu[i] * bx[j] * f.coef[i * f.nx + j];
    }
  }
  return v;
}

/// edge part: coef[i * nc + c] of B_i(d) (last d basis dropped) x hat_c(drift), or (single) coef[c] of hat_c(drift)
struct EdgeFit {
  bool single = false;
  std::vector<double> coef;
  std::vector<double> centres;
  int nd = 0;
  double chi2 = 0.;
  int nPoints = 0;
  bool valid() const { return !coef.empty(); }
};

struct EdgePoint {
  double d;     // measured d [pads]
  double drift; // drift-bin centre of the simulation [cm]
  double y;     // edge part [um]
  double err;   // [um]
};

inline EdgeFit fitEdge(const std::vector<EdgePoint>& pts, const std::vector<double>& centres, const std::vector<double>& dKnots, double lam, bool single)
{
  EdgeFit f;
  f.single = single;
  f.centres = centres;
  const int nc = int(centres.size());
  f.nd = single ? 0 : int(dKnots.size()) + 1;
  const int nPar = single ? nc : f.nd * nc;
  std::vector<SparseRow> rows;
  std::vector<double> y;
  std::vector<double> e;
  for (const auto& p : pts) {
    const auto h = driftHats(p.drift, centres);
    SparseRow r;
    if (single) {
      r = nonZero(h);
    } else {
      auto bd = bsplineBasis(p.d, dKnots);
      bd.pop_back(); // E = 0 at the last knot
      for (int i = 0; i < f.nd; ++i) {
        for (int c = 0; c < nc && bd[i] != 0.; ++c) {
          if (h[c] != 0.) {
            r.emplace_back(i * nc + c, bd[i] * h[c]);
          }
        }
      }
    }
    rows.push_back(r);
    y.push_back(p.y);
    e.push_back(p.err);
  }
  const Matrix pen = single ? Matrix{} : kronIdentity(secondDifferencePenalty(f.nd), nc, true);
  f.coef = penalisedLeastSquares(rows, y, e, pen, single ? 0. : lam, nPar);
  if (f.valid()) {
    for (size_t n = 0; n < rows.size(); ++n) {
      double v = 0.;
      for (const auto& [i, a] : rows[n]) {
        v += a * f.coef[i];
      }
      f.chi2 += (y[n] - v) * (y[n] - v) / (e[n] * e[n]);
    }
    f.nPoints = int(rows.size());
  }
  return f;
}

/// edge part at d and drift: constant below dKnots.front(), 0 from dKnots.back() on (single: the fitted value below singleMax)
inline double evalEdge(const EdgeFit& f, const std::vector<double>& dKnots, double d, double drift, double singleMax = 1.)
{
  const auto h = driftHats(drift, f.centres);
  const int nc = int(f.centres.size());
  if (f.single) {
    if (d >= singleMax) {
      return 0.;
    }
    double v = 0.;
    for (int c = 0; c < nc; ++c) {
      v += f.coef[c] * h[c];
    }
    return v;
  }
  if (d >= dKnots.back()) {
    return 0.;
  }
  auto bd = bsplineBasis(std::clamp(d, dKnots.front(), dKnots.back() - 1e-9), dKnots);
  bd.pop_back();
  double v = 0.;
  for (int i = 0; i < f.nd; ++i) {
    for (int c = 0; c < nc; ++c) {
      v += bd[i] * h[c] * f.coef[i * nc + c];
    }
  }
  return v;
}

} // namespace o2::tpc::cluster_bias::smooth

#endif
