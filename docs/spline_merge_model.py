# Copyright 2019-2026 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

"""1D model of O2's compact (Hermite) spline as used by TPCFastSpaceChargeCorrection, to estimate
(A) the error of TPCFastSpaceChargeCorrectionHelper::addCorrections when the added map has other knots
    or another range than the main map, and
(B) the edge error of TPCFastSpaceChargeCorrectionHelper::initInverse (real grid := measured grid, clamping).
The formulas follow GPU/Utils/Spline1DSpec.h (getSderivativesOverParsAtU) and Spline2DHelper::approximateDataPoints."""
import numpy as np


class Spline1D:
    """compact cubic Hermite spline on integer knots u_k in [0, umax]; params (f_k, d_k = df/du)"""

    def __init__(self, knots_u, x0, x1, fade=None):
        self.u = np.asarray(knots_u, float)
        self.umax = self.u[-1]
        self.x0, self.x1 = x0, x1
        self.scale = self.umax / (x1 - x0)
        self.f = np.zeros(len(self.u))
        self.d = np.zeros(len(self.u))
        self.fade = fade  # (zOut, zReadout) or None

    def to_u(self, x):
        return np.clip((np.asarray(x, float) - self.x0) * self.scale, 0.0, self.umax)

    def basis(self, u):
        i = np.clip(np.searchsorted(self.u, u, side="right") - 1, 0, len(self.u) - 2)
        h = self.u[i + 1] - self.u[i]
        v = (u - self.u[i]) / h
        sr = v * v * (3 - 2 * v)
        return i, 1 - sr, sr, h * v * (1 - v) ** 2, h * v * v * (v - 1)

    def spline_scale(self, x):
        if self.fade is None:
            return np.ones_like(np.asarray(x, float))
        z_out, z_ro = self.fade
        k = 1.0 / (z_ro - z_out) if abs(z_ro - z_out) > 1 else 0.0
        return 1 - np.clip((np.asarray(x, float) - z_out) * k, 0, 1)

    def __call__(self, x):
        u = self.to_u(x)
        i, a, b, c, e = self.basis(u)
        return self.spline_scale(x) * (a * self.f[i] + b * self.f[i + 1] + c * self.d[i] + e * self.d[i + 1])

    def value_and_du(self, x):
        """value and d/du at x (u clamped), without the z fade (as interpolateParametersAtU)"""
        u = self.to_u(x)
        i, *_ = self.basis(u)
        h = self.u[i + 1] - self.u[i]
        v = (u - self.u[i]) / h
        a, b, c, e = 1 - v * v * (3 - 2 * v), v * v * (3 - 2 * v), h * v * (1 - v) ** 2, h * v * v * (v - 1)
        val = a * self.f[i] + b * self.f[i + 1] + c * self.d[i] + e * self.d[i + 1]
        da, db = -6 * v * (1 - v) / h, 6 * v * (1 - v) / h
        dc, de = (1 - v) * (1 - 3 * v), v * (3 * v - 2)
        dval = da * self.f[i] + db * self.f[i + 1] + dc * self.d[i] + de * self.d[i + 1]
        return val, dval

    def fit(self, x, y, w=None):
        """least squares like Spline2DHelper::approximateDataPoints (tiny smoothness ignored: 1e-8)"""
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        w = np.ones_like(x) if w is None else w
        n = len(self.u)
        A = np.zeros((2 * n, 2 * n))
        B = np.zeros(2 * n)
        u = (x - self.x0) * self.scale  # data already inside or clamped by the caller
        u = np.clip(u, 0, self.umax)
        i, a, b, c, e = self.basis(u)
        for k in range(len(x)):
            idx = [2 * i[k], 2 * i[k] + 2, 2 * i[k] + 1, 2 * i[k] + 3]
            co = [a[k], b[k], c[k], e[k]]
            for p in range(4):
                B[idx[p]] += w[k] * co[p] * y[k]
                for q in range(4):
                    A[idx[p], idx[q]] += w[k] * co[p] * co[q]
        A += 1e-8 * np.eye(2 * n)
        sol = np.linalg.solve(A, B)
        self.f, self.d = sol[0::2], sol[1::2]
        return self


def residual_knots(n_bins, centre, row_x):
    """createFromTrackResiduals: a knot at every second bin centre, boundary bins included; integer grid with pitch 0.1*min spacing"""
    idx = []
    i, j = 0, n_bins - 1
    while i <= j:
        idx.append(i)
        if j >= i + 1:
            idx.append(j)
        i += 2
        j -= 2
    pos = np.sort([centre(k) for k in idx])
    pitch = 0.1 * np.min(np.diff(pos))
    ku = np.rint((pos - pos[0]) / pitch).astype(int)
    return ku, pos[0] * row_x, pos[-1] * row_x


def merged(main_like, extra, x):
    """addCorrections for one extra map onto the main knots: values and chain-rule derivatives at the main knots,
    evaluated with the main spline (clamping and z fade of the main map)"""
    m = Spline1D(main_like.u, main_like.x0, main_like.x1, main_like.fade)
    knot_x = m.x0 + m.u / m.scale
    val, du_e = extra.value_and_du(knot_x)
    ls = extra.spline_scale(knot_x)  # extra map's own z fade at the knot
    m.f = ls * val
    m.d = ls * du_e * (extra.scale / m.scale)
    return m(x)


def report(name, x, diff, masks):
    out = [f"{name}:"]
    for label, mk in masks:
        if mk.any():
            out.append(f"  {label:42s} max |diff| = {1e4 * np.max(np.abs(diff[mk])):8.1f} um")
    print("\n".join(out))


L = 250.0
print("=== (A1) M-shape-like map (4 uniform knots over z 0..250, no fade) merged onto a residual-map z grid ===")
for row_x, nz, maxz2x in [(85.0, 28, 1.4), (85.0, 5, 1.0), (200.0, 28, 1.4)]:
    ku, z_lo, z_hi = residual_knots(nz, lambda k: (k + 0.5) * maxz2x / nz, row_x)
    z_out = min(z_hi, L)
    main = Spline1D(ku, z_lo, z_hi, fade=(z_out, L))
    zs = np.linspace(0, L, 1001)
    mshape = Spline1D([0, 1, 2, 3], 0.0, L)
    mshape.fit(zs, 0.1 * np.sin(np.pi * zs / L) ** 2 * np.cos(np.pi * zs / L * 1.5))  # ~1 mm, M-like profile
    z = np.linspace(0.01, L - 0.01, 20001)
    diff = merged(main, mshape, z) - mshape(z)
    report(f"row x = {row_x:.0f} cm, z/x binning {nz} bins up to {maxz2x} (main grid z {z_lo:.1f}..{z_hi:.1f} cm)", z, diff,
           [("inside the main grid", (z >= z_lo) & (z <= min(z_hi, L))),
            ("z below the main grid (towards the CE)", z < z_lo),
            ("z beyond the main grid (fade to readout)", z > z_out)])

print("\n=== (A2) sector-edge-like map (finer grid over the full pad row) merged onto a residual-map y grid ===")
for row_x, ny, n_pads, pad_w in [(85.0, 20, 66, 0.4166), (85.0, 15, 66, 0.4166), (245.0, 20, 138, 0.6)]:
    max_y2x = np.tan(np.pi / 18) - 1.5 / row_x
    ku, y_lo, y_hi = residual_knots(ny, lambda k: -max_y2x + (k + 0.5) * 2 * max_y2x / ny, row_x)
    main = Spline1D(ku, y_lo, y_hi)
    edge = n_pads * pad_w / 2
    ys = np.linspace(-edge, edge, 2001)
    sec = Spline1D(np.arange(41), -edge, edge)
    lam = 1.0  # cm, assumed decay length of the sector-edge distortion
    sec.fit(ys, 0.1 * (np.exp(-(edge - ys) / lam) - np.exp(-(edge + ys) / lam)))  # 1 mm at the edges
    y = np.linspace(-edge, edge, 20001)
    diff = merged(main, sec, y) - sec(y)
    report(f"row x = {row_x:.0f} cm, y/x binning {ny} bins (main grid |y| <= {y_hi:.2f} cm, pad edge {edge:.2f} cm)", y, diff,
           [("|y| <= main grid edge - 3 cm", np.abs(y) <= y_hi - 3),
            ("main grid edge - 3 cm < |y| <= main grid edge", (np.abs(y) > y_hi - 3) & (np.abs(y) <= y_hi)),
            ("beyond the main grid (up to the pad edge)", np.abs(y) > y_hi)])

print("\n=== (B) initInverse edge effect: real grid := measured grid, real points outside clamped onto the edge ===")
row_x, ny = 85.0, 20
max_y2x = np.tan(np.pi / 18) - 1.5 / row_x
ku, y_lo, y_hi = residual_knots(ny, lambda k: -max_y2x + (k + 0.5) * 2 * max_y2x / ny, row_x)
for shift, edge_amp, lam in [(0.2, 0.0, 1.0), (0.5, 0.0, 1.0), (0.2, 0.3, 2.0), (0.0, 0.5, 2.0), (-0.5, 0.0, 1.0)]:
    def d(m):  # direct correction (measured -> real = m + d): a shift, a 2 % gradient and an edge bump at +y
        return shift + 0.02 * m + edge_amp * np.exp((m - y_hi) / lam)
    direct = Spline1D(ku, y_lo, y_hi)
    m_fit = np.linspace(y_lo, y_hi, 4001)
    direct.fit(m_fit, d(m_fit))
    # initInverse: sample the measured grid (10 points per knot interval), real = m + d(m), clamp to the grid, fit d
    us = np.concatenate([np.linspace(direct.u[k], direct.u[k + 1], 11)[:-1] for k in range(len(direct.u) - 1)] + [[direct.umax]])
    m_s = direct.x0 + us / direct.scale
    r_s = m_s + direct(m_s)
    inverse = Spline1D(ku, y_lo, y_hi).fit(np.clip(r_s, y_lo, y_hi), direct(m_s))
    # truth: for real r, the measured m solves m + d(m) = r; inverse value = d(m)
    m_true = np.linspace(y_lo - 2, y_hi + 2, 40001)
    r_true = m_true + direct(m_true)
    sel = (r_true >= y_lo) & (r_true <= y_hi)
    err = inverse(r_true[sel]) - direct(m_true[sel])
    r = r_true[sel]
    knot_pitch = 2 * 2 * max_y2x / ny * row_x
    print(f"  d = {shift:+.1f} cm + 2% gradient + {edge_amp:.1f} cm edge bump (lambda {lam:.0f} cm): "
          f"max |err| interior {1e4 * np.max(np.abs(err[np.abs(r) < y_hi - knot_pitch])):7.1f} um, "
          f"last knot interval {1e4 * np.max(np.abs(err[np.abs(r) >= y_hi - knot_pitch])):7.1f} um")

print("\n=== (B2) the same, split by what TPCFastSpaceChargeCorrection::testInverse checks (measured and real inside the grid) ===")
for shift, edge_amp, lam in [(0.5, 0.0, 1.0), (0.2, 0.3, 2.0), (0.0, 0.5, 2.0), (-0.5, 0.0, 1.0)]:
    def d(m):
        return shift + 0.02 * m + edge_amp * np.exp((m - y_hi) / lam)
    direct = Spline1D(ku, y_lo, y_hi)
    m_fit = np.linspace(y_lo, y_hi, 4001)
    direct.fit(m_fit, d(m_fit))
    us = np.concatenate([np.linspace(direct.u[k], direct.u[k + 1], 11)[:-1] for k in range(len(direct.u) - 1)] + [[direct.umax]])
    m_s = direct.x0 + us / direct.scale
    inverse = Spline1D(ku, y_lo, y_hi).fit(np.clip(m_s + direct(m_s), y_lo, y_hi), direct(m_s))
    m_true = np.linspace(y_lo - 2, y_hi + 2, 40001)
    r_true = m_true + direct(m_true)
    real_in = (r_true >= y_lo) & (r_true <= y_hi)
    meas_in = (m_true >= y_lo) & (m_true <= y_hi)
    def max_err(sel):
        return 1e4 * np.max(np.abs(inverse(r_true[sel]) - direct(m_true[sel]))) if sel.any() else float("nan")
    print(f"  shift {shift:+.1f} cm, edge bump {edge_amp:.1f} cm: max |err| seen by testInverse {max_err(real_in & meas_in):6.1f} um, "
          f"real inside but measured outside {max_err(real_in & ~meas_in):6.1f} um")
