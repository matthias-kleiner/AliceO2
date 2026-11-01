# Review: TPC correction maps per TF (`TPCScalerSpec` → `TPCFastTransformPOD` → reconstruction)

This reviews upstream `dev` at `1f4567d` (2026-10-31). I read the code only; nothing was built or run. The review starts from
`Detectors/TPC/workflow/src/TPCScalerSpec.cxx` and follows the map through:

- `CorrectionMapsLoader` / `CorrectionMapsHelper`, which pick the maps and compute the lumi scale;
- `TPCFastSpaceChargeCorrectionHelper::addCorrections`, which merges the maps;
- `TPCFastTransformPOD`, with its creation, correction and inverse-correction functions and the zero-copy spline path;
- the consumers: GPU TPC tracking and the merger, ITS-TPC matching, TOF, TRD, SVertexer, the TPC refitter and alignment.

**Legend.**

- **certain**: follows from the code as written; I read every code path involved.
- **plausible**: depends on configuration, CCDB content or intent.

Severity is given as high / medium / low. Findings are ranked within each category.

## How the map is built and used (for orientation)

1. **Per TF, in `TPCScalerSpec::run`:**
   - the IDC is the mean of the A and C scalers at the TF time (l. 152-171);
   - the loader updates the lumi scale (l. 173);
   - `buildMap` (l. 185-256) clones the main map (or, in `StaticMapOnly`, the reference map) and merges in the lumi-scaled reference
     or derivative map, the M-shape map and the sector-edge map with `addCorrections`. It then converts the result into a
     `TPCFastTransformPOD` with the VDrift and T0 written in, and publishes it as `TPC/TPCCORRMAP`.
2. **`addCorrections`** (`TPCFastSpaceChargeCorrectionHelper.cxx:1033-1180`) adds spline knot parameters. The values and first and
   mixed derivatives of each extra map are evaluated at the main map's knots. This is done for the direct map and, separately, for
   the two inverse maps.
3. **Consumers** read the POD straight from the message (`TPCFastTransformPOD::get(raw)`). The GPU workflow copies it every TF
   (`GPUWorkflowTPC.cxx:345-355`). They use `Transform` / `TransformXYZ` (cluster correction), and `InverseTransformYZtoX` /
   `InverseTransformYZtoNominalYZ` (seeding, attaching clusters, edge rejection, dE/dx dead-pad checks).

---

## Summary

### Certain bugs

| # | Where | Bug | Severity |
|---|-------|-----|----------|
| C1 | `TPCScalerSpec::buildMap` + `CorrectionMapsHelper::updateLumiScale` | With `lumi-type 0` ("no scaling, use map as is") and `corrmap-lumi-mode 0` (Linear) — **both are the option defaults** — the lumi scale is 0, and `buildMap` builds `0·(mean − ref) + ref` = **the reference map only**. The mean map is dropped. The legacy `TPCFastTransform::TransformLocal` treats scale 0 as "no scaling" and applies the mean map. The same happens in Linear mode whenever `meanLumi` is 0, e.g. a map without IDC / lumi info. | high (if this configuration is used) |
| C2 | `CorrectionMapsLoader::extractCCDBInputs` | The IDC → CTP fallback does not check that CTP lumi is available. With `--disable-ctp-lumi-request` and an invalid scaler, it switches the type to CTP and the mean lumis to CTP values, but the instantaneous lumi is never updated: it keeps the last IDC value. The lumi scale then mixes IDC and CTP units for the whole fallback period. | medium (only in that configuration) |
| C3 | `TPCFastTransformPOD::InverseTransformXYZtoNominalXYZ` (also legacy `TPCFastTransform`) | The interpolation weights read the **output** parameter `nx` before it is assigned. The result depends on whatever the caller passed in. No caller found. | low (dead code) |
| C4 | `TPCFastTransformPOD::TransformInTimeFrame(sector, time, z, maxTimeBin)` (also legacy) | Off by ±TPC length compared with `convTimeToZinTimeFrame` and the pad/row overload: it gives z = L instead of 0 at `time = T0 + maxTimeBin`. No caller found. | low (dead code) |
| C5 | `TPCFastTransformPOD::create` | `mApplyCorrection` is always set to `true`. `setApplyCorrectionOff()` on the source (the `NoCorrection` mode, `TPCScalerSpec.cxx:200`) is lost. Numerically harmless, because the dummy splines are zero, but `isCorrectionApplied()` reports the wrong value and every cluster still evaluates three splines. | low |
| C6 | `CorrectionMapsLoader::extractCCDBInputs` | `dumRep` is a local variable reset on every call, so the "warning n of 5" limiter never limits. The alarm fires every TF, always as "1 of 5". | low |
| C7 | `TPCMShapeCorrection::getBoundaryPotential` | Edge handling: a TF at or before the first stored time, or after the last one, gets no M-shape even when it is within `mMaxDeltaTimeMS` of it. `idx == 0` and `idx == size` are never accepted. | low |

### Plausible issues

- **P-C1:** once computed, the M-shape map stays applied indefinitely, even for TFs with no stored M-shape nearby.
- **P-C2:** the sector-edge correction is on by default, including for MC, although the code comment says "only for data".
- **P-C3:** a scaler that is invalid on one side only is not detected: (A + (−1))/2 is accepted as the IDC.
- **P-C4:** linear combination of inverse maps (an approximation, unchanged from the legacy code), and edge effects of the derivative
  merge.

### Follow-up (section 4)

- **Merging maps with different spline grids.** The merge pins every extra map to the main map's grid. Outside that grid the extra
  map's correction is frozen at the edge value, and in z it is faded out with the main map's fade. For the M-shape map this is certain;
  for the sector-edge map it depends on that map's grid.
- **`initInverse`.** It resets the inverse's (real-coordinate) grid to the measured grid and clamps points onto the edges. The error
  concentrates in the last knot interval of each row.
- **How to check.** The macro `docs/checkCorrMapMerge.C` runs these checks on the real CCDB maps. The numbers in section 4 come from a 1D
  model (`docs/spline_merge_model.py`).

### Performance

- **PERF1:** with M-shape enabled, the Poisson solve and the spline fit run **every TF**.
- **PERF2:** the full map is rebuilt every TF, with two passes in Linear mode where one is enough.
- **PERF3:** the GPU workflow copies and re-uploads the map every TF.

### Checked and found correct

- **Formulas.** The Linear (`(mean − ref)·s + ref`) and derivative (`mean + s·deriv`) combinations in `buildMap` match the legacy
  runtime formulas in `TPCFastTransform::TransformLocal`, `InverseTransformYZtoX` and `InverseTransformYZtoNominalYZ`, for the direct
  and the inverse maps.
- **`addCorrections` knot-parameter order.** It is {S, ∂S/∂v, ∂S/∂u, ∂²S/∂u∂v}. This is consistent between
  `interpolateParametersAtU` (`Spline2DSpec.h:399-403`), the knot layout read by `interpolateAtU` / `interpolateAtUZeroCopy` (a[i]
  coefficients), and the chain-rule factors `{ls, ls·scaleV, ls·scaleU, ls·scaleU·scaleV}`. The scale factors are the correct ratios
  `yScale_corr / yScale_main` and `zScale_corr / zScale_main`. The inverse maps share the direct map's knots, so `knotIndex` is right
  for all three.
- **POD layout.** `estimateSize` and `create` follow the same offset and alignment sequence, and `create` checks the result.
- **Zero-copy reads.** `interpolateAtUZeroCopy` reproduces the grid flat-buffer layout of `Spline2DContainerBase::setActualBufferAddress`
  and `Spline1DContainer::setActualBufferAddress`: knots, then the u → knot map, then grid 2 at
  `alignSize(g1size, g2align)`. It only reads plain-int members of the slim spline objects. The slim splines' internal pointers stay
  null, and no consumer uses `getSplineForRow` or `getKnot` on the POD.
- **Reinterpreting the 3D spline as 1D / 2D** for the inverses is safe. The spec-1 `getNumberOfParameters()` uses the compile-time
  `YdimT`, and the grids are shared.
- **The POD functions match the legacy `TPCFastTransform` versions one for one.** This includes `convZOffsetToVertexTime`,
  `convVertexTimeToZOffset`, the in-time-frame conversions and `getMaxDriftTime`. So no behaviour changed in the POD switch itself,
  apart from C5.
- **Consumers' inverse calls.** `GPUTPCTrackletConstructor.cxx:291` passes the same variables as input and output, which is safe
  because the inputs are taken by value. Its `do { … } while (false)` runs once, so the inverse is not applied twice.

---

## 1. Correctness (ranked)

### C1 [certain; high if used] "No scaling" plus Linear mode builds the reference map, not the mean map

**Where:** `TPCScalerSpec.cxx:207-218`; `CorrectionMapsHelper.cxx:49-50`; `CorrectionMapsLoader.cxx:165-176`; defaults in
`CorrectionMapsOptions.cxx` (`lumi-type 0`, `corrmap-lumi-mode 0`); meaning of `NoScaling` in `CorrectionMapsTypes.h:23`.

- **Why the scale is 0.** With `LumiScaleType::NoScaling`:
  - the loader never sets `mMeanLumi` (`mapMeanRate` is only filled for CTP or IDC, l. 166-170);
  - it never sets `mInstLumi` (only for the CTP or IDC types).

  So, with the default `CorrMapParam` (`lumiMean = 0`, no override), `mLumiScale = mMeanLumi ? mInstLumi / mMeanLumi : 0` = **0**.
- **What `buildMap` does with it.** In Linear mode it computes `final = mean − ref`, then `final = 0·final + ref`, i.e. **ref**.
- **What the legacy code did.** `TPCFastTransform::TransformLocal` (`TPCFastTransform.h:414-432`) only scales when `scale > 0`. With
  scale 0 it applies the mean map unchanged, which matches the enum's "no scaling, use map as is".
- **The same TF-level effect** occurs for the Linear mode with any type whenever `meanLumi` is 0 (a map without IDC / lumi and no
  override).
- **Context.** Scale 0 now acts as a real scale: zero luminosity gives the reference map. That is physically right when the luminosity
  really is 0. It is wrong when 0 only means "no scaling information".
- **Fix.** In `buildMap`, when `getLumiScaleType() == NoScaling`, or when no valid mean lumi exists, use the mean map unchanged (skip
  the Linear steps), or fail loudly.
- **Reproducer of the scale** (ROOT, O2 environment):

  ```cpp
  o2::gpu::CorrectionMapsHelper h;
  h.setLumiScaleType(o2::tpc::LumiScaleType::NoScaling);
  h.setLumiScaleMode(o2::tpc::LumiScaleMode::Linear);
  h.getLumiScale()   // 0  -> buildMap: 0 * (mean - ref) + ref = ref
  ```

### C2 [certain; medium, configuration-dependent] The IDC → CTP fallback without CTP lumi mixes units

**Where:** `CorrectionMapsLoader.cxx:50-65, 80-93`.

- **What happens.** With an invalid scaler (`tpcScaler == -1`), `canUseCTPScaling` only checks that both maps carry IDC and lumi. It
  then switches to the maps' CTP mean lumis and the CTP type.
- **Where it breaks.** The instantaneous CTP lumi is only read when `getLumiCTPAvailable()` is true (l. 80). With
  `--disable-ctp-lumi-request`, `mInstLumi` keeps the previous TF's IDC value, so the lumi scale is `(IDC_inst − CTP_mean) / CTP_meanRef`
  for every TF until the scaler is valid again.
- **Fix.** Add `getLumiCTPAvailable()` to `canUseCTPScaling`. Otherwise use the "stored IDC" branch, which gives scale 0 in derivative
  mode and 1 in Linear mode.

### C3 [certain; low, dead code] `InverseTransformXYZtoNominalXYZ` uses `nx` before it is set

**Where:** `TPCFastTransformPOD.h:705-723`; also `TPCFastTransform.h:757-775`.

- `c1 = (nx2 - nx) / (nx2 - nx1)` and `c2 = (nx - nx1) / (nx2 - nx1)` read the output reference `nx`. Only afterwards is `nx = x`
  assigned. It should use `x`.
- Since c1 + c2 = 1, the result is an extrapolation with arbitrary weights. It only reduces to the correct answer for a zero
  correction.
- No caller exists in the repository. Fix it or remove it before someone uses it.

### C4 [certain; low, dead code] 3-argument `TransformInTimeFrame` is off by ±L

**Where:** `TPCFastTransformPOD.h:503-507`; also `TPCFastTransform.h:555-559`.

- It computes `l = (time − T0 − maxTimeBin)·vd` and then `z = convDriftLengthToZ(sector, l)`, which gives L − l on side A. At
  `time = T0 + maxTimeBin` that is z = L.
- The "Z(maxTimeBin) == 0" convention of `convTimeToZinTimeFrame` and of the pad/row overload (`convPadTimeToLocalInTimeFrame` adds L
  to l) gives 0 instead.
- No caller exists. Fix it (`l = L + …`) or remove it.

### C5 [certain; low] The POD ignores the source's apply-correction flag

**Where:** `TPCFastTransformPOD.cxx:108` sets `mApplyCorrection = true`, and the overload that takes a `TPCFastTransform` (l. 189-203)
does not copy the flag.

- `TPCScalerSpec.cxx:200` calls `finalMap.setApplyCorrectionOff()` for `NoCorrection`, and that has no effect.
- The corrections stay numerically zero, because `setNoCorrection` zeroes the data.
- But every `TransformLocal` call still evaluates the spline, and `isCorrectionApplied()` returns true.
- **Fix.** Copy `src.isCorrectionApplied()`, or call `pod->setApplyCorrectionOff()` in `buildMap`.
- **Related:** `updateCalibration(*pod, 0, …)` (`TPCScalerSpec.cxx:254`) sets the POD timestamp to 0. The POD's `mIDC` / `mLumi` are
  the CCDB map's mean values, not the TF's. These are informational fields, and no consumer reads them (checked), so it is low.

### C6 [certain; low] The CTP dummy-input alarm limiter never limits

**Where:** `CorrectionMapsLoader.cxx:45-46, 84-86`.

`int dumRep = 0` is local, so `dumRep < maxDumRep` is always true. Make it a member, or `static`, like `lumiPrev`.

### C7 [certain; low] `getBoundaryPotential` loses the window at both ends

**Where:** `TPCMShapeCorrection.cxx:85-106`.

- `lower_bound` returns `idx = 0` for `t <= time[0]`, and `idx = size` for `t > time.back()`. Both are rejected by
  `(idx > 0) && (idx < size)`.
- So a TF exactly at, or shortly before, the first entry, or shortly after the last entry, gets no M-shape although it is within
  `mMaxDeltaTimeMS`.
- **Fix.** Check `idx == 0` against `time[0]`, and `idx == size` against `time.back()`.

### P-C1 [plausible, depends on intent; medium] Once computed, the M-shape correction is never removed

**Where:** `TPCScalerSpec.cxx:110-140, 224-227`.

- `mCorrMapMShape` is only replaced when a non-empty boundary potential is found. When `getBoundaryPotential` returns empty (no stored
  M-shape within 50 ms: gaps, or after the last entry, see C7), the previous M-shape map is still added in `buildMap`, for the rest of
  the run.
- If "empty" means "no M-shape" (which the 50 ms window suggests), the correction should be reset (`setCorrMapMShape(nullptr)`). If
  "hold the last value" is intended, add a comment and a time limit.

### P-C2 [plausible; medium] The sector-edge correction is on for MC

**Where:** `CorrectionMapsTypes.h` (`enableSecEdgeFlucCorrection = true` by default); `CorrectionMapsLoader.cxx:120-124` (comment
"load … only for data").

- Nothing checks for MC: neither the `DerivativeMapMC` mode nor an MC flag disables it.
- In anchored MC, `getSectorsAtTime(run, t)` would return the data run's intervals and add a correction for a distortion that was
  never simulated, unless every MC workflow passes `--disable-sec-edge-fluc-correction`.
- Please check O2DPG, or disable it automatically for `DerivativeMapMC`.

### P-C3 [plausible; low-medium] A scaler that is invalid on one side is not detected

**Where:** `TPCScalerSpec.cxx:159-162`; `CorrectionMapsLoader.cxx:52`.

- `getMeanScaler` returns −1 per side when the data are missing (`TPCScaler.cxx:121-123`), and 0 when the weights sum to 0.
- The fallback only triggers on an average of exactly −1, i.e. both sides invalid. One invalid side gives (A − 1)/2, which is accepted
  as the IDC.
- Times beyond the object's data range get the clamped edge values, also without notice.
- **Fix.** Validate each side (> 0) before averaging, and treat any invalid side as −1.

### P-C4 [plausible, known approximation; low-medium] How the inverse maps and derivatives are merged

**Where:** `TPCFastSpaceChargeCorrectionHelper.cxx:1132-1156`.

- **Inverse maps.** These are combined linearly: `inv_final(r) = s_main·inv_main(r) + Σ s_i·inv_i(r)`. The exact inverse of
  `A + sB` differs by roughly `∇A·sB + s∇B·((1 − s)B − A)`: the static map's gradient times the fluctuation correction. Near the
  inner field cage this can reach about 0.1 mm for cm-scale static maps.
  - This is the **same** approximation as the legacy runtime composition (`TPCFastTransform.h:677-691, 717-738`), so it is not a
    regression.
  - Now that the merged map is built once per TF on the CPU, the inverse could be re-fitted from the merged direct map
    (`initInverse`) for the M-shape and sector-edge parts, or at least monitored with `testInverse`.
- **Edges of the extra maps.** Where the main map's knots lie outside an extra map's grid (clamped `lu` / `lv`), or in its z fade
  region (`ls < 1`), its edge derivatives are copied onto the main knots. This creates a spurious slope beyond the edge. Small, but
  worth a check when the grids differ (the M-shape map has 4 × 4 knots).

---

## 2. Performance (ranked)

### PERF1 [certain mechanism; high when M-shape is enabled] The M-shape map is recomputed every TF

**Where:** `TPCScalerSpec.cxx:110-140`.

- **The cost.** Whenever a stored boundary potential is within 50 ms (`mMaxDeltaTimeMS`, which is also the sampling time), the device
  runs a full `SpaceCharge` Poisson solve (129 × 129 × 5 grid), the field and global-correction integration, and
  `createFromGlobalCorrection`, with a spline fit and the inverse. It does this for **every TF**.
- **Why that is unnecessary.** With 50 ms sampling, consecutive TFs mostly select the **same** stored entry.
- **Fix.** Cache the selected entry index (or its time) and rebuild only when it changes. This also makes P-C1 easy to handle.

### PERF2 [certain; medium] The full map is rebuilt every TF

**Where:** `TPCScalerSpec.cxx:185-256`.

- **Per TF:**
  - `cloneFromObject` copies the whole map.
  - Then one to three `addCorrections` passes run. Each pass touches every knot of every sector and row, for the direct map and both
    inverse maps, with one `interpolateParametersAtU` per knot per extra map. It also spawns `nThreads` threads 36 times (once per
    sector).
  - Then the POD is created (another copy), and the GPU workflow copies it again (`GPUWorkflowTPC.cxx:351-352`).
- **Linear mode wastes a pass.** It needs two passes for `s·(mean − ref) + ref`. A single
  `addCorrections(final, s, {{ref, 1 − s}, M-shape, sector-edge})` gives the same result in one pass. In derivative mode everything is
  already one pass.
- **Further options.**
  - Use a persistent thread pool instead of spawning threads per sector.
  - Skip the rebuild when the scale changes by less than a tolerance and the maps, VDrift and the M-shape / sector-edge state are
    unchanged.
  - Measure the time spent in `buildMap` against the TF rate on the EPN.

### PERF3 [certain; low-medium] The GPU calibration is updated every TF

**Where:** `GPUWorkflowTPC.cxx:345-355`.

`mustUpdate = true` every TF, together with a `memcpy` of the whole POD, means the whole map is uploaded to the GPU every TF. In the
legacy flow this only happened on map updates. It follows from PERF2, but it becomes avoidable once the producer marks unchanged maps.

### PERF4 [certain; low] Info-level logging every TF

Several messages log at info level on every TF: `reportScaling()`, "Adding M-shape correction …", "Checking for sector edge
fluctuation", "fast space charge correction helper: Add corrections", and "Merge of corrections tooks: …s". Use debug level, or log
only on change.

---

## 3. Improvements

1. **Simplify `buildMap`:** one `addCorrections` call per mode (PERF2), and an explicit branch for `NoScaling` (C1).
2. **Cache the `NoCorrection` dummy map.** `TPCFastTransformHelperO2::create(0)` is called every TF (`TPCScalerSpec.cxx:198`).
3. **Replace the dummy-map heuristic.** `isCorrMapMShapeDummy()` treats any map with 4 knots in row 0 as a dummy
   (`CorrectionMapsHelper.h:124-131`). A genuine 2 × 2 M-shape map (`--n-knots-y 2 --n-knots-z 2`) would be silently dropped. Use an
   explicit validity flag instead.
4. **Fix misleading or outdated messages and comments.**
   - "Adding M-shape correction to the final map with scaling factor {}" prints `mMShapeScalingFac`, but the map is added with scale
     1, because the factor was already applied to the potential.
   - "Correction mode unknown! Choose either 0 (default) or 1 …" (`CorrectionMapsLoader.cxx:117`): there are five modes now.
   - "load sector edge fluctuation correction only for data" (P-C2).
5. **`CorrectionMapsHelper::clear()`** does not reset `mCorrMapSecEdgeFluc`.
6. **Remove or fix the unused legacy functions** in both `TPCFastTransform` and the POD (C3, C4). Also add a POD unit test:
   `TPCFastTransformPOD::test` exists, but nothing calls it except `TPCFastTransformInit*.C` and `testGPUCATracking`. A test of
   `addCorrections` against the legacy runtime composition (direct and inverse, random points) would cover the merge math.
7. **Make the POD's metadata meaningful.** Set `mTimeStamp` to the TF time, and `mIDC` / `mLumi` to the instantaneous values actually
   used (or rename them `mMeanIDC` / `mMeanLumi`). Initialise `mLumi` / `mIDC` in the overload that takes only a correction
   (`TPCFastTransformPOD.cxx:117-121`); they rely on the caller's buffer being zeroed.
8. **Message alignment.** Consumers keep a pointer into the DPL message (`TPCFastTransformPOD::get(raw)`). The POD needs 8-byte
   alignment, and its spline parameters assume 16 (`getParameterAlignmentBytes`) relative to the buffer start. Add an assertion on
   the message alignment in `get()`, or document the guarantee.

---

## 4. Follow-up: `addCorrections` with different spline grids, and `initInverse`

### 4.1 Which grids the maps actually have [certain, from the creation code]

- **Maps made from residuals** (`TPCFastSpaceChargeCorrectionHelper::createFromTrackResiduals`, l. 399-509): main, reference and
  derivative.
  - Knots sit at every second y/x and z/x bin centre, with both boundary bins always included.
  - Per row, the y grid spans `rowX · y2x(first bin centre) … rowX · y2x(last bin centre)`, and the z grid
    `rowX · z2x(first centre) … rowX · z2x(last centre)`.
  - The correction fades to 0 between `zOut = min(zMax, 250 cm)` and the readout. `gridReal` is the same as `gridMeasured`.
  - So these maps share a grid exactly when they were made with the same `TrackResiduals` binning (y/x and z/x bins, `maxZ2X`).
- **M-shape map** (`createFromGlobalCorrection` → `createFromLocalCorrection`, with `n-knots-y` / `n-knots-z` = 4 by default): uniform
  knots over the full pad row and over z = 0 … 250 cm. It keeps the default grid of `setNoCorrection`, so it has no fade. Its grid is
  therefore **always different** from that of a residual map.
- **Sector-edge map** (`TPC/Calib/CorrectionMapSecEdgeFluc`): this repository has no code that creates it, so its grid is unknown.

### 4.2 What the merge does when the grids differ [certain]

`addCorrections` evaluates the extra map's value, first derivatives and mixed derivative at the **main map's knots**, and the main map's
Hermite splines then interpolate between them. The legacy runtime composition evaluated every map with **its own** spline. The two
agree exactly only where both of these hold:

- the main map's cell contains no knot line of the extra map (inside one cell of each, both are bicubic);
- the point lies inside the main grid, before the fade.

Otherwise:

1. **A main cell crosses a knot line of the extra map.** The extra map is a compact spline, continuous only in its first derivatives,
   so it cannot be reproduced exactly within that cell.
2. **The point is outside the main grid.**
   - In y, this is the strip between the outermost bin centre and the pad edge: 0.9 cm in row 0 and 1.8 cm in row 151 with 20 y/x
     bins.
   - In z, it covers z/x beyond the last bin centre (for example z > 117 cm in row 0 with 28 bins up to 1.4) and the region near the
     central electrode.
   - There the main map clamps u and v, so the **extra map is frozen at the main grid's edge value**. In z it is also multiplied by the
     main map's fade, `s_main(z)`, instead of following its own z profile.
3. **Main knots outside the extra map's grid.** There the extra map's slope at its own edge is copied onto the main knots, which
   creates a spurious slope.

### 4.3 Size, from the 1D model (`docs/spline_merge_model.py`, illustrative profiles)

| Extra map | Region | max \|merged − legacy\| |
|---|---|---|
| Derivative map with the main map's binning | everywhere | 0 (exact) |
| M-shape-like, 1 mm amplitude, 4 knots; main map with 28 z/x bins up to 1.4 | inside the main grid (row x = 85 / 200 cm) | 0.3 µm / 20 µm |
| same | row x = 85 cm, z beyond the main grid (z > 117 cm, z/x > 1.375) | **450 µm** |
| same, main map with 5 z/x bins up to 1.0 | row x = 85 cm, z > 77 cm | **930 µm** |
| Sector-edge-like, 1 mm at the pad edge, 1 cm decay, 40 knots over the full row; main map with 20 y/x bins | within 3 cm inside the main grid edge (x = 85 / 245 cm) | 16 µm / 85 µm |
| same | between the outermost bin centre and the pad edge | **600-820 µm** (most of the edge correction is lost) |
| same | interior, x = 245 cm (main knots 8 cm apart) | up to 90 µm of ripple |

**Conclusions.**

- **M-shape (certain).** Merging changes the M-shape correction for inner-row clusters at large |z|/x, i.e. forward tracks: there it is
  frozen at the main grid's edge and faded towards the readout. Inside the main grid the merge is essentially exact.
- **Sector-edge correction (plausible).** If its CCDB map has a finer grid, or one that reaches the pad edges, the merged per-TF map
  loses most of the correction in the last 1-2 cm of each pad row, which is exactly where it is meant to act. If it was made with the
  main map's binning, the merge is exact.
- **How to tell.** `docs/checkCorrMapMerge.C` prints the knots of both maps and compares merged against legacy over the whole area, so
  one run on the real objects decides it.
- **Possible fixes.**
  - Merge onto a grid that covers the union of the inputs, in coverage and resolution.
  - Or keep maps with incompatible grids as a separate runtime term, as the legacy `ref2` path did.
  - Or at least give the residual maps grids that reach the pad edges and the readout.
  - Refitting the sum on the main grid would fix only point 1, not the coverage.

### 4.4 `initInverse` (`TPCFastSpaceChargeCorrectionHelper.cxx:886-1014`)

**I1 [certain mechanism, plausible size; medium] The real grid is the measured grid, and points outside are clamped.**

- **What the code does.**
  - l. 905-908 sets `gridReal = gridMeasured` for every row.
  - The inverse is sampled on the *measured* grid (ten points per knot interval). Each sample is placed at
    `real = measured + d` and clamped into the grid (`convRealLocalToGrid`).
- **The consequences.**
  - Where the distortion points outward at an edge, those samples pile up on the edge line with full weight and with their own values
    of `d`.
  - Where it points inward, the strip of the real grid that no sample reaches is fixed only by the `1e-8` smoothness of
    `approximateDataPoints`.
  - The `gridReal` field exists to describe the real-coordinate domain, but this setting throws that away.
- **Size in the 1D model** (row x = 85 cm, 20 y/x bins):
  - a 2-5 mm shift plus a 2 % gradient: 45-70 µm in the last knot interval, ≤ 6 µm further inside;
  - a 3-5 mm edge bump with 2 cm decay (a sector-edge-like distortion): **470-690 µm** in the last knot interval and 40-55 µm further
    inside.
- **Visibility.** `testInverse` sees most of this, because the large errors occur where both the measured and the real point lie inside
  the grid. The creation macros print its maximum, so the map producers can check their logs. It does skip real points whose measured
  partner lies outside the grid (≈ 50 µm in the model).
- **Fix.**
  - Set `gridReal` per row from the actual range of `measured + d`.
  - Also sample measured points beyond the measured grid (the direct map is constant there), so that the real grid is fully covered.
  - Drop points that fall outside instead of clamping them.

**I2 [certain; low] The fade factor is applied twice in the inverse.** Each inverse data point is stored as
`scale · d` (l. 974-976), with `scale` the fade factor at the real z, and `getCorrectionXatRealYZ` / `getCorrectionYZatRealYZ` multiply
by the fade factor again. This only affects samples whose real z passes `zOut`, at the ≲ 0.5 % level. Store `d`.

**I3 [certain; low] Every interior knot line is sampled twice.** `for (double s = 1.; s > 0.; s -= 0.1)` (l. 931, 944) runs **11** times:
the last value is s = 1.4·10⁻¹⁶ (checked with IEEE doubles). Harmless double weight; use an integer loop.

**I4 [certain; low] A related dead check in `fillSpaceChargeCorrectionFromMap`** (l. 157-160). It is used for the direct fit and for the
inverse from an inverse voxel tree. `if (scale - 1.f > 1.e-6) continue; // point is outside the grid` can never be true, because the fade
scale is ≤ 1. Points outside the grid are therefore not skipped but clamped onto the edges. The intended test is probably
`!correction.isLocalInsideGrid(...)`, or `1.f - scale > 1.e-6` for the fade region.

**I5 [low]** Inside the sampling loop, `float gridU, gridV` (l. 970) shadow the `std::vector<double> gridU, gridV` of the enclosing
scope. This is legal, but confusing.
