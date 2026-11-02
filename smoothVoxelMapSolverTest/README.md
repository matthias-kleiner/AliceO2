# smoothVoxelMap.C solver test (TEST ONLY, branch `tpc-fine-map-gbl-solver`, not for upstream)

`Detectors/TPC/calibration/SpacePoints/macro/smoothVoxelMap.C` smooths a TPC space-charge voxel map by a penalised fit per sector:
one node (dX, dY, dZ) per voxel of the grid (152 pad rows × y/x bins × z/x bins), the voxels with data entering with their χ² as a
3 × 3 quadratic form at their position (bilinear between the 4 nodes of their pad row in y/x and z/x), plus second-derivative penalties
along the pad rows, y/x and z/x (`Stencil`: λ (Σ_k c_k x_{first + k·stride})², 3 nodes each, none across the stack boundaries), and an
optional ridge. `solveSector` solves the normal equations by conjugate gradients (CG) with a block-Jacobi preconditioner (the inverse of
each node's 3 × 3 diagonal block), tolerance 1e-8 on the relative residual, at most 20000 iterations.

**The problem:** a stiffer prior (the taus ÷ 30, "t30" = 33 / 67 / 33 µm/cm² instead of 1000 / 2000 / 1000) gives better maps on tracks
and is the candidate default, but needs 3–4× the CG iterations (t1: 4700–6500 per sector, t30: 17500–20000, 7 of 36 sectors hit the
cap at a residual of 2.7e-8). The block-Jacobi preconditioner ignores the coupling between nodes, which the stiff penalty makes dominant.

## Files

- `solverTest.cxx`: includes the macro itself with `SMOOTHVOXELMAP_SOLVER_ONLY` (only the solver part, no ROOT / O2), reads dumped
  sector systems, solves each from zero with the macro's `solveSector` and compares with the macro's solution stored in the file:
  iterations, residual, wall time, rms / max |x − x_ref| per component in µm, and the objective F(x) at both solutions
  (F(new) − F(ref) ≤ ~0: at least as converged).
- `data/<tag>_sectorNN.bin`: sector systems + solutions dumped by the macro (`SMOOTHVOXELMAP_DUMP=<prefix>`; format in the macro,
  `writeSectorSystem`) from a real PbPb 2025 map (run 568041, one 5-min window), for t1 and t30, in the sectors with the most iterations.

```
c++ -std=c++20 -O2 -o solverTest solverTest.cxx
./solverTest data/*.bin
```

## Task: a better preconditioner (or solver) in `solveSector`

Goal: at t30 far fewer iterations and less wall time per sector (target: ≤ the t1 cost today), the same solution to within the CG
tolerance, at t1 not slower than today.

Constraints:
- Change only the solver inside `smoothVoxelMap.C` (`solveSector` and helpers in the `SMOOTHVOXELMAP_SOLVER_ONLY` part). No new
  arguments, knobs or parameters of the macro; the tolerance 1e-8 and the output unchanged; std C++ only (no Eigen or other libraries:
  ACLiC on the GRID has none); O2 style (one declaration per line, descriptive names, `clang-format` with the repo's `.clang-format`).
- `solveSector` runs concurrently for different sectors (one thread per sector, 8 threads over 36 sectors): no shared mutable state.
  Single-threaded within a sector is fine.
- The preconditioner must stay symmetric positive definite (CG). Memory per sector: a few times the current is fine.
- The robust cycles (`robust > 0`, Tukey reweighting, warm start from the previous cycle) must keep working.

Hints (not prescriptions): the nodes are indexed `row + nX · (y/x bin + nY2X · z/x bin)` (`strides` = {1, nX, nX · nY2X}); the data
couple nodes only within a pad row (bilinear in y/x and z/x), the row penalty couples consecutive rows of a column, and the row spacing
(~0.75–1.5 cm) is much smaller than the y/x and z/x spacing, so the row-penalty terms are the stiffest. Candidates: a line (column)
preconditioner — per (y/x, z/x) column the exact block-banded (3 × 3 blocks, bandwidth 2 rows) solve of the data diagonal blocks + the
row penalty; symmetric line Gauss–Seidel / ADI over the three directions; a geometric multigrid V-cycle as preconditioner. Measure,
don't guess: report iterations and time per sector for t1 and t30 before / after in the test.

Acceptance (in `solverTest`, all files): finite; final residual ≤ 1e-8 where the reference converged; rms |x − x_ref| ≤ 5 µm per
component (the reference itself is converged to 1e-8; 1e-8 vs 1e-12 changed the map by 1–4 µm rms); F(new) − F(ref) ≤ 1e-9 · |F(ref)|.

Deliverable: commit on a new branch `tpc-fine-map-gbl-solver-precond` from `tpc-fine-map-gbl-solver`, the solver change in its own
commit(s) touching only `smoothVoxelMap.C` (harness changes, if any, separate), with the before / after table of the test in the commit
message. It is validated afterwards on the full map on GSI (the whole map identical within the tolerance, wall time on 8 threads).
