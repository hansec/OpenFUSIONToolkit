# Follow-ups: xMHD quadrature-point cache (PR #457)

Open items, limitations and ideas collected while developing and validating the quadrature-point cache.
Option names are those of the merged code (`xmhd_qcache=0/1`, `xmhd_jac_posmap`, `xmhd_jac_posmap_maxmb`,
`xmhd_coarse_quad`). Numbers refer to `RESULTS.md`.

## Coverage of the cache and the position map

- **Linear runs.** `xmhd_massmatrix_apply` (used by `xmhd_lin_run`) still evaluates the basis directly. It is
  applied once or twice per time step and linear runs are dominated by the assembled-Jacobian solve, so the gain
  would be small. The pattern is the same as the residual: a shared pointwise kernel plus the cell skeleton of
  `xmhd_resid_qcache`.
- **PETSc matrices** are assembled with the generic `fem_mat_add_local` (per-field locks, column search). A
  position map would need access to PETSc's CSR layout (eg. `MatSeqAIJGetArray` and the local-to-global
  mapping); the cache path itself already works with PETSc (validated, see `RESULTS.md`).
- **Geometric MG levels** (coarser meshes) use the direct path with cell-centered fields. They are cheap compared
  with the finest mesh, so this was not worth a cache of their own.
- **Position map memory.** `4*sum(nce_s*nce_t)` bytes per cell and level (25 KB/cell for tets p3, 275 KB/cell for
  hexes p3). The default cap is 1 GB per level and rank (`xmhd_jac_posmap_maxmb=1024`); with several p-levels
  the total can reach the cap times the number of levels. A block CSR layout over the field blocks would shrink
  both the matrix index traffic and the map.

## Experimental coarse-level quadrature (`xmhd_coarse_quad`)

Coarse p-level Jacobian builds are 1.6-1.9x faster with the level's own (lower) quadrature order, and the one run
case tested (`tet_p3_mg_run`) had identical iteration counts at every step. Before considering it as a default,
check p-MG convergence on harder problems: strong advection, large `dt`, strong anisotropic conduction, curved
and hex meshes, and p=4 hierarchies. It only changes the preconditioner's coarse operators, never the solution.

## Performance

- **Jacobian build** (tet p3, 1 thread, after this work): scatter with the position map 27% (memory bound, about
  34k entries per cell into a 27 MB value array), trial-function contraction `xmhd_jac_wform` 26%, `dgemm` 17%,
  `xmhd_jac_qp` 8%, BCs/scaling 6%, field reconstruction 6%, basis rebuild 5%. On hexes p3 `dgemm` is 57%.
  Options:
  - compact the coefficient tensor `D` (at most 43 x 43 per qp, about a quarter non-zero for the default
    two-temperature model) to its static non-zero pattern, so `xmhd_jac_wform` loops over non-zeros only;
  - skip zero rows (Dirichlet rows) in the scatter;
  - reorder cells for locality of the matrix rows they touch.
- **NL residual** (tet p3, cache): about 75% of the time is the two small dense contractions per cell and
  quantity (reconstruction and test-function accumulation, eg. 20 x 105 x 6 for Lagrange p3), which run at about
  2 GFLOP/s per core in the generic `-O2` (SSE2) build. `-march`-level vectorization of the Fortran kernels
  (OpenBLAS already dispatches AVX-512 at run time), or batching several cells into one GEMM, would help.
- **Time stepping** is now dominated by the linear solve and the preconditioner (the Jacobian dropped from about
  8% to 2.5% of a tet p3 p-MG step). Candidates:
  - reuse the symbolic ILU/LU factorization between preconditioner updates (the sparsity pattern is fixed);
  - retune `xmhd_prefreq`/`nl_update` now that Jacobian rebuilds are cheap relative to factorizations;
  - re-measure "assembled Jacobian for the matvec, matrix-free only for the residual" in MFNK;
  - deeper p hierarchies (p=4 -> 3 -> 2 -> 1) or more coarse-level smoothing;
  - block CSR (BSR) storage over the field blocks for matvecs and triangular solves.

## GPU porting

The cell loops already have the structure of batched GPU kernels:
- per cell: indexed gathers through the DOF tables, cell-independent reference tables (they fit in shared memory
  up to p=3), per-qp 3x3 geometry, a pointwise kernel without allocation or type-bound calls, two batched GEMMs per
  row space (eg. tet p3: Lagrange `20 x 35 x ~600` and `20 x 105 x ~240`, H(Curl) `45 x 105 x ~400`), and a
  scatter with precomputed positions (atomic adds, no search);
- needed first: collect the module settings read by `xmhd_resid_qp`/`xmhd_jac_qp` into one constant structure,
  replace the `res_profile` procedure pointer by a tabulated profile, precompute the shaping factors per
  (cell, qp), and provide device versions of `brag_*`, `elec_ion_therm_rate`, `cross_product`, `magnitude`;
- orientation handling (signed permutations, dense blocks, the per-orientation H(Curl) face blocks on tets) and
  boundary conditions are per cell and branchy: bin cells by orientation signature and by kind (fluid/solid), and
  apply BCs as a separate pass over flagged rows;
- the `dt<0` diagnostics and floor flags are per-cell reductions.

## BLAS threading

The Jacobian calls `dgemm` from inside OpenMP parallel regions. With a pthreads-threaded OpenBLAS (eg. Ubuntu's
default `libopenblas0-pthread`) the build stops scaling unless `OPENBLAS_NUM_THREADS=1` is set (tet p4: 2.33 s at
1 thread vs 2.21 s at 4 threads; 0.71 s at 4 threads with the variable set). This is documented (xMHD options
table, MUG documentation page, main page, README, `oft-build` skill), but not detected at run time. A portable,
cheap detection was not found: `openblas_get_parallel` exists only in OpenBLAS, and calling it would break
linking against MKL, Accelerate or reference BLAS. A CMake-time check of the BLAS variant, or a weak symbol in
the C layer, would be possible.

## Profiling and benchmarking tooling

The profiling counters (`xmhd_prof_*`, `xmhd_prof_report`, the `XMHD_PROF` output lines, configurable
repetitions of profile mode) and the benchmark harness were removed from the PR before merging; snapshots are in
`harness/`. They need the instrumentation of the pre-cleanup PR head (see `README.md`) to run.

- If performance work on xMHD continues, consider a permanent, lightweight timing facility (per-level Jacobian
  build times, residual and preconditioner-update totals, cache and position-map memory) that the harness could
  parse, instead of carrying the instrumentation on a development branch.
- Profile mode (`xmhd_run(..., profile_only=.TRUE.)`) on `main` uses fixed 10 x 10 x 10 repetitions; the
  development branch made them configurable (`xmhd_prof_reps`) and printed checksums of `A(u)` and `J u`
  (`XMHD_PROF nlop_norm`/`matvec_norm`) for comparing runs. Both are useful for any future A/B work.

## Basis choice (from `basis_study.md`)

The orientation handling on tets (dense H^1 face transforms at order >= 4, own per-orientation table blocks for
H(Curl) faces at p >= 2) comes from OFT's face function families, which are not invariant under the face's
symmetry group. The study shows the remaining obstruction is structural (a single 2x2 block per face cannot be
removed by any hierarchical basis), but that a symmetric, hierarchical choice of face functions would make every
transform an own-entity signed permutation plus one 2x2 block per face. Ranked recommendations from the study:
1. a symmetric hierarchical face family as an **opt-in** basis variant for H(Curl) (p >= 2) and H^1 (order >= 4);
2. without a basis change, cross-space transforms (H(Curl) faces expressed with Grad(H^1)) for composite B fields;
3. a mesh-side alternative (sort tet vertices by global index with a parity swap, leaving two orientation
   classes) only if the per-face 2x2 blocks matter on GPUs;
4. no wholesale adoption of a basis from the literature.

It also found that OFT's H(Curl) **cell-interior** block is not invariant under relabeling a cell's vertices, so
the standalone H(Curl) space depends on each cell's local vertex order (the composite H(Curl) + Grad(H^1) space
does not). This is a caveat for any mesh-reordering option and is fixable by projecting the interior levels off
the H^1 interior gradients (section 3.7 of the study). Changing the basis would change basis-dependent test
reference values and restart files, so none of this was done in the PR.

## Smaller items

- `test_MUG_qcache.py::test_qcache_curved_p3` skips the Jacobian comparison because the p=3 Jacobian of that
  mesh needs about 8 GB for the check; `test_qcache_curved_p2` covers the curved Jacobian.
- On curved cells the cache and the direct Jacobian differ by up to 2.4e-12 of the row maximum at p=3 (7e-13 at
  p=2), in the electron temperature / H(Curl)-curl couplings: the cache rebuilds physical curls from folded
  logical components with the per-qp `cgop`, and the evaluator combines them in a different order with
  cancellation. Harmless, but the reason the test tolerance is `rel_row <= 1e-11` rather than 1e-13.
- No canonical-frame variant of the Jacobian (`fem_qcache_mat_from_ref`): the basis rebuild is only 5-11% of a
  build, and the canonical frame does not apply to H(Curl) faces on tets.
