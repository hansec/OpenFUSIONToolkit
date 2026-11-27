# PR #457: MUG quadrature-point cache for the xMHD residual and Jacobian

- **PR:** [OpenFUSIONToolkit/OpenFUSIONToolkit#457](https://github.com/OpenFUSIONToolkit/OpenFUSIONToolkit/pull/457),
  "MUG: Quadrature-point cache for xMHD residual and Jacobian assembly"
- **Head branch:** `claude/xmhd-qcache-final` on the fork (`hansec/OpenFUSIONToolkit`)

## Summary

The 3D xMHD module (`src/physics/xmhd.F90`) evaluated every basis function, gradient and curl at every
quadrature point of every cell in each nonlinear residual and Jacobian build. The PR adds a quadrature-point
cache (`src/fem/fem_quad_cache.F90`): cell-independent reference tables of the logical components of the basis
functions, indexed by entity orientation, plus per-cell geometry. Cell-batched kernels reconstruct fields and
contract test functions from these tables. The pointwise physics of the residual and the Jacobian was factored
into single routines (`xmhd_resid_qp`, `xmhd_jac_qp`) shared by the direct and cached paths. Local Jacobian
matrices are formed by batched BLAS `dgemm` contractions, and assembly into native CSR matrices uses
precomputed entry positions instead of column searches. A residual scatter without critical sections and
lock-free Jacobian adds for cells away from thread-partition boundaries were also added.

On a 4-core Xeon VM the NL residual is 1.4-2.0x faster on tets and about 4x on hexes, and the Jacobian build
3-8x faster at one thread (2.7-7.3x at four), with identical linear/nonlinear iteration histories and
round-off-level agreement in every tested configuration. The new options are `xmhd_qcache` (0 direct, 1 cache,
default), `xmhd_jac_posmap`, `xmhd_jac_posmap_maxmb` and the experimental `xmhd_coarse_quad`. Multi-threaded
runs need a BLAS that runs single-threaded inside OpenMP regions (documented in the PR).

## Files

| file | content |
|---|---|
| `RESULTS.md` | Final report: options, design, validation (native, PETSc, geometric MG, solid regions), performance tables (residual, Jacobian at 1/4 threads, time stepping, memory, coarse-level quadrature) with hardware and libraries |
| `FOLLOW_UPS.md` | Open items and ideas: coverage gaps (linear runs, PETSc position map), performance hotspots, time stepping, GPU porting, BLAS threading detection, profiling tooling, basis choice |
| `harness/` | Snapshot of the benchmark harness removed from the PR before merging: `run_bench.py`, its `README.md` (how to run, output fields, A/B procedure on other platforms) and the Fortran driver `xmhd_bench.F90` |
| `basis_study.md`, `basis_study/` | Analysis study of the basis/orientation structure behind the cache's orientation handling, with standalone Python scripts and their outputs (not adopted; see `FOLLOW_UPS.md`) |

## Branches and commits (fork `hansec/OpenFUSIONToolkit`)

The work was done in a chain of agent sessions, one phase per branch, each starting from the previous one:

| branch | head | content |
|---|---|---|
| `claude/xmhd-qcache-phase0` | `a3c6a74` | benchmark harness and driver, baseline profile of `main` (`phase0_baseline.md`, raw results, perf reports) |
| `claude/xmhd-qcache-phase1` | `b10a39c` | conditional `log2phys`, partition-aware lock-free residual scatter and Jacobian adds (`phase1_results.md`) |
| `claude/xmhd-qcache-phase2` | `5328aab` | `fem_quad_cache` module with FULL and REFERENCE modes, kernel benchmark (`phase2_results.md`) |
| `claude/xmhd-qcache-basis-study` | `d344ffe` | basis study (copied here as `basis_study*`) |
| `claude/xmhd-qcache-phase3` | `b6f2c8b` | NL residual on the cache (`phase3_results.md`, equivalence logs, A/B results) |
| `claude/xmhd-qcache-phase4` | `d36f2b2` | Jacobian on the cache, CSR position map, cache as default (`phase4_results.md`, valgrind, A/B results) |
| `claude/xmhd-qcache-final` | | PR head: FULL mode removed, options simplified to 0/1, documentation, PETSc validation, merges of `main` |
| `claude/xmhd-qcache-nla-compat` | `5bdbfcd` | scratch merges of the PR branch with in-progress native linear-algebra work (compatibility testing only) |

Commits on the PR branch referenced from these notes:

| commit | meaning |
|---|---|
| `c466aae` | last PR head **with** the development tooling: `benchmarks/xmhd_qcache/` (harness, `RESULTS.md`), the driver `src/tests/physics/MUG/xmhd_bench.F90`, and the profiling instrumentation in `xmhd.F90` (`xmhd_prof_*`, `xmhd_prof_report`, `XMHD_PROF` output) |
| `dd58222` | development tooling removed for merging; the equivalence driver became `src/tests/physics/MUG/test_MUG_qcache.F90` |

The per-phase reports and raw outputs (per-run JSON/CSV, perf reports, equivalence logs, valgrind output) are
only on the phase branches. Their option values differ from the final ones (`xmhd_qcache=2` on the Phase 3/4
branches is `xmhd_qcache=1` in the merged code, and the FULL mode `xmhd_qcache=1` there no longer exists);
`RESULTS.md` uses the final values.

## Reproducing

- **Correctness** (merged code): `src/tests/fem/test_fem_qcache.py` (basis cache against the direct evaluators)
  and `src/tests/physics/MUG/test_MUG_qcache.py` (residual and Jacobians of every level, `xmhd_qcache` 0 vs 1, with
  and without the position map, native and PETSc), run from `<build_dir>/tests` as usual.
- **Benchmarks:** check out `c466aae` and follow `benchmarks/xmhd_qcache/README.md` there (identical to
  `harness/README.md`). The harness needs the profiling instrumentation of that commit. To benchmark a later
  version of the code, restore the instrumentation (or an equivalent timing facility, see `FOLLOW_UPS.md`) on a
  development branch and copy `harness/` back to `benchmarks/xmhd_qcache/` and
  `src/tests/physics/MUG/xmhd_bench.F90`.
- Record the CPU, core count, compiler, BLAS variant and MPI with any new numbers, and compare only runs from the
  same machine (see `harness/README.md`, "Testing on another platform").
