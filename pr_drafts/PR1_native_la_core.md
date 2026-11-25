**Title:** Linear Algebra: Native LA core refactor (access modes, kernel layer, shared patterns, scatter plans)

**Branch:** `hansec/OpenFUSIONToolkit:claude/native-la-core` (b46866a) → `OpenFUSIONToolkit/OpenFUSIONToolkit:main`

## Naming your pull request
This change is not specific to one physics tool, so the title starts with `Linear Algebra:`.

## Primary changes
Refactor of the native linear algebra backend (`src/lin_alg/native_*`, `src/base/oft_stitching.F90`, new `src/base/oft_scatter.F90`), in four parts. Results are unchanged except where noted (bitwise identical in most paths; see "Behaviour changes").

**1. Local performance and correctness fixes**
- Non-zero balanced contiguous row ranges per thread in CSR products; cached transposed (CSC) index for `applyt`; cached diagonal locations for assembly.
- Native norms and mixed real/complex dots without temporary copies; threaded value copy in `update_slice`; fused sweeps in the symmetric Jacobi smoother.
- Block-Jacobi keeps persistent local vectors instead of allocating per apply.
- Fixed an orphaned worksharing loop in `vec_get_local` (block path) and sizing of connection flags in `native_matrix_setup_full`.
- Global stitching keeps persistent work arrays and posts MPI sends from the calling thread instead of inside a worksharing loop under `!$omp critical`.

**2. Access modes, `copy`/`maxpy`/`mdot`, bitwise-consistent stitching**
- Abstract vectors get `borrow_local`/`release_local` with access modes (`OFT_ACCESS_READ`, `_WRITE`, `_READWRITE`, `_ACCUMULATE`). Native vectors return a pointer into their own storage (no copy); other backends get a copying default. `release_local(...,sync)` lets callers skip the MPI exchange when shared values are already consistent (checked in `OFT_DEBUG` builds by `native_check_shared`).
- `copy` and single-pass `maxpy` (bitwise equal to successive `add` calls) for real and complex vectors; single-pass, single-reduction `mdot`. Solvers (GMRES orthogonalization and updates, CG, smoothers) and operators use them instead of `add(0,1,x)` loops.
- BCs, multigrid base transfers, block preconditioners (FEM, TokaMaker, ThinCurr, block-Jacobi) borrow instead of copying.
- Summation stitching is bitwise consistent across all copies of an element (elements with three or more contributions, eg. periodic images plus remote copies, are summed in ascending order).

**3. Kernel layer and shared pattern object**
- `native_kernels.F90`: all inner loops of the native vector/matrix/solver classes as plain-array kernels with templated value types (real, complex and mixed), so threading/vectorization is applied in one place.
- `oft_native_pattern`: reference-counted CSR non-zero pattern (with cached linkage, transposed index and diagonal locations) shared between matrices. `create_matrix` takes an optional `pattern` cache; FE representations, the ThinCurr resistance matrix and the TokaMaker linearized LHS/RHS share patterns.

**4. Scatter abstraction and accumulate access**
- `oft_scatter`: scatter plans (send/receive lists by neighbour plus a weighted local part, flat persistent buffers, own MPI requests and tag), forward/reverse, insert/add, split into begin/end. Global stitching runs on a plan built from the seam; slice gathers, block-Jacobi block gather/scatter and submatrix value gathers use plans.
- `OFT_ACCESS_ACCUMULATE`: zeroed accumulation storage whose contributions are summed across processes once, on the completing release. Used by element-assembled results (MUG matrix-free operators and solid-node flags, Lagrange/H1/H(Curl)/H(Curl)+Grad projections and matrix-free transposes).

**Merge with `main`** (23 commits since the base): conflicts were reconciled with #450 (block-Jacobi: no parallel region for non-thread-safe sub-solvers or inside an active region; order-independent initial guesses with overlap), #456 (owned-only reductions, extended here to the native vector kernels using the same ownership weights; `clamp_nonneg` kept), #458/#461/#422 (accumulate borrows re-applied on the new `xmhd_lag`/`xmhd_2d` code, keeping main's periodic repack and `get_field_weights`), #467 (documentation updated for `MPI_THREAD_SERIALIZED`).

## Secondary changes
- Fixed latent bugs in unused native matrix paths: `keep_diag` was inverted in `zero_rows`/`zero_cols` (CSR and dense, real and complex), the diagonal check when zeroing rows of a block compared with the block-local row, and zeroing columns of a block used the row offset. No caller passes `keep_diag` or zeroes block columns, so results are unchanged.
- `OFT_DEBUG` is defined for Debug builds. It enables `native_check_shared`, which verifies on `release_local(...,sync=.FALSE.)` that shared values really are consistent across copies (no MPI for process-local seams).
- New regression test `lin_alg/test_native_la` (serial, 2 and 4 processes, periodic meshes): transpose products, norms, complex/mixed products, access modes, `copy`/`maxpy`/`mdot`, BC operators, stitching consistency, accumulate access, shared patterns, row/column zeroing.
- PETSc variants of the gyroviscosity, Hall, viscosity and Taylor–Green MUG tests; generic (non-native) accumulate borrows no longer read the vector; MUG examples only swap the block-Jacobi sub-solver for native backends.

## API changes
No input file or XML changes. Fortran API additions/changes:
- New constants `OFT_ACCESS_READ`, `OFT_ACCESS_WRITE`, `OFT_ACCESS_READWRITE`, `OFT_ACCESS_ACCUMULATE` (`oft_la_base`).
- New type-bound procedures on `oft_vector`/`oft_cvector`: `borrow_local`, `release_local`, `copy` (complex: `copy_real`/`copy_complex` behind generic `copy`), `maxpy`; native overrides of these and of `mdot`.
- New `oft_native_pattern` type and `native_pattern_create/retain/release`, `native_matrix_set_pattern`, `native_matrix_adopt_pattern`, `native_matrix_share`, `native_pattern_live` (`oft_native_la`); `create_matrix` has an optional `pattern` argument.
- New module `oft_scatter` (`oft_scatter_plan`, `scatter_plan_setup[_local]`, `scatter_begin/end/wait_recv/unpack/wait_send/local`, `scatter_plan_check`, `OFT_SCATTER_INSERT/ADD`, `OFT_TAG_STITCH/SCATTER`); new module `oft_native_kernels`.
- Removed `send`, `recv`, `csend`, `crecv` components of `oft_seam` (stitching buffers now live in the seam's scatter plan). The ownership weights from main (`own`, `nown`, `seam_build_own`) are now also used by the native reductions.
- Behaviour changes: native `set_values` aborts on entries outside the pattern (as `add_values` already did; it previously wrote to a wrong location) and the block form locates entries like `add_values`; `add` with zero scale on the destination no longer reads it (non-finite old values do not propagate); complex native `mdot` skips the MPI reduction for complete local meshes like the real one.

## Testing
Ubuntu, GCC 13, OpenMPI, OpenBLAS, SuperLU, UMFPACK, ARPACK, HDF5, METIS (no PETSc, so PETSc variants are skipped), `OMP_NUM_THREADS=2`, 4 cores.

| Suite | Branch | `main` (4cb430f, same sandbox) |
|---|---|---|
| Release `make test` | 3 failed, 1163 passed, 394 skipped (44:07, at b46866a) | 3 failed, 1154 passed, 378 skipped (44:53) |
| Checked build (`-DOFT_DEBUG`) `make test` | 3 failed, 1163 passed, 394 skipped (43:01); no shared-value check failures | – |
| Slow MUG tests, checked build | 110 passed, 73 skipped, 0 failed (`test_alfven` 23/12 in 54:49, `test_alfven_lag` 23/12 in 27:02, `test_sound` 32/24 in 22:02, `test_sound_lag` 32/25 in 18:03); no shared-value check failures | – |

The 3 failures are the same on `main` with the same values: `lin_alg/test_native_bjacobi.py::test_ILU_mpi4_nlocal2[0]` and `..._part[0]` (iteration counts 77 vs 73 and 76 vs 70 expected) and `physics/test_tracing.py::test_poincare[1-True-4]` (tracer timeout).

**Numerics vs `main`** (`lin_alg`, `fem`, `physics/MUG`, result files compared per test, branch run twice to measure run-to-run variation from threaded reductions):
- Most iteration-count differences (25, all ±1–5) are within the branch's own run-to-run spread. Eight `fem` cases (mostly periodic or multigrid) differ reproducibly by 1–3 iterations because of the changed summation order (owned-only reductions, bitwise-consistent stitching); converged values agree to ≤1e-11.
- MUG values agree to ≤2e-5, except the matrix-free cubic nonlinear `test_sound_2d` cases, where `main` differs by up to 1.3e-3 and the branch matches both the assembled-matrix result and the pre-merge branch (~1e-5). All are well inside the tests' 5% tolerance.

**Performance** (536,769 DOF Lagrange problem, range over 3 runs, branch vs `main`):

| Operation | 1 process, 2 threads | 2 processes, 2 threads each |
|---|---|---|
| SpMV (mean per call over 200) | 14.9–15.9 ms vs 16.8–19.4 ms | 16.2–16.4 ms vs 14.7–17.0 ms |
| GMRES(20)+Jacobi, 400 its | 12.4–12.6 s vs 15.7–16.9 s | 12.3–12.8 s vs 13.3–14.4 s |
| CG + multigrid (27 / 102 its, both) | 2.6–2.9 s vs 4.2–4.4 s | 10.1–10.2 s vs 12.4–13.0 s |

Solution norms agree to ~1e-13 relative.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_01NKQvsAfyx8o9bePd9xRb4f
