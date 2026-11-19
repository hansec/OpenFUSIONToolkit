# Plan: Algebraic Multigrid (HYPRE now, AMGX later) for OFT

Status: design proposal only. No code has been changed. All paths are relative to `src/`.

---

## 1. What the code has today

### 1.1 Solver framework

- **Base classes.** `oft_solver` / `oft_csolver` (`lin_alg/solver_base.F90`) provide `apply`, `update(new_pattern)`, `setup_from_xml`, `view`, `check_thread` and `delete`. Preconditioners are solvers chained through `%pre`.
- **XML factories** (`lin_alg/solver_utils.F90`):
  - `create_solver_xml` (L488) accepts `cg | gmres | sjacobi | block_jacobi | add_schwarz | lu | ilu`.
  - `create_pre_xml` (L691) accepts `jacobi | block_jacobi | add_schwarz | lu | ilu`.
  - `create_native_pre` (L585) and `create_petsc_pre` (L606) dispatch to the native or PETSc class.
  - `create_ml_xml` (L769) builds the geometric/p-multigrid preconditioner from `<smoother direction=...>` and an optional `<coarse><solver .../></coarse>`.
- **Pattern for a new external backend.** `oft_lusolver` (`lin_alg/oft_lu.F90`) does it this way:
  - a `<package>` element chooses the backend;
  - a per-package options group (`<strumpack>`) holds backend options;
  - a C/C++ bridge (`strumpack_bridge.cxx`) is wrapped in `#ifdef HAVE_X`;
  - setup is lazy, on the first `apply`.
- **PETSc preconditioners** extend `oft_petsc_precond` (`lin_alg/petsc_solvers.F90:72`) and only override `setup(pc,dist)`. No PETSc wrapper calls `PCSetFromOptions`/`KSPSetFromOptions`, so `-pc_type hypre` cannot be passed in today.
- **No AMG exists anywhere.** There are no matches for hypre, amgx, gamg or boomer.

### 1.2 Native matrix layout (matters for the bridge)

- **Local storage.** `oft_native_matrix` (`lin_alg/native_la.F90:123`) is a local CSR (`kr`, `lc`, `M`) over all local DOFs, including DOFs shared with other ranks. A matrix-vector product is `y = stitch(A_loc x)`, so shared rows are only partially assembled on each rank.
- **Globally summed rows already exist.** `native_matrix_setup_full` (public) plus `%update_slice` (L2768) fill `Mfull`, the row values summed across ranks. Redundant periodic rows become identity rows. Block-Jacobi already uses this.
- **Ownership.** Owned rows per field block are `map(i)%nslice` / `map(i)%slice`. `setup_petsc_vec` (`lin_alg/lin_alg_utils.F90:168-208`) shows how to build a contiguous global numbering, owned rows first, and push it to ghost copies with `oft_global_stitch`. HYPRE IJ and AMGX both need exactly this numbering.
- **Multi-field vectors are field-blocked, not interleaved.** Fields are stored in blocks at `map(i)%offset`.
- **Dirichlet rows are zeroed and given a 1 on the diagonal, but their columns are kept.** For example, `grad_shaf.F90:5855` and the `fem/*_operators.F90` `get*op` routines do this. The assembled matrices are therefore structurally non-symmetric even for SPD operators.

### 1.3 Linear systems across OFT, and the right AMG for each

| # | System (where) | Operator | Fit | Recommended AMG |
|---|---|---|---|---|
| 1 | Lagrange / H1 3D Laplacians (`oft_lag_getlop`, `oft_h1_getlop`). Used in Marklin divergence cleaning and vacuum fields, MUG `divout` (H1(p+1), `grnd`), MUG examples (heat), `tests/lin_alg/test_solver_xml` | SPD scalar, orders 1-5 | **Excellent** | BoomerAMG. For p ≥ 2, the best option is to use AMG as the `<coarse>` solver of OFT's existing p-multigrid. |
| 2 | Marklin H(curl) curl-curl "WOP" (`oft_hcurl_getwop`, `taylor_hmodes`, `taylor_vac_curr`) | SPSD with no mass term; the kernel is the gradients | **Excellent** | **HYPRE AMS** with β = 0 (`HYPRE_AMSSetBetaPoissonMatrix(NULL)`), used as a CG preconditioner alongside the existing `hcurl_divout` cleaner. |
| 3 | MUG xMHD Jacobian blocks (`xmhd.F90`, `xmhd_lag.F90`): velocity (3-field vector Lagrange), density n, temperature T, magnetic field B (H(curl) `Bc` plus grad H1 `Bg`) | Non-symmetric, 7-10 fields | **Good per field, poor monolithic** | Field split. Per-field BoomerAMG: systems AMG for velocity (`num_functions=3`), scalar AMG for n and T. AMS for the `Bc` block, which has the form M/dt + η curl-curl (AMS α/β form). |
| 4 | ThinCurr resistance matrix R (`tw_compute_Rmat_exp`). R⁻¹ is needed in the ARPACK eigen-solve (`thin_wall_solvers.F90:150`), with a CG+Jacobi fallback | SPD surface P1 Laplacian, plus hole and coil rows | **Moderate** | BoomerAMG. Use it in place of the CG+Jacobi fallback, or for large R. The near-dense hole rows may hurt coarsening, so this is experimental. |
| 5 | Anisotropic parallel diffusion (`oft_lag_getpdop`, `marklin_compute_pardiff`) | SPD, strongly anisotropic along field lines | **Moderate / hard** | BoomerAMG with a high `strong_threshold` (0.5-0.9). Good as a stress test; expect iteration counts to degrade. |
| 6 | All mass matrices (Lagrange, vector Lagrange, H(curl), `hcurl_grad`, TokaMaker projections) | SPD, well conditioned | **Don't** | Jacobi-PCG already converges in O(1) iterations. |
| 7 | TokaMaker `dels` (free boundary). Fixed-boundary and time-dependent variants use LU | 2D; dense BEM boundary block; coil rows | **Don't** | Dense rows break coarsening, and 2D serial LU is faster at these sizes. The fixed-boundary variant is AMG-friendly but gains nothing over LU. |
| 8 | ThinCurr L + dt·R, and iωL + R | Dense or HODLR; complex | **Don't** | AMG does not apply to dense operators. Complex HYPRE is a separate library build that cannot be linked alongside the real build. |
| 9 | Marklin JMLB (∫curl·curl − λ u·curl) | Indefinite above the first eigenvalue | **Experimental only** | AMS on the WOP part as a GMRES preconditioner, for small λ only. |
| 10 | 2D MHD (`xmhd_2d`, `mugtok_td`) | 2D; saddle point (incompressible); dense BEM/coil blocks | **Don't (for now)** | LU remains appropriate. |

ADS needs an H(div) space, which OFT does not have, so ADS is out of scope.

**H(curl) details for AMS (verified in code):**
- OFT tetrahedral H(curl) has exactly one DOF per edge at every order (`fem/hcurl_basis.F90:303-313`). Higher orders add only rotational face and cell DOFs; the missing gradients live in the separate grad(H1(p+1)) field of `hcurl_grad`.
- `hcurl_grad` (`fem/hcurl_operators.F90:344`) applies G matrix-free. Its stencil is `±1` (sign of `mesh%global%le`) from the two vertices of each edge into that edge's first DOF.
- So for **p = 1**, G and the vertex coordinates (`mesh%r`) are all AMS needs.
- For **p > 1**, the same G (with zero rows for face/cell DOFs) is probably still the correct kernel. However, AMS's coordinate-built Nédélec interpolation Π only covers the edge DOFs, so convergence will be sub-optimal.
- Recommendation for p > 1: p-coarsen with OFT's existing ML hierarchy down to p = 1 and use AMS as the `<coarse>` solver. A full high-order Π (`HYPRE_AMSSetInterpolations`) can follow later if needed.

---

## 2. Architecture

The goal is one backend-neutral "AMG preconditioner" layer. HYPRE goes in now; AMGX later reuses the same matrix-conversion code.

```
oft_solver
 └─ oft_amg_solver            (abstract, new: lin_alg/oft_amg.F90)
     ├─ oft_hypre_boomeramg   (native, HAVE_HYPRE)
     ├─ oft_hypre_ams         (native, HAVE_HYPRE; needs G + coordinates)
     └─ oft_amgx_solver       (future, HAVE_AMGX)
oft_petsc_precond
 └─ oft_petsc_amgprecond      (PETSc backend: PCHYPRE boomeramg/ams, PCGAMG, PCAMGX later)
```

### 2.1 Shared layer `oft_amg_solver` (Fortran, backend-neutral)

It holds what every AMG package needs, so that backends see only a "distributed owned-row CSR".

- **Global numbering.**
  - Compute a contiguous numbering with owned rows first, using an MPI exclusive scan of the summed `map(:)%nslice`, then push it to shared and ghost copies with `oft_global_stitch`. This is the same algorithm as `setup_petsc_vec`.
  - Store it as `INTEGER(i8) gidx(n)` plus `(ilower, iupper)`.
  - Options:
    - *Field-blocked* numbering (the default; matches OFT storage).
    - *Interleaved* numbering when every block shares an identical map, as with vector Lagrange. Interleaving costs nothing because we own the numbering, and it lets BoomerAMG nodal coarsening and AMGX block formats work.
- **Owned-row CSR extraction.**
  - Call `native_matrix_setup_full` and `update_slice` to get `Mfull`.
  - Keep the owned rows only and map their columns through `gidx`.
  - Inserting owned rows only means HYPRE needs no off-process `AddToValues` communication.
  - The pattern (row pointers and global columns) is cached, so a value-only refresh just re-gathers values.
- **`symmetrize_bc` (default .TRUE.).** Detect identity rows (Dirichlet or redundant periodic) and drop their columns from other rows. This restores the symmetry that PCG and AMG theory assume. It is safe because OFT's BC objects already zero those residual entries.
- **Vectors.** Use `get_slice` / `restore_slice` to move owned values into the backend and back, followed by `stitch` to refresh shared copies. This works for any `oft_vector`.
- **Local mode.** When `A%nr == A%nrg`, for example inside `oft_bjprecond` blocks, use a single-rank communicator (`MPI_COMM_SELF`) with no stitching. This gives a per-rank, per-field AMG sub-solver through the existing `block_jacobi nlocal=-1` with no other plumbing.
- **Setup reuse.**
  - `update(new_pattern=.FALSE.)` marks the values as stale. The next `apply` re-gathers values and re-runs setup only every `rebuild_every` updates; between rebuilds it lags the hierarchy, which suits Newton-Krylov in MUG.
  - `new_pattern=.TRUE.` rebuilds the cached pattern.
- **Laziness and `check_thread`.** Setup is lazy, on the first `apply`. This is required because `oft_bjprecond` clones sub-solvers with `ALLOCATE(..., source=pre)` before they are used. `check_thread` returns `.FALSE.` so that HYPRE or AMGX owns the threading.
- **Deferred interface.** Backends implement `backend_setup(csr)`, `backend_update_values(vals)`, `backend_apply(x_owned, b_owned)` and `backend_delete`. AMGX maps cleanly onto this: `AMGX_matrix_upload_all_global` takes the same owned-row CSR with int64 global columns, and `AMGX_matrix_replace_coefficients` is the value refresh.

### 2.2 Native HYPRE backend

- **New `lin_alg/hypre_bridge.c`.** Plain C; its body is wrapped in `#ifdef HAVE_HYPRE`. An opaque handle struct holds the `HYPRE_IJMatrix`, `HYPRE_ParCSRMatrix`, the x/b `HYPRE_IJVector`s, the solver, and the type. Entry points:
  - `oft_hypre_create(handle, fcomm, type)` converts the communicator with `MPI_Comm_f2c`; in a sequential HYPRE build it passes `hypre_MPI_COMM_WORLD`/0.
  - `oft_hypre_set_matrix(handle, ilower, iupper, nrows, rowptr, gcols(int64), vals, values_only)`. On the C side, `gcols` is converted to `HYPRE_BigInt`, with an overflow check unless HYPRE was built with mixed-int.
  - `oft_hypre_boomeramg_set_opts(handle, opts*)` takes a `BIND(C)` struct where `-1` means the HYPRE default (the STRUMPACK convention), plus `oft_hypre_set_dof_func(handle, n, dof_func)`.
  - `oft_hypre_ams_set_aux(handle, G..., x, y, z, dim, beta_zero)`, `oft_hypre_setup(handle)`, `oft_hypre_apply(handle, b, x)` and `oft_hypre_destroy(handle)`.
  - `oft_hypre_init()` / `oft_hypre_finalize()` are called lazily. They are guarded by `HYPRE_Initialized()` so they don't clash with PETSc's own HYPRE initialization. A finalize hook goes in `oft_finalize` (`base/oft_base.F90:365`).
- **GPU readiness (written now, enabled later).**
  - At init, when `HYPRE_USING_GPU` is defined, call `HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE)` and `HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE)`, and select the device from the node-local rank before `HYPRE_Initialize`.
  - Assemble with `HYPRE_IJMatrixInitialize_v2(ij, HYPRE_MEMORY_HOST)` and then `HYPRE_IJMatrixMigrate(ij, HYPRE_MEMORY_DEVICE)`, so the Fortran side always passes host arrays.
  - Validate the requested options against the GPU-supported subset (§3.3): warn, then remap.
- **API compatibility.** Avoid HYPRE 3.x-only APIs. The bridge must also compile against whichever hypre version the pinned PETSc (3.20-3.25) downloads; check this for each PETSc version.
- **New Fortran module `lin_alg/oft_hypre.F90`** (or a section of `oft_amg.F90`):
  - `oft_hypre_boomeramg` and `oft_hypre_ams`, each with `setup_from_xml` and `view`.
  - The `BIND(C)` interfaces, modelled on the STRUMPACK block in `oft_lu.F90:290`.
- **AMS auxiliary data.** AMS cannot be built from XML alone, because the discretization must supply G and the coordinates.
  - Add `oft_hcurl_get_gmat(hcurl_fe, lag_fe/h1_fe, G)` in `fem/hcurl_operators.F90`. It assembles the matrix-free `hcurl_grad` stencil as an `oft_matrix` from vertex DOFs to H(curl), with zero rows for face/cell DOFs and **no BCs applied**, as HYPRE requires.
  - Add `hcurl_get_vertex_coords(...)`, which returns three vertex-space `oft_vector`s from `mesh%r`.
  - Add `amg_attach_ams_data(solver, G, coords)` in `solver_utils.F90`. It walks a solver tree (`%pre`, and each `oft_ml_precond` level) and attaches the per-level data to any `oft_hypre_ams` it finds.
  - `hcurl_getwop_pre`, `hcurl_getjmlb_pre`, Marklin and MUG call it once they have created the preconditioner from XML.

### 2.3 PETSc backend: `oft_petsc_amgprecond`

This is a small addition because PETSc already wraps HYPRE.

- `setup(pc,dist)` does the following:
  - Calls `PCSetType(pc, PCHYPRE)` and `PCHYPRESetType(pc,'boomeramg'|'ams')`, or `PCSetType(pc, PCGAMG)`; later `PCAMGX`.
  - Translates the OFT-canonical options (§3) into `-<prefix>pc_hypre_boomeramg_*` entries with `PetscOptionsSetValue`. The prefix is unique per object (for example `oft_amg<N>_`) so that the global options database is not polluted.
  - Then calls `PCSetOptionsPrefix` and `PCSetFromOptions`.
  - An optional `<petsc_options>` string is an escape hatch for options OFT does not translate.
- AMS needs extra calls:
  - `PCHYPRESetDiscreteGradient(pc, G%M)`, with G created as an `oft_petsc_matrix` using the same `oft_hcurl_get_gmat`;
  - `PCSetCoordinates(pc, 3, nloc, coords)`;
  - `PCHYPRESetBetaPoissonMatrix(pc, NULL)` when β = 0. Check that PETSc treats NULL as "beta is zero".
- GAMG comes almost for free and is a useful baseline for vector problems. It needs `MatSetNearNullSpace` for velocity rigid-body modes, which can be built from `oft_lag_npos`.

### 2.4 Factory and XML wiring (`solver_utils.F90`)

- Add `"amg"` and `"ams"` to `create_pre_xml`, `create_solver_xml` (so they work as `<coarse>` and block-Jacobi sub-solvers), `create_native_pre` and `create_petsc_pre`.
- Backend choice follows the `lu` precedent:
  - `use_petsc` selects `oft_petsc_amgprecond`;
  - otherwise the native class is used, chosen by `<package>`;
  - `native="true"` on a `<solver>` still forces the native path.
- If the package is not compiled in, `oft_abort("... requires HYPRE")`.
- Add a programmatic helper `create_amg_pre(pre, package)`, like `create_diag_pre`/`create_bjacobi_pre`, for physics codes that do not read XML (TokaMaker, ThinCurr).

---

## 3. XML input specification

### 3.1 BoomerAMG / generic AMG: `<pre type="amg">` (also valid as `<solver type="amg">`)

```xml
<solver type="cg">
  <its>-3</its>
  <pre type="amg">
    <package>hypre</package>          <!-- hypre | gamg (PETSc only) | amgx (future) -->
    <cycles>1</cycles>                <!-- AMG cycles per preconditioner application -->
    <symmetrize_bc>true</symmetrize_bc>
    <rebuild_every>1</rebuild_every>  <!-- re-run setup every N matrix updates (lagging) -->
    <systems>none</systems>           <!-- none | unknown | nodal (multi-field operators) -->
    <ordering>blocked</ordering>      <!-- blocked | interleaved (needs identical field maps) -->
    <verbose>false</verbose>
    <hypre>
      <coarsen_type>hmis</coarsen_type>      <!-- falgout | pmis | hmis | int -->
      <interp_type>ext+i</interp_type>       <!-- classical | direct | ext+i | ext | ext+e | ext+i-mm | int -->
      <strong_threshold>0.5</strong_threshold>
      <max_row_sum>0.9</max_row_sum>
      <pmax>4</pmax>                         <!-- max interp elements/row -->
      <trunc_factor>0.0</trunc_factor>
      <agg_levels>1</agg_levels>             <!-- aggressive coarsening levels -->
      <agg_interp_type>multipass</agg_interp_type>
      <relax_type>l1-gs</relax_type>         <!-- jacobi | l1-jacobi | hgs | shgs | l1-gs | 2s-gs | cheby | int -->
      <relax_sweeps>1</relax_sweeps>
      <relax_weight>1.0</relax_weight>
      <cycle_type>v</cycle_type>             <!-- v | w -->
      <max_levels>25</max_levels>
      <max_coarse_size>-1</max_coarse_size>
      <num_functions>-1</num_functions>      <!-- -1: number of OFT field blocks when systems!=none -->
      <print_level>0</print_level>
    </hypre>
    <petsc_options>-pc_hypre_boomeramg_print_statistics 1</petsc_options>  <!-- PETSc escape hatch -->
  </pre>
</solver>
```

| Group | Description | Type | Notes |
|---|---|---|---|
| `package` | AMG package | str | Default: `hypre` if compiled in. `gamg` requires PETSc. |
| `cycles` | Cycles per application (HYPRE `MaxIter`, with `Tol=0`) | int | Default 1 |
| `symmetrize_bc` | Remove Dirichlet columns before handing the matrix to AMG | bool | Default true |
| `rebuild_every` | Hierarchy rebuild interval on value-only updates | int | Default 1 |
| `systems` | Multi-field treatment | str | `unknown` uses `dof_func` (any ordering). `nodal` needs `interleaved`. |
| `ordering` | Global numbering of fields | str | Default `blocked` |
| `hypre/*` | BoomerAMG parameters | – | `-1` or omitted means the HYPRE default. Each accepts a name or the raw HYPRE integer. |

**Canonical names and their translations.** One table maps each OFT name to both the native HYPRE integer code and the PETSc option string, so XML files are portable between backends.

| OFT name | HYPRE code | PETSc option |
|---|---|---|
| `falgout` / `pmis` / `hmis` | 6 / 8 / 10 | `-pc_hypre_boomeramg_coarsen_type Falgout\|PMIS\|HMIS` |
| `classical` / `direct` / `ext+i` / `ext` / `ext+e` | 0 / 3 / 6 / 14 / 18 | `-pc_hypre_boomeramg_interp_type classical\|direct\|ext+i\|ext\|ext+e-mm` |
| `jacobi` / `hgs` / `shgs` / `l1-gs` / `cheby` / `l1-jacobi` | 0 / 3 / 6 / 8 (13/14 directional) / 16 / 18 | `-pc_hypre_boomeramg_relax_type_all Jacobi\|SOR/Jacobi\|symmetric-SOR/Jacobi\|l1scaled-SOR/Jacobi\|Chebyshev\|l1scaled-Jacobi` |
| `pmax`, `agg_levels`, `strong_threshold`, `max_levels`, `num_functions`, `nodal`, `cycles` | – | `_P_max`, `_agg_nl`, `_strong_threshold`, `_max_levels`, `_numfunctions`, `_nodal_coarsen`, `_max_iter` |

Before implementing, verify the exact PETSc spellings against the documentation for the pinned PETSc version.

### 3.2 AMS: `<pre type="ams">`

```xml
<pre type="ams">
  <package>hypre</package>
  <cycles>1</cycles>
  <beta_zero>auto</beta_zero>      <!-- auto | true | false ; auto = true for WOP (no mass term) -->
  <symmetrize_bc>true</symmetrize_bc>
  <hypre>
    <cycle_type>1</cycle_type>     <!-- 1 (default) or 13 recommended; additive 2,4,6,8,12,14 OK as preconditioner only -->
    <relax_type>2</relax_type>     <!-- AMS fine-grid smoother (l1 options) -->
    <relax_sweeps>1</relax_sweeps>
    <relax_weight>1.0</relax_weight>
    <omega>1.0</omega>
    <alpha_amg> <coarsen_type>hmis</coarsen_type> <agg_levels>1</agg_levels>
                <relax_type>shgs</relax_type> <strong_threshold>0.25</strong_threshold>
                <interp_type>classical</interp_type> <pmax>0</pmax> </alpha_amg>
    <beta_amg>  ...same keys... </beta_amg>
    <print_level>0</print_level>
  </hypre>
</pre>
```

- The calling FE code must supply the discrete gradient and the vertex coordinates through `amg_attach_ams_data`. If they are missing at the first `apply`, abort with a clear message.
- **Marklin WOP (p-MG with AMS at the coarse level):**
  `<hcurl><wop><pre type="mg">…<coarse><solver type="ams">…</solver></coarse></pre></wop></hcurl>`
- **MUG:** a per-field block-Jacobi (`nlocal=-1`) with an `amg` sub-solver works right away, but only rank-locally. A true global field split is the Phase 3 item below.

### 3.3 GPU-safe defaults (documented, and enforced when HYPRE is built with CUDA)

| Setting | GPU-supported values |
|---|---|
| Coarsening | PMIS (8) |
| Interpolation | 3, 6, 14, 15, 18 |
| Relaxation | 3, 4, 6, 7, 11, 12, 16, 18 |
| Aggressive coarsening interpolation | 5, 6, 7 |
| Relaxation order | 0 (no C/F ordering) |

On CPU, OFT uses HYPRE's own defaults: HMIS with ext+i, `pmax` 4, and l1-Gauss-Seidel. Recommended starting point for 3D and higher-order operators: `strong_threshold` 0.5-0.6 and one level of aggressive coarsening.

### 3.4 AMGX (future): `<pre type="amg"><package>amgx</package>`

```xml
<amgx>
  <config_file>amgx_classical_pmis.json</config_file>   <!-- AMGX JSON config -->
  <config>max_iters=1, selector=PMIS, smoother=MULTICOLOR_GS</config>  <!-- AMGX "key=value" overrides -->
  <mode>dDDI</mode>
</amgx>
```

The shared keys (`cycles`, `systems`, `ordering`, `rebuild_every`, `symmetrize_bc`) carry over unchanged. Mapping:

| OFT key | AMGX equivalent |
|---|---|
| `systems=nodal` | `block_dimx/y` |
| `rebuild_every` | `AMGX_matrix_replace_coefficients` + `AMGX_solver_resetup`, or `structure_reuse_levels` |

---

## 4. Building HYPRE: `src/utilities/build_libs.py`

### 4.1 New `HYPRE(package)` class

Model it on `STRUMPACK` (L1748), and on `SUPERLU_DIST` (L1578) for the MPI handling.

- **Source.**
  - `url`: `https://github.com/hypre-space/hypre/archive/refs/tags/v3.2.1.tar.gz` (the latest tag on the GitHub releases page). Pin it.
  - `build_dir`: `hypre-3.2.1`. HYPRE's CMake project lives in `src/`, so configure with `{CMAKE} ... ../src`.
  - HYPRE requires CMake ≥ 3.21; check this the same way STRUMPACK does.
- **`setup()`.**
  - Call `setup_root_struct()` and set `config_dict['HYPRE_LIBS'] = '-lHYPRE'`.
  - `install_chk_files = [HYPRE_LIB/libHYPRE.a]`, or the `DYN_EXT` file if built shared.
- **CMake options.**

```
-DCMAKE_INSTALL_PREFIX:PATH={HYPRE_ROOT} -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_BUILD_TYPE=Release
-DBUILD_SHARED_LIBS:BOOL=OFF -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON   # static+PIC (Python shared libs, packaging)
-DHYPRE_ENABLE_MPI:BOOL=ON|OFF            # follows use_mpi; OFF => sequential HYPRE for OpenMP-only OFT builds
-DMPI_C_COMPILER={MPI_CC} -DMPI_CXX_COMPILER={MPI_CXX}   # when wrappers given; MPI_HOME for lib/include style
-DHYPRE_ENABLE_OPENMP:BOOL=ON|OFF         # follows --oft_use_openmp (as STRUMPACK); forced OFF with --hypre_cuda
-DHYPRE_ENABLE_BIGINT:BOOL=OFF            # never: OFT local indices are i4
-DHYPRE_ENABLE_MIXEDINT:BOOL=OFF|ON       # --hypre_mixedint: 64-bit global ids (>2^31 DOFs), 32-bit local
-DHYPRE_ENABLE_COMPLEX:BOOL=OFF -DHYPRE_ENABLE_SINGLE:BOOL=OFF -DHYPRE_ENABLE_MIXED_PRECISION:BOOL=OFF
-DHYPRE_ENABLE_FORTRAN:BOOL=OFF           # OFT talks through its own C bridge
-DHYPRE_ENABLE_HYPRE_BLAS:BOOL=ON -DHYPRE_ENABLE_HYPRE_LAPACK:BOOL=ON   # internal (hypre_-prefixed, no clashes);
                                          #   optional later: OFF + TPL_BLAS_LIBRARIES={BLAS_LIB_PATH}, TPL_LAPACK_LIBRARIES
-DHYPRE_ENABLE_SUPERLU:BOOL=OFF -DHYPRE_ENABLE_DSUPERLU:BOOL=OFF -DHYPRE_ENABLE_MAGMA:BOOL=OFF
-DHYPRE_ENABLE_UMPIRE:BOOL=OFF -DHYPRE_ENABLE_CALIPER:BOOL=OFF
-DHYPRE_BUILD_TESTS:BOOL=OFF -DHYPRE_BUILD_EXAMPLES:BOOL=OFF
-DCMAKE_OSX_SYSROOT={MACOS_SDK_PATH}      # when set
# GPU (future, --hypre_cuda):
-DHYPRE_ENABLE_CUDA:BOOL=ON -DCMAKE_CUDA_ARCHITECTURES={CUDA_ARCH}
-DHYPRE_ENABLE_UNIFIED_MEMORY:BOOL=OFF -DHYPRE_ENABLE_GPU_AWARE_MPI:BOOL=OFF|ON
(cuSPARSE/cuBLAS/cuRAND/cuSOLVER default ON)
```

- **Build lines.** Use the STRUMPACK sequence: `export CC/CXX/FC` → `{CMAKE} ... ../src` → `make -j{MAKE_THREADS}` → `make install`.
- **MPI without wrappers** (the `MPI_LIBS`/`MPI_INCLUDE` case): follow ARPACK/HDF5. Note that `SUPERLU_DIST` does not handle this case today.

### 4.2 Command-line flags (new `"HYPRE"` argparse group)

| Flag | Type | Meaning |
|---|---|---|
| `--build_hypre` | int {0,1} | Build HYPRE. With `--build_petsc=1` this means PETSc `--download-hypre` instead (see below). |
| `--hypre_mixedint` | store_true | Use 64-bit global indices |
| `--hypre_cuda` | store_true | **Future.** Build HYPRE with CUDA. Requires `--cuda_arch`, and rejected together with PETSc until PETSc is built with `--with-cuda`. |
| `--cuda_arch` | str | **Future.** Shared by HYPRE and AMGX, e.g. `80;90` |
| `--build_amgx` | int {0,1} | **Future.** Build NVIDIA AMGX v2.5.0 (CUDA ≥ 12, SM70+, CMake). Uses `CMAKE_CUDA_ARCHITECTURES`; needs MPI for the distributed version. |

### 4.3 Registration and guards (`L2194-2286`)

- **Without PETSc:** append `HYPRE(...)` after STRUMPACK. HYPRE works with or without MPI, so unlike SuperLU-DIST there is **no MPI requirement**.
- **With PETSc:** add a `with_hypre` argument to `PETSC(...)`:
  - it appends `--download-hypre`, and sets `need_cxx=True` if PETSc's hypre recipe needs C++ (check for the pinned PETSc version);
  - the native bridge then links PETSc's libHYPRE, so there are never two copies;
  - error out on `--hypre_mixedint` unless PETSc is also built with 64-bit indices, or simply disallow the combination.
- **`build_cmake_script` (L520):** add `if "HYPRE_ROOT" in mydict: -DOFT_HYPRE_ROOT:PATH=...` next to STRUMPACK, outside the PETSc `else`.
- **Separate fix found during the survey:** in `PETSC.setup`, `install_chk_files` is appended to and then overwritten (L1963-1975), and the appends also mutate the class-level list.

### 4.4 CMake

- **New `cmake/FindOFT_HYPRE.cmake`** (based on `FindOFT_STRUMPACK.cmake`):
  - `find_library(HYPRE)` and `find_path(HYPRE.h)`, reading `HYPRE_config.h` directly rather than HYPRE's CMake package config;
  - **fatal** on `HYPRE_BIGINT` and `HYPRE_COMPLEX`;
  - **fatal** on a `HYPRE_SEQUENTIAL` versus `OFT_USE_MPI` mismatch, and on `HYPRE_USING_OPENMP` without `OFT_USE_OpenMP`;
  - when `HYPRE_MIXEDINT` is set, define `OFT_HYPRE_MIXEDINT`;
  - when `HYPRE_USING_CUDA` is set, call `find_package(CUDAToolkit)` and add `CUDA::cudart cusparse cublas curand cusolver` to `OFT_HYPRE_LIBRARIES`.
- **`CMakeLists.txt` (L400-473):** add a STRUMPACK-style block: `add_definitions(-DHAVE_HYPRE)`, `include_directories`, and prepend `EXT_LIBS`.
- **`FindOFT_PETSc.cmake`:** detect `libHYPRE` in PETSc's lib directory. When found, set `OFT_HYPRE_FOUND` and add `-DHAVE_HYPRE -DHAVE_PETSC_HYPRE` (the PETSc class needs this to enable `PCHYPRE`).
- **`lin_alg/CMakeLists.txt`:** add `oft_amg.F90`, `oft_hypre.F90` and `hypre_bridge.c`.
- **Packaging:** nothing needed while HYPRE is static. If it is ever shared, or once CUDA/AMGX runtimes exist, add them to `filter_list` in `cmake/patch_package.py`.
- **Spack:** add a `+hypre` variant: `depends_on("hypre~int64~complex+mpi", when="+hypre+mpi")`, the `~mpi` equivalent, and `depends_on("petsc+hypre", when="+hypre+petsc")`.

---

## 5. Tests (`src/tests`)

- **`tests/lin_alg/test_solver_xml.py`.** Add `test_cg_amg[petsc_flag]` (BoomerAMG on the existing 4th-order Lagrange Laplacian), `test_cg_amg_coarse` (p-MG with an AMG `<coarse>` level), and `test_cg_gamg[petsc]`.
  - `oft_testing.run_OFT` needs a mapping from "requires HYPRE" to `pytest.skip`, matching the existing PETSc case.
- **New `tests/lin_alg/test_hypre.F90/.py`** (skip with `SKIP TEST` when `HAVE_HYPRE` is absent):
  1. p = 1 Lagrange Laplacian on refined cube levels. Check that iterations stay roughly constant, i.e. AMG scalability.
  2. 1, 2 and 4 MPI ranks give the same iteration counts within tolerance. This validates the global-numbering and ghost-row path.
  3. A periodic mesh, which exercises the redundant-row handling. Periodic native-mesh tests now exist.
  4. A value-only `update(.FALSE.)` refresh, and lagging with `rebuild_every`.
  5. Vector Lagrange with `systems=unknown` versus `nodal` + `interleaved`.
  6. Block-Jacobi `nlocal=-1` with an `amg` sub-solver.
- **`tests/fem/test_hcurl.py` and `tests/physics/Marklin/test_taylor.py`:** AMS on WOP, p = 1 direct and p = 2 as the coarse level. Compare iteration counts against the current Jacobi-smoothed MG.
- **CI** (`.github/workflows/ci_build.yaml`):
  - Add `--build_hypre=1` to both `build_flags` and `mpi_build_flags` of the Ubuntu 24.04 GCC 14 entry. This covers sequential and MPI HYPRE.
  - Add `--build_hypre=1` (via PETSc) to the `mpich_petsc` include and to `cov_build.yaml`.
  - `copilot-setup-steps.yml` can stay as it is.

---

## 6. Phasing

| Phase | Scope | Rough size |
|---|---|---|
| **0: Build** | `HYPRE` class and flags, PETSc `--download-hypre`, `FindOFT_HYPRE`, `HAVE_HYPRE`, CI flags, init/finalize hooks | Small to medium |
| **1: BoomerAMG** | `oft_amg_solver` (numbering, CSR extraction, BC symmetrization, reuse), HYPRE bridge, `oft_hypre_boomeramg`, `oft_petsc_amgprecond` (hypre and gamg), factory and XML wiring, `test_solver_xml` plus `test_hypre` cases 1-4, a `doc_solvers.md` "Algebraic Multigrid" section | Medium. Most of the work is the conversion layer. |
| **2: Use in physics** | AMG as `<coarse>` in Lagrange/H1 MG; `create_amg_pre` for MUG divergence cleaning and the ThinCurr R⁻¹ fallback; systems AMG (cases 5-6) | Small to medium |
| **3: AMS and field split** | `oft_hcurl_get_gmat`, vertex coordinates, `amg_attach_ams_data`, `oft_hypre_ams` and PETSc AMS; Marklin WOP. MUG: a global field-split preconditioner (native sub-matrix extraction across ranks, or exposing `PCFIELDSPLIT` on the PETSc side) with AMG/AMS per field | Medium to large |
| **4: GPU** | `--hypre_cuda` (device selection, IJ migrate, enforcing GPU-supported options), then the `oft_amgx_solver` backend and `--build_amgx` / PETSc `PCAMGX`, all on the shared CSR layer | Medium |

### Documentation to update with the code

- `docs/doc_solvers.md`: an AMG/AMS section with the tables from §3.
- `docs/doc_petsc_usage.md`: `<petsc_options>`.
- `.agents/skills/oft-build/SKILL.md`: the flag table.
- `.agents/skills/oft-test/SKILL.md`: the skip reasons.
- `AGENTS.md`: the `lin_alg` library list.

### Open questions

1. Should `amg` without `<package>` default to HYPRE, or to GAMG when `use_petsc` is set and HYPRE is absent?
2. Is a lagged hierarchy (`rebuild_every > 1`) acceptable as the default for MUG Newton-Krylov runs, or should every Jacobian update rebuild?
3. Is the Phase 3 global field split worth doing natively, or should MUG large runs simply rely on the PETSc backend plus `PCFIELDSPLIT`?
