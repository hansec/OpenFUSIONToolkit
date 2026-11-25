# Structure-aware smoothers for MUG p-multigrid: study report

Branch `claude/pmg-structure-aware-smoothers`, implementation commit `2c4919e`, based on `5b4252d`.
Status: **development stopped**. This report records what was built, what was measured (including negative
results), and the open problems, for future work.

---

## 1. Summary

* **Problem.** MUG's p-multigrid V-cycle smooths with a few FGMRES iterations preconditioned by block-Jacobi with
  local LU (`oft_bjprecond`). As a stationary iteration that smoother diverges on MUG Jacobians. It fails for two
  separate reasons, one at small Alfvén CFL and one at large CFL.
* **What was built.** A new native preconditioner/smoother `type="patch"` (`oft_patch_precond`,
  `src/lin_alg/oft_patch_solvers.F90`):
  * overlapping patches solved exactly with dense LU;
  * additive, multiplicative and symmetric modes;
  * opt-in single-precision factors;
  * patches on process boundaries assembled and solved in full across MPI ranks.

  A generic builder (`comp_fem_vertex_patches`, `src/fem/fem_utils.F90`) attaches **vertex-star** patches to the
  matrix: every DOF of every field on a vertex and on the edges, faces and cells containing it. MUG (`xmhd` and
  `xmhd_lag`) attaches them on every multigrid level.
* **Headline results.**
  * **Small CFL (≲ 0.1): stars work as theory predicts.**
    * Stationary convergence factor 0.175 per symmetric sweep (two-grid 0.013). Block-Jacobi diverges at 1.8–2.2.
    * Same for both B representations.
    * In MUG: about 2 outer iterations per step versus 3 for block-Jacobi with overlap, about 4× faster.
    * Insensitive to rank count (2–3 iterations on 1–4 ranks). Block-Jacobi+overlap degrades to 8 and 14.
  * **CFL ≳ 1: stars alone fail (negative result).**
    * The coupled B–V star blocks become nearly singular (scaled condition about 4e8).
    * Stationary factor 1.6 at dt = 1e-5; up to 2.5e6 at dt = 1e-4.
    * The p-coarse level does not rescue it.
  * **A hybrid helps at CFL ≈ 1.** Star down smoother plus overlapped block-Jacobi up smoother: 4 iterations
    serial and 9–10 on 2/4 ranks, against 6 and 29–31 for block-Jacobi+overlap alone (dt = 1e-5 cube).
  * **Spheromak at CFL ≈ 20.** No configuration was satisfactory. Stars run in 2.6 GB but need dt cut 4×. Every
    LU-based block-Jacobi variant exceeded the 14 GB sandbox. ILU variants need 52–122 iterations per step.
* **Open problem.** Robust smoothing for CFL ≳ 1 needs subdomains or a preconditioner structure that captures
  Alfvén-wave coupling over more than one vertex star (Section 11).

---

## 2. Background and problem statement

MUG (`src/physics/xmhd.F90`, `src/physics/xmhd_lag.F90`):

* **Time stepping:** implicit, with Newton–Krylov and an outer FGMRES.
* **Preconditioner:** a V-cycle (`oft_ml_precond`, `native_solvers.F90`; built in `create_ml_xml`,
  `solver_utils.F90`) on an assembled approximate Jacobian, rebuilt per level by `xmhd_set_ops`/`xmhd_build_ops`.
* **Levels:** polynomial orders p, …, 1 on the finest mesh, optionally with geometric levels below.
* **Fields:**
  * `xmhd`: B in H(curl)_p (Bc) + grad H1_{p+1} (Bg, vertex-gradient DOFs constrained to identity rows); V as 3
    Lagrange blocks; then n, T.
  * `xmhd_lag`: all Lagrange (B with divergence cleaning).

### Two failure mechanisms

1. **Small CFL: B-basis conditioning.**
   * The hierarchical H(curl) + grad(H1) basis is badly conditioned. Block-Jacobi on the B mass block alone has a
     stationary factor of 1.07 / 1.48 / 1.60 for p = 1 / 2 / 3.
   * Non-overlapping subdomains do not contain the gradients of vertex hat functions, so the near-kernel of curl is
     not resolved locally. That is exactly the Lee–Wu–Xu–Zikatanov condition (Section 3).
2. **Large CFL: B–V Alfvén coupling.** The Schur coupling scales like CFL². Errors that couple B and V over
   wavelengths longer than a subdomain are not damped by local solves.

### Baselines before this work (measured by the parent session)

| Case | One-level | p-MG |
|---|---|---|
| Spheromak, dt = 1e-6, block-Jacobi | 74 | 39–67 |
| Spheromak, dt = 1e-6, block-Jacobi + `<loverlap>` (one graph layer, RAS-like) | 14 | 11–13 |
| Spheromak, dt = 1e-5 | nothing converged | nothing converged |

The parent session also found that MPI interfaces cost a factor of 5–10. On the refined non-periodic Alfvén cube
(p = 3, about 40k DOF, dt = 1e-5, 16 blocks), outer iterations at step 1 were:

| Configuration | One-level | p-MG |
|---|---|---|
| serial, `loverlap` | 13 | 12 |
| 2 ranks, `loverlap` | 114 | 66 |
| 2 ranks, `loverlap` + `boverlap` | 84 | 42 |
| 4 ranks, `loverlap` | 164 | 75 |

This finding drove the cross-rank seam-patch design (Section 4.4).

---

## 3. Literature and other codes

Full texts on arxiv.org and the Tufts author pages could not be fetched (network policy). Items marked
*(unverified)* rest on abstracts or memory and should be checked.

**Proven theory**

* **Arnold, Falk & Winther, Numer. Math. 85 (2000).** For (u,v) + τ(curl u, curl v) in H(curl), multigrid with
  Schwarz smoothers on **vertex-star patches** converges uniformly in h and τ. Stars are needed because they
  contain grad φ_v.
* **Hiptmair, SINUM (1998).** A hybrid smoother (edge Gauss–Seidel plus Gauss–Seidel on the nodal potential via the
  discrete gradient) is h-uniform for lowest-order curl-curl + mass.
* **Lee, Wu, Xu & Zikatanov, M3AS (2007).** For A0 + εM, a subspace correction is ε-robust iff
  ker A0 = Σ_i (V_i ∩ ker A0).
  * For MUG's B space, every kernel function is either the gradient of a vertex hat (in the edge DOFs of star(v)) or
    an explicit Bg DOF on an edge, face or cell.
  * So **vertex stars satisfy the condition; METIS blocks and one-layer graph overlap do not.**
* **Pavarino (1994); Schöberl, Melenk, Pechstein & Zaglmayr (2008).** Overlapping vertex-star Schwarz is p-robust
  for H1. The H(curl) p-robust result (Zaglmayr thesis 2006) relies on an exact-sequence basis like MUG's
  *(unverified detail)*.

**Empirical (MHD and solver frameworks)**

* **Adler, Benson, Cyr, Farrell, MacLachlan & Tuminaro, SISC 43(5) (2021), "Monolithic multigrid methods for
  magnetohydrodynamics"** (arXiv 2006.15700). Vanka-type monolithic MG for resistive MHD, up to 1.7e8 DOF; exact
  patch definitions *(unverified)*.
* **Laakmann, Farrell & Mitchell, SISC 44(4) (2022)** (arXiv 2104.14855). Augmented-Lagrangian outer solver with
  star-patch MG on the (E,B) block, implemented in Firedrake + PCPATCH. Note that they wrap the patches in a
  physics-based outer block structure rather than using monolithic stars alone, which is consistent with our
  large-CFL finding.
* **Farrell, Knepley, Mitchell & Wechsung, ACM TOMS (2021), PCPATCH.** Topological star and Vanka patches,
  additive or multiplicative, dense patch LU. This was the design template here.
* **Brubeck & Farrell (2024).** Star patches with p-MG for the de Rham complex. Dense patch solves at high p are the
  cost bottleneck; they sparsify with FDM bases.
* **Shadid et al. (2016), Drekar.** Fully coupled AMG for nodal equal-order MHD with nodal-block and additive
  Schwarz/ILU smoothers that keep all fields at a node together. Relevant to `xmhd_lag`.
* **hypre AMS / MFEM.** Auxiliary-space solver for the H(curl) block only. It is not in OFT's native stack, and
  using it would require a physics-based block (Schur) preconditioner.

**Inference made at plan time, and what happened**

* Expected: stars fix the small-CFL divergence. **Confirmed.**
* Expected: a big improvement at moderate CFL, but no CFL-robustness from smoothing alone. **Worse than expected:**
  monolithic stars fail outright at CFL ≈ 1 rather than degrading gracefully.

---

## 4. Design

### 4.1 Patch definition

For each vertex v, the patch holds all DOFs of all fields on:

* v itself;
* the edges containing v;
* the faces of the star cells that contain v;
* the cells of star(v) (`mesh%kpc/lpc`).

DOF numbering is entity-ordered by `gstruct = (g1,g2,g3,g4)`:

| Entity | Local DOF index |
|---|---|
| vertex | `(v-1)*g1+j` |
| edge | `g1*np+(e-1)*g2+j` |
| face | `g1*np+g2*ne+(f-1)*g3+j` |
| cell | `g1*np+g2*ne+g3*nf+(c-1)*g4+j` |

Composite fields are concatenated with offsets equal to the sum of the preceding `fe%ne`. Patch lists are sorted
and unique. The builder is generic over any `oft_fem_comp_type`, so it serves `xmhd`, `xmhd_lag` and single-field
Lagrange (used in the unit test).

Typical tetrahedral patch sizes (B+V+n+T): about 35 DOF at p = 1, 230 at p = 2, 770 at p = 3. The measured mean/max
was 196/316 on the p = 2 periodic cube and 636/1044 at p = 3 (H(curl) B); for Lagrange B at p = 3 it was 390/630.

### 4.2 Periodic meshes

Local vertices are grouped by `mesh%global%lp`, so periodic copies of one vertex form one patch. Two conventions
are needed:

* **Matrix entries:** taken from `Mfull` using non-redundant ("rep") rows and rep columns only.
* **Interpolation P** (offline analysis only): rows of all periodic copies are summed, because `apply` does
  `stitch(1)`.

Both conventions were validated against Fortran mat-vecs to round-off.

### 4.3 Colouring and modes

* **Colouring:** a greedy distance-1 colouring of patches, so that no two same-coloured patches share a cell.
  Same-coloured patches are therefore independent and are processed in parallel with OpenMP.
* **additive:** u += df · Σ_i w_i R_iᵀ A_i⁻¹ R_i r, with w = 1/(global patch count of each DOF). The count is made
  consistent across ranks.
* **multiplicative:** a coloured sweep over interior patches with local residual updates. Then one global residual
  (`A%apply`), then an additive correction from the seam patches.
* **symmetric:** multiplicative followed by a backward interior sweep.
* **Semantics:** when `its > 0` the object performs `its` damped Richardson iterations, accepts a nonzero initial u
  and returns the residual in g, so it can be a multigrid smoother directly. When `its <= 0` it acts as a
  preconditioner, for example under GMRES (the configuration used in all MUG runs).

### 4.4 Cross-rank seam patches

A patch whose vertex is shared between ranks is a seam patch:

* **Ownership.** The owner is the minimum rank among the vertex sharers (from `mesh%pstitch`: `proc_con`, `kle/lle`,
  `lbp`).
* **Assembly.** Each contributor sends the global DOF ids (i8, from `i_map%lge`/offsets) of its part of the star,
  plus the dense block of its **unassembled** local matrix `M` on those DOFs.
  * Every cell containing a pair of star DOFs contains the vertex, and each cell lives on exactly one rank. So the
    sum over contributors is exactly the assembled global block, including Dirichlet rows, which only the owner
    sets.
  * The owner sums the blocks into a global-id-indexed block and factors it.
* **Apply.** Contributors send residual values. The owner solves and returns each contributor's part of the
  correction. Every holder of a DOF sees every patch containing it, so counts and corrections are consistent. Shared
  rows are combined with owner-wins `restore_local`.
* **Communication.** Point-to-point (`MPI_IRECV`/`MPI_ISEND`, tag 31) over neighbour lists built once in
  `patch_precond_structure`, with send and receive lists ordered by a global key.
* **Consequence.** Additive mode gives the same result on any number of ranks (verified: identical GMRES iteration
  counts on 1, 2 and 4 ranks). The multiplicative modes treat seam patches additively, which weakens them as the
  number of ranks grows (Section 7.4).

### 4.5 Storage

* Dense LU per patch (`dgetrf/dgetrs`). Memory is Σ n_i².
* `single = T` stores the factors in r4 (`sgetrf/sgetrs`), halving memory, with the residual kept in r8.
* Estimates: about 0.4 MB per patch at p = 2 and about 4.7 MB at p = 3. Measured: 2.6 GB peak for the Spheromak
  p = 2 run (about 2.85e5 DOF), whole process.

---

## 5. Implementation map (commit 2c4919e, 17 files, +1882 lines)

| File | Change |
|---|---|
| `src/lin_alg/oft_patch_solvers.F90` (new, 1237 lines) | Module `oft_patch_solvers`, type `oft_patch_precond` (extends `oft_solver`). See the routine list below. |
| `src/lin_alg/lin_alg_base.F90` | New public type `oft_patch_info` (`npatch`, `ncolors`, `patches(:)`, `color(:)`, `gid(:)` i8 global vertex key, `owner(:)`, `share(:)`). New field `oft_matrix%patch_info`. New `matrix_delete_patches`. |
| `src/lin_alg/native_la.F90` | `mat_delete` frees `patch_info`. |
| `src/lin_alg/solver_utils.F90` | `CASE("patch")` in `create_solver_xml`, `create_native_pre` and `create_pre_xml`. Aborts under PETSc. |
| `src/lin_alg/CMakeLists.txt` | Adds `oft_patch_solvers.F90` to `SOLVER_OBJS` (after `oft_lu.F90`). |
| `src/fem/fem_utils.F90` | New `comp_fem_vertex_patches(self_comp,pinfo)`: star construction, DOF formulas, periodic merge, colouring, sharers/owner. Internal helpers take explicit arguments, so they are OpenMP-safe. |
| `src/physics/xmhd.F90`, `xmhd_lag.F90` | In `xmhd_alloc_ops`, after the existing colouring: `CALL comp_fem_vertex_patches(xmhd_rep,ops%J%patch_info)` on every level. This costs only integer memory. |
| `src/docs/doc_solvers.md` | New section `doc_solvers_basic_patch` (XML spec, memory note, CFL limitation, hybrid). Adds the missing `loverlap` row to the block_jacobi table. |
| `src/tests/lin_alg/test_native_patch.{F90,py}` + CMake | Unit test (Section 9). |
| `src/tests/physics/MUG/oft_in_xmhd_mg_patch.xml`, `mug_testing.py` (`mg_xml` argument), `test_alfven.py`, `test_alfven_lag.py`, CMake | MUG regression tests with the patch smoother in p-MG. |

**Routines in `oft_patch_solvers.F90`:**

* `nbr_exchange`, `nbr_exchange_sized`, `free_bufs`: neighbour communication.
* `patch_precond_structure`: builds the structure.
  * rep map via sorted `lge`;
  * global isolated-row flags via a vector restore-add;
  * local and seam patch lists, neighbour lists, owned-patch gid unions and positions;
  * colours, overlap check, weights.
* `free_structure`.
* `patch_precond_factor`: interior blocks from `Mfull` (`fill_patch`); seam blocks from the unassembled `M`
  (`partial_block`, `partial_block_owned`), exchanged and summed; then LU.
* `patch_solve`, `patch_precond_seam`, `patch_precond_sweep` (with `interior_sweep` and `solve_patch`).
* Solver interface: `patch_precond_apply`, `_update`, `_setup_xml`, `_view`, `_check_thread` (returns `.FALSE.`),
  `_delete`.

**Limitations.**

* Native backend only.
* The matrix must carry `patch_info`. Currently only MUG attaches it, and the unit test attaches it by hand.
* Not thread-safe as a sub-solver; it is threaded internally.
* `update()` refactors on every call.

---

## 6. Usage

**XML options** (in `<pre type="patch">` or `<solver type="patch">`):

| Option | Meaning |
|---|---|
| `mode` | `additive` (default), `multiplicative` or `symmetric` |
| `its` | Richardson iterations, per-level list; `<= 0` (the default, -1) means preconditioner mode |
| `df` | damping factor, per-level list; default 1 |
| `single` | single-precision factors; default false |

**Stars in p-MG** (small CFL; this is the regression-test configuration):

```xml
<oft><xmhd><pre type="mg">
  <smoother direction="up"><solver type="gmres"><its>2</its><nrits>2</nrits>
    <pre type="patch"><mode>symmetric</mode></pre>
  </solver></smoother>
</pre></xmhd></oft>
```

**Hybrid** (CFL ≈ 1: star down smoother, overlapped block-Jacobi up smoother):

```xml
<oft><xmhd><pre type="mg">
  <smoother direction="down"><solver type="gmres"><its>2</its><nrits>2</nrits>
    <pre type="patch"><mode>symmetric</mode></pre>
  </solver></smoother>
  <smoother direction="up"><solver type="gmres"><its>2</its><nrits>2</nrits>
    <pre type="block_jacobi"><nlocal>16</nlocal><loverlap>T</loverlap><solver type="lu"></solver></pre>
  </solver></smoother>
</pre></xmhd></oft>
```

**One-level star preconditioner:**

```xml
<oft><xmhd><pre type="gmres"><its>2</its><nrits>2</nrits>
  <pre type="patch"><mode>symmetric</mode></pre>
</pre></xmhd></oft>
```

Guidance:

* Use `symmetric` for serial runs and small rank counts. Use `additive` where results must be independent of the
  rank count.
* At p = 3 with large meshes consider `single`.
* Do not use stars alone when the Alfvén CFL number is ≳ 1.

---

## 7. Results

All runs used `OMP_NUM_THREADS=2` unless stated otherwise. MPI timings with 4 ranks × 2 threads were
oversubscribed on the sandbox, so wall times at 4 ranks are not meaningful.

### 7.1 Phase 0: offline stationary convergence factors (scipy, dumped MUG operators)

Method:

* Per-level J, P and patch lists were dumped from MUG by a temporary hook (Appendix B).
* They were reduced to non-redundant rows; identity Dirichlet/constraint rows were dropped.
* rho(S) was estimated by power iteration on the error-propagation operator.
* Dumps reproduce Fortran mat-vecs to round-off, and p-level Galerkin products are exact (relative difference
  about 1e-16).

Periodic Alfvén cube. "tiny" is p = 2; "small" is p = 3, 58,880 DOF (H(curl) B) or 41,472 (Lagrange B); "medium"
is p = 2, 167,936 DOF.

| Case | block-Jacobi (8 parts) | RAS, 1 layer | star additive (ω = 0.7) | star mult. | star sym. | two-grid, star sym. (pre+post) |
|---|---|---|---|---|---|---|
| tiny, p = 2, dt = 8e-7, B block only | 1.81 | 0.50 | 0.62 | 0.11 | 0.08 | – |
| tiny, p = 2, dt = 8e-7, all fields | 1.81 | 0.50 | 1.00 (0.66) | 0.26 | 0.26 | 0.04 (bj 1.04, ras1 0.19, star_add 0.45) |
| medium, p = 2, H(curl), dt = 8e-7 | 1.89 | – (out of memory) | – | – | – | – |
| medium, p = 2, H(curl), dt = 1e-4 | – | – | 14.3 (10.0) | – | 4.5e10 | – |
| small, p = 3, H(curl), dt = 8e-7 | 2.20 | 0.47 | 0.62 (0.64) | – | **0.175** | **0.013** (post-only 0.079) |
| small, p = 3, Lagrange, dt = 8e-7 | – | – | 0.60 (0.59) | – | **0.175** | **0.013** (post-only 0.079) |
| small, p = 3, H(curl), dt = 1e-4 | 13.6 | 2.37 | 10.1 (7.1) | – | 2.5e6 | 8.7e12 |
| small, p = 3, Lagrange, dt = 1e-4 | – | – | 5.7 (4.0) | – | 6.9e3 | 4.5e7 |

Notes:

* On the tiny mesh each of the 8 RAS parts is about 1/8 of the domain, so that comparison is not size-matched.
* The coarse-only two-grid factor is 1.000 everywhere: the p-coarse correction alone does not reduce the error, as
  expected.

**Non-periodic p = 3 cube, dt = 1e-5 (CFL ≈ 1): field-block diagnosis**

| Block | n | star sym. rho | scaled patch condition (max / median) |
|---|---|---|---|
| B only | 20,568 | 0.080 | 1.9e3 / 1.7e3 |
| V only | 8,979 | 0.110 | about 8e0 |
| B + V coupled | 29,547 | **1.632** | **4.7e8 / 3.5e8** |

### 7.2 MUG, non-periodic Alfvén cube, serial

Setup:

* `test_alfven` (H(curl) B), p = 3, linear, `ref_per=F`, `nlevels=2`, `nbase=2`.
* About 4.08e4 DOF: B 2.08e4, V 1.20e4, n 4.0e3, T 4.0e3.
* p-MG uses `minlev=2` (p = 3, 2, 1). One-level runs also use `nlevels=2`, `nbase=2`.
* Each cell gives outer FGMRES iterations for steps 1/2/3, then the run time.

| Smoother / preconditioner | dt = 1e-6 | dt = 1e-5 |
|---|---|---|
| 1-level, BJ16 + loverlap (LU) | 3/3/2, 29.3 s | 13 (parent baseline) |
| 1-level, star symmetric | **2/2/2, 7.4 s** | 74–149 (fails) |
| 1-level, star additive | 4/4/4, 7.7 s | – |
| p-MG, BJ16 + loverlap | 3/2/2, 33.4 s | 6/6/6, 47.0 s |
| p-MG, star symmetric | **2/2/2, 8.5 s** | 136/100/96, 154.7 s (fails) |
| p-MG hybrid (star down, BJ16 + loverlap up) | 1/1/1, 32.3 s | **4/4/4, 44.1 s** |

### 7.3 MUG, same cube, MPI

* MPI runs add one geometric level (`nlevels=3`, `nbase=2`).
* dt = 1e-6 runs use the one-level preconditioner. All runs used 2 threads per rank, so 4-rank wall times are
  oversubscribed and omitted.
* dt = 1e-5 runs use p-MG.

| Configuration | 1 rank | 2 ranks | 4 ranks |
|---|---|---|---|
| dt = 1e-6, star symmetric | 2/2/2 | 2/2/2, 5.0 s | 3/2/2 |
| dt = 1e-6, star additive | 4/4/4, 7.7 s | 4/4/4, 4.4 s | 4/4/4 |
| dt = 1e-6, BJ + loverlap | 3/3/2 | 8/5/6, 19.3 s | 14/10/7 |
| dt = 1e-5, BJ16 + loverlap (p-MG) | 6/6/6 | 29/31/29, 55.5 s | 31/30/29 |
| dt = 1e-5, hybrid (p-MG) | 4/4/4 | **9/9/9, 32.6 s** | **10/9/9** |

### 7.4 Unit test: quadratic Lagrange Laplacian on the unit cube

GMRES iterations to convergence:

| Mode | serial | 2 ranks | 4 ranks |
|---|---|---|---|
| additive | 33 | 33 | 33 |
| multiplicative | 18 | – | – |
| symmetric | 15 | 25 | 34 |

* Stationary factors: 0.893 for additive and 0.623 for symmetric.
* A single patch covering all rows reproduces the exact solve. Symmetric mode with single-precision factors gives
  the same iteration count (15) as double precision.
* **Negative result:** the multiplicative and symmetric modes degrade with rank count, because the seam phase is
  additive; at 4 ranks symmetric is no better than additive. The degradation did not show in MUG at small CFL.

### 7.5 Spheromak tilt, p = 2 (about 2.85e5 DOF), dt = 1e-6 (Alfvén CFL ≈ 20), serial, 3 steps

Setup:

* Hand-built `MUG_sph_tilt`, matrix-free Newton, `xmhd_prefreq=20`, `ittarget=40`.
* For each step: outer iterations @ simulated time. MUG halves dt automatically when a solve fails.
* Peak is whole-machine memory used.

| Configuration | Result | Peak memory |
|---|---|---|
| 1-level star symmetric | dt halved twice; then 9 / 15 / 14 iterations (dt ≈ 2.5e-7 → 2.9e-7) | 2.4 GB |
| p-MG star symmetric | dt halved twice; 9@2.5e-7, 15@5.2e-7, 15@8.1e-7 (959 s) | 2.6 GB |
| 1-level BJ16 + loverlap, LU | **out of memory** (killed at 14 GB) | > 14 GB |
| 1-level BJ64 + loverlap, LU | **out of memory** | 14.2 GB |
| p-MG hybrid, BJ64 LU up | **out of memory** | 14.1 GB |
| 1-level BJ16 + loverlap, ILU | dt halved once; 92@5e-7, 122@9.4e-7, 60@1.3e-6 (410 s) | 2.9 GB |
| p-MG hybrid, BJ16 + loverlap ILU up | dt halved once; 78@5e-7, 79@9.6e-7, 52@1.4e-6 (1116 s) | 3.9 GB |

The parent session's LU-based 14 / 11–13 iteration baseline could not be reproduced within the 14 GB sandbox.

---

## 8. Analysis

* **Why stars work at small CFL.** The operator is mass-dominated. The only hard part is the curl near-kernel of the
  hierarchical B basis, which vertex stars contain (LWXZ condition). The results match the AFW theory: a factor of
  about 0.18 per sweep, independent of the B representation, and a two-grid factor of 0.013.
* **Why stars fail at CFL ≳ 1.**
  * Each field block alone is fine: B-only 0.08, V-only 0.11, with well-conditioned patches.
  * The coupled B+V star block is nearly singular: scaled condition about 4e8, against about 2e3 for B alone. The
    local Alfvén operator restricted to a star with homogeneous Dirichlet data on the star boundary has near-zero
    eigenvalues (local wave resonances), so the exact local solves amplify error instead of damping it.
  * The symmetric sweep compounds this (rho up to 1e10 at dt = 1e-4); the additive mode fails less violently.
  * This is a structural limitation of monolithic one-vertex patches for hyperbolic-like coupling at large CFL, not
    an implementation bug: the same code is accurate at small CFL and in the single-field blocks.
* **The p-coarse level does not rescue it.** The coarse p = 1 level lives on the same mesh, so it cannot represent
  the long-wavelength Alfvén modes any better. The coarse-only two-grid factor is 1.0.
* **Why the hybrid helps.** The large overlapped blocks damp the coupled long-range error, and the stars clean up the
  B-basis conditioning. Each smoother covers the other's weakness. But the block-Jacobi half still pays the
  rank-interface penalty (no cross-rank overlap), and its LU memory is prohibitive for the Spheromak.
* **Multiplicative modes across ranks.** Interior patches are swept multiplicatively, but seam patches are only
  applied additively after one global residual. As the seam fraction grows with rank count, the benefit goes away.
  A per-colour seam sweep would fix this at the cost of one mat-vec and exchange per seam colour.

---

## 9. Testing status

**New tests**

* `lin_alg/test_native_patch.py`, 11 cases:
  * additive, multiplicative, symmetric, whole-matrix (exact) and single-precision runs, serial;
  * Richardson factors for additive (`RHO[1] = 0.8932600`) and symmetric (`RHO[3] = 0.6228771`);
  * MPI runs on 2 and 4 ranks (`ITS_MPI`: additive 33/33, symmetric 25/34);
  * `converged_norm 0.79845123052314`; serial `ITS_ADD 33`, `ITS_MULT 18`, `ITS_SYM 15`.
* `physics/MUG/test_alfven.py` and `test_alfven_lag.py`: `test_nl_r1_p2_patch` (serial) and
  `test_nl_r1_p2_patch_mpi` (2 ranks), with the same expected errors as `test_nl_r1_p2`. These are the first MUG
  tests that use a non-LU multigrid smoother.

**Full-suite runs on the final build** (`OMP_NUM_THREADS=2`):

| Suite | Passed | Failed | Skipped |
|---|---|---|---|
| `lin_alg` | 33 | 2 | 84 |
| MUG `test_alfven.py` | 38 | 0 | 19 |
| MUG `test_alfven_lag.py` | 38 | 0 | 19 |
| Other 12 MUG files | 103 | 0 | 18 |

* The other MUG files are gyro, hall, viscous, revolved and solid (H(curl) and Lagrange), plus 2D Alfvén, Jacobian,
  sound and Taylor–Green.
* The two failures are `test_native_bjacobi.py::test_ILU_mpi4_nlocal2` and `..._part` (76 iterations against 70
  expected). They are believed environment-related and unrelated to this change, since block-Jacobi and ILU are
  untouched, but this was **not** confirmed by re-running on the base commit.
* **Not run:** the 3D sound tests (`test_sound.py`, `test_sound_lag.py`, 169 tests), for time. Every MUG run now
  builds patch lists, so CI should cover them.

---

## 10. Known issues and loose ends

1. **`oft_bjprecond` XML bugs** (found while reading; not fixed):
   * `slice_group=read_tmp(1)` (around `native_solvers.F90:2445`) copies the first group id into every entry.
   * `sflag` is uninitialised for non-owned rows when `boverlap` and `loverlap` are both enabled.
2. `<noverlap>N</noverlap>` (multi-layer graph overlap for block-Jacobi) was planned as a baseline but not done.
3. **Patch memory at p = 3** is about 12 GB for Spheromak-sized meshes. Use `single`, small meshes, or sparse or
   matrix-free patch solves (future).
4. The multiplicative modes apply seam patches additively (Section 8).
5. The patch smoother is unavailable under PETSc, which aborts.
6. Patches are attached only by MUG. Other physics modules would need a `comp_fem_vertex_patches` call on their
   operator.
7. `update()` always refactors; there is no reuse of the patch structure beyond a pattern change.

---

## 11. Future directions (ranked)

1. **Physics-based block preconditioning around the stars.** For example, an Alfvén Schur complement or
   augmented-Lagrangian splitting, with star-patch MG (or AMS) on the B block and stars/AMG on V. This is the route
   taken by Laakmann–Farrell–Mitchell and the most likely path to CFL robustness. Large effort.
2. **Larger or field-aligned patches.** Two-layer (macro) stars, or patches elongated along B₀, capture more of the
   Alfvén coupling. Cheap to prototype offline with the Phase 0 scripts by changing the patch lists. Memory grows
   quickly.
3. **Cross-rank overlap for block-Jacobi**, to make the hybrid scale. The seam-assembly machinery in
   `oft_patch_solvers` (global ids, partial-M exchange) can be reused. Medium effort.
4. **Sparse or matrix-free patch solves** (FDM/sparse LU) to cut memory at p ≥ 3.
5. **Per-colour seam sweep** for proper multiplicative behaviour across ranks.
6. The two small `oft_bjprecond` fixes and `<noverlap>N`.

---

## 12. Related findings from the parallel basis-conditioning study

A parallel session studied the B-basis conditioning directly; its summary was relayed by the parent session. Its
smoother experiments (cell-overlap Schwarz) come from branch `claude/schwartz-cell-adjacency-overlaps-15f6zn`. These
are its reported results; they were **not re-measured here**.

* **Changing the basis does not help.** A Schöberl–Zaglmayr basis spanning the same H(curl)+grad space gives the
  same block-smoother conditioning, within a few percent. Exact local solves are invariant to basis changes inside
  entity blocks.
* **Partitioning is not the cause.** On the B mass block at p = 1/2/3, entity-aligned blocks (`xmhd_nparts`) give
  rho 1.02/1.80/1.99 and METIS partitions give 1.21/1.85/1.70.
* **The bad modes at CFL ≈ 0 are global, smooth and clustered.** Gradient edge functions nearly cancel face
  rotationals (p ≥ 2) or Whitney functions (p = 1); each mode spans about 1e3 entities.
  * Two-level κ (p = 1 coarse + blocks) for H(curl) vs Lagrange B: 52 vs 5.7 (refined cube, p = 2) and 347 vs 8.6
    (tiny mesh, p = 3).
  * Two-grid rho with cell-overlap-1 RAS: 0.97 vs 0.61.
* **Practical impact in MUG p-MG is modest:** at most about 2× the iterations of Lagrange B at small and moderate dt.
  At large dt both B representations are limited by Alfvén coupling, consistent with Section 8.
* **Spheromak p = 2, dt = 1e-6, iterations per step:**
  * plain block-Jacobi: 74–119;
  * any overlap (`loverlap`, or cell-Schwarz overlap 0): about 24–27;
  * cell overlap 1 (`<pre type="add_schwarz"><overlap_mode>cell</overlap_mode><overlap>1</overlap>`): best on the
    Alfvén cube (dt = 1e-4, p = 2: 11.5 iterations vs 26 with `loverlap`), but runs out of memory or fails in
    UMFPACK at p = 3 and on the spheromak.
* **MPI:** 2 ranks × 8 blocks with `loverlap` take 91–114 iterations against 10–13 serial; `boverlap` recovers
  about half. This agrees with the motivation for the seam patches in Section 4.4.

**Open question: the two studies need reconciling.** That study concludes local smoothers cannot fix the CFL ≈ 0
modes. Here, vertex stars gave a stationary factor of 0.08 on the B block and 0.175 on all fields, with two-grid
0.013, at dt = 8e-7 (Section 7.1). A likely explanation is that the near-cancelling gradient/rotational
combinations lie inside vertex stars (LWXZ condition) but not inside entity, METIS or one-layer-overlap blocks. If
so, stars already deliver what a Hiptmair–Xu nodal auxiliary-space correction for the B block would target, at
least at small CFL. To check: apply both diagnostics (global-mode spectra, and the star rho) to the same dumped
operator.

---

## Appendix A: reproduction

**Environment**

* Ubuntu sandbox, GCC, OpenMP + MPI (OpenMPI, run as root), OpenBLAS, UMFPACK, SuperLU, ARPACK, HDF5, METIS.
* About 14–15 GB RAM.
* Single flat build tree `builds/build_release` (see the `oft-build` skill).

**Tests**

```bash
cd builds/build_release/tests
export OMP_NUM_THREADS=2 OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
./run_test.sh lin_alg/test_native_patch.py
./run_test.sh physics/MUG/test_alfven.py -k patch
./run_test.sh physics/MUG/test_alfven_lag.py -k patch
```

The full Alfvén files take about 86 and 32 minutes.

**Alfvén cube experiments**

These ran the test executables directly. The input below is what the run script generated. `PER=F` and `LIN=T`
for all Section 7 cube runs; for MPI, `nlevels` and `minlev` increase by one and `nbase` stays 2.

```
&runtime_options ppn=1 debug=0 test_run=T /
&mesh_options meshname='cube' cad_type=92 nlevels=2 nbase=2 grid_order=1 /
&cube_options mesh_type=1 ni=1,2,4 rscale=1.,1.,2. shift=-0.5,-0.5,-1. ref_per=F,F,F /
&test_alfven_options order=3 minlev=2 delta=1.e-3 linear=T /
&xmhd_options xmhd_jcb=T xmhd_advec=F xmhd_hall=F xmhd_adv_temp=F xmhd_adv_den=F
  bbc='bc' nu_visc=1.E-3   ! test_alfven_lag additionally: k_db=1.
  vbc='all' dt=1e-5 eta=1.e-3 nsteps=3 rst_freq=1000 lin_tol=1.E-9 nl_tol=1.E-5
  nu_xmhd=-1 nclean=5000 rst_ind=0 maxextrap=2 ittarget=100 /
&xmhd_plot_options t0=.9e-5 dt=1.e-5 rst_start=0 rst_end=1 /
```

* Command: `mpirun --oversubscribe -np N builds/build_release/tests/physics/MUG/test_alfven oft.in in.xml`.
* `minlev` = `nlevels` (that is, p = 3 only) gives one level.
* The XML files are those in Section 6. `BJ16 + loverlap` is
  `<pre type="block_jacobi"><nlocal>16</nlocal><loverlap>T</loverlap><solver type="lu"/></pre>` under GMRES
  (its 2, nrits 2). The ILU variant replaces the LU solver with ILU.

**Spheromak**

* A standalone driver built from the `MUG_sph_tilt` example (`cyl_tilt.h5` mesh), with `minlev` added to its
  namelist. It is linked by hand against the installed OFT libraries; the object file must come before the
  libraries on the link line.
* Serial runs need `nlevels=1`, `nbase=1`; `nbase=2` gives 2.2e6 DOF and runs out of memory.
* Key `xmhd_options`: `xmhd_adv_den=F xmhd_adv_temp=F bbc='bc' vbc='all' dt=1e-6 eta=1. nu_visc=10. lin_tol=1.E-9
  nl_tol=1.E-5 xmhd_mfnk=T ittarget=40 mu_ion=2. xmhd_prefreq=20`.
* Key `sph_tilt_options`: `order=2 linear=F b0_scale=1.e-1 b1_scale=1.e-4 n0=1.E19 t0=3.`.

## Appendix B: offline analysis procedure (not committed)

* **Temporary hook.** `xmhd_dump_ops`, included in `xmhd.F90`/`xmhd_lag.F90` behind the `OFT_DUMP_DIR` environment
  variable and called once after the first operator build.
  * For each level it calls `native_matrix_setup_full`, then writes as text: `nr nnz ni nred`, `kr`, `lc`, the
    `Mfull` values, `redind`, each field's `i_map` (offset, n, `offsetg`, `lge`), the patch lists and colours, an
    8-part `comp_fem_partition` (for the block-Jacobi/RAS cases), mesh edge and periodic-vertex maps, and the
    interpolation matrix P (sparse, with vectors from `xmhd_ml_vecspace%interp` for checking).
  * `OFT_DUMP_STOP=1` aborts after the dump.
  * Gotcha: read `nr` values per row buffer and take `[:nc]`; the vals buffer is reused at the larger size.
* **Python.** Run with `numpy`/`scipy` only.
  * `oftload.py`: load the dump and reduce to non-redundant rows and columns, summing periodic-copy rows of P.
  * `smooth.py`: block-Jacobi, RAS with k layers, star additive/multiplicative/symmetric as dense-LU patch
    operators; `stationary_rho` by power iteration; `two_grid_rho` with a Galerkin coarse level.
  * `analyze.py`: the driver.
  * `diag2.py`: the per-field-block diagnosis in Section 7.1.
  * The scripts are short (under 450 lines in total) and can be rewritten from this description. Use the lean
    (star-only) cases on meshes above about 1e5 DOF: the RAS cases ran out of memory in 14 GB.
