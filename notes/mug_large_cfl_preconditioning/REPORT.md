# Preconditioning MUG (`xmhd.F90`) at large Alfvén CFL: findings report

Branch: `claude/pmg-large-cfl-preconditioning`. Code commits: `9f3b4c6` (physics-based preconditioner types 1–3 and the K-cycle MG option) and `4f23692` (B-Schur, type 4). The branch starts from `claude/busy-einstein-4qu62p` @ `5b4252d`.
Status: **paused**. This document records what was tried and what was learned, so the work can be picked up later.
**Update (2026-11-27):** the workstation-scale benchmark campaign from [`HANDOFF_BENCHMARKS.md`](HANDOFF_BENCHMARKS.md) has been run on an Apple M2 Max. See [`RESULTS_M2Max.md`](RESULTS_M2Max.md) for the accuracy cost of `xmhd_si` on the spheromak, scaling in p and h, and timing against production-style preconditioners.

---

## 1. Summary

| Question | Answer |
|---|---|
| Does Krylov-accelerated MG (K-cycle, coarse-correction scaling) fix the large-CFL failure? | **No.** It made no measurable difference on the coupled system or inside the block preconditioner. The coupled two-grid cycle itself diverges, and Krylov acceleration cannot fix that. |
| Does the textbook "parabolized velocity Schur" preconditioner (NIMROD / M3D-C1 / PIXIE3D style) work for `xmhd.F90`? | **No.** Its iteration count stalls at about 35 per step at dt=1e-4 on the cube. The continuum operator over-estimates the discrete Schur complement (about 2.2× in norm). Its smallest eigenvalue relative to the exact Schur complement drops to 3e-5 at dt=1e-3. Cause: B lives in H(curl), and the continuum form ignores the discrete projection. |
| Is there a Schur approximation that is dt-robust? | **Yes: eliminate V instead of B ("B-Schur", type 4).** S_X = J_XX − J_XV·diag(J_VV)⁻¹·J_VX. With exact inner solves on the cube, dense GMRES needs 9 / 26 / 30 iterations at dt 1e-5 / 1e-4 / 1e-3, against 13 / 54 / 500 (iteration cap) for the best velocity-Schur variant. In MUG with LU inner solves: 5 / 14–15 / 19–40 outer iterations per step. |
| Is it practical yet? | **Not for production-size problems.** S_X has 2.6× the nonzeros of J. Cheap inner solvers lose much of the gain (block-Jacobi 16–92 its; p-MG fails at dt=1e-3), and ILU(0) hits a zero pivot on S_X for the spheromak. On the spheromak at dt=1e-5 with LU on 1024 blocks, type 4 held dt=2.5e-6 against 1.25e-6 for the baseline, but used 6.3 GB and was too slow (one step in 50 min). |
| Would a NIMROD-style semi-implicit term (in the time discretization) help? | **Conditional yes (§7a, offline, cube only).** Stable with scale factor s ≥ 1.25 up to CFL ~1900; the continuum operator matches the discrete one within 8–17% on the physical Alfvén modes (the large over-estimates are on spurious modes); resolved-mode phase accuracy is within a few % of CN. Non-uniform equilibria and nonlinear terms are untested. |
| Can the SI operator live inside the normal CN/BDF2 step? | **Yes; implemented (§7b, `xmhd_si`).** It is the increment s·(dt/2)²·W·(Vⁿ⁺¹−Vⁿ). With it, PB type 1 is dt-independent in linear cube runs to CFL ~190 (8–14 its vs 13–39), and the spheromak runs at **dt = 1e-5 without dt cuts** (115 its, 762 s per 3 steps), where s=0 is cut to 2.5e-6. Cost: a time error about 2.5× CN's at s=1 (cube). Nonlinear CFL ≳ 100 is then limited by MUG's lagged Newton Jacobian; spheromak accuracy is not yet measured. |
| What does it cost and buy at workstation scale? ([`RESULTS_M2Max.md`](RESULTS_M2Max.md)) | **Accuracy:** tilt growth-rate error Δγ/γ ≈ −1.5 %·s·(dt/1 µs)². The nonlinear tilt is within 2 % at dt=1e-6 with s=1 and 8–21 % low at 3e-6, so accuracy limits s=1 to dt ≲ 1–2 µs. **Speed:** PB + s=1 at dt 1e-6 is 1.4–2.4× faster than the example ILU at dt 2e-7, but only 1.06–1.18× faster than full-Jacobian block-Jacobi LU at s=0, dt 1e-6 (which is accurate). It uses about half the memory. SI with full-Jacobian LU looks faster still (being confirmed). **Bottleneck:** with an exact LU of P, outer iterations are 5–14 at every p, h and CFL up to 300, so the inner solve on P is where the remaining cost is. **Railgun (rotational periodicity):** correct; PB uses 2.7× fewer iterations and 35 % less memory. |
| Recommended next step | Use the **SI increment term (§7b)** with PB type 1; validate its accuracy cost on the spheromak and improve the inner solver on P (an **auxiliary-space (Hiptmair–Xu / AMS-type) solver** for the B block), which is curl–curl-like with a large gradient near-null space. Alternative: a matrix-free exact-Schur variant (planned "type 5"). |

---

## 2. Problem statement and prior evidence

MUG (`src/physics/xmhd.F90`) solves each time step with Newton–Krylov: outer FGMRES on the assembled Jacobian J, with Crank–Nicolson centering (`jac_dt = dt/2`).
- Fields: B = B_c (H(curl), field 1) + B_g (grad H1, field 2) (+ `j2` for hyper-resistivity); V = vector Lagrange (fields 3–5); n (6); T (7); plus optional Te and n2.
- J is a single composite native CSR matrix built from `xmhd_mat_mask`.

Findings from the earlier investigation (branch `claude/busy-einstein-4qu62p`):
- At Alfvén CFL ≳ 1, no smoother or p-MG configuration converged.
- Non-overlapping block-Jacobi(LU) diverges (spectral radius up to 19).
- The energy-norm B–V Schur coupling grows like CFL²: 0.06 / 9 / 900 at dt 8e-7 / 1e-5 / 1e-4 on the cube.
- The Galerkin coarse correction is an oblique projection that *amplifies* error: ρ(S·C) = 3.16 > ρ(S) = 2.04.
- Transfers and coarse operators were verified correct.
- On the spheromak-tilt example at dt=1e-5, nothing converged.
- Conclusion: the B–V hyperbolic coupling is not elliptic, so smoothing plus coarse correction cannot work on the fully coupled matrix. The coupling has to be removed algebraically first.

Test problems used here:
- **Alfvén-wave cube** (`src/tests/physics/MUG/test_alfven`, `mug_testing.ALFVEN_TEMPLATE`): r1, p=2, `nlevels=1` unless noted, `nu_visc=1e-3`, `bbc='bc'`, 3 steps. Global sizes after folding periodic copies: B 1856, V 768, n/T 512. Alfvén CFL ≈ 0.25 / 3 / 30 / 300 at dt 8e-7 / 1e-5 / 1e-4 / 1e-3.
- **Spheromak tilt** (`src/examples/MUG/Spheromak_tilt`): p=2, `nlevels=1` (serial), fixed dt, 3 steps. DOFs: B 1.66e5, V 7.1e4, n 2.4e4, T 2.4e4.

---

## 3. Literature and other codes

Tags: **[P]** proven; **[E]** empirical (in the literature or other codes); **[I]** our inference.

**Krylov-accelerated MG**
- Notay & Vassilevski, *NLAA* 15 (2008), K-cycle; Notay, *NLAA* 17 (2010), nonsymmetric two-grid analysis.
  - [P] If the two-grid method converges uniformly, the K-cycle converges uniformly at V-cycle-like cost.
  - [I] Here the two-grid method diverges, and the number of amplified modes (k·v_A·Δt ≳ 1) grows with CFL and with refinement. A small Krylov space can only remove a bounded number of outliers. **Confirmed empirically (§5.5).**

**Physics-based / approximate block factorization**
- Chacón, Knoll & Finn, *JCP* 178 (2002); Chacón, *Phys. Plasmas* 15, 056103 (2008) (PIXIE3D).
  - [E] Eliminate B, ρ and T and parabolize the velocity Schur complement: S_V = ρ/Δt + θ²Δt·W(B0, p0). Then use MG on S_V.
  - [E] Iteration counts are nearly flat up to CFL ~10³.
  - **For `xmhd.F90` this works only after reversing the elimination order (§6.1).**
- Chacón, arXiv:2407.07031 (2024/25), Hall MHD. [E] The same factorization works with a vector potential and electron-MHD block treatment. Relevant for future Hall support.
- Knoll & Keyes, *JCP* 193 (2004), JFNK review. [E] Recommends physics-based preconditioning as the standard remedy for stiff waves.
- Sandia Drekar/Teko: Shadid et al., *CMAME* 2016; Cyr et al., *SISC* 2013; Phillips, Elman, Cyr, Shadid & Pawlowski, *SISC* 2014/2016.
  - [E] Fully coupled AMG with block smoothers works at moderate coupling.
  - [E] At high CFL or Lundquist number they switch to approximate block factorization with physics-based Schur complements.
- Ma, Hu, Hu & Xu, *JCP* 2016. [P] Robust block preconditioners, but they require structure-preserving (exact-sequence) spaces, which Lagrange V is not.
- Laakmann, Farrell & Mitchell, *SISC* 44 (2022), augmented Lagrangian (arXiv:2104.14855). [P/E] Very robust, but needs exactly divergence-free spaces and star-patch MG. Not applicable without changing the discretization.

**Production codes**
- **NIMROD** (Sovinec et al., *JCP* 2004). [E] Semi-implicit operator ρ − Δt²·L_ideal(B0, p0) in the velocity advance, solved with SuperLU_DIST per Fourier block. **B is nodal in NIMROD**, so the continuum operator matches the discrete one.
- **M3D-C1** (Jardin et al.). [E] Split-implicit velocity solve with the same parabolization, using direct solves per plane.
- **JOREK** (Hölzl et al.). [E] GMRES with block-Jacobi over toroidal modes, using direct solves per block (PaStiX / MUMPS / STRUMPACK).
- **BOUT++** (Dudson et al.). [E] Shear-Alfvén Schur parabolization in reduced MHD, inverted with Laplacian solvers.
- **Consensus [E]:** none of these smooth the coupled hyperbolic system directly. All of them either parabolize a velocity operator or rely on direct solves of large blocks.

---

## 4. What was implemented

All new options are opt-in through XML. With no new XML, the defaults are unchanged; the non-slow `physics/MUG` regression suite passed after the refactor (166 passed, 56 skipped).

### 4.1 M1 — K-cycle and coarse-correction scaling in generic MG
- `src/lin_alg/native_solvers.F90`:
  - New fields on `oft_ml_precond` (l.229–230): `kcycle` (number of flexible-GMRES steps on the coarse correction; 0 means a plain V-cycle) and `cscale` (logical).
  - `ml_precond_apply` (~l.1843 and ~l.1856) calls the internal routines `kcycle_correction` (l.1890) and `scale_correction` (l.1988).
  - The K-cycle takes the coarse operator from `base_solve%A`, else `smooth_up%A`, else `smooth_down%A`. It is skipped for `oft_ml_trans` levels. It uses FGMRES with modified Gram–Schmidt and Givens rotations; work storage is allocated on first use and freed in `ml_precond_delete`.
  - `cscale` scales the correction by α = (r, Ap)/(Ap, Ap).
- `src/lin_alg/solver_utils.F90` `create_ml_xml` (l.832–880): reads `kcycle="N"` and `cscale="T"` on the MG `<pre>` node, applies them to every level, and sets `A` on each coarse level.

```xml
<pre type="mg" kcycle="2" cscale="T"> <smoother direction="up"> ... </smoother> </pre>
```

### 4.2 M2 — Physics-based block preconditioner in `xmhd.F90`
- **XML:** a `<pb_pre>` element inside `<xmhd>` (read in `xmhd_read_settings`, l.542–559). Attributes:
  - `type`: 1 = full block factorization, 2 = lower triangular, 3 = block diagonal, 4 = B-Schur;
  - `alfven`, `sound`: scale factors on the parabolization terms (types 1–3);
  - `dump="T"`: write J, P and the field layout to HDF5 for offline analysis.
  - Content: a `<solver>` element for the inner solve, or MG children (`<smoother>`, `kcycle`/`cscale` attributes) when `nlevels>1`.
- **Type:** `oft_xmhd_pb_precond` (l.253). Procedures:
  - `xmhd_pb_apply` (l.5163) computes the residual, then predictor → middle solve → corrector, using J matvecs on group-masked vectors (`xmhd_pb_zero_groups`, l.5150);
  - `xmhd_pb_update` (l.5254) also handles the dump;
  - `xmhd_pb_delete` (l.5282).
  - Note: the doc comment on `pb_type` lists only types 1–3; type 4 is missing there.
- **Groups:** B = {B_c, B_g, j2}; V = {Vx, Vy, Vz}; S = everything else.
- **Matrix P** (`xmhd_ops%P`, l.143), one per level, in array `ml_P`.
  - Types 1–3: P is block-diagonal, from mask `xmhd_pb_mask` (l.2881–2888). Its B and S blocks are copies of J's blocks; its V block is J_VV plus the parabolization.
    - The parabolization integrand is assembled in `xmhd_build_ops` (~l.2579, coefficients `pb_ca`, `pb_cs`). Form: S_V(φ,v) = J_VV + jac_dt²/(n·m_i)·[(1/μ0)·Q(φ)·Q(v) − (1/μ0)·φ·(J0×Q(v)) + (∇·φ)(γ·p0·∇·v + v·∇p0)], with Q(v) = (B0·∇)v − B0(∇·v).
  - Type 4: P holds J_VV on the V rows and S_X = J_XX − J_XV·diag(J_VV)⁻¹·J_VX on the X rows (X = all non-V fields). It is built by `xmhd_pb_build_bschur` (l.4881), called from `xmhd_set_ops` after each `xmhd_build_ops` (l.4779, l.4813).
    - The local J is folded into global CSR through `lge` (duplicate periodic or shared copies summed); the global diagonal of V is taken; rows are accumulated in two passes (`q_row`, `acc_add`).
    - On the first call, per-block `oft_graph` objects are built and passed to `create_matrix(..., native=.TRUE.)`. Entries go on the representative (first) local copy, with zero diagonal placeholders on the redundant copies.
    - Type-4 apply: x_V = J_VV⁻¹ r_V (with P's V block), then an S_X solve on r_X − J_XV·x_V, then a V corrector.
- **Refactor:** the BC handling of the local matrices was moved out of `xmhd_build_ops` into `xmhd_local_bcs` (l.2674) and `xmhd_bc_diag` (l.2750), so J and P get identical BCs. In solid regions, P gets the B blocks only.
- **Wiring:** `xmhd_create_pb` (l.4832) is called from `xmhd_run` (l.912) and `xmhd_lin_run` (l.1477) when `<pb_pre>` is present. It prints "Using physics-based block preconditioner".

Example (type 4, overlapped block-Jacobi inner):
```xml
<oft><xmhd>
  <pb_pre type="4">
    <solver type="gmres"><its>5</its><nrits>5</nrits>
      <pre type="block_jacobi"><nlocal>8</nlocal><loverlap>T</loverlap><solver type="lu"></solver></pre>
    </solver>
  </pb_pre>
</xmhd></oft>
```

**Limitations / known gaps**
- Type 4 is serial and native-backend only (it folds via the global `lge` map). MPI is not implemented.
- `dump` triggers in `update`, and the first preconditioner update happens at step 2, so dumping needs `nsteps ≥ 2`.
- The type-4 matrix keeps zero-diagonal placeholder rows for redundant local copies. This is harmless for LU but suspect for the ILU failure (§5.7).
- No regression tests or user documentation were added (planned milestone M3, not done).
- Hall terms are not treated; there is no robustness claim with Hall on.

---

## 5. Results

Cube runs: outer FGMRES iterations per time step for 3 steps (Newton iterations are 1 unless noted). "dt cut" means MUG's Newton failed and it silently reduced dt; the final dt is given. Raw data: `results/`.

### 5.1 Stage (i): continuum velocity Schur (types 1–3) with LU blocks, cube

| dt (Alfvén CFL) | LU on full J | 1-level overlap BJ(LU), nlocal 8 | PB type 1 (full) | PB type 2 (lower) | PB type 3 (diag) |
|---|---|---|---|---|---|
| 8e-7 (0.25) | 1,1,1 | 5,4,4 | 3,2,2 | 3,2,2 | 4,4,4 |
| 1e-5 (3) | 1,1,2 | 7,6,6 | 12,11,11 | 12,12,12 | 18,18,18 |
| 1e-4 (30) | 1,1,4 | 29,30,29 | 35,38,39 | 37,38,43 | 86,86,97 |
| 1e-3 (300) | 3,2,17 (Newton 2,1,2) | 97,120,213, **dt cut → 3.2e-4** | 86,53,173 | 92,54,288 | 210,162,190, **dt cut → 5e-4** |

Errors from types 1–3 agree with the LU reference wherever dt was not cut. Memory: PB about 93 MB vs about 120 MB for full LU and 210 MB for overlap BJ.

**Alfvén-scale sweep** (type 1, dt=1e-4): `alfven` = 0 → 65,72; 0.25 → 34,35; 1 → 35,38; 4 → 46,46; −1 → 87,94. The sign and scale are right, but the operator stalls around 35, so the *form* of the approximation is the problem, not the coefficient.

### 5.2 Offline Schur analysis, cube (`scripts/schur_analysis.py`, full output in `results/offline_schur_analysis.txt`)
J and P are dumped with `dump="T"`, folded to global numbering, and compared against the exact Schur complements using dense linear algebra.

Dense right-preconditioned GMRES iterations (tol 1e-8, random RHS, exact inner block inverses):

| dt | block-diag(J) | V-Schur, continuum S_V (types 1–3) | V-Schur, diag(J_XX) | **B-Schur, diag(J_VV) (type 4)** | exact Schur |
|---|---|---|---|---|---|
| 1e-5 | 22 | 24 | 13 | **9** | 1 |
| 1e-4 | 500 (max) | 166 | 54 | **26** | 1 |
| 1e-3 | 500 (max) | 500 | 500 | **30** | 2 |

These counts come from re-running the committed script on the saved dumps when the report was written. An earlier interactive run gave 159 / 157 / 55 for the first three columns at dt=1e-4 (probably a different dump or version of the script); the conclusions are unchanged.

Spectra:
- eig(S_V,cont⁻¹·S_V,exact), smallest |λ|: 0.22 / 3e-3 / 3e-5 at dt 1e-5 / 1e-4 / 1e-3.
- ‖δW_continuum‖ / ‖δW_exact‖ = 2.24 at dt=1e-5, where δW = S − J_VV.
- eig(S_B,diag⁻¹·S_B,exact) range: [0.76, 1.22] / [0.32, 2.83] / [0.26, 3.74]. This range stays bounded as dt grows 100×.
- Row-sum lumping |J_VV| in place of diag(J_VV) gives the same counts.

### 5.3 Type 4 (B-Schur) in MUG with LU inner solve, cube

| dt | type 4 + LU | (LU on full J) | (overlap BJ) |
|---|---|---|---|
| 8e-7 | 2,1,1 | 1,1,1 | 5,4,4 |
| 1e-5 | 5,5,5 | 1,1,2 | 7,6,6 |
| 1e-4 | 14,15,15 | 1,1,4 | 29,30,29 |
| 1e-3 | 31,19,40 (Newton 2,1,2), no dt cut | 3,2,17 | 97–213, dt cut |

Errors match the LU reference.

### 5.4 Stage (ii): type 4 with scalable inner solvers, cube

| configuration | dt=1e-5 | dt=1e-4 | dt=1e-3 | RSS |
|---|---|---|---|---|
| type 4 + GMRES(2)/overlap BJ(LU, nlocal 8) | 6,6,5 | 23,23,24 | 118,79,130 | ~390 MB |
| type 4 + GMRES(5)/overlap BJ(LU) | 5,5,5 | 16,17,17 | 61,44,92 (no dt cut) | ~390 MB |
| type 4 + p-MG (p2→p1), overlap BJ smoother | 6,5,5 | 23,24,24 | 122,179,37, **dt cut → 1.6e-4** | ~430 MB |
| type 4 + p-MG + K-cycle(2) | 6,5,5 | 23,24,24 | 96,159,37, **dt cut → 1.6e-4** | ~430 MB |
| p-MG on full J (reference) | 6,5,5 | 28,28,28 | 36,44,80, **dt cut → 1.6e-4** | ~230 MB |
| p-MG on full J + K-cycle(2) | 6,5,5 | 28,28,28 | 37,46,80, **dt cut → 1.6e-4** | ~230 MB |

### 5.5 K-cycle / coarse-correction scaling
- No measurable effect in any configuration (§5.4). This was a 2-level hierarchy, which is weak evidence.
- It is consistent with the expectation in §3: the cycle being accelerated diverges on many modes, not a few.
- `cscale` was implemented but not studied systematically.

### 5.6 Memory and fill
On the cube, S_X has nnz ≈ 1.25M vs nnz(J) ≈ 0.48M (2.6×). The reason is that eliminating V through its diagonal couples every pair of B DOFs that share a velocity node: a two-element stencil.

### 5.7 Spheromak tilt (serial, p=2)

| configuration | dt requested | outer its / step | Newton | dt actually taken | wall | peak RSS |
|---|---|---|---|---|---|---|
| example preconditioner: GMRES(8)/BJ(nlocal −1)/ILU on J | 1e-6 | 136,117,118 | 2,2,2 | 1e-6 | — | 1.36 GB |
| same | 1e-5 | 149,179,208 | 2,2,2 | **cut → 1.25e-6 … 1.58e-6** | — | — |
| GMRES(8)/BJ(LU, nlocal 1024) on J | 1e-5 | 132,178,205 | 2,2,2 | **cut → 1.25e-6 … 1.58e-6** | 1040 s | 3.1 GB |
| type 4 + GMRES(8)/BJ(LU, nlocal 1024) | 1e-5 | 159 (1 step) | 2 | **cut → 2.5e-6** | timeout at 3000 s | 6.3 GB |
| type 4 + GMRES(8)/BJ(nlocal −1)/ILU | 1e-5 | — | — | **ILU(0) zero pivot** ("Factorization failed", `ilusolver_apply`) | — | ~4.2 GB at setup |
| GMRES(8)/overlap BJ(LU), nlocal 8 and 64, on J | 1e-6 | — | — | **out of memory** (~13.6 GB) during factorization | — | >13 GB |

Conclusions:
- Type 4 doubles the time step that Newton accepts on the spheromak.
- With block LU inner solves it costs about 2× the memory and much more time per step.
- dt=1e-5 was not reached by any configuration.
- With the sandbox's ~14 GB, larger LU blocks are not an option.

---

## 6. Analysis

### 6.1 Why the continuum velocity Schur fails in `xmhd.F90`
- The discrete Schur complement on V is C_VB·M_B⁻¹·C_BV: the induction coupling goes *through the conforming H(curl)+grad space*, and only the part of curl(v×B0) representable there contributes.
- The continuum form ∫Q(φ)·Q(v) evaluates curl(v×B0) pointwise. That acts like a broken (element-wise) L² projection, which is larger, especially at high wavenumber. Hence ‖δW_cont‖ ≈ 2.2·‖δW_exact‖.
- More importantly, there are directions where the exact δW is far smaller than the continuum δW. Along those, S_V,cont⁻¹·S_exact has eigenvalues of order 1/CFL² (3e-5 at CFL 300). GMRES cannot handle such a cluster cheaply.
- B0 is uniform and η and ν are small in the cube, so the dropped ∇B0, resistive and viscous terms cannot explain the stall. The `alfven` sweep confirms that no rescaling helps.
- Codes where this works (NIMROD, M3D-C1, PIXIE3D) use nodal or finite-difference B, where the continuum and discrete operators agree. **Expected to work for `xmhd_lag.F90` (nodal B); not tested.**

### 6.2 Why the B-Schur is dt-robust
- Eliminating V uses J_VV ≈ M_V/dt + ν-terms. Here M_V is a *Lagrange* mass matrix, and its diagonal is spectrally equivalent to it with dt-independent constants.
- So S_X,approx and S_X,exact differ only through that spectral equivalence, and the eigenvalue range stays bounded ([0.26, 3.74] at CFL 300).
- In contrast, the diagonal of the *hierarchical H(curl)* mass is a poor approximation (diagonally scaled condition number about 3271 at p=3), which is why V-Schur with diag(J_XX) also fails.

### 6.3 Why cheap inner solves on S_X degrade
- S_X ≈ M_B/dt + (dt/ρ)·curl·(B0 × ·)ᵀ(·)curl-type terms. It is a curl–curl-like operator with a large near-null space (gradients and field-aligned components) whose stiffness grows with CFL.
- Standard (Lagrange-style) smoothers and p-MG do not handle curl–curl near-null spaces. This is the textbook reason for Hiptmair smoothers and auxiliary-space (HX/AMS) methods.
- ILU(0) zero pivot: probably from the placeholder rows for redundant copies, or from the indefinite or ill-conditioned S_X pattern. Not diagnosed.

### 6.4 Why the K-cycle cannot help
The K-cycle accelerates a convergent or nearly convergent cycle. On the coupled J the two-grid cycle amplifies every mode with k·v_A·Δt ≳ 1, and their number grows with CFL and resolution. Inside type 4 the bottleneck is the smoother on curl–curl modes, not the coarse correction.

---

## 7. Issues found that predate this work
1. **Newton failure and silent dt reduction:** in the test setup, MUG's Newton fails at dt=8e-7 even with a full-J LU solve and silently halves dt (2e-7 is the first accepted step in `smoke_*` runs). Worth investigating separately; it affects how test results are interpreted.
2. MUG's adaptive dt cut, combined with `ittarget`, can hide preconditioner failures as smaller time steps. The experiments set `ittarget=100000` to keep dt fixed where possible.

---

## 7a. Semi-implicit (NIMROD-style) feasibility study — Phase A, offline

**Question:** instead of using the parabolized operator inside the *solver* (where §5.1 showed it fails), should MUG use it in the *time discretization*, as NIMROD and M3D-C1 do? In a semi-implicit (SI) advance the mismatch between the continuum operator and the exact discrete coupling becomes a time-discretization error rather than extra Krylov iterations, and every linear solve is near-elliptic: a V solve with M_V + s·jdt²·W (which is exactly the V block of the type-1 P matrix, already assembled) followed by B and n/T solves that are mass and diffusion dominated. Script: `scripts/si_analysis.py`; full output: `results/semi_implicit_analysis.txt`.

**Method.**
- Inputs: three cube dumps (type 1, `dump="T"`) at dt = 1e-5 / 2e-5 / 4e-5, none with dt cuts.
- MUG's Jacobian is J(dt) = M + jdt·L with jdt = dt/2, so two dumps give the mass matrix M and the linearized operator L; the third checks the affine fit (residual 3.5e-4, the size of the state dependence between runs). The SI operator is W = (P_VV − J_VV)/jdt² (consistent across dumps to 6e-4). M_VX is set to 0, since the fitted value is fit residue.
- Constrained (identity) rows are dropped (544 of 3136, mostly the frozen n/T fields: this test has `xmhd_adv_den/temp=F`, so there are no sound waves).
- Exact discrete ideal operator in V space: Kr = −L_VX·M_X⁻¹·L_XV, so that M_V·V'' = −Kr·V.
- SI scheme modelled (staggered, leapfrog-like; V at half steps, X = all non-V fields, explicit in the V equation):
  - (M_V + jdt·L_VV + s·jdt²·W)·V⁺ = (M_V − jdt·L_VV + s·jdt²·W)·V⁻ − dt·L_VX·X
  - (M_X + jdt·L_XX)·X⁺ = (M_X − jdt·L_XX)·X − dt·L_XV·V⁺
  - s is the SI scale factor (the existing `alfven` attribute).
- Compared against MUG's actual Crank–Nicolson map (M + jdt·L)⁻¹(M − jdt·L).

**Results.**
1. **Stability (full one-step map, including dissipation).** Spectral radius:

   | dt (max CFL) | CN | SI s=0 | s=0.5 | s=1 | s=1.25 | s=2 | s=4 |
   |---|---|---|---|---|---|---|---|
   | 1e-5 (1.9) | 1.00011 | 1.00011 | 1.00011 | 1.00011 | 1.00011 | 1.00011 | 1.00011 |
   | 1e-4 (19) | 1.0011 | 345 | 4.95 | 1.0011 | 1.0011 | 1.0011 | 1.0011 |
   | 1e-3 (186) | 1.011 | 3.4e4 | 416 | 1.011 | 1.011 | 1.011 | 1.011 |
   | 1e-2 (1863) | 1.115 | 3.2e6 | 7.3e3 | **1.45** | 1.103 | 1.102 | 1.099 |

   The excess over 1 also appears for CN. It is a weak growing semi-discrete mode (max Re μ = 11 s⁻¹, against Alfvén frequencies of about 3e4–2e5 s⁻¹), probably a fit artifact or a weak instability of the perturbed state. **With s ≥ 1.25 the SI map is never less stable than CN, up to CFL ~1900; s = 1 fails only at CFL ~1900.** Without the SI term (s=0) the explicit coupling is violently unstable, as expected.

2. **Where the continuum operator over-estimates.** Per-mode ratio s_eff = (yᵀWy)/(yᵀKr·y) on the discrete ideal normal modes (Kr·y = ω²·M_V·y):
   - The **physical** lowest Alfvén modes sit exactly at ω = π·v_A = 3.14e4 s⁻¹ (k∥ = 2π/L_z), with **s_eff = 1.08–1.17**. 42 modes have s_eff < 1.5, including the next shear and fast branches.
   - Modes with **ω below π·v_A**, which cannot exist physically in a uniform periodic box, have s_eff = 6–426. These are spurious discrete modes: k∥ ≈ 0 flows that the H(curl)/Lagrange coupling feels weakly. The large over-estimates in §5.2 and §6.1 are concentrated there and at grid scale. Slowing them is harmless and arguably beneficial.
   - 273 V modes are in the discrete null space (ω < 1e-3·ω_max), for example parallel flows.
   - Global norms: ‖W‖/‖Kr‖ = 2.24, the same ratio as in §5.2.

3. **Accuracy of physical modes** (Rayleigh estimate cos(Ω·dt) = 1 − θ/2, θ = dt²·yᵀKr·y / yᵀ(M_V + s·jdt²·W)·y; ratio of numerical to true frequency):

   | mode | ω·dt | SI s=1 | SI s=1.25 | SI s=2 | CN |
   |---|---|---|---|---|---|
   | lowest Alfvén (s_eff 1.08) | 0.31 | 0.991 | 0.988 | 0.978 | 0.992 |
   | next branch (s_eff 1.4–1.5) | 0.63 | 0.95 | 0.93 | 0.89 | 0.97 |
   | s_eff 1.33 | 1.38 | 0.83 | 0.79 | 0.69 | 0.88 |
   | lowest Alfvén | 3.14 | 0.61 | 0.54 | 0.43 | 0.64 |
   | lowest Alfvén | 31 | 0.082 | 0.066 | 0.047 | 0.096 |

   For resolved modes (ω·dt ≲ 1), SI with s ≈ 1–1.25 is within a few percent of CN. For under-resolved modes both schemes are badly wrong: CN is already at 0.64 for ω·dt = 3. SI is somewhat more dispersive, by a factor of roughly s·s_eff in the error.

**Verdict: conditional GO.** All three Phase-A criteria are met on the cube:
- stable with a modest scale factor (s = 1.25);
- the over-stabilized modes are spurious or unresolved;
- resolved-mode accuracy is comparable to CN.

Practical attraction: the SI V operator is SPD and near-elliptic, and is **already assembled** (P's V block); the B and n/T solves are mass/diffusion dominated. So there would be no coupled hyperbolic solve at all, which is precisely what defeated every preconditioner above.

**Caveats and what Phase A does not show.**
- The cube has uniform B0, no flow, no pressure dynamics (no sound waves), and is linear. W also contains J0×Q and ∇p0 terms, which are sign-indefinite in sheared or pressurized equilibria (spheromak, tokamak). The condition s·W ≥ Kr has to be re-checked there, possibly with a larger s, as NIMROD does.
- A nonlinear SI needs an extra (isotropic) operator for nonlinear-amplitude stability, and flow advection keeps its own CFL limit unless it is treated implicitly.
- The modelled scheme is a staggered leapfrog with first-order time centering of the coupling. It loses MUG's exact CN energy behaviour, and the Hall/two-fluid terms are not covered.
- The spheromak matrices are too large for this dense analysis. The real test is a prototype (Phase B): staggered V/X advance in `xmhd_lin_run` first, reusing `xmhd_ops%P` and the `<pb_pre>` inner-solver infrastructure, checked on spheromak tilt growth rates vs dt against the fully implicit runs.

---

## 7b. Semi-implicit increment term inside the CN/BDF2 step — Phase B, implemented

**Idea (user's suggestion).** Rather than a separate staggered integrator, put the SI operator *inside* MUG's existing time step, as a stabilizing term for the B–V coupling of under-resolved modes while every other term keeps its normal treatment (as is already done when dissipation is moved to backward Euler). MUG's step is A(+dt/2)·uⁿ⁺¹ = A(−dt/2)·uⁿ, so a momentum term proportional to `self%dt**2` has the same sign at both time levels and automatically becomes the increment **s·(dt/2)²·W·(Vⁿ⁺¹ − Vⁿ)**. That addition is O(dt²)-consistent. The B–V coupling stays implicit.

**Implementation** (`src/physics/xmhd.F90`, commits `9c001a8` and the follow-up consistency commit):
- New `xmhd_options`: `xmhd_si` (s, magnetic part), `xmhd_si_sound` (pressure part), `xmhd_si_full` (include the sign-indefinite J0× and ∇p0 terms; default F keeps W positive semidefinite). All default off, and the non-slow MUG suite is unchanged (166 passed, 56 skipped).
- `xmhd_si_local` (new) assembles the parabolized operator into a local matrix. It is shared by the `pb_pre` P matrix and the Jacobian.
- **Nonlinear runs (CN):**
  - The residual (`xmhd_errmatrix_apply`) adds s·self%dt²·W_ref·V.
  - W_ref is evaluated at the start-of-step state (`up`, via the module pointer `xmhd_si_ref`), so the term is exactly linear in V and an exact increment.
  - The Jacobian adds s·jac_dt²·W_ref, built at the same reference state on the finest level (coarse preconditioner levels use the linearization state).
  - A temporary finite-difference check showed the residual and the assembled Jacobian agree to about 1e-12 along velocity directions.
- **Linear runs (BDF2/BE):** the term is a modified velocity mass M_V + s·jac_dt²·W(u_eq), added to both the Jacobian and `xmhd_massmatrix_apply` (BDF2 right-hand side).

**Accuracy** (cube, regression-test setup, fixed final time T = 1e-4, exact LU solves; relative L2 waveform error of B; spatial floor about 4.27e-2; `results/si_accuracy.txt`):

| dt (ω·dt of the test wave) | CN s=0 | CN s=0.5 | CN s=1 | CN s=2 | BDF2 s=0 | BDF2 s=1 |
|---|---|---|---|---|---|---|
| 8e-7 (0.025) | 4.267e-2 | 4.267e-2 | 4.267e-2 | 4.266e-2 | 4.267e-2 | 4.266e-2 |
| 4e-6 (0.13) | 4.282e-2 | 4.318e-2 | 4.383e-2 | 4.577e-2 | 4.668e-2 | 5.126e-2 |
| 1e-5 (0.31) | 4.933e-2 | 6.123e-2 | 7.577e-2 | 0.1088 | 0.1180 | 0.1687 |
| 2e-5 (0.63) | 0.1057 | 0.1742 | 0.2441 | 0.3792 | 0.3270 | 0.5244 |

The time-error part grows by about 2.5× at s=1 (CN, ω·dt = 0.31), close to the predicted (1/12 + 1.1/8)/(1/12) ≈ 2.65. Below ω·dt ≈ 0.1 the term is invisible.

**Preconditioner robustness** (cube, `pb_pre` type 1 with LU blocks, outer iterations per step; `results/si_pb1_iterations_v1.txt`, `results/si_nonlinear.txt`):

| dt (max Alfvén CFL) | linear BDF2 s=0 / 1 / 2 | nonlinear CN s=0 / 1 / 2 |
|---|---|---|
| 1e-5 (1.9) | 13–17 / 8–9 / 6–7 | 11–12 / 8 / 6–7 |
| 1e-4 (19) | 29–33 / **11** / **9** | 35–39 / **11–12** / **9** |
| 1e-3 (186) | 36–39 / **13–14** / **10–11** | see below |

With the SI term the type-1 preconditioner is essentially **dt-independent in linear runs up to CFL ~190**. That matches the offline bound eig(P_V⁻¹·S) ≥ s/(1+s). The continuum S_V failed as a pure preconditioner (§5.1); it works once the same operator is part of the physics.

**Nonlinear runs at CFL ~190: limited by Newton, not by the linear preconditioner.**
- At dt = 1e-3 (6 steps), all settings show Newton trouble, s=0 included (dt cut to 5e-4).
- With s ≥ 1 and MUG's default assembled-Jacobian Newton, Newton needs 3–5 iterations and dt is cut.
- With matrix-free Newton (`xmhd_mfnk=T`), s=1 holds dt = 1e-3 (Newton 2–3), but with 52–152 Krylov iterations per step.
- Diagnosis by temporary finite-difference checks (since removed):
  - The SI term is exactly consistent.
  - The mismatch is entirely in the pre-existing Lorentz block ∂R_V/∂B, and is about 100% in the V∥B0 rows. Those rows hold only O(δ) terms ((J×b)_z + (curl b × δB)_z) that depend on the *current* B.
  - MUG builds the Jacobian at the start of the step, and at CFL ≳ 100 the perturbation changes by order one within one step, so those entries are stale. This affects s=0 as well; the SI term exposes it because Newton takes larger corrections.
  - With the Jacobian refreshed at every Newton iteration (temporary test), s=1 held dt = 1e-3 with 2 Newton iterations per step.
- Practical options: `xmhd_mfnk=T`, or an input option for the Newton Jacobian-refresh interval (`nksolver%up_freq` is hard-coded to 4 in `xmhd_run`).

**Spheromak tilt** (serial, p=2, matrix-free Newton as in the example, `pb_pre` type 1 with inner GMRES(8) + block-Jacobi LU (1024 blocks) on P, 3 steps; `results/si_spheromak.txt`):

| configuration | dt requested | outer its/step | Newton | dt taken | wall (3 steps) | peak RSS |
|---|---|---|---|---|---|---|
| example ILU preconditioner (reference, §5.7) | 1e-6 | 117–136 | 2 | 1e-6 | — | 1.36 GB |
| PB type 1, s=0 | 1e-6 | 91–104 | 2 | 1e-6 | 622 s | 2.3 GB |
| PB type 1, **s=1** | 1e-6 | **15** | 1 | 1e-6 | **130 s** | 2.3 GB |
| PB type 1, **s=1** | 3e-6 | **31–32** | 1 | 3e-6 | 232 s | 2.3 GB |
| PB type 1, s=0 | 1e-5 | 199–261 | 2 | **cut → 2.5e-6 … 3.2e-6** | 2376 s | 2.7 GB |
| PB type 1, **s=1** | 1e-5 | **115–117** | 2 | **1e-5 (no cut)** | **762 s** | 2.3 GB |
| PB type 1, s=1, inner GMRES(20) on P | 1e-5 | 66–69 | 2 | 1e-5 | 1065 s | 2.4 GB |

**For the first time in this work, the spheromak runs at dt = 1e-5 without dt reduction.**
- Per step it is about 3× cheaper than the s=0 run, which itself only manages about 2.8e-6, so about 10× more physical time per wall-clock second.
- At 1e-6 it is about 4.8× faster than s=0 with the same preconditioner.
- Iterations still grow with dt (15 → 32 → 115). That is consistent with the weak inner solve on P: block-Jacobi LU on a stiff, anisotropic V block and the H(curl) B block. A stronger inner solve (GMRES(20) instead of GMRES(8) on P) cuts the outer iterations to 66–69, but raises the wall time to 1065 s: the inner solver on P is now the bottleneck.

**Not shown.** Spheromak accuracy: a tilt growth-rate comparison needs about 0.3 ms of physical time (2+ hours for the cheapest s=1 run, and a converged small-dt reference is out of reach in this sandbox). The cube gives the accuracy evidence; the spheromak growth rate vs dt and s is the most important open check.

**Verdict.** The SI increment term is the most effective large-CFL measure found in this work. It is cheap, opt-in, consistent, and it makes the existing type-1 preconditioner dt-robust. The remaining limits are (i) the nonlinear Newton Jacobian lag at CFL ≳ 100 (pre-existing), (ii) the inner solver on P, and (iii) the unmeasured accuracy cost on realistic problems.

---

## 8. Recommendations / future work (ranked)
0. **Semi-implicit increment term (§7b, implemented):** the accuracy cost on the spheromak tilt is now measured ([`RESULTS_M2Max.md`](RESULTS_M2Max.md) §6). See RESULTS §10 for the ranked follow-ups from that campaign (SI with full-Jacobian LU, a scalable inner solve on P, accuracy-based dt guidance, MUMPS for large blocks). Tune the inner solver on P at large dt (stronger V-block solve; auxiliary-space solver for the B block, item 1). Add an input option for the Newton Jacobian-refresh interval, or document `xmhd_mfnk=T`, for nonlinear runs at CFL ≳ 100. Add a regression test (linear cube, `xmhd_si=1`, PB type 1, iteration cap).
1. **Auxiliary-space preconditioner for S_X** (Hiptmair–Xu / AMS-style). Hiptmair smoothing on the H(curl) block, plus nodal vector-Lagrange and scalar-Lagrange auxiliary corrections. OFT already has H1 and Lagrange hierarchies and discrete gradient operators. The B_g (grad H1) part of the composite B already exposes the gradient subspace, which may simplify this.
2. **Type 5: matrix-free exact Schur.** The V solve is FGMRES on the exact S_V applied matrix-free (inner D_X solves), preconditioned by S_X or by continuum S_V. Robust, but the inner–outer cost needs measuring.
3. **Port the continuum S_V path to `xmhd_lag.F90`** (nodal B, as in NIMROD), where §6.1 suggests it should work.
4. **Reduce S_X fill:** drop small entries, lump J_XV·diag⁻¹·J_VX element-wise, or keep only the Alfvénic (B0·∇) part. Check spectral equivalence with `schur_analysis.py`.
5. **MPI for type 4:** replace the global `lge` fold with distributed row ownership.
6. **Tests and docs (M3):** `test_alfven.py` cases for `pb_pre` types 1 and 4 with an iteration cap; document `<pb_pre>`, `kcycle` and `cscale`; add type 4 to the `pb_type` doc comment.
7. **Hall/whistler:** follow Chacón 2024 (electron-MHD block treatment).
8. Diagnose the ILU(0) zero pivot on S_X: try removing the placeholder rows, or ILU with a diagonal shift.

---

## 9. Reproduction
- **Build:** see the `oft-build` skill. The system-library build used here: MPI + OpenMP, SuperLU, ARPACK, in `builds/build_release`. Configure with `-DOFT_BUILD_EXAMPLES=ON`, then `make MUG_sph_tilt` for the spheromak executable. Always use `OMP_NUM_THREADS=2`.
- **Scripts** (`scripts/`). Paths resolve relative to the repo; set `OFT_BUILD_DIR` for another build tree. Each script writes into `./runs/<tag>` or `./runs_sph/<tag>` under the current directory.
  - Cube: `python3 run_alfven.py --tag T --kind {lu|bj|mg|pb_lu|pb_bj|pb_mg} [--pbtype 1-4] [--dt 1.E-4] [--nsteps 3] [--nlocal 8 --overlap] [--its 2] [--inner_its 5] [--inner lu|ilu] [--kcycle 2] [--cscale] [--alfven 1.0] [--nlevels 2 --minlev 1]`
  - Spheromak: `python3 run_sph.py --tag T --kind {bj|pb_bj} --pbtype 4 --nlocal 1024 --inner lu --its 8 --inner_its 8 --dt 1.e-5 --nsteps 3`. Run it alongside `./watchdog.sh 12000000` to kill a run that exceeds 12 GB.
  - Matrix dump plus offline analysis: `python3 run_alfven.py --tag T --kind pb_lu --pbtype 1 --dump --nsteps 2 --dt 1.E-5`, then `python3 schur_analysis.py runs/T`.
  - Semi-implicit study (§7a): dump at dt = 1.E-5, 2.E-5, 4.E-5 (as above, tags `sidump_<dt>`), then `python3 si_analysis.py runs/sidump_1.E-5 runs/sidump_2.E-5 runs/sidump_4.E-5 [--scan_s 1,1.25]`. It takes about 10–20 minutes on 2 cores (dense eigensolves of size about 2600).
  - SI increment term (§7b): add `--si 1` to `run_alfven.py` or `run_sph.py` (and `--linear` / `--mfnk` as needed), for example `python3 run_alfven.py --tag T --kind pb_lu --pbtype 1 --dt 1.E-4 --nsteps 6 --si 1`. In an input file: `xmhd_si=1.` in `xmhd_options`, together with `<pb_pre type="1">` in the XML.
  - `python3 summarize.py runs runs_sph` prints the per-run summary lines (`results/runs_summary.txt`).
- **Resources:** the sandbox had about 14 GB RAM and 2 cores. The spheromak LU block-Jacobi with nlocal ≤ 64 does not fit.
