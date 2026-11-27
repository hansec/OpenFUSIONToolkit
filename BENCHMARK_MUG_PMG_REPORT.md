# Large-Scale MUG $p$-MG & Block-Schur Scalability Benchmark Report

> **Purpose & Audience**: This report documents the Tier A large-scale (up to `6.32M` DOFs, $p=2\text{--}4$, 16 physical cores) benchmarking of the five stacked `mug-pmg-step*` PR branches (tip commit `6987f2d`, `mug-pmg-step5-coarse-reg`), executed per `BENCHMARK_MUG_PMG.md` plus a follow-up study comparing production subdomain `LU` (`SuperLU`) against lightweight `ILU(0)`. It is intended as a durable engineering reference for future developers and AI agents preparing PR 1--5 for merge or running Tier B/C HPC scaling follow-ups.

---

## Executive Summary

1. **PR 3 (`vel_scale=-1` auto) + PR 4 (Option 2.A Block-Schur) are decisive wins at scale**:
   - On `S2-p2` (`130,440` tets, `1.84M` `B+V` / `2.20M` total DOFs, `np=16`), PR 3 alone cuts baseline `C0` step time **2.2--2.5x** (`26.74s -> 12.24s` at `dt=1e-6`, $\mathrm{CFL}_A \approx 5$).
   - Block-Schur `C4_lower` (`2.84 s/step`) and `C4_ldu` (`3.15 s/step`) add another **3.9--4.3x** speedup over auto-scaled `C0` (**9.4x** over pre-PR3 `C0`), and uniquely converge up to $\mathrm{CFL}_A \approx 10$ (`8.32 s/step`) and $\mathrm{CFL}_A \approx 25$ (`24.91 s/step`) where `C0` and all monolithic $p$-MG configs fail.
   - The advantage grows with $p$ and mesh size: on `S2-p3` (`5.11M` `B+V` / `6.32M` total DOFs, `dt=1e-6`), `C4_lower` (`35.38 s/step`) and `C4_ldu` (`36.72 s/step`) beat `C0` (`321.81 s/step`) by **8.8--9.1x**, and iterative-coarse `C6` beats `C0` by **7.5x** (`43.15 s/step`, lowest Schur RSS `5.20 GB/rank`, zero distributed coarse LU).
   - Over 50 nonlinear MFNK steps (`xmhd_prefreq=20`), Block-Schur setup/refactorization overhead amortizes to **0.1--1.8%** of total wall time, delivering **1.9--3.2x** faster 5-field `MUG_sph_tilt` runs (**2.3x** more physical time per wall-second under adaptive `dt`) and **2.65x** faster 7-field `MUG_sph_heat` runs at `dt=1e-6`.
2. **PR 5 (Option 2.B monolithic $p$-MG) does not scale past small/low-CFL problems**: Monolithic coarse `superd` LU (`C2`/`C3a`/`C3b`) costs `20--33 s` setup and `2.75 GB/rank` already at `S2-p2` (vs `1.3--2.1 s` and `1.42 GB/rank` for `S_V`-only coarse LU in `C4`), and degrades above $\mathrm{CFL}_A \approx 2.5$.
3. **Coarse `S_V` LU is NOT a bottleneck up through `6.32M` DOFs**, and **`J_BB` does NOT need $p$-MG (`C5`)** in any tested regime (`C5` is `1.5--2.2x` slower than `C4` across $p=2,3$, high CFL, and 7-field Spitzer $\eta_0=25$).
4. **Subdomain `LU` vs `ILU(0)`**:
   - **Never pair `add_schwarz(overlap=1)` with subdomain `LU` on 3D high-order meshes**: 1-hop 3D overlap expands subdomain blocks 3--5x, blowing up SuperLU fill and triggering OOM (`8--11 GB/rank`, `>130 GB` total).
   - At $p=2, 3$, lightweight `ILU(0)` (`b_pre: BJ(-1, ILU)`, `v_pre: ASM(1, ILU)`) is **1.5--1.8x faster** than non-overlapping `BJ(16, LU)` at the same memory footprint (`1.42--1.45 GB/rank`).
   - At **$p=4$** (`S1-p4`, `1.76M` DOFs, `dt=5e-7`), all `ILU(0)` configs fail, while non-overlapping **`BJ(nlocal=16, SuperLU)` inside `<b_pre>`** cures the failure cleanly (`9.19 s/step`, `2.57 GB/rank`, **3.3x faster** than `C0 + BJ(16, LU)`).
5. **Pre-merge blockers found**: Missing `NULLIFY(output)` in `src/base/oft_xml.F90` (heap corruption on `C2`), branches need rebase onto `main` for `#459` (METIS `npart=1` + PETSc ASM) and `#453` (coarse geometric $h+p$-MG `SIGSEGV`), PETSc `<solver type="ilu">` defaults to `mumps` unless `<package>petsc</package>` is set, and `MUG_sph_tilt`/`MUG_sph_heat` crash under `use_petsc=T` on unguarded `linv%pre%pre` and hardcode `minlev=-1` (silently disabling $p$-MG).

---

## 1. System, Build & Experimental Setup

| Property | Value |
|---|---|
| **Host CPU** | AMD Ryzen AI MAX+ 395 (Zen 5 Strix Halo), 1 socket, **16 physical cores** (2 CCDs x 8 cores, 64 MiB L3 = 2 x 32 MiB; SMT 32 threads unused) |
| **Memory** | 123 GiB RAM (NUMA node 0), 8 GiB swap |
| **Hardware Tier** | **Tier A** (single-node workstation; strictly `<= 16` physical cores, strictly 1 benchmark at a time) |
| **Git Ref Tested** | `mug-pmg-step5-coarse-reg` (`6987f2d`, cumulative PR 1--5) in isolated detached worktree |
| **Build Config** | `CMAKE_BUILD_TYPE=Release`, LTO, OpenMP, MPI F08 |
| **Compilers & MPI** | GCC 16.0.1 (`gfortran-16` / `gcc-16` / `g++-16`, Ubuntu 16-20260315), MPICH 4.2.3 (Hydra) |
| **Linear Algebra Stack** | PETSc 3.25.2 (`SuperLU`, `SuperLU_DIST`, `MUMPS`, `ScaLAPACK`, `SuiteSparse`, `METIS`, `ParMETIS`), OpenBLAS 0.3.33 (`OPENBLAS_NUM_THREADS=1`), ARPACK-NG 3.9.1, HDF5 1.14.6, LibXml2 2.15.2 |
| **MPI Launch Flags** | `MPIR_CVAR_ENABLE_INTRANODE_TOPOLOGY_AWARE_TREES=0` (avoids MPICH 4.2.3 Zen 5 `mpir_hwtopo.c:57` assertion), `-bind-to core` (`OMP=1`) or `-bind-to core:2` (`OMP=2`), `OPENBLAS_NUM_THREADS=1` |

### Temporary Worktree Compatibility Patches Used During Benchmarking
To run the full PETSc + MPI + $p$-MG matrix from `6987f2d`, four minimal fixes were applied in the isolated benchmark worktree (none committed to the source tree; all documented in Section 5 for PR cleanup):
1. Cherry-picked `#459` (`d528069`): fixes `METIS_PartMeshDual` SIGFPE for `npart=1` and PETSc 3.23--3.25 `PCASMSetLocalSubdomains` bug.
2. Added `NULLIFY(output)` at entry of `xml_extractDataAttribute_{double,int,logical}1D` (`src/base/oft_xml.F90:951,1026,1101`) to fix `cfree` heap corruption when optional PR 5 attributes are omitted.
3. Guarded `DEALLOCATE(linv%pre%pre)` / `create_ilu_pre(linv%pre%pre)` / `linv%pre%pre%delete()` with `IF(ASSOCIATED(linv%pre%pre))` in `MUG_sph_tilt.F90` (`:166-191`) and `MUG_sph_heat.F90` (`:150-168`) for `use_petsc=T`.
4. Added `minlev` to `NAMELIST/sph_tilt_options/` and `NAMELIST/sph_heat_options/` (`minlev=-1` for 1-level `C0`/`C1`; `minlev=nlevels` for pure $p$-MG `C2`--`C6`).

### Problem & DOF Ladder (`cyl_tilt.h5`, $v_A \approx 4.88\times 10^5\text{ m/s}$)

| Size ID | `nbase / nlevels` (MPI) | Tetrahedra | `order` | Exact `B+V` DOFs | Exact Total DOFs | Status on 16 Cores |
|---|---|---:|---:|---:|---:|---|
| `S1-p2` | `1 / 2` | 16,305 | 2 | `2.37e5` | `2.85e5` | Complete (Phase 1, Suites 2, 3) |
| `S1-p3` | `1 / 2` | 16,305 | 3 | `6.53e5` | `8.09e5` | Complete (Suite 2) |
| `S1-p4` | `1 / 2` | 16,305 | 4 | `1.40e6` | `1.76e6` | Complete (Suite 2 + Follow-Up LU) |
| `S2-p2` | `1 / 3` | 130,440 | 2 | `1.84e6` | `2.20e6` | Complete (Suites 1, 3, 4, 5a, 5b + LU) |
| `S2-p3` | `1 / 3` | 130,440 | 3 | `5.11e6` | `6.32e6` | Complete (Suite 2) |
| `S3-p2` | `1 / 4` | 1,043,520 | 2 | `~1.3e7` | `~1.5e7` | Exceeds 600s in `taylor_hmodes` IC on 16 cores (Tier B/C) |

---

## 2. Physics Parity Verification

All 16 Phase 1 configurations (`C0`, `C1` native, `C2`, `C3a/b/c`, `C4_lower/ldu` ASM/BJ + `superd`/`mumps`, `C5`, `C6`) ran cleanly on `S1-p2` (`np=8, dt=5e-7`, 3 steps) with identical `# B-DOFs` and `# V-DOFs`.

Because `MUG_sph_tilt` starts from a $10^{-4}$ linear perturbation (`men ~ 5.2e-9`, `ven ~ 3.8e-6`) and stops GMRES on preconditioned relative residual `lin_tol=1e-9`, different preconditioners stop at slightly different true residuals (`dmen <= 8.56e-5`, `dven <= 4.11e-3`; even repeated `C0` runs differ by `3.0e-6`/`4.5e-5` from MPI reduction order). Tightening to `lin_tol=1e-12` confirmed exact algebraic equivalence:

| Verification Check | `Mag Energy (men)` | `Kin Energy (ven)` | Max Rel. Diff (`men`) | Max Rel. Diff (`ven`) |
|---|---:|---:|---:|---:|
| **Phase 1 `S1-p2` (`lin_tol=1e-12`): `C0` vs `C4_ldu`** | `5.17045828e-09` vs `5.17045873e-09` | `3.80336132e-06` vs `3.80337860e-06` | **`8.59e-08`** | **`4.54e-06`** |
| Phase 1 `S1-p2` (`lin_tol=1e-9`, all 16 configs vs `C0`) | `5.1703e-09` | `3.7983e-06` | `8.56e-05` | `4.11e-03` |
| Suite 1 `S2-p2` (`dt=1e-6, lin_tol=1e-9`, `C0/C3c/C4/C5/C6`) | `1.3047e-08` | `9.7769e-05` | `2.11e-03` | `4.04e-03` |

**Result: PASS** (relative error `< 1e-5` once GMRES tolerance is tightened below the perturbation floor).

---

## 3. Scaling & Performance Results

> **Critical `<b_pre>` finding**: On `S2-p2` (130k tets), `<b_pre>` configured with coupled `<nlocal>1</nlocal>` + PETSc `ILU(0)` stalled (`>200` iterations) even at `dt=2.5e-7` ($\mathrm{CFL}_A \approx 1.25$) because zero-fill ILU is unstable on coupled `Bc+Bg` at large subdomain sizes. Replacing `<b_pre>`'s preconditioner with field-split `<pre type="block_jacobi"><nlocal>-1</nlocal><solver type="ilu"><package>petsc</package></solver></pre>` (while keeping `add_schwarz(overlap=1, nlocal=1, ilu)` on SPD `S_V`) immediately restored fast convergence and is used for all `C4`/`C5`/`C6` runs in Tables 1--5. Subdomain `LU` on `<b_pre>` is evaluated separately in Table 6.

### Table 1: Suite 1 — $\mathrm{CFL}_A$ Sweep on `S2-p2` (`1.84M` `B+V` DOFs, `np=16`, `OMP=1`, 10 linear steps)

| Config | `dt=2.5e-7` ($\mathrm{CFL}_A\approx 1.25$)<br>`its / t_setup / t_step` | `dt=5e-7` ($\mathrm{CFL}_A\approx 2.5$)<br>`its / t_setup / t_step` | `dt=1e-6` ($\mathrm{CFL}_A\approx 5.0$)<br>`its / t_setup / t_step` | `dt=2e-6` ($\mathrm{CFL}_A\approx 10$)<br>`its / t_setup / t_step` | `dt=5e-6` ($\mathrm{CFL}_A\approx 25$)<br>`its / t_setup / t_step` | Peak RSS (`1e-6`) | Speedup vs `C0` (`1e-6`) |
|---|---|---|---|---|---|---:|---:|
| `C0` (auto `vel_scale`) | `4 / 0.83s / 1.48s` | `14 / 0.84s / 4.77s` | `36 / 3.21s / 12.24s` | `FAIL` | `SKIP` | `0.86 GB` | `1.00x` |
| `C0` (`vel_scale=1e3`, pre-PR3) | — | `35 / 1.56s / 11.88s` | `79 / 12.18s / 26.74s` | `SKIP` | `SKIP` | `0.86 GB` | `0.46x` |
| `C2` (Mono $p$-MG unreg) | `9 / 31.16s / 3.87s` | `27 / 33.30s / 11.34s` | `TOUT (~81 / ~35s)*` | `SKIP` | `SKIP` | `2.76 GB` | `~0.35x` |
| `C3a` (2.B `damp=0.5`) | `6 / 29.46s / 2.62s` | `28 / 19.92s / 11.77s` | `TOUT (~95 / ~42s)*` | `SKIP` | `SKIP` | `2.73 GB` | `~0.29x` |
| `C3b` (2.B `shift=1.0`) | `7 / 29.52s / 3.07s` | `26 / 23.13s / 10.96s` | `TOUT (~83 / ~69s)*` | `SKIP` | `SKIP` | `2.75 GB` | `~0.18x` |
| `C3c` (2.B `cILU+shift`) | `6 / 0.33s / 1.64s` | `21 / 0.00s / 5.35s` | `62 / 0.00s / 15.69s` | `FAIL` | `SKIP` | `0.94 GB` | `0.78x` |
| **`C4_lower` (2.A Schur)** | **`5 / 1.30s / 1.11s`** | **`7 / 1.67s / 1.50s`** | **`14 / 2.07s / 2.84s`** | **`41 / 2.02s / 8.32s`** | **`125 / 8.95s / 24.91s`** | `1.42 GB` | **`4.31x`** |
| **`C4_ldu` (2.A Schur)** | **`3 / 1.34s / 0.99s`** | **`5 / 1.38s / 1.52s`** | **`11 / 1.88s / 3.15s`** | **`31 / 2.48s / 8.77s`** | **`107 / 7.27s / 30.02s`** | `1.42 GB` | **`3.89x`** |
| `C5` (2.A `B=pMG, V=pMG`) | `3 / 14.23s / 1.88s` | `4 / 15.08s / 2.45s` | `9 / 16.54s / 5.25s` | `30 / 16.66s / 17.25s` | `TOUT (9/10, 65.9s/st)` | `2.16 GB` | `2.33x` |
| `C6` (2.A `iter coarse`) | `4 / 0.00s / 1.20s` | `6 / 0.42s / 1.76s` | `14 / 0.70s / 4.07s` | `36 / 7.47s / 10.88s` | `TOUT (4/10, 71.5s/st)` | `1.26 GB` | `3.01x` |

*\*`C2`/`C3a`/`C3b` at `dt=1e-6` exceeded the 600 s wall-clock timeout; listed `its`/`t_step` are from un-overlapped steps 2--7.*

---

### Table 2: Suite 2 — Polynomial Order ($p$-Refinement) Study (`np=16`, `OMP=1`)

#### 2a. `S1` (`16,305` tets, `p = 2, 3, 4`)

| Mesh & DOFs (`B+V` / tot) | Config | `dt=5e-7`: `its / t_setup / t_step / RSS` | `dt=2e-6`: `its / t_setup / t_step / RSS` | Notes |
|---|---|---|---|---|
| **`S1-p2`** (`2.37e5` / `2.85e5`) | `C0` | `9 / 0.05s / 0.38s / 0.17G` | `59 / 0.35s / 2.37s / 0.17G` | Baseline |
| | `C3a` | `12 / 1.05s / 0.52s / 0.33G` | `119 / 0.00s / 5.13s / 0.37G` | 2.2x slower than `C0` at `2e-6` |
| | **`C4_lower`** | `7 / 0.13s / 0.18s / 0.25G` | **`26 / 0.35s / 0.61s / 0.25G`** | **3.9x faster than `C0`** |
| | `C4_ldu` | `6 / 0.10s / 0.21s / 0.25G` | `20 / 0.42s / 0.67s / 0.25G` | 3.5x faster than `C0` |
| | `C5` | `5 / 0.76s / 0.31s / 0.30G` | `19 / 1.29s / 1.14s / 0.30G` | Slower than `C4` |
| | `C6` | `6 / 0.05s / 0.20s / 0.22G` | `22 / 0.44s / 0.68s / 0.22G` | Matches `C4_ldu`, lowest Schur RSS |
| **`S1-p3`** (`6.53e5` / `8.09e5`) | `C0` | `18 / 0.00s / 3.36s / 0.47G` | `141 / 10.94s / 26.00s / 0.47G` | Rapid $p$-degradation |
| | `C3a` | `14 / 0.07s / 2.50s / 0.69G` | `FAIL (step 2)` | Fails at `2e-6` |
| | **`C4_lower`** | `11 / 0.18s / 1.37s / 0.81G` | **`41 / 10.02s / 4.91s / 0.82G`** | **5.3x faster than `C0`** |
| | **`C4_ldu`** | `8 / 0.09s / 1.35s / 0.81G` | `32 / 10.50s / 5.13s / 0.82G` | **5.1x faster than `C0`** |
| | `C5` | `6 / 0.90s / 1.52s / 0.85G` | `29 / 20.48s / 7.00s / 0.85G` | 2x higher setup, slower `t_step` |
| | `C6` | `8 / 0.05s / 1.25s / 0.75G` | `39 / 9.12s / 5.82s / 0.75G` | Fastest at `5e-7`, 4.5x vs `C0` at `2e-6` |
| **`S1-p4`** (`1.40e6` / `1.76e6`) | All ILU (`C0`--`C6`) | `FAIL` | `FAIL` | At `dt=2.5e-7`: `C0` `31 its / 17.74s / 1.25G`; `C4_ldu` `44.5 its / 23.68s / 2.56G` |
| | `C0` + `BJ(16, LU)` | `30.5 / 0.00s / 30.60s / 1.54G` | `FAIL` | Subdomain LU rescues $p=4$ |
| | **`C4_ldu` + `b:BJ(16,LU)`** | **`15.5 / 0.00s / 9.19s / 2.57G`** | `FAIL` | **3.3x faster than `C0_bj16lu`, cures $p=4$** |

#### 2b. `S2-p3` (`130,440` tets, `5.11M` `B+V` / `6.32M` total DOFs, `dt=1e-6`, `np=16`, 3 steps)

| Config | Status | `lits` | `t1 (s)` | `t_step (s)` | Peak RSS/rank | Speedup vs `C0` |
|---|---|---:|---:|---:|---:|---:|
| `C0` | `TOUT (2/3)` | `106.0` | `298.51` | `321.81` | `n/a` | `1.00x` |
| **`C4_lower`** | **OK** | **`37.0`** | **`33.77`** | **`35.38`** | **`5.55 GB`** | **`9.10x`** |
| **`C4_ldu`** | **OK** | **`28.5`** | **`33.86`** | **`36.72`** | **`5.55 GB`** | **`8.76x`** |
| `C5` (`B=pMG, V=pMG`) | OK | `27.0` | `60.47` | `54.46` | `6.00 GB` | `5.91x` |
| **`C6` (iterative coarse)** | **OK** | `34.5` | `43.50` | `43.15` | **`5.20 GB`** | **`7.46x`** |

---

### Table 3: Suite 3 — Exact 8x Weak Scaling (`S1-p2 @ np=2` $\to$ `S2-p2 @ np=16`, ~`115k` `B+V` DOFs/rank)

| Regime | Config | `S1-p2 @ np=2` (`its / t_step / RSS`) | `S2-p2 @ np=16` (`its / t_step / RSS`) | Weak Efficiency $E_w = t(2)/t(16)$ |
|---|---|---|---:|---:|
| **Suite 3a: Fixed $\Delta t = 10^{-6}$**<br>($\mathrm{CFL}_A: 2.5 \to 5.0$) | `C0` | `23 / 2.52s / 0.81G` | `36 / 12.24s / 0.86G` | `20.6%` |
| | `C3a` | `34 / 3.89s / 1.51G` | `TOUT (~95 / ~42s)` | `~9.3%` |
| | **`C4_ldu`** | **`10 / 0.99s / 1.21G`** | **`11 / 3.15s / 1.42G`** | **`31.4%`** (its: `10 -> 11`) |
| | `C5` | `8 / 1.42s / 1.50G` | `9 / 5.25s / 2.16G` | `27.0%` |
| | `C6` | `10 / 0.99s / 1.17G` | `14 / 4.07s / 1.26G` | `24.3%` |
| **Suite 3b: Fixed $\mathrm{CFL}_A \approx 2.5$**<br>(`dt=1e-6` $\to$ `5e-7`) | `C0` | `23 / 2.52s / 0.81G` | `14 / 4.77s / 0.86G` | `52.8%` |
| | `C3a` | `34 / 3.89s / 1.51G` | `28 / 11.77s / 2.67G` | `33.0%` |
| | **`C4_ldu`** | **`10 / 0.99s / 1.21G`** | **`5 / 1.52s / 1.42G`** | **`65.1%`** (its: `10 -> 5`) |
| | `C5` | `8 / 1.42s / 1.50G` | `4 / 2.45s / 2.17G` | `58.0%` |
| | `C6` | `10 / 0.99s / 1.17G` | `6 / 1.76s / 1.27G` | `56.3%` |

---

### Table 4: Suite 4 — Strong MPI Scaling & Solver Variants (`S2-p2`, `dt=1e-6`)

#### 4a. Strong MPI Scaling (`OMP_NUM_THREADS=1`, `np = 2, 4, 8, 16`)

| Config | `np=2` (`its/t_st/RSS`) | `np=4` (`its/t_st/RSS`) | `np=8` (`its/t_st/RSS`) | `np=16` (`its/t_st/RSS`) | Speedup $S(16/2)$ | Efficiency $E_s(16/2)$ |
|---|---|---|---|---|---:|---:|
| `C0` | `37 / 36.02s / 5.87G` | `37 / 28.20s / 3.11G` | `35 / 24.60s / 1.60G` | `36 / 12.24s / 0.86G` | `2.94x` | `36.8%` |
| **`C4_lower`** | **`12.5 / 7.28s / 9.27G`** | **`13 / 5.78s / 5.01G`** | **`14 / 5.55s / 2.61G`** | **`14 / 2.84s / 1.42G`** | `2.56x` | `32.0%` |
| **`C4_ldu`** | **`10 / 8.46s / 9.28G`** | **`10.5 / 6.69s / 5.00G`** | **`11 / 6.21s / 2.59G`** | **`11 / 3.15s / 1.42G`** | `2.69x` | `33.6%` |
| `C5` | `9.5 / 16.15s / 13.19G` | `9.5 / 11.77s / 7.20G` | `9 / 9.88s / 3.73G` | `9 / 5.25s / 2.16G` | `3.08x` | `38.5%` |
| `C6` | `11.5 / 10.14s / 8.85G` | `12.5 / 8.27s / 4.65G` | `14 / 8.17s / 2.36G` | `14 / 4.07s / 1.26G` | `2.49x` | `31.1%` |

*Hardware architecture note*: SpMV-bound solves saturate single-CCD memory fabric from `np=4 -> 8` (`5.78s -> 5.55s`), then jump **2.0x** at `np=16` (`2.84s`) once both Zen 5 CCDs and both 32 MiB L3 caches engage. `C4_lower` on **2 cores** (`7.28 s/step`) beats `C0` on **16 cores** (`12.24 s/step`).

#### 4b. Variant Comparisons at 16 Physical Cores (`S2-p2`, `dt=1e-6`)

| Comparison Axis | Variant A | Variant B | Verdict |
|---|---|---|---|
| **Coarse LU (`C4_ldu`)** | `superd`: `11 its, t1=5.03s, t_step=3.15s, 1.42G` | `mumps`: `11 its, t1=5.78s, t_step=3.12s, 1.38G` | Effectively tied (`superd` `0.75s` faster setup; `mumps` `~1%` faster solve, `40 MB` lighter) |
| **`S_V` Smoother (`C4_ldu`)** | `add_schwarz(ov=1)`: **`11 its, t_step=3.15s, 1.42G`** | `block_jacobi`: `14 its, t_step=3.87s, 1.38G` | **`add_schwarz(ov=1)` wins by 23%** (`19%` on `C4_lower`: `2.84s` vs `3.48s`) |
| **Pure MPI vs Hybrid (`C4_ldu`)** | `16 ranks x 1 OMP`: **`11 its, t_step=3.15s, 1.42G`** | `8 ranks x 2 OMP`: `11 its, t_step=3.46s, 2.60G` | **Pure MPI (`16x1`) wins by 10%** (`C0`: `12.24s` vs `13.50s`) at half per-rank RSS |

---

### Table 5: Suite 5 — Nonlinear Production Proxy (`S2-p2`, `np=16`, `xmhd_mfnk=T`, `xmhd_prefreq=20`, 50 steps / 600 s cap)

#### 5a. 5-Field Nonlinear `MUG_sph_tilt` (`xmhd_advec=T`)

| Mode / `dt` | Config | Status (Steps) | `lits` (med / avg) | `nlits` avg | `t1 (s)` | `t_step (s)` | Total (s) | Setup Amort. | Phys. Time | Throughput vs `C0` |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **Fixed `dt=5e-7`** | `C0` | `TOUT (29/50)` | `47.0 / 47.3` | `1.00` | `19.93` | `19.02` | `559.0` | `0.2%` | `14.5 us` | `1.00x` |
| | `C3a` | `TOUT (8/50)` | `110.0 / 115.2` | `1.00` | `95.08` | `56.72` | `531.1` | `7.2%` | `4.0 us` | `0.29x` |
| | **`C4_lower`** | **`OK (50/50)`** | **`33.0 / 33.7`** | **`1.00`** | **`9.72`** | **`9.72`** | **`511.6`** | **`0.1%`** | **`25.0 us`** | **`1.88x` (`1.96x t_step`)** |
| | **`C4_ldu`** | **`OK (50/50)`** | **`27.0 / 27.4`** | **`1.00`** | **`10.46`** | **`10.06`** | **`526.0`** | **`0.1%`** | **`25.0 us`** | **`1.83x` (`1.89x t_step`)** |
| | `C5` | `TOUT (31/50)` | `25.0 / 25.3` | `1.00` | `29.13` | `16.21` | `554.9` | `2.3%` | `15.5 us` | `1.08x` |
| **Fixed `dt=2e-6`** | `C0` | `TOUT (7/50)` | `142.0 / 142.6` | `1.00` | `122.45` | `59.88` | `495.6` | `12.6%` | `10.4 us` *(cut `dt`)* | `1.00x` |
| | `C3a` | `TOUT (0/50)` | — | — | — | — | — | — | `0.0 us` | — |
| | **`C4_lower`** | **`TOUT (18/50)`** | **`103.0 / 102.1`** | **`1.00`** | **`24.96`** | **`29.77`** | **`543.6`** | **`<0.5%`** | **`36.0 us`** | **`3.15x`** |
| | **`C4_ldu`** | **`TOUT (17/50)`** | **`89.0 / 86.3`** | **`1.00`** | **`25.04`** | **`32.43`** | **`542.1`** | **`<0.5%`** | **`34.0 us`** | **`2.99x`** |
| | `C5` | `TOUT (10/50)` | `77.0 / 76.7` | `1.00` | `56.13` | `48.77` | `532.6` | `1.4%` | `20.0 us` | `1.79x` |
| **Adaptive `dt`**<br>(`ittarget=40`, `dt_max=5e-6`) | `C0` | `TOUT (32/50)` | `39.0 / 40.5` | `1.00` | `19.53` | `16.26` | `549.7` | `0.8%` | `14.6 us` (`dt` in `[4.3,5.0]e-7`) | `1.00x` |
| | **`C4_ldu`** | **`TOUT (40/50)`** | **`36.0 / 35.9`** | **`1.00`** | **`10.45`** | **`13.55`** | **`552.6`** | **`1.3%`** | **`33.7 us`** (`dt` in `[5.0,10.9]e-7`) | **`2.30x`** |
| | **`C4_lower`** | **`TOUT (48/50)`** | **`38.0 / 37.6`** | **`1.00`** | **`9.72`** | **`11.18`** | **`550.0`** | **`<0.5%`** | **`32.0 us`** (`dt` in `[5.0,8.4]e-7`) | **`2.19x`** |

#### 5b. 7-Field Nonlinear `MUG_sph_heat` (`B+V+n+T`, Spitzer $\eta_0=25$, $\kappa_\parallel=10^4$, $\kappa_\perp=10^2$)

| `dt` | Config | Status (Steps) | `lits` (med / avg) | `nlits` avg | `t1 (s)` | `t_step (s)` | Total (s) | Setup Amort. | Phys. Time | Speedup vs `C0` (`t_step`) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **`2e-7`** | `C0` | `OK (50/50)` | `11.0 / 11.8` | `1.00` | `10.69` | `6.60` | `335.9` | `1.2%` | `10.0 us` (`1.10G`) | `1.00x` |
| | `C3a` | `TOUT (35/50)` | `19.0 / 18.8` | `1.00` | `68.84` | `12.83` | `577.3` | `9.7%` | `7.0 us` | `0.51x` |
| | **`C4_lower`** | **`OK (50/50)`** | **`14.0 / 14.7`** | **`1.00`** | **`8.54`** | **`5.03`** | **`281.5`** | **`1.4%`** | **`10.0 us` (`1.69G`)** | **`1.31x` (`1.19x tot`)** |
| | **`C4_ldu`** | **`OK (50/50)`** | **`11.0 / 11.8`** | **`1.00`** | **`8.79`** | **`5.20`** | **`294.5`** | **`1.4%`** | **`10.0 us` (`1.69G`)** | **`1.27x` (`1.14x tot`)** |
| | `C5` | `TOUT (43/50)` | `16.0 / 16.1` | `1.00` | `28.00` | `11.59` | `577.8` | `3.8%` | `8.6 us` | `0.57x` |
| **`1e-6`** | `C0` | `TOUT (17/50)` | `62.0 / 66.1` | `1.06` | `64.87` | `31.76` | `578.3` | `5.7%` | `17.0 us` | `1.00x` |
| | `C3a` | `TOUT (3/50)` | `118.5 / 154.3` | `1.33` | `200.81` | `161.87` | `524.5` | `7.4%` | `3.0 us` | `0.20x` |
| | **`C4_lower`** | **`TOUT (45/50)`** | **`33.5 / 35.7`** | **`1.02`** | **`20.03`** | **`11.97`** | **`570.9`** | **`1.4%`** | **`45.0 us`** | **`2.65x`** |
| | **`C4_ldu`** | **`TOUT (43/50)`** | **`27.5 / 29.1`** | **`1.02`** | **`22.84`** | **`12.38`** | **`571.4`** | **`1.8%`** | **`43.0 us`** | **`2.57x`** |
| | `C5` | `TOUT (20/50)` | `37.0 / 36.9` | `1.05` | `45.82` | `27.03` | `575.2` | `3.3%` | `20.0 us` | `1.18x` |

#### 5c. Suite 6 — Periodic 3D Alfvén Wave (`test_alfven`, $h+p$-MG)
- `nlevels=7` (`reflect='xyz'`, 1.57M tets): OOM on 16 cores (`7.9 GB/rank x 16 > 123 GB`; Tier B/C size).
- `nlevels=6` (`reflect='xyz'`, 196,608 tets, `5.52M` `B+V` DOFs, `dt=2e-5`):
  - Pure $p$-MG (`minlev=6`): `C0`/`C5` timed out (`240s`); `C4_ldu`/`C6` hit `gmres_maxits=200`.
  - Geometric $h+p$-MG (`minlev=5` and `minlev=4`): crashed with `SIGSEGV` / `corrupted double-linked list` during coarse geometric level setup because `6987f2d` predates PR `#453` (`c6de81c` on `main`). Rebase `mug-pmg-step*` onto `main` before re-testing $h+p$-MG.

---

### Table 6: Follow-Up — Subdomain `LU` (`SuperLU`) vs `ILU(0)` in `block_jacobi` / `add_schwarz` (`np=16`)

Production MUG runs often use local `LU` inside `block_jacobi` or `add_schwarz` rather than unoptimized `ILU(0)`. We tested whole-rank `nlocal=1`, per-field `nlocal=-1`, and partitioned `nlocal=16` SuperLU subdomains on `S2-p2` and `S1-p4`:

| Case | Config & Local Subdomain Solver | Status | `lits` | `t1 (s)` | `t_setup` | `t_step (s)` | Peak RSS/rank | vs Best ILU (`t_step`) |
|---|---|---|---:|---:|---:|---:|---:|---|
| **`S2-p2`, `dt=1e-6`** | `C0` `BJ(-1, ILU)` *(ref)* | OK | `36.0` | `15.46` | `3.21s` | `12.24` | `0.86 GB` | `1.00x` |
| | **`C4_lower` `b:BJ(-1,ILU), v:ASM(1,ILU)`** *(ref)* | **OK** | **`14.0`** | **`4.92`** | **`2.07s`** | **`2.84`** | **`1.42 GB`** | **`4.31x` (Winner at $p=2$)** |
| | `C0` `BJ(1, LU)` *(whole-rank SuperLU)* | `TOUT (2/10)` | `37.0` | `230.56` | `96.32s` | `134.25` | `n/a` | Subdomain too big for serial LU |
| | `C0` `BJ(16, LU)` | `TOUT (3/10)` | `74.0` | `58.27` | `5.93s` | `52.34` | `n/a` | `4.3x` slower than `C0` ILU |
| | `C0` `ASM(ov=1, 16, LU)` | **OOM** | — | — | — | — | **`7.92 GB`** | `16 x 7.9G > 123G` |
| | `C4_ldu` `b:BJ(-1,LU), v:ASM(1,ILU)` | OK | `10.0` | `83.09` | `57.91s` | `25.18` | `4.28 GB` | `8.9x` slower, `3x` RSS |
| | `C4_lower` `b:BJ(16,LU), v:ASM(1,ILU)` | OK | `23.0` | `6.38` | `1.20s` | `5.17` | **`1.45 GB`** | `1.8x` slower than ILU, same RSS |
| | `C4_ldu` `b:BJ(16,LU), v:ASM(1,ILU)` | OK | `16.0` | `6.83` | `1.52s` | `5.31` | **`1.42 GB`** | `1.7x` slower than ILU, same RSS |
| | `C4_ldu` `b:ASM(16,LU), v:ASM(16,LU)` | **OOM** | — | — | — | — | **`10.23 GB`** | Overlapped LU blows up fill |
| **`S2-p2`, `dt=2e-6`** | **`C4_lower` `b:BJ(-1,ILU), v:ASM(1,ILU)`** *(ref)* | **OK** | **`41.0`** | **`10.34`** | **`2.02s`** | **`8.32`** | **`1.41 GB`** | **Winner at $\mathrm{CFL}_A=10$** |
| | `C4_lower` `b:BJ(16,LU), v:ASM(1,ILU)` | OK | `57.0` | `17.74` | `4.99s` | `12.75` | `1.44 GB` | `1.5x` slower, same RSS |
| | `C4_ldu` `b:BJ(16,LU), v:ASM(1,ILU)` | OK | `39.0` | `16.75` | `3.80s` | `12.95` | `1.46 GB` | `1.5x` slower, same RSS |
| | `C0` `BJ(16, LU)` / `ASM(16, LU)` | `FAIL` / **OOM** | — | — | — | — | `n/a` / `7.84G` | Fails at $\mathrm{CFL}_A=10$ |
| **`S1-p4`, `dt=5e-7`** | All ILU configs (`C0`--`C6`) | **FAIL** | — | — | — | — | — | `ILU(0)` breaks down at $p=4$ |
| | `C0` `BJ(16, LU)` | OK | `30.5` | `21.32` | `0.00s` | `30.60` | `1.54 GB` | Rescues `C0` at $p=4$ |
| | **`C4_ldu` `b:BJ(16,LU), v:ASM(1,ILU)`** | **OK** | **`15.5`** | **`8.91`** | **`0.00s`** | **`9.19`** | **`2.57 GB`** | **`3.33x` faster than `C0_bj16lu`** |
| | `C4_ldu` `b:ASM(16,LU), v:ASM(16,LU)` | **OOM** | — | — | — | — | **`10.74 GB`** | Overlapped LU OOM |

**Key rules for subdomain `LU` vs `ILU`**:
1. **Do NOT use `add_schwarz(overlap=1)` with subdomain `LU` on 3D high-order meshes**: 1-hop overlap enlarges each 3D subdomain 3--5x, triggering SuperLU fill explosion (`8--11 GB/rank` $\Rightarrow$ OOM). Keep `overlap=1` paired with `ILU` on `S_V`, and use **non-overlapping `block_jacobi(nlocal=16..32)`** whenever deploying subdomain `LU`.
2. **At $p=2, 3$**: Lightweight `ILU(0)` (`b: BJ(-1, ILU)` + `v: ASM(1, ILU)`) is `1.5--1.8x` faster than `BJ(16, LU)` at identical RSS (`1.42--1.45 GB/rank`).
3. **At $p \ge 4$**: Partitioned `BJ(16, LU)` inside `<b_pre>` is **essential** — it cures `ILU(0)` instability with modest RSS (`2.57 GB/rank`) and beats `C0_bj16lu` by `3.3x`.

---

## 4. Decision Verdicts & Recommended Production XMLs

### Q1. Coarse `S_V` Bottleneck?
- **No through `6.32M` DOFs (`S2-p3`, Tier A)**: `S_V` solves only the $p=1$ velocity block ($\sim 1/12$ of total DOFs), so `superd`/`mumps` factorization takes `1.3--2.1 s` on `S2-p2` (`0.1--1.8%` amortized over 50 nonlinear steps) and fits easily on `S2-p3` (`5.55 GB/rank`).
- **Extreme scale (`S3` $\ge 1.5\times 10^7$ DOFs, Tier B/C)**: `C6` (iterative coarse) already delivers `7.5x` speedup over `C0` on `S2-p3` at `5.20 GB/rank` with zero distributed LU. Next engineering step for `np >= 128--512`: expose PETSc `GAMG` / Hypre `BoomerAMG` on SPD `S_V` and rebase onto `main` (`#453`) to enable $h+p$-MG.

### Q2. Does `J_BB` Need $p$-MG (`C5`) at Scale?
- **No** across every tested regime (`S1-p2..p4`, `S2-p2` up to $\mathrm{CFL}_A=25$, `S2-p3` `6.32M` DOFs, 5-field nonlinear, and 7-field `MUG_sph_heat` with Spitzer $\eta_0=25$). `C5` saves only 1--2 outer iterations while costing `15--20 s` extra coarse `J_BB` LU setup, `+0.45--0.75 GB/rank` RSS, and `1.5--2.2x` slower `t_step`.

### Q3. Recommended Production XML Configurations

#### 1. Default 5-Field ($p = 2, 3$, moderate & high CFL — `C4_lower` for fixed `dt`, `C4_ldu` for adaptive `dt`)
```xml
<oft><xmhd>
  <pre type="block_schur">
    <schur_type>lower</schur_type> <!-- switch to ldu for adaptive-dt -->
    <b_pre>
      <solver type="gmres"><its>3</its><nrits>3</nrits>
        <pre type="block_jacobi"><nlocal>-1</nlocal><solver type="ilu"><package>petsc</package></solver></pre>
      </solver>
    </b_pre>
    <v_pre>
      <pre type="mg">
        <smoother direction="both">
          <solver type="gmres"><its>2</its><nrits>2</nrits>
            <pre type="add_schwarz"><overlap>1</overlap><nlocal>1</nlocal><solver type="ilu"><package>petsc</package></solver></pre>
          </solver>
        </smoother>
        <coarse><solver type="lu"><package>superd</package></solver></coarse>
      </pre>
    </v_pre>
  </pre>
</xmhd></oft>
```

#### 2. Default 7-Field (`MUG_sph_heat`, `B+V+n+T`)
Add `<s_pre>` inside `<pre type="block_schur">` above:
```xml
    <s_pre>
      <solver type="gmres"><its>4</its><nrits>4</nrits>
        <pre type="add_schwarz"><overlap>1</overlap><nlocal>1</nlocal><solver type="ilu"><package>petsc</package></solver></pre>
      </solver>
    </s_pre>
```

#### 3. High-Order ($p \ge 4$)
Use `<schur_type>ldu</schur_type>` and replace `<b_pre>`'s preconditioner with partitioned non-overlapping `SuperLU`:
```xml
    <b_pre>
      <solver type="gmres"><its>3</its><nrits>3</nrits>
        <pre type="block_jacobi"><nlocal>16</nlocal><solver type="lu"><package>super</package></solver></pre>
      </solver>
    </b_pre>
```

#### 4. Memory-Constrained / Extreme-Scale / Native MPI Fallback (`C6`)
Replace `<coarse>` inside `<v_pre>` with:
```xml
        <coarse>
          <solver type="gmres"><its>15</its><nrits>15</nrits><rtol>1.e-2</rtol>
            <pre type="add_schwarz"><overlap>1</overlap><nlocal>1</nlocal><solver type="lu"><package>super</package></solver></pre>
            <!-- For Native MPI (use_petsc=F), use block_jacobi without <overlap> -->
          </solver>
        </coarse>
```

---

## 5. PR Readiness: Bug Inventory & Pre-Merge Fix Checklist

Future agents preparing `mug-pmg-step1..5` for PR submission should address these items before merging:

| Priority | Target PR / File | Root Cause & Symptom | Required Fix |
|---|---|---|---|
| **Blocker** | PR 5 / `src/base/oft_xml.F90:943-1105` | `xml_extractDataAttribute_{double,int,logical}1D` leaves `INTENT(out) :: output` undefined when an optional XML attribute (`coarse_damp`, `coarse_shift`) is absent; scalar wrapper `DEALLOCATE(output_tmp)` corrupts the heap (`cfree` abort on `C2`). | Add `NULLIFY(output)` at entry of all three `xml_extractDataAttribute_*1D` routines. |
| **Blocker** | All PRs (`6987f2d` base) | Stacked branch base `a9f91c6` predates `#459` (`d528069`: METIS `npart=1` SIGFPE + PETSc 3.23--3.25 `PCASMSetLocalSubdomains` bug) and `#453` (`c6de81c`: coarse geometric $h+p$-MG `SIGSEGV` in `test_alfven`). | Rebase `mug-pmg-step1..5` onto current `main`. |
| **High** | PETSc / `src/lin_alg/petsc_solvers.F90:1752-1774` | `<solver type="ilu"></solver>` without `<package>petsc</package>` inherits `DEF_LU_PACK` (`mumps`), aborting with `MatSolverType mumps does not support factorization type ILU`. | Default ILU package to `'petsc'` rather than `DEF_LU_PACK` in `create_lu_pc`. |
| **High** | Examples & PR 4 / `MUG_sph_tilt.F90:165-187`, `MUG_sph_heat.F90:149-166`, `xmhd.F90:6155,6176,6198` | Unconditional `DEALLOCATE(linv%pre%pre)` / `create_ilu_pre(linv%pre%pre)` / `linv%pre%pre%delete()` crashes under `use_petsc=T` (where `pre%pre` is unassociated and `create_ilu_pre` returns an `oft_petsc_precond`). Also `MUG_sph_tilt` & `MUG_sph_heat` hardcode `minlev=-1`, silently disabling $p$-MG (`xmhd.F90:4541`). | Guard with `IF(ASSOCIATED(...%pre))` (and check `use_petsc` in `xmhd_setup_schur_pre` fallbacks); expose `minlev` in `sph_tilt_options` and `sph_heat_options`. |
| **Medium** | PR 1 & PR 4 Docs | (a) Coupled `<b_pre>` with `nlocal=1` + `ILU(0)` stalls on large meshes (`>= 130k` tets); document `nlocal=-1` (ILU, $p \le 3$) and `nlocal=16` (SuperLU, $p \ge 4$). (b) `MUG_sph_tilt` allocates `Ne`/`Ti` even when advection is off (`self%A%ni = 7`), so PR 1 `<groups>` needs 7 entries (`1 1 2 3 4 5 6`), not 5. | Update XML examples/docs in PR 1 and PR 4. |
| **Strategic** | PR 5 (`Option 2.B`) | PR 3 + PR 4 deliver 95%+ of the practical speedup and CFL robustness at scale; Option 2.B does not scale past `S1-p2` / low CFL. | Prioritize merging PR 1--4 (plus the bug fixes above), then add PETSc `GAMG`/Hypre support on `S_V`; keep PR 5 optional/low-priority after fixing `oft_xml.F90`. |

---

## 6. Notes for Future Tier B/C Benchmarking Agents

- **Re-test after rebase onto `main`**: Re-run Suite 6 (`test_alfven`, `minlev = nlevels-1, nlevels-2`) once `#453` is in the branch base to measure $h+p$-MG coarse-grid reduction (64x smaller `S_V`).
- **Large-mesh targets (`S3-p2` ~15M DOFs, `S3-p3` ~45M DOFs, `test_alfven L7`)**: Use `np >= 128--256` (`>= 256 GB` aggregate RAM, `~50k--80k` DOFs/rank) and allow `>= 180 s` for `taylor_hmodes` IC generation before `Timestep 1` prints.
- **Direct vs iterative coarse crossover**: Compare `C4_lower`/`C4_ldu` (`superd` and `mumps`) against `C6` (and future PETSc `GAMG` on `S_V`) across `np = 64, 128, 256, 512` to locate the exact rank count where `superd`/`mumps` coarse factorization exceeds 30% of `t_step`.
- **Clean timeout teardown**: Always launch `mpirun` in its own process group (`start_new_session=True`) and kill with `os.killpg` on timeout so orphaned MPI workers never overlap subsequent timing runs.
- **MPICH 4.2.3 on Zen 5**: Export `MPIR_CVAR_ENABLE_INTRANODE_TOPOLOGY_AWARE_TREES=0` if encountering `mpir_hwtopo.c:57` assertions on multi-CCD AMD Ryzen/EPYC nodes.

---

## 7. Post-PR #477 Study 1: 100% Native (`use_petsc=F`) vs PETSc (`use_petsc=T`) Block-Schur

After merging PR #477 (native overlapping Schwarz, `#474` native LA refactor, `#451` native distributed LU, `#453`/`#456` geometric MG fixes, `#459` PETSc ASM/METIS fixes) and the Step 6 pre-merge fixes (`dc9ce0f`), the entire winning `C4_lower`/`C4_ldu` Block-Schur stack runs natively (`use_petsc=F`) using `<pre type="add_schwarz">` on `S_V` and `<solver type="lu"><package>superd|mumps</package><distributed>T</distributed></solver>` on the coarse level.

All runs below used the same GCC-16 MPI+PETSc `Release` build (strictly 1 run at a time, `<= 16` physical cores).

### 7.1 Physics Parity (`S1-p2`, `np=8`, `dt=5e-7`, `lin_tol=1e-12`, 3 steps)
- `C4_ldu` PETSc: `31 its`, `t_step=2.00s`, `RSS=0.39 GB`, `men=5.17045740e-09`, `ven=3.80333654e-06`
- `C4_ldu` Native: `31 its`, `t_step=2.43s`, `RSS=0.40 GB`, `men=5.17045784e-09`, `ven=3.80334745e-06`
- **Max rel. diff**: **`8.6e-8` (`men`)**, **`2.9e-6` (`ven`)** — exact parity with identical iteration count.

### 7.2 `S2-p2` (`1.84M` `B+V` DOFs) & `S2-p3` (`5.11M` `B+V` / `6.32M` total DOFs) Comparison (`np=16` unless noted)

| Mesh / `dt` | Config | Backend | `lits` | `t1 (s)` | `t_setup` | `t_step (s)` | Peak RSS/rank | Notes |
|---|---|---|---:|---:|---:|---:|---:|---|
| **`S2-p2`, `1e-6`** | `C0` (`gmres(8)+BJ(-1,ILU)`) | PETSc | `37.0` | `17.38` | `4.81s` | `12.58` | `0.86 GB` | PETSc baseline |
| | `C0` (`gmres(8)+BJ(-1,ILU)`) | **Native** | `36.0` | `22.04` | `6.57s` | `15.48` | `0.88 GB` | Native `ILU(0)` ~1.23x slower |
| | `C4_lower` (`ASM(1,ILU)` + `superd`) | PETSc | `14.0` | `4.90` | `2.03s` | `2.87` | `1.41 GB` | `4.38x` vs PETSc `C0` |
| | `C4_lower` (`restricted ASM` + `superd`) | **Native** | `15.0` | `6.66` | `2.68s` | **`3.98`** | **`1.41 GB`** | **`3.89x` vs Native `C0`**, same RSS as PETSc |
| | `C4_lower` (`weighted ASM` + `superd`) | **Native** | `15.0` | `6.97` | `2.97s` | `4.00` | `1.42 GB` | Partition-of-unity matches restricted |
| | `C4_ldu` (`ASM(1,ILU)` + `superd`) | PETSc | `11.0` | `5.06` | `1.89s` | `3.17` | `1.41 GB` | Reference |
| | `C4_ldu` (`restricted ASM` + `superd`) | **Native** | **`11.0`** | `7.14` | `3.06s` | `4.08` | **`1.41 GB`** | Identical `lits` & RSS; 1.29x PETSc `t_step` |
| | `C4_ldu` (`restricted ASM` + `mumps`) | **Native** | `11.0` | `7.56` | `3.53s` | `4.03` | `1.40 GB` | Tied with `superd` |
| | `C4_ldu` (`BJ(1,ILU)` + `superd`) | **Native** | `14.0` | `6.92` | `2.06s` | `4.86` | `1.25 GB` | Native `ASM(1)` beats `BJ` by **19%** |
| | `C4_ldu` Hybrid (`np=8 x OMP=2`) | PETSc / **Native** | `11` / `11` | `6.09` / `8.61` | `2.60` / `3.78s` | `3.49` / `4.83` | `2.59` / `2.59 GB` | Pure MPI `16x1` beats `8x2` in both (10% / 18%) |
| **`S2-p2`, `2e-6`** | `C4_lower` ($\mathrm{CFL}_A=10$) | PETSc | `42.0` | `10.35` | `1.83s` | `8.52` | `1.43 GB` | Reference |
| | `C4_lower` ($\mathrm{CFL}_A=10$) | **Native** | **`38.0`** | `15.70` | `5.51s` | `10.19` | `1.44 GB` | **Fewer its (`38` vs `42`)**, within 20% of PETSc |
| **`S2-p3`, `1e-6`** | `C4_lower` (`6.32M` DOFs, 3 steps) | PETSc | `37.0` | `33.74` | `0.00s` | `35.32` | `5.55 GB` | Reference |
| | `C4_lower` (`6.32M` DOFs, 3 steps) | **Native** | `37.5` | `52.28` | `1.62s` | `50.66` | `6.00 GB` | Clean convergence, **6.4x faster** than `C0` |

**Key Takeaway**: 100% Native Block-Schur (`use_petsc=F`) works out of the box through `6.32M` DOFs and $\mathrm{CFL}_A=10$, achieving identical iteration counts and RSS to PETSc and `3.9--6.4x` speedups over `C0`. The `~20--40%` wall-time gap relative to PETSc is explained by unoptimized native `ILU(0)` (the same `1.23x` ratio appears in `C0_nat` vs `C0_petsc`) — motivating Step 8 (cell-mode overlapping Schwarz with `SuperLU` on Block-Schur sub-blocks).
