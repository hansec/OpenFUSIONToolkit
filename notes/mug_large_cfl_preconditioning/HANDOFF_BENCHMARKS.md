# Handoff: benchmarking the semi-implicit (`xmhd_si`) preconditioning approach at scale

> **Status (2026-11-27):** Tier 0 and the T1 (workstation) tier have been run on an Apple M2 Max (12 cores, 64 GB). Results are in [`RESULTS_M2Max.md`](RESULTS_M2Max.md) (raw tables in `results/m2max/`, campaign scripts in `scripts/m2max_campaign/`). §2 below is updated with what that campaign settled. Still open: the T2/T3 tiers (R2, p=3 at R1, strong scaling beyond 8 cores) and the follow-ups in RESULTS §10.
>
> **Script status:** `bench_case.py` and `analyze_runs.py` have been exercised on the full M2 Max campaign (577 case directories) and now handle BSD/macOS `time -l`. Caveats: `analyze_runs.py`'s automatic growth-rate window differs from a fixed 100–300 µs fit, and its `wall_per_simtime` for restarted runs includes MUG's dt ramp after a restart (RESULTS §2).

Branch `claude/pmg-large-cfl-preconditioning`. Background and all small-scale results: [`REPORT.md`](REPORT.md) (§7b is the starting point). Scripts: [`scripts/`](scripts/).

This document is a test plan for machines larger than the 4-core / ~14 GB VM used so far. It covers three questions:
- scaling with mesh resolution;
- scaling with finite-element order;
- the accuracy (physics) cost of large time steps, including runs long enough for a physics comparison.

---

## 1. What is being tested

`xmhd_si = s` adds the term **s·(dt/2)²·W·(Vⁿ⁺¹ − Vⁿ)** to MUG's momentum equation. W is the parabolized (linearized ideal-MHD) velocity operator, so the term is a NIMROD-style semi-implicit operator placed inside MUG's existing CN step (or BDF2/BE step for linear runs).
- It is O(dt²)-consistent.
- It keeps the step unconditionally stable when W is positive semidefinite (the default).
- It leaves all other terms unchanged.

Used together with the physics-based block preconditioner `<pb_pre type="1">`, it makes the velocity Schur approximation match the operator being solved. Offline analysis bounds eig(P_V⁻¹·S) ≥ s/(1+s), independent of dt.

Results so far (serial, small meshes; REPORT §7b):

| problem | result |
|---|---|
| Cube Alfvén wave, linear runs, PB type 1 + LU blocks | 8–14 outer its/step from CFL 2 to 190 with s = 1–2, against 13–39 with s = 0 |
| Cube, nonlinear CN | flat to CFL ~19 (11–12 its, against 38); at CFL ~190 limited by MUG's Newton (§2) |
| Spheromak tilt, base mesh (R0), p=2, PB type 1 + GMRES(8)/BJ-LU(1024) on P | dt=1e-6: 15 its (s=1) vs 91–104 (s=0). dt=3e-6: 31–32 its. **dt=1e-5: 115 its, no dt cuts** (s=0 is cut to 2.5e-6), about 254 s/step at 2.3 GB |
| Spheromak R0 on 2 MPI ranks | 16 its/step at dt=1e-6, s=1 (serial: 15); first MPI data point |
| Accuracy (cube, fixed final time) | invisible below ω·dt ≈ 0.1; time error about 2.5× CN's at s=1 (prediction: (1/12 + s·1.1/8)/(1/12)) |

## 2. Known vs. unknown

| item | status |
|---|---|
| `xmhd_si` with defaults off | regression suite unchanged (non-slow MUG: 166 passed, 56 skipped) |
| Residual and Jacobian consistency of the SI term | verified with a finite-difference check (agreement to about 1e-12) |
| MPI | **M2 Max:** cube and spheromak on 2–8 ranks; iteration counts within ±2 of serial; restarts are decomposition-independent. Coax railgun (rotationally periodic) on 3 ranks. Not tested: PETSc backend (not supported by `pb_pre`), > 8 ranks |
| Refined meshes (R ≥ 1) | **M2 Max:** cube R1 (A1/A1x/A2) and spheromak R1 (B): iterations depend on CFL, not h. R2 untested |
| FE order p=3, 4 | **M2 Max:** cube p=2–4 and spheromak p=3 run correctly with SI. With the default ~300-DOF inner blocks iterations grow with p and CFL. With an exact LU of P they stay at 5–14 |
| MG hierarchy on P (`pb_pre` with `nlevels>1` FE levels) | runs on the cube (Tier 0.5: 17–21 its), not benchmarked further |
| Spheromak accuracy (growth rate, saturation) | **M2 Max:** measured. Δγ/γ ≈ −1.5 %·s·(dt/1 µs)²; nonlinear tilt within 2 % at dt=1e-6 (s=1), 8–21 % low at 3e-6 (RESULTS §6) |
| Nonlinear Newton at CFL ≳ 100 | **known issue (pre-existing):** MUG's assembled Jacobian goes stale within a step (Lorentz-force block, parallel-velocity rows), and Newton stalls once SI allows big steps. Use `xmhd_mfnk=T` (the spheromak example already does). A Newton refresh-interval option does not exist yet (`nksolver%up_freq=4` is hard-coded in `xmhd_run`) |
| Inner solver on P | **known limit, confirmed on M2 Max:** iterations grow with dt (15 → 32 → 115 on the spheromak), and with p and CFL on the cube. An exact LU of P keeps outer iterations at 5–14 up to CFL 300, so this is the main remaining cost (RESULTS §5.2) |
| Production-style baselines | **M2 Max:** full-Jacobian block-Jacobi LU at 64–256 blocks is 1.2–2× faster than the example ILU and accurate at dt=1e-6 with s=0. PB + s=1 is only 1.06–1.18× faster than that, but uses about half the memory. SI + full-Jacobian LU looks fastest (from t=0) and is being confirmed (RESULTS §7) |

## 3. Build and setup

```bash
git fetch origin claude/pmg-large-cfl-preconditioning && git checkout claude/pmg-large-cfl-preconditioning
# Build per the oft-build skill / build_libs.py; enable examples so MUG_sph_tilt is built:
#   cmake ... -DOFT_BUILD_EXAMPLES=ON
make -j MUG_sph_tilt test_alfven      # in <build_dir>
export OFT_BUILD_DIR=<build_dir>      # used by the scripts (default: <repo>/builds/build_release)
```

- **Python:** numpy, h5py (and optionally the built `OpenFUSIONToolkit` package on `PYTHONPATH`, e.g. `<build_dir>/python`, for the history reader; `analyze_runs.py` has a fallback parser).
- **Threads:** use `OMP_NUM_THREADS=1–2` per MPI rank (`bench_case.py --threads`).
- **Containers running as root:** set `OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1`.

**Input snippets.** `bench_case.py` writes these for you.
```fortran
&xmhd_options
 ...
 xmhd_si=1.          ! semi-implicit scale s (0 = off). Optional: xmhd_si_sound, xmhd_si_full
 xmhd_mfnk=T         ! recommended with SI at Alfven CFL >~ 100
 ittarget=100000     ! benchmarking: keeps dt fixed unless Newton fails (MUG otherwise adapts dt)
/
```
Preconditioner presets (`--pre`):

| preset | XML |
|---|---|
| `ilu` (P0, the example's default) | `<pre type="gmres"><its>8</its><nrits>8</nrits><pre type="block_jacobi"><nlocal>-1</nlocal><solver type="ilu"/></pre></pre>` |
| `pb1` (P1, recommended) | `<pb_pre type="1"><solver type="gmres"><its>8</its><nrits>8</nrits><pre type="block_jacobi"><nlocal>N</nlocal><solver type="lu"/></pre></solver></pb_pre>` |
| `pb1_lu` | `pb_pre` type 1 with a direct LU of P (small cases) |
| `pb1_mg` (P2, cube only) | `<pb_pre type="1"><smoother direction="up">…GMRES/BJ…</smoother></pb_pre>` with `minlev` below the top FE level |
| `bj`, `lu`, `mg` | the same solvers applied to the full Jacobian (baselines) |

`nlocal` is the number of block-Jacobi blocks **per rank**. Keep about 250–300 DOFs per block, as at R0: 2.85e5 DOFs / 1024 blocks. Smaller `nlocal` means larger LU blocks and can blow up memory (`nlocal` ≤ 64 ran out of memory at R0).

**Mesh levels.** `nlevels` counts all levels including the MPI transfer level.
- Serial: `nbase=nlevels`.
- MPI: `nlevels = nbase + 1 + (refinements after decomposition)`.
- Refinements R = `nlevels−1` (serial) or `nlevels−2` (MPI), each giving 8× the cells. For example, on MPI, R0 = `nbase=1 nlevels=2`, R1 = `nbase=1 nlevels=3` (refine after decomposition; cheaper on memory), R2 = `nbase=1 nlevels=4`.

## 4. Tier 0: pre-flight smoke tests (minutes; run first on the target machine)

| # | test | command sketch | pass criterion |
|---|---|---|---|
| 0.1 | regression, SI off | `cd <build_dir>/tests && ./run_test.sh physics/MUG -m "not slow"` | all pass |
| 0.2 | cube, serial reference | `bench_case.py --geom cube --out t02 --dt 1.E-4 --nsteps 3 --si 1 --pre pb1_lu` | 11 its/step, 1 Newton it |
| 0.3 | cube MPI (R0 distributed) | same plus `--nbase 1 --nlevels 2 --ranks 2 --pre pb1 --nlocal 4` | no crash; its within ±2 of the serial `pb1` run with the same `nlocal` per rank |
| 0.4 | cube p=3 and p=4 | `--order 3` / `--order 4`, `--pre pb1_lu`, dt=1e-4 and 1e-3, `--linear` | its bounded in dt (compare s=0 vs s=1) |
| 0.5 | MG on P | cube `--order 3 --minlev 1 --pre pb1_mg`, s=1, `--linear`, dt=1e-4 | runs; its reported (no reference yet) |
| 0.6 | spheromak serial reference | `--geom sph --dt 1.e-6 --nsteps 2 --si 1 --pre pb1 --nlocal 1024` | 15 its/step |
| 0.7 | spheromak MPI | same at `--nbase 1 --nlevels 2 --ranks 2 --nlocal 512` (done on the VM: 16 its) and at 8–32 ranks | its within a few of serial |

Report every failure: these combinations are new.

## 5. Benchmark A: cube Alfvén wave (analytic accuracy and clean scaling)

A traveling Alfvén wave in a triply periodic 1×1×2 box (`test_alfven`); B0 = ẑ, v_A = 1e4 m/s. The lowest wave has ω = π·v_A = 3.14e4 s⁻¹, and the fastest discrete wave at R0 p=2 is about 1.9e5 s⁻¹.

**Matrix:**
- resolution R0 → R3 (R0 = the regression mesh; 192 cells, B 1856 / V 768 DOFs at p=2);
- order p = 2, 3, 4 (DOFs relative to p=2: B ×2.6 / ×5.3, V ×3.4 / ×8);
- dt from a CFL ladder, max Alfvén CFL ∈ {0.3, 3, 30, 300}. The fastest wave grows about ∝ p/h, so choose dt per (R, p) to hit the CFL values: measure ω_max once per mesh with `scripts/si_analysis.py`, or scale from 1.9e5·(p/2)·2^R;
- s ∈ {0, 0.5, 1, 2};
- linear (BDF2) and nonlinear (CN, both `--mfnk F` and `T`).

**Runs:**
- A1, scaling: 6–10 steps per case.
- A2, accuracy: fixed final time T = 1e-4 (the regression setup: `nsteps = T/dt`), reading `alfven.results` (relative L2 errors of B and V).

**Metrics:**
- outer its/step, Newton its, dt reductions, time/step, peak RSS;
- B/V waveform error.

Split the error into a spatial floor (the smallest dt) and a time error (the difference). Compare the time error to the prediction ≈ (1/12 + s·s_eff/8)(ω·dt)², with s_eff ≈ 1.1 measured at R0 p=2.

**What "good" looks like:**
- its/step bounded in CFL for s ≥ 1 at every (R, p);
- mild growth with R and p, at most about √2 per refinement;
- time error ratio vs s=0 close to the prediction.

**Optional A3:** a nonlinear long run, 20–50 wave periods at CFL 3 and 30 with s = 0 / 0.5 / 1, to check phase and amplitude drift.

## 6. Benchmark B: spheromak tilt scaling (short runs)

The example `src/examples/MUG/Spheromak_tilt`: an ideal tilt instability in a tall cylindrical flux conserver. No density/temperature advance, η=1, ν=10, matrix-free Newton.

**Matrix:**
- resolution R0, R1, R2 (`nbase=1`, `nlevels` = 2, 3, 4 on MPI). For different geometric resolutions, finer meshes can also be built from `cyl_tilt.geo` with gmsh and `OFT_convert_gmsh.py`;
- order p = 2, 3;
- dt ∈ {1e-6, 3e-6, 1e-5, 3e-5};
- s ∈ {0, 1, 2};
- preconditioner `pb1` (BJ-LU inner, `--inner_its` 8 and 20) against `ilu` (the example default) and `bj` on J.

Run 5–10 steps per case (steps 1–2 include setup).

**Record:**
- outer its/step, Newton its/step, dt reductions;
- wall time per step and per µs of simulated time (`analyze_runs.py` column `wall_per_simtime_s_per_us`);
- peak RSS per node (`time.log`);
- setup/factorization cost vs solve cost. Compare the first step with later steps, and steps at which the preconditioner is rebuilt (`xmhd_prefreq`).

**Scaling:**
- **strong:** the best configuration at R1 on 8 → 64 ranks;
- **weak:** R0/2 ranks → R1/16 → R2/128, keeping DOFs per rank and per block fixed.

**What "good" looks like:**
- the SI advantage (no dt cuts at 1e-5; its far below s=0) survives refinement;
- its grow at most mildly with R and p at fixed dt.

Expect higher CFL at fixed dt on finer meshes. Present results both at fixed dt and at fixed CFL.

## 7. Benchmark C: accuracy and physics (long runs; the main open question)

**C1. Linear tilt growth rate.**
- Setup: `--linear` (BDF2), with s ∈ {0, 0.5, 1, 2} and dt ∈ {2e-7, 1e-6, 3e-6, 1e-5}, at R0 and R1, for p=2 and p=3.
- Run length: until the kinetic energy has grown by at least 10⁴ in a clean exponential phase. Check the length on the reference first.
- Measure: γ = ½·d ln(KE)/dt (`analyze_runs.py`, auto window or `--fit t0,t1`; check `gamma_r2` > 0.999).
- Reference: s=0 at the smallest dt, with a dt-convergence check (2e-7 vs 1e-7).
- Acceptance: |γ − γ_ref|/γ_ref ≤ 1–2% at the dt of interest. Report, for each s, the largest dt meeting 1% and 5%.
- **Note:** the tilt is an *Alfvénic* instability, so γ·dt (not the grid CFL) controls the error. An SI benefit is expected where γ·dt ≪ 1 while the grid CFL ≫ 1.

**C2. Nonlinear tilt to saturation (physics comparison).**
- Setup: the full example, 0.6 ms of simulated time.
  - Reference: s=0, dt=2e-7 (the example's own setting, 3000 steps; use the `ilu` or `pb1` preconditioner).
  - Test runs: s=1 at dt = 1e-6 (600 steps), 3e-6 (200), 1e-5 (60). Optional: s=0.5 at 3e-6.
  - Resolution: R0 first, then R1 if affordable.
  - Use `--rst_freq 50` so long jobs can restart.
- Compare (`analyze_runs.py --traces traces.csv`; plot `ven` and `tcurr` against time for all cases on one figure):
  - time and value of the kinetic-energy peak (`t_ke_peak`, `ke_peak`);
  - time of the toroidal-current collapse (`t_curr_half`: current below 50% of its initial value);
  - final/initial magnetic energy (`men_final_over_initial`);
  - visual agreement of the linear-growth slope.
- Acceptance (suggested): peak and collapse times within 3–5%, peak amplitude within ~5–10%, saturated magnetic energy within a few %. Report deviations as they are; the trade-off is the point of the study.

**C3 (optional).** Repeat C2 with `xmhd_si_full=T` (adds the J0 term; the spheromak is force-free, so J0 ∥ B0). This shows whether the positive-semidefinite simplification costs accuracy or iterations.

## 8. Run matrix by tier

DOF and memory figures are **estimates** extrapolated from measured R0 values: R0 p=2 spheromak has 2.85e5 DOFs (B 1.66e5, V 7.1e4, n/T 2.4e4 each) on 16,305 cells. DOFs scale by 8× per refinement and about 2.6–3.4× for p=3. Measured R0 cost (s=1, `pb1`, 2 threads, 2.3 GB): 50 s/step at dt=1e-6, 77 s at 3e-6, 254 s at 1e-5.

| tier | machine | cases | sizes (est.) | est. cost |
|---|---|---|---|---|
| T1 workstation | 8–16 cores, 32–64 GB | Tier 0; A at R0–R1, p=2–4; B at R0 p=2/3; C1 at R0 p=2; C2 at R0 for s=1 at 1e-5 and 3e-6 | spheromak R0: 2.9e5 DOFs (p=2), ~9e5 (p=3, ~8–12 GB) | C2 s=1 at 1e-5: ~4 h; at 3e-6: ~4.5 h; C1 sweep: ~1 day |
| T2 node | 32–128 cores, 256–512 GB | A at R2–R3; B at R1 p=2/3 + strong scaling; C1 at R0/R1 p=2/3; C2 reference (s=0, dt=2e-7) at R0; C2 s=1 at R1 | spheromak R1: 2.3e6 DOFs (p=2, ~20–30 GB), 7e6 (p=3) | C2 reference at R0: 3000 steps, likely 1–2 days on a node; size it from the first 20 steps |
| T3 multi-node | 100–1000 ranks | B at R2 + weak scaling; C2 at R1–R2 (reference and s=1) | spheromak R2: 1.8e7 DOFs (p=2) | dominated by the C2 reference; may need restarts |

Do T1 first and only scale up once the T0/T1 results look sane.

## 9. Data to collect and reporting

- **One directory per case** (`bench_case.py --out <name>`). Suggested name: `<geom>_R<r>_p<p>_<pre>_s<s>_dt<dt>[_lin]_np<ranks>`. Keep `log.txt`, `xmhd.hist`, `case.json`, `time.log`, `alfven.results`. Restart and plot files are optional.
- **Summary:** `python3 analyze_runs.py <dirs> --csv summary.csv [--traces traces.csv]`.
  - Columns: case parameters, `steps_done`, `lits_mean/max`, `nlits_mean`, `dt_min`, `dt_reduced_steps`, `wall_per_step_s`, `wall_per_simtime_s_per_us`, `peak_rss_GB`, `berr/verr`, `gamma`, `gamma_r2`, `ke_peak`, `t_ke_peak`, `t_curr_half`, `men_final_over_initial`, `failed`.
  - `peak_rss_GB` comes from GNU time on the launcher. For `mpirun` it reflects the largest single process it waited on, not the node total; use the scheduler's accounting for node memory.
- **Plots:**
  - its/step against CFL, one curve per (R, p, s);
  - wall per simulated µs against dt;
  - γ error against dt per s;
  - `ven`/`tcurr` traces for C2.
- **Write-up:** add a section to `REPORT.md` (or a new `RESULTS_<machine>.md`) using the same table layout as REPORT §7b, and commit the `summary.csv`/`traces.csv` files (not the raw run directories).

## 10. Decision criteria

| observation | implication |
|---|---|
| its/step bounded in dt but growing ≳ 2× per refinement or per order | the inner solver on P is the bottleneck. Prioritize a better V-block solve and an auxiliary-space (Hiptmair–Xu) solver for the H(curl) B block (REPORT §8 item 1) |
| γ error at s=1 > 2% where dt is otherwise attractive | try s=0.5 (the bound s/(1+s) = 1/3 still holds), or accept a smaller dt; quantify the trade-off curve |
| C2 saturation physics changes at large dt even with good γ | SI slows nonlinear/fast dynamics too much; restrict to linear or slowly evolving phases, or adapt s |
| Newton failures / dt cuts at CFL ≳ 100 with `xmhd_mfnk=T` | implement a Newton Jacobian-refresh option, or a refresh at the iterate with the SI block kept at the start-of-step reference (already the case in the code) |
| MPI its much higher than serial at the same blocks per rank | rank-interface effect (seen in sibling work with block-Jacobi); test `boverlap`/`loverlap` options on the inner solver |
| s ≥ 1 worse than s=0 anywhere | report; unexpected (possible J0/pressure effects in non-uniform equilibria; try `xmhd_si_full=T`) |

## 11. Pitfalls

- **Silent dt changes:** MUG halves dt on Newton failure and adapts dt toward `ittarget`. Use a large `ittarget` (the default in `bench_case.py`) and always check `dt_reduced_steps`. The cube's first steps at dt=8e-7 have been seen to cut dt in some settings.
- **`nlevels` vs MPI:** serial needs `nbase=nlevels`; MPI needs `nlevels>nbase`. `nlevels=2, nbase=1` on MPI is the *same* mesh as serial R0.
- **xMHD multigrid is polynomial only;** h-levels cannot be MG levels. The spheromak example's preconditioner is single-level.
- **`pb_pre` type 4** (B-Schur) is serial only and expensive. It is not part of this study.
- **LU memory:** keep `nlocal` per rank so blocks stay at about 300 DOFs.
- **The `dump` attribute** (matrix dumps for offline analysis) only triggers at the first preconditioner update, so it needs `nsteps ≥ 2`.
- **Conservation diagnostics:** the SI term changes the effective inertia, so kinetic/magnetic energy budgets at large dt include its O(dt²) contribution.
- **Restarting with a different `xmhd_si` or dt mid-run:** allowed, but treat the restart point as a new initial condition.
