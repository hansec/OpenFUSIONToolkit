# Benchmark results: `xmhd_si` + `pb_pre` on an Apple M2 Max workstation (Tier 0 / T1)

Branch `claude/pmg-large-cfl-preconditioning` @ `5e3b5c7` (code unchanged by this campaign), following
[`HANDOFF_BENCHMARKS.md`](HANDOFF_BENCHMARKS.md). Also covers the coax railgun example (rotational periodicity) from
`claude/mug-railgun-example-437435`. Background: [`REPORT.md`](REPORT.md) §7b.

**Status (2026-11-27):** Tier 0, A1, A1x, A2, C1, C2, D and B at R0 are **complete**. Still running: B at R1 (weak and
strong scaling), the railgun timing set, and a follow-up timing set T2 (§9). Numbers marked *pending* will be filled in
when those finish.

Raw tables: [`results/m2max/`](results/m2max/) (`summary.csv` from `analyze_runs.py` for all 577 finished case
directories, `traces_C1.csv`/`traces_C2.csv`, and one text table per study). Campaign scripts:
[`scripts/m2max_campaign/`](scripts/m2max_campaign/).

---

## 1. Headline findings

| # | Finding | Evidence |
|---|---|---|
| 1 | **The SI term and `pb_pre` type 1 are correct in every setting tried**: MPI (2–8 ranks), p = 2–4, R0/R1, nonlinear, and the railgun's rotationally periodic velocity packing. | Tier 0 all pass (§3). Railgun: PB with s=0 reproduces the example to ~1e-5, with s=1 to ~1e-3, with identical sheet speed (§8). |
| 2 | **The SI error is predictable: Δγ/γ ≈ −1.5 %·s·(dt/1 µs)² on the tilt**, linear in s and quadratic in dt (saturating above ~3 µs). On the cube the SI/CN time-error ratios match (1/12 + 1.1·s/8)/(1/12) within ~10 %, and the error depends on ω·dt of the physical wave, not on CFL, mesh or order. | C1 (§6.1), A2 (§5.3) |
| 3 | **Accuracy limits s=1 to dt ≲ 0.8 µs (1 % in γ) or ≲ 1.8 µs (5 %) on the spheromak.** At dt = 3e-6 the nonlinear KE peak is 8 % (s=0.5) to 21 % (s=1) low; `xmhd_si_full` fixes the peak *time* but not the amplitude. dt = 1e-5 is qualitatively wrong (−50 % peak, +60 % peak time). | C1, C2 (§6) |
| 4 | **P is an excellent preconditioner; the inner solve on P is the bottleneck.** With an exact LU of P, outer iterations stay at 5–14 for every p (2–4), mesh (R0/R1) and CFL up to 300. With the default ~300-DOF block-Jacobi inner solve they reach 73–312 or fail. | A1x (§5.2) |
| 5 | **Outer iterations depend on CFL, not on h** (R0 and R1 agree at fixed CFL). At fixed dt, each refinement doubles the CFL. | A1 (§5.1) |
| 6 | **Replacing the example's ILU with block-Jacobi LU is the simplest robust speedup**: 1.2–1.36× at the same dt, and it allows dt = 1e-6 at s=0 with *no* measurable accuracy loss (KE peak −0.1 %, γ −0.002 %), for 1.24–2.0× over ILU at dt = 2e-7. Memory 1.7 → 3.6–6 GB. | D (§7.1), C2 |
| 7 | **PB + s=1 at dt = 1e-6 is 1.43–2.40× faster than ILU, but only 1.06–1.18× faster than block-Jacobi LU at s=0, dt = 1e-6**, at a cost of −1.5 % γ / −1.8 % KE peak. Its clear advantage is memory (~3 GB vs 6–21 GB for full-Jacobian LU of similar speed). | D (§7.1) |
| 8 | **SI combined with *full-Jacobian* block-Jacobi LU is the fastest configuration seen so far** (from t=0: 1316 µs/h at dt = 1e-6, s=1, vs 324 at s=0 and 781 for PB s=1). This disagrees in magnitude with the D restart segments (parity with PB at 16 blocks). **Being confirmed by T2.** | B (§7.2) |
| 9 | Large dt (≥ 3e-6) does not pay off on the tilt: iterations triple, Newton iterations double, throughput is no better than at 1e-6, and accuracy is clearly degraded. | D, C2 |
| 10 | Railgun: the −4 % sheet-speed loss at 2–4× dt comes from the time step itself (s=0 and s=1 agree, and the example's own solver gives the same at 6e-8). **The shipped dt = 3e-8 is itself ~1–2 % from converged.** | §8 |

---

## 2. Machine, build and protocol

- **Machine:** Apple M2 Max, 12 cores (8 performance + 4 efficiency), 64 GB, macOS. Concurrent jobs capped at ≤ 54 GB total RSS
  (watchdog in `runq.py`); no swap was used at any point.
- **Build:** GCC 16 (gfortran) + MPICH 4.2.3 + PETSc 3.23 (UMFPACK/MUMPS/SuperLU), OFT native solvers at runtime
  (`use_petsc=F`). Needed the PR #479 CMake fix for SuiteSparse include paths (local only, not committed here).
- **Two lanes:**
  - *Accuracy lane* (A1, A1x, A2, C1, C2, railgun correctness): several jobs at once on shared cores. Used **only** for
    iteration counts and accuracy. Wall times from this lane are indicative only.
  - *Timing lane* (D, B, railgun timing, T2): **one case at a time on an otherwise idle machine, always 4 MPI ranks ×
    2 OpenMP threads** (8 performance cores), so mixed core types and contention do not distort timings.
- **Metric:** simulated µs per wall-clock hour (µs/h) from per-step wall times. D segments restart from the C2
  reference; MUG starts a restart at the file's dt (2e-7) and ramps up over ~7 steps, so only full-dt steps are counted.
  These include one Jacobian refactorization per segment, which slightly overstates it relative to once per
  `xmhd_prefreq` = 20 steps.
- **Script fixes made for macOS** (committed): `bench_case.py` detects GNU `time -v` vs BSD `time -l` and sets
  `GFORTRAN_UNBUFFERED_PRECONNECTED=y` (buffered output looked like stalled runs); `analyze_runs.py` parses the BSD
  peak-RSS format (bytes). `analyze_runs.py` output was spot-checked against the hand-computed tables below (iterations, Newton counts, peak
  RSS, C2 metrics). Its automatic growth-rate window differs from the fixed 100–300 µs window used in §6.1, so use
  `fit_C1.py` (or `--fit 1e-4,3e-4`) for γ. Its `wall_per_simtime` for D includes the restart dt ramp; use
  `tab_D.py` for D throughput.

---

## 3. Tier 0 (all pass)

| # | Test | Result |
|---|---|---|
| 0.1 | Non-slow MUG regression suite | 102 passed, 0 failed, 120 skipped (skips are by design in a PETSc build) |
| 0.2 | Cube serial, s=1, `pb1_lu` | 11 its, 1 Newton (matches reference) |
| 0.3 | Cube on 2 ranks, `pb1` | 13–14 its vs 13 serial |
| 0.4 | Cube p = 3/4, exact LU of P | s=1: 9–13 its; s=0: 46–70 |
| 0.5 | p-MG on P (`pb1_mg`) | runs, 17–21 its |
| 0.6 | Spheromak serial, s=1, dt=1e-6 | 15 its (matches reference) |
| 0.7 | Spheromak on 2/4/8 ranks | 16 / 15 / 16 its; ~13–16 s/step serial (about 3× faster than the VM) |

---

## 4. Test problems

| Problem | Size | Reference / metric |
|---|---|---|
| Cube Alfvén wave (`test_alfven`) | R0: 192 cells, ~3e3 DOFs at p=2 (×8 per refinement; ×2.6–8 for p = 3–4) | Analytic wave; `alfven.results` B error at T = 1e-4 |
| Spheromak tilt (`MUG_sph_tilt`), R0 p=2 | 16,305 cells, 2.85e5 DOFs | Linear γ vs s=0, dt=1e-7 (γ_ref = 2.5953e4 s⁻¹, fit 100–300 µs); nonlinear to 600 µs vs ILU s=0, dt=2e-7 (KE peak 24.42 at 217 µs, current half-time 224 µs, men_f/i 0.9488) |
| Coax railgun (`MUG-railgun`, rotational periodicity) | 23,648 cells, ~4.0e5 DOFs (B 2.34e5, V 9.9e4, N/T 3.3e4) | Shipped example at dt=3e-8 on 3 ranks: sheet speed 88.3 / 89.0 / 88.4 km/s; relative L2 field differences |

---

## 5. Benchmark A: cube

### 5.1 A1 — outer iterations per step (`pb1`, ~300-DOF block-Jacobi LU inner, 6 steps)
Full tables: `results/m2max/A1_iterations.txt` (R0 and R1; linear, CN, CN + matrix-free Newton).

Linear (BDF2), s = 0 / 0.5 / 1 / 2:

| p | CFL | R0 | R1 |
|---|---|---|---|
| 2 | 3 | 27 / 12 / 10 / 9 | 23 / 12 / 10 / 9 |
| 2 | 30 | F / 20 / 20 / 22 | 215 / 22 / 23 / 26 |
| 2 | 300 | F / 79 / 73 / 62 | F / 150 / 144 / 141 |
| 3 | 30 | F / 31 / 34 / 44 | F / 28 / 29 / 33 |
| 3 | 300 | F / 301 / 312 / 345 | F / 265 / F / F |
| 4 | 30 | F / 49 / 55 / 60 | F / 53 / 53 / 54 |
| 4 | 300 | F / F / F / F | F / F / F / F |

- s > 0 removes the s=0 breakdown at CFL ≥ 30, and iteration counts are nearly the same at R0 and R1 for fixed CFL.
- They are **not** CFL-bounded with this inner solver: they grow with CFL and p. §5.2 shows why.
- Nonlinear CN at CFL 300 needs 2.5–5 Newton iterations and dt cuts; matrix-free Newton roughly halves the Newton count.

### 5.2 A1x — inner solver on P (linear, s=1; outer its/step)
`results/m2max/A1x_inner_solver.txt`

| case | default (~300-DOF blocks) | exact LU of P | ~3000-DOF blocks | 20 inner its |
|---|---|---|---|---|
| R0 p2 CFL 300 | 73 | **14** | 14 | 27 |
| R0 p3 CFL 300 | 312 | **13** | 23 | 85 |
| R0 p4 CFL 30 | 55 | **9** | 31 | 19 |
| R0 p4 CFL 300 | F | **11** | 191 | 141 |
| R1 p2 CFL 300 | 144 | **13** | 53 | 52 |
| R1 p3 CFL 300 | F | **11** | 183 | 101 |
| R1 p4 CFL 30 | 53 | UMFPACK OOM | 36 | 20 |

With an exact solve of P, the outer iteration count is essentially independent of dt, p and h. The three R1 p4
exact-LU cases failed with "UMFPACK: Out of Memory in factorization" at 5.7 GB RSS: this is the 2 GiB factor cap of
the 32-bit `umfpack_di` interface, not machine memory. MUMPS would get past it.

### 5.3 A2 — accuracy at fixed final time T = 1e-4 (CN; relative B error)
`results/m2max/A2_accuracy.txt`. Example rows (the spatial floor is negligible for R0 p4, R1 p3, R1 p4):

| N steps (dt) | s=0 | s=0.5 | s=1 | s=2 |
|---|---|---|---|---|
| 125 (8e-7) | 2.1e-4 | 3.3e-4 | 4.5e-4 | 6.9e-4 |
| 25 (4e-6) | 4.13e-3 | 7.22e-3 | 1.03e-2 | 1.65e-2 |
| 10 (1e-5) | 2.55e-2 | 4.45e-2 | 6.34e-2 | 0.101 |
| 5 (2e-5) | 9.76e-2 | 0.170 | 0.240 | 0.376 |

- The time error at fixed dt is **identical across R0 p4, R1 p3 and R1 p4** (CFL ranging 1.5–15): it depends on ω·dt of
  the physical wave, not the grid CFL.
- Ratios to s=0 at N=25: 1.75 / 2.49 / 3.98 against the predicted 1.83 / 2.65 / 4.30.
- Oddity: at p=4 the N=500 error (8.0e-4 at R0, 8.7e-4 at R1) is *larger* than at N=125. Most likely the
  per-step Newton tolerance (`nl_tol` = 1e-5) accumulating over 500 steps; not a time-step effect.

---

## 6. Benchmark C: spheromak accuracy

### 6.1 C1 — linear growth rate (BDF2, PB type 1 with 1024 blocks for all s; fit 100–300 µs)
`results/m2max/C1_growth_rate.txt`, `traces_C1.csv`. Error relative to γ_ref = 2.5953e4 s⁻¹ (s=0, dt=1e-7):

| dt | s=0 | s=0.5 | s=1 | s=2 | outer its (s=0 / s=1) |
|---|---|---|---|---|---|
| 2e-7 | +0.001 % | −0.03 % | −0.07 % | −0.15 % | 11 / 6 |
| 1e-6 | −0.002 % | −0.76 % | −1.51 % | −2.94 % | 46 / 15 |
| 3e-6 | +0.17 % | −6.0 % | −10.8 % | −18.1 % | 147 / 39 |
| 1e-5 | fails (GMRES max its) | −33 % | −45 % | −57 % | — / 102 |

- **Model:** Δγ/γ ≈ −1.5 %·s·(dt/1 µs)² (linear in s, quadratic in dt; saturates above ~3 µs).
- **The time discretization itself costs nothing up to 3e-6 at s=0**: essentially all of the error is the SI term.
- **Largest dt for a given γ error:** 1 % → ~1.15 µs (s=0.5), ~0.8 µs (s=1), ~0.6 µs (s=2); 5 % → ~2.6 µs (s=0.5),
  ~1.8 µs (s=1).

### 6.2 C2 — nonlinear tilt to 600 µs (serial, 2 threads)
`results/m2max/C2_nonlinear.txt`, `traces_C2.csv`

| case | KE peak | peak time | current half-time | men_f/i | outer its | Newton | peak RSS |
|---|---|---|---|---|---|---|---|
| reference: ILU, s=0, dt 2e-7 | 24.42 | 217 µs | 224 µs | 0.9488 | 12 | 1.00 | 1.4 GB |
| **block-Jacobi LU (16 blocks), s=0, dt 1e-6** | 24.39 (−0.1 %) | 217 µs | 224 µs | 0.9488 | 24 | 1.55 | 7.3 GB |
| PB, s=1, dt 1e-6 | 23.99 (−1.8 %) | 219 µs (+0.9 %) | 226 µs (+1.0 %) | 0.9487 | 23 | 1.68 | 2.8 GB |
| PB, s=0.5, dt 3e-6 | 22.38 (−8 %) | 225 µs (+3.7 %) | 234 µs (+4.6 %) | 0.9486 | 94 | 2.53 | 2.6 GB |
| PB, s=1, dt 3e-6 | 19.16 (−21 %) | 234 µs (+7.8 %) | 243 µs (+8.6 %) | 0.9485 | 88 | 2.41 | 2.6 GB |
| PB, s=1 + `xmhd_si_full`, dt 3e-6 | 19.24 (−21 %) | 219 µs (+0.9 %) | 231 µs (+3.1 %) | 0.9484 | 92 | 2.48 | 2.6 GB |
| PB, s=1, dt 1e-5 | 12.13 (−50 %) | 347 µs (+60 %) | 355 µs (+58 %) | 0.9524 | 310 | 3.10 | 2.6 GB |

- Against the suggested acceptance (times within 3–5 %, amplitude within 5–10 %), only dt = 1e-6 (s=1) passes; s=0.5
  at 3e-6 is borderline.
- The final magnetic energy is insensitive to dt and s, apart from the dt = 1e-5 run.
- `xmhd_si_full` (adding the J0 term) corrects the timing but not the amplitude. Iterations are unchanged.

---

## 7. Timing (exclusive, 4 ranks × 2 threads)

### 7.1 D — production-style preconditioners over the tilt (20-step restart segments from the C2 reference)
`results/m2max/D_timing.txt`. Throughput in simulated µs per wall hour (× = relative to ILU at dt 2e-7 in the same
phase); phases restart at t = 100 µs (linear growth), 220 µs (KE peak), 400 µs (decay).

| config | dt | linear | KE peak | late | outer its (lin/peak/late) | peak RSS | accuracy (C1/C2) |
|---|---|---|---|---|---|---|---|
| ILU (example), s=0 | 2e-7 | 340 (1.00) | 302 (1.00) | 509 (1.00) | 16 / 18 / 10 | 1.7 GB | reference |
| BJ LU 256 blocks, s=0 | 2e-7 | 464 (1.36) | 411 (1.36) | 614 (1.21) | 9 / 10 / 6 | 3.6 GB | reference |
| BJ LU 64, s=0 | 2e-7 | 364 (1.07) | 348 (1.15) | 530 (1.04) | 8 / 9 / 5 | 5.5–5.9 GB | reference |
| BJ LU 16, s=0 | 2e-7 | 215 (0.63) | 202 (0.67) | 311 (0.61) | 8 / 8 / 5 | 9.3 GB | reference |
| BJ LU 4 (1/rank), s=0 | 2e-7 | 87 (0.25) | 80 (0.27) | 135 (0.26) | 7 / 8 / 4 | 16–18 GB | reference |
| ASM LU 16 (overlap), s=0 | 2e-7 | 263 (0.77) | 206 (0.68) | 363 (0.71) | 4 / 5 / 2 | 12.5–13.6 GB | reference |
| **BJ LU 64, s=0** | 1e-6 | **693 (2.04)** | **374 (1.24)** | **852 (1.68)** | 23 / 41 / 19 | 5.7–6.4 GB | ~exact |
| BJ LU 16, s=0 | 1e-6 | 458 (1.35) | 230 (0.76) | 524 (1.03) | 20 / 39 / 18 | 13–20 GB | ~exact |
| ASM LU 16, s=0 | 1e-6 | 317 (0.93) | 192 (0.64) | 280 (0.55) | 20 / 30 / 22 | 20–26 GB | ~exact |
| PB 1024, s=1 | 1e-6 | 810 (2.38) | 404 (1.34) | 961 (1.89) | 16 / 30 / 13 | 2.8–3.0 GB | −1.5 % γ, −1.8 % KE |
| **PB 256, s=1** | 1e-6 | **816 (2.40)** | **431 (1.43)** | **902 (1.77)** | 14 / 25 / 12 | 3.0–3.3 GB | −1.5 % γ, −1.8 % KE |
| PB 64, s=1 | 1e-6 | 705 (2.07) | 332 (1.10) | 756 (1.49) | 12 / 24 / 11 | 4.6–5.0 GB | same |
| PB 16, s=1 | 1e-6 | 415 (1.22) | 193 (0.64) | 455 (0.89) | 11 / 23 / 10 | 8–10 GB | same |
| BJ LU 16, s=1 | 1e-6 | 890 (2.61) | 375 (1.24) | 930 (1.83) | 9 / 20 / 9 | 17–21 GB | same |
| PB 256, s=0.5 | 3e-6 | 648 (1.91) | 323 (1.07) | 749 (1.47) | 52 / 106 / 45 | 2.9–3.2 GB | −6 % γ, −8 % KE |
| BJ LU 16, s=0.5 | 3e-6 | 1151 (3.38) | 425 (1.41) | 1158 (2.28) | 19 / 62 / 20 | 16–21 GB | same |
| PB 64, s=1 | 3e-6 | 720 (2.12) | 294 (0.97) | 633 (1.24) | 34 / 85 / 39 | 4.8–4.9 GB | −11 % γ, −21 % KE |

Observations:
- **Block size matters more than the method.** Full-Jacobian block-Jacobi LU is best at 64–256 blocks in total
  (~1–4.5k DOFs per block). Fewer, larger blocks cut iterations a little but cost far more in factorization time and
  memory. Additive Schwarz halves iterations but loses in wall time and memory.
- **The KE peak is the expensive phase for every method.** Newton iterations double (1 → 2–3), and the large-dt cases
  suffer most.
- **PB with s=1 vs a well-chosen conventional setup at the same dt** (BJ LU 64, s=0, dt 1e-6): 1.18× / 1.15× / 1.06×
  faster, with about half the memory and a 1.5–1.8 % accuracy cost.
- **Smaller PB inner blocks are faster.** 256–1024 blocks beat 64 and 16: the LU cost grows faster than the
  outer-iteration savings.
- **SI also helps full-Jacobian LU.** BJ LU 16 at dt 1e-6 goes from 20 to 9 its/step with s=1. The missing
  BJ LU 64/256 + s>0 combinations are queued as T2.

### 7.2 B — spheromak from t=0 (8 steps; includes setup and the first factorization)
`results/m2max/B_timing.txt`. R0 p=2, µs/h (outer its/step):

| config | dt 1e-6 | dt 3e-6 | dt 1e-5 | dt 3e-5 | peak RSS |
|---|---|---|---|---|---|
| ILU, s=0 | 230 (122) | 176 (494) | 135 (406, dt cut) | fails | 1.7 GB |
| BJ LU 64, s=0 | 324 (48) | 385 (127) | 586 (265) | 281 (468, dt cut) | 4.4–7.0 GB |
| **BJ LU 64, s=1** | **1316 (10)** | **3080 (13)** | **3050 (47)** | **6697 (66)** | 4.2–5.7 GB |
| PB 1024, s=0 | 138 (100) | 159 (264) | 110 (417, cut) | 97 (416, cut) | 2.8 GB |
| PB 1024, s=1 | 781 (16) | 1178 (33) | 891 (146) | 1989 (203) | 2.8–3.0 GB |
| PB 1024, s=2 | 833 (15) | 1175 (33) | 1326 (99) | 2290 (174) | 2.8–3.0 GB |
| PB 1024 + 20 inner its, s=1 | 459 (11) | 874 (18) | 708 (76) | 1091 (151) | 2.9 GB |
| PB 64, s=1 | 663 (12) | 1029 (24) | 1042 (83) | 1498 (170) | 3.7–3.9 GB |

R0 p=3 (8.6e5 DOFs), µs/h (its): PB 2800, s=1: 163 (23) / 142 (83) / 126 (318) / 72 (280, cut), 11 GB;
PB 2800, s=0: 52 / 54 / 33 (cut) / 27 (cut), 11 GB; BJ LU 256, s=0: 105 (54) / 36 (223, cut) / 38 (265, cut) / fails,
16.5–18.9 GB.

- B starts from the equilibrium (the early linear phase) and includes the setup and first factorization in only 8
  steps. Combined with the accuracy limits of §6, only the dt ≤ 1e-6 to 2e-6 columns are usable for physics with s > 0.
- The large BJ LU 64 + s=1 advantage here (4× BJ LU s=0, 1.7× PB s=1 at dt 1e-6) is larger than D suggests. T2
  repeats this combination in the D restart segments.

**R1 (2.3e6 DOFs) weak scaling, same 4 × 2 layout, PB 8192 blocks** (partly *pending*):

| config | dt 1e-6 | dt 1e-5 | peak RSS |
|---|---|---|---|
| s=0 | 101 its, 181 s/step, 19.9 µs/h | *pending* | 20.4 GB |
| s=1 | 26 its, 52 s/step, 69.2 µs/h | *pending* | 20.3 GB |

At dt 1e-6 with s=1, R0 → R1 raises its/step from 16 to 26 (the CFL doubles) and the step time by 11× for 8× the
DOFs. s=1 is 3.5× faster than s=0 at R1, compared with 5.7× at R0. R1 strong scaling at 8 cores (2×4, 8×1) is
*pending*.

---

## 8. Coax railgun (rotational periodicity)

`results/m2max/railgun_correctness.txt`. Correctness runs: 3 ranks × 1 thread on shared cores. Wall times are
indicative only; the exclusive timing set is *pending* (§9).

| case | dt | its/step | Newton | sheet speed r = 0.11/0.15/0.18 m (km/s) | dN / dB vs shipped at ~10.5 µs | peak RSS |
|---|---|---|---|---|---|---|
| example solver (GMRES(2) + BJ LU, 1 block/rank), shipped | 3e-8 | 16.1 | 1.00 | 88.3 / 89.0 / 88.4 | — | 7.6 GB |
| example solver | 1.5e-8 | 11.9 | 1.00 | 89.5 / 90.1 / 89.4 | 0.11 / 0.028 | 7.3 GB |
| example solver | 6e-8 | 24.0 | 1.00 | 85.3 / 85.4 / 85.4 | 0.34 / 0.074 | 7.6 GB |
| PB, s=0 | 3e-8 | 5.9 | 1.00 | 88.3 / 89.0 / 88.4 | 8e-6 / 8e-6 | 4.9 GB |
| PB, s=1 | 3e-8 | 6.0 | 1.00 | 88.3 / 89.0 / 88.4 | 1.4e-3 / 3.6e-4 | 4.9 GB |
| PB, s=1 | 6e-8 | 9.2 | 1.00 | 85.4 / 85.4 / 85.4 | 0.34 / 0.074 | 4.9 GB |
| PB, s=0 | 1.2e-7 | 41.1 | 1.66 | 85.0 / 84.8 / 85.1 | 0.42 / 0.097 (9.6 µs) | 4.9 GB |
| PB, s=1 | 1.2e-7 | 40.0 | 1.67 | 84.6 / 84.6 / 85.4 | 0.43 / 0.105 | 4.9 GB |

- **Correctness: pass.** `pb_pre` type 1 (with `xmhd_local_bcs` and the revolved-periodic velocity packing) reproduces
  the example to ~1e-5 with s=0. The SI term adds ~1e-3 at the shipped dt, and the sheet speed is unchanged.
- **The PB/SI solution at a given dt matches the example's own solver at that dt** (6e-8: 85.4 vs 85.3–85.4 km/s,
  identical field differences). All the error at larger dt is time-step error, present with either solver.
- **The shipped dt is not fully converged.** Halving it raises the sheet speed by ~1.2 km/s (+1.3 %) and changes the
  density by ~11 % (L2, from a slightly shifted front). Speed changes per halving (3.1, then 1.2 km/s) suggest order
  ~1.4 and a zero-dt limit near 90 km/s. Above 6e-8 the speed barely changes (~85 km/s), so the large-dt end is not in
  the asymptotic regime.
- **PB s=0 at 1.2e-7** reduced dt on 15 of 88 steps (recovering automatically) and reported negative-density cells.
  s=1 held dt with the same iteration count. SI makes large dt stable here but not more accurate.
- **Cost:** PB cuts iterations 2.5–2.7× (16 → 6 at 3e-8; 24 → 9 at 6e-8) and memory from 7.6 to 4.9 GB. Wall time on
  shared cores was not better (2001–2266 s vs 1591 s at 3e-8); exclusive timing is *pending*.

---

## 9. Remaining runs (queued, run automatically one at a time at 4 × 2)

| set | cases | purpose | status |
|---|---|---|---|
| B R1 weak scaling | `R1_p2_pb8192_s{0,1}_dt1e-5` | weak scaling at large dt | running |
| B R1 strong scaling | `R1_p2_pb8192_s1_dt1e-5_np{2t4,8t1}` | ranks vs threads at fixed 8 cores | queued |
| Railgun timing (`Trail.list`) | example solver with 1 and 4 blocks/rank at 3e-8; PB 330/rank with s=0 and s=1 at 3e-8; PB s=1 at 6e-8 and 1.2e-7; PB 32/rank with s=1 at 3e-8 and 1.2e-7 | exclusive wall time for the railgun | queued |
| T2 (`T2.list`) | D segments at all three phases for BJ LU 64 + s=1 (dt 1e-6), BJ LU 256 + s=1 (1e-6), BJ LU 256 + s=0 (1e-6), BJ LU 64 + s=0.5 (3e-6); B: R0 p3 BJ LU 256 + s=1 at 4 dts | confirm finding 8 (SI + full-Jacobian LU) and fill the missing p=3 comparison | queued |

Expected to finish within about 10 h of this writing. Generate the tables with the `tab_*.py` scripts and update §7.2,
§8 and this section.

---

## 10. Suggested follow-ups (ranked)

1. **Settle SI + full-Jacobian block-Jacobi LU (T2).** If B's advantage holds in the D phases, the practical
   recommendation becomes "`xmhd_si` with a moderate-block full-Jacobian LU", and `pb_pre` is mainly a memory-saving
   option.
2. **A stronger, scalable inner solve on P.** With an exact LU of P the outer iteration count is 5–14 at every
   p, h and CFL (A1x), so this is where the remaining cost is. Options: larger P blocks with MUMPS (no 2 GiB cap),
   p-multigrid on P (`pb1_mg` ran in Tier 0 but was not benchmarked), and an auxiliary-space (Hiptmair–Xu/AMS) solver for
   the H(curl) B block (REPORT §8 item 1).
3. **Choose dt from accuracy, not stability.** On the tilt, s=1 needs dt ≲ 0.8 µs for 1 % in γ; s=0.5 gives about
   1.4× more dt at the same error. Consider documenting Δγ/γ ≈ −1.5 %·s·(dt/1 µs)² style guidance (or (1/12 + 1.1·s/8)(ω·dt)²
   in general) next to `xmhd_si`.
4. **Use MUMPS (or `umfpack_dl`) for large LU blocks.** `umfpack_di` fails above 2 GiB of factor memory (R1 p4 exact
   LU of P), well below the machine limit.
5. **Railgun example:** consider halving the shipped dt (or noting the ~1–2 % sheet-speed error at 3e-8), and offering
   the `pb_pre` type 1 XML as a lower-memory option (7.6 → 4.9 GB, 2.7× fewer iterations), pending the exclusive timing.
6. **MUG restart behavior:** a restart starts from the restart file's dt and ramps up to the input dt over ~7 steps.
   Benchmarks have to account for this, and users changing dt at a restart may not expect it. Consider an option to
   start at the input dt.
7. **The A2 p=4 small-dt anomaly:** rerun N=500 with a tighter `nl_tol` to confirm that the per-step Newton tolerance
   accumulates.
8. **Larger machines (HANDOFF tiers T2/T3):** R2 meshes, p=3 at R1, the C2 reference at R1, and strong scaling beyond 8
   cores were out of reach here (A1 R2 was not run: its memory estimate was uncertain).
9. **A regression test for `xmhd_si`** (linear cube, `pb_pre` type 1, iteration cap), as in REPORT §8 item 0.
