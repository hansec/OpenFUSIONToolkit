# Notes: Apple M5 Max, macOS 26.7.3

## Machine and run conditions

- M5 Max (Mac17,7), 18 cores, 128 GB, on AC power. The benchmark ran with no other builds
  or test runs active.
- **Core labels**: the M5 Max has no efficiency cores. `hw.perflevel0` is 6 "Super" cores and
  `hw.perflevel1` is 12 "Performance" cores. `sysinfo.txt` labels these `perf_cores: 6` and
  `eff_cores: 12` because `run_bench.sh` assumes perflevel1 means efficiency cores.
  `NTHREADS` used the default of 6 (the Super cores), so the threaded cases leave the 12
  Performance cores idle.
- The threaded OpenBLAS (step 4) was skipped, so there is no `OpenBLAS-OMP Nt` column. All
  OpenBLAS numbers are OFT's default serial build.
- No `timeout`/`gtimeout` command was installed, so the sparse LU runs had no time limit.
  None took longer than 11 s.

## Build notes

- Homebrew's default `python3` is 3.15, which has no h5py wheel yet (a source build failed),
  so the venv used `python3.14`.
- `--use_accelerate` built cleanly with SDK 26.5. No link errors or missing legacy symbols.

## Correctness

- Micro-benchmark checks are all below 3.2e-12 (dgemm at N=4000).
- Dense ThinCurr eigenvalues and frequency-response norms agree with OpenBLAS to 1.3e-13.
- HODLR eigenvalues agree to 7.2e-5.
- All sparse LU cases report `ok` with 1 GMRES iteration and norm differences of 1.3e-13 or
  less, except UMFPACK at 228k unknowns. That case fails in both trees with
  `UMFPACK: Out of Memory in factorization`, which is the expected 32-bit interface limit.
- No Accelerate sparse LU crashes. The two workarounds in `accelerate_bridge.c` (one RHS at a
  time for transposed solves, serialized symbolic factorization) were not separately
  retested on macOS 26.7.

## Differences from the M2 Max reference

- **The dgemv regression does not reproduce.** Accelerate dgemv is 8 to 10x faster than serial
  OpenBLAS at N ≤ 1000 and 2.5x faster at N ≥ 2000 (44 vs 18 GFLOP/s). On the M2 Max it was
  0.44 to 0.64x at N ≥ 2000. The drop from about 175 to 44 GFLOP/s between N=1000 and 2000
  is still there, but it stays well above OpenBLAS.
- **Accelerate threading has almost no effect.** `Accel Nt` (6 threads) matches `Accel 1t`
  for N ≥ 1000 in every kernel. In ThinCurr, `VECLIB_MAXIMUM_THREADS=1` gives the same times
  as 6 threads (freq_resp_direct 2.80 s in both cases; on the M2 Max it was 4.1 vs 2.0 s). A single
  thread appears to saturate the matrix unit on this chip.
- ThinCurr direct eigenvalues: 35.0 s → 15.1 s (2.3x; the M2 Max got 1.35x), probably from
  faster dgeev/dsyevd. Frequency response: 15.1 s → 2.8 s (5.4x). Time-dependent: 2.66 →
  1.23 s.
- Sparse LU, first solve at 228k unknowns: Accelerate 6.8 to 7.3 s, STRUMPACK 3.7 s with
  Accelerate BLAS and 7.1 s with OpenBLAS. With Accelerate BLAS, STRUMPACK is still about
  2x faster than Accelerate's sparse LU. Accelerate BLAS speeds up STRUMPACK 1.9x here (1.4x
  on the M2 Max). Supernodal pivoting with METIS is the fastest Accelerate variant at 228k,
  but only by a few percent.
