## System
```
date: 2026-11-22T22:02Z
host_model: Mac17,7
cpu: Apple M5 Max
arch: arm64
ncpu: 18
perf_cores: 6
eff_cores: 12
mem_GB: 128
macos: 26.7.3
sdk: 26.5
fc: /opt/homebrew/bin/gfortran-16 (GNU Fortran (Homebrew GCC 16.2.0) 16.2.0)
python: Python 3.14.8
oft_rev: f6f1287e (build-accelerate)
openblas_lib: /Users/hansec/Claude/Repos/OpenFUSIONToolkit/builds/GCC/MPI_openblas/OpenBLAS-0_3_33/lib/libopenblas.a
omp_openblas_lib: none
nthreads: 6
--- Accelerate liboftpy BLAS linkage (expect liboftaccelerate) ---
	/Users/hansec/Claude/Repos/OpenFUSIONToolkit/builds/GCC/MPI_accelerate/accelerate_shim/lib/liboftaccelerate.dylib (compatibility version 0.0.0, current version 0.0.0)
	/System/Library/Frameworks/Accelerate.framework/Versions/A/Accelerate (compatibility version 1.0.0, current version 4.0.0)
```

## BLAS/LAPACK kernels

Time in seconds (dgemv: GFLOP/s, higher is better). "Speedup" is OpenBLAS 1t / Accel 1t. Max accuracy check across columns.

| kernel | N | OpenBLAS 1t | Accel 1t | Accel Nt | speedup 1t | max check |
|---|---|---|---|---|---|---|
| dgemm | 500 | 0.0063 | 0.0008 | 0.0005 | 7.88x | 2.3e-13 |
| dgetrf | 500 | 0.0022 | 0.0010 | 0.0006 | 2.20x | 2.4e-15 |
| dgetri | 500 | 0.0041 | 0.0011 | 0.0007 | 3.73x | 8.9e-16 |
| dsyevd | 500 | 0.0223 | 0.0149 | 0.0096 | 1.50x | 1.8e-13 |
| dgeev | 500 | 0.0677 | 0.0559 | 0.0509 | 1.21x | 1.3e-15 |
| zgemm | 500 | 0.0187 | 0.0020 | 0.0020 | 9.35x | 1.3e-15 |
| dgemv | 500 | 17.9 | 145.3 | 168.9 | 8.12x | 0.0e+00 |
| dgemm | 1000 | 0.0334 | 0.0037 | 0.0038 | 9.03x | 9.1e-13 |
| dgetrf | 1000 | 0.0124 | 0.0026 | 0.0026 | 4.77x | 3.7e-15 |
| dgetri | 1000 | 0.0250 | 0.0039 | 0.0039 | 6.41x | 1.7e-15 |
| dsyevd | 1000 | 0.1267 | 0.0425 | 0.0427 | 2.98x | 3.9e-13 |
| dgeev | 1000 | 0.4070 | 0.2388 | 0.2405 | 1.70x | 1.3e-15 |
| zgemm | 1000 | 0.1496 | 0.0176 | 0.0176 | 8.50x | 5.2e-16 |
| dgemv | 1000 | 17.8 | 175.1 | 186.9 | 9.83x | 0.0e+00 |
| dgemm | 2000 | 0.2669 | 0.0349 | 0.0350 | 7.65x | 6.8e-13 |
| dgetrf | 2000 | 0.0938 | 0.0181 | 0.0181 | 5.18x | 5.7e-15 |
| dgetri | 2000 | 0.1948 | 0.0323 | 0.0326 | 6.03x | 1.2e-15 |
| dsyevd | 2000 | 0.9228 | 0.2949 | 0.2978 | 3.13x | 1.6e-12 |
| dgeev | 2000 | 2.6647 | 1.1223 | 1.1681 | 2.37x | 4.4e-15 |
| zgemm | 2000 | 1.1924 | 0.1469 | 0.1460 | 8.12x | 3.1e-15 |
| dgemv | 2000 | 18.0 | 44.4 | 44.6 | 2.47x | 0.0e+00 |
| dgemm | 4000 | 2.1482 | 0.2983 | 0.3124 | 7.20x | 3.2e-12 |
| dgetrf | 4000 | 0.7387 | 0.1370 | 0.1395 | 5.39x | 1.0e-14 |
| dgetri | 4000 | 1.5627 | 0.2661 | 0.2674 | 5.87x | 2.4e-15 |
| zgemm | 4000 | 9.6081 | 1.2350 | 1.2209 | 7.78x | 7.9e-15 |
| dgemv | 4000 | 17.7 | 43.8 | 43.7 | 2.47x | 0.0e+00 |

## ThinCurr dense (direct) solvers

Torus 48x96 (4609 unknowns), time in seconds.

| case | Lmat | eigs_direct | freq_resp_direct | td_direct |
|---|---|---|---|---|
| OpenBLAS, OFT N threads | 1.565 | 35.043 | 15.140 | 2.664 |
| Accelerate, OFT N threads | 1.592 | 15.134 | 2.803 | 1.230 |
| Accelerate, OFT N threads, VECLIB_MAXIMUM_THREADS=1 | 1.642 | 15.048 | 2.805 | 1.230 |
| OpenBLAS, 1 thread | 8.176 | 35.226 | 15.146 | 3.799 |
| Accelerate, 1 thread | 8.171 | 15.212 | 2.802 | 2.420 |

Agreement with OpenBLAS (max relative difference):

- OpenBLAS, OFT N threads: eigenvalues 0.0e+00, freq. response norm 0.0e+00
- Accelerate, OFT N threads: eigenvalues 1.1e-13, freq. response norm 5.0e-14
- Accelerate, OFT N threads, VECLIB_MAXIMUM_THREADS=1: eigenvalues 1.1e-13, freq. response norm 5.0e-14
- OpenBLAS, 1 thread: eigenvalues 9.3e-14, freq. response norm 3.7e-15
- Accelerate, 1 thread: eigenvalues 1.3e-13, freq. response norm 5.2e-14

## ThinCurr HODLR (compressed) path, OFT N threads

| unknowns | Lmat_hodlr OpenBLAS | Lmat_hodlr Accel | eigs_iterative OpenBLAS | eigs_iterative Accel | eig rel. diff |
|---|---|---|---|---|---|
| 4609 | 1.525 | 1.540 | 0.057 | 0.048 | 4.9e-05 |
| 18433 | 7.857 | 8.066 | 0.258 | 0.184 | 7.2e-05 |

## Sparse direct (LU) solvers, OFT N threads

Order-3 Lagrange Laplacian on a unit cube (`lin_alg/test_native_bjacobi`, GMRES with a full LU preconditioner). "solve" is the first solve including factorization, in seconds; "wall" also includes mesh and matrix setup. "norm diff" is the relative difference of the solution norm from the first successful case.

### 29449 unknowns (nlevels=4)

| package | status | solve | wall | GMRES its | norm diff |
|---|---|---|---|---|---|
| Accelerate sparse LU (default options) | ok | 0.40 | 0.59 | 1 | 0.0e+00 |
| Accelerate sparse LU, METIS ordering | ok | 0.43 | 0.62 | 1 | 6.4e-16 |
| Accelerate sparse LU, supernodal pivoting + METIS | ok | 0.42 | 0.61 | 1 | 6.4e-16 |
| STRUMPACK (Accelerate BLAS) | ok | 0.24 | 0.43 | 1 | 3.2e-15 |
| STRUMPACK (OpenBLAS) | ok | 0.28 | 0.47 | 1 | 4.1e-15 |
| UMFPACK (Accelerate BLAS) | ok | 0.59 | 0.78 | 1 | 4.8e-14 |
| UMFPACK (OpenBLAS) | ok | 1.14 | 1.33 | 1 | 3.2e-15 |
| SuperLU (Accelerate BLAS) | ok | 10.55 | 10.74 | 1 | 1.3e-14 |

### 228241 unknowns (nlevels=5)

| package | status | solve | wall | GMRES its | norm diff |
|---|---|---|---|---|---|
| Accelerate sparse LU (default options) | ok | 6.94 | 7.75 | 1 | 0.0e+00 |
| Accelerate sparse LU, METIS ordering | ok | 7.30 | 8.12 | 1 | 2.7e-14 |
| Accelerate sparse LU, supernodal pivoting + METIS | ok | 6.81 | 7.62 | 1 | 2.7e-14 |
| STRUMPACK (Accelerate BLAS) | ok | 3.67 | 4.52 | 1 | 1.3e-13 |
| STRUMPACK (OpenBLAS) | ok | 7.11 | 7.95 | 1 | 1.9e-13 |
| UMFPACK (Accelerate BLAS) | failed | — | 5.12 | — | — |
| UMFPACK (OpenBLAS) | failed | — | 11.03 | — | — |

