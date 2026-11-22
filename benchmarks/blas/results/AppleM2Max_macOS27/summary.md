## System
```
date: 2026-11-22T13:04Z
host_model: Mac14,6
cpu: Apple M2 Max
arch: arm64
ncpu: 12
perf_cores: 8
eff_cores: 4
mem_GB: 64
macos: 27.2.1
sdk: 27.2
fc: /opt/homebrew/bin/gfortran-16 (GNU Fortran (Homebrew GCC 16.2.0) 16.2.0)
python: Python 3.14.10
oft_rev: 29e9f58 (build-accelerate)
openblas_lib: /Users/hansec/Agents/Claude/OpenFUSIONToolkit/builds/GCC/MPI_openblas/OpenBLAS-0_3_33/lib/libopenblas.a
omp_openblas_lib: /private/tmp/claude-501/-Users-hansec-Agents-Claude-OpenFUSIONToolkit/4c7413a7-e06a-4de9-854d-f69c3abd0001/scratchpad/bench/OpenBLAS-0.3.33/libopenblas.a
nthreads: 8
--- Accelerate liboftpy BLAS linkage (expect liboftaccelerate) ---
	/Users/hansec/Agents/Claude/OpenFUSIONToolkit/builds/GCC/MPI_accelerate/accelerate_shim/lib/liboftaccelerate.dylib (compatibility version 0.0.0, current version 0.0.0)
```

## BLAS/LAPACK kernels

Time in seconds (dgemv: GFLOP/s, higher is better). "Speedup" is OpenBLAS 1t / Accel 1t. Max accuracy check across columns.

| kernel | N | OpenBLAS 1t | Accel 1t | OpenBLAS-OMP Nt | Accel Nt | speedup 1t | max check |
|---|---|---|---|---|---|---|---|
| dgemm | 500 | 0.0061 | 0.0010 | 0.0014 | 0.0004 | 6.10x | 0.0e+00 |
| dgetrf | 500 | 0.0025 | 0.0012 | 0.0027 | 0.0009 | 2.08x | 2.6e-15 |
| dgetri | 500 | 0.0071 | 0.0012 | 0.0040 | 0.0010 | 5.92x | 8.9e-16 |
| dsyevd | 500 | 0.0213 | 0.0183 | 0.0328 | 0.0158 | 1.16x | 1.4e-13 |
| dgeev | 500 | 0.0855 | 0.0733 | 0.1042 | 0.0740 | 1.17x | 1.0e-15 |
| zgemm | 500 | 0.0228 | 0.0027 | 0.0036 | 0.0014 | 8.44x | 1.1e-15 |
| dgemv | 500 | 12.6 | 36.2 | 12.6 | 35.4 | 2.88x | 0.0e+00 |
| dgemm | 1000 | 0.0371 | 0.0055 | 0.0061 | 0.0028 | 6.75x | 6.8e-13 |
| dgetrf | 1000 | 0.0146 | 0.0034 | 0.0069 | 0.0029 | 4.29x | 3.7e-15 |
| dgetri | 1000 | 0.0349 | 0.0050 | 0.0095 | 0.0054 | 6.98x | 1.7e-15 |
| dsyevd | 1000 | 0.1304 | 0.0747 | 0.1036 | 0.0748 | 1.75x | 4.0e-13 |
| dgeev | 1000 | 0.4973 | 0.3482 | 0.4201 | 0.3667 | 1.43x | 1.7e-15 |
| zgemm | 1000 | 0.1814 | 0.0228 | 0.0297 | 0.0114 | 7.96x | 5.8e-16 |
| dgemv | 1000 | 12.5 | 37.3 | 11.8 | 37.1 | 2.98x | 0.0e+00 |
| dgemm | 2000 | 0.2927 | 0.0456 | 0.0441 | 0.0238 | 6.42x | 4.5e-13 |
| dgetrf | 2000 | 0.1070 | 0.0256 | 0.0281 | 0.0144 | 4.18x | 5.7e-15 |
| dgetri | 2000 | 0.2319 | 0.0441 | 0.0531 | 0.0337 | 5.26x | 1.2e-15 |
| dsyevd | 2000 | 0.9146 | 0.4555 | 0.3865 | 0.4416 | 2.01x | 8.9e-13 |
| dgeev | 2000 | 3.1871 | 1.8738 | 1.8401 | 1.8997 | 1.70x | 4.0e-15 |
| zgemm | 2000 | 1.4477 | 0.1949 | 0.1993 | 0.0965 | 7.43x | 3.5e-15 |
| dgemv | 2000 | 12.4 | 8.0 | 15.3 | 8.3 | 0.64x | 0.0e+00 |
| dgemm | 4000 | 2.3470 | 0.3935 | 0.3549 | 0.1917 | 5.96x | 2.7e-12 |
| dgetrf | 4000 | 0.8178 | 0.2094 | 0.1560 | 0.1023 | 3.91x | 1.0e-14 |
| dgetri | 4000 | 1.7921 | 0.3663 | 0.4768 | 0.2489 | 4.89x | 2.4e-15 |
| zgemm | 4000 | 11.6701 | 1.8049 | 1.6189 | 0.7868 | 6.47x | 7.9e-15 |
| dgemv | 4000 | 11.9 | 5.2 | 13.7 | 5.3 | 0.44x | 0.0e+00 |

## ThinCurr dense (direct) solvers

Torus 48x96 (4609 unknowns), time in seconds.

| case | Lmat | eigs_direct | freq_resp_direct | td_direct |
|---|---|---|---|---|
| OpenBLAS, OFT N threads | 1.936 | 43.875 | 18.651 | 3.935 |
| Accelerate, OFT N threads | 1.947 | 32.586 | 2.003 | 2.396 |
| Accelerate, OFT N threads, VECLIB_MAXIMUM_THREADS=1 | 1.935 | 36.343 | 4.146 | 2.994 |
| OpenBLAS, 1 thread | 13.386 | 44.385 | 19.157 | 5.350 |
| Accelerate, 1 thread | 13.389 | 36.324 | 4.168 | 4.083 |

Agreement with OpenBLAS (max relative difference):

- OpenBLAS, OFT N threads: eigenvalues 0.0e+00, freq. response norm 0.0e+00
- Accelerate, OFT N threads: eigenvalues 2.6e-13, freq. response norm 4.5e-14
- Accelerate, OFT N threads, VECLIB_MAXIMUM_THREADS=1: eigenvalues 2.6e-13, freq. response norm 4.7e-14
- OpenBLAS, 1 thread: eigenvalues 8.0e-14, freq. response norm 7.3e-15
- Accelerate, 1 thread: eigenvalues 2.6e-13, freq. response norm 4.9e-14

## ThinCurr HODLR (compressed) path, OFT N threads

| unknowns | Lmat_hodlr OpenBLAS | Lmat_hodlr Accel | eigs_iterative OpenBLAS | eigs_iterative Accel | eig rel. diff |
|---|---|---|---|---|---|
| 4609 | 1.809 | 1.831 | 0.090 | 0.073 | 4.6e-05 |
| 18433 | 8.925 | 9.026 | 0.447 | 0.307 | 2.0e-04 |

