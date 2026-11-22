# Regression tests: Apple M2 Max, macOS 27.2

Compilers: Homebrew GCC 16.2.0 (`gcc-16`, `g++-16`, `gfortran-16`).

`build_libs.py` flags, identical in both trees except for `--use_accelerate`:

```
--nthread=6 --build_mpich=1 --build_arpack=1 --build_superlu=1 --build_superlu_dist=1
--build_umfpack=1 --build_strumpack=1 --oft_build_tests=1 --no_dl_progress [--use_accelerate]
```

The OpenBLAS tree uses the default serial, static OpenBLAS 0.3.33. Neither tree has PETSc.

`OMP_NUM_THREADS=2 make test`, with the two trees run one after the other:

| tree | result |
|---|---|
| Accelerate | 1148 passed, 349 skipped, 229 deselected, 0 failed in 1491 s (24:51) |
| OpenBLAS | 1148 passed, 349 skipped, 229 deselected, 0 failed in 1505 s (25:05) |

- The 349 skips are PETSc-only tests.
- One timing difference reproduced across repeats: `physics/MUG/test_alfven.py::test_nl_r2_p2_gmg` took 54.6 s with OpenBLAS and 43.4 s with Accelerate.
