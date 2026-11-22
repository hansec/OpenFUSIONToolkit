# Regression tests: Apple M5 Max, macOS 26.7.3

Compilers: Homebrew GCC 16.2.0 (`gcc-16`, `g++-16`, `gfortran-16`). Python 3.14.8
(Homebrew) venv with pytest, numpy, scipy, h5py, netcdf4, matplotlib. CMake 4.4.4,
macOS SDK 26.5 (Command Line Tools).

`build_libs.py` flags, identical in both trees except for `--use_accelerate`:

```
--nthread=8 --build_mpich=1 --build_arpack=1 --build_superlu=1 --build_superlu_dist=1
--build_umfpack=1 --build_strumpack=1 --oft_build_tests=1 --no_dl_progress [--use_accelerate]
```

The OpenBLAS tree uses the default serial, static OpenBLAS 0.3.33. Neither tree has PETSc.
CMake reported `Using Accelerate sparse solvers` in the Accelerate tree, so the
`accelerate` sparse LU package was built and tested.

`OMP_NUM_THREADS=2 make test`, with the two trees run one after the other:

| tree | result |
|---|---|
| Accelerate | 1175 passed, 349 skipped, 229 deselected, 0 failed in 945 s (15:45) |
| OpenBLAS | 1148 passed, 376 skipped, 229 deselected, 0 failed in 936 s (15:36) |

- The 349 skips common to both trees are PETSc-only tests.
- The 27 extra skips in the OpenBLAS tree are the Accelerate sparse LU tests
  (`lin_alg/test_accelerate.py`, 18, and the `[7]` cases of `lin_alg/test_native_bjacobi.py`, 9),
  which all pass in the Accelerate tree.
- No failures in either tree.
