# Accelerate vs OpenBLAS benchmark for OFT (macOS)

This directory holds a repeatable procedure for comparing an OFT build linked to
Apple's Accelerate BLAS/LAPACK (`build_libs.py --use_accelerate`) with the
default OpenBLAS build, on any Mac. It is written so that a coding agent can run
it start to finish on a new machine and commit the results here.

Contents:

| File | Purpose |
|---|---|
| `run_bench.sh` | Driver: builds and runs the kernel micro-benchmark, runs the ThinCurr cases, writes `summary.md` |
| `micro.f90` | Kernel micro-benchmark (dgemm, dgetrf, dgetri, dsyevd, dgeev, zgemm, dgemv) with accuracy checks |
| `bench_thincurr.py` | ThinCurr benchmark on a generated torus mesh (dense direct and HODLR paths) |
| `summarize.py` | Turns a results directory into Markdown tables (called by `run_bench.sh`) |
| `results/<machine>/` | Committed results, one directory per machine |

## Background (read before debugging)

- `--use_accelerate` needs macOS 13.3 or later (arm64 or x86_64). It generates
  `<config_dir>/accelerate_shim/lib/liboftaccelerate.dylib`, which maps the
  standard Fortran names (`dgemm_`, ...) to Accelerate's LAPACK 3.11 entry points
  (`dgemm$NEWLAPACK`). Accelerate's plain `dgemm_` symbols are an old LAPACK
  (3.2.1) with f2c calling conventions that gfortran does not match, so OFT
  never links them directly.
- All external libraries (ARPACK, SuperLU, SuperLU-DIST, UMFPACK, STRUMPACK, ...)
  and OFT itself link against the shim. OFT gets it through the CMake variable
  `OFT_LAPACK_LIBRARIES`.
- The default OFT OpenBLAS is built **serial** (no threads). Accelerate is threaded
  and controlled by `VECLIB_MAXIMUM_THREADS` (default: all cores). The benchmark
  pins both explicitly. For a fair threaded comparison, also build a threaded
  OpenBLAS (step 4, optional).

## Procedure

Use the repository's agent instructions (`AGENTS.md` and the `oft-build` and
`oft-test` skills) for general build and test conventions. In particular:
ask the user which Python environment to use if they have not said, and do not
start the 20 to 30 minute library builds without confirming with the user.

### 1. Prerequisites

- macOS ≥ 13.3, Xcode Command Line Tools (`xcrun --show-sdk-path` works).
- GCC with gfortran (eg. Homebrew `gcc`). `gcc` on macOS is Apple clang, so pass
  the versioned names explicitly (`gcc-15`, `gfortran-15`, ... whatever is installed).
- CMake ≥ 3.27.
- A Python environment with `pytest numpy scipy h5py netcdf4 matplotlib`.

### 2. Get the branch

```bash
git fetch https://github.com/hansec/OpenFUSIONToolkit.git build-accelerate
git checkout -b build-accelerate FETCH_HEAD    # or: git worktree add ...
```

The `--use_accelerate` option lives in `src/utilities/build_libs.py` on this
branch. Use this checkout's `build_libs.py` for **both** trees, so the only
difference between them is BLAS.

### 3. Build two configurations

Create two sibling directories (anywhere gitignored, eg. `builds/GCC/`) and run
`build_libs.py` in each, with the Python environment active. Substitute the
installed GCC version and a sensible `--nthread` (number of cores minus a few).
The two builds can run at the same time.

```bash
SRC=<repo_root>/src
export CC=gcc-NN CXX=g++-NN FC=gfortran-NN
COMMON="--nthread=6 --build_mpich=1 --build_arpack=1 --build_superlu=1 \
  --build_superlu_dist=1 --build_umfpack=1 --build_strumpack=1 \
  --oft_build_tests=1 --no_dl_progress"

mkdir -p builds/GCC/MPI_accelerate builds/GCC/MPI_openblas
(cd builds/GCC/MPI_accelerate && python $SRC/utilities/build_libs.py $COMMON --use_accelerate > build_libs.log 2>&1) &
(cd builds/GCC/MPI_openblas   && python $SRC/utilities/build_libs.py $COMMON                  > build_libs.log 2>&1) &
wait
```

Then, in each directory:

```bash
bash config_cmake.sh > cmake_config.log 2>&1
cd build_release && make -j$(getconf _NPROCESSORS_ONLN)
```

Check the Accelerate tree before going further:

```bash
grep OFT_LAPACK_LIBRARIES builds/GCC/MPI_accelerate/config_cmake.sh          # points at liboftaccelerate.dylib
grep 'Using specified BLAS/LAPACK' builds/GCC/MPI_accelerate/cmake_config.log
otool -L builds/GCC/MPI_accelerate/build_release/bin/liboftpy.dylib | grep liboftaccelerate
```

If the shim step of `build_libs.py` fails, its log is copied to
`builds/GCC/MPI_accelerate/build/build_error.log` (eg. a missing SDK, or an SDK
whose `vecLib.framework/Versions/A/libBLAS.tbd` has no `$NEWLAPACK` symbols).
An undefined `_xxx_` symbol at link time means a routine Accelerate only exports
under its old name. Report it rather than working around it: the fix belongs in
the `legacy_funcs` table of the `ACCELERATE` class in `build_libs.py`.

### 4. (Optional) threaded OpenBLAS for the micro-benchmark

Only needed for the `OpenBLAS-OMP Nt` column. Build it outside the repository:

```bash
curl -L -o OpenBLAS-0.3.33.tar.gz https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.33/OpenBLAS-0.3.33.tar.gz
tar xzf OpenBLAS-0.3.33.tar.gz && cd OpenBLAS-0.3.33
make CC=gcc-NN FC=gfortran-NN NO_CBLAS=1 NO_LAPACKE=1 NO_SHARED=1 USE_THREAD=1 USE_OPENMP=1 -j8 libs netlib
# result: OpenBLAS-0.3.33/libopenblas.a  -> pass as OMP_OPENBLAS_LIB below
```

### 5. Run the regression tests in both trees

Per `AGENTS.md`, always set `OMP_NUM_THREADS=2`. Each run takes about 25
minutes on an M2 Max. The two runs can share the machine, but they then compete
for cores, so don't use the timings from simultaneous runs.

```bash
cd builds/GCC/MPI_accelerate/build_release && OMP_NUM_THREADS=2 make test > test.log 2>&1
cd builds/GCC/MPI_openblas/build_release   && OMP_NUM_THREADS=2 make test > test.log 2>&1
tail -3 */build_release/test.log      # pytest summary line
```

Expected: identical pass/skip counts in both trees and **no failures** in the
Accelerate tree beyond any that also fail with OpenBLAS. Without PETSc, about
349 tests skip (`Not compiled with PETSc`). For any failure that appears only
with Accelerate, rerun it alone (`cd tests && ./run_test.sh <file> -k <name>`)
and record the output.

### 6. Run the benchmark

The machine should be otherwise idle: close other builds and test runs, and keep
it on AC power.

```bash
cd <repo_root>/benchmarks/blas
PY=$(which python) NTHREADS=<perf cores> OMP_OPENBLAS_LIB=<path>/libopenblas.a \
  ./run_bench.sh <abs>/builds/GCC/MPI_accelerate <abs>/builds/GCC/MPI_openblas \
  results/<Chip>_macOS<NN>
```

- Name the results directory after the chip and OS, eg. `AppleM4Pro_macOS26`
  or `Intel_i9_macOS15`.
- `NTHREADS` defaults to the number of performance cores
  (`sysctl hw.perflevel0.physicalcpu`), which is usually the right choice.
- Omit `OMP_OPENBLAS_LIB` if you skipped step 4.
- Set `SKIP_LARGE=1` on machines with less than 16 GB of memory.
- The run takes about 5 minutes on an M2 Max. Expect longer on slower machines,
  mostly from the serial OpenBLAS cases.

What it runs:

| Output | Description |
|---|---|
| `sysinfo.txt` | Machine, OS, SDK, compiler, OFT revision, and the Accelerate tree's BLAS linkage |
| `micro_*.txt` | Kernels at N = 500 to 4000, best of 3 (1 run at N=4000). Columns: time, GFLOP/s, accuracy check |
| `tc_<case>.json/.log` | ThinCurr on a 48x96 torus (4609 unknowns): `Lmat` (no BLAS), `eigs_direct` (dgemm + dgeev), `freq_resp_direct` (complex LU/inverse), `td_direct` (200 time steps). Also the HODLR path at 4609 and 18433 unknowns |
| `summary.md` | Tables produced by `summarize.py` |

Cases: `*_tN` use NTHREADS OpenMP threads in OFT. `acc_tN_veclib1` keeps OFT
threaded but limits Accelerate to 1 thread. `*_t1` are fully serial.

### 7. Check, record and commit

1. Read `summary.md` and check that:
   - All micro-benchmark `max check` values are small: below about 1e-10. The dgemm, dgetrf and dsyevd checks are
     absolute residuals, so they grow with N (the M2 Max reference reaches 2.7e-12 for dgemm at N=4000).
   - Dense ThinCurr eigenvalues and frequency-response norms agree with OpenBLAS
     to better than about 1e-10 relative.
   - HODLR eigenvalues agree to about 1e-3, which is within the compression tolerance.

   Anything larger is a correctness problem. Report it.
2. Add a `tests.md` to the results directory containing:
   - both trees' pytest summary lines (passed/skipped/failed and wall time),
   - any failures, with the failing test names,
   - the exact `build_libs.py` flags and compiler versions used.
3. Add a short `notes.md` for anything unusual: build problems, workarounds,
   thermal throttling, or another load on the machine.
4. Commit only the results directory (not build trees or scratch files), then push:

```bash
git add benchmarks/blas/results/<machine>
git commit -m "Benchmarks: Accelerate vs OpenBLAS results for <machine>"
git push origin build-accelerate
```

## Reference results

`results/AppleM2Max_macOS27/` (M2 Max, 8 performance + 4 efficiency cores,
64 GB, gcc 16):

- Accelerate is 4 to 8x faster than OFT's serial OpenBLAS on dgemm, zgemm, dgetrf
  and dgetri, even on one thread.
- With 8 threads it is about 1.5 to 2x faster than threaded OpenBLAS, with
  similar dsyevd and dgeev times.
- dgemv is 3x faster for N ≤ 1000 but 1.5 to 2.3x **slower** at N ≥ 2000.
- ThinCurr direct frequency response: 19.2 s → 2.1 s. Direct eigenvalues:
  44 s → 33 s.
- The HODLR path shows no meaningful difference.
- `make test`: 1148 passed and 349 skipped in both trees, with no failures.

## Interpreting results on other machines

- **Intel Macs**: Accelerate on x86_64 is a different code path, and OpenBLAS is
  much better tuned for x86 (AVX2/AVX-512). Expect a much smaller advantage,
  possibly none.
- **Newer Apple Silicon (M3/M4/M5)**: Accelerate uses the AMX/SME matrix units,
  so the BLAS-3 advantage should hold or grow.
- **dgemv regression**: check whether the large-N dgemv slowdown reproduces. It
  matters for ThinCurr time stepping with a dense L matrix.
