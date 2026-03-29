# Open FUSION Toolkit (OFT) – Copilot Agent Instructions

Trust these instructions. Only search the codebase if something is incomplete or found to be incorrect.

## What This Repo Does
OFT is a finite-element simulation suite for plasma/fusion research. It provides four component tools:
**TokaMaker** (axisymmetric MHD equilibria), **ThinCurr** (thin-wall eddy currents), **MUG** (extended MHD dynamics), and **Marklin** (3-D force-free equilibria). The framework is built on variable-order FEM on unstructured tetrahedral/hexahedral grids.

## Languages & Frameworks
- **Fortran** (primary, ~118 `.F90` files, ~98k lines) – physics, FEM framework
- **C / C++** – low-level utilities, external-library wrappers
- **Python** (~77 files) – high-level interfaces, build tooling, tests
- **CMake** – build system (min 3.12); always driven through `build_libs.py`
- Supported compilers: GCC 12–15, Intel oneAPI (icx/ifx); OpenMP enabled by default

## Repository Layout
```
/                   – repo root: README.md, CONTRIBUTING.md, .gitignore, setup_env.sh
src/                – all source code (≈91 MB)
  CMakeLists.txt    – master CMake configuration
  base/             – Fortran runtime utilities (I/O, sorting, stitching)
  grid/             – unstructured grid management
  lin_alg/          – linear-algebra layer (solvers, ARPACK, PETSc, SuperLU, UMFPACK)
  fem/              – finite-element infrastructure (H1, HCurl, H(div), Lagrange)
  physics/          – physics modules: TokaMaker, ThinCurr, MUG, Marklin
  python/           – Python interfaces
    OpenFUSIONToolkit/  – importable package (TokaMaker, ThinCurr, Marklin, MUG APIs)
    pyproject.toml  – ruff linting config (target-version py37, ignores E722/F403/F405)
    wrappers/       – Fortran→Python shared-library wrappers
  tests/            – pytest test suite (run via CMake `make test`)
    conftest.py     – markers: slow, coverage, mpi, linear
    oft_testing.py  – shared test helpers
    base/ fem/ grid/ lin_alg/ physics/  – per-module test directories
  utilities/        – build scripts and standalone tools
    build_libs.py   – downloads and builds all external dependencies + configures OFT
    generate_stack.py – validates internal Fortran debug-stack entries
  examples/         – example cases (TokaMaker, ThinCurr, MUG, Marklin)
  bin/              – auxiliary executables
  ext_libs/         – CMake wrappers for external libraries
  include/          – generated Fortran include files
  cmake/            – custom CMake find-modules
  docs/             – Doxygen documentation sources
.github/
  workflows/
    ci_build.yaml         – CI: build + test on Ubuntu 22/24 + macOS (GCC & Intel)
    lint.yaml             – Lint: ruff + generate_stack.py check
    cov_build.yaml        – Coverage build with GCC --coverage
    cd_nightly.yaml       – Nightly packaging/release builds
    copilot-setup-steps.yml – Minimal Ubuntu 24.04/GCC-14 setup for Copilot agents
```

## Environment Setup (do this first, every time)
```bash
# Activate the Python venv (always required before any Python or build step)
source setup_env.sh
# Install Python dependencies into the venv
pip install pytest numpy scipy h5py matplotlib xarray jupyter nbconvert pyvista
```
`setup_env.sh` contains `source /path/to/oft_venv/bin/activate`. The venv lives in `oft_venv/` (gitignored).

## Building OFT

All build steps run from a **`libs/`** directory adjacent to `src/`. `build_libs.py` downloads, compiles external libraries, and generates `libs/config_cmake.sh` which configures OFT's CMake build.

```bash
mkdir -p libs && cd libs

# 1. Build external libraries + configure OFT (≈15–30 min first time, cached thereafter)
python ../src/utilities/build_libs.py \
    --oblas_dynamic_arch --build_umfpack=1 --build_superlu=1 \
    --build_arpack=1 --oft_build_tests=1 --no_dl_progress --nthread=4

# 2. Configure OFT via the generated script
bash config_cmake.sh

# 3. Build OFT (from libs/build_release/)
cd build_release
make -j4

# 4. Install OFT
make install
```

For MPI builds, add `--build_mpich=1` (or `--build_openmpi=1`) to step 1.  
The cache key for external libraries is based on `src/utilities/build_libs.py` and `compiler_id.txt`.

## Running Tests
```bash
cd libs/build_release
source ../../setup_env.sh   # venv must be active
make test                   # runs pytest with "-m 'not slow'" marker
make test_full              # runs all tests including slow ones
make test_examples          # runs example-based tests
```
Tests are pytest-based, orchestrated by CMake. JUnit XML output: `tests/OFT.junit.xml`.

## Linting (runs on every PR via lint.yaml)
Always run both checks before submitting a PR:

```bash
# 1. Check Fortran debug-stack entries (run from src/)
source setup_env.sh
cd src
python utilities/generate_stack.py -l

# 2. Check Python code with ruff (run from src/python/)
cd src/python
python -m ruff check
```
Ruff config: `src/python/pyproject.toml` – targets Python 3.7, ignores E722/F403/F405.

## CI Checks (must all pass)
| Workflow | Trigger | What it checks |
|---|---|---|
| `lint.yaml` | PR / push to main | ruff + generate_stack.py |
| `ci_build.yaml` | PR / push to main | Build + test matrix (Ubuntu 22/24 GCC, Intel, macOS GCC) |
| `cov_build.yaml` | PR / push to main | Coverage build with GCC |

## Key Notes
- **Never** commit `oft_venv/`, `libs/`, `build*/`, `install*/`, or `compiler_id.txt` – all are gitignored.
- Fortran sources use free-form format (`.F90`). Line length is unlimited (`-ffree-line-length-none`).  
  Keep lines ≤132 chars to avoid `-Werror=line-truncation` failures.
- Python code in `src/python/` must remain compatible with **Python 3.7** (ruff enforces this).
- Documentation is generated by Doxygen; add docstrings to new Fortran subroutines/functions using the `!>` syntax.
- Adding a new Fortran subroutine/function: if `OFT_DEBUG_STACK` is enabled, run `generate_stack.py` to register it; the lint check will fail if any stack entries are missing or stale.
