#!/usr/bin/env bash
# Benchmark an Accelerate-based OFT build against an OpenBLAS-based one.
#
# Usage: run_bench.sh <accelerate_config_dir> <openblas_config_dir> <results_dir>
#
# <*_config_dir> are build_libs.py output directories (holding config_cmake.sh and
# build_release/), identical except for --use_accelerate. Results (text/JSON/logs and
# summary.md) are written to <results_dir>; scratch files go to $WORK_DIR.
#
# Environment variables (all optional):
#   PY               Python with numpy/h5py used to build OFT (default: python3)
#   FC               Fortran compiler (default: the one in config_cmake.sh)
#   NTHREADS         Thread count for the multi-threaded runs (default: performance cores)
#   OMP_OPENBLAS_LIB Threaded (OpenMP) libopenblas.a for an extra micro-benchmark column
#   SKIP_MICRO=1     Skip the BLAS/LAPACK kernel micro-benchmark
#   SKIP_THINCURR=1  Skip the ThinCurr benchmarks
#   SKIP_LARGE=1     Skip the large (18k unknown) HODLR case
#   WORK_DIR         Scratch directory (default: a new temporary directory)
set -euo pipefail

if [ $# -ne 3 ]; then
    sed -n '2,20p' "$0"
    exit 1
fi
HERE=$(cd "$(dirname "$0")" && pwd)
ACC=$(cd "$1" && pwd)
OBL=$(cd "$2" && pwd)
mkdir -p "$3"
OUT=$(cd "$3" && pwd)
WORK=${WORK_DIR:-$(mktemp -d -t oft_blas_bench)}
mkdir -p "$WORK"
PY=${PY:-python3}
FC=${FC:-$(sed -n 's/.*-DCMAKE_Fortran_COMPILER:FILEPATH=\([^ ]*\).*/\1/p' "$ACC/config_cmake.sh")}
NT=${NTHREADS:-$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || getconf _NPROCESSORS_ONLN)}

ACC_LIB="$ACC/accelerate_shim/lib/liboftaccelerate.dylib"
OBL_LIB=$(ls "$OBL"/OpenBLAS-*/lib/libopenblas.a 2>/dev/null | head -1 || true)
[ -f "$ACC_LIB" ] || { echo "Missing $ACC_LIB (was $ACC built with --use_accelerate?)" >&2; exit 1; }
[ -n "$OBL_LIB" ] || { echo "No OpenBLAS-*/lib/libopenblas.a in $OBL" >&2; exit 1; }
for d in "$ACC" "$OBL"; do
    [ -d "$d/build_release/python/OpenFUSIONToolkit" ] || { echo "No built Python package in $d/build_release" >&2; exit 1; }
done
echo "Accelerate tree: $ACC"
echo "OpenBLAS tree:   $OBL"
echo "Results:         $OUT"
echo "Scratch:         $WORK"
echo "Threads:         $NT"

# ---- System and build information ----
{
    echo "date: $(date -u +%Y-%m-%dT%H:%MZ)"
    echo "host_model: $(sysctl -n hw.model 2>/dev/null)"
    echo "cpu: $(sysctl -n machdep.cpu.brand_string 2>/dev/null)"
    echo "arch: $(uname -m)"
    echo "ncpu: $(sysctl -n hw.ncpu 2>/dev/null)"
    echo "perf_cores: $(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || echo n/a)"
    echo "eff_cores: $(sysctl -n hw.perflevel1.physicalcpu 2>/dev/null || echo n/a)"
    echo "mem_GB: $(( $(sysctl -n hw.memsize) / 1073741824 ))"
    echo "macos: $(sw_vers -productVersion)"
    echo "sdk: $(xcrun --show-sdk-version 2>/dev/null)"
    echo "fc: $FC ($("$FC" --version | head -1))"
    echo "python: $("$PY" --version 2>&1)"
    echo "oft_rev: $(git -C "$HERE" rev-parse --short HEAD 2>/dev/null) ($(git -C "$HERE" rev-parse --abbrev-ref HEAD 2>/dev/null))"
    echo "openblas_lib: $OBL_LIB"
    echo "omp_openblas_lib: ${OMP_OPENBLAS_LIB:-none}"
    echo "nthreads: $NT"
    echo "--- Accelerate liboftpy BLAS linkage (expect liboftaccelerate) ---"
    otool -L "$ACC"/build_release/bin/liboftpy.dylib 2>/dev/null | grep -E 'liboftaccelerate|openblas|Accelerate\.framework|vecLib' || echo "(no dynamic BLAS reference found)"
} > "$OUT/sysinfo.txt"
cat "$OUT/sysinfo.txt"

# ---- BLAS/LAPACK kernel micro-benchmark ----
if [ "${SKIP_MICRO:-0}" != "1" ]; then
    echo "== Building micro-benchmark"
    "$FC" -O2 "$HERE/micro.f90" "$OBL_LIB" -o "$WORK/micro_oblas"
    "$FC" -O2 "$HERE/micro.f90" "$ACC_LIB" -o "$WORK/micro_acc"
    echo "== micro: OpenBLAS (OFT default, serial)"
    OMP_NUM_THREADS=1 "$WORK/micro_oblas" > "$OUT/micro_oblas_t1.txt"
    echo "== micro: Accelerate, 1 thread"
    VECLIB_MAXIMUM_THREADS=1 "$WORK/micro_acc" > "$OUT/micro_acc_t1.txt"
    echo "== micro: Accelerate, $NT threads"
    VECLIB_MAXIMUM_THREADS=$NT "$WORK/micro_acc" > "$OUT/micro_acc_tN.txt"
    if [ -n "${OMP_OPENBLAS_LIB:-}" ]; then
        "$FC" -O2 -fopenmp "$HERE/micro.f90" "$OMP_OPENBLAS_LIB" -o "$WORK/micro_oblas_omp"
        echo "== micro: threaded OpenBLAS, $NT threads"
        OMP_NUM_THREADS=$NT "$WORK/micro_oblas_omp" > "$OUT/micro_oblasomp_tN.txt"
    fi
fi

# ---- ThinCurr benchmarks ----
# run_tc <tag> <config_dir> <omp_threads> <veclib_threads> <ntheta> <nphi> [hodlr]
run_tc() {
    local tag=$1 tree=$2 nt=$3 vt=$4 nth=$5 nph=$6 mode=${7:-}
    echo "== ThinCurr: $tag"
    mkdir -p "$WORK/run_$tag"
    (cd "$WORK/run_$tag" && VECLIB_MAXIMUM_THREADS=$vt PYTHONPATH="$tree/build_release/python" \
        "$PY" "$HERE/bench_thincurr.py" "$WORK/mesh" "$nth" "$nph" "$nt" "$OUT/tc_$tag.json" $mode \
        > "$OUT/tc_$tag.log" 2>&1) || { echo "   FAILED (see $OUT/tc_$tag.log)"; return 0; }
    grep '^BENCH' "$OUT/tc_$tag.log" | sed 's/^/   /'
}
if [ "${SKIP_THINCURR:-0}" != "1" ]; then
    # Dense direct solvers: OFT threads = NT
    run_tc oblas_tN       "$OBL" "$NT" "$NT" 48 96
    run_tc acc_tN         "$ACC" "$NT" "$NT" 48 96
    run_tc acc_tN_veclib1 "$ACC" "$NT" 1     48 96
    # Fully serial
    run_tc oblas_t1       "$OBL" 1 1 48 96
    run_tc acc_t1         "$ACC" 1 1 48 96
    # HODLR (compressed) path
    run_tc oblas_tN_hodlr "$OBL" "$NT" "$NT" 48 96 hodlr
    run_tc acc_tN_hodlr   "$ACC" "$NT" "$NT" 48 96 hodlr
    if [ "${SKIP_LARGE:-0}" != "1" ]; then
        run_tc oblas_tN_hodlrL "$OBL" "$NT" "$NT" 96 192 hodlr
        run_tc acc_tN_hodlrL   "$ACC" "$NT" "$NT" 96 192 hodlr
    fi
fi

"$PY" "$HERE/summarize.py" "$OUT" > "$OUT/summary.md"
echo "== Summary written to $OUT/summary.md"
cat "$OUT/summary.md"
