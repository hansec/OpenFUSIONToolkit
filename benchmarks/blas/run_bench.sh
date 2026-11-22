#!/usr/bin/env bash
# Benchmark an Accelerate-based OFT build against an OpenBLAS-based one, including
# the sparse direct (LU) solver packages (UMFPACK, STRUMPACK, SuperLU, Accelerate).
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
#   SKIP_SPARSE=1    Skip the sparse direct (LU) solver benchmarks
#   SKIP_LARGE=1     Skip the large cases (18k unknown HODLR, 228k unknown sparse LU)
#   SPARSE_TIMEOUT   Time limit in seconds for each sparse LU run (default: 900)
#   WORK_DIR         Scratch directory (default: a new temporary directory)
set -euo pipefail

if [ $# -ne 3 ]; then
    sed -n '2,23p' "$0"
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
# Keep an existing sysinfo.txt (eg. when adding results with SKIP_* to an existing directory)
SYSINFO="$OUT/sysinfo.txt"
[ -f "$SYSINFO" ] && SYSINFO="$OUT/sysinfo_$(date -u +%Y%m%dT%H%MZ).txt"
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
} > "$SYSINFO"
cat "$SYSINFO"

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

# ---- Sparse direct (LU) solvers ----
# Order-3 Lagrange Laplacian on a unit cube, solved with GMRES preconditioned by a full LU
# (lin_alg/test_native_bjacobi). The time of the first solve (factorization + GMRES, usually
# 1 iteration) is recorded, along with total wall time. Results go to sparse.txt:
#   <tag> <nlevels> <unknowns> <status> <wall_s> <solve_s> <iterations> <solution_norm>
# run_sparse <tag> <config_dir> <nlevels> <sol_type> [<factorization> <ordering>]
run_sparse() {
    local tag=$1 tree=$2 nl=$3 st=$4 fact=${5:-} ord=${6:-}
    local exe="$tree/build_release/tests/lin_alg/test_native_bjacobi"
    local dir="$WORK/sparse_${tag}_nl$nl" use_xml=F args="oft.in"
    echo "== Sparse LU: $tag, nlevels=$nl"
    [ -x "$exe" ] || { echo "   Missing $exe (build with --oft_build_tests=1)"; return 0; }
    mkdir -p "$dir"
    if [ -n "$fact" ]; then
        use_xml=T
        args="oft.in oft_in.xml"
        printf '<oft><lu><package>accelerate</package><accelerate><factorization>%s</factorization><ordering>%s</ordering></accelerate></lu></oft>\n' \
            "$fact" "$ord" > "$dir/oft_in.xml"
    fi
    cat > "$dir/oft.in" <<EOF
&runtime_options
 ppn=1
 debug=1
 test_run=T
/
&mesh_options
 meshname='cube'
 cad_type=92
 nlevels=$nl
 nbase=$nl
 grid_order=1
/
&cube_options
 mesh_type=1
/
&test_bj_options
 nlocal=1
 sol_type=$st
 use_xml=$use_xml
/
EOF
    local log="$OUT/sparse_${tag}_nl$nl.log" status=ok rc=0 t0 t1
    rm -f "$dir/bjacobi.results"
    t0=$("$PY" -c 'import time; print(time.time())')
    (cd "$dir" && OMP_NUM_THREADS=$NT VECLIB_MAXIMUM_THREADS=$NT ${TIMEOUT_CMD:+$TIMEOUT_CMD "${SPARSE_TIMEOUT:-900}"} "$exe" $args > "$log" 2>&1) || rc=$?
    t1=$("$PY" -c 'import time; print(time.time())')
    if grep -q 'SKIP TEST' "$log"; then
        status=skipped
    elif [ $rc -eq 124 ]; then
        status=timeout
    elif [ $rc -ne 0 ] || [ ! -f "$dir/bjacobi.results" ]; then
        status=failed
    fi
    # Order-3 Lagrange unknowns = points + 2*edges + faces on the finest mesh level
    local ndofs
    ndofs=$(awk '/# of points  /{p=$NF} /# of edges  /{e=$NF} /# of faces  /{f=$NF} END{print p+2*e+f}' "$log")
    local res="- - -"
    [ "$status" = ok ] && res=$(awk 'NR==1{c=$1} NR==2{u=$1} NR==3{t=$1} END{print t, c, u}' "$dir/bjacobi.results")
    local wall
    wall=$("$PY" -c "print('{0:.2f}'.format($t1-$t0))")
    echo "$tag $nl $ndofs $status $wall $res" >> "$OUT/sparse.txt"
    echo "   $status: wall=${wall} s, solve/iterations/norm = $res"
}
if [ "${SKIP_SPARSE:-0}" != "1" ]; then
    # macOS has no timeout command; use GNU coreutils (Homebrew) if available
    TIMEOUT_CMD=$(command -v timeout || command -v gtimeout || true)
    [ -n "$TIMEOUT_CMD" ] || echo "Warning: no timeout/gtimeout command found, sparse LU runs are not time limited"
    rm -f "$OUT/sparse.txt"
    SPARSE_LEVELS="4"
    [ "${SKIP_LARGE:-0}" != "1" ] && SPARSE_LEVELS="4 5"
    for nl in $SPARSE_LEVELS; do
        run_sparse acc_accelerate           "$ACC" "$nl" 7
        run_sparse acc_accelerate_metis     "$ACC" "$nl" 7 lu metis
        run_sparse acc_accelerate_spp_metis "$ACC" "$nl" 7 lu_spp metis
        run_sparse acc_strumpack            "$ACC" "$nl" 6
        run_sparse oblas_strumpack          "$OBL" "$nl" 6
        run_sparse acc_umfpack              "$ACC" "$nl" 4
        run_sparse oblas_umfpack            "$OBL" "$nl" 4
        # SuperLU takes more than 30 minutes at nlevels=5
        [ "$nl" = 4 ] && run_sparse acc_superlu "$ACC" "$nl" 1
    done
fi

"$PY" "$HERE/summarize.py" "$OUT" > "$OUT/summary.md"
echo "== Summary written to $OUT/summary.md"
cat "$OUT/summary.md"
