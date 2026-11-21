import os
import sys
import subprocess
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT, env_flag

# Basic template for input file
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels={0}
 nbase={3}
 grid_order=1
/

&cube_options
 mesh_type=1
/

&test_ilu_options
 case_id={1}
 order={2}
 ilu_pkg='{4}'
/
"""

# ILU packages run on the operator tests: native ILU(0) and SuperLU ILUTP (skipped by the
# test program when OFT is built without SuperLU)
ILU_packages = ('native', 'super')

def write_input(nprocs, case_id, order=2, package='native'):
    # Distributed runs decompose the base mesh and refine once more
    nlevels = 3
    if nprocs > 1:
        nlevels = 4
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nlevels, case_id, order, 3, package))

def run_ilu(nprocs, case_id, nthreads=None, package='native'):
    write_input(nprocs, case_id, package=package)
    # Thread count for this run only (other tests use the caller's OMP_NUM_THREADS)
    omp_save = os.environ.get('OMP_NUM_THREADS')
    if nthreads is not None:
        os.environ['OMP_NUM_THREADS'] = str(nthreads)
    try:
        if not run_OFT("./test_native_ilu", nprocs, 60):
            return None
    finally:
        if omp_save is None:
            os.environ.pop('OMP_NUM_THREADS', None)
        else:
            os.environ['OMP_NUM_THREADS'] = omp_save
    with open('ilu.results', 'r') as fid:
        nfail = int(fid.readline())
        its = [int(val) for val in fid.readline().split()]
    if nfail != 0:
        print("FAILED: {0} ILU checks failed".format(nfail))
        return None
    return its

@pytest.mark.coverage
def test_ilu_local():
    '''ILU(0) on matrices built directly: one application solves banded systems with no dropped fill
    to round-off (ILU(0) is then the exact LU factorization), refactoring after value changes and after
    switching to matrices with a different pattern or size (detected automatically, or requested with
    update(new_pattern=.TRUE.)), sign-preserving replacement of small and zero pivots (checked against
    the solution with the replaced pivot), the diagonal shift (exact solve of the shifted matrix) and
    deep copies of a factored solver.'''
    assert run_ilu(1, 1) is not None

@pytest.mark.mpi
def test_ilu_local_mpi2():
    '''Local ILU(0) checks run independently on each of 2 processes.'''
    assert run_ilu(2, 1) is not None

@pytest.mark.parametrize("package", ILU_packages)
def test_ilu_poisson(package):
    '''GMRES on the Poisson equation preconditioned by ILU (native ILU(0) or SuperLU ILUTP) applied
    directly to the operator and by block-Jacobi with boundary overlap (ILU of the full local operator)
    on a single process; both must converge (true residual checked).'''
    its = run_ilu(1, 2, package=package)
    assert its is not None
    print("Iterations: ", its)
    assert min(its) > 0

@pytest.mark.mpi
@pytest.mark.parametrize("package", ILU_packages)
@pytest.mark.parametrize("nprocs", (2, 4))
def test_ilu_poisson_mpi(nprocs, package):
    '''GMRES on the Poisson equation preconditioned by block-Jacobi with boundary overlap: the ILU
    sub-solver factors the full local matrix of the distributed operator, which block-Jacobi requests
    explicitly (true residual checked).'''
    its = run_ilu(nprocs, 2, package=package)
    assert its is not None
    print("Iterations: ", its)
    assert its[1] > 0

@pytest.mark.parametrize("package", ILU_packages)
def test_ilu_direct_serial(package):
    '''ILU applied directly to the operator on a single process (not distributed) is allowed.'''
    its = run_ilu(1, 3, package=package)
    assert its is not None
    assert its[0] > 0

@pytest.mark.mpi
@pytest.mark.parametrize("package", ILU_packages)
def test_ilu_direct_mpi2_abort(package):
    '''ILU applied directly to an operator distributed across 2 processes aborts with a message
    asking for block-Jacobi or Schwarz, instead of factoring local values that lack contributions
    from other processes.'''
    if not env_flag('OFT_HAVE_MPI'):
        pytest.skip("Not compiled with MPI")
    write_input(2, 3, package=package)
    command = "mpirun --map-by :OVERSUBSCRIBE -np 2 ./test_native_ilu"
    timeout = 60
    if env_flag('OFT_DEBUG_TEST'):
        timeout *= 4
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            timeout=timeout)
    std_out = result.stdout.decode()
    std_err = result.stderr.decode()
    print(std_out)
    print(std_err)
    if std_out.find('SKIP TEST') > -1:
        pytest.skip("Configuration not supported")
    assert result.returncode != 0
    assert std_err.find('distributed across processes') > -1

@pytest.mark.coverage
@pytest.mark.parametrize("nthreads", (1, 2, 4))
def test_ilu_threaded(nthreads):
    '''Level-scheduled (OpenMP threaded) ILU(0) solve on a non-symmetric variant of the quadratic
    Lagrange Laplacian, run with 1, 2 and 4 threads: the schedules group rows by dependency level
    (checked against the matrix structure), the scheduled solve is used exactly when more than one
    thread is available, and its solutions are bitwise identical to the serial solve for real and
    complex values, with and without replaced pivots and a diagonal shift, and 2 to 4 threads.
    Inside a parallel region the solve runs serially, patterns without a diagonal entry in every
    row are never threaded, and the threaded ILU(0) solver matches the serial one bitwise.'''
    its = run_ilu(1, 5, nthreads)
    assert its is not None

@pytest.mark.mpi
def test_ilu_schwarz_mpi2():
    '''GMRES on the Poisson equation preconditioned by Schwarz (one subdomain per process, overlap 1)
    with ILU(0) on 2 processes and 2 threads, with and without the threaded solve: iteration
    counts and solutions must be identical. Also reports wall times (end-to-end timing harness;
    larger cases are run by hand by raising nlevels and order in oft.in).'''
    its = run_ilu(2, 6, 2)
    assert its is not None
    print("Iterations: ", its)
    assert its[0] > 0 and its[0] == its[1]

def test_ilu_timing():
    '''Timing harness for the ILU(0) factorization and solve (Laplacian, quadratic Lagrange
    elements, one process) on a small mesh; checks that it runs and reports the timings, including
    the factorization and the serial and level-scheduled solves for 1, 2, 4, ... threads (up to
    OMP_NUM_THREADS).
    Larger cases are run by hand by raising nlevels and nbase (equal, for one process) and order
    in oft.in.'''
    write_input(1, 4)
    assert run_OFT("./test_native_ilu", 1, 60)
    with open('ilu_timing.results', 'r') as fid:
        nr, nnz = [int(val) for val in fid.readline().split()[:2]]
        tfactor, tsolve = [float(val) for val in fid.readline().split()[1:]]
    print("nr = {0}, nnz = {1}: factor {2:.3e} s, solve {3:.3e} s".format(nr, nnz, tfactor, tsolve))
    assert nr > 0 and nnz > 0
    assert tsolve > 0.0

@pytest.mark.coverage
def test_ilu_superlu():
    '''SuperLU ILUTP ("super" package) on matrices built directly: without dropping it is a complete
    LU (residual at round-off where ILU(0) is not), also for several right hand sides and matrices
    without a pattern, and does not modify the matrix values; refactoring after value changes and
    after changes of the numerical options reuses the column permutation and elimination tree and
    matches a new factorization bitwise; pattern and size changes, update(new_pattern=.TRUE.) and
    column permutation changes redo the symbolic setup; zero pivots are replaced and counted; copies
    of a factored solver factor anew and are deleted independently (no shared SuperLU factors).
    Skipped without SuperLU.'''
    assert run_ilu(1, 7) is not None

def test_ilu_compare():
    '''Comparison harness for native ILU(0) and SuperLU ILUTP at drop tolerances 1E-2, 1E-3 and 1E-4
    (the last also with a larger fill bound) and MILU variants on the quadratic Lagrange Laplacian
    (one process): all preconditioned GMRES solves converge, ILUTP with the default drop tolerance
    (1E-3) needs fewer iterations than ILU(0), a larger drop tolerance does not need fewer, and a
    smaller one with a larger fill bound does not need more. Reports factorization and solve times
    and the fill of the factors; larger cases are run by hand by raising nlevels, nbase and order in
    oft.in. Skipped without SuperLU.'''
    write_input(1, 8)
    assert run_OFT("./test_native_ilu", 1, 120)
    with open('ilu.results', 'r') as fid:
        assert int(fid.readline()) == 0
    with open('ilu_compare.results', 'r') as fid:
        nr, nnz = [int(val) for val in fid.readline().split()]
        results = [fid.readline().split() for _ in range(8)]
    its = [int(row[5]) for row in results]
    print("nr = {0}, nnz = {1}, iterations: {2}".format(nr, nnz, its))
    assert min(its) > 0
    # Configurations 1-5: native ILU(0), then ILUTP (SILU) with drop_tol 1E-2, 1E-3, 1E-4 and 1E-4 with
    # fill_factor 30 (with the default bound of 10, 1E-4 can reach the bound on larger meshes, where
    # entries are then dropped by size alone and the preconditioner weakens)
    assert its[2] < its[0]
    assert its[4] <= its[2] <= its[1]
