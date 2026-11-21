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
/
"""

def write_input(nprocs, case_id, order=2):
    # Distributed runs decompose the base mesh and refine once more
    nlevels = 3
    if nprocs > 1:
        nlevels = 4
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nlevels, case_id, order, 3))

def run_ilu(nprocs, case_id, nthreads=None):
    write_input(nprocs, case_id)
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

def test_ilu_poisson():
    '''GMRES on the Poisson equation preconditioned by ILU(0) applied directly to the operator and by
    block-Jacobi with boundary overlap (ILU(0) of the full local operator) on a single process; both
    must converge (true residual checked).'''
    its = run_ilu(1, 2)
    assert its is not None
    print("Iterations: ", its)
    assert min(its) > 0

@pytest.mark.mpi
@pytest.mark.parametrize("nprocs", (2, 4))
def test_ilu_poisson_mpi(nprocs):
    '''GMRES on the Poisson equation preconditioned by block-Jacobi with boundary overlap: the ILU(0)
    sub-solver factors the full local matrix of the distributed operator, which block-Jacobi requests
    explicitly (true residual checked).'''
    its = run_ilu(nprocs, 2)
    assert its is not None
    print("Iterations: ", its)
    assert its[1] > 0

def test_ilu_direct_serial():
    '''ILU(0) applied directly to the operator on a single process (not distributed) is allowed.'''
    its = run_ilu(1, 3)
    assert its is not None
    assert its[0] > 0

@pytest.mark.mpi
def test_ilu_direct_mpi2_abort():
    '''ILU(0) applied directly to an operator distributed across 2 processes aborts with a message
    asking for block-Jacobi or Schwarz, instead of factoring local values that lack contributions
    from other processes.'''
    if not env_flag('OFT_HAVE_MPI'):
        pytest.skip("Not compiled with MPI")
    write_input(2, 3)
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
