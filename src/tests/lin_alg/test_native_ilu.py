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
 nbase=3
 grid_order=1
/

&cube_options
 mesh_type=1
/

&test_ilu_options
 case_id={1}
/
"""

def write_input(nprocs, case_id):
    # Distributed runs decompose the base mesh and refine once more
    nlevels = 3
    if nprocs > 1:
        nlevels = 4
    os.chdir(test_dir)
    with open('oft.in', 'w+') as fid:
        fid.write(oft_in_template.format(nlevels, case_id))

def run_ilu(nprocs, case_id):
    write_input(nprocs, case_id)
    if not run_OFT("./test_native_ilu", nprocs, 60):
        return None
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
