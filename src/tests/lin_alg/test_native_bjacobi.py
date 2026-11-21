import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT

LU_solvers = (1,2,3,4,5)
# ILU packages: 0 native ILU(0), 1 SuperLU ILUTP (default options), 5 MKL ILU(0)
ILU_solvers = (0,1,5)

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

&test_bj_options
 nlocal={1}
 sol_type={2}
 use_ilu={3}
/
"""

# Common setup function and process handling
def native_bjacobi_setup(nprocs=1, nlocal=1, sol_type=1, use_ilu='F'):
    nlevels = 3
    if nprocs > 1:
        nlevels = 4
    #
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nlevels, nlocal, sol_type, use_ilu))
    return run_OFT("./test_native_bjacobi", nprocs, 20)

def validate_result(iteration_count,converged_error):
    #
    fid = open('bjacobi.results','r')
    its_test = int(fid.readline())
    if iteration_count is not None:
        if its_test>iteration_count*1.05:
            print("FAILED: Iteration count incorrect!")
            print("  Expected = {0}".format(iteration_count))
            print("  Actual =   {0}".format(its_test))
            return False
    error_test = float(fid.readline())
    if abs(converged_error-error_test) > 1.E-6:
        print("FAILED: Residual error incorrect!")
        print("  Expected = {0}".format(converged_error))
        print("  Actual =   {0}".format(error_test))
        return False
    return True

def ilu_count(solver_type, ilu0_count, ilut_count):
    # Expected iterations for an ILU package: the ILU(0) packages share counts
    if solver_type == 1:
        return ilut_count
    return ilu0_count

#============================================================================
#
def check_solver(nprocs,nlocal,sol_type,iteration_count,converged_error,use_ilu='F'):
    assert native_bjacobi_setup(nprocs,nlocal,sol_type,use_ilu=use_ilu)
    assert validate_result(iteration_count,converged_error)

# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("solver_type", LU_solvers)
def test_LU_base(solver_type):
    iteration_count = 1
    converged_error = 1.4688646289584524
    check_solver(1, 1, solver_type, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.parametrize("solver_type", LU_solvers)
def test_LU_base_nlocal2(solver_type):
    iteration_count = 41
    converged_error = 1.4688646289584524
    check_solver(1, 2, solver_type, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("solver_type", LU_solvers)
def test_LU_base_nlocal2_part(solver_type):
    iteration_count = 46
    converged_error = 1.4688646289584524
    check_solver(1, -2, solver_type, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("solver_type", LU_solvers)
def test_LU_mpi2(solver_type):
    iteration_count = 50
    converged_error = 1.4688646289584524
    check_solver(2, 1, solver_type, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.parametrize("solver_type", LU_solvers)
def test_LU_mpi2_nlocal2(solver_type):
    iteration_count = 56
    converged_error = 1.4688646289584524
    check_solver(2, 2, solver_type, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("solver_type", LU_solvers)
def test_LU_mpi2_nlocal2_part(solver_type):
    iteration_count = 57
    converged_error = 1.4688646289584524
    check_solver(2, -2, solver_type, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.parametrize("solver_type", LU_solvers)
def test_LU_mpi4(solver_type):
    iteration_count = 60
    converged_error = 1.4688646289584524
    check_solver(4, 1, solver_type, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.parametrize("solver_type", LU_solvers)
def test_LU_mpi4_nlocal2(solver_type):
    iteration_count = 65
    converged_error = 1.4688646289584524
    check_solver(4, 2, solver_type, iteration_count, converged_error)
# Test runner for base test case
@pytest.mark.parametrize("solver_type", LU_solvers)
def test_LU_mpi4_nlocal2_part(solver_type):
    iteration_count = 65
    converged_error = 1.4688646289584524
    check_solver(4, -2, solver_type, iteration_count, converged_error)


# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("solver_type", ILU_solvers)
def test_ILU_base(solver_type):
    iteration_count = ilu_count(solver_type, 40, 8)
    converged_error = 1.4688646289584524
    check_solver(1, 1, solver_type, iteration_count, converged_error, use_ilu='T')
# Test runner for base test case
@pytest.mark.parametrize("solver_type", ILU_solvers)
def test_ILU_base_nlocal2(solver_type):
    iteration_count = ilu_count(solver_type, 59, 40)
    converged_error = 1.4688646289584524
    check_solver(1, 2, solver_type, iteration_count, converged_error, use_ilu='T')
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("solver_type", ILU_solvers)
def test_ILU_base_nlocal2_part(solver_type):
    iteration_count = ilu_count(solver_type, 61, 46)
    converged_error = 1.4688646289584524
    check_solver(1, -2, solver_type, iteration_count, converged_error, use_ilu='T')
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("solver_type", ILU_solvers)
def test_ILU_mpi2(solver_type):
    iteration_count = ilu_count(solver_type, 64, 49)
    converged_error = 1.4688646289584524
    check_solver(2, 1, solver_type, iteration_count, converged_error, use_ilu='T')
# Test runner for base test case
@pytest.mark.parametrize("solver_type", ILU_solvers)
def test_ILU_mpi2_nlocal2(solver_type):
    iteration_count = ilu_count(solver_type, 69, 55)
    converged_error = 1.4688646289584524
    check_solver(2, 2, solver_type, iteration_count, converged_error, use_ilu='T')
# Test runner for base test case
@pytest.mark.coverage
@pytest.mark.parametrize("solver_type", ILU_solvers)
def test_ILU_mpi2_nlocal2_part(solver_type):
    iteration_count = ilu_count(solver_type, 67, 58)
    converged_error = 1.4688646289584524
    check_solver(2, -2, solver_type, iteration_count, converged_error, use_ilu='T')
# Test runner for base test case
@pytest.mark.parametrize("solver_type", ILU_solvers)
def test_ILU_mpi4(solver_type):
    iteration_count = ilu_count(solver_type, 68, 57)
    converged_error = 1.4688646289584524
    check_solver(4, 1, solver_type, iteration_count, converged_error, use_ilu='T')
# Test runner for base test case
@pytest.mark.parametrize("solver_type", ILU_solvers)
def test_ILU_mpi4_nlocal2(solver_type):
    iteration_count = ilu_count(solver_type, 73, 58)
    converged_error = 1.4688646289584524
    check_solver(4, 2, solver_type, iteration_count, converged_error, use_ilu='T')
# Test runner for base test case
@pytest.mark.parametrize("solver_type", ILU_solvers)
def test_ILU_mpi4_nlocal2_part(solver_type):
    iteration_count = ilu_count(solver_type, 70, 64)
    converged_error = 1.4688646289584524
    check_solver(4, -2, solver_type, iteration_count, converged_error, use_ilu='T')
