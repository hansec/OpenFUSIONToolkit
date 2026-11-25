#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Tests for the overlapping vertex-star patch smoother (oft_patch_solvers) on a quadratic Lagrange Laplacian
import os
import sys
import pytest
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..')))
from oft_testing import run_OFT

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

&test_patch_options
 mode={1}
 rich_its={2}
 whole={3}
 single={4}
/
"""

converged_norm = 0.79845123052314
# GMRES iteration counts for additive, multiplicative and symmetric smoothers (serial)
ITS_ADD = 33
ITS_MULT = 18
ITS_SYM = 15
# Average residual reduction per Richardson iteration
RHO = {1: 0.8932600, 3: 0.6228771}
# Distributed runs: additive seam patches reproduce the serial operator exactly
ITS_MPI = {(1, 2): 33, (1, 4): 33, (3, 2): 25, (3, 4): 34}

def patch_setup(nprocs=1, mode=1, rich_its=0, whole='F', single='F'):
    '''Write input and run test_native_patch (one extra (distributed) mesh level for MPI runs)'''
    nlevels = 3
    if nprocs > 1:
        nlevels = 4
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(nlevels, mode, rich_its, whole, single))
    return run_OFT("./test_native_patch", nprocs, 60)

def validate_result(iteration_count, value, tol=1.E-8):
    '''Check iteration count (at most 5% above expected) and solution norm (or reduction factor)'''
    with open('patch.results','r') as fid:
        its_test = int(fid.readline())
        val_test = float(fid.readline())
    retval = True
    if its_test > iteration_count*1.05:
        print("FAILED: Iteration count too high! Expected = {0}, Actual = {1}".format(iteration_count, its_test))
        retval = False
    if abs(value-val_test) > tol:
        print("FAILED: Result incorrect! Expected = {0}, Actual = {1}".format(value, val_test))
        retval = False
    return retval

def test_whole():
    '''A single patch containing every row is an exact solve'''
    assert patch_setup(whole='T')
    assert validate_result(1, converged_norm)

@pytest.mark.coverage
@pytest.mark.parametrize("mode,its", ((1, ITS_ADD), (2, ITS_MULT), (3, ITS_SYM)))
def test_gmres(mode, its):
    '''GMRES preconditioned by additive, multiplicative and symmetric vertex-star smoothers'''
    assert patch_setup(mode=mode)
    assert validate_result(its, converged_norm)

def test_gmres_single():
    '''Symmetric smoother with single precision patch factors'''
    assert patch_setup(mode=3, single='T')
    assert validate_result(ITS_SYM, converged_norm)

@pytest.mark.parametrize("mode", (1, 3))
def test_richardson(mode):
    '''Stationary iteration: the average residual reduction factor over 10 iterations'''
    assert patch_setup(mode=mode, rich_its=10)
    assert validate_result(10, RHO[mode], tol=0.02)

@pytest.mark.mpi
@pytest.mark.parametrize("nprocs", (2, 4))
@pytest.mark.parametrize("mode", (1, 3))
def test_gmres_mpi(mode, nprocs):
    '''Distributed runs: patches at process boundaries are completed across processes'''
    assert patch_setup(nprocs=nprocs, mode=mode)
    assert validate_result(ITS_MPI[(mode, nprocs)], converged_norm)
